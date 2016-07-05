// Bistro: main.C
// Bret Larget and Claudia Solis-Lemus
// Copyright 2016
// May 20, 2016

#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <string>
#include <random>

#include "parameter.h"
#include "sequence.h"
#include "alias.h"
#include "tree.h"
#include "model.h"
#include "ccdprobs.h"

#include "Eigen/Core"
#include "Eigen/Eigenvalues"

using namespace std;
using namespace Eigen;

int main(int argc, char* argv[])
{
  // Read command line and process parameters
  cerr << "Processing command line ...";
  Parameter parameters;
  parameters.processCommandLine(argc,argv);
  cerr << " done." << endl;

  // Read in sequences from FASTA file
  cerr << "Reading alignment from FASTA file ...";
  Alignment alignment(parameters.getSequenceFileName());
  cerr << " done." << endl;
  cerr << endl << "Alignment summary:" << endl;
  alignment.summarize(cout);
  cout << endl;

  // Find Jukes-Cantor pairwise distances
  cerr << "Finding initial Jukes-Cantor pairwise distances ...";
  MatrixXd jcDistanceMatrix(alignment.getNumTaxa(),alignment.getNumTaxa());
  alignment.calculateJCDistances(jcDistanceMatrix);
  cerr << " done." << endl;

  cerr << endl << jcDistanceMatrix << endl << endl;
  
  // Find Initial Neighbor-joining tree
  cerr << "Finding initial neighbor-joining tree ...";
  MatrixXd jcDistanceMatrixCopy(alignment.getNumTaxa(),alignment.getNumTaxa());
  jcDistanceMatrixCopy = jcDistanceMatrix;
  Tree jctree(jcDistanceMatrixCopy);
  jctree.reroot(1);
  jctree.sortCanonical();
  cerr << " done." << endl;
  cerr << endl << "Tree topology:" << endl;
  cerr << jctree.makeTopologyNumbers() << endl << endl;

  // Initialize random number generator
  cerr << "Initializing random number generator ...";
  if ( parameters.getSeed() == 0 )
  {
    random_device rd;
    parameters.setSeed(rd());
  }
  mt19937_64 rng(parameters.getSeed());
  cerr << " done." << endl;
  cerr << "Seed = " << parameters.getSeed() << endl;

  // Run MCMC on tree to estimate Q matrix parameters
  //   initial Q matrix
  
  vector<double> p_init(4,0.25);
  vector<double> s_init(6,0.1);
  s_init[1] = 0.3;
  s_init[4] = 0.3;

  QMatrix q_init(p_init,s_init);

  // burnin
  q_init.mcmc(alignment,jctree,1000,alignment.getNumSites(),rng);

  // mcmc to get final Q
  q_init.mcmc(alignment,jctree,1000,alignment.getNumSites(),rng);

//  QMatrix model(parameters.getStationaryP(),parameters.getSymmetricQP());
  QMatrix model(q_init.getStationaryP(),q_init.getSymmetricQP());

  // Recalculate pairwise distances using estimated Q matrix (TODO: add site rate heterogeneity)
  cerr << "Finding initial GTR pairwise distances ...";
  MatrixXd gtrDistanceMatrix(alignment.getNumTaxa(),alignment.getNumTaxa());
  alignment.calculateGTRDistances(model,jcDistanceMatrix,gtrDistanceMatrix);
  cerr << " done." << endl;

  cerr << endl << gtrDistanceMatrix << endl << endl;
  
  // Do bootstrap with new pairwise distances
  cerr << "Beginning " << parameters.getNumBootstrap() << " bootstrap replicates ..." << endl;
  map<string,int> topologyToCountMap;
  MatrixXd bootDistanceMatrix(alignment.getNumTaxa(),alignment.getNumTaxa());
  vector<int> weights(alignment.getNumSites());
  cerr << '|';
  for ( int b=0; b<parameters.getNumBootstrap(); ++b )
  {
    if ( (b+1) % (parameters.getNumBootstrap() / 100) == 0 )
      cerr << '*';
    if ( (b+1) % (parameters.getNumBootstrap() / 10) == 0 )
      cerr << '|';
    alignment.setBootstrapWeights(weights,rng);
    alignment.calculateGTRDistancesUsingWeights(weights,model,gtrDistanceMatrix,bootDistanceMatrix);
//    alignment.calculateJCDistancesUsingWeights(weights,bootDistanceMatrix);
    Tree bootTree(bootDistanceMatrix);
    bootTree.reroot(1);
    bootTree.sortCanonical();
    topologyToCountMap[ bootTree.makeTopologyNumbers() ]++;
  }
  cerr << endl << "done." << endl;

  cerr << endl << "Topology counts:" << endl;
  for ( map<string,int>::iterator m=topologyToCountMap.begin(); m != topologyToCountMap.end(); ++m )
    cerr << (*m).first << " " << setw(5) << (*m).second << endl;
  cerr << endl;

  vector<int> taxaNumbers;
  vector<string> taxaNames;
  alignment.getTaxaNumbersAndNames(taxaNumbers,taxaNames);
  CCDProbs ccd(topologyToCountMap,taxaNumbers,taxaNames);

  ofstream f(parameters.getOutFileRoot().c_str());

  f << "tree logl logTop logProp logPrior logWt" << endl;

  int numRandom = parameters.getNumRandom();

  vector<double> logwt(numRandom,0);
  double maxLogWeight;

  if ( parameters.getNumBootstrap() > 0 )
  {
    cerr << "Generating " << numRandom << " random trees:" << endl << '|';
    for ( int k=0; k<numRandom; ++k )
    {
      if ( (k+1) % (numRandom / 100) == 0 )
	cerr << '*';
      if ( (k+1) % (numRandom / 10) == 0 )
	cerr << '|';
      double logTopologyProbability=0;
      string treeString = ccd.randomTree(rng,logTopologyProbability);
      Tree tree(treeString);
      tree.relabel(alignment);
      tree.unroot();
      MatrixXd jcDistanceMatrixCopy(alignment.getNumTaxa(),alignment.getNumTaxa());
      jcDistanceMatrixCopy = jcDistanceMatrix;
      tree.setNJDistances(jcDistanceMatrixCopy,rng);
      tree.sortCanonical();
      tree.randomize(rng);
      double logProposalDensity = 0;
      for ( int i=0; i<parameters.getNumMLE(); ++i )
      {
	tree.randomEdges(alignment,model,rng,logProposalDensity,true);
      }
      tree.randomEdges(alignment,model,rng,logProposalDensity,false);
      double logBranchLengthPriorDensity = tree.logPriorExp(0.1);
      double logLik = tree.calculate(alignment, model);
      double logWeight = logTopologyProbability + logBranchLengthPriorDensity + logLik - logProposalDensity;
      tree.reroot(1);
      tree.sortCanonical();
      f << tree.makeTopologyNumbers() << " " << logLik << " " << logTopologyProbability << " " << logProposalDensity << " " << logBranchLengthPriorDensity << " " << logWeight << endl;
      logwt[k] = logWeight;
      if ( k==0 || logWeight > maxLogWeight )
	maxLogWeight = logWeight;
    }
    cerr << endl << "done." << endl;
    vector<double> wt(numRandom,0);
    double sum=0;
    for ( int k=0; k<numRandom; ++k )
    {
      wt[k] = exp(logwt[k] - maxLogWeight);
      sum += wt[k];
    }
    double essInverse=0;
    for ( int k=0; k<numRandom; ++k )
    {
      wt[k] /= sum;
      essInverse += wt[k]*wt[k];
    }
    cout << "ESS = " << fixed << setprecision(2) << 1.0/essInverse << ", or "
	 << setprecision(2) << 100.0 / essInverse / numRandom << " percent." << endl;
  }
  return 0;
}
