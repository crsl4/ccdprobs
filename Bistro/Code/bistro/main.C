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
#include "random.h"

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
  cout << endl << "Alignment summary:" << endl;
  alignment.summarize(cout);
  cout << endl;

  // Find Jukes-Cantor pairwise distances
  cerr << alignment.getNumTaxa() << endl;
  cerr << "Finding initial Jukes-Cantor pairwise distances ...";
  MatrixXd jcDistanceMatrix(alignment.getNumTaxa(),alignment.getNumTaxa());
  alignment.calculateJCDistances(jcDistanceMatrix);
  cerr << " done." << endl;

  cout << "Jukes-Cantor Distance Matrix:" << endl;
  cout << endl << jcDistanceMatrix << endl << endl;

  // Find Initial Neighbor-joining tree
  cerr << "Finding initial neighbor-joining tree ...";
  MatrixXd jcDistanceMatrixCopy(alignment.getNumTaxa(),alignment.getNumTaxa());
  jcDistanceMatrixCopy = jcDistanceMatrix;
  Tree jctree(jcDistanceMatrixCopy);
  jctree.reroot(1); //warning: if 1 changed, need to change makeBinary if called after
  jctree.sortCanonical();
  cerr << " done." << endl;
  cout << endl << "Tree topology:" << endl;
  cout << jctree.makeTopologyNumbers() << endl << endl;
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
  cout << "Seed = " << parameters.getSeed() << endl;

  cerr << "Running MCMC to estimate Q matrix ..." << endl;
  // Run MCMC on tree to estimate Q matrix parameters
  //   initial Q matrix

  vector<double> p_init(4,0.25);
  vector<double> s_init(6,0.1);
  s_init[1] = 0.3;
  s_init[4] = 0.3;

  QMatrix q_init(p_init,s_init);

  // burnin
  cerr << "burn-in:" << endl;
  q_init.mcmc(alignment,jctree,1000,alignment.getNumSites(),rng);
  cerr << endl << " done." << endl;

  // mcmc to get final Q
  cerr << "sampling:" << endl;
  q_init.mcmc(alignment,jctree,10000,alignment.getNumSites(),rng);
  cerr << endl << " done." << endl;

//  QMatrix model(parameters.getStationaryP(),parameters.getSymmetricQP());
  QMatrix model(q_init.getStationaryP(),q_init.getSymmetricQP());

  cerr << endl << " done." << endl;

  // Recalculate pairwise distances using estimated Q matrix (TODO: add site rate heterogeneity)
  cerr << "Finding initial GTR pairwise distances ...";
  MatrixXd gtrDistanceMatrix(alignment.getNumTaxa(),alignment.getNumTaxa());
  alignment.calculateGTRDistances(model,jcDistanceMatrix,gtrDistanceMatrix);
  cerr << " done." << endl;

  cout << "GTR Distance Matrix:" << endl;
  cout << endl << gtrDistanceMatrix << endl << endl;

  // Do bootstrap with new pairwise distances
  cerr << "Beginning " << parameters.getNumBootstrap() << " bootstrap replicates ..." << endl;
  map<string,int> topologyToCountMap;
  map<string,int> topologyToParsimonyScoreMap;
  int minimumParsimonyScore = -1;
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
    Tree bootTree(bootDistanceMatrix);
    bootTree.reroot(1); //warning: if 1 changes, need to change makeBinary if called after
    bootTree.makeBinary();
    bootTree.sortCanonical();

    string top = bootTree.makeTopologyNumbers();
    // add score to map if it is not already there
    // update minimum parsimony score if new minimum found
    if ( topologyToParsimonyScoreMap.find(top) == topologyToParsimonyScoreMap.end() )
    {
      int score = bootTree.parsimonyScore(alignment);
      topologyToParsimonyScoreMap[ top ] = score;
      if ( score < minimumParsimonyScore || b==0 )
	minimumParsimonyScore = score;
    }
    topologyToCountMap[ top ]++;
  }
  cerr << endl << "done." << endl;

  map<string,double> topologyToWeightMap; // not normalized; normalization taken care of during alias creation
  for ( map<string,int>::iterator m=topologyToCountMap.begin(); m != topologyToCountMap.end(); ++m )
  {
    topologyToWeightMap[ (*m).first ] =
      (*m).second * exp( parameters.getParsimonyScale() *
			 (minimumParsimonyScore - topologyToParsimonyScoreMap[ (*m).first ]) );
  }

  cout << endl << "Topology counts:" << endl;
  {
    map<string,double>::iterator wm=topologyToWeightMap.begin();
    for ( map<string,int>::iterator cm=topologyToCountMap.begin(); cm != topologyToCountMap.end(); ++cm )
    {
      cout << (*cm).first << " " << setw(5) << (*cm).second << " " << setw(10) << setprecision(4) << fixed << (*wm).second;
      cout << " " << setw(5) << topologyToParsimonyScoreMap[ (*cm).first ] << " " << setw(4) << minimumParsimonyScore - topologyToParsimonyScoreMap[ (*cm).first ] << endl;
      ++wm;
    }
    cout << endl;
  }

  vector<int> taxaNumbers;
  vector<string> taxaNames;
  alignment.getTaxaNumbersAndNames(taxaNumbers,taxaNames);

  // here do parsimony score (keep the option to do without parsimony)
  CCDProbs<int> ccd(topologyToCountMap,taxaNumbers,taxaNames);
  CCDProbs<double> ccdParsimony(topologyToWeightMap,taxaNumbers,taxaNames);
  // write map out to temp files to check
  ofstream smap("temp.smap");
  ccd.writeCladeCount(smap);
  smap.close();
  ofstream tmap("temp.tmap");
  ccd.writePairCount(tmap);
  tmap.close();
  smap.open("pars.smap");
  ccdParsimony.writeCladeCount(smap);
  smap.close();
  tmap.open("pars.tmap");
  ccdParsimony.writePairCount(tmap);
  tmap.close();

  //   CCDProbs<double> ccd(topologyToCountMap,taxaNumbers,taxaNames);
  //CCDProbs ccd(topologyToCountMap,taxaNumbers,taxaNames);

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
      if ( numRandom > 99 && (k+1) % (numRandom / 100) == 0 )
	cerr << '*';
      if ( numRandom > 9 && (k+1) % (numRandom / 10) == 0 )
	cerr << '|';
      double logTopologyProbability=0;
      //string treeString = ccd.randomTree(rng,logTopologyProbability);
      string treeString = ccdParsimony.randomTree(rng,logTopologyProbability);
      Tree tree(treeString);
      tree.relabel(alignment);
      tree.unroot();
      MatrixXd gtrDistanceMatrixCopy(alignment.getNumTaxa(),alignment.getNumTaxa());
      gtrDistanceMatrixCopy = gtrDistanceMatrix;
      tree.setNJDistances(gtrDistanceMatrixCopy,rng);
      tree.randomize(rng);
      tree.print(cout);
      cout << tree.makeTreeNumbers() << endl;
      double logProposalDensity = 0;
      for ( int i=0; i<parameters.getNumMLE(); ++i )
      {
	tree.randomEdges(alignment,model,rng,logProposalDensity,true);
	cout << tree.makeTreeNumbers() << endl;
      }
      if( parameters.getIndependent() )
	{
	  cout << "Branch lengths sampled independently" << endl;
	  tree.randomEdges(alignment,model,rng,logProposalDensity,false);
	}
      else
	{
	  cout << "Branch lengths sampled jointly in 2D" << endl;
	  tree.generateBranchLengths(alignment,model,rng, logProposalDensity, false);
	}
      cout << tree.makeTreeNumbers() << endl;
      double logBranchLengthPriorDensity = tree.logPriorExp(0.1);
      double logLik = tree.calculate(alignment, model);
      double logWeight = logTopologyProbability + logBranchLengthPriorDensity + logLik - logProposalDensity;
      tree.reroot(1); //warning: if 1 changes, need to change makeBinary if called after
      tree.sortCanonical();
//      int score = tree.parsimonyScore(alignment);
//      f << tree.makeTopologyNumbers() << " " << logLik << " " << logTopologyProbability << " " << logProposalDensity << " " << logBranchLengthPriorDensity << " " << logWeight << " " << score << endl;
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
    cerr << "ESS = " << fixed << setprecision(2) << 1.0/essInverse << ", or "
	 << setprecision(2) << 100.0 / essInverse / numRandom << " percent." << endl;
  }
  return 0;
}
