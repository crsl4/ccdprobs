// Bistro: main.C
// Bret Larget and Claudia Solis-Lemus
// Copyright 2016
// May 20, 2016
// will only run ccdprobs on bootstrap sample (with parsimony and lik weights)

#define MCMC_Q_BURN 100
#define MCMC_Q_SAMPLE 1000


#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <string>
#include <random>
#include <thread>
#include <chrono>

#include "parameter.h"
#include "sequence.h"
#include "alias.h"
#include "tree.h"
#include "model.h"
#include "ccdprobs.h"
#include "random.h"

#include "Eigen/Core"
#include "Eigen/Eigenvalues"
#include "Eigen/SVD"

using namespace std;
using namespace std::chrono;
using namespace Eigen;

bool comparePairStringDouble(const pair<string, double>  &p1, const pair<string, double> &p2)
{
  return p1.second < p2.second;
}

VectorXd dirichletProposalDensity(VectorXd x,double scale,double& logProposalDensity,mt19937_64& rng)
{
  VectorXd alpha(x.size());
  VectorXd y(x.size());
  double sum = 0;
  double sumalpha = 0;
  for ( int i=0; i<x.size(); ++i )
  {
    alpha(i) = scale*x(i) + 1;
    gamma_distribution<double> rgamma(alpha(i),1.0);
    y(i) = rgamma(rng);
    sum += y(i);
    sumalpha += alpha(i);
  }
  for ( int i=0; i<x.size(); ++i )
  {
    y(i) /= sum;
  }
  for ( int i=0; i<x.size(); ++i )
  {
    logProposalDensity += ( - lgamma(alpha(i)) + scale*x(i)*log(y(i)) ) ;
  }
  logProposalDensity += lgamma(sumalpha);
  return y;
}


int main(int argc, char* argv[])
{
  milliseconds ms0 = duration_cast< milliseconds >( system_clock::now().time_since_epoch() );
  // Read command line and process parameters
  cerr << "Processing command line ...";
  Parameter parameters;
  parameters.processCommandLine(argc,argv);
  cerr << " done." << endl;

  milliseconds ms1 = duration_cast< milliseconds >( system_clock::now().time_since_epoch() );
  // Read in sequences from FASTA file
  cerr << "Reading alignment from FASTA file ...";
  Alignment alignment(parameters.getSequenceFileName());
  cerr << " done." << endl;
  cout << endl << "Alignment summary:" << endl;
  alignment.summarize(cout);
  cout << endl;

  milliseconds ms2 = duration_cast< milliseconds >( system_clock::now().time_since_epoch() );

  // Find Jukes-Cantor pairwise distances
//  cerr << alignment.getNumTaxa() << endl;
  cerr << "Finding initial Jukes-Cantor pairwise distances ...";
  MatrixXd jcDistanceMatrix(alignment.getNumTaxa(),alignment.getNumTaxa());
  alignment.calculateJCDistances(jcDistanceMatrix);
  cerr << " done." << endl;

  cout << "Jukes-Cantor Distance Matrix:" << endl;
  cout << endl << jcDistanceMatrix << endl << endl;

  milliseconds ms3 = duration_cast< milliseconds >( system_clock::now().time_since_epoch() );

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

  // for ( int i=0; i<4; ++i )
  // {
  //   double logFoo;
  //   jctree.randomEdges(alignment,model,rng,logFoo,true);
  // }

  milliseconds ms4 = duration_cast< milliseconds >( system_clock::now().time_since_epoch() );

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

  milliseconds ms5 = duration_cast< milliseconds >( system_clock::now().time_since_epoch() );

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
  jctree.mcmc(q_init,alignment,(MCMC_Q_BURN),alignment.getNumSites(),rng, true);
  cerr << endl << " done." << endl;
  cout << "JC tree: " << endl;
  cout << jctree.makeTreeNumbers() << endl;

  // mcmc to get final Q
  cerr << "sampling:" << endl;
  jctree.mcmc(q_init,alignment,(MCMC_Q_SAMPLE),alignment.getNumSites(),rng, false);
  cerr << endl << " done." << endl;
  cout << "JC tree: " << endl;
  cout << jctree.makeTreeNumbers() << endl;


//  QMatrix model(parameters.getStationaryP(),parameters.getSymmetricQP());
  QMatrix model(q_init.getStationaryP(),q_init.getSymmetricQP());
  model.setMcmcVarP(q_init.getMcmcVarP());
  model.setMcmcVarQP(q_init.getMcmcVarQP());
  cerr << endl << " done." << endl;

  milliseconds ms6 = duration_cast< milliseconds >( system_clock::now().time_since_epoch() );


  // Recalculate pairwise distances using estimated Q matrix (TODO: add site rate heterogeneity)
  cerr << "Finding initial GTR pairwise distances ...";
  MatrixXd gtrDistanceMatrix(alignment.getNumTaxa(),alignment.getNumTaxa());
  alignment.calculateGTRDistances(model,jcDistanceMatrix,gtrDistanceMatrix);
  cerr << " done." << endl;

  cout << "GTR Distance Matrix:" << endl;
  cout << endl << gtrDistanceMatrix << endl << endl;

  milliseconds ms7 = duration_cast< milliseconds >( system_clock::now().time_since_epoch() );

  MatrixXd gtrDistanceMatrixCopy = gtrDistanceMatrix;
  Tree gtrtree(gtrDistanceMatrixCopy);
  gtrtree.reroot(1);
  gtrtree.sortCanonical();

  cout << "GTR Tree" << endl;
  cout << gtrtree.makeTreeNumbers() << endl;

  for ( int i=0; i<4; ++i )
  {
    double logFoo;
    gtrtree.randomEdges(alignment,model,rng,logFoo,true);
  }

  cout << "MLE Branch Lengths" << endl;
  cout << gtrtree.makeTreeNumbers() << endl;

  // Do bootstrap with new pairwise distances
  if( parameters.getNumBootstrap() == 0)
    parameters.setNumBootstrap(1000);
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
    // add parsimony score to map if it is not already there
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

  milliseconds ms8 = duration_cast< milliseconds >( system_clock::now().time_since_epoch() );

  map<string,double> topologyToWeightMap; // not normalized; normalization taken care of during alias creation
  map<string,double> topologyToLogLikWeightMap; // not normalized; normalization taken care of during alias creation
  for ( map<string,int>::iterator m=topologyToCountMap.begin(); m != topologyToCountMap.end(); ++m )
  {
    topologyToWeightMap[ (*m).first ] =
      (*m).second * exp( parameters.getParsimonyScale() *
			 (minimumParsimonyScore - topologyToParsimonyScoreMap[ (*m).first ]) );
  }

  cout << endl << "Topology counts to file" << endl;
  {
    string topologyCountsFile = parameters.getOutFileRoot() + ".topCounts";
    ofstream topCounts(topologyCountsFile.c_str());
    map<string,double>::iterator wm=topologyToWeightMap.begin();
    topCounts << "tree count parsimonyWt parsimonyScore parsimonyDiff" << endl;
    for ( map<string,int>::iterator cm=topologyToCountMap.begin(); cm != topologyToCountMap.end(); ++cm )
    {
      Tree t((*cm).first);
      t.relabel(alignment);
      t.unroot();
      t.reroot(1);
      t.sortCanonical();
      string top = t.makeTopologyNumbers();
      topCounts << top << " " << setw(5) << (*cm).second << " ";
      topCounts << setw(10) << setprecision(4) << fixed << (*wm).second;
      topCounts << " " << setw(5) << topologyToParsimonyScoreMap[ (*cm).first ] << " " << setw(4) << minimumParsimonyScore - topologyToParsimonyScoreMap[ (*cm).first ];
      topCounts << " " << endl;
      ++wm;
    }
    topCounts << endl;
    topCounts.close();
  }

  vector<int> taxaNumbers;
  vector<string> taxaNames;
  alignment.getTaxaNumbersAndNames(taxaNumbers,taxaNames);

  milliseconds ms9 = duration_cast< milliseconds >( system_clock::now().time_since_epoch() );

  CCDProbs<int> ccd(topologyToCountMap,taxaNumbers,taxaNames);
  CCDProbs<double> ccdParsimony(topologyToWeightMap,taxaNumbers,taxaNames);
  CCDProbs<double> ccdLogLik(topologyToLogLikWeightMap,taxaNumbers,taxaNames);
  // write map out to temp files to check
  string originalSmapFile = parameters.getOutFileRoot() + "-nopars.smap";
  string originalTmapFile = parameters.getOutFileRoot() + "-nopars.tmap";
  string parsimonySmapFile = parameters.getOutFileRoot() + "-pars.smap";
  string parsimonyTmapFile = parameters.getOutFileRoot() + "-pars.tmap";
  ofstream smap(originalSmapFile.c_str());
  ccd.writeCladeCount(smap);
  smap.close();
  ofstream tmap(originalTmapFile.c_str());
  ccd.writePairCount(tmap);
  tmap.close();
  smap.open(parsimonySmapFile);
  ccdParsimony.writeCladeCount(smap);
  smap.close();
  tmap.open(parsimonyTmapFile);
  ccdParsimony.writePairCount(tmap);
  tmap.close();

  return 0;
}
