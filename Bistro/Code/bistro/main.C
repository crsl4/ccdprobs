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

  // Find Initial Neighbor-joining tree
  cerr << "Finding initial neighbor-joining tree ... ";
  Tree tree(jcDistanceMatrix);
  tree.reroot(1);
  tree.sortCanonical();
  cerr << " done." << endl;
  cerr << endl << "Tree topology:" << endl;
  cerr << tree.makeTopologyNumbers() << endl << endl;

  // Run MCMC on tree to estimate Q matrix parameters

  // set up Q matrix
  
  // QMatrix model(parameters.getStationaryP(),parameters.getSymmetricQP());

  // Recalculate pairwise distances using estimated Q matrix (TODO: add site rate heterogeneity)

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

  // Do bootstrap with new pairwise distances
  cerr << "Beginning " << parameters.getNumBootstrap() << " bootstrap replicates ..." << endl;
  map<string,int> topologyToCountMap;
  MatrixXd bootDistanceMatrix(alignment.getNumTaxa(),alignment.getNumTaxa());
  vector<int> weights(alignment.getNumSites());
  for ( int b=0; b<parameters.getNumBootstrap(); ++b )
  {
    if ( (b+1) % (parameters.getNumBootstrap() / 100) == 0 )
      cerr << '*';
    alignment.setBootstrapWeights(weights,rng);
    alignment.calculateJCDistancesUsingWeights(weights,bootDistanceMatrix);
    Tree bootTree(bootDistanceMatrix);
    bootTree.reroot(1);
    bootTree.sortCanonical();
    topologyToCountMap[ bootTree.makeTopologyNumbers() ]++;
  }
  cerr << '*' << " ... done." << endl;

  cerr << endl << "Topology counts:" << endl;
  for ( map<string,int>::iterator m=topologyToCountMap.begin(); m != topologyToCountMap.end(); ++m )
    cout << (*m).first << " " << setw(5) << (*m).second << endl;
  cerr << endl;

  vector<int> taxaNumbers;
  vector<string> taxaNames;
  alignment.getTaxaNumbersAndNames(taxaNumbers,taxaNames);
  CCDProbs ccd(topologyToCountMap,taxaNumbers,taxaNames);

  int numRandom = 10000;
  cerr << "Generating " << numRandom << " random trees:" << endl;
  for ( int k=0; k<numRandom; ++k )
  {
    cout << ccd.randomTree(rng) << endl;
  }
  
  return 0;
}
