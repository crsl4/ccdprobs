// Bret Larget and Claudia Solis-Lemus
// distances -f fasta --mbfile mb.tre --bistrofile bistro.treeBL -s seed

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
#include "cladegraph.h"

#include "Eigen/Core"
#include "Eigen/Eigenvalues"
#include "Eigen/SVD"

using namespace std;
using namespace std::chrono;
using namespace Eigen;

void distanceMatrix(MatrixXd& mat,Tree mtree, vector<string> trees, Alignment alignment)
{
  int badTrees = 0;
  for ( int i = 0; i<trees.size(); ++i)
  {
    string fooi = trees[i];
    if ( fooi.find("nan") != string::npos )
    {
      ++badTrees;
      continue;
    }
    Tree* treei = new Tree(fooi, alignment);
    for ( int j = i+1; j<trees.size(); ++j)
    {
      string fooj = trees[j];
      if ( fooj.find("nan") != string::npos )
      {
	++badTrees;
	continue;
      }
      Tree* treej = new Tree(fooj, alignment);
      double dd = treei->distance(treej);
      mat(i,j) = dd;
      mat(j,i) = dd;
    }
    double d = mtree.distance(treei);
    mat(i,trees.size()) = d;
    mat(trees.size(),i) = d;
  }
  cerr << "Found " << badTrees << " for the pairwise distance matrix" << endl;
  cerr << "Tree distance matrix: " << endl << endl;
  cerr << mat << endl;
}

vector<string> sampleTrees(vector<string> trees, int sampleSize, int seed)
{
  random_shuffle ( trees.begin(), trees.end(), seed );
  vector<string>::const_iterator first = trees.begin();
  vector<string>::const_iterator last = trees.begin() + sampleSize - 1;
  vector<string> newtrees(first, last);
  return newtrees;
}

int main(int argc, char* argv[])
{
  // Read command line and process parameters
  cerr << "Processing command line ...";
  Parameter parameters;
  parameters.processCommandLine(argc,argv);
  cerr << " done." << endl;

  // fixit: we are reading the sequences because we need them for the Tree constructor
  // Read in sequences from FASTA file
  cerr << "Reading alignment from FASTA file ...";
  Alignment alignment(parameters.getSequenceFileName());
  cerr << " done." << endl;
  cout << endl << "Alignment summary:" << endl;
  alignment.summarize(cout);
  cout << endl;
  if ( parameters.getSeed() == 0 )
  {
    random_device rd;
    parameters.setSeed(rd());
  }


  // here we need to read the trees in the two files
  vector<string> mbtrees;
  vector<string> bistrotrees;

  // mean tree
  CladeGraph MBcladeGraph;
  CladeGraph BistrocladeGraph;
  MBcladeGraph.findMeanTree(mbtrees,alignment); //warning: this will fail is fasta and nexus files have different translate table
  BistrocladeGraph.findMeanTree(bistrotrees,alignment);
  string MBmeanTree = MBcladeGraph.getMeanTree();
  string BistromeanTree = BistrocladeGraph.getMeanTree();
  cerr << "MB meanTree = " << MBmeanTree << endl;
  cerr << "Bistro meanTree = " << BistromeanTree << endl;

  // now we need to choose a sample of 100 from the mbtrees and bistrotrees
  int sampleSize = 100;
  vector<string> newBistroTrees;
  if(bistrotrees.size() > sampleSize)
    newBistroTrees = sampleTrees(bistrotrees, sampleSize, parameters.getSeed());
  else
    newBistroTrees = bistrotrees;
  vector<string> newMBTrees;
  if(mbtrees.size() > sampleSize)
    newMBTrees = sampleTrees(mbtrees, sampleSize, parmeters.getSeed());
  else
    newMBTrees = mbtrees;

  Tree mbtree(MBmeanTree,alignment);
  Tree bistrotree(BistromeanTree,alignment);
  MatrixXd MBtreeDistanceMatrix = MatrixXd::Zero(sampleSize+1,sampleSize+1);
  MatrixXd BistrotreeDistanceMatrix = MatrixXd::Zero(sampleSize+1,sampleSize+1);
  distanceMatrix(MBtreeDistanceMatrix,mbtree,newMbtrees,alignment);
  distanceMatrix(BistrotreeDistanceMatrix,bistrotree,newBistrotrees,alignment);
  string mbdistances = parameters.getMBfile() + ".distances";
  ofstream mbstream(mbdistances.c_str());
  mbstream << MBtreeDistanceMatrix << endl;
  mbstream.close();
  string bistrodistances = parameters.getBistrofile() + ".distances";
  ofstream bistrostream(bistrodistances.c_str());
  bistrostream << BistrotreeDistanceMatrix << endl;
  bistrostream.close();
}


