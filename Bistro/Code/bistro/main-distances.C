// Bret Larget and Claudia Solis-Lemus
// distances -f fasta --mbfile mb.tre --bistrofile bistro.treeBL -s seed
// for now we assume that the fasta file (used in bistro) and the nexus file
// (used in mrbayes) have the same translate table

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

vector<string> readFile(istream& f, const string& name, int skip)
{
  string str;
  vector<string> trees;
  for(int i=0;i<skip;i++)
    if(!getline(f,str)) {
      cerr << "Error: " << (name=="-" ? "Standard input" : ("Input file " + name))
	   << " contains " << (i==skip-1 ? "exactly" : "fewer than ") << skip << " lines." << endl;
      exit(1);
    }
  if(getline(f,str)) {
    do {
      trees.push_back(str);
    }
    while(getline(f,str));
  }
  return trees;
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

// need to add here the parameter for skip

  // here we read the trees in the two files
  string mbname = parameters.getMBfile();
  ifstream f(mbname.c_str());
  if(!f) {
    cerr << "Error: Could not open " << mbname << endl;
    exit(1);
  }
  f.close();
  vector<string> mbtrees = readFile(f,mbname,skip);

  string bistroname = parameters.getBistrofile();
  ifstream f(bistroname.c_str());
  if(!f) {
    cerr << "Error: Could not open " << bistroname << endl;
    exit(1);
  }
  f.close();
  vector<string> bistrotrees = readFile(f,bistroname,skip);

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


