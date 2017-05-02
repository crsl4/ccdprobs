// Bret Larget and Claudia Solis-Lemus
// distances -f fasta --mbfile mb.tre --bistrofile bistro---0-249.treeBL -s seed --skip int
// fixit: we want to give the root, not each filename
// the code will combine mb.run1.tre and mb.run2.tre (gotten after mb2badger),
// and it will also combine bistro---*.treeBL
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
  cerr << "Found " << badTrees << " bad trees for the computation of pairwise distance matrix" << endl;
  // cerr << "Tree distance matrix: " << endl << endl;
  // cerr << mat << endl;
}

vector<string> sampleTrees(vector<string> trees, int sampleSize, int seed,mt19937_64& rng)
{
  random_shuffle( trees.begin(), trees.end()); //fixit: not using seed
  vector<string>::const_iterator first = trees.begin();
  vector<string>::const_iterator last = trees.begin() + sampleSize - 1;
  vector<string> newtrees(first, last);
  return newtrees;
}

vector<string> readFile(istream& f, const string& name, int skip, bool skipTrees)
{
  string str;
  vector<string> trees;
  if(skipTrees)
  {
    for(int i=0;i<skip;i++)
      if(!getline(f,str)) {
	cerr << "Error: " << (name=="-" ? "Standard input" : ("Input file " + name))
	     << " contains " << (i==skip-1 ? "exactly" : "fewer than ") << skip << " lines." << endl;
	exit(1);
      }
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
  mt19937_64 rng(parameters.getSeed());
  int skip = parameters.getSkip();

  // here we read the trees in the two files
  string mbname = parameters.getMBfile();
  ifstream f(mbname.c_str());
  if(!f) {
    cerr << "Error: Could not open " << mbname << endl;
    exit(1);
  }
  vector<string> mbtrees = readFile(f,mbname,skip,true);
  cerr << "Read MrBayes file: " << mbname << endl;
  cerr << "read " << mbtrees.size() << " trees" << endl;
  f.close();


  string bistroname = parameters.getBistrofile();
  ifstream f2(bistroname.c_str());
  if(!f2) {
    cerr << "Error: Could not open " << bistroname << endl;
    exit(1);
  }
  vector<string> bistrotrees = readFile(f2,bistroname,skip,false);
  cerr << "Read Bistro file: " << bistroname << endl;
  cerr << "read " << bistrotrees.size() << " trees" << endl;
  f2.close();


  // mean tree
  CladeGraph MBcladeGraph;
  CladeGraph BistrocladeGraph;
  cerr << "Computing MrBayes mean tree" << endl;
  MBcladeGraph.findMeanTree(mbtrees,alignment); //warning: this will fail is fasta and nexus files have different translate table
  cerr << "done." << endl;
  cerr << "Computing Bistro mean tree" << endl;
  BistrocladeGraph.findMeanTree(bistrotrees,alignment);
  cerr << "done." << endl;
  string MBmeanTree = MBcladeGraph.getMeanTree();
  string BistromeanTree = BistrocladeGraph.getMeanTree();
  cerr << "MB meanTree = " << MBmeanTree << endl;
  cerr << "Bistro meanTree = " << BistromeanTree << endl;

  // now we need to choose a sample of 100 from the mbtrees and bistrotrees
  int sampleSize = 100;
  cerr << "Getting random sample of " << sampleSize << " trees" << endl;
  vector<string> newBistrotrees;
  if(bistrotrees.size() > sampleSize)
    newBistrotrees = sampleTrees(bistrotrees, sampleSize, parameters.getSeed(),rng);
  else
    newBistrotrees = bistrotrees;
  vector<string> newMBtrees;
  if(mbtrees.size() > sampleSize)
    newMBtrees = sampleTrees(mbtrees, sampleSize, parameters.getSeed(),rng);
  else
    newMBtrees = mbtrees;

  Tree mbtree(MBmeanTree,alignment);
  Tree bistrotree(BistromeanTree,alignment);
  MatrixXd MBtreeDistanceMatrix = MatrixXd::Zero(sampleSize+1,sampleSize+1);
  MatrixXd BistrotreeDistanceMatrix = MatrixXd::Zero(sampleSize+1,sampleSize+1);
  cerr << "Distance for MrBayes trees" << endl;
  distanceMatrix(MBtreeDistanceMatrix,mbtree,newMBtrees,alignment);
  cerr << "done." << endl;
  cerr << "Distance for Bistro trees" << endl;
  distanceMatrix(BistrotreeDistanceMatrix,bistrotree,newBistrotrees,alignment);
  cerr << "done." << endl;

  string mbdistances = parameters.getMBfile() + ".distances";
  cerr << "Writing MrBayes matrix to file " << mbdistances << endl;
  ofstream mbstream(mbdistances.c_str());
  mbstream << MBtreeDistanceMatrix << endl;
  mbstream.close();
  cerr << "done." << endl;
  string bistrodistances = parameters.getBistrofile() + ".distances";
  cerr << "Writing Bistro matrix to file " << bistrodistances << endl;
  ofstream bistrostream(bistrodistances.c_str());
  bistrostream << BistrotreeDistanceMatrix << endl;
  bistrostream.close();
  cerr << "done." << endl;
}


