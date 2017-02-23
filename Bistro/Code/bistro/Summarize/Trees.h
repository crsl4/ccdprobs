#ifndef TREESH
#define TREESH

using namespace std;

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <math.h>
#include <cstdlib>
#include <algorithm>

#include "CladeTable.h"
#include "SplitTable.h"
#include "Tree.h"
#include "CladeTree.h"
class Trees {

public:

  double Found, NotFound, CompFound, CompNotFound;

  Trees() : maxTrees(10000), numTrees(0), totalTrees(0) {};
    
  inline double minimum(double a, double b) { return a<b ? a : b; }

  inline double sqr(double x) { return x * x; }

  void readFiles(const vector<string>& names, int skip);

  void readTrees(const vector<string>& names);

  void findClades(double threshold);

  void findNamedClades(double threshold, int maxTopologies);

  void printMeanTree(ostream &c);

  void printNamedClades(ostream& c);

  void printClades(ostream& c, double threshold);

  void printTreeTopologies(ostream& c, int m);

  void printCladeTreeTopologies(ostream& c, int m);

  void printTrans(ostream& c, int maxTopologies);

  void printProbableTreeClades(ostream& c, double threshold) ;

  void printDistanceMatrices(ostream &c);

private:

  int maxTrees;				           // maximum number of trees.
  int numTrees;					   // number of unique trees of two or more taxa.
  int totalTrees;				   // number of trees read in (including duplicates).
  int maxLen;				           // length of a whole tree, i.e., number of taxa.
  bool hasLength;				   // whether or not the input contains branch lengths.
  double blSquared;				   // Sum of the squares of the branchlengths in 
  						   // all the trees.
  vector<TreePtr> allTrees;		           // list of trees read in.
  vector< vector<TreePtr> > trees;                 // trees[i] is the list of unique trees of size i.
  vector< vector<CladePtr> > clades;	           // clades[i] is the list of clades of size i that 
  						   // have sufficient counts.
  vector<CladePtr> namedClades;                    // list of named clades.
  Set cladeSet;					   // Temporary cladeSet (so we don't create it each time).
  vector<CladeTable> cladeTable;		   // Table of clades of each size.
  SplitTable splitTable;			   // Table of the splits.
  vector< vector<double> > leafDist;		   // Matrix of distances between leaves.
  vector< vector<double> > leafDistSqr;		   // Matrix of squared distances between leaves.

  void findSplits();

  void addEdge(CladePtr parent, CladePtr first, CladePtr second, int count, vector<int> offsets);

  void readFile(istream& f, const string& name, int skip, vector<TreePtr>& sindex, bool& first, 
		vector<TreePtr>& taxa, bool& hasLength, bool& hasSemicolon, double& blSquared);

  void readTree(string str, vector<TreePtr>& sindex, bool& first, 
		vector<TreePtr>& taxa, bool& hasLength, bool& hasSemicolon, double& blSquared);

  void syntaxError(const string& msg, const string& line, string::const_iterator pos);

  void initializeTrees(string str, vector<TreePtr>& taxa, bool& hasLength, bool& hasSemicolon);

  void skipLength(string::const_iterator& top, const string& str, bool& hasLength);

  string::const_iterator getTaxaTop(string::const_iterator top, const string& str, vector<int>& taxaNames,
				    bool& hasLength, bool& hasSemicolon);

  string::const_iterator getTaxa(string::const_iterator top, const string& str, vector<int>& taxaNames,
				 bool& hasLength);

  bool isNumeric(char ch);

  void readLength(string::const_iterator& top, const string& str, double& bl);

  string::const_iterator storeTopWithLen(string::const_iterator top, const string& str, 
					 TreePtr& tree, vector<TreePtr>& sindex, 
					 const vector<TreePtr>& taxa, bool hasSemicolon, 
					 double& blSquared);

  string::const_iterator storeWithLen(string::const_iterator top, const string& str, 
				      TreePtr& tree, unsigned int& hash, 
				      vector<TreePtr>& sindex, const vector<TreePtr>& taxa,
				      double &bl, double& blSquared,
				      vector<double> &dist, vector<int> &leaves);

  string::const_iterator storeTop(string::const_iterator top, const string& str, 
				  TreePtr& tree, vector<TreePtr>& sindex, 
				  const vector<TreePtr>& taxa, bool hasSemicolon, double& blSquared);

  string::const_iterator store(string::const_iterator top, const string& str, 
			       TreePtr& tree, unsigned int& hash, int& len,
			       vector<TreePtr>& sindex, const vector<TreePtr>& taxa);

  TreePtr add(unsigned int &hash, TreePtr ltree, TreePtr rtree, 
	      vector<TreePtr>& sindex, int llen, int rlen, double branchLength);

  inline unsigned int hashfn(unsigned int a, unsigned int b, size_t s) const { return (((a * 257) + b) % s); }

  void printTotalDist(ostream& c, CladePtr j);

  void printMeanTree2(ostream &c, CladePtr j, int topTree);
	
  void prettyPrintMeanTree(ostream &c, CladePtr j, int topTree, int indent, int stayOnLine);
	    
  void addTrans(TreePtr tree, int time, vector<int>& lastTime, vector<int>& lastTopNum, 
		vector< vector< vector<int> > >& trans);

  void collectClades(TreePtr tree, int votesNeeded, vector<CladePtr>& treeClades);

};

#endif
