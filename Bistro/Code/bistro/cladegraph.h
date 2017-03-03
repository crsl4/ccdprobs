#ifndef __CLADEGRAPH_H
#define __CLADEGRAPH_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <algorithm>
#include <sstream>

#include "tree.h"

#include <boost/dynamic_bitset.hpp>

using namespace boost;
using namespace std;

class CladeNode
{
private:
  dynamic_bitset<unsigned char> clade; // taxa numbers - 1
  bool leaf; // true if clade has only one taxon
  set<pair<dynamic_bitset<unsigned char>,dynamic_bitset<unsigned char> > > subclades; // list of all subclade pairs
  CladeNode* left; // the best left subclade
  CladeNode* right; // the other clade to the best left subclade
  string subtree; // string of subtree with branch lengths
  double value; // for a leaf, meanLength^2; otherwise, sum of meanLength^2 + values of left and right
  double sumOfLengths; // sum of all split lengths from input trees
  double meanLength; // sumOfLengths divided by total number of trees
public:
  CladeNode(dynamic_bitset<unsigned char>);
  dynamic_bitset<unsigned char> getClade() { return clade; }
  bool getLeaf() { return leaf; }
  dynamic_bitset<unsigned char> otherSubClade(dynamic_bitset<unsigned char> x) { return clade - x; }
  CladeNode* getLeft() { return left; }
  void setLeft(CladeNode* x) { left = x; }
  CladeNode* getRight() { return right; }
  void setRight(CladeNode* x) { right = x; }
  string getSubtree() { return subtree; }
  void setSubtree(bool);
  double getValue() { return value; }
  void setValue(map<dynamic_bitset<unsigned char>,CladeNode*>);
  double getSumOfLengths() { return sumOfLengths; }
  void addLength(double x) { sumOfLengths += x; }
  double getMeanLength() { return meanLength; }
  void setMeanLength(double denom) { meanLength = sumOfLengths / denom; }
};

class CladeGraph
{
private:
  int numTaxa; // number of taxa
  int numTrees; // number of trees summarized in the graph
  CladeNode* root; // cladeNode corresponding to all taxa
  map<dynamic_bitset<unsigned char>,CladeNode*> cladeMap;
  multimap<int,CladeNode*> sizeMap;
  string meanTree;
public:
  CladeGraph() {};
  void processTrees(vector<string>);
  void processTree(Tree*,map<dynamic_bitset<unsigned char>,CladeNode*>&,multimap<int,CladeNode*>&);
  void setMeanLengths();
  void setValues();
  string getMeanTree() { return meanTree; }
  void findMeanTree(vector<string>);
};

#endif
