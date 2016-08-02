#ifndef __CCDPROBS_H
#define __CCDPROBS_H

// Bret Larget
// University of Wisconsin, Madison
//
// July 24, 2012
// February 20, 2013
// May 23, 2016
//
// Copyright (C) 2012, 2013, 2016

// ccdprobs adapted for Bistro

#include <iostream>
#include <iomanip>
#include <vector>
#include <set>
#include <string>
#include <sstream>
#include <fstream>
#include <cctype>
#include <map>
#include <cmath>
#include <ctime>
#include <random>

#include "alias.h"
#include "tree.h"

using namespace std;

// Program uses the BOOST dynamic bitset library in order to
// store clades more efficiently and without a limit on the size of the tree
#include <boost/dynamic_bitset.hpp>

using namespace boost;

// ************************************************************
// Clade and Tree stuff from tree string
// 
// The data in a Clade is a dynamic bitset of size numTaxa
// where clade[i-1] = 1 if taxon i is in the clade.
// ************************************************************
class Clade {
private:
  dynamic_bitset<unsigned char> clade; // one bit per taxa 0/1
public:
  Clade() {}
  Clade(int n) { clade.resize(n); }
  ~Clade() { clade.clear(); }
  Clade(dynamic_bitset<unsigned char> x) { clade=x; }
  dynamic_bitset<unsigned char> get() const { return clade; }
  int count() const { return clade.count(); } // number of 1s
  int size() const { return clade.size(); } // number of bits
  void clear() { clade.clear(); }
  void set(dynamic_bitset<unsigned char> x) { clade=x; }
  void resize(int n) { clade.resize(n); }
  void add(int x) { clade[size() - x] = 1; } // taxa #x stored at size-x (1 is at size-1, numTaxa is at 0).
//  void add(int x) { clade[x-1] = 1; }
  void add(Clade c) { clade |= c.get(); }
//  void subtract(int x) { clade[x-1] = 0; }
  void print(ostream&) const ;
  friend bool operator< (const Clade&,const Clade&);
  friend bool operator> (const Clade&,const Clade&);
  friend bool operator== (const Clade&,const Clade&);
  friend Clade operator- (const Clade&,const Clade&);
  string randomTree(multimap<Clade,pair<Clade,int> >&,
		    map<Clade,Alias<dynamic_bitset<unsigned char> >* >&,
		    map<pair<dynamic_bitset<unsigned char>,dynamic_bitset<unsigned char> >,double>&,
		    mt19937_64&,
		    double&);
};

// ############################################################

class CCDTree {
private:
  string top;
  int numTaxa;
  vector<Clade> clades;
public:
  CCDTree(string,int);
  ~CCDTree() { 
    for ( vector<Clade>::iterator p=clades.begin(); p != clades.end(); p++ )
      (*p).clear();
    clades.clear();
  }
  void print(ostream& f) { f << top; }
  void printClades(ostream& f)
  {
    vector<Clade>::iterator p = clades.begin();
    (*p).print(f);
    ++p;
    for (; p != clades.end(); ++p)
      (*p).print(f);
    f << endl;
  }    
  string getTop() { return top; }
  int getNumTaxa() { return numTaxa; }
  int getNumClades() { return clades.size(); }
  Clade getClade(int i) { return clades[i]; }
};

// ************************************************************
// CladePair stuff
// 
// A clade pair will be two clades where the second clade is a strict subset of the first clade
// and the second clade will precede the clade which is the set difference between the first and second clades.
// ************************************************************
class CladePair {
private:
  Clade clade1,clade2;
public:
  CladePair(Clade c1,Clade c2) {
    clade1 = c1;
    clade2 = c2;
  }
  Clade getClade1() const { return clade1; }
  Clade getClade2() const { return clade2; }
  void print(ostream& f) const {
    f << "[[";
    clade1.print(f);
    f << ",";
    clade2.print(f);
    f <<  "]]";
  }

  friend bool operator< (const CladePair&,const CladePair&);
};

// ************************************************************
// Rooted Binary Tree stuff
// ************************************************************

class RootedNode {
private:
  RootedNode* parent;
  RootedNode* left;
  RootedNode* right;
  bool leaf;
  bool root;
  Clade clade;
public:
  RootedNode() { left = NULL; right = NULL; parent = NULL; leaf=false; root=false; }
  RootedNode(RootedNode* p) { left = NULL; right = NULL; parent = p; }
  RootedNode* getParent() const { return parent; }
  RootedNode* getLeft() const { return left; }
  RootedNode* getRight() const { return right; }
  bool getLeaf() const { return leaf; }
  bool getRoot() const { return root; }
  Clade getClade() const { return clade; }
  void setParent(RootedNode* p) { parent = p; }
  void setLeft(RootedNode* p) { left = p; }
  void setRight(RootedNode* p) { right = p; }
  void setLeaf(bool b) { leaf = b; }
  void setRoot(bool b) { root = b; }
  void addChild(RootedNode* c) {
    if ( left==NULL )
      left = c;
    else if ( right==NULL )
      right = c;
    else {
      cerr << "Error: input tree is not binary." << endl;
      exit(1);
    }
  }
  void add(int x) { clade.add(x); }
  void add(Clade c) { clade.add(c); }
  void copyClade(Clade c) { clade = c; }
  void setNumTaxa(int n) { clade.resize(n); }
  string randomTree(multimap<Clade,pair<Clade,int> >&,map<Clade,Alias<dynamic_bitset<unsigned char> >* >&,mt19937_64&,double&);
};

class RootedTree
{
private:
  string top;
  string binaryTop;
  int numTaxa;
  RootedNode* root;
  vector<RootedNode*> nodes;
  vector<Clade> clades;
public:
  RootedTree() {};
  RootedTree(string,int);
  RootedNode* getRoot() const { return root; }
  int getNumTaxa() const { return numTaxa; }
  string getTop() const { return top; }
  string getBinaryTop() const { return binaryTop; }
  int getNumClades() { return clades.size(); }
  Clade getClade(int i) { return clades[i]; }
  void print(ostream&);
  void printClades(ostream&);
  void count(int,map<Clade,int>&,map<CladePair,int>&);
  double estimateProbability(int, map<Clade,int>&, map<CladePair,int>&);
};

class CCDProbs
{
private:
  int sampleSize;
  int numTaxa;
  map<Clade,int> cladeCount;
  map<CladePair,int> pairCount;
  vector<int> taxaNumbers;
  vector<string> taxaNames;
  Clade all;
public:
  multimap<Clade,pair<Clade,int> > mm; //one clade (key) has many pairs clade,int (values): ways to split
  map<Clade,Alias<dynamic_bitset<unsigned char> >* > am;
  // add a map from pair<Clade,Clade> to double to store the log of the probability of selecting the second clade from the first
  map<pair<dynamic_bitset<unsigned char>,dynamic_bitset<unsigned char> >,double> cladeLogProbMap;
public:
  CCDProbs(map<string,int>&,vector<int>&,vector<string>&);
  void writeTranslateTable(ostream&);
  void writeCladeCount(ostream&);
  void writePairCount(ostream&);
  string randomTree(mt19937_64&,double&);
//  Tree rtop(mt19937_64&);
};


#endif
