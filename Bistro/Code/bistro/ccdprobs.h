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

// question: do we need to put T in all classes?
// for Clade, only randomTree uses it, still need to define Clade<T>?
// it seems that a member template is enough: http://www.cplusplus.com/forum/general/46654/
// http://en.cppreference.com/w/cpp/language/member_template
// http://www.tutorialspoint.com/cplusplus/cpp_templates.htm
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
  void add(int x) { clade[x-1] = 1; }
  //void add(int x) { clade[size() - x] = 1; } // taxa #x stored at size-x (1 is at size-1, numTaxa is at 0).
//  void add(int x) { clade[x-1] = 1; }
  void add(Clade c) { clade |= c.get(); }
//  void subtract(int x) { clade[x-1] = 0; }
  void print(ostream&) const ;
  friend bool operator< (const Clade&,const Clade&);
  friend bool operator> (const Clade&,const Clade&);
  friend bool operator== (const Clade&,const Clade&);
  friend Clade operator- (const Clade&,const Clade&);
  template<typename T>
  string randomTree(multimap<Clade,pair<Clade,T> >&, // T instead of int
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
  //template<typename T> with T instead of int: not implemented function
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
  template<typename T> //with T instead of all ints
  void count(T,map<Clade,T>&,map<CladePair,T>&);
  double estimateProbability(int, map<Clade,int>&, map<CladePair,int>&); //not implemented
};

template <typename T>
class CCDProbs
{
private:
  T sampleSize;
  int numTaxa;
  map<Clade,T> cladeCount; //T instead of int
  map<CladePair,T> pairCount; //T instead of int
  vector<int> taxaNumbers;
  vector<string> taxaNames;
  Clade all;
public: //T instead of int in mm
  multimap<Clade,pair<Clade,T> > mm; //one clade (key) has many pairs clade,int (values): ways to split
  map<Clade,Alias<dynamic_bitset<unsigned char> >* > am;
  // add a map from pair<Clade,Clade> to double to store the log of the probability of selecting the second clade from the first
  map<pair<dynamic_bitset<unsigned char>,dynamic_bitset<unsigned char> >,double> cladeLogProbMap;
public:
  CCDProbs(map<string,T>&,vector<int>&,vector<string>&); //T instead of 1st int in map
  void writeTranslateTable(ostream&);
  void writeCladeCount(ostream&);
  void writePairCount(ostream&);
  string randomTree(mt19937_64&,double&);
//  Tree rtop(mt19937_64&);
};

// to add n to count of clade
template<typename T> //and change int for T
void RootedTree::count(T n,map<Clade,T>& cladeCount,map<CladePair,T>& pairCount)
{
  print(cout); // print tree info
  for ( vector<RootedNode*>::iterator p = nodes.begin(); p != nodes.end(); p++ ) {
    if ( (*p)->getLeaf() )
      continue;
    Clade z = (*p)->getClade();
    z.print(cout);
    cout << "added count " << n << endl;
    cladeCount[z] += n;
    Clade x = (*p)->getLeft()->getClade();
    Clade y = (*p)->getRight()->getClade();
    if ( x > y )
      x = y;
    CladePair c(z,x);
    pairCount[c] += n;
  }
}

template<typename T> //change int for T and put CCDProbs<T>
CCDProbs<T>::CCDProbs(map<string,T>& topologyToCountMap,vector<int>& taxaNumbers,vector<string>& taxaNames)
{
  sampleSize = 0;
  numTaxa = taxaNames.size();
  for ( typename map<string,T>::iterator m=topologyToCountMap.begin(); m != topologyToCountMap.end(); ++m )
  {
    RootedTree rt(m->first,numTaxa);
    cout << "Rooted tree top and binarytop" << endl;
    cout << rt.getTop() << endl;
    cout << rt.getBinaryTop() << endl;
    rt.count<T>(m->second,cladeCount,pairCount);
    sampleSize += m->second;
  }

  for ( typename map<CladePair,T>::iterator p=pairCount.begin(); p!=pairCount.end(); ++p ) {
    Clade parent=(p->first).getClade1();
    Clade child=(p->first).getClade2();
    mm.insert( pair<Clade,pair<Clade,T> >(parent,make_pair(child, p->second)) );
  }

  all.resize(taxaNames.size());
  for ( int k=1; k<=numTaxa; ++k )
    all.add(k);
}

template<typename T>
void CCDProbs<T>::writeTranslateTable(ostream& f)
{
  f << "translate" << endl;
  for ( int i=0; i < taxaNumbers.size(); ++i ) {
    f << setw(8) << taxaNumbers[i] << " " << taxaNames[i];
    if ( i < taxaNumbers.size()-1 )
      f << "," << endl;
    else
      f << ";" << endl;
  }
}

template<typename T>
void CCDProbs<T>::writeCladeCount(ostream& f)
{
  writeTranslateTable(f);
  f.setf(ios::fixed,ios::floatfield);
  for ( typename map<Clade,T>::iterator p=cladeCount.begin(); p != cladeCount.end(); p++ ) {
    f << setw(10) << setprecision(8) << p->second / (double) sampleSize << " ";
    f << setw(10) << p->second << " ";
    p->first.print(f);
    f << endl;
  }
}

template<typename T>
void CCDProbs<T>::writePairCount(ostream& f)
{
  writeTranslateTable(f);
  f.setf(ios::fixed,ios::floatfield);
  for ( typename map<CladePair,T>::iterator p=pairCount.begin(); p != pairCount.end(); p++ ) {
    f << setw(10) << setprecision(8) << p->second / (double) sampleSize << " ";
    f << setw(10) << p->second << " ";
    f << "[ ";
    (p->first).getClade1().print(f);
    f << " , ";
    (p->first).getClade2().print(f);
    f << " ]" << endl;
  }
}

template<typename T> //with T instead of int
string Clade::randomTree(multimap<Clade,pair<Clade,T> >& mm,
			 map<Clade,Alias<dynamic_bitset<unsigned char> >* >& am,
			 map<pair<dynamic_bitset<unsigned char>,dynamic_bitset<unsigned char> >,double>& cladeLogProbMap,
			 mt19937_64& rng,
			 double& logTopologyProbability)
{
  if ( count()==1 ) { // one leaf
    stringstream ss;
    dynamic_bitset<unsigned char>::size_type first = clade.find_first(); //find first 1
    ss << first + 1;
    // back when subsets were stored right to left;
    //ss << size() - first;
    return ss.str();
  }
  if ( am.find(*this) == am.end() ) { // need to initialize alias for this clade
    pair< typename multimap<Clade,pair<Clade,T> >::iterator,typename multimap<Clade,pair<Clade,T> >::iterator> ret = mm.equal_range(*this);
    vector<double> probs;
    vector<dynamic_bitset<unsigned char> > indices;
    int index = 0;
    T total = 0;
    for ( typename multimap<Clade,pair<Clade,T> >::iterator p=ret.first; p!= ret.second; ++p, ++index ) {
      T counter = (p->second).second; //p->second: pair clade,int
      probs.push_back((double)counter);
      indices.push_back( (p->second).first.get() ); //get(): get dynamic bit set
      total += counter;
    }  
    for (int i = 0; i < probs.size(); ++i) {
      probs[i] = probs[i]/(double)total;
      cladeLogProbMap[ pair<dynamic_bitset<unsigned char>,dynamic_bitset<unsigned char> >(get(),indices[i]) ] = log(probs[i]);
    }
    am[*this] = new Alias<dynamic_bitset<unsigned char> >(probs,indices);
  }
  Clade c1( (am[*this])->pick(rng) );
  Clade c2(clade - c1.get());
  string s1 = c1.randomTree<T>(mm,am,cladeLogProbMap,rng,logTopologyProbability);
  string s2 = c2.randomTree<T>(mm,am,cladeLogProbMap,rng,logTopologyProbability);
  string out;
  if ( c1 > c2 )
    out = '(' + s1 + ',' + s2 + ')';
  else
    out = '(' + s2 + ',' + s1 + ')';
  logTopologyProbability += cladeLogProbMap[ pair<dynamic_bitset<unsigned char>,dynamic_bitset<unsigned char> >(get(),c1.get()) ];
  return out;
}

template<typename T>
string CCDProbs<T>::randomTree(mt19937_64& rng,double& logTopologyProbability)
{
  return all.randomTree<T>(mm,am,cladeLogProbMap,rng,logTopologyProbability) + ';';
}

#endif
