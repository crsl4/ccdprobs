#include <iostream>
#include <iomanip>
#include <vector>
#include <list>
#include <fstream>
#include <sstream>
#include <string>
#include <locale> // isdigit(), tolower()
#include <cctype>
#include <random>

//#include "alias.h"
#include "ccdprobs.h"

using namespace std;

#include <Boost/dynamic_bitset.hpp>

using namespace boost;

bool operator< (const Clade& c1,const Clade& c2)
{
  return c1.get() < c2.get();
}

bool operator> (const Clade& c1,const Clade& c2)
{
  return c1.get() > c2.get();
}

bool operator== (const Clade& c1,const Clade& c2)
{
  return c1.get() == c2.get();
}

bool operator< (const CladePair& c1,const CladePair& c2)
{
  return ( c1.getClade1() < c2.getClade1() ? true : (c1.getClade1() == c2.getClade1() ? (c1.getClade2() < c2.getClade2()) : false ) );
}

void Clade::print(ostream& f) const {
    
  string comma = "";

  dynamic_bitset<unsigned char>::size_type first = clade.find_first();
  dynamic_bitset<unsigned char>::size_type last = first;
  dynamic_bitset<unsigned char>::size_type i = first;

  f << '{';
  while ( i != dynamic_bitset<unsigned char>::npos ) {
    i = clade.find_next(i);
    if ( i > last+1 ) { // print the previous range and set first to i
      f << comma << first+1;
      comma = ",";
      if ( last > first )
	f << '-' << last+1;
      first = i;
    }
    last = i;
  }
  f << '}';
}

CCDTree::CCDTree(string x,int n)
{
  top = x;
  numTaxa = n;
  stringstream s(x);
  vector<Clade*> w;
  while ( true ) {
    char c = s.peek();
    if ( c=='(' ) { // begin a new clade
      w.push_back( new Clade(numTaxa) );
      s >> c;
    }
    else if ( c==')' ) { // complete a clade
      clades.push_back(w.back()->get());
      w.pop_back();
      s >> c;
    }
    else if ( isdigit(c) ) { // add taxon to each working clade
      int tax;
      s >> tax;
      for ( vector<Clade*>::iterator p=w.begin(); p!= w.end(); p++ )
	(*p)->add(tax);
    }
    else if ( c==',' ) {
      s >> c;
      continue;
    }
    else if ( c==';' ) {
      s >> c;
      break;
    }
  }
}

void makeBinary(string& x)
{
  stringstream s(x);
  vector<string> subtrees;
  char c;
  // read opening '('
  s >> c;
  int level = 1;
  string subtree;
  while ( level > 0 && !s.fail() ) {
    s >> c;
    if ( c == ',' && level==1 ) {
      subtrees.push_back(subtree);
      subtree.clear();
      continue;
    }
    else if ( c == '(' )
      level += 1;
    else if ( c == ')' ) {
      level -= 1;
      if ( level == 0 ) {
	subtrees.push_back(subtree);
	subtree.clear();
      }
    }
    if ( level > 0 )
      subtree += c;
  }

  if ( subtrees.size() == 2 )
    return;
  if ( subtrees.size() == 3 ) {
    x = "(" + subtrees[0] + ",(" + subtrees[1] + "," + subtrees[2] + "));";
    return;
  }
  else {
    cerr << "Error: input tree " << x << " has more than three subtrees or only one subtree at the root level." << endl;
    exit(1);
  }
}

// In version 1.0, code assumed that the outermost level of the tree was binary,
// as was the case in older versions of mbsum.  However, the distributed version
// of mbsum has three subtrees at the top level.  This induces a bug when the code
// expected this top level to be binary.
//
// This version of ccdprobs corrects this and assumes that the top level of the tree has three children.
// It assumes that the third subtree is the same for all trees and it adds a node to combine the first two subtrees.
RootedTree::RootedTree(string x,int n)
{
  top = x;
  makeBinary(x);
  binaryTop = x;
  numTaxa = n;
  stringstream s(x);
  bool first=true;
  vector<RootedNode*> w; // working vector of nodes
  while ( true ) {
    char c = s.peek();
    if ( c=='(' ) { // add new internal node
      RootedNode* node = new RootedNode();
      if ( first ) { // root of the tree
	node->setRoot(true);
	first = false;
      }
      else {
	RootedNode* parent = w.back();
	parent->addChild(node);
	node->setParent(parent);
      }
      w.push_back( node );
      s >> c;
    }
    else if ( c==')' ) { // complete a subtree
      nodes.push_back(w.back());
      nodes.back()->copyClade(nodes.back()->getLeft()->getClade());
      nodes.back()->add(nodes.back()->getRight()->getClade());
      w.pop_back();
      clades.push_back(nodes.back()->getClade());
      s >> c;
    }
    else if ( isdigit(c) ) { // add leaf node
      int tax;
      s >> tax;
      RootedNode* node = new RootedNode(w.back());
      node->setLeaf(true);
      node->setNumTaxa(numTaxa);
      node->add(tax);
      w.back()->addChild(node);
      nodes.push_back(node);
    }
    else if ( c==',' ) {
      s >> c;
      continue;
    }
    else if ( c==';' ) {
      s >> c;
      break;
    }
  }
}

void RootedTree::print(ostream& f)
{
  f << top << endl;
  for ( vector<RootedNode*>::iterator p=nodes.begin(); p != nodes.end(); p++ ) {
    (*p)->getClade().print(f);
    if ( !((*p)->getLeaf()) ) {
      f << " :";
      (*p)->getLeft()->getClade().print(f);
      f << " ";
      (*p)->getRight()->getClade().print(f);
    }
    f << endl;
  }
  f << endl;
}

void RootedTree::printClades(ostream& f)
{
  for ( vector<Clade>::iterator p = clades.begin(); p != clades.end(); ++p) {
    (*p).print(f); 
    f << " ";
  }
}    

// to add n to count of clade
template<typename T> //and change int for T
void RootedTree::count(T n,map<Clade,T>& cladeCount,map<CladePair,T>& pairCount)
{
  for ( vector<RootedNode*>::iterator p = nodes.begin(); p != nodes.end(); p++ ) {
    if ( (*p)->getLeaf() )
      continue;
    Clade z = (*p)->getClade();
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
    rt.count<T>(m->second,cladeCount,pairCount);
    sampleSize += m->second;
  }

  for ( map<CladePair,int>::iterator p=pairCount.begin(); p!=pairCount.end(); ++p ) {
    Clade parent=(p->first).getClade1();
    Clade child=(p->first).getClade2();
    mm.insert( pair<Clade,pair<Clade,int> >(parent,make_pair(child, p->second)) );
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
  for ( map<Clade,int>::iterator p=cladeCount.begin(); p != cladeCount.end(); p++ ) {
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
  for ( map<CladePair,int>::iterator p=pairCount.begin(); p != pairCount.end(); p++ ) {
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
    ss << size() - first;
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
  return all.randomTree(mm,am,cladeLogProbMap,rng,logTopologyProbability) + ';';
}

/*
Tree CCDProbs::rtop(mt19937_64& rng)
{
  
}
*/
