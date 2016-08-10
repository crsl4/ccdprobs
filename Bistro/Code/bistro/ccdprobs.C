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
      Clade clad = node->getClade();
      cout << "Clade for read internal node: ";
      clad.print(cout);
      cout << endl;
    }
    else if ( c==')' ) { // complete a subtree
      nodes.push_back(w.back());
      nodes.back()->copyClade(nodes.back()->getLeft()->getClade()); //using clades already in RootedNode
      nodes.back()->add(nodes.back()->getRight()->getClade());
      w.pop_back();
      clades.push_back(nodes.back()->getClade());
      s >> c;
    }
    else if ( isdigit(c) ) { // add leaf node
      int tax;
      s >> tax;
      RootedNode* node = new RootedNode(w.back()); //why this?
      node->setLeaf(true);
      node->setNumTaxa(numTaxa);
      Clade clad = node->getClade();
      cout << "Clade for read leaf node: " << tax << " ";
      clad.print(cout);
      cout << endl;
      node->add(tax); //this is the problem!
      clad = node->getClade();
      cout << "Clade for read leaf node: " << tax << " ";
      clad.print(cout);
      cout << endl;
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
  printClades(cout);
  cout << endl;
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

/*
Tree CCDProbs::rtop(mt19937_64& rng)
{
  
}
*/
