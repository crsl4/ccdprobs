// Bret Larget
// University of Wisconsin, Madison
//
// July 24, 2012
// February 20, 2013
//
// Copyright (C) 2012, 2013

// Use mbsum to summarize the tree file from a MrBayes analysis.
// Then use ccdprob on the output from mbsum to reestimate tree probabilities
//   using conditional clade probability distributions.


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
///////////////////////////////////////////////////////
#include "rand.h";
#include "alias.h"
///////////////////////////////////////////////////////////
using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////
string VERSION="1.0.5";

// Version history
// 1.0.1 --- 23 July 2012 --- added counts to the output map files
// 1.0.2 --- 24 July 2012 --- fixed bug when input trees have three subtrees at the root

////////////////////////////////////////////////////////////////////////////////////////////////////

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


class CladeInt;

class Clade {
private:
  dynamic_bitset<unsigned char> clade;
public:
  Clade() {}
  Clade(int n) { clade.resize(n); }
  ~Clade() { clade.clear(); }
  Clade(dynamic_bitset<unsigned char> x) { clade=x; }
  dynamic_bitset<unsigned char> get() const { return clade; }
  int count() const { return clade.count(); }
  int size() const { return clade.size(); }
  void clear() { clade.clear(); }
  void set(dynamic_bitset<unsigned char> x) { clade=x; }
  void resize(int n) { clade.resize(n); }
//  void add(int x) { clade[x-1] = 1; }  // original
  void add(int x) { clade[size() - x] = 1; } // taxa #x stored at size-x (1 is at size-1, numTaxa is at 0).
  void add(Clade c) { clade |= c.get(); }
//  void subtract(int x) { clade[x-1] = 0; }
  void print(ostream&) const ;
  friend bool operator< (const Clade&,const Clade&);
  friend bool operator> (const Clade&,const Clade&);
  friend bool operator== (const Clade&,const Clade&);
  friend Clade operator- (const Clade&,const Clade&);
  string randomTree(Rand&,multimap<Clade,CladeInt>&,map<Clade,Alias<dynamic_bitset<unsigned char> >* >&);
};

////////////////////////////////////////////////////////////////////////////////
// CladeInt stuff
//
// A clade int pair will be a clade and an int paired together.
class CladeInt {
private:
  Clade clade;
  int count;
public:
  CladeInt(Clade c,int x) {
    clade = c;
    count = x;
  }
  Clade getClade() const { return clade; }
  int getInt() const { return count; }
};

string Clade::randomTree(Rand &rand, multimap<Clade,CladeInt>& mm, map<Clade,Alias<dynamic_bitset<unsigned char> >* >& am)
{
  if ( count()==1 ) {
    stringstream ss;
    dynamic_bitset<unsigned char>::size_type first = clade.find_first();
    ss << size() - first;
    return ss.str();
  }
  if ( am.find(*this) == am.end() ) { // need to initialize alias for this clade
    pair<multimap<Clade,CladeInt>::iterator,multimap<Clade,CladeInt>::iterator> ret = mm.equal_range(*this);
    vector<double> probs;
    vector<dynamic_bitset<unsigned char> > indices;
    int index = 0;
    int total = 0;
    for ( multimap<Clade,CladeInt>::iterator p=ret.first; p!= ret.second; ++p, ++index ) {
      int counter = (p->second).getInt();
      probs.push_back((double)counter);
      indices.push_back( (p->second).getClade().get() );
      total += counter;
    }  
    for (int i = 0; i < probs.size(); ++i) {
      probs[i] = probs[i]/(double)total;
    }
    am[*this] = new Alias<dynamic_bitset<unsigned char> >(probs,indices);
  }
  Clade c1( (am[*this])->pick(rand) );
  Clade c2(clade - c1.get());
  string s1 = c1.randomTree(rand,mm,am);
  string s2 = c2.randomTree(rand,mm,am);
  string out;
  if ( c1 > c2 )
    out = '(' + s1 + ',' + s2 + ')';
  else
    out = '(' + s2 + ',' + s1 + ')';
  return out;
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

Clade operator- (const Clade& c1,const Clade& c2)
{
  return Clade(c1.get() - c2.get());
}

// ############################################################

class Tree {
private:
  string top;
  int numTaxa;
  vector<Clade> clades;
public:
  Tree(string,int);
  ~Tree() { 
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

Tree::Tree(string x,int n)
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

bool operator< (const CladePair& c1,const CladePair& c2)
{
  return ( c1.getClade1() < c2.getClade1() ? true : (c1.getClade1() == c2.getClade1() ? (c1.getClade2() < c2.getClade2()) : false ) );
}



// ************************************************************
// Rooted Binary Tree stuff
// ************************************************************

class Node {
private:
  Node* parent;
  Node* left;
  Node* right;
  bool leaf;
  bool root;
  Clade clade;
public:
  Node() { left = NULL; right = NULL; parent = NULL; leaf=false; root=false; }
  Node(Node* p) { left = NULL; right = NULL; parent = p; }
  Node* getParent() const { return parent; }
  Node* getLeft() const { return left; }
  Node* getRight() const { return right; }
  bool getLeaf() const { return leaf; }
  bool getRoot() const { return root; }
  Clade getClade() const { return clade; }
  void setParent(Node* p) { parent = p; }
  void setLeft(Node* p) { left = p; }
  void setRight(Node* p) { right = p; }
  void setLeaf(bool b) { leaf = b; }
  void setRoot(bool b) { root = b; }
  void addChild(Node* c) {
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
};

class RootedTree {
private:
  string top;
  string binaryTop;
  int numTaxa;
  Node* root;
  vector<Node*> nodes;
  vector<Clade> clades;
public:
  RootedTree() {};
  RootedTree(string,int);

  Node* getRoot() const { return root; }
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

// reads in a string of form "(subtree1,subtree2,subtree3);"
// and changes it to "((subtree1,subtree2),subtree3);"
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
    x = "((" + subtrees[0] + "," + subtrees[1] + ")," + subtrees[2] + ");";
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
  vector<Node*> w; // working vector of nodes
  while ( true ) {
    char c = s.peek();
    if ( c=='(' ) { // add new internal node
      Node* node = new Node();
      if ( first ) { // root of the tree
	node->setRoot(true);
	first = false;
      }
      else {
	Node* parent = w.back();
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
      Node* node = new Node(w.back());
      node->setLeaf(true);
      node->setNumTaxa(numTaxa);
      node->add(tax);
      w.back()->addChild(node);
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
  for ( vector<Node*>::iterator p=nodes.begin(); p != nodes.end(); p++ ) {
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

void RootedTree::count(int n,map<Clade,int>& cladeCount,map<CladePair,int>& pairCount)
{
  for ( vector<Node*>::iterator p = nodes.begin(); p != nodes.end(); p++ ) {
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

double RootedTree::estimateProbability(int sampleSize, map<Clade,int>& cladeCount, map<CladePair,int>& pairCount)
{
  double prob = 1.0;
  for ( vector<Node*>::iterator p = nodes.begin(); p != nodes.end(); p++ ) {
    if ( (*p)->getLeaf() )
      continue;
    Clade z = (*p)->getClade();
    Clade x = (*p)->getLeft()->getClade();
    Clade y = (*p)->getRight()->getClade();
    if ( x > y )
      x = y;
    CladePair c(z,x);
    prob *= pairCount[c] / (double)cladeCount[z];
  }
  return prob;
}

// ************************************************************

// Simple class to contain all input and output file streams

class Files
{
public:
  ifstream in; // BUCKy format file beginning with a translate table and followed by lines that consist of a tree, space, and a count.
  ofstream log; // output log file
  ofstream out; // output with table of trees, input sample probability, and clade estimated probability
  ofstream smap; // output with clade map for single clades
  ofstream tmap; // output with clade map for calde triples
  ifstream tin; // file with translate table and lines of trees for which the clade estimated probabilities will be computed; like BUCKy format, but no counts.
  ofstream tout; // output with extra trees and their clade estimated probabilities
////////////////////////////////////////////////////////////////////////////
  ofstream rtrees; // output with randomly generated trees
  ofstream bin; // rewriting input file, but trees written as binary
////////////////////////////////////////////////////////////////////////////
};

// ************************************************************

void usage()
{
/////////////////////////////////////////////////////////////////////////////////
  cerr << "Usage: ccdprobs <--out|-o> root-name-for-output-files [--extra] extra_tree_file_input <--in|-i> infile <--s1|-s1> seed1 <--s2|-s2> seed2 <--n|-n> number-of-random-trees" << endl << endl;
////////////////////////////////////////////////////////////////////////////////////
  cerr << "ccdprobs VERSION " << VERSION << endl << endl;
  cerr << "There is a required root name for all output files preceded by the flag --out or -o ." << endl;
  cerr << "There is a required input files preceded by the flag --in or -i ." << endl; 
/////////////////////////////////////////////////////////////////////////////////////
  cerr << "There is a required integer preceded by the flag --s1 or -s1 ." << endl;
  cerr << "There is a required integer preceded by the flag --s2 or -s2 ." << endl;
  cerr << "There is a required integer preceded by the flag --n or -n ." << endl;
//////////////////////////////////////////////////////////////////////////////////////
  cerr << "There is an optional extra file preceded by the flag --extra . " << endl;
  cerr << "These output files are always produced:" << endl << endl;
  cerr << "  .out   --- one line per sampled tree with the tree followed by four numbers:" << endl;
  cerr << "             (1) simple relative frequency of tree (SRF)" << endl;
  cerr << "             (2) conditional clade distribution (CCD) estimate of tree probability" << endl;
  cerr << "             (3) log10 of SRF" << endl;
  cerr << "             (4) log10 of CCD" << endl << endl;
  cerr << "  .smap  --- one line for each clade in some sampled tree preceded by its observed relative frequency" << endl;
  cerr << "  .tmap  --- one line for each triple of clades found in some sampled tree preceded by its observed relative frequency" << endl << endl;
  cerr << "  .log   --- log file that contains copy of the input command and the output to the terminal" << endl << endl;
  cerr << "This output file is only produced if the --extra flag is included with a file of extra trees" << endl << endl;
  cerr << "  .extra --- one line per input extra tree with the tree and its CCD estimated probability" << endl << endl;
//////////////////////////////////////////////////////////////////////////////////////
  cerr << "  .rtrees --- one line per tree randomly generated " << endl << endl;
///////////////////////////////////////////////////////////////////////////////////////
  cerr << "Note that the trees in the extra input file are not used to compute the maps for the CCD probability calculations" << endl << endl;
  cerr << "Example beginning with output from MrBayes:" << endl << endl;
  cerr << "   mbsum --skip 1001 --out example.in example.nex.run1.t" << endl;
  cerr << "   ccdprobs --out example --extra extras --in example.in" << endl << endl;
  cerr << "produces files:" << endl;
  cerr << "   example.out" << endl;
  cerr << "   example.smap" << endl;
  cerr << "   example.tmap" << endl;
  cerr << "   example.log" << endl;
  cerr << "   example.extra" << endl;
////////////////////////////////////////////////////////////////////////////
  cerr << "   example.rtrees" << endl << endl;
////////////////////////////////////////////////////////////////////////////
  exit(1);
}

void checkFlag(string key,int i,int argc) {
  if ( i==argc ) {
    cerr << "Error: flag " << key << " is not followed by required argument." << endl;
    usage();
  }
}
//////////////////////////////////////////////////////////////////////////////
void initialize(Files& f,unsigned int& seed1, unsigned int& seed2, unsigned int& randTreeNumbers, int argc,char* argv[])
///////////////////////////////////////////////////////////////////////////
{
  if ( argc==1 )
    usage();

  int i=1;
  string rootFile;
  bool infileGiven=false;
  bool rootFileGiven=false;
  bool extraFileGiven=false;
////////////////////////////////////////
  bool randTreeNumbersGiven=false;
////////////////////////////////////////
  while ( i < argc ) {
    string key=argv[i++];

    if ( key=="--help" || key=="=h" )
      usage();

    if ( key=="--version" ) {
      cerr << "ccdprob VERSION " << VERSION << endl;
      exit(2);
    }

    if ( key=="-i" || key=="--in" ) {
      checkFlag(key,i,argc);
      f.in.open(argv[i]);
      if ( f.in.fail() || f.in.bad() ) {
	cerr << "Error: cannot read file " << argv[i] << endl;
	exit(1);
      }
      infileGiven = true;
    }
    else if ( key=="--extra" ) {
      checkFlag(key,i,argc);
      f.tin.open(argv[i]);
      if ( f.tin.fail() || f.tin.bad() ) {
	cerr << "Error: cannot read file " << argv[i] << endl;
	exit(1);
      }
      extraFileGiven = true;
    }
    else if ( key=="-o" || key=="--out" ) {
      checkFlag(key,i,argc);
      rootFile = (string)(argv[i]);
      rootFileGiven = true;
    }    
//////////////////////////////////////////////////////////
    else if (key == "-s1" || key == "--seed1") {
      checkFlag(key,i,argc);
      string num = argv[i];
      istringstream f(num);
      if (!(f >> seed1)) {
        cerr << "Error in seed1: --seed1 integer " << endl;    
        exit(1);   
      }
      if (seed1 == 0) {
        cerr << "Warning: parameter seed1 must be at least one. Ignorning command --seed1 " << seed1 << "." << endl;
        exit(1);
      }      
    } 
    else if (key == "-s2" || key == "--seed2") {
      checkFlag(key,i,argc);
      string num = argv[i];
      istringstream f(num);
      if (!(f >> seed2)) {
        cerr << "Error in seed2: --seed2 integer " << endl;    
        exit(1);   
      }
      if (seed2 == 0) {
        cerr << "Warning: parameter seed2 must be at least one. Ignorning command --seed2 " << seed2 << "." << endl;
        exit(1);
      } 
    }  
    else if (key == "-n") {
      checkFlag(key,i,argc);
      string num = argv[i];
      istringstream f(num);
      if (!(f >> randTreeNumbers)) {
        cerr << "Error in number of trees: -n integer " << endl;    
        exit(1);   
      }
      if (randTreeNumbers == 0) {
        cerr << "Warning: parameter randTreeNumbers must be at least one. Ignorning command -n " << randTreeNumbers << "." << endl;
        exit(1);
      }       
      randTreeNumbersGiven=true;
    }   
////////////////////////////////////////////////////////////////
    else {
      cerr << "Error: Did not recognize argument flag " << key << endl;
      usage();
    }
    ++i;
  }
  if ( !infileGiven )
    cerr << "Error: no flag --in or -i followed by input file." << endl;
  if ( !rootFileGiven )
    cerr << "Error: no flag --out or -o followed by output root file name." << endl;
//////////////////////////////////////////////////////////////////////////
  if ( !randTreeNumbersGiven )
    cerr << "Error: no flag -n followed by number of random trees." << endl;  
  if ( !infileGiven || !rootFileGiven || !randTreeNumbersGiven )
    usage();
//////////////////////////////////////////////////////////////////////////////
  // Open required output files
  string outString = rootFile + ".out";
  f.out.open(outString.c_str());
  if ( f.out.fail() || f.out.bad() ) {
    cerr << "Error: cannot open file " << outString << " for output." << endl;
    exit(1);
  }

  string smapString = rootFile + ".smap";
  f.smap.open(smapString.c_str());
  if ( f.smap.fail() || f.smap.bad() ) {
    cerr << "Error: cannot open file " << smapString << " for single clade map output." << endl;
    exit(1);
  }

  string tmapString = rootFile + ".tmap";
  f.tmap.open(tmapString.c_str());
  if ( f.tmap.fail() || f.tmap.bad() ) {
    cerr << "Error: cannot open file " << tmapString << " for triple clade map output." << endl;
    exit(1);
  }

  string logString = rootFile + ".log";
  f.log.open(logString.c_str());
  if ( f.log.fail() || f.log.bad() ) {
    cerr << "Error: cannot open file " << logString << " for output log file." << endl;
    exit(1);
  }
///////////////////////////////////////////////////////////////////////////////////////
  string rtreesString = rootFile + ".rtrees";
  f.rtrees.open(rtreesString.c_str());
  if ( f.rtrees.fail() || f.rtrees.bad() ) {
    cerr << "Error: cannot open file " << rtreesString << " for output rtrees file." << endl;
    exit(1);
  }

  string binaryInputTreeString = rootFile + ".bintree";
  f.bin.open(binaryInputTreeString.c_str());
  if ( f.bin.fail() || f.bin.bad() ) {
    cerr << "Error: cannot open file " << binaryInputTreeString << " for output of input trees in binary format file." << endl;
    exit(1);
  }
///////////////////////////////////////////////////////////////////////////////////////
  if ( extraFileGiven ) {
    string extraString = rootFile + ".extra";
    f.tout.open(extraString.c_str());
    if ( f.tout.fail() || f.tout.bad() ) {
      cerr << "Error: cannot open file " << extraString << " for output file for extra trees." << endl;
      exit(1);
    }
  }

  // Write call to log file
  f.log << "ccdprobs version " << VERSION << endl << endl;
  time_t rawtime;
  time( &rawtime );
  
  f.log << "  produced at " << ctime( &rawtime ) << endl;
  f.log << "Call:";
  for ( int i=0; i<argc; i++ )
    f.log << " " << argv[i];
  f.log << endl << endl;


  return;
}

void writeTranslateTable(ofstream& f, vector<int>& taxaNumbers,vector<string>& taxaNames) {
  f << "translate" << endl;
  for ( int i=0; i < taxaNumbers.size(); ++i ) {
    f << setw(8) << taxaNumbers[i] << " " << taxaNames[i];
    if ( i < taxaNumbers.size()-1 )
      f << "," << endl;
    else
      f << ";" << endl;
  }
}

void readTranslateTable(ifstream& f,vector<int>& taxaNumbers,vector<string>& taxaNames,char& eol)
{
  // determine the eol character
  stringstream s;
  do {
    eol = f.get();
    s << eol;
  } while( (eol != '\n') && (eol != '\r') && !f.fail() );
  if(eol == '\r') {
    eol = f.peek();
    if(eol != '\n')
      eol = '\r';
  }
  // now use eol as EOL character
  string line;
  s >> line;
  if ( line != "translate" ) {
    cerr << "Error: expected first line to be \"translate\"; instead, it is " << line << endl;
    exit(1);
  }
  while ( getline(f,line,eol) ) {
    int n;
    string name;
    stringstream ss(line);
    ss >> n;
    taxaNumbers.push_back(n);
    // assume here no space in the name; also, last character is either , or ;
    char c;
    bool done=false;
    while ( ss >> c ) {
      if ( c==',' ) {
	taxaNames.push_back(name);
	break;
      }
      else if ( c==';' ) {
	taxaNames.push_back(name);
	done=true;
	break;
      }
      else
	name += c;
    }
    if ( done )
      break;
  }
}

// readTrees
//
// f is an input file stream that should be pointed at a line that contain a Newick tree followed by space and an integer count
// out is the file stream for output produced by the function
// eol is the end of line character previously determined
// numTaxa is the total number of taxa in each tree, previously determined from the translate table
// sampleSize is the total number of sampled trees, and will be the sum of the integer counts on each line
// cladeCount is a map from the unsigned integer representation of the clade to the count of the total number of trees with the clade
// pairCount is a map from a clade pair (A,B) that counts the number of trees with clade A at an internal node and clade B at a child node of the node with clade A
// 

void readTrees(ifstream& f,ofstream& out,char eol,int numTaxa,int& sampleSize,map<Clade,int>& cladeCount,map<CladePair,int>& pairCount,vector<int>& counts,vector<RootedTree>& trees)
{
  string line;
  while ( getline(f,line,eol) ) {
    stringstream s(line);
    string tree;
    int count;
    s >> tree >> count;
    sampleSize += count;

    RootedTree rt(tree,numTaxa);
    rt.count(count,cladeCount,pairCount);

    counts.push_back(count);
    trees.push_back(rt);
  }
}

void calculateExtraTrees(ifstream& tin, ofstream& tout, int sampleSize, map<Clade,int>& cladeCount, map<CladePair,int>& pairCount)
{
  // read translate table; for now, do not check to see if it matches
  vector<int> taxaNumbers;
  vector<string> taxaNames;
  char eol;
  readTranslateTable(tin,taxaNumbers,taxaNames,eol);
  int numTaxa = taxaNames.size();
    
  string line;
  while ( getline(tin,line,eol) ) {
    stringstream s(line);
    string tree;
    s >> tree;
    RootedTree rt(tree,numTaxa);

    tout.setf(ios::fixed,ios::floatfield);
    double p = rt.estimateProbability(sampleSize,cladeCount,pairCount);
    tout << tree << setw(11) << setprecision(8) << p << setw(13) << setprecision(8) << log10(p) << endl;
  }
  tin.close();
  tout.close();
}

void writeOutput(Files& f, vector<RootedTree>& trees, vector<int>& counts, int sampleSize, map<Clade,int>& cladeCount, map<CladePair,int>& pairCount, vector<int>& taxaNumbers, vector<string>& taxaNames, unsigned int randTreeNumbers, Clade& all, Rand &rand, multimap<Clade,CladeInt>& mm, map<Clade,Alias<dynamic_bitset<unsigned char> >* >& am)
{
  f.out.setf(ios::fixed,ios::floatfield);

  double total = 0.0,min=1.0,max=0.0;

  for ( int i=0; i< trees.size(); i++ ) {
    double p = counts[i] / (double)(sampleSize);
    double est = trees[i].estimateProbability(sampleSize,cladeCount,pairCount);
    total += est;
    if ( est > max )
      max = est;
    if( est < min )
      min = est;
    f.out << trees[i].getTop();
    f.out << setw(11) << setprecision(8) << p << setw(11) << setprecision(8) << est << setw(13) << setprecision(8) << log10(p) << setw(13) << setprecision(8) << log10(est) << endl;
  }

  cerr.setf(ios::fixed,ios::floatfield);
  cerr << endl << "sample size = " << sampleSize << " total trees." << endl;
  cerr << endl << "total   estimated probability = " << setw(10) << setprecision(8) << total << endl;
  cerr         << "maximum estimated probability = " << setw(10) << setprecision(8) << max << endl;
  cerr         << "minimum estimated probability = " << setw(10) << setprecision(8) << min << endl;

  f.log.setf(ios::fixed,ios::floatfield);
  f.log << endl << "sample size = " << sampleSize << " total trees." << endl;
  f.log << endl << "total   estimated probability = " << setw(10) << setprecision(8) << total << endl;
  f.log         << "maximum estimated probability = " << setw(10) << setprecision(8) << max << endl;
  f.log         << "minimum estimated probability = " << setw(10) << setprecision(8) << min << endl;

  writeTranslateTable(f.smap,taxaNumbers,taxaNames);
  f.smap.setf(ios::fixed,ios::floatfield);
  for ( map<Clade,int>::iterator p=cladeCount.begin(); p != cladeCount.end(); p++ ) {
    f.smap << setw(10) << setprecision(8) << p->second / (double) sampleSize << " ";
    f.smap << setw(10) << p->second << " ";
    (p->first).print(f.smap);
    f.smap << endl;
  }

  writeTranslateTable(f.tmap,taxaNumbers,taxaNames);
  f.tmap.setf(ios::fixed,ios::floatfield);
  for ( map<CladePair,int>::iterator p=pairCount.begin(); p != pairCount.end(); p++ ) {
    f.tmap << setw(10) << setprecision(8) << p->second / (double) sampleSize << " ";
    f.tmap << setw(10) << p->second << " ";
    f.tmap << "[ ";
    (p->first).getClade1().print(f.tmap);
    f.tmap << " , ";
    (p->first).getClade2().print(f.tmap);
    f.tmap << " ]" << endl;
  }

  writeTranslateTable(f.bin,taxaNumbers,taxaNames);
  int index=0;
  for ( vector<RootedTree>::iterator p=trees.begin(); p!= trees.end(); ++p )
    f.bin << (*p).getBinaryTop() << " " << counts[index++] << endl;

  for ( int k=0; k< randTreeNumbers; ++k ) {
    f.rtrees << all.randomTree(rand,mm,am) << ';' << endl;
  }
  
}

int main(int argc, char* argv[])
{
  Files f;
///////////////////////////////////////////////////////////////////////////////
  unsigned int seed1, seed2;
  unsigned int randTreeNumbers;
  initialize(f,seed1, seed2, randTreeNumbers, argc,argv);

///////////////////////////////////////////////////////////////////////////////
  vector<int> taxaNumbers;
  vector<string> taxaNames;
  char eol;
  readTranslateTable(f.in,taxaNumbers,taxaNames,eol);
  int numTaxa = taxaNames.size();

  map<Clade,int> cladeCount;
  map<CladePair,int> pairCount;

  int sampleSize=0; // sum of counts
  vector<int> counts; // integer count for each tree
  vector<RootedTree> trees;

  readTrees(f.in,f.out,eol,numTaxa,sampleSize,cladeCount,pairCount,counts,trees);

////////////////////////////////////////////////////////////////////////////////
  // Produce a random tree
  // Open the random number generator.
/////////////////////////////////////////////////////////////////////////////

  Rand rand; //default  
  if ( seed1 != 0 )
    rand.setSeed1(seed1);
  if ( seed2 != 0)
    rand.setSeed2(seed2);

////////////////////////////////////////////////////////////////////////////////

  multimap<Clade,CladeInt> mm;
  for ( map<CladePair,int>::iterator p=pairCount.begin(); p!=pairCount.end(); ++p ) {
    Clade parent=(p->first).getClade1();
    Clade child=(p->first).getClade2();
    mm.insert( pair<Clade,CladeInt>(parent,CladeInt(child, p->second)) );
  }

  Clade all(numTaxa);
  for ( int k=1; k<=numTaxa; ++k )
    all.add(k);

  // calculate for extra trees if given
//  if ( f.tin.is_open() )
//    calculateExtraTrees(f.tin,f.tout,sampleSize,cladeCount,pairCount);
/////////////////////////////////////////////////////////////////////////////////////
  // Write output files
  map<Clade,Alias<dynamic_bitset<unsigned char> >* > am;
  writeOutput(f,trees,counts,sampleSize,cladeCount,pairCount,taxaNumbers,taxaNames, randTreeNumbers, all, rand, mm, am);
////////////////////////////////////////////////////////////////////////////////////////
  return 0;
}

