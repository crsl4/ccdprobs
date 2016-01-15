#ifndef TREE_H_
#define TREE_H_

#include <iostream>
#include <iomanip>
#include <vector>
//#include <set>
#include <string>
#include <sstream>
#include <fstream>
#include <cctype>
//#include <map>
//#include <cmath>
//#include <ctime>

#include "distance.h"

using namespace std;

class Node;

class Edge {
private:
  int number;
  double length;
  Node* nodes[2];
public:
  Edge() {}
  ~Edge() {
    for(int i=0;i<2;i++)
      nodes[i] = NULL;
  }
  Edge(int n,double x) : number(n), length(x) {}
  int getNumber() { return number; }
  void setNumber(int x) { number=x; }
  void setLength(double x) { length=x; }
  double getLength() { return length; }
  Node* getOtherNode(Node *n) { return (nodes[0] == n ? nodes[1] : nodes[0]); }
  void setNode(Node *n,int i) { nodes[i]=n; }
  void setNodes(Node *m,Node *n) { nodes[0]=m; nodes[1]=n; }
  Node* getNode(int n) { return nodes[n]; }
  void swapNodes() {
    Node *temp=nodes[0];
    nodes[0] = nodes[1];
    nodes[1] = temp;
  }
};

class Node {
private:
  int number;
  bool leaf;
  vector<Edge*> edges;
public:
  Node() { number = -1; }
  ~Node() { edges.clear(); }
  Node(int n,bool l) : number(n), leaf(l) {}
  int getNumber() const { return number; }
  void setNumber(int x) { number = x; }
  bool isLeaf() { return leaf; }
  void setLeaf(bool x) { leaf = x; }
  int getNumEdges() { return edges.size(); }
  Edge* getEdge(int n) { return edges[n]; }
  void addEdge(Edge* e) { edges.push_back(e); }
  void setEdge(Edge *e,int n) { edges[n] = e; }
  void clearEdges() { edges.clear(); }
  void setNumbers(int&,Edge*);
  void setDistances(DistanceMatrix&,vector<double>&,vector<int>&,Edge*);
};

class Tree {
private:
  int numTaxa;
  int numNodes;
  int numEdges;
  vector<Node*> nodes;
  vector<Edge*> edges;
  Node* root;
public:
  Tree(string,int);
  ~Tree() {
    for(vector<Edge*>::iterator e=edges.begin();e!=edges.end();e++)
      delete *e;
    for(vector<Node*>::iterator n=nodes.begin();n!=nodes.end();n++)
      delete *n;
    nodes.clear();
    edges.clear();
  }
  int getNumTaxa() { return numTaxa; }
  int getNumNodes() { return numNodes; }
  int getNumEdges() { return numEdges; }
  Node* getNode(int n) { return nodes[n]; }
  Edge* getEdge(int n) { return edges[n]; }
  void readSubtree(istringstream&,Node*,int);
  void setNumbers(Node*);
  void setDistances(DistanceMatrix&);
};

void Node::setNumbers(int& nextNumber,Edge* parentEdge) {
  // Set numbers for internal nodes for subtrees along non-parent edges.
  // Numbers for leaves should already be set.
  // Set number of parent edge to match number of the node

  if ( !leaf ) {
    for ( vector<Edge*>::iterator e=edges.begin(); e!=edges.end(); e++ )
      if ( *e != parentEdge )
	(*e)->getOtherNode(this)->setNumbers(nextNumber,*e);
    setNumber(nextNumber++);
  }
  parentEdge->setNumber(number);
}

void Tree::setNumbers(Node* root) {
  // Internal nodes are not numbered when tree is read.
  // After reading the tree, we can give numbers to the internal nodes.
  int nextNumber = numTaxa + 1;
  for ( int i=0; i<root->getNumEdges(); ++i) {
    Edge* e = root->getEdge(i);
    e->getOtherNode(root)->setNumbers(nextNumber,e);
  }
  root->setNumber(nextNumber);
}

void advance(istringstream& s,char& c,int lineNumber) {
  c = s.peek();
  if ( s.eof() || s.fail() ) {
    cerr << "Error: Tree on line " << lineNumber << " ends prematurely." << endl;
    string a;
    s >> a;
    cerr << "Remainder of line is /" << a << "/." << endl;
    exit(1);
  }
  s >> c;
}

void checkCharacter(char c,char key,istringstream& s,int lineNumber) {
  if ( c != key ) {
    cerr << "Error: Tree does not match expected format." << endl;
    cerr << "Read " << c << " when expecting to read " << key << " on line " << lineNumber << endl;
    string a;
    s >> a;
    cerr << "Remainder of line is /" << a << "/." << endl;
    exit(1);
  }
}

void Tree::readSubtree(istringstream& s,Node* parent,int lineNumber)
{
  Node* n;
  char c = s.peek();
  if(isdigit(c)) {
    int num;
    s >> num;
    n = new Node();
    n->setNumber(num);
    n->setLeaf(true);
    numTaxa++;
  }
  else if(c=='(') {
    n = new Node();
    n->setLeaf(false);
    s >> c;
    readSubtree(s,n,lineNumber);
    advance(s,c,lineNumber);
    checkCharacter(c,',',s,lineNumber);
    while ( c == ',' ) {
      readSubtree(s,n,lineNumber);
      advance(s,c,lineNumber);
    }
    checkCharacter(c,')',s,lineNumber);
  }
  else {
    cerr << "Error: Expected beginning of subtree, but read '" << c << "' on line " << lineNumber << "." << endl;
    string a;
    s >> a;
    cerr << "Remainder of line is /" << a << "/." << endl;
    exit(1);
  }
  Edge* e = new Edge();
  e->setNodes(n,parent);
  n->addEdge(e);
  parent->addEdge(e);
  nodes.push_back(n);
  numNodes++;
  edges.push_back(e);
  numEdges++;

  c = s.peek();
  double length = 1.0;
  if ( c == ':' ) {
    s >> c;
    s >> length;
    if(s.fail()) {
      cerr << "Error: Tree on line " << lineNumber << ".  Expected an edge length after reading colon." << endl;
      string a;
      s >> a;
      cerr << "Remainder of line is /" << a << "/." << endl;
      exit(1);
    }
  }
  e->setLength(length);
}

Tree::Tree(string line,int lineNumber)
{
  // read in the tree from parenthetic representation.
  // Create new nodes and edges on the fly.
  // Then renumber.
  // Then, reorder so that leaves come first.

  numEdges = 0;
  numTaxa = 0;
  numNodes = 0;

  istringstream s(line);
  char c = s.peek();
  if(c =='(') {
    Node *n = new Node(); // root of the tree
    // root is not a leaf
    n->setLeaf(false);
    nodes.push_back(n);
    numNodes++;
    s >> c; // read in left parenthesis
    readSubtree(s,n,lineNumber);
    advance(s,c,lineNumber);
    checkCharacter(c,',',s,lineNumber);
    while ( c == ',' ) {
      readSubtree(s,n,lineNumber);
      advance(s,c,lineNumber);
    }
    checkCharacter(c,')',s,lineNumber);
    advance(s,c,lineNumber);
    checkCharacter(c,';',s,lineNumber);
  }
  else if(isdigit(c)) { // base tree may just be a single node
    int n;
    s >> n;
    nodes.push_back(new Node()); // root of the tree, thinks it is a leaf....
    numNodes = 1;
    nodes[0]->setLeaf(true);
    nodes[0]->setNumber(n);
    numTaxa++;
    advance(s,c,lineNumber);
    checkCharacter(c,';',s,lineNumber);
  }
  root = nodes[0];
  // give numbers to internal nodes, and picks root to be last taxon's neighbor
  setNumbers(root);
}

#endif
