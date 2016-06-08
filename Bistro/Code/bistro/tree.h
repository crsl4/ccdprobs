#ifndef __TREE_H
#define __TREE_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <list>
#include <string>
#include <random>

#include "Eigen/Core"
#include "Eigen/Eigenvalues"

#include "sequence.h"
#include "model.h"

using namespace std;
using namespace Eigen;

class Node;

class Edge
{
private:
  int number;
  double length;
  Node* nodes[2];
  // things for likelihood calculation
  Matrix4d transitionMatrix;
public:
  Edge() {}
  ~Edge() {
    nodes[0] = nodes[1] = NULL;
  }
  Edge(int n,double x) : number(n), length(x) {}
  int getNumber() { return number; }
  void setNumber(int x) { number=x; }
  double getLength() { return length; }
  void setLength(double x) { length=x; }
  Node* getOtherNode(Node *n) { return (nodes[0] == n ? nodes[1] : nodes[0]); }
  void setNode(Node *n,int i) { nodes[i]=n; }
  void setNodes(Node *m,Node *n) { nodes[0]=m; nodes[1]=n; }
  Node* getNode(int n) { return nodes[n]; }
  void swapNodes() {
    Node *temp=nodes[0];
    nodes[0] = nodes[1];
    nodes[1] = temp;
  }
  void print(ostream&);
  void calculate(QMatrix& qmatrix) { transitionMatrix = qmatrix.getTransitionMatrix(length); }
  Matrix4d getTransitionMatrix() { return transitionMatrix; }
};

class Node
{
private:
  int number;
  string name;
  vector<Edge*> edges;
  bool leaf;
  // things for likelihood calculation
  string pattern;
  int level;
  vector<Node*> activeChildren;
  Node* nodeParent;
  Edge* edgeParent;
  int minNumber; // smallest number in subtree rooted at node
public:
  map<string,pair<double,Vector4d> > patternToProbMap;
  Node() { number = -1; level = 0; }
  ~Node() { edges.clear(); }
  Node(int n,bool l) : number(n), leaf(l) { level = 0;}
  Node(bool l) : leaf(l) { level = 0; }
  int getNumber() const { return number; }
  void setNumber(int x) { number = x; }
  string getName() const { return name; }
  void setName(string n) { name = n; }
  bool getLeaf() { return leaf; }
  void setLeaf(bool x) { leaf = x; }
  int getNumEdges() { return edges.size(); }
  Edge* getEdge(int n) { return edges[n]; }
  void addEdge(Edge* e) { edges.push_back(e); }
  void setEdge(Edge *e,int n) { edges[n] = e; }
  void clearEdges() { edges.clear(); }
  Node* getNeighbor(Edge* e) { return e->getOtherNode(this); }
  Node* getNeighbor(int i) { return edges[i]->getOtherNode(this); }
  Edge* findEdge(Node* n) {
    Edge* e=NULL;
    for ( vector<Edge*>::iterator p=edges.begin(); p!=edges.end(); ++p )
      if ( (*p)->getOtherNode(this) == n ) {
	e = *p;
	break;
      }
    return e;
  }
  void print(ostream&);
  void makeTopology(stringstream&,Edge*,bool);
  string getPattern() const { return pattern; }
  int getLevel() const { return level; }
  void setLevel(Edge*);
  void calculate(int,const Alignment&,Edge*,bool);
  pair<double,Vector4d> getProb();
  void randomize(mt19937_64&,Edge*);
  void clearProbMaps(Edge*);
  void depthFirstNodeList(list<Node*>&,Edge*);
  Node* getNodeParent() const { return nodeParent; }
  Node* getActiveChild(int i) const { return activeChildren[i]; }
  int getActiveChildrenSize() const { return activeChildren.size(); }
  void deactivateChild(int i) { activeChildren.erase(activeChildren.begin() + i); }
  void setActiveChildrenAndNodeParents(Edge*);
  Node* closeRelative();
  Edge* getEdgeParent() const { return edgeParent; }
  void calculateEdges(QMatrix&);
  int getMinNumber() const { return minNumber; }
  int setMinNumber(Edge*);
  void sortCanonical(Edge*);
};

class Tree
{
private:
  string treeString;
  int numTaxa;
  vector<Node*> nodes;
  vector<Edge*> edges;
  Node* root;
public:
  Tree(string);
  Tree(MatrixXd&);
  ~Tree() {
    for(vector<Edge*>::iterator e=edges.begin();e!=edges.end();e++)
      delete *e;
    for(vector<Node*>::iterator n=nodes.begin();n!=nodes.end();n++)
      delete *n;
    nodes.clear();
    edges.clear();
  }
  int getNumTaxa() { return numTaxa; }
  int getNumNodes() { return nodes.size(); }
  int getNumEdges() { return edges.size(); }
  Node* getNode(int n) { return nodes[n]; }
  Edge* getEdge(int n) { return edges[n]; }
  string getTreeString() { return treeString; }
  void readSubtree(istringstream&,Node*,int);
  void setNumbers(Node*);
  void print(ostream&);
  void readSubtree(istringstream&,Node*,vector<Node*>&,vector<Node*>&);
  void relabel(Alignment&);
  string makeTopology(bool);
  string makeTopologyNames() { return makeTopology(true); }
  string makeTopologyNumbers() { return makeTopology(false); }
  void randomBranchLengths(mt19937_64&,double);
  double calculate(const Alignment&,QMatrix&);
  void randomize(mt19937_64&);
  void clearProbMaps();
  void partialPathCalculations(double,Alignment&,Node*,Edge*,Node*,Edge*,QMatrix&,double&,double&,double&,bool);
  double pathLogLikelihood(double,Alignment&,Node*,Edge*,Node*,Edge*,QMatrix&,bool);
  double pathLogLikelihoodDerivative(double,Alignment&,Node*,Edge*,Node*,Edge*,QMatrix&,bool);
  double pathLogLikelihoodSecondDerivative(double,Alignment&,Node*,Edge*,Node*,Edge*,QMatrix&,bool);
  double mleDistance(Alignment&,Node*,Edge*,Node*,Edge*,QMatrix&);
  void setNodeLevels();
  void depthFirstNodeList(list<Node*>&);
  void setActiveChildrenAndNodeParents();
  void generateBranchLengths(Alignment&,QMatrix&);
  void reroot(int);
  void setMinNumber() { root->setMinNumber(NULL); }
  void sortCanonical()
  {
    root->setMinNumber(NULL);
    root->sortCanonical(NULL);
  }
};

#endif
