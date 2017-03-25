#ifndef __TREE_H
#define __TREE_H

#define MIN_EDGE_LENGTH 0.00001
#define MAX_EDGE_LENGTH 10
#define TOL 1.0e-6
#define VERBOSE false
#define PRIOR_MEAN 0.1
#define LAMBDA 0.19

#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <list>
#include <string>
#include <random>
#include <algorithm>

#include "Eigen/Core"
#include "Eigen/Eigenvalues"

#include "sequence.h"
#include "model.h"
#include "ccdprobs.h"

using namespace std;
using namespace Eigen;

class CladeGraph;

class Clade;

class Node;

class Edge
{
private:
  int number;
  int current; // 0 or 1 for index into which length is current for MCMC purposes
  double length[2];
  Node* nodes[2];
  // things for likelihood calculation
  Matrix4d transitionMatrix;
public:
  Edge()
  {
    length[0] = 0;
    length[1] = 0;
    current = 0;
  }
  ~Edge()
  {
    nodes[0] = nodes[1] = NULL;
  }
  Edge(int n,double x) : number(n)
  {
    current = 0;
    length[0] = length[1] = x;
  }
  int getNumber() { return number; }
  void setNumber(int x) { number=x; }
  int getCurrent() { return current; }
  double getLength() { return length[current]; }
  void setLength(double x) { length[current]=x; }
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
  void calculate(QMatrix& qmatrix)
  {
    transitionMatrix = qmatrix.getTransitionMatrix(length[current]);
  }
  Matrix4d getTransitionMatrix() { return transitionMatrix; }
  void mleError(bool&);
  void calculate(double,Alignment&,QMatrix&,double&,double&,double&);
  double mleLength(Alignment&,QMatrix&,bool&);
  void randomLength(Alignment&,QMatrix&,mt19937_64&,double&,bool);
  bool isTerminal();
  void switchCurrent() { current = 1 - current; }
};

class Node
{
private:
  int number;
  int current; // 0 or 1
  string name;
  vector<Edge*> edges;
  bool leaf;
  // things for likelihood calculation
  string pattern;
  int level;
  vector<Node*> activeChildren;
  Node* nodeParent;
  Edge* edgeParent;
  Edge* mapParent; // pointer to parent edge used when creating patternToProbMap, Clau: I think we need this separate from parent
  int minNumber; // smallest number in subtree rooted at node
public:
  map<string,pair<double,Vector4d> > patternToProbMap[2];
  Node() { number = -1; level = 0; current = 0; }
  ~Node()
  {
    edges.clear();
    activeChildren.clear();
  }
  Node(int n,bool l) : number(n), leaf(l)
  {
    level = 0;
    current = 0;
  }
  Node(bool l) : leaf(l)
  {
    level = 0;
    current = 0;
  }
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
  void deleteEdge(Edge* e)
  {
    vector<Edge*>::iterator p = find(edges.begin(),edges.end(),e);
    if ( p != edges.end() )
      edges.erase(p);
  }
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
  void makeTreeNumbers(stringstream&,Edge*);
  string getPattern() const { return pattern; }
  int getLevel() const { return level; }
  void setLevel(Edge*);
  void setPattern(int,const Alignment& ,Edge* );
  void calculateAfterPattern(int,const Alignment& ,Edge*);
  void calculate(int,const Alignment&,Edge*);//,bool);
  pair<double,Vector4d> getProb();
  void randomize(mt19937_64&,Edge*);
  void clearProbMaps(Edge*);
  void clearProbMapsSmart(Edge*);
  void depthFirstNodeList(list<Node*>&,Edge*);
  void postorderCherryNodeList(list<Node*>&,Edge*);
  bool isPrunedLeaf(Edge*);
  Node* getNodeParent() const { return nodeParent; }
  Node* getActiveChild(int i) const { return activeChildren[i]; }
  int getActiveChildrenSize() const { return activeChildren.size(); }
  void deactivateChild(int i) { activeChildren.erase(activeChildren.begin() + i); }
  void setActiveChildrenAndNodeParents(Edge*);
  void setMapParent(Edge*);
  Edge* getMapParent() { return mapParent; }
  Node* closeRelative();
  Edge* getEdgeParent() const { return edgeParent; }
  void calculateEdges(QMatrix&);
  int getMinNumber() const { return minNumber; }
  int setMinNumber(Edge*);
  void sortCanonical(Edge*);
  void randomEdges(Alignment&,QMatrix&,mt19937_64&,Edge*,double&,bool);
  void parsimonyScore(Alignment&,Edge*,int,int&,int&);
  void clearMapParent();
  void distance(map<dynamic_bitset<unsigned char>,double>&,Clade&,Edge*);
  void processTree(CladeGraph*,dynamic_bitset<unsigned char>&,Edge* parent);
};

class Tree
{
private:
  string treeString;
  int numTaxa;
  vector<Node*> nodes;
  vector<Edge*> edges;
  Node* root;
  map<dynamic_bitset<unsigned char>,double> cladeToLengthMap;
public:
  Tree(string);
  Tree(string,Alignment);
  Tree(MatrixXd&);
  ~Tree() {
    for(vector<Edge*>::iterator e=edges.begin();e!=edges.end();e++)
      delete *e;
    for(vector<Node*>::iterator n=nodes.begin();n!=nodes.end();n++)
      delete *n;
    nodes.clear();
    edges.clear();
  }
  void treeInit(string);
  int getNumTaxa() { return numTaxa; }
  int getNumNodes() { return nodes.size(); }
  int getNumEdges() { return edges.size(); }
  Node* getNode(int n) { return nodes[n]; }
  Edge* getEdge(int n) { return edges[n]; }
  string getTreeString() { return treeString; }
  vector<Edge*> getEdges() {return edges; }
  Node* getRoot() { return root; }
  void readSubtree(istringstream&,Node*,int);
  void setNumbers(Node*);
  void print(ostream&);
  void readSubtree(istringstream&,Node*,vector<Node*>&,vector<Node*>&);
  void relabel(Alignment&);
  map<dynamic_bitset<unsigned char>,double> getCladeToLengthMap() {return cladeToLengthMap;}
  void setCladeToLengthMap( map<dynamic_bitset<unsigned char>,double> c ){ cladeToLengthMap = c;}
  string makeTopology(bool);
  string makeTopologyNames() { return makeTopology(true); }
  string makeTopologyNumbers() { return makeTopology(false); }
  string makeTreeNumbers();
  void randomBranchLengths(mt19937_64&,double);
  double calculate(const Alignment&,QMatrix&);
  void randomize(mt19937_64&);
  void randomizeBL(mt19937_64&);
  void clearProbMaps();
  void clearProbMapsSmart();
  /* void partialPathCalculations(double,Alignment&,Node*,Edge*,Node*,Edge*,QMatrix&,double&,double&,double&,bool); */
  /* double pathLogLikelihood(double,Alignment&,Node*,Edge*,Node*,Edge*,QMatrix&,bool); */
  /* double pathLogLikelihoodDerivative(double,Alignment&,Node*,Edge*,Node*,Edge*,QMatrix&,bool); */
  /* double pathLogLikelihoodSecondDerivative(double,Alignment&,Node*,Edge*,Node*,Edge*,QMatrix&,bool); */
  double mleDistance(Alignment&,Node*,Edge*,Node*,Edge*,QMatrix&,double);
  void setNodeLevels();
  void depthFirstNodeList(list<Node*>&);
  void postorderCherryNodeList(list<Node*>&);
  void setActiveChildrenAndNodeParents();
  void generateBranchLengths(Alignment&,QMatrix&,mt19937_64&,double&, bool, double,bool);
  void reroot(int);
  void setMinNumber() { root->setMinNumber(NULL); }
  void sortCanonical()
  {
    root->setMinNumber(NULL);
    root->sortCanonical(NULL);
  }
  void setNJDistances(MatrixXd&,mt19937_64&);
  void unroot();
  void randomEdges(Alignment&,QMatrix&,mt19937_64&,double&,bool);
  double logPriorExp(double);
  Vector2d mleLength2D(Alignment&,Node*,Edge*,Node*,Edge*,Node*,Edge*,QMatrix&,bool&);
  void partialPathCalculations2D(Vector2d ,double,Alignment& ,Node* ,Edge* ,Node* ,Edge* ,Node* ,Edge* ,QMatrix& ,double& ,Vector2d& ,Matrix2d&);// ,bool);
  Vector3d mleLength3D(Alignment&,Node*,Edge*,Node*,Edge*,Node*,Edge*,QMatrix&,bool&);
  void partialPathCalculations3D(Vector3d ,Alignment& ,Node* ,Edge* ,Node* ,Edge* ,Node* ,Edge* ,QMatrix& ,double& ,Vector3d& ,Matrix3d&);// ,bool);
  double vectorProduct(vector<Vector4d>);
  double vectorProduct4D(Vector4d, Vector4d, Vector4d, Vector4d);
  void mleError(bool&);
  void makeBinary();
  int parsimonyScore(Alignment&);
  double logLikelihoodScore(Alignment&, QMatrix&);
  void clearMapParent();
  Edge* whichMaxBranch(); //finds edge with max length
  void mcmc(QMatrix&, Alignment&, unsigned int, double, mt19937_64&, bool);
  void mcmc(QMatrix&, Alignment&, unsigned int, double, mt19937_64&, ofstream&, ofstream&, bool);
  void mcmc(QMatrix&, Alignment&, unsigned int, double, mt19937_64&, ofstream&, ofstream&, bool, bool);
  void setInitialEdgeLengths(double); // set all edge lengths to x
  double distance(Tree*);
  void processTree(CladeGraph*);
  void mcmcUpdateQ(int,MCMCStats&,QMatrix&,Alignment&,double,mt19937_64&);
  void mcmcUpdatePi(int,MCMCStats&,QMatrix&,Alignment&,double,mt19937_64&);
  void mcmcUpdateS(int,MCMCStats&,QMatrix&,Alignment&,double,mt19937_64&);
  void mcmcUpdateEdges(int,MCMCStats&,QMatrix&,Alignment&,mt19937_64&);
};

class MCMCStats
{
private:
  double currLogLikelihood;
  double sumAcceptP;
  double sumAcceptS;
  double sumAcceptBL;
  Vector4d avgP;
  Vector4d avgPold;
  Vector4d sP;
  VectorXd avgS;
  VectorXd avgSold;
  VectorXd sS;
  VectorXd avgBL;
public:
  MCMCStats(int,double);
  double getCurrLogLikelihood() { return currLogLikelihood; }
  void setCurrLogLikelihood(double x) { currLogLikelihood = x; }
  double getSumAcceptP() { return sumAcceptP; }
  void addSumAcceptP(double x) { sumAcceptP += x; }
  double getSumAcceptS() { return sumAcceptS; }
  void addSumAcceptS(double x) { sumAcceptS += x; }
  double getSumAcceptBL() { return sumAcceptBL; }
  void addSumAcceptBL(double x) { sumAcceptBL += x; }
  Vector4d getAvgP() { return avgP; }
  void setAvgP(Vector4d x) { avgP = x; }
  double getAvgP(int i) { return avgP(i); }
  double getAvgPold(int i) { return avgPold(i); }
  Vector4d getAvgPold() { return avgPold; }
  void getAvgPold(Vector4d x) { avgPold = x; }
  void saveAvgP() { avgPold = avgP; }
  Vector4d getSP() { return sP; }
  double getSP(int i) { return sP(i); }
  void addSP(Vector4d x) { sP += x; }
  VectorXd getAvgS() { return avgS; }
  void setAvgS(VectorXd x) { avgS = x; }
  double getAvgS(int i) { return avgS(i); }
  double getAvgSold(int i) { return avgSold(i); }
  VectorXd getAvgSold() { return avgSold; }
  void saveAvgS() { avgSold = avgS; }
  VectorXd getSS() { return sS; }
  double getSS(int i) { return sS(i); }
  void addSS(VectorXd x) { sS += x; }
  double getAvgBL(int i) { return avgBL(i); }
  void addAvgBL(int i,double x) { avgBL(i) += x; }
  void printMCMCSummary(ostream&,QMatrix&,int,unsigned int);
};

#endif
