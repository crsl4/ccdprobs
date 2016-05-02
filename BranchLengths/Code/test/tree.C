#include <iostream>
#include <iomanip>
#include <vector>
#include <list>
#include <fstream>
#include <sstream>
#include <string>
#include <locale> // isdigit(), tolower()
#include <cctype>
#include <algorithm> // sort()
#include <random>

#include "Eigen/Core"
#include "Eigen/Eigenvalues"

#include "tree.h"

using namespace std;
using namespace Eigen;

void Node::print(ostream& f)
{
  f << setw(3) << number << ": " << name << ": " << level << ": ";
  f << "Edges:";
  for ( vector<Edge*>::iterator e=edges.begin(); e!=edges.end(); ++e )
  {
    if ( e != edges.begin() )
      f << ",";
    f << " e[" << (*e)->getNumber() << "] -> " << "n[" << getNeighbor(*e)->getNumber() << "]";
  }
  f << endl;
}

void Edge::print(ostream& f)
{
  f << setw(3) << number << ":";
  f << " n[" << (nodes[0])->getNumber() << "] <-> " << "n[" << (nodes[1])->getNumber() << "]";
  f << " length = " << setprecision(8) << length << endl;
}

void Tree::print(ostream& f)
{
  f << makeTopologyNumbers() << endl;
//  f << treeString << endl;
  f << "numTaxa = " << numTaxa << ", numNodes = " << numNodes << ", numEdges = " << numEdges << endl;
  f << "Nodes:" << endl;
  for(int i=0;i<numNodes;i++)
    nodes[i]->print(f);
  f << "Edges:" << endl;
  for(int i=0;i<numEdges;i++)
    edges[i]->print(f);
  f << endl;
}

////////////////////////////////////////////////////////////////////////////////
//
// Functions to read in a tree
//
////////////////////////////////////////////////////////////////////////////////

char readOneTreeCharacter(istringstream& s) {
  char c = s.peek();
  if ( s.eof() || s.fail() )
  {
    cerr << "Error: Tree ends prematurely." << endl;
    exit(1);
  }
  s >> c;
  return c;
}

string readName(istringstream& s) {
  string name;
  while ( !s.fail() )
  {
    char c = s.peek();
    if ( (c==',') || (c==':') || (c==')') || (c==';') )
      break;
    s >> c;
    if ( s.fail() )
    {
      cerr << "Error while reading name." << endl;
      cerr << "Failed when reading " << c << endl;
      cerr << "name = " << name << endl;
      exit(1);
    }
    name += c;
  }
  return(name);
}

double readNumber(istringstream& s)
{
  double x;
  char c = s.peek();
  if ( s.fail() || !isdigit(c) )
  {
    cerr << "Error while reading number." << endl;
    exit(1);
  }
  s >> x;
  return(x);
}

void printRest(istringstream& s)
{
  string a;
  s >> a;
  cerr << "Remainder is /" << a << "/." << endl;
}

void expectedCharacter(char expected,char actual)
{
  cerr << "Error: expected '" << expected << "', but read '" << actual << "'." << endl;
}

void checkTreeEnd(istringstream& s) {
  char c = readOneTreeCharacter(s);
  if ( c != ';' )
  {
    expectedCharacter(';',c);
    printRest(s);
    exit(1);
  }
}

void Tree::readSubtree(istringstream& s,Node* parent,vector<Node*>& leafNodes,vector<Node*>& internalNodes)
{
  Node* n = new Node(); // root of this subtree
  Edge* e = new Edge(); // edge connecting n to parent
  edges.push_back(e);

  e->setNodes(n,parent);
  n->addEdge(e);
  parent->addEdge(e);

  char c = s.peek();

  if ( c=='(' )
  {
    n->setLeaf(false);
    internalNodes.push_back(n);
    s >> c; // read the '('
    readSubtree(s,n,leafNodes,internalNodes); // first subtree must be followed by a comma
    c = readOneTreeCharacter(s);
    if(c != ',') {
      expectedCharacter(',',c);
      printRest(s);
      exit(1);
    }
    while ( true )
    {
      readSubtree(s,n,leafNodes,internalNodes);
      c = readOneTreeCharacter(s);
      if(c==',')
        continue;
      else
        break;
    }
    if ( c != ')' )
    {
      expectedCharacter(')',c);
      printRest(s);
      exit(1);
    }
  }
  else
  {
    string name = readName(s);
    leafNodes.push_back(n);
    n->setLeaf(true);
    n->setName(name);
    numTaxa++;
  }
  return;
}

// This function numbers leaves in the order they are read in.
// Names are treated as strings.
// The function relabelTree() will change node numbers and names to match the sequences.

Tree::Tree(string line)
{
  // Create the tree from parenthetic representation in line.
  // Create new nodes and edges on the fly.
  // Leaf node and internal node pointers placed in separate vectors.  Combine to vector nodes at end.

  numTaxa = 0;
  treeString = line;

  istringstream s(line);

  char c = s.peek();
  if(c != '(') { // trivial case where tree is a single node
    nodes.push_back(new Node(0,true));
    numTaxa++;
    nodes[0]->setName(readName(s));
    root = nodes[0];
    checkTreeEnd(s);
    return;
  }

  // string begins with '(', so tree is not trivial
  vector<Node*> leafNodes;
  vector<Node*> internalNodes;

  internalNodes.push_back(new Node(false));
  root = internalNodes[0];
  s >> c; // read in left parenthesis

  while ( true ) {
    readSubtree(s,root,leafNodes,internalNodes);
    c = readOneTreeCharacter(s);
    if ( c==',' ) // root has another subtree
      continue;
    else if( c==')' ) { // end of base subtree
      checkTreeEnd(s);
      break;
    }
    else { // error
      cerr << "Error: expected ',' or ')', but read ";
      printRest(s);
      exit(1);
    }
  }
  // Now, load up the nodes vector and set node numbers.
  // leaf nodes are numbered in the order read in, which may not match the name
  for ( vector<Node*>::iterator p=leafNodes.begin(); p!=leafNodes.end(); ++p )
    nodes.push_back(*p);
  for ( vector<Node*>::reverse_iterator p=internalNodes.rbegin(); p!=internalNodes.rend(); ++p )
    nodes.push_back(*p);
  for ( int i=0; i<nodes.size(); ++i )
    nodes[i]->setNumber(i+1);
  for ( int i=0; i<edges.size(); ++i )
    edges[i]->setNumber(i+1);

  numNodes = nodes.size();
  numEdges = edges.size();

  setNodeLevels();

  leafNodes.resize(0);
  internalNodes.resize(0);
}

// The tree constructor numbers leafs in the order they appear with the topology
// and uses the names as written, either names or numbers.
// This function will relabel and renumber the nodes so that the names/numbers match the input sequence file.
// It also reorders the vector nodes in tree so that the order matches the numbers.
void Tree::relabel(Alignment& alignment)
{
  for ( vector<Node*>::iterator p=nodes.begin(); p!=nodes.end(); ++p )
  {
    if ( !(*p)->getLeaf() )
      continue;

    string name = (*p)->getName();
    if ( isdigit( name[0] ) )
    {
      stringstream s(name);
      int number;
      s >> number;
      map<int,string>::iterator m=alignment.numberToNameMap.find(number);
      if ( m==alignment.numberToNameMap.end() )
      {
        cerr << "Error: could not find " << number << " in map." << endl;
        exit(1);
      }
      (*p)->setName(m->second);
      (*p)->setNumber(m->first);
    }
    else
    {
      map<string,int>::iterator m=alignment.nameToNumberMap.find(name);
      if ( m==alignment.nameToNumberMap.end() )
      {
        cerr << "Error: could not find " << name << " in map." << endl;
        exit(1);
      }
      (*p)->setNumber(m->second);

    }
  }
  // Now, sort the nodes
  vector<pair<int,Node*> > v;
  for ( vector<Node*>::iterator p=nodes.begin(); p!=nodes.end(); ++p )
  {
    if ( !(*p)->getLeaf() )
      continue;
    v.push_back( pair<int,Node*> ((*p)->getNumber(),*p) );
  }
  sort(v.begin(),v.end());

  vector<pair<int,Node*> >::iterator p=v.begin();
  for ( int k=0; k<numNodes; ++k )
  {
    if ( !(nodes[k])->getLeaf() )
      continue;
    nodes[k] = p->second;
    p++;
  }
}

void Node::makeTopology(stringstream& s,Edge* parent,bool useName)
{
  if ( leaf )
  {
    if ( useName )
      s << name;
    else
      s << number;
    return;
  }
  // internal node
  s << '(';
  bool first=true;
  for ( vector<Edge*>::iterator e=edges.begin(); e !=edges.end(); ++e )
  {
    if ( (*e)==parent )
      continue;
    if ( first )
      first = false;
    else
      s << ',';
    getNeighbor(*e)->makeTopology(s,*e,useName);
  }
  s << ')';
}

string Tree::makeTopology(bool useName)
{
  stringstream s;
  root->makeTopology(s,NULL,useName);
  s << ';';
  string top;
  s >> top;
  return top;
}

void Tree::randomBranchLengths(mt19937_64& rng,double lambda)
{
  exponential_distribution<double> rexp(lambda);
  for ( vector<Edge*>::iterator e=edges.begin(); e!=edges.end(); ++e )
  {
    (*e)->setLength( rexp(rng) );
  }
}

Vector4d translate(char base)
{
  Vector4d condProbs (0,0,0,0);
  switch ( base )
  {
  case 'a' :
    condProbs(0) = 1;
    break;
  case 'c' :
    condProbs(1) = 1;
    break;
  case 'g' :
    condProbs(2) = 1;
    break;
  case 't' :
  case 'u' :
    condProbs(3) = 1;
    break;
  case 'r' :
    condProbs(0) = condProbs(2) = 1;
    break;
  case 'y' :
    condProbs(1) = condProbs(3) = 1;
    break;
  case 'k' :
    condProbs(2) = condProbs(3) = 1;
    break;
  case 'm' :
    condProbs(0) = condProbs(1) = 1;
    break;
  case 's' :
    condProbs(1) = condProbs(2) = 1;
    break;
  case 'w' :
    condProbs(0) = condProbs(4) = 1;
    break;
  case 'b' :
    condProbs(1) = condProbs(2) = condProbs(3) = 1;
    break;
  case 'v' :
    condProbs(0) = condProbs(1) = condProbs(2) = 1;
    break;
  case 'h' :
    condProbs(0) = condProbs(1) = condProbs(3) = 1;
    break;
  case 'd' :
    condProbs(0) = condProbs(2) = condProbs(3) = 1;
    break;
  case 'x' :
  case 'n' :
  case '-' :
    condProbs(0) = condProbs(1) = condProbs(2) = condProbs(3) = 1;
    break;
  }
  return condProbs;
}

pair<double,Vector4d> Node::getProb() // pattern must be set
{
  map<string,pair<double,Vector4d> >::iterator m=patternToProbMap.find(pattern);
  if ( m != patternToProbMap.end() )
    return m->second;
  else
  {
    cerr << "Error: node " << number << " has no probability for pattern " << pattern << endl;
    exit(1);
  }
}

void Node::calculateEdges(QMatrix& qmatrix)
{
  for ( vector<Edge*>::iterator e=edges.begin(); e != edges.end(); ++e )
  {
    if ( *e != edgeParent )
      (*e)->calculate(qmatrix);
  }
}

void Node::calculate(int site,const Alignment& alignment,Edge* parent,bool recurse)
{
  pattern.clear();
  map<string,pair<double,Vector4d> >::iterator m;
  if ( leaf )
  {
    char base = alignment.getBase(number,site);
    base = tolower(base);
    pattern = base;
    m = patternToProbMap.find(pattern);
    if ( m == patternToProbMap.end() ) // first time this pattern calculated
    {
      patternToProbMap[ pattern ] = pair<double,Vector4d> (0,translate(base)); //clau: translate->4 probs
    }
    return;
  }
  // internal node
  for ( vector<Edge*>::iterator e=edges.begin(); e!=edges.end(); ++e )
  {
    if ( (*e) != parent )
    {
      if ( recurse )
        getNeighbor(*e)->calculate(site,alignment,*e,recurse);
      pattern += getNeighbor(*e)->getPattern(); //clau: pattern of int node is the concatenation of pattern of children
    }
  }
  m = patternToProbMap.find(pattern);
  if ( m == patternToProbMap.end() ) // first time this pattern calculated
  {
    double scale=0;
    Vector4d tempProb(1,1,1,1);

    for ( vector<Edge*>::iterator e=edges.begin(); e!=edges.end(); ++e )
    {
      if ( (*e) != parent )
      {
        pair<double,Vector4d> condProbPair = getNeighbor(*e)->getProb();
        scale += condProbPair.first;
        Vector4d partProb = (*e)->getTransitionMatrix() * condProbPair.second;
        for ( int i=0; i<4; ++i )
          tempProb[i] *= partProb(i);
      }
    }
    double vmax = tempProb.maxCoeff();
    scale += log( vmax );
    tempProb /= vmax;
    patternToProbMap[ pattern ] = pair<double,Vector4d> (scale,tempProb);
  }
}

// must have accurate patternToProbMaps prior to calling
// call clearProbMaps() on either whole tree or nodes for which
// this may have changed due to rerooting or branch length changes
// this function computes transition probabilities for all edges
// this could be changed to only do so for edges for which this is not accurate
double Tree::calculate(const Alignment& alignment,QMatrix& qmatrix)
{
  // initialize vector to store site log likelihoods
  VectorXd logLikelihood;
  logLikelihood.resize(alignment.getNumSites());
  // calculate transition probabilities for all edges
  for ( vector<Edge*>::iterator e=edges.begin(); e!=edges.end(); ++e )
    (*e)->calculate(qmatrix);

  // loop over sites and calculate log-likelihood while traversing the tree
  for ( int k=0; k < alignment.getNumSites(); ++k )
  {
    root->calculate(k,alignment,NULL,true);
    pair<double,Vector4d> condProbPair = root->getProb();
    logLikelihood(k) = condProbPair.first + log( qmatrix.getStationaryP().dot(condProbPair.second) );
//    cout << logLikelihood(k) << endl;
  }
  return logLikelihood.sum();
}

void Node::clearProbMaps(Edge* parent)
{
  for ( vector<Edge*>::iterator e=edges.begin(); e!=edges.end(); ++e )
  {
    if ( (*e) != parent )
      getNeighbor(*e)->clearProbMaps(*e);
  }
  patternToProbMap.clear();
}

void Tree::clearProbMaps()
{
  root->clearProbMaps(NULL);
}

// put edges in random order
void Node::randomize(mt19937_64& rng,Edge* parent)
{
  if ( leaf )
    return;
  for ( vector<Edge*>::iterator e=edges.begin(); e!=edges.end(); ++e )
  {
    if ( (*e) != parent )
      getNeighbor(*e)->randomize(rng,*e);
  }
  int numEdges = edges.size();
  for ( int i=0; i<numEdges-1; ++i )
  {
    uniform_int_distribution<> rint(i,numEdges-1);
    int j = rint(rng);
    if ( j != i )
    {
      Edge* temp = edges[i];
      edges[i] = edges[j];
      edges[j] = temp;
    }
  }
}

// pick a random internal node for the root
// randomly reorder all children of all nodes
// depends on node order with all leaves first and internal nodes later
void Tree::randomize(mt19937_64& rng)
{
  uniform_int_distribution<> rint(numTaxa,numNodes-1);
  root = nodes[ rint(rng) ];
  root->randomize(rng,NULL);
  setNodeLevels();
}

void Tree::partialPathCalculations(double t,Alignment& alignment,Node* na,Edge* ea,Node* nb,Edge* eb,QMatrix& qmatrix,double& logl,double& dlogl,double& ddlogl,bool recurse)
{
  Matrix4d P = qmatrix.getTransitionMatrix(t);
  Matrix4d QP = qmatrix.getQP(t);
  Matrix4d QQP = qmatrix.getQQP(t);
  int numSites = alignment.getNumSites();

  logl = 0;
  dlogl = 0;
  ddlogl = 0;
  for ( int k=0; k<numSites; ++k )
  {
    na->calculate(k,alignment,ea,recurse); // set pattern and put probability in map if not already there
    nb->calculate(k,alignment,eb,recurse);
    pair<double,Vector4d> pa = na->patternToProbMap[na->getPattern()];
    pair<double,Vector4d> pb = nb->patternToProbMap[nb->getPattern()];
    Vector4d va = pa.second;
    Vector4d vq = qmatrix.getStationaryP();
    for ( int i=0; i<4; ++i )
      va(i) *= vq(i);
    Vector4d vb = pb.second;
    double f0 = (va.asDiagonal() * P * vb.asDiagonal()).sum();
    double f1 = (va.asDiagonal() * QP * vb.asDiagonal()).sum();
    double f2 = (va.asDiagonal() * QQP * vb.asDiagonal()).sum();
    logl += pa.first + pb.first + log( f0 );
    dlogl += f1/f0;
    ddlogl += (f0*f2 - f1*f1)/(f0*f0);
  }
}

double Tree::pathLogLikelihood(double t,Alignment& alignment,Node* na,Edge* ea,Node* nb,Edge* eb,QMatrix& qmatrix,bool recurse)
{
  double logl,dlogl,ddlogl;
  partialPathCalculations(t,alignment,na,ea,nb,eb,qmatrix,logl,dlogl,ddlogl,recurse);
  return logl;
}

double Tree::pathLogLikelihoodDerivative(double t,Alignment& alignment,Node* na,Edge* ea,Node* nb,Edge* eb,QMatrix& qmatrix,bool recurse)
{
  double logl,dlogl,ddlogl;
  partialPathCalculations(t,alignment,na,ea,nb,eb,qmatrix,logl,dlogl,ddlogl,recurse);
  return dlogl;

}

double Tree::pathLogLikelihoodSecondDerivative(double t,Alignment& alignment,Node* na,Edge* ea,Node* nb,Edge* eb,QMatrix& qmatrix, bool recurse)
{
  double logl,dlogl,ddlogl;
  partialPathCalculations(t,alignment,na,ea,nb,eb,qmatrix,logl,dlogl,ddlogl,recurse);
  return ddlogl;
}

void mleError(Node* na,Node* nb,double curr,double prop,double curr_dlogl,double prop_dlogl)
{
  cerr << "Error: too many iterations in mleDistance." << endl;
  cerr << "Nodes " << na->getNumber() << " and " << nb->getNumber() << endl;
  cerr << "Derivative = " << curr_dlogl << " and " << prop_dlogl << endl;
  exit(1);
}
// Find mle distance from node a to b through a path that uses edges ea and eb.
// Conditon on data in subtrees through other edges.
// Assumes that edge lengths in these subtrees exist
// and that the patternToProbMaps are accurate if edges ea and eb head toward the root.
double Tree::mleDistance(Alignment& alignment,Node* na,Edge* ea,Node* nb,Edge* eb,QMatrix& qmatrix)
{
  bool recurse=true;
  int iter=0;
  double curr = 0.05;
  // get a decent starting point
  double curr_logl,curr_dlogl,curr_ddlogl;
  partialPathCalculations(curr,alignment,na,ea,nb,eb,qmatrix,curr_logl,curr_dlogl,curr_ddlogl,true);
  double prop = curr;
  double prop_logl = curr_logl;
  double prop_dlogl = curr_dlogl;
  double prop_ddlogl = curr_ddlogl;
  if ( curr_dlogl > 0 )
  {
    do
    {
      curr = prop;
      curr_logl = prop_logl;
      curr_dlogl = prop_dlogl;
      curr_ddlogl = prop_ddlogl;
      prop = 2*curr;
      partialPathCalculations(prop,alignment,na,ea,nb,eb,qmatrix,prop_logl,prop_dlogl,prop_ddlogl,recurse);
      if ( ++iter > 100 )
        mleError(na,nb,curr,prop,curr_dlogl,prop_dlogl);
    } while ( prop_dlogl > 0);
  }
  else
  {
    do
    {
      curr = prop;
      curr_logl = prop_logl;
      curr_dlogl = prop_dlogl;
      curr_ddlogl = prop_ddlogl;
      prop = 0.5*curr;
      partialPathCalculations(prop,alignment,na,ea,nb,eb,qmatrix,prop_logl,prop_dlogl,prop_ddlogl,recurse);
      if ( ++iter > 100 )
        mleError(na,nb,curr,prop,curr_dlogl,prop_dlogl);
    } while ( prop_dlogl < 0 );
  }
  // switch to protected Newton-Raphson
  prop = curr - curr_dlogl * (prop - curr) / (prop_dlogl - curr_dlogl);
  do
  {
    curr = prop;
    partialPathCalculations(curr,alignment,na,ea,nb,eb,qmatrix,curr_logl,curr_dlogl,curr_ddlogl,recurse);
    if ( ++iter > 100 )
      mleError(na,nb,curr,prop,curr_dlogl,prop_dlogl);
    double delta = -curr_dlogl / curr_ddlogl;
    prop = curr + delta;
    while ( prop < 0 )
    {
      delta = 0.5*delta;
      prop = curr + delta;
    }
    partialPathCalculations(prop,alignment,na,ea,nb,eb,qmatrix,prop_logl,prop_dlogl,prop_ddlogl,recurse);
    if ( ++iter > 100 )
      mleError(na,nb,curr,prop,curr_dlogl,prop_dlogl);
    while ( ( fabs(prop_dlogl) > fabs(curr_dlogl) ) && fabs(curr - prop) >1.0e-8 )
    {
      delta = 0.5*delta;
      prop = curr + delta;
      partialPathCalculations(prop,alignment,na,ea,nb,eb,qmatrix,prop_logl,prop_dlogl,prop_ddlogl,recurse);
      if ( ++iter > 100 )
        mleError(na,nb,curr,prop,curr_dlogl,prop_dlogl);
    }
  } while ( fabs(curr - prop) > 1.0e-8);
  return prop;
}

void Node::setLevel(Edge* parent)
{
  if ( parent==NULL )
    level = 0;
  else
    level = getNeighbor(parent)->getLevel() + 1;
  if ( leaf )
    return;
  for ( vector<Edge*>::iterator e=edges.begin(); e != edges.end(); ++e )
    if ( (*e) != parent )
      getNeighbor(*e)->setLevel(*e);
}

void Tree::setNodeLevels()
{
  root->setLevel(NULL);
}

void Node::depthFirstNodeList(list<Node*>& nodeList,Edge* parent)
{
  list<Node*>::iterator p=nodeList.begin();
  while ( p != nodeList.end() && level <= (*p)->getLevel() )
    ++p;
  nodeList.insert(p,this);
  for ( vector<Edge*>::iterator e=edges.begin(); e!=edges.end(); ++e )
  {
    if ( *e != parent )
      getNeighbor(*e)->depthFirstNodeList(nodeList,*e);
  }
}

void Tree::depthFirstNodeList(list<Node*>& nodeList)
{
  root->depthFirstNodeList(nodeList,NULL);
}

// clau: set attribute nodeParent, edgeParent for Node
void Node::setActiveChildrenAndNodeParents(Edge* parent)
{
  if ( parent != NULL )
  {
    nodeParent = getNeighbor(parent);
    edgeParent = parent;
  }
  else
  {
    nodeParent = NULL;
    edgeParent = NULL;
  }
  if ( leaf )
    return;
  activeChildren.clear(); // clau: activeChildren is vector of nodes of Node
  for ( vector<Edge*>::iterator e=edges.begin(); e!= edges.end(); ++e ) //clau: edges is vector of edges of Node
  {
    if ( (*e) != parent ) // clau: go to all children
    {
      Node* n = getNeighbor(*e); //clau: *e = pointer to e
      activeChildren.push_back(n); //clau: make n an active child
      n->setActiveChildrenAndNodeParents(*e); //clau: n->bla ~ n.bla (conceptually)
    }
  }
}

void Tree::setActiveChildrenAndNodeParents()
{
  root->setActiveChildrenAndNodeParents(NULL);
}

Node* Node::closeRelative()
{
  int i=0;
  Node* relative;
  do
  {
    relative = nodeParent->getActiveChild(i++);
  } while ( relative == this && i < nodeParent->getActiveChildrenSize() );
  if ( relative == this )
  {
    cerr << "Error: closeRelative failed." << endl;
    exit(1);
  }
  if ( relative->getActiveChildrenSize() == 0 )
    return relative;
  return relative->getActiveChild(0);
}

pair<int,int> getPair(int x,int y)
{
  if ( x > y )
  {
    int temp;
    temp=x; x=y; y=temp;
  }
  return pair<int,int> (x,y);
}

void Tree::generateBranchLengths(Alignment& alignment,QMatrix& qmatrix)
{
  map<pair<int,int>,double> distanceMap;
  list<Node*> nodeList;
  depthFirstNodeList(nodeList);
  setActiveChildrenAndNodeParents();
  cout << "Node List:";
  for ( list<Node*>::iterator p=nodeList.begin(); p!= nodeList.end(); ++p )
    cout << " " << (*p)->getNumber();
  cout << endl;

  list<Node*>::iterator p=nodeList.begin();
  while ( true )
  {
    Node* x;
    Node* y;
    Node* z;
    Node* par; //clau: parent of x,y,z
    if ( (*p)->getNodeParent() == root ) //clau: if p parent is root
    {
      if ( root->getActiveChildrenSize() != 3) //clau: need root to have 3 children
      {
        cerr << "yeah, write the general code...." << root->getActiveChildrenSize() << endl;
        cerr << root->getActiveChild(0)->getNumber() << endl;
        cerr << (*p)->getNumber() << endl;
        exit(1);
      }
      par = root;
      x = *p++;
      y = *p++;
      z = *p;
    }
    else //clau: p parent not root
    {
      // cout << "*p is " << (*p)->getNumber();
      // cout << endl;
      x = *p++; //question: dont understand here, seems like *p==*p++, I think it means take *p and move
      y = *p++;
      par = x->getNodeParent();
      z = par->closeRelative();
      // cout << "*p is " << (*p)->getNumber();
      // cout << endl;
    }
    // cout << "x is " << x->getNumber();
    // cout << ", y is " << y->getNumber();
    // cout << ", z is " << z->getNumber();
    // cout << ", par is " << par->getNumber();
    // cout << endl;
    // do calculations
    x->calculateEdges(qmatrix); //question: why do we do this? aren't bl meaningless at this point?
    y->calculateEdges(qmatrix);
    z->calculateEdges(qmatrix);
    for ( int k=0; k<alignment.getNumSites(); ++k )
    {
      x->calculate(k,alignment,x->getEdgeParent(),false);
      y->calculate(k,alignment,y->getEdgeParent(),false);
      z->calculate(k,alignment,z->getEdgeParent(),false);
    } //clau: after this, we have in each node the patternProb map, but we dont know how many different patterns, do we? (question)
    // aqui voy: to do, finish studying generateBL function
//    cerr << x->getNumber() << " " << y->getNumber() << " " << z->getNumber() << " " << par->getNumber() << endl;
    map<pair<int,int>,double>::iterator m;
    double dxy, dxz, dyz;
    m = distanceMap.find( getPair(x->getNumber(),y->getNumber()) );
    if ( m == distanceMap.end() )
      dxy = mleDistance(alignment,x,x->getEdgeParent(),y,y->getEdgeParent(),qmatrix);
    else
      dxy = m->second;
    m = distanceMap.find( getPair(x->getNumber(),z->getNumber()) );
    if ( m == distanceMap.end() )
      dxz = mleDistance(alignment,x,x->getEdgeParent(),z,z->getEdgeParent(),qmatrix);
    else
      dxz = m->second;
    m = distanceMap.find( getPair(y->getNumber(),z->getNumber()) );
    if ( m == distanceMap.end() )
      dyz = mleDistance(alignment,y,y->getEdgeParent(),z,z->getEdgeParent(),qmatrix);
    else
      dyz = m->second;
    double lengthX = (dxy + dxz - dyz)*0.5;
    if ( lengthX < 0 )
    {
      cerr << "Warning: fix negative edge length." << endl;
      lengthX = 0;
    }
    double lengthY = (dxy + dyz - dxz)*0.5;
    if ( lengthY < 0 )
    {
      cerr << "Warning: fix negative edge length." << endl;
      lengthY = 0;
    }
    x->getEdgeParent()->setLength( lengthX );
    y->getEdgeParent()->setLength( lengthY );
    double lengthZ = (dxz + dyz - dxy)*0.5;
    if ( lengthZ < 0 )
    {
      cerr << "Warning: fix negative edge length." << endl;
      lengthZ = 0;
    }
    if ( par==root )
    {
      z->getEdgeParent()->setLength( lengthZ );
      break;
    }
    distanceMap[ getPair(z->getNumber(),par->getNumber()) ] = lengthZ;
    par->deactivateChild(1);
    par->deactivateChild(0);
  }

  for ( map<pair<int,int>,double>::iterator m=distanceMap.begin(); m!= distanceMap.end(); ++m )
    cout << (*m).first.first << " " << (*m).first.second << " --> " << (*m).second << endl;
}
