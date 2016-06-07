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
#include "random.h"

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
  f << " n[" << (nodes[0])->getNumber() << "] <-> " << "n[" << (nodes[1])->getNumber() << "] ,";
  f << " length = " << setprecision(8) << length << endl;
}

void Edge::printTable(ostream& f)
{
  f << (nodes[0])->getNumber() << "," << (nodes[1])->getNumber() << ",";
  f << setprecision(8) << length;
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
  logdensity = 0;
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
//      cerr << "calculate called on node " << number << " for site " << site;
//      cerr << " prob for " << pattern << " set to " << patternToProbMap[ pattern ].first << " " << patternToProbMap[ pattern ].second.transpose() << endl;
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
//    cerr << "calculate called on node " << number << " for site " << site;
//    cerr << " prob for " << pattern << " set to " << patternToProbMap[ pattern ].first << " " << patternToProbMap[ pattern ].second.transpose() << endl;
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
    nb->calculate(k,alignment,eb,recurse); //clau: i think we are doing this twice (before calling mleDistance, and here)
    pair<double,Vector4d> pa = na->patternToProbMap[na->getPattern()];
    pair<double,Vector4d> pb = nb->patternToProbMap[nb->getPattern()];
    Vector4d va = pa.second;
    Vector4d vq = qmatrix.getStationaryP(); //fixit: do we need to do this inside for each site?
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

// May 18, 2016
// Assumes Exponential( lambda ) prior on t
// void Tree::partialMaxPosteriorPathCalculations(double t,Alignment& alignment,Node* na,Edge* ea,Node* nb,Edge* eb,QMatrix& qmatrix,
// 					       double& logl,double& dlogl,double& ddlogl,double lambda,bool recurse)
// {
//   partialPathCalculations(t,alignment,na,ea,nb,eb,qmatrix,logl,dlogl,ddlogl,recurse);
//   logl += ( log(lambda) - lambda*t);
//   dlogl -= lambda;
// }

void Tree::partialPathCalculations3D(Vector3d t,Alignment& alignment,Node* nx,Edge* ex,Node* ny,Edge* ey,Node* nz,Edge* ez,QMatrix& qmatrix,double& logl,Vector3d& gradient,Matrix3d& hessian,bool recurse)
{
  Matrix4d P1 = qmatrix.getTransitionMatrix(t[0]);
  Matrix4d QP1 = qmatrix.getQP(t[0]);
  Matrix4d QQP1 = qmatrix.getQQP(t[0]);
  Matrix4d P2 = qmatrix.getTransitionMatrix(t[1]);
  Matrix4d QP2 = qmatrix.getQP(t[1]);
  Matrix4d QQP2 = qmatrix.getQQP(t[1]);
  Matrix4d P3 = qmatrix.getTransitionMatrix(t[2]);
  Matrix4d QP3 = qmatrix.getQP(t[2]);
  Matrix4d QQP3 = qmatrix.getQQP(t[2]);
  Vector4d vq = qmatrix.getStationaryP();
  Vector4d ones(1,1,1,1);

  int numSites = alignment.getNumSites();

  logl = 0;
  double dll1 = 0;
  double dll2 = 0;
  double dll3 = 0;
  double d2ll_11 = 0;
  double d2ll_12 = 0;
  double d2ll_13 = 0;
  double d2ll_22 = 0;
  double d2ll_23 = 0;
  double d2ll_33 = 0;

  for ( int k=0; k<numSites; ++k )
  {
    nx->calculate(k,alignment,ex,recurse); // set pattern and put probability in map if not already there
    ny->calculate(k,alignment,ey,recurse);
    nz->calculate(k,alignment,ez,recurse);
    pair<double,Vector4d> px = nx->patternToProbMap[nx->getPattern()];
    pair<double,Vector4d> py = ny->patternToProbMap[ny->getPattern()];
    pair<double,Vector4d> pz = nz->patternToProbMap[nz->getPattern()];

    Vector4d S1 = P1 * px.second.asDiagonal() * ones;
    Vector4d S2 = P2 * py.second.asDiagonal() * ones;
    Vector4d S3 = P3 * pz.second.asDiagonal() * ones;
    Vector4d S1pr = QP1 * px.second.asDiagonal() * ones;
    Vector4d S2pr = QP2 * py.second.asDiagonal() * ones;
    Vector4d S3pr = QP3 * pz.second.asDiagonal() * ones;
    Vector4d S1doublepr = QQP1 * px.second.asDiagonal() * ones;
    Vector4d S2doublepr = QQP2 * py.second.asDiagonal() * ones;
    Vector4d S3doublepr = QQP3 * pz.second.asDiagonal() * ones;

    //fixt: make a function of this
    // question: how to make a function to return many things?
    double fk = vectorProduct4D(vq,S1,S2,S3);
    double fkpr1 = vectorProduct4D(vq,S1pr,S2,S3);
    double fkpr2 = vectorProduct4D(vq,S1,S2pr,S3);
    double fkpr3 = vectorProduct4D(vq,S1,S2,S3pr);
    double fkdoublepr11 = vectorProduct4D(vq,S1doublepr,S2,S3);
    double fkdoublepr22 = vectorProduct4D(vq,S1,S2doublepr,S3);
    double fkdoublepr33 = vectorProduct4D(vq,S1,S2,S3doublepr);
    double fkdoublepr12 = vectorProduct4D(vq,S1pr,S2pr,S3);
    double fkdoublepr13 = vectorProduct4D(vq,S1pr,S2,S3pr);
    double fkdoublepr23 = vectorProduct4D(vq,S1,S2pr,S3pr);

    logl += px.first + py.first + pz.first + log( fk );
    dll1 += fkpr1/fk;
    dll2 += fkpr2/fk;
    dll3 += fkpr3/fk;
    d2ll_11 += (fk*fkdoublepr11 - fkpr1*fkpr1)/(fk*fk);
    d2ll_12 += (fk*fkdoublepr12 - fkpr1*fkpr2)/(fk*fk);
    d2ll_13 += (fk*fkdoublepr13 - fkpr1*fkpr3)/(fk*fk);
    d2ll_22 += (fk*fkdoublepr22 - fkpr2*fkpr2)/(fk*fk);
    d2ll_23 += (fk*fkdoublepr23 - fkpr2*fkpr3)/(fk*fk);
    d2ll_33 += (fk*fkdoublepr33 - fkpr3*fkpr3)/(fk*fk);
  }
  gradient << dll1, dll2, dll3;
  hessian << d2ll_11, d2ll_12, d2ll_13, d2ll_12, d2ll_22, d2ll_23, d2ll_13, d2ll_23, d2ll_33; //row-wise
  //MatrixXd L( hessian.llt().matrixL() );
  //cerr << "LLT of Hessian 3D: " << endl << L << endl;
}

// similar to partialPathCalculations3D, t=(t1,t2), sum: t3=sum-t1
// warning: need to be careful in the order of t
// gradient and hessian are 2d now
void Tree::partialPathCalculations2D(Vector2d t, double sum,Alignment& alignment,Node* nx,Edge* ex,Node* ny,Edge* ey,Node* nz,Edge* ez,QMatrix& qmatrix,double& logl,Vector2d& gradient,Matrix2d& hessian,bool recurse)
{
  if (sum < t[0])
    {
      cerr << "Sum smaller than summand in partialPathCalculations2D" << endl;
      cerr << "t: " << t.transpose() << endl;
      cerr << "sum: " << sum << endl;
      //exit(1);
      throw 20;
    }
  Matrix4d P1 = qmatrix.getTransitionMatrix(t[0]);
  Matrix4d QP1 = qmatrix.getQP(t[0]);
  Matrix4d QQP1 = qmatrix.getQQP(t[0]);
  Matrix4d P2 = qmatrix.getTransitionMatrix(t[1]);
  Matrix4d QP2 = qmatrix.getQP(t[1]);
  Matrix4d QQP2 = qmatrix.getQQP(t[1]);
  Matrix4d P3 = qmatrix.getTransitionMatrix(sum-t[0]);
  Matrix4d QP3 = qmatrix.getQP(sum-t[0]);
  Matrix4d QQP3 = qmatrix.getQQP(sum-t[0]);
  Vector4d vq = qmatrix.getStationaryP();
  Vector4d ones(1,1,1,1);

  int numSites = alignment.getNumSites();

  logl = 0;
  double dll1 = 0;
  double dll2 = 0;
  double d2ll_11 = 0;
  double d2ll_12 = 0;
  double d2ll_22 = 0;

  for ( int k=0; k<numSites; ++k )
  {
    nx->calculate(k,alignment,ex,recurse); // set pattern and put probability in map if not already there
    ny->calculate(k,alignment,ey,recurse);
    nz->calculate(k,alignment,ez,recurse);
    pair<double,Vector4d> px = nx->patternToProbMap[nx->getPattern()];
    pair<double,Vector4d> py = ny->patternToProbMap[ny->getPattern()];
    pair<double,Vector4d> pz = nz->patternToProbMap[nz->getPattern()];

    // cout << "Inside partialPathCalc2D" << endl;
    // cout << "site " << k << ", recurse " << recurse << endl;
    // cout << "Nodes " << nx->getNumber() << ", " << ny->getNumber() << ", " << nz->getNumber() << endl;
    // cout << "Probs" << px.second << ", " << py.second << ", " << pz.second << endl;

    Vector4d S1 = P1 * px.second.asDiagonal() * ones;
    Vector4d S2 = P2 * py.second.asDiagonal() * ones;
    Vector4d S3 = P3 * pz.second.asDiagonal() * ones;
    Vector4d S1pr = QP1 * px.second.asDiagonal() * ones;
    Vector4d S2pr = QP2 * py.second.asDiagonal() * ones;
    Vector4d S3pr = (-1) * QP3 * pz.second.asDiagonal() * ones;
    Vector4d S1doublepr = QQP1 * px.second.asDiagonal() * ones;
    Vector4d S2doublepr = QQP2 * py.second.asDiagonal() * ones;
    Vector4d S3doublepr = QQP3 * pz.second.asDiagonal() * ones;

    //    cout << "vq as diagonal " << vq.asDiagonal() << ", S1 as diagonal " << S1.asDiagonal() << endl;
    // question: how to make a function here?
    double fk = vectorProduct4D(vq,S1,S2,S3);
    double fkpr1 = vectorProduct4D(vq,S1pr,S2,S3);
    double fkpr2 = vectorProduct4D(vq,S1,S2pr,S3);
    double fkpr3 = vectorProduct4D(vq,S1,S2,S3pr);
    double fkdoublepr11 = vectorProduct4D(vq,S1doublepr,S2,S3);
    double fkdoublepr22 = vectorProduct4D(vq,S1,S2doublepr,S3);
    double fkdoublepr33 = vectorProduct4D(vq,S1,S2,S3doublepr);
    double fkdoublepr12 = vectorProduct4D(vq,S1pr,S2pr,S3);
    double fkdoublepr13 = vectorProduct4D(vq,S1pr,S2,S3pr);
    double fkdoublepr23 = vectorProduct4D(vq,S1,S2pr,S3pr);


    logl += px.first + py.first + pz.first + log( fk );
    dll1 += (fkpr1+fkpr3)/fk;
    dll2 += fkpr2/fk;
    d2ll_11 += (fk*(fkdoublepr11+2*fkdoublepr13+fkdoublepr33) - (fkpr1+fkpr3)*(fkpr1+fkpr3))/(fk*fk);
    d2ll_12 += (fk*(fkdoublepr12+fkdoublepr23) - fkpr2*(fkpr1+fkpr3))/(fk*fk);
    d2ll_22 += (fk*fkdoublepr22 - fkpr2*fkpr2)/(fk*fk);
  }
  gradient << dll1, dll2;
  hessian << d2ll_11, d2ll_12, d2ll_12, d2ll_22; //row-wise
}

// similar to partialPathCalculations2D, t=t1, sum1,sum2
  void Tree::partialPathCalculations1D(double t1, double sum1, double sum2, Alignment& alignment,Node* nx,Edge* ex,Node* ny,Edge* ey,Node* nz,Edge* ez,QMatrix& qmatrix,double& logl,double& dlogl,double& ddlogl,bool recurse)
{
  if (sum1 < t1 || sum2 < t1)
    {
      cerr << "Sum " << sum1 << " or " << sum2 << "smaller than summand " << t1 << " in partialPathCalculations1D" << endl;
      //exit(1);
      throw 20;
    }
  Matrix4d P1 = qmatrix.getTransitionMatrix(t1);
  Matrix4d QP1 = qmatrix.getQP(t1);
  Matrix4d QQP1 = qmatrix.getQQP(t1);
  Matrix4d P2 = qmatrix.getTransitionMatrix(sum1-t1);
  Matrix4d QP2 = qmatrix.getQP(sum1-t1);
  Matrix4d QQP2 = qmatrix.getQQP(sum1-t1);
  Matrix4d P3 = qmatrix.getTransitionMatrix(sum2-t1);
  Matrix4d QP3 = qmatrix.getQP(sum2-t1);
  Matrix4d QQP3 = qmatrix.getQQP(sum2-t1);
  Vector4d vq = qmatrix.getStationaryP();
  Vector4d ones(1,1,1,1);

  int numSites = alignment.getNumSites();

  logl = 0;
  dlogl = 0;
  ddlogl = 0;

  for ( int k=0; k<numSites; ++k )
  {
    nx->calculate(k,alignment,ex,recurse); // set pattern and put probability in map if not already there
    ny->calculate(k,alignment,ey,recurse);
    nz->calculate(k,alignment,ez,recurse);
    pair<double,Vector4d> px = nx->patternToProbMap[nx->getPattern()];
    pair<double,Vector4d> py = ny->patternToProbMap[ny->getPattern()];
    pair<double,Vector4d> pz = nz->patternToProbMap[nz->getPattern()];

    Vector4d S1 = P1 * px.second.asDiagonal() * ones;
    Vector4d S2 = P2 * py.second.asDiagonal() * ones;
    Vector4d S3 = P3 * pz.second.asDiagonal() * ones;
    Vector4d S1pr = QP1 * px.second.asDiagonal() * ones;
    Vector4d S2pr = (-1) * QP2 * py.second.asDiagonal() * ones;
    Vector4d S3pr = (-1) * QP3 * pz.second.asDiagonal() * ones;
    Vector4d S1doublepr = QQP1 * px.second.asDiagonal() * ones;
    Vector4d S2doublepr = QQP2 * py.second.asDiagonal() * ones;
    Vector4d S3doublepr = QQP3 * pz.second.asDiagonal() * ones;

    //    cout << "vq as diagonal " << vq.asDiagonal() << ", S1 as diagonal " << S1.asDiagonal() << endl;
    // question: how to make a function here?
    double fk = vectorProduct4D(vq,S1,S2,S3);
    double fkpr1 = vectorProduct4D(vq,S1pr,S2,S3);
    double fkpr2 = vectorProduct4D(vq,S1,S2pr,S3);
    double fkpr3 = vectorProduct4D(vq,S1,S2,S3pr);
    double fkdoublepr11 = vectorProduct4D(vq,S1doublepr,S2,S3);
    double fkdoublepr22 = vectorProduct4D(vq,S1,S2doublepr,S3);
    double fkdoublepr33 = vectorProduct4D(vq,S1,S2,S3doublepr);
    double fkdoublepr12 = vectorProduct4D(vq,S1pr,S2pr,S3);
    double fkdoublepr13 = vectorProduct4D(vq,S1pr,S2,S3pr);
    double fkdoublepr23 = vectorProduct4D(vq,S1,S2pr,S3pr);

    logl += px.first + py.first + pz.first + log( fk );
    dlogl += (fkpr1+fkpr2+fkpr3)/fk;
    ddlogl += (fk*(fkdoublepr11+2*fkdoublepr12+2*fkdoublepr13+2*fkdoublepr23+fkdoublepr22+fkdoublepr33) - (fkpr1+fkpr2+fkpr3)*(fkpr1+fkpr2+fkpr3))/(fk*fk);
  }
}

//clau: I dont think we need these 3 functions
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
  //exit(1);
  throw 20;
}

void mleErrorJoint(Node* nx,Node* ny,Node* nz)
{
  cerr << "Error: too many iterations in mleDistanceJoint." << endl;
  cerr << "Nodes " << nx->getNumber() << ", " << ny->getNumber() << ", and " << nz->getNumber() << endl;
  //exit(1); 
  throw 20;
}

// Find mle distance from node a to b through a path that uses edges ea and eb.
// Conditon on data in subtrees through other edges.
// Assumes that edge lengths in these subtrees exist
// and that the patternToProbMaps are accurate if edges ea and eb head toward the root.
double Tree::mleDistance(Alignment& alignment,Node* na,Edge* ea,Node* nb,Edge* eb,QMatrix& qmatrix, bool verbose)
{
  if(verbose)
    cout << "mleDistance for nodes " << na->getNumber() << ", " << nb->getNumber() << endl;
  bool recurse=true; //warning: need to keep true
  int iter=0;
  double curr = 0.05;
  // if(!(na->getLeaf()) || !(nb->getLeaf())) //fixit: does not compile
  //   gridPlotFile(alignment,na,ea,nb,eb,qmatrix);
  // get a decent starting point
  double curr_logl,curr_dlogl,curr_ddlogl;
  partialPathCalculations(curr,alignment,na,ea,nb,eb,qmatrix,curr_logl,curr_dlogl,curr_ddlogl,recurse);
  if(verbose)
    cout << "Starting at curr="<<curr << ", logl, dlogl and ddlogl: " << curr_logl << ", " << curr_dlogl << ", " << curr_ddlogl << endl;
  double prop = curr;
  double prop_logl = curr_logl;
  double prop_dlogl = curr_dlogl;
  double prop_ddlogl = curr_ddlogl;
  if ( curr_dlogl > 0 ) // clau: this is for a good starting point
  {
    do
    {
      curr = prop;
      curr_logl = prop_logl;
      curr_dlogl = prop_dlogl;
      curr_ddlogl = prop_ddlogl;
      // cout << "curr in mleDistance for curr_dlogl>0: " << curr << endl;
      // cout << "curr_logl in mleDistance for curr_dlogl>0: " << curr_logl << endl;
      // cout << "curr_dlogl in mleDistance for curr_dlogl>0: " << curr_dlogl << endl;
      // cout << "curr_ddlogl in mleDistance for curr_dlogl>0: " << curr_ddlogl << endl;
      prop = 2*curr;
      if(verbose)
	cout << "double starting point in mleDistance because curr_dlogl>0" << endl;
      partialPathCalculations(prop,alignment,na,ea,nb,eb,qmatrix,prop_logl,prop_dlogl,prop_ddlogl,recurse);
      if(verbose)
	{
	  cout << "prop in mleDistance for curr_dlogl>0: " << prop << endl;
	  cout << "prop_logl in mleDistance for curr_dlogl>0: " << prop_logl << endl;
	  cout << "prop_dlogl in mleDistance for curr_dlogl>0: " << prop_dlogl << endl;
	  cout << "prop_ddlogl in mleDistance for curr_dlogl>0: " << prop_ddlogl << endl;
	}
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
      // cout << "curr in mleDistance for curr_dlogl<0: " << curr << endl;
      // cout << "curr_logl in mleDistance for curr_dlogl<0: " << curr_logl << endl;
      // cout << "curr_dlogl in mleDistance for curr_dlogl<0: " << curr_dlogl << endl;
      // cout << "curr_ddlogl in mleDistance for curr_dlogl<0: " << curr_ddlogl << endl;
      prop = 0.5*curr;
      if(verbose)
	cout << "half starting point in mleDistance because curr_dlogl<0" << endl;
      partialPathCalculations(prop,alignment,na,ea,nb,eb,qmatrix,prop_logl,prop_dlogl,prop_ddlogl,recurse);
      if(verbose)
	{
	  cout << "prop in mleDistance for curr_dlogl<0: " << prop << endl;
	  cout << "prop_logl in mleDistance for curr_dlogl<0: " << prop_logl << endl;
	  cout << "prop_dlogl in mleDistance for curr_dlogl<0: " << prop_dlogl << endl;
	  cout << "prop_ddlogl in mleDistance for curr_dlogl<0: " << prop_ddlogl << endl;
	}
      if ( ++iter > 100 )
        mleError(na,nb,curr,prop,curr_dlogl,prop_dlogl);
    } while ( prop_dlogl < 0 );
  }
  // switch to protected Newton-Raphson
  if(verbose)
    {
      cout << "Two points to interpolate: " << curr << ", " << prop << endl;
      cout << "Derivatives: " << curr_dlogl << ", " << prop_dlogl << endl;
    }
  prop = curr - curr_dlogl * (prop - curr) / (prop_dlogl - curr_dlogl);
  partialPathCalculations(prop,alignment,na,ea,nb,eb,qmatrix,prop_logl,prop_dlogl,prop_ddlogl,recurse);
  if(verbose)
    {
      cout << "Starting prop: " << prop << endl;
      cout << "prop_logl, prop_dlogl, prop_ddlogl" << prop_logl << ", " << prop_dlogl << ", " << prop_ddlogl << endl;
    }
  do
  {
    if(verbose)
      cout << "entered while" << endl;
    curr = prop;
    partialPathCalculations(curr,alignment,na,ea,nb,eb,qmatrix,curr_logl,curr_dlogl,curr_ddlogl,recurse);
    // cout << "curr in mleDistance: " << curr << endl;
    if ( ++iter > 100 )
      mleError(na,nb,curr,prop,curr_dlogl,prop_dlogl);
    double delta = -curr_dlogl / curr_ddlogl;
    // cout << "Delta in mleDistance: " << delta << endl;
    prop = curr + delta;
    // cout << "prop in mleDistance: " << prop << endl;
    while ( prop < 0 )
    {
      // cout << "found negative" << endl;
      delta = 0.5*delta;
      // cout << "Delta: " << delta << endl;
      prop = curr + delta;
    }
    partialPathCalculations(prop,alignment,na,ea,nb,eb,qmatrix,prop_logl,prop_dlogl,prop_ddlogl,recurse);
    if ( ++iter > 100 )
      mleError(na,nb,curr,prop,curr_dlogl,prop_dlogl);
    while ( ( fabs(prop_dlogl) > fabs(curr_dlogl) ) && fabs(curr - prop) > (TOL*TOL) )
    {
      // cout << "found big jump" << endl;
      delta = 0.5*delta;
      // cout << "Delta: " << delta << endl;
      prop = curr + delta;
      partialPathCalculations(prop,alignment,na,ea,nb,eb,qmatrix,prop_logl,prop_dlogl,prop_ddlogl,recurse);
      if ( ++iter > 100 )
        mleError(na,nb,curr,prop,curr_dlogl,prop_dlogl);
    }
  } while ( fabs(curr - prop) > (TOL*TOL) && fabs(prop_dlogl) > (TOL*TOL)); //clau: added stopping rule dlogl~0g
  if(verbose)
    {
      cout << "Finally converged to prop: " << prop << " with logl, dlogl, ddlogl: " << prop_logl << ", " << prop_dlogl << ", " << prop_ddlogl << endl;
      cout << "-----" << endl;
    }
  return prop;
}

// double Tree::maxPosteriorDistance(Alignment& alignment,Node* na,Edge* ea,Node* nb,Edge* eb,QMatrix& qmatrix,double lambda)
// {
//   cerr << "maxPosteriorDistance for nodes " << na->getNumber() << ", " << nb->getNumber() << endl;
//   bool recurse=false;
//   int iter=0;
//   double curr = 0.05;
//   // get a decent starting point
//   double curr_logl,curr_dlogl,curr_ddlogl;
//   partialMaxPosteriorPathCalculations(curr,alignment,na,ea,nb,eb,qmatrix,curr_logl,curr_dlogl,curr_ddlogl,lambda,false);
//   cerr << "Starting at curr="<<curr << ", logl, dlogl and ddlogl: " << curr_logl << ", " << curr_dlogl << ", " << curr_ddlogl << endl;
//   double prop = curr;
//   double prop_logl = curr_logl;
//   double prop_dlogl = curr_dlogl;
//   double prop_ddlogl = curr_ddlogl;
//   if ( curr_dlogl > 0 ) // clau: this is for a good starting point
//   {
//     do
//     {
//       curr = prop;
//       curr_logl = prop_logl;
//       curr_dlogl = prop_dlogl;
//       curr_ddlogl = prop_ddlogl;
//       prop = 2*curr;
//       cerr << "double starting point in mleDistance because curr_dlogl>0" << endl;
//       partialMaxPosteriorPathCalculations(prop,alignment,na,ea,nb,eb,qmatrix,prop_logl,prop_dlogl,prop_ddlogl,lambda,recurse);
//       cerr << "prop in mleDistance for curr_dlogl>0: " << prop << endl;
//       cerr << "prop_logl in mleDistance for curr_dlogl>0: " << prop_logl << endl;
//       cerr << "prop_dlogl in mleDistance for curr_dlogl>0: " << prop_dlogl << endl;
//       cerr << "prop_ddlogl in mleDistance for curr_dlogl>0: " << prop_ddlogl << endl;
//       if ( ++iter > 100 )
//         mleError(na,nb,curr,prop,curr_dlogl,prop_dlogl);
//     } while ( prop_dlogl > 0);
//   }
//   else
//   {
//     do
//     {
//       curr = prop;
//       curr_logl = prop_logl;
//       curr_dlogl = prop_dlogl;
//       curr_ddlogl = prop_ddlogl;
//       prop = 0.5*curr;
//       cerr << "half starting point in mleDistance because curr_dlogl<0" << endl;
//       partialMaxPosteriorPathCalculations(prop,alignment,na,ea,nb,eb,qmatrix,prop_logl,prop_dlogl,prop_ddlogl,lambda,recurse);
//       cerr << "prop in mleDistance for curr_dlogl<0: " << prop << endl;
//       cerr << "prop_logl in mleDistance for curr_dlogl<0: " << prop_logl << endl;
//       cerr << "prop_dlogl in mleDistance for curr_dlogl<0: " << prop_dlogl << endl;
//       cerr << "prop_ddlogl in mleDistance for curr_dlogl<0: " << prop_ddlogl << endl;
//       if ( ++iter > 100 )
//         mleError(na,nb,curr,prop,curr_dlogl,prop_dlogl);
//     } while ( prop_dlogl < 0 );
//   }
//   // switch to protected Newton-Raphson
//   cerr << "Two points to interpolate: " << curr << ", " << prop << endl;
//   cerr << "Derivatives: " << curr_dlogl << ", " << prop_dlogl << endl;
//   prop = curr - curr_dlogl * (prop - curr) / (prop_dlogl - curr_dlogl);
//   partialMaxPosteriorPathCalculations(prop,alignment,na,ea,nb,eb,qmatrix,prop_logl,prop_dlogl,prop_ddlogl,lambda,recurse);
//   cerr << "Starting prop: " << prop << endl;
//   cerr << "prop_logl, prop_dlogl, prop_ddlogl" << prop_logl << ", " << prop_logl << ", " << prop_dlogl << ", " << prop_ddlogl << endl;
//   do
//   {
//     cerr << "entered while" << endl;
//     curr = prop;
//     partialMaxPosteriorPathCalculations(curr,alignment,na,ea,nb,eb,qmatrix,curr_logl,curr_dlogl,curr_ddlogl,lambda,recurse);

//     // cout << "curr in mleDistance: " << curr << endl;
//     if ( ++iter > 100 )
//       mleError(na,nb,curr,prop,curr_dlogl,prop_dlogl);
//     double delta = -curr_dlogl / curr_ddlogl;
//     // cout << "Delta in mleDistance: " << delta << endl;
//     prop = curr + delta;
//     // cout << "prop in mleDistance: " << prop << endl;
//     while ( prop < 0 )
//     {
//       // cout << "found negative" << endl;
//       delta = 0.5*delta;
//       // cout << "Delta: " << delta << endl;
//       prop = curr + delta;
//     }
//     partialMaxPosteriorPathCalculations(prop,alignment,na,ea,nb,eb,qmatrix,prop_logl,prop_dlogl,prop_ddlogl,lambda,recurse);
//     if ( ++iter > 100 )
//       mleError(na,nb,curr,prop,curr_dlogl,prop_dlogl);
//     while ( ( fabs(prop_dlogl) > fabs(curr_dlogl) ) && fabs(curr - prop) >1.0e-8 )
//     {
//       // cout << "found big jump" << endl;
//       delta = 0.5*delta;
//       // cout << "Delta: " << delta << endl;
//       prop = curr + delta;
//       partialMaxPosteriorPathCalculations(prop,alignment,na,ea,nb,eb,qmatrix,prop_logl,prop_dlogl,prop_ddlogl,lambda,recurse);
//       if ( ++iter > 100 )
//         mleError(na,nb,curr,prop,curr_dlogl,prop_dlogl);
//     }
//   } while ( fabs(curr - prop) > 1.0e-8 && fabs(prop_dlogl) > 1.0e-8); //clau: added stopping rule dlogl~0g
//   cout << "Finally converged to prop: " << prop << " with logl, dlogl, ddlogl: " << prop_logl << ", " << prop_dlogl << ", " << prop_ddlogl << endl;
//   cout << "-----" << endl;
//   return prop;
// }

// Find mle distance between nodes nx,ny,nz with edges ex,ey,ez
// conditional on two sums: t1,s1-t1,s2-t1
// warning: input t1,t2,t3,s1,s2
// Conditon on data in subtrees through other edges.
// Assumes that edge lengths in these subtrees exist
// and that the patternToProbMaps are accurate if edges ea and eb head toward the root.
void Tree::mleDistance1D(Alignment& alignment,Node* nx,Edge* ex,Node* ny,Edge* ey,Node* nz,Edge* ez,QMatrix& qmatrix, double& t1, double& t2, double& t3, double& sum1, double& sum2, mt19937_64& rng, bool verbose, bool mvnormal)
{
  if(sum1 < t1 || sum2 < t1)
    {
      cerr << "Sum smaller than summand in mle distance1D" << endl;
      //exit(1);
      throw 20;
    }
  if(verbose)
    cout << "Entering mleDistance1D" << endl;
  bool recurse=true; //warning: need to keep true
  int iter=0;
  double curr = t1;
  double curr_logl,curr_dlogl,curr_ddlogl;
  partialPathCalculations1D(curr,sum1,sum2,alignment,nx,ex,ny,ey,nz,ez,qmatrix,curr_logl,curr_dlogl,curr_ddlogl,recurse); //true un original mleDist
  double prop = curr;
  double prop_logl = curr_logl;
  double prop_dlogl = curr_dlogl;
  double prop_ddlogl = curr_ddlogl;
  if ( curr_dlogl > 0 ) //to find a good starting point
  {
    do
    {
      if(verbose)
	cout << "mleDistance1D Newton-Raphson curr: " << curr << endl;
      curr = prop;
      curr_logl = prop_logl;
      curr_dlogl = prop_dlogl;
      curr_ddlogl = prop_ddlogl;
      prop = 2*curr;
      while ( prop > min(sum1,sum2) )
	prop = (prop+curr)/2;
      partialPathCalculations1D(prop,sum1,sum2,alignment,nx,ex,ny,ey,nz,ez,qmatrix,prop_logl,prop_dlogl,prop_ddlogl,recurse);
      if ( ++iter > 100 )
        mleErrorJoint(nx,ny,nz);
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
      partialPathCalculations1D(prop,sum1,sum2,alignment,nx,ex,ny,ey,nz,ez,qmatrix,prop_logl,prop_dlogl,prop_ddlogl,recurse);
      if ( ++iter > 100 )
        mleErrorJoint(nx,ny,nz);
    } while ( prop_dlogl < 0 );
  }
  // switch to protected Newton-Raphson
  if(verbose)
    {
      cout << "Two points to interpolate: " << curr << ", " << prop << endl;
      cout << "Derivatives: " << curr_dlogl << ", " << prop_dlogl << endl;
    }
  prop = curr - curr_dlogl * (prop - curr) / (prop_dlogl - curr_dlogl);
  partialPathCalculations1D(prop,sum1,sum2,alignment,nx,ex,ny,ey,nz,ez,qmatrix,prop_logl,prop_dlogl,prop_ddlogl,recurse);
  if(verbose)
    {
      cout << "Starting prop: " << prop << endl;
      cout << "prop_logl, prop_dlogl, prop_ddlogl" << prop_logl << ", " << prop_dlogl << ", " << prop_ddlogl << endl;
    }
  do
  {
    curr = prop;
    partialPathCalculations1D(curr,sum1,sum2,alignment,nx,ex,ny,ey,nz,ez,qmatrix,curr_logl,curr_dlogl,curr_ddlogl,recurse);
    if(verbose)
      {
	cout << "mleDistance1D Newton-Raphson curr: " << curr << endl;
	cout << "mleDistance1D Newton-Raphson dlogl: " << curr_dlogl << endl;
	cout << "mleDistance1D Newton-Raphson ddlogl: " << curr_ddlogl << endl;
      }
    if ( ++iter > 100 )
      mleErrorJoint(nx,ny,nz);
    double delta = curr_dlogl / curr_ddlogl;
    prop = curr - delta;
    while ( prop < 0 )
    {
      if(verbose)
	cerr << "found negative" << endl;
      delta = 0.5*delta;
      prop = curr - delta;
    }
    while ( prop > min(sum1,sum2))
      {
	if(verbose)
	  cerr << "found jump bigger than sum" << endl;
	delta = 0.5*delta;
	prop = curr - delta;
      }
    //cerr << "Delta " << delta << endl;
    partialPathCalculations1D(prop,sum1,sum2,alignment,nx,ex,ny,ey,nz,ez,qmatrix,prop_logl,prop_dlogl,prop_ddlogl,recurse);
    if ( ++iter > 100 )
      mleErrorJoint(nx,ny,nz);
    while ( ( fabs(prop_dlogl) > fabs(curr_dlogl) ) && fabs(curr - prop) > (TOL*TOL) )
    {
      if(verbose)
	cerr << "found bigger step" << endl;
      delta = 0.5*delta;
      prop = curr - delta;
      partialPathCalculations1D(prop,sum1,sum2,alignment,nx,ex,ny,ey,nz,ez,qmatrix,prop_logl,prop_dlogl,prop_ddlogl,recurse);
      if ( ++iter > 100 )
        mleErrorJoint(nx,ny,nz);
    }
  } while ( fabs(curr - prop) > (TOL*TOL) && fabs(prop_dlogl) > (TOL*TOL));
  if(verbose)
    {
      cout << "Finally converged to prop: " << prop << " with logl, dlogl, ddlogl: " << prop_logl << ", " << prop_dlogl << ", " << prop_ddlogl << endl;
      cout << "-----" << endl;
    }
  double mu = prop;
  double var = -1/prop_ddlogl;
  if(mvnormal)
    {
      normal_distribution<double> rnorm;
      double z = rnorm(rng);
      t1 = mu + sqrt(var)*z;
      t2 = sum1 - t1;
      t3 = sum2 - t1;
      logdensity += -0.5*z*z;
    }
  else
    {
      Vector3d bl = multivariateGamma1D(mu,var,sum1,sum2,rng, logdensity, verbose);
      if(verbose)
	cout << "bl after multivariateGamma1D: " << bl.transpose() << endl;
      t1 = bl[0];
      t2 = bl[1];
      t3 = bl[2];
    }
  if(t1<0 || t2<0 || t3<0)
    {
      cerr << "Sampled negative branches in mleDistance1D" << endl;
      //exit(1);
      throw 20;
    }
}

// void Tree::maxPosteriorDistance1D(Alignment& alignment,Node* nx,Edge* ex,Node* ny,Edge* ey,Node* nz,Edge* ez,QMatrix& qmatrix, double lambda,
// 				  double& t1, double& t2, double& t3, double& sum1, double& sum2, mt19937_64& rng)
// {
//   if (sum1 < t1 || sum2 < t1)
//     {
//       cerr << "Sum smaller than summand in mle distance1D" << endl;
//       exit(1);
//     }
//   cout << "Entering mleDistance1D" << endl;
//   bool recurse=false;
//   int iter=0;
//   double curr = t1;
//   double curr_logl,curr_dlogl,curr_ddlogl;
//   partialPathCalculations1D(curr,sum1,sum2,alignment,nx,ex,ny,ey,nz,ez,qmatrix,curr_logl,curr_dlogl,curr_ddlogl,true); //true un original mleDist
//   double prop = curr;
//   double prop_logl = curr_logl;
//   double prop_dlogl = curr_dlogl;
//   double prop_ddlogl = curr_ddlogl;
//   if ( curr_dlogl > 0 ) //to find a good starting point
//   {
//     do
//     {
//       cout << "mleDistance1D Newton-Raphson curr: " << curr << endl;
//       curr = prop;
//       curr_logl = prop_logl;
//       curr_dlogl = prop_dlogl;
//       curr_ddlogl = prop_ddlogl;
//       prop = 2*curr;
//       partialPathCalculations1D(prop,sum1,sum2,alignment,nx,ex,ny,ey,nz,ez,qmatrix,prop_logl,prop_dlogl,prop_ddlogl,recurse);
//       if ( ++iter > 100 )
//         mleErrorJoint(nx,ny,nz);
//     } while ( prop_dlogl > 0);
//   }
//   else
//   {
//     do
//     {
//       curr = prop;
//       curr_logl = prop_logl;
//       curr_dlogl = prop_dlogl;
//       curr_ddlogl = prop_ddlogl;
//       prop = 0.5*curr;
//       partialPathCalculations1D(prop,sum1,sum2,alignment,nx,ex,ny,ey,nz,ez,qmatrix,prop_logl,prop_dlogl,prop_ddlogl,recurse);
//       if ( ++iter > 100 )
//         mleErrorJoint(nx,ny,nz);
//     } while ( prop_dlogl < 0 );
//   }
//   // switch to protected Newton-Raphson
//   cout << "Two points to interpolate: " << curr << ", " << prop << endl;
//   cout << "Derivatives: " << curr_dlogl << ", " << prop_dlogl << endl;
//   prop = curr - curr_dlogl * (prop - curr) / (prop_dlogl - curr_dlogl);
//   partialPathCalculations1D(prop,sum1,sum2,alignment,nx,ex,ny,ey,nz,ez,qmatrix,prop_logl,prop_dlogl,prop_ddlogl,recurse);
//   cout << "Starting prop: " << prop << endl;
//   cout << "prop_logl, prop_dlogl, prop_ddlogl" << prop_logl << ", " << prop_logl << ", " << prop_dlogl << ", " << prop_ddlogl << endl;
//   do
//   {
//     curr = prop;
//     partialPathCalculations1D(curr,sum1,sum2,alignment,nx,ex,ny,ey,nz,ez,qmatrix,curr_logl,curr_dlogl,curr_ddlogl,recurse);
//     cout << "mleDistance1D Newton-Raphson curr: " << curr << endl;
//     cout << "mleDistance1D Newton-Raphson dlogl: " << curr_dlogl << endl;
//     cout << "mleDistance1D Newton-Raphson ddlogl: " << endl << curr_ddlogl << endl;
//     if ( ++iter > 100 )
//       mleErrorJoint(nx,ny,nz);
//     double delta = curr_dlogl / curr_ddlogl;
//     prop = curr - delta;
//     while ( prop < 0 )
//     {
//       cerr << "found negative" << endl;
//       delta = 0.5*delta;
//       prop = curr - delta;
//     }
//     //cerr << "Delta " << delta << endl;
//     partialPathCalculations1D(prop,sum1,sum2,alignment,nx,ex,ny,ey,nz,ez,qmatrix,prop_logl,prop_dlogl,prop_ddlogl,recurse);
//     if ( ++iter > 100 )
//       mleErrorJoint(nx,ny,nz);
//     while ( ( fabs(prop_dlogl) > fabs(curr_dlogl) ) && fabs(curr - prop) >1.0e-8 )
//     {
//       cerr << "found bigger step" << endl;
//       delta = 0.5*delta;
//       prop = curr - delta;
//       partialPathCalculations1D(prop,sum1,sum2,alignment,nx,ex,ny,ey,nz,ez,qmatrix,prop_logl,prop_dlogl,prop_ddlogl,recurse);
//       if ( ++iter > 100 )
//         mleErrorJoint(nx,ny,nz);
//     }
//   } while ( fabs(curr - prop) > 1.0e-8 && fabs(prop_dlogl) > 1.0e-8);
//   cout << "Finally converged to prop: " << prop << " with logl, dlogl, ddlogl: " << prop_logl << ", " << prop_dlogl << ", " << prop_ddlogl << endl;
//   cout << "-----" << endl;
//   double mu = prop;
//   double var = -1/prop_ddlogl;
//   double s = min(sum1,sum2);
//   double part1 = (mu * mu * (s-mu)) / (s * var);
//   double a =  part1 - mu / s;
//   double b = part1 - (s - mu) / s;
//   t1 = beta(a,b,rng);
//   //t1 = a/b; //temporarily while we create beta generator r.v.
//   cout << "1D mean: " << mu << ", variance: " << var << endl;
//   cout << "Sample 1D bl: " << t1 << endl;
//   t2 = sum1-t1;
//   t3 = sum2-t1;
// }

// Find joint mle distance between nodes nx,ny,nz with edges ex,ey,ez
// Conditon on data in subtrees through other edges.
// Assumes that edge lengths in these subtrees exist
// and that the patternToProbMaps are accurate if edges ex and ey head toward the root.
// it calls different functions depending on the number of sums known
// double& to modify inside? yes
// warning: this would break if a sum is in fact 0, but i don't think this would happen
void Tree::mleDistanceJoint(Alignment& alignment,Node* nx,Edge* ex,Node* ny,Edge* ey,Node* nz,Edge* ez,QMatrix& qmatrix, double& lx, double& ly, double& lz, double sxy, double sxz, double syz,
			    bool foundXY, bool foundXZ, bool foundYZ, mt19937_64& rng, bool verbose, bool mvnormal, ofstream& par3D, ofstream& par2D)
{ //todo: add ofstream to each mledistance3d,2d and write needed things to file
  if ( verbose )
  {
    cout << "Entering mleDistanceJoint" << endl;
    cout << "lx,ly,lz: " << lx << ", " << ly << ", " << lz << endl;
    cout << "sxy,sxz,syz: " << sxy << ", " << sxz << ", " << syz << endl;
  }
//  if(sxy > 0 && sxz > 0 && syz > 0)
  if ( foundXY && foundXZ && foundYZ )
  {
    cerr << "Trying to condition al all three sums, impossible!" << endl;
    //exit(1);
    throw 20;
  }
//  if ( sxy + sxz + syz == 0) //not conditional
  if ( !foundXY && !foundXZ && !foundYZ ) //not conditional
  {
    mleDistance3D(alignment,nx,ex,ny,ey,nz,ez,qmatrix, lx, ly, lz,rng, verbose,mvnormal, par3D);
    return;
  }
//  if( sxy > 0 && sxz + syz == 0) //conditional on sxy
  if ( foundXY && !foundXZ && !foundYZ ) //conditional on sxy
  {
    mleDistance2D(alignment, nx,ex,nz,ez,ny,ey,qmatrix,lx,lz,ly,sxy,rng, verbose,mvnormal, par2D); //warning in order branches
    return;
  }
//  if( sxz > 0 && sxy + syz == 0) //conditional on sxz
  if ( !foundXY && foundXZ && !foundYZ ) //conditional on sxz
  {
    mleDistance2D(alignment, nx,ex,ny,ey,nz,ez,qmatrix,lx,ly,lz,sxz,rng, verbose,mvnormal, par2D); //warning in order branches
    return;
  }
//  if( syz > 0 && sxz + sxy == 0) //conditional on syz
  if ( !foundXY && !foundXZ && foundYZ ) //conditional on syz
  {
    mleDistance2D(alignment, ny,ey,nx,ex,nz,ez,qmatrix,ly,lx,lz,syz,rng, verbose,mvnormal, par2D); //warning in order branches
    return;
  }
//  if( sxy > 0 && sxz > 0) //conditional on sxy,sxz
  if ( foundXY && foundXZ ) //conditional on sxy,sxz
  {
    mleDistance1D(alignment,nx,ex,ny,ey,nz,ez,qmatrix, lx, ly, lz, sxy, sxz, rng, verbose,mvnormal);
    return;
  }
//  if( sxy > 0 && syz > 0) //conditional on sxy,syz
  if ( foundXY && foundYZ ) //conditional on sxy,syz
  {
    mleDistance1D(alignment,ny,ey,nx,ex,nz,ez,qmatrix, ly, lx, lz, sxy, syz, rng, verbose,mvnormal);
    return;
  }
//  if( sxz > 0 && syz > 0) //conditional on sxz,syz
  if ( foundXZ && foundYZ ) //conditional on sxz,syz
  {
    mleDistance1D(alignment,nz,ez,nx,ex,ny,ey,qmatrix, lz, lx, ly, sxz, syz, rng, verbose,mvnormal);
    return;
  }
  cerr << "Error: called mleDistanceJoint and no case found" << endl;
  throw 20;
}

// Find joint mle distance between nodes nx,ny,nz with edges ex,ey,ez
// Conditon on data in subtrees through other edges.
// Assumes that edge lengths in these subtrees exist
// and that the patternToProbMaps are accurate if edges ex and ey head toward the root.
void Tree::mleDistance3D(Alignment& alignment,Node* nx,Edge* ex,Node* ny,Edge* ey,Node* nz,Edge* ez,QMatrix& qmatrix, double& lx, double& ly, double& lz, mt19937_64& rng, bool verbose, bool mvnormal, ofstream& par3D)
{
  if(verbose)
    cout << "Entering mleDistance3D" << endl;
  bool recurse=true; //warning: need to keep true
  int iter=0;
  Vector3d curr(lx,ly,lz);
  double curr_logl;
  Vector3d curr_gradient;
  Matrix3d curr_hessian;
  partialPathCalculations3D(curr,alignment,nx,ex,ny,ey,nz,ez,qmatrix,curr_logl,curr_gradient,curr_hessian,recurse); //true in original mleDist
  Vector3d prop = curr;
  double prop_logl = curr_logl;
  Vector3d prop_gradient = curr_gradient;
  Matrix3d prop_hessian = curr_hessian;
  Vector3d delta = curr - prop;
  bool keepZero1 = false; //whether to keep that entry at 0
  bool keepZero2 = false;
  bool keepZero3 = false;

  // find starting point: we want a good prop
  do
  {
    curr = prop;
    partialPathCalculations3D(curr,alignment,nx,ex,ny,ey,nz,ez,qmatrix,curr_logl,curr_gradient,curr_hessian,recurse);
    if(verbose)
      {
	cout << "mleDistance3D Newton-Raphson curr: " << curr.transpose() << endl;
	cout << "mleDistance3D Newton-Raphson gradient: " << curr_gradient.transpose() << endl;
	cout << "mleDistance3D Newton-Raphson inverse hessian: " << endl << curr_hessian.inverse() << endl;
      }
    if ( ++iter > 100 )
      mleErrorJoint(nx,ny,nz);
    delta = curr_hessian.inverse() * curr_gradient;
    if(keepZero1)
      delta[0] = 0;
    if(keepZero2)
      delta[1] = 0;
    if(keepZero3)
      delta[2] = 0;
    prop = curr - delta;
    if(curr[0] > TOL && curr[1] > TOL && curr[2] > TOL)
      {
	while ( prop[0] < 0 || prop[1] < 0 || prop[2] < 0)
	  {
	    if(verbose)
	      cerr << "found negative with big curr, will shrink delta" << endl;
	    if(prop[0] < 0)
	      delta[0] = 0.5* delta[0];
	    if(prop[1] < 0)
	      delta[1] = 0.5* delta[1];
	    if(prop[2] < 0)
	      delta[2] = 0.5* delta[2];
	    prop = curr - delta;
	  }
      }
    //    cerr << "Delta " << delta.transpose() << endl;
    if(prop[0] < 0 && curr[0] < TOL)
      {
	if(verbose)
	  cerr << "found negative for 1st element with curr small, will set to zero" << endl;
    	prop[0] = 0;
	keepZero1 = true;
      }
    if(prop[1] < 0 && curr[1] < TOL)
      {
	if(verbose)
	  cerr << "found negative for 2nd element with curr small, will set to zero" << endl;
    	prop[1] = 0;
	keepZero2 = true;
      }
    if(prop[2] < 0 && curr[2] < TOL)
      {
	if(verbose)
	  cerr << "found negative for 3rd element with curr small, will set to zero" << endl;
    	prop[2] = 0;
	keepZero3 = true;
      }
    if(verbose)
      cout << "prop befofe partialPathCalculations3D: " << prop.transpose() << endl;

    partialPathCalculations3D(prop,alignment,nx,ex,ny,ey,nz,ez,qmatrix,prop_logl,prop_gradient,prop_hessian,recurse);
    if ( ++iter > 100 )
      mleErrorJoint(nx,ny,nz);
    if(!keepZero1 && !keepZero2 && !keepZero3)
      {
	while ( prop_gradient.squaredNorm() > curr_gradient.squaredNorm() && delta.squaredNorm() > (TOL*TOL) )
	  {
	    if(verbose)
	      cerr << "found bigger step" << endl;
	    if(keepZero1)
		delta[0] = 0;
	    if(keepZero2)
		delta[1] = 0;
	    if(keepZero3)
		delta[2] = 0;
	    delta = 0.5 *delta;
	    prop = curr - delta;
	    partialPathCalculations3D(prop,alignment,nx,ex,ny,ey,nz,ez,qmatrix,prop_logl,prop_gradient,prop_hessian,recurse);
	    if ( ++iter > 100 )
	      mleErrorJoint(nx,ny,nz);
	  }
      }
  } while ( delta.squaredNorm() > (TOL*TOL) && prop_gradient.squaredNorm() > (TOL*TOL));
  Matrix3d cov = (-1) * prop_hessian.inverse();
  par3D << prop[0] << "," << prop[1] << "," << prop[2] << "," << cov(0,0) << "," << cov(0,1) << "," << cov(0,2) << "," << cov(1,0) << "," << cov(1,1) << "," << cov(1,2) << "," << cov(2,0) << "," << cov(2,1) << "," << cov(2,2) << ",";
  if(verbose)
    {
      cout << "Finally converged to" << endl;
      cout << "Gradient " << endl << prop_gradient.transpose() << endl;
      cout << "3D mean: " << prop.transpose() << endl << ", cov matrix: " << endl << cov << endl;
    }
  Vector3d bl;
  if(mvnormal)
    {
      bl = multivariateNormal(prop,cov,rng, logdensity);
      par3D << "0,0,0,0,0,0" <<endl;
    }
  else
    bl = multivariateGamma3D(prop,cov,rng, logdensity, verbose,par3D);
  if(verbose)
    cout << "Sample bl 3D: "<< bl.transpose() << endl;
  //Vector3d bl = prop;
  lx = bl[0];
  ly = bl[1];
  lz = bl[2];
  if( lx<0 || ly<0 || lz<0)
    {
      cerr << "Sampled negative bl in mleDistance3D" << endl;
      //exit(1);
      throw 20;
    }
}

// Find joint mle distance between nodes nx,ny,nz with edges ex,ey,ez
// conditional in one sum: t1,t2,s-t1
// warning: input t1,t2,t3,sum
// warning: nx is associated with t1, ny with t2, and nz with t3, careful in mleDistanceJoint
// Conditon on data in subtrees through other edges.
// Assumes that edge lengths in these subtrees exist
// and that the patternToProbMaps are accurate if edges ex and ey head toward the root.
void Tree::mleDistance2D(Alignment& alignment,Node* nx,Edge* ex,Node* ny,Edge* ey,Node* nz,Edge* ez,QMatrix& qmatrix, double& t1, double& t2, double& t3, double& sum, mt19937_64& rng, bool verbose, bool mvnormal, ofstream& par2D)
{
  if(verbose)
    cout << "Entering mleDistance2D" << endl;
  if (sum < t1)
    {
      cerr << "Sum smaller than summand in mleDistance2D" << endl;
      //exit(1);
      throw 20;
    }
  bool recurse=true; //warning: need to keep true
  int iter=0;
  Vector2d curr(t1, t2);
  double curr_logl;
  Vector2d curr_gradient;
  Matrix2d curr_hessian;
  partialPathCalculations2D(curr,sum,alignment,nx,ex,ny,ey,nz,ez,qmatrix,curr_logl,curr_gradient,curr_hessian,recurse); //true in original mleDist
  if(verbose)
    {
      cout << "Starting point in mleDistance2D" << endl;
      cout << "mleDistance2D Newton-Raphson curr: " << curr.transpose() << endl;
      cout << "mleDistance2D Newton-Raphson gradient: " << curr_gradient.transpose() << endl;
      cout << "mleDistance2D Newton-Raphson inverse hessian: " << endl << curr_hessian.inverse() << endl;
    }
  Vector2d prop = curr;
  double prop_logl = curr_logl;
  Vector2d prop_gradient = curr_gradient;
  Matrix2d prop_hessian = curr_hessian;
  Vector2d delta = curr - prop;
  bool keepZero1 = false; //whether to keep that entry at 0
  bool keepZero2 = false;
  // still need to find a starting point
  do
  {
    if(verbose)
      cout << "entered while" << endl;
    curr = prop;
    partialPathCalculations2D(curr,sum,alignment,nx,ex,ny,ey,nz,ez,qmatrix,curr_logl,curr_gradient,curr_hessian,recurse);
    if(verbose)
      {
	cout << "mleDistance2D Newton-Raphson curr: " << curr.transpose() << endl;
	cout << "mleDistance2D Newton-Raphson gradient: " << curr_gradient.transpose() << endl;
	cout << "mleDistance2D Newton-Raphson inverse hessian: " << endl << curr_hessian.inverse() << endl;
      }
    if ( ++iter > 100 )
      mleErrorJoint(nx,ny,nz);
    delta = curr_hessian.inverse() * curr_gradient;
    if(keepZero1)
      delta[0] = 0;
    if(keepZero2)
      delta[1] = 0;
    if(verbose)
      cout << "First delta: " << delta.transpose() << endl;
    prop = curr - delta;
    if(verbose)
      {
	cout << "New proposed: " << prop.transpose() << endl;
	cout << "with sum: " << sum << endl;
      }
    if(curr[0] > TOL && curr[1] > TOL)
      {
	while ( prop[0] < 0 || prop[1] < 0)
	  {
	    if(verbose)
	      cerr << "found negative with curr big, will shrink delta" << endl;
	    // if(prop[0] < 0)
	    //   delta[0] = 0.5* delta[0];
	    // if(prop[1] < 0)
	    //   delta[1] = 0.5* delta[1];
	    delta = 0.5* delta;
	    if(verbose)
	      cerr << "new delta: " << delta.transpose() << endl;
	    prop = curr - delta;
	  }
      }
    //cerr << "Delta " << delta.transpose() << endl;
    if(prop[0] < 0 && curr[0] < TOL)
      {
	if(verbose)
	  cerr << "found negative for 1st element with curr small, will set to zero" << endl;
    	prop[0] = 0;
	keepZero1 = true;
      }
    if(prop[1] < 0 && curr[1] < TOL)
      {
	if(verbose)
	  cerr << "found negative for 2nd element with curr small, will set to zero" << endl;
    	prop[1] = 0;
	keepZero2 = true;
      }
    if(verbose)
      cout << "prop befofe partialPathCalculations2D: " << prop.transpose() << endl;
    partialPathCalculations2D(prop,sum,alignment,nx,ex,ny,ey,nz,ez,qmatrix,prop_logl,prop_gradient,prop_hessian,recurse);
    if ( ++iter > 100 )
      mleErrorJoint(nx,ny,nz);
    // if(!keepZero1 && !keepZero2)
    //   {
    // 	while ( prop_gradient.squaredNorm() > curr_gradient.squaredNorm() && delta.squaredNorm() > (TOL*TOL) )
    // 	  {
    // 	    if(verbose)
    // 	      cerr << "found bigger step" << endl;
    // 	    if(keepZero1)
    // 		delta[0] = 0;
    // 	    if(keepZero2)
    // 		delta[1] = 0;
    // 	    delta = delta*0.5;
    // 	    if(verbose)
    // 	      cerr << "new delta: " << delta.transpose() << endl;
    // 	    prop = curr - delta;
    // 	    if(verbose)
    // 	      cout << "prop befofe partialPathCalculations2D inside check for bigger step: " << prop.transpose() << endl;
    // 	    partialPathCalculations2D(prop,sum,alignment,nx,ex,ny,ey,nz,ez,qmatrix,prop_logl,prop_gradient,prop_hessian,recurse);
    // 	    if ( ++iter > 100 )
    // 	      mleErrorJoint(nx,ny,nz);
    // 	  }
    //   }
    } while ( delta.squaredNorm() > (TOL*TOL) && prop_gradient.squaredNorm() > (TOL*TOL) );
  if(verbose)
    cout << "Finally converged. delta squared norm: " << delta.squaredNorm() << " and prop gradient squared norm " << prop_gradient.squaredNorm() << "with TOL " << TOL << endl;
  if(prop[0] < 0)
    {
      if(verbose)
	cout << "after while, we need to fix one negative" << endl;
      prop[0] = 0;
      partialPathCalculations2D(prop,sum,alignment,nx,ex,ny,ey,nz,ez,qmatrix,prop_logl,prop_gradient,prop_hessian,recurse);
    }
  if(prop[1] < 0)
    {
      if(verbose)
	cout << "after while, we need to fix one negative" << endl;
      prop[1] = 0 ;
      partialPathCalculations2D(prop,sum,alignment,nx,ex,ny,ey,nz,ez,qmatrix,prop_logl,prop_gradient,prop_hessian,recurse);
    }

  Matrix2d cov = (-1) * prop_hessian.inverse();
  par2D << prop[0] << "," << prop[1] << "," << cov(0,0) << "," << cov(0,1) << "," << cov(1,0) << "," << cov(1,1) << ",";
  if(verbose)
    {
      cout << "Gradient " << endl << prop_gradient.transpose() << endl;
      cout << "2D mean: " << prop.transpose() << endl << ", cov matrix: " << endl << cov << endl;
    }
  //  double logdensity = 0; //fixit: needs to be the logdensity saved until this point
  Vector3d bl;
  if(mvnormal)
    {
      Vector2d bl12 = multivariateNormal(prop,cov,rng, logdensity); //t3=sum-t1
      par2D << "0,0,0,0" << endl;
      bl[0] = bl12[0];
      bl[1] = bl12[1];
      bl[2] = sum - bl[0];
    }
  else
    bl = multivariateGamma2D(prop,cov,sum,rng, logdensity, verbose, par2D); //t3=sum-t1
  if(verbose)
    cout << "Sample 2D bl: " << bl.transpose() << endl;
  t1 = bl[0];
  t2 = bl[1];
  t3 = bl[2];
  if(t1<0 || t2<0 || t3<0)
    {
      cerr << "Sampled negative bl in mleDistance2D" << endl;
      //exit(1);
      throw 20;
    }
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
    //exit(1);
    throw 20;
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

void Tree::generateBranchLengths(Alignment& alignment,QMatrix& qmatrix, mt19937_64& rng, bool verbose, bool mvnormal, ofstream& par3D, ofstream& par2D) //clau: added seed, verbose, mvnormal
{
  //cout << "Starting generateBL with verbose: " << verbose << endl;

  map<pair<int,int>,double> distanceMap;
  list<Node*> nodeList;
  depthFirstNodeList(nodeList);
  setActiveChildrenAndNodeParents();
  logdensity = 0;
  if(verbose)
    {
      cout << "Node List:";
      for ( list<Node*>::iterator p=nodeList.begin(); p!= nodeList.end(); ++p )
	cout << " " << (*p)->getNumber();
      cout << endl;
    }

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
        //exit(1);
	throw 20;
      }
      par = root;
      x = *p++;
      y = *p++;
      z = *p;
    }
    else //clau: p parent not root
    {
      x = *p++; //it means take *p and move right
      y = *p++;
      par = x->getNodeParent();
      z = par->closeRelative();
    }
    if(verbose)
      {
	cout << "------------- Sample branch lengths for nodes:------------------" << endl;
	cout << x->getNumber() << " " << y->getNumber() << " " << z->getNumber() << " " << par->getNumber() << endl;
      }
    // do calculations
    x->calculateEdges(qmatrix); //ignores the parent edge, only for children edges
    y->calculateEdges(qmatrix);
    z->calculateEdges(qmatrix);
    for ( int k=0; k<alignment.getNumSites(); ++k )
    {
      x->calculate(k,alignment,x->getEdgeParent(),false); //clau: recurse=false
      y->calculate(k,alignment,y->getEdgeParent(),false);
      z->calculate(k,alignment,z->getEdgeParent(),false);
    } //clau: after this, we have in each node the patternProb map, but we dont know how many different patterns, do we? no, we search the pattern in mleDistance


    map<pair<int,int>,double>::iterator m;
    double dxy, dxz, dyz;
    bool foundxy = false; //clau: depending on which bl are found it is condition on 1 or 2 sums
    bool foundyz = false;
    bool foundxz = false;
    m = distanceMap.find( getPair(x->getNumber(),y->getNumber()) );
    if ( m == distanceMap.end() ) //clau: did not find x,y; when do you put them inside distanceMap? at the end here
      dxy = mleDistance(alignment,x,x->getEdgeParent(),y,y->getEdgeParent(),qmatrix, verbose);
    else
      {
	dxy = m->second;
	foundxy = true;
      }
    m = distanceMap.find( getPair(x->getNumber(),z->getNumber()) );
    if ( m == distanceMap.end() )
      dxz = mleDistance(alignment,x,x->getEdgeParent(),z,z->getEdgeParent(),qmatrix, verbose);
    else
      {
	dxz = m->second;
	foundxz = true;
      }
    m = distanceMap.find( getPair(y->getNumber(),z->getNumber()) );
    if ( m == distanceMap.end() )
      dyz = mleDistance(alignment,y,y->getEdgeParent(),z,z->getEdgeParent(),qmatrix, verbose);
    else
      {
	dyz = m->second;
	foundyz = true;
      }
    if ( foundxy + foundxz + foundyz == 3 ) //all found, error
      {
	cerr << "Error: three sums found, cannot condition on three sums" << endl;
	//exit(1);
	throw 20;
      }
    if(verbose)
      cout << "foundxy " << foundxy << ", foundyz " << foundyz << ", foundxz " << foundxz << endl;
    double sxy = 0;
    double sxz = 0;
    double syz = 0;
    if ( foundxy )
	sxy = dxy;
    if ( foundxz )
      sxz = dxz;
    if ( foundyz )
      syz = dyz;

    if(verbose)
      cout << "after mleDistance or map, dxy,dxz,dyz " << dxy << ", " << dxz << ", " << dyz << endl;
    //clau: lengthX0, lengthY0, lengthZ0 (from NJ) are starting points for joint N-R:
    //fixit: need to make a function of this
    double lengthX0 = (dxy + dxz - dyz)*0.5;
    double lengthY0 = (dxy + dyz - dxz)*0.5;
    double lengthZ0 = (dxz + dyz - dxy)*0.5;
    if(verbose)
      {
	cout << "lengthX0: " << lengthX0 << endl;
	cout << "lengthY0: " << lengthY0 << endl;
	cout << "lengthZ0: " << lengthZ0 << endl;
      }
    // question: could it be that two or more are negative?
    if ( lengthX0 < 0 )
    {
      if(verbose)
	cerr << "Warning: fix negative edge length in starting point." << endl;
      lengthX0 = 0.0001; //clau: does not work very well starting point of 0
      if (foundxy)
	lengthY0 = dxy - lengthX0;
      if (foundxz)
      	lengthZ0 = dxz - lengthX0;
    }
    if ( lengthY0 < 0 )
    {
      if(verbose)
	cerr << "Warning: fix negative edge length in starting point." << endl;
      lengthY0 = 0.0001;
      if (foundxy)
	lengthX0 = dxy - lengthY0;
      if (foundyz)
      	lengthZ0 = dyz - lengthY0;
    }
    if ( lengthZ0 < 0 )
    {
      if(verbose)
	cerr << "Warning: fix negative edge length in starting point." << endl;
      lengthZ0 = 0.0001;
      if (foundxz)
	lengthX0 = dxz - lengthZ0;
      if (foundyz)
      	lengthY0 = dyz - lengthZ0;
    }
    if(lengthX0 < 0 || lengthY0 < 0 || lengthZ0 < 0)
      {
	cerr << "after fixing negative bl, still found negative" << endl;
	//exit(1);
	throw 20;
      }

    double lx = lengthX0;
    double ly = lengthY0;
    double lz = lengthZ0;

    if(verbose)
      {
	cout << "dxy " << dxy << ", dxz " << dxz << ", dyz " << dyz << endl;
	cout << "sxy " << sxy << ", sxz " << sxz << ", syz " << syz << endl;
	cout << "Starting points: " << endl;
	cout << "lx " << lx << ", ly " << ly << ", lz " << lz << endl;
      }

    mleDistanceJoint(alignment, x, x->getEdgeParent(), y, y->getEdgeParent(), z, z->getEdgeParent(), qmatrix, lx,ly,lz, sxy,sxz,syz,
		     foundxy,foundxz,foundyz,rng, verbose, mvnormal, par3D, par2D);
    //change to pass the parent
    if(verbose)
      {
	cout << "after mleDistanceJoint: " << endl;
	cout << "lx " << lx << ", ly " << ly << ", lz " << lz << endl;
      }

    if ( lx < 0 ) //fixit: make a function of this
    {
      if(verbose)
	cerr << "Warning: fix negative edge length." << endl;
      lx = 0.0;
    }
    if ( ly < 0 )
    {
      if(verbose)
	cerr << "Warning: fix negative edge length." << endl;
      ly = 0.0;
    }
    if ( lz < 0 )
    {
      if(verbose)
	cerr << "Warning: fix negative edge length." << endl;
      lz = 0.0;
    }


    x->getEdgeParent()->setLength( lx );
    y->getEdgeParent()->setLength( ly );

    if ( par==root )
    {
      z->getEdgeParent()->setLength( lz );
      break;
    }
    distanceMap[ getPair(z->getNumber(),par->getNumber()) ] = lz; //why it seems you are only putting lengthZ in distanceMap? dont we need to input X,Y?
                                                                        // is it because dist X,par and Y,par is always one edge only? will it be?
                                                                        // clau: I think yes, we only need lengthZ because that is more than 1 edge
    // for ( int k=0; k<alignment.getNumSites(); ++k )
    // {
    //   if ( par == root)
    // 	par->calculate(k,alignment,NULL,false);
    //   else
    // 	par->calculate(k,alignment,par->getEdgeParent(),false);
    // }

    par->deactivateChild(1);
    par->deactivateChild(0);
  }
  for ( map<pair<int,int>,double>::iterator m=distanceMap.begin(); m!= distanceMap.end(); ++m )
    cout << (*m).first.first << " " << (*m).first.second << " --> " << (*m).second << endl;
}


double Tree::logPriorExp(double mean,ofstream& logwfile)
{
  double logprior = 0;
  for(vector<Edge*>::iterator e=edges.begin();e!=edges.end();e++)
    {
      double length = (*e)->getLength();
      logprior += length;
      //      logwfile << length << ",";
      (*e)->printTable(logwfile);
      logwfile << ",";
    }
  logprior = (1/mean) * logprior;
  return logprior;
}

double Tree::calculateWeight(const Alignment& alignment,QMatrix& qmatrix, double mean, bool verbose, ofstream& logwfile)
{
  double logprior = logPriorExp(mean,logwfile);
  double loglik = calculate(alignment,qmatrix);
  double logdens = getLogdensity();
  if(verbose)
    {
      cout << "Loglik for tree: " << loglik << endl;
      cout << "LogDensity for tree: " <<  logdens << endl;
      cout << "LogPrior for tree: " << logprior << endl;
    }
  double weight = logprior + loglik - logdens;
  if(verbose)
    cout << "LogWeight: " << weight << endl;
  logwfile << loglik << "," << logprior << "," << logdens << "," << weight << endl;
  return weight;
}


double vectorProduct(vector<Vector4d> v)
{
  double sum = 0;
  for ( int i=0; i<4; ++i )
  {
    double product = 1;
    for(  vector<Vector4d>::iterator m = v.begin(); m!=v.end(); m++)
    {
      product *= (*m)(i);
    }
    sum += product;
  }
  return sum;
}

double vectorProduct4D(Vector4d v1, Vector4d v2, Vector4d v3, Vector4d v4)
{
  vector<Vector4d> v;
  v.push_back(v1);
  v.push_back(v2);
  v.push_back(v3);
  v.push_back(v4);
  return vectorProduct(v);
}

// // function just to write grid for logl, dlogl, ddlogl to plot it
// fixit: does not work
// void gridPlotFile(Alignment& alignment,Node* na,Edge* ea,Node* nb,Edge* eb,QMatrix& qmatrix)
// {
//   ofstream currfile;
//   double curr = 0.0001;
//   double curr_logl,curr_dlogl,curr_ddlogl;
//   currfile.open ("currfile.txt");
//   while(curr < 1.0)
//     {
//       partialPathCalculations(curr,alignment,na,ea,nb,eb,qmatrix,curr_logl,curr_dlogl,curr_ddlogl,true);
//       currfile << curr << " " << curr_logl << " " << curr_dlogl << " " << curr_ddlogl << endl;
//       curr = 2 * curr;
//     }
//   currfile.close();
// }

// void gridPlotScreen(Alignment& alignment,Node* na,Edge* ea,Node* nb,Edge* eb,QMatrix& qmatrix)
// {
//   double curr;
//   double curr_logl,curr_dlogl,curr_ddlogl;
//   while(curr < 1.0)
//     {
//       partialPathCalculations(curr,alignment,na,ea,nb,eb,qmatrix,curr_logl,curr_dlogl,curr_ddlogl,true);
//       cout << curr << " " << curr_logl << " " << curr_dlogl << " " << curr_ddlogl << endl;
//       curr = 2 * curr;
//     }
// }
