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
    nb->calculate(k,alignment,eb,recurse); //clau: i think we are doing this twice (before calling mleDistance, and here)
    pair<double,Vector4d> pa = na->patternToProbMap[na->getPattern()];
    pair<double,Vector4d> pb = nb->patternToProbMap[nb->getPattern()];
    Vector4d va = pa.second;
    Vector4d vq = qmatrix.getStationaryP(); //question: do we need to do this inside for each site?
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

    Vector4d S1 = P1 * px.second.asDiagonal() * ones; //question: do i need to save px.second as a variable?
    Vector4d S2 = P2 * py.second.asDiagonal() * ones;
    Vector4d S3 = P3 * pz.second.asDiagonal() * ones;
    Vector4d S1pr = QP1 * px.second.asDiagonal() * ones;
    Vector4d S2pr = QP2 * py.second.asDiagonal() * ones;
    Vector4d S3pr = QP3 * pz.second.asDiagonal() * ones;
    Vector4d S1doublepr = QQP1 * px.second.asDiagonal() * ones;
    Vector4d S2doublepr = QQP2 * py.second.asDiagonal() * ones;
    Vector4d S3doublepr = QQP3 * pz.second.asDiagonal() * ones;

    //    cout << "vq as diagonal " << vq.asDiagonal() << ", S1 as diagonal " << S1.asDiagonal() << endl;

    //question: matrix product error
    // double fk = (vq.asDiagonal() * S1.asDiagonal() * S2.asDiagonal() * S3.asDiagonal()).sum();
    // double fkpr1 = (vq.asDiagonal() * S1pr.asDiagonal() * S2.asDiagonal() * S3.asDiagonal()).sum();
    // double fkpr2 = (vq.asDiagonal() * S1.asDiagonal() * S2pr.asDiagonal() * S3.asDiagonal()).sum();
    // double fkpr3 = (vq.asDiagonal() * S1.asDiagonal() * S2.asDiagonal() * S3pr.asDiagonal()).sum();
    // double fkdoublepr11 = (vq.asDiagonal() * S1doublepr.asDiagonal() * S2.asDiagonal() * S3.asDiagonal()).sum();
    // double fkdoublepr22 = (vq.asDiagonal() * S1.asDiagonal() * S2doublepr.asDiagonal() * S3.asDiagonal()).sum();
    // double fkdoublepr33 = (vq.asDiagonal() * S1.asDiagonal() * S2.asDiagonal() * S3doublepr.asDiagonal()).sum();
    // double fkdoublepr12 = (vq.asDiagonal() * S1pr.asDiagonal() * S2pr.asDiagonal() * S3.asDiagonal()).sum();
    // double fkdoublepr13 = (vq.asDiagonal() * S1pr.asDiagonal() * S2.asDiagonal() * S3pr.asDiagonal()).sum();
    // double fkdoublepr23 = (vq.asDiagonal() * S1.asDiagonal() * S2pr.asDiagonal() * S3pr.asDiagonal()).sum();

    double fk = 0.1 ;
    double fkpr1 = 0.1;
    double fkpr2 = 0.1;
    double fkpr3 = 0.1;
    double fkdoublepr11 = 0.1;
    double fkdoublepr22 = 0.1;
    double fkdoublepr33 = 0.1;
    double fkdoublepr12 = 0.1;
    double fkdoublepr13 = 0.1;
    double fkdoublepr23 = 0.1;


    logl += px.first + py.first + pz.first + log( fk ); //question: scaling correct?
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
  gradient << dll1, dll2, dll3; //question: set vector and matrix like this?
  hessian << d2ll_11, d2ll_12, d2ll_13, d2ll_12, d2ll_22, d2ll_23, d2ll_13, d2ll_23, d2ll_33;
}

// similar to partialPathCalculations3D, t=(t1,t2), sum=s-t1
// warning: need to be careful in the order of t
// gradient and hessian are 2d now
void Tree::partialPathCalculations2D(Vector2d t, double sum,Alignment& alignment,Node* nx,Edge* ex,Node* ny,Edge* ey,Node* nz,Edge* ez,QMatrix& qmatrix,double& logl,Vector2d& gradient,Matrix2d& hessian,bool recurse)
{
  if (sum < t[0])
    {
      cerr << "Sum smaller than summand in partialPathCalculations2D" << endl;
      exit(1);
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

    Vector4d S1 = P1 * px.second.asDiagonal() * ones; //question: do i need to save px.second as a variable?
    Vector4d S2 = P2 * py.second.asDiagonal() * ones;
    Vector4d S3 = P3 * pz.second.asDiagonal() * ones;
    Vector4d S1pr = QP1 * px.second.asDiagonal() * ones;
    Vector4d S2pr = QP2 * py.second.asDiagonal() * ones;
    Vector4d S3pr = (-1) * QP3 * pz.second.asDiagonal() * ones; //question: can we multiply by constant?
    Vector4d S1doublepr = QQP1 * px.second.asDiagonal() * ones;
    Vector4d S2doublepr = QQP2 * py.second.asDiagonal() * ones;
    Vector4d S3doublepr = QQP3 * pz.second.asDiagonal() * ones;

    //    cout << "vq as diagonal " << vq.asDiagonal() << ", S1 as diagonal " << S1.asDiagonal() << endl;

    //question: matrix product error
    // double fk = (vq.asDiagonal() * S1.asDiagonal() * S2.asDiagonal() * S3.asDiagonal()).sum();
    // double fkpr1 = (vq.asDiagonal() * S1pr.asDiagonal() * S2.asDiagonal() * S3.asDiagonal()).sum();
    // double fkpr2 = (vq.asDiagonal() * S1.asDiagonal() * S2pr.asDiagonal() * S3.asDiagonal()).sum();
    // double fkpr3 = (vq.asDiagonal() * S1.asDiagonal() * S2.asDiagonal() * S3pr.asDiagonal()).sum();
    // double fkdoublepr11 = (vq.asDiagonal() * S1doublepr.asDiagonal() * S2.asDiagonal() * S3.asDiagonal()).sum();
    // double fkdoublepr22 = (vq.asDiagonal() * S1.asDiagonal() * S2doublepr.asDiagonal() * S3.asDiagonal()).sum();
    // double fkdoublepr33 = (vq.asDiagonal() * S1.asDiagonal() * S2.asDiagonal() * S3doublepr.asDiagonal()).sum();
    // double fkdoublepr12 = (vq.asDiagonal() * S1pr.asDiagonal() * S2pr.asDiagonal() * S3.asDiagonal()).sum();
    // double fkdoublepr13 = (vq.asDiagonal() * S1pr.asDiagonal() * S2.asDiagonal() * S3pr.asDiagonal()).sum();
    // double fkdoublepr23 = (vq.asDiagonal() * S1.asDiagonal() * S2pr.asDiagonal() * S3pr.asDiagonal()).sum();

    double fk = 0.1 ;
    double fkpr1 = 0.1;
    double fkpr2 = 0.1;
    double fkpr3 = 0.1;
    double fkdoublepr11 = 0.1;
    double fkdoublepr22 = 0.1;
    double fkdoublepr33 = 0.1;
    double fkdoublepr12 = 0.1;
    double fkdoublepr13 = 0.1;
    double fkdoublepr23 = 0.1;


    logl += px.first + py.first + pz.first + log( fk ); //question: scaling correct?
    dll1 += (fkpr1+fkpr3)/fk;
    dll2 += fkpr2/fk;
    d2ll_11 += (fk*(fkdoublepr11+2*fkdoublepr13+fkdoublepr33) - (fkpr1+fkpr3)*(fkpr1+fkpr3))/(fk*fk);
    d2ll_12 += (fk*(fkdoublepr12+fkdoublepr23) - fkpr2*(fkpr1+fkpr3))/(fk*fk);
    d2ll_22 += (fk*fkdoublepr22 - fkpr2*fkpr2)/(fk*fk);
  }
  gradient << dll1, dll2; //question: set vector and matrix like this?
  hessian << d2ll_11, d2ll_12, d2ll_12, d2ll_22;
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
  exit(1);
}
void mleErrorJoint(Node* nx,Node* ny,Node* nz)
{
  cerr << "Error: too many iterations in mleDistanceJoint." << endl;
  cerr << "Nodes " << nx->getNumber() << ", " << ny->getNumber() << ", and " << nz->getNumber() << endl;
  exit(1);
}
// Find mle distance from node a to b through a path that uses edges ea and eb.
// Conditon on data in subtrees through other edges.
// Assumes that edge lengths in these subtrees exist
// and that the patternToProbMaps are accurate if edges ea and eb head toward the root.
double Tree::mleDistance(Alignment& alignment,Node* na,Edge* ea,Node* nb,Edge* eb,QMatrix& qmatrix)
{
  bool recurse=true; //question: not sure this should be true, need to think
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

// Find joint mle distance between nodes nx,ny,nz with edges ex,ey,ez
// Conditon on data in subtrees through other edges.
// Assumes that edge lengths in these subtrees exist
// and that the patternToProbMaps are accurate if edges ex and ey head toward the root.
// it calls different functions depending on the number of sums known
// question: double& to modify inside?
void Tree::mleDistanceJoint(Alignment& alignment,Node* nx,Edge* ex,Node* ny,Edge* ey,Node* nz,Edge* ez,QMatrix& qmatrix, double& lx, double& ly, double& lz, double sxy, double sxz, double syz,mt19937_64& rng)
{
  cout << "Entering mleDistanceJoint" << endl;
  if ( sxy + sxz + syz == 0) //not conditional
    {
      mleDistance3D(alignment,nx,ex,ny,ey,nz,ez,qmatrix, lx, ly, lz,rng);
      return;
    }
  if( sxy > 0 && sxz + syz == 0) //conditional on sxy
    {
      mleDistance2D(alignment, nx,ex,ny,ey,nz,ez,qmatrix,lx,lz,ly,sxy,rng); //warning in order branches
      return;
    }
  if( sxz > 0 && sxy + syz == 0) //conditional on sxz
    {
      mleDistance2D(alignment, nx,ex,ny,ey,nz,ez,qmatrix,lx,ly,lz,sxz,rng); //warning in order branches
      return;
    }
  if( syz > 0 && sxz + sxy == 0) //conditional on syz
    {
      mleDistance2D(alignment, nx,ex,ny,ey,nz,ez,qmatrix,ly,lx,lz,syz,rng); //warning in order branches
      return;
    }
}

// Find joint mle distance between nodes nx,ny,nz with edges ex,ey,ez
// Conditon on data in subtrees through other edges.
// Assumes that edge lengths in these subtrees exist
// and that the patternToProbMaps are accurate if edges ex and ey head toward the root.
void Tree::mleDistance3D(Alignment& alignment,Node* nx,Edge* ex,Node* ny,Edge* ey,Node* nz,Edge* ez,QMatrix& qmatrix, double& lx, double& ly, double& lz, mt19937_64& rng)
{
  cout << "Entering mleDistance3D" << endl;
  bool recurse=true; //question: not sure this should be true, need to think
  int iter=0;
  Vector3d curr(lx,ly,lz);
  double curr_logl;
  Vector3d curr_gradient; //question: declare like this?
  Matrix3d curr_hessian;
  partialPathCalculations3D(curr,alignment,nx,ex,ny,ey,nz,ez,qmatrix,curr_logl,curr_gradient,curr_hessian,true);
  Vector3d prop = curr; //question: how to make equal to curr?
  double prop_logl = curr_logl;
  Vector3d prop_gradient = curr_gradient; //question: how to make vectors/matrices equal?
  Matrix3d prop_hessian = curr_hessian;
  do
  {
    curr = prop;
    partialPathCalculations3D(curr,alignment,nx,ex,ny,ey,nz,ez,qmatrix,curr_logl,curr_gradient,curr_hessian,recurse);
    if ( ++iter > 100 )
      mleErrorJoint(nx,ny,nz);
    //    Vector3d delta = curr_hessian.llt() * curr_gradient; //question: inverse(hessian) * gradient?
    Vector3d delta(0.01,0.01,0.01);
    prop = curr - delta;
    while ( prop[0] < 0 || prop[1] < 0 || prop[2] < 0) //question: better way?
    {
      delta = 0.5*delta;
      prop = curr - delta;
    }
    partialPathCalculations3D(prop,alignment,nx,ex,ny,ey,nz,ez,qmatrix,prop_logl,prop_gradient,prop_hessian,recurse);
    if ( ++iter > 100 )
      mleErrorJoint(nx,ny,nz);
  } while ( fabs(curr[0] - prop[0]) >1.0e-8 || fabs(curr[1] - prop[1]) >1.0e-8 || fabs(curr[2] - prop[2]) >1.0e-8 ); //question: better way?

  //Matrix3d cov = -prop_hessian.llt(); //question: how to take inverse? error here
  //Vector3d bl = multivariateNormal(prop,cov,rng);
  Vector3d bl = prop;
  lx = bl[0]; //question: better input as vector?
  ly = bl[1];
  lz = bl[2];
}

// Find joint mle distance between nodes nx,ny,nz with edges ex,ey,ez
// conditional in one sum: t1,t2,s-t1
// warning: input t1,t2,t3,sum
// Conditon on data in subtrees through other edges.
// Assumes that edge lengths in these subtrees exist
// and that the patternToProbMaps are accurate if edges ex and ey head toward the root.
void Tree::mleDistance2D(Alignment& alignment,Node* nx,Edge* ex,Node* ny,Edge* ey,Node* nz,Edge* ez,QMatrix& qmatrix, double& t1, double& t2, double& t3, double& sum, mt19937_64& rng)
{
  cout << "Entering mleDistance2D" << endl;
  if (sum < t1)
    {
      cerr << "Sum smaller than summand" << endl;
      exit(1);
    }
  bool recurse=true; //question: not sure this should be true, need to think
  int iter=0;
  Vector2d curr(t1, t2);
  double curr_logl;
  Vector2d curr_gradient; //question: declare like this?
  Matrix2d curr_hessian;
  partialPathCalculations2D(curr,sum,alignment,nx,ex,ny,ey,nz,ez,qmatrix,curr_logl,curr_gradient,curr_hessian,true);
  Vector2d prop = curr; //question: how to make equal to curr?
  double prop_logl = curr_logl;
  Vector2d prop_gradient = curr_gradient; //question: how to make vectors/matrices equal?
  Matrix2d prop_hessian = curr_hessian;
  do
  {
    curr = prop;
    partialPathCalculations2D(curr,sum,alignment,nx,ex,ny,ey,nz,ez,qmatrix,curr_logl,curr_gradient,curr_hessian,recurse);
    if ( ++iter > 100 )
      mleErrorJoint(nx,ny,nz);
    //    Vector2d delta = curr_hessian.llt() * curr_gradient; //question: inverse(hessian) * gradient?
    Vector2d delta(0.01,0.01);
    prop = curr - delta;
    while ( prop[0] < 0 || prop[1] < 0) //question: better way?
    {
      delta = 0.5*delta;
      prop = curr - delta;
    }
    partialPathCalculations2D(prop,sum,alignment,nx,ex,ny,ey,nz,ez,qmatrix,prop_logl,prop_gradient,prop_hessian,recurse);
    if ( ++iter > 100 )
      mleErrorJoint(nx,ny,nz);
  } while ( fabs(curr[0] - prop[0]) >1.0e-8 || fabs(curr[1] - prop[1]) >1.0e-8); //question: better way?

  //Matrix2d cov = -prop_hessian.llt(); //question: how to take inverse? error here
  //Vector2d bl = multivariateNormal(prop,cov,rng);
  Vector2d bl = prop;
  t1 = bl[0]; //question: better input as vector?
  t2 = bl[1];
  t3 = sum-bl[0];
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

void Tree::generateBranchLengths(Alignment& alignment,QMatrix& qmatrix, mt19937_64& rng) //clau: added seed
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
      x = *p++; //it means take *p and move right
      y = *p++;
      par = x->getNodeParent();
      z = par->closeRelative();
    }
    // do calculations
    x->calculateEdges(qmatrix); //question: why do we do this? aren't bl meaningless at this point?
    y->calculateEdges(qmatrix);
    z->calculateEdges(qmatrix);
    for ( int k=0; k<alignment.getNumSites(); ++k )
    {
      x->calculate(k,alignment,x->getEdgeParent(),false); //clau: recurse=false
      y->calculate(k,alignment,y->getEdgeParent(),false);
      z->calculate(k,alignment,z->getEdgeParent(),false);
    } //clau: after this, we have in each node the patternProb map, but we dont know how many different patterns, do we? no, we search the pattern in mleDistance

//    cerr << x->getNumber() << " " << y->getNumber() << " " << z->getNumber() << " " << par->getNumber() << endl;
    map<pair<int,int>,double>::iterator m;
    double dxy, dxz, dyz;
    bool foundxy = false; //clau: depending on which bl are found it is condition on 1 or 2 sums
    bool foundyz = false;
    bool foundxz = false;
    m = distanceMap.find( getPair(x->getNumber(),y->getNumber()) );
    if ( m == distanceMap.end() ) //clau: did not find x,y; when do you put them inside distanceMap? at the end here
      dxy = mleDistance(alignment,x,x->getEdgeParent(),y,y->getEdgeParent(),qmatrix);
    else
      {
	dxy = m->second;
	foundxy = true;
      }
    m = distanceMap.find( getPair(x->getNumber(),z->getNumber()) );
    if ( m == distanceMap.end() )
      dxz = mleDistance(alignment,x,x->getEdgeParent(),z,z->getEdgeParent(),qmatrix);
    else
      {
	dxz = m->second;
	foundxz = true;
      }
    m = distanceMap.find( getPair(y->getNumber(),z->getNumber()) );
    if ( m == distanceMap.end() )
      dyz = mleDistance(alignment,y,y->getEdgeParent(),z,z->getEdgeParent(),qmatrix);
    else
      {
	dyz = m->second;
	foundyz = true;
      }
    //clau: lengthX0, lengthY0, lengthZ0 (from NJ) are starting points for joint N-R:
    double lengthX0 = (dxy + dxz - dyz)*0.5;
    if ( lengthX0 < 0 )
    {
      cerr << "Warning: fix negative edge length." << endl;
      lengthX0 = 0.0001; //clau: does not work very well starting point of 0
    }
    double lengthY0 = (dxy + dyz - dxz)*0.5;
    if ( lengthY0 < 0 )
    {
      cerr << "Warning: fix negative edge length." << endl;
      lengthY0 = 0.0001;
    }
    double lengthZ0 = (dxz + dyz - dxy)*0.5;
    if ( lengthZ0 < 0 )
    {
      cerr << "Warning: fix negative edge length." << endl;
      lengthZ0 = 0.0001;
    }
    //cout << "foundxy " << foundxy << ", foundyz " << foundyz << ", foundxz " << foundxz << endl;

    double lx = lengthX0;
    double ly = lengthY0;
    double lz = lengthZ0;
    double sxy = 0;
    double sxz = 0;
    double syz = 0;
    if ( foundxy + foundxz + foundyz == 3 ) //all found, error
      {
	cerr << "Error: three sums found, cannot condition on three sums" << endl;
	exit(1);
      }
    if ( foundxy )
      sxy = dxy;
    if ( foundxz )
      sxz = dxz;
    if ( foundyz )
      syz = dyz;
    cout << "sxy " << sxy << ", sxz " << sxz << ", syz " << syz << endl;
    mleDistanceJoint(alignment, x, x->getEdgeParent(), y, y->getEdgeParent(), z, z->getEdgeParent(), qmatrix, lx,ly,lz, sxy,sxz,syz, rng);

    x->getEdgeParent()->setLength( lx );
    y->getEdgeParent()->setLength( ly );

    if ( par==root )
    {
      z->getEdgeParent()->setLength( lz );
      break;
    }
    distanceMap[ getPair(z->getNumber(),par->getNumber()) ] = lz; //question: why it seems you are only putting lengthZ in distanceMap? dont we need to input X,Y?
                                                                        // is it because dist X,par and Y,par is always one edge only? will it be?
                                                                        // clau: I think yes, we only need lengthZ because that is more than 1 edge
    par->deactivateChild(1);
    par->deactivateChild(0);
  }

  for ( map<pair<int,int>,double>::iterator m=distanceMap.begin(); m!= distanceMap.end(); ++m )
    cout << (*m).first.first << " " << (*m).first.second << " --> " << (*m).second << endl;
}
