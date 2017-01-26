#include <iostream>
#include <iomanip>
#include <vector>
#include <list>
#include <fstream>
#include <sstream>
#include <string>
#include <locale> // isdigit(), tolower()
#include <cctype>
#include <algorithm> // sort(), find for vector<>
#include <random>

#include "Eigen/Core"
#include "Eigen/Eigenvalues"

#include "tree.h"
#include "random.h"
#include "model.h"

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
  f << "numTaxa = " << numTaxa << ", numNodes = " << getNumNodes() << ", numEdges = " << getNumEdges() << endl;
  f << "Nodes:" << endl;
  for ( int i=0; i<getNumNodes(); ++i )
    nodes[i]->print(f);
  f << "Edges:" << endl;
  for ( int i=0; i<getNumEdges(); ++i )
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
// The function relabel() will change node numbers and names to match the sequences.

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
  {
    edges[i]->setNumber(i+1);
    edges[i]->setLength(0);
  }

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
  for ( int k=0; k<nodes.size(); ++k )
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

void Node::makeTreeNumbers(stringstream& s,Edge* parent)
{
  if ( leaf )
  {
    s << number;
  }
  else // internal node
  {
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
      getNeighbor(*e)->makeTreeNumbers(s,*e);
    }
    s << ')';
  }
  if ( parent != NULL )
    s << ":" << setprecision(7) << fixed << parent->getLength();
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

string Tree::makeTreeNumbers()
{
  stringstream s;
  root->makeTreeNumbers(s,NULL);
  s << ";";
  string tree;
  s >> tree;
  return tree;
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
    condProbs(0) = condProbs(3) = 1;
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


void Node::setPattern(int site,const Alignment& alignment,Edge* parent)
{
  pattern.clear();
  if ( leaf )
  {
    char base = alignment.getBase(number,site);
    base = tolower(base);
    pattern = base;
    return;
  }
  // internal node
  for ( vector<Edge*>::iterator e=edges.begin(); e!=edges.end(); ++e )
  {
    if ( (*e) != parent )
    {
      getNeighbor(*e)->setPattern(site,alignment,*e);
      pattern += getNeighbor(*e)->getPattern();
    }
  }
}


void Node::calculateAfterPattern(int site,const Alignment& alignment,Edge* parent)
{
  map<string,pair<double,Vector4d> >::iterator m;

  if ( leaf )
  {
    char base = alignment.getBase(number,site);
    base = tolower(base);
    m = patternToProbMap.find(pattern);
    if ( m == patternToProbMap.end() ) // first time this pattern calculated
    {
      patternToProbMap[ pattern ] = pair<double,Vector4d> (0,translate( base ));
//      cerr << "  " << pattern << " --> " << 0 << "," << translate(base).transpose() << endl;
    }
    return;
  }
  // internal node
  for ( vector<Edge*>::iterator e=edges.begin(); e!=edges.end(); ++e )
  {
    if ( (*e) != parent )
    {
      Node* neighbor = getNeighbor(*e);
      if( (*e) != neighbor->getMapParent() ) //only calculate for children if *e is not mapParent
	neighbor->calculateAfterPattern(site,alignment,*e);//,recurse);
    }
  }
  // if( VERBOSE)
  //   cout << "Calculate after pattern for node: " << number << " for site " << site << endl << flush;
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
	// if(VERBOSE)
	//   {
	//     cout << "pattern not found in map" << endl;
	//     cout << "Neighbor is " << getNeighbor(*e)->getNumber() << endl;
	//     cout << "with condProbPair: " << condProbPair.first << ", " << condProbPair.second.transpose() << endl;
	//   }
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
    // if(VERBOSE)
    //   cout << "  " << pattern << " --> " << scale << "," << tempProb.transpose() << endl;
  }
}


// assumes patternToProbMap is current
void Node::calculate(int site,const Alignment& alignment,Edge* parent)//,bool recurse)
{
  setPattern(site,alignment,parent); // this will clear the pattern of this node, and call setPattern on all the children
  calculateAfterPattern(site,alignment,parent);
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
    root->calculate(k,alignment,NULL);//,true);
    pair<double,Vector4d> condProbPair = root->getProb();
    logLikelihood(k) = condProbPair.first + log( qmatrix.getStationaryP().dot(condProbPair.second) );

    // if(VERBOSE)
    //   {
    // 	Node* x = root->getEdge(0)->getOtherNode(root);
    // 	Node* y = root->getEdge(1)->getOtherNode(root);
    // 	Node* z = root->getEdge(2)->getOtherNode(root);
    // 	pair<double,Vector4d> condProbPairx = x->getProb();
    // 	pair<double,Vector4d> condProbPairy = y->getProb();
    // 	pair<double,Vector4d> condProbPairz = z->getProb();
    // 	cout << "Site: " << k << endl;
    // 	cout << "root getProb: " << condProbPair.first << ", " << condProbPair.second.transpose() << endl;
    // 	cout << "x getProb: " << condProbPairx.first << ", " << condProbPairx.second.transpose() << endl;
    // 	cout << "y getProb: " << condProbPairy.first << ", " << condProbPairy.second.transpose() << endl;
    // 	cout << "z getProb: " << condProbPairz.first << ", " << condProbPairz.second.transpose() << endl;
    // 	cout << logLikelihood(k) << endl;
    //   }
  }
  return logLikelihood.sum();
}

// same function as always, traverse all the tree and does not worry about mapParent Edge
void Node::clearProbMaps(Edge* parent)
{
  for ( vector<Edge*>::iterator e=edges.begin(); e!=edges.end(); ++e )
    {
      if ( (*e) != parent )
	getNeighbor(*e)->clearProbMaps(*e);
    }
//  cout << "clearing prob maps and mapParent for node: " << number << endl << flush;
  patternToProbMap.clear();
  mapParent = NULL;
}

// same function as always, traverse all the tree and does not worry about mapParent Edge
void Tree::clearProbMaps()
{
  root->clearProbMaps(NULL);
}

// checks if mapParent is the same used to create the map, if so, it does not clear the map
void Node::clearProbMapsSmart(Edge* parent)
{
  for ( vector<Edge*>::iterator e=edges.begin(); e!=edges.end(); ++e )
    {
      if ( (*e) != parent )
	{
	  Node* neighbor = getNeighbor(*e);
	  if( (*e) != neighbor->getMapParent() ) //only keep going for children if *e is not mapParent
	    neighbor->clearProbMapsSmart(*e);
	}
    }
  if( mapParent != parent)
    {
//      cout << "clearing prob maps and mapParent for node: " << number << endl << flush;
      patternToProbMap.clear();
      mapParent = NULL;
    }
}

void Tree::clearProbMapsSmart()
{
  root->clearProbMapsSmart(NULL);
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
  uniform_int_distribution<> rint(numTaxa,getNumNodes()-1);
  root = nodes[ rint(rng) ];
  root->randomize(rng,NULL);
  setNodeLevels();
}

// similar to randomize, but uses edges' lengths to
// choose the root
void Tree::randomizeBL(mt19937_64& rng)
{
  // if tree has a long branch, choose one of the nodes attached to it
  Edge* maxEdge = whichMaxBranch();
  if( maxEdge->isTerminal() )
    {
      if( maxEdge->getNode(0)->getLeaf() )
	root = maxEdge->getNode(1);
      else
	root = maxEdge->getNode(0);
    }
  else
    {
      uniform_real_distribution<double> rint(0.0,1.0);
      if(rint(rng) < 0.5)
	      //if(true)
	root = maxEdge->getNode(0);
      else
	root = maxEdge->getNode(1);
    }
  root->randomize(rng,NULL);
  setNodeLevels();
}

Edge* Tree::whichMaxBranch()
{
  double maxbl = 0;
  Edge* maxEdge = *edges.begin();
  // int indmax = 0;
  // int i = 0;
  for ( vector<Edge*>::iterator e=edges.begin(); e!=edges.end(); ++e )
    {
      if( (*e)->getLength() > maxbl )
	{
	  maxbl = (*e)->getLength();
	  //	  indmax = i;
	  maxEdge = *e;
	}
      //i++;
    }
  return maxEdge;
}
// void Tree::partialPathCalculations(double t,Alignment& alignment,Node* na,Edge* ea,Node* nb,Edge* eb,QMatrix& qmatrix,double& logl,double& dlogl,double& ddlogl,bool recurse)
// {
//   Matrix4d P = qmatrix.getTransitionMatrix(t);
//   Matrix4d QP = qmatrix.getQP(t);
//   Matrix4d QQP = qmatrix.getQQP(t);
//   int numSites = alignment.getNumSites();

//   logl = 0;
//   dlogl = 0;
//   ddlogl = 0;
//   for ( int k=0; k<numSites; ++k )
//   {
//     na->calculate(k,alignment,ea,recurse); // set pattern and put probability in map if not already there
//     nb->calculate(k,alignment,eb,recurse);
//     pair<double,Vector4d> pa = na->patternToProbMap[na->getPattern()];
//     pair<double,Vector4d> pb = nb->patternToProbMap[nb->getPattern()];
//     Vector4d va = pa.second;
//     Vector4d vq = qmatrix.getStationaryP();
//     for ( int i=0; i<4; ++i )
//       va(i) *= vq(i);
//     Vector4d vb = pb.second;
//     double f0 = (va.asDiagonal() * P * vb.asDiagonal()).sum();
//     double f1 = (va.asDiagonal() * QP * vb.asDiagonal()).sum();
//     double f2 = (va.asDiagonal() * QQP * vb.asDiagonal()).sum();
//     logl += pa.first + pb.first + log( f0 );
//     dlogl += f1/f0;
//     ddlogl += (f0*f2 - f1*f1)/(f0*f0);
//   }
// }

// double Tree::pathLogLikelihood(double t,Alignment& alignment,Node* na,Edge* ea,Node* nb,Edge* eb,QMatrix& qmatrix,bool recurse)
// {
//   double logl,dlogl,ddlogl;
//   partialPathCalculations(t,alignment,na,ea,nb,eb,qmatrix,logl,dlogl,ddlogl,recurse);
//   return logl;
// }

// double Tree::pathLogLikelihoodDerivative(double t,Alignment& alignment,Node* na,Edge* ea,Node* nb,Edge* eb,QMatrix& qmatrix,bool recurse)
// {
//   double logl,dlogl,ddlogl;
//   partialPathCalculations(t,alignment,na,ea,nb,eb,qmatrix,logl,dlogl,ddlogl,recurse);
//   return dlogl;

// }

// double Tree::pathLogLikelihoodSecondDerivative(double t,Alignment& alignment,Node* na,Edge* ea,Node* nb,Edge* eb,QMatrix& qmatrix, bool recurse)
// {
//   double logl,dlogl,ddlogl;
//   partialPathCalculations(t,alignment,na,ea,nb,eb,qmatrix,logl,dlogl,ddlogl,recurse);
//   return ddlogl;
// }

// void mleError(Node* na,Node* nb,double curr,double prop,double curr_dlogl,double prop_dlogl)
// {
//   cerr << "Error: too many iterations in mleDistance." << endl;
//   cerr << "Nodes " << na->getNumber() << " and " << nb->getNumber() << endl;
//   cerr << "Derivative = " << curr_dlogl << " and " << prop_dlogl << endl;
//   exit(1);
// }

void Edge::mleError(bool& converge)
{
  cerr << "Warning: too many iterations in Edge::mleLength()" << endl;
  converge = false;
}

void Tree::mleError(bool& converge)
{
  cerr << "Warning: too many iterations in Tree::mleLength3D()" << endl;
  converge = false;
}

// Find mle distance from node a to b through a path that uses edges ea and eb.
// Conditon on data in subtrees through other edges.
// Assumes that edge lengths in these subtrees exist
// and that the patternToProbMaps are accurate if edges ea and eb head toward the root.
// double Tree::mleDistance(Alignment& alignment,Node* na,Edge* ea,Node* nb,Edge* eb,QMatrix& qmatrix,double initialLength)
// {
//   bool recurse=true;
//   int iter=0;
// //  double curr = 0.05;
//   double curr = initialLength;
//   // get a decent starting point
//   double curr_logl,curr_dlogl,curr_ddlogl;
//   partialPathCalculations(curr,alignment,na,ea,nb,eb,qmatrix,curr_logl,curr_dlogl,curr_ddlogl,true);
//   double prop = curr;
//   double prop_logl = curr_logl;
//   double prop_dlogl = curr_dlogl;
//   double prop_ddlogl = curr_ddlogl;
//   if ( curr_dlogl > 0 )
//   {
//     do
//     {
//       curr = prop;
//       curr_logl = prop_logl;
//       curr_dlogl = prop_dlogl;
//       curr_ddlogl = prop_ddlogl;
//       prop = 2*curr;
//       partialPathCalculations(prop,alignment,na,ea,nb,eb,qmatrix,prop_logl,prop_dlogl,prop_ddlogl,recurse);
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
//       partialPathCalculations(prop,alignment,na,ea,nb,eb,qmatrix,prop_logl,prop_dlogl,prop_ddlogl,recurse);
//       if ( ++iter > 100 )
//         mleError(na,nb,curr,prop,curr_dlogl,prop_dlogl);
//     } while ( prop_dlogl < 0 );
//   }
//   // switch to protected Newton-Raphson
//   prop = curr - curr_dlogl * (prop - curr) / (prop_dlogl - curr_dlogl);
//   do
//   {
//     curr = prop;
//     partialPathCalculations(curr,alignment,na,ea,nb,eb,qmatrix,curr_logl,curr_dlogl,curr_ddlogl,recurse);
//     if ( ++iter > 100 )
//       mleError(na,nb,curr,prop,curr_dlogl,prop_dlogl);
//     double delta = -curr_dlogl / curr_ddlogl;
//     prop = curr + delta;
//     while ( prop < 0 )
//     {
//       delta = 0.5*delta;
//       prop = curr + delta;
//     }
//     partialPathCalculations(prop,alignment,na,ea,nb,eb,qmatrix,prop_logl,prop_dlogl,prop_ddlogl,recurse);
//     if ( ++iter > 100 )
//       mleError(na,nb,curr,prop,curr_dlogl,prop_dlogl);
//     while ( ( fabs(prop_dlogl) > fabs(curr_dlogl) ) && fabs(curr - prop) >1.0e-8 )
//     {
//       delta = 0.5*delta;
//       prop = curr + delta;
//       partialPathCalculations(prop,alignment,na,ea,nb,eb,qmatrix,prop_logl,prop_dlogl,prop_ddlogl,recurse);
//       if ( ++iter > 100 )
//         mleError(na,nb,curr,prop,curr_dlogl,prop_dlogl);
//     }
//   } while ( fabs(curr - prop) > 1.0e-8);
//   return prop;
// }

void Edge::calculate(double t,Alignment& alignment,QMatrix& qmatrix,double& logl,double& dlogl,double& ddlogl)
{
  Matrix4d P = qmatrix.getTransitionMatrix( t );
  Matrix4d QP = qmatrix.getQP( t );
  Matrix4d QQP = qmatrix.getQQP( t );
  int numSites = alignment.getNumSites();

  if(false)
    {
      cout << "Edge:: calculate on edge " << number << " between nodes " << nodes[0]->getNumber() << " and " << nodes[1]->getNumber() << endl << endl << flush;
      cout << "P =" << endl << P << endl << endl;
      cout << "QP =" << endl << QP << endl << endl;
      cout << "QQP =" << endl << QQP << endl << endl;
    }

  logl = 0;
  dlogl = 0;
  ddlogl = 0;
  for ( int k=0; k<numSites; ++k )
  {
    if(false)
      cout << "k=" << k << endl;
    nodes[0]->calculate(k,alignment,this);//,true); // sets pattern for this site
    nodes[1]->calculate(k,alignment,this);//,true);
    if(false)
      {
	cout << "Nodes: " << nodes[0]->getNumber() << ", " << nodes[1]->getNumber() << endl;
	// cout << nodes[0]->getPattern() << endl;
	// cout << nodes[1]->getPattern() << endl;
      }
    pair<double,Vector4d> pa = nodes[0]->patternToProbMap[nodes[0]->getPattern()];
    pair<double,Vector4d> pb = nodes[1]->patternToProbMap[nodes[1]->getPattern()];
    if(false)
      cout << pa.second.transpose() << " // " << pb.second.transpose() << endl;
    Vector4d va = pa.second;
    Vector4d vq = qmatrix.getStationaryP();
    for ( int i=0; i<4; ++i )
      va(i) *= vq(i);
    Vector4d vb = pb.second;
//    cerr << va.transpose() << " // " << vb.transpose() << endl;
    double f0 = (va.asDiagonal() * P * vb.asDiagonal()).sum();
    double f1 = (va.asDiagonal() * QP * vb.asDiagonal()).sum();
    double f2 = (va.asDiagonal() * QQP * vb.asDiagonal()).sum();
//    cerr << "f0,f1,f2 = " << f0 << ", " << f1 << ", " << f2 << endl;
    logl += pa.first + pb.first + log( f0 );
    dlogl += f1/f0;
    ddlogl += (f0*f2 - f1*f1)/(f0*f0);
    if(false)
      cout << "logl: " << logl << endl;
  }
  nodes[0] -> setMapParent(this); //here we are resetting and traversing every time, maybe we could avoid this
  nodes[1] -> setMapParent(this);
}

double Edge::mleLength(Alignment& alignment,QMatrix& qmatrix,bool& converge)
{
  if(VERBOSE)
    cout << "starting mleLength" << endl;
  int iter=0;
  double curr = length;
  // trying to avoid really bad large initial values
  if ( curr > 1 )
    curr = 1.0;
  
//  cerr << "Edge::mleLength, initial length = " << length << endl;;

  /*
  if ( number == 2 )
  {
    ofstream f("edge2.txt");
    f << "t logl dlogl ddlogl" << endl;
    for ( int i=1; i<=2000; ++i )
    {
      double t,logl,dlogl,ddlogl;
      t = i* 0.001;
      calculate(t,alignment,qmatrix,logl,dlogl,ddlogl);
      f << t << " " << logl << " " << dlogl << " " << ddlogl << endl;
    }
    f.close();
  }
  */
  // try to find two spanning points
  double curr_logl,curr_dlogl,curr_ddlogl;
  calculate(curr,alignment,qmatrix,curr_logl,curr_dlogl,curr_ddlogl);

//  cerr << "curr,logl,dlogl,ddlogl = " << curr << ", " << curr_logl << ", " << curr_dlogl << ", " << curr_ddlogl << endl;

  double prop = curr;
  double prop_logl = curr_logl;
  double prop_dlogl = curr_dlogl;
  double prop_ddlogl = curr_ddlogl;
  if(-TOL < curr_dlogl && curr_dlogl < TOL)
  {
    double prop1 = 0.5*curr;
    double prop1_logl = curr_logl;
    double prop1_dlogl = curr_dlogl;
    double prop1_ddlogl = curr_ddlogl;
    calculate(prop1,alignment,qmatrix,prop1_logl,prop1_dlogl,prop1_ddlogl);
    double prop2 = 2*curr;
    double prop2_logl = curr_logl;
    double prop2_dlogl = curr_dlogl;
    double prop2_ddlogl = curr_ddlogl;
    calculate(prop2,alignment,qmatrix,prop2_logl,prop2_dlogl,prop2_ddlogl);
    do
    {
      if(prop1_dlogl < 0)
      {
	prop1 = 0.5*curr;
	calculate(prop1,alignment,qmatrix,prop1_logl,prop1_dlogl,prop1_ddlogl);
      }
      if(prop2_dlogl > 0)
      {
	prop2 = 2*curr;
	calculate(prop2,alignment,qmatrix,prop2_logl,prop2_dlogl,prop2_ddlogl);
      }
      if ( ++iter > 100 )
      {
        mleError(converge);
	cerr << "Can't find MLE?" << endl;
	prop = 0.5*(prop1+prop2);
	return prop;
      }
    } while ( prop1_dlogl < 0 || prop2_dlogl > 0);
    // here we want to change prop1, prop2, to prop and curr
    prop = prop2;
    prop_logl = prop2_logl;
    prop_dlogl = prop2_dlogl;
    prop_ddlogl = prop2_ddlogl;
    curr = prop1;
    curr_logl = prop1_logl;
    curr_dlogl = prop1_dlogl;
    curr_ddlogl = prop1_ddlogl;
  }
  else if ( curr_dlogl > TOL )
  {
    do
    {
      if(VERBOSE)
	cout << "Positive derivative, curr: " << curr << endl;
      curr = prop;
      curr_logl = prop_logl;
      curr_dlogl = prop_dlogl;
      curr_ddlogl = prop_ddlogl;
      prop = 2*curr;
      calculate(prop,alignment,qmatrix,prop_logl,prop_dlogl,prop_ddlogl);
      if ( ++iter > 100 )
      {
        mleError(converge);
	cerr << "Infinite MLE?" << endl;
	prop = MAX_EDGE_LENGTH;
	return prop;
      }
    } while ( prop_dlogl > TOL && (prop_logl > curr_logl) );
  }
  else
  {
    do
    {
      if(VERBOSE)
	cout << "Negative derivative, curr: " << curr << endl;
      curr = prop;
      curr_logl = prop_logl;
      curr_dlogl = prop_dlogl;
      curr_ddlogl = prop_ddlogl;
      prop = 0.5*curr;
      calculate(prop,alignment,qmatrix,prop_logl,prop_dlogl,prop_ddlogl);
      if ( ++iter > 100 )
      {
        mleError(converge);
	cerr << "Warning: edge length set to minimum." << endl;
	prop = MIN_EDGE_LENGTH;
	return prop;
      }
    } while ( prop_dlogl < TOL && (prop_logl > curr_logl) && prop > MIN_EDGE_LENGTH );
    if ( prop < MIN_EDGE_LENGTH )
    {
      prop = MIN_EDGE_LENGTH;
      return prop;
    }
  }
  double lowerlimit;
  double upperlimit;
  if(prop_dlogl > TOL && curr_dlogl < -TOL)
  {
    lowerlimit = prop;
    upperlimit = curr;
  }
  else if(prop_dlogl < -TOL && curr_dlogl > TOL)
  {
    lowerlimit = curr;
    upperlimit = prop;
  }
  else
  {
    cerr << "Warning: something wrong here. prop and curr should have different sign in derivative" << endl;
    cerr << "prop and dlogl: " << prop << ", " << prop_dlogl << endl;
    cerr << "curr and dlogl: " << curr << ", " << curr_dlogl << endl;
  }
  
  // switch to protected Newton-Raphson
//  cerr << "protected" << endl;
// keep the positive and negative points as limits for the NR, change these limits as NR goes, with the previous point (if the new point is bigger than the previous point, then set the lower limit to be the previous point, and viceversa)
// if trying to jump outside with NR, we do a parabolic interpolation to get a new point
// not doing this now, let's see if the modified NR works
  prop = curr - curr_dlogl * (prop - curr) / (prop_dlogl - curr_dlogl);
  do
  {
    curr = prop;
    calculate(curr,alignment,qmatrix,curr_logl,curr_dlogl,curr_ddlogl);
//    cerr << "curr,logl,dlogl,ddlogl = " << curr << ", " << curr_logl << ", " << curr_dlogl << ", " << curr_ddlogl << endl;

    if ( ++iter > 100 )
    {
      cerr << "Warning: too many iterations in protected Newton-Raphson";
      length = curr;
      mleError(converge);
    }
    double delta = -curr_dlogl / curr_ddlogl;
    prop = curr + delta;
    while ( prop < lowerlimit || prop > upperlimit ) // we don't want to go outside the (lower,upper)
    {
      delta = 0.5*delta;
      prop = curr + delta;
    }
    calculate(prop,alignment,qmatrix,prop_logl,prop_dlogl,prop_ddlogl);
    if ( ++iter > 100 )
      mleError(converge);
    while ( ( fabs(prop_dlogl) > fabs(curr_dlogl) ) && fabs(curr - prop) >1.0e-8 ) // we don't want to go far from maximum
    {
      delta = 0.5*delta;
      prop = curr + delta;
      calculate(prop,alignment,qmatrix,prop_logl,prop_dlogl,prop_ddlogl);
      if ( ++iter > 100 )
        mleError(converge);
    }
    if(prop_dlogl > TOL && curr_dlogl < -TOL) //if close to 0, don't update the limits anymore
    {
      lowerlimit = prop;
      upperlimit = curr;
    }
    else if(prop_dlogl < -TOL && curr_dlogl > TOL)
    {
      lowerlimit = curr;
      upperlimit = prop;
    }
  } while ( fabs(curr - prop) > 1.0e-8);
  return prop;
}

// clear prob maps of both nodes
// do all conditional calculations for both subtrees
// find MLE of edge length conditional on rest of tree
// gamma distributed random length
void Edge::randomLength(Alignment& alignment,QMatrix& qmatrix,mt19937_64& rng,double& logProposalDensity, bool onlyMLE)
{
  nodes[0]->clearProbMapsSmart(this);
  nodes[1]->clearProbMapsSmart(this);

  bool converge;
  // set length to MLE distance
  length = mleLength(alignment,qmatrix,converge); // this is the Edge attribute length
//  cout << "Setting length from mleLength to edge: " << number << endl << flush;
  if ( !onlyMLE )
  {
    if ( length < MIN_EDGE_LENGTH + 1.0e-08 ) // generate from exponential
    {
      double lambda = 1.0 / (MIN_EDGE_LENGTH) ;
      exponential_distribution<double> rexp(lambda);
      length = rexp(rng);
      logProposalDensity += log(lambda) - lambda*length;
    }
    else // generate from gamma
    {
      double logl,dlogl,ddlogl;
      calculate(length,alignment,qmatrix,logl,dlogl,ddlogl);
      double lambda = -1 * length * ddlogl;
      double alpha = length*lambda;
      gamma_distribution<double> rgamma(alpha,1.0 / lambda);
      length = rgamma(rng);
      // if(isnan(length))
      // 	cerr << "found nan bl with alpha: " << alpha << " and lambda: " << lambda << endl;
      logProposalDensity += alpha * log(lambda) - lgamma(alpha) + (alpha-1)*log(length) - lambda*length;
    }
  }
//  cerr << "Setting edge " << number << " length to " << length << endl;
  // still to do: generate a random length

  // recalculate the transition matrices
  calculate(qmatrix);
}

// assumes tree already has some decent starting values for all edges
// can improve efficiency by not calculating with full recursion when not needed, but worry about that later
//
void Node::randomEdges(Alignment& alignment,QMatrix& qmatrix,mt19937_64& rng,Edge* parent,double& logProposalDensity,bool onlyMLE)
{
//  cerr << "node::randomEdges called on" << number << endl;
  // first call on all children
  for ( vector<Edge*>::iterator e=edges.begin(); e != edges.end(); ++e )
  {
    if ( *e != parent )
      getNeighbor(*e)->randomEdges(alignment,qmatrix,rng,*e, logProposalDensity,onlyMLE);
  }
  //  get random length for parent edge
  if ( parent != NULL ) // call edge command on parent
  {
    parent->randomLength(alignment,qmatrix,rng,logProposalDensity,onlyMLE);
  }
}

void Tree::randomEdges(Alignment& alignment,QMatrix& qmatrix,mt19937_64& rng,double& logProposalDensity,bool onlyMLE)
{
  // clear probability maps from all nodes for fresh calculation, and clear mapParent for every node
  clearProbMaps(); //we need this here, otherwise it does not work
  // compute transition matrices for all edges using provisional edge lengths
  for ( vector<Edge*>::iterator e=edges.begin(); e!=edges.end(); ++e )
    (*e)->calculate(qmatrix);
//  cerr << "root = " << root->getNumber() << endl;
  // traverse tree and find MLE's for each edge
  root->randomEdges(alignment,qmatrix,rng,NULL,logProposalDensity,onlyMLE);
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

// will put in nodeList a traversal of left cherries, right cherries and then current node,
// if it is a prunedLeaf
void Node::postorderCherryNodeList(list<Node*>& nodeList,Edge* parent)
{
  for ( vector<Edge*>::iterator e=edges.begin(); e!=edges.end(); ++e )
    if ( *e != parent )
      getNeighbor(*e)->postorderCherryNodeList(nodeList,*e);
  if( isPrunedLeaf(parent) || parent == NULL )
    {
      for ( vector<Edge*>::iterator e=edges.begin(); e!=edges.end(); ++e )
	  if ( *e != parent )
	    {
	      Node* neighbor = getNeighbor(*e);
	      if( neighbor->getLeaf() )
		nodeList.push_back(neighbor);
	    }
      nodeList.push_back(this);
    }
}

void Tree::postorderCherryNodeList(list<Node*>& nodeList)
{
  root->postorderCherryNodeList(nodeList,NULL);
}

void Node::setMapParent(Edge* edge)
{
  if( leaf )
    {
      mapParent = edge;
      return;
    }
  for ( vector<Edge*>::iterator e=edges.begin(); e!=edges.end(); ++e )
  {
    if ( *e != edge )
      getNeighbor(*e)->setMapParent(*e);
  }
  mapParent = edge;
}

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
  activeChildren.clear();
  for ( vector<Edge*>::iterator e=edges.begin(); e!= edges.end(); ++e )
  {
    if ( (*e) != parent )
    {
      Node* n = getNeighbor(*e);
      activeChildren.push_back(n);
      n->setActiveChildrenAndNodeParents(*e);
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

void njCalculateAdjustedMatrix(MatrixXd& dist,list<pair<int,Node*> >& active,MatrixXd& adj,VectorXd& sums)
{
  // copy active part of dist to adj
  int i=0;
  int j=0;
  for ( list<pair<int,Node*> >::iterator pi = active.begin(); pi!= active.end(); ++pi )
  {
    j = 0;
    for ( list<pair<int,Node*> >::iterator pj = active.begin(); pj!= active.end(); ++pj )
    {
      adj(i,j) = dist(pi->first,pj->first);
      ++j;
    }
    ++i;
  }
  // calculate adj
//  VectorXd sums = adj.colwise().sum();
  sums = adj.colwise().sum();
  adj *= (active.size() - 2);
  adj.colwise() -= sums;
  adj.rowwise() -= sums.transpose();
}

void njFindIndicesAtMinimum(MatrixXd& adj,int& mini,int& minj)
{
  mini = 0;
  minj = 1;
  double minAdj = adj(0,1);
  for ( int i=0; i<adj.rows()-1; ++i )
  {
    for ( int j=i+1; j<adj.rows(); ++j)
    {
      if ( adj(i,j) < minAdj )
      {
	mini = i;
	minj = j;
	minAdj = adj(i,j);
      }
    }
  }
}

// rows and columns with *pi with new data and inactivate *pj
void njUpdateDistanceMatrix(MatrixXd& dist,list<pair<int,Node*> >& active,list<pair<int,Node*> >::iterator pi,list<pair<int,Node*> >::iterator pj,Node* n)
{
  double dij = dist(pi->first,pj->first);
  for ( list<pair<int,Node*> >::iterator px=active.begin(); px != active.end(); ++px )
  {
    if ( px == pi || px == pj)
      continue;
    double d = 0.5 * (dist(px->first,pi->first) + dist(px->first,pj->first) - dij);
    if ( d < 0 )
      d = 0;
    dist(px->first,pi->first) = dist(pi->first,px->first) = d;
  }
  *pi = make_pair(pi->first,n);
  active.erase(pj);
}

// make a tree using neighbor-joining algorithm from the distance matrix
Tree::Tree(MatrixXd& dist)
{
  numTaxa = dist.rows();
  int nextNodeIndex = numTaxa + 1;
  int nextEdgeIndex = 0;
  list<pair<int,Node*> > active;

  for ( int i=0; i<numTaxa; ++i )
  {
    nodes.push_back( new Node(i+1,true) );
    active.push_back( make_pair(i,nodes[i]) );
  }
  root = nodes[0];

  while ( active.size() > 3 )
  {
    MatrixXd adj(active.size(),active.size());
    VectorXd sums(active.size());
    njCalculateAdjustedMatrix(dist,active,adj,sums);

    // find indices where adj is at minimum
    int i,j;
    njFindIndicesAtMinimum(adj,i,j);

    // create new node and edges connecting this new node to the two nodes found
    list<pair<int,Node*> >::iterator pi = active.begin();
    for ( int k=0; k<i; ++k )
      ++pi;
    list<pair<int,Node*> >::iterator pj = active.begin();
    for ( int k=0; k<j; ++k )
      ++pj;

    Node* n = new Node(nextNodeIndex++,false);
    Edge* ei = new Edge();
    Edge* ej = new Edge();
    ei->setNumber(nextEdgeIndex++);
    ej->setNumber(nextEdgeIndex++);
    n->addEdge(ei);
    pi->second->addEdge(ei);
    ei->setNodes(pi->second,n);
    n->addEdge(ej);
    pj->second->addEdge(ej);
    ej->setNodes(pj->second,n);
    // set edge lengths
    double d = 0.5*dist(pi->first,pj->first);
    double foo = 0.5*(sums(i) - sums(j))/(active.size()-2);
    ei->setLength( d+foo > 0 ? d+foo : 0.00001 );
    ej->setLength( d-foo > 0 ? d-foo : 0.00001 );
    // add these edges and the nodes to the tree
    nodes.push_back(n);
    edges.push_back(ei);
    edges.push_back(ej);
    // update the distance matrix and active list
    njUpdateDistanceMatrix(dist,active,pi,pj,n);
  }

  // join last three nodes and set lengths
  vector<int> x;
  for ( list<pair<int,Node*> >::iterator p = active.begin(); p != active.end(); ++p )
  {
    x.push_back(p->first);
  }
  vector<double> d(3);
  d[0] = 0.5*(dist(x[0],x[1]) + dist(x[0],x[2]) - dist(x[1],x[2]));
  d[1] = 0.5*(dist(x[0],x[1]) + dist(x[1],x[2]) - dist(x[0],x[2]));
  d[2] = 0.5*(dist(x[0],x[2]) + dist(x[1],x[2]) - dist(x[0],x[1]));
  for ( int i=0; i<3; ++i )
    if ( d[i] < 0 )
      d[i] = 0.00001;

  Node* n = new Node(nextNodeIndex++,false);
  int ii=0;
  for ( list<pair<int,Node*> >::iterator p = active.begin(); p != active.end(); ++p )
  {
    Edge* e = new Edge();
    e->setNumber(nextEdgeIndex++);
    n->addEdge(e);
    p->second->addEdge(e);
    e->setNodes(n,p->second);
    e->setLength(d[ii++]);
    edges.push_back(e);
  }
  nodes.push_back(n);
  root = n;
}

void Tree::reroot(int num)
{
  // find taxon with number int and make its adjacent node the root

  vector<Node*>::iterator p=nodes.begin();
  do
  {
    if ( (*p)->getNumber() == num )
      break;
  } while ( p != nodes.end() );
  if ( p == nodes.end() )
  {
    cerr << "Error: reroot did not find node with number " << num << endl;
    exit(1);
  }
  if ( (*p)->getNumEdges() > 0 )
    root = (*p)->getNeighbor(0);
}

int Node::setMinNumber(Edge* parent)
{
  if ( leaf )
    minNumber = number;
  // recurse on all children
  else
  {
    bool first = true;
    for ( vector<Edge*>::iterator p=edges.begin(); p!= edges.end(); ++p )
    {
      if ( *p != parent )
      {
	int m = getNeighbor(*p)->setMinNumber(*p);
	if ( first )
	{
	  minNumber = m;
	  first = false;
	}
	else
	{
	  if ( m < minNumber )
	    minNumber = m;
	}
      }
    }
  }
  return minNumber;
}

void Node::sortCanonical(Edge* parent)
{
  if ( leaf )
    return;
  vector<pair<int,Edge*> > v;
  if ( parent != NULL )
    v.push_back( make_pair(0,parent) );
  for ( vector<Edge*>::iterator p=edges.begin(); p!= edges.end(); ++p )
  {
    Edge* e = *p;
    if ( e != parent )
    {
      getNeighbor(e)->sortCanonical(e);
      v.push_back( make_pair(getNeighbor(e)->getMinNumber(),e) );
    }
  }
  sort(v.begin(),v.end());
  for ( int i=0; i<edges.size(); ++i )
    edges[i] = (v[i]).second;
}

// use NJ algorithm to set edge lengths cherry by cherry to a tree topology
void Tree::setNJDistances(MatrixXd& dist,mt19937_64& rng)
{
  randomize(rng);
  map<Node*,int> nodeToDistIndexMap;
  // this assumes nodes[0]..nodes[numTaxa-1] are the leaf nodes in order
  for ( int i=0; i<numTaxa; ++i )
  {
    nodeToDistIndexMap[nodes[i]] = i;
  }
  list<Node*> nodeList; // list of all nodes in tree in depth first order
  depthFirstNodeList(nodeList);
  setActiveChildrenAndNodeParents();
  list<Node*>::iterator p=nodeList.begin();

  while ( nodeToDistIndexMap.size() > 3)
  {
    Node* x = *p++;
    Node* y = *p++;
    int xi = nodeToDistIndexMap[x];
    int yi = nodeToDistIndexMap[y];
    double sumx=0;
    double sumy=0;

    for ( map<Node*,int>::iterator m=nodeToDistIndexMap.begin(); m!=nodeToDistIndexMap.end(); ++m)
    {
      sumx += dist(xi,m->second);
      sumy += dist(yi,m->second);
    }
    // set edge lengths
    double dxy = dist(xi,yi);
    double d = 0.5*dxy;
    double foo = 0.5*(sumx - sumy)/(nodeToDistIndexMap.size()-2);
    x->getEdgeParent()->setLength( d+foo > 0 ? d+foo : 0.00001 );
    y->getEdgeParent()->setLength( d-foo > 0 ? d-foo : 0.00001 );
    // update the distance matrix and map
    for ( map<Node*,int>::iterator m=nodeToDistIndexMap.begin(); m!=nodeToDistIndexMap.end(); ++m)
    {
      if ( m->first == x || m->first == y )
	continue;
      d = 0.5 * ( dist(xi,m->second) + dist(yi,m->second) - dxy );
      if ( d < 0 )
	d = 0;
      dist(xi,m->second) = dist(m->second,xi) = d;
    }
    nodeToDistIndexMap[x->getNodeParent()] = xi;
    nodeToDistIndexMap.erase(x);
    nodeToDistIndexMap.erase(y);
  }
  // now, root and three children of root left
  Node* x = *p++;
  Node* y = *p++;
  Node* z = *p++;
  int xi = nodeToDistIndexMap[x];
  int yi = nodeToDistIndexMap[y];
  int zi = nodeToDistIndexMap[z];
  double d = 0.5 * (dist(xi,yi) + dist(xi,yi) - dist(yi,zi));
  x->getEdgeParent()->setLength( d > 0 ? d : 0.00001 );
  d = 0.5 * (dist(xi,yi) + dist(yi,zi) - dist(xi,zi));
  y->getEdgeParent()->setLength( d > 0 ? d : 0.00001 );
  d = 0.5 * (dist(xi,zi) + dist(yi,zi) - dist(xi,yi));
  z->getEdgeParent()->setLength( d > 0 ? d : 0.00001 );
}

// if root node has degree two, remove the node and one of the edges from the tree
// use remaining edge to attach the children.
// set one of these children (a non-leaf) to be the root of the tree

void Tree::unroot()
{
  if ( root == NULL ) // no node is set to be the root
    return;
  if ( root->getNumEdges() < 2 ) // root is a leaf!
    return;
  if ( root->getNumEdges() > 2 ) // already unrooted
  {
    return;
  }
  // root has exactly two children
  Edge* ex = root->getEdge(0);
  Node* x = root->getNeighbor(ex);
  Edge* ey = root->getEdge(1);
  Node* y = root->getNeighbor(ey);
  ex->setNodes(x,y);
  y->deleteEdge(ey);
  y->addEdge(ex);
  vector<Edge*>::iterator p = find(edges.begin(),edges.end(),ey);
  if ( p != edges.end() )
    edges.erase(p);
  vector<Node*>::iterator n = find(nodes.begin(),nodes.end(),root);
  if ( n != nodes.end() )
    nodes.erase(n);

  int deletedNumber = ey->getNumber();

  for ( p=edges.begin(); p!=edges.end(); ++p )
  {
    int num = (*p)->getNumber();
    if ( num > deletedNumber )
      (*p)->setNumber(num-1);
  }

  delete ey;
  delete root;
  if ( !x->getLeaf() )
    root = x;
  else
    root = y;
}

double Tree::logPriorExp(double mean)
{
  double sum = 0;
  double lambda = 1/mean;
  for ( vector<Edge*>::iterator e=edges.begin();e!=edges.end();e++ )
  {
    sum += (*e)->getLength();
  }
  return edges.size()*log(lambda) - lambda*sum;
}


// based on generateBranchLengths in BranchLengths/Code/test/tree.C
void Tree::generateBranchLengths(Alignment& alignment,QMatrix& qmatrix, mt19937_64& rng, double& logdensity, bool jointMLE, double eta, bool weightMean)
{
  // clear probability maps from all nodes for fresh calculation
  clearProbMaps();
  // compute transition matrices for all edges using provisional edge lengths
  for ( vector<Edge*>::iterator e=edges.begin(); e!=edges.end(); ++e )
    (*e)->calculate(qmatrix);
  list<Node*> nodeList;
  postorderCherryNodeList(nodeList);
  //depthFirstNodeList(nodeList);
  setActiveChildrenAndNodeParents(); //need to keep this to have getParentEdge ok
  if(VERBOSE)
    {
      cout << "Starting generateBL. Postorder Node List:";
      for ( list<Node*>::iterator p=nodeList.begin(); p!= nodeList.end(); ++p )
	cout << " " << (*p)->getNumber();
      cout << endl << flush;
    }

  list<Node*>::iterator p=nodeList.begin();
  while ( true )
  {
    Node* x;
    Node* y;
    Node* z;
    Node* par; //clau: parent of x,y
    if( (*p) == root )
    //if ( (*p)->getNodeParent() == root )
      {
	if ( root->getActiveChildrenSize() != 3) //clau: need root to have 3 children
	  {
	    cerr << "yeah, write the general code...." << root->getActiveChildrenSize() << endl;
	    cerr << root->getActiveChild(0)->getNumber() << endl;
	    cerr << (*p)->getNumber() << endl;
	    throw 20;
	  }
	// par = root;
	// x = *p++;
	// y = *p++;
	// z = *p;
	par = root;
	x = par->getEdge(0)->getOtherNode(par);
	y = par->getEdge(1)->getOtherNode(par);
	z = par->getEdge(2)->getOtherNode(par);
      }
    else if ( (*p)->getNodeParent() == root ) //clau: if p parent is root, just skip
      {
    	x = *p++;
    	continue;
      }
    else //clau: p parent not root, and p not root
    {
      x = *p++; //it means take *p and move right
      y = *p++;
      par = x->getNodeParent();
      z = par->getNodeParent();
    }

    // clear prob maps recursively through entire tree
    // there is a smarter way to do this for only part of the tree, depending on order of edges
    // worry about increased efficiency later
    //clearProbMaps(); //we don't need this with new approach

    if ( par==root )
    {
      if(VERBOSE)
	cout << "Setting branch lengths for x,y,z " << x->getNumber() << " " << y->getNumber() << " " << z->getNumber() << endl;
      x-> clearProbMapsSmart(x->getEdgeParent());
      y-> clearProbMapsSmart(y->getEdgeParent());
      z-> clearProbMapsSmart(z->getEdgeParent()); //edge parent of par to go in opposite direction

      // this is computed inside mleLength3D, get rid of this
      // compute probabilities at subtrees x,y,z
      // for ( int k=0; k<alignment.getNumSites(); ++k )
      // 	{
      // 	  x->calculate(k,alignment,x->getEdgeParent());//,true);
      // 	  y->calculate(k,alignment,y->getEdgeParent());//,true);
      // 	  z->calculate(k,alignment,z->getEdgeParent());//,true);
      // 	}

      bool converge = true;
      Vector3d t(x->getEdgeParent()->getLength(),y->getEdgeParent()->getLength(),z->getEdgeParent()->getLength());
      // if(isnan(t[0]) || isnan(t[1]) || isnan(t[2]))
      // 	cerr << "found nan bl in generateBranchLengths, initial BL" << endl;
      if(jointMLE)
	t = mleLength3D(alignment,x,x->getEdgeParent(), y, y->getEdgeParent(), z, z->getEdgeParent(), qmatrix, converge);

      if(!converge)
	cout << "Newton-Raphson did not converge" << endl;
      double prop_logl;
      Vector3d prop_gradient;
      Matrix3d prop_hessian;
      partialPathCalculations3D(t,alignment,x,x->getEdgeParent(),y,y->getEdgeParent(),z,z->getEdgeParent(),qmatrix,prop_logl,prop_gradient,prop_hessian);//,true);
      Vector3d mu = t;
      Matrix3d cov = (-1) * prop_hessian.inverse();
      	// XXX replace mle mu with weighted average of mle and prior mean
	// use PRIOR_MEAN as the mean and sd of the exponential prior edge length
	// use weights proportional to 1 / variance of each estimate

      if(weightMean)
	{
	  double priorVariance = (PRIOR_MEAN * PRIOR_MEAN);
	  mu(0) = ( (PRIOR_MEAN) * cov(0,0) + mu(0)*priorVariance ) / (priorVariance + cov(0,0));
	  mu(1) = ( (PRIOR_MEAN) * cov(1,1) + mu(1)*priorVariance ) / (priorVariance + cov(1,1));
	  mu(2) = ( (PRIOR_MEAN) * cov(2,2) + mu(2)*priorVariance ) / (priorVariance + cov(2,2));
	}
      t = multivariateGamma3D(mu,cov,rng, logdensity,eta);
      //t = multivariateNormal(mu,cov,rng, logdensity, eta);
      // if(isnan(t[0]) || isnan(t[1]) || isnan(t[2]))
      // 	cerr << "found nan bl in generateBranchLengths, after multivariate gamma BL" << endl;

      x->getEdgeParent()->setLength( t[0] );
      y->getEdgeParent()->setLength( t[1] );
      z->getEdgeParent()->setLength( t[2] );
      x->getEdgeParent()->calculate(qmatrix);
      y->getEdgeParent()->calculate(qmatrix);
      z->getEdgeParent()->calculate(qmatrix);
      break;
    }
    else //par not root
      {
	if(VERBOSE)
	  cout << "Setting branch lengths for x,y " << x->getNumber() << " " << y->getNumber() << endl;
	x-> clearProbMapsSmart(x->getEdgeParent());
	y-> clearProbMapsSmart(y->getEdgeParent());
	z-> clearProbMapsSmart(par->getEdgeParent()); //edge parent of par to go in opposite direction

	// this is done inside mleLength2D, get rid of this now
	// compute probabilities at subtrees x,y,z
	// for ( int k=0; k<alignment.getNumSites(); ++k )
	//   {
	//     x->calculate(k,alignment,x->getEdgeParent());//,true);
	//     y->calculate(k,alignment,y->getEdgeParent());//,true);
	//     z->calculate(k,alignment,par->getEdgeParent());//,true); //edge parent of par to go in opposite direction
	//   }

	Vector2d t(x->getEdgeParent()->getLength(),y->getEdgeParent()->getLength());
	bool converge;
	if(jointMLE)
	  t = mleLength2D(alignment,x,x->getEdgeParent(), y, y->getEdgeParent(), z, par->getEdgeParent(), qmatrix, converge);

	double prop_logl;
	Vector2d prop_gradient;
	Matrix2d prop_hessian;
	double lz = par->getEdgeParent()->getLength(); //kept fixed throughout
	partialPathCalculations2D(t,lz,alignment,x,x->getEdgeParent(),y,y->getEdgeParent(),z,par->getEdgeParent(),qmatrix,prop_logl,prop_gradient,prop_hessian);//,true);
	Vector2d mu = t;
	Matrix2d cov = (-1) * prop_hessian.inverse();
	// XXX replace mle mu with weighted average of mle and prior mean
	// use PRIOR_MEAN as the mean and sd of the exponential prior edge length
	// use weights proportional to 1 / variance of each estimate
	if(weightMean)
	  {
	    double priorVariance = (PRIOR_MEAN * PRIOR_MEAN);
	    mu(0) = ( (PRIOR_MEAN) * cov(0,0) + mu(0)*priorVariance ) / (priorVariance + cov(0,0));
	    mu(1) = ( (PRIOR_MEAN) * cov(1,1) + mu(1)*priorVariance ) / (priorVariance + cov(1,1));
	  }
	t = multivariateGamma2D(mu,cov,rng, logdensity,eta);
	//	t = multivariateNormal(mu,cov,rng, logdensity,eta);
	x->getEdgeParent()->setLength( t[0] );
	y->getEdgeParent()->setLength( t[1] );
	x->getEdgeParent()->calculate(qmatrix);
	y->getEdgeParent()->calculate(qmatrix);
	par->clearMapParent();
      }
  }
}

// modified to do 2D MLE instead of 3D
Vector2d Tree::mleLength2D(Alignment& alignment,Node* nx,Edge* ex,Node* ny,Edge* ey,Node* nz,Edge* ez,QMatrix& qmatrix,bool& converge)
{
//  cout << "Entering mleLength2D";
//  cout << " x=" << nx->getNumber() << " y=" << ny->getNumber() << " z=" << nz->getNumber() << endl;
  int iter=0;
  Vector2d curr(ex->getLength(),ey->getLength());
  double lz = ez->getLength(); //kept fixed throughout
  //  cout << "Initial bl curr: " << curr.transpose() << endl;
  double curr_logl;
  Vector2d curr_gradient;
  Matrix2d curr_hessian;
  partialPathCalculations2D(curr,lz,alignment,nx,ex,ny,ey,nz,ez,qmatrix,curr_logl,curr_gradient,curr_hessian);//,true);
  Vector2d prop = curr;
  double prop_logl = curr_logl;
  Vector2d prop_gradient = curr_gradient;
  Matrix2d prop_hessian = curr_hessian;
  Vector2d delta = curr - prop;
  bool keepZero1 = false; //whether to keep that entry at 0
  bool keepZero2 = false;

  // find starting point: we want a good prop
  do
  {
    curr = prop;
    partialPathCalculations2D(curr,lz,alignment,nx,ex,ny,ey,nz,ez,qmatrix,curr_logl,curr_gradient,curr_hessian);//,true);
    // cout << "mleDistance3D Newton-Raphson curr: " << curr.transpose() << endl;
    // cout << "mleDistance3D Newton-Raphson gradient: " << curr_gradient.transpose() << endl;
    // cout << "mleDistance3D Newton-Raphson inverse hessian: " << endl << curr_hessian.inverse() << endl;
    if ( ++iter > 100 )
      mleError(converge);
    delta = curr_hessian.inverse() * curr_gradient;
    if(keepZero1)
      delta[0] = 0;
    if(keepZero2)
      delta[1] = 0;
    prop = curr - delta;
    if(curr[0] > TOL && curr[1] > TOL)
      {
	while ( prop[0] < 0 || prop[1] < 0)
	  {
//	    cout << "found negative with big curr, will shrink delta" << endl;
	    if(prop[0] < 0)
	      delta[0] = 0.5* delta[0];
	    if(prop[1] < 0)
	      delta[1] = 0.5* delta[1];
	    prop = curr - delta;
	  }
      }
    //    cerr << "Delta " << delta.transpose() << endl;
    if(prop[0] < 0 && curr[0] < TOL)
      {
//	cout << "found negative for 1st element with curr small, will set to zero" << endl;
    	prop[0] = MIN_EDGE_LENGTH;
	keepZero1 = true;
      }
    if(prop[1] < 0 && curr[1] < TOL)
      {
//	cout << "found negative for 2nd element with curr small, will set to zero" << endl;
    	prop[1] = MIN_EDGE_LENGTH;
	keepZero2 = true;
      }
//    cout << "prop befofe partialPathCalculations2D: " << prop.transpose() << endl;

    partialPathCalculations2D(prop,lz,alignment,nx,ex,ny,ey,nz,ez,qmatrix,prop_logl,prop_gradient,prop_hessian);//,true);
    if ( ++iter > 100 )
      mleError(converge);
    if(!keepZero1 && !keepZero2)
      {
	while ( prop_gradient.squaredNorm() > curr_gradient.squaredNorm() && delta.squaredNorm() > (TOL*TOL) )
	  {
//	    cout << "found bigger step" << endl;
	    if(keepZero1)
		delta[0] = 0;
	    if(keepZero2)
		delta[1] = 0;
	    delta = 0.5 *delta;
	    prop = curr - delta;
	    partialPathCalculations2D(prop,lz,alignment,nx,ex,ny,ey,nz,ez,qmatrix,prop_logl,prop_gradient,prop_hessian);//,true);
	    if ( ++iter > 100 )
	      mleError(converge);
	  }
      }
  } while ( delta.squaredNorm() > (TOL*TOL) && prop_gradient.squaredNorm() > (TOL*TOL));
//  cout << "Finally converged to" << endl;
//  cout << "Gradient " << endl << prop_gradient.transpose() << endl;
//  cout << "prop: " << prop.transpose() << endl;
  return prop;
}


Vector3d Tree::mleLength3D(Alignment& alignment,Node* nx,Edge* ex,Node* ny,Edge* ey,Node* nz,Edge* ez,QMatrix& qmatrix,bool& converge)
{
//  cout << "Entering mleLength3D";
//  cout << " x=" << nx->getNumber() << " y=" << ny->getNumber() << " z=" << nz->getNumber() << endl;
  int iter=0;
  Vector3d curr(ex->getLength(),ey->getLength(),ez->getLength());
  //  cout << "Initial bl curr: " << curr.transpose() << endl;
  double curr_logl;
  Vector3d curr_gradient;
  Matrix3d curr_hessian;
  partialPathCalculations3D(curr,alignment,nx,ex,ny,ey,nz,ez,qmatrix,curr_logl,curr_gradient,curr_hessian);//,true);
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
    partialPathCalculations3D(curr,alignment,nx,ex,ny,ey,nz,ez,qmatrix,curr_logl,curr_gradient,curr_hessian);//,true);
    // cout << "mleDistance3D Newton-Raphson curr: " << curr.transpose() << endl;
    // cout << "mleDistance3D Newton-Raphson gradient: " << curr_gradient.transpose() << endl;
    // cout << "mleDistance3D Newton-Raphson inverse hessian: " << endl << curr_hessian.inverse() << endl;
    if ( ++iter > 100 )
      mleError(converge);
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
//	    cout << "found negative with big curr, will shrink delta" << endl;
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
//	cout << "found negative for 1st element with curr small, will set to zero" << endl;
    	prop[0] = MIN_EDGE_LENGTH;
	keepZero1 = true;
      }
    if(prop[1] < 0 && curr[1] < TOL)
      {
//	cout << "found negative for 2nd element with curr small, will set to zero" << endl;
    	prop[1] = MIN_EDGE_LENGTH;
	keepZero2 = true;
      }
    if(prop[2] < 0 && curr[2] < TOL)
      {
//	cout << "found negative for 3rd element with curr small, will set to zero" << endl;
    	prop[2] = MIN_EDGE_LENGTH;
	keepZero3 = true;
      }
//    cout << "prop befofe partialPathCalculations3D: " << prop.transpose() << endl;

    partialPathCalculations3D(prop,alignment,nx,ex,ny,ey,nz,ez,qmatrix,prop_logl,prop_gradient,prop_hessian);//,true);
    if ( ++iter > 100 )
      mleError(converge);
    if(!keepZero1 && !keepZero2 && !keepZero3)
      {
	while ( prop_gradient.squaredNorm() > curr_gradient.squaredNorm() && delta.squaredNorm() > (TOL*TOL) )
	  {
//	    cout << "found bigger step" << endl;
	    if(keepZero1)
		delta[0] = 0;
	    if(keepZero2)
		delta[1] = 0;
	    if(keepZero3)
		delta[2] = 0;
	    delta = 0.5 *delta;
	    prop = curr - delta;
	    partialPathCalculations3D(prop,alignment,nx,ex,ny,ey,nz,ez,qmatrix,prop_logl,prop_gradient,prop_hessian);//,true);
	    if ( ++iter > 100 )
	      mleError(converge);
	  }
      }
  } while ( delta.squaredNorm() > (TOL*TOL) && prop_gradient.squaredNorm() > (TOL*TOL));
//  cout << "Finally converged to" << endl;
//  cout << "Gradient " << endl << prop_gradient.transpose() << endl;
//  cout << "prop: " << prop.transpose() << endl;
  return prop;
}

// modified to make t a Vector2d (instead of Vector3d) to keep the length to node z fixed
// gradient and hessian are 2d as well
void Tree::partialPathCalculations2D(Vector2d t, double lz,Alignment& alignment,Node* nx,Edge* ex,Node* ny,Edge* ey,Node* nz,Edge* ez,QMatrix& qmatrix,double& logl,Vector2d& gradient,Matrix2d& hessian)//,bool recurse)
{
  Matrix4d P1 = qmatrix.getTransitionMatrix(t[0]);
  Matrix4d QP1 = qmatrix.getQP(t[0]);
  Matrix4d QQP1 = qmatrix.getQQP(t[0]);
  Matrix4d P2 = qmatrix.getTransitionMatrix(t[1]);
  Matrix4d QP2 = qmatrix.getQP(t[1]);
  Matrix4d QQP2 = qmatrix.getQQP(t[1]);
  Matrix4d P3 = qmatrix.getTransitionMatrix(lz);
  Matrix4d QP3 = qmatrix.getQP(lz);
  Matrix4d QQP3 = qmatrix.getQQP(lz);
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
    nx->calculate(k,alignment,ex);//,recurse); // set pattern and put probability in map if not already there
    ny->calculate(k,alignment,ey);//,recurse);
    nz->calculate(k,alignment,ez);//,recurse);
    pair<double,Vector4d> px = nx->patternToProbMap[nx->getPattern()];
    pair<double,Vector4d> py = ny->patternToProbMap[ny->getPattern()];
    pair<double,Vector4d> pz = nz->patternToProbMap[nz->getPattern()];

    Vector4d S1 = P1 * px.second.asDiagonal() * ones;
    Vector4d S2 = P2 * py.second.asDiagonal() * ones;
    Vector4d S3 = P3 * pz.second.asDiagonal() * ones;
    Vector4d S1pr = QP1 * px.second.asDiagonal() * ones;
    Vector4d S2pr = QP2 * py.second.asDiagonal() * ones;
    Vector4d S1doublepr = QQP1 * px.second.asDiagonal() * ones;
    Vector4d S2doublepr = QQP2 * py.second.asDiagonal() * ones;

    //fixt: make a function of this
    // question: how to make a function to return many things?
    double fk = vectorProduct4D(vq,S1,S2,S3);
    double fkpr1 = vectorProduct4D(vq,S1pr,S2,S3);
    double fkpr2 = vectorProduct4D(vq,S1,S2pr,S3);
    double fkdoublepr11 = vectorProduct4D(vq,S1doublepr,S2,S3);
    double fkdoublepr22 = vectorProduct4D(vq,S1,S2doublepr,S3);
    double fkdoublepr12 = vectorProduct4D(vq,S1pr,S2pr,S3);

    logl += px.first + py.first + pz.first + log( fk );
    dll1 += fkpr1/fk;
    dll2 += fkpr2/fk;
    d2ll_11 += (fk*fkdoublepr11 - fkpr1*fkpr1)/(fk*fk);
    d2ll_12 += (fk*fkdoublepr12 - fkpr1*fkpr2)/(fk*fk);
    d2ll_22 += (fk*fkdoublepr22 - fkpr2*fkpr2)/(fk*fk);
  }
  gradient << dll1, dll2;
  hessian << d2ll_11, d2ll_12, d2ll_12, d2ll_22; //row-wise

  nx -> setMapParent(ex); //here we are resetting and traversing every time, maybe we could avoid this
  ny -> setMapParent(ey);
  nz -> setMapParent(ez);
}


void Tree::partialPathCalculations3D(Vector3d t,Alignment& alignment,Node* nx,Edge* ex,Node* ny,Edge* ey,Node* nz,Edge* ez,QMatrix& qmatrix,double& logl,Vector3d& gradient,Matrix3d& hessian)//,bool recurse)
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
    nx->calculate(k,alignment,ex);//,recurse); // set pattern and put probability in map if not already there
    ny->calculate(k,alignment,ey);//,recurse);
    nz->calculate(k,alignment,ez);//,recurse);
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

  nx -> setMapParent(ex); //here we are resetting and traversing every time, maybe we could avoid this
  ny -> setMapParent(ey);
  nz -> setMapParent(ez);
}


double Tree::vectorProduct(vector<Vector4d> v)
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

double Tree::vectorProduct4D(Vector4d v1, Vector4d v2, Vector4d v3, Vector4d v4)
{
  vector<Vector4d> v;
  v.push_back(v1);
  v.push_back(v2);
  v.push_back(v3);
  v.push_back(v4);
  return vectorProduct(v);
}

// resolves the root so that it has 2 children
// needs to be called after reroot(1)
// will assume the outgroup will be 1
void Tree::makeBinary()
{
  if ( root == NULL ) // no node is set to be the root
    return;
  if ( root->getNumEdges() < 2 ) // root is a leaf!
    return;
  if ( root->getNumEdges() == 2 ) // already rooted
    {
      return;
    }
  if( root -> getNumEdges() > 3 ) // cannot handle more than three children
    exit(1);
  // root has 3 children
  int i;
  Edge* ex;
  Node* x;
  for ( i=0; i<3; ++i )
    {
      ex = root->getEdge(i);
      x = root->getNeighbor(ex);
      if( x->getNumber() == 1 ) //want to find node number 1
	break;
    }
  if(x->getNumber() != 1) //not found
    {
      cerr << "not found node number 1 for resolve root" << endl;
      exit(1);
    }
  root -> deleteEdge(ex);
  Node* newnode = new Node(); //number -1
  newnode -> setLeaf(false);
  Edge* newedge = new Edge(); //no number
  newedge -> setNodes(root,newnode);
  ex -> setNodes(newnode,x);
  newnode -> addEdge(newedge);
  newnode -> addEdge(ex);
  root -> addEdge(newedge);
  edges.push_back(newedge);
  nodes.push_back(newnode);
  root = newnode;
}


// recurse tree for given site and add to the parsimony score
// bases is an integer representative of a set of DNA bases, A=1, C=2, G=4, T=8, the set is the sum of these
// use | for union and & for intersection
void Node::parsimonyScore(Alignment& alignment,Edge* parent,int site,int& score,int& bases)
{
  if ( leaf )
  {
    char b = alignment.getBase(number,site);
    Vector4d cond = translate(tolower(b));
    bases = 0;
    int foo = 1;
    for ( int i=0; i<4; ++i )
    {
      if ( cond(i) > 0 ) // actually, should be a double equal to one
	bases += foo;
      foo *= 2;
    }
    score = 0;
    return;
  }
  // now code for internal nodes
  // assume fully resolved tree
  // case where node is not the root
  if ( parent != NULL )
  {
    Edge* left = NULL;
    Edge* right = NULL;
    for ( vector<Edge*>::iterator e=edges.begin(); e!=edges.end(); ++e )
    {
      if ( *e != parent )
      {
	if ( left == NULL )
	{
	  left = *e;
	}
	else
	{
	  if ( right == NULL )
	  {
	    right = *e;
	  }
	  else
	  {
	    cerr << "Error: found another nonparent edge but left and right already set." << endl;
	  }
	}
      }
    }

    int leftScore = 0,rightScore=0,leftBases=0,rightBases=0;
    getNeighbor(left)->parsimonyScore(alignment,left,site,leftScore,leftBases);
    getNeighbor(right)->parsimonyScore(alignment,right,site,rightScore,rightBases);
    int intersection = leftBases & rightBases;
    score = leftScore + rightScore;
    if ( intersection == 0 )
    {
      bases = leftBases | rightBases;
      ++score;
    }
    else
    {
      bases = intersection;
    }
    return;
  }
  // now case where node is root of the tree and has three children
  //// changed code 8 aug 2016 so it works for root with 2 or more children
  // if ( getNumEdges() != 3 )
  // {
  //   cerr << "Error: was expecting root of the tree with three children.";
  //   exit(1);
  // }
  int score0=0,bases0=0;
  getNeighbor(edges[0])->parsimonyScore(alignment,edges[0],site,score0,bases0);
  score += score0;
  for ( int k=1; k<getNumEdges(); ++k )
  {
    int score1=0,bases1=0;
    getNeighbor(edges[k])->parsimonyScore(alignment,edges[k],site,score1,bases1);
    score += score1;
    int intersection = bases0 & bases1;
    if ( intersection == 0 )
    {
      ++ score;
      bases0 = bases0 | bases1;
    }
    else
    {
      bases0 = intersection;
    }
  }
//  int score0=0,score1=0,score2=0,bases0=0,bases1=0,bases2=0;
//  getNeighbor(edges[0])->parsimonyScore(alignment,edges[0],site,score0,bases0);
//  getNeighbor(edges[1])->parsimonyScore(alignment,edges[1],site,score1,bases1);
//  getNeighbor(edges[2])->parsimonyScore(alignment,edges[2],site,score2,bases2);
//  score = score0 + score1 + score2;
//  int intersection = bases0 & bases1;
  // if ( intersection == 0 )
  // {
  //   ++score;
  //   bases1 = bases0 | bases1;
  // }
  // else
  // {
  //   bases1 = intersection;
  // }
  // intersection = bases1 & bases2;
  // if ( intersection == 0 )
  // {
  //   ++score;
  // }
  // bases = 0;
}

// return parsimony score of the tree
// tree must be rooted
int Tree::parsimonyScore(Alignment& alignment)
{
  int score = 0;
  for ( int site=0; site<alignment.getNumSites(); ++site )
  {
    int siteScore = 0;
    int bases = 0;
    root->parsimonyScore(alignment,NULL,site,siteScore,bases);
    score += siteScore;
  }
  return score;
}

// do 2 MLE passes, and calculate likelihood of tree
// assumes the tree is unrooted, and has set NJ branch lengths
double Tree::logLikelihoodScore(Alignment& alignment, QMatrix& model)
{
  mt19937_64 rng(1234); //just need one as input, never used when onlyMLE=true
  double logBL = 0; //just need one as input, never used when onlyMLE=true
  for ( int i=0; i<2; ++i )
    {
      randomEdges(alignment,model,rng,logBL,true);
    }
  double logLik = calculate(alignment, model);
  return logLik;
}



void Tree::clearMapParent()
{
  for ( vector<Node*>::iterator n=nodes.begin(); n!=nodes.end(); ++n )
    (*n)->clearMapParent();
}

void Node::clearMapParent()
{
  mapParent = NULL;
}


bool Node::isPrunedLeaf(Edge* parent)
{
  if( leaf )
    return false;
  bool prunedLeaf = true;
  for ( vector<Edge*>::iterator e=edges.begin(); e!=edges.end(); ++e )
    if( (*e) != parent )
      {
	Node* neighbor = getNeighbor(*e);
	prunedLeaf &= ( neighbor->getLeaf() || neighbor->isPrunedLeaf(*e) );
      }
  return prunedLeaf;
}

bool Edge::isTerminal()
{
  if( nodes[0]->getLeaf() || nodes[1]->getLeaf() )
    return true;
  else
    return false;
}

void Tree::mcmc(QMatrix& Q,Alignment& alignment,int numGenerations,double scale,mt19937_64& rng, bool burnin)
{
  ofstream treeStream;
  ofstream parStream;
  mcmc(Q,alignment,numGenerations,scale,rng,treeStream,parStream,false, burnin);
}

void Tree::mcmc(QMatrix& Q,Alignment& alignment,int numGenerations,double scale,mt19937_64& rng, ofstream& treeStream, ofstream& parStream, bool burnin)
{
  mcmc(Q,alignment,numGenerations,scale,rng,treeStream,parStream,true, burnin);
}

// for the variance: https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
void Tree::mcmc(QMatrix& Q,Alignment& alignment,int numGenerations,double scale,mt19937_64& rng, ofstream& treeStream, ofstream& parStream, bool printOutput, bool burnin)
{
  clearProbMaps();
  double currLogLikelihood = calculate(alignment,Q);
  double sumAcceptP = 0;
  double sumAcceptS = 0;
  double sumAcceptBL = 0;
  Vector4d avgP = Vector4d::Zero();
  Vector4d avgPold = Vector4d::Zero();
  Vector4d sP = Vector4d::Zero();
  Vector4d prod1 = Vector4d::Zero();
  Vector4d prod2 = Vector4d::Zero();
  Vector4d prod12 = Vector4d::Zero();
  VectorXd avgS = VectorXd::Zero(6);
  VectorXd avgSold = VectorXd::Zero(6);
  VectorXd sS = VectorXd::Zero(6);
  VectorXd prod3 = VectorXd::Zero(6);
  VectorXd prod4 = VectorXd::Zero(6);
  VectorXd prod34 = VectorXd::Zero(6);
  VectorXd avgBL = VectorXd::Zero(getNumEdges());

  cerr << '|';
  for ( int i=0; i<numGenerations; ++i )
  {
    if ( (i+1) % (numGenerations / 100) == 0 )
      cerr << '*';
    if ( (i+1) % (numGenerations / 10) == 0 )
      cerr << '|';
    for( int kk=0; kk<10; ++kk ) //we want to update Q 10 times per once for BL
    {
      Vector4d x = Q.getStationaryP();
      double logProposalRatio = 0;
      Vector4d y = dirichletProposal(x,scale,logProposalRatio,rng);
      QMatrix propQ(y,Q.getSymmetricQP());
      clearProbMaps();
      double propLogLikelihood = calculate(alignment,propQ);
      double acceptProbability = exp(propLogLikelihood - currLogLikelihood + logProposalRatio);
      if ( acceptProbability > 1 )
	acceptProbability = 1;
      sumAcceptP += acceptProbability;
      uniform_real_distribution<double> runif(0,1);
      if ( runif(rng) < acceptProbability )
      {
	Q.resetStationaryP(y);
	currLogLikelihood = propLogLikelihood;
      }
      prod1 = Q.getStationaryP() - avgPold;
      avgP = avgPold + prod1/(i+1);
      prod2 = Q.getStationaryP() - avgP;
      for ( int j=0; j<4; ++j )
	prod12(j) = prod1(j)*prod2(j);
      sP += prod12;
      avgPold = avgP;
      VectorXd xx = Q.getSymmetricQP();
      logProposalRatio = 0;
      VectorXd yy(6);
      yy = dirichletProposal(xx,scale,logProposalRatio,rng);
      propQ.reset(Q.getStationaryP(),yy);
      clearProbMaps();
      propLogLikelihood = calculate(alignment,propQ);
      acceptProbability = exp(propLogLikelihood - currLogLikelihood + logProposalRatio);
      if ( acceptProbability > 1 )
	acceptProbability = 1;
      sumAcceptS += acceptProbability;
      if ( runif(rng) < acceptProbability )
      {
	Q.resetSymmetricQP(yy);
	currLogLikelihood = propLogLikelihood;
      }
      prod3 = Q.getSymmetricQP() - avgSold;
      avgS = avgSold + prod3/(i+1);
      prod4 = Q.getSymmetricQP() - avgS;
      for ( int j=0; j<6; ++j )
	prod34(j) = prod3(j)*prod4(j);
      sS += prod34;
      avgSold = avgS;
    }
    // now with Q, we want to sample branch lengths
    clearProbMaps();
    int j = 0;
    //make tree function
    for ( vector<Edge*>::iterator e=edges.begin(); e!=edges.end(); ++e )
    {
      double xxx = (*e)->getLength();
      logProposalRatio = 0;
      uniform_real_distribution<double> runif(0,1);
      double r = exp(LAMBDA*(runif(rng)-0.5));
      double yyy = xxx * r;
      (*e)->setLength(yyy);
      clearProbMaps();
      propLogLikelihood = calculate(alignment,Q); //make faster later
      acceptProbability = exp(propLogLikelihood - currLogLikelihood + log(r));
      if ( acceptProbability > 1 )
	acceptProbability = 1;
      sumAcceptBL += acceptProbability; // total for all branches
      if ( runif(rng) < acceptProbability )
      {
	(*e)->setLength(yyy);
	currLogLikelihood = propLogLikelihood;
      }
      else
      {
	(*e)->setLength(xxx);
      }
      avgBL(j) += (*e)->getLength();
      j++;
    }
    if ( printOutput )
    {
      parStream << currLogLikelihood << " " << Q.getStationaryP().transpose() << " " << Q.getSymmetricQP().transpose() << endl;
      sortCanonical();
      treeStream << makeTreeNumbers() << endl;
    }
  }
  cout << "stationary acceptance: " << sumAcceptP / numGenerations << endl;
  cout << "symmetric acceptance: " << sumAcceptS / numGenerations << endl;
  cout << "branch length acceptance: " << sumAcceptBL / (getNumEdges() * numGenerations) << endl;
  sP /= numGenerations;
  sS /= numGenerations;
  cout << "avgP: " << avgP.transpose() << endl;
  cout << "sP: " << sP.transpose() << endl;
  cout << "avgS: " << avgS.transpose() << endl;
  cout << "sS: " << sS.transpose() << endl;
  if ( !burnin )
  {
    Q.reset(avgP,avgS);
    Q.setMcmcVarP(sP);
    Q.setMcmcVarQP(sS);
  }
  cout << Q.getStationaryP().transpose() << endl;
  cout << Q.getSymmetricQP().transpose() << endl;
  cout << endl << Q.getQ() << endl << endl;
  if ( !burnin )
  {
    int i = 0;
    for ( vector<Edge*>::iterator e=edges.begin(); e!=edges.end(); ++e )
    {
      (*e)->setLength(avgBL(i)/numGenerations);
      i++;
    }
  }
}

void Tree::setInitialEdgeLengths(double x)
{
  for ( vector<Edge*>::iterator e=edges.begin(); e!=edges.end(); ++e )
  {
    // cerr << "edge number: " << (*e)->getNumber() << endl;
    // cerr << "edge length: " << (*e)->getLength() << endl;
    (*e)->setLength(x);
  }
  // getEdge(0)->setLength(0.03751); // for CCLT
  // getEdge(1)->setLength(0.1185);
  // getEdge(2)->setLength(0.095308);
  // getEdge(3)->setLength(0.052919);
  // getEdge(4)->setLength(0.070793);
}
