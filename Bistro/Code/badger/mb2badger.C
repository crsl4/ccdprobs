// mb2badger.C

// Bret Larget, September 6, 2006

// Input: a MrBayes .t file
// Output: a .top file and a .tre file
//         where each tree is binary and in sorted subtree order
// 
// 25 October 2016
//   Changed code so that parent edge is not always last
//     and do not print binary tree if a node has degree greater than 3

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>

using namespace std; // put this after the #includes

class Node;

class Edge {
public:
  Edge() {}
  Edge(int n,double x) :
    number(n),
    length(x)
  {}
  int getNumber() { return number; }
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
private:
  int number;
  double length;
  Node* nodes[2];
};

class Node {
 public:
  Node() {}
  Node(int n,bool l) : number(n), leaf(l) {}
  int getNumber() const { return number; }
  void setNumber(int x) { number = x; }
  bool isLeaf() { return leaf; }
  void setLeaf(bool x) { leaf = x; }
  int getNumEdges() { return edges.size(); }
  Edge* getEdge(int n) { return edges[n]; }
  void setEdge(Edge* e) { edges.push_back(e); }
  void setEdge(Edge *e,int n) { edges[n] = e; }
  void reorder(Edge*);
  void rootReorder();
  int getMinTaxa() { return minTaxa; }
  int setMinTaxa(int x) { minTaxa = x; return( minTaxa ); }
  int setMinTaxa(Edge*);
  void setNumbers(int&,Edge*);
  void printInfo(ostream&);
  void printTop(ostream&,Edge*);
  void printTree(ostream&,Edge*);
 private:
  int number;
  bool leaf;
  vector<Edge*> edges;
  int minTaxa;
};

class EdgeCompare : public binary_function<Edge*,Edge*,bool> {
public:
  EdgeCompare(Edge* e,Node* n) : parent(e),node(n) {}

  bool operator()(Edge* x,Edge* y) {
    if(x==parent)
      return false;
    if(y==parent)
      return true;
    return (x->getOtherNode(node)->getMinTaxa() < y->getOtherNode(node)->getMinTaxa());
  }

private:
  const Edge* parent;
  Node* node;
};

class RootEdgeCompare : public binary_function<Edge*,Edge*,bool> {
public:
  RootEdgeCompare(Node* r) : root(r) {}

  bool operator()(Edge* x,Edge* y) {
    return (x->getOtherNode(root)->getMinTaxa() < y->getOtherNode(root)->getMinTaxa());
  }

private:
  Node* root;
};

void Node::rootReorder() {
  // only call on the root!
  if(leaf)
    return;

  if(edges.size()>2)
    sort(edges.begin(),edges.end(),RootEdgeCompare(this)); // XXX changed this line to sort all edges of the root
}

void Node::reorder(Edge* par) {
  // sort edges by order of minTaxa of other nodes with parent edge coming last
  // recurse to children
  if(leaf)
    return;
  for(int i=0;i<edges.size();i++) {
    Edge* e = edges[i];
    if(e != par)
      e->getOtherNode(this)->reorder(e);
  }

  sort(edges.begin(),edges.end(),EdgeCompare(par,this));
}

class Tree {
public:
  Tree(string,int);
  ~Tree() {
    nodes.clear();
    edges.clear();
  }
  int getNumTaxa() { return numTaxa; }
  int getNumNodes() { return numNodes; }
  int getNumEdges() { return numEdges; }
  Node* getNode(int n) { return nodes[n]; }
  Edge* getEdge(int n) { return edges[n]; }
  void printTop(ostream&);
  void printTree(ostream&);
  void printInfo(ostream&);
  void reorder(Node*);
  void countTaxa(string,int);
  void readSubtree(istringstream&,Node*,int&,int);
  void setNumbers();
  void setMinTaxa(Node*);
private:
  int numTaxa;
  int numNodes;
  int numEdges;
  vector<Node*> nodes;
  vector<Edge*> edges;
};

bool cmpNodes(Node* x,Node* y) { return( x->getNumber() < y->getNumber() ); }

void Node::setNumbers(int& nextNumber,Edge* parentEdge) {
  // Set numbers for internal nodes for subtrees along non-parent edges.
  // Numbers for leaves should already be set.
  if(leaf)
    return;
  for(vector<Edge*>::iterator e=edges.begin();e!=edges.end();e++)
    if(*e != parentEdge)
      (*e)->getOtherNode(this)->setNumbers(nextNumber,*e);
  setNumber(nextNumber++);
}

void Tree::setNumbers() {
  // Internal nodes are not numbered when tree is read.
  // After reading the tree, we can give numbers to the internal nodes.
  int nextNumber = numTaxa+1;
  for(int i=0;i<nodes[0]->getNumEdges();i++) {
    Edge* e = nodes[0]->getEdge(i);
    e->getOtherNode(nodes[0])->setNumbers(nextNumber,e);
  }
  nodes[0]->setNumber(nextNumber);
}

int Node::setMinTaxa(Edge* parentEdge) {
  // assumes that number is set for node and that it is larger than leaf numbers
  if(leaf) {
    setMinTaxa(number);
    return number;
  }
  int m = number,n;
  for(vector<Edge*>::iterator e=edges.begin();e!=edges.end();e++)
    if(*e != parentEdge) {
      n = (*e)->getOtherNode(this)->setMinTaxa(*e);
      if(n < m)
	m = n;
    }
  setMinTaxa(m);
  return m;
}

void Tree::setMinTaxa(Node* root) {
  // assumes that root has the largest number
  int m = root->getNumber(),n;
  for(int i=0;i<root->getNumEdges();i++) {
    Edge* e = root->getEdge(i);
    n = e->getOtherNode(root)->setMinTaxa(e);
    if(n < m)
      m = n;
  }
  root->setMinTaxa(m);
}

void Tree::reorder(Node* root) {
  // reorder edges in all nodes, last edge for root is still the last
  for(int i=0;i<root->getNumEdges();i++) {
    Edge* e = root->getEdge(i);
    e->getOtherNode(root)->reorder(e);
  }
  root->rootReorder();
}

void Tree::printInfo(ostream& f ) {
  f << "numTaxa = " << numTaxa << ", numNodes = " << numNodes << " , numEdges = " << numEdges << endl;
}

void Tree::countTaxa(string line,int lineNumber)
{
  int left=0,right=0,colon=0,semi=0,comma=0;
  istringstream input(line);
  char c;
  while(input >> c)
    switch(c) {
    case '(': left++;break;
    case ')': right++;break;
    case ':': colon++;break;
    case ';': semi++;break;
    case ',': comma++;break;
    default: break;
    }
  numTaxa = comma+1;
  if(left!=right) {
    cerr << "Error: Parentheses unbalanced on line "<< lineNumber << endl;
    cerr << "Tree string is /" << line << "/" << endl;
    exit(1);
  }
}

void usageError(char *name)
{
  cerr << "Usage: " << name << " <filename>" << endl;
  exit(1);
}

void advance(istringstream& s,char& c,int numLeft,int lineNumber) {
  c = s.peek();
  if(s.eof() || s.fail()) {
    cerr << "Error: Tree on line " << lineNumber
	 << " ends prematurely while reading subtree after left parenthesis " << numLeft << "." << endl;
    exit(1);
  }
  s >> c;
}

void Tree::readSubtree(istringstream& s,Node* parent,int& numLeft,int numLine)
{
  Node* n = new Node();
  nodes.push_back(n);
  numNodes++;
  Edge* e = new Edge();
  edges.push_back(e);
  numEdges++;
  e->setNodes(n,parent);
  n->setEdge(e);
  parent->setEdge(e);
  char c = s.peek();
  if(c=='(') {
    n->setLeaf(false);
    s >> c;
    numLeft++;
    readSubtree(s,n,numLeft,numLine);
    advance(s,c,numLeft,numLine);
    if(c == ',')
      readSubtree(s,n,numLeft,numLine);
    else {
      cerr << "Error: Tree on line " << numLine << " after left parenthesis " << numLeft
	   << " expected comma, but read " << c << "." << endl;
      string a;
      s >> a;
      cerr << "Remainder of line is /" << a << "/." << endl;
    }
    advance(s,c,numLeft,numLine);
    if(c != ')') {
      cerr << "Error: Tree on line " << numLine << " after left parenthesis " << numLeft
	   << " expected right parenthesis, but read " << c << "." << endl;
      string a;
      s >> a;
      cerr << "Remainder of line is /" << a << "/." << endl;
      exit(1);
    }
  }
  else if(isdigit(c)) {
    int num;
    s >> num;
    n->setNumber(num);
    n->setLeaf(true);
    numTaxa++;
  }
  else {
    cerr << "Error: Expected beginning of subtree, but read '" << c << "' on line " << numLine << "." << endl;
    exit(1);
  }
  c = s.peek();
  double length = 1.0;
  if(c == ':') {
    s >> c;
    s >> length;
    if(s.fail()) {
      cerr << "Error: Tree on line " << numLine << " after left parenthesis " << numLeft
	   << " expected an edge length after reading colon." << endl;
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

  nodes.push_back(new Node()); // root of the tree, thinks it is a leaf....
  numNodes = 1;

  int numLeft = 0;
  istringstream s(line);
  char c = s.peek();
  if(c =='(') {
    numLeft++;
    // root is not a leaf
    nodes[0]->setLeaf(false);
    s >> c; // read in left parenthesis
    while(true) {
      readSubtree(s,nodes[0],numLeft,lineNumber);
      s >> c;
      if(s.eof() || s.fail()) {
	cerr << "Error: Tree on line " << lineNumber << " ended while reading in subtree beginning with left parenthesis number "
	     << numLeft << "." << endl;
	exit(1);
      }
      else if(c==',') // read in another subtree to the root.
	continue;
      else if(c==')') { // base subtree is finished.
	c = s.peek();
	if(c==';')
	  s >> c;
	break;
      }
    }
  }
  else if(isdigit(c)) { // base tree may just be a single node
    int n;
    s >> n;
    nodes[0]->setLeaf(true);
    nodes[0]->setNumber(n);
    numTaxa++;
    if(c==';')
      s >> c;
    s >> c;
    if(!(s.eof())) {
      cerr << "Error: Tree on line  " << lineNumber << " begins with a number, but is not a one-taxon tree." << endl;
      exit(1);
    }
  }
  // give numbers to internal nodes
  setNumbers();
  // set minTaxa
  setMinTaxa(nodes[0]);
  // sort nodes
  sort(nodes.begin(),nodes.end(),cmpNodes);
  // reorder edges in nodes
  reorder(nodes[numNodes-1]);
}

void Node::printInfo(ostream& f) {
  int degree = getNumEdges();
  f << setw(3) << number << setw(4) << minTaxa << setw(7) << degree << " ";
  for(int i=0;i<degree;i++) {
    Edge* e = getEdge(i);
    f << e->getOtherNode(this)->getNumber();
    if( i < degree - 1)
      f << ", ";
  }
  f << endl;
}

void Node::printTop(ostream& f,Edge* parent) {
  if(isLeaf()) {
    f << number;
    return;
  }
  int remainingCommas = getNumEdges()-2;
  f << "(";
  for(vector<Edge*>::iterator e=edges.begin();e!=edges.end();e++)
    if(*e != parent) {
      (*e)->getOtherNode(this)->printTop(f,*e);
      if(remainingCommas-- > 0)
	f << ",";
    }
  f << ")";
}

void Tree::printTop(ostream& f) {
  Node* root = nodes[numNodes-1];
  if(numNodes==1) {
    f << root->getNumber() << ";" << endl;
    return;
  }
  f << "(";
  int degree = root->getNumEdges();
//  if(degree==2) { XXX always do this case
    for(int i=0;i<degree;i++) {
      Edge* e = root->getEdge(i);
      e->getOtherNode(root)->printTop(f,e);
      if(i < degree-1)
	f << ",";
    }
//  }
//  else { // degree is 3 or more
    // combine first degree-1 subtrees XXX do not combine subtrees
//    f << "(";
//    for(int i=0;i<degree;i++) { // XXX changed this line
//      Edge* e = root->getEdge(i);
//      e->getOtherNode(root)->printTop(f,e);
//      if(i < degree-1) // XXX changed this line
//	f << ",";
//    }
//    f << "),";
//    Edge* last = root->getEdge(degree-1);
//    last->getOtherNode(root)->printTop(f,last);
//  }
  f << ");" << endl;
}

void Node::printTree(ostream& f,Edge* parent) {
  if(isLeaf()) {
    f << number << ":" << parent->getLength();
    return;
  }
  int remainingCommas = getNumEdges()-2;
  f << "(";
  for(vector<Edge*>::iterator e=edges.begin();e!=edges.end();e++)
    if(*e != parent) {
      (*e)->getOtherNode(this)->printTree(f,*e);
      if(remainingCommas-- > 0)
	f << ",";
    }
  f << ")" << ":" << parent->getLength();
}

void Tree::printTree(ostream& f) {
  Node* root = nodes[numNodes-1];
  if(numNodes==1) {
    f << root->getNumber() << ";" << endl;
    return;
  }
  f << "(";
  int degree = root->getNumEdges();
//  if(degree==2) {
    for(int i=0;i<degree;i++) {
      Edge* e = root->getEdge(i);
      e->getOtherNode(root)->printTree(f,e);
      if(i < degree-1)
	f << ",";
    }
//  }
//  else { // degree is 3 or more
//    // combine first degree-1 subtrees 
//    f << "(";
//    for(int i=0;i<degree-1;i++) {
//      Edge* e = root->getEdge(i);
//      e->getOtherNode(root)->printTree(f,e);
//      if(i < degree-2)
//	f << ",";
//    }
//    f << "):0,";
//    Edge* last = root->getEdge(degree-1);
//    last->getOtherNode(root)->printTree(f,last);
//  }
  f << ");" << endl;
}

int main(int argc, char *argv[])
{
  if(argc != 2)
    usageError(argv[0]);

  char *fileName;
  fileName = argv[1];

  int numTrees=0,numTaxa=0;

  ifstream f(fileName);
  if(f.fail()) {
    cerr << "Error: Could not open file " << fileName << " ." << endl;
    exit(1);
  }

  string fileRoot = (string)(fileName);
  int ext = fileRoot.rfind('.');

  if(ext!=string::npos) // '.' in filename, erase from last '.'
    fileRoot.erase(ext,fileRoot.length() - ext);
    
  string topFile = fileRoot + ".top";
  string treeFile = fileRoot + ".tre";

  ofstream topOut(topFile.c_str());
  ofstream treeOut(treeFile.c_str());

  string line;
  int lineNumber=0;
  while(getline(f,line)) {
    lineNumber++;
    // skip if line is not in format "  tree name = treeRep"
    istringstream s(line);
    string keyTree,name,equalSign;
    s >> keyTree;
    if(keyTree != "tree")
      continue;
    s >> name >> equalSign;
    if(equalSign != "=")
      continue;
    // rest should be a parenthetic representation of a tree.  assume no spaces!
    string treeString;
    s >> treeString;
    Tree tree(treeString,lineNumber);
    numTrees++;
    tree.printTop(topOut);
    tree.printTree(treeOut);
  }
  cout << "Read " << numTrees << " trees." << endl;

  //  EdgeCompare ec(nodes(numNodes
  return(0);
}
