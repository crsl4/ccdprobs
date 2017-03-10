#include <sstream>
#include <string>

#include "cladegraph.h"
#include "tree.h"

using namespace std;

CladeNode::CladeNode(dynamic_bitset<unsigned char> taxa)
{
  clade = taxa;
  sumOfLengths = 0;
  if ( clade.count() == 1 )
    leaf = true;
}

void CladeNode::setSubtree(bool isRoot) // set the string
{
  stringstream s;
  if ( leaf )
  {
    s << clade.find_first() + 1 << ':' << meanLength;
  }
  else
  {
    if ( isRoot )
      one->setSubtree(false);
    left->setSubtree(false);
    right->setSubtree(false);
    s << '(';
    if ( isRoot )
      s << one->getSubtree() << ',';
    s << left->getSubtree() << ',' << right->getSubtree() << ')';
    if ( isRoot )
      s << ';';
    else
      s << ':' << meanLength;
  }
  s >> subtree;
}

void CladeNode::setValue(map<dynamic_bitset<unsigned char>,CladeNode*>& cladeMap)
{
  if ( leaf )
  {
    value = meanLength*meanLength;
    return;
  }
  // find best subclades
  value = 0;

  cerr << "CladeNode::setValue: clade = " << clade << endl;
  for( set<pair<dynamic_bitset<unsigned char>,dynamic_bitset<unsigned char> > >::iterator p=subclades.begin(); p!= subclades.end(); ++p )
  {
    CladeNode* n1 = cladeMap[ p->first ];
    CladeNode* n2 = cladeMap[ p->second ];
    double valueSum = n1->getValue() + n2->getValue();
    cerr << "<" << p->first << ": " << n1->getValue() << ", " << p->second << ": " << n2->getValue() << endl;
    if ( valueSum > value )
    {
      left = n1;
      right = n2;
      value = valueSum;
      cerr << "changed left and right" << endl;
    }
  }
  value += meanLength*meanLength;
  cerr << "setValue() clade = " << clade << ", value = " << value << endl;
}

void CladeNode::addPairToSet(dynamic_bitset<unsigned char> clade1,dynamic_bitset<unsigned char> clade2)
{
  cerr << "adding <" << clade1 << "," << clade2 << "> to set for clade " << clade << endl;
  subclades.insert( pair<dynamic_bitset<unsigned char>,dynamic_bitset<unsigned char> >(clade1,clade2) );
}

void Node::processTree(CladeGraph* graph,dynamic_bitset<unsigned char>& clade,Edge* parent)
{
  if ( parent != NULL )
    cerr << parent->getLength() << endl;
  else
    cerr << "root" << endl;

  cerr << "Node::processTree clade.size() = " << clade.size() << endl;
  cerr << "number = " << number << endl;
  
  if ( leaf )
    clade[number-1] = 1;

  cerr << clade << endl;

  bool left = true;
  dynamic_bitset<unsigned char> leftClade(clade.size());
  dynamic_bitset<unsigned char> rightClade(clade.size());
  dynamic_bitset<unsigned char> oneClade(clade.size());

  for ( vector<Edge*>::iterator e=edges.begin(); e!= edges.end(); ++ e)
  {
    if ( (*e) != parent )
    {
      // root
      if ( parent == NULL )
      {
	Node* n = getNeighbor(*e);
	// edge to taxon 1
	if ( n->getLeaf() && (n->getNumber() == 1) )
	{
	  getNeighbor(*e)->processTree(graph,oneClade,*e);
	  clade |= oneClade;
	  continue;
	}
      }
	
      if ( left )
      {
	getNeighbor(*e)->processTree(graph,leftClade,*e);
	clade |= leftClade;
	left = false;
      }
      else
      {
	getNeighbor(*e)->processTree(graph,rightClade,*e);
	clade |= rightClade;
      }
    }
  }

  graph->addCladeToMaps(clade);
  if ( !leaf )
    graph->addPairToCladeNodeSet(clade,leftClade,rightClade);
  if ( parent != NULL )
  {
    graph->addLength(clade,parent->getLength());
  }
}

void Tree::processTree(CladeGraph* graph)
{
  dynamic_bitset<unsigned char> clade(numTaxa);
  root->processTree(graph,clade,NULL);
}

void CladeGraph::processTrees(vector<string> trees)
{
  numTrees = trees.size();
  cerr << "numTrees = " << numTrees << endl;
  if ( numTrees == 0 )
    return;
  Tree* firstTree = new Tree(trees[0]);
  numTaxa = firstTree->getNumTaxa();
  delete firstTree;
  
  for ( vector<string>::iterator p=trees.begin(); p!=trees.end(); ++p )
  {
    cerr << *p << endl;
    Tree* t = new Tree(*p);
    t->reroot(1);
    t->processTree(this);
    delete t;
  }
  // need to set root of the CladeGraph
  dynamic_bitset<unsigned char> allTaxa(numTaxa);
  for ( int i=0; i<numTaxa; ++i )
    allTaxa[i] = 1;
  root = cladeMap[ allTaxa ];
}

void CladeGraph::setMeanLengths()
{
  cerr << "setMeanLengths(), numTrees = " << numTrees << endl;
  for ( map<dynamic_bitset<unsigned char>,CladeNode*>::iterator p=cladeMap.begin(); p!=cladeMap.end(); ++p )
  {
    p->second->setMeanLength( (double)(numTrees) );
    cerr << p->first << " " << p->second->getSumOfLengths() << " " << p->second->getMeanLength() << endl;
  }
}

// iterate through cladeNodes, making sure to do smaller ones before larger ones
void CladeGraph::setValues()
{
  for ( multimap<int,CladeNode*>::iterator p=sizeMap.begin(); p!=sizeMap.end(); ++p )
  {
    p->second->setValue(cladeMap);
//    cerr << "size = " << p->first << ", clade = " << p->second->getClade() << ", value = " << p->second->getValue() << endl;
  }
}

void CladeGraph::findMeanTree(vector<string> trees)
{
  int foo = 0;
  cerr << "Checkpoint " << ++foo << endl;
  processTrees(trees);
  cerr << "Checkpoint " << ++foo << endl;
  setMeanLengths();
  cerr << "Checkpoint " << ++foo << endl;
  setValues();
  cerr << "Checkpoint " << ++foo << endl;
  root->setSubtree(true);
  cerr << "Checkpoint " << ++foo << endl;
  meanTree = root->getSubtree();
  cerr << "Checkpoint " << ++foo << endl;
}

void CladeGraph::addCladeToMaps(dynamic_bitset<unsigned char> clade)
{
  map<dynamic_bitset<unsigned char>,CladeNode*>::iterator p=cladeMap.find(clade);
  if ( p == cladeMap.end() )
  {
    cladeMap[clade] = new CladeNode(clade);
    sizeMap.insert( pair<int,CladeNode*> (clade.count(),cladeMap[clade]) );
  }
}

void CladeGraph::addPairToCladeNodeSet(dynamic_bitset<unsigned char> clade,
				       dynamic_bitset<unsigned char> left,
				       dynamic_bitset<unsigned char> right)
{
  cladeMap[clade]->addPairToSet(left,right);
}

void CladeGraph::addLength(dynamic_bitset<unsigned char> clade,double x)
{
  cerr << "adding " << x << " to " << clade << endl;
  cladeMap[clade]->addLength(x);
}
