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
  else
    leaf = false;
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

void CladeNode::setValue(map<dynamic_bitset<unsigned char>,CladeNode*>& cladeMap,bool isRoot)
{
  if ( isRoot)
  {
    dynamic_bitset<unsigned char> taxonOne(clade.size());
    taxonOne[0] = 1;
    one = cladeMap[ taxonOne ];
  }
  
  if ( leaf )
  {
    value = meanLength*meanLength;
    return;
  }
  // find best subclades
  value = 0;

  for( set<pair<dynamic_bitset<unsigned char>,dynamic_bitset<unsigned char> > >::iterator p=subclades.begin(); p!= subclades.end(); ++p )
  {
    CladeNode* n1 = cladeMap[ p->first ];
    CladeNode* n2 = cladeMap[ p->second ];
    double valueSum = n1->getValue() + n2->getValue();
    if ( valueSum > value )
    {
      left = n1;
      right = n2;
      value = valueSum;
    }
  }
  value += meanLength*meanLength;
}

void CladeNode::addPairToSet(dynamic_bitset<unsigned char> clade1,dynamic_bitset<unsigned char> clade2)
{
  subclades.insert( pair<dynamic_bitset<unsigned char>,dynamic_bitset<unsigned char> >(clade1,clade2) );
}

void Node::processTree(CladeGraph* graph,dynamic_bitset<unsigned char>& clade,Edge* parent)
{
  if ( leaf )
    clade[number-1] = 1;

  bool leftIsNext = true;
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
	
      if ( leftIsNext )
      {
	getNeighbor(*e)->processTree(graph,leftClade,*e);
	clade |= leftClade;
	leftIsNext = false;
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

void CladeGraph::processTrees(vector<string> trees,Alignment& alignment)
{
  numTrees = trees.size();
  if ( numTrees == 0 )
    return;
  Tree* firstTree = new Tree(trees[0],alignment);
  numTaxa = firstTree->getNumTaxa();
  delete firstTree;
  
  for ( vector<string>::iterator p=trees.begin(); p!=trees.end(); ++p )
  {
    Tree* t = new Tree(*p,alignment);
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
  for ( map<dynamic_bitset<unsigned char>,CladeNode*>::iterator p=cladeMap.begin(); p!=cladeMap.end(); ++p )
  {
    p->second->setMeanLength( (double)(numTrees) );
  }
}

// iterate through cladeNodes, making sure to do smaller ones before larger ones
void CladeGraph::setValues()
{
  for ( multimap<int,CladeNode*>::iterator p=sizeMap.begin(); p!=sizeMap.end(); ++p )
  {
    if ( p->first == numTaxa )
      p->second->setValue(cladeMap,true);
    else
      p->second->setValue(cladeMap,false);
  }
}

void CladeGraph::printCladeMap(ostream& f)
{
    for ( map<dynamic_bitset<unsigned char>,CladeNode*>::iterator p=cladeMap.begin(); p!=cladeMap.end(); ++p )
    {
      f << p->first;
      if ( !p->second->getLeaf() )
	f << " " << p->second->getLeft()->getClade() << " " << p->second->getRight()->getClade();
      f << " " << p->second->getValue() << endl;
    }
}

void CladeGraph::printSizeMap(ostream& f)
{
    for ( multimap<int,CladeNode*>::iterator p=sizeMap.begin(); p!=sizeMap.end(); ++p )
      f << setw(2) << p->first << " " << p->second->getClade() << endl;
}

void CladeGraph::printMaps(ostream& f)
{
  f << "Clade Map:" << endl;
  printCladeMap(f);
  f << endl;
  f << "Size Map:" << endl;
  printSizeMap(f);
  f << endl;
}

void CladeGraph::findMeanTree(vector<string> trees,Alignment& alignment)
{
  int foo = 0;
  processTrees(trees,alignment);
  setMeanLengths();
  setValues();
//  printMaps(cerr);
  root->setSubtree(true);
  meanTree = root->getSubtree();
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
  cladeMap[clade]->addLength(x);
}
