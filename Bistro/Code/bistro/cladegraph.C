#include <sstream>
#include <string>

#include "cladegraph.h"
#include "tree.h"

using namespace std;

CladeNode::CladeNode(dynamic_bitset<unsigned char> taxa)
{
  clade = taxa;
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
    left->setSubtree(false);
    right->setSubtree(false);
    s << '(' << left->getSubtree() << ',' << right->getSubtree() << ')';
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
  for( set<pair<dynamic_bitset<unsigned char>,dynamic_bitset<unsigned char> > >::iterator p=subclades.begin(); p!= subclades.end(); ++p )
  {
    CladeNode* n1 = cladeMap[ (*p).first];
    CladeNode* n2 = cladeMap[ (*p).second];
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

void CladeGraph::processTree(Tree* tree,map<dynamic_bitset<unsigned char>,CladeNode*>& cladeMap,multimap<int,CladeNode*>& sizeMap)
{
}
    
void CladeGraph::processTrees(vector<string> trees)
{
  numTrees = trees.size();
  if ( numTrees == 0 )
    return;
  Tree* firstTree = new Tree(trees[0]);
  numTaxa = firstTree->getNumTaxa();
  delete firstTree;
  
  for ( vector<string>::iterator p=trees.begin(); p!=trees.end(); ++p )
  {
    Tree* t = new Tree(*p);
    processTree(t,cladeMap,sizeMap);
    delete t;
  }
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
    p->second->setValue(cladeMap);
  }
}

void CladeGraph::findMeanTree(vector<string> trees)
{
  processTrees(trees);
  setMeanLengths();
  setValues();
  root->setSubtree(true);
}

