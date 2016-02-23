#ifndef TREEH
#define TREEH

using namespace std;

#include <iostream>
#include <iomanip>

#include "Clade.h"

#if (__GNUC__>=3)
#include <sstream>
#define ISTRSTREAM istringstream
#define OSTRSTREAM ostringstream
#else
#include <strstream.h>
#define ISTRSTREAM istrstream
#define OSTRSTREAM ostrstream
#endif

class Tree;

typedef Tree* TreePtr;

class Tree {
 public:

  Tree(TreePtr l, TreePtr r, int c, int n=0) : left(l), right(r), count(c), clade(NULL), num(n) {}

  TreePtr left,right;    // left and right subtrees.
  int count;		 // number of times tree appears, for tree > 1 taxa. Taxon number for 
    		         // trees of 1 taxa.
  CladePtr clade;	 // clade corresponding to the tree.
  int num;   		 // for a topology of a named clade, the number of that topology, also 
    			 // used temporarily as a hash value.

  void print(ostream& c);

  void printTopology(ostream& c);

  void printWithClade(ostream& c);

  inline friend int cmpTrees(const TreePtr& a, const TreePtr& b) {
    if(a==b) 
      return 0;
    else if(a->left==NULL)
      if(b->left==NULL)
	return a->count - b->count;
      else
	return -1;
    else if(b->left==NULL)
      return 1;
    else {
      int c=cmpTrees(a->left,b->left);
      return (c !=0 ? c : cmpTrees(a->right,b->right));
    }
  }

  inline int getCount() const { return count; }
};

#endif
