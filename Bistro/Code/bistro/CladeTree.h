#ifndef CLADETREEH
#define CLADETREEH

using namespace std;

#include "Tree-summarize.h"

class CladeTree {
public:

  CladeTree(TreePtr t, int c) : tree(t), count(c) {}

  TreePtr getTree() const { return tree; }

  int getCount() const { return count; }

  void setCount(int c) { count = c; }

private:
  TreePtr tree;
  int count;
};

#endif
