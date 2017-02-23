#include "Tree-summarize.h"

void Tree2::print(ostream& c) {
  if(left) {
    c << "(";
    left->print(c);
    c << ",";
    right->print(c);
    c << ")";
  }
  else
    c << count;
}

void Tree2::printTopology(ostream& c) {
  if(left)
//    if (num!=0)
    if (false)
      c << topologyName(clade->name,num);
    else {
      c << "(";
      left->printTopology(c);
      c << ",";
      right->printTopology(c);
      c << ")";
    }
  else 
    c << count;
}

void Tree2::printWithClade(ostream& c) {
  if(left)
    if (num!=0)
      c << cladeName(clade->name);
    else {
      c << "(";
      left->printWithClade(c);
      c << ",";
      right->printWithClade(c);
      c << ")";
    }
  else 
    c << count;
}

