#ifndef CLADETABLEH
#define CLADETABLEH

#include <vector>
#include "Prime.h"
#include "Clade.h"

using namespace std;

class CladeTable {
  int maxClades;
  vector<int> idx;
  vector<CladePtr> els;

public:

  CladeTable(int m) : maxClades(m), idx(Prime::nextPrime(maxClades*10/6),0), els(1) 
  {} // add dummy first el to els

  CladeTable() : maxClades(1000), idx(Prime::nextPrime(maxClades*10/6),0), els(1) 
  {} // add dummy first el to els

  CladePtr operator[](int  y) { return els[y]; }

  int size() const { return els.size(); }

  vector<CladePtr> &getElements() { return els; }

  inline void push_back(CladePtr cp) { 
    // Bypasses the index. Should only be used for cladeTable[1] (the leaves).
    els.push_back(cp);
  }
    
  bool find(Set &set, unsigned int &k, CladePtr &clade);

  void add(unsigned int k, CladePtr cp);
};

#endif
