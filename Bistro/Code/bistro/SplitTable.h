#ifndef SPLITTABLEH
#define SPLITTABLEH

#include <vector>
#include "Prime.h"
#include "Clade.h"

using namespace std;

class SplitTable {

private:

  double Found2, NotFound2, CompFound2, CompNotFound2;
    
  class SplitEntry {
  public:
    CladePtr first;
    CladePtr second;
    int count;
    int nextSplit;

    SplitEntry() : first(NULL), second(NULL), count(0), nextSplit(0) {}

    SplitEntry(CladePtr f, CladePtr s, int c, int n) : first(f), second(s), count(c), nextSplit(n) {}
  };

  int maxSplits;
  vector<int> idx;
  vector<SplitEntry> els;

public:

  // add a dummy element to els, so idx=0 means empty.
  SplitTable() : maxSplits(10000), idx(Prime::nextPrime(maxSplits*10/6),0), els(1) {}

  void resetCounts() { Found2 = NotFound2 = CompFound2 = CompNotFound2 = 0; }

  void printCounts(int count) {
    cerr << setw(10) << count << ": ave. # probes (found) = " << setw(8) 
	 << setprecision(6) << CompFound2/Found2 << ", Ave. # probes (not found) = " 
	 << setw(8) << setprecision(6) << CompNotFound2/NotFound2 << endl;
  }

  //unsigned int hash(CladePtr first, CladePtr second, vector<int> offsets) { return (first->num * (second->num + idx.size())) % idx.size(); }

  inline unsigned int hash(CladePtr first, CladePtr second, vector<int> offsets) { 
    //unsigned int h = ((second->num+1) * 257 + (first->num+1)) % idx.size();
    //unsigned int h = (first->num << 10 + second->num) % idx.size();
    //unsigned int h = (((first->num+1) << 8) * (second->num+1)) % idx.size();
    //unsigned int h = ((idx.size() - first->num) * (second->num+1)) % idx.size();
    //unsigned int h = ((first->num * nthprime(first->len)) << 8 + second->num * nthprime(second->len)) % idx.size();
    //unsigned int h = ((first->num * nthprime(first->len)) << 12 + second->num * nthprime(second->len)) % idx.size();
    //unsigned int h = (first->hashval << 12 + second->hashval) % idx.size();
    //unsigned int h = (first->hashval << 8 + second->hashval) % idx.size();
    //unsigned int h = (first->hashval << 16 + second->hashval) % idx.size();
    //unsigned int h = ((first->num + offsets[first->len]) << 8 + (second->num * offsets[second->len])) % idx.size();
    //unsigned int h = ((second->num + offsets[second->len]) << 8 + (first->num * offsets[first->len])) % idx.size();
    //unsigned int h = ((first->num + offsets[first->len]) << 12 + (second->num * offsets[second->len])) % idx.size();
    //unsigned int h = ((first->num + offsets[first->len]) << 16 + (second->num * offsets[second->len])) % idx.size();
    //unsigned int h = ((first->num + offsets[first->len]) * 257 + (second->num * offsets[second->len])) % idx.size();
    //unsigned int h = ((first->num + offsets[first->len]) * 513 + (second->num * offsets[second->len])) % idx.size();
    //unsigned int h = ((second->num + offsets[second->len]) * 513 + (first->num * offsets[first->len])) % idx.size();
    //unsigned int h = ((second->num + offsets[second->len]) * 257 + (first->num * offsets[first->len])) % idx.size();
    //unsigned int h = ((second->num + offsets[second->len]) * 67 + (first->num * offsets[first->len])) % idx.size();
    //unsigned int h = ((second->num + offsets[second->len]) * 1023 + (first->num * offsets[first->len])) % idx.size();
    //unsigned int h = ((second->num + offsets[second->len]) * 3421 + (first->num * offsets[first->len])) % idx.size();
    //unsigned int h = ((second->num + offsets[second->len]) * 4311 + (first->num * offsets[first->len])) % idx.size();
    //unsigned int h = ((second->num + offsets[second->len]) << 12 + (first->num * offsets[first->len])) % idx.size();
    //***unsigned int h = ((second->num + offsets[second->len]) * 3421 + (first->num * offsets[first->len])) % idx.size();
    //unsigned int h = ((second->num + offsets[second->len]) * 3881 + (first->num * offsets[first->len])) % idx.size();
    //unsigned int h = ((second->num + offsets[second->len]) * 7829 + (first->num * offsets[first->len])) % idx.size();
    //unsigned int h = ((second->num + offsets[second->len]) * 7829 + (first->num * offsets[first->len]) * 3881) % idx.size();
    //***unsigned int h = ((second->num + offsets[second->len]) * 7829 + (first->num * offsets[first->len]) * 313) % idx.size();
    unsigned int h = ((second->num + offsets[second->len]) * offsets[offsets.size()-1] + (first->num * offsets[first->len]))
      % idx.size();

    //cerr << first->num << " " << second->num << " " << h << endl;
    return h;
  }

  inline CladePtr first(int n) { return els[n].first; }

  inline CladePtr second(int n) { return els[n].second; }

  inline int count(int n) { return els[n].count; }

  inline int nextSplit(int n) { return els[n].nextSplit; }

  bool find(CladePtr first, CladePtr second, int &t, unsigned int &k, vector<int> offsets);

  inline void addCount(int t, int count) { els[t].count += count; }

  void add(unsigned int k, CladePtr first, CladePtr second, int count, CladePtr parent, 
	   vector<int> offsets);
};

#endif
