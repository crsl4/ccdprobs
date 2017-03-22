using namespace std;

#include "SplitTable.h"

bool SplitTable::find(CladePtr first, CladePtr second, int &t, unsigned int &k, vector<int> offsets) { 
  // returns t, the entry in the SplitTable, and k, the entry in the SplitTable's index.
  k = hash(first, second, offsets);
  double comps = 0;
  while(comps++,t=idx[k]) // assignment, not equality.
    if(els[t].first == first && els[t].second == second) {
      Found2++;
      CompFound2 += comps;
      return true;
    }
    else 
      k = (k+1) % idx.size();
  NotFound2++;
  CompNotFound2 += comps;
  return false;
}

void SplitTable::add(unsigned int k, CladePtr first, CladePtr second, int count, CladePtr parent, 
		     vector<int> offsets) {
  // k is the position in the SplitTable's index where the new entry should go.
  if(els.size() >=  maxSplits) { // Re-hash all the splits.

    maxSplits *= 2;
    idx.clear();
    idx.resize(Prime::nextPrime(maxSplits*10/6),0);

    cerr << "Re-hashing splits, new size = " << idx.size() << endl;
      
    int t = 1;
    for(vector<SplitEntry>::iterator p=els.begin()+1;p!=els.end();p++,t++) { // skip dummy element
      k = hash(p->first, p->second, offsets);
      while(idx[k]) // not empty
	k = (k + 1) % idx.size();
      idx[k] = t;
    }
    // Find k in the new table.
    k = hash(first, second, offsets);
    while(idx[k])
      k = (k + 1) % idx.size();
  }
  // Actually add the entry to the table.
  idx[k] = els.size();
  els.push_back(SplitEntry(first,second,count,parent->firstSplit));
  parent->firstSplit = idx[k];
}
