#include "CladeTable.h"

bool CladeTable::find(Set &set, unsigned int &k, CladePtr &clade) {
  // k will be the position in the CladeTable's index, clade will be the pointer to the clade.
  k = set.hash() % idx.size();
  while(idx[k]) // not empty
    if(set == els[idx[k]]->set) {
      clade = els[idx[k]];
      return true;
    }
    else 
      k = (k + 1) % idx.size();
  return false;
}

void CladeTable::add(unsigned int k, CladePtr cp) {
  // k is the position in the CladeTable's index where the new entry should go.

  if(els.size() >= maxClades) { // Re-hash all the clades.
    maxClades *=2;
    idx.clear();
    idx.resize(Prime::nextPrime(maxClades*10/6),0);

    int t=1;
    for(vector<CladePtr>::iterator p=els.begin()+1;p!=els.end();p++,t++) {
      //k = (*p)->set.hash() % idx.size();
      k = (*p)->hashval % idx.size();
      while(idx[k]) // not empty.
	k = (k + 1) % idx.size();
      idx[k] = t;
    }

    // Reset k to be the position in the new table.
    //k = cp->set.hash() % idx.size();
    k = cp->hashval % idx.size();
    while(idx[k]) // not empty.
      k = (k + 1) % idx.size();
  }

  cp->num = idx[k] = els.size();
  els.push_back(cp);
}
