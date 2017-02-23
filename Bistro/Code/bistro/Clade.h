#ifndef CLADEH
#define CLADEH

using namespace std;

#include <vector>
#include "Set.h"

#if (__GNUC__>=3)
#include <sstream>
#define ISTRSTREAM istringstream
#define OSTRSTREAM ostringstream
#else
#include <strstream.h>
#define ISTRSTREAM istrstream
#define OSTRSTREAM ostrstream
#endif

class Tree2;
typedef Tree2* TreePtr;

class Clade2;
typedef Clade2* CladePtr;

string cladeNameInt(int n, string& str);

string cladeName(int name);

string topologyName(int name, int num);

class Clade2 {
 public:
  Clade2(Set s, int c, int t, double bl, int ln, int n) : set(s), count(c), name(0), numTopologies(t), 
    sumBL(bl), len(ln), firstSplit(0), num(n), 
    hashval(set.hash()) {}
  
    Set set;			 // set of taxa.
    int count;		         // number of times clade appears.
    int name;			 // name of the clade, also used before clade is named to denoted 
    				 // if clade should be named.
    int numTopologies;	         // number of subtree topologies. We count first before collecting them, 
    				 // so topologies.size can be different.
    double sumBL;		 	 // sum of the branch lengths.
    int len;			 // length of the clade.
    int firstSplit;		 // first split in a linked list of splits.
    int num;			 // position of the clade in the clade table.
    vector<TreePtr> topologies;    // list of subtree topologies.
    unsigned int hashval;	         // the hash value of the set.
    double value;		         // least squares distance

    inline Set getSet() const { return set; }

    inline int getCount() const { return count; }
};

#endif
