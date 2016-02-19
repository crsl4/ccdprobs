#define VERSION 8

/* Summarize for BADGER, Version 8.
 * (c) Copyright 1998, 1999, 2004 by Donald Simon & Bret Larget
 * Department of Mathematics/Computer Science, Duquesne University
 * 5/8/04
 */
using namespace std;

#define VERSION 8

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <utility>
#include <math.h>

#if (__GNUC__>=3)
#include <sstream>
#define ISTRSTREAM istringstream
#define OSTRSTREAM ostringstream
#else
#include <strstream.h>
#define ISTRSTREAM istrstream
#define OSTRSTREAM ostrstream
#endif

#include "prime.h"
#include "votelist.h"

class Set {
  // Simple set class, storing ints in a bit vector. Assumes 32-bit words.
public:
  static void init(int s) { size = (s+30)/31; } // must be called before creating any sets.

  Set() : els(size) {}

  Set& operator= (const Set &a) { 
    for(int i=0;i<size;i++) 
      els[i] = a.els[i]; 
    return *this;
  }

  friend int operator== (const Set &a, const Set &b) { 
    for(int i=0;i<size;i++) 
      if(a.els[i] != b.els[i])
	return 0;
    return 1;
  }

  friend int operator!= (const Set &a, const Set &b) { 
    for(int i=0;i<size;i++) 
      if(a.els[i] != b.els[i])
	return 1;
    return 0;
  }

  friend int operator<(const Set &a, const Set &b) { return a.cmp(b)<0; }

  friend int operator>(const Set &a, const Set &b) { return a.cmp(b)>0; }

  void clear() { 
    for(int i=0;i<size;i++) 
      els[i] = 0; 
  }

  void singleton(const int n) {
    clear();
    els[(n-1)/31] = (1 << (30 - ((n-1) % 31)));
  }

  unsigned int hash() const {
    unsigned int n = els[0];
    for(int i=1;i<size;i++)
      n += els[i];
    return n;
  }

  void myunion(const Set &a, const Set &b) { 
    for(int i=0;i<size;i++) 
      els[i] = a.els[i] | b.els[i];
  }

  void difference(const Set &a, const Set &b) { 
    for(int i=0;i<size;i++) 
      els[i] = a.els[i] & ~b.els[i];
  }

  int subset(const Set &a) const {
    for(int i=0;i<size;i++)
      if(els[i] & ~a.els[i])
	return 0;
    return 1;
  }

  void print(ostream &out) const {
    static unsigned int bit = (1 << 30);
    int first = 1;
    unsigned int n;
    int low = 0, high;

    out << "{";
    for(int i=0;i<size;i++) {
      n = els[i];
      for(int j=1;j<32;j++,n<<=1)
	if(n & bit) {
	  int k = j+31*i;
	  if(low==0)
	    low = high = k;
	  else if (high==k-1)
	    high = k;
	  else {
	    printRange(out,low,high,first);
	    low = high = k;
	    first = 0;
	  }
	}
    }
    if(low!=0)
      printRange(out,low,high,first);
    out << "}";
  }

  void printElements(ostream& c) {
    static unsigned int bit = (1 << 30);
    for(int i=0;i<size;i++) {
      unsigned int n = els[i];
      for(int j=1;j<32;j++,n<<=1)
	if(n & bit)
	  c << setw(12) << " " << setw(3) << j+31*i << endl;
    }
  }

private:
  static int size;
  vector<unsigned int> els;

  int cmp(const Set &a) const {
    int c;
    for(int i=0;i<size;i++) 
      if((c=els[i]-a.els[i])!=0)
	return c;
    return 0;
  }

  void printRange(ostream &out, int low, int high, int first) const
  {
    if(!first)
      out << ",";
    out << low;
    if(high-low>0)
      out << "-" << high;
  }
};

int Set::size;

class Trees {
public:

  Trees() : maxTrees(1000), numTrees(0), totalTrees(0) {}

  void readFiles(const vector<string>& names, int skip) {
    // Read in tree topologies from the files listed in names, skipping the first skip lines of each.

    // Create hash table used to find unique trees (and subtrees).
    int size = FindPrime::nextPrime(maxTrees*10/6);
    vector<TreePtr> sindex(size,(TreePtr)NULL); 

    vector<TreePtr> taxa; // Maps taxon numbers to the corresponding trees.
    int first=1;
    for(vector<string>::const_iterator name=names.begin();name!=names.end();name++) {
      if(*name=="-")
	readFile(cin,*name,skip,sindex,first,taxa);
      else {
	ifstream f((*name).c_str());
	if(!f) {
	  cerr << "Error: Could not open " << *name << endl;
	  exit(1);
	}
	readFile(f,*name,skip,sindex,first,taxa);
	f.close();
      }
    }
    totalTrees = allTrees.size();
  }

  void findClades(double threshold) {
    // Find all the sets that occur in at least threshold (in (0,1]) fraction of the trees.
    // threshold should be the minimum of the threshold for clades to print and the threshold 
    // for named clades.

    int votesNeeded=(int)(ceil(threshold*totalTrees));
    Set set;
    // Create singleton clades for the leaves.
    for(vector<TreePtr>::iterator p=trees[1].begin();p!=trees[1].end();p++) {
      set.singleton((*p)->count);
      (*p)->clade = new Clade(set,0);
      clades[1].push_back((*p)->clade);
    }

    for(int i=2;i<=maxLen;i++) {
      cerr << "Number of trees of size " << i << " = " << trees[i].size() << endl;
      int numCandidates=(maxLen/i*totalTrees)/votesNeeded;
      VoteList<Set> voteList(numCandidates,clades[1][0]->set);
      // Pass 1: Find candidates
      for(vector<TreePtr>::iterator p=trees[i].begin();p!=trees[i].end();p++)
	voteList.voteForCandidate(getClade(*p),(*p)->count);
      // Pass 2: Tally votes for candidates.
      for(vector<TreePtr>::iterator p=trees[i].begin();p!=trees[i].end();p++)
	voteList.checkCandidate(getClade(*p),(*p)->count);
      // Pass 3: Collect winners and link corresponding trees to them.
      vector<CladePtr> id2clade(numCandidates+1); // Convert contenders' ids to clades.
      for(VoteList<Set>::iterator p=voteList.begin();p!=voteList.end();p++) {
	if(p->votes>=votesNeeded) {
	  id2clade[p->id] = new Clade(p->item,p->votes);
	  clades[i].push_back(id2clade[p->id]);
	}
      }
      // Set clade field of trees to the clade if saved, NULL otherwise.
      Contender<Set> q;
      for(vector<TreePtr>::iterator p=trees[i].begin();p!=trees[i].end();p++)
	(*p)->clade = voteList.getContender(getClade(*p),q) && q.votes>=votesNeeded ? id2clade[q.id] : (CladePtr)NULL;
    }
  }

  void findNamedClades(double threshold, int maxTopologies) {
    // Find all named clades, i.e., clades that are maximal, occur in at least threshold fraction of trees,
    // and have <= maxTopologies topologies.
    // findClades must be called first.

    // Eliminate clades with too many topologies.
    for(int i=2;i<=maxLen;i++)
      for(vector<TreePtr>::iterator p=trees[i].begin();p!=trees[i].end();p++)
	if((*p)->clade && ++((*p)->clade->numTopologies)>maxTopologies)
	  (*p)->clade->name = 0;

    // Eliminate clades with too few votes and non-maximal clades. Collect named clades.
    int votesNeeded=(int)(ceil(threshold*totalTrees));
    for(int i=maxLen;i>=2;i--)
      for(vector<CladePtr>::iterator p=clades[i].begin();p!=clades[i].end();p++) {
	if((*p)->count<votesNeeded)
	  (*p)->name = 0;
	if((*p)->name!=0) {
	  namedClades.push_back(*p);
	  for(int j=2;j<i;j++)
	    for(vector<CladePtr>::iterator q=clades[j].begin();q!=clades[j].end();q++)
	      if((*q)->name!=0 && (*q)->set.subset((*p)->set))
		(*q)->name = 0;
	}
      }

    // Reset tree nums to 0 and collect topologies for each clade.
    for(int i=2;i<=maxLen;i++) 
      for(vector<TreePtr>::iterator p=trees[i].begin();p!=trees[i].end();p++) {
	(*p)->num = 0;
	if((*p)->clade && (*p)->clade->name!=0)
	  (*p)->clade->topologies.push_back(*p);
      }

    // Sort named clades, name them, and sort topologies within named clades and name them as well.
    int name=1;
    sort(namedClades.begin(),namedClades.end(),namedCladeLt);
    for(vector<CladePtr>::iterator p=namedClades.begin();p!=namedClades.end();p++,name++) {
      (*p)->name = name;
      sort((*p)->topologies.begin(),(*p)->topologies.end(),treeLt);
      int k=1;
      for(vector<TreePtr>::iterator q=(*p)->topologies.begin();q!=(*p)->topologies.end();q++,k++) 
	(*q)->num = k;
    }
  }
	    
  void printNamedClades(ostream& c) {
    c << "******************** Named clades ********************" << endl << endl;
    Set others=clades[maxLen][0]->set;
    for(vector<CladePtr>::iterator p=namedClades.begin();p!=namedClades.end();p++) {
      others.difference(others,(*p)->set);
      c << setw(6) << (*p)->count << "  " << setw(4) << left << cladeName((*p)->name) << right << " ";
      (*p)->set.print(c);
      c << endl;
      for(vector<TreePtr>::iterator q=(*p)->topologies.begin();q!=(*p)->topologies.end();q++) {
	c << setw(6) << (*q)->count << "  " << setw(4) << left << topologyName((*p)->name,(*q)->num) << right << " ";
	(*q)->print(c);
	c << endl;
      }
      c << endl;
    }
    others.printElements(c);
    c << endl;
  }

  void printClades(ostream& c, double threshold) {
    c << "************************ Common Clades ************************" << endl << endl;
    int votesNeeded = (int)(ceil(threshold*totalTrees));
    vector<CladePtr> pclades;
    for(int i=2;i<=maxLen;i++)
      for(vector<CladePtr>::iterator j=clades[i].begin();j!=clades[i].end();j++)
	if((*j)->count>=votesNeeded)
	  pclades.push_back(*j);
    sort(pclades.begin(),pclades.end(),namedCladeLt);
    for(vector<CladePtr>::iterator j=pclades.begin();j!=pclades.end();j++) {
      c << setw(10) << (*j)->count << " " << setw(6) << setprecision(3)
	<< (double)(*j)->count/(double)totalTrees << " ";
      (*j)->set.print(c);
      c << endl;
    }
    c << endl;
  }

  void printTreeTopologies(ostream& c, int m) {
    // Print out the top m (or all if m=0) tree topologies by count in order from largest to smallest,
    // with named clade-topology number substituted for corresponding named clade topology in the tree.
    // findNamedClades must be called before this.

    c <<  "******************** Tree topologies " << "********************"  << endl << endl;
    c << "Count  Prob.  Cum.  Tree topology" << endl;

    // Find the top m trees and sort them.
    if(m==0 || m>=trees[maxLen].size())
      m = trees[maxLen].size();
    else
      nth_element(trees[maxLen].begin(),trees[maxLen].begin()+m,trees[maxLen].end(),treeLt);
    sort(trees[maxLen].begin(),trees[maxLen].begin()+m,treeLt); 

    // Print out the top m trees.
    int sum=0;
    for(vector<TreePtr>::iterator p=trees[maxLen].begin();p!=trees[maxLen].begin()+m;p++) {
      int ct = (*p)->count;
      sum += ct;
      c << setw(5) << ct << " " << setw(6) << setprecision(3) << ct/(double)totalTrees
	<< " " << setw(6) << sum/(double)totalTrees << " ";
      (*p)->printTopology(c);
      c << endl;
    }
    c << endl;
  }

  void printCladeTreeTopologies(ostream& c, int m) {
    // Print out the top m (or all if m=0) tree topologies by count in order from largest to smallest,
    // with a named clade substituted for the corresponding named clade in the tree.
    // findNamedClades must be called before this.

    c << "******************** Clade tree topologies " << "********************" << endl << endl;
    c << "Count  Prob.  Cum.  Tree topology" << endl;

    // Find the counts of all trees with the same topology when named clades are considered the same despite topology.
    vector<CladeTree> cladeTrees;
    for(vector<TreePtr>::iterator q=trees[maxLen].begin();q!=trees[maxLen].end();q++)
      cladeTrees.push_back(CladeTree(*q,(*q)->count));

    sort(cladeTrees.begin(),cladeTrees.end(),cladeTreeLt);

    vector<CladeTree>::iterator last=cladeTrees.begin();
    for(vector<CladeTree>::iterator q=last+1;q!=cladeTrees.end();q++) {
      if(cladeTreeEq(*last,*q))
	last->setCount(last->getCount()+q->getCount());
      else 
	*(++last) = *q;
    }
    last++;
    int size=last-cladeTrees.begin();

    // Find the top m trees by size and sort.
    if(m==0 || m>=size)
      m = size;
    else
      nth_element(cladeTrees.begin(),cladeTrees.begin()+m,last,cladeTreeLess);
    sort(cladeTrees.begin(),cladeTrees.begin()+m,cladeTreeLess); 

    // Print out the top m trees.
    int sum=0;
    for(vector<CladeTree>::iterator p=cladeTrees.begin();p!=cladeTrees.begin()+m;p++) {
      int ct = p->getCount();
      sum += ct;
      c << setw(5) << ct << " " << setw(6) << setprecision(3) << ct/(double)totalTrees
	<< " " << setw(6) << sum/(double)totalTrees << " ";
      p->getTree()->printWithClade(c);
      c << endl;
    }
    c << endl;
  }

  void printTrans(ostream& c, int maxTopologies) {
    // Print out the named clade transition matrices.
    // findNamedClades must be called first.

    c << "******************** Clade transition matrices ********************" << endl << endl;

    // Set time of last transition for each clade and create transition matrices.
    vector<int> lastTime(namedClades.size()+1,-1); // last time (number of tree in the input) that named clade i was present.
    vector<int> lastTopNum(namedClades.size()+1);  // number of the topology of clade i last time it was present.
    vector< vector< vector<int> > > trans(namedClades.size()+1); // trans[i][j][k] is the number of transitions for clade i from
								 // topology j to topology k, with 0 meaning not present.

    for(vector<CladePtr>::iterator p=namedClades.begin();p!=namedClades.end();p++) {
      int n=(*p)->numTopologies;
      trans[(*p)->name].resize(n+1);
      for(int k=0;k<=n;k++)
	trans[(*p)->name][k].resize(n+1,0);
    }

    // Count all transitions.
    int time=0;
    for(vector<TreePtr>::iterator p=allTrees.begin();p!=allTrees.end();p++,time++)
      addTrans(*p,time,lastTime,lastTopNum,trans);

    // Count last transitions.
    for(vector<CladePtr>::iterator p=namedClades.begin();p!=namedClades.end();p++) {
      int name=(*p)->name;
      if(lastTime[name]!=totalTrees-1)
	trans[name][lastTopNum[name]][0]++;

      // Print out the transition tables.
      int n=(*p)->numTopologies, zeroCol=0, zeroRow=0, sum=0;
      c << "     | ";
      for(int k=1;k<n+1;k++) {
	c << "  " << setw(4) << left << topologyName(name,k) << right;
	sum += trans[name][k][0] + trans[name][0][k];
	for(int j=1;j<n+1;j++)
	  sum += trans[name][k][j];
      }

      trans[name][0][0] = totalTrees - 1 - sum;

      for(int k=0;k<n+1;k++) {
	zeroCol = zeroCol || trans[name][k][0];
	zeroRow = zeroRow || trans[name][0][k];
      }

      if(zeroCol)
	c << "   -  ";
      c << endl << "-----+";
      for(int j=1;j<n;j++)
	c << "------";
      if(zeroCol)
	c << "------";
      c << "-----";
      c << endl;
      for(int k=1;k<n+1;k++) {
	c << " " << setw(4) << left << topologyName(name,k) << right << "|";
	for(int j=1;j<n+1;j++)
	  c << setw(5) << trans[name][k][j] << " ";
	if(zeroCol)
	  c << setw(5) << trans[name][k][0] << " ";
	c << endl;
      }
      if(zeroRow) {
	c << "  -  |";
	for(int j=1;j<n+1;j++)
	  c << setw(5) << trans[name][0][j] << " ";
	if(zeroCol)
	  c << setw(5) << trans[name][0][0] << " ";
	c << endl;
      }
      c << endl << endl;
    }
  }

  void printProbableTreeClades(ostream& c, double threshold) {
    c << "***** Posterior probabilities of clades in most probable tree topology *****" << endl << endl;
    c << "     Count  Prob. Tree topology" << endl;

    TreePtr tree = *max_element(trees[maxLen].begin(),trees[maxLen].end(),treeLt);
    int votesNeeded = (int)(ceil(threshold*totalTrees));
    vector<CladePtr> treeClades;
    collectClades(tree,votesNeeded,treeClades);
    sort(treeClades.begin(),treeClades.end(),cladeLt);
    for(vector<CladePtr>::iterator j=treeClades.begin();j!=treeClades.end();j++) {
      c << setw(10) << (*j)->count << " " << setw(6) << setprecision(3) << (double)(*j)->count/(double)totalTrees << " ";
      (*j)->set.print(c);
      c << endl;
    }
    c << endl;
  }

private:
  class Tree;
  typedef Tree *TreePtr;

  class Clade;
  typedef Clade *CladePtr;

  class Tree {
  public:
    Tree(TreePtr l, TreePtr r, int c, int n=0) : left(l), right(r), count(c), clade(NULL), num(n) {}
    TreePtr left,right;  // left and right subtrees.
    int count;		 // number of times tree appears, for tree > 1 taxa. Taxon number for trees of 1 taxa.
    CladePtr clade;	 // clade corresponding to the tree, if saved. Null otherwise.
    unsigned int num;	 // for a topology of a named clade, the number of that topology. Also used when reading in
    			 // trees to store the trees' hash values, hence the unsigned int.

    void print(ostream& c) {
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

    void printTopology(ostream& c) {
      if(left)
	if (num!=0)
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

    void printWithClade(ostream& c) {
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
  };

  class Clade {
  public:
    Clade() {}
    Clade(Set s, int c) : set(s), count(c), name(1), numTopologies(0) {}
  
    Set set;			 // set of taxa.
    int count;			 // number of times clade appears.
    int name;			 // name of the clade, 0 means unnamed.
    int numTopologies;		 // number of distinct subtree topologies.
    vector<TreePtr> topologies;  // list of subtree topologies.
  };

  int maxTrees;				           // maximum number of trees.
  int numTrees;					   // number of unique trees of two or more taxa.
  int totalTrees;				   // number of trees read in (including duplicates).
  int maxLen;				           // length of a whole tree, i.e., number of taxa.
  vector<TreePtr> allTrees;		           // list of trees read in.
  vector< vector<TreePtr> > trees;                 // trees[i] is the list of unique trees of size i.
  vector< vector<CladePtr> > clades;	           // clades[i] is the list of clades of size i that have sufficient counts.
  vector<CladePtr> namedClades;                    // list of named clades.

  class CladeTree {
  public:

    CladeTree(TreePtr t, int c) : tree(t), count(c) {}

    TreePtr getTree() const { return tree; }

    int getCount() const { return count; }

    void setCount(int c) { count = c; }

    friend bool cladeTreeLt(const CladeTree& a, const CladeTree& b) { return cmpCladeTree(a.tree,b.tree)<0; }
  
    friend bool cladeTreeEq(const CladeTree& a, const CladeTree& b) { return cmpCladeTree(a.tree,b.tree)==0; }

    friend bool cladeTreeLess(const CladeTree& a, const CladeTree& b) { 
      int c=b.count-a.count;
      if(c==0)
	c = cmpTrees(a.tree,b.tree);
      return c<0;
    }

    friend int cmpCladeTree(const TreePtr a, const TreePtr b) {
      if(a==b) 
	return 0;
      if(a->left==NULL)
	if(b->left==NULL)
	  return a->count - b->count;
	else
	  return -1;
      else if(b->left==NULL)
	return 1;
      else if (a->num!=0)
	if(b->num!=0)
	  return (a->clade->name - b->clade->name);
	else 
	  return -1;
      else if(b->num!=0)
	return 1;
      else {
	int c=cmpCladeTree(a->left,b->left);
	return (c !=0 ? c : cmpCladeTree(a->right,b->right));
      }
    }

  private:
    TreePtr tree;
    int count;
  };

  friend bool cladeLt(CladePtr c1, CladePtr c2) { return c1->count > c2->count; }

  friend bool namedCladeLt(CladePtr c1, CladePtr c2) { return c1->set > c2->set; }

  void readFile(istream& f, const string& name, int skip, vector<TreePtr>& sindex, int& first, vector<TreePtr>& taxa) {
    string str;
    for(int i=0;i<skip;i++)
      if(!getline(f,str)) {
	cerr << "Error: " << (name=="-" ? "Standard input" : ("Input file " + name)) 
	     << " contains " << (i==skip-1 ? "exactly" : "fewer than ") << skip << " lines." << endl;
	exit(1);
      }
    while(getline(f,str)) {
      if(first) {
	initializeTrees(str,taxa);
	first = 0;
      }
      string::const_iterator top=str.begin();
      TreePtr tree;
      int len;
      unsigned int hash;
      store(top,str,tree,hash,len,sindex,taxa);
      allTrees.push_back(tree);
    }
  }

  void syntaxError(const string& msg, const string& line, string::const_iterator pos) {
    cerr << "Error: " << msg << " - /" << string(line.begin(),pos) << "** " 
	 << *pos << " **" << string(pos+1,line.end()) << "/" << endl;
    exit(1);
  }

  void initializeTrees(string str, vector<TreePtr>& taxa) {
    trees.resize(2);
    vector<int> taxaNames;
    string::const_iterator top=str.begin();
    getTaxa(top,str,taxaNames);
    maxLen = trees[1].size();
    Set::init(maxLen);
    trees.resize(maxLen+1);
    clades.resize(maxLen+1);
    int maxTaxa = *(max_element(taxaNames.begin(),taxaNames.end()));
    taxa.resize(maxTaxa+1);
    for(int i=0;i<taxaNames.size();i++)
      taxa[taxaNames[i]] = trees[1][i];
  }

  string::const_iterator getTaxa(string::const_iterator top, const string& str, vector<int>& taxaNames) {
    if(*top == '(') {
      top = getTaxa(top+1,str,taxaNames);
      if(*top != ',')
	syntaxError("Missing comma",str,top);
      top = getTaxa(top+1,str,taxaNames);
      if(*top != ')') 
	syntaxError("Missing right parenthesis",str,top);
      return top + 1;
    }
    else {
      int n = 0;
      for(;isdigit(*top);top++)
	n = 10 * n + (*top-'0');
      if(n==0)
	syntaxError("Taxon cannot be zero",str,top);
      taxaNames.push_back(n);
      trees[1].push_back(new Tree(NULL,NULL,n,n));
      return top;
    }
  }

  string::const_iterator store(string::const_iterator top, const string& str, 
			       TreePtr& tree, unsigned int& hash, int& len,
			       vector<TreePtr>& sindex, const vector<TreePtr>& taxa) {
    TreePtr ltree,rtree;
    int lnew,rnew,llen,rlen;
    unsigned int lhash,rhash;

    if(*top == '(') {
      top = store(top+1,str,ltree,lhash,llen,sindex,taxa);
      if(*top != ',')
	syntaxError("Missing comma",str,top);
      top = store(top+1,str,rtree,rhash,rlen,sindex,taxa);
      if(*top != ')') 
	syntaxError("Missing right parenthesis",str,top);
      len = llen + rlen;
      hash = hashfn(lhash,rhash,sindex.size());
      tree = add(hash,ltree,rtree,len,sindex);
      return top + 1;
    }
    else {
      hash = 0;
      for(;isdigit(*top);top++)
	hash = 10 * hash + (*top-'0');
      if(hash==0)
	syntaxError("Taxon cannot be zero",str,top);
      tree = taxa[hash];
      len = 1;
      return top;
    }
  }

  TreePtr add(unsigned int &hash, TreePtr ltree, TreePtr rtree, int len, vector<TreePtr>& sindex) {
    TreePtr tree;
    int hash1 = hash;
    while(tree=sindex[hash1])
      if(tree->left == ltree && tree->right == rtree) {
	tree->count++;
	return tree;
      }
      else
	hash1 = (hash1 + 1) % sindex.size();
    if(numTrees < maxTrees) {
      tree = sindex[hash1] = new Tree(ltree,rtree,1);
      trees[len].push_back(tree);
      numTrees++;
      return tree;
    }
    else { // Re-hash all the trees.
      maxTrees *=2;
      sindex.clear();
      sindex.resize(FindPrime::nextPrime(maxTrees*10/6),NULL);

      for(int i=2;i<=maxLen;i++)
	for(vector<TreePtr>::iterator p=trees[i].begin();p!=trees[i].end();p++) {
	  hash1 = (*p)->num = hashfn((*p)->left->num,(*p)->right->num,sindex.size());
	  while(sindex[hash1])
	    hash1 = (hash1 + 1) % sindex.size();
	  sindex[hash1] = *p;
	}

      hash = hashfn(ltree->num,rtree->num,sindex.size());
      return add(hash,ltree,rtree,len,sindex);
    }
  }

  unsigned int hashfn(unsigned int a, unsigned int b, size_t s) const { return (((a << 6) + b) % s); }
      
  friend bool treeLt(const TreePtr& a, const TreePtr& b)  { 
    int c=b->count-a->count;
    if(c==0)
      c = cmpTrees(a,b);
    return c<0;
  }

  friend int cmpTrees(const TreePtr& a, const TreePtr& b) {
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

  Set getClade(const TreePtr& p) {
    Set lclade = (p->left->clade ? p->left->clade->set : getClade(p->left));
    Set rclade = (p->right->clade ? p->right->clade->set : getClade(p->right));
    Set clade;
    clade.myunion(lclade,rclade);
    return clade;
  }

  friend string topologyName(int name, int num) {
    string buf;
    OSTRSTREAM f(buf);
    f << num;
    return cladeName(name) + f.str();
  }

  friend string cladeName(int name) { 
    string str="";
    return cladeNameInt(name,str);
  }

  friend string cladeNameInt(int n, string& str) { 
    if(n>26)
      cladeNameInt(n/26,str);
    str.push_back(char('A'+n-1));
    return str;
  }

  void addTrans(TreePtr tree, int time, vector<int>& lastTime, vector<int>& lastTopNum, vector< vector< vector<int> > >& trans) {
    int name;
    if(tree->left) {
      if(tree->num!=0) {
	int name=tree->clade->name;
	if(lastTime[name] == -1) {
	  if(time>0)
	    trans[name][0][tree->num]++;
	}
	else if(lastTime[name] == time-1) {
	  trans[name][lastTopNum[name]][tree->num]++;
	}
	else {
	  trans[name][lastTopNum[name]][0]++;
	  trans[name][0][tree->num]++;
	}
	lastTime[name] = time;
	lastTopNum[name] = tree->num;
      }
      addTrans(tree->left,time,lastTime,lastTopNum,trans);
      addTrans(tree->right,time,lastTime,lastTopNum,trans);
    }
  }

  void collectClades(TreePtr tree, int votesNeeded, vector<CladePtr>& treeClades) {
    if(tree->left) {
      CladePtr clade=tree->clade;
      if(clade) // && clade->count>=votesNeeded)
	treeClades.push_back(clade);
      collectClades(tree->left,votesNeeded,treeClades);
      collectClades(tree->right,votesNeeded,treeClades);
    }
  }
};

const char* names[]= {"skip", "trees", "cthreshold", "nthreshold", "maxtopologies"};
const char* desc[] = {"skipped_lines", "number_of_trees_to_print", "threshold_for_clades", "threshold_for_named_clades",
                      "max_tree_topologies" };
const char* defaults[] = {"0", "100", ".01", ".80", "100"};
const int skipField=0, numTreesField=1, cthresholdField=2, nthresholdField=3, maxTopsField=4;

void usageError(ostream& c, char *name)
{
  c << "Usage: " << name;
  int n=sizeof(names)/sizeof(names[0]);
  c << " [-h] [--help]" << endl;
  for(int i=0;i<n;i++)
    if(defaults[i]==string(""))
      c << " {-" << names[i][0] << "|--" <<  names[i] << "} " << desc[i] << endl;
    else
      c << " [-" << names[i][0] << " " << desc[i] << "] [--" <<  names[i] << " " << desc[i] << "]" << endl;
  c << " <file1> <file2> ..." << endl;
  exit(1);
}

void summarize(int skip, int numTrees, double cthreshold, double nthreshold, int maxTops,
	       const vector<string>& filenames, int argc, char *argv[])
{
  Trees trees;

  cout << "******************** BADGER Summarize Version " << VERSION << " ********************" << endl << endl;
  cout << "Invocation: ";
  for(int i=0;i<argc;i++)
    cout << argv[i] << " ";
  cout << endl << endl;

  cout.setf(ios::fixed, ios::floatfield);
  cout.setf(ios::showpoint);

  trees.readFiles(filenames,skip);
  trees.findClades(min(cthreshold,nthreshold));
  trees.findNamedClades(nthreshold,maxTops);

  trees.printNamedClades(cout);
  trees.printTreeTopologies(cout,numTrees);
  trees.printProbableTreeClades(cout,nthreshold);
  trees.printTrans(cout,maxTops);
  trees.printCladeTreeTopologies(cout,numTrees);
  trees.printClades(cout,cthreshold);
  cout << endl;
}

int main(int argc, char *argv[])
{
  int n=sizeof(names)/sizeof(names[0]);
  int okay=1;
  vector<string> val(n,"");
  vector<string> filenames;

  for(int i=1;i<argc;i++) {
    int found = 0;
    for(int j=0;j<n;j++)
      if(argv[i]==string("-h") || argv[i]==string("--help")) {
	usageError(cout,argv[0]);
	exit(0);
      }
    else if(string(argv[i]) == string("-") + string(names[j])[0] || string(argv[i]) == string("--") + string(names[j])) {
	found = 1;
	if(i+1<argc) {
	    val[j] = string(argv[i+1]);
	    i++;
	}
	else
	  okay = 0;
	break;
      }
    if(!found)
      filenames.push_back(argv[i]);
  }

  for(int j=0;j<n;j++) {
    if(val[j]=="")
      val[j] = defaults[j];
    if(val[j]=="")
      okay = 0;
  }

  if(!okay)
    usageError(cerr,argv[0]);


  int skip,numTrees,maxTops;
  double cthreshold,nthreshold;

  ISTRSTREAM f1(val[skipField]);
  if(!(f1 >> skip))
    usageError(cerr,argv[0]);

  ISTRSTREAM f2(val[numTreesField]);
  if(!(f2 >> numTrees))
    usageError(cerr,argv[0]);

  ISTRSTREAM f3(val[cthresholdField]);
  if(!(f3 >> cthreshold))
    usageError(cerr,argv[0]);
  if(cthreshold < 0.0 || cthreshold > 1.0) {
    cerr << "Error: named clade threshold value (" << setprecision(3) 
	 << cthreshold << ") must be between 0.0 and 1.0." << endl;
    exit(1);
  }

  ISTRSTREAM f4(val[nthresholdField]);
  if(!(f4 >> nthreshold))
    usageError(cerr,argv[0]);
  if(nthreshold < 0.5 || nthreshold > 1.0) {
    cerr << "Error: named clade threshold value (" << setprecision(3) 
	 << nthreshold << ") must be between 0.5 and 1.0." << endl;
    exit(1);
  }

  ISTRSTREAM f5(val[maxTopsField]);
  if(!(f5 >> maxTops))
    usageError(cerr,argv[0]);
  if(maxTops < 1) {
    cerr << "Error: Number of topologies (" << maxTops 
	 << ") of a named clade must be >= 1" << "." << endl;
    exit(1);
  }

  if(filenames.size()==0)
    filenames.push_back("-");
  summarize(skip,numTrees,cthreshold,nthreshold,maxTops,filenames,argc,argv);
  return 0;
}

