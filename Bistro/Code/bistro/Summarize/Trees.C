using namespace std;

#include "Trees.h"

void Trees::readFiles(const vector<string>& names, int skip) {
  //cerr << "sizeof(Tree) = " << sizeof(Tree) << endl;
  // Read in tree topologies from the files listed in names, skipping the first skip lines of each.

  // Create hash table used to find unique trees (and subtrees).
  int size = Prime::nextPrime(maxTrees*10/6);
  vector<TreePtr> sindex(size,(TreePtr)NULL); 

  vector<TreePtr> taxa; // Maps taxon numbers to the corresponding trees.

  bool first = true;
  bool hasSemicolon = false;

  hasLength = false;
  blSquared = 0.0;
  Found = NotFound = CompFound = CompNotFound = 0;

  for(vector<string>::const_iterator name=names.begin();name!=names.end();name++) {
    if(*name=="-")
      readFile(cin,*name,skip,sindex,first,taxa,hasLength,hasSemicolon,blSquared);
    else {
      ifstream f((*name).c_str());
      if(!f) {
	cerr << "Error: Could not open " << *name << endl;
	exit(1);
      }
      readFile(f,*name,skip,sindex,first,taxa,hasLength,hasSemicolon,blSquared);
      f.close();
    }
  }
  totalTrees = allTrees.size();
}

void Trees::findClades(double threshold) {
  // Find all the sets that occur in at least threshold (in (0,1]) fraction of the trees.
  // threshold should be the minimum of the threshold for clades to print and the threshold 
  // for named clades.

  // Collect the clades with sufficent votes.
  int totalClades = 0;
  int votesNeeded=(int)(ceil(threshold*totalTrees));
  for(int i=2;i<=maxLen;i++) {
    //cerr << "Number of clades of size " << i << " = " << cladeTable[i].size()-1 << endl;
    totalClades += cladeTable[i].size()-1;
    clades[i].clear();
    for(vector<CladePtr>::iterator p=cladeTable[i].getElements().begin()+1; 
	p!=cladeTable[i].getElements().end(); 
	p++)
      if((*p)->count >= votesNeeded)
	clades[i].push_back(*p);
  }
  //cerr << "Total number of clades = " << totalClades << endl;
}

bool cladeLt(CladePtr c1, CladePtr c2) { return c1->getCount() > c2->getCount(); }

bool namedCladeLt(CladePtr c1, CladePtr c2) { return c1->getSet() > c2->getSet(); }

bool treeLt(const TreePtr& a, const TreePtr& b)  { 
  int c=b->getCount() - a->getCount();
  if(c==0)
    c = cmpTrees(a,b);
  return c<0;
}

void Trees::findNamedClades(double threshold, int maxTopologies) {
  // Find all named clades, i.e., clades that are maximal, occur in at least threshold fraction of trees,
  // and have <= maxTopologies topologies.
  // findClades must be called first.

  vector< vector<CladePtr> > clades2(maxLen+1);

  // Eliminate clades with too few votes or too many topologies.
  int votesNeeded=(int)(ceil(threshold*totalTrees));
  for(int i=2;i<=maxLen;i++) {
    clades2[i].clear();
    for(vector<CladePtr>::iterator p=cladeTable[i].getElements().begin()+1; 
	p!=cladeTable[i].getElements().end(); 
	p++)
      if((*p)->count >= votesNeeded && (*p)->numTopologies <= maxTopologies) {
	(*p)->name = 1;
	clades2[i].push_back(*p);
      }
  }
  // Eliminate non-maximal clades. Collect named clades.
  for(int i=maxLen;i>=2;i--)
    if(clades2[i].size() > 0) {
      for(vector<CladePtr>::iterator p=clades2[i].begin();p!=clades2[i].end();p++) {
	if((*p)->name!=0) {
	  namedClades.push_back(*p);
	  for(int j=2;j<i;j++)
	    for(vector<CladePtr>::iterator q=clades2[j].begin();q!=clades2[j].end();q++)
	      if((*q)->name!=0 && (*q)->set.subset((*p)->set))
		(*q)->name = 0;
	}
      }
    }
    
  // Reset tree nums to 0 and collect topologies for each clade.
  for(int i=2;i<=maxLen;i++) {
    for(vector<TreePtr>::iterator p=trees[i].begin();p!=trees[i].end();p++) {
      (*p)->num = 0;
      if((*p)->clade && (*p)->clade->name!=0)
	(*p)->clade->topologies.push_back(*p);
    }
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

void Trees::printMeanTree(ostream &c) { //not only prints the mean, but computes it
  findSplits();

  for(int i=0;i<maxLen;i++) // going over all leaves (external edges)
    trees[1][i]->clade->value = sqr(trees[1][i]->clade->sumBL/(double)totalTrees);

      
  for(vector<CladePtr>::iterator j=clades[2].begin();j!=clades[2].end();j++)
    (*j)->value = sqr((*j)->sumBL/(double)totalTrees) + splitTable.first((*j)->firstSplit)->value 
      + splitTable.second((*j)->firstSplit)->value;

  for(int i=3;i<=maxLen;i++) {
    for(vector<CladePtr>::iterator j=clades[i].begin();j!=clades[i].end();j++) {
      double val=0;
      for(int k=(*j)->firstSplit;k!=0;k=splitTable.nextSplit(k)) { //finding the max value here
	double val2 = splitTable.first(k)->value + splitTable.second(k)->value; //splitTable must keep track of the best first and second of the subtree in the clade
	if (val2>val)
	  val = val2;
      }
      (*j)->value = sqr((*j)->sumBL/(double)totalTrees) + val;
    }
  }

  c << "********************* Mean Tree *********************" << endl << endl;

  printTotalDist(c,trees[maxLen][0]->clade);
  printMeanTree2(c,trees[maxLen][0]->clade,1);
  c << endl << endl;
  //prettyPrintMeanTree(c,trees[maxLen][0]->clade,1,0,0);
  //c << endl << endl;
}

void Trees::printNamedClades(ostream& c) {
  c << "******************** Named clades ********************" << endl << endl;
  Set others=clades[maxLen][0]->set;
  for(vector<CladePtr>::iterator p=namedClades.begin();p!=namedClades.end();p++) {
    others.difference(others,(*p)->set);
    c << setw(6) << (*p)->count << "  " << setw(4) << left << cladeName((*p)->name) << right << " ";
    (*p)->set.print(c);
    c << endl;
    for(vector<TreePtr>::iterator q=(*p)->topologies.begin();q!=(*p)->topologies.end();q++) {
      c << setw(8) << (*q)->count << "  " << setw(4) << left << topologyName((*p)->name,1) 
	<< right << " ";
      (*q)->print(c);
      c << endl;
    }
    c << endl;
  }
  others.printElements(c);
  c << endl;
}

void Trees::printClades(ostream& c, double threshold) {
  c << "************************ Common Clades ************************" << endl << endl;
  int votesNeeded = (int)(ceil(threshold*totalTrees));
  vector<CladePtr> pclades;
  for(int i=2;i<=maxLen;i++)
    for(vector<CladePtr>::iterator j=clades[i].begin();j!=clades[i].end();j++) 
      if((*j)->count>=votesNeeded)
	pclades.push_back(*j);
  sort(pclades.begin(),pclades.end(),namedCladeLt);
  int i = 0;
  for(vector<CladePtr>::iterator j=pclades.begin();j!=pclades.end();j++, i++) {
    c << setw(10) << (*j)->count << " " << setw(6) << setprecision(3)
      << (double)(*j)->count/(double)totalTrees << " ";
    (*j)->set.print(c);
    c << endl;
  }
  c << endl;
}

void Trees::printTreeTopologies(ostream& c, int m) {
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


int cmpCladeTree(const TreePtr a, const TreePtr b) {
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

bool cladeTreeLt(const CladeTree& a, const CladeTree& b) { return cmpCladeTree(a.getTree(),b.getTree())<0; }
  
bool cladeTreeEq(const CladeTree& a, const CladeTree& b) { return cmpCladeTree(a.getTree(),b.getTree())==0; }

bool cladeTreeLess(const CladeTree& a, const CladeTree& b) { 
  int c=b.getCount()-a.getCount();
  if(c==0)
    c = cmpTrees(a.getTree(),b.getTree());
  return c<0;
}

void Trees::printCladeTreeTopologies(ostream& c, int m) {
  // Print out the top m (or all if m=0) tree topologies by count in order from largest to smallest,
  // with a named clade substituted for the corresponding named clade in the tree.
  // findNamedClades must be called before this.

  c << "******************** Clade tree topologies " << "********************" << endl << endl;
  c << "Count  Prob.  Cum.  Tree topology" << endl;

  // Find the counts of all trees with the same topology when named clades are considered 
  // the same despite topology.
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

void Trees::printTrans(ostream& c, int maxTopologies) {
  // Print out the named clade transition matrices.
  // findNamedClades must be called first.

  c << "******************** Clade transition matrices ********************" << endl << endl;

  // Set time of last transition for each clade and create transition matrices.
  vector<int> lastTime(namedClades.size()+1,-1); // last time (number of tree in the input) that 
  // named clade i was present.
  vector<int> lastTopNum(namedClades.size()+1);  // number of the topology of clade i last time 
  // it was present.
  vector< vector< vector<int> > > trans(namedClades.size()+1); // trans[i][j][k] is the number of 
  // transitions for clade i from
  // topology j to topology k, 
  // with 0 meaning not present.

  for(vector<CladePtr>::iterator p=namedClades.begin();p!=namedClades.end();p++) {
    int n=(*p)->topologies.size();
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
    int n=(*p)->topologies.size(), zeroCol=0, zeroRow=0, sum=0;
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

void Trees::printProbableTreeClades(ostream& c, double threshold) {
  c << "***** Posterior probabilities of clades in most probable tree topology *****" << endl << endl;
  c << "     Count  Prob. Tree topology" << endl;

  TreePtr tree = *max_element(trees[maxLen].begin(),trees[maxLen].end(),treeLt);
  int votesNeeded = (int)(ceil(threshold*totalTrees));
  vector<CladePtr> treeClades;
  collectClades(tree,votesNeeded,treeClades);
  sort(treeClades.begin(),treeClades.end(),cladeLt);
  for(vector<CladePtr>::iterator j=treeClades.begin();j!=treeClades.end();j++) {
    c << setw(10) << (*j)->count << " " << setw(6) << setprecision(3) 
      << (double)(*j)->count/(double)totalTrees << " ";
    (*j)->set.print(c);
    c << endl;
  }
  c << endl;
}

void Trees::printDistanceMatrices(ostream &c) {
  if(hasLength) {
    c << "*********************** Distance Mean ";
    c << "***********************" << endl << endl;
    for(int i1=0;i1<maxLen;i1++) {
      for(int i2=0;i2<maxLen;i2++) 
	c << setprecision(6) << setw(12) << leafDist[i1][i2]/totalTrees << " ";
      c << endl;
    }
    c << endl;
      
    c << "************************ Variance ";
    c << "************************" << endl << endl;
    for(int i1=0;i1<maxLen;i1++) {
      for(int i2=0;i2<maxLen;i2++) 
	c << setprecision(6) << setw(12) 
	  << leafDistSqr[i1][i2]/totalTrees - leafDist[i1][i2]/sqr(totalTrees) << " ";
      c << endl;
    }
  }
}

void Trees::findSplits() {
  vector<int> offsets(maxLen+2);
  offsets[0] = 0;
  for(int i=0;i<=maxLen;i++)
    offsets[i+1] = offsets[i] + cladeTable[i].size() - 1;
  for(int i=2;i<=maxLen;i++) {
    //cerr << i << ":" << endl;
    splitTable.resetCounts();
    int count = 0;
    for(vector<TreePtr>::iterator p=trees[i].begin();p!=trees[i].end();p++) {
      CladePtr left=(*p)->left->clade, right=(*p)->right->clade;
      int lnum = left->num, llen = left->len, rnum = right->num, rlen = right->len;
      // We use the num test as we don't have the least element on hand.
      if(llen > rlen || (llen == rlen && lnum < rnum))
	addEdge((*p)->clade, left, right, (*p)->count, offsets);
      else
	addEdge((*p)->clade, right, left, (*p)->count, offsets);

      //if(++count % 10000 == 0)
      //splitTable.printCounts(count);
    }
    //cerr << endl;
  }
}

void Trees::addEdge(CladePtr parent, CladePtr first, CladePtr second, int count, vector<int> offsets) {
  int t; 	    // entry in the splitTable.
  unsigned int k; // entry in the splitTable's index.

  if(splitTable.find(first, second, t, k, offsets))
    splitTable.addCount(t,count);
  else 
    splitTable.add(k, first, second, count, parent, offsets);
}

void Trees::readFile(istream& f, const string& name, int skip, vector<TreePtr>& sindex, bool& first, 
		     vector<TreePtr>& taxa, bool& hasLength, bool& hasSemicolon, double& blSquared) {
  string str;
  for(int i=0;i<skip;i++)
    if(!getline(f,str)) {
      cerr << "Error: " << (name=="-" ? "Standard input" : ("Input file " + name)) 
	   << " contains " << (i==skip-1 ? "exactly" : "fewer than ") << skip << " lines." << endl;
      exit(1);
    }
  if(getline(f,str)) {
    if(first) {
      initializeTrees(str,taxa,hasLength,hasSemicolon);
      first = false;
    }

    int count = 0;
    TreePtr tree;
    if(hasLength)
      do {
	string::const_iterator top=str.begin();
	storeTopWithLen(top,str,tree,sindex,taxa,hasSemicolon,blSquared);
	allTrees.push_back(tree);
      }
      while(getline(f,str));
    else 
      do {
	//if(++count % 1000 == 0)
	//cerr << ".";
	if(++count % 10000 == 0)
	  cerr << setw(10) << count << " (" << setw(15) << numTrees << "): ave. # probes (found) = " 
	       << setw(8) << setprecision(6) << CompFound/Found << ", Ave. # probes (not found) = " 
	       << setw(8) << setprecision(6) << CompNotFound/NotFound << endl;
	string::const_iterator top=str.begin();
	storeTop(top,str,tree,sindex,taxa,hasSemicolon,blSquared);
	allTrees.push_back(tree);
      }
      while(getline(f,str));
    cerr << endl;
    //cerr << "Done reading in file." << endl;
  }
}

void Trees::syntaxError(const string& msg, const string& line, string::const_iterator pos) {
  cerr << "Error: " << msg << " - /" << string(line.begin(),pos) << "** " 
       << *pos << " **" << string(pos+1,line.end()) << "/" << endl;
  exit(1);
}

void Trees::initializeTrees(string str, vector<TreePtr>& taxa, bool& hasLength, bool& hasSemicolon) {
  trees.resize(2);
  vector<int> taxaNames;
  string::const_iterator top=str.begin();
  hasLength = false;
  getTaxaTop(top,str,taxaNames,hasLength,hasSemicolon);
  maxLen = trees[1].size();
  trees.resize(maxLen+1);
  clades.resize(maxLen+1);
  int maxTaxa = *(max_element(taxaNames.begin(),taxaNames.end()));
  taxa.resize(maxTaxa+1);
  cladeSet.init(maxLen); 
  cladeTable.resize(maxLen+1);
  leafDist.resize(maxLen);
  leafDistSqr.resize(maxLen);
  for(int i=0;i<maxLen;i++) {
    leafDist[i].resize(maxLen);
    leafDistSqr[i].resize(maxLen);
  }

  for(int i=0;i<taxaNames.size();i++) {
    taxa[taxaNames[i]] = trees[1][i];
    // Create singleton clades for the leaves.
    cladeSet.singleton(trees[1][i]->count,maxLen);
    trees[1][i]->clade = new Clade(cladeSet,1,1,0,1,i);
    cladeTable[1].push_back(trees[1][i]->clade);
    clades[1].push_back(trees[1][i]->clade);
  }
}

void Trees::skipLength(string::const_iterator& top, const string& str, bool& hasLength) {
  if(top != str.end() && *top == ':') {
    top++;
    while (top!=str.end() && *top != ',' && *top != ')')
      top++;
    hasLength = true;
  }
}

string::const_iterator Trees::getTaxaTop(string::const_iterator top, const string& str, vector<int>& taxaNames,
					 bool& hasLength, bool& hasSemicolon) {
  if(*top == '(') {
    top = getTaxa(top+1,str,taxaNames,hasLength);
    if(*top != ',')
      syntaxError("Missing comma",str,top);
    top = getTaxa(top+1,str,taxaNames,hasLength);
    if(*top != ')') 
      syntaxError("Missing right parenthesis",str,top);
    top++;
    if(top != str.end() && *top == ';') {
      hasSemicolon = true;
      top++;
    }
    return top;
  }
  else {
    int n = 0;
    for(;isdigit(*top);top++)
      n = 10 * n + (*top-'0');
    if(n==0)
      syntaxError("Taxon cannot be zero",str,top);
    if(top != str.end() && *top == ';') {
      hasSemicolon = true;
      top++;
    }
    taxaNames.push_back(n);
    trees[1].push_back(new Tree(NULL,NULL,n,n));
    return top;
  }
}

string::const_iterator Trees::getTaxa(string::const_iterator top, const string& str, vector<int>& taxaNames,
				      bool& hasLength) {
  if(*top == '(') {
    top = getTaxa(top+1,str,taxaNames,hasLength);
    if(*top != ',')
      syntaxError("Missing comma",str,top);
    top = getTaxa(top+1,str,taxaNames,hasLength);
    if(*top != ')') 
      syntaxError("Missing right parenthesis",str,top);
    top++;
    skipLength(top,str,hasLength);

    return top;
  }
  else {
    int n = 0;
    for(;isdigit(*top);top++)
      n = 10 * n + (*top-'0');
    if(n==0)
      syntaxError("Taxon cannot be zero",str,top);
    skipLength(top,str,hasLength);
    taxaNames.push_back(n);
    trees[1].push_back(new Tree(NULL,NULL,n,n));
    return top;
  }
}

bool Trees::isNumeric(char ch) {
  switch(ch) {
  case '0': case '1': case '2': case '3': case '4': case '5': case '6': case '7': case '8': case '9': 
  case '.': case '+': case '-': case 'E': case 'e': return true;
  default: return false;
  }
}

void Trees::readLength(string::const_iterator& top, const string& str, double& bl) {
  if(top == str.end() || *top != ':')
    syntaxError("Missing colon",str,top);
  top++;
  string::const_iterator top0 = top;
  while(top != str.end() && isNumeric(*top))
    top++;
  if(top == top0)
    syntaxError("Missing branch length",str,top0);
  string number = string(top0,top);
  ISTRSTREAM f1(number.c_str());
  if(!((f1 >> bl)))
    syntaxError("Ill-formed branch length",str,top0);
}

string::const_iterator Trees::storeTopWithLen(string::const_iterator top, const string& str, 
					      TreePtr& tree, vector<TreePtr>& sindex, 
					      const vector<TreePtr>& taxa, bool hasSemicolon, 
					      double& blSquared) {
  int len;
  unsigned int hash;
  TreePtr ltree,rtree;
  int lnew,rnew,llen,rlen;
  unsigned int lhash,rhash;
  double bl1,bl2;
  vector<double> dist,rdist;
  vector<int> leaves,rleaves;

  if(*top == '(') {
    top = storeWithLen(top+1,str,ltree,lhash,sindex,taxa,bl1,blSquared,dist,leaves);
    if(*top != ',')
      syntaxError("Missing comma",str,top);
    top = storeWithLen(top+1,str,rtree,rhash,sindex,taxa,bl2,blSquared,rdist,rleaves);
    if(*top != ')') 
      syntaxError("Missing right parenthesis",str,top);
    len = llen + rlen;
    blSquared += sqr(bl1+bl2);
    hash = hashfn(lhash,rhash,sindex.size());
    tree = add(hash,ltree,rtree,sindex,leaves.size(),rleaves.size(),0);
    top++;
    if(hasSemicolon)
      if(top == str.end() || *top != ';')
	syntaxError("Missing semicolon",str,top);

    vector<int>::iterator i1, i2;
    vector<double>::iterator d1, d2;
    for(i1=leaves.begin(),d1=dist.begin();i1!=leaves.end();i1++,d1++)
      for(i2=rleaves.begin(),d2=rdist.begin();i2!=rleaves.end();i2++,d2++)
	if(*i1 < *i2) {
	  leafDist[*i1][*i2] += *d1 + *d2;
	  leafDistSqr[*i1][*i2] += sqr(*d1 + *d2);
	}
	else {
	  leafDist[*i2][*i1] += *d1 + *d2;
	  leafDistSqr[*i2][*i1] += sqr(*d1 + *d2);
	}

    return top;
  }
  else {
    hash = 0;
    for(;isdigit(*top);top++)
      hash = 10 * hash + (*top-'0');
    if(hash==0)
      syntaxError("Taxon cannot be zero",str,top);
    if(hasSemicolon)
      if(top == str.end() || *top != ';')
	syntaxError("Missing semicolon",str,top);
    tree = taxa[hash];
    len = 1;
    return top;
  }
}

string::const_iterator Trees::storeWithLen(string::const_iterator top, const string& str, 
					   TreePtr& tree, unsigned int& hash, 
					   vector<TreePtr>& sindex, const vector<TreePtr>& taxa,
					   double &bl, double& blSquared,
					   vector<double> &dist, vector<int> &leaves) {
  TreePtr ltree,rtree;
  unsigned int lhash,rhash;
  double bl1,bl2;
  vector<double> rdist;
  vector<int> rleaves;

  if(*top == '(') {
    top = storeWithLen(top+1,str,ltree,lhash,sindex,taxa,bl1,blSquared,dist,leaves);
    if(*top != ',')
      syntaxError("Missing comma",str,top);
    top = storeWithLen(top+1,str,rtree,rhash,sindex,taxa,bl2,blSquared,rdist,rleaves);
    if(*top != ')') 
      syntaxError("Missing right parenthesis",str,top);
    top++;
    readLength(top,str,bl);
    blSquared += sqr(bl1) + sqr(bl2);
    hash = hashfn(lhash,rhash,sindex.size());
    tree = add(hash,ltree,rtree,sindex,leaves.size(),rleaves.size(),bl);

    vector<int>::iterator i1, i2;
    vector<double>::iterator d1, d2;
    for(i1=leaves.begin(),d1=dist.begin();i1!=leaves.end();i1++,d1++)
      for(i2=rleaves.begin(),d2=rdist.begin();i2!=rleaves.end();i2++,d2++)
	if(*i1 < *i2) {
	  leafDist[*i1][*i2] += *d1 + *d2;
	  leafDistSqr[*i1][*i2] += sqr(*d1 + *d2);
	}
	else {
	  leafDist[*i2][*i1] += *d1 + *d2;
	  leafDistSqr[*i2][*i1] += sqr(*d1 + *d2);
	}
    leaves.insert(leaves.end(),rleaves.begin(),rleaves.end());
    dist.insert(dist.end(),rdist.begin(),rdist.end());
    for(vector<double>::iterator p=dist.begin();p!=dist.end();p++)
      (*p) += bl;
    return top;
  }
  else {
    hash = 0;
    for(;isdigit(*top);top++)
      hash = 10 * hash + (*top-'0');
    if(hash==0)
      syntaxError("Taxon cannot be zero",str,top);
    readLength(top,str,bl);
    tree = taxa[hash];
    tree->clade->sumBL += bl;
    leaves.push_back(hash - 1);
    dist.push_back(bl);
    return top;
  }
}

string::const_iterator Trees::storeTop(string::const_iterator top, const string& str, 
				       TreePtr& tree, vector<TreePtr>& sindex, 
				       const vector<TreePtr>& taxa, bool hasSemicolon, double& blSquared) {
  int len;
  unsigned int hash;
  TreePtr ltree,rtree;
  int lnew,rnew,llen,rlen;
  unsigned int lhash,rhash;

  blSquared += 2*maxLen - 2;

  if(*top == '(') {
    top = store(top+1,str,ltree,lhash,llen,sindex,taxa);
    if(*top != ',')
      syntaxError("Missing comma",str,top);
    top = store(top+1,str,rtree,rhash,rlen,sindex,taxa);
    if(*top != ')') 
      syntaxError("Missing right parenthesis",str,top);
    len = llen + rlen;
    hash = hashfn(lhash,rhash,sindex.size());
    tree = add(hash,ltree,rtree,sindex,llen,rlen,0);
    top++;
    if(hasSemicolon)
      if(top == str.end() || *top != ';')
	syntaxError("Missing semicolon",str,top);
    return top;
  }
  else {
    hash = 0;
    for(;isdigit(*top);top++)
      hash = 10 * hash + (*top-'0');
    if(hash==0)
      syntaxError("Taxon cannot be zero",str,top);
    if(hasSemicolon)
      if(top == str.end() || *top != ';')
	syntaxError("Missing semicolon",str,top);
    tree = taxa[hash];
    len = 1;
    return top;
  }
}

string::const_iterator Trees::store(string::const_iterator top, const string& str, 
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
    tree = add(hash,ltree,rtree,sindex,llen,rlen,1);
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
    tree->clade->sumBL += 1;
    return top;
  }
}

TreePtr Trees::add(unsigned int &hash, TreePtr ltree, TreePtr rtree, 
		   vector<TreePtr>& sindex, int llen, int rlen, double branchLength) {
  TreePtr tree;
  int hash1 = hash;
  int comps = 0;
  int len = llen + rlen;

  while(tree=sindex[hash1])
    if((comps++, tree->left == ltree && tree->right == rtree)) {
      tree->count++;
      tree->clade->count++;
      tree->clade->sumBL += branchLength;
      Found++;
      CompFound += comps;
      return tree;
    }
    else
      hash1 = (hash1 + 1) % sindex.size();

  if(numTrees < maxTrees) {
    tree = sindex[hash1] = new Tree(ltree,rtree,1);
    trees[len].push_back(tree);
    numTrees++;

    NotFound++;
    CompNotFound += comps;

    cladeSet.myunion(ltree->clade->set,rtree->clade->set);

    unsigned int k; // position in cladeTable's index
    if(cladeTable[len].find(cladeSet,k,tree->clade)) {
      tree->clade->count++;
      tree->clade->sumBL += branchLength;
      tree->clade->numTopologies++;
    }
    else {
      tree->clade = new Clade(cladeSet,1,1,branchLength,len,0); // num will be set by add method.
      cladeTable[len].add(k, tree->clade);
    }
    return tree;
  }
  else { // Re-hash all the trees.
    maxTrees *=2;
    sindex.clear();
    sindex.resize(Prime::nextPrime(maxTrees*10/6),NULL);
    cerr << "Re-hashing trees, new size = " << sindex.size() << ", space = " 
	 << sizeof(TreePtr)*sindex.size() << endl;

    for(int i=2;i<=maxLen;i++) {
      for(vector<TreePtr>::iterator p=trees[i].begin();p!=trees[i].end();p++) {
	hash1 = (*p)->num = hashfn((*p)->left->num,(*p)->right->num,sindex.size());
	while(sindex[hash1])
	  hash1 = (hash1 + 1) % sindex.size();
	sindex[hash1] = *p;
      }
    }
    //cerr << "Re-hashed trees" << endl;

    hash = hashfn(ltree->num,rtree->num,sindex.size());
    return add(hash,ltree,rtree,sindex,llen,rlen,branchLength);
  }
}

void Trees::printTotalDist(ostream& c, CladePtr j) {
  double val=0;

  CladePtr first,second;
  if(j->len == 2) {
    first = splitTable.first(j->firstSplit);
    second = splitTable.second(j->firstSplit);
  }
  else {
    for(int k=j->firstSplit;k!=0;k=splitTable.nextSplit(k)) {
      double val2 = splitTable.first(k)->value + splitTable.second(k)->value;
      if (val2>val) {
	val = val2;
	first = splitTable.first(k);
	second = splitTable.second(k);
      }
    }
  }
  double firstBL,secondBL;
  firstBL = first->sumBL/double(totalTrees);
  secondBL = second->sumBL/double(totalTrees);
  c << "Total Distance = " 
    << blSquared 
    - (totalTrees * j->value - sqr(firstBL) - sqr(secondBL) + sqr(firstBL + secondBL)) << endl;
}


void Trees::printMeanTree2(ostream &c, CladePtr j, int topTree) {
  if(j->len == 1)
    c << trees[1][j->num]->count << ":" << j->sumBL/(double)totalTrees;
  else {
    CladePtr first,second;
    c << "(";
    if(j->len == 2) {
      first = splitTable.first(j->firstSplit);
      second = splitTable.second(j->firstSplit);
      c << trees[1][first->num]->count << ":" << first->sumBL/(double)totalTrees << ","
	<< trees[1][second->num]->count << ":" << second->sumBL/(double)totalTrees;
    }
    else {
      double val=0;
      for(int k=j->firstSplit;k!=0;k=splitTable.nextSplit(k)) { //going over all sub-splits of clade j
	double val2 = splitTable.first(k)->value + splitTable.second(k)->value;
	if (val2>val) { //searching for the maximum
	  val = val2;
	  first = splitTable.first(k);
	  second = splitTable.second(k);
	}
      }
      if(first->set > second->set) {
	printMeanTree2(c,first,0);
	c << ",";
	printMeanTree2(c,second,0);
      }
      else {
	printMeanTree2(c,second,0);
	c << ",";
	printMeanTree2(c,first,0);
      }
    }
    c << ")";
    if(topTree)
      c << ";";
    else {
      c << ":" << j->sumBL/(double)totalTrees;
    }
  }
}
	
void Trees::prettyPrintMeanTree(ostream &c, CladePtr j, int topTree, int indent, int stayOnLine) {
  for(int i=0;i<indent;i++)
    c << " ";
  indent++;
  if(j->len == 1)
    c << trees[1][j->num]->count << ":" << j->sumBL/(double)totalTrees;
  else {
    CladePtr first,second;
    c << "(";
    if(j->len==2) {
      first = splitTable.first(j->firstSplit);
      second = splitTable.second(j->firstSplit);
      c << trees[1][first->num]->count << ":" << first->sumBL/(double)totalTrees << ","
	<< trees[1][second->num]->count << ":" << second->sumBL/(double)totalTrees;
    }
    else {
      double val=0;
      for(int k=j->firstSplit;k!=0;k=splitTable.nextSplit(k)) {
	double val2 = splitTable.first(k)->value + splitTable.second(k)->value;
	if (val2>val) {
	  val = val2;
	  first = splitTable.first(k);
	  second = splitTable.second(k);
	}
      }
      if(first->set > second->set) {
	c << endl;
	prettyPrintMeanTree(c,first,0,indent,1);
	c << "," << endl;
	prettyPrintMeanTree(c,second,0,indent,0);
      }
      else {
	if(second->len == 1)
	  prettyPrintMeanTree(c,second,0,0,1);
	else {
	  c << endl;
	  prettyPrintMeanTree(c,second,0,indent,1);
	}
	c << "," << endl;
	prettyPrintMeanTree(c,first,0,indent,0);
      }
    }
    c << ")";
    if(topTree)
      c << ";";
    else {
      c << ":" << j->sumBL/(double)totalTrees;
      if(!stayOnLine) {
	c << endl;
	for(int i=0;i<indent-2;i++)
	  c << " ";
      }
    }
  }
}
	    
void Trees::addTrans(TreePtr tree, int time, vector<int>& lastTime, vector<int>& lastTopNum, 
		    vector< vector< vector<int> > >& trans) {
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

void Trees::collectClades(TreePtr tree, int votesNeeded, vector<CladePtr>& treeClades) {
  if(tree->left) {
    CladePtr clade=tree->clade;
    if(clade) // && clade->count>=votesNeeded)
      treeClades.push_back(clade);
    collectClades(tree->left,votesNeeded,treeClades);
    collectClades(tree->right,votesNeeded,treeClades);
  }
}
