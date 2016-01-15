#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
//#include <algorithm>
//#include <functional>
//#include <map>
#include "tree.h"

using namespace std;

string VERSION="1.0.0";

// dist is a numTaxa by numTaxa matrix
//   taxon number i (1 to n) is stored in row/column i-1.
//   dist[i][j] will contain the distance from taxon i-1 to taxon j-1.
//   dist[i][j] is set by the mrca of these taxa in the rooted tree.
// depths is a vector of length numTaxa
//   Each node sets it to be the distance from itself to the leaves in its subtree.
//   The node does not touch the distances for nodes not in its subtree as these may
//     have been set elsewhere.
// subtreeTaxa is the set of taxa numbers contained in the subtree rooted at this node.
//
void Node::setDistances(DistanceMatrix& dist,vector<double>& depths,vector<int>& subtreeTaxa,Edge* parent)
{
  subtreeTaxa.clear();
  // if leaf
  if ( leaf ) {
    dist.set(number-1,number-1,0.0);
    depths[number-1] = 0.0;
    subtreeTaxa.push_back(number);
    return;
  }
  // internal nodes:
  // if called on the root, then parent==NULL
  // children[i] is the vector of taxa numbers in the ith subtree of this node
  vector<vector<int> > children;
  for ( vector<Edge*>::iterator p=edges.begin(); p != edges.end(); ++p ) {
    if ( (*p) != parent ) {
      vector<int> taxa;
      taxa.resize(0);
      children.push_back( taxa );
      (*p)->getOtherNode(this)->setDistances(dist,depths,children.back(),*p);
      // add length of edge *p to the depths of all children in the corresponding subtree
      // add child taxon number of subtreeTaxon return vector
      for ( vector<int>::iterator c=children.back().begin(); c != children.back().end(); ++c ) {
	depths[*c-1] += (*p)->getLength();
	subtreeTaxa.push_back(*c);
      }
    }
  }
  
  // set distances between all pairs of taxa in separate subtrees
  for ( vector<vector<int> >::iterator subtree1=children.begin(); subtree1 != children.end() - 1; ++subtree1 )
    for ( vector<vector<int> >::iterator subtree2=subtree1 + 1; subtree2 != children.end(); ++subtree2 )
      for ( vector<int>::iterator i=(*subtree1).begin(); i != (*subtree1).end(); ++i )
	for ( vector<int>::iterator j=(*subtree2).begin(); j != (*subtree2).end(); ++j ) {
	  double d = depths[*i-1] + depths[*j-1];
	  dist.set(*i-1,*j-1,d);
	  dist.set(*j-1,*i-1,d);
	}
}

void Tree::setDistances(DistanceMatrix& dist)
{
  dist.empty();
  vector<double> depths(numTaxa);
  vector<int> taxa;
  Edge* parent=NULL;
  root->setDistances(dist,depths,taxa,parent);
}

void usage()
{
  cerr << "Usage: distance [--help] [--version] [--skip number-of-trees-to-skip] [--print-distance-matrices] --tree treefile --out outfile" << endl;
  exit(1);
}

void checkFlag(string key,int i,int argc) {
  if ( i==argc ) {
    cerr << "Error: flag " << key << " is not followed by required argument." << endl;
    usage();
  }
}

void processCommandLine(int argc, char* argv[], ifstream& tin, ofstream& out, int& numSkip, bool& printDistanceMatrices, char& eol)
{
  if ( argc == 1 )
    usage();

  bool treeFileGiven=false;
  bool outFileGiven=false;
  int i=1;
  while ( i < argc ) {
    string key=argv[i++];

    if ( key=="--help" || key=="=h" )
      usage();

    if ( key=="--version" ) {
      cerr << "distance VERSION " << VERSION << endl;
      exit(2);
    }

    if ( key=="--tree" ) {
      checkFlag(key,i,argc);
      tin.open(argv[i]);
      if ( tin.fail() || tin.bad() ) {
	cerr << "Error: cannot read file " << argv[i] << endl;
	exit(1);
      }
      treeFileGiven = true;
    }
    else if ( key=="--out" ) {
      checkFlag(key,i,argc);
      out.open(argv[i]);
      if ( out.fail() || out.bad() ) {
	cerr << "Error: cannot open output file " << argv[i] << endl;
	exit(1);
      }
      outFileGiven = true;
    }
    else if ( key == "--skip" ) {
      checkFlag(key,i,argc);
      istringstream f((string)(argv[i]));
      f >> numSkip;
      if ( numSkip < 0 ) {
	cerr << "Error: Cannot skip a negative number of trees." << endl;
	usage();
      }
    }
    else if ( key=="print-distance-matrices" ) {
      printDistanceMatrices = true;
    }
    else {
      cerr << "Error: Did not recognize argument flag " << key << endl;
      usage();
    }
    ++i;
  }
  // Error if tree file or outfile are missing
  if ( !treeFileGiven )
    cerr << "Error: no flag --tree followed by tree input file." << endl;
  if ( !outFileGiven )
    cerr << "Error: no flag --out followed by output file name." << endl;
  if ( !treeFileGiven || !outFileGiven )
    usage();

  // determine EOL character
  do {
    eol = tin.get();
  } while ( (eol != '\n') && (eol != '\r') && !tin.fail() );
  if ( eol == '\r' ) {
    eol = tin.peek();
    if ( eol != '\n' )
      eol = '\r';
  }
}

void provisionalMeans(DistanceMatrix& means, DistanceMatrix& sumsq, const DistanceMatrix& dist, int k)
{
  int n = dist.size();
  for ( int i=0; i<n; ++i )
    for ( int j=i; j<n; ++j ) {
      double x = dist.get(i,j) - means.get(i,j); // dist - old mean
      means.set(i,j, means.get(i,j) + x/k);
      means.set(j,i, means.get(i,j));
      sumsq.set(i,j, sumsq.get(i,j) + x*(dist.get(i,j)-means.get(i,j)));
      sumsq.set(j,i, sumsq.get(i,j));
    }
}

int main(int argc, char *argv[])
{

  ifstream tin;
  ofstream out;
  int numSkip=0;
  bool printDistanceMatrices=false;
  char eol;
  processCommandLine(argc,argv,tin,out,numSkip,printDistanceMatrices,eol);
  int numActuallySkipped=0;
  int skipInFile = numSkip;


  // read in file one line at a time
  // if it is a tree, process it
  // if it is not a tree, just skip it
  string line;
  int lineNumber=0;
  int numTrees=0;
  DistanceMatrix means;
  DistanceMatrix sumsq;

  while ( getline(tin,line,eol) ) {
      lineNumber++;
      istringstream s(line);
      string keyTree,name,equalSign;

      s >> keyTree;
      if ( keyTree != "tree" )
	continue;

      if(skipInFile-- > 0) {
	numActuallySkipped++;
	continue;
      }

      s >> name >> equalSign;

      if ( equalSign != "=" )
	continue;

      // Skip spurious characters before the first '('
      while ( !s.fail() ) {
	char c;
	c = s.get();
	if ( c == '(' ) {
	  s.putback(c);
	  break;
	}
      }
      
      // The rest should be a parenthetic representation of a tree.  Assume no spaces!

      string treeString;
      s >> treeString;
      Tree tree(treeString,lineNumber);
      numTrees++;
      DistanceMatrix dist(tree.getNumTaxa());
      tree.setDistances(dist);
      if ( printDistanceMatrices ) {
	out << "Tree " << numTrees << endl;
	dist.print(out);
      }
      if ( means.size() == 0 ) { // first non-skipped tree
	means.reset(tree.getNumTaxa());
	sumsq.reset(tree.getNumTaxa());
      }
      provisionalMeans(means,sumsq,dist,numTrees);
  }
  out << "Mean Distances:" << endl;
  means.print(out);
  out << "Standard Deviations:" << endl;
  sumsq.divide(numTrees); // divide sum of squares by numTrees for each element
  sumsq.sqrt(); // take square root of each element
  sumsq.print(out);

  cout << "Skipped " << numActuallySkipped << " trees." << endl;
  cout << "Read " << numTrees << " trees." << endl;
  return 0;
}
