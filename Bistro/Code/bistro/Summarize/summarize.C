#define VERSION "CVS"

/* Summarize for BADGER 
 * (c) Copyright 1998, 1999, 2004, 2005, 2006, 2007  by Donald Simon & Bret Larget
 * Department of Mathematics/Computer Science, Duquesne University
 * 4/12/2007
 */

/* Changelist: Incorporates mean tree code, average tree code, reading
 * in different trees with branchlengths or no branchlengths, either
 * int or double
 *
 * Still to do: Read in unrooted trees, turn on/off parts of the output,
 * read in NEXUS blocks.
 */

using namespace std;

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <assert.h>
#include <sstream>

#include "Trees.h"

bool debug = false;

const char* names[]= {"skip", "trees", "cthreshold", "nthreshold", "maxtopologies"};
const char* desc[] = {"skipped_lines", "number_of_trees_to_print", "threshold_for_clades", 
		      "threshold_for_named_clades",
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

void test()
{
  cerr << "test starting" << endl;
  Trees trees;
  const char *args[] = {"(1:1,(2:1,(3:1,4:1):1):1);", "((1:1,2:1):1,(3:1,4:1):1);"};
  vector<string> names(args, std::end(args));
  cerr << "succesful creation of names vector" << endl;
  trees.readTrees(names);
  cerr << "successful reading of trees" << endl;
  ostringstream c;
  string meanTree = trees.printMeanTree(c);
  cerr << "successful mean of trees" << endl << meanTree << endl;
}

void summarize(int skip, int numTrees, double cthreshold, double nthreshold, int maxTops,
	       const vector<string>& filenames, int argc, char *argv[])
{
  Trees trees;

  cout << "******************** BADGER Summarize Version " << VERSION << " ********************" 
       << endl << endl;
  cout << "Invocation: ";
  for(int i=0;i<argc;i++)
    cout << argv[i] << " ";
  cout << endl << endl;

  cout.setf(ios::fixed, ios::floatfield);
  cout.setf(ios::showpoint);

  trees.readFiles(filenames,skip);
  trees.findClades(min(cthreshold,nthreshold));
  trees.findNamedClades(nthreshold,maxTops);

  trees.printMeanTree(cout);
  // trees.printNamedClades(cout);
  // trees.printTreeTopologies(cout,numTrees);
  // trees.printProbableTreeClades(cout,nthreshold);
  // trees.printTrans(cout,maxTops);
  // trees.printCladeTreeTopologies(cout,numTrees);
  // trees.printClades(cout,cthreshold);
  // trees.printDistanceMatrices(cout);
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
    else if(string(argv[i]) == string("-") + string(names[j])[0] 
	    || string(argv[i]) == string("--") + string(names[j])) {
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
//  summarize(skip,numTrees,cthreshold,nthreshold,maxTops,filenames,argc,argv);
  test();
  return 0;
}

