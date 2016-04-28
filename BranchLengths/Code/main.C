#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <string>
#include <random>

#include "parameter.h"
#include "sequence.h"
#include "tree.h"
#include "model.h"

using namespace std;

void checkDistances(Tree& tree,Alignment& alignment,QMatrix& model)
{
  cout << tree.makeTopologyNumbers() << endl;
  for ( int i=0; i<11; ++i )
  {
    for ( int j=i+1; j<12; ++j )
    {
      Node* na = tree.getNode(i);
      Edge* ea = na->getEdge(i);
      Node* nb = tree.getNode(j);
      Edge* eb = nb->getEdge(j);
      double t = tree.mleDistance(alignment,na,ea,nb,eb,model);
      double logl,dlogl,ddlogl;
      tree.partialPathCalculations(t,alignment,na,ea,nb,eb,model,logl,dlogl,ddlogl,true);
      cout << setw(3) << i
           << setw(3) << j << " "
           << setprecision(6) << fixed << t
           << " " << setprecision(6) << logl
           << " " << setprecision(6) << setw(12) << dlogl
           << " " << setprecision(6) << setw(18) << ddlogl
           << " " << setprecision(6) << setw(10) << -1/ddlogl
           << " " << setprecision(4) << setw(8) << -t*t*ddlogl
           << " " << setprecision(4) << setw(10) << -t*ddlogl << endl;
    }
  }
}

int main(int argc, char* argv[])
{
  // Read command line and process parameters
  Parameter parameters;
  parameters.processCommandLine(argc,argv);

  // Read in sequences from FASTA file
  Alignment alignment(parameters.getSequenceFileName());
  alignment.summarize(cout);

  // Initialize random number generator
  if ( parameters.getSeed() == 0 )
  {
    random_device rd;
    parameters.setSeed(rd());
  }
  mt19937_64 rng(parameters.getSeed());

  // set up Q matrix
  QMatrix model(parameters.getStationaryP(),parameters.getSymmetricQP());

  // create tree
  if ( !parameters.getTopology().empty() )
  {
    Tree tree(parameters.getTopology());
    tree.relabel(alignment);
    tree.randomBranchLengths(rng,20); // only here to test likelihood calculations
//    tree.print(cout);
//    tree.setSitePatterns(sequences);

    // log-likelihood calculation
    checkDistances(tree,alignment,model);
    tree.randomize(rng);
    cout << tree.makeTopologyNumbers() << endl;
    tree.generateBranchLengths(alignment,model);
    tree.print(cout);
  }
  return 0;
}
