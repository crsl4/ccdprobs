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
#include "random.h"
#include "Eigen/Dense"

using namespace std;
using namespace Eigen;


void checkDistances(Tree& tree,Alignment& alignment,QMatrix& model)
{
  cout << tree.makeTopologyNumbers() << endl;
  int numTaxa = alignment.getNumTaxa();
  for ( int i=0; i < numTaxa-1; ++i )
  {
    for ( int j=i+1; j<numTaxa; ++j )
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

  cerr << "Seed: " << parameters.getSeed() << endl;

  // set up Q matrix
  QMatrix model(parameters.getStationaryP(),parameters.getSymmetricQP());

  // create vector of weights
  int sampleSize = 1000; //fixit: make it an argument
  VectorXd logw(sampleSize);

  // create tree
  if ( !parameters.getTopology().empty() )
  {
    Tree tree(parameters.getTopology());
    tree.relabel(alignment);
//    tree.randomBranchLengths(rng,20); // only here to test likelihood calculations
//    tree.print(cout);
//    tree.setSitePatterns(sequences);

    // log-likelihood calculation
    //checkDistances(tree,alignment,model);
    tree.randomize(rng);
    cout << tree.makeTopologyNumbers() << endl;
    int errors = 0;
    for(int i = 0; i < sampleSize; i++)
      {
	try
	  {
	    tree.generateBranchLengths(alignment,model,rng);
	    double weight = tree.calculateWeight(alignment, model, 0.05);
	    logw(i) = weight;
	    tree.print(cout);
	  }
	catch ( int e )
	  {
	    cerr << "Found error in replicate number: " << i << endl;
	    errors++;
	  }
      }
    // normalize vector
    double suma = 0;
    for(int i = 0; i < sampleSize; i++)
      {
	suma += log(i);
      }

    logw = logw / suma;
    double ess = 0;
    for(int i = 0; i < sampleSize; i++)
      {
	ess += log(i)*logw(i);
      }
    ess = 1/ess;
    ess = ess/sampleSize;
    cout << "Effective sample size: " << ess << endl;
    cout << "Errors: " << errors << endl;
  }
  return 0;
}
