#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <string>
#include <random>
#include <fstream>
#include <sstream>

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
  bool verbose = false;
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
      double t = tree.mleDistance(alignment,na,ea,nb,eb,model, verbose);
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

  cout << "Seed: " << parameters.getSeed() << endl;

  // set up Q matrix
  QMatrix model(parameters.getStationaryP(),parameters.getSymmetricQP());

  // create vector of weights
  unsigned int sampleSize;
  if ( parameters.getSampleSize() == 0 )
    sampleSize = 1000;
  else
    sampleSize = parameters.getSampleSize();

  bool verbose = parameters.getVerbose();
  bool mvnormal = parameters.getMvnormal();
  //cout << "Verbose: " << verbose << endl;

  cout << "Sample size: " << sampleSize << endl;

  VectorXd logw(sampleSize);

  // create tree
  if ( !parameters.getTopology().empty() )
  {
    Tree tree(parameters.getTopology());
    tree.relabel(alignment);
//    tree.randomBranchLengths(rng,20); // only here to test likelihood calculations
    tree.print(cout);
//    tree.setSitePatterns(sequences);

    // log-likelihood calculation
    if ( verbose )
    {
      cerr << "Check distances" << endl;
      checkDistances(tree,alignment,model);
      cerr << "----------------------------------------------------------------" << endl;
    }
	
    if(verbose)
      cout << tree.makeTopologyNumbers() << endl;
    int errors = 0;
    ofstream logwfile;
    logwfile.open("logw.txt"); //warning: only works for 4taxa tree
    //    logwfile << "n1,n2,bl1,n1,n2,bl2,n1,n2,bl3,n1,n2,bl4,n1,n2,bl5,loglik,logprior,logdens,logweight" << endl;
    // files to study form of lik:
    ofstream table3D;
    ofstream table2D;
    ofstream par3D;
    ofstream par2D;
    //table3D.open("table3D.txt");
    //table2D.open("table2D.txt");
    par3D.open("par3D.txt");
    par2D.open("par2D.txt");

    par3D << "mu1,mu2,mu3,s1,s2,s3,s4,s5,s6,s7,s8,s9,a1,b1,a2,b2,a3,b3" << endl;
    par2D << "mu1,mu2,s1,s2,s3,s4,a1,b1,a2,b2" << endl;

    for(int i = 0; i < sampleSize; i++)
      {
	cout << "------------- rep " << i << " -------------" << endl;
	try
	  {
	    tree.randomize(rng);
	    if ( verbose )
	      tree.print(cerr);
	    tree.clearProbMaps();
	    tree.generateBranchLengths(alignment,model,rng, verbose, mvnormal, par3D, par2D);
	    double weight = tree.calculateWeight(alignment, model, 0.05, verbose, logwfile);
	    logw(i) = weight;
	    tree.print(cout);
	  }
	catch ( ... )
	  {
	    cerr << "Found error in replicate number: " << i << endl;
	    errors++;
	    logw(i) = 0;
	  }
      }
    logwfile.close();
    par3D.close();
    par2D.close();

    // double average = 0;
    // for(int i = 0; i < sampleSize; i++)
    //   {
    // 	average += logw(i);
    //   }
    // average = average / sampleSize;
    // cout << "Average logw: " << average << endl;

    // // normalize vector
    // for(int i = 0; i < sampleSize; i++)
    //   {
    // 	logw(i) = logw(i) - average;
    // 	logw(i) = exp(logw(i));
    //   }

    // double suma = 0;
    // for(int i = 0; i < sampleSize; i++)
    //   {
    // 	suma += logw(i);
    //   }

    // logw = logw / suma;
    // double ess = 0;
    // for(int i = 0; i < sampleSize; i++)
    //   {
    // 	ess += log(i)*logw(i);
    //   }
    // ess = 1/ess;
    // ess = ess/sampleSize;
    // cout << "Effective sample size: " << ess << endl;
    cout << "Errors: " << errors << endl;
  }
  return 0;
}
