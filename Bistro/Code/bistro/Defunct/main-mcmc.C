// Bistro: main.C
// Bret Larget and Claudia Solis-Lemus
// Copyright 2016
// May 20, 2016

#define MCMC_Q_BURN 100
#define MCMC_Q_SAMPLE 1000


#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <string>
#include <random>
#include <thread>
#include <chrono>

#include "parameter.h"
#include "sequence.h"
#include "alias.h"
#include "tree.h"
#include "model.h"
#include "ccdprobs.h"
#include "random.h"

#include "Eigen/Core"
#include "Eigen/Eigenvalues"
#include "Eigen/SVD"

using namespace std;
using namespace std::chrono;
using namespace Eigen;

bool comparePairStringDouble(const pair<string, double>  &p1, const pair<string, double> &p2)
{
  return p1.second < p2.second;
}

VectorXd dirichletProposalDensity(VectorXd x,double scale,double& logProposalDensity,mt19937_64& rng)
{
  VectorXd alpha(x.size());
  VectorXd y(x.size());
  double sum = 0;
  double sumalpha = 0;
  for ( int i=0; i<x.size(); ++i )
  {
    alpha(i) = scale*x(i) + 1;
    gamma_distribution<double> rgamma(alpha(i),1.0);
    y(i) = rgamma(rng);
    sum += y(i);
    sumalpha += alpha(i);
  }
  for ( int i=0; i<x.size(); ++i )
  {
    y(i) /= sum;
  }
  for ( int i=0; i<x.size(); ++i )
  {
    logProposalDensity += ( - lgamma(alpha(i)) + scale*x(i)*log(y(i)) ) ;
  }
  logProposalDensity += lgamma(sumalpha);
  return y;
}

// will use x as the mean to the dirichlet, and v as the variance for each entry
// we calculate the scale inside
// E(yi) = alpha_i/sum{alpha_i}=x_i (because sum{x_i}=1)
// Var(yi) = alpha_i*(sum{alpha_i}-alpha_i)/(sum{alpha_i})^2*(sum{alpha_i}+1)
// if alpha_i=scale*x_i => Var(yi) = xi*(1-xi)/(scale+1) = vi => scale = xi(1-xi)/vi 
VectorXd dirichletProposalDensityScale(VectorXd x,VectorXd v,double& logProposalDensity,mt19937_64& rng)
{
  //  cout << "vector x: " << x.transpose() << endl;
  //cout << "vector v: " << v.transpose() << endl;
  VectorXd alpha(x.size());
  VectorXd y(x.size());
  double sum = 0;
  double sumalpha = 0;
  double scale;
  for ( int i=0; i<x.size(); ++i )
    scale += x(i)*(1-x(i))/v(i);
  scale /= x.size(); //we need to use the same scale or we create bias
  cout << "Scale: " << scale << endl;  
  for ( int i=0; i<x.size(); ++i )
  {
    alpha(i) = scale*x(i) + 1;
    gamma_distribution<double> rgamma(alpha(i),1.0);
    y(i) = rgamma(rng);
    sum += y(i);
    sumalpha += alpha(i);
  }
  for ( int i=0; i<x.size(); ++i )
  {
    y(i) /= sum;
  }
  for ( int i=0; i<x.size(); ++i )
  {
    logProposalDensity += ( - lgamma(alpha(i)) + scale*x(i)*log(y(i)) ) ;
  }
  logProposalDensity += lgamma(sumalpha);
  return y;
}

// arguments made reference to avoid copying in each thread
template<typename T>
void randomTrees(int coreID, int indStart, int indEnd, vector<double>& logwt, double& maxLogWeight, CCDProbs<T> ccd, mt19937_64& rng, Alignment& alignment, MatrixXd& gtrDistanceMatrix, QMatrix& q_init, Parameter& parameters, multimap<string,double>& topologyToLogweightMMap)
{
//     cerr << "Random trees from " << indStart << " to " << indEnd << endl;

  string outFile = parameters.getOutFileRoot() + "---" + to_string(indStart) + "-" + to_string(indEnd-1) + ".out";
  ofstream f(outFile.c_str());
  string treeBLFile = parameters.getOutFileRoot() + "---" + to_string(indStart) + "-" + to_string(indEnd-1) + ".treeBL";
  ofstream treebl(treeBLFile.c_str());
  string treeBLFileSort = parameters.getOutFileRoot() + "---" + to_string(indStart) + "-" + to_string(indEnd-1) + ".treeBLsort";
  ofstream treeblsort(treeBLFileSort.c_str());
   f << "tree logl logTop logBL logPrior logQ logWt pi1 pi2 pi3 pi4 s1 s2 s3 s4 s5 s6" << endl;

   for ( int k=indStart; k<indEnd; ++k )
     {
       if(indStart == 0)
	 {
	   if ( indEnd > 99 && (k+1) % (indEnd / 100) == 0 )
	     cerr << '*';
	   if ( indEnd > 9 && (k+1) % (indEnd / 10) == 0 )
	     cerr << '|';
	 }
       // here we need to sample a Q
       double logQ = 0;
       double scale = 1;
       if( parameters.getFixedQ() )
	 scale = 10000000;
       VectorXd p_star = dirichletProposalDensityScale(q_init.getStationaryP(), q_init.getMcmcVarP(), logQ, rng);
       VectorXd s_star = dirichletProposalDensityScale(q_init.getSymmetricQP(), q_init.getMcmcVarQP(), logQ, rng);

       if( parameters.getFixedQ() )
	 logQ = 0;
       QMatrix model(p_star,s_star);

      double logTopologyProbability=0;
      string treeString;
      treeString = ccd.randomTree(rng,logTopologyProbability);

      Tree tree(treeString);
      tree.relabel(alignment);
      tree.unroot();
      MatrixXd gtrDistanceMatrixCopy(alignment.getNumTaxa(),alignment.getNumTaxa());
      gtrDistanceMatrixCopy = gtrDistanceMatrix;
      tree.setNJDistances(gtrDistanceMatrixCopy,rng);
      if( parameters.getRootFix() )
	tree.randomizeBL(rng);
      else
	tree.randomize(rng);

      if(VERBOSE)
	{
	  cout << coreID << " " << k << "th tree" << endl << flush;
	  tree.print(cout);
	  cout << tree.makeTreeNumbers() << endl << flush;
	  list<Node*> nodeList;
	  tree.postorderCherryNodeList(nodeList);
	  cout << "Postorder Node List:";
	  for ( list<Node*>::iterator p=nodeList.begin(); p!= nodeList.end(); ++p )
	    cout << " " << (*p)->getNumber();
	  cout << endl << flush;
	}
      double logBL = 0;

      if( parameters.getNumMLE() == 0 )
	parameters.setNumMLE(2);

      for ( int i=0; i<parameters.getNumMLE(); ++i )
	{
	  tree.randomEdges(alignment,model,rng,logBL,true);
//	  cout << tree.makeTreeNumbers() << endl;
	}
      if( parameters.getIndependent() )
	{
//	  cout << "Branch lengths sampled independently" << endl;
	  tree.randomEdges(alignment,model,rng,logBL,false);
	}
      else
	{
//	  cout << "Branch lengths sampled jointly in 2D" << endl;
	  tree.generateBranchLengths(alignment,model,rng, logBL, parameters.getJointMLE(), parameters.getEta(), parameters.getWeightMean());
	}
//      cout << tree.makeTreeNumbers() << endl;
      treebl << tree.makeTreeNumbers() << endl;
      double logBranchLengthPriorDensity = tree.logPriorExp( (PRIOR_MEAN) );
      // if(VERBOSE)
      // 	cout << "calculating the final loglik now before clearing map" << endl;
      // double logLik0 = tree.calculate(alignment, model);
      // tree.clearProbMaps(); //added just to see if this was missing
      if(VERBOSE)
	cout << "calculating the final loglik now without clearing map" << endl;
      double logLik = tree.calculate(alignment, model);
      double logWeight = logBranchLengthPriorDensity + logLik - logBL - logTopologyProbability - logQ;
      // string top0 = tree.makeTopologyNumbers();
      tree.reroot(1); //warning: if 1 changes, need to change makeBinary if called after
      tree.sortCanonical();
      treeblsort << tree.makeTreeNumbers() << endl;
      string top = tree.makeTopologyNumbers();
      f << top << " " << logLik << " " << " " << logTopologyProbability << " " << logBL << " " << logBranchLengthPriorDensity << " " << logQ << " " << logWeight << " ";
      f << model.getStationaryP().transpose() << " " << model.getSymmetricQP().transpose() << endl;
      logwt.push_back( logWeight );
      //logwt[k] = logWeight;
      if ( k==indStart || logWeight > maxLogWeight )
	maxLogWeight = logWeight;
      //add to the multimap for the topology the logweight
      topologyToLogweightMMap.insert( pair<string,double>(top,logWeight) ) ;
    }
  treebl.close();
  treeblsort.close();
}


int main(int argc, char* argv[])
{
  milliseconds ms0 = duration_cast< milliseconds >( system_clock::now().time_since_epoch() );
  // Read command line and process parameters
  cerr << "Processing command line ...";
  Parameter parameters;
  parameters.processCommandLine(argc,argv);
  cerr << " done." << endl;

  milliseconds ms1 = duration_cast< milliseconds >( system_clock::now().time_since_epoch() );
  // Read in sequences from FASTA file
  cerr << "Reading alignment from FASTA file ...";
  Alignment alignment(parameters.getSequenceFileName());
  cerr << " done." << endl;
  cout << endl << "Alignment summary:" << endl;
  alignment.summarize(cout);
  cout << endl;

// Read in initial tree topology
  Tree tree(parameters.getTopology());
  tree.setInitialEdgeLengths(0.1);
  
  milliseconds ms4 = duration_cast< milliseconds >( system_clock::now().time_since_epoch() );

  // Initialize random number generator
  cerr << "Initializing random number generator ...";
  if ( parameters.getSeed() == 0 )
  {
    random_device rd;
    parameters.setSeed(rd());
  }
  mt19937_64 rng(parameters.getSeed());
  cerr << " done." << endl;
  cout << "Seed = " << parameters.getSeed() << endl;

  milliseconds ms5 = duration_cast< milliseconds >( system_clock::now().time_since_epoch() );

  // Initialize Q matrix

  vector<double> p_init(4,0.25);
  vector<double> s_init(6,0.1);
  s_init[1] = 0.3;
  s_init[4] = 0.3;

  QMatrix q_init(p_init,s_init);

  // use parameters.getNumRandom() as the MCMC sample size
  // use parameters.getNumBootstrap() as the burnin size
  // file name for trees is root + .tre
  // file name for pi and s is root + .par

  string treeFile = parameters.getOutFileRoot() + ".tre";
  string parFile = parameters.getOutFileRoot() + ".par";
  ofstream treeStream(treeFile.c_str());
  ofstream parStream(parFile.c_str());
    
  cerr << "burn-in:" << endl;
  tree.mcmc(q_init,alignment,parameters.getNumBootstrap(),alignment.getNumSites(),rng,treeStream,parStream, true);
  cerr << endl << " done." << endl;
  cerr << "tree after burn-in: " << endl;
  cerr << tree.makeTreeNumbers() << endl;
  cerr << "average pi and s after burn-in" << endl;
  cerr << q_init.getStationaryP().transpose() << " " << q_init.getSymmetricQP().transpose() << endl;

  // mcmc
  cerr << "sampling:" << endl;
  tree.mcmc(q_init,alignment,parameters.getNumRandom(),alignment.getNumSites(),rng,treeStream,parStream, false);
  cerr << endl << " done." << endl;
  cerr << "average tree: " << endl;
  cerr << tree.makeTreeNumbers() << endl;
  cerr << "average pi and s after sampling" << endl;
  cerr << q_init.getStationaryP().transpose() << " " << q_init.getSymmetricQP().transpose() << endl;

  treeStream.close();
  parStream.close();
  
  return 0;
}
