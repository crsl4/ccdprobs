
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
#include "cladegraph.h"

#include "Eigen/Core"
#include "Eigen/Eigenvalues"
#include "Eigen/SVD"

using namespace std;
using namespace std::chrono;
using namespace Eigen;

// summarize stuff:
//#include "Trees.h"
//bool debug = false;
//const char* names[]= {"skip", "trees", "cthreshold", "nthreshold", "maxtopologies"};
//const char* desc[] = {"skipped_lines", "number_of_trees_to_print", "threshold_for_clades",
//		      "threshold_for_named_clades",
//                      "max_tree_topologies" };
//// changing default values for tree, cthreshold, and maxtopologies
//const char* defaults[] = {"0", "100", "0.00001", ".80", "100"};
//const int skipField=0, numTreesField=1, cthresholdField=2, nthresholdField=3, maxTopsField=4;

bool comparePairStringDouble(const pair<string, double>  &p1, const pair<string, double> &p2)
{
  return p1.second < p2.second;
}

// will use x as the mean to the dirichlet, and v as the variance for each entry
// we calculate the scale inside
// E(yi) = alpha_i/sum{alpha_i}=x_i (because sum{x_i}=1)
// Var(yi) = alpha_i*(sum{alpha_i}-alpha_i)/(sum{alpha_i})^2*(sum{alpha_i}+1)
// if alpha_i=scale*x_i => Var(yi) = xi*(1-xi)/(scale+1) = vi => scale = xi(1-xi)/vi
// this function is not used in mcmc, that one used is in model.C
VectorXd dirichletProposalDensityScale(VectorXd x,double scale,double& logProposalDensity,mt19937_64& rng)
{
//  cout << "vector x: " << x.transpose() << endl;
  VectorXd alpha(x.size());
  VectorXd y(x.size());
  double sum = 0;
  double sumalpha = 0;
  for ( int i=0; i<x.size(); ++i )
  {
    alpha(i) = scale*x(i);
    if (VERBOSE)
      cout << "alpha: " << alpha(i) << endl;
    gamma_distribution<double> rgamma(alpha(i),1.0);
    y(i) = rgamma(rng);
    if (VERBOSE)
      cout << "proposed y: " << y(i) << endl;
    sum += y(i);
    sumalpha += alpha(i);
  }
  if (VERBOSE)
    cout << "sum: " << sum << endl;
  for ( int i=0; i<x.size(); ++i )
  {
    y(i) /= sum;
    if (VERBOSE)
      cout << "proposed y: " << y(i) << endl;
  }
  for ( int i=0; i<x.size(); ++i )
  {
    logProposalDensity += ( - lgamma(alpha(i)) + (alpha(i)-1)*log(y(i)) ) ;
  }
  logProposalDensity += lgamma(sumalpha);
  return y;
}

// arguments made reference to avoid copying in each thread
template<typename T>
void randomTrees(int coreID, int indStart, int indEnd, vector<double>& logwt, double& maxLogWeight, CCDProbs<T> ccd, mt19937_64& rng, Alignment& alignment, MatrixXd& gtrDistanceMatrix, QMatrix& q_init, Parameter& parameters, multimap<string,double>& topologyToLogweightMMap)
{
  string outFile = parameters.getOutFileRoot() + "---" + to_string(indStart) + "-" + to_string(indEnd-1) + ".out";
  ofstream f(outFile.c_str());
  string treeBLFile = parameters.getOutFileRoot() + "---" + to_string(indStart) + "-" + to_string(indEnd-1) + ".treeBL";
  ofstream treebl(treeBLFile.c_str());
  string treeBLFileSort = parameters.getOutFileRoot() + "---" + to_string(indStart) + "-" + to_string(indEnd-1) + ".treeBLsort";
  ofstream treeblsort(treeBLFileSort.c_str());
   f << "tree logl logTop logBL logPrior logQ logWt pi1 pi2 pi3 pi4 s1 s2 s3 s4 s5 s6" << endl;

   // debug variable
   int count = 0;

// calculate the scale for P and QP: fixit, this is done twice, also in main
   double scaleP = 100000000;
   double temp;
   VectorXd x = q_init.getStationaryP();
   VectorXd v = q_init.getMcmcVarP();
   if(VERBOSE)
   {
     cout << "p in q_init: " << q_init.getStationaryP() << endl;
     cout << "var(p) in q_init: " << q_init.getMcmcVarP() << endl;
   }
   for ( int i=0; i<x.size(); ++i )
   {
     temp = x(i)*(1-x(i))/v(i);
     if(temp < scaleP)
       scaleP = temp;
   }

   double scaleQP = 100000000;
   x = q_init.getSymmetricQP();
   v = q_init.getMcmcVarQP();
   if(VERBOSE)
   {
     cout << "s in q_init: " << q_init.getSymmetricQP() << endl;
     cout << "var(s) in q_init: " << q_init.getMcmcVarQP() << endl;
   }
   for ( int i=0; i<x.size(); ++i )
   {
     temp = x(i)*(1-x(i))/v(i);
     if( temp < scaleQP )
       scaleQP = temp;
   }

   if(VERBOSE)
     cout << "inside randomTrees, scaleP: " << scaleP << " scaleQP: " << scaleQP << endl;

   for ( int k=indStart; k<indEnd; ++k )
   {
     if(indStart == 0)
     {
       if ( indEnd > 99 && (k+1) % (indEnd / 100) == 0 )
	 cerr << '*';
       if ( indEnd > 9 && (k+1) % (indEnd / 10) == 0 )
	 cerr << '|';
     }
     if(VERBOSE)
       cout << "------------------" << endl;
     // here we need to sample a Q
     double logQ = 0;
     double otherscale = 1;
     if( parameters.getFixedQ() )
       otherscale = 10000000;
     if ( VERBOSE )
     {
       cout << "logQ before p_star = " << logQ << endl;
       cout << "eta: " << parameters.getEta() << endl;
     }

     VectorXd p_star = dirichletProposalDensityScale(q_init.getStationaryP(), parameters.getEta()*otherscale*scaleP, logQ, rng);
     if ( VERBOSE )
     {
       cout << "logQ after p_star = " << logQ << endl;
       cout << "p_star: " << p_star.transpose() << endl;
       cout << "scaleP: " << scaleP << endl;
     }
     VectorXd s_star = dirichletProposalDensityScale(q_init.getSymmetricQP(), parameters.getEta()*otherscale*scaleQP, logQ, rng);
     if ( VERBOSE )
     {
       cout << "logQ after s_star = " << logQ << endl;
       cout << "s_star: " << s_star.transpose() << endl;
       cout << "scaleQP: " << scaleQP << endl;
     }
     if( parameters.getFixedQ() ) // not used anymore
       logQ = 0;
     QMatrix model(p_star,s_star);

     double logTopologyProbability=0;
     string treeString = ccd.randomTree(rng,logTopologyProbability);
     Tree tree(treeString,alignment);
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
       // list<Node*> nodeList;
       // tree.postorderCherryNodeList(nodeList);
       // cout << "Postorder Node List:";
       // for ( list<Node*>::iterator p=nodeList.begin(); p!= nodeList.end(); ++p )
       // 	 cout << " " << (*p)->getNumber();
       cout << endl << flush;
     }

     for ( int i=0; i<parameters.getNumMLE(); ++i )
     {
       tree.mleLengths(alignment,model);
     }

     if(VERBOSE)
     {
       cout << "After MLE passes: " << endl << flush;
       tree.print(cout);
       tree.printDerivatives(cout,alignment,model);
       cout << tree.makeTreeNumbers() << endl << flush;
       cout << endl << flush;
     }

     double logBL = 0;

     if( parameters.getIndependent() )
     {//	  cout << "Branch lengths sampled independently" << endl;
       tree.randomEdges(alignment,model,rng,logBL);
     }
     else
     {
//	  cout << "Branch lengths sampled jointly in 2D" << endl;
       tree.generateBranchLengths(alignment,model,rng, logBL, parameters.getEta());
     }

     if(VERBOSE)
     {
       cout << "After generating branch lengths: " << endl;
       cout << tree.makeTreeNumbers() << endl;
     }
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
     if ( VERBOSE )
     {
       cout << "logBranchLengthPriorDensity = " << logBranchLengthPriorDensity << endl;
       cout << "logLik = " << logLik << endl;
       cout << "logBL = " << logBL << endl;
       cout << "logTopologyProbability = " << logTopologyProbability << endl;
       cout << "logQ = " << logQ << endl;
     }

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

  milliseconds ms2 = duration_cast< milliseconds >( system_clock::now().time_since_epoch() );

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

  // -------------- Find Jukes-Cantor pairwise distances and JC tree -----------------------------
  cerr << "Finding initial Jukes-Cantor pairwise distances ...";
  MatrixXd jcDistanceMatrix(alignment.getNumTaxa(),alignment.getNumTaxa());
  alignment.calculateJCDistances(jcDistanceMatrix);
  cerr << " done." << endl;
  cout << "Jukes-Cantor Distance Matrix:" << endl;
  cout << endl << jcDistanceMatrix << endl << endl;
  milliseconds ms3 = duration_cast< milliseconds >( system_clock::now().time_since_epoch() );

  // Find Initial Neighbor-joining tree
  cerr << "Finding initial neighbor-joining tree ...";
  MatrixXd jcDistanceMatrixCopy(alignment.getNumTaxa(),alignment.getNumTaxa());
  jcDistanceMatrixCopy = jcDistanceMatrix;
  Tree jctree(jcDistanceMatrixCopy);
  jctree.reroot(1); //warning: if 1 changed, need to change makeBinary if called after
  jctree.sortCanonical();
  cerr << " done." << endl;
  cout << endl << "JC Tree topology:" << endl;
  cout << jctree.makeTopologyNumbers() << endl << endl;
  cerr << endl << "JC Tree topology:" << endl;
  cerr << jctree.makeTopologyNumbers() << endl << endl;
  milliseconds ms4 = duration_cast< milliseconds >( system_clock::now().time_since_epoch() );

  // --------------------------- Use JC tree or input tree, and NJ branch lengths  ---------------
  string treetext; // the tree to use in MCMC will depend on whether a tree was input or not
  if ( !parameters.getTopology().empty() )
    treetext = parameters.getTopology();
  else
    treetext = jctree.makeTopologyNumbers();
  cerr << "treetext: " << treetext << endl;
  Tree starttree(treetext,alignment);
  starttree.unroot();
  cerr << "Start tree: " << starttree.makeTopologyNumbers() << endl;
  jcDistanceMatrixCopy = jcDistanceMatrix;
  starttree.setNJDistances(jcDistanceMatrixCopy,rng);
  cerr << "Setting NJ distances to tree: " << endl;
  cerr << starttree.makeTreeNumbers() << endl;

  milliseconds ms5 = duration_cast< milliseconds >( system_clock::now().time_since_epoch() );

  // ---------------------- Estimate Q with base frequencies and pairwise counts -----------------
  vector<double> p_init = alignment.baseFrequencies();
  cerr << "Estimated base frequencies: " << p_init.at(0) << "," << p_init.at(1) << "," << p_init.at(2) << "," << p_init.at(3) << endl;

  VectorXd s_pairwise(6);
  s_pairwise = alignment.averagePairwiseS();
  cerr << "Estimated rates: " << endl;
  cerr << s_pairwise.transpose() << endl;

  // checking to see if p and s were initialized at command line

  Vector4d p0;
  VectorXd s0(6);

  if ( parameters.getEnteredP() )
    p0 = convert(parameters.getStationaryP());
  else
    p0 = convert(p_init);
  if ( parameters.getEnteredS() )
    s0 = convert(parameters.getSymmetricQP());
  else
    s0 = s_pairwise;

  QMatrix q_init(p0,s0);

// initial mcmc var (just so that there are not zeros)
  Vector4d vp0;
  for ( int j=0; j<4; ++j )
    vp0(j) = p0(j)*(1-p0(j))/alignment.getNumSites();

  VectorXd vs0(6);
  for ( int j=0; j<6; ++j )
    vs0(j) = s0(j)*(1-s0(j))/alignment.getNumSites();

  q_init.setMcmcVarP(vp0);
  q_init.setMcmcVarQP(vs0);

// ------------------------------ Put MLE branch lengths to tree -------------------------
  for ( int i=0; i<4; ++i )
  {
    starttree.mleLengths(alignment,q_init);
  }
  cout << "MLE Branch Lengths" << endl;
  cout << starttree.makeTreeNumbers() << endl;

// MCMC
  if ( parameters.getDoMCMC() )
  {
    cerr << "Running MCMC to estimate Q matrix ..." << endl;

    // burnin
    cerr << "burn-in:" << endl;
    unsigned int mcmcGenerations = parameters.getNumMCMC();
    unsigned int mcmcBurn =  mcmcGenerations / 10;
    starttree.mcmc(q_init,alignment,mcmcBurn,alignment.getNumSites(),rng,true);
    cerr << endl << " done." << endl;

    // mcmc to get final Q
    cerr << "sampling:" << endl;
    starttree.mcmc(q_init,alignment,mcmcGenerations,alignment.getNumSites(),rng,false);
    cerr << endl << " done." << endl;
  }
  cerr << "After MCMC block" << endl;

// calculate the scale for P and QP (just to write it down, because it is calculated in randomTrees, but threads cannot
// save to file)
  double scaleP = 100000000;
  double temp;
  VectorXd x = q_init.getStationaryP();
  VectorXd v = q_init.getMcmcVarP();
  cout << "mean p: " << x.transpose() << endl;
  cout << "var p: " << v.transpose() << endl;

  for ( int i=0; i<x.size(); ++i )
  {
    temp = x(i)*(1-x(i))/v(i);
    cout << "Scale for p" << i << " is " << temp << endl;
    if(temp < scaleP)
      scaleP = temp;
  }
  cout << "Dirichlet for P scale: " << scaleP << endl;
  double scaleQP = 100000000;
  x = q_init.getSymmetricQP();
  v = q_init.getMcmcVarQP();
  cout << "mean s: " << x.transpose() << endl;
  cout << "var s: " << v.transpose() << endl;

  for ( int i=0; i<x.size(); ++i )
  {
    temp = x(i)*(1-x(i))/v(i);
    cout << "Scale for s" << i << " is " << temp << endl;
    if( temp < scaleQP )
      scaleQP = temp;
  }
  cout << "Dirichlet for QP scale: " << scaleQP << endl;

  // Initial Q, either from naive estimate or MCMC on fixed tree with MLE BL
  //QMatrix model_init(parameters.getStationaryP(),parameters.getSymmetricQP());
  QMatrix model_init(q_init.getStationaryP(),q_init.getSymmetricQP());
  model_init.setMcmcVarP(q_init.getMcmcVarP());
  model_init.setMcmcVarQP(q_init.getMcmcVarQP());
  cerr << endl << " done." << endl;

// will use the starting estimate of Q for the bootstrap:
//  QMatrix model_init(convert(p_init),s_pairwise);

  milliseconds ms6 = duration_cast< milliseconds >( system_clock::now().time_since_epoch() );

  // Recalculate pairwise distances using estimated Q matrix (TODO: add site rate heterogeneity)
  cerr << "Finding initial GTR pairwise distances ...";
  MatrixXd gtrDistanceMatrix(alignment.getNumTaxa(),alignment.getNumTaxa());
  alignment.calculateGTRDistances(model_init,jcDistanceMatrix,gtrDistanceMatrix);
  cerr << " done." << endl;

  cout << "GTR Distance Matrix:" << endl;
  cout << endl << gtrDistanceMatrix << endl << endl;

  milliseconds ms7 = duration_cast< milliseconds >( system_clock::now().time_since_epoch() );

// not normalized maps; normalization taken care of during alias creation
  map<string,int> topologyToCountMap;
  map<string,double> topologyToWeightMap;
  map<string,double> topologyToDistanceWeightMap;

// --------------------------------------- Bootstrap of topologies from ccdprobs -----------------------
  if( parameters.getTopology().empty() ) //only do bootstrap if not fixed topology as input
  {
    cerr << "Do boostrap with new pairwise distances"<< endl;
   // Do bootstrap with new pairwise distances
    if( parameters.getNumBootstrap() == 0)
      parameters.setNumBootstrap(1000);
    cerr << "Beginning " << parameters.getNumBootstrap() << " bootstrap replicates ..." << endl;
    map<string,int> topologyToParsimonyScoreMap;
    int minimumParsimonyScore = -1;
    MatrixXd bootDistanceMatrix(alignment.getNumTaxa(),alignment.getNumTaxa());
    vector<int> weights(alignment.getNumSites());
    cerr << '|';
    // output files
    string bootstrapTreesFile = parameters.getOutFileRoot() + ".bootstrap";
    ofstream bootstrapTrees(bootstrapTreesFile.c_str());
    string bootstrapTreesFileBL = parameters.getOutFileRoot() + ".bootstrapBL";
    ofstream bootstrapTreesBL(bootstrapTreesFileBL.c_str());
    string bootstrapTreesDistances = parameters.getOutFileRoot() + ".distances";
    ofstream bootstrapTreesDist(bootstrapTreesDistances.c_str());
    // vector of strings trees
    // change to using push_back() to avoid keeping bad trees
    vector<string> bootstrapStrings;
    string meanTree;
    { // block to contain badTrees
      int badTrees=0;
      for ( int b=0; b<parameters.getNumBootstrap(); ++b )
      {
	if ( parameters.getNumBootstrap() > 99 && (b+1) % (parameters.getNumBootstrap() / 100) == 0 )
	  cerr << '*';
	if ( parameters.getNumBootstrap() > 9 && (b+1) % (parameters.getNumBootstrap() / 10) == 0 )
	  cerr << '|';
	alignment.setBootstrapWeights(weights,rng);
	alignment.calculateGTRDistancesUsingWeights(weights,model_init,gtrDistanceMatrix,bootDistanceMatrix);
	Tree bootTree(bootDistanceMatrix);
	bootTree.reroot(1); //warning: if 1 changes, need to change makeBinary if called after
	bootTree.sortCanonical();
	string unrootedTreeString = bootTree.makeTreeNumbers(); // for the new distance method
	bootTree.makeBinary();
	bootTree.sortCanonical();

	string top = bootTree.makeTopologyNumbers();
	string topBL = bootTree.makeTreeNumbers();
	// write bootstrap tree to file
	if ( topBL.find("nan") != string::npos )
	{
	  ++badTrees;
	  continue;
	}
	bootstrapTrees << top << endl;
	bootstrapTreesBL << topBL << endl;
	bootstrapStrings.push_back(unrootedTreeString);

	// add parsimony score to map if it is not already there
	// update minimum parsimony score if new minimum found
	if ( topologyToParsimonyScoreMap.find(top) == topologyToParsimonyScoreMap.end() )
	{
	  int score = bootTree.parsimonyScore(alignment);
	  topologyToParsimonyScoreMap[ top ] = score;
	  if ( score < minimumParsimonyScore || b==0 )
	    minimumParsimonyScore = score;
	}
	topologyToCountMap[ top ]++;
      } // end of bootstrap
      cerr << endl;
      if ( badTrees > 0 )
	cerr << "Warning: found " << badTrees << " bootstrap tree" << (badTrees > 1 ? "s" : "") << "with nan edge lengths." << endl;
    }
    bootstrapTrees.close();
    bootstrapTreesBL.close();

    cerr << endl << "    ... done." << endl;

    // mean tree
    CladeGraph cladeGraph;
    cladeGraph.findMeanTree(bootstrapStrings,alignment);
    meanTree = cladeGraph.getMeanTree();
    cerr << "meanTree = " << meanTree << endl;
    cout << "meanTree = " << meanTree << endl;

    string meanTreeFile = parameters.getOutFileRoot() + ".meanTree";
    ofstream meanTreeFileStream(meanTreeFile.c_str());
    meanTreeFileStream << meanTree << endl;
    meanTreeFileStream.close();
//    cerr << "written mean tree to a file ok" << endl;

    // calculate distance from trees to mean tree
    Tree mtree(meanTree,alignment);
    cerr << "mean tree topology = " << mtree.makeTopologyNumbers() << endl;
    int badTrees = 0;
    for ( vector<string>::iterator t = bootstrapStrings.begin(); t!=bootstrapStrings.end(); ++t )
    {
//      cerr << "bootstrap tree: " << (*t) << endl;
      // check if tree string contains nan
      string foo = *t;
      if ( foo.find("nan") != string::npos )
      {
	++badTrees;
	continue;
      }
      Tree* boottree = new Tree(*t, alignment);
      double d = mtree.distance(boottree);
      cout << boottree->makeTopologyNumbers() << " " << d << endl;
      boottree->reroot(1); //warning: if 1 changes, need to change makeBinary if called after
      boottree->makeBinary();
      boottree->sortCanonical();
      string topDist = boottree->makeTopologyNumbers();
      bootstrapTreesDist << topDist << " " << d << endl;
      if ( topologyToDistanceWeightMap.find(topDist) == topologyToDistanceWeightMap.end() )
	topologyToDistanceWeightMap[ topDist ] = exp(-parameters.getWeightScale() * d);
      else
	topologyToDistanceWeightMap[ topDist ] += exp(-parameters.getWeightScale() * d);
      delete boottree;
    }
    bootstrapTreesDist.close();

//    for ( map<string,double>::iterator m=topologyToDistanceWeightMap.begin(); m != topologyToDistanceWeightMap.end(); ++m )
//      cerr << "topology to distance weight, tree: " << (*m).first << " distance weight: " << (*m).second << endl;

    if ( badTrees > 0 )
      cerr << "Warning: found " << badTrees << " bootstrap trees with nan edge lengths." << endl;

    for ( map<string,int>::iterator m=topologyToCountMap.begin(); m != topologyToCountMap.end(); ++m )
    {
      topologyToWeightMap[ (*m).first ] =
	(*m).second * exp( 0.5 *
			   (minimumParsimonyScore - topologyToParsimonyScoreMap[ (*m).first ]) );
    }
    cout << endl << "Topology counts to file" << endl;
    cerr << endl << "Topology counts to file" << endl;
    {
      string topologyCountsFile = parameters.getOutFileRoot() + ".topCounts";
      ofstream topCounts(topologyCountsFile.c_str());
      map<string,double>::iterator wm=topologyToWeightMap.begin();
      map<string,double>::iterator distwm=topologyToDistanceWeightMap.begin();
      topCounts << "tree count parsimonyWt parsimonyScore parsimonyDiff";
      topCounts << " distanceWt" << endl;
      for ( map<string,int>::iterator cm=topologyToCountMap.begin(); cm != topologyToCountMap.end(); ++cm )
      {
	Tree t((*cm).first,alignment);
	t.unroot();
	t.reroot(1);
	t.sortCanonical();
	string top = t.makeTopologyNumbers();
	topCounts << top << " " << setw(5) << (*cm).second << " " << flush;
	topCounts << setw(10) << setprecision(4) << fixed << (*wm).second;
	topCounts << " " << setw(5) << topologyToParsimonyScoreMap[ (*cm).first ] << " " << setw(4) << minimumParsimonyScore - topologyToParsimonyScoreMap[ (*cm).first ];
	topCounts << " " << flush;
	topCounts << setw(10) << setprecision(4) << fixed << (*distwm).second << flush;
	++distwm;
	topCounts << endl;
	++wm;
      }
      topCounts << endl;
      topCounts.close();
    }
    cerr << "after writing to files" << endl;

    // computing tree distance matrix, only if there are fewer than 100 bootstrap trees (fixit: sample)
    if( parameters.getNumBootstrap() < 101 )
    {
      cerr << "Fewer than 100 bootstrap trees, so we will compute the tree distance matrix" << endl;
      Tree mtree(meanTree,alignment);
      MatrixXd treeDistanceMatrix = MatrixXd::Zero(parameters.getNumBootstrap()+1,parameters.getNumBootstrap()+1);
      for ( int i = 0; i<bootstrapStrings.size(); ++i)
      {
	string fooi = bootstrapStrings[i];
	if ( fooi.find("nan") != string::npos )
	{
	  ++badTrees;
	  continue;
	}
	Tree* treei = new Tree(fooi, alignment);
	for ( int j = i+1; j<bootstrapStrings.size(); ++j)
	{
	  string fooj = bootstrapStrings[j];
	  if ( fooj.find("nan") != string::npos )
	  {
	    ++badTrees;
	    continue;
	  }
	  Tree* treej = new Tree(fooj, alignment);
	  double dd = treei->distance(treej);
	  treeDistanceMatrix(i,j) = dd;
	  treeDistanceMatrix(j,i) = dd;
	}
	double d = mtree.distance(treei);
	treeDistanceMatrix(i,bootstrapStrings.size()) = d;
	treeDistanceMatrix(bootstrapStrings.size(),i) = d;
      }
      cerr << "Found " << badTrees << " for the pairwise distance matrix" << endl;
//      cerr << treeDistanceMatrix << endl;
      cout << "Tree distance matrix: " << endl << endl;
      cout << treeDistanceMatrix << endl;
    }
  }
  else // topology input, so only one tree in the map
  {

    Tree tree(parameters.getTopology(),alignment);
    string t = tree.makeTopologyNumbers();
    cerr << "Will not do bootstrap, using input tree instead: " << t << endl;
    topologyToCountMap[ t ]++;
    topologyToWeightMap[ t ] = 1.0;
    topologyToDistanceWeightMap[ t ] = 1.0;
  }


  milliseconds ms8 = duration_cast< milliseconds >( system_clock::now().time_since_epoch() );
  vector<int> taxaNumbers;
  vector<string> taxaNames;
  alignment.getTaxaNumbersAndNames(taxaNumbers,taxaNames);
  milliseconds ms9 = duration_cast< milliseconds >( system_clock::now().time_since_epoch() );

// --------------------- Clade distribution from bootstrap sample ------------------------------------
  CCDProbs<int> ccd(topologyToCountMap,taxaNumbers,taxaNames);
  CCDProbs<double> ccdParsimony(topologyToWeightMap,taxaNumbers,taxaNames);
  CCDProbs<double> ccdDist(topologyToDistanceWeightMap,taxaNumbers,taxaNames);
  // write map out to temp files to check
  string originalSmapFile = parameters.getOutFileRoot() + "-nopars.smap";
  string originalTmapFile = parameters.getOutFileRoot() + "-nopars.tmap";
  string parsimonySmapFile = parameters.getOutFileRoot() + "-pars.smap";
  string parsimonyTmapFile = parameters.getOutFileRoot() + "-pars.tmap";
  string distSmapFile = parameters.getOutFileRoot() + "-dist.smap";
  string distTmapFile = parameters.getOutFileRoot() + "-dist.tmap";
  ofstream smap(originalSmapFile.c_str());
  ccd.writeCladeCount(smap);
  smap.close();
  ofstream tmap(originalTmapFile.c_str());
  ccd.writePairCount(tmap);
  tmap.close();
  smap.open(parsimonySmapFile);
  ccdParsimony.writeCladeCount(smap);
  smap.close();
  tmap.open(parsimonyTmapFile);
  ccdParsimony.writePairCount(tmap);
  tmap.close();
  smap.open(distSmapFile);
  ccdDist.writeCladeCount(smap);
  smap.close();
  tmap.open(distTmapFile);
  ccdDist.writePairCount(tmap);
  tmap.close();
  milliseconds ms10 = duration_cast< milliseconds >( system_clock::now().time_since_epoch() );

  if(parameters.getOnlyBootstrap())
    return 0;

  // --------------------- Random sample of trees ------------------------------------
  // using the same model_init (either from initial MCMC or naive estimate)
  QMatrix model(model_init.getStationaryP(),model_init.getSymmetricQP());
  model.setMcmcVarP(model_init.getMcmcVarP());
  model.setMcmcVarQP(model_init.getMcmcVarQP());

  if ( parameters.getNumRandom() > 0 )
  {
    int numRandom = parameters.getNumRandom();

    unsigned int cores;
    if( parameters.getNumThreads() == 0 )
    {
      cores = thread::hardware_concurrency();
      if( cores == 0 ) //hardware_concurrency could return 0
	cores = 1;
    }
    else
      cores = parameters.getNumThreads();
    if( numRandom < cores )
      cores = numRandom; // do not create more threads than trees
    // i want to check that cores not > than hardware_concurrency?
    cerr << "Generating " << numRandom << " random trees in " << cores << " cores:" << endl;
    int div = numRandom / cores;
    int remainder = numRandom % cores;
    vector<int> startTreeNumber(cores);
    for ( int i=0; i<cores; ++i )
      startTreeNumber[i] = i*div + (i < remainder ? i : remainder);
    startTreeNumber.push_back(numRandom);

    // generating the seeds for each core
    uniform_int_distribution<> rint_orig(0,4294967295);
    unsigned int initial_seed = rint_orig(rng);
    cout << "Seed = " << initial_seed << endl;
    minstd_rand seed_rng(initial_seed);
    uniform_int_distribution<> rint(0,4294967295);
    vector<unsigned int> seeds(cores);
    vector<mt19937_64*> vrng;
    for ( vector<unsigned int>::iterator p=seeds.begin(); p!=seeds.end(); ++p )
      {
    	*p = rint(seed_rng);
    	vrng.push_back( new mt19937_64(*p) );
      }
    cout << "Seeds per core: " << endl;
    for ( vector<unsigned int>::iterator p=seeds.begin(); p!=seeds.end(); ++p )
      cout << setw(15) << *p << endl;

    vector<thread> threads;
    vector< multimap<string,double> > topologymm(cores); //vector of multimaps
    vector< vector<double> > logwt0(cores);
    vector<double> maxLogW(cores); //vector of maxlogweight
    cerr << "jointMLE " << parameters.getJointMLE() << endl;
    cerr << "fixedQ " << parameters.getFixedQ() << endl;
    cerr << "eta " << parameters.getEta() << endl;

    for ( int i=0; i<cores; ++i )
    {
      cerr << "core = " << i << endl;
      if ( !parameters.getReweight() )
	threads.push_back(thread(randomTrees<double>,i,startTreeNumber[i], startTreeNumber[i+1], ref(logwt0[i]), ref(maxLogW[i]), ccdParsimony, ref(*(vrng[i])), ref(alignment), ref(gtrDistanceMatrix), ref(model), ref(parameters), ref(topologymm[i])));
      else
	threads.push_back(thread(randomTrees<double>,i,startTreeNumber[i], startTreeNumber[i+1], ref(logwt0[i]), ref(maxLogW[i]), ccdDist, ref(*(vrng[i])), ref(alignment), ref(gtrDistanceMatrix), ref(model), ref(parameters), ref(topologymm[i])));
    }

    for(auto &t : threads){
      t.join();
    }
    // for ( int i=0; i<cores; ++i )
    //   threads.at(i).join();
    cerr << endl << "done." << endl;

    // combine all maxlogweight
    double maxLogWeight = maxLogW[0];
    for(int i=1; i<cores; ++i)
      {
	if(maxLogW[i] > maxLogWeight)
	  maxLogWeight = maxLogW[i];
      }
    cerr << "maxLogWeight: " << maxLogWeight << endl;

    // combine vector of all logwt
    vector<double> logwt(numRandom,0);
    int k = 0;
    for( vector< vector<double> >::iterator p=logwt0.begin(); p != logwt0.end(); ++p)
      {
	for(vector<double>::iterator q=(*p).begin(); q != (*p).end(); ++q)
	  {
	    logwt[k] = *q;
	    //cerr << "logwt[" << k << "]= " << logwt[k] << endl;
	    ++k;
	  }
      }

    // combine topology maps to one
    multimap<string,double> topologyToLogweightMMap;
    for ( vector< multimap<string,double>>::iterator p=topologymm.begin(); p!=topologymm.end(); ++p) //for every multimap in topologymm
      {
	for ( multimap<string,double >::iterator q=(*p).begin(); q!= (*p).end(); ++q) // traverse this multimap (*p)
	  {
	    topologyToLogweightMMap.insert( pair<string,double>(q->first,q->second) ) ;
	  }
      }

    vector<double> wt(numRandom,0);
    double sum=0;
    for ( int k=0; k<numRandom; ++k )
    {
      //      cerr << "logwt[k]: " << logwt[k] << "for k= " << k << endl;
      wt[k] = exp(logwt[k] - maxLogWeight);
      //cerr << "wt[k]: " << wt[k] << endl;
      sum += wt[k];
    }
    double essInverse=0;
    for ( int k=0; k<numRandom; ++k )
    {
      wt[k] /= sum;
      essInverse += wt[k]*wt[k];
    }
    cout << "ESS = " << fixed << setprecision(2) << 1.0/essInverse << ", or "
	 << setprecision(2) << 100.0 / essInverse / numRandom << " percent." << endl;
    cerr << "ESS = " << fixed << setprecision(2) << 1.0/essInverse << ", or "
	 << setprecision(2) << 100.0 / essInverse / numRandom << " percent." << endl;

    // substract max from multimap
    for ( multimap<string,double >::iterator p=topologyToLogweightMMap.begin(); p!= topologyToLogweightMMap.end(); ++p) {
      p->second -= maxLogWeight; // need to test this
    }

    // map from topology to unnormalized weight
    map<string,double> topologyToUnnormalizedWeightMap;
    for ( multimap<string,double >::iterator p=topologyToLogweightMMap.begin(); p!= topologyToLogweightMMap.end(); ++p) {
      pair< multimap<string,double >::iterator,multimap<string,double >::iterator> ret = topologyToLogweightMMap.equal_range(p->first);
      double sumWt = 0;
      for ( multimap<string,double >::iterator q=ret.first; q!= ret.second; ++q) {
	sumWt += exp(q->second);
      }
      //cout << "topology: " << (p->first) << " and sum of weights " << sumWt << endl;
      topologyToUnnormalizedWeightMap[ p->first ] = sumWt;
    }

    CCDProbs<double> splits(topologyToUnnormalizedWeightMap,taxaNumbers,taxaNames);
    string splitsWeightsFile = parameters.getOutFileRoot() + ".splits";
    ofstream splitsWt(splitsWeightsFile.c_str());
    splits.writeCladeCountOrdered(splitsWt, essInverse);
    splitsWt.close();

    string topologyPPFile = parameters.getOutFileRoot() + ".topPP";
    ofstream topPP(topologyPPFile.c_str());

    vector< pair<string,double> > v;
    copy(topologyToUnnormalizedWeightMap.begin(), topologyToUnnormalizedWeightMap.end(), back_inserter(v));
    sort(v.begin(), v.end(), comparePairStringDouble);
    double prob;
    double se;
    topPP << "tree prob crudeSE(ESS=" << fixed << setprecision(2) << 1.0/essInverse << ")" << endl;

    for (vector<pair<string,double>>::reverse_iterator i = v.rbegin(); i != v.rend(); ++i ) {
      prob = i->second / sum;
      se = sqrt(prob * (1-prob) * essInverse);
      topPP << (i->first) << " " << fixed << setprecision(4) << prob << " " << se << endl;
    }

    // for ( map<string,double >::iterator p=topologyToUnnormalizedWeightMap.begin(); p!= topologyToUnnormalizedWeightMap.end(); ++p) {
    //   prob = (p->second)/sum;
    //   se = sqrt(prob * (1-prob) * essInverse);
    //   topPP << (p->first) << fixed << setprecision(4) << prob << " " << se << endl;
    // }
  }

  milliseconds ms11 = duration_cast< milliseconds >( system_clock::now().time_since_epoch() );

  cout << "Times: " << endl;
  // cout << "Processing command-line arguments: " << (ms1.count() - ms0.count())/1000 << "seconds" << endl;
  // cout << "Reading FASTA: " << (ms2.count() - ms1.count())/1000 << "seconds" << endl;
  // cout << "Jukes-Cantor distances: " << (ms3.count() - ms2.count())/1000 << "seconds" << endl;
  // cout << "Getting initial JC NJ tree: " << (ms4.count() - ms3.count())/1000 << "seconds" << endl;
  // cout << "Setting NJ branch lengths on tree: " << (ms5.count() - ms4.count())/1000 << "seconds" << endl;
  cout << "MCMC: " << (ms6.count() - ms5.count())/1000 << "seconds" << endl;
//  cout << "GTR distance matrix: " << (ms7.count() - ms6.count())/1000 << "seconds" << endl;
  cout << "Bootstrap: " << (ms8.count() - ms7.count())/1000 << "seconds" << endl;
//  cout << ms9.count() - ms8.count() << endl;
  cout << "Estimate ccdprobs from bootstrap: " << (ms10.count() - ms9.count())/1000 << "seconds" << endl;
  cout << "Importance sampling: " << (ms11.count() - ms10.count())/1000 << "seconds" << endl;

  return 0;
}
