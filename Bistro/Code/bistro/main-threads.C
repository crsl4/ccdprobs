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

  milliseconds ms2 = duration_cast< milliseconds >( system_clock::now().time_since_epoch() );

  // Find Jukes-Cantor pairwise distances
//  cerr << alignment.getNumTaxa() << endl;
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
  cout << endl << "Tree topology:" << endl;
  cout << jctree.makeTopologyNumbers() << endl << endl;
  cerr << endl << "Tree topology:" << endl;
  cerr << jctree.makeTopologyNumbers() << endl << endl;

  // for ( int i=0; i<4; ++i )
  // {
  //   double logFoo;
  //   jctree.randomEdges(alignment,model,rng,logFoo,true);
  // }

  
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

  cerr << "Running MCMC to estimate Q matrix ..." << endl;
  // Run MCMC on tree to estimate Q matrix parameters
  //   initial Q matrix

  vector<double> p_init(4,0.25);
  vector<double> s_init(6,0.1);
  s_init[1] = 0.3;
  s_init[4] = 0.3;

  QMatrix q_init(p_init,s_init);

  // burnin
  cerr << "burn-in:" << endl;
  //  q_init.mcmc(alignment,jctree,(MCMC_Q_BURN),alignment.getNumSites(),rng);
  jctree.mcmc(q_init,alignment,(MCMC_Q_BURN),alignment.getNumSites(),rng);
  cerr << endl << " done." << endl;

  // mcmc to get final Q
  cerr << "sampling:" << endl;
  //  q_init.mcmc(alignment,jctree,(MCMC_Q_SAMPLE),alignment.getNumSites(),rng);
  jctree.mcmc(q_init,alignment,(MCMC_Q_SAMPLE),alignment.getNumSites(),rng);
  cerr << endl << " done." << endl;

//  QMatrix model(parameters.getStationaryP(),parameters.getSymmetricQP());
  QMatrix model(q_init.getStationaryP(),q_init.getSymmetricQP());
  model.setMcmcVarP(q_init.getMcmcVarP());
  model.setMcmcVarQP(q_init.getMcmcVarQP());
  cerr << endl << " done." << endl;

  milliseconds ms6 = duration_cast< milliseconds >( system_clock::now().time_since_epoch() );


  // Recalculate pairwise distances using estimated Q matrix (TODO: add site rate heterogeneity)
  cerr << "Finding initial GTR pairwise distances ...";
  MatrixXd gtrDistanceMatrix(alignment.getNumTaxa(),alignment.getNumTaxa());
  alignment.calculateGTRDistances(model,jcDistanceMatrix,gtrDistanceMatrix);
  cerr << " done." << endl;

  cout << "GTR Distance Matrix:" << endl;
  cout << endl << gtrDistanceMatrix << endl << endl;

  milliseconds ms7 = duration_cast< milliseconds >( system_clock::now().time_since_epoch() );

  MatrixXd gtrDistanceMatrixCopy = gtrDistanceMatrix;
  Tree gtrtree(gtrDistanceMatrixCopy);
  gtrtree.reroot(1);
  gtrtree.sortCanonical();

  cout << "GTR Tree" << endl;
  cout << gtrtree.makeTreeNumbers() << endl;

  for ( int i=0; i<4; ++i )
  {
    double logFoo;
    gtrtree.randomEdges(alignment,model,rng,logFoo,true);
  }

  cout << "MLE Branch Lengths" << endl;
  cout << gtrtree.makeTreeNumbers() << endl;

  // Do bootstrap with new pairwise distances
  if( parameters.getNumBootstrap() == 0)
    parameters.setNumBootstrap(1000);
  cerr << "Beginning " << parameters.getNumBootstrap() << " bootstrap replicates ..." << endl;
  map<string,int> topologyToCountMap;
  map<string,int> topologyToParsimonyScoreMap;
  map<string,double> topologyToLogLikScoreMap;
  int minimumParsimonyScore = -1;
  int maximumLogLikScore = -1;
  MatrixXd bootDistanceMatrix(alignment.getNumTaxa(),alignment.getNumTaxa());
  vector<int> weights(alignment.getNumSites());
  cerr << '|';
  for ( int b=0; b<parameters.getNumBootstrap(); ++b )
  {
    if ( (b+1) % (parameters.getNumBootstrap() / 100) == 0 )
      cerr << '*';
    if ( (b+1) % (parameters.getNumBootstrap() / 10) == 0 )
      cerr << '|';
    alignment.setBootstrapWeights(weights,rng);
    alignment.calculateGTRDistancesUsingWeights(weights,model,gtrDistanceMatrix,bootDistanceMatrix);
    Tree bootTree(bootDistanceMatrix);
    bootTree.reroot(1); //warning: if 1 changes, need to change makeBinary if called after
    bootTree.makeBinary();
    bootTree.sortCanonical();

    string top = bootTree.makeTopologyNumbers();
    // add parsimony score to map if it is not already there
    // update minimum parsimony score if new minimum found
    if ( topologyToParsimonyScoreMap.find(top) == topologyToParsimonyScoreMap.end() )
    {
      int score = bootTree.parsimonyScore(alignment);
      topologyToParsimonyScoreMap[ top ] = score;
      if ( score < minimumParsimonyScore || b==0 )
	minimumParsimonyScore = score;
    }
    // add likelihood score to map if not there already
    // update minimum loglik
    if ( topologyToLogLikScoreMap.find(top) == topologyToLogLikScoreMap.end() )
    {
      bootTree.unroot();
      MatrixXd gtrDistanceMatrixCopy2(alignment.getNumTaxa(),alignment.getNumTaxa());
      gtrDistanceMatrixCopy2 = gtrDistanceMatrix;
      bootTree.setNJDistances(gtrDistanceMatrixCopy2,rng);
      int score = bootTree.logLikelihoodScore(alignment, model);
      topologyToLogLikScoreMap[ top ] = score;
      if ( score > maximumLogLikScore || b==0 )
	maximumLogLikScore = score;
    }
    topologyToCountMap[ top ]++;
  }
  cerr << endl << "done." << endl;

  milliseconds ms8 = duration_cast< milliseconds >( system_clock::now().time_since_epoch() );

  map<string,double> topologyToWeightMap; // not normalized; normalization taken care of during alias creation
  map<string,double> topologyToLogLikWeightMap; // not normalized; normalization taken care of during alias creation
  for ( map<string,int>::iterator m=topologyToCountMap.begin(); m != topologyToCountMap.end(); ++m )
  {
    topologyToWeightMap[ (*m).first ] =
      (*m).second * exp( parameters.getParsimonyScale() *
			 (minimumParsimonyScore - topologyToParsimonyScoreMap[ (*m).first ]) );
    topologyToLogLikWeightMap[ (*m).first ] =
      (*m).second * exp( parameters.getParsimonyScale() *
			 (topologyToLogLikScoreMap[ (*m).first ] - maximumLogLikScore) );
  }

  cout << endl << "Topology counts to file" << endl;
  {
    string topologyCountsFile = parameters.getOutFileRoot() + ".topCounts";
    ofstream topCounts(topologyCountsFile.c_str());
    map<string,double>::iterator wm=topologyToWeightMap.begin();
    map<string,double>::iterator loglwm=topologyToLogLikWeightMap.begin();
    topCounts << "tree count parsimonyWt parsimonyScore parsimonyDiff loglikWt loglikScore loglikDiff" << endl;
    for ( map<string,int>::iterator cm=topologyToCountMap.begin(); cm != topologyToCountMap.end(); ++cm )
    {
      Tree t((*cm).first);
      t.relabel(alignment);
      t.unroot();
      t.reroot(1);
      t.sortCanonical();
      string top = t.makeTopologyNumbers();
      topCounts << top << " " << setw(5) << (*cm).second << " ";
      topCounts << setw(10) << setprecision(4) << fixed << (*wm).second;
      topCounts << " " << setw(5) << topologyToParsimonyScoreMap[ (*cm).first ] << " " << setw(4) << minimumParsimonyScore - topologyToParsimonyScoreMap[ (*cm).first ];
      topCounts << " " << setw(10) << setprecision(4) << fixed << (*loglwm).second;
      topCounts << " " << setw(5) << topologyToLogLikScoreMap[ (*cm).first ] << " " << setw(4) << topologyToLogLikScoreMap[ (*cm).first ] - maximumLogLikScore << endl;
      ++wm;
      ++loglwm;
    }
    topCounts << endl;
    topCounts.close();
  }

  vector<int> taxaNumbers;
  vector<string> taxaNames;
  alignment.getTaxaNumbersAndNames(taxaNumbers,taxaNames);

  milliseconds ms9 = duration_cast< milliseconds >( system_clock::now().time_since_epoch() );

  CCDProbs<int> ccd(topologyToCountMap,taxaNumbers,taxaNames);
  CCDProbs<double> ccdParsimony(topologyToWeightMap,taxaNumbers,taxaNames);
  CCDProbs<double> ccdLogLik(topologyToLogLikWeightMap,taxaNumbers,taxaNames);
  // write map out to temp files to check
  string originalSmapFile = parameters.getOutFileRoot() + "-nopars.smap";
  string originalTmapFile = parameters.getOutFileRoot() + "-nopars.tmap";
  string parsimonySmapFile = parameters.getOutFileRoot() + "-pars.smap";
  string parsimonyTmapFile = parameters.getOutFileRoot() + "-pars.tmap";
  string loglikSmapFile = parameters.getOutFileRoot() + "-loglik.smap";
  string loglikTmapFile = parameters.getOutFileRoot() + "-loglik.tmap";
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
  smap.open(loglikSmapFile);
  ccdLogLik.writeCladeCount(smap);
  smap.close();
  tmap.open(loglikTmapFile);
  ccdLogLik.writePairCount(tmap);
  tmap.close();
  milliseconds ms10 = duration_cast< milliseconds >( system_clock::now().time_since_epoch() );

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
    cerr << "loglik weights " << parameters.getLoglikWt() << endl;
    cerr << "eta " << parameters.getEta() << endl;

    for ( int i=0; i<cores; ++i )
    {
      cerr << "core = " << i << endl;
      if ( !parameters.getReweight() )
	threads.push_back(thread(randomTrees<int>,i,startTreeNumber[i], startTreeNumber[i+1], ref(logwt0[i]), ref(maxLogW[i]), ccd, ref(*(vrng[i])), ref(alignment), ref(gtrDistanceMatrix), ref(model), ref(parameters), ref(topologymm[i])));
      else
	if( parameters.getLoglikWt() )
	  threads.push_back(thread(randomTrees<double>,i,startTreeNumber[i], startTreeNumber[i+1], ref(logwt0[i]), ref(maxLogW[i]), ccdLogLik, ref(*(vrng[i])), ref(alignment), ref(gtrDistanceMatrix), ref(model), ref(parameters), ref(topologymm[i])));
	else
	  threads.push_back(thread(randomTrees<double>,i,startTreeNumber[i], startTreeNumber[i+1], ref(logwt0[i]), ref(maxLogW[i]), ccdParsimony, ref(*(vrng[i])), ref(alignment), ref(gtrDistanceMatrix), ref(model), ref(parameters), ref(topologymm[i])));
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
  cout << ms1.count() - ms0.count() << endl;
  cout << ms2.count() - ms1.count() << endl;
  cout << ms3.count() - ms2.count() << endl;
  cout << ms4.count() - ms3.count() << endl;
  cout << ms5.count() - ms4.count() << endl;
  cout << ms6.count() - ms5.count() << endl;
  cout << ms7.count() - ms6.count() << endl;
  cout << ms8.count() - ms7.count() << endl;
  cout << ms9.count() - ms8.count() << endl;
  cout << ms10.count() - ms9.count() << endl;
  cout << ms11.count() - ms10.count() << endl;

  return 0;
}
