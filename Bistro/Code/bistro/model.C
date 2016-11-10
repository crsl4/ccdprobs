#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <cmath>
#include <random>

#include "Eigen/Core"
#include "Eigen/Eigenvalues"

#include "sequence.h"
#include "parameter.h"
#include "model.h"
#include "tree.h"

using namespace std;
using namespace Eigen;

void QMatrix::completeConstruction()
{
    // map to get correct element of s from (i,j) index
  map<pair<int,int>,int> ijMap;
  ijMap[ pair<int,int> (0,1) ] = 0;
  ijMap[ pair<int,int> (0,2) ] = 1;
  ijMap[ pair<int,int> (0,3) ] = 2;
  ijMap[ pair<int,int> (1,2) ] = 3;
  ijMap[ pair<int,int> (1,3) ] = 4;
  ijMap[ pair<int,int> (2,3) ] = 5;

  // Q matrix
  Matrix4d qmatrix;

  for ( int i=0; i<4; ++i )
  {
    qmatrix(i,i) = 0;
    for (int j=0; j<4; ++j )
    {
      if ( i < j )
      {
        qmatrix(i,j) = symmetricQP( ijMap[ pair<int,int> (i,j) ] ) / (2.0 * stationaryP(i));
        qmatrix(i,i) -= qmatrix(i,j);
      }
      else if ( i > j ) {
        qmatrix(i,j) = symmetricQP( ijMap[ pair<int,int> (j,i) ] ) / (2.0 * stationaryP(i));
        qmatrix(i,i) -= qmatrix(i,j);
      }
    }
  }

  Vector4d psqrt;
  Vector4d psqrtinv;
  for ( int i=0; i<4; ++i )
  {
    psqrt(i) = sqrt( stationaryP(i) );
    psqrtinv(i) = 1 / psqrt(i);
  }

  DiagonalMatrix<double,4> pSqrtDiag = psqrt.asDiagonal();
  DiagonalMatrix<double,4> pSqrtDiagInv = psqrtinv.asDiagonal();

  Matrix4d symmat = pSqrtDiag * qmatrix * pSqrtDiagInv;
  SelfAdjointEigenSolver<Matrix4d> es(symmat);
  lambda = es.eigenvalues();
  V = pSqrtDiagInv * es.eigenvectors();
  Vinv = es.eigenvectors().transpose() * pSqrtDiag;
}

QMatrix::QMatrix(vector<double> p,vector<double> s) // p and s, assumes sum to 1 and positive
{
  for ( int i=0; i<4; ++i )
    stationaryP(i) = p[i];
  for ( int i=0; i<6; ++i )
    symmetricQP(i) = s[i];

  completeConstruction();
}

QMatrix::QMatrix(Vector4d p,VectorXd s)
{
  stationaryP = p;
  symmetricQP = s;
  completeConstruction();
}

Vector4d QMatrix::vectorExp(double t)
{
  Vector4d x;
  for ( int i=0; i<4; ++i )
    x(i) = exp(lambda(i) * t);
  return x;
}

Matrix4d QMatrix::getTransitionMatrix(double t) // assumes t>=0
{
  Vector4d x = vectorExp(t);
  return V * x.asDiagonal() * Vinv;
}

Matrix4d QMatrix::getQ()
{
  return V * lambda.asDiagonal() * Vinv;
}

Matrix4d QMatrix::getQP(double t)
{
  Vector4d x = vectorExp(t);
  return V * lambda.asDiagonal() * x.asDiagonal() * Vinv;
}

Matrix4d QMatrix::getQQP(double t)
{
  Vector4d x = vectorExp(t);
  return V * lambda.asDiagonal() * lambda.asDiagonal() * x.asDiagonal() * Vinv;
}

VectorXd dirichletProposal(VectorXd x,double scale,double& logProposalRatio,mt19937_64& rng)
{
  VectorXd alpha(x.size());
  VectorXd y(x.size());
  double sum = 0;
  for ( int i=0; i<x.size(); ++i )
  {
    alpha(i) = scale*x(i) + 1;
    gamma_distribution<double> rgamma(alpha(i),1.0);
    y(i) = rgamma(rng);
    sum += y(i);
  }
  for ( int i=0; i<x.size(); ++i )
  {
    y(i) /= sum;
  }
  for ( int i=0; i<x.size(); ++i )
  {
    logProposalRatio += ( lgamma(alpha(i)) - lgamma(scale*y(i) + 1) + scale*(y(i)*log(x(i)) - x(i)*log(y(i))) );
  }
  return y;
}

// for the variance: https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
void QMatrix::mcmc(Alignment& alignment,Tree& tree,int numGenerations,double scale,mt19937_64& rng)
{
  tree.clearProbMaps();
  double currLogLikelihood = tree.calculate(alignment,*this);
  double sumAcceptP = 0;
  double sumAcceptS = 0;
  Vector4d avgP;
  Vector4d avgPold;
  Vector4d sP;
  Vector4d prod1;
  Vector4d prod2;
  Vector4d prod12;
  avgP << 0,0,0,0;
  avgPold << 0,0,0,0;
  sP << 0,0,0,0;
  prod1 << 0,0,0,0;
  prod2 << 0,0,0,0;
  prod12 << 0,0,0,0;
  VectorXd avgS(6);
  VectorXd avgSold(6);
  VectorXd sS(6);
  VectorXd prod3(6);
  VectorXd prod4(6);
  VectorXd prod34(6);
  avgS << 0,0,0,0,0,0;
  avgSold << 0,0,0,0,0,0;
  sS << 0,0,0,0,0,0;
  prod3 << 0,0,0,0,0,0;
  prod4 << 0,0,0,0,0,0;
  prod34 << 0,0,0,0,0,0;
  cerr << '|';
  for ( int i=0; i<numGenerations; ++i )
  {
    if ( (i+1) % (numGenerations / 100) == 0 )
      cerr << '*';
    if ( (i+1) % (numGenerations / 10) == 0 )
      cerr << '|';
    Vector4d x = getStationaryP();
    double logProposalRatio = 0;
    Vector4d y = dirichletProposal(x,scale,logProposalRatio,rng);
    QMatrix propQ(y,getSymmetricQP());
    tree.clearProbMaps();
    double propLogLikelihood = tree.calculate(alignment,propQ);
    double acceptProbability = exp(propLogLikelihood - currLogLikelihood + logProposalRatio);
    if ( acceptProbability > 1 )
      acceptProbability = 1;
    sumAcceptP += acceptProbability;
    uniform_real_distribution<double> runif(0,1);
    if ( runif(rng) < acceptProbability )
    {
      resetStationaryP(y);
      currLogLikelihood = propLogLikelihood;
    }
    //    avgP += getStationaryP();
    prod1 = getStationaryP() - avgPold;
    avgP = avgPold + prod1/(i+1);
    prod2 = getStationaryP() - avgP;
    for ( int i=0; i<4; ++i )
	prod12(i) = prod1(i)*prod2(i);
    sP += prod12;
    avgPold = avgP;
    VectorXd xx = getSymmetricQP();
    logProposalRatio = 0;
    VectorXd yy(6);
    yy = dirichletProposal(xx,scale,logProposalRatio,rng);
    propQ.reset(getStationaryP(),yy);
    tree.clearProbMaps();
    propLogLikelihood = tree.calculate(alignment,propQ);
    acceptProbability = exp(propLogLikelihood - currLogLikelihood + logProposalRatio);
    if ( acceptProbability > 1 )
      acceptProbability = 1;
    sumAcceptS += acceptProbability;
    if ( runif(rng) < acceptProbability )
    {
      resetSymmetricQP(yy);
      currLogLikelihood = propLogLikelihood;
    }
    //    avgS += getSymmetricQP();
    prod3 = getSymmetricQP() - avgSold;
    avgS = avgSold + prod3/(i+1);
    prod4 = getSymmetricQP() - avgS;
    for ( int i=0; i<6; ++i )
	prod34(i) = prod3(i)*prod4(i);
    sS += prod34;
    avgSold = avgS;
  }
  cout << "stationary acceptance: " << sumAcceptP / numGenerations << endl;
  cout << "symmetric acceptance: " << sumAcceptS / numGenerations << endl;
  //  avgP /= numGenerations;
  //  avgS /= numGenerations;
  sP /= (numGenerations-1);
  sS /= (numGenerations-1);
  cout << "avgP: " << avgP.transpose() << endl;
  cout << "sP: " << sP.transpose() << endl;
  cout << "avgS: " << avgS.transpose() << endl;
  cout << "sS: " << sS.transpose() << endl;
  reset(avgP,avgS);
  setMcmcVarP(sP);
  setMcmcVarQP(sS);
  cout << stationaryP.transpose() << endl;
  cout << symmetricQP.transpose() << endl;
  cout << endl << getQ() << endl << endl;
}

// need to double check
double gammaProposal(double x,double ddlogl,double& logProposalRatio,mt19937_64& rng)
{
  double lambday = -1 * x * ddlogl;
  double alphay = x*lambday;
  gamma_distribution<double> rgamma(alphay,1.0 / lambday);
  double y = rgamma(rng);
  double lambdax = -1 * y * ddlogl; //?? really ddlogl here too? maybe scale?
  double alphax = y*lambdax;
  logProposalRatio += alphax*log(lambdax) - lgamma(alphax) + (alphax-1)*log(x) - lambdax*x - alphax*log(lambdax) + lgamma(alphax) - (alphax-1)*log(x) + lambdax*x;
  return y;
}

void mcmc(QMatrix& Q,Alignment& alignment,Tree& tree,int numGenerations,double scale,mt19937_64& rng)
{
  tree.clearProbMaps();
  double currLogLikelihood = tree.calculate(alignment,Q);
  double sumAcceptP = 0;
  double sumAcceptS = 0;
  Vector4d avgP;
  Vector4d avgPold;
  Vector4d sP;
  Vector4d prod1;
  Vector4d prod2;
  Vector4d prod12;
  avgP << 0,0,0,0;
  avgPold << 0,0,0,0;
  sP << 0,0,0,0;
  prod1 << 0,0,0,0;
  prod2 << 0,0,0,0;
  prod12 << 0,0,0,0;
  VectorXd avgS(6);
  VectorXd avgSold(6);
  VectorXd sS(6);
  VectorXd prod3(6);
  VectorXd prod4(6);
  VectorXd prod34(6);
  avgS << 0,0,0,0,0,0;
  avgSold << 0,0,0,0,0,0;
  sS << 0,0,0,0,0,0;
  prod3 << 0,0,0,0,0,0;
  prod4 << 0,0,0,0,0,0;
  prod34 << 0,0,0,0,0,0;
  VectorXd avgBL(tree.getNumEdges());
  cerr << '|';
  for ( int i=0; i<numGenerations; ++i )
  {
    if ( (i+1) % (numGenerations / 100) == 0 )
      cerr << '*';
    if ( (i+1) % (numGenerations / 10) == 0 )
      cerr << '|';
    Vector4d x = Q.getStationaryP();
    double logProposalRatio = 0;
    Vector4d y = dirichletProposal(x,scale,logProposalRatio,rng);
    QMatrix propQ(y,Q.getSymmetricQP());
    tree.clearProbMaps();
    double propLogLikelihood = tree.calculate(alignment,propQ);
    double acceptProbability = exp(propLogLikelihood - currLogLikelihood + logProposalRatio);
    if ( acceptProbability > 1 )
      acceptProbability = 1;
    sumAcceptP += acceptProbability;
    uniform_real_distribution<double> runif(0,1);
    if ( runif(rng) < acceptProbability )
    {
      Q.resetStationaryP(y);
      currLogLikelihood = propLogLikelihood;
    }
    //    avgP += Q.getStationaryP();
    prod1 = Q.getStationaryP() - avgPold;
    avgP = avgPold + prod1/(i+1);
    prod2 = Q.getStationaryP() - avgP;
    for ( int i=0; i<4; ++i )
	prod12(i) = prod1(i)*prod2(i);
    sP += prod12;
    avgPold = avgP;
    VectorXd xx = Q.getSymmetricQP();
    logProposalRatio = 0;
    VectorXd yy(6);
    yy = dirichletProposal(xx,scale,logProposalRatio,rng);
    propQ.reset(Q.getStationaryP(),yy);
    tree.clearProbMaps();
    propLogLikelihood = tree.calculate(alignment,propQ);
    acceptProbability = exp(propLogLikelihood - currLogLikelihood + logProposalRatio);
    if ( acceptProbability > 1 )
      acceptProbability = 1;
    sumAcceptS += acceptProbability;
    if ( runif(rng) < acceptProbability )
    {
      Q.resetSymmetricQP(yy);
      currLogLikelihood = propLogLikelihood;
    }
    prod3 = Q.getSymmetricQP() - avgSold;
    avgS = avgSold + prod3/(i+1);
    prod4 = Q.getSymmetricQP() - avgS;
    for ( int i=0; i<6; ++i )
	prod34(i) = prod3(i)*prod4(i);
    sS += prod34;
    avgSold = avgS;
    // now with Q, we want to sample branch lengths
    tree.clearProbMaps();
    double xxx;
    double yyy;
    int j = 0;
    double logl,dlogl,ddlogl; //??? or scale??
    for ( vector<Edge*>::iterator e=tree.getEdges().begin(); e!=tree.getEdges().end(); ++e )
      {
	xxx = (*e)->getLength();
	logProposalRatio = 0;
	(*e)->calculate(xxx,alignment,Q,logl,dlogl,ddlogl);
	yyy = gammaProposal(xxx,ddlogl,logProposalRatio,rng);
	tree.clearProbMaps();
	propLogLikelihood = tree.calculate(alignment,Q); //for the whole tree?
	acceptProbability = exp(propLogLikelihood - currLogLikelihood + logProposalRatio);
	if ( acceptProbability > 1 )
	  acceptProbability = 1;
	//	sumAcceptBL += acceptProbability; // ??
	uniform_real_distribution<double> runif(0,1);
	if ( runif(rng) < acceptProbability )
	  {
	    (*e)->setLength(yyy);
	    currLogLikelihood = propLogLikelihood;
	  }
	avgBL(j) += (*e)->getLength();
	j++;
      }
  }
  cout << "stationary acceptance: " << sumAcceptP / numGenerations << endl;
  cout << "symmetric acceptance: " << sumAcceptS / numGenerations << endl;
  sP /= (numGenerations-1);
  sS /= (numGenerations-1);
  cout << "avgP: " << avgP.transpose() << endl;
  cout << "sP: " << sP.transpose() << endl;
  cout << "avgS: " << avgS.transpose() << endl;
  cout << "sS: " << sS.transpose() << endl;
  Q.reset(avgP,avgS);
  Q.setMcmcVarP(sP);
  Q.setMcmcVarQP(sS);
  cout << Q.getStationaryP().transpose() << endl;
  cout << Q.getSymmetricQP().transpose() << endl;
  cout << endl << Q.getQ() << endl << endl;
  int i = 0;
  for ( vector<Edge*>::iterator e=tree.getEdges().begin(); e!=tree.getEdges().end(); ++e )
    {
      (*e)->setLength(avgBL(i));
      i++;
    }
}

