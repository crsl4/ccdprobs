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

// old function: new one in tree.C
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
    for ( int j=0; j<4; ++j )
	prod12(j) = prod1(j)*prod2(j);
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
    for ( int j=0; j<6; ++j )
	prod34(j) = prod3(j)*prod4(j);
    sS += prod34;
    avgSold = avgS;
  }
  cout << "stationary acceptance: " << sumAcceptP / numGenerations << endl;
  cout << "symmetric acceptance: " << sumAcceptS / numGenerations << endl;
  //  avgP /= numGenerations;
  //  avgS /= numGenerations;
  sP /= (numGenerations-1);
  sS /= (numGenerations-1);
  reset(avgP,avgS);
  setMcmcVarP(sP);
  setMcmcVarQP(sS);
  cout << "avgP: " << avgP.transpose() << endl;
  cout << "sP: " << sqrt(sP(0)) << "," << sqrt(sP(1)) << "," << sqrt(sP(2)) << "," << sqrt(sP(3)) << endl;
  cout << "avgS: " << avgS.transpose() << endl;
  cout << "sS: " << sqrt(sS(0)) << "," << sqrt(sS(1)) << "," << sqrt(sS(2)) << "," << sqrt(sS(3)) << "," << sqrt(sS(4)) << "," << sqrt(sS(5))  << endl;
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

VectorXd convert(vector<double> x)
{
  VectorXd y(x.size());
  for ( int i=0; i<x.size(); ++i )
    y(i) = x[i];
  return y;
}

vector<double> convert(VectorXd x)
{
  vector<double> y(x.size());
  for ( int i=0; i<x.size(); ++i )
    y[i] = x(i);
  return y;
}

void QMatrix::resetAfterMCMC(MCMCStats& stats,unsigned int numGenerations)
{
  reset(stats.getAvgP(),stats.getAvgS());
  setMcmcVarP(stats.getSP()/numGenerations);
  setMcmcVarQP(stats.getSS()/numGenerations);
}

void QMatrix::calculateAlphaLambdaForGenDirichlet()
{
  // pi
  {
    Vector4d mu = getStationaryP();
    Vector4d v = getMcmcVarP();
    double total = 0;
    int k=mu.size();
    for ( int i=0; i<k; ++i )
    {
      lambdaForGenDirichletPi(i) = mu(i)*(1-mu(i))/v(i);
      total += lambdaForGenDirichletPi(i);
      alphaForGenDirichletPi(i) = mu(i)*lambdaForGenDirichletPi(i);
    }
    lambdaForGenDirichletPi *= (k/total);
  }
  // s
  {
    Vector6d mu = getSymmetricQP();
    Vector6d v = getMcmcVarQP();
    double total = 0;
    int k=mu.size();
    for ( int i=0; i<k; ++i )
    {
      lambdaForGenDirichletS(i) = mu(i)*(1-mu(i))/v(i);
      total += lambdaForGenDirichletS(i);
      alphaForGenDirichletS(i) = mu(i)*lambdaForGenDirichletS(i);
    }
    lambdaForGenDirichletS *= (k/total);
  }
}

void QMatrix::genDirichletProposal(double& logQ,mt19937_64& rng,Vector4d& p_star,Vector6d& s_star)
{
  // pi
  double sum=0;
  double sumAlpha=0;
  double sumLambdaX=0;
  double sumAlphaLogLambda=0;
  double sumLogGammaAlpha=0;
  double sumAlphaX=0;

  for ( int i=0; i<4; ++i )
  {
    gamma_distribution<double> rgamma(alphaForGenDirichletPi(i),1.0/lambdaForGenDirichletPi(i));
    p_star(i) = rgamma(rng);
    sum += p_star(i);
  }
  p_star /= sum;

  for ( int i=0; i<4; ++i )
  {
    sumAlpha += alphaForGenDirichletPi(i);
    sumLambdaX += lambdaForGenDirichletPi(i)*p_star(i);
    sumAlphaLogLambda += alphaForGenDirichletPi(i)*log( lambdaForGenDirichletPi(i) );
    sumLogGammaAlpha += lgamma(alphaForGenDirichletPi(i));
    sumAlphaX += (alphaForGenDirichletPi(i)-1) * log(p_star(i));
  }
  logQ += lgamma(sumAlpha) + sumAlphaLogLambda - sumLogGammaAlpha + sumAlphaX - sumAlpha*log(sumLambdaX);

  // s
  sum=0;
  sumAlpha=0;
  sumLambdaX=0;
  sumAlphaLogLambda=0;
  sumLogGammaAlpha=0;
  sumAlphaX=0;
  
  for ( int i=0; i<6; ++i )
  {
    gamma_distribution<double> rgamma(alphaForGenDirichletS(i),1.0/lambdaForGenDirichletS(i));
    s_star(i) = rgamma(rng);
    sum += s_star(i);
  }
  s_star /= sum;

  for ( int i=0; i<6; ++i )
  {
    sumAlpha += alphaForGenDirichletS(i);
    sumLambdaX += lambdaForGenDirichletS(i)*s_star(i);
    sumAlphaLogLambda += alphaForGenDirichletS(i)*log( lambdaForGenDirichletS(i) );
    sumLogGammaAlpha += lgamma(alphaForGenDirichletS(i));
    sumAlphaX += (alphaForGenDirichletS(i)-1) * log(s_star(i));
  }
  logQ += lgamma(sumAlpha) + sumAlphaLogLambda - sumLogGammaAlpha + sumAlphaX - sumAlpha*log(sumLambdaX);
}

void QMatrix::copyAlphaLambda(QMatrix& x)
{
  alphaForGenDirichletPi = x.getAlphaPi();
  lambdaForGenDirichletPi = x.getLambdaPi();
  alphaForGenDirichletS = x.getAlphaS();
  lambdaForGenDirichletS = x.getLambdaS();
}

