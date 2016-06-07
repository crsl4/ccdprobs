#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <cmath>
#include <random>

#include "Eigen/Core"
#include "Eigen/Eigenvalues"

#include "parameter.h"
#include "model.h"

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

// x is current state; assume 0 < x[i] < 1 and sum x[i] = 1
// y is the returned proposed state with equivalent constraints
// proposal is Dirichlet distribution with alpha[i] = scale*x[i] + 0.5
// this helps prevent the state from getting sucked into a corner where x[i] = 0 for some i
//void mcmcProposeDirichlet(vector<double>& x,vector<double>& y,double scale,mt19937_64& rng)
//{
//  
//}
