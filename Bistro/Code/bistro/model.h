#ifndef __MODEL_H
#define __MODEL_H

#include <iostream>
#include <iomanip>
#include <vector>

#include "Eigen/Core"
#include "Eigen/Eigenvalues"

#include "parameter.h"

using namespace std;
using namespace Eigen;

typedef Matrix<double, 6, 1> Vector6d;

// Q matrix for GTR Model
// Use parameterization q(i,j) = s(i,j) / (2 * p(i)) for i != j
// S = diag(sqrt(p)) * Q * diag( 1/sqrt(p) ) is similar symmetric matrix; use for eigenvalue decomposition
// Q = V L Vinv where L is diagonal of eigenvalues and V is eigenvectors
// S = U L U^T where U are eigenvectors of S and U^T is both the transpose and inverse of U
// V = diag( 1/sqrt(p) ) * U, Vinv = U^T * diag(sqrt(p))
// P(t) = V * exp(L*t) * Vinv = diag( 1/sqrt(p) ) * U * exp(L*t) * U * diag ( sqrt(p) )
class QMatrix
{
private:
  Vector6d symmetricQP;
  Vector4d stationaryP;
  Vector4d lambda;
  Matrix4d V;
  Matrix4d Vinv;
  Matrix4d Q;
public:
  void completeConstruction();
  QMatrix(vector<double>,vector<double>); // p and s, assumes sum to 1 and positive
  QMatrix(Vector4d,VectorXd); // p and s, assumes sum to 1 and positive
  Vector4d vectorExp(double t);
  Matrix4d getTransitionMatrix(double);
  Matrix4d getQ();
  Matrix4d getQP(double);
  Matrix4d getQQP(double);
  Vector4d getStationaryP() const { return stationaryP; }
  Vector6d getSymmetricQP() const { return symmetricQP; }
  Vector4d getStatP() const { return stationaryP; }
};
  
#endif
