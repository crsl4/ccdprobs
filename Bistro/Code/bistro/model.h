#ifndef __MODEL_H
#define __MODEL_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <random>

#include "Eigen/Core"
#include "Eigen/Eigenvalues"

#include "parameter.h"
#include "sequence.h"

class Tree;
class MCMCStats;

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
  Vector6d symmetricQP; //AC,AG,AT,CG,CT,GT?
  Vector4d stationaryP; //A,C,G,T?
  Vector4d lambda;
  Matrix4d V;
  Matrix4d Vinv;
  Matrix4d Q;
  Vector6d mcmcVarQP;
  Vector4d mcmcVarP;
  Vector4d alphaForGenDirichletPi;
  Vector4d lambdaForGenDirichletPi;
  Vector6d alphaForGenDirichletS;
  Vector6d lambdaForGenDirichletS;
public:
  void completeConstruction();
  QMatrix(vector<double>,vector<double>); // p and s, assumes sum to 1 and positive
  QMatrix(Vector4d,VectorXd); // p and s, assumes sum to 1 and positive
  QMatrix(Vector4d,VectorXd, Vector4d, VectorXd);
  Vector4d vectorExp(double t);
  Matrix4d getTransitionMatrix(double);
  Matrix4d getQ();
  Matrix4d getQP(double);
  Matrix4d getQQP(double);
  Vector4d getStationaryP() const { return stationaryP; }
  Vector6d getSymmetricQP() const { return symmetricQP; }
  Vector4d getStatP() const { return stationaryP; } //fixit: why this?
  void resetStationaryP(Vector4d p)
  {
    stationaryP = p;
    completeConstruction();
  }
  void resetSymmetricQP(VectorXd s)
  {
    symmetricQP = s;
    completeConstruction();
  }
  void reset(Vector4d p,VectorXd s)
  {
    stationaryP = p;
    symmetricQP = s;
    completeConstruction();
  }
  void mcmc(Alignment&,Tree&,int,double,mt19937_64&);
  Vector6d getMcmcVarQP() const { return mcmcVarQP; }
  Vector4d getMcmcVarP() const { return mcmcVarP; }
  void setMcmcVarQP(VectorXd v) { mcmcVarQP = v; }
  void setMcmcVarP(Vector4d v) { mcmcVarP = v; }
  void resetAfterMCMC(MCMCStats&,unsigned int);
  void calculateAlphaLambdaForGenDirichlet();
  void genDirichletProposal(double&,mt19937_64&,Vector4d&,Vector6d&);
  void copyAlphaLambda(QMatrix&);
  Vector4d getAlphaPi() { return alphaForGenDirichletPi; }
  Vector4d getLambdaPi() { return lambdaForGenDirichletPi; }
  Vector6d getAlphaS() { return alphaForGenDirichletS; }
  Vector6d getLambdaS() { return lambdaForGenDirichletS; }
};

VectorXd dirichletProposal(VectorXd ,double ,double& ,mt19937_64& );

VectorXd convert(vector<double>);

vector<double> convert(VectorXd);

#endif
