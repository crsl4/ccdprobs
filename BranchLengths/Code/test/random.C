#include <random>
#include <iostream>

#include "Eigen/Core"
#include "Eigen/Eigenvalues"

#include "random.h"

using namespace std;
using namespace Eigen;


// code to generate multivariate normal random variables
// mu is vector of means
// vc is symmetric square positive definite variance covariance matrix
// if length of mu=n, then vc is n by n.
// returns n instances
// only compute cholesky decomposition once
// I think we only need it once ever, but think about it
// Idea is if we want variance-covariance matrix to be S,
//   then if S = LL^t where L is lower triangular
//   and if z is a vector of independent standard normals,
//   then if y = Lz,
//   E (y y^t) = E( Lz z^t L^t) = L E(z z^t) L^t = LL^t = S
// For the final implementation if we really will only use 2by2 and 3by3 matrices,
//   we may want different functions for each case
//   because the small matrices are optimized in the Eigen library
VectorXd multivariateNormal(VectorXd mu,MatrixXd vc,mt19937_64& rng)
{
  int r = vc.rows();
  VectorXd v( r );
  MatrixXd L( vc.llt().matrixL() );

  normal_distribution<double> rnorm;
  for ( int i=0; i<r; ++i )
    v( i ) = rnorm( rng );

  return mu + L * v;
}

double normal(double mu, double vc, mt19937_64& rng)
{
  normal_distribution<double> rnorm(mu,sqrt(vc));
  return rnorm( rng );
}

double gamma(double alpha, double beta, mt19937_64& rng)
{
  gamma_distribution<double> gam1(alpha,beta); 
  return gam1(rng);
}


Vector3d multivariateGamma3D(Vector3d mu,Matrix3d vc,mt19937_64& rng)
{
  Vector3d bl;
  Matrix3d L( vc.llt().matrixL() );
  // ------------ T1 ------------------
  double alpha1 = mu[0]*mu[0] / (L(0,0) * L(0,0));
  double lambda1 = mu[0] / (L(0,0) * L(0,0));
  bl[0] = gamma(alpha1,1/lambda1,rng); //c++ gamma has different parametrization
  // ------------ T2 ------------------
  double z1 = (bl[0] - mu[0]) / L(0,0);
  if( z1 < 0)
    {
      cerr << "z1 in multivariateGamma3D is negative" << endl;
      exit(1);
    }
  double alpha2 = (mu[1] + L(1,0)*z1)*(mu[1] + L(1,0)*z1) / (L(1,1) * L(1,1));
  double lambda2 = (mu[1] + L(1,0)*z1) / (L(1,1) * L(1,1));
  bl[1] = gamma(alpha2,1/lambda2,rng); //c++ gamma has different parametrization
  // ------------ T3 ------------------
  double z2 = (bl[1] - mu[1] - L(1,0)*z1) / L(1,1);
  if( z2 < 0)
    {
      cerr << "z2 in multivariateGamma3D is negative" << endl;
      exit(1);
    }
  double alpha3 = (mu[2] + L(2,0)*z1 + L(2,1)*z2)*(mu[2] + L(2,0)*z1 + L(2,1)*z2) / (L(2,2) * L(2,2));
  double lambda3 = (mu[2] + L(2,0)*z1 + L(2,1)*z2) / (L(2,2) * L(2,2));
  bl[2] = gamma(alpha3,1/lambda3,rng); //c++ gamma has different parametrization
  return bl;
}
