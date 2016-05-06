#include <random>

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
