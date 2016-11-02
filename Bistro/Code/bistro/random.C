#include <random>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "Eigen/Core"
#include "Eigen/Eigenvalues"

#include "random.h"
#include <cmath> // cbrt, tgamma

using namespace std;
using namespace Eigen;

void calculateAlphaLambda(double mu, double v, double eta, double& alpha, double& lambda)
{
  alpha = mu*mu / (v * v);
  lambda = mu / (v * v);
  if(alpha<1)
    {
      //cout << "Entering alpha<1 case" << endl;
      lambda = sqrt( (alpha+1.0)/alpha ) * lambda;
      alpha = (alpha + 1.0);
    }
  else
    {
      alpha = eta*alpha;
      lambda = eta*lambda;
    }
}

double gamma(double alpha, double b, mt19937_64& rng)
{
  gamma_distribution<double> gam1(alpha,b);
  return gam1(rng);
}

// eta is dividing the variance: var/eta
Vector3d multivariateGamma3D(Vector3d mu,Matrix3d vc,mt19937_64& rng, double& logdensity, double eta)
{
  Vector3d bl;
  Matrix3d L( vc.llt().matrixL() );
  // ------------ T1 ------------------
  double alpha1;
  double lambda1;
  double mu0;
  if( mu[0] < MIN_EDGE_LENGTH + TOL)
      mu0 = MIN_EDGE_LENGTH;
  else
      mu0 = mu[0];
  calculateAlphaLambda(mu0,L(0,0),eta,alpha1,lambda1);
  bl[0] = gamma(alpha1,1.0 / lambda1,rng); //c++ gamma has different parametrization

  // ------------ T2 ------------------
  double z1 = (bl[0] - (alpha1/lambda1)) / L(0,0);
  //  double z1 = (bl[0] - mu[0]) / L(0,0);
  //double z1 = (bl[0] - mu0) / L(0,0);
  double num = mu[1] + L(1,0)*z1;
  double alpha2;
  double lambda2;
  if( num < MIN_EDGE_LENGTH + TOL)
      num = MIN_EDGE_LENGTH;
  calculateAlphaLambda(num,L(1,1),eta,alpha2,lambda2);
  bl[1] = gamma(alpha2,1.0 / lambda2,rng); //c++ gamma has different parametrization

  // ------------ T3 ------------------
  double z2 = (bl[1] - (alpha2/lambda2)) / L(1,1);
  //  double z2 = (bl[1] - mu[1] - L(1,0)*z1) / L(1,1);
  //double z2 = (bl[1] - num) / L(1,1);
  num = mu[2] + L(2,0)*z1 + L(2,1)*z2;
  double alpha3;
  double lambda3;
  if( num < MIN_EDGE_LENGTH + TOL)
      num = MIN_EDGE_LENGTH;
  calculateAlphaLambda(num,L(2,2),eta,alpha3,lambda3);
  bl[2] = gamma(alpha3,1.0 / lambda3,rng); //c++ gamma has different parametrization

  logdensity += alpha1*log(lambda1) + alpha2*log(lambda2) + alpha3*log(lambda3) + (alpha1-1)*log(bl[0])-lambda1*bl[0]+(alpha2-1)*log(bl[1])-lambda2*bl[1]+(alpha3-1)*log(bl[2])-lambda3*bl[2] - lgamma(alpha1) - lgamma(alpha2) - lgamma(alpha3);
  return bl;
}

// eta divides the variance: var/eta
Vector2d multivariateGamma2D(Vector2d mu,Matrix2d vc,mt19937_64& rng, double& logdensity, double eta)
{
  Vector2d bl;
  Matrix2d L( vc.llt().matrixL() );
  // ------------ T1 ------------------
  double alpha1;
  double lambda1;
  double mu0;
  if( mu[0] < MIN_EDGE_LENGTH + TOL)
      mu0 = MIN_EDGE_LENGTH;
  else
      mu0 = mu[0];
  calculateAlphaLambda(mu0,L(0,0),eta,alpha1,lambda1);
  bl[0] = gamma(alpha1,1.0 / lambda1,rng); //c++ gamma has different parametrization

  // ------------ T2 ------------------
  double z1 = (bl[0] - (alpha1/lambda1)) / L(0,0);
  //  double z1 = (bl[0] - mu[0]) / L(0,0);
  //double z1 = (bl[0] - mu0) / L(0,0);
  double num = mu[1] + L(1,0)*z1;
  double alpha2;
  double lambda2;
  if( num < MIN_EDGE_LENGTH + TOL)
      num = MIN_EDGE_LENGTH;
  calculateAlphaLambda(num,L(1,1),eta,alpha2,lambda2);
  bl[1] = gamma(alpha2,1.0 / lambda2,rng); //c++ gamma has different parametrization

  logdensity += alpha1*log(lambda1) + alpha2*log(lambda2) + (alpha1-1)*log(bl[0])-lambda1*bl[0]+(alpha2-1)*log(bl[1])-lambda2*bl[1] - lgamma(alpha1) - lgamma(alpha2);
  return bl;
}

VectorXd multivariateNormal(VectorXd mu,MatrixXd vc,mt19937_64& rng, double& logdensity, double eta)
{
  int r = vc.rows();
  VectorXd v( r );
  MatrixXd L( vc.llt().matrixL() );

  normal_distribution<double> rnorm;
  for ( int i=0; i<r; ++i )
    v( i ) = rnorm( rng );

  logdensity -= v.dot(v);

  return mu + L * v;
}
