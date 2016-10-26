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
  double alpha1;
  double lambda1;
  double mu0;
  // if(mu[0] < 0)
  //   {
  //     cerr << "Error with mu1<0 in multivariateGamma3D" << endl;
  //     throw 20;
  //   }
  if( mu[0] < MIN_EDGE_LENGTH + TOL)
    {
//      cout << "Case mu1==0 (or negative )in multivariateGamma3D (will use exponential): " << mu[0] << endl;
      // alpha1 = 1.0;
      // lambda1 = 1.0 / (MIN_EDGE_LENGTH);
      mu0 = MIN_EDGE_LENGTH;
    }
  else
    {
      // alpha1 = eta*mu[0]*mu[0] / (L(0,0) * L(0,0));
      // lambda1 = eta*mu[0] / (L(0,0) * L(0,0));
      // if(alpha1 < 1) // avoid alpha1<1
      // 	{
      // 	  lambda1 = lambda1 / alpha1;
      // 	  alpha1 = 1.0;
      // 	}
      mu0 = mu[0];
    }
  alpha1 = eta*mu0*mu0 / (L(0,0) * L(0,0));
  lambda1 = eta*mu0 / (L(0,0) * L(0,0));

  if(alpha1<1)
    {
      alpha1 = alpha1 + 1.0;
      lambda1 = lambda1 + 1.0;
    }
  // ------------ T1 ------------------
//  cout << "Random gamma: alpha = " << alpha1 << ", lambda = " << lambda1 << ", mu = " << alpha1 / lambda1 << ", sigma = " << sqrt(alpha1) / lambda1 << endl;

  bl[0] = gamma(alpha1,1.0 / lambda1,rng); //c++ gamma has different parametrization
//  cout << "T1: " << bl[0] << endl;
  // ------------ T2 ------------------
  //  double z1 = (bl[0] - mu[0]) / L(0,0);
  double z1 = (bl[0] - mu0) / L(0,0);
  double num = mu[1] + L(1,0)*z1;
  double alpha2;
  double lambda2;
  // if( num < 0)
  //   {
  //     cerr << "mu2 + L21*z1 in multivariateGamma3D is negative, mu2,L21,z1: " << mu[1] << ","<< L(1,0) << "," << z1 << "," << num << endl;
  //     throw 20;
  //   }
  if( num < MIN_EDGE_LENGTH + TOL)
    {
//      cout << "mu2 + L21*z1 in multivariateGamma3D is zero or negative (will use exponential): " << num << endl;
      // alpha2 = 1.0;
      // lambda2 = 1.0 / (MIN_EDGE_LENGTH);
      num = MIN_EDGE_LENGTH;
    }
  // else
  //   {
      // alpha2 = eta * (num * num) / (L(1,1) * L(1,1));
      // lambda2 = eta * num / (L(1,1) * L(1,1));
      // if(alpha2 < 1)
      // 	{
      // 	  lambda2 = lambda2 / alpha2;
      // 	  alpha2 = 1.0;
      // 	}
    // }
  alpha2 = eta * (num * num) / (L(1,1) * L(1,1));
  lambda2 = eta * num / (L(1,1) * L(1,1));

//  cout << "Random gamma: alpha = " << alpha2 << ", lambda = " << lambda2 << ", mu = " << alpha2 / lambda2 << ", sigma = " << sqrt(alpha2) / lambda2 << endl;
  if(alpha2<1)
    {
      alpha2 = alpha2 + 1.0;
      lambda2 = lambda2 + 1.0;
    }

  bl[1] = gamma(alpha2,1.0 / lambda2,rng); //c++ gamma has different parametrization
//  cout << "T2: " << bl[1] << endl;
  // ------------ T3 ------------------
  //  double z2 = (bl[1] - mu[1] - L(1,0)*z1) / L(1,1);
  double z2 = (bl[1] - num) / L(1,1);
  num = mu[2] + L(2,0)*z1 + L(2,1)*z2;
  double alpha3;
  double lambda3;
  // if( num < 0)
  //   {
  //     cerr << "mu3+L31*z1+L32*z2 in multivariateGamma3D is negative" << endl;
  //     throw 20;
  //   }
  if( num < MIN_EDGE_LENGTH + TOL)
    {
//      cout << "mu3+L31*z1+L32*z2 in multivariateGamma3D is zero or negative (will use exponential): " << num << endl;
      // alpha3 = 1.0;
      // lambda3 = 1.0 / (MIN_EDGE_LENGTH);
      num = MIN_EDGE_LENGTH;
    }
  // else
  //   {
  //     alpha3 = eta * (num * num) / (L(2,2) * L(2,2));
  //     lambda3 = eta * num / (L(2,2) * L(2,2));
  //     if(alpha3 < 1)
  // 	{
  // 	  lambda3 = lambda3 / alpha3;
  // 	  alpha3 = 1.0;
  // 	}
  //   }
  alpha3 = eta * (num * num) / (L(2,2) * L(2,2));
  lambda3 = eta * num / (L(2,2) * L(2,2));
//  cout << "Random gamma: alpha = " << alpha3 << ", lambda = " << lambda3 << ", mu = " << alpha3 / lambda3 << ", sigma = " << sqrt(alpha3) / lambda3 << endl;
  if(alpha3<1)
    {
      alpha3 = alpha3 + 1.0;
      lambda3 = lambda3 + 1.0;
    }

  bl[2] = gamma(alpha3,1.0 / lambda3,rng); //c++ gamma has different parametrization
//  cout << "T3: " << bl[2] << endl;
  if(logdensity > 1000000000)
    {
      cerr << "before computing new: " << logdensity << endl;
      exit(1);
    }
  logdensity += alpha1*log(lambda1) + alpha2*log(lambda2) + alpha3*log(lambda3) + (alpha1-1)*log(bl[0])-lambda1*bl[0]+(alpha2-1)*log(bl[1])-lambda2*bl[1]+(alpha3-1)*log(bl[2])-lambda3*bl[2] - lgamma(alpha1) - lgamma(alpha2) - lgamma(alpha3);
  if(logdensity > 1000000000)
    {
      cerr << alpha1 << endl;
      cerr << alpha2 << endl;
      cerr << alpha3 << endl;
      cerr << lambda1 << endl;
      cerr << lambda2 << endl;
      cerr << lambda3 << endl;
      cerr << bl.transpose() << endl;
      cerr << logdensity << endl;
      exit(1);
    }
  return bl;
}

// eta divides the variance: var/eta
Vector2d multivariateGamma2D(Vector2d mu,Matrix2d vc,mt19937_64& rng, double& logdensity, double eta)
{
  Vector2d bl;
  Matrix2d L( vc.llt().matrixL() );
  double alpha1;
  double lambda1;
  double mu0;
  // if(mu[0] < 0)
  //   {
  //     cerr << "Error with mu1<0 in multivariateGamma2D" << endl;
  //     throw 20;
  //   }
  if( mu[0] < MIN_EDGE_LENGTH + TOL)
    {
//      cout << "Case mu1==0 (or negative) in multivariateGamma2D (will use exponential): " << mu[0] << endl;
      // alpha1 = 1.0;
      // lambda1 = 1.0 / (MIN_EDGE_LENGTH);
      mu0 = MIN_EDGE_LENGTH;
    }
  else
    {
      // alpha1 = eta*mu[0]*mu[0] / (L(0,0) * L(0,0));
      // lambda1 = eta*mu[0] / (L(0,0) * L(0,0));
      // if(alpha1 < 1)
      // 	{
      // 	  lambda1 = lambda1 / alpha1;
      // 	  alpha1 = 1.0;
      // 	}
      mu0 = mu[0];
    }
  alpha1 = eta*mu0*mu0 / (L(0,0) * L(0,0));
  lambda1 = eta*mu0 / (L(0,0) * L(0,0));
  // ------------ T1 ------------------
//  cout << "Random gamma: alpha = " << alpha1 << ", lambda = " << lambda1 << ", mu = " << alpha1 / lambda1 << ", sigma = " << sqrt(alpha1) / lambda1 << endl;
  if(alpha1<1)
    {
      alpha1 = alpha1 + 1.0;
      lambda1 = lambda1 + 1.0;
    }

  bl[0] = gamma(alpha1,1.0 / lambda1,rng); //c++ gamma has different parametrization
//  cout << "T1: " << bl[0] << endl;
  // ------------ T2 ------------------
  //  double z1 = (bl[0] - mu[0]) / L(0,0);
  double z1 = (bl[0] - mu0) / L(0,0);
  double num = mu[1] + L(1,0)*z1;
  double alpha2;
  double lambda2;
  // if( num < 0)
  //   {
  //     cerr << "mu2 + L21*z1 in multivariateGamma2D is negative" << endl;
  //     throw 20;
  //   }
  if( num < MIN_EDGE_LENGTH + TOL)
    {
//      cout << "mu2 + L21*z1 in multivariateGamma2D is zero or negative (will use exponential): " << num << endl;
      // alpha2 = 1.0;
      // lambda2 = 1.0 / (MIN_EDGE_LENGTH);
      num = MIN_EDGE_LENGTH;
    }
  // else
  //   {
  //     alpha2 = eta*(num * num) / (L(1,1) * L(1,1));
  //     lambda2 = eta*num / (L(1,1) * L(1,1));
  //     if(alpha2 < 1)
  // 	{
  // 	  lambda2 = lambda2 / alpha2;
  // 	  alpha2 = 1.0;
  // 	}
  //   }

//  cout << "Random gamma: alpha = " << alpha2 << ", lambda = " << lambda2 << ", mu = " << alpha2 / lambda2 << ", sigma = " << sqrt(alpha2) / lambda2 << endl;
  alpha2 = eta*(num * num) / (L(1,1) * L(1,1));
  lambda2 = eta*num / (L(1,1) * L(1,1));
  if(alpha2<1)
    {
      alpha2 = alpha2 + 1.0;
      lambda2 = lambda2 + 1.0;
    }

  bl[1] = gamma(alpha2,1.0 / lambda2,rng); //c++ gamma has different parametrization
//  cout << "T2: " << bl[1] << endl;
  logdensity += alpha1*log(lambda1) + alpha2*log(lambda2) + (alpha1-1)*log(bl[0])-lambda1*bl[0]+(alpha2-1)*log(bl[1])-lambda2*bl[1] - lgamma(alpha1) - lgamma(alpha2);
  if(logdensity > 1000000000)
    {
      cerr << "Random gamma: alpha = " << alpha2 << ", lambda = " << lambda2 << ", mu = " << alpha2 / lambda2 << ", sigma = " << sqrt(alpha2) / lambda2 << endl;
      cerr << alpha1 << endl;
      cerr << lambda1 << endl;
      cerr << alpha2 << endl;
      cerr << lambda2 << endl;
      cerr << bl.transpose() << endl;
      cerr << logdensity << endl;
      cerr << num << endl;
      exit(1);
    }

  return bl;
}
