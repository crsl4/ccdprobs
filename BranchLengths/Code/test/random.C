#include <random>
#include <iostream>
#include <iomanip>

#include "Eigen/Core"
#include "Eigen/Eigenvalues"

#include "random.h"
#include <cmath> // cbrt

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
VectorXd multivariateNormal(VectorXd mu,MatrixXd vc,mt19937_64& rng, double& logdensity)
{
  int r = vc.rows();
  VectorXd v( r );
  MatrixXd L( vc.llt().matrixL() );

  normal_distribution<double> rnorm;
  double rss = 0;
  for ( int i=0; i<r; ++i )
    {
      v( i ) = rnorm( rng );
      rss += v(i)*v(i);
    }

  logdensity = -0.5*rss;

  return mu + L * v;
}

double normal(double mu, double vc, mt19937_64& rng)
{
  normal_distribution<double> rnorm(mu,sqrt(vc));
  return rnorm( rng );
}

double gamma(double alpha, double b, mt19937_64& rng)
{
  gamma_distribution<double> gam1(alpha,b);
  return gam1(rng);
}


double beta(double alpha,double b,mt19937_64& rng)
{
  if(alpha<0 || b< 0)
    {
      cerr << "Trying to generate beta rv with negative arguments" << endl;
      exit(1);
    }
  gamma_distribution<double> rgammaA(alpha,1);
  gamma_distribution<double> rgammaB(b,1);
  double x1 = rgammaA(rng);
  double x2 = rgammaB(rng);
  return( x1/(x1+x2) );
 }

Vector3d multivariateGamma3D(Vector3d mu,Matrix3d vc,mt19937_64& rng, double& logdensity, bool verbose)
{
  Vector3d bl;
  Matrix3d L( vc.llt().matrixL() );
  double alpha1;
  double lambda1;
  if(mu[0] < 0)
    {
      cerr << "Error with mu1<0 in multivariateGamma3D" << endl;
      exit(1);
    }
  if( mu[0] < TOL)
    {
      if(verbose)
	cerr << "Case mu1==0 in multivariateGamma3D" << endl;
      //exit(1);
      alpha1 = 1;
      lambda1 = 1/L(0,0);
    }
  if( mu[0] > TOL)
    {
      alpha1 = mu[0]*mu[0] / (L(0,0) * L(0,0));
      lambda1 = mu[0] / (L(0,0) * L(0,0));
    }
  // ------------ T1 ------------------
  bl[0] = gamma(alpha1,1/lambda1,rng); //c++ gamma has different parametrization
  if(verbose)
    cout << "T1: " << bl[0] << endl;
  // ------------ T2 ------------------
  double z1 = (bl[0] - mu[0]) / L(0,0);
  double num = mu[1] + L(1,0)*z1;
  double alpha2;
  double lambda2;
  if( num < 0)
    {
      cerr << "mu2 + L21*z1 in multivariateGamma3D is negative" << endl;
      exit(1);
    }
  if( num < TOL)
    {
      if(verbose)
	cerr << "mu2 + L21*z1 in multivariateGamma3D is zero" << endl;
      //exit(1);
      alpha2 = 1;
      lambda2 = 1/L(1,1);
    }
  if( num > TOL)
    {
      alpha2 = (num * num) / (L(1,1) * L(1,1));
      lambda2 = num / (L(1,1) * L(1,1));
    }
  bl[1] = gamma(alpha2,1/lambda2,rng); //c++ gamma has different parametrization
  if(verbose)
    cout << "T2: " << bl[1] << endl;
  // ------------ T3 ------------------
  double z2 = (bl[1] - mu[1] - L(1,0)*z1) / L(1,1);
  num = mu[2] + L(2,0)*z1 + L(2,1)*z2;
  double alpha3;
  double lambda3;
  if( num < 0)
    {
      cerr << "mu3+L31*z1+L32*z2 in multivariateGamma3D is negative" << endl;
      exit(1);
    }
  if( num < TOL)
    {
      if(verbose)
	cerr << "mu3+L31*z1+L32*z2 in multivariateGamma3D is zero" << endl;
      //exit(1);
      alpha3 = 1;
      lambda3 = 1/L(1,1);
    }
  if( num > TOL)
    {
      alpha3 = (num * num) / (L(2,2) * L(2,2));
      lambda3 = num / (L(2,2) * L(2,2));
    }
  bl[2] = gamma(alpha3,1/lambda3,rng); //c++ gamma has different parametrization
  if(verbose)
    cout << "T3: " << bl[2] << endl;
  logdensity += (alpha1-1)*log(bl[0])-lambda1*bl[0]+(alpha2-1)*log(bl[1])-lambda2*bl[1]+(alpha3-1)*log(bl[2])-lambda3*bl[2];
  return bl;
 }


  // t1,t2,sum-t1; so mu, vc are for t1,t2
Vector3d multivariateGamma2D(Vector2d mu,Matrix2d vc,double sum, mt19937_64& rng, double& logdensity, bool verbose)
 {
   if(verbose)
     cout << "multivariateGamma2D with mu: " << mu.transpose() << endl;
   double alpha1;
   double beta1;
   Vector3d bl;
   Matrix2d L( vc.llt().matrixL() );
   if( mu[0] < 0)
     {
       cerr << "Error with mu1<0 in multivariateGamma2D" << endl;
       exit(1);
     }
   if( mu[0] < TOL)
     {
       if(verbose)
	 cerr << "Case mu1==0 in multivariateGamma2D" << endl;
       //exit(1);
       pathologicalBetaPar(L(0,0),alpha1,beta1);
     }
   if( mu[0] > sum - TOL)
     {
       if(verbose)
	 cerr << "Case mu1==sum in multivariateGamma2D" << endl;
       //exit(1);
       pathologicalBetaPar(L(0,0),beta1,alpha1);
     }
   if( mu[0] > TOL && mu[0] <= sum - TOL)
     {
       if( mu[0]*(sum-mu[0]) <= (L(0,0)*L(0,0)))
	 {
	   cerr << "mu(s-mu)<sigma2 in multivariateGamma2D" << endl;
	   exit(1);
	 }
       double part1 = (mu[0] * mu[0] * (sum-mu[0])) / (sum * L(0,0) * L(0,0));
       double part2 = (mu[0] * (sum-mu[0])*(sum-mu[0])) / (sum * L(0,0) * L(0,0));
       alpha1 = part1 - (mu[0] / sum);
       beta1 = part2 - ((sum - mu[0]) / sum);
     }
   if(alpha1 < TOL || beta1 < TOL )
     {
       cerr << "Alpha or beta negative or too close to zero in multivariate2D: " << alpha1<< ", " << beta1 << endl;
       exit(1);
     }
   //------------------ B ------------------------
   double b = beta(alpha1,beta1,rng);
   bl[0] = sum * b; //t1
   bl[2] = sum *(1-b); //t3
   // ------------ t2 ------------------
   double z1 = (bl[0] - mu[0]) / L(0,0);
   double num = mu[1] + L(1,0)*z1;
   double alpha2;
   double lambda2;
   if( num < 0)
     {
       cerr << "mu2+L21*z1 in multivariateGamma2D is negative" << endl;
       exit(1);
     }
   if( num < TOL)
     {
       if(verbose)
	 cerr << "mu2+L21*z1 in multivariateGamma2D is zero" << endl;
       //exit(1);
       alpha2 = 1;
       lambda2 = 1/L(1,1);
     }
   if( num > TOL)
     {
       alpha2 = (num * num) / (L(1,1) * L(1,1));
       lambda2 = num / (L(1,1) * L(1,1));
     }
   bl[1] = gamma(alpha2,1/lambda2,rng); //c++ gamma has different parametrization
   logdensity += (alpha1-1)*log(bl[0])+(beta1-1)*log(sum-bl[0])+(alpha2-1)*log(bl[1])-lambda2*bl[1];
   return bl;
 }


  // t1,sum1-t1,sum2-t1
Vector3d multivariateGamma1D(double mu,double var,double sum1, double sum2, mt19937_64& rng, double& logdensity, bool verbose)
 {
   if(verbose)
     cout << "multivariateGamma1D with mu: " << mu << endl;
   double a;
   double b;
   double s = min(sum1,sum2);
   Vector3d bl;
   if( mu < 0)
     {
       cerr << "Error with mu<0 in multivariateGamma1D" << endl;
       exit(1);
     }
   if( mu < TOL)
     {
       if(verbose)
	 cerr << "Case mu==0 in multivariateGamma1D" << endl;
       //exit(1);
       pathologicalBetaPar(var,a,b);
     }
   if( mu > s - TOL)
     {
       if(verbose)
	 cerr << "Case mu==sum in multivariateGamma1D" << endl;
       //exit(1);
       pathologicalBetaPar(var,b,a);
     }
   if( mu > TOL && mu <= s - TOL)
     {
       if( mu*(s-mu) <= var)
	 {
	   cerr << "mu(s-mu)<sigma2 in multivariateGamma1D" << endl;
	   exit(1);
	 }
       double part1 = (mu * mu * (s-mu)) / (s * var);
       double part2 = (mu * (s-mu)*(s-mu)) / (s * var);
       a =  part1 - mu / s;
       b = part2 - (s - mu) / s;
       if(verbose)
	 cout << "a,b in multivariateGamma1D" << a << " , " << b << endl;
     }
   if(a < TOL || b < TOL )
     {
       cerr << "Alpha or beta negative or too close to zero in multivariate1D: " << a<< ", " << b << endl;
       exit(1);
     }
   if(verbose)
     {
       cout << "a,b in multivariateGamma1D" << a << " , " << b << endl;
       cout << "1D mean: " << mu << ", variance: " << var << endl;
       cout << "sum1, sum2 " << sum1 << ", " << sum2 << endl;
       cout << "min sum: " << s << endl;
     }
   double rbeta = beta(a,b,rng);
   if(verbose)
     cout << "rbeta: " << rbeta << " with a,b " << a << " , " << b << endl;
   bl[0] = rbeta * s;
   bl[1] = sum1-bl[0];
   bl[2] = sum2-bl[0];
   logdensity += (a-1)*log(bl[0])+(b-1)*log(s-bl[0]);
   if(verbose)
     cout << "after mleDistance1D, bl: " << bl.transpose() << endl;
   return bl;
 }

// set alpha=1, and beta=long expression
void pathologicalBetaPar(double v,double& alpha,double& beta)
{
  alpha = 1;
  double ct = v*v*v*v*v + 11*v*v*v*v - v*v*v;
  if( ct < TOL)
    ct = 0;
  double cb = -v*v*v - 18*v*v +3*sqrt(3)*sqrt(ct);
  if( cb < 0)
    {
      //cerr << "made inside of cbrt positive in pathologicalBetaPar: " << cb << endl;
      cb = fabs(cb);
    }
  //cerr << "thing inside sqrt: " << ct << endl;
  //cerr << "thing inside cbrt: " << cb  << endl;
  double den = cbrt(cb);
  if( den < TOL)
    {
      cerr << "denominator for pathologicalBetaPar is too close to zero: " << den << endl;
      exit(1);
    }
  beta = (v+3)/(3*den) + den/(3*v) - (4/3.0); // made positive always
}


