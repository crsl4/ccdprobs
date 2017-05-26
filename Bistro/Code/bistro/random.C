#include <random>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "Eigen/Core"
#include "Eigen/Eigenvalues"

#include "random.h"
#include <cmath> // cbrt, tgamma
#include <boost/math/distributions/normal.hpp>

using namespace std;
using namespace Eigen;

void checkMuV(double mu, double v)
{
  if ( mu < 0 || v < 0 )
  {
    cerr << "Error with either mu or v negative: " << mu << ", " << v << endl;
    exit(1);
  }
  else if ( mu < TOL || v < TOL )
    cerr << "Warning with either mu or v smaller than TOL: " << mu << ", " << v << endl;
}

// In original code, v was the sd
// Now, v is the variance
void calculateAlphaLambda(double mu, double v, double eta, double& alpha, double& lambda)
{
  checkMuV(mu,v);
  alpha = mu*mu / v;
  lambda = mu / v;
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
  //cout << "almv, " << alpha << ", " << lambda << ", " << mu << ", " << v << endl;
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
  calculateAlphaLambda(mu0,L(0,0)*L(0,0),eta,alpha1,lambda1);
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
  calculateAlphaLambda(num,L(1,1)*L(1,1),eta,alpha2,lambda2);
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
  calculateAlphaLambda(num,L(2,2)*L(2,2),eta,alpha3,lambda3);
  bl[2] = gamma(alpha3,1.0 / lambda3,rng); //c++ gamma has different parametrization
  logdensity += alpha1*log(lambda1) + alpha2*log(lambda2) + alpha3*log(lambda3) + (alpha1-1)*log(bl[0])-lambda1*bl[0]+(alpha2-1)*log(bl[1])-lambda2*bl[1]+(alpha3-1)*log(bl[2])-lambda3*bl[2] - lgamma(alpha1) - lgamma(alpha2) - lgamma(alpha3);
  if( logdensity > 1000000 )
  {
    cerr << "found inf logBL" << endl;
    cout << "found inf logBL" << endl;
    cout << alpha1 << "," << alpha2 << "," << alpha3 << "," << lambda1 << "," << lambda2 << "," << lambda3 << endl;
    exit(1);
  }
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
  calculateAlphaLambda(mu0,L(0,0)*L(0,0),eta,alpha1,lambda1);
  bl[0] = gamma(alpha1,1.0 / lambda1,rng); //c++ gamma has different parametrization

  // ------------ T2 ------------------
  double z1 = (bl[0] - (alpha1/lambda1)) / L(0,0);
  //double z1 = (bl[0] - mu[0]) / L(0,0);
  //double z1 = (bl[0] - mu0) / L(0,0);
  double num = mu[1] + L(1,0)*z1;
  double alpha2;
  double lambda2;
  if( num < MIN_EDGE_LENGTH + TOL)
    num = MIN_EDGE_LENGTH;
  calculateAlphaLambda(num,L(1,1)*L(1,1),eta,alpha2,lambda2);
  bl[1] = gamma(alpha2,1.0 / lambda2,rng); //c++ gamma has different parametrization

  logdensity += alpha1*log(lambda1) + alpha2*log(lambda2) + (alpha1-1)*log(bl[0])-lambda1*bl[0]+(alpha2-1)*log(bl[1])-lambda2*bl[1] - lgamma(alpha1) - lgamma(alpha2);
  if( logdensity > 1000000 )
  {
    cerr << "found inf logBL" << endl;
    cout << "found inf logBL" << endl;
    cout << alpha1 << "," << alpha2 << "," << lambda1 << "," << lambda2 << endl;
    exit(1);
  }
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

void calculateMuSigma(double t0, double t1, double t2, double y0, double y1, double y2, double& mu, double& sigma2)
{
  if(VERBOSE)
    cerr << "t0,t1,t2,y0,y1,y2: " << t0 << "," << t1 << "," << t2 << "," << y0 << "," << y1 << "," << y2 << endl;
  double u01 = y0-y1;
  double v10 = t1*t1-t0*t0;
  double w10 = t1-t0;
  double u12 = y1-y2;
  double v21 = t2*t2-t1*t1;
  double w21 = t2-t1;
  mu = (v10*u12-v21*u01)/(2*(w10*u12-w21*u01));
  sigma2 = w10*w21*(t2-t0)/(2*(u12*w10-u01*w21));
  // if ( sigma2 < 0 )
  // {
  //   cerr << "found sigma2 negative in calculate mu and sigma from three points" << endl;
  //   cerr << mu << "," << sigma2 << endl;
  //   exit(1);
  // }
  // sigma = sqrt(sigma2);
}

double normalTailArea(double z)
{
  return 0.5 * erfc(z/sqrt(2.0));
}

double randomHalfNormal(double mu, double sigma, double x0,double& logdensity,mt19937_64& rng)
{
  double z0 = (x0-mu)/sigma;
  uniform_real_distribution<> runif(0, 1);
  double tailArea = normalTailArea(z0);
  double a = 1 - (1-runif(rng))*tailArea;
  if(VERBOSE)
  {
    cerr << "in random half normal, with mu: " << mu << ", sigma " << sigma << endl;
    cerr << "z0: " << z0 << endl;
    cerr << "erfc(z0): " << erfc(z0) << endl;
    cerr << "a= " << a << endl;
    cerr << "log tail= " << log(tailArea) << endl;
  }
  boost::math::normal rnorm(0.0, 1.0);
  double q = quantile(rnorm, a);
  logdensity += -log(tailArea) - log(sigma) - (0.5)*log(2*M_PI) - (0.5)*q*q;
  if( logdensity > 1000000 )
  {
    cerr << "found inf logBL" << endl;
    cout << "found inf logBL" << endl;
    cout << mu << "," << sigma << "," << tailArea << endl;
    exit(1);
  }
  return mu + sigma*q;
}

// half normal or gamma depending on the length
// using exponential when z0 > 8
double halfNormalGamma(Edge* e, Alignment& alignment, QMatrix& qmatrix, double& logdensity,mt19937_64& rng)
{
  double t = e->getLength();
  if( t < MIN_EDGE_LENGTH + TOL ) //half normal
  {
    double x1 = 0.001;
    double x2 = 0.002;
    double x3 = 0.003;
    double y1,y2,y3;
    double dlogl,ddlogl;
    e->calculate(x1,alignment,qmatrix,y1,dlogl,ddlogl);
    e->calculate(x2,alignment,qmatrix,y2,dlogl,ddlogl);
    e->calculate(x3,alignment,qmatrix,y3,dlogl,ddlogl);
    double mu;
    double sigma2;
//    cerr << "t = " << t << endl;
//    e->printLikMinLength(cerr,alignment,qmatrix);
    calculateMuSigma(x1,x2,x3,y1,y2,y3,mu,sigma2);
    if(VERBOSE)
      cerr << "mu: " << mu << ",. sigma2: " << sigma2 << endl;
    double sigma;
    double z0;
    double s;
    if ( sigma2 > 0 )
    {
      sigma = sqrt(sigma2);
      z0 = -mu/sigma;
    }

    if ( sigma2 > 0 && z0 < 7)// && mu < 0 )
    {
      s = randomHalfNormal(mu,sigma,0.0,logdensity,rng);
      e->addSampler(0,1,0);
      if(VERBOSE)
	cerr << "found half normal case, t: " << t << ", mu: " << mu << ", sigma: " << sigma << ", sampled: " << s << endl;
      if ( mu > 0 )
      {
	cout << "t = " << t << ", mu = " << mu << ", sigma2 = " << sigma2 << endl;
	e->printLikMinLength(cout,alignment,qmatrix);
//	exit(1);
      }
    }
    else
    {
      double lambda = (y3-y1) / (x3-x1);
      lambda = (-1)*lambda;
      s = gamma(1.0, 1.0 / lambda ,rng);
      e->addSampler(0,0,1);
      if(VERBOSE)
	cerr << "found exponential case, t: " << t << ", lambda= " << lambda << ", sampled: " << s << endl;
      logdensity += log(lambda) - lambda*s;
      if( logdensity > 1000000 )
      {
	cerr << "found inf logBL" << endl;
	cout << "found inf logBL" << endl;
	cout << lambda << endl;
	exit(1);
      }
    }
    return s;
  }
  else //gamma
  {
    double logl,dlogl,ddlogl;
    e->calculate(t,alignment,qmatrix,logl,dlogl,ddlogl);
    double lambda = -1 * t * ddlogl;
    double alpha = t*lambda;
    gamma_distribution<double> rgamma(alpha,1.0 / lambda);
    t = rgamma(rng);
    e->addSampler(1,0,0);
    if(VERBOSE)
      cerr << "found gamma case, t: " << t << ", alpha: " << alpha << ", lambda: " << lambda << ", sampled: " << t << endl;
    logdensity += alpha * log(lambda) - lgamma(alpha) + (alpha-1)*log(t) - lambda*t;
    return t;
  }
}
