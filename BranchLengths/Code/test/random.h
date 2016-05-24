#ifndef __RANDOM_H
#define __RANDOM_H

#define TOL 1.0e-4

#include <random>
#include <iostream>

#include "Eigen/Core"
#include "Eigen/Eigenvalues"

using namespace std;
using namespace Eigen;

VectorXd multivariateNormal(VectorXd,MatrixXd,mt19937_64&);
double normal(double, double, mt19937_64&);
double gamma(double, double, mt19937_64&);
double beta(double, double, mt19937_64&);
Vector3d multivariateGamma3D(Vector3d,Matrix3d,mt19937_64&, double&);
Vector3d multivariateGamma2D(Vector2d,Matrix2d,double,mt19937_64&, double&);
Vector3d multivariateGamma1D(double,double,double,double,mt19937_64&, double&);
void pathologicalBetaPar(double,double&,double&);

#endif
