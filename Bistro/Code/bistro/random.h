#ifndef __RANDOM_H
#define __RANDOM_H

#include <random>
#include <iostream>

#include "Eigen/Core"
#include "Eigen/Eigenvalues"

#include "tree.h"

using namespace std;
using namespace Eigen;

double gamma(double, double, mt19937_64&);
Vector3d multivariateGamma3D(Vector3d,Matrix3d,mt19937_64&, double&, double);
Vector2d multivariateGamma2D(Vector2d,Matrix2d,mt19937_64&, double&, double);

#endif
