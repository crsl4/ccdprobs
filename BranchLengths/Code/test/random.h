#ifndef __RANDOM_H
#define __RANDOM_H

#include <random>

#include "Eigen/Core"
#include "Eigen/Eigenvalues"

using namespace std;
using namespace Eigen;

VectorXd multivariateNormal(VectorXd,MatrixXd,mt19937_64&);
double normal(double, double, mt19937_64&);

#endif
