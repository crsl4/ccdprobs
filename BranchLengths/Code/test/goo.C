#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <random>

#include "Eigen/Core"
#include "Eigen/Eigenvalues"
#include "Eigen/LU"

using namespace std;
using namespace Eigen;

int main()
{
  Matrix4d qmat = MatrixXd::Random(4,4);

  cout << qmat << endl << endl;

  cout << qmat.inverse() << endl;
  
}
