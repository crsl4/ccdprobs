#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <random>

#include "Eigen/Core"
#include "Eigen/Eigenvalues"

using namespace std;
using namespace Eigen;

int main()
{
  random_device rd;
  unsigned int seed = rd();

  mt19937_64 rng(seed);
  uniform_real_distribution<double> runif(0.0,1.0);

  Vector4d p;
//  vector<double> p(4);
//  double sum=0;
  for ( int i=0; i<4; ++i )
  {
//    p[i] = runif(rng);
//    sum += p[i];
    p(i) = runif(rng);
  }
//  for ( int i=0; i<4;++i )
//    p[i] /= sum;
  p = p/p.sum();
  
//  vector<double> s(6);
//  sum = 0;
//  for ( int i=0; i<6; ++i )
//  {
//    s[i] = runif(rng);
//    sum += s[i];
//  }
//  for ( int i=0; i<6; ++i )
//    s[i] /= sum;
  VectorXd s(6);
  for ( int i=0; i<6; ++i )
    s(i) = runif(rng);
  
//  cout << "p:";
//  for ( int i=0; i<4; ++i)
//    cout << " " << p[i];
//  cout << endl;
//
//  cout << "s:";
//  for ( int i=0; i<6; ++i)
//    cout << " " << s[i];
//  cout << endl;

  cout << p.transpose() << endl << endl;
  cout << s.transpose() << endl << endl;
  
  map<pair<int,int>,int> ijMap;
  ijMap[ pair<int,int> (0,1) ] = 0;
  ijMap[ pair<int,int> (0,2) ] = 1;
  ijMap[ pair<int,int> (0,3) ] = 2;
  ijMap[ pair<int,int> (1,2) ] = 3;
  ijMap[ pair<int,int> (1,3) ] = 4;
  ijMap[ pair<int,int> (2,3) ] = 5;

  Matrix4d qmat;
  
  for ( int i=0; i<4; ++i )
  {
    qmat(i,i) = 0;
    for (int j=0; j<4; ++j )
    {
      if ( i < j )
      {
        qmat(i,j) = s( ijMap[ pair<int,int> (i,j) ] ) / (2.0 * p(i));
        qmat(i,i) -= qmat(i,j);
      }
      else if ( i > j ) {
        qmat(i,j) = s( ijMap[ pair<int,int> (j,i) ] ) / (2.0 * p(i));
        qmat(i,i) -= qmat(i,j);
      }
    }
  }

  cout << setprecision(9) << qmat << endl << endl;

  Vector4d psqrt;
  Vector4d psqrtinv;
  for ( int i=0; i<4; ++i )
  {
    psqrt(i) = sqrt( p(i) );
    psqrtinv(i) = 1 / psqrt(i);
  }
  
  DiagonalMatrix<double,4> pSqrtDiag = psqrt.asDiagonal();
  DiagonalMatrix<double,4> pSqrtDiagInv = psqrtinv.asDiagonal();
  
  Matrix4d symmat = pSqrtDiag * qmat * pSqrtDiagInv;
  
  cout << symmat << endl << endl;

  SelfAdjointEigenSolver<Matrix4d> es(symmat);
  cout << es.eigenvalues() << endl << endl;
  cout << es.eigenvectors() << endl << endl;
  
  
}
