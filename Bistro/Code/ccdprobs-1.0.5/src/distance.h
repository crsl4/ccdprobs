#ifndef DISTANCE_H_
#define DISTANCE_H_

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <cctype>
#include <ctime>
#include <cmath>

using namespace std;

class DistanceMatrix
{
 private:
  int numTaxa;
  vector<vector<double> > dist;
 public:
  DistanceMatrix() {};
  void reset(int n) {
    numTaxa = n;
    dist.resize(n);
    for ( int i=0; i<n; i++ )
      dist[i].resize(n,0);
  }
  DistanceMatrix(int n) {
    reset(n);
  }
  double get(int i,int j) const { return dist[i][j]; }
  void set(int i, int j,double x) { dist[i][j] = x; }
  void empty() {
    for ( int i = 0; i<numTaxa; i++ )
      for ( int j = 0; j< numTaxa; j++ )
	dist[i][j] = 0;
  }
  int size() const { return dist.size(); }
  void print(ostream& f) {
      f.setf(ios::fixed,ios::floatfield);
      f << endl;
      for ( int i=0; i<numTaxa; i++ ) {
	for ( int j=0; j<numTaxa; j++ )
	  f << setprecision(6) << setw(9) << dist[i][j];
	f << endl;
      }
      f << endl;
  }
  void divide(double x) {
    for ( int i=0; i<numTaxa; ++i )
      for ( int j=0; j<numTaxa; ++j )
	dist[i][j] /= x;
  }
  void sqrt() {
    for ( int i=0; i<numTaxa; ++i )
      for ( int j=0; j<numTaxa; ++j )
	dist[i][j] = std::sqrt(dist[i][j]);
  }
};

#endif
