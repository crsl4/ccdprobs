#include <iostream>
#include <fstream>
#include <random>

using namespace std;

int main() {
  mt19937 rng(12345);
  
//  mt19937 rng(12345);
//  mt19937_64 rng(12345);

  normal_distribution<double> rnorm;

  for ( int i=0; i<100; ++i )
    double foo = rnorm(rng);

  ofstream fout("save.txt");
  fout << rng;
  fout.close();

  cout << "Sample 1:";
  for ( int i=0; i<5; ++i )
    cout << " " << rnorm(rng);
  cout << endl;
  
// reset the rng    
  ifstream fin("save.txt");
  if ( fin.fail() || fin.bad() )
    cerr << "fin failed." << endl;
  fin >> rng;
  fin.close();
  cout << "Sample 2:";
  for ( int i=0; i<5; ++i )
    cout << " " << rnorm(rng);
  cout << endl;

// new rng
  mt19937_64 rng2;
// new stream
  ifstream fin2;
//  mt19937_64 rng2;
  fin2.open("save.txt");
  if ( fin2.fail() || fin2.bad() )
    cerr << "fin2 failed." << endl;
  
  fin2 >> rng2;
  cout << "Sample 3:";
  for ( int i=0; i<5; ++i )
    cout << " " << rnorm(rng2);
  cout << endl;
  

  return 0;
}
