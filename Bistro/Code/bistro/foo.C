#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <string>
#include <sstream>

#include <boost/math/distributions/normal.hpp>

using namespace std;

int main()
{
  for ( int i=0; i< 30; ++i )
    cout << setw(2) << i << " " << erfc(double(i)) << endl;
  return 0;
}
