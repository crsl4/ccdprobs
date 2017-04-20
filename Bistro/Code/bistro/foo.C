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
  {
    double tail = 0.5 * erfc(i/sqrt(2.0));
    boost::math::normal rnorm(0.0, 1.0);
    double q = quantile(rnorm, 1-tail);
    cout << setw(2) << i << " " << 1-tail << " " << q << endl;
  }
  return 0;
}
