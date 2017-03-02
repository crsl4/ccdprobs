#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <string>
#include <sstream>

using namespace std;

int main()
{
  stringstream s;
  s << "thisIsAString" << endl;
  string foo;
  s >> foo;
  s << foo;

  cout << foo << endl;

  foo.clear();

  cout << foo << endl;

  s >> foo;

  cout << foo << endl;
  
  return 0;
}
