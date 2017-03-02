#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <string>

using namespace std;

int main()
{
  string tree = "(1:0.0577584,(2:0.0572218,((3:0.0677952,((4:0.0375919,5:0.0385211):0.0105002,6:0.0574576):0.0098151):0.0140263,(7:0.0000100,(8:0.0000100,(9:0.0000100,(10:0.0000100,(11:nan,12:nan):nan):0.0000100):0.0000100):0.0000100):0.0000100):0.0096315):0.0000000);";
  if ( tree.find("nan") != std::string::npos )
    cerr << "Found nan" << endl;
  return 0;
}
