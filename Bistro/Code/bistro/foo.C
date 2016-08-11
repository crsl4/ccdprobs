#include <iostream>
#include <iomanip>
#include <boost/dynamic_bitset.hpp>
#include <string>
#include <sstream>

using namespace std;
using namespace boost;

class Clade {
private:
  dynamic_bitset<unsigned char> clade; // one bit per taxa 0/1
public:
  Clade() {}
  Clade(int n) { clade.resize(n); }
  ~Clade() { clade.clear(); }
  void add(int x) { clade[x-1] = 1; }
  void print(ostream&) const ;
};

void Clade::print(ostream& f) const {
    
  string comma = "";

  dynamic_bitset<unsigned char>::size_type first = clade.find_first();
  dynamic_bitset<unsigned char>::size_type last = first;
  dynamic_bitset<unsigned char>::size_type i = first;

  f << '{';
  while ( i != dynamic_bitset<unsigned char>::npos ) {
    i = clade.find_next(i);
    if ( i > last+1 ) { // print the previous range and set first to i
      f << comma << first+1;
      comma = ",";
      if ( last > first )
	f << '-' << last+1;
      first = i;
    }
    last = i;
  }
  f << '}';
}

int main()
{
  Clade clade(10);
  clade.add(4);
  clade.add(5);
  clade.add(8);
  
  clade.print(cout);
  cout << endl;
  
  return 0;
}
