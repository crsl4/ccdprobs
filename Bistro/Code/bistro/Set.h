#ifndef SETH
#define SETH

using namespace std;

#include <iostream>
#include <iomanip>

class Set {
  // Simple set class, storing ints in a bit vector. Assumes 32-bit words.
public:

  Set() : els(0) {}

  Set(int len) : els((len+30)/31) {}

  int size() const { return els.size(); }

  void init(int len) { els.resize((len+30)/31); }

  Set& operator= (const Set &a) { 
    els.resize(size());
    for(int i=0;i<size();i++) 
      els[i] = a.els[i]; 
    return *this;
  }

  friend int operator== (const Set &a, const Set &b) { 
    for(int i=0;i<a.size();i++) 
      if(a.els[i] != b.els[i])
	return 0;
    return 1;
  }

  friend int operator!= (const Set &a, const Set &b) { 
    for(int i=0;i<a.size();i++) 
      if(a.els[i] != b.els[i])
	return 1;
    return 0;
  }

  friend int operator<(const Set &a, const Set &b) { return a.cmp(b)<0; }

  friend int operator>(const Set &a, const Set &b) { return a.cmp(b)>0; }

  void clear() { 
    for(int i=0;i<size();i++) 
      els[i] = 0; 
  }

  void singleton(const int n, int len) {
    els.clear();
    els.resize((len+30)/31,0);
    els[(n-1)/31] = (1 << (30 - ((n-1) % 31)));
  }

  unsigned int hash() const {
    unsigned int n = els[0];
    for(int i=1;i<size();i++) {
      n += els[i];
      //n ^= els[i];
    }
    return n;
  }

  void myunion(const Set &a, const Set &b) { 
    els.resize(a.size());
    for(int i=0;i<a.size();i++) 
      els[i] = a.els[i] | b.els[i];
  }

  void difference(const Set &a, const Set &b) { 
    els.resize(a.size());
    for(int i=0;i<a.size();i++) 
      els[i] = a.els[i] & ~b.els[i];
  }

  int subset(const Set &a) const {
    for(int i=0;i<size();i++)
      if(els[i] & ~a.els[i])
	return 0;
    return 1;
  }

  void print(ostream &out) const {
    static unsigned int bit = (1 << 30);
    int first = 1;
    unsigned int n;
    int low = 0, high;

    out << "{";
    for(int i=0;i<size();i++) {
      n = els[i];
      for(int j=1;j<32;j++,n<<=1)
	if(n & bit) {
	  int k = j+31*i;
	  if(low==0)
	    low = high = k;
	  else if (high==k-1)
	    high = k;
	  else {
	    printRange(out,low,high,first);
	    low = high = k;
	    first = 0;
	  }
	}
    }
    if(low!=0)
      printRange(out,low,high,first);
    out << "}";
  }

  void printElements(ostream& c) {
    static unsigned int bit = (1 << 30);
    for(int i=0;i<size();i++) {
      unsigned int n = els[i];
      for(int j=1;j<32;j++,n<<=1)
	if(n & bit)
	  c << setw(12) << " " << setw(3) << j+31*i << endl;
    }
  }

  void printHex(ostream& c) {
    for(int i=0;i<size();i++)
      c << hex << els[i] << " ";
    c << dec << endl;
  }

private:
  vector<unsigned int> els;

  int cmp(const Set &a) const {
    int c;
    for(int i=0;i<size();i++) 
      if((c=els[i]-a.els[i])!=0)
	return c;
    return 0;
  }

  void printRange(ostream &out, int low, int high, int first) const
  {
    if(!first)
      out << ",";
    out << low;
    if(high-low>0)
      out << "-" << high;
  }
};

#endif
