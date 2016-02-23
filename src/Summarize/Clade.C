#include "Clade.h"

string cladeNameInt(int n, string& str) { 
  if(n>26)
    cladeNameInt(n/26,str);
  str.push_back(char('A'+n-1));
  return str;
}

string cladeName(int name) { 
  string str="";
  return cladeNameInt(name,str);
}

string topologyName(int name, int num) {
  string buf;
  OSTRSTREAM f(buf);
  f << num;
  return cladeName(name) + f.str();
}

