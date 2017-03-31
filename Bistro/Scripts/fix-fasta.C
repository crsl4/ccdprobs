#include <iostream>
#include <iomanip>
#include <map>
#include <string>
#include <sstream>
#include <fstream>

using namespace std;

string trim(const string& str)
{
    size_t first = str.find_first_not_of(' ');
    if (string::npos == first)
    {
        return str;
    }
    size_t last = str.find_last_not_of(' ');
    return str.substr(first, (last - first + 1));
}

int readFastaSequence(istream& f,int& name,string& sequence)
{
  char c;
  f >> c;
  
  if ( c != '>' )
  {
    return 1;
  }

  string line;
  getline(f,line);
  string namestring = trim(line);
  stringstream n;
  n << namestring;
  n >> name;

  stringstream s;
  while ( f.good() )
  {
    f >> std::ws; // read leading white space
    c = f.peek();
    if ( c =='>' || c==EOF)
      break;
    else {
      f >> c;
      s << c;
    }
  }
  s >> sequence;
  return 0;
}

int main(int argc, char* argv[])
{
  if ( argc != 2 )
  {
    cerr << "Error: usage fix-fasta file" << endl;
    exit(1);
  }

  stringstream s;
  s << argv[1];
  string filename;
  s >> filename;
  ifstream f(filename);
  map<int,string> nameToSequenceMap;
  int name;
  string sequence;
  while ( readFastaSequence(f,name,sequence) == 0 )
    nameToSequenceMap[name] = sequence;
  for ( map<int,string>::iterator p=nameToSequenceMap.begin(); p!= nameToSequenceMap.end(); ++p )
  {
    cout << "> " << p->first << endl << p->second << endl;
  }

  return 0;
}
