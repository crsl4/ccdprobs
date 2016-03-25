#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>

#include "sequence.h"

using namespace std;

int Sequence::readFastaSequence(istream& f)
{
  char c;
  f >> c;
  
  if ( c != '>' )
  {
    return 1;
  }
  
  string line;
  getline(f,line);
  name = line;

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

void Sequence::print(ostream& f)

{
  f << '>' << name << endl << sequence << endl;
}

void Alignment::addSequence(Sequence* seq)
{
  if ( sequences.size() == 0 ) // first sequence
  {
    numSites = seq->getSize();    
  }
  else
  {
    if ( seq->getSize() != numSites )
    {
      cerr << "Error: sequence " << numTaxa+1 << ", " << seq->getName() << ", has length " << seq->getSize();
      cerr << " when " << numSites << " was expected." << endl;
      exit(1);
    }
  }
  ++numTaxa;
  sequences.push_back(seq);
}

void Alignment::makeMaps()
{
  for ( int k=0; k < sequences.size(); )
  {
    string name = sequences[k]->getName();
    nameToNumberMap[ name ] = ++k;
    numberToNameMap[ k ] = name;
  }
}

Alignment::Alignment(string filename)
{
  numTaxa=0;
  ifstream f(filename.c_str());
  while ( f.good() )
  {
    Sequence seq;
    if (seq.readFastaSequence(f) >0 )
      break;
    addSequence( new Sequence(seq.getName(),seq.getSequence()) );
  }
  makeMaps();
}

void Alignment::summarize(ostream& f )
{
  for ( vector<Sequence*>::iterator p=sequences.begin(); p!= sequences.end(); ++p )
  {
    map<char,int> counts;
    for ( int k=0; k<(*p)->getSize(); ++k )
      counts[(*p)->getBase(k)]++;
    f << (*p)->getName() << ":";
    for ( map<char,int>::iterator m=counts.begin(); m !=counts.end(); ++m )
      f << " " << m->first << "=" << m->second;
    f << endl;
  }
}

