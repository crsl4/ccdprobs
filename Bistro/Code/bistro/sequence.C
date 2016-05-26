#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <locale> // isdigit(), tolower()
#include <random>

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

bool baseNotUncertain(char a)
{
  return ( (a=='a') || (a=='c') || (a=='g') || (a=='t') );
}

// Assumes that matrix is correctly sized before
// Even though matrix is symmetric, store distance in both upper and lower triangle to faciliate
//   calculation of row sums
// For now, ignore any sites with one or more uncertain characters
void Alignment::calculateJCDistancesUsingWeights(vector<int>& weights,MatrixXd& jc)
{
  for ( int i=0; i<numTaxa; ++i ) // when i==numTaxa-1, j loop not executed
  {
    jc(i,i) = 0;
    for ( int j=i+1; j<numTaxa; ++j )
    {
      int numDifferent = 0;
      int numTotal = 0;
      for ( int k=0; k<numSites; ++k )
      {
	if ( weights[k] > 0 )
	{
	  char a = tolower( getBase(i+1,k) );
	  char b = tolower( getBase(j+1,k) );
	  if ( baseNotUncertain(a) && baseNotUncertain(b) )
	  {
	    numTotal += weights[k];
	    if ( a != b )
	      numDifferent += weights[k];
	  }
	}
      }
      if ( numTotal == 0 )
      {
	cerr << "Error: sequences " << i+1 << " and " << j+1 << " have no sites in common with certain bases." << endl;
	exit(1);
      }
      double p = numDifferent / (double)(numTotal);
      if ( p >= 0.75 )
      {
	cerr << "Warning: sequences " << i+1 << " and " << j+1 << " have infinite MLE distance." << endl;
      }
      jc(i,j) = -3.0 * log( 1 - 4.0*p/3.0 ) / 4.0;
      jc(j,i) = jc(i,j);
    }
  }
}

void Alignment::calculateJCDistances(MatrixXd& jc)
{
  vector<int> weights(numSites,1);
  calculateJCDistancesUsingWeights(weights,jc);
}

void Alignment::setBootstrapWeights(vector<int>& weights,mt19937_64& rng)
{
  for ( int i=0; i<weights.size(); ++i )
    weights[i] = 0;
  uniform_int_distribution<int> rint(0,numSites-1);
  for ( int i=0; i<numSites; ++i )
  {
    ++weights[rint(rng)];
  }
}

void Alignment::getTaxaNumbersAndNames(vector<int>& taxaNumbers,vector<string>& taxaNames)
{
  for ( map<int,string>::iterator m=numberToNameMap.begin(); m!= numberToNameMap.end(); ++m )
  {
    taxaNumbers.push_back(m->first);
    taxaNames.push_back(m->second);
  }
}

