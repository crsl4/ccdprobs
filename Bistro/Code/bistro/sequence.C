#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <locale> // isdigit(), tolower()
#include <random>

#include "sequence.h"
#include "model.h"

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
  name = trim(line);

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

// fixit: Change baseFreqencies to return a Vector4d
vector<double> Alignment::baseFrequencies()
{
  map<char,int> counts;
  vector<double> pi(4);
  for ( vector<Sequence*>::iterator p=sequences.begin(); p!= sequences.end(); ++p )
    for ( int k=0; k<(*p)->getSize(); ++k )
      counts[(*p)->getBase(k)]++;
  int i = 0;
  cerr << "read all counts" << endl;
  for ( map<char,int>::iterator m=counts.begin(); m !=counts.end(); ++m )
  {
    pi[i] = (m->second)/(counts.size()*(double)numSites);
    i++;
  }
  return pi;
}

bool baseNotUncertain(char a)
{
  return ( (a=='a') || (a=='c') || (a=='g') || (a=='t') );
}

int baseToInt(char a)
{
  switch ( a )
  {
  case 'a' : return 0;
  case 'c' : return 1;
  case 'g' : return 2;
  case 't' : return 3;
  default :
    cerr << "Error: baseToInt called on character that is not 'a', 'c', 'g', or 't'." << endl;
    exit(1);
  }
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


// ignore constant from log stationary prob in log likelihood
void gtrLogLike(Matrix4i& dnaCounts,QMatrix& qmatrix,double t,double& logl,double& dlogl,double& ddlogl)
{
  Matrix4d P = qmatrix.getTransitionMatrix(t);
  Matrix4d QP = qmatrix.getQP(t);
  Matrix4d QQP = qmatrix.getQQP(t);
  logl=0;
  dlogl=0;
  ddlogl=0;
  for ( int i=0; i<4; ++i)
  {
    for ( int j=0; j<4; ++j )
    {
      logl += dnaCounts(i,j) * log( P(i,j) );
      dlogl += dnaCounts(i,j) * QP(i,j) / P(i,j);
      ddlogl += dnaCounts(i,j) * (P(i,j)*QQP(i,j) - QP(i,j)*QP(i,j)) / (P(i,j)*P(i,j));
    }
  }
}

double gtrMLE(Matrix4i& dnaCounts,QMatrix& model,double t0)
{
  double curr,prop;
  double curr_logl,prop_logl;
  double curr_dlogl,prop_dlogl;
  double curr_ddlogl,prop_ddlogl;
  prop = t0;
  gtrLogLike(dnaCounts,model,prop,prop_logl,prop_dlogl,prop_ddlogl);
  do
  {
    curr = prop;
    curr_logl = prop_logl;
    curr_dlogl = prop_dlogl;
    curr_ddlogl = prop_ddlogl;
    double delta = -curr_dlogl / curr_ddlogl;
    prop = curr + delta;
    while ( prop < 0 )
    {
      delta = 0.5*delta;
      prop = curr + delta;
    }
    gtrLogLike(dnaCounts,model,prop,prop_logl,prop_dlogl,prop_ddlogl);
    while ( ( fabs(prop_dlogl) > fabs(curr_dlogl) ) && fabs(curr - prop) >1.0e-8 )
    {
      delta = 0.5*delta;
      prop = curr + delta;
      gtrLogLike(dnaCounts,model,prop,prop_logl,prop_dlogl,prop_ddlogl);
    }
  } while ( fabs(curr - prop) > 1.0e-8);
  return prop;
}

// use init for initial lengths
// store mle distances in gtr
void Alignment::calculateGTRDistancesUsingWeights(vector<int>& weights,QMatrix model,MatrixXd& init,MatrixXd& gtr)
{
  for ( int i=0; i<numTaxa; ++i ) // when i==numTaxa-1, j loop not executed
  {
    gtr(i,i) = 0;
    for ( int j=i+1; j<numTaxa; ++j )
    {
      int numTotal = 0;
      Matrix4i dnaCounts;
      dnaCounts << 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
      for ( int k=0; k<numSites; ++k )
      {
	if ( weights[k] > 0 )
	{
	  char a = tolower( getBase(i+1,k) );
	  char b = tolower( getBase(j+1,k) );
	  if ( baseNotUncertain(a) && baseNotUncertain(b) )
	  {
	    numTotal += weights[k];
	    dnaCounts( baseToInt(a), baseToInt(b) ) += weights[k];
	  }
	}
      }
      if ( numTotal == 0 )
      {
	cerr << "Error: sequences " << i+1 << " and " << j+1 << " have no sites in common with certain bases." << endl;
	exit(1);
      }
      gtr(i,j) = gtr(j,i) = gtrMLE(dnaCounts,model,init(i,j));
    }
  }
}

void Alignment::calculateGTRDistances(QMatrix model,MatrixXd& init,MatrixXd& gtr)
{
  vector<int> weights(numSites,1);
  calculateGTRDistancesUsingWeights(weights,model,init,gtr);
}

// matrix has been initialized to zero outside
MatrixXi Alignment::calculatePairwiseCounts(int i, int j)
{
  MatrixXi jc=MatrixXi::Zero(4,4);
  for ( int k=0; k<numSites; ++k )
  {
    char a = tolower( getBase(i+1,k) );
    char b = tolower( getBase(j+1,k) );
    if ( baseNotUncertain(a) && baseNotUncertain(b) )
      jc(baseToInt(a),baseToInt(b))++;
  }
  return jc;
}

VectorXd Alignment::averagePairwiseS()
{
  int total=0;
  VectorXd s = VectorXd::Zero(6);
  for ( int i=0; i<numTaxa; ++i ) // when i==numTaxa-1, j loop not executed
  {
    for ( int j=i+1; j<numTaxa; ++j )
    {
      MatrixXi mat = calculatePairwiseCounts(i,j);
      VectorXd foo(6);
      foo(0) = mat(0,1) + mat(1,0);
      foo(1) = mat(0,2) + mat(2,0);
      foo(2) = mat(0,3) + mat(3,0);
      foo(3) = mat(1,2) + mat(2,1);
      foo(4) = mat(1,3) + mat(3,1);
      foo(5) = mat(2,3) + mat(3,2);
      foo /= foo.sum();
      ++total;
      s += foo;
      cerr <<  i << " " << j << ": " << foo.transpose() << endl;
    }
  }
  s /= total;
  return s;
}
