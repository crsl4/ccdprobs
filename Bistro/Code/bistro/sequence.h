#ifndef __SEQUENCE_H
#define __SEQUENCE_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <string>
#include <random>

#include "Eigen/Core"
#include "Eigen/Eigenvalues"

class QMatrix;

using namespace std;
using namespace Eigen;

class Sequence
{
private:
  string name;
  string sequence;
public:
  Sequence() {}
  Sequence(string n,string s) : name(n), sequence(s) {}
  int readFastaSequence(istream&);
  string getName() const { return name; }
  int getSize() const { return sequence.size(); }
  string getSequence() const { return sequence; }
  char getBase(int i) const { return sequence[i]; }
  void print(ostream&);
};

class Alignment
{
private:
  int numTaxa;
  int numSites;
  vector<Sequence*> sequences;
public: // data in the class, but inconvenient to make private
  map<string,int> nameToNumberMap;
  map<int,string> numberToNameMap;
public:
  Alignment(string);
  void addSequence(Sequence*);
  void makeMaps();
  int getNumTaxa() const { return numTaxa; }
  int getNumSites() const { return numSites; }
  char getBase(int taxon,int site) const { return sequences[taxon-1]->getBase(site); } // taxon number is one-based
  void summarize(ostream&);
  void calculateJCDistancesUsingWeights(vector<int>&,MatrixXd&);
  void calculateJCDistances(MatrixXd&);
  void calculateGTRDistancesUsingWeights(vector<int>&,QMatrix,MatrixXd&,MatrixXd&);
  void calculateGTRDistances(QMatrix,MatrixXd&,MatrixXd&);
  void setBootstrapWeights(vector<int>&,mt19937_64&);
  void getTaxaNumbersAndNames(vector<int>&,vector<string>&);
  vector<double> baseFrequencies();
  void calculatePairwiseCounts(int,int,MatrixXd&);
};
#endif
