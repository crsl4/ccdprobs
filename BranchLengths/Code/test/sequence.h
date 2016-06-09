#ifndef __SEQUENCE_H
#define __SEQUENCE_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <string>
#include <cstdio>
#include <cstring>
#include <cstdlib>

using namespace std;

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
public:
  map<string,int> nameToNumberMap;
  map<int,string> numberToNameMap;
  Alignment(string);
  void addSequence(Sequence*);
  void makeMaps();
  int getNumTaxa() const { return numTaxa; }
  int getNumSites() const { return numSites; }
  char getBase(int taxon,int site) const { return sequences[taxon-1]->getBase(site); } // taxon number is one-based
  void summarize(ostream&);
};
#endif
