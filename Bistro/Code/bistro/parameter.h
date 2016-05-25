#ifndef __PARAMETER_H
#define __PARAMETER_H

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

using namespace std;

class Parameter
{
private:
  string sequenceFileName;
  string topology;
  vector<double> stationaryP;
  vector<double> symmetricQP;
  unsigned int seed;
  int numBootstrap;
public:
  Parameter()
  {
    stationaryP.resize(4,0.25);
    symmetricQP.resize(6,1.0/6.0);
    seed = 0;
    numBootstrap = 0;
  }
  string getSequenceFileName() const { return sequenceFileName; }
  void setSequenceFileName(string x) { sequenceFileName = x; }
  string getTopology() const { return topology; }
  void setTopology(string x) { topology = x; }
  vector<double> getStationaryP() const { return stationaryP; }
  void setStationaryP(vector<double> p) { stationaryP = p; }
  vector<double> getSymmetricQP() const { return symmetricQP; }
  void setSymmetricQP(vector<double> qp) { symmetricQP = qp; }
  unsigned int getSeed() const { return seed; }
  void setSeed(unsigned int x) { seed = x; }
  int getNumBootstrap() { return numBootstrap; }
  void setNumBootstrap(int n) { numBootstrap = n; }

  void processCommandLine(int,char* []);
};

#endif
