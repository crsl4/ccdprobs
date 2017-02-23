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
  int numRandom;
  int numMLE;
  string outFileRoot;
  bool independent;
  bool jointMLE;
  bool reweight;
  double parsimonyScale;
  unsigned int numThreads;
  bool fixedQ;
  bool loglikwt;
  double eta;
  bool rootFix;
  bool weightMean;
  bool onlyBootstrap;
public:
  Parameter()
  {
    stationaryP.resize(4,0.25);
    symmetricQP.resize(6,0.1);
    symmetricQP[1] = symmetricQP[4] = 0.3;
    seed = 0;
    numBootstrap = 0;
    numRandom = 0;
    numMLE = 0;
    outFileRoot = (string)("run1");
    independent = false;
    jointMLE = false;
    reweight = true;
    parsimonyScale = 0.5;
    numThreads = 0;
    fixedQ = false;
    loglikwt = false;
    eta = 1.0;
    rootFix = false;
    weightMean = false;
    onlyBootstrap = false;
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
  int getNumBootstrap() const { return numBootstrap; }
  void setNumBootstrap(int n) { numBootstrap = n; }
  int getNumRandom() const { return numRandom; }
  void setNumRandom(int n) { numRandom = n; }
  int getNumMLE() const { return numMLE; }
  void setNumMLE(int n) { numMLE = n; }
  bool getIndependent() const { return independent; }
  void setIndependent(bool b) { independent = b; }
  bool getJointMLE() const { return jointMLE; }
  void setJointMLE(bool b) { jointMLE = b; }
  string getOutFileRoot() const { return outFileRoot; }
  void setOutFileRoot(string name) { outFileRoot = name; }
  bool getReweight() const { return reweight; }
  void setReweight(bool b) { reweight = b; }
  double getParsimonyScale() const { return parsimonyScale; }
  void setParsimonyScale(double x) { parsimonyScale = x; }
  void processCommandLine(int,char* []);
  void setNumThreads(unsigned int b) { numThreads = b; }
  unsigned int getNumThreads() const { return numThreads; }
  bool getFixedQ() const { return fixedQ; }
  void setFixedQ(bool b) { fixedQ = b; }
  bool getLoglikWt() const { return loglikwt; }
  void setLoglikWt(bool b) { loglikwt = b; }
  double getEta() const { return eta; }
  void setEta(double x) { eta = x; }
  bool getRootFix() const { return rootFix; }
  void setRootFix(bool b) { rootFix = b; }
  bool getWeightMean() const { return weightMean; }
  void setWeightMean(bool b) { weightMean = b; }
  bool getOnlyBootstrap() const { return onlyBootstrap; }
  void setOnlyBootstrap(bool b) { onlyBootstrap = b; }
};

#endif
