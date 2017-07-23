#ifndef __PARAMETER_H
#define __PARAMETER_H

#include <iostream>
#include <iomanip>
#include <string>
#include <cstring>
#include <vector>
#include <sstream>

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
  double weightScale;
  unsigned int numThreads;
  double eta;
  bool rootFix;
  bool weightMean;
  bool onlyBootstrap;
  bool enteredP;
  bool enteredS;
  bool doMCMC;
  int numMCMC;
  string mcmcfile;
  string mbfile; //for main-distances.C
  string bistrofile; //for main-distances.C
  bool onlyMCMC;
  bool onlyMean; //for main-distances.C
  int skip; //for main-distances.C
  bool noMCMCThreads;
public:
  Parameter()
  {
    stationaryP.resize(4,0.25);
    symmetricQP.resize(6,0.1);
    symmetricQP[1] = symmetricQP[4] = 0.3;
    seed = 0;
    numBootstrap = 0;
    numRandom = 0;
    numMLE = 2;
    outFileRoot = (string)("run1");
    independent = false;
    jointMLE = false;
    reweight = true;
    weightScale = 0;
    numThreads = 0;
    eta = 1.0;
    rootFix = false;
    weightMean = false;
    onlyBootstrap = false;
    enteredP = false;
    enteredS = false;
    doMCMC = true;
    numMCMC = 1000;
    onlyMCMC = false;
    onlyMean = false;
    skip=0;
    noMCMCThreads = false;
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
  double getWeightScale() const { return weightScale; }
  void setWeightScale(double x) { weightScale = x; }
  void processCommandLine(int,char* []);
  void setNumThreads(unsigned int b) { numThreads = b; }
  unsigned int getNumThreads() const { return numThreads; }
  double getEta() const { return eta; }
  void setEta(double x) { eta = x; }
  bool getRootFix() const { return rootFix; }
  void setRootFix(bool b) { rootFix = b; }
  bool getWeightMean() const { return weightMean; }
  void setWeightMean(bool b) { weightMean = b; }
  bool getOnlyBootstrap() const { return onlyBootstrap; }
  void setOnlyBootstrap(bool b) { onlyBootstrap = b; }
  bool getEnteredP() const { return enteredP; }
  void setEnteredP(bool b) { enteredP = b; }
  bool getEnteredS() const { return enteredS; }
  void setEnteredS(bool b) { enteredS = b; }
  bool getDoMCMC() const { return doMCMC; }
  void setDoMCMC(bool b) { doMCMC = b; }
  int getNumMCMC() { return numMCMC; }
  void setNumMCMC(int x) { numMCMC = x; }
  string getMCMCfile() const { return mcmcfile; }
  void setMCMCfile(string x) { mcmcfile = x; }
  string getMBfile() const { return mbfile; }
  void setMBfile(string x) { mbfile = x; }
  string getBistrofile() const { return bistrofile; }
  void setBistrofile(string x) { bistrofile = x; }
  bool getOnlyMCMC() const { return onlyMCMC; }
  void setOnlyMCMC(bool b) { onlyMCMC = b; }
  bool getOnlyMean() const { return onlyMean; }
  void setOnlyMean(bool b) { onlyMean = b; }
  int getSkip() const { return skip; }
  void setSkip(int n) { skip = n; }
  bool getNoMCMCThreads() const { return noMCMCThreads; }
  void setNoMCMCThreads(bool b) { noMCMCThreads = b; }

};


#endif
