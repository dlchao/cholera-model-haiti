/*
 * Population.h
 * 11/2010
 */
#ifndef __POPULATION_H
#define __POPULATION_H

#include <math.h> 
#include <string.h> 
#include <gsl/gsl_rng.h>
#include "Community.h"

class Population {
 public:
  Population();
  virtual ~Population();
  int loadGrid(gsl_rng *rng, int nGridSize, int nCommunitySize);
  int loadLandScan(gsl_rng *rng, int nGridSizeX, int nGridSizeY, 
		   double workingfraction,
		   double tau1, double tau2, double rho,
		   const char *landscanfilename,
		   int nOriginXDeg, int nOriginXMin, 
		   int nOriginYDeg, int nOriginYMin,
		   bool bR0run=false);
  int loadRivers(const char *riverfilename);
  int loadHighways(const char *highwayfilename);
  
  int infect(gsl_rng *rng, int num, double x, double y, bool bMakeSymptomatic=false);
  int infect(gsl_rng *rng, int gridindex, double frac, bool bMakeSymptomatic=false);
  int infectOne(gsl_rng *rng, bool bMakeSymptomatic=false);
  int drive(gsl_rng *rng);
  virtual int step(gsl_rng *rng);
  int getDay() { return _nDay; }
  int getNumCells() { return _grid->getSize(); }
  int getNumCommunities() { return _nNumCommunities; }
  int getNumInfectious(int gridindex) {
    int total=0;
    for (int i=_nCommunityStart[gridindex]; i<_nCommunityStart[gridindex+1]; i++)
      total += _community[i].getNumInfectious();
    return total;
  }
  int getNumSusceptible(int gridindex) {
    int total=0;
    for (int i=_nCommunityStart[gridindex]; i<_nCommunityStart[gridindex+1]; i++)
      total += _community[i].getNumSusceptible();
    return total;
  }
  int getNumImmune(int gridindex) {
    int total=0;
    for (int i=_nCommunityStart[gridindex]; i<_nCommunityStart[gridindex+1]; i++)
      total += _community[i].getNumImmune();
    return total;
  }
  int getNumSymptomatic(int gridindex) {
    int total=0;
    for (int i=_nCommunityStart[gridindex]; i<_nCommunityStart[gridindex+1]; i++)
      total += _community[i].getNumSymptomatic();
    return total;
  }
  int getNumNewSymptomatic(int gridindex) {
    int total=0;
    for (int i=_nCommunityStart[gridindex]; i<_nCommunityStart[gridindex+1]; i++)
      total += _community[i].getNumNewSymptomatic();
    return total;
  }
  int getNumResidents(int gridindex) {
    int total=0;
    for (int i=_nCommunityStart[gridindex]; i<_nCommunityStart[gridindex+1]; i++)
      total += _community[i].getNumResidents();
    return total;
  }
  int getCumulativeSymptomatic(int gridindex) {
    int total=0;
    for (int i=_nCommunityStart[gridindex]; i<_nCommunityStart[gridindex+1]; i++)
      total += _community[i].getCumulativeSymptomatic();
    return total;
  }
  int getNumVaccinated(int gridindex) {
    int total=0;
    for (int i=_nCommunityStart[gridindex]; i<_nCommunityStart[gridindex+1]; i++)
      total += _community[i].getNumVaccinated();
    return total;
  }
  double getVibrioLevel(int gridindex) {
    double total=0.0;
    for (int i=_nCommunityStart[gridindex]; i<_nCommunityStart[gridindex+1]; i++)
      total += _community[i].getVibrioLevel(); 
    return total;
  }
  double getHyperVibrioLevel(int gridindex) {
    double total=0.0;
    for (int i=_nCommunityStart[gridindex]; i<_nCommunityStart[gridindex+1]; i++)
      total += _community[i].getHyperVibrioLevel(); 
    return total;
  }
  double getRiverVibrioLevel(int gridindex) { return _riverVibrio[gridindex]; }
  void setRiverVibrioLevel(int gridindex, double f) { if (_riverDist && _riverDist[gridindex]>0) _riverVibrio[gridindex] = f; }
  double getRiverHyperVibrioLevel(int gridindex) { return _riverHyperVibrio[gridindex]; }
  int getRiver(int gridindex) { if (_riverDist) return _riverDist[gridindex]; else return -1; }
  bool isHighway(int gridindex) { if (_bHighway) return _bHighway[gridindex]; else return false; }
  int getX(int gridindex) { return _grid->getX(gridindex); }
  int getY(int gridindex) { return _grid->getY(gridindex); }
  const char *getLabel(int gridindex) { return _szLabels[gridindex].c_str(); }
  const char *getLabel2(int gridindex) { return _szLabels2[gridindex].c_str(); }
  int getRiverShedCycles() { return _nRiverShedCycles; }
  void setRiverShedCycles(int n) { _nRiverShedCycles=n; }
  double getRiverShedFraction() { return _fRiverShedFraction; }
  void setRiverShedFraction(double f) { _fRiverShedFraction=f; }
  double getRiverFlowDelta() { return _fRiverFlowDelta; }
  void setRiverFlowDelta(double f) { _fRiverFlowDelta=f; }
  double getTravelProb() { return _fTravelProb; }
  void setTravelProb(double f) { _fTravelProb=f; }
  double getDriveProb() { return _fDriveProb; }
  void setDriveProb(double f) { _fDriveProb=f; }
  double getVaccinationTarget() { return _fVaccinationTarget; }
  void setHygieneTarget(double f) { _fHygieneTarget = f; }
  double getHygieneTarget() { return _fHygieneTarget; }
  void setVaccinationTarget(double f) { _fVaccinationTarget=f; }
  int getVaccinationThreshold() { return _nVaccinationThreshold; }
  void setVaccinationThreshold(int n) { _nVaccinationThreshold=n; }
  int getVaccinationDelay() { return _nVaccinationDelay; }
  void setVaccinationDelay(int n) { _nVaccinationDelay=n; }
  void prevaccinate(gsl_rng *rng);
  void setNumVaccinesAvailable(int n) { _nNumVaccinesAvailable=n; }
  int getNumVaccinesAvailable() { return _nNumVaccinesAvailable; }
  void prioritizeCell(int gridindex);
  bool isVaccinated(int gridindex) { return _bGridVaccinated[gridindex]; }
  bool wantVaccine(int gridindex) { return _nVaccinationDay[gridindex]<1000; }

 protected:
  int getRadius(GridCells *g, int center, int radius, int *buf, int bufsize);
  int convertDegToInt(int deg, int min, int sec);
  int initGridPopulation(gsl_rng *rng, int nCommunitySize);
  int assignWorkPlaces(gsl_rng *rng, int nSourceID, double fWorkFraction, double tau1, double tau2, double rho);
  int computedrivelength(gsl_rng *rng, int start);

  int _nGridSizeX;       // how many cells wide is the grid?
  int _nGridSizeY;       // how many cells tall is the grid?
  int _nXOrigin;         // for LandScan reading
  int _nYOrigin;
  double _fWorkingFraction; // what fraction of the total population works
  double _fTravelProb;   // daily probability that someone will travel
  GridCells *_grid;  // a single 1km^2 area, which can have multiple communities
  Community *_community; // communities are well-mixed populations
  int *_nCommunityStart; // the first community that belongs to each grid cell
  int _nNumCommunities;  // total number of communities
  string *_szLabels;     // department names for each cell
  string *_szLabels2;    // commune names for each cell
  int _nTargetGroupSize; // ideal size for a single community
  int _nDay;             // current simulation day

  int *_riverDist;           // how far is this from the river source?
  int *_riverNeighbors;      // how many cells are downstream of here?
  double *_riverVibrio;      // contains level of vibrio in each cell
  double *_riverHyperVibrio; // contains level of fresh vibrio in each cell
  double *_riverVibrio2;     // temporary data storage
  double *_riverHyperVibrio2;// temporary data storage
  double *_newVibrio;        // temporary data storage
  double _fRiverShedFraction;     // amount fresh vibrio that goes to the river (as a fraction of fresh vibrio)
  double _fRiverFlowDelta; // how much vibrio is lost per cell traveled
  int _nRiverShedCycles;   // number of times to let vibrio flow downstream
  int _nMaxCommunitiesOnRiver; // largest number of communties in a cell that can be on a river

  double _fDriveProb;   // daily probability that someone will drive on highway
  bool *_bHighway;      // is this cell on a highway?
  int *_nNumHighwayDests;   // number of driving destinations for each cell
  int **_nHighwayDests;     // list of driving destinations for each cell
  int *_nHighwayDepth;      // temporary data storage

  double _fVaccinationTarget;  // desired fraction to vaccinate in a community
  bool *_bGridVaccinated;      // has this cell been vaccinated?
  int *_nVaccinationDay;       // when to vaccinate this cell?
  int _nNumVaccinesUsed;       // how many people vaccinated so far
  int _nNumVaccinesAvailable;  // how many more people can be vaccinated
  int _nVaccinationThreshold;  // how many cases before reactive vaccination
  int _nVaccinationDelay;      // how many days before reactive vaccination
  double _fHygieneTarget;      // 
};
#endif
