/*
 * Community.h
 * 11/2010
 */
#ifndef __COMMUNITY_H
#define __COMMUNITY_H

#include <math.h> 
#include <string.h> 
#include <gsl/gsl_rng.h>

using namespace std;

class Person;
class Community;

class Person {
 public:
  Person() {
    _nID = _nNextID++;
    _nVaccinationDay=-1;
    _nIncubationCountdown=-1;
    _nInfectiousCountdown=-1;
    _nInfectiousCountup=-1;
    _fInfectiousness=0.0;
    _fSusceptibility=1.0;
    _fBaseInfectiousness=0.0;
    _fBaseSusceptibility=1.0;
    _bSymptomatic=false;
    _bCase=false;
    _nFamilySize=-1;
    _nFamilyID=-1;
    _nHomeCommunity=_nWorkCommunity=-1;
    _fHygiene = 0.0;
  }
  ~Person() {}
  int getID() const { return _nID; }
  double getInfectiousness() const { return _fInfectiousness; }
  double getBaseInfectiousness() const { return _fBaseInfectiousness; }
  double getSusceptibility() const { return _fSusceptibility; }
  double getBaseSusceptibility() const { return _fBaseSusceptibility; }
  double getHygiene() const { return _fHygiene; }
  void setHygiene(double f) { _fHygiene=f; }
  void setHomeCommunity(int n) { _nHomeCommunity=n; }
  void setWorkCommunity(int n) { _nWorkCommunity=n; }
  void setFamily(int n) { _nFamilyID=n; }
  int getFamily(void) const { return _nFamilyID; }
  int getHomeCommunity() const { return _nHomeCommunity; }
  int getWorkCommunity() const { return _nWorkCommunity; }
  void setFamilySize(int n) { _nFamilySize = n; }
  int getFamilySize() { return _nFamilySize; }
  void infect(gsl_rng *rng);
  int makesymptomatic(gsl_rng *rng);
  void makeImmune();
  bool step(gsl_rng *rng);
  bool isSymptomatic() const { return _bSymptomatic; }
  bool isCase() const { return _bCase; }
  int getInfectiousCountdown() const { return _nInfectiousCountdown; }
  int getInfectiousCountup() const { return _nInfectiousCountup; }
  int getIncubationCountdown() const { return _nIncubationCountdown; }
  void vaccinate();
  void prevaccinate();
  bool isVaccinated() { if (_nVaccinationDay<0) return false; else return true; }
  int getVaccinationDay() const { return _nVaccinationDay; }
  int getDaysInfectious() const { return _nInfectiousCountup; }
  void copyInfectionStatus(const Person &p);

  static int getNextPersonID() { assert(_nUsedID+1<MAXPERSONS); return _nUsedID++; }
  static int getLastPersonID() { return _nUsedID; }
  const static int MAXPERSONS;
  static Person personArray[];

 protected:
  int _nID;
  int _nHomeCommunity;
  int _nWorkCommunity;
  int _nVaccinationDay;
  int _nIncubationCountdown;
  int _nInfectiousCountdown;
  int _nInfectiousCountup;
  double _fInfectiousness;
  double _fSusceptibility;
  double _fBaseInfectiousness; // baseline infectiousness (without vaccination)
  double _fBaseSusceptibility; // baseline susceptibility (without vaccination)
  double _fHygiene;    // how hygienic is this person (0 is baseline, 1 is perfect)
  bool _bSymptomatic;  // is currently symptomatic
  bool _bCase;         // is or was ever symptomatic
  int _nFamilySize;
  int _nFamilyID;
  static int _nNextID;
  static int _nUsedID;
};

class Community {
 public:
  Community() {
    _nID = _nNextID++;
    _nNumResidents=0;
    _nNumVisitors=0;
    _nMaxVisitors=0;
    _fVibrio=_fHyperVibrio=_fRiverVibrio=_fRiverHyperVibrio=0.0;
    _bRiver = false;
    _residents = _visitors = NULL;
  }
  virtual ~Community() {
    if (_residents)
      delete [] _residents;
    if (_visitors)
      delete [] _visitors;
  }
  void addSusceptibles(int n);
  int populate(gsl_rng *rng, int nNumPeople);
  int employ(gsl_rng *rng, double fWorkFraction, int *neighbors, double *weights, int nNeighborLength);
  void addVisitor(int id);
  int infect(gsl_rng *rng, int n, bool bMakeSymptomatic=false);
  double poop(gsl_rng *rng);
  virtual void drink(gsl_rng *rng);
  double getVibrioLevel() { return _fVibrio; }
  double getHyperVibrioLevel() { return _fHyperVibrio; }
  double getRiverVibrioLevel() { return _fRiverVibrio; }
  double getRiverHyperVibrioLevel() { return _fRiverHyperVibrio; }
  void setVibrioLevel(double f) { _fVibrio=f; }
  void setHyperVibrioLevel(double f) { _fHyperVibrio=f; }
  void setRiverVibrioLevel(double f) { _fRiverVibrio=f; }
  void setRiverHyperVibrioLevel(double f) { _fRiverHyperVibrio=f; }
  void addRiverVibrioLevel(double f) { _fRiverVibrio+=f; }
  void addRiverHyperVibrioLevel(double f) { _fRiverHyperVibrio+=f; }
  int getNumSusceptible();
  int getNumImmune();
  int getNumInfectious();
  int getNumSymptomatic();
  int getNumNewSymptomatic();
  int getCumulativeSymptomatic();
  int getNumResidents() { return _nNumResidents; }
  int getNumVaccinated();
  int getResidentID(int i) { return _residents[i]; }
  int vaccinate(gsl_rng *rng, double f);
  int prevaccinate(gsl_rng *rng, double f);
  int setHygiene(double f);
  int tick(gsl_rng *rng);
  void setRiver(bool b) { _bRiver=b; }
  bool getRiver() { return _bRiver; }
  static void setSymptomaticFraction(double f) { _fSymptomaticFraction=f; }
  static double getSymptomaticFraction() { return _fSymptomaticFraction; }
  static void setHyperVibrioMultiplier(double f) { _fHyperVibrioMultiplier=f; }
  static void setVibrioDecay(double f) { _fVibrioDecay=f; }
  static double getVibrioDecay() { return _fVibrioDecay; }
  static void setVibrio50(double f) { _fVibrio50=f; }
  static void setBeta(double f) { _fVibrioBeta=f; }
  static double getBeta() { return _fVibrioBeta; }
  static bool setFamilySizeCDF(double *f, int n);
  static void setHouseholdContactProbability(double f) {_fHouseholdContactProbability=f;}
  static void setAsymptomaticInfectiousnessMultiplier(double f) { _fAsymptomaticInfectiousnessMultiplier=f; }
  static double getAsymptomaticInfectiousnessMultiplier() { return _fAsymptomaticInfectiousnessMultiplier; }
  static double getMaxVaccineDays() { return _nMaxVaccineDays; }

  static void setMaxRuralPop(int i) { _nMaxRuralPop=i; }
  static void setRuralVibrio50Multiplier(double f) { _fRuralVibrio50Multiplier = f; }
  static void setRuralSheddingMultiplier(double f) { _fRuralSheddingMultiplier = f; }
  static double getRiverBetaMultiplier() { return _fRiverBetaMultiplier; }
  static void setRiverBetaMultiplier(double f) { _fRiverBetaMultiplier = f; }
  static double getRiverSheddingMultiplier() { return _fRiverSheddingMultiplier; }
  static void setRiverSheddingMultiplier(double f) { _fRiverSheddingMultiplier = f; }
  static void setVEI(double f) { _fVEI = f; }
  static double getVEI() { return _fVEI; }
  static void setVES(double f) { _fVES = f; }
  static double getVES() { return _fVES; }
  static void setVEP(double f) { _fVEP = f; }
  static double getVEP() { return _fVEP; }

 public:
  static double _fVaccineBuildup[];
  static double _fVES, _fVEI, _fVEP; // vaccine efficacy

 protected:
  int _nID;
  int *_residents; // those who live here
  int *_visitors;  // those from elsewhere who work here
  int _nNumResidents;
  int _nNumVisitors;
  int _nMaxVisitors;
  double _fVibrio;
  double _fHyperVibrio;
  double _fRiverVibrio;
  double _fRiverHyperVibrio;
  bool _bRiver;
  static double _fHouseholdContactProbability;
  //  static double _fSheddingRate;
  static double _fHyperVibrioMultiplier;
  static double _fVibrioDecay;
  static double _fVibrio50;
  static double _fVibrioBeta;
  static int _nNextID;
  static double _fWorkTimeFraction;
  static const int MAXFAMILYSIZE;  // largest family size
  static int familyArray[];
  static double _fFamiltySizeCDF[];
  static double _fAsymptomaticInfectiousnessMultiplier;
  static double _fSymptomaticFraction;
  static const int _nMaxVaccineDays; // time for vaccine to reach max efficacy
  static int _nMaxRuralPop;
  static double _fRuralVibrio50Multiplier;
  static double _fRuralSheddingMultiplier;
  static double _fRiverBetaMultiplier;
  static double _fRiverSheddingMultiplier;
};
#endif
