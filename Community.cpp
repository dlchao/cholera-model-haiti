/*
 * Community.cpp
 * 11/2010
 * Contains classes for Person and Community (well-mixed collection of people)
 */

#include <iostream>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "Community.h"

using namespace std;

const int Person::MAXPERSONS = 10000000;
Person Person::personArray[MAXPERSONS];
int Person::_nNextID = 0;
int Person::_nUsedID = 0;
const int Community::MAXFAMILYSIZE = 10;
double Community::_fFamiltySizeCDF[Community::MAXFAMILYSIZE] = 
  {1.0,1.0,1.0,1.0,1.0,
   1.0,1.0,1.0,1.0,1.0}; // cumulative fraction of families of size n+1
int Community::familyArray[Person::MAXPERSONS];
double Community::_fWorkTimeFraction = 0.3;
double Community::_fHouseholdContactProbability = 0.0;
double Community::_fVibrioDecay = 1.0-1.0/30.0;
double Community::_fVibrio50 = 650000;
double Community::_fHyperVibrioMultiplier=700.0;
double Community::_fVibrioBeta=1.0;
int Community::_nNextID = 0;
double Community::_fSymptomaticFraction = 0.1;
double Community::_fAsymptomaticInfectiousnessMultiplier = 0.1;
int Community::_nMaxRuralPop=2505;
double Community::_fRuralVibrio50Multiplier=1.0;
double Community::_fRuralSheddingMultiplier=1.0; // not used yet
double Community::_fRiverBetaMultiplier=1.0;
double Community::_fRiverSheddingMultiplier=1.0;

int nMinInfectiousDays = 7;
int nMaxInfectiousDays = 14;
double fWithdrawFraction = 0.75;
double Community::_fVES = 0.0;
double Community::_fVEI = 0.5;
double Community::_fVEP = 0.64;
const int _nMaxIncubationDays = 5;
double _fIncubationCDF[_nMaxIncubationDays] = {0.4,0.4,0.07,0.07,0.06};
const int Community::_nMaxVaccineDays = 21;
double Community::_fVaccineBuildup[Community::_nMaxVaccineDays] = 
  {0.001274,0.008382,0.02524,0.05516,0.1012,0.1661,0.2525,0.363,0.5,0.5,0.5,0.5,
   0.5,0.5,0.5025,0.5166,0.55,0.6092,0.7003,0.8288,1.0};
int _fIncubationCDF100[100] = 
  {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
   1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
   2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
   2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
   3,3,3,3,3,3,3,4,4,4,4,4,4,4,5,5,5,5,5,5};

void Person::vaccinate() { 
  if (_nVaccinationDay<0) 
    _nVaccinationDay=0; 
}

void Person::prevaccinate() { 
  if (_nVaccinationDay<0) {
    _nVaccinationDay=Community::getMaxVaccineDays()-1; 
    _fInfectiousness = _fBaseInfectiousness * (1.0-Community::_fVEI*Community::_fVaccineBuildup[_nVaccinationDay]);
    _fSusceptibility = _fBaseSusceptibility * (1.0-Community::_fVES*Community::_fVaccineBuildup[_nVaccinationDay]);
  }
}

void Person::infect(gsl_rng *rng) {
  if (getSusceptibility()>0.0) {
    int firstday = _fIncubationCDF100[gsl_rng_uniform_int(rng, 100)];
    _nIncubationCountdown=firstday;
    _fBaseSusceptibility = _fSusceptibility = 0.0;
  }
}

// makesymptomatic - turns this (possibly uninfected) person into a symptomatic case
int Person::makesymptomatic(gsl_rng *rng) {
  if (_bSymptomatic)
    return 0;
  _nInfectiousCountdown = nMinInfectiousDays+gsl_rng_uniform_int(rng, nMaxInfectiousDays-nMinInfectiousDays);
  _nInfectiousCountup=0;
  _fBaseSusceptibility=_fSusceptibility=0.0;
  _bSymptomatic = true;
  _bCase = true;
  _fBaseInfectiousness=_fInfectiousness=1.0;
  return 1;
}

// makeImmune - makes this person immune to infection
void Person::makeImmune() {
  _fBaseSusceptibility=_fSusceptibility=0.0;
}

// copyInfectionStatus - copies infection and vaccination data from p
void Person::copyInfectionStatus(const Person &p) {
  _nVaccinationDay = p.getVaccinationDay();
  _nIncubationCountdown = p.getIncubationCountdown();
  _nInfectiousCountdown = p.getInfectiousCountdown();
  _nInfectiousCountup = p.getInfectiousCountup();
  _fInfectiousness = p.getInfectiousness();
  _fSusceptibility = p.getSusceptibility();
  _fBaseInfectiousness = p.getBaseInfectiousness();
  _fBaseSusceptibility = p.getBaseSusceptibility();
  _bSymptomatic = p.isSymptomatic();
  _bCase = p.isCase();
}

// step - update timers and internal state for this person
bool Person::step(gsl_rng *rng) {
  if (_nInfectiousCountdown>0) {
    _nInfectiousCountdown--;
    _nInfectiousCountup++;
  } else if (_nInfectiousCountdown==0) {
    _fBaseInfectiousness=_fInfectiousness=0.0;
    _bSymptomatic = false;
  }
  if (_nIncubationCountdown>0)
    _nIncubationCountdown--;
  else if (_nIncubationCountdown==0) {
    _nIncubationCountdown=-1;
    if ((_nVaccinationDay<0 && 
	 gsl_rng_uniform(rng)<Community::getSymptomaticFraction()) ||
	(_nVaccinationDay>=0 && 
	 gsl_rng_uniform(rng)<(1.0-Community::_fVEP*Community::_fVaccineBuildup[_nVaccinationDay])*Community::getSymptomaticFraction())) {
      _bSymptomatic = true;
      _bCase = true;
      _fBaseInfectiousness=_fInfectiousness=1.0;
    } else {
      _bSymptomatic = false;
      _bCase = false;
      _fBaseInfectiousness=_fInfectiousness=Community::getAsymptomaticInfectiousnessMultiplier();
    }
    _nInfectiousCountdown = nMinInfectiousDays+gsl_rng_uniform_int(rng, nMaxInfectiousDays-nMinInfectiousDays);
    _nInfectiousCountup=0;
  }
  if (_nVaccinationDay>=0) {
    if (_nInfectiousCountdown>0)
      _fInfectiousness = _fBaseInfectiousness * (1.0-Community::_fVEI*Community::_fVaccineBuildup[_nVaccinationDay]);
    else if (_fBaseSusceptibility>0.0)
      _fSusceptibility = _fBaseSusceptibility * (1.0-Community::_fVES*Community::_fVaccineBuildup[_nVaccinationDay]);
  }  
  if (_nVaccinationDay>=0 && _nVaccinationDay<Community::getMaxVaccineDays()-1)
    _nVaccinationDay++;
  return true;
}

// setFamilySizeCDF
// set cumulative density function of family sizes, with the first
// entry being fraction of families of size 1, the second being the
// fraction of size 2 or less, etc.
// n is the size of the array.
bool Community::setFamilySizeCDF(double *f, int n) {
  if (n>MAXFAMILYSIZE)
    return false;
  for (int i=0; i<n; i++)
    _fFamiltySizeCDF[i] = f[i];
  for (int i=n; i<MAXFAMILYSIZE; i++)
    _fFamiltySizeCDF[i] = 1.0;
  return true;
}

int Community::populate(gsl_rng *rng, int nNumPeople) {
  int familycount=0;
  int leftinfamily=0;
  int familysize=0;

  _nNumResidents=0;
  _nNumVisitors=0;
  _residents = new int[nNumPeople];
  _nMaxVisitors=nNumPeople/2+1;
  _visitors = new int[_nMaxVisitors];
  
  for (int i=nNumPeople; i>0; i--) {
    int id = Person::getNextPersonID();
    _residents[_nNumResidents++] = id;
    Person &person = Person::personArray[id];
    person.setHomeCommunity(_nID);
    person.setWorkCommunity(_nID);
    if (leftinfamily<=0) {
      double r = gsl_rng_uniform(rng);
      for (; leftinfamily < MAXFAMILYSIZE && _fFamiltySizeCDF[leftinfamily]<r; leftinfamily++)
	;
      familysize = leftinfamily+1;
      familycount++;
    } else
      leftinfamily--;
    familyArray[id] = familycount;
    person.setFamilySize(familysize);
    person.setFamily(familycount);
  }
  return 1;
}

void Community::addVisitor(int id) {
  if (_nNumVisitors+1>=_nMaxVisitors) {
    int *temp = new int[_nMaxVisitors*2];
    memcpy(temp, _visitors, sizeof(int)*_nNumVisitors);
    _nMaxVisitors*=2;
    delete [] _visitors;
    _visitors=temp;
  }
  _visitors[_nNumVisitors++] = id;
}

int Community::getNumSusceptible() {
  int total = 0;
  for (int i=0; i<_nNumResidents; i++) {
    Person &p = Person::personArray[_residents[i]];
    if (p.getSusceptibility()>0.0)
      total++;
  }
  return total;
}

int Community::getNumImmune() {
  int total = 0;
  for (int i=0; i<_nNumResidents; i++) {
    Person &p = Person::personArray[_residents[i]];
    if (p.getSusceptibility()==0.0)
      total++;
  }
  return total;
}

int Community::getNumInfectious() {
  int total = 0;
  for (int i=0; i<_nNumResidents; i++) {
    Person &p = Person::personArray[_residents[i]];
    if (p.getInfectiousness()>0.0)
      total++;
  }
  return total;
}

int Community::getNumSymptomatic() {
  int total = 0;
  for (int i=0; i<_nNumResidents; i++) {
    Person &p = Person::personArray[_residents[i]];
    if (p.isSymptomatic())
      total++;
  }
  return total;
}

int Community::getNumNewSymptomatic() {
  int total = 0;
  for (int i=0; i<_nNumResidents; i++) {
    Person &p = Person::personArray[_residents[i]];
    if (p.isSymptomatic() && p.getInfectiousCountup()==0)
      total++;
  }
  return total;
}

int Community::getCumulativeSymptomatic() {
  int total = 0;
  for (int i=0; i<_nNumResidents; i++) {
    Person &p = Person::personArray[_residents[i]];
    if (p.isCase())
      total++;
  }
  return total;
}

int Community::getNumVaccinated() {
  int total = 0;
  for (int i=0; i<_nNumResidents; i++) {
    Person &p = Person::personArray[_residents[i]];
    if (p.isVaccinated())
      total++;
  }
  return total;
}

// vaccinate - vaccinates fraction of the population
// returns number vaccinated
int Community::vaccinate(gsl_rng *rng, double f) {
  int nNumVaccinate=0;
  for (int i=0; i<_nNumResidents; i++) {
    if (gsl_rng_uniform(rng)<f) {
      Person &p = Person::personArray[_residents[i]];
      if (!p.isVaccinated()) {
	p.vaccinate();
	nNumVaccinate++;
      }
    }
  }
  return nNumVaccinate;
}

int Community::prevaccinate(gsl_rng *rng, double f) {
  int nNumVaccinate=0;
  for (int i=0; i<_nNumResidents; i++) {
    if (gsl_rng_uniform(rng)<f) {
      Person &p = Person::personArray[_residents[i]];
      if (!p.isVaccinated()) {
	p.prevaccinate();
	nNumVaccinate++;
      }
    }
  }
  return nNumVaccinate;
}

// setHygiene - sets the hygiene levels of the entire population
// returns population size
int Community::setHygiene(double f) {
  for (int i=0; i<_nNumResidents; i++) {
    Person &p = Person::personArray[_residents[i]];
    p.setHygiene(f);
  }
  return _nNumResidents;
}

// infect
// n is the total number of new infections in an unvaccinated population.
// vaccinated people have lower chance of being infected.
// returns actual number of people infected.
int Community::infect(gsl_rng *rng, int n, bool bMakeSymptomatic) {
  int count=0;
  if (_nNumResidents>0) {
    for (int i=0; i<n; i++) {
      int id = _residents[gsl_rng_uniform_int(rng, _nNumResidents)]; // person to try infecting
      Person &p = Person::personArray[id]; 
      if (bMakeSymptomatic) {
	p.makesymptomatic(rng);
	count++;
      } else {
	double s = p.getSusceptibility();
	if (s==1.0 || (s>0.0 && gsl_rng_uniform(rng)<s)) {
	  p.infect(rng);
	  count++;
	}
      }
    }
  }
  return count;
}

// poop
// update amount of Vibrio in the environment
// returns the amount of new poop
double Community::poop(gsl_rng *rng) {
  //  cerr << "poop " << _nNumResidents << "," << _nID << endl;
  double hf = (1.0-_fWorkTimeFraction); // time spent at home
  _fVibrio*=_fVibrioDecay;              // level of Vibrio goes down
  _fHyperVibrio=0.0;                    // hyperinfectious vibrio disappears each day
  double newvibrio = 0.0;
  for (int i=0; i<_nNumResidents; i++) {
    assert(_residents[i]<Person::getLastPersonID());
    Person &p = Person::personArray[_residents[i]];
    //    cerr << "comm " << _nID << ", person " << _residents[i] << ": " << i <<  " / " << _nNumResidents << endl;
    if (p.getInfectiousness()>0.0) {
      if (p.getWorkCommunity()==_nID) { // are they here all day?
	newvibrio += p.getInfectiousness();
	_fHyperVibrio += p.getInfectiousness();
      } else {
	newvibrio += hf*p.getInfectiousness();
	_fHyperVibrio += hf*p.getInfectiousness();
      }
    }
  }

  for (int i=0; i<_nNumVisitors; i++) {
    assert(_visitors[i]<Person::getLastPersonID());
    Person &p = Person::personArray[_visitors[i]];
    if (p.getWorkCommunity()==_nID) { // make sure they are still working here
      newvibrio += _fWorkTimeFraction*p.getInfectiousness();
      _fHyperVibrio += _fWorkTimeFraction*p.getInfectiousness();
    }
  }
  if (_bRiver) {
    newvibrio *= _fRiverSheddingMultiplier;
    _fHyperVibrio *= _fRiverSheddingMultiplier;
  }
  _fVibrio += newvibrio;
  return newvibrio;
}

// drink
// infect people (from environment and households)
void Community::drink(gsl_rng *rng) {
  double v = (_fVibrio + _fHyperVibrioMultiplier*_fHyperVibrio)/_nNumResidents + _fRiverVibrio + _fHyperVibrioMultiplier*_fRiverHyperVibrio;
  assert(_fRiverVibrio>=0.0);
  assert(_fRiverHyperVibrio>=0.0);
  ///  if (_nNumResidents<=_nMaxRuralPop) // rural or urban population?
  //    p = _fVibrioBeta * v/(_fVibrio50*_fRuralVibrio50Multiplier + v);
  //  else
  double p = _fVibrioBeta * v/(_fVibrio50 + v);
  //  if (_bRiver)
  //    p *= _fRiverBetaMultiplier * log(_nNumResidents);

  // transmission from environment
  if (p>0.0) {
    double ph = p*(1.0-_fWorkTimeFraction);
    double pw = p*_fWorkTimeFraction;
    for (int i=0; i<_nNumResidents; i++) {
      Person &person = Person::personArray[_residents[i]];
      if (person.getSusceptibility()>0.0 && 
	  gsl_rng_uniform(rng)<(person.getWorkCommunity()==_nID?p:ph)*person.getSusceptibility()*(1.0-person.getHygiene()))
	person.infect(rng);
    }
    for (int i=0; i<_nNumVisitors; i++) {
      Person &person = Person::personArray[_visitors[i]];
      if (person.getWorkCommunity()==_nID &&
	person.getSusceptibility()>0.0 && 
	  gsl_rng_uniform(rng)<pw*person.getSusceptibility()*(1.0-person.getHygiene()))
	person.infect(rng);
    }
  }

  // transmission within families
  if (Community::_fHouseholdContactProbability>0.0) {
    for (int i=0; i<_nNumResidents; i++) {
      Person &person = Person::personArray[_residents[i]];
      // find infected person
      if (person.getInfectiousness()>0.0) {
	// find other infected people in the family
	double escapeprob = 1.0;
	int family = person.getFamily();
	int familylast; // id+1 of last person in family
	int familyfirst;// id of first person in family
	for (familylast=i; familylast<_nNumResidents && familylast<i+person.getFamilySize(); familylast++) {
	  Person &person2 = Person::personArray[_residents[familylast]];
	  if (person2.getFamily()!=family)
	    break;
	  if (person2.getInfectiousness()>0.0)
	    escapeprob *= (1.0-Community::_fHouseholdContactProbability*person2.getInfectiousness());
	}
	for (familyfirst=i-1; familyfirst>=0 && familyfirst>i-person.getFamilySize(); familyfirst--) {
	    Person &person2 = Person::personArray[_residents[familyfirst]];
	    if (person2.getFamily()!=family)
	      break;
	}
	familyfirst++;
	// expose all susceptibles in family
	double infectprob = 1.0-escapeprob;
	//	cerr << " person: " << person.getID() << ", family size: " << person.getFamilySize() << ",  p=" << Community::_fHouseholdContactProbability << ", infprob=" << infectprob << endl;
	if (infectprob>0.0) {
	  for (int j=familyfirst; j<familylast; j++) {
	    Person &person2 = Person::personArray[_residents[j]];
	    if (person2.getSusceptibility()>0.0 && 
		gsl_rng_uniform(rng)<infectprob*person2.getSusceptibility()*(1.0-person2.getHygiene()))
	      person2.infect(rng);
	  }
	}
	// go to next family
	i = familylast;
      }
    }
  }
}

    //    escape[0] *= 1.0-fContactProbability*_fWorkTimeFraction*Person::personArray[_residents[i]].getInfectiousness();
    //    escape[1] *= 1.0-fContactProbability*(1.0-_fWorkTimeFraction)*Person::personArray[_residents[i]].getInfectiousness();

  /*  // probability of escaping infection from residents
  for (int i=0; i<_nNumResidents; i++) {
    escape[0] *= 1.0-fContactProbability*_fWorkTimeFraction*Person::personArray[_residents[i]].getInfectiousness();
    escape[1] *= 1.0-fContactProbability*(1.0-_fWorkTimeFraction)*Person::personArray[_residents[i]].getInfectiousness();
  }
  // probability of escaping infection from residents who work elsewhere
  for (int i=0; i<_nNumVisitors; i++)
    escape[2] *= 1.0-fContactProbability*(1.0-_fWorkTimeFraction)*Person::personArray[_visitors[i]].getInfectiousness();
  // probability of escaping infection from visiting visitors
  for (int i=0; i<_nNumVisitingVisitors; i++)
    escape[3] *= 1.0-fContactProbability*_fWorkTimeFraction*Person::personArray[_visitingvisitors[i]].getInfectiousness();
*/

  /*
  // "nighttime" infection
  FOI = 1.0-escape[1]*escape[2];
  if (FOI>0.0) {
    for (int i=0; i<_nNumResidents; i++) {
      Person &p = Person::personArray[_residents[i]];
      if (p.getSusceptibility()>0.0 && gsl_rng_uniform(rng)<FOI*p.getSusceptibility())
	p.infect(rng);
    }
    for (int i=0; i<_nNumVisitors; i++) {
      Person &p = Person::personArray[_visitors[i]];
      if (p.getSusceptibility()>0.0 && gsl_rng_uniform(rng)<FOI*p.getSusceptibility())
	p.infect(rng);
    }
  }
  for (int i=0; i<_nNumResidents; i++)
    Person::personArray[_residents[i]].step(rng);
  for (int i=0; i<_nNumVisitors; i++)
    Person::personArray[_visitors[i]].step(rng);
   */

// tick - advance counters for each resident
int Community::tick(gsl_rng *rng) {
  //  cerr << "tick" << endl;
  int count=0;
  for (int i=0; i<_nNumResidents; i++) {
    Person &p = Person::personArray[_residents[i]];
    p.step(rng);
    if (p.getInfectiousCountup()==0) {
      count++;
    } else if (p.getDaysInfectious()==0 && p.isSymptomatic() &&
	       p.getHomeCommunity()!=p.getWorkCommunity() &&
	       gsl_rng_uniform(rng)<fWithdrawFraction) {
      // this worker withdraws
      // equivalent to turning the worker into a non-worker
      p.setWorkCommunity(p.getHomeCommunity());
    }
  }
  return count;
}

/*
int main(void) {
  gsl_rng * rng = gsl_rng_alloc(gsl_rng_taus2);
  int nNumPops = 10;
  Community *pop = new Community[nNumPops];
  for (int i=0; i<nNumPops; i++)
    pop[i].addSusceptibles(1000);

  for (int i=0; i<35; i++) {
    int id = Person::getNextPersonID();
    pop[0].addWorker(id);
    pop[1].addVisitor(id);
  }
  for (int i=0; i<50; i++) {
    int id = Person::getNextPersonID();
    pop[1].addWorker(id);
    pop[0].addVisitor(id);
  }

  // seed infected
  for (int i=0; i<20; i++)
    Person::personArray[i].infect(rng);

  // main loop
  for (int day=0; day<50; day++) {
    for (int i=0; i<nNumPops; i++)
      pop[i].step(rng);
    for (int i=0; i<Person::getLastPersonID(); i++)
      Person::personArray[i].step(rng);

    // output
    for (int i=0; i<nNumPops; i++) {
      int p = pop[i].getNumInfectious();
      if (p>0)
	cout << day << "," << i << "," << p << endl;
    }
  }
  delete [] pop;
  delete rng;
  return 0;
}
*/
