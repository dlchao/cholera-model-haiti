/*
 * driver.cpp
 *
 * Each grid cell has a single cholera model.
 * Each pair of grid cells (i,j) for which residents of i work in j also has a cholera model.
 */

#include <math.h>
#include <assert.h>
#include <fstream>
#include <iostream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "GridCells.h"
#include "Community.h"
#include "Population.h"

using namespace std;

int getRadius(GridCells *g, int center, int radius, int *buf, int bufsize) {
  const int maxneighbors = 6;
  int searched[1200];
  int searchedsize=1;
  int oldsearchedsize=1;
  searched[0] = center;
  // breadth first search
  for (int r=0; r<radius; r++) {
    oldsearchedsize=searchedsize;
    for (int i=searchedsize-1; i>=0; i--) {
      int neighbors[maxneighbors];
      int numneighbors = g->getNeighbors(searched[i], neighbors, maxneighbors);
      for (int j=0; j<numneighbors; j++) {
	bool bFound=false;
	for (int k=0; k<searchedsize && !bFound; k++) {
	  if (searched[k]==neighbors[j])
	    bFound=true;
	}
	if (!bFound)
	  searched[searchedsize++] = neighbors[j];
	assert(searchedsize<1200);
      }
    }
  }
  for (int i=0; i<searchedsize-oldsearchedsize; i++)
    buf[i] = searched[oldsearchedsize+i];
  return searchedsize-oldsearchedsize;
}

int main(int argc, char *argv[]) {
  gsl_rng * rng = gsl_rng_alloc(gsl_rng_taus2);
  double familycdf[10] = {0.093,0.238,0.383,0.520,0.658,0.795,0.851,0.907,963,1.0};
  int randomseed = 5489;
  GridCells *g=NULL;
  Community *model=NULL;
  Population *pop=NULL;
  int communitysize = 1000;
  bool bMatlab = false; // use Matlab grid?
  bool bHaiti = false; // use Haiti LandScan data?
  bool bNoTransmit = false; // is cholera transmissible?
  int nStartCell = -1; // where was the first case?
  int nGridSize = -1; // use grid grid?
  double kappa = 250000;
  double hyper = 0.0;
  double beta = 1.0;
  double rho = 3.8;
  int runlength = 300;
  double fVaccinateFraction = 0.0;
  int nVaccinateThreshold = 0; // set to -1 for pre-vaccination
  int nVaccinateDelay = -1;    // set to >=0 for reactive vaccination
  int nVaccineStockpile = 1000000000;  // amount of vaccine available (after first day)
  int nVaccineCap = 1000000000;  // amount of vaccine ever available (maximum supply)
  int nVaccineFirstDay = -1;// when vaccines arrive
  int nVaccinePerDay = 0;  // amount of vaccine that arrives per day after first day
  double fHygiene = 0.0;
  int nOutputInterval = 1;
  double fWorkingFraction=0.35; // fraction of people who work
  double fRiverShedFraction=0.0;
  double fRiverFlowDelta = 0.0;
  int nRiverShedCycles=0;
  bool bPrioritizeRiver = false;
  double fTravelProb = 0.0;
  double fDriveProb = 0.0;
  string szTractFile=""; //grid-cells.csv";
  string szOutputFile="grid-infections.csv";
  string szPeopleFile="";
  string szSummaryFile="";

  if (argc>1) {
    for (int i=1; i<argc; i++) {
      char **end = NULL;
      if (strcmp(argv[i], "-p")==0) {
	double p = strtof(argv[i+1],end);
	cerr << "p = " << p << endl;
	Community::setHouseholdContactProbability(p);
	i++;
      } else if (strcmp(argv[i], "-size")==0) {
	communitysize = strtol(argv[i+1],end,10);
	cerr << "community size = " << communitysize << endl;
	i++;
      } else if (strcmp(argv[i], "-k")==0 || 
		 strcmp(argv[i], "-kappa")==0) {
	kappa = strtof(argv[i+1],end);
	cerr << "kappa = " << kappa << endl;
	i++;
      } else if (strcmp(argv[i], "-hyper")==0) {
	hyper = strtof(argv[i+1],end);
	cerr << "hyper = " << hyper << endl;
	i++;
      } else if (strcmp(argv[i], "-beta")==0) {
	beta = strtof(argv[i+1],end);
	cerr << "beta = " << beta << endl;
	i++;
      } else if (strcmp(argv[i], "-rho")==0) {
	rho = strtof(argv[i+1],end);
	cerr << "rho = " << rho << endl;
	i++;
      } else if (strcmp(argv[i], "-asymptomaticmultiplier")==0 ||
		 strcmp(argv[i], "-m")==0) {
	double shed = strtof(argv[i+1],end);
	Community::setAsymptomaticInfectiousnessMultiplier(shed);
	cerr << "asymptomatic infectiousness multiplier = " << shed << endl;
	i++;
      } else if (strcmp(argv[i], "-symptomaticfraction")==0) {
	double f = strtof(argv[i+1],end);
	Community::setSymptomaticFraction(f);
	cerr << "symptomatic fraction = " << f << endl;
	i++;
      } else if (strcmp(argv[i], "-randomseed")==0) {
	randomseed = strtol(argv[i+1],end,10);
	cerr << "random seed = " << randomseed << endl;
	i++;
      } else if (strcmp(argv[i], "-vaccinatefraction")==0) {
	fVaccinateFraction = strtof(argv[i+1],end);
	cerr << "vaccinate fraction = " << fVaccinateFraction << endl;
	i++;
      } else if (strcmp(argv[i], "-vaccinatethreshold")==0) {
	nVaccinateThreshold = strtol(argv[i+1],end,10);
	cerr << "vaccination threshold = " << nVaccinateThreshold << endl;
	i++;
      } else if (strcmp(argv[i], "-vaccinatedelay")==0) {
	nVaccinateDelay = strtol(argv[i+1],end,10);
	cerr << "vaccination delay = " << nVaccinateDelay << endl;
	i++;
      } else if (strcmp(argv[i], "-vaccinefirstday")==0) {
	nVaccineFirstDay = strtol(argv[i+1],end,10);
	cerr << "first day of vaccination = " << nVaccineFirstDay << endl;
	i++;
      } else if (strcmp(argv[i], "-vaccineperday")==0) {
	nVaccinePerDay = strtol(argv[i+1],end,10);
	cerr << "new vaccine available per day = " << nVaccinePerDay << endl;
	i++;
      } else if (strcmp(argv[i], "-workingfraction")==0) {
	fWorkingFraction = strtof(argv[i+1],end);
	cerr << "working fraction = " << fWorkingFraction << endl;
	i++;
	/*      } else if (strcmp(argv[i], "-ruralkmultiplier")==0) {
	double temp = strtof(argv[i+1],end);
	Community::setRuralVibrio50Multiplier(temp);
	cerr << "rural kappa multiplier = " << temp << endl;
	i++; */
      } else if (strcmp(argv[i], "-rivershedfraction")==0) {
	fRiverShedFraction=strtof(argv[i+1],end);
	cerr << "fraction shed into the river = " << fRiverShedFraction << endl;
	i++;
      } else if (strcmp(argv[i], "-rivershedcycles")==0) {
	nRiverShedCycles=strtol(argv[i+1],end,10);
	cerr << "river shed cycles per day = " << nRiverShedCycles << endl;
	i++;
      } else if (strcmp(argv[i], "-riverflowdelta")==0) {
	fRiverFlowDelta=strtof(argv[i+1],end);
	cerr << "fraction that disappears per distance travelled in river = " << fRiverFlowDelta << endl;
	i++;
      } else if (strcmp(argv[i], "-prioritizeriver")==0) {
	bPrioritizeRiver = true;
	cerr << "Prioritizing rivers" << endl;
	/*      } else if (strcmp(argv[i], "-rivershedmultiplier")==0) {
	double temp = strtof(argv[i+1],end);
	Community::setRiverSheddingMultiplier(temp);
	cerr << "river shedding multiplier = " << temp << endl;
	i++;
      } else if (strcmp(argv[i], "-riverbetamultiplier")==0) {
	double temp = strtof(argv[i+1],end);
	Community::setRiverBetaMultiplier(temp);
	cerr << "river beta multiplier = " << temp << endl;
	i++; */
      } else if (strcmp(argv[i], "-travelprob")==0) {
	fTravelProb = strtof(argv[i+1],end);
	cerr << "daily travel probability = " << fTravelProb << endl;
	i++;
      } else if (strcmp(argv[i], "-driveprob")==0) {
	fDriveProb = strtof(argv[i+1],end);
	cerr << "daily highway driving probability = " << fDriveProb << endl;
	i++;
      } else if (strcmp(argv[i], "-ves")==0 ||
		 strcmp(argv[i], "-VES")==0) {
	double temp = strtof(argv[i+1],end);
	cerr << "VE_S = " << temp << endl;
	Community::setVES(temp);
	i++;
      } else if (strcmp(argv[i], "-vep")==0 ||
		 strcmp(argv[i], "-VEP")==0) {
	double temp = strtof(argv[i+1],end);
	cerr << "VE_P = " << temp << endl;
	Community::setVEP(temp);
	i++;
      } else if (strcmp(argv[i], "-vei")==0 ||
		 strcmp(argv[i], "-VEI")==0) {
	double temp = strtof(argv[i+1],end);
	cerr << "VE_I = " << temp << endl;
	Community::setVEI(temp);
	i++;
      } else if (strcmp(argv[i], "-days")==0) {
	runlength = strtol(argv[i+1],end,10);
	cerr << "runlength = " << runlength << endl;
	i++;
      } else if (strcmp(argv[i], "-outputinterval")==0) {
	nOutputInterval = strtol(argv[i+1],end,10);
	cerr << "output interval (in days) = " << nOutputInterval << endl;
	i++;
      } else if (strcmp(argv[i], "-vaccinestockpile")==0) {
	nVaccineStockpile = strtol(argv[i+1],end,10);
	cerr << "vaccine stockpile size = " << nVaccineStockpile << endl;
	i++;
      } else if (strcmp(argv[i], "-vaccinecap")==0) {
	nVaccineCap = strtol(argv[i+1],end,10);
	cerr << "total vaccines = " << nVaccineCap << endl;
	i++;
      } else if (strcmp(argv[i], "-hygiene")==0) {
	fHygiene = strtof(argv[i+1],end);
	cerr << "hygiene = " << fHygiene << endl;
	i++;
      } else if (strcmp(argv[i], "-cellfile")==0) {
	szTractFile = argv[i+1];
	cerr << "cell file = " << szTractFile << endl;
	i++;
      } else if (strcmp(argv[i], "-outputfile")==0) {
	szOutputFile = argv[i+1];
	if (szOutputFile.length()>0)
	  cerr << "output file = " << szOutputFile << endl;
	else 
	  cerr << "no output file" << endl;
	i++;
      } else if (strcmp(argv[i], "-summaryfile")==0) {
	szSummaryFile = argv[i+1];
	if (szSummaryFile.length()>0)
	  cerr << "summary file = " << szSummaryFile << endl;
	else 
	  cerr << "no summary file" << endl;
	i++;
      } else if (strcmp(argv[i], "-peoplefile")==0) {
	szPeopleFile = argv[i+1];
	cerr << "people file = " << szPeopleFile << endl;
	i++;
      } else if (strcmp(argv[i], "-grid")==0) {
	nGridSize = strtol(argv[i+1],end,10);
	bMatlab=false;
	cerr << "grid dimension = " << nGridSize << endl;
	i++;
      } else if (strcmp(argv[i], "-matlab")==0) {
	bMatlab = true;
	nGridSize = -1;
	cerr << "Matlab population" << endl;
      } else if (strcmp(argv[i], "-haiti")==0) {
	bHaiti = true;
	nGridSize = -1;
	cerr << "Haiti population" << endl;
      } else if (strcmp(argv[i], "-notransmit")==0) {
	bNoTransmit=true;
	cerr << "no secondary transmission" << endl;
      } else {
	cerr << "****Unknown parameter: " << argv[i] << endl;
	exit(-1);
      }
    }
  }  
  gsl_rng_set(rng, randomseed);
  Community::setFamilySizeCDF(familycdf, 10);
  Community::setVibrioDecay(1.0-1.0/30.0);
  Community::setVibrio50(kappa);
  Community::setHyperVibrioMultiplier(hyper);
  Community::setBeta(beta);

  if (bHaiti) { // Haiti
    pop=new Population();
    pop->loadLandScan(rng, 342, 248,
		      fWorkingFraction,
		      0.3, 0.64, rho,
		      "landscan_population_dms.txt",
		      71,37,18,1,
		      bNoTransmit);
    pop->loadRivers("rivers.txt");
    pop->loadHighways("highways.txt");
    pop->setRiverShedFraction(fRiverShedFraction);
    pop->setRiverShedCycles(nRiverShedCycles);
    pop->setRiverFlowDelta(fRiverFlowDelta);
    pop->setTravelProb(fTravelProb);
    pop->setDriveProb(fDriveProb);
    pop->setVaccinationTarget(fVaccinateFraction);
    if (nVaccineFirstDay<0)
      pop->setNumVaccinesAvailable(nVaccineStockpile);
    else
      pop->setNumVaccinesAvailable(0);
    pop->setHygieneTarget(fHygiene);
    if (bPrioritizeRiver) {
      // send vaccine to river when available
      for (int i=0; i < pop->getNumCells(); i++)
	if (pop->getRiver(i)>0)
	  pop->prioritizeCell(i);
      pop->setVaccinationThreshold(100000000);
      pop->setVaccinationDelay(-1);
    } else {
      pop->setVaccinationThreshold(nVaccinateThreshold);
      pop->setVaccinationDelay(nVaccinateDelay);
    }

    // map goes from 71 37', 18 1'  to  74 28', 20 5'

    // start epidemic in Mirebalais 18°50′0″N 72°6′19″W
    //    int source = pop->infect(rng, 100, 57, 98);
    if (!bNoTransmit) {
      // seed in PETITE_RIVIERE_DE_L'ARTI, SAINT-MARC, VERRETTES, LASCAHOBAS, MIREBALAIS
      // pops of 144467, 215390, 112986, 28507, 83375
      // The first days of cas vu in Artibonite: 1111 1840 1833 2138 1675
      // First days of cas vu in Centre: 61  92 244 165
      //      int count1 = 0;
      //      int count2 = 0;
      for (int i=0; i < pop->getNumCells(); i++) {
	const char *lab2 = pop->getLabel2(i);
	if (pop->getRiver(i)>0) {
	  if (strcmp(lab2, "PETITE_RIVIERE_DE_L'ARTI")==0 || 
	      strcmp(lab2, "SAINT-MARC")==0 || 
	      strcmp(lab2, "VERRETTES")==0) { // Artibonite province
	    pop->setRiverVibrioLevel(i, 20);
	    //	    count1 += pop->infect(rng, i, 0.04, false);
	  } else if (//strcmp(lab2, "LASCAHOBAS")==0 || 
		     strcmp(lab2, "MIREBALAIS")==0) { // Centre province
	    pop->setRiverVibrioLevel(i, 1);
	    //	    count2 += pop->infect(rng, i, 0.007, false);
	  }
	}
      }
      //      cerr << "Infected " << count1 << " in Artibonite." << endl;
      //      cerr << "Infected " << count2 << " in Centre." << endl;
      /*
      int source = pop->infect(rng, 10, 57, 90);
      cerr << "infected 10/" << pop->getNumResidents(source) << " residents in community " << source << endl;
      //    source = pop->infect(rng, 100, 58, 98);
      source = pop->infect(rng, 10, 62, 89);
      cerr << "infected 10/" << pop->getNumResidents(source) << " residents in community " << source << endl;
      */
    } else {
      // R0 test run
      nStartCell = pop->infectOne(rng);
      cerr << "infected one person" << endl;
    }
    //points(x=fromdegminsec(72,28,15), y=fromdegminsec(19,3,15))
    //    source = pop->infect(rng, 30, 102,124);
    //    cerr << "infected 30/" << pop->getNumResidents(source) << " residents in community " << source << endl;
  } else if (bMatlab) { // Matlab
    const int xmax=8, ymax=8; // source data grid size
    const int pops[xmax*ymax] = {0,0,119,6002,1481,0,0,0,
				 675,1753,2969,7570,5705,3891,707,0,
				 7341,9643,4419,6629,10326,6389,3162,0,
				 3506,3755,4301,8017,5421,4588,7135,3421,
				 0,0,2650,1217,6171,5330,6597,2437,
				 0,0,0,377,5414,6608,5717,0,
				 0,0,0,2356,3451,5308,1821,4195,
				 0,0,0,4931,88,233,0,0}; // population of each grid square (from Matlab pop, cholerapop.dat)
    if (communitysize>1) { // matlab on 6km grid
      g = new SquareGridCells(xmax,ymax,false,false);
      double fWorkingFraction=0.35; // fraction of people who work
      //GridCells *g = new HexGridCells(10,10,false,false);

      // initialize grid and non-working residents
      model = new Community[g->getSize()]; // populations of people

      for (int i=0; i<g->getSize(); i++) {
	//	int x = g->getX(i);
	//	int y = g->getY(i);

	int neighbors[101];
	double probs[101];
	// figure out where people in this cell work
	if (fWorkingFraction>0.0) {
	  neighbors[0] = i;
	  int numneighbors1 = getRadius(g, i, 1, neighbors+1, 100);
	  // remove neighbors with no residents
	  for (int j=numneighbors1-1; j>=0; j--) {
	    if (model[neighbors[1+j]].getNumResidents()<=0) {
	      neighbors[1+j] = neighbors[1+numneighbors1-1]; // copy last valid neighbor to this spot
	      numneighbors1--;
	    }
	  }
	  int numneighbors2 = getRadius(g, i, 2, neighbors+numneighbors1+1, 100-numneighbors1);
	  // remove neighbors with no residents
	  for (int j=numneighbors2-1; j>=0; j--) {
	    if (model[neighbors[numneighbors1+1+j]].getNumResidents()<=0) {
	      neighbors[numneighbors1+1+j] = neighbors[numneighbors1+1+numneighbors2-1]; // copy last valid neighbor to this spot
	      numneighbors2--;
	    }
	  }
	  assert(numneighbors1+numneighbors2<100);

	  int sum1=0,
	    sum2=0;
	  for (int j=1; j<numneighbors1+1; j++)
	    sum1+=neighbors[j];
	  for (int j=numneighbors1+1; j<numneighbors2+1; j++)
	    sum2+=neighbors[j];
	  probs[0] = 0.51; // 51% of workers stay in their home region
	  for (int j=0; j<numneighbors1; j++)
	    probs[j] = neighbors[j]*(0.8*(1.0-probs[0]))/sum1;
	  for (int j=numneighbors1; j<numneighbors2; j++)
	    probs[j] = neighbors[j]*(0.2*(1.0-probs[0]))/sum2;
	  //	  model[i].populate(rng, pops[x+y*xmax], fWorkingFraction, neighbors, probs, 1+numneighbors1+numneighbors2);
	} else {
	  //	  model[i].populate(rng, pops[x+y*xmax], 0.0, NULL, NULL, 0);
	}
      }
      // hook up workers to their communities
      if (fWorkingFraction>0.0) {
	for (int i=0; i<=Person::getLastPersonID(); i++) {
	  Person &p = Person::personArray[i];
	  if (p.getHomeCommunity()!=p.getWorkCommunity())
	    model[p.getWorkCommunity()].addVisitor(p.getID());
	}
      }

      // seed infected (0.5% of population)
      for (int i=0; i<g->getSize(); i++) {
	int y = g->getY(i);
	if (y<3) {
	  Community &m = model[i];
	  int temp = m.infect(rng, round(m.getNumResidents()*0.005));
	}
      }
    } else { // Matlab on 1km grid
      const int xmax2=48, ymax2=48; // grid size
      g = new SquareGridCells(xmax2,ymax2,false,false);
      double fWorkingFraction=0.35; // fraction of people who work

      // initialize grid and non-working residents
      model = new Community[g->getSize()]; // populations of people who don't leave their communities (e.g., non-workers)

      for (int i=0; i<g->getSize(); i++) {
	//	cerr << "i=" << i << endl;
	int x = g->getX(i);
	int y = g->getY(i);

	int neighbors[1000];
	double probs[1000];
	int popsize = round(pops[x/6+y/6*xmax]/36.0);
	// figure out where people in this cell work
	if (popsize>0 && fWorkingFraction>0.0) {
	  neighbors[0] = i;
	  probs[0] = 1.0/3.95991;
	  int totalneighbors = 1;
	  for (int radius=1; radius<15; radius++) {
	    int numneighbors = getRadius(g, i, radius, neighbors+totalneighbors, 1000);
	    // remove neighbors with no residents
	    for (int j=numneighbors-1; j>=0; j--) {
	      if (model[neighbors[totalneighbors+j]].getNumResidents()<=0) {
		neighbors[totalneighbors+j] = neighbors[totalneighbors+numneighbors-1]; // copy last valid neighbor to this spot
		numneighbors--;
	      }
	    }

	    for (int j=0; j<numneighbors; j++)
	      probs[totalneighbors+j] = pow(0.75,radius)/3.95991/numneighbors;
	    //	    cerr << radius << " -- " << pow(0.75,radius) << " -- " << (pow(0.75,radius)/3.95991/numneighbors) << " -- " << numneighbors << endl;
	    totalneighbors+=numneighbors;
	    assert(totalneighbors<1000);
	    //	    cerr << "totalneighbors = " << totalneighbors << endl;
	  }
	  //	  model[i].populate(rng, popsize, fWorkingFraction, neighbors, probs, totalneighbors);
	  //	  cerr << round(pops[x/6+y/6*xmax]/36.0) << " residents" << endl;
	} else {
	  //	  model[i].populate(rng, popsize, 0.0, NULL, NULL, 0);
	}
      }
      // hook up workers to their communities
      if (fWorkingFraction>0.0) {
	for (int i=0; i<Person::getLastPersonID(); i++) {
	  Person &p = Person::personArray[i];
	  if (p.getHomeCommunity()!=p.getWorkCommunity()) {
	    assert(p.getWorkCommunity()<g->getSize());
	    assert(p.getWorkCommunity()>=0);
	    model[p.getWorkCommunity()].addVisitor(p.getID());
	  }
	}
      }

      // seed infected (0.5% of population)
      for (int i=0; i<g->getSize(); i++) {
	int y = g->getY(i);
	if (y<3*6) {
	  Community &m = model[i];
	  int temp = m.infect(rng, round(m.getNumResidents()*0.005));
	}
      }
    }
  } else if (nGridSize>0) { // hex grid
    pop=new Population();
    pop->loadGrid(rng,nGridSize,200);
  } else { // one community
    cerr << "one community" << endl;
    pop=new Population();
    pop->loadGrid(rng,1,communitysize);
    int source = pop->infect(rng, 1, 0, 0, true); // put one symptomatic person here
    cerr << "start symptomatic " << (pop->getNumResidents(source)-pop->getNumSusceptible(source)) << " / " << pop->getNumResidents(source) << " residents in community " << source << endl;

    //    g = new SquareGridCells(1,1,false,false);
    //    model = new Community[g->getSize()]; // populations of people who don't leave their communities (e.g., non-workers)
    //    model[0].populate(rng, communitysize);
    //    model[0].infect(rng, 20);    // seed infected
  }

  // pre-vaccinate?
  if (fVaccinateFraction>0.0 && nVaccinateThreshold<0)
    pop->prevaccinate(rng);

  // output tract file
  if (pop && szTractFile.length()>0) {
    ofstream tractFile;
    tractFile.open(szTractFile.c_str());
    if(tractFile.fail()) {
      cerr << "ERROR: Cell file '" << szTractFile << "' cannot be open for writing." << endl;
      return false;
    }
    cerr << "outputing cell information to " << szTractFile << endl;
    tractFile << "id,x,y,pop,river,highway,label,label2" << endl;
    for (int i=0; i < pop->getNumCells(); i++) {
      double popsize=pop->getNumResidents(i);
      int river = pop->getRiver(i);
      int highway = (pop->isHighway(i)?1:0);
      if (popsize>0 || river>0) {
	double x = pop->getX(i);
	double y = pop->getY(i);
	tractFile << i << "," << x << "," << y << "," << popsize << "," << river << "," << highway << "," << pop->getLabel(i) <<  "," << pop->getLabel2(i) << endl;
      }
    }
    tractFile.close();
  }

  // output people file
  if (szPeopleFile.length()>0) {
    cerr << "outputing people information to " << szPeopleFile << endl;
    ofstream peopleFile;
    peopleFile.open(szPeopleFile.c_str());
    if(peopleFile.fail()) {
      cerr << "ERROR: People file '" << szPeopleFile << "' cannot be open for writing." << endl;
      return false;
    }
    peopleFile << "id,home,work" << endl;
    for (int i=0; i <= Person::getLastPersonID(); i++) {
      peopleFile << Person::personArray[i].getID() << "," << Person::personArray[i].getHomeCommunity() << "," << Person::personArray[i].getWorkCommunity() << endl;
    }
    peopleFile.close();
  }

  // set up main output file
  ofstream outputFile;
  if (szOutputFile.length()>0) {
    outputFile.open(szOutputFile.c_str());
    if(outputFile.fail()) {
      cerr << "ERROR: Output file '" << szOutputFile << "' cannot be open for writing." << endl;
      return false;
    }
    cerr << "outputing information to " << szOutputFile << endl;
    outputFile << "time,id,susceptible,infectious,symptomatic,newsymptomatic,cumulativesymptomatic,vibrio,hypervibrio,rivervibrio,riverhypervibrio,vaccines" << endl;
  } else {
    cerr << "not outputing information" << endl;
  }

  // main loop
  for (int t=0; t<runlength; t++) {
    if (t%25==0)
      cerr << "time " << t << endl;
    if (nVaccineFirstDay==t)
      pop->setNumVaccinesAvailable(nVaccineStockpile);
    if (nVaccineFirstDay<=t)
      pop->setNumVaccinesAvailable(pop->getNumVaccinesAvailable()+nVaccinePerDay);
    if (pop->getNumVaccinesAvailable()>nVaccineCap)
      pop->setNumVaccinesAvailable(nVaccineCap);
    if (pop) {
      pop->step(rng);
      if (bPrioritizeRiver) {
	// check to see if river is all vaccinated
	int count = 0;
	for (int i=0; i < pop->getNumCells(); i++)
	  if (pop->wantVaccine(i) && !pop->isVaccinated(i))
	    count++;
	if (count==0) {
	  // done vaccinating river
	  bPrioritizeRiver = false;
	  pop->setVaccinationThreshold(nVaccinateThreshold);
	  pop->setVaccinationDelay(nVaccinateDelay);
	}
      }
    } else {
      cerr << " no pop!" << endl;
      for (int i=0; i < g->getSize(); i++)
	model[i].poop(rng);
      for (int i=0; i < g->getSize(); i++)
	model[i].drink(rng);
      for (int i=0; i < g->getSize(); i++)
	model[i].tick(rng);
    }
    // print day's results
    if (nOutputInterval>0 && t%nOutputInterval==0 && outputFile.is_open()) {
      for (int i=0; i < (pop?pop->getNumCells():g->getSize()); i++) {
	int infectious, s,news,sus,cases,res,vac;
	double vlevel, hvlevel, rvlevel, rhvlevel;
	if (pop) {
	  res=pop->getNumResidents(i);
	  infectious=pop->getNumInfectious(i);
	  s=pop->getNumSymptomatic(i);
	  news=pop->getNumNewSymptomatic(i);
	  sus=pop->getNumSusceptible(i);
	  cases=pop->getCumulativeSymptomatic(i);
	  vac=pop->getNumVaccinated(i);
	  vlevel = pop->getVibrioLevel(i);
	  hvlevel = pop->getHyperVibrioLevel(i);
	  rvlevel = pop->getRiverVibrioLevel(i);
	  rhvlevel = pop->getRiverHyperVibrioLevel(i);
	} else {
	  res=model[i].getNumResidents();
	  infectious=model[i].getNumInfectious();
	  s=model[i].getNumSymptomatic();
	  news=model[i].getNumNewSymptomatic();
	  sus=model[i].getNumSusceptible();
	  cases=model[i].getCumulativeSymptomatic();
	  vac=pop->getNumVaccinated(i);
	  vlevel = model[i].getVibrioLevel();
	  hvlevel = model[i].getHyperVibrioLevel();
	  rvlevel = 0.0;
	  rhvlevel = 0.0;
	}
	if (res>0)
	  outputFile << t << "," << i << "," << sus << "," << infectious << "," << s << "," << news << "," <<  cases << "," << vlevel << "," << hvlevel << "," << rvlevel << "," << rhvlevel << "," << vac << endl;
      }
    }
  }
  if (outputFile.is_open())
    outputFile.close();

  // final output
  // family attack rates
  int family = -1;
  int numcase = 0;
  int numinfected = 0;
  int numpeople = 0;
  int numcasetotal = 0;
  int numinfectedtotal = 0;
  int numpeopletotal = 0;
  int numcasetotalminus = 0;
  int numinfectedtotalminus = 0;
  int numpeopletotalminus = 0;
  int count = 0;
  int countcase = 0;
  for (int i=0; i<Person::getLastPersonID(); i++) {
    Person &person = Person::personArray[i];
    if (person.getFamily()!=family) {
      if (numcase>0) {
	numcasetotal+=numcase;
	numinfectedtotal+=numinfected;
	numpeopletotal+=numpeople;
	numcasetotalminus+=numcase-1;
	numinfectedtotalminus+=numinfected-1;
	numpeopletotalminus+=numpeople-1;
      }
      family = person.getFamily();
      numcase = numinfected = numpeople = 0;
    }
    if (person.getSusceptibility()==0.0) {
      numinfected++;
      count++;
    }
    if (person.isCase()) {
      numcase++;
      countcase++;
    }
    numpeople++;
  }
  cerr << "total cases/infected/pop: " << countcase << "/" << count << "/" << Person::getLastPersonID() << endl;
  cerr << "total family infected with case: " << numcasetotal << "," << numinfectedtotal << "/" << numpeopletotal << endl;
  cerr << "total family infected with case (minus one infected person): " << numinfectedtotalminus << "/" << numpeopletotalminus << endl;

  if (szSummaryFile.length()>0) {
    ofstream summaryFile;
    summaryFile.open(szSummaryFile.c_str());
    if(summaryFile.fail()) {
      cerr << "ERROR: Summary file '" << szSummaryFile << "' cannot be open for writing." << endl;
    } else {
      int nSumCases = 0;
      int nSumRecovered = 0;
      summaryFile << "Started in cell, " << nStartCell << endl;
      summaryFile << "Immune at origin, " << pop->getNumImmune(nStartCell) << endl;
      for (int i=0; i < pop->getNumCells(); i++) {
	nSumCases += pop->getCumulativeSymptomatic(i);
	nSumRecovered += pop->getNumImmune(i);
      }
      summaryFile << "Cases total, " << nSumCases << endl;
      summaryFile << "Immune total, " << nSumRecovered << endl;
      summaryFile.close();
    }
  }
  
  // clean up and exit
  delete [] model;
  delete g;
  delete pop;
  gsl_rng_free(rng);
  return 0;
}
