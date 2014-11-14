/*
 * Population.cpp
 *
 * A network of "cells".  Each cell can have multiple Communities.
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <assert.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "GridCells.h"
#include "Population.h"
#include "R0Community.h"

using namespace std;

Population::Population() {
  _community = NULL;
  _nCommunityStart = NULL;
  _nNumCommunities = 0;
  _grid = NULL;
  _nGridSizeX=0;
  _nGridSizeY=0;
  _riverDist = NULL;
  _riverNeighbors = NULL;
  _newVibrio = NULL;
  _riverVibrio = _riverHyperVibrio = _riverVibrio2 = _riverHyperVibrio2 = NULL;
  _szLabels = _szLabels2 = NULL;
  _fWorkingFraction = 0.0;
  _nXOrigin = _nYOrigin = 0;
  _nRiverShedCycles = 0;
  _fRiverShedFraction = 0.0;
  _fRiverFlowDelta = 0.0;
  _fTravelProb = 0.0;
  _fVaccinationTarget = 0.0;
  _nVaccinationThreshold = 0;
  _nVaccinationDelay = -1;
  _bGridVaccinated = NULL;
  _nVaccinationDay = NULL;
  _nNumVaccinesUsed = 0;
  _nNumVaccinesAvailable = 100000000;
  _fHygieneTarget = 0.0;
  _nTargetGroupSize = 500;
  _nMaxCommunitiesOnRiver = 6;
  _fDriveProb = 0.00001;
  _bHighway = NULL;
  _nHighwayDests = NULL;
  _nNumHighwayDests = NULL;
  _nHighwayDepth = NULL;
  _nDay = -1;
}

Population::~Population() {
  if (_nHighwayDests) {
    for (int i=0; i<_grid->getSize(); i++)
      delete [] _nHighwayDests[i];
    delete [] _nHighwayDests;
  }
  if (_community)
    delete [] _community;
  if (_nCommunityStart)
    delete [] _nCommunityStart;
  if (_grid)
    delete _grid;
  if (_riverDist)
    delete [] _riverDist;
  if (_riverVibrio)
    delete [] _riverVibrio;
  if (_riverHyperVibrio)
    delete [] _riverHyperVibrio;
  if (_riverVibrio2)
    delete [] _riverVibrio2;
  if (_riverHyperVibrio2)
    delete [] _riverHyperVibrio2;
  if (_riverNeighbors)
    delete [] _riverNeighbors;
  if (_newVibrio)
    delete [] _newVibrio;
  if (_bHighway)
    delete [] _bHighway;
  if (_nHighwayDepth)
    delete [] _nHighwayDepth;
  if (_nNumHighwayDests)
    delete [] _nNumHighwayDests;
  if (_bGridVaccinated)
    delete [] _bGridVaccinated;
  if (_nVaccinationDay)
    delete [] _nVaccinationDay;
  if (_szLabels)
    delete [] _szLabels;
  if (_szLabels2)
    delete [] _szLabels2;
}

int Population::loadGrid(gsl_rng *rng, int nGridSize, int nCommunitySize) {
  //  _nGridSize = nGridSize;
  _nGridSizeX = _nGridSizeY = nGridSize;
  _grid = new HexGridCells(_nGridSizeX,_nGridSizeY,false,false);
  _community = new Community[_grid->getSize()]; // populations of people
  _newVibrio = new double[_grid->getSize()];

  // initialize grid and non-working residents
  initGridPopulation(rng, nCommunitySize);

  // seed infected (infect 2 people in each cell on one edge)
  /*  for (int i=0; i<_nGridSizeX; i++) {
    int temp = gsl_ran_binomial(rng, 
				0.005,
				_community[i].getNumResidents());
    _community[i].infect(rng, temp);
    }*/
  return 1;
}

// Initialize a population based on LandScan data
// The population will be on a rectangular grid
// The data from the provided landscan file will be used, starting at
// the specified origin
int Population::loadLandScan(gsl_rng *rng, int nGridSizeX, int nGridSizeY, 
			     double workingfraction,
			     double tau1, double tau2, double rho,
			     const char *landscanfilename,
			     int nOriginXDeg, int nOriginXMin, 
			     int nOriginYDeg, int nOriginYMin,
			     bool bR0run) {
  cerr << "Loading " << landscanfilename << endl;
  _fWorkingFraction = workingfraction;
  ostringstream oss;
  if (landscanfilename)
    oss.str(landscanfilename);
  else
    return -1;

  _nGridSizeX=nGridSizeX;
  _nGridSizeY=nGridSizeY;
  _grid = new SquareGridCells(_nGridSizeX,_nGridSizeY,false,false);
  _nCommunityStart = new int[_grid->getSize()+1]; // index of first community per landscan cell
  _szLabels = new string[_grid->getSize()];
  _szLabels2 = new string[_grid->getSize()];
  _newVibrio = new double[_grid->getSize()];

  for (int i=0; i<=_grid->getSize(); i++)
    _nCommunityStart[i] = 0;
  _nNumCommunities = 0;

  _nXOrigin = convertDegToInt(nOriginXDeg, nOriginXMin, 0);
  _nYOrigin = convertDegToInt(nOriginYDeg, nOriginYMin, 0);

  // count communities and people in LandScan data
  ifstream iss(oss.str().c_str());
  if (!iss) {
    cerr << "ERROR: " << landscanfilename << " not found." << endl;
    return -1;
  }
  string line;
  getline(iss, line); // throw away header line
  while (getline(iss, line)) {
    istringstream iss;
    iss.str(line);
    int pop,xdeg,xmin,xsec,ydeg,ymin,ysec;
    string xdir, ydir, commune, department;
    iss >> pop >> xdeg >> xmin >> xsec >> xdir >> ydeg >> ymin >> ysec >> ydir >> commune >> department;
    int xloc = convertDegToInt(xdeg, xmin, xsec)-_nXOrigin;
    int yloc = convertDegToInt(ydeg, ymin, ysec)-_nYOrigin;
    if (pop>0 &&
	xloc>=0 && xloc<_nGridSizeX &&
	yloc>=0 && yloc<_nGridSizeY) {
      int gridindex = (int)(round(yloc*_nGridSizeX+xloc)); // assumes square lattice - should be fixed!
      assert (gridindex<_grid->getSize());
      _nCommunityStart[gridindex]=round(pop/_nTargetGroupSize);
      if (_nCommunityStart[gridindex]==0)
	_nCommunityStart[gridindex] = 1;
      _nNumCommunities+=_nCommunityStart[gridindex];
    }
  }
  iss.close();
  for (int i=_grid->getSize()-1; i>=0; i--)
    _nCommunityStart[i+1] = _nCommunityStart[i];
  _nCommunityStart[0] = 0;
  for (int i=1; i<=_grid->getSize(); i++)
    _nCommunityStart[i] += _nCommunityStart[i-1]; // cumulative
  cerr <<       _nNumCommunities << " communities" << endl;

  // load LandScan data and populate communities
  if (bR0run)
    _community = new R0Community[_nNumCommunities]; // well-mixed groups of people
  else 
    _community = new Community[_nNumCommunities]; // well-mixed groups of people
  iss.open(oss.str().c_str());
  if (!iss) {
    cerr << "ERROR: " << landscanfilename << " not found." << endl;
    return -1;
  }
  getline(iss, line); // throw away header line
  while (getline(iss, line)) {
    istringstream iss;
    iss.str(line);
    int pop,xdeg,xmin,xsec,ydeg,ymin,ysec;
    string xdir, ydir, commune, department;
    iss >> pop >> xdeg >> xmin >> xsec >> xdir >> ydeg >> ymin >> ysec >> ydir >> commune >> department;
    int xloc = convertDegToInt(xdeg, xmin, xsec)-_nXOrigin;
    int yloc = convertDegToInt(ydeg, ymin, ysec)-_nYOrigin;

    //    cout << (xloc-_nXOrigin) << ": " << xdeg << "," << xmin  << "," << xsec  << ".  " << (yloc-nOriginYLoc) << ": " << ydeg << "," << ymin  << "," << ysec << endl;
    if (pop>0 &&
	xloc>=0 && xloc<_nGridSizeX &&
	yloc>=0 && yloc<_nGridSizeY) {
      //     cerr << "x,y loc: " << xloc << "," << yloc << "," << pop << endl;
      int gridindex = (int)(round(yloc*_nGridSizeX+xloc)); // assumes square lattice - should be fixed!
      assert (gridindex<_grid->getSize());
      _szLabels[gridindex] = department;
      _szLabels2[gridindex] = commune;

      // divide cell population among communities
      int numcomms = _nCommunityStart[gridindex+1]-_nCommunityStart[gridindex];
      int popleft = pop;
      for (int i=0; i<numcomms; i++) {
	//	cout << "grid id = " << gridindex << " : " << _nCommunityStart[gridindex]+i << "(" << i << ") : " << round(popleft/(numcomms-i)) << "/" << pop << endl;
	Community &comm = _community[_nCommunityStart[gridindex]+i];
	comm.populate(rng, round(popleft/(numcomms-i)));
	popleft-=round(popleft/(numcomms-i));
      }
    }
  }
  iss.close();

  _bGridVaccinated = new bool[_grid->getSize()];
  for (int i=0; i<_grid->getSize(); i++)
    _bGridVaccinated[i] = false;
  _nVaccinationDay = new int[_grid->getSize()];
  for (int i=0; i<_grid->getSize(); i++)
    _nVaccinationDay[i] = 10000000;

  // assign people workplaces
  for (int sourceindex=0; sourceindex<_grid->getSize(); sourceindex++)
    assignWorkPlaces(rng, sourceindex, _fWorkingFraction,
		     tau1, tau2, rho);
  return 1;
}

// assignWorkPlaces - assigns workplaces for residents of cell nSourceID
// returns number of employed people
int Population::assignWorkPlaces(gsl_rng *rng, int nSourceID, 
				 double fWorkFraction,
				 double tau1, double tau2, double rho) {
  const int maxradius = 40;
  if (fWorkFraction<=0.0)
    return 0;
  int P1 = getNumResidents(nSourceID); // population in this cell
  if (P1==0)
    return 0;
  double sourcex = _grid->getX(nSourceID);
  double sourcey = _grid->getY(nSourceID);

  // calculate probabilities of working in other neighborhoods
  int neighbors[maxradius*maxradius*4]; // id of neighbors
  double probs[maxradius*maxradius*4];  // probabilities based on gravity model
  int totalneighbors = 1;
  for (int destindex=0; destindex<_grid->getSize(); destindex++) {
    double destx = _grid->getX(destindex);
    double desty = _grid->getY(destindex);
    if (destindex!=nSourceID && fabs(destx-sourcex)<maxradius && fabs(desty-sourcey)<maxradius) {
      double dist = sqrt((sourcex-destx)*(sourcex-destx) + (sourcey-desty)*(sourcey-desty));
      int P2 = getNumResidents(destindex);
      if (P2>0 && dist<maxradius) {
	neighbors[totalneighbors] = destindex;
	probs[totalneighbors] = pow(P1, tau1) * pow(P2, tau2) / pow(dist, rho); // do we need theta * ???
	//	      cerr << "dist=" << dist << "," << probs[totalneighbors] <<  "," << 		pow(P1, tau1)  << ","  << pow(P2, tau2) << ","  << (pow(dist, rho)) << endl;
	totalneighbors++;
	assert(totalneighbors<maxradius*maxradius*4);
      }
    }
  }
  if (totalneighbors==1) {
    return 0;
  } else {
    int count = 0;
    neighbors[0] = nSourceID;
    probs[0] = 0.3; // 30% chance of working in home community
    double sum=0.0;
    for (int i=1; i<totalneighbors; i++)
      sum += probs[i];
    sum/=(1.0-probs[0]); // normalize over neighbors >=1 cell away
    for (int i=1; i<totalneighbors; i++)
      probs[i] = probs[i-1]+probs[i]/sum;
    assert(fabs(1.0-probs[totalneighbors-1])<0.000001);
	   
    for (int commindex=_nCommunityStart[nSourceID]; commindex<_nCommunityStart[nSourceID+1]; commindex++) {
      Community &sourcecomm = _community[commindex];
      for (int i=0; i<sourcecomm.getNumResidents(); i++) {
	int id = sourcecomm.getResidentID(i);
	Person &person = Person::personArray[id];
	assert(person.getHomeCommunity()==commindex);
	if (gsl_rng_uniform(rng)<fWorkFraction) { // this person works
	  double r = gsl_rng_uniform(rng);
	  for (int j=0; j<totalneighbors; j++)
	    if (probs[j]>=r) {
	      int dest = _nCommunityStart[neighbors[j]];
	      int diff = _nCommunityStart[neighbors[j]+1]-_nCommunityStart[neighbors[j]];
	      if (diff>1)
		dest += gsl_rng_uniform_int(rng, diff);
	      person.setWorkCommunity(dest);
	      _community[dest].addVisitor(id);
	      count++;
	      break;
	    }
	}
      }
    }
    return count;
  }
}

int Population::loadRivers(const char *riverfilename) {
  ostringstream oss;
  if (riverfilename)
    oss.str(riverfilename);
  else
    return -1;
  ifstream iss(oss.str().c_str());
  if (!iss) {
    cerr << "ERROR: " << riverfilename << " not found." << endl;
    return -1;
  }
  _riverDist = new int[_grid->getSize()];
  _riverNeighbors = new int[_grid->getSize()];
  _riverVibrio = new double[_grid->getSize()];
  _riverHyperVibrio = new double[_grid->getSize()];
  _riverVibrio2 = new double[_grid->getSize()];
  _riverHyperVibrio2 = new double[_grid->getSize()];
  for (int i=0; i<_grid->getSize(); i++) {
    _riverDist[i] = -1;
    _riverVibrio[i] = 0.0;
    _riverHyperVibrio[i] = 0.0;
  }

  cerr << "loading rivers " << _nGridSizeX << " x " << _nGridSizeY << endl;
  string line;
  getline(iss, line); // throw away header line
  int nNumLinesRead = 0;
  while (getline(iss, line)) {
    istringstream iss;
    iss.str(line);
    int xdeg,xmin,xsec,ydeg,ymin,ysec,dist;
    iss >> ydeg >> ymin >> ysec >> xdeg >> xmin >> xsec >> dist;
    nNumLinesRead++;
    int xloc = convertDegToInt(xdeg, xmin, xsec)-_nXOrigin;
    int yloc = convertDegToInt(ydeg, ymin, ysec)-_nYOrigin;
    if (dist>0 &&
	xloc>=0 && xloc<_nGridSizeX &&
	yloc>=0 && yloc<_nGridSizeY) {
      //     cerr << "river x,y loc: " << xloc << "," << yloc << "," << dist << endl;
     int gridindex = xloc + yloc*_nGridSizeX; // assumes rectangular lattice. should be fixed!
     assert (gridindex<_grid->getSize());
     _riverDist[gridindex] = dist;
    for (int commindex=_nCommunityStart[gridindex]; commindex<_nCommunityStart[gridindex+1]; commindex++)
      _community[commindex].setRiver(true); // !!!!!!!!!!!!!1
    }
  }
  // count cells that are downstream from each cell
  for (int i=0; i < _grid->getSize(); i++) {
    _riverNeighbors[i] = 0;
    if (_riverDist[i]>0) {
      for (int x=-1; x<=1; x++)
	for (int y=-1; y<=1; y++)
	  if ((x!=0 || y!=0) &&
	      _riverDist[i+x+y*_nGridSizeX] > _riverDist[i])
	    _riverNeighbors[i]++;
    }
  }

  return nNumLinesRead;
}

int Population::loadHighways(const char *highwayfilename) {
  ostringstream oss;
  if (highwayfilename)
    oss.str(highwayfilename);
  else
    return -1;
  ifstream iss(oss.str().c_str());
  if (!iss) {
    cerr << "ERROR: " << highwayfilename << " not found." << endl;
    return -1;
  }
  _bHighway = new bool[_grid->getSize()];
  for (int i=0; i<_grid->getSize(); i++)
    _bHighway[i] = false;

  cerr << "loading highways " << _nGridSizeX << " x " << _nGridSizeY << endl;
  string line;
  getline(iss, line); // throw away header line
  int nNumLinesRead = 0;
  while (getline(iss, line)) {
    istringstream iss;
    iss.str(line);
    int xdeg,xmin,xsec,ydeg,ymin,ysec,dist;
    iss >> ydeg >> ymin >> ysec >> xdeg >> xmin >> xsec >> dist;
    nNumLinesRead++;
    int xloc = convertDegToInt(xdeg, xmin, xsec)-_nXOrigin;
    int yloc = convertDegToInt(ydeg, ymin, ysec)-_nYOrigin;
    if (dist>0 &&
	xloc>=0 && xloc<_nGridSizeX &&
	yloc>=0 && yloc<_nGridSizeY) {
      int gridindex = xloc + yloc*_nGridSizeX; // assumes rectangular lattice. should be fixed!
      assert (gridindex<_grid->getSize());
      _bHighway[gridindex] = true;
    }
  }
  return nNumLinesRead;
}

// population is on a uniform lattice
int Population::initGridPopulation(gsl_rng *rng, int nCommunitySize) {
  // initialize grid and non-working residents
  int maxradius = min(15,_nGridSizeX);
  double nnorm = 0.0;
  double fade = 0.75;
  for (int i=0; i<maxradius; i++)
    nnorm += pow(fade,i);

  for (int i=0; i<_grid->getSize(); i++) {
    int neighbors[1000];
    double probs[1000];
    int pop = nCommunitySize;
    // figure out where people in this cell work
    if (pop>0 && _fWorkingFraction>0.0) {
      neighbors[0] = i;
      probs[0] = 1.0/nnorm;
      int totalneighbors = 1;
      for (int radius=1; radius<15; radius++) {
	int numneighbors = getRadius(_grid, i, radius, neighbors+totalneighbors, 1000);
	// remove neighbors with no residents
	for (int j=numneighbors-1; j>=0; j--) {
	  if (getNumResidents(neighbors[totalneighbors+j])<=0) {
	    neighbors[totalneighbors+j] = neighbors[totalneighbors+numneighbors-1]; // copy last valid neighbor to this spot
	    numneighbors--;
	  }
	}
	for (int j=0; j<numneighbors; j++)
	  probs[totalneighbors+j] = pow(fade,radius)/nnorm/numneighbors;
	totalneighbors+=numneighbors;
	assert(totalneighbors<1000);
      }
      for (int j=1; j<totalneighbors; j++)
	probs[j] = probs[j-1]+probs[j];
      _community[i].populate(rng, pop);
    } else {
      _community[i].populate(rng, pop);
    }
  }

  for (int sourceindex=0; sourceindex<_grid->getSize(); sourceindex++)
    assignWorkPlaces(rng, sourceindex, _fWorkingFraction,
		     0.0,0.0,2); //tau1, tau2, rho); // fix this!!!!!!!

  return 1;
}

int Population::getRadius(GridCells *g, int center, int radius, int *buf, int bufsize) {
  const int maxneighbors = 6;
  int searched[4000];
  int searchedsize=1;
  int oldsearchedsize=1;
  searched[0] = center;
  // breadth first search
  for (int r=0; r<radius; r++) {
    oldsearchedsize=searchedsize;
    for (int i=searchedsize-1; i>=0; i--) {
      int neighbors[maxneighbors];
      int numneighbors = _grid->getNeighbors(searched[i], neighbors, maxneighbors);
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

// each degree is broken into 120 intervals
int Population::convertDegToInt(int deg, int min, int sec) {
  return deg*120+min*2+(sec>30?1:0);
}

// infect - infect frac of the people in the cell gridindex
int Population::infect(gsl_rng *rng, int gridindex, double frac, bool bMakeSymptomatic) {
  int count = 0;
  for (int commindex=_nCommunityStart[gridindex];
       commindex < _nCommunityStart[gridindex+1];
       commindex++) {
	 if (_community[commindex].getNumResidents()>0) {
	int nNumInfect = gsl_ran_binomial(rng, 
					   frac,
					  _community[commindex].getNumResidents());
	count += _community[commindex].infect(rng, nNumInfect, bMakeSymptomatic);
	 }
       }
  return count;
}

// infect - infect num people in the cell closest to x,y
int Population::infect(gsl_rng *rng, int num, double x, double y, bool bMakeSymptomatic) {
  int index=0;
  double xindex, yindex;
  double distindex=1000000.0;
  for (int i=0; i<_grid->getSize(); i++) {
    double newx = _grid->getX(i);
    double newy = _grid->getY(i);
    double newdist = (newx-x)*(newx-x)+(newy-y)*(newy-y);
    if (newdist<distindex && getNumResidents(i)>0) {
      index = i;
      xindex = newx;
      yindex = newy;
      distindex = newdist;
    }
  }
  if (_community[_nCommunityStart[index]].getNumResidents()>0)
    _community[_nCommunityStart[index]].infect(rng, num, bMakeSymptomatic);
  return index;
}

// infectOne - infects a single person in a random community
int Population::infectOne(gsl_rng *rng, bool bMakeSymptomatic) {
  int commindex = gsl_rng_uniform_int(rng, _nNumCommunities);
  _community[commindex].infect(rng, 1, bMakeSymptomatic);
  for (int i=0; i<_grid->getSize(); i++)
    if (commindex<_nCommunityStart[i+1])
      return i;
  return -1;
}

// prevaccinate
// we do not check how much vaccine is available
void Population::prevaccinate(gsl_rng *rng) {
  if (_fVaccinationTarget>0.0 && _bGridVaccinated) {
    for (int i=0; i < _grid->getSize(); i++)
      if (!_bGridVaccinated[i]) {
	for (int j=_nCommunityStart[i]; j<_nCommunityStart[i+1]; j++)
	  _nNumVaccinesUsed += _community[j].prevaccinate(rng, _fVaccinationTarget);
	_bGridVaccinated[i] = true;
	_nVaccinationDay[i] = -1;
      }
  }
}

void Population::prioritizeCell(int gridindex) {
  if (!_bGridVaccinated[gridindex])
    _nVaccinationDay[gridindex] = _nDay;
}

int Population::computedrivelength(gsl_rng *rng, int start) {
  const int MAXDIST=200;
  const int MAXQUEUESIZE=MAXDIST*30;
  int gridqueue[MAXQUEUESIZE+1];
  if (_nHighwayDepth==NULL)
    _nHighwayDepth = new int[_grid->getSize()];
  for (int i=0; i<_grid->getSize(); i++)
    _nHighwayDepth[i] = 0;
  int head = 0;
  int tail = 1;
  gridqueue[head] = start;
  _nHighwayDepth[start] = 1;
  //  cerr << "start at " << start << ": " << (start%_nGridSizeX) << "," << (start/_nGridSizeX) << endl;
  while(head<tail && tail<MAXQUEUESIZE-1) {
    int loc = gridqueue[head++];
    for (int x=-1; x<=1; x++)
      for (int y=-1; y<=1; y++) {
	int pos = loc+_nGridSizeX*y+x;
	if (pos>=0 && pos<_grid->getSize())
	  if ((x!=0 || y!=0) &&
	      _bHighway[pos] &&
	      _nHighwayDepth[pos]<=0) {
	    //	  cerr << "  " << loc << ": " << (loc%_nGridSizeX)+x << "," << (loc/_nGridSizeX)+y << "," << (_nHighwayDepth[loc]+1) << endl;
	    if (_nHighwayDepth[loc]<MAXDIST) {
	      assert(tail<MAXQUEUESIZE);
	      if (tail>=MAXQUEUESIZE)
		break;
	      gridqueue[tail++] = loc+_nGridSizeX*y+x;
	      _nHighwayDepth[pos] = _nHighwayDepth[loc]+1;
	    }
	  }
      }
  }
  return start;
}

int Population::drive(gsl_rng *rng) {
  if (_fDriveProb>0.0)
    for (int i=0; i < _grid->getSize(); i++) {
      if (_bHighway[i] && (_nNumHighwayDests==NULL || _nNumHighwayDests[i]!=0) && getNumResidents(i)>0) {
	int nNumDrivers = gsl_ran_binomial(rng, 
					   _fDriveProb,
					   getNumResidents(i));
	if (nNumDrivers>0) {
	  if (!_nNumHighwayDests) {
	    _nNumHighwayDests = new int[_grid->getSize()];
	    for (int j=0; j < _grid->getSize(); j++)
	      _nNumHighwayDests[j] = -1;
	    _nHighwayDests = new int *[_grid->getSize()];
	  }
	  if (_nNumHighwayDests[i]<0) {
	    // store possible destinations
	    _nNumHighwayDests[i] = 0;
	    computedrivelength(rng, i);
	    _nHighwayDests[i] = new int[500];
	    for (int j=0; j<_grid->getSize(); j++)
	      if (_nHighwayDepth[j]>30 && getNumResidents(j)>2000) {
		_nHighwayDests[i][_nNumHighwayDests[i]++] = j;
		assert(_nNumHighwayDests[i]<500);
	      }
	  }
	  if (_nNumHighwayDests[i]>0 && _nHighwayDests[i])
	    while (--nNumDrivers >= 0) {
	      // swap non-symptomatic person at source with dest
	      int destnum = gsl_rng_uniform_int(rng, _nNumHighwayDests[i]);
	      int commsource = _nCommunityStart[i];
	      if (_nCommunityStart[i+1]>_nCommunityStart[i]+1)
		commsource += gsl_rng_uniform_int(rng, _nCommunityStart[i+1]-_nCommunityStart[i]);
	      int commdest = _nCommunityStart[_nHighwayDests[i][destnum]];
	      if (_nCommunityStart[_nHighwayDests[i][destnum]+1]>_nCommunityStart[_nHighwayDests[i][destnum]]+1)
		commdest += gsl_rng_uniform_int(rng, _nCommunityStart[_nHighwayDests[i][destnum]+1]-_nCommunityStart[_nHighwayDests[i][destnum]]);
	      int person1id = _community[commsource].getResidentID(gsl_rng_uniform_int(rng, _community[commsource].getNumResidents()));
	      int person2id = _community[commdest].getResidentID(gsl_rng_uniform_int(rng, _community[commdest].getNumResidents()));
	      Person &p1 = Person::personArray[person1id];
	      Person &p2 = Person::personArray[person2id];
	      if ((!p1.isSymptomatic() || p1.getDaysInfectious()==0) &&
		  (!p2.isSymptomatic() || p2.getDaysInfectious()==0)) { // don't travel if symptomatic for more than one day
		Person temp;
		temp.copyInfectionStatus(p1);
		p1.copyInfectionStatus(p2);
		p2.copyInfectionStatus(temp);
	      }
	      //	    cout << (start%_nGridSizeX) << "," << (start/_nGridSizeX) << "," << _nHighwayDepth[start] << endl;
	    }
	}
      }
    }
  return 0;
}

int Population::step(gsl_rng *rng) {
  _nDay++;
  // shed vibrio and record amount shed today
  for (int i=0; i < _grid->getSize(); i++) {
    _newVibrio[i] = 0;
    for (int j=0; j<_nCommunityStart[i+1]-_nCommunityStart[i]; j++)
      if (j<_nMaxCommunitiesOnRiver)
	_newVibrio[i] += _community[_nCommunityStart[i]+j].poop(rng);
      else 
	_community[_nCommunityStart[i]+j].poop(rng);
  }
  if (_riverDist && _fRiverShedFraction>0.0 && _nRiverShedCycles>0) {
    // decay of vibrio in rivers
    for (int i=0; i < _grid->getSize(); i++)
      if (_riverDist[i]>0) {
	_riverVibrio[i]*=Community::getVibrioDecay();     // level of Vibrio goes down
	_riverHyperVibrio[i] = 0.0;          // hyperinfectious vibrio disappears each day
      }
    // send vibrio down river (diffusion)
    for (int cycle=0; cycle<_nRiverShedCycles; cycle++) {
      for (int i=0; i < _grid->getSize(); i++)
	if (_riverDist[i]>0) {
	  _riverVibrio2[i] = 0.0;
	  _riverHyperVibrio2[i] = 0.0;
	}
      // copy vibrio to one cell downstream
      for (int i=0; i < _grid->getSize(); i++) 
	if (_riverDist[i]>0 &&
	    (_riverVibrio[i]>0.0 || _riverHyperVibrio[i]>0.0) &&
	    _riverNeighbors[i]>0) 
	  for (int x=-1; x<=1; x++)
	    for (int y=-1; y<=1; y++) {
	      int index = i+x+y*_nGridSizeX;
	      if ((x!=0 || y!=0) &&
		  _riverDist[index] > _riverDist[i]) {
		_riverVibrio2[index] += (1.0-_fRiverFlowDelta) * _riverVibrio[i]/_riverNeighbors[i];
		_riverHyperVibrio2[index] += (1.0-_fRiverFlowDelta) * _riverHyperVibrio[i]/_riverNeighbors[i];
	      }
	    }
      // add new vibrio from original sources
      for (int i=0; i < _grid->getSize(); i++)
	if (_riverDist[i]>0) {
	  double temp = _newVibrio[i]*_fRiverShedFraction/(_nRiverShedCycles+1);
	  _riverVibrio2[i] += temp;
	  _riverHyperVibrio2[i] += temp;
	}
      double *temp = _riverVibrio;
      _riverVibrio = _riverVibrio2;
      _riverVibrio2 = temp;
      temp = _riverHyperVibrio;
      _riverHyperVibrio = _riverHyperVibrio2;
      _riverHyperVibrio2 = temp;
    }
    // update river vibrio levels in communities
    for (int i=0; i < _grid->getSize(); i++) {
      if (_riverDist[i]>0) {
	int diff = _nCommunityStart[i+1]-_nCommunityStart[i];
	if (diff>_nMaxCommunitiesOnRiver)
	  diff=_nMaxCommunitiesOnRiver;
	for (int j=0; j<diff; j++) {
	  _community[_nCommunityStart[i]+j].setRiverVibrioLevel(_riverVibrio[i]);
	  _community[_nCommunityStart[i]+j].setRiverHyperVibrioLevel(_riverHyperVibrio[i]);
	}
      }
    }
  }

  // infect susceptibles
  for (int i=0; i < _nNumCommunities; i++)
    _community[i].drink(rng);
  for (int i=0; i < _nNumCommunities; i++)
    _community[i].tick(rng);

  // reactive vaccination if appropriate
  if (_fVaccinationTarget>0.0) {
    if (_nVaccinationThreshold>=0 && _nVaccinationDelay>=0) {
      // trigger (schedule) reactive vaccination
      for (int i=0; i < _grid->getSize(); i++) {
	if (!_bGridVaccinated[i] &&
	    _nVaccinationDay[i]>_nDay+_nVaccinationDelay && 
	    getCumulativeSymptomatic(i)>_nVaccinationThreshold) {
	  _nVaccinationDay[i] = _nDay+_nVaccinationDelay;
	}
      }
    }

    // vaccinate if ready
    if (_nNumVaccinesAvailable>_nNumVaccinesUsed) {
      int nNumPeopleWant = 0;
      int nNumCellsWant = 0;
      for (int i=0; i < _grid->getSize(); i++)
	if (!_bGridVaccinated[i] &&
	    _nVaccinationDay[i]<=_nDay &&
	    getNumResidents(i)>0) {
	  nNumPeopleWant+=getNumResidents(i)*_fVaccinationTarget;
	  nNumCellsWant++;
	}
      if (nNumPeopleWant<=_nNumVaccinesAvailable-_nNumVaccinesUsed) {
	// enough vaccine for everyone who wants it
	for (int i=0; i < _grid->getSize(); i++) {
	  if (!_bGridVaccinated[i] &&
	      _nVaccinationDay[i]<=_nDay &&
	      _nNumVaccinesAvailable>_nNumVaccinesUsed) {
	    for (int j=_nCommunityStart[i]; j<_nCommunityStart[i+1]; j++) {
	      _nNumVaccinesUsed += _community[j].vaccinate(rng, _fVaccinationTarget);
	      if (_fHygieneTarget>0.0)
		_community[j].setHygiene(_fHygieneTarget);
	    }
	    _bGridVaccinated[i] = true;
	  }
	}
      } else if (nNumPeopleWant>0) {
	// not enough vaccine - randomly choose cells to vaccinate until we run out
	int miss = 0;
	while (_nNumVaccinesAvailable>_nNumVaccinesUsed && nNumCellsWant>0 && miss<5) {
	  int gridcount = 0;
	  int gridid = 0;
	  gridcount = gsl_rng_uniform_int(rng, nNumCellsWant); // choose the nth unvaccinated cell
	  for (gridid=0; gridid<_grid->getSize() && gridcount>=0; gridid++)
	    if (!_bGridVaccinated[gridid] &&
		_nVaccinationDay[gridid]<=_nDay &&
		getNumResidents(gridid)>0)
	      gridcount--;

	  if (getNumResidents(gridid)*_fVaccinationTarget*0.75>_nNumVaccinesAvailable-_nNumVaccinesUsed) {
	    // not enough vaccine for this cell. try again
	    miss++;
	  } else {
	    // vaccinate this cell
	    for (int j=_nCommunityStart[gridid]; j<_nCommunityStart[gridid+1]; j++)
	      _nNumVaccinesUsed += _community[j].vaccinate(rng, _fVaccinationTarget);
	    _bGridVaccinated[gridid] = true;
	    nNumCellsWant--;
	  }
	}
      }
    }
  }
 
  // some people drive on the highways
  drive(rng);

  // "travel" by swapping infection and vaccination status of 2 random (healthy) people
  if (_fTravelProb>0.0) {
    int nNumTravelers = gsl_ran_binomial(rng, 
					 _fTravelProb,
					 Person::getLastPersonID()+1);
    for (int i=0; i<nNumTravelers; i++) {
      // should we check for symptomatics?????
      int person1id = gsl_rng_uniform_int(rng, Person::getLastPersonID()+1);
      int person2id = gsl_rng_uniform_int(rng, Person::getLastPersonID()+1);
      Person &p1 = Person::personArray[person1id];
      Person &p2 = Person::personArray[person2id];
      if ((!p1.isSymptomatic() || p1.getDaysInfectious()==0) &&
	  (!p2.isSymptomatic() || p2.getDaysInfectious()==0)) { // don't travel if symptomatic for more than one day
	Person temp;
	temp.copyInfectionStatus(p1);
	p1.copyInfectionStatus(p2);
	p2.copyInfectionStatus(temp);
      }
    }
  }
  return 1;
}
