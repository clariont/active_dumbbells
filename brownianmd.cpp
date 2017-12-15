

///////////////////////////////////////////////////////////////////////////////////////////////////////
//  Date:   14 May 2014
//  Description:    Performs Brownian MD
//
//  Usage Syntax:
//
///////////////////////////////////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <string>
#include "genarray.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;

// Global Variables
const double PI = 3.14159265359;
gsl_rng * mrRand;
double mySeed = 1;
int natoms = 20;
double volfrac = 0.3;
double boxrad = sqrt(natoms/4.0/volfrac);
double TEMP = 1;
double actforce = 0;
double r_cut = 2.5;
double skin = 2.0;
double mass = 1.0;
double myGammaPar = 1.0;
double myGammaPerp = 2*myGammaPar;
double myGammaRot = 2/3.0*myGammaPar;
int totalSteps = 10000000;
int writeEvery = 1000;
double timeStep = 0.0001;
bool yesRods;
int yesRestart;
string restartName;
int yesEquil;


// Potential Parameters
const double ljEps = 2.5;
const double ljSig = 1.0;

const double irSig = 1.0;
const double irN = 12;
const double irA = 0.018;
const double irRa = irSig;
const double irEa = 1.0;
const double irRr = 2*irSig;
const double irEr = 1.0;

double ir_Temp = 0;		    // Following the Reatto convention - their temperature is in the range of [0,1.0]
double ir_invTemp = 0;		    // This scales the potential so that we can set T = 1.0.

double rodSpacing = irSig;



// Function Declarations
double distsq(double p1x, double p1y, double p2x, double p2y);
    //
double calcEnergy_lj(genarray<double> &atomPositions, genarray<int> &verlNums, genarray<int> &verlLists); 
double calcEnergy_lj_noverl(genarray<double> &atomPositions, genarray<int> &verlNums, genarray<int> &verlLists); 
    //
double calcForces_lj(genarray<double> &atomPositions, genarray<int> &verlNums, genarray<int> &verlLists, genarray<double> &forces); 
    //
void integrate(genarray<double> &atomPositions, genarray<int> &verlNums, genarray<int> &verlLists, genarray<double> &verlPos, genarray<double> &forces, double dt);
    //
void setVerlet(genarray<double> &atomPositions, genarray<int> &verlNums, genarray<int> &verlLists, genarray<double> &verlPos);
    //
void checkVerlet(genarray<double> &atomPositions, genarray<int> &verlNums, genarray<int> &verlLists, genarray<double> &verlPos); 
    // 
void initConfig(genarray<double> &atomPositions);
    //
void writeConfig(genarray<double> &atomPositions, string outName);
    //

double calcEnergy_ir(genarray<double> &atomPositions, genarray<int> &verlNums, genarray<int> &verlLists); 
    //
double calcForces_ir(genarray<double> &atomPositions, genarray<int> &verlNums, genarray<int> &verlLists, genarray<double> &forces); 
    //
void calcRodForces(genarray<double> &rodList, genarray<double> &rodForces, genarray<double> &forces);
    //
void paramReader (string fileName);
    //
void integrateRods(genarray<double> &atomPositions, genarray<int> &verlNums, genarray<int> &verlLists, genarray<double> &verlPos, genarray<double> &forces, genarray<double> &rodForces, double dt);
    //
int mySign(double x);
    //
double calcForces_ljrep(genarray<double> &atomPositions, genarray<int> &verlNums, genarray<int> &verlLists, genarray<double> &forces); 

void readRestart (string restartName, genarray<double> &atomPositions);

bool checkExplosion(genarray<double> &atomPositions); 

// Main
int main(int argv, char *argc[]) {

    // Read parameter file
    string paramFile(argc[1]);
    cout << "reading parameters from: " << paramFile << endl;
    paramReader(paramFile);
    boxrad = sqrt(natoms/4.0/volfrac);

    // Seed
    mrRand = gsl_rng_alloc (gsl_rng_taus);
    gsl_rng_set(mrRand, 1);

    // Initialize variables
    genarray<double> atomPositions(natoms*3);	    // Elements: 0 - x, 1 - y, 2 - theta
    genarray<double> atomForces(natoms*2);	    // Elements: 0 - Fx, 1 - Fy

    genarray<int> verlNums(natoms);			    // In the future we can make a verlet class to encapsulate these arrays for cleanliness
						    // - or, wrap all the arrays in a whole class...
    genarray<double> verlPos(natoms*2);			    // Elements: 0 - x, 1 - y
    genarray<int> verlLists(natoms*natoms);			    // Elements: 0..(natoms-1) is the verlet list for the first particle, natoms..2*natoms-1 is
						    // the list for the second particle, etc...


    cout << "the temperature is: " << TEMP << endl;
    if (TEMP < 0.5) cout << "it is chilly.\n";
    else if (TEMP >= 0.5 and TEMP < 1.5) cout << "it is nice out.\n";
    else cout << "it's hot!\n";
    string movieName("hello.lammpstrj");
    ofstream myOut;
    myOut.open(movieName.c_str(), ios::out);
    myOut.close();

    if (yesRestart == 1) { readRestart(restartName, atomPositions); }
    else { initConfig(atomPositions); }

    writeConfig(atomPositions, movieName);
    setVerlet(atomPositions, verlNums, verlLists, verlPos);
    

    int equilSteps = 50000;
    double dt = timeStep;

    r_cut = 10.0*irSig;

    cout << "volume fraction: " << natoms/4.0/boxrad/boxrad << endl;
    cout << "box length: " << 2*boxrad << endl;

    // Calculate U0:
    //	(this allows us to rescale the potential and set temperature=1.0)
    double u0 = irA*pow(irSig, irN) - irEa*irSig*irSig/irRa/irRa*exp(-1.0/irRa) + irEr*irSig*irSig/irRr/irRr*exp(-1.0/irRr);
    u0 = fabs(u0);
    cout << "U_0 is: " << u0 << endl;	    // the value of the potential at sigma.
    ir_invTemp = 1/ir_Temp/u0;		    // this factor rescales the potential.

    // Print out the well depth (to get an idea of the active force magnitude )
    double wellMin = 0.974936140452;	    // this is the value of 'r' at which the I-R potential has a minimum.
    double myWellDepth = ir_invTemp*( irA*pow(irSig/wellMin, irN) - irEa*irSig*irSig/irRa/irRa*exp(-wellMin/irRa) + irEr*irSig*irSig/irRr/irRr*exp(-wellMin/irRr) );
    cout << "Well depth of the potential is: " << myWellDepth << endl;
    cout << "starting energy: " << calcEnergy_ir(atomPositions, verlNums, verlLists) << endl;


    // Run
    TEMP = 1.0;
    r_cut = 10.0*irSig;
    dt = timeStep;


    genarray<double> dummyList(1);
    genarray<double> rodForces;
    rodForces.resize(int(natoms/2.0)*2);
    
    double saveinv = ir_invTemp;
    double saveFORCE = actforce;

    // Equilibration:
    if (yesEquil == 1) {
	double saveTEMP = TEMP;
	TEMP = 1.5;
	actforce = 0;

	cout << "\n ***  equilibrating...   ***\n";
	cout << "starting lj rep potential" << endl;
	for (int i = 0; i < 40000; i++) {
	    calcForces_ljrep(atomPositions, verlNums, verlLists, atomForces);
	    calcRodForces(dummyList, rodForces, atomForces);
	    integrateRods(atomPositions, verlNums, verlLists, verlPos, atomForces, rodForces, dt*0.1);
	    if (i%writeEvery == 0) {
    //	    cout << "timestep, energy verl: " << i << " " << calcEnergy_lj(atomPositions, verlNums, verlLists) << endl;
    //	    cout << "timestep, energy: " << i << " " << calcEnergy_ir(atomPositions, verlNums, verlLists) << endl;
		cout << "timestep: " << i << endl;
		writeConfig(atomPositions, movieName);
		if (checkExplosion(atomPositions) == true) {
		    cout << "*** SIMULATION EXPLODED ***" << endl;
		    cout << "*** u suck lol ***" << endl;
		    return 0;
		}
	    }
	}


	cout << "starting real pot (small timestep)" << endl;
	TEMP = saveTEMP;
	ir_invTemp = ir_invTemp;

	// Next do an annealing in timestep:
	double dt_equil = 5e-9;
	for (int i = 0; i < 10000; i++) {
	    calcForces_ir(atomPositions, verlNums, verlLists, atomForces);
	    calcRodForces(dummyList, rodForces, atomForces);
	    integrateRods(atomPositions, verlNums, verlLists, verlPos, atomForces, rodForces, dt_equil);
	    if (i%writeEvery == 0) {
    //	    cout << "timestep, energy verl: " << i << " " << calcEnergy_lj(atomPositions, verlNums, verlLists) << endl;
//		dt_equil += dt_adder;
		cout << "timestep, dt, energy: " << i << " " << dt_equil << " " << calcEnergy_ir(atomPositions, verlNums, verlLists) << endl;
		writeConfig(atomPositions, movieName);
		if (checkExplosion(atomPositions) == true) {
		    cout << "*** SIMULATION EXPLODED ***" << endl;
		    cout << "*** u suck lol ***" << endl;
		    return 0;
		}
	    }
	}
	double dt0 = 1e-6;
	double myu = 5*(gsl_rng_uniform(mrRand)+1);
	dt_equil = dt0;
	for (int i = 0; i < equilSteps; i++) {
	    calcForces_ir(atomPositions, verlNums, verlLists, atomForces);
	    calcRodForces(dummyList, rodForces, atomForces);
	    integrateRods(atomPositions, verlNums, verlLists, verlPos, atomForces, rodForces, dt_equil*myu);
	    if (i%writeEvery == 0) {
		cout << "timestep, dt, energy: " << i << " " << dt_equil*myu << " " << calcEnergy_ir(atomPositions, verlNums, verlLists) << endl;
		writeConfig(atomPositions, movieName);
		if (checkExplosion(atomPositions) == true) {
		    cout << "*** SIMULATION EXPLODED ***" << endl;
		    cout << "*** u suck lol ***" << endl;
		    return 0;
		}
	    }
	    myu = 5*(gsl_rng_uniform(mrRand)+1);
	}
    }

    dt = timeStep;
    ir_invTemp = saveinv;
    actforce = saveFORCE;

    cout << "***************************************************************\n";
    cout << "***************************************************************\n";
    cout << "***            starting real simulation                     ***\n";
    cout << "***                  8)     8)    8)                        ***\n";
    cout << "***************************************************************\n";
    for (int i = 0; i < totalSteps; i++) {
	calcForces_ir(atomPositions, verlNums, verlLists, atomForces);
	calcRodForces(dummyList, rodForces, atomForces);
	integrateRods(atomPositions, verlNums, verlLists, verlPos, atomForces, rodForces, dt);
	if (i%writeEvery == 0) {
//	    cout << "timestep, energy verl: " << i << " " << calcEnergy_lj(atomPositions, verlNums, verlLists) << endl;
	    cout << "timestep, energy: " << i << " " << calcEnergy_ir(atomPositions, verlNums, verlLists) << endl;
	    writeConfig(atomPositions, movieName);
		if (checkExplosion(atomPositions) == true) {
		    cout << "*** SIMULATION EXPLODED ***" << endl;
		    cout << "*** u suck lol ***" << endl;
		    return 0;
		}
	}
    }

    cout << "finished" << endl;

}




// Function Definitions
double distsq(double p1x, double p1y, double p2x, double p2y) {
    double dx = p1x - p2x;
    double dy = p1y - p2y;
    if (dx > boxrad) dx -= 2*boxrad;
    if (dy > boxrad) dy -= 2*boxrad;
    if (dx < -boxrad) dx += 2*boxrad;
    if (dy < -boxrad) dy += 2*boxrad;
    return (dx*dx + dy*dy);
}

void setVerlet(genarray<double> &atomPositions, genarray<int> &verlNums, genarray<int> &verlLists, genarray<double> &verlPos) {
// Set verlet lists.
    double drsq;
    double r_vsq = (r_cut + skin)*(r_cut + skin);
    
    // Reset verlNum, verlPos
    for (int i = 0; i < natoms; i++) {
	verlNums(i) = 0;
	verlPos(i*2) = atomPositions(i*3);
	verlPos(i*2+1) = atomPositions(i*3+1);
    }

//    cout << "in set Verlet " << endl;
    // For every pair of particles, check if they are in each other's verlet lists.
    for (int i = 0; i < natoms; i++) {
	for (int j = (i+1); j < natoms; j++) {
	    drsq = distsq(atomPositions(i*3), atomPositions(i*3+1), atomPositions(j*3), atomPositions(j*3+1));
//	    cout << "dr: " << sqrt(drsq) << endl;
	    if (drsq < r_vsq) {
		verlLists(i*natoms+verlNums(i)) = j;
//		cout << "adding pair i,j: " << i << " " << j << endl;
//		cout << "cout j, verlNums(i): " << j << " " << verlNums(i) << endl;
//		verlLists(verlNums(j)) = i;		    // For MD, we only need to count one particle of every pair for energy/force calculations
		verlNums(i) = verlNums(i) + 1;
//		verlNums(j) = verlNums(j) + 1;
	    }
	}
    }
}


void checkVerlet(genarray<double> &atomPositions, genarray<int> &verlNums, genarray<int> &verlLists, genarray<double> &verlPos) { 
// Check if a particle has moved beyond its skin distance.  If so, recalculate Verlet lists.
    double drsq;
    double skinsq = skin*skin;

    for (int i = 0; i < natoms; i++) {
	drsq = distsq(atomPositions(i*3), atomPositions(i*3+1), verlPos(i*2), verlPos(i*2+1));
	if (drsq < skin) {
	    setVerlet(atomPositions, verlNums, verlLists, verlPos);
	    i = natoms;
	}
    }
}

double calcEnergy_lj(genarray<double> &atomPositions, genarray<int> &verlNums, genarray<int> &verlLists) { 
    int verlMax = 0;
    int j_index;
    double dx, dy, dr, dr6inv;
    double sig6 = pow(ljSig, 6.0);
    double en = 0;
    double en_shift = 4*ljEps*sig6*(sig6/pow(r_cut, 12.0) - 1/pow(r_cut, 6.0));
    for (int i = 0; i < natoms; i++) {
	verlMax = verlNums(i);
	for (int j = 0; j < verlMax; j++) {
	    j_index = verlLists(i*natoms+j);
	    dx = atomPositions(i*3) - atomPositions(j_index*3);
	    dy = atomPositions(i*3+1) - atomPositions(j_index*3+1);
//	    cout << "i, j, ix, iy, jx, jy: " << i << " " << j_index << " " << atomPositions(i*3) << " " << atomPositions(i*3+1) << " ";
//	    cout << atomPositions(j_index*3) << " " << atomPositions(j_index*3+1) << endl;
	    if (dx > boxrad) dx -= 2*boxrad;
	    if (dy > boxrad) dy -= 2*boxrad;
	    if (dx < -boxrad) dx += 2*boxrad;
	    if (dy < -boxrad) dy += 2*boxrad;
	    dr = sqrt(dx*dx + dy*dy);
	    if (dr < r_cut) {
		dr6inv = 1/pow(dr, 6.0);
		en += 4*ljEps*sig6*dr6inv*(sig6*dr6inv - 1) - en_shift;
	    }
	}
    }
    return en;
}

double calcForces_lj(genarray<double> &atomPositions, genarray<int> &verlNums, genarray<int> &verlLists, genarray<double> &forces) { 
    int verlMax = 0;
    int j_index = 0;
    double dx, dy, dr, dr7inv, sig6;
    double fx, fy, myForce;
    for (int i = 0; i < forces.length(); i++) {
	forces(i) = 0;
    }

    if (yesRods == 1) {
	for (int i = 0; i < natoms; i++) {
	    verlMax = verlNums(i);
	    for (int j = 0; j < verlMax; j++) {
		if (!( (max(i,j_index)%2 == 1) and (fabs(j_index - i) == 1) )) {
		    j_index = verlLists(i*natoms+j);
		    dx = atomPositions(i*3) - atomPositions(j_index*3);
		    dy = atomPositions(i*3+1) - atomPositions(j_index*3+1);
		    if (dx > boxrad) dx -= 2*boxrad;
		    if (dy > boxrad) dy -= 2*boxrad;
		    if (dx < -boxrad) dx += 2*boxrad;
		    if (dy < -boxrad) dy += 2*boxrad;
		    dr = sqrt(dx*dx + dy*dy);
		    dr7inv = 1.0/pow(dr, 7);
		    sig6 = pow(ljSig, 6);
		    myForce = 4.0*ljEps*sig6*(12*dr7inv*dr7inv*dr*sig6 - 6*dr7inv);
		    fx = dx/dr*myForce;
		    fy = dy/dr*myForce;
		    forces(i*2) += fx;
		    forces(i*2+1) += fy;
		    forces(j_index*2) -= fx;
		    forces(j_index*2+1) -= fy;
		}
	    }
	}
    }
    else {
	for (int i = 0; i < natoms; i++) {
	    verlMax = verlNums(i);
	    for (int j = 0; j < verlMax; j++) {
		j_index = verlLists(i*natoms+j);
		dx = atomPositions(i*3) - atomPositions(j_index*3);
		dy = atomPositions(i*3+1) - atomPositions(j_index*3+1);
		if (dx > boxrad) dx -= 2*boxrad;
		if (dy > boxrad) dy -= 2*boxrad;
		if (dx < -boxrad) dx += 2*boxrad;
		if (dy < -boxrad) dy += 2*boxrad;
		dr = sqrt(dx*dx + dy*dy);
		dr7inv = 1.0/pow(dr, 7);
		sig6 = pow(ljSig, 6);
		myForce = 4.0*ljEps*sig6*(12*dr7inv*dr7inv*dr*sig6 - 6*dr7inv);
		fx = dx/dr*myForce;
		fy = dy/dr*myForce;
		forces(i*2) += fx;
		forces(i*2+1) += fy;
		forces(j_index*2) -= fx;
		forces(j_index*2+1) -= fy;
	    }
	}
    }
}

void integrate(genarray<double> &atomPositions, genarray<int> &verlNums, genarray<int> &verlLists, genarray<double> &verlPos, genarray<double> &forces, double dt) {
    double ran_constPar = sqrt(2*myGammaPar*TEMP/dt);
    double ran_constRot = sqrt(2*myGammaRot*TEMP/dt);
//    double mass_inv = 1/mass/gamma;
    double theta, vx, vy;
    double xn, yn;
//    cout << "in integrate" << endl;

    for (int i = 0; i < natoms; i++) {
	
//	// Langevin Part
	vx = forces(i*2) + ran_constPar*gsl_ran_gaussian(mrRand, 1.0);
	vy = forces(i*2+1) + ran_constPar*gsl_ran_gaussian(mrRand, 1.0);

//	// Diffusion check: turn off forces (ideal particles)
//	vx = ran_const*gsl_ran_gaussian(mrRand, 1.0);
//	vy = ran_const*gsl_ran_gaussian(mrRand, 1.0);

	// Active Force
	theta = atomPositions(i*3+2) + ran_constRot*gsl_ran_gaussian(mrRand, 1.0)*dt;
	if (theta > PI) theta -= 2*PI;
	if (theta < -PI) theta += 2*PI;
	vx += actforce*cos(theta);
	vy += actforce*sin(theta);

//	cout << "\t after active: vx, vy: " << vx << " " << vy << endl;

	// Update
	xn = atomPositions(i*3) + vx*dt;
	yn = atomPositions(i*3+1) + vy*dt; 
	if (xn > boxrad) xn -= 2*boxrad;
	if (yn > boxrad) yn -= 2*boxrad;
	if (xn < -boxrad) xn += 2*boxrad;
	if (yn < -boxrad) yn += 2*boxrad;
	atomPositions(i*3) = xn;
	atomPositions(i*3+1) = yn;
	atomPositions(i*3+2) = theta;
    }
    checkVerlet(atomPositions, verlNums, verlLists, verlPos);

}


void writeConfig(genarray<double> &atomPositions, string outName) {
    double length_inv = 0.5/boxrad;
    ofstream myOut;
    myOut.open(outName.c_str(), ios::app);
    myOut << "ITEM: TIMESTEP\n" << "0" << "\nITEM: NUMBER OF ATOMS\n";
    myOut << natoms << "\n" << "ITEM: BOX BOUNDS pp pp pp\n";
    myOut << -boxrad << " " << boxrad << "\n" << -boxrad << " " << boxrad << "\n" << -boxrad << " " << boxrad << "\n";
    myOut << "ITEM: ATOMS id type xs ys zs q\n";
    for (int i = 0; i < natoms; i++) {
	myOut << i+1 << " 1 " << (atomPositions(i*3)+boxrad)*length_inv << " " << (atomPositions(i*3+1)+boxrad)*length_inv;
	myOut << " 0 " << atomPositions(i*3+2) << "\n";	    // z-coordinate and angle
    }
    myOut.close();
}



void initConfig(genarray<double> &atomPositions) {
// Initial configuration is a square lattice
//    double d = 1.0;
    double d = (2*boxrad-2.0)/sqrt(0.52*natoms);
    double mytheta = 0;
    double x0 = -boxrad + 1.0;
    double y0 = -boxrad + 1.0;
    cout << "natoms is: " << natoms << endl;

    if (yesRods == 1) {
	int nrods = natoms/2.0;
	for (int i = 0; i < nrods; i++) {
	    // The rod is made of two overlapping disks with centers one radius apart.
	    // Particle 1:
	    atomPositions(i*6) = x0;	
	    atomPositions(i*6+1) = y0;	
	    mytheta = 2*PI*(gsl_rng_uniform(mrRand)-0.5);
//	    atomPositions(i*6+2) = mytheta;
	    atomPositions(i*6+2) = PI*0.5;

	    // Particle 2:
	    atomPositions(i*6+3) = x0;	
	    atomPositions(i*6+4) = y0 + rodSpacing;
//	    atomPositions(i*6+5) = mytheta;
	    atomPositions(i*6+5) = PI*0.5;
	    
	    if ((x0+d) > boxrad) {
		x0 = -boxrad + 1.0;
		y0 = y0 + d;
	    }
	    else {
		x0 = x0 + d;
	    }
	}

////	// Rod 1 -    |__
//	atomPositions(0) = 0;
//	atomPositions(1) = 0;
//	atomPositions(2) = 1.570796;
//	atomPositions(3) = 0;
//	atomPositions(4) = 0.5;
//	atomPositions(5) = 1.570796;
//	// Rod 2
//	atomPositions(6) = 1.2;
//	atomPositions(7) = 0;
//	atomPositions(8) = 0;
//	atomPositions(9) = 1.7;
//	atomPositions(10) = 0;
//	atomPositions(11) = 0;
	
////	// Rod 1 - vertical
//	atomPositions(0) = 0;
//	atomPositions(1) = 0;
//	atomPositions(2) = 1.570796;
//	atomPositions(3) = 0;
//	atomPositions(4) = 0.5;
//	atomPositions(5) = 1.570796;
//	// Rod 2
//	atomPositions(6) = 1.2;
//	atomPositions(7) = 0;
//	atomPositions(8) = 1.570796;
//	atomPositions(9) = 1.2;
//	atomPositions(10) = 0.5;
//	atomPositions(11) = 1.570796;

//	// Rod 1 - horizontal
//	atomPositions(0) = 0;
//	atomPositions(1) = 0;
//	atomPositions(2) = 0;
//	atomPositions(3) = 0.5;
//	atomPositions(4) = 0;
//	atomPositions(5) = 0;
//	// Rod 2
//	atomPositions(6) = 0;
//	atomPositions(7) = 1.2;
//	atomPositions(8) = 0;
//	atomPositions(9) = 0.5;
//	atomPositions(10) = 1.2;
//	atomPositions(11) = 0;

    }
    else {
	for (int i = 0; i < natoms; i++) {
	    atomPositions(i*3) = x0;	
	    atomPositions(i*3+1) = y0;	
	    atomPositions(i*3+2) = 2*PI*(gsl_rng_uniform(mrRand)-0.5);
	    if ((x0+d) > boxrad) {
		x0 = -boxrad + 1.0;
		y0 = y0 + d;
	    }
	    else {
		x0 = x0 + d;
	    }
	}
    }
    cout << "finished" << endl;

}


double calcEnergy_lj_noverl(genarray<double> &atomPositions, genarray<int> &verlNums, genarray<int> &verlLists) {
    // calculates energy without using verlet lists - to debug

    double dx, dy, dr, dr6inv;
    double sig6 = pow(ljSig, 6.0);
    double en = 0;
    double en_shift = 4*ljEps*sig6*(sig6/pow(r_cut, 12.0) - 1/pow(r_cut, 6.0));
    for (int i = 0; i < natoms; i++) {
	for (int j = (i+1); j < natoms; j++) {
	    dx = atomPositions(i*3) - atomPositions(j*3);
	    dy = atomPositions(i*3+1) - atomPositions(j*3+1);
//	    cout << "i, j, ix, iy, jx, jy: " << i << " " << j << " " << atomPositions(i*3) << " " << atomPositions(i*3+1) << " ";
//	    cout << atomPositions(j*3) << " " << atomPositions(j*3+1) << endl;
	    if (dx > boxrad) dx -= 2*boxrad;
	    if (dy > boxrad) dy -= 2*boxrad;
	    if (dx < -boxrad) dx += 2*boxrad;
	    if (dy < -boxrad) dy += 2*boxrad;
	    dr = sqrt(dx*dx + dy*dy);
	    if (dr < r_cut) {
		dr6inv = 1/pow(dr, 6.0);
		en += 4*ljEps*sig6*dr6inv*(sig6*dr6inv - 1) - en_shift;
	    }
	    	
	}
    }
    return en;
} 


double calcEnergy_ir(genarray<double> &atomPositions, genarray<int> &verlNums, genarray<int> &verlLists) { 
// Potential is from Imperio, Reatto, Zeppari, PRE 78, 2008
// I am adding in an additional factor that scales the interaction strength (instead of raising/lower the temperature).
    int j_index;
    double dx, dy, dr, dr6inv;
    double sig6 = pow(irSig, 6.0);
    double en = 0;
    double en_shift = irA*pow(irSig/r_cut, irN) - irEa*irSig*irSig/irRa/irRa*exp(-r_cut/irRa) + irEr*irSig*irSig/irRr/irRr*exp(-r_cut/irRr);

    // debug
    int ctr = 0;
    double r_avg = 0;

    // end debug


    double sig12 = pow(irSig, 12.0);

    for (int i = 0; i < natoms; i++) {
	for (int j = 0; j < verlNums(i); j++) {
	    j_index = verlLists(i*natoms+j);
//	    if (fabs(j_index - i) != 1) {
	    if (!( (max(i,j_index)%2 == 1) and (fabs(j_index - i) == 1) )) {
		dx = atomPositions(i*3) - atomPositions(j_index*3);
		dy = atomPositions(i*3+1) - atomPositions(j_index*3+1);
    //	    cout << "i, j, ix, iy, jx, jy: " << i << " " << j_index << " " << atomPositions(i*3) << " " << atomPositions(i*3+1) << " ";
    //	    cout << atomPositions(j_index*3) << " " << atomPositions(j_index*3+1) << endl;
		if (dx > boxrad) dx -= 2*boxrad;
		if (dy > boxrad) dy -= 2*boxrad;
		if (dx < -boxrad) dx += 2*boxrad;
		if (dy < -boxrad) dy += 2*boxrad;
		dr = sqrt(dx*dx + dy*dy);
		if (dr < r_cut) {
		    en += ir_invTemp*( irA*pow(irSig/dr, irN) - irEa*irSig*irSig/irRa/irRa*exp(-dr/irRa) + irEr*irSig*irSig/irRr/irRr*exp(-dr/irRr) );
    //		cout << "\tdr: " << dr << endl;
//		    ctr++;
//		    r_avg += dr;
		}
	    }
	}
    }
//    cout << "average dr: " << r_avg/ctr << endl;
    return en;
}



double calcForces_ir(genarray<double> &atomPositions, genarray<int> &verlNums, genarray<int> &verlLists, genarray<double> &forces) { 
// I am adding in an additional factor that scales the interaction strength (instead of raising/lower the temperature).
    int j_index = 0;
    double dx, dy, dr, dr7inv, sig6;
    double fx, fy, myForce, myConst1, myConst2;
    for (int i = 0; i < forces.length(); i++) {
	forces(i) = 0;
    }

//    const double irSig = 1.0;
//    const double irN = 12;
//    const double irA = 0.018;
//    const double irRa = irSig;
//    const double irEa = 1.0;
//    const double irRr = 2*irSig;
//    const double irEr = 1.0;
    if (yesRods == 1) {
	myConst1 = irEa*irSig*irSig/irRa*irRa*irRa;
	myConst2 = irEr*irSig*irSig/irRr/irRr/irRr;
	for (int i = 0; i < natoms; i++) {
	    for (int j = 0; j < verlNums(i); j++) {
		j_index = verlLists(i*natoms+j);
		// Ignore bonded particles:
		if (!( (max(i,j_index)%2 == 1) and (fabs(j_index - i) == 1) )) {
//		    cout << "\tforce pair: " << i << "\t" << j_index << endl;
		    dx = atomPositions(i*3) - atomPositions(j_index*3);
		    dy = atomPositions(i*3+1) - atomPositions(j_index*3+1);
//		    cout << "\tdx, dy: " << dx << " " << dy << endl;
		    if (dx > boxrad) dx -= 2*boxrad;
		    if (dy > boxrad) dy -= 2*boxrad;
		    if (dx < -boxrad) dx += 2*boxrad;
		    if (dy < -boxrad) dy += 2*boxrad;
		    // Can calculate some of the prefactors outside the loop! 
		    dr = sqrt(dx*dx + dy*dy);
//		    myForce = ir_invTemp*( irA*irN*pow(irSig, irN)*pow(dr, (-irN-1)) - irEa*irSig*irSig/irRa*irRa*irRa*exp(-dr/irRa) + irEr*irSig*irSig/irRr/irRr/irRr*exp(-dr/irRr) );
		    myForce = ir_invTemp*( irA*irN*pow(irSig, irN)*pow(dr, (-irN-1)) - myConst1*exp(-dr/irRa) + myConst2*exp(-dr/irRr) );
//		    if (myForce > 1000) {
//			cout<< "\t bigForce.  pair, dr: " << i << " " << j_index << " " << dr << endl; 
//			cout << "x1, y1, x2, y2: " <<  atomPositions(i*3) << " " << atomPositions(i*3+1) << ", " << atomPositions(j_index*3) << " " << atomPositions(j_index*3+1) << endl;
//		    }
		    fx = dx/dr*myForce;
		    fy = dy/dr*myForce;
		    forces(i*2) += fx;
		    forces(i*2+1) += fy;
		    forces(j_index*2) -= fx;
		    forces(j_index*2+1) -= fy;
		}
	    }
	}
    }
    else {
	for (int i = 0; i < natoms; i++) {
	    for (int j = 0; j < verlNums(i); j++) {
		j_index = verlLists(i*natoms+j);
		dx = atomPositions(i*3) - atomPositions(j_index*3);
		dy = atomPositions(i*3+1) - atomPositions(j_index*3+1);
		if (dx > boxrad) dx -= 2*boxrad;
		if (dy > boxrad) dy -= 2*boxrad;
		if (dx < -boxrad) dx += 2*boxrad;
		if (dy < -boxrad) dy += 2*boxrad;
		// Can calculate some of the prefactors outside the loop! 
		dr = sqrt(dx*dx + dy*dy);
		myForce = ir_invTemp*( irA*irN*pow(irSig, irN)*pow(dr, (-irN-1)) - irEa*irSig*irSig/irRa*irRa*irRa*exp(-dr/irRa) + irEr*irSig*irSig/irRr/irRr/irRr*exp(-dr/irRr) );
		fx = dx/dr*myForce;
		fy = dy/dr*myForce;
		forces(i*2) += fx;
		forces(i*2+1) += fy;
		forces(j_index*2) -= fx;
		forces(j_index*2+1) -= fy;
	    }
	}
    }
}


void paramReader (string fileName)
{
// Read in Parameters
    ifstream in;
    in.open(fileName.c_str(), ios::in);
    string junk1;
    in >> junk1 >> natoms;
    in >> junk1 >> volfrac;
    in >> junk1 >> TEMP;
    in >> junk1 >> timeStep;
    in >> junk1 >> totalSteps;
    in >> junk1 >> writeEvery;
    in >> junk1 >> mySeed;
    in >> junk1 >> actforce;
    in >> junk1 >> ir_Temp;
    yesRods = 0;
    in >> junk1 >> yesRods;
    in >> junk1 >> rodSpacing;
    if (yesRods == 1) {
	cout << "Simulating rods..." << endl;
    }
    else {
	cout << "Simulating disks..." << endl;
    }
    in >> junk1 >> yesRestart;
    in >> junk1 >> restartName;
    in >> junk1 >> yesEquil;
    
    in.close();

}


void calcRodForces(genarray<double> &rodList, genarray<double> &rodForces, genarray<double> &forces) {
    // rodForces should be half the size of forces (each rod is composed of two particles)
    // rodList...structure?
    // 1-2, 3-4, 5-6, 7-8...?  Let's do that for now.  You can specify different pairs later.
    int fourxi = 0;
    if (natoms%2 != 0)
	cout << "Problems! odd number of atoms for a rod system." << endl;
    // Add in exception handling...catch/throw later?

    double fcmx, fcmy;
    int nrods = natoms/2.0;
//    rodForces.resize(natoms);
    
    for (int i = 0; i < nrods; i++) {
	fourxi = i*4;
	fcmx = forces(fourxi)+forces(fourxi+2);
	fcmy = forces(fourxi+1)+forces(fourxi+3);
	rodForces(i*2) = fcmx;
	rodForces(i*2+1) = fcmy;
    }
//    cout << "finished calcRodForces\n";

}


//void integrateRods(genarray<double> &atomPositions, genarray<int> &verlNums, genarray<int> &verlLists, genarray<double> &verlPos, genarray<double> &forces, genarray<double> &rodForces, double dt) {
//    double ran_const = sqrt(2*myGamma*TEMP/dt);
//    double ran_constR = sqrt(2*myGammaR*TEMP/dt);
//    double randx, randy;
//    double theta, vx, vy, dtheta, theta_old;
//    double xn, yn;
////    cout << "in integrate" << endl;
//
//    int nrods = natoms/2.0;
//    int rodIndex = 0;
//    double irRad = irSig*0.5;
//
//    for (int i = 0; i < nrods; i++) {
//    // The sequence:
//    //	1) Rotational part: Draw a random number for angular displacement, rotate both rod particles.
//    //	2) Get the pairwise potential part of the force.
//    //	3) Add the random noise.
//    //	4) Add the active force along the new angle.
//    //	5) Integrate forwards in time xn = xo + v*dt.
//	
////	// Langevin Part
////	vx = forces(i*2) + ran_const*gsl_ran_gaussian(mrRand, 1.0);
////	vy = forces(i*2+1) + ran_const*gsl_ran_gaussian(mrRand, 1.0);
//
//	rodIndex = int(i/2);		// rod particle pairs are 0-1, 2-3, 4-5, etc.
//	if (i%2 == 0) {
//	    // Both of the rod particles get the same random force:
//	    // (Do we need to divide by 2?  or change the variance to account for the rod mass?)
//	    randx = ran_const*gsl_ran_gaussian(mrRand, 1.0);
//	    randy = ran_const*gsl_ran_gaussian(mrRand, 1.0);
//	}
////	vx = rodForces(rodIndex*2)/2.0 + randx;
////	vy = rodForces(rodIndex*2+1)/2.0 + randy; 
//	vx = randx;
//	vy = randy;
//
//	// The rod rotates about its center of mass due to bath collisions.
//	// Rotations: Update the thermostatted angle, then rotate the vector about the axis.  2D is easy: 
//	if (i%2 == 0) {
//	    dtheta = ran_constR*gsl_ran_gaussian(mrRand, 1.0)*dt;
//	}
//	theta_old = atomPositions(i*3+2);
//	theta = theta_old + dtheta;
//	if (theta > PI) theta -= 2*PI;
//	if (theta < -PI) theta += 2*PI;
//	// Now rotate:
//	// what about the second particle...it needs to rotate the other way.
//	if (i%2 == 0) {
//	    // The even particle is anti (the force directs in the opposite direction.
//	    xn = atomPositions(i*3) + irRad*(cos(theta) - cos(theta_old));
//	    yn = atomPositions(i*3+1) + irRad*(sin(theta) - sin(theta_old));
//	}
//	else {
//	    xn = atomPositions(i*3) + irRad*(cos(theta) - cos(theta_old));
//	    yn = atomPositions(i*3+1) + irRad*(sin(theta) - sin(theta_old));
//	}
//	
//
//	// Apply the Active Force
//	vx += actforce*cos(theta);
//	vy += actforce*sin(theta);
//
////	cout << "\t after active: vx, vy: " << vx << " " << vy << endl;
//
//	// Update
//	xn = xn + vx*dt;
//	yn = yn + vy*dt; 
//	if (xn > boxrad) xn -= 2*boxrad;
//	if (yn > boxrad) yn -= 2*boxrad;
//	if (xn < -boxrad) xn += 2*boxrad;
//	if (yn < -boxrad) yn += 2*boxrad;
//	atomPositions(i*3) = xn;
//	atomPositions(i*3+1) = yn;
//	atomPositions(i*3+2) = theta;
//    }
//
//    checkVerlet(atomPositions, verlNums, verlLists, verlPos);
//
////}

void integrateRods(genarray<double> &atomPositions, genarray<int> &verlNums, genarray<int> &verlLists, genarray<double> &verlPos, genarray<double> &forces, genarray<double> &rodForces, double dt) {
    double diffusePerp = TEMP/myGammaPerp;
    double diffusePar = TEMP/myGammaPar;
    double diffuseRot = TEMP/myGammaRot;
    double ran_constPerp = sqrt(2*diffusePerp*TEMP/dt);		    // D = kT/gamma
    double ran_constPar = sqrt(2*diffusePar*TEMP/dt);		    // D = kT/gamma
    double ran_constR = sqrt(2*diffuseRot*TEMP/dt);
//    cout << "diffusion par, perp, rot: " << diffusePar << ", " << diffusePerp << ", " << diffuseRot << endl;
    double randx, randy;
    double theta, vx, vy, dtheta, theta_old;
    double xn, yn, mycos, mysin;
    double xodd_n, yodd_n, xeven_n, yeven_n, xcm, ycm;
    double xe, ye, xo, yo, dx, dy;

    // 27 Aug 2015:
    //		    Let's redo the dumbbells properly following (Heptner and Dzubiella, Moleular Physics 2015).  Need to thermostat the
    //		    parallel, perpendicular, and rotational motion of the rods separately.

    // New variables:
    double rpari, rparj, rperpi, rperpj, ni, nj, nperpi, nperpj;	// i is the first component, j is the second component of the vector
    double fxo, fyo, fxe, fye;
    double fpari, fparj, fperpi, fperpj;
    double rpari_n, rparj_n, rperpi_n, rperpj_n;
    double randpar, randperp, randrot;

    int nrods = natoms/2.0;
    int rodIndex = 0;
    int rodxthree = 0;
    double irRad = irSig*0.5;
    double rX, rY, myTorque = 0;
    double randAngle = 0;
    double torqE, torqO;

//    cout << "rodForces.length(): " << rodForces.length() << endl;
    for (int i = 0; i < nrods; i++) {
	
	rodIndex = int(i*2);		// rod particle pairs are 0-1, 2-3, 4-5, etc.
	rodxthree = rodIndex*3;		// used for accessing array elements.

	// Calculate center of mass:
	xe = atomPositions(rodxthree);
	ye = atomPositions(rodxthree+1);
	xo = atomPositions(rodxthree+3);
	yo = atomPositions(rodxthree+4);
	dx = xo - xe;
	dy = yo - ye;
	if (fabs(dx) > boxrad) xe = xe + mySign(dx)*2*boxrad;		// Checks if a rod lies across a boundary.
	if (fabs(dy) > boxrad) ye = ye + mySign(dy)*2*boxrad;
	xcm = 0.5*(xe + xo);
	ycm = 0.5*(ye + yo);
//	cout << "xcm, ycm: " << xcm << ", " << ycm << endl;

	theta_old = atomPositions(rodxthree+2);
	ni = cos(theta_old);
	nj = sin(theta_old);
//	cout << "ni, nj: " << ni << ", " << nj << endl;

	// parallel component along the longitudinal axis:  rpar = (n dot r)*n
	rpari = (ni*xcm + nj*ycm)*ni;
	rparj = (ni*xcm + nj*ycm)*nj;
	fpari = (ni*rodForces(i*2) + nj*rodForces(i*2+1))*ni;
	fparj = (ni*rodForces(i*2) + nj*rodForces(i*2+1))*nj;
//	cout << "rpari, rparj: " << rpari << ", " << rparj << endl;
	
	// perpendicular component to the longitudinal axis: rperp = r - rpar
	rperpi = xcm - rpari;
	rperpj = ycm - rparj;
	fperpi = rodForces(i*2) - fpari;
	fperpj = rodForces(i*2+1) - fparj;
//	cout << "rperpi, rperpj: " << rperpi << ", " << rperpj << endl;


	// Evolve rperp, rpar:
	randpar = ran_constPar*gsl_ran_gaussian(mrRand, 1.0);
	randperp = ran_constPerp*gsl_ran_gaussian(mrRand, 1.0);

	    // (let's take mass of the dumbbell to be 1, otherwise need to divide Force by the mass)
	rpari_n = rpari + dt*( diffusePar*(fpari + actforce*ni) + randpar*ni );
	rparj_n = rparj + dt*( diffusePar*(fparj + actforce*nj) + randpar*nj );		    // active force only acts longitudinally

	rperpi_n = rperpi + dt*( diffusePerp*fperpi + randperp*(-nj) );
	rperpj_n = rperpj + dt*( diffusePerp*fperpj + randperp*(ni) );

//	cout << "diffusepar*fpari+actforce: " << diffusePar*(fpari+actforce) << ", randpar*ni: " << randpar*ni << endl;
//	cout << "diffusepar*fparj+actforce: " << diffusePar*(fparj+actforce) << ", randpar*nj: " << randpar*nj << endl;
//	cout << "diffusePerp*fperpi: " << diffusePerp*fperpi << ", randperp*ni: " << randperp*ni << endl;
//	cout << "diffusePerp*fperpj: " << diffusePerp*fperpj << ", randperp*nj: " << randperp*nj << endl;
//	cout << "rpar old, new: " << rpari << ", " << rparj << "    " << rpari_n << ", " << rparj_n << endl;
//	cout << "rperp old, new: " << rperpi << ", " << rperpj << "    " << rperpi_n << ", " << rperpj_n << endl;

	//
	// The active force vector points from the even particle to the odd particle:
	//	e.g.   0 ---> 1
	//
	
	// Calculate Torque:  n * F, where F = F_odd - F_even
	fxe = forces(rodIndex*2);
	fye = forces(rodIndex*2+1);
	fxo = forces(rodIndex*2+2);
	fyo = forces(rodIndex*2+3);
	myTorque = 2*irSig/2.0 * ( -ni*(fyo - fye) + nj*(fxo - fxe) );

	// Evolve Torque:
	randrot = ran_constR*gsl_ran_gaussian(mrRand, 1.0);
	theta = theta_old + dt*(diffuseRot*myTorque + randrot);

//	theta = theta_old;

	
	// Convert from dumbbell coords back to Cartesian and update positions:
	//	    r = rperp + rpar

	xcm = rpari_n + rperpi_n;
	ycm = rparj_n + rperpj_n;
//	cout << "xcm, ycm new: " << xcm << ", " << ycm << endl;
	mycos = rodSpacing*0.5*cos(theta);
	mysin = rodSpacing*0.5*sin(theta);

//	cout << "\n";

	xodd_n = xcm + mycos;
	yodd_n = ycm + mysin;
	xeven_n = xcm - mycos;
	yeven_n = ycm - mysin;
	
	// PBC: 
	if (xodd_n > boxrad)  { xodd_n -= 2*boxrad; } 
	if (yodd_n > boxrad)  { yodd_n -= 2*boxrad; } 
	if (xodd_n < -boxrad) { xodd_n += 2*boxrad; } 
	if (yodd_n < -boxrad) { yodd_n += 2*boxrad; } 
	if (xeven_n > boxrad) { xeven_n -= 2*boxrad; }
	if (yeven_n > boxrad) { yeven_n -= 2*boxrad; } 
	if (xeven_n < -boxrad){ xeven_n += 2*boxrad; } 
	if (yeven_n < -boxrad){ yeven_n += 2*boxrad; } 
	//

	atomPositions(rodxthree) = xodd_n;
	atomPositions(rodxthree+1) = yodd_n;
	atomPositions(rodxthree+2) = theta;

	atomPositions(rodxthree+3) = xeven_n;
	atomPositions(rodxthree+4) = yeven_n;
	atomPositions(rodxthree+5) = theta;
    }
	
    checkVerlet(atomPositions, verlNums, verlLists, verlPos);

}

int mySign(double x) {
    if (x > 0) return 1;
    if (x < 0) return -1;
    return 0;
}


double calcForces_ljrep(genarray<double> &atomPositions, genarray<int> &verlNums, genarray<int> &verlLists, genarray<double> &forces) { 
// lj repulsive - to help equilibrating rods
// V(r) = 4*epsilon*( (sigma/r)^12 + (sigma/r)^6 )
    int verlMax = 0;
    int j_index = 0;
    double dx, dy, dr, dr7inv, sig6;
    double fx, fy, myForce;
    for (int i = 0; i < forces.length(); i++) {
	forces(i) = 0;
    }

    if (yesRods == 1) {
	for (int i = 0; i < natoms; i++) {
	    verlMax = verlNums(i);
	    for (int j = 0; j < verlMax; j++) {
		if (!( (max(i,j_index)%2 == 1) and (fabs(j_index - i) == 1) )) {
		    j_index = verlLists(i*natoms+j);
		    dx = atomPositions(i*3) - atomPositions(j_index*3);
		    dy = atomPositions(i*3+1) - atomPositions(j_index*3+1);
		    if (dx > boxrad) dx -= 2*boxrad;
		    if (dy > boxrad) dy -= 2*boxrad;
		    if (dx < -boxrad) dx += 2*boxrad;
		    if (dy < -boxrad) dy += 2*boxrad;
		    dr = sqrt(dx*dx + dy*dy);
		    dr7inv = 1.0/pow(dr, 7);
		    sig6 = pow(ljSig, 6);
		    myForce = 4.0*ljEps*sig6*(12*dr7inv*dr7inv*dr*sig6 + 6*dr7inv);
		    fx = dx/dr*myForce;
		    fy = dy/dr*myForce;
		    forces(i*2) += fx;
		    forces(i*2+1) += fy;
		    forces(j_index*2) -= fx;
		    forces(j_index*2+1) -= fy;
		}
	    }
	}
    }
    else {
	for (int i = 0; i < natoms; i++) {
	    verlMax = verlNums(i);
	    for (int j = 0; j < verlMax; j++) {
		j_index = verlLists(i*natoms+j);
		dx = atomPositions(i*3) - atomPositions(j_index*3);
		dy = atomPositions(i*3+1) - atomPositions(j_index*3+1);
		if (dx > boxrad) dx -= 2*boxrad;
		if (dy > boxrad) dy -= 2*boxrad;
		if (dx < -boxrad) dx += 2*boxrad;
		if (dy < -boxrad) dy += 2*boxrad;
		dr = sqrt(dx*dx + dy*dy);
		dr7inv = 1.0/pow(dr, 7);
		sig6 = pow(ljSig, 6);
		myForce = 4.0*ljEps*sig6*(12*dr7inv*dr7inv*dr*sig6 + 6*dr7inv);
		fx = dx/dr*myForce;
		fy = dy/dr*myForce;
		forces(i*2) += fx;
		forces(i*2+1) += fy;
		forces(j_index*2) -= fx;
		forces(j_index*2+1) -= fy;
	    }
	}
    }
}

void readRestart (string restartName, genarray<double> &atomPositions)
{
    int ctr = 0;
    int atomctr = 0;
    double xlo, xhi, ylo, yhi, zlo, zhi;
    double x,y,z;
    double theta;
    cout << "reading file " << restartName << " for restart." << endl;
    ifstream in(restartName.c_str());
    string junk, myLine;

    while (getline(in, myLine)) {
	if (myLine != "") {
	    stringstream ss(myLine);
	    if (ctr == 5) {
		ss >> xlo >> xhi;
	    }
	    if (ctr == 6) {
		ss >> ylo >> yhi;
	    }
	    if (ctr == 7) {
		ss >> zlo >> zhi;
	    }
	    if (ctr%(natoms+9) > 8) {
		ss >> junk >> junk >> x >> y >> z >> theta;
		// Each frame has natoms*3 elements.
		atomPositions(atomctr*3) = x*(xhi-xlo) + xlo;
		atomPositions(atomctr*3+1) = y*(yhi-ylo) + ylo;
		atomPositions(atomctr*3+2) = theta;
		atomctr++;
	    }
	    ctr++;
	}
    }

    cout << "finished reading:" << endl;
    for (int i = 0; i < natoms; i++) {
	cout << "\t" << atomPositions(i*3) << " " << atomPositions(i*3+1) << " " << atomPositions(i*3+2) << "\n";
    }
    cout << "finished reading:" << endl;



}


bool checkExplosion(genarray<double> &atomPositions) {


    double myX, myY;
    for (int i = 0; i < natoms; i++) {
	myX = atomPositions(i*3);
	myY = atomPositions(i*3+1);
	if (myX > boxrad)
	    return true;
	if (myX < -boxrad)
	    return true;
	if (myY > boxrad)
	    return true;
	if (myY < -boxrad)
	    return true;
	if (isnan(myX) || isnan(myY))
	    return true;
    }
    return false;
}


