#ifndef DILEP_INPUT
#define DILEP_INPUT

#include <TLorentzVector.h>
#include <TRandom.h>
#include <fstream>
#include <vector>
#include "../../LipMiniAnalysis/LipMiniAnalysis.h"

using namespace std;

class DilepInput {

	TLorentzVector *z_lep, *c_lep;
	TLorentzVector *z_bj, *c_bj;
	TLorentzVector z_bl, c_bl;

	TLorentzVectorWFlags *z_bjWFlags, *c_bjWFlags;
	TLorentzVectorWFlags *z_lepWFlags, *c_lepWFlags;

	double in_mpx[2];
	double in_mpy[2];
	double in_mpz[2];

	double t_mass[2];
	double w_mass[2];

	double MissPx, MissPy;

public:
	DilepInput (TLorentzVector _z_lep, TLorentzVector _c_lep, TLorentzVector _z_bj, TLorentzVector _c_bj, 
				TLorentzVectorWFlags _z_bjWFlags, TLorentzVectorWFlags _c_bjWFlags, TLorentzVectorWFlags _z_lepWFlags,
				TLorentzVectorWFlags _c_lepWFlags, double _in_mpx[], double _in_mpy[], double _in_mpz[], double _MissPx,
				double _MissPy, double _t_mass[], double _w_mass[]);

	//DilepInput (DilepInput &other);

	void print (char* filename);
	void applyVariance (float res, int seed);

	// Getters
	TLorentzVector getZlep (void);
	TLorentzVector getClep (void);
	TLorentzVector getZbj (void);
	TLorentzVector getCbj (void);
	TLorentzVector getZbl (void);
	TLorentzVector getCbl (void);
	TLorentzVectorWFlags getZlepW (void);
	TLorentzVectorWFlags getClepW (void);
	TLorentzVectorWFlags getZbjW (void);
	TLorentzVectorWFlags getCbjW (void);
	double getMissPx (void);
	double getMissPy (void);
	double getInMpx (int);
	double getInMpy (int);
	double getInMpz (int);
	double getTmass (int);
	double getWmass (int);
};

//vector<DilepInput> applyVariance (DilepInput di, float res, int amount, int seed);

#endif
