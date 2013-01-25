#ifndef DILEP_INPUT
#define DILEP_INPUT

#include <TLorentzVector.h>
#include <TRandom.h>
#include <fstream>
#include <vector>
#include "../../LipMiniAnalysis/LipMiniAnalysis.h"

using namespace std;

class DilepInput {

	TLorentzVector z_lep, c_lep;
	TLorentzVector z_bj, c_bj;
	TLorentzVector z_bl, c_bl;

	TLorentzVectorWFlags z_bjWFlags, c_bjWFlags;
	TLorentzVectorWFlags z_lepWFlags, c_lepWFlags;
	TLorentzVectorWFlags jet1_HiggsWFlags, jet2_HiggsWFlags;

	double in_mpx[2];
	double in_mpy[2];
	double in_mpz[2];

	double t_mass[2];
	double w_mass[2];

	double MissPx, MissPy;

public:
	DilepInput (TLorentzVector _z_lep, TLorentzVector _c_lep, TLorentzVector _z_bj, TLorentzVector _c_bj, 
				TLorentzVectorWFlags _z_bjWFlags, TLorentzVectorWFlags _c_bjWFlags, TLorentzVectorWFlags _z_lepWFlags,
				TLorentzVectorWFlags _c_lepWFlags, TLorentzVectorWFlags _jet1_HiggsWFlags, TLorentzVectorWFlags _jet2_HiggsWFlags, 
				double _in_mpx[], double _in_mpy[], double _in_mpz[], double _MissPx,
				double _MissPy, double _t_mass[], double _w_mass[]);

	DilepInput (const DilepInput &other);

	void print (char* filename);
	void applyVariance (float res, int seed);

	// Getters
	TLorentzVector getZlep (void) const;
	TLorentzVector getClep (void) const;
	TLorentzVector getZbj (void) const;
	TLorentzVector getCbj (void) const;
	TLorentzVector getZbl (void) const;
	TLorentzVector getCbl (void) const;
	TLorentzVectorWFlags getZlepW (void) const;
	TLorentzVectorWFlags getClepW (void) const;
	TLorentzVectorWFlags getZbjW (void) const;
	TLorentzVectorWFlags getCbjW (void) const;
	TLorentzVectorWFlags getJet1HiggsW (void) const;
	TLorentzVectorWFlags getJet2HiggsW (void) const;
	double getMissPx (void) const;
	double getMissPy (void) const;
	double getInMpx (int) const;
	double getInMpy (int) const;
	double getInMpz (int) const;
	double getTmass (int) const;
	double getWmass (int) const;
};

vector<DilepInput> applyVariance (DilepInput di, float res, int amount, int seed);

#endif
