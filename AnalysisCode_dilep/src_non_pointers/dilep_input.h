#ifndef DILEP_INPUT
#define DILEP_INPUT

#include <TLorentzVector.h>
#include <TRandom.h>
#include <fstream>
#include <vector>
#include "myvector.h"
#include "utilities.h"
#include "../../LipMiniAnalysis/LipMiniAnalysis.h"

#define SEED 123456.0

extern TRandom3 *t_rnd;

using namespace std;


class DilepInput {
	
	//TRandom3 t_rnd;

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

	int hasSolution;

	vector<myvector> result;

public:

	DilepInput (TLorentzVector _z_lep, TLorentzVector _c_lep, TLorentzVector _z_bj, TLorentzVector _c_bj, 
				TLorentzVectorWFlags _z_bjWFlags, TLorentzVectorWFlags _c_bjWFlags, TLorentzVectorWFlags _z_lepWFlags,
				TLorentzVectorWFlags _c_lepWFlags, TLorentzVectorWFlags _jet1_HiggsWFlags, TLorentzVectorWFlags _jet2_HiggsWFlags, 
				double _in_mpx[], double _in_mpy[], double _in_mpz[], double _MissPx,
				double _MissPy, double _t_mass[], double _w_mass[]);

	DilepInput (TLorentzVector _z_lep, TLorentzVector _c_lep, TLorentzVector _z_bj, TLorentzVector _c_bj, 
				TLorentzVectorWFlags _z_bjWFlags, TLorentzVectorWFlags _c_bjWFlags, TLorentzVectorWFlags _z_lepWFlags,
				TLorentzVectorWFlags _c_lepWFlags, TLorentzVectorWFlags _jet1_HiggsWFlags, TLorentzVectorWFlags _jet2_HiggsWFlags, 
				double _in_mpx[], double _in_mpy[], double _in_mpz[], double _MissPx,
				double _MissPy, double _t_mass[], double _w_mass[], TRandom3 &_t_rnd);
	DilepInput (const DilepInput &other);
	DilepInput (void);

	void print (char* filename);
	void applyVariance (float res);
	void applyVariance (float res, int seed);

	// Getters
	__attribute__((target(mic))) TLorentzVector getZlep (void) const;
	__attribute__((target(mic))) TLorentzVector getClep (void) const;
	__attribute__((target(mic))) TLorentzVector getZbj (void) const;
	__attribute__((target(mic))) TLorentzVector getCbj (void) const;
	__attribute__((target(mic))) TLorentzVector getZbl (void) const;
	__attribute__((target(mic))) TLorentzVector getCbl (void) const;
	__attribute__((target(mic))) TLorentzVectorWFlags getZlepW (void) const;
	__attribute__((target(mic))) TLorentzVectorWFlags getClepW (void) const;
	__attribute__((target(mic))) TLorentzVectorWFlags getZbjW (void) const;
	__attribute__((target(mic))) TLorentzVectorWFlags getCbjW (void) const;
	__attribute__((target(mic))) int getZlepWisb (void) const;
	__attribute__((target(mic))) int getClepWisb (void) const;
	__attribute__((target(mic))) int getZbjWisb (void) const;
	__attribute__((target(mic))) int getCbjWisb (void) const;
	__attribute__((target(mic))) TLorentzVectorWFlags getJet1HiggsW (void) const;
	__attribute__((target(mic))) TLorentzVectorWFlags getJet2HiggsW (void) const;
	__attribute__((target(mic))) double getMissPx (void) const;
	__attribute__((target(mic))) double getMissPy (void) const;
	__attribute__((target(mic))) double getInMpx (int) const;
	__attribute__((target(mic))) double getInMpy (int) const;
	__attribute__((target(mic))) double getInMpz (int) const;
	__attribute__((target(mic))) double getTmass (int) const;
	__attribute__((target(mic))) double getWmass (int) const;
	__attribute__((target(mic))) int getHasSol (void) const;
	__attribute__((target(mic))) vector<myvector> getResult (void) const;

	// Setters
	__attribute__((target(mic))) void setHasSol (int);
	__attribute__((target(mic))) void setResult (vector<myvector> *);
	__attribute__((target(mic))) void setZblCbl (void);
};

void applyVariance (vector<DilepInput> &vdi, float res, int amount);
vector<DilepInput> applyVariance (vector<DilepInput> &vdi, float res, int amount, int seed);

#endif
