
#include "ttDKF_Best_Sol.h"

ttDKF_Best_Sol::ttDKF_Best_Sol (double _probability, double _mHiggsJet1_ttDKF, double _mHiggsJet2_ttDKF,
								TLorentzVectorWFlags _n1_ttDKF, TLorentzVectorWFlags _n2_ttDKF,
								TLorentzVectorWFlags _b1_ttDKF, TLorentzVectorWFlags _b2_ttDKF, TLorentzVectorWFlags _l1_ttDKF,
								TLorentzVectorWFlags _l2_ttDKF, TLorentzVectorWFlags _W1_ttDKF, TLorentzVectorWFlags _W2_ttDKF,
								TLorentzVectorWFlags _t1_ttDKF, TLorentzVectorWFlags _t2_ttDKF, TLorentzVectorWFlags _ttbar_ttDKF,
								TLorentzVectorWFlags _b1_Higgs_ttDKF, TLorentzVectorWFlags _b2_Higgs_ttDKF, TLorentzVectorWFlags _Higgs_ttDKF) {

	probability = _probability;
	mHiggsJet1_ttDKF = _mHiggsJet1_ttDKF;
	mHiggsJet2_ttDKF = _mHiggsJet2_ttDKF;
	n1_ttDKF = _n1_ttDKF;
	n2_ttDKF = _n2_ttDKF;
	b1_ttDKF = _b1_ttDKF;
	b2_ttDKF = _b2_ttDKF;
	l1_ttDKF = _l1_ttDKF;
	l2_ttDKF = _l2_ttDKF;
	W1_ttDKF = _W1_ttDKF;
	W2_ttDKF = _W2_ttDKF;
	t1_ttDKF = _t1_ttDKF;
	t2_ttDKF = _t2_ttDKF;
	ttbar_ttDKF = _ttbar_ttDKF;
	b1_Higgs_ttDKF = _b1_Higgs_ttDKF;
	b2_Higgs_ttDKF = _b2_Higgs_ttDKF;
	Higgs_ttDKF = _Higgs_ttDKF;
}

ttDKF_Best_Sol::ttDKF_Best_Sol (void) {
	probability = 0.0;
	mHiggsJet1_ttDKF = 0.0;
	mHiggsJet2_ttDKF = 0.0;
}

// Getters
double ttDKF_Best_Sol::getProb (void) const {
	return probability;
}

double ttDKF_Best_Sol::getMHiggsJet (int x) const {
	return (x == 1) ? mHiggsJet1_ttDKF : mHiggsJet2_ttDKF;
}

TLorentzVectorWFlags ttDKF_Best_Sol::getN (int x) const {
	return (x == 1) ? n1_ttDKF : n2_ttDKF;
}

TLorentzVectorWFlags ttDKF_Best_Sol::getB (int x) const {
	return (x == 1) ? b1_ttDKF : b2_ttDKF;
}

TLorentzVectorWFlags ttDKF_Best_Sol::getL (int x) const {
	return (x == 1) ? l1_ttDKF : l2_ttDKF;
}

TLorentzVectorWFlags ttDKF_Best_Sol::getW (int x) const {
	return (x == 1) ? W1_ttDKF : W2_ttDKF;
}

TLorentzVectorWFlags ttDKF_Best_Sol::getT (int x) const {
	return (x == 1) ? t1_ttDKF : t2_ttDKF;
}

TLorentzVectorWFlags ttDKF_Best_Sol::getTTbar (void) const {
	return ttbar_ttDKF;
}

TLorentzVectorWFlags ttDKF_Best_Sol::getBHiggs (int x) const {
	return (x == 1) ? b1_Higgs_ttDKF : b2_Higgs_ttDKF;
}

TLorentzVectorWFlags ttDKF_Best_Sol::getHiggs (void) const {
	return Higgs_ttDKF;
}

inline bool operator> (const ttDKF_Best_Sol &t1, const ttDKF_Best_Sol &t2) {
	return (t1.getProb() > t2.getProb()) ? true : false;
}

inline bool operator< (const ttDKF_Best_Sol &t1, const ttDKF_Best_Sol &t2) {
	return operator> (t2, t1);
}

inline bool operator<= (const ttDKF_Best_Sol &t1, const ttDKF_Best_Sol &t2) {
	return !operator> (t1, t2);
}

inline bool operator>= (const ttDKF_Best_Sol &t1, const ttDKF_Best_Sol &t2) {
	return !operator< (t1, t2);
}

inline bool operator== (const ttDKF_Best_Sol &t1, const ttDKF_Best_Sol &t2) {
	return (t1.getProb() == t2.getProb()) ? true : false;
}

inline bool operator!= (const ttDKF_Best_Sol &t1, const ttDKF_Best_Sol &t2) {
	return !operator== (t1, t2);
}
