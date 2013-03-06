
#include "../../LipMiniAnalysis/TLorentzVectorWFlags.h"
#include "../../LipMiniAnalysis/LipMiniAnalysis.h"

class ttDKF_Best_Sol {

	double probability;
	double mHiggsJet1_ttDKF;
	double mHiggsJet2_ttDKF;
	TLorentzVectorWFlags n1_ttDKF;
	TLorentzVectorWFlags n2_ttDKF;
	TLorentzVectorWFlags b1_ttDKF;
	TLorentzVectorWFlags b2_ttDKF;
	TLorentzVectorWFlags l1_ttDKF;
	TLorentzVectorWFlags l2_ttDKF;
	TLorentzVectorWFlags W1_ttDKF;
	TLorentzVectorWFlags W2_ttDKF;
	TLorentzVectorWFlags t1_ttDKF;
	TLorentzVectorWFlags t2_ttDKF;
	TLorentzVectorWFlags ttbar_ttDKF;
	TLorentzVectorWFlags b1_Higgs_ttDKF;
	TLorentzVectorWFlags b2_Higgs_ttDKF;
	TLorentzVectorWFlags Higgs_ttDKF;

public:
	ttDKF_Best_Sol (double, double, double, TLorentzVectorWFlags, TLorentzVectorWFlags, TLorentzVectorWFlags, TLorentzVectorWFlags,
					TLorentzVectorWFlags, TLorentzVectorWFlags, TLorentzVectorWFlags, TLorentzVectorWFlags, TLorentzVectorWFlags,
					TLorentzVectorWFlags, TLorentzVectorWFlags, TLorentzVectorWFlags, TLorentzVectorWFlags, TLorentzVectorWFlags);

	ttDKF_Best_Sol (void);

	// Getters
	double getProb (void) const;
	double getMHiggsJet (int) const;
	TLorentzVectorWFlags getN (int) const;
	TLorentzVectorWFlags getB (int) const;
	TLorentzVectorWFlags getL (int) const;
	TLorentzVectorWFlags getW (int) const;
	TLorentzVectorWFlags getT (int) const;
	TLorentzVectorWFlags getTTbar (void) const;
	TLorentzVectorWFlags getBHiggs (int) const;
	TLorentzVectorWFlags getHiggs (void) const;


};

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

inline ttDKF_Best_Sol& operator= (const ttDKF_Best_Sol& rhs) {
	probability = rhs.getProb();
	mHiggsJet1_ttDKF = rhs.getMHiggsJet(1);
	mHiggsJet2_ttDKF = rhs.getMHiggsJet(2);
	n1_ttDKF = rhs.getN(1);
	n2_ttDKF = rhs.getN(2);
	b1_ttDKF = rhs.getB(1);
	b2_ttDKF = rhs.getB(2);
	l1_ttDKF = rhs.getL(1);
	l2_ttDKF = rhs.getL(2);
	W1_ttDKF = rhs.getW(1);
	W2_ttDKF = rhs.getW(2);
	t1_ttDKF = rhs.getT(1);
	t2_ttDKF = rhs.getT(2);
	ttbar_ttDKF = rhs.getTTbar();
	b1_Higgs_ttDKF = rhs.getBHiggs(1);
	b2_Higgs_ttDKF = rhs.getBHiggs(2);
	Higgs_ttDKF = rhs.getHiggs();

	return *this;
}
