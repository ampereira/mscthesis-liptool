
#include "../../LipMiniAnalysis/TLorentzVectorWFlags.h"

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

inline bool operator> (const ttDKF_Best_Sol&, const ttDKF_Best_Sol&);
inline bool operator< (const ttDKF_Best_Sol&, const ttDKF_Best_Sol&);
inline bool operator<= (const ttDKF_Best_Sol&, const ttDKF_Best_Sol&);
inline bool operator>= (const ttDKF_Best_Sol&, const ttDKF_Best_Sol&);
inline bool operator== (const ttDKF_Best_Sol&, const ttDKF_Best_Sol&);
inline bool operator!= (const ttDKF_Best_Sol&, const ttDKF_Best_Sol&);
