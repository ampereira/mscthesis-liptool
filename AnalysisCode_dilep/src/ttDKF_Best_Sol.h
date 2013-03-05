
#include "../../LipMiniAnalysis/LipMiniAnalysis.h"

class ttDKF_Best_Sol {

	double probability;
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
	double mHiggsJet1_ttDKF;
	double mHiggsJet2_ttDKF;

public:
	ttDKF_Best_Sol (double, TLorentzVectorWFlags, TLorentzVectorWFlags, TLorentzVectorWFlags, TLorentzVectorWFlags,
					TLorentzVectorWFlags, TLorentzVectorWFlags, TLorentzVectorWFlags, TLorentzVectorWFlags, TLorentzVectorWFlags,
					TLorentzVectorWFlags, TLorentzVectorWFlags, TLorentzVectorWFlags, TLorentzVectorWFlags, TLorentzVectorWFlags,
					double, double);

};

inline bool operator> (const ttbar_ttDKF&, const ttbar_ttDKF&);
inline bool operator< (const ttbar_ttDKF&, const ttbar_ttDKF&);
inline bool operator<= (const ttbar_ttDKF&, const ttbar_ttDKF&);
inline bool operator>= (const ttbar_ttDKF&, const ttbar_ttDKF&);
inline bool operator== (const ttbar_ttDKF&, const ttbar_ttDKF&);
inline bool operator!= (const ttbar_ttDKF&, const ttbar_ttDKF&);
