#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>
#include <TROOT.h>
#include <TLorentzVector.h>
#include <iomanip>
#include <math.h>
#include <sys/time.h>
#include <string>
#include <fstream>
#include <sstream>
#include "../myvector.h"
#include "../dilep_input.h"
#include "../utilities.h"

namespace Dilep {
	namespace CPU {
		void Csqrt(double _ar, double _ai, double _my[]);
		void cubic(double a[], double rr[], double ri[]);
		void toz(double k[], double l[], double g[]);
		void my_qu( double my_in[], double my_val[]);
		double calcMass(double x, double y, double z, double e);

		std::vector<myvector>* dilep(double t_mass[], double w_mass[], 
									double in_mpx[], double in_mpy[], double in_mpz[],
									TLorentzVector* lep_a, TLorentzVector* lep_b, 
									TLorentzVector* bl_a, TLorentzVector* bl_b);

		void dilep (vector<DilepInput> &vdi);

		void dilep (DilepInput &di);

		std::vector<myvector>* calc_dilep(double t_mass[], double w_mass[], 
										double in_mpx[], double in_mpy[], double in_mpz[],
										TLorentzVector& lep_a, TLorentzVector& lep_b, 
										TLorentzVector& bl_a, TLorentzVector& bl_b);
	}
}
