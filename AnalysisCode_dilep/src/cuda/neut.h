#ifndef NEUT
#define NEUT

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
	namespace GPU {
		__device__ __host__ void Csqrt(double _ar, double _ai, double _my[]);
		__device__ __host__ void cubic(double a[], double rr[], double ri[]);
		__device__ __host__ void toz(double k[], double l[], double g[]);
		__device__ __host__ void my_qu( double my_in[], double my_val[]);
		__device__ __host__ double calcMass(double x, double y, double z, double e);

		std::vector<myvector>* dilep(unsigned iterations, double t_mass[], double w_mass[], 
									double in_mpx[], double in_mpy[], double in_mpz[],
									TLorentzVector* lep_a, TLorentzVector* lep_b, 
									TLorentzVector* bl_a, TLorentzVector* bl_b, int *hasSol);

		void dilep (vector<DilepInput> &vdi);

		void dilep (DilepInput &di);
		void dilep (DilepInput &di, int EveNumber);

		std::vector<myvector>* calc_dilep(double t_mass[], double w_mass[], 
										double in_mpx[], double in_mpy[], double in_mpz[],
										TLorentzVector* lep_a, TLorentzVector* lep_b, 
										TLorentzVector* bl_a, TLorentzVector* bl_b);

		__global__ void calc_dilep(double t_mass[], double w_mass[], 
						double in_mpx[], double in_mpy[], double in_mpz[], double lep_a[], 
						double lep_b[], double bl_a[], double bl_b[], 
						double nc[], int a[]);
	}
}

#endif
