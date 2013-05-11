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

#define TPI 3.14159265358979312
#define NUM_THREADS 2
#define TO1D(nc,tid,sol,did) nc[tid*16+sol*4+did]
#define STRIDE2(a,i) a[tid * 2 + i]
#define STRIDE5(a,i) a[tid * 5 + i]

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

		//void dilep (DilepInput &di);
		void dilep (vector<DilepInput> &di);
		void dilep (vector<DilepInput> &vdi, int x);


		std::vector<myvector>* calc_dilep(double t_mass[], double w_mass[], 
										double in_mpx[], double in_mpy[], double in_mpz[],
										TLorentzVector* lep_a, TLorentzVector* lep_b, 
										TLorentzVector* bl_a, TLorentzVector* bl_b);

		__global__ void calc_dilep(double t_mass[], double w_mass[], 
						double in_mpx[], double in_mpy[], double _lep_a[], 
						double _lep_b[], double _bl_a[], double _bl_b[], 
						double nc[], int a[]);

		__host__ void calc_dilep(double t_mass[], double w_mass[], 
						double in_mpx[], double in_mpy[], double _lep_a[], 
						double _lep_b[], double _bl_a[], double _bl_b[], 
						double nc[], int a[], unsigned _tid);
	}
}

#endif
