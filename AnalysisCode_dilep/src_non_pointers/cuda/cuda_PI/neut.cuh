#ifndef NEUT_CUDA
#define NEUT_CUDA

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

//#include <cuda.h>
//#include <cuda_runtime.h>
//#include <curand_kernel.h>


#define PI 3.14159265358979312
// Transpose the 2D abstraction of the array to 1D
#define TO1D(nc,tid,sol,did)	nc[tid*16+sol*4+did]

// Normalize the random number generated
//#define cuRand(x) (((double) curand(x))/((double) UINT_MAX))


namespace Dilep {
	namespace GPU {

		//__device__ 
		//double gaus_kernel (double mean, double sigma, curandStateMtgp32 *state);
		//__device__ 
		void applyVariance ();
		//__host__ 
		void gpu_init (int blocks, int threads);
		//__device__ void Csqrt (double _ar, double _ai, double _my[]);
		//__device__ void cubic (double a[], double rr[], double ri[]);
		//__device__ void toz (double k[], double l[], double g[]);
		//__device__ void my_qu (double my_in[], double my_val[]);
		//__device__ double calcMass (double x, double y, double z, double e);
		void Csqrt (double _ar, double _ai, double _my[]);
		void cubic (double a[], double rr[], double ri[]);
		void toz (double k[], double l[], double g[]);
		void my_qu (double my_in[], double my_val[]);
		double calcMass (double x, double y, double z, double e);

		std::vector<myvector>* dilep(unsigned iterations, double t_mass[], double w_mass[], 
									double in_mpx[], double in_mpy[], double in_mpz[],
									TLorentzVector* lep_a, TLorentzVector* lep_b, 
									TLorentzVector* bl_a, TLorentzVector* bl_b, int *hasSol);

		void dilep (vector<DilepInput> &vdi, int);

		void dilep (DilepInput &di);

		void calc_dilep(double t_mass[], double w_mass[], double in_mpx[], 
		//__global__ void calc_dilep(double t_mass[], double w_mass[], double in_mpx[], 
										double in_mpy[], double in_mpz[], double lep_a[], double lep_b[], 
										double bl_a[], double bl_b[], double nc[], int a[]);
	}
}
#endif
