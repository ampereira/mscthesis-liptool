#ifndef NEUT
#define NEUT

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <TRandom3.h>
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
#include "../../../LipMiniAnalysis/TLorentzVectorWFlags.h"

#define TPI 3.14159265358979312
#define TO1D(nc,tid,sol,did) nc[tid*16+sol*4+did]
#define STRIDE2(a,i) a[tid * 2 + i]
#define STRIDE5(a,i) a[tid * 5 + i]
#define CUDA_THREAD_STRIDE2(a,i) a[(tid / size) * 2 + i]
#define CUDA_THREAD_STRIDE5(a,i) a[(tid / size) * 5 + i]
#define RESOLUTION 0.02
//#define cuRand(x) (((double) curand(x))/((double) UINT_MAX))

//#define CUDA_CALL(x) if((x) != cudaSuccess) { \
//	printf("Error at %s:%d - %d\n",__FILE__,__LINE__,x); \
//	exit (EXIT_FAILURE);}
#define CUDA_CALL(x) x

extern int dilep_iterations;


namespace Dilep {
	namespace GPU {
		//static TRandom3 _t_rnd_ (123456.0);

		__device__ void Csqrt(double _ar, double _ai, double _my[]);
		__device__ void cubic(double a[], double rr[], double ri[]);
		__device__ void toz(double k[], double l[], double g[]);
		__device__ void my_qu( double my_in[], double my_val[]);
		__device__ double calcMass(double x, double y, double z, double e);
		__device__ void calcMass (double *array);

		//void dilep (DilepInput &di);

		__host__ void dilep (vector<DilepInput> &di);
		//__host__ void dilep (vector<DilepInput> &di);
		//void dilep (vector<DilepInput> &vdi, int x);

		__global__ void dilep_kernel (double _in_mpx[], double _in_mpy[], double _z_lepWFlags[], double _c_lepWFlags[],
			double _z_bjWFlags[], double _c_bjWFlags[], double _z_lep[], double _c_lep[], double _z_bj[], double _c_bj[],
			double *_MissPx, double *_MissPy, double _t_mass[], double _w_mass[], double nc[], int a[], unsigned *size);

		__device__ void calc_dilep(double t_mass[], double w_mass[], 
						double in_mpx[], double in_mpy[], double _lep_a[], 
						double _lep_b[], double _bl_a[], double _bl_b[], 
						double nc[], int a[]);

	}
}

#endif
