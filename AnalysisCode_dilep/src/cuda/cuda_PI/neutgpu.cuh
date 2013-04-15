
#include <iostream>
#include <vector>
#include <cmath>
#include <cuda.h>
#include "myvector.cuh"
#include <cuda_runtime.h>
#include <TLorentzVector.h>

using std::vector;

__host__
void gpu_init(int gridSize, int blockSize);

__host__
void runDilep(double t_mass[], double w_mass[], 
			double in_mpx[], double in_mpy[], double in_mpz[],
			TLorentzVector* lep_a, TLorentzVector* lep_b, 
			TLorentzVector* bl_a, TLorentzVector* bl_b,
			vector<myvector>* &resultfinal);

__device__
double calcMass(double x, double y, double z, double e);

__global__
void dilep_kernel(double t_mass[], double w_mass[], double in_mpx[],
		double in_mpy[], double in_mpz[], double lep_a[5],
		double lep_b[5], double bl_a[5], double bl_b[5],
		double nc[], int *a);

__device__
void toz_kernel(double k[], double l[], double g[]);

__device__
void my_qu_kernel(double my_in[], double my_val[]);

__device__
void Csqrt_kernel(double _ar, double _ai, double _my[]);

__device__
void cubic_kernel(double a[], double rr[], double ri[]);

