
#include "neutgpu.cuh"

/*
 *    Using the CUDA Occupancy Calculator it was determined that the optimum
 *    block size is 128.
 *    The best grid size, according to our profiling, is 4096.
 *    This creates a total of 524288 threads, i.e., dilep executions.
 *
 *    The input variance is not implemented. Further profiling is required to
 *    better identify the bottlenecks with the variance.
 *
 *    Note: Avoid float point comparations with absolute values, such as == 0.0,
 *    or f1 == f2. Add a threshold instead.
 *
 *    Note 2: Most of the functions suffered modifications to reduce the amount
 *    of local variables used - reducing register spilling and improving performance - 
 *    at the cost of ilegibility.
 */

unsigned int GRID_SIZE;		// number of blocks per grid
unsigned int BLOCK_SIZE;	// number of threads per block

unsigned int NUM_THREADS;

// Initializes the Grid and Block sizes
__host__ 
void gpu_init(int blocks, int threads)	{
	GRID_SIZE = blocks;
	BLOCK_SIZE = threads;
	NUM_THREADS = blocks*threads;
}

__device__
double calcMass(double x, double y, double z, double e){
	double mm, mass;

	mm = e * e - (x * x + y * y + z * z);

	if(mm < 0.0)
		mass = -sqrt(-mm);
	else
		mass = sqrt(mm);

	return mass;
}

// NEUTRINO SOLUTIONS

// dilep kernel - called from the host (CPU) and executed on the device (GPU)
// Many of the inputs were altered from the CPU version
	__global__
void dilep_kernel(double t_mass[], double w_mass[], double in_mpx[], 
		double in_mpy[], double in_mpz[], double lep_a[], 
		double lep_b[], double bl_a[], double bl_b[], 
		double nc[], int a[])
{

	double G_1, G_3;
	double WMass_a, WMass_b, tMass_a, tMass_b;


	WMass_a = w_mass[0];
	tMass_a = t_mass[0];
	WMass_b = w_mass[1];
	tMass_b = t_mass[1];  
	
	G_1 = (WMass_a - lep_a[4]) * (WMass_a + lep_a[4]);
	G_3 = (WMass_b - lep_b[4]) * (WMass_b + lep_b[4]);

	double G_5,G_6,G_7,G_8,G_9,G_10,G_11,G_12;
	G_5 = ( bl_a[0]/bl_a[3] - lep_a[0]/lep_a[3] );
	G_6 = ( bl_a[1]/bl_a[3] - lep_a[1]/lep_a[3] );
	G_7 = ( bl_a[2]/bl_a[3] - lep_a[2]/lep_a[3] );
	G_8 = ( G_1/lep_a[3] - ((tMass_a - bl_a[4]) * (tMass_a + bl_a[4]))/bl_a[3] )/2.;

	G_9 =	( bl_b[0]/bl_b[3] - lep_b[0]/lep_b[3] );
	G_10 =	( bl_b[1]/bl_b[3] - lep_b[1]/lep_b[3] );
	G_11 =	( bl_b[2]/bl_b[3] - lep_b[2]/lep_b[3] );
	G_12 =	( G_3/lep_b[3] - ((tMass_b - bl_b[4]) * (tMass_b + bl_b[4]))/bl_b[3] )/2.;

	///////////////////////////////////////////////////////////////////
	//// 	G_5 *x1 + G_6*y1 + G_7*z1 = G8;  		(6)
	////  	G_9 *x2 + G_10*y2 + G_11*z2 = G12; 		(7)
	////  	2*El_1*sqrt() - 2*(ax1+by1+cz1) = G_1;  	(8)
	////  	2*El_2*sqrt() - 2*(Ax2+By2+Cz2) = G_3;		(9)
	////  	x1+x2 = S;					(10)
	////  	y1+y2 = T;					(11)
	////  	bring z1 and z2 (from 6/7) to 7 and 8
	///////////////////////////////////////////////////////////////////

	//// 1st top decay product /////
	
	double in_a[5],out_a[6];
	in_a[0] = G_8/G_7;
	in_a[1] = -1.0*G_5/G_7;
	in_a[2] = -1.0*G_6/G_7;
	in_a[3] = lep_a[3];
	in_a[4] = G_1;
	toz_kernel(in_a, lep_a, out_a);

	double in_c[5],out_c[6];
	in_c[0] = G_12/G_11;
	in_c[1] = -1*G_9/G_11;
	in_c[2] = -1*G_10/G_11;
	in_c[3] = lep_b[3];
	in_c[4] = G_3;
	toz_kernel(in_c, lep_b, out_c);
	/////////////////////////////////////////////////////
	//////change x2 y2 equation to x1 and y1 by using
	////// 		x1+x2 = S = in_mpx[0]
	////// 		y1+y2 = T = mpy
	/////////////////////////////////////////////////////
	double out_e[6];
	out_e[0] = out_c[0];
	out_e[1] = out_c[1];
	out_e[2] = -1*( out_c[0]*in_mpx[0] + out_c[2] + out_c[4]*in_mpy[0]);
	out_e[3] = -1*( out_c[1]*in_mpy[0] + out_c[3] + out_c[4]*in_mpx[0]);
	out_e[4] = out_c[4]; 
	out_e[5] =( out_c[0]*in_mpx[0]*in_mpx[0] + out_c[1]*in_mpy[0]*in_mpy[0] + 2*out_c[2]*in_mpx[0] + 2*out_c[3]*in_mpy[0] + out_c[5] + 2*out_c[4]*in_mpx[0]*in_mpy[0]);

	///////////////////////////////////////////////////
	///  solve 
	/// {ax2+by2+2dx+2ey+2fxy+g=0		(12)
	/// {Ax2+By2+2Dx+2Ey+2Fxy+G=0		(13)
	/// out_a[6]: 0   1    2    3    4     5
	/// out_a[6]: a   b    d    e    f     g
	/// out_e[6]: A   B    D    E    F     G
	/// if a!=0, everything is OK.
	///
	/// if a==0, then we can get x2 = f(x,y) from (13)
	/// (12) --> [x2 - f(x,y)] + by2 + ... = 0
	///////////////////////////////////////////////////
	
	double fx_1, fx_2, fx_3, fx_4, fx_5;
	double k_1, k_2, k_3, k_4, k_5;

	// bad organization of the code; diminished register spilling
	if ( out_a[0]!=0  ){
		fx_1 = 2.*(out_e[0]*out_a[2] - out_a[0]*out_e[2]);
		fx_2 = 2.*(out_e[0]*out_a[4] - out_a[0]*out_e[4]);
		fx_3 = out_a[0]*out_e[1] - out_e[0]*out_a[1];
		fx_4 = 2.*(out_e[3]*out_a[0] - out_e[0]*out_a[3]);
		fx_5 = out_a[0]*out_e[5] - out_a[5]*out_e[0];

		k_1 = ( out_a[4]*out_a[4] - out_a[0]*out_a[1] )/out_a[0]/out_a[0];
		k_2 = ( 2.*out_a[2]*out_a[4] - 2.*out_a[0]*out_a[3] )/out_a[0]/out_a[0];
		k_3 = ( out_a[2]*out_a[2]-out_a[0]*out_a[5] )/out_a[0]/out_a[0];
		k_4 = -out_a[2]/out_a[0];
		k_5 = -out_a[4]/out_a[0];
	} else {
		if (out_a[0]==0 && out_e[0]!=0 ) {
			fx_1 = 2.*(out_e[0]*(out_a[2] + out_e[2]/out_e[0]) - out_e[2]);
			fx_2 = 2.*(out_e[0]*(out_a[4] + out_e[4]/out_e[0]) - out_e[4]);
			fx_3 = out_e[1] - out_e[0]*(out_a[1] + out_e[1]/out_e[0]);
			fx_4 = 2.*(out_e[3] - out_e[0]*(out_a[3] + out_e[3]/out_e[0]));
			fx_5 = out_e[5] - (out_a[5] + out_e[5]/out_e[0])*out_e[0];

			k_1 = ( (out_a[4] + out_e[4]/out_e[0])*(out_a[4] + out_e[4]/out_e[0]) - (out_a[1] + out_e[1]/out_e[0]) );
			k_2 = ( 2.*(out_a[2] + out_e[2]/out_e[0])*(out_a[4] + out_e[4]/out_e[0]) - 2.*(out_a[3] + out_e[3]/out_e[0]) );
			k_3 = ( (out_a[2] + out_e[2]/out_e[0])*(out_a[2] + out_e[2]/out_e[0])-(out_a[5] + out_e[5]/out_e[0]) );
			k_4 = -(out_a[2] + out_e[2]/out_e[0]);
			k_5 = -(out_a[4] + out_e[4]/out_e[0]);
		}
	}

	if ( out_a[0]==0 && out_e[0]==0){
		return;
	}


	/////
	///// the part above is 
	///// x = (fx3*y**2 + fx4*y + fx5)/(fx1 + fx2*y)
	///// used to get x value once y is known
	/////
	//// if fx1 + fx2*y == 0, then x is 
	//// x = +/-sqrt(k1*y**2 + k2*y + k3) + (k4 + k5*y)
	////


	double g_1 = 4.*out_e[0]*out_e[0]*k_5*k_5 + 4.*out_e[4]*out_e[4] + 8.*out_e[0]*out_e[4]*k_5;
	double m_1 = g_1*k_1;
	double g_2 = 8.*out_e[0]*out_e[0]*k_4*k_5 + 8.*out_e[0]*out_e[2]*k_5 + 8.*out_e[0]*out_e[4]*k_4 + 8.*out_e[2]*out_e[4];
	double g_3 = 4.*out_e[0]*out_e[0]*k_4*k_4 + 4.*out_e[2]*out_e[2] + 8.*out_e[0]*out_e[2]*k_4;
	double g_4 = out_e[0]*k_1 + out_e[0]*k_5*k_5;
	double g_5 = out_e[0]*k_2 + 2.*out_e[0]*k_4*k_5 + 2.*out_e[2]*k_5;
	double g_6 = out_e[0]*k_3 + out_e[0]*k_4*k_4 + 2.*out_e[2]*k_4 + out_e[5];

	double m_2 = g_1*k_2 + g_2*k_1;
	double m_3 = g_1*k_3 + g_2*k_2 + g_3*k_1;
	double m_4 = g_2*k_3 + g_3*k_2;
	double m_5 = g_3*k_3;

	double m_6  = out_e[1]*out_e[1] + 4.*out_e[4]*out_e[4]*k_5*k_5 + 4.*out_e[1]*out_e[4]*k_5;
	double m_7  = 4.*out_e[1]*out_e[3] + 8.*out_e[4]*out_e[4]*k_4*k_5 + 4.*out_e[1]*out_e[4]*k_4 + 8.*out_e[3]*out_e[4]*k_5;
	double m_8  = 4.*out_e[3]*out_e[3] + 4.*out_e[4]*out_e[4]*k_4*k_4 + 8.*out_e[3]*out_e[4]*k_4;
	double m_80 = pow(g_4,2);
	double m_81 = 2*g_4*g_5;
	double m_9  = pow(g_5,2) + 2.*g_4*g_6;
	double m_10 = 2.*g_5*g_6;
	double m_11 = g_6*g_6;

	double m_12 = 	2.*out_e[0]*out_e[1]*k_1 + 2.*out_e[0]*out_e[1]*k_5*k_5 + 4.*out_e[0]*out_e[4]*k_1*k_5 + 4.*out_e[0]*out_e[4]*pow(k_5,3);
	double m_13 = 	2.*out_e[0]*out_e[1]*k_2 + 4.*out_e[0]*out_e[1]*k_4*k_5 + 4.*out_e[1]*out_e[2]*k_5 + 
		4.*out_e[0]*(out_e[3]*k_1 + out_e[3]*k_5*k_5 + out_e[4]*k_1*k_4 + out_e[4]*k_2*k_5) + 
		12.*out_e[0]*out_e[4]*k_4*k_5*k_5 + 8.*out_e[2]*out_e[4]*k_5*k_5;
	double m_14 = 	2.*out_e[0]*out_e[1]*k_3 + 2.*out_e[0]*out_e[1]*k_4*k_4 + 4.*out_e[2]*out_e[1]*k_4 + 2.*out_e[1]*out_e[5] + 4.*out_e[0]*out_e[3]*k_2 + 
		8.*out_e[0]*out_e[3]*k_4*k_5 + 8.*out_e[3]*out_e[2]*k_5 + 4.*out_e[0]*out_e[4]*k_2*k_4 + 4.*out_e[0]*out_e[4]*k_3*k_5 + 
		12.*out_e[0]*out_e[4]*k_4*k_4*k_5 + 16.*out_e[2]*out_e[4]*k_4*k_5 + 4.*out_e[4]*out_e[5]*k_5;
	double m_15 = 	4.*out_e[0]*out_e[3]*(k_3 + k_4*k_4) + 8.*out_e[3]*out_e[2]*k_4 + 4.*out_e[3]*out_e[5] + 4.*out_e[0]*out_e[4]*(k_3*k_4 + pow(k_4,3)) + 
		8.*out_e[2]*out_e[4]*k_4*k_4 + 4.*out_e[4]*out_e[5]*k_4;

	double  re[5];
	re[0] = m_1 - m_6 - m_12 - m_80;
	re[1] = m_2 - m_7 - m_13 - m_81;
	re[2] = m_3 - m_8 - m_9 - m_14;
	re[3] = m_4 - m_10 - m_15;
	re[4] = m_5 - m_11;  



	double output[8];
	my_qu_kernel(re,output);

	int ncand(0);

	double rec_x1, rec_y1, rec_z1, rec_e1, rec_x2, rec_y2, rec_z2, rec_e2;

	unsigned long tid = threadIdx.x + blockIdx.x * blockDim.x;

	for (int j=0; j<8; j+=2){
		double delta = k_1*output[j]*output[j] + k_2*output[j] + k_3;
		if ( output[j+1]==0 && delta >=0) {
			if ( (fx_1 + fx_2*output[j])!=0 ) {
				rec_x1 = (fx_3*pow(output[j],2) + fx_4*output[j] + fx_5)/(fx_1 + fx_2*output[j]);
			} else {
				rec_x1 = sqrt(delta)+k_4+k_5*output[j];
			}  

			rec_y1 = output[j];
			rec_z1 = G_8/G_7 - G_5*rec_x1/G_7 - G_6*rec_y1/G_7;
			rec_e1 = sqrt(rec_x1*rec_x1 + rec_y1*rec_y1 + rec_z1*rec_z1);
			rec_x2 = in_mpx[0] - rec_x1;
			rec_y2 = in_mpy[0] - rec_y1;
			rec_z2 = G_12/G_11 - G_9*rec_x2/G_11 - G_10*rec_y2/G_11;
			rec_e2 = sqrt(rec_x2*rec_x2 + rec_y2*rec_y2 + rec_z2*rec_z2);
			
			// self-consistence check and control of the solutions

			double m_w11 = calcMass(rec_x1+lep_a[0], rec_y1+lep_a[1], rec_z1+lep_a[2], rec_e1+lep_a[3]);
			double m_w12 = calcMass(rec_x2+lep_b[0], rec_y2+lep_b[1], rec_z2+lep_b[2], rec_e2+lep_b[3]);
			double m_t11 = calcMass(rec_x1+ bl_a[0], rec_y1+ bl_a[1], rec_z1+ bl_a[2], rec_e1+ bl_a[3]);
			double m_t12 = calcMass(rec_x2+ bl_b[0], rec_y2+ bl_b[1], rec_z2+ bl_b[2], rec_e2+ bl_b[3]);

			// m_delta_mass is 1000.0
			bool m_good_eq1 = ( fabs(in_mpx[0] -(rec_x1+rec_x2)) <= 0.01 ) * true + 
							  ( fabs(in_mpx[0] -(rec_x1+rec_x2)) > 0.01 ) * false;
			bool m_good_eq2 = ( fabs(in_mpy[0] -(rec_y1+rec_y2)) <= 0.01 ) * true +
							  ( fabs(in_mpy[0] -(rec_y1+rec_y2)) > 0.01 ) * false;
			bool m_good_eq3 = ( fabs(m_w11 - w_mass[0]) <= 1000.0 ) * true + 
							  ( fabs(m_w11 - w_mass[0]) > 1000.0 ) * false;
			bool m_good_eq4 = ( fabs(m_w12 - w_mass[1]) <= 1000.0 ) * true +
							  ( fabs(m_w12 - w_mass[1]) > 1000.0 ) * false;
			bool m_good_eq5 = ( fabs(m_t11 - t_mass[0]) <= 1000.0 ) * true +
							  ( fabs(m_t11 - t_mass[0]) > 1000.0 ) * false;
			bool m_good_eq6 = ( fabs(m_t12 - t_mass[1]) <= 1000.0 ) * true +
							  ( fabs(m_t12 - t_mass[1]) <= 1000.0 ) * false;

			bool cond = m_good_eq1 && m_good_eq2 && m_good_eq3 && m_good_eq4 && m_good_eq5 && m_good_eq6;
			syncthreads();
			nc[tid * 16 + 2*j] = cond * rec_x1;
			nc[tid * 16 + 2*j + 1] = cond * rec_y1;
			nc[tid * 16 + 2*j + 2] = cond * rec_z1;
			nc[tid * 16 + 2*j + 3] = cond * rec_z2;
			ncand += cond * 1;
		}
	}

	// indicates the number of solutions that this thread found
	a[tid] = ncand;
}

#define TO1D(nc,tid,sol,did)	nc[tid*16+sol*4+did]

// =========================================
// function that runs the dilep_kernel
// =========================================

	__host__
void runDilep(double t_mass[], double w_mass[], 
		double in_mpx[], double in_mpy[], double in_mpz[],
		TLorentzVector* lep_a, TLorentzVector* lep_b, 
		TLorentzVector* bl_a, TLorentzVector* bl_b, vector<myvector>* &resultfinal)
{
	float time;
	double *dev_t_mass, *dev_w_mass, *dev_in_mpx, *dev_in_mpy, *dev_in_mpz;

	double a[5],b[5],c[5],d[5];
	
	a[0] = lep_a->Px();
	a[1] = lep_a->Py();
	a[2] = lep_a->Pz();
	a[3] = lep_a->E();
	a[4] = lep_a->M();

	b[0] = lep_b->Px();
	b[1] = lep_b->Py();
	b[2] = lep_b->Pz();
	b[3] = lep_b->E();
	b[4] = lep_b->M();

	c[0] = bl_a->Px();
	c[1] = bl_a->Py();
	c[2] = bl_a->Pz();
	c[3] = bl_a->E();
	c[4] = bl_a->M();
     	
	d[0] = bl_b->Px();
	d[1] = bl_b->Py();
	d[2] = bl_b->Pz();
	d[3] = bl_b->E();
	d[4] = bl_b->M();

	double *dev_lep_a, *dev_lep_b, *dev_bl_a, *dev_bl_b;
	double nc[16*NUM_THREADS];
	double *dev_nc;
	int count[NUM_THREADS], *dev_count;

	// GPU memory allocation of the inputs and outputs of the dilep kernel
	cudaMalloc(&dev_t_mass, 2*sizeof(double));
	cudaMalloc(&dev_w_mass, 2*sizeof(double));
	cudaMalloc(&dev_in_mpx, 2*sizeof(double));
	cudaMalloc(&dev_in_mpy, 2*sizeof(double));
	cudaMalloc(&dev_in_mpz, 2*sizeof(double));

	cudaMalloc(&dev_lep_a, sizeof(a));
	cudaMalloc(&dev_lep_b, sizeof(b));
	cudaMalloc(&dev_bl_a, sizeof(c));
	cudaMalloc(&dev_bl_b, sizeof(d));
	
	cudaMalloc(&dev_nc, 16*NUM_THREADS*sizeof(double));
	cudaMalloc(&dev_count, NUM_THREADS*sizeof(int));

	// transfer the inputs to GPU memory
	cudaMemcpy(dev_t_mass, t_mass, 2*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_w_mass, w_mass, 2*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_in_mpx, in_mpx, 2*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_in_mpy, in_mpy, 2*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_in_mpz, in_mpz, 2*sizeof(double), cudaMemcpyHostToDevice);

	cudaMemcpy(dev_lep_a, &a,	5*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_lep_b, &b,	5*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_bl_a, &c,	5*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_bl_b, &d,	5*sizeof(double), cudaMemcpyHostToDevice);

	//cudaMemcpy(dev_nc, nc, 16*NUM_THREADS*sizeof(double), cudaMemcpyHostToDevice);
	//cudaMemcpy(dev_count, count,NUM_THREADS*sizeof(int), cudaMemcpyHostToDevice);

	// define the dimensions of the grid and blocks
	// i.e. the number of times dilep is executed
	dim3 dimGrid(GRID_SIZE, 1);
	dim3 dimBlock(BLOCK_SIZE, 1);

	// dilep kernel call
	dilep_kernel<<<dimGrid,dimBlock>>>(
			dev_t_mass, dev_w_mass, dev_in_mpx, dev_in_mpy, dev_in_mpz, 
			dev_lep_a, dev_lep_b, dev_bl_a, dev_bl_b, dev_nc, dev_count);

	// memory transfer of the results from the GPU
	cudaMemcpy(nc, dev_nc, 16*NUM_THREADS*sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(count, dev_count,NUM_THREADS*sizeof(int), cudaMemcpyDeviceToHost);

	unsigned int aux_size = GRID_SIZE*BLOCK_SIZE*16, is;

	// reconstruction of the normal output of dilep
	for(int thread = 0 ; thread < GRID_SIZE*BLOCK_SIZE ; thread++)	{
		for( int sol = 0 ; sol < count[thread] && sol<4 ; sol++)	{
			
			myvector *mv = new myvector( 
				TO1D(nc,thread,sol,0),
				TO1D(nc,thread,sol,1),
				TO1D(nc,thread,sol,2),
				TO1D(nc,thread,sol,3) );
			
			resultfinal->push_back(*mv);
		}
	}
	
	// frees the memory allocated on GPU
	cudaFree(dev_t_mass);
	cudaFree(dev_w_mass);
	cudaFree(dev_in_mpx);
	cudaFree(dev_in_mpy);
	cudaFree(dev_in_mpz);

	cudaFree(dev_lep_a);
	cudaFree(dev_lep_b);
	cudaFree(dev_bl_a);
	cudaFree(dev_bl_b);

	cudaFree(dev_count);
	cudaFree(dev_nc);

}

__device__
void toz_kernel(double k[], double l[], double g[]){
	//// checked !!
	///////////////////////////////////////////////////////////////////////////
	///// bring z=A+Bx+Cy to 2*D*sqrt(x**2+y**2+z**2)-2(ax+by+dz) = E
	///// simplify it to g1*x^2 + g2*y^2 + 2*g3*x + 2*g4*y + 2*g5*xy + g6 = 0
	///////////////////////////////////////////////////////////////////////////
	g[0] = 4*pow(k[3],2)*( 1 + pow(k[1],2)) - 4*pow(l[0],2) -4*pow(l[2],2)*pow(k[1],2) - 8*l[0]*l[2]*k[1];
	if ( g[0]!=0 ) { 
		g[1] = ( 4*pow(k[3],2)*( 1 + pow(k[2],2)) - 4*pow(l[1],2) -4*pow(l[2],2)*pow(k[2],2) - 8*l[1]*l[2]*k[2] )/g[0] ;
		g[2] = ( (4*pow(k[3],2)-4*l[2]*l[2])*k[0]*k[1] - 4*l[0]*l[2]*k[0] - 2*k[4]*l[0] - 2*k[4]*l[2]*k[1] )/g[0];
		g[3] = ( (4*pow(k[3],2)-4*l[2]*l[2])*k[0]*k[2] - 4*l[1]*l[2]*k[0] - 2*k[4]*l[1] - 2*k[4]*l[2]*k[2] )/g[0];
		g[4] = ( (4*pow(k[3],2)-4*l[2]*l[2])*k[1]*k[2] - 4*l[0]*l[1]   - 4*l[0]*l[2]*k[2] - 4*l[1]*l[2]*k[1] )/g[0];
		g[5] = ( (4*pow(k[3],2)-4*l[2]*l[2])*k[0]*k[0] - k[4]*k[4] - 4*k[4]*l[2]*k[0] )/g[0];
		g[0] = 1.0; 
	} else {
		g[1] = ( 4*pow(k[3],2)*( 1 + pow(k[2],2)) - 4*pow(l[1],2) -4*pow(l[2],2)*pow(k[2],2) - 8*l[1]*l[2]*k[2] ) ;
		g[2] = ( (4*pow(k[3],2)-4*l[2]*l[2])*k[0]*k[1] - 4*l[0]*l[2]*k[0] - 2*k[4]*l[0] - 2*k[4]*l[2]*k[1] );
		g[3] = ( (4*pow(k[3],2)-4*l[2]*l[2])*k[0]*k[2] - 4*l[1]*l[2]*k[0] - 2*k[4]*l[1] - 2*k[4]*l[2]*k[2] );
		g[4] = ( (4*pow(k[3],2)-4*l[2]*l[2])*k[1]*k[2] - 4*l[0]*l[1]   - 4*l[0]*l[2]*k[2] - 4*l[1]*l[2]*k[1] );
		g[5] = ( (4*pow(k[3],2)-4*l[2]*l[2])*k[0]*k[0] - k[4]*k[4] - 4*k[4]*l[2]*k[0] );
		g[0] = 0.;  
	}	 
	return;
}

	__device__
void my_qu_kernel( double my_in[], double my_val[])
{

	///////////////////////////////////////////
	/////  find the solution /////////////////
	/////  ax^4+bx^3+cx^2+dx+e=0
	//////////////////////////////////////////
	double real[3]={0,0,0};
	double img[3]={0,0,0};
	double x1_r = 0; double x1_i = 0;
	double x2_r = 0; double x2_i = 0;
	double x3_r = 0; double x3_i = 0;
	double x4_r = 0; double x4_i = 0;


	/////////////////////////////////////////////
	///// in case of a==0, simplify to cubic
	///// bx^3+cx^2+dx+e=0
	/////////////////////////////////////////////

	if ( my_in[0] ==0 && my_in[1]!=0 ){
		double input[4]={my_in[1],my_in[2],my_in[3],my_in[4]};
		cubic_kernel(input,real,img);
		x1_r = real[0];	x1_i = img[0];
		x2_r = real[1];	x2_i = img[1];
		x3_r = real[2];	x3_i = img[2];
		my_val[0] = x1_r;
		my_val[1] = x1_i;
		my_val[2] = x2_r;
		my_val[3] = x2_i;
		my_val[4] = x3_r;
		my_val[5] = x3_i;
		my_val[6] = x4_r;
		my_val[7] = x4_i;

		return;
	}


	/////////////////////////////////////////////
	///// in case of a==0 && b==0, 
	///// simplify to quadratic
	///// cx*2 + d*x + e =0
	/////////////////////////////////////////////

	if ( my_in[0] ==0 && my_in[1]==0 && my_in[2]!=0){
		double alpha = pow(my_in[3],2)-4*my_in[2]*my_in[4];
		if (alpha>=0) {
			x1_r = (-1*my_in[3] + sqrt(alpha))/2/my_in[2];	x1_i = 0;
			x2_r = (-1*my_in[3] - sqrt(alpha))/2/my_in[2];	x2_i = 0;
		} else {
			x1_r = (-1*my_in[3] )/2/my_in[2];	x1_i = sqrt(-alpha)/2/my_in[2];
			x2_r = (-1*my_in[3] )/2/my_in[2];	x2_i = -1*sqrt(-alpha)/2/my_in[2];  	
		}

		my_val[0] = x1_r;
		my_val[1] = x1_i;
		my_val[2] = x2_r;
		my_val[3] = x2_i;
		my_val[4] = x3_r;
		my_val[5] = x3_i;
		my_val[6] = x4_r;
		my_val[7] = x4_i;

		return;
	}


	/////////////////////////////////////////////
	///// in case of a==0 && b==0 && c==0, 
	///// simplify to linear equation dx + e =0
	/////////////////////////////////////////////

	if (my_in[0] ==0 && my_in[1]==0 && my_in[2]==0 && my_in[3]!=0){
		x1_r = -my_in[4]/my_in[3];	x1_i = 0;
		my_val[0] = x1_r;
		my_val[1] = x1_i;
		my_val[2] = x2_r;
		my_val[3] = x2_i;
		my_val[4] = x3_r;
		my_val[5] = x3_i;
		my_val[6] = x4_r;
		my_val[7] = x4_i;

		return;	
	}


	//////////////////////////////////////////////
	////  (1)the common expression
	///////////////////////////////////////////////

	//////////////////////////////////////////////
	////  (2) the equation changes to 
	////	x^4 + bb*x^3 + cc*x^2 + dd*x + ee=0
	////
	////  (3) substitude x=y-aa/4, then we get
	////    y^4 + f*y^2 + g*y + h =0; where
	///////////////////////////////////////////

	double _f = my_in[2]/my_in[0] - 3*pow(my_in[1]/my_in[0],2)/8;
	double _g = my_in[3]/my_in[0] + (pow(my_in[1]/my_in[0],3)/8) - (my_in[1]/my_in[0]*my_in[2]/my_in[0]/2);
	double _h = my_in[4]/my_in[0] - (3*pow(my_in[1]/my_in[0],4)/256) + (pow(my_in[1]/my_in[0],2)*my_in[2]/my_in[0]/16) - (my_in[1]/my_in[0]*my_in[3]/my_in[0]/4);

	////////////////////////////////////////////////////////////////////////////
	///  (4) the normal situation is f, g and h are non-zero; then		////
	//// 	the related cubic equation is					////
	//// 	z^3 + (f/2) z^2 + ( (f^2-4h)/16 )*z -g^2/64 =0;			////
	////	 it has three "squared" roots, for example, p,q and l, 		////
	//// 	then p^2, q^2 and l^2 are the root of equation above		////
	//// 	set r=-g/8(pq), then the four roots of the original quartic are	////
	//// 		x = p + q + r -bb/4;					////
	//// 		x = p - q - r -bb/4;					////
	//// 		x = -p + q - r -bb/4;					////
	//// 		x = -p - q + r -bb/4;					////
	////////////////////////////////////////////////////////////////////////////

	double input[4]={1.0, _f/2, (pow(_f,2)-4*_h)/16., -1*pow(_g,2)/64.};
	cubic_kernel(input, real, img);

	////////////////////////////////////////////////
	////// (5) sqrt root of the cubic equation solutions
	////////////////////////////////////////////////

	double out_r[3]={0,0,0};
	double out_i[3]={0,0,0};
	bool _img[3]={false,false,false};
	int img_index =0;
	double my[2];

	for (int ii=0; ii<3; ii++){
		Csqrt_kernel(real[ii],img[ii],my);
		out_r[ii] = my[0];
		out_i[ii] = my[1];

		if ( my[1]!=0 ) {
			_img[ii]=true;
			img_index++;
		}
	}    


	///////////////////////////////////////////////
	////(x + yi)(u + vi) = (xu - yv) + (xv + yu)i
	//// calculating r = -g/(8pq)
	////////////////////////////////////////////////
	double r_r;	double r_i;
	double p_r;	double p_i;
	double q_r;	double q_i;

	if (_img[0]==_img[1] && out_r[0]*out_r[1]>0 ){
		r_r	= out_r[0]*out_r[1] - out_i[0]*out_i[1];
		r_i	= out_r[0]*out_i[1] + out_r[1]*out_i[0];
		p_r	= out_r[0]; p_i = out_i[0];
		q_r	= out_r[1]; q_i = out_i[1];
	} else {	
		for (int kk=0; kk<2; kk++) {
			for (int k=kk+1; k<3; k++) {
				if (_img[kk]==_img[k]){ 
					r_r	= out_r[kk]*out_r[k] - out_i[kk]*out_i[k];
					r_i	= out_r[kk]*out_i[k] + out_r[k]*out_i[kk];
					p_r	= out_r[kk]; p_i = out_i[kk];
					q_r	= out_r[k]; q_i = out_i[k];
				}
			}
		}
	}


	if( r_r !=0 )	{ 
		r_r 	= -1.*_g/8/r_r;
	} else {
		r_r 	= 0;
	}

	if( r_i >= 0.0000001 && r_i <= -0.0000001 )	{ 
		r_r 	= -1.*_g/8/r_i;
	} else {
		r_i 	= 0;
	}	

	x1_r = p_r + q_r + r_r -my_in[1]/my_in[0]/4;
	x1_i = p_i + q_i + r_i;
	x2_r = p_r - q_r - r_r -my_in[1]/my_in[0]/4;

	x1_r = p_r + q_r + r_r -my_in[1]/my_in[0]/4;
	x1_i = p_i + q_i + r_i;
	x2_r = p_r - q_r - r_r -my_in[1]/my_in[0]/4;
	x2_i = p_i - q_i - r_i;
	x3_r = -1*p_r + q_r - r_r -my_in[1]/my_in[0]/4;
	x3_i = -1*p_i + q_i - r_i;
	x4_r = -1*p_r - q_r + r_r -my_in[1]/my_in[0]/4;
	x4_i = -1*p_i - q_i + r_i;

	my_val[0] = x1_r;
	my_val[1] = x1_i;
	my_val[2] = x2_r;
	my_val[3] = x2_i;
	my_val[4] = x3_r;
	my_val[5] = x3_i;
	my_val[6] = x4_r;
	my_val[7] = x4_i;
}

	__device__
void Csqrt_kernel(double _ar, double _ai, double _my[])
{
	///// complex sqrt
	double x,y,r,w;
	if  (  (_ar == 0.0) && (_ai == 0.0) ) {
		_my[0]=0.0;
		_my[1]=0.0;

		return;
	} else {	
		x=fabs(_ar);
		y=fabs(_ai);

		if (x >= y) {
			r=y/x;
			w=sqrt(x)*sqrt(0.5*(1.0+sqrt(1.0+r*r)));
		} else {
			r=x/y;
			w=sqrt(y)*sqrt(0.5*(r+sqrt(1.0+r*r)));
		}

		if (_ar>= 0.0) {
			_my[0]=w;
			_my[1]=_ai/(2.0*w);
		} else {
			_my[1]=(_ai >= 0) ? w : -w;
			_my[0]=_ai/(2.0*_my[1]);
		}

		return;
	}
}

//////////////////////////////////////////////////////////////////
/// cubic /// a[0]x^3+a[1]x^2+a[2]x+a[3]=0
//////////////////////////////////////////////////////////////////
#define pi 3.1415926535897932384626433832795

__device__
void cubic_kernel(double a[], double rr[], double ri[])
{
	int i;
	double a1, a2;
	double g, y1, sh, theta;
	double z1, z2, z3, z4;
	//// initialize the results
	for (i = 0; i < 3; i ++)
	{
		rr[i] = 0.0;
		ri[i] = 0.0;
	}

	a1 = a[1]/3.0;
	a2 = a[2]/3.0;

	g = (a[0] * a[0]) * a[3] - 3.0 * a[0] * a1 * a2 + 2.0 * pow(a1, 3);
	y1 = g * g + 4.0 * pow(a[0] * a2 - a1 * a1, 3);

	if (y1 < 0.0){
		sh = sqrt(-(a[0] * a2 - a1 * a1));
		theta = acos(g / (2.0 * (a[0] * a2 - a1 * a1) * sh)) / 3.0;
		
		rr[0] = (2.0 * sh * cos(theta) - a1) / a[0];
		rr[1] = (2.0 * sh * cos(theta + (2.0 * pi / 3.0)) - a1) / a[0];
		rr[2] = (2.0 * sh * cos(theta + (4.0 * pi / 3.0)) - a1) / a[0];
		return;
	} else {
		z1 = (g + sqrt(y1)) / 2.0;
		z2 = (g - sqrt(y1)) / 2.0;
		if (z1 < 0.0){
			z3 = pow(-z1, 1.0/3.0);
			z3 = -z3;
		} else  z3 = pow(z1, 1.0/3.0);
		if (z2 < 0.0){
			z4 = pow(-z2, 1.0/3.0);
			z4 = - z4;
		}
		else  z4 = pow(z2, 1.0/3.0);

		rr[0] = -(a1 + z3 + z4) / a[0];
		rr[1] = (-2.0 * a1 + z3 + z4) / (2.0 * a[0]);
		ri[1] = sqrt(3.0) * (z4 - z3) / (2.0 * a[0]);
		rr[2] = rr[1];
		ri[2] = -ri[1];

		return;

	}
}
