
#include "neutgpu.cuh"
#include <cmath>
// defines para o CUDA

unsigned int GRID_SIZE;		// num de blocos por grid, vai dar 10^6
unsigned int BLOCK_SIZE;	// num threads por bloco

float total_M1, total_M2, total_K;

unsigned int NUM_THREADS; // = (1*8192*128);

__host__ 
void gpu_init(int blocks,int threads)	{
	GRID_SIZE=blocks;
	BLOCK_SIZE=threads;
	NUM_THREADS=blocks*threads;
}

__host__
void printIterations()	{
	printf("%d",NUM_THREADS);
}
/**
 *		Comparacoes de numeros em virgula flutuante sao problematicas, pelo que todas terao associadas
 * um erro relativo de 0.00000
 */

__device__
bool equal ( const double _this, const double &other) 	{
		double tmp = abs(_this-other);
		return ((tmp/abs(_this))< 0.000001);
}

__device__
bool nequal(const double _this, const double &other)	{
	return (!equal(_this,other));
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

// FILE* flog;

// NEUTRINO SOLUTIONS


	__global__
void dilep_kernel(double t_mass[], double w_mass[], double in_mpx[], 
//		double in_mpy[], double in_mpz[], LorentzVector *lep_a, 
		double in_mpy[], double in_mpz[], double lep_a[], 
//		LorentzVector *lep_b, LorentzVector *bl_a, LorentzVector *bl_b, 
		double lep_b[], double bl_a[], double bl_b[], 
		double nc[], int a[])
{

	///////////////////////////////////////////////////////////////////////////////////////////////////////
	// Neste codigo esta-se a assumir que LorentzVector.m contem o valor precalculado de LorentzVector->M()
	///////////////////////////////////////////////////////////////////////////////////////////////////////
	/*_LorentzVectorM(*lep_a);
	_LorentzVectorM(*lep_b);
	_LorentzVectorM(*bl_a);
	_LorentzVectorM(*bl_b);
*/
	unsigned tid = threadIdx.x + blockIdx.x * blockDim.x;
	double mpx, mpy, G_1, G_2, G_3, G_4;
	double WMass_a, WMass_b, tMass_a, tMass_b;

	mpx = in_mpx[0];
	mpy = in_mpy[0];

	WMass_a = w_mass[0];
	tMass_a = t_mass[0];
	WMass_b = w_mass[1];
	tMass_b = t_mass[1];  
	
	// Cancelamento subractivo
	/*
	G_1 = WMass_a*WMass_a - ( lep_a[4] )*( lep_a[4] );
	G_3 = WMass_b*WMass_b - ( lep_b[4] )*( lep_b[4] );
	G_2 = tMass_a*tMass_a - ( bl_a[4]  )*( bl_a[4]  ); 
	G_4 = tMass_b*tMass_b - ( bl_b[4]  )*( bl_b[4]  );
	*/
	
	G_1 = (WMass_a - lep_a[4]) * (WMass_a + lep_a[4]);
	G_3 = (WMass_b - lep_b[4]) * (WMass_b + lep_b[4]);
	G_2 = (tMass_a - bl_a[4]) * (tMass_a + bl_a[4]);
	G_4 = (tMass_b - bl_b[4]) * (tMass_b + bl_b[4]);

	double S=mpx;
	double T=mpy;

	double G_5,G_6,G_7,G_8,G_9,G_10,G_11,G_12;
	G_5 = ( bl_a[0]/bl_a[3] - lep_a[0]/lep_a[3] );
	G_6 = ( bl_a[1]/bl_a[3] - lep_a[1]/lep_a[3] );
	G_7 = ( bl_a[2]/bl_a[3] - lep_a[2]/lep_a[3] );
	G_8 = ( G_1/lep_a[3] - G_2/bl_a[3] )/2.;

	G_9 =	( bl_b[0]/bl_b[3] - lep_b[0]/lep_b[3] );
	G_10 =	( bl_b[1]/bl_b[3] - lep_b[1]/lep_b[3] );
	G_11 =	( bl_b[2]/bl_b[3] - lep_b[2]/lep_b[3] );
	G_12 =	( G_3/lep_b[3] - G_4/bl_b[3] )/2.;

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
	
	double in_a[5],on_a[3],out_a[6];
	in_a[0] = G_8/G_7;
	in_a[1] = -1.0*G_5/G_7;
	in_a[2] = -1.0*G_6/G_7;
	in_a[3] = lep_a[3];
	in_a[4] = G_1;

	on_a[0] = lep_a[0];
	on_a[1] = lep_a[1];
	on_a[2] = lep_a[2];
	toz_kernel(in_a, on_a, out_a);

	double in_c[5],on_c[3],out_c[6];
	in_c[0] = G_12/G_11;
	in_c[1] = -1*G_9/G_11;
	in_c[2] = -1*G_10/G_11;
	in_c[3] = lep_b[3];
	in_c[4] = G_3;

	on_c[0] = lep_b[0];
	on_c[1] = lep_b[1];
	on_c[2] = lep_b[2];
	toz_kernel(in_c, on_c, out_c);
	/////////////////////////////////////////////////////
	//////change x2 y2 equation to x1 and y1 by using
	////// 		x1+x2 = S = mpx
	////// 		y1+y2 = T = mpy
	/////////////////////////////////////////////////////
	double out_e[6];
	out_e[0] = out_c[0];
	out_e[1] = out_c[1];
	out_e[2] = -1*( out_c[0]*mpx + out_c[2] + out_c[4]*mpy);
	out_e[3] = -1*( out_c[1]*mpy + out_c[3] + out_c[4]*mpx);
	out_e[4] = out_c[4]; 
	out_e[5] =( out_c[0]*mpx*mpx + out_c[1]*mpy*mpy + 2*out_c[2]*mpx + 2*out_c[3]*mpy + out_c[5] + 2*out_c[4]*mpx*mpy);

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
	
	double _A = out_e[0];
	double _B = out_e[1];
	double _D = out_e[2];
	double _E = out_e[3];
	double _F = out_e[4];
	double _G = out_e[5];
	
	double _d, _za, _f, _b, _e, _g;
	double fx_1, fx_2, fx_3, fx_4, fx_5;
	double k_1, k_2, k_3, k_4, k_5;
/*
	//if ( out_a[0]!=0  ){
	_a = ( out_a[0]!=0  ) * out_a[0];
	_b = ( out_a[0]!=0  ) * out_a[1];
	_d = ( out_a[0]!=0  ) * out_a[2];
	_e = ( out_a[0]!=0  ) * out_a[3];
	_f = ( out_a[0]!=0  ) * out_a[4];
	_g = ( out_a[0]!=0  ) * out_a[5];

	//if (out_a[0]==0 && _A!=0 ) {
	_a = (out_a[0]==0 && _A!=0) * 1.;
	_b = (out_a[0]==0 && _A!=0) * (out_a[1] + _B/_A);
	_d = (out_a[0]==0 && _A!=0) * (out_a[2] + _D/_A);
	_e = (out_a[0]==0 && _A!=0) * (out_a[3] + _E/_A);
	_f = (out_a[0]==0 && _A!=0) * (out_a[4] + _F/_A);
	_g = (out_a[0]==0 && _A!=0) * (out_a[5] + _G/_A);*/

	if ( out_a[0]!=0  ){
		_za = out_a[0];
		_b = out_a[1];
		_d = out_a[2];
		_e = out_a[3];
		_f = out_a[4];
		_g = out_a[5];

		/* badalhoquice!!! */
		fx_1 = 2.*(_A*_d - _za*_D);
		fx_2 = 2.*(_A*_f - _za*_F);
		fx_3 = _za*_B - _A*_b;
		fx_4 = 2.*(_E*_za - _A*_e);
		fx_5 = _za*_G - _g*_A;

		k_1 = ( _f*_f - _za*_b )/_za/_za;
		k_2 = ( 2.*_d*_f - 2.*_za*_e )/_za/_za;
		k_3 = ( _d*_d-_za*_g )/_za/_za;
		k_4 = -1.*_d/_za;
		k_5 = -1.*_f/_za;
	} else {
		if (out_a[0]==0 && _A!=0 ) {
			_za = 1.;
			_b = out_a[1] + _B/_A;
			_d = out_a[2] + _D/_A;
			_e = out_a[3] + _E/_A;
			_f = out_a[4] + _F/_A;
			_g = out_a[5] + _G/_A;

		/* badalhoquice!!! */
			fx_1 = 2.*(_A*_d - _za*_D);
			fx_2 = 2.*(_A*_f - _za*_F);
			fx_3 = _za*_B - _A*_b;
			fx_4 = 2.*(_E*_za - _A*_e);
			fx_5 = _za*_G - _g*_A;

			k_1 = ( _f*_f - _za*_b )/_za/_za;
			k_2 = ( 2.*_d*_f - 2.*_za*_e )/_za/_za;
			k_3 = ( _d*_d-_za*_g )/_za/_za;
			k_4 = -1.*_d/_za;
			k_5 = -1.*_f/_za;
		}
	}

	if ( out_a[0]==0 && _A==0){
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


	double g_1 = 4.*_A*_A*k_5*k_5 + 4.*_F*_F + 8.*_A*_F*k_5;
	double m_1 = g_1*k_1;
	double g_2 = 8.*_A*_A*k_4*k_5 + 8.*_A*_D*k_5 + 8.*_A*_F*k_4 + 8.*_D*_F;
	double g_3 = 4.*_A*_A*k_4*k_4 + 4.*_D*_D + 8.*_A*_D*k_4;
	double g_4 = _A*k_1 + _A*k_5*k_5;
	double g_5 = _A*k_2 + 2.*_A*k_4*k_5 + 2.*_D*k_5;
	double g_6 = _A*k_3 + _A*k_4*k_4 + 2.*_D*k_4 + _G;

	//double m_1 = g_1*k_1;
	double m_2 = g_1*k_2 + g_2*k_1;
	double m_3 = g_1*k_3 + g_2*k_2 + g_3*k_1;
	double m_4 = g_2*k_3 + g_3*k_2;
	double m_5 = g_3*k_3;

	double m_6  = _B*_B + 4.*_F*_F*k_5*k_5 + 4.*_B*_F*k_5;
	double m_7  = 4.*_B*_E + 8.*_F*_F*k_4*k_5 + 4.*_B*_F*k_4 + 8.*_E*_F*k_5;
	double m_8  = 4.*_E*_E + 4.*_F*_F*k_4*k_4 + 8.*_E*_F*k_4;
	double m_80 = pow(g_4,2);
	double m_81 = 2*g_4*g_5;
	double m_9  = pow(g_5,2) + 2.*g_4*g_6;
	double m_10 = 2.*g_5*g_6;
	double m_11 = g_6*g_6;

	double m_12 = 	2.*_A*_B*k_1 + 2.*_A*_B*k_5*k_5 + 4.*_A*_F*k_1*k_5 + 4.*_A*_F*pow(k_5,3);
	double m_13 = 	2.*_A*_B*k_2 + 4.*_A*_B*k_4*k_5 + 4.*_B*_D*k_5 + 
		4.*_A*(_E*k_1 + _E*k_5*k_5 + _F*k_1*k_4 + _F*k_2*k_5) + 
		12.*_A*_F*k_4*k_5*k_5 + 8.*_D*_F*k_5*k_5;
	double m_14 = 	2.*_A*_B*k_3 + 2.*_A*_B*k_4*k_4 + 4.*_D*_B*k_4 + 2.*_B*_G + 4.*_A*_E*k_2 + 
		8.*_A*_E*k_4*k_5 + 8.*_E*_D*k_5 + 4.*_A*_F*k_2*k_4 + 4.*_A*_F*k_3*k_5 + 
		12.*_A*_F*k_4*k_4*k_5 + 16.*_D*_F*k_4*k_5 + 4.*_F*_G*k_5;
	double m_15 = 	4.*_A*_E*(k_3 + k_4*k_4) + 8.*_E*_D*k_4 + 4.*_E*_G + 4.*_A*_F*(k_3*k_4 + pow(k_4,3)) + 
		8.*_D*_F*k_4*k_4 + 4.*_F*_G*k_4;

	double  re[5];
	re[0] = m_1 - m_6 - m_12 - m_80;
	re[1] = m_2 - m_7 - m_13 - m_81;
	re[2] = m_3 - m_8 - m_9 - m_14;
	re[3] = m_4 - m_10 - m_15;
	re[4] = m_5 - m_11;  



	/////////////////////////////////////////////////
	double output[8];
	my_qu_kernel(re,output);

	//int sign[2] = {-1, +1}; // NEVER USED
	int ncand(0);

	// deve haver uma solucao mais decente... (devia ser neut..= new std::..
	// mas o new nao funciona em CUDA)
	//std::vector<LorentzVector> n_a;
	//std::vector<LorentzVector> n_b;
	// nao e permitido fazer isto pelo device...
	//vector<myvector> n_c ();
	//vector<LorentzVector> *neutrinoContainer_a = &n_a;
	//vector<LorentzVector> *neutrinoContainer_b = &n_b;

	//////////////////////////////////////////////////
	//// please replace your code starting from here 
	//// ////////////////////////////////////////////


	double rec_x1, rec_y1, rec_z1, rec_e1, rec_x2, rec_y2, rec_z2, rec_e2;
	double m_delta_mass;

	for (int j=0; j<8; j+=2){
		double delta = k_1*output[j]*output[j] + k_2*output[j] + k_3;
//		if ( equal(output[j+1],0) && delta>=0) {
		if ( output[j+1]==0 && delta >=0) {
		//if ( output[j+1]==0 && (k_1*output[j]*output[j] + k_2*output[j] + k_3) >=0) {
			if ( (fx_1 + fx_2*output[j])!=0 ) {
//			if ( nequal((fx_1 + fx_2*output[j]),0) ) {
				rec_x1 = (fx_3*pow(output[j],2) + fx_4*output[j] + fx_5)/(fx_1 + fx_2*output[j]);
			} else {
				rec_x1 = sqrt(delta)+k_4+k_5*output[j];
			}  

			rec_y1 = output[j];
			rec_z1 = G_8/G_7 - G_5*rec_x1/G_7 - G_6*rec_y1/G_7;
			rec_e1 = sqrt(rec_x1*rec_x1 + rec_y1*rec_y1 + rec_z1*rec_z1);
			rec_x2 = S - rec_x1;
			rec_y2 = T - rec_y1;
			rec_z2 = G_12/G_11 - G_9*rec_x2/G_11 - G_10*rec_y2/G_11;
			rec_e2 = sqrt(rec_x2*rec_x2 + rec_y2*rec_y2 + rec_z2*rec_z2);
			
			// self-consistence check and control of the solutions
			double m_w11 = calcMass(rec_x1+lep_a[0], rec_y1+lep_a[1], rec_z1+lep_a[2], rec_e1+lep_a[3]);
			double m_w12 = calcMass(rec_x2+lep_b[0], rec_y2+lep_b[1], rec_z2+lep_b[2], rec_e2+lep_b[3]);
			double m_t11 = calcMass(rec_x1+ bl_a[0], rec_y1+ bl_a[1], rec_z1+ bl_a[2], rec_e1+ bl_a[3]);
			double m_t12 = calcMass(rec_x2+ bl_b[0], rec_y2+ bl_b[1], rec_z2+ bl_b[2], rec_e2+ bl_b[3]);

			m_delta_mass = 1000.0; // allow mass variation range for reco W and tops..
			
			// Avoid thread divergence
			bool m_good_eq1 = ( fabs(S -(rec_x1+rec_x2)) <= 0.01 ) * true + 
							  ( fabs(S -(rec_x1+rec_x2)) > 0.01 ) * false;
			bool m_good_eq2 = ( fabs(T -(rec_y1+rec_y2)) <= 0.01 ) * true +
							  ( fabs(T -(rec_y1+rec_y2)) > 0.01 ) * false;
			bool m_good_eq3 = ( fabs(m_w11 - w_mass[0]) <= m_delta_mass ) * true + 
							  ( fabs(m_w11 - w_mass[0]) > m_delta_mass ) * false;
			bool m_good_eq4 = ( fabs(m_w12 - w_mass[1]) <= m_delta_mass ) * true +
							  ( fabs(m_w12 - w_mass[1]) > m_delta_mass ) * false;
			bool m_good_eq5 = ( fabs(m_t11 - t_mass[0]) <= m_delta_mass ) * true +
							  ( fabs(m_t11 - t_mass[0]) > m_delta_mass ) * false;
			bool m_good_eq6 = ( fabs(m_t12 - t_mass[1]) <= m_delta_mass ) * true +
							  ( fabs(m_t12 - t_mass[1]) <= m_delta_mass ) * false;

			int aux = 2*j;
			bool cond = m_good_eq1 && m_good_eq2 && m_good_eq3 && m_good_eq4 && m_good_eq5 && m_good_eq6;
			syncthreads();
			nc[tid * 16 + aux] = cond * rec_x1;
			nc[tid * 16 + aux + 1] = cond * rec_y1;
			nc[tid * 16 + aux + 2] = cond * rec_z1;
			nc[tid * 16 + aux + 3] = cond * rec_z2;
			ncand += cond * 1;
		}
	}/*
	nc[0] = output[0];
	nc[1] = output[1];
	nc[2] = output[2];
	nc[3] = output[3];
	nc[4] = output[4];
	nc[5] = output[5];
	nc[6] = output[6];
	nc[7] = output[7];*/
	a[tid] = ncand;
	//a[0] = -5;
}  //// end of function

//__host__
//double inline TO1D (double nc[], int tid, int sol, int did)	{
//	return nc[tid*16+sol*4+did];
//}

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
	// aloca a memoria no GPU (falta tratar os erros e verificar a capacidade do GPU)
	double *dev_t_mass, *dev_w_mass, *dev_in_mpx, *dev_in_mpy, *dev_in_mpz;
	
//	double ax = lep_a->Px();
//	double ay = lep_a->Py();
//	double az = lep_a->Pz();
//	double ae = lep_a->E();
//	double am = lep_a->M();

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
	//resultfinal = new std::vector<myvector>();

	//choseDevice();

// Prepare
	cudaEvent_t start, stop;
	float elapsedTime;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	// Start record
	cudaEventRecord(start, 0);

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

	// transfere os dados para a memoria do GPU
	cudaMemcpy(dev_t_mass, t_mass, 2*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_w_mass, w_mass, 2*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_in_mpx, in_mpx, 2*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_in_mpy, in_mpy, 2*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_in_mpz, in_mpz, 2*sizeof(double), cudaMemcpyHostToDevice);

	cudaMemcpy(dev_lep_a, &a,	5*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_lep_b, &b,	5*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_bl_a, &c,	5*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_bl_b, &d,	5*sizeof(double), cudaMemcpyHostToDevice);

	cudaMemcpy(dev_nc, nc, 16*NUM_THREADS*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_count, count,NUM_THREADS*sizeof(int), cudaMemcpyHostToDevice);

	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&elapsedTime, start, stop);
	total_M1 += elapsedTime;
	//	cout << "Grid " << GRID_SIZE << endl << "Block " << BLOCK_SIZE << endl; 
	// executa a dilep_kernel - 1 dimensao para ser mais facil trabalhar com o array
	dim3 dimGrid(GRID_SIZE, 1);
	dim3 dimBlock(BLOCK_SIZE, 1);

	cudaEventRecord(start, 0);
	// todos os blocos fazem o mesmo... 
	dilep_kernel<<<dimGrid,dimBlock>>>(
			dev_t_mass, dev_w_mass, dev_in_mpx, dev_in_mpy, dev_in_mpz, 
			dev_lep_a, dev_lep_b, dev_bl_a, dev_bl_b, dev_nc, dev_count);


	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&elapsedTime, start, stop);
	total_K += elapsedTime;

	cudaEventRecord(start, 0);
	// recebe o resultado, nao sei se o nc tambem tem de ser apenas 
	cudaMemcpy(nc, dev_nc, 16*NUM_THREADS*sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(count, dev_count,NUM_THREADS*sizeof(int), cudaMemcpyDeviceToHost);

	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&elapsedTime, start, stop);
	total_M2 += elapsedTime;
//	cout << "LOL " << count[0] << endl;
	unsigned int aux_size = GRID_SIZE*BLOCK_SIZE*16, is;

	cudaEventDestroy(start);
	cudaEventDestroy(stop);
//	cout << resultfinal->size() << " MY_WEIRD_SIZE = " << aux_size;

//	fprintf(flog, "aux_size: %u\n", aux_size);
	for(int thread = 0 ; thread < GRID_SIZE*BLOCK_SIZE ; thread++)	{
		for( int sol = 0 ; sol < count[thread] && sol<4 ; sol++)	{
			


			myvector *mv = new myvector( 
				TO1D(nc,thread,sol,0),
				TO1D(nc,thread,sol,1),
				TO1D(nc,thread,sol,2),
				TO1D(nc,thread,sol,3) );
			
//			cout << endl << "rec_x1: "   << TO1D(nc,thread,sol,0)
//				<< "\trec_y1: " << TO1D(nc,thread,sol,1)
//				<< "\trec_z1: " << TO1D(nc,thread,sol,2)
//				<< "\trec_z2: " << TO1D(nc,thread,sol,3)  << endl;
			
			resultfinal->push_back(*mv);
		}
	}
	
	// i aponta sempre para o inicio do resultado de cada kernel
	// count[i] indica o numero de solucoes desse kernel, cada uma com 4 doubles
/*	for(is = 0; is < aux_size; is += 16){
		for(int aux_j = 0; aux_j < count[is / 16]*4; aux_j += 4){
			
			cout << endl << "rec_x1:" << nc[is+aux_j] << "\trec_y1:" << nc[is+aux_j+1] << " \trec_z1:" << nc[is+aux_j+2] << "\trec_z2" << nc[is+aux_j+3] << endl;
			
			myvector *mv = new myvector(nc[is+aux_j], nc[is+aux_j+1],
					nc[is+aux_j+2], nc[is+aux_j+3]);
			resultfinal->push_back(*mv);
		}
		fprintf(flog,"count: %d\n", count[is/16]);
	}
	fprintf(flog, "sols: %d\n", resultfinal->size());
*/	// faz free das variaveis
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

// CUDA version must be higher than 2.0 (nao funciona na versao das 8800gt)
// Enough memory to run the kernels
// Highest number of SMPs
// Maximum of 8 GPUs

__host__
int choseDevice(){
	int count, device = -1,  numSMP = 0, sel[8], flag = 0, flag2;
	cudaDeviceProp properties[8];
	// assumir doubles - ainda nao sei aquilo dos floats
	long int globMem = (30 + 16 * BLOCK_SIZE*GRID_SIZE + BLOCK_SIZE*GRID_SIZE)*8;

	cudaGetDeviceCount(&count);
	
	for(int i = 0; i < count && i < 8; i++){

		cudaGetDeviceProperties( &properties[i], i);
		
			printf("DEV %s\n", properties[i].name);
		if(strcmp(properties[i].name,"Tesla C2050")==0)	{
			cudaSetDevice(i);
			return 0;
		}

		continue;
		
		if(properties[i].major >= 2){
			sel[i] = 1;
			flag++;
		}
		else
			sel[i] = 0;
		
		if(flag == 0){
	//		fprintf(flog, "There is no GPUs capable of running this program on your system.\n");
			return device;
		}

	//	fprintf(flog, "\nNumber of CUDA capable devices: %d\n", count);

		flag = 1;
		for(i = 0; i < count && i < 8; i++){
			if(sel[i]){
				if(properties[i].totalGlobalMem > globMem && properties[i].totalGlobalMem < 8000000000){
					sel[i] += flag;
					flag++;
					globMem = properties[i].totalGlobalMem;
				}
			}
		}
		
		flag2 = 1;

		for(i = 0; i < count && i < 8; i++){
			if(sel[i] == flag){
				if(properties[i].multiProcessorCount > numSMP){
					sel[i] += flag2;
					flag2++;
					numSMP = properties[i].multiProcessorCount;
				}
			}
		}

		for(i = 0; i < count && i < 8; i++){
			if(sel[i] == flag + flag2 - 1){
				cudaChooseDevice(&i, &properties[i]);
				cudaSetDevice(i);
			//	fprintf(flog, "\n%s GPU chosen\n", properties[i].name);
			//	fprintf(flog, "ID: %d\n", i);
			//	fprintf(flog, "Global Memory: %ld bytes\n", globMem);
			//	fprintf(flog, "Number of SMP: %d\n\n", numSMP);
				return i;
			}
		}
	}
	// Se nao houver nenhuma "Tesla C2050" retorna 1
	return 1;

	return device;
}

//////////////////////////////////////
__device__
void toz_kernel(double k[], double l[], double g[]){
	//// checked !!
	///////////////////////////////////////////////////////////////////////////
	///// bring z=A+Bx+Cy to 2*D*sqrt(x**2+y**2+z**2)-2(ax+by+dz) = E
	///// simplify it to g1*x^2 + g2*y^2 + 2*g3*x + 2*g4*y + 2*g5*xy + g6 = 0
	///////////////////////////////////////////////////////////////////////////
	double A = k[0];
	double B = k[1];
	double C = k[2];
	double D = k[3];
	double E = k[4];
	double a = l[0];
	double b = l[1];
	double d = l[2];
	g[0] = 4*pow(D,2)*( 1 + pow(B,2)) - 4*pow(a,2) -4*pow(d,2)*pow(B,2) - 8*a*d*B;
	if ( g[0]!=0 ) { 
//	if ( nequal(g[0],0) ) { 
		g[1] = ( 4*pow(D,2)*( 1 + pow(C,2)) - 4*pow(b,2) -4*pow(d,2)*pow(C,2) - 8*b*d*C )/g[0] ;
		g[2] = ( (4*pow(D,2)-4*d*d)*A*B - 4*a*d*A - 2*E*a - 2*E*d*B )/g[0];
		g[3] = ( (4*pow(D,2)-4*d*d)*A*C - 4*b*d*A - 2*E*b - 2*E*d*C )/g[0];
		g[4] = ( (4*pow(D,2)-4*d*d)*B*C - 4*a*b   - 4*a*d*C - 4*b*d*B )/g[0];
		g[5] = ( (4*pow(D,2)-4*d*d)*A*A - E*E - 4*E*d*A )/g[0];
		g[0] = 1.0; 
	} else {
		g[1] = ( 4*pow(D,2)*( 1 + pow(C,2)) - 4*pow(b,2) -4*pow(d,2)*pow(C,2) - 8*b*d*C ) ;
		g[2] = ( (4*pow(D,2)-4*d*d)*A*B - 4*a*d*A - 2*E*a - 2*E*d*B );
		g[3] = ( (4*pow(D,2)-4*d*d)*A*C - 4*b*d*A - 2*E*b - 2*E*d*C );
		g[4] = ( (4*pow(D,2)-4*d*d)*B*C - 4*a*b   - 4*a*d*C - 4*b*d*B );
		g[5] = ( (4*pow(D,2)-4*d*d)*A*A - E*E - 4*E*d*A );
		g[0] = 0.;  
	}	 
	return;
}


///////////////////////////////////////////
	__device__
void my_qu_kernel( double my_in[], double my_val[])
{

	///////////////////////////////////////////
	/////  find the solution /////////////////
	/////  ax^4+bx^3+cx^2+dx+e=0
	//////////////////////////////////////////
	double a=my_in[0];
	double b=my_in[1];
	double c=my_in[2];
	double d=my_in[3];
	double e=my_in[4];

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

	//void cubic(double a[4], double rr[3], double ri[3]);
	if ( a ==0 && b!=0 ){
//	if ( equal(a, 0) && nequal(b,0) ){
		double input[4]={b,c,d,e};
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
		// printf(" the 1st cubic root is: %3.5f %5c %3.5f \n", x1_r, "+i*", x1_i);
		// printf(" the 2nd cubic root is: %3.5f %5c %3.5f \n", x2_r, "+i*", x2_i);
		// printf(" the 3rd cubic root is: %3.5f %5c %3.5f \n", x3_r, "+i*", x3_i);
		return;
	}


	/////////////////////////////////////////////
	///// in case of a==0 && b==0, 
	///// simplify to quadratic
	///// cx*2 + d*x + e =0
	/////////////////////////////////////////////

	if ( a ==0 && b==0 && c!=0){
//	if ( equal(a, 0) && equal(b,0) && nequal(c,0)){
		//double input[3]={c,d,e}; // NEVER USED
		double alpha = pow(d,2)-4*c*e;
		if (alpha>=0) {
			x1_r = (-1*d + sqrt(alpha))/2/c;	x1_i = 0;
			x2_r = (-1*d - sqrt(alpha))/2/c;	x2_i = 0;
		} else {
			x1_r = (-1*d )/2/c;	x1_i = sqrt(-alpha)/2/c;
			x2_r = (-1*d )/2/c;	x2_i = -1*sqrt(-alpha)/2/c;  	
		}
		// printf(" the 1st quadratic root is: %3.5f %5c %3.5f \n", x1_r, "+i*", x1_i);
		// printf(" the 2nd quadratic root is: %3.5f %5c %3.5f \n", x2_r, "+i*", x2_i);
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

	if (a ==0 && b==0 && c==0 && d!=0){
//	if (equal(a,0) && equal(b,0) && equal(c,0) && nequal(d,0)){
		x1_r = -e/d;	x1_i = 0;
		my_val[0] = x1_r;
		my_val[1] = x1_i;
		my_val[2] = x2_r;
		my_val[3] = x2_i;
		my_val[4] = x3_r;
		my_val[5] = x3_i;
		my_val[6] = x4_r;
		my_val[7] = x4_i;
		// printf(" the answer to linear equation is: %3.5f %5c %3.5f \n", x1_r, "+i*", x1_i);
		return;	
	}


	//////////////////////////////////////////////
	////  (1)the common expression
	///////////////////////////////////////////////
	// double aa=1; // NEVER USED
	double bb=b/a; double cc=c/a;
	double dd=d/a; double ee=e/a;

	//////////////////////////////////////////////
	////  (2) the equation changes to 
	////	x^4 + bb*x^3 + cc*x^2 + dd*x + ee=0
	////
	////  (3) substitude x=y-aa/4, then we get
	////    y^4 + f*y^2 + g*y + h =0; where
	///////////////////////////////////////////

	double _f = cc - 3*pow(bb,2)/8;
	double _g = dd + (pow(bb,3)/8) - (bb*cc/2);
	double _h = ee - (3*pow(bb,4)/256) + (pow(bb,2)*cc/16) - (bb*dd/4);
	//printf(" _f , _g and _h: %3.5f %3.5f %3.5f \n", _f, _g, _h);

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

	double c_1 = 1.;
	double c_2 = _f/2;
	double c_3 = (pow(_f,2)-4*_h)/16.;
	double c_4 = -1*pow(_g,2)/64.;
	//printf(" c_2, c_3 and c_4: %3.5f %3.5f %3.5f\n\n", c_2,c_3,c_4); 

	double input[4]={c_1,c_2,c_3,c_4};
	cubic_kernel(input,real,img);
	//printf(" real part, 1/2/3: %3.5f %3.5f %3.5f \n",real[0],real[1],real[2]);
	//printf(" img part, 1/2/3: %3.5f %3.5f %3.5f \n",img[0],img[1],img[2]);

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
//		if ( nequal(my[1],0) ) {
			_img[ii]=true;
			img_index++;
		}
		//    cout <<" my "<<out_r[ii]<<" "<<out_i[ii]<< " "<< _img[ii]<<endl;
	}    


	///////////////////////////////////////////////
	////(x + yi)(u + vi) = (xu - yv) + (xv + yu)i
	//// calculating r = -g/(8pq)
	////////////////////////////////////////////////
	double r_r;	double r_i;
	double p_r;	double p_i;
	double q_r;	double q_i;

	if (_img[0]==_img[1] && out_r[0]*out_r[1]>0 ){
//	if (equal(_img[0],_img[1]) && out_r[0]*out_r[1]>0 ){
		r_r	= out_r[0]*out_r[1] - out_i[0]*out_i[1];
		r_i	= out_r[0]*out_i[1] + out_r[1]*out_i[0];
		p_r	= out_r[0]; p_i = out_i[0];
		q_r	= out_r[1]; q_i = out_i[1];
	} else {	
		for (int kk=0; kk<2; kk++) {
			for (int k=kk+1; k<3; k++) {
				if (_img[kk]==_img[k]){ 
//				if (equal(_img[kk],_img[k])){ 
					r_r	= out_r[kk]*out_r[k] - out_i[kk]*out_i[k];
					r_i	= out_r[kk]*out_i[k] + out_r[k]*out_i[kk];
					p_r	= out_r[kk]; p_i = out_i[kk];
					q_r	= out_r[k]; q_i = out_i[k];
					//	cout <<  "taking "<< kk+1 << " and "<< k+1 <<endl;
					//	cout << " kk/k "<< q_r<< " " <<q_i<<" "<<p_r << " "<<p_i<<endl;
				}
			}
		}
	}


	if( r_r !=0 )	{ 
//	if( nequal(r_r, 0))	{ 
		r_r 	= -1.*_g/8/r_r;
	} else {
		r_r 	= 0;
	}
//	if( r_i !=0 )	{ 
//	if( nequal(r_i, 0))	{ 
	if( r_i >= 0.0000001 && r_i <= -0.0000001 )	{ 
		r_r 	= -1.*_g/8/r_i;
	} else {
		r_i 	= 0;
	}	

	//cout << "  rr "<<r_r <<" " <<r_i<<endl;

	x1_r = p_r + q_r + r_r -b/a/4;
	x1_i = p_i + q_i + r_i;
	x2_r = p_r - q_r - r_r -b/a/4;

	x1_r = p_r + q_r + r_r -b/a/4;
	x1_i = p_i + q_i + r_i;
	x2_r = p_r - q_r - r_r -b/a/4;
	x2_i = p_i - q_i - r_i;
	x3_r = -1*p_r + q_r - r_r -b/a/4;
	x3_i = -1*p_i + q_i - r_i;
	x4_r = -1*p_r - q_r + r_r -b/a/4;
	x4_i = -1*p_i - q_i + r_i;

	my_val[0] = x1_r;
	my_val[1] = x1_i;
	my_val[2] = x2_r;
	my_val[3] = x2_i;
	my_val[4] = x3_r;
	my_val[5] = x3_i;
	my_val[6] = x4_r;
	my_val[7] = x4_i;
	/*
	   cout << " 1st root is: "<<x1_r <<" + i*"<<x1_i<<endl;
	   cout << " 2nd root is: "<<x2_r <<" + i*"<<x2_i<<endl;
	   cout << " 3rd root is: "<<x3_r <<" + i*"<<x3_i<<endl;
	   cout << " 4th root is: "<<x4_r <<" + i*"<<x4_i<<endl;
	 */
}
////////////////////end of main
///////////////////////////////////////////////////////////////
////+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	__device__
void Csqrt_kernel(double _ar, double _ai, double _my[])
{
	///// complex sqrt
	//cout << " coming here "<<_ar<<" " <<_ai<<endl; 
	double x,y,r,w;
	if  (  (_ar == 0.0) && (_ai == 0.0) ) {
//	if  (  (equal(_ar , 0.0)) && (equal(_ai , 0.0)) ) {
		_my[0]=0.0;
		_my[1]=0.0;
		return;
	} else {	
		x=fabs(_ar);
		y=fabs(_ai);
		//		cout << x <<" "<<y<<endl;
		if (x >= y) {
			r=y/x;
			w=sqrt(x)*sqrt(0.5*(1.0+sqrt(1.0+r*r)));
		} else {
			r=x/y;
			w=sqrt(y)*sqrt(0.5*(r+sqrt(1.0+r*r)));
		}
		//	cout <<"r/w "<<r<<" "<<w<<endl;
		if (_ar>= 0.0) {
			_my[0]=w;
			_my[1]=_ai/(2.0*w);
		} else {
			_my[1]=(_ai >= 0) ? w : -w;
			_my[0]=_ai/(2.0*_my[1]);
		}
		//		cout << _my[0] <<" " <<_my[1]<<endl;
		return;
	}
}

//////////////////////////////////////////////////////////////////
/// cubic /// a[0]x^3+a[1]x^2+a[2]x+a[3]=0
//////////////////////////////////////////////////////////////////
	__device__
void cubic_kernel(double a[], double rr[], double ri[])
{
	int i;
	double a0, a1, a2, a3;
	double g, h, y1, sh, theta, pi, xy1, xy2, xy3;
	double y2, z1, z2, z3, z4;
	//// initialize the results
	for (i = 0; i < 3; i ++)
	{
		rr[i] = 0.0;
		ri[i] = 0.0;
	}

	a0 = a[0];
	a1 = a[1]/3.0;
	a2 = a[2]/3.0;
	a3 = a[3];

	g = (a0 * a0) * a3 - 3.0 * a0 * a1 * a2 + 2.0 * pow(a1, 3);
	h = a0 * a2 - a1 * a1;
	y1 = g * g + 4.0 * pow(h, 3);

	if (y1 < 0.0){
		sh = sqrt(-h);
	//	xk = 2.0 * sh;
		theta = acos(g / (2.0 * h * sh)) / 3.0;
		xy1 = 2.0 * sh * cos(theta);
		//                pi in double precision
		pi = 3.1415926535897932384626433832795;
		xy2 = 2.0 * sh * cos(theta + (2.0 * pi / 3.0));
		xy3 = 2.0 * sh * cos(theta + (4.0 * pi / 3.0));
		rr[0] = (xy1 - a1) / a0;
		rr[1] = (xy2 - a1) / a0;
		rr[2] = (xy3 - a1) / a0;
		return;
	} else {
		y2 = sqrt(y1);
		z1 = (g + y2) / 2.0;
		z2 = (g - y2) / 2.0;
		if (z1 < 0.0){
			z3 = pow(-z1, 1.0/3.0);
			z3 = -z3;
		} else  z3 = pow(z1, 1.0/3.0);
		if (z2 < 0.0){
			z4 = pow(-z2, 1.0/3.0);
			z4 = - z4;
		}
		else  z4 = pow(z2, 1.0/3.0);

		rr[0] = -(a1 + z3 + z4) / a0;
		rr[1] = (-2.0 * a1 + z3 + z4) / (2.0 * a0);
		ri[1] = sqrt(3.0) * (z4 - z3) / (2.0 * a0);
		rr[2] = rr[1];
		ri[2] = -ri[1];

		return;

	}
}

__device__
 //void _LorentzVectorM(LorentzVector &v);
 double _LorentzVectorM(double x, double y, double z, double e)  {
     double mm=e*e-(x*x+y*y+z*z);
     double res;
	res = sqrt(-mm);

	res *= (mm < 0.0) * -1.0;

	return res;
}

