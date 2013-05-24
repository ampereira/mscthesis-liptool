#include "neut.h"

using std::vector;
using namespace std;

namespace Dilep {
	namespace GPU {

		unsigned _tid = 0, size;
	
		double calcMass (double x, double y, double z, double e) {
			double mm, mass;

			mm = e*e - (x*x + y*y + z*z);

			if (mm < 0.0)
				mass = -sqrt(-mm);
			else
				mass = sqrt(mm);

			return mass;
		}

		double calcMass (double array[]) {
			double mm = array[3]*array[3] - (array[0]*array[0] + array[1]*array[1] + array[2]*array[2]);

			if (mm < 0.0)
				array[4] = -sqrt(-mm);
			else
				array[4] = sqrt(mm);
		}

		// Wrapper for the dilep calculation using the input class
		void dilep (DilepInput &di) {
			std::vector<myvector> *result = new std::vector<myvector> ();
			int hasSolution = 0;

			// time measurement
			#ifdef MEASURE_DILEP
			long long int time = startTimer();
			#endif

			double in_mpx[2], in_mpy[2], in_mpz[2], t_mass[2], w_mass[2];
			TLorentzVector lep_a, lep_b, bl_a, bl_b;

			// Transferring the inputs to local variables
			in_mpx[0] = di.getInMpx(0);
			in_mpx[1] = di.getInMpx(1);
			in_mpy[0] = di.getInMpy(0);
			in_mpy[1] = di.getInMpy(1);
			in_mpz[0] = di.getInMpz(0);
			in_mpz[1] = di.getInMpz(1);
			t_mass[0] = di.getTmass(0);
			t_mass[1] = di.getTmass(1);
			w_mass[0] = di.getWmass(0);
			w_mass[1] = di.getWmass(1);

			lep_a = di.getZlep();
			lep_b = di.getClep();
			bl_a = di.getZbl();
			bl_b = di.getCbl();

			result = calc_dilep(t_mass, w_mass, in_mpx, in_mpy, in_mpz, &lep_a, &lep_b, &bl_a, &bl_b);

			// Check if there is any solutions for this reconstruction
			if (result->size())
				++hasSolution;  // increment solution counter

			di.setHasSol(hasSolution);
			di.setResult(result);

			// time measurement
			#ifdef MEASURE_DILEP
			stopTimer(time);
			#endif
		}

		// TLorentzs/Flags
		// [0] -> x
		// [1] -> y
		// [2] -> z
		// [3] -> E/isb
		// [4] -> M

		void applyVariance (double _in_mpx[], double _in_mpy[], double _z_lepWFlags[], double _c_lepWFlags[],
			double _z_bjWFlags[], double _c_bjWFlags[], double _z_lep[], double _c_lep[], double _z_bj[], double _c_bj[],
			double _z_bl[], double _c_bl[], double _MissPx, double _MissPy) {

			unsigned tid = _tid;
			//unsigned tid = threadIdx.x + blockIdx.x * blockDim.x;;

			// Using pointers for better code readbility - does it affect the performance in the kernel?
			double *in_mpx		= &STRIDE2(_in_mpx, 0);
			double *in_mpy		= &STRIDE2(_in_mpy, 0);
			double *z_lepWFlags = &STRIDE5(_z_lepWFlags, 0);
			double *c_lepWFlags = &STRIDE5(_c_lepWFlags, 0);
			double *z_bjWFlags	= &STRIDE5(_z_bjWFlags, 0);
			double *c_bjWFlags  = &STRIDE5(_c_bjWFlags, 0);
			double *z_lep 		= &STRIDE5(_z_lep, 0);
			double *c_lep 		= &STRIDE5(_c_lep, 0);
			double *z_bj 		= &STRIDE5(_z_bj, 0);
			double *c_bj 		= &STRIDE5(_c_bj, 0);
			double *z_bl 		= &STRIDE5(_z_bl, 0);
			double *c_bl 		= &STRIDE5(_c_bl, 0);


			// new four-vectors	
			double n_Px, n_Py, n_Pz, n_Pt, n_E;	
			double delPx, delPy;

			// Vary!

			// _______________________________
			// _______z_lep___________________
			// _______________________________
			if (  abs(  z_lepWFlags[3]  )  ==  11  ){ //___electrons____
				n_Px = z_lepWFlags[0] * ( 1. + _t_rnd_.Gaus( 0., RESOLUTION ) );
				n_Py = z_lepWFlags[1] * ( 1. + _t_rnd_.Gaus( 0., RESOLUTION ) );
				n_Pz = z_lepWFlags[2] * ( 1. + _t_rnd_.Gaus( 0., RESOLUTION ) );
			} else if (  abs(z_lepWFlags[3]) == 13 ){ //_____muons______
				n_Px = z_lepWFlags[0] * ( 1. + _t_rnd_.Gaus( 0., RESOLUTION ) );
				n_Py = z_lepWFlags[1] * ( 1. + _t_rnd_.Gaus( 0., RESOLUTION ) );
				n_Pz = z_lepWFlags[2] * ( 1. + _t_rnd_.Gaus( 0., RESOLUTION ) );
			}
			// Recalculate z_lep
			n_E = sqrt ( n_Px*n_Px + n_Py*n_Py + n_Pz*n_Pz + z_lepWFlags[4]*z_lepWFlags[4] );
			z_lep[0] = n_Px;	// Change Px 				
			z_lep[1] = n_Py; 	// Change Py 	
			z_lep[2] = n_Pz; 	// Change Pz 
			z_lep[3] =  n_E; 	// Change E 
			// Propagate to MissPx and MissPy
			delPx = z_lepWFlags[0] - n_Px; 
			delPy = z_lepWFlags[1] - n_Py;			
			in_mpx[0] = _MissPx + delPx; in_mpx[1] = _MissPx + delPx; // initialize miss(Px,Py) neutrino 1
			in_mpy[0] = _MissPy + delPy; in_mpy[1] = _MissPy + delPy; // initialize miss(Px,Py) neutrino 2

			// _______________________________
			// _______c_lep___________________
			// _______________________________
			if (  abs(  c_lepWFlags[3]  )  ==  11  ){ //___electrons____
				n_Px = c_lepWFlags[0] * ( 1. + _t_rnd_.Gaus( 0., RESOLUTION ) );
				n_Py = c_lepWFlags[1] * ( 1. + _t_rnd_.Gaus( 0., RESOLUTION ) );
				n_Pz = c_lepWFlags[2] * ( 1. + _t_rnd_.Gaus( 0., RESOLUTION ) );
			} else if (  abs(c_lepWFlags[3]) == 13 ){ //_____muons______
				n_Px = c_lepWFlags[0] * ( 1. + _t_rnd_.Gaus( 0., RESOLUTION ) );
				n_Py = c_lepWFlags[1] * ( 1. + _t_rnd_.Gaus( 0., RESOLUTION ) );
				n_Pz = c_lepWFlags[2] * ( 1. + _t_rnd_.Gaus( 0., RESOLUTION ) );
			}
			// Recalculate c_lep
			n_E = sqrt ( n_Px*n_Px + n_Py*n_Py + n_Pz*n_Pz + c_lepWFlags[4]*c_lepWFlags[4] );
			c_lep[0] = n_Px;	// Change Px 				
			c_lep[1] = n_Py; 	// Change Py 	
			c_lep[2] = n_Pz; 	// Change Pz 
			c_lep[3] = n_E ; 	// Change E 
			// Propagate to MissPx and MissPy
			delPx = c_lepWFlags[0] - n_Px; 
			delPy = c_lepWFlags[1] - n_Py;			
			in_mpx[0] += delPx; in_mpx[1] += delPx; // correct miss(Px,Py) neutrino 1
			in_mpy[0] += delPy; in_mpy[1] += delPy; // correct miss(Px,Py) neutrino 2

			// _______________________________
			// _______z_bj____________________
			// _______________________________
			n_Px = z_bjWFlags[0] * ( 1. + _t_rnd_.Gaus( 0., RESOLUTION ) );
			n_Py = z_bjWFlags[1] * ( 1. + _t_rnd_.Gaus( 0., RESOLUTION ) );
			n_Pz = z_bjWFlags[2] * ( 1. + _t_rnd_.Gaus( 0., RESOLUTION ) );
			// Recalculate z_bj
			n_E = sqrt ( n_Px*n_Px + n_Py*n_Py + n_Pz*n_Pz + z_bjWFlags[4]*z_bjWFlags[4] );
			z_bj[0] = n_Px;	// Change Px 				
			z_bj[1] = n_Py; 	// Change Py 	
			z_bj[2] = n_Pz; 	// Change Pz 
			z_bj[3] = n_E ; 	// Change E 
			// Propagate to MissPx and MissPy
			delPx = z_bjWFlags[0] - n_Px; 
			delPy = z_bjWFlags[1] - n_Py;			
			in_mpx[0] += delPx; in_mpx[1] += delPx; // correct miss(Px,Py) neutrino 1
			in_mpy[0] += delPy; in_mpy[1] += delPy; // correct miss(Px,Py) neutrino 2


			// _______________________________
			// _______c_bj____________________
			// _______________________________
			n_Px = c_bjWFlags[0] * ( 1. + _t_rnd_.Gaus( 0., RESOLUTION ) );
			n_Py = c_bjWFlags[1] * ( 1. + _t_rnd_.Gaus( 0., RESOLUTION ) );
			n_Pz = c_bjWFlags[2] * ( 1. + _t_rnd_.Gaus( 0., RESOLUTION ) );
		//	n_Pt = c_bjWFlags.Pt() * ( 1. + _t_rnd_.Gaus( 0., St_j ) );
		//	n_E  = c_bjWFlags.E()  * ( 1. + _t_rnd_.Gaus( 0., Se_j ) );
			// Recalculate c_bj
			n_E = sqrt ( n_Px*n_Px + n_Py*n_Py + n_Pz*n_Pz + c_bjWFlags[4]*c_bjWFlags[4] );
			c_bj[0] = n_Px;	// Change Px 				
			c_bj[1] = n_Py; 	// Change Py 	
			c_bj[2] = n_Pz; 	// Change Pz 
			c_bj[3] = n_E ; 	// Change E 
			// Propagate to MissPx and MissPy
			delPx = c_bjWFlags[0] - n_Px; 
			delPy = c_bjWFlags[1] - n_Py;			
			in_mpx[0] += delPx; in_mpx[1] += delPx; // correct miss(Px,Py) neutrino 1
			in_mpy[0] += delPy; in_mpy[1] += delPy; // correct miss(Px,Py) neutrino 2

			// ---------------------------------------
			// Define TLorentzVectors for (b,l) system
			// ---------------------------------------
			//z_bl = z_bj + z_lep;
			//c_bl = c_bj + c_lep;

			z_bl[0] = z_bj[0] + z_lep[0];
			z_bl[1] = z_bj[1] + z_lep[1];
			z_bl[2] = z_bj[2] + z_lep[2];
			z_bl[3] = z_bj[3] + z_lep[3];

			c_bl[0] = c_bj[0] + c_lep[0];
			c_bl[1] = c_bj[1] + c_lep[1];
			c_bl[2] = c_bj[2] + c_lep[2];
			c_bl[3] = c_bj[3] + c_lep[3];

			// Re-calculate the masses
			calcMass(z_bl);
			calcMass(c_bl);
		}

		//__global__
		void dilep_kernel (double _in_mpx[], double _in_mpy[], double _z_lepWFlags[], double _c_lepWFlags[],
			double _z_bjWFlags[], double _c_bjWFlags[], double _z_lep[], double _c_lep[], double _z_bj[], double _c_bj[],
			double _MissPx, double _MissPy, double _t_mass[], double _w_mass[], double nc[], int a[]) {

			// CPU version
			double _z_bl[5 * size], _c_bl[5 * size];

			for (_tid = 0; _tid < size; ++_tid)
				applyVariance(_in_mpx, _in_mpy, _z_lepWFlags, _c_lepWFlags, _z_bjWFlags, _c_bjWFlags,
					_z_lep, _c_lep, _z_bj, _c_bj, _z_bl, _c_bl, _MissPx, _MissPy);
			
			for (unsigned tid = 0; tid < size; ++tid)
				calc_dilep(_t_mass, _w_mass, _in_mpx, _in_mpy, 
							_z_lep, _c_lep, _z_bl, _c_bl, nc, a, tid);

			// GPU version
			//double _z_bl[5], _c_bl[5];

			//applyVariance(_in_mpx, _in_mpy, _z_lepWFlags, _c_lepWFlags, _z_bjWFlags, _c_bjWFlags,
			//		_z_lep, _c_lep, _z_bj, _c_bj, _z_bl, _c_bl, _MissPx, _MissPy);

			//calc_dilep(_t_mass, _w_mass, _in_mpx, _in_mpy, 
			//				_z_lep, _c_lep, _z_bl, _c_bl, nc, a);
		}

		void dilep (vector<DilepInput> &di) {

			size = di.size();
			
			double in_mpx[2 * size], in_mpy[2 * size], in_mpz[2 * size], 
				   t_mass[2 * size], w_mass[2 * size];
			
			double *dev_t_mass, *dev_w_mass, *dev_in_mpx, *dev_in_mpy;
			double a[5 * size], b[5 * size], c[5 * size], d[5 * size], e[5 * size], f[5 * size];	// e and f are the z/c_bl
			double aFlags[5 * size], bFlags[5 * size], cFlags[5 * size], dFlags[5 * size];
		
			double *dev_lep_a, *dev_lep_b, *dev_bj_a, *dev_bj_b;
			double *dev_lep_aFlags, *dev_lep_bFlags, *dev_bj_aFlags, *dev_bj_bFlags;
			double nc[16*size];
			double *dev_nc, *dev_MissPx, *dev_MissPy;
			int count[size], *dev_count;
			int hasSolution = 0;

			double _misspx = di[0].getMissPx(), _misspy = di[0].getMissPy();

			// time measurement
			#ifdef MEASURE_DILEP
			long long int time = startTimer();
			#endif

			for (unsigned i = 0; i < size; ++i) {

				in_mpx[i * 2]		= di[i].getInMpx(0);
				in_mpx[(i * 2) + 1] = di[i].getInMpx(1);
				in_mpy[i * 2]		= di[i].getInMpy(0);
				in_mpy[(i * 2) + 1] = di[i].getInMpy(1);
				t_mass[i * 2]		= di[i].getTmass(0);
				t_mass[(i * 2) + 1] = di[i].getTmass(1);
				w_mass[i * 2]		= di[i].getWmass(0);
				w_mass[(i * 2) + 1] = di[i].getWmass(1);
				
				// z_lep
				a[i * 5]	   = di[i].getZlep().Px();
				a[(i * 5) + 1] = di[i].getZlep().Py();
				a[(i * 5) + 2] = di[i].getZlep().Pz();
				a[(i * 5) + 3] = di[i].getZlep().E();
				a[(i * 5) + 4] = di[i].getZlep().M();

				// z_lepWFlags
				aFlags[i * 5]	    = di[i].getZlepW().Px();
				aFlags[(i * 5) + 1] = di[i].getZlepW().Py();
				aFlags[(i * 5) + 2] = di[i].getZlepW().Pz();
				aFlags[(i * 5) + 3] = di[i].getZlepW().isb;
				aFlags[(i * 5) + 4] = di[i].getZlepW().M();

				// c_lep
				b[i * 5]	   = di[i].getClep().Px();
				b[(i * 5) + 1] = di[i].getClep().Py();
				b[(i * 5) + 2] = di[i].getClep().Pz();
				b[(i * 5) + 3] = di[i].getClep().E();
				b[(i * 5) + 4] = di[i].getClep().M();

				// c_lepWFlags
				bFlags[i * 5]	    = di[i].getClepW().Px();
				bFlags[(i * 5) + 1] = di[i].getClepW().Py();
				bFlags[(i * 5) + 2] = di[i].getClepW().Pz();
				bFlags[(i * 5) + 3] = di[i].getClepW().isb;
				bFlags[(i * 5) + 4] = di[i].getClepW().M();

				// z_bj
				c[i * 5]	   = di[i].getZbj().Px();
				c[(i * 5) + 1] = di[i].getZbj().Py();
				c[(i * 5) + 2] = di[i].getZbj().Pz();
				c[(i * 5) + 3] = di[i].getZbj().E();
				c[(i * 5) + 4] = di[i].getZbj().M();

				// z_bjWFlags
				cFlags[i * 5]	    = di[i].getZbjW().Px();
				cFlags[(i * 5) + 1] = di[i].getZbjW().Py();
				cFlags[(i * 5) + 2] = di[i].getZbjW().Pz();
				cFlags[(i * 5) + 3] = di[i].getZbjW().isb;
				cFlags[(i * 5) + 4] = di[i].getZbjW().M();

				// c_bj
				d[i * 5]	   = di[i].getCbj().Px();
				d[(i * 5) + 1] = di[i].getCbj().Py();
				d[(i * 5) + 2] = di[i].getCbj().Pz();
				d[(i * 5) + 3] = di[i].getCbj().E();
				d[(i * 5) + 4] = di[i].getCbj().M();

				// c_bjWFlags
				dFlags[i * 5]	    = di[i].getCbjW().Px();
				dFlags[(i * 5) + 1] = di[i].getCbjW().Py();
				dFlags[(i * 5) + 2] = di[i].getCbjW().Pz();
				dFlags[(i * 5) + 3] = di[i].getCbjW().isb;
				dFlags[(i * 5) + 4] = di[i].getCbjW().M();
			}

			// GPU memory allocation of the inputs and outputs of the dilep kernel
			//cudaMalloc(&dev_t_mass, size*2*sizeof(double));
			//cudaMalloc(&dev_w_mass, size*2*sizeof(double));
			//cudaMalloc(&dev_in_mpx, size*2*sizeof(double));
			//cudaMalloc(&dev_in_mpy, size*2*sizeof(double));

			//cudaMalloc(&dev_lep_a, sizeof(a));
			//cudaMalloc(&dev_lep_b, sizeof(b));
			//cudaMalloc(&dev_bj_a, sizeof(c));
			//cudaMalloc(&dev_bj_b, sizeof(d));

			//cudaMalloc(&dev_lep_aFlags, sizeof(aFlags));
			//cudaMalloc(&dev_lep_bFlags, sizeof(bFlags));
			//cudaMalloc(&dev_bj_aFlags, sizeof(cFlags));
			//cudaMalloc(&dev_bj_bFlags, sizeof(dFlags));

			//cudaMalloc(&dev_MissPx, sizeof(double));
			//cudaMalloc(&dev_MissPy, sizeof(double));

			// allocation of the results
			//cudaMalloc(&dev_nc, size*16*sizeof(double));
			//cudaMalloc(&dev_count, size*sizeof(int));

			// transfer the inputs to GPU memory
			//cudaMemcpy(dev_t_mass, t_mass, sizeof(t_mass), cudaMemcpyHostToDevice);
			//cudaMemcpy(dev_w_mass, w_mass, sizeof(w_mass), cudaMemcpyHostToDevice);
			//cudaMemcpy(dev_in_mpx, in_mpx, sizeof(in_mpx), cudaMemcpyHostToDevice);
			//cudaMemcpy(dev_in_mpy, in_mpy, sizeof(in_mpy), cudaMemcpyHostToDevice);

			//cudaMemcpy(dev_lep_a, a, sizeof(a), cudaMemcpyHostToDevice);
			//cudaMemcpy(dev_lep_b, b, sizeof(b), cudaMemcpyHostToDevice);
			//cudaMemcpy(dev_bj_a, c, sizeof(c), cudaMemcpyHostToDevice);
			//cudaMemcpy(dev_bj_b, d, sizeof(d), cudaMemcpyHostToDevice);

			//cudaMemcpy(dev_lep_aFlags, aFlags, sizeof(aFlags), cudaMemcpyHostToDevice);
			//cudaMemcpy(dev_lep_bFlags, bFlags, sizeof(bFlags), cudaMemcpyHostToDevice);
			//cudaMemcpy(dev_bj_aFlags, cFlags, sizeof(cFlags), cudaMemcpyHostToDevice);
			//cudaMemcpy(dev_bj_bFlags, dFlags, sizeof(dFlags), cudaMemcpyHostToDevice);

			//cudaMemcpy(dev_MissPx, _misspx, sizeof(double), cudaMemcpyHostToDevice);
			//cudaMemcpy(dev_MissPy, _misspy, sizeof(double), cudaMemcpyHostToDevice);
			

			dilep_kernel(in_mpx, in_mpy, aFlags, bFlags, cFlags, dFlags,
					a, b, c, d, _misspx, _misspy, t_mass, w_mass, nc, count);


			//dilep_kernel <<< 1, size >>> (dev_in_mpx, dev_in_mpy, dev_lep_aFlags, dev_lep_bFlags, dev_bj_aFlags, dev_bj_bFlags,
			//		dev_lep_a, dev_lep_b, dev_bj_a, dev_bj_b, dev_MissPx, dev_MissPy, dev_t_mass, dev_w_mass, dev_nc, dev_count);
			
			
			// memory transfer of the results from the GPU
			//cudaMemcpy(nc, dev_nc, 16*size*sizeof(double), cudaMemcpyDeviceToHost);
			//cudaMemcpy(count, dev_count, size*sizeof(int), cudaMemcpyDeviceToHost);

			// reconstruction of the normal output of dilep
			// o num de combs*vars e o num de threads

			for (unsigned comb = 0; comb < size; ++comb) {
				vector<myvector> result;

				for (int sol = 0 ; sol < count[comb] && sol<4 ; sol++) {
					myvector *mv = new myvector( 
						TO1D(nc,comb,sol,0),
						TO1D(nc,comb,sol,1),
						TO1D(nc,comb,sol,2),
						TO1D(nc,comb,sol,3) );
					
					result.push_back(*mv);
				}
				if(result.size())
					++hasSolution;

				di[comb].setHasSol(hasSolution);
				di[comb].setResult(&result);
			}

			// time measurement
			#ifdef MEASURE_DILEP
			stopTimer(time);
			#endif
			
		}

		//__global__
		void calc_dilep(double t_mass[], double w_mass[], 
				double in_mpx[], double in_mpy[], double _lep_a[], 
				double _lep_b[], double _bl_a[], double _bl_b[], 
				double nc[], int a[])
		{

			//unsigned tid = threadIdx.x + blockIdx.x * blockDim.x;
			unsigned tid = 1;
			double G_1, G_3;
			double WMass_a, WMass_b, tMass_a, tMass_b, lep_a[5], lep_b[5], bl_a[5], bl_b[5];
			double in_mpz[2] = {0.0, 0.0};


			WMass_a = STRIDE2(w_mass, 0);
			tMass_a = STRIDE2(t_mass, 0);
			WMass_b = STRIDE2(w_mass, 1);
			tMass_b = STRIDE2(t_mass, 1);

			for (unsigned i = 0; i < 5; ++i) {
				lep_a[i] = STRIDE5(_lep_a, i);
				lep_b[i] = STRIDE5(_lep_b, i);

				bl_a[i] = STRIDE5(_bl_a, i);
				bl_b[i] = STRIDE5(_bl_b, i);
			}
			

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
			toz(in_a, lep_a, out_a);

			double in_c[5],out_c[6];
			in_c[0] = G_12/G_11;
			in_c[1] = -1*G_9/G_11;
			in_c[2] = -1*G_10/G_11;
			in_c[3] = lep_b[3];
			in_c[4] = G_3;
			toz(in_c, lep_b, out_c);
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
			my_qu(re, output);

			int ncand(0);

			double rec_x1, rec_y1, rec_z1, rec_e1, rec_x2, rec_y2, rec_z2, rec_e2;

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
					
					// aqui podem nao chegar as threads todas
					//__syncthreads();
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

		//__host__
		void calc_dilep(double t_mass[], double w_mass[], 
				double in_mpx[], double in_mpy[], double _lep_a[], 
				double _lep_b[], double _bl_a[], double _bl_b[], 
				double nc[], int a[], unsigned _tid)
		{

			unsigned tid = _tid;
			double G_1, G_3;
			double WMass_a, WMass_b, tMass_a, tMass_b, lep_a[5], lep_b[5], bl_a[5], bl_b[5];
			double in_mpz[2] = {0.0, 0.0};


			WMass_a = STRIDE2(w_mass, 0);
			tMass_a = STRIDE2(t_mass, 0);
			WMass_b = STRIDE2(w_mass, 1);
			tMass_b = STRIDE2(t_mass, 1);

			for (unsigned i = 0; i < 5; ++i) {
				lep_a[i] = STRIDE5(_lep_a, i);
				lep_b[i] = STRIDE5(_lep_b, i);

				bl_a[i] = STRIDE5(_bl_a, i);
				bl_b[i] = STRIDE5(_bl_b, i);
			}

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
			toz(in_a, lep_a, out_a);

			double in_c[5],out_c[6];
			in_c[0] = G_12/G_11;
			in_c[1] = -1*G_9/G_11;
			in_c[2] = -1*G_10/G_11;
			in_c[3] = lep_b[3];
			in_c[4] = G_3;
			toz(in_c, lep_b, out_c);
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
			my_qu(re, output);

			int ncand(0);

			double rec_x1, rec_y1, rec_z1, rec_e1, rec_x2, rec_y2, rec_z2, rec_e2;

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
					
					// aqui podem nao chegar as threads todas
					//__syncthreads();
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

		// NEUTRINO SOLUTIONS
		vector<myvector>* calc_dilep(double t_mass[], double w_mass[], 
				double in_mpx[], double in_mpy[], double in_mpz[],
				TLorentzVector* lep_a, TLorentzVector* lep_b, 
				TLorentzVector* bl_a, TLorentzVector* bl_b)
		{

			double mpx, mpy, G_1, G_2, G_3, G_4, _d, _a, _f, _b, _e, _g;

			/////////// input part

			mpx = in_mpx[0];
			mpy = in_mpy[0];

			//  printf( " mpx, mpy are: %3.3f  %3.3f\n", mpx, mpy);
			const double WMass_a = w_mass[0];
			const double tMass_a = t_mass[0];
			const double WMass_b = w_mass[1];
			const double tMass_b = t_mass[1];  
			/////   
			G_1 = WMass_a*WMass_a - ( lep_a->M() )*( lep_a->M() );
			G_3 = WMass_b*WMass_b - ( lep_b->M() )*( lep_b->M() );
			G_2 = tMass_a*tMass_a - ( bl_a->M()  )*( bl_a->M()  ); 
			G_4 = tMass_b*tMass_b - ( bl_b->M()  )*( bl_b->M()  );

			double S=mpx;
			double T=mpy;

			///////////////////////////////////////////////////////////////////

			double G_5 =	( bl_a->Px()/bl_a->E() - lep_a->Px()/lep_a->E() );
			double G_6 =	( bl_a->Py()/bl_a->E() - lep_a->Py()/lep_a->E() );
			double G_7 =	( bl_a->Pz()/bl_a->E() - lep_a->Pz()/lep_a->E() );
			double G_8 =	( G_1/lep_a->E() - G_2/bl_a->E() )/2.;

			double G_9 =	( bl_b->Px()/bl_b->E() - lep_b->Px()/lep_b->E() );
			double G_10 =	( bl_b->Py()/bl_b->E() - lep_b->Py()/lep_b->E() );
			double G_11 =	( bl_b->Pz()/bl_b->E() - lep_b->Pz()/lep_b->E() );
			double G_12 =	( G_3/lep_b->E() - G_4/bl_b->E() )/2.;

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
			in_a[3] = lep_a->E();
			in_a[4] = G_1;

			on_a[0] = lep_a->Px();
			on_a[1] = lep_a->Py();
			on_a[2] = lep_a->Pz();
			toz(in_a, on_a, out_a);

			///////// 2nd ////////////////////
			double in_c[5],on_c[3],out_c[6];
			in_c[0] = G_12/G_11;
			in_c[1] = -1*G_9/G_11;
			in_c[2] = -1*G_10/G_11;
			in_c[3] = lep_b->E();
			in_c[4] = G_3;

			on_c[0] = lep_b->Px();
			on_c[1] = lep_b->Py();
			on_c[2] = lep_b->Pz();
			toz(in_c, on_c, out_c);

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

			if ( out_a[0]!=0  ){
				_a = out_a[0];
				_b = out_a[1];
				_d = out_a[2];
				_e = out_a[3];
				_f = out_a[4];
				_g = out_a[5];
			} 
			if (out_a[0]==0 && _A!=0 ) {
				_a = 1.;
				_b = out_a[1] + _B/_A;
				_d = out_a[2] + _D/_A;
				_e = out_a[3] + _E/_A;
				_f = out_a[4] + _F/_A;
				_g = out_a[5] + _G/_A;
			}
			if ( out_a[0]==0 && _A==0){
				return NULL;
			}

			double fx_1 = 2.*(_A*_d - _a*_D);
			double fx_2 = 2.*(_A*_f - _a*_F);
			double fx_3 = _a*_B - _A*_b;
			double fx_4 = 2.*(_E*_a - _A*_e);
			double fx_5 = _a*_G - _g*_A;

			double k_1 = ( _f*_f - _a*_b )/_a/_a;
			double k_2 = ( 2.*_d*_f - 2.*_a*_e )/_a/_a;
			double k_3 = ( _d*_d-_a*_g )/_a/_a;
			double k_4 = -1.*_d/_a;
			double k_5 = -1.*_f/_a;

			/////
			///// the part above is 
			///// x = (fx3*y**2 + fx4*y + fx5)/(fx1 + fx2*y)
			///// used to get x value once y is known
			/////
			//// if fx1 + fx2*y == 0, then x is 
			//// x = +/-sqrt(k1*y**2 + k2*y + k3) + (k4 + k5*y)
			////


			double g_1 = 4.*_A*_A*k_5*k_5 + 4.*_F*_F + 8.*_A*_F*k_5;
			double g_2 = 8.*_A*_A*k_4*k_5 + 8.*_A*_D*k_5 + 8.*_A*_F*k_4 + 8.*_D*_F;
			double g_3 = 4.*_A*_A*k_4*k_4 + 4.*_D*_D + 8.*_A*_D*k_4;
			double g_4 = _A*k_1 + _A*k_5*k_5;
			double g_5 = _A*k_2 + 2.*_A*k_4*k_5 + 2.*_D*k_5;
			double g_6 = _A*k_3 + _A*k_4*k_4 + 2.*_D*k_4 + _G;

			double m_1 = g_1*k_1;
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
			my_qu(re,output);

			double delta;
			int ncand(0);

			vector<myvector> *neutrinoContainer = new std::vector<myvector>; //// must use new

			//////////////////////////////////////////////////
			//// please replace your code starting from here 
			//// ////////////////////////////////////////////


			double rec_x1;
			double rec_y1;
			double rec_z1;
			double rec_e1;
			double rec_x2;
			double rec_y2;
			double rec_z2;
			double rec_e2;

			double m_delta_mass;

			for (int j=0; j<8; j+=2){
				delta = k_1*output[j]*output[j] + k_2*output[j] + k_3;
				if ( output[j+1]==0 && delta>=0) {
					if ( (fx_1 + fx_2*output[j])!=0 ) {
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
					double m_w1 = calcMass(rec_x1+lep_a->Px(), rec_y1+lep_a->Py(), rec_z1+lep_a->Pz(), rec_e1+lep_a->E());
					double m_w2 = calcMass(rec_x2+lep_b->Px(), rec_y2+lep_b->Py(), rec_z2+lep_b->Pz(), rec_e2+lep_b->E());
					double m_t1 = calcMass(rec_x1+ bl_a->Px(), rec_y1+ bl_a->Py(), rec_z1+ bl_a->Pz(), rec_e1+ bl_a->E());
					double m_t2 = calcMass(rec_x2+ bl_b->Px(), rec_y2+ bl_b->Py(), rec_z2+ bl_b->Pz(), rec_e2+ bl_b->E());

					m_delta_mass = 1000.0; // allow mass variation range for reco W and tops..
					bool m_good_eq1 = ( fabs(S -(rec_x1+rec_x2)) 	  <=0.01 )?true:false;
					bool m_good_eq2 = ( fabs(T -(rec_y1+rec_y2)) 	  <=0.01 )?true:false;
					bool m_good_eq3 = ( fabs(m_w1 - w_mass[0]) <= m_delta_mass )?true:false;
					bool m_good_eq4 = ( fabs(m_w2 - w_mass[1]) <= m_delta_mass )?true:false;
					bool m_good_eq5 = ( fabs(m_t1 - t_mass[0]) <= m_delta_mass )?true:false;
					bool m_good_eq6 = ( fabs(m_t2 - t_mass[1]) <= m_delta_mass )?true:false;

					if ( m_good_eq1 && m_good_eq2 && m_good_eq3 && m_good_eq4 && m_good_eq5 && m_good_eq6 ) {
						neutrinoContainer->push_back( myvector(rec_x1, rec_y1, rec_z1, rec_z2) );
						ncand++;
					} /*else {
						cout<<" One solution is found to be out of range: "<< m_w1<<" "<<m_w2<<" "<< m_t1 <<" "<< m_t2<<endl;
						}*/  	  
				}
			}

			return neutrinoContainer;
		}

		//////////////////////////////////////
		void toz(double k[], double l[], double g[]){
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
		void my_qu( double my_in[], double my_val[])
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

			if ( a ==0 && b!=0 ){
				double input[4]={b,c,d,e};
				cubic(input,real,img);
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

			if ( a ==0 && b==0 && c!=0){
				double alpha = pow(d,2)-4*c*e;
				if (alpha>=0) {
					x1_r = (-1*d + sqrt(alpha))/2/c;	x1_i = 0;
					x2_r = (-1*d - sqrt(alpha))/2/c;	x2_i = 0;
				} else {
					x1_r = (-1*d )/2/c;	x1_i = sqrt(-alpha)/2/c;
					x2_r = (-1*d )/2/c;	x2_i = -1*sqrt(-alpha)/2/c;  	
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

			if (a ==0 && b==0 && c==0 && d!=0){
				x1_r = -e/d;	x1_i = 0;
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

			double input[4]={c_1,c_2,c_3,c_4};
			cubic(input,real,img);

			////////////////////////////////////////////////
			////// (5) sqrt root of the cubic equation solutions
			////////////////////////////////////////////////

			double out_r[3]={0,0,0};
			double out_i[3]={0,0,0};
			bool _img[3]={false,false,false};
			int img_index =0;
			double my[2];

			for (int ii=0; ii<3; ii++){
				Csqrt(real[ii],img[ii],my);
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
			if( r_i !=0 )	{ 
				r_r 	= -1.*_g/8/r_i;
			} else {
				r_i 	= 0;
			}	

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
		}
		////////////////////end of main
		///////////////////////////////////////////////////////////////
		////+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		void Csqrt(double _ar, double _ai, double _my[])
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
		void cubic(double a[], double rr[], double ri[])
		{
			int i;
			double a0, a1, a2, a3;
			double g, h, y1, sh, theta, xy1, xy2, xy3;
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
				theta = acos(g / (2.0 * h * sh)) / 3.0;
				xy1 = 2.0 * sh * cos(theta);
				xy2 = 2.0 * sh * cos(theta + (2.0 * TPI / 3.0));
				xy3 = 2.0 * sh * cos(theta + (4.0 * TPI / 3.0));
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
	}
}

