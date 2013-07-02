#include "neut.h"

using std::vector;
using namespace std;

namespace Dilep {
	namespace Phi {
		

		// Wrapper for the dilep calculation using a vector of the input class
		// vdi vector with DilepInput varied for a jet combo
		void dilep (std::vector<DilepInput> &vdi) {
			
			// time measurement
			#ifdef MEASURE_DILEP
			long long int time = startTimer();
			#endif
			unsigned size = vdi.size();
			int hasSolution = 0;
			unsigned size2  = 2  * size;
			unsigned size5  = 5  * size;	
			unsigned size16 = 16 * size;


			double in_mpx[size2], in_mpy[size2], t_mass[size2], w_mass[size2];
			double lep_a[size5], lep_b[size5], bl_a[size5], bl_b[size5];
			
			double __attribute__((target(mic))) nc[size16];
			int __attribute__((target(mic))) count[size];

			for (unsigned i = 0; i < size; ++i) {

				in_mpx[(i * 2) + 0] = vdi[i].getInMpx(0);
				in_mpx[(i * 2) + 1] = vdi[i].getInMpx(1);
				in_mpy[(i * 2) + 0] = vdi[i].getInMpy(0);
				in_mpy[(i * 2) + 1] = vdi[i].getInMpy(1);
				t_mass[(i * 2) + 0] = vdi[i].getTmass(0);
				t_mass[(i * 2) + 1] = vdi[i].getTmass(1);
				w_mass[(i * 2) + 0] = vdi[i].getWmass(0);
				w_mass[(i * 2) + 1] = vdi[i].getWmass(1);

				lep_a[(i * 5) + 0] = vdi[i].getZlep().Px();
				lep_a[(i * 5) + 1] = vdi[i].getZlep().Py();
				lep_a[(i * 5) + 2] = vdi[i].getZlep().Pz();
				lep_a[(i * 5) + 3] = vdi[i].getZlep().E(); 
				lep_a[(i * 5) + 4] = vdi[i].getZlep().M();

				lep_b[(i * 5) + 0] = vdi[i].getClep().Px();
				lep_b[(i * 5) + 1] = vdi[i].getClep().Py();
				lep_b[(i * 5) + 2] = vdi[i].getClep().Pz();
				lep_b[(i * 5) + 3] = vdi[i].getClep().E(); 
				lep_b[(i * 5) + 4] = vdi[i].getClep().M();

				bl_a[(i * 5) + 0] = vdi[i].getZbl().Px();
				bl_a[(i * 5) + 1] = vdi[i].getZbl().Py();
				bl_a[(i * 5) + 2] = vdi[i].getZbl().Pz();
				bl_a[(i * 5) + 3] = vdi[i].getZbl().E(); 
				bl_a[(i * 5) + 4] = vdi[i].getZbl().M();

				bl_b[(i * 5) + 0] = vdi[i].getCbl().Px();
				bl_b[(i * 5) + 1] = vdi[i].getCbl().Py();
				bl_b[(i * 5) + 2] = vdi[i].getCbl().Pz();
				bl_b[(i * 5) + 3] = vdi[i].getCbl().E(); 
				bl_b[(i * 5) + 4] = vdi[i].getCbl().M();
			}
					

			
			calc_dilep_mic(t_mass, w_mass, in_mpx, in_mpy, lep_a, 
						lep_b, bl_a, bl_b, nc, count, size);
					

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

				vdi[comb].setHasSol(hasSolution);
				vdi[comb].setResult(&result);
			}		

			// time measurement
			#ifdef MEASURE_DILEP
			stopTimer(time);
			#endif
			
		}
}
}