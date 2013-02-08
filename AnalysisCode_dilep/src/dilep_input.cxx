#include "dilep_input.h"

// Constructor
DilepInput::DilepInput (TLorentzVector _z_lep, TLorentzVector _c_lep, TLorentzVector _z_bj, TLorentzVector _c_bj,
						TLorentzVectorWFlags _z_bjWFlags, TLorentzVectorWFlags _c_bjWFlags, TLorentzVectorWFlags _z_lepWFlags,
						TLorentzVectorWFlags _c_lepWFlags, TLorentzVectorWFlags _jet1_HiggsWFlags, TLorentzVectorWFlags _jet2_HiggsWFlags,
						double _in_mpx[], double _in_mpy[], double _in_mpz[], double _MissPx,
						double _MissPy, double _t_mass[], double _w_mass[]) {

	z_lep = _z_lep;
	c_lep = _c_lep;
	z_bj = _z_bj;
	c_bj = _c_bj;

	z_bjWFlags = _z_bjWFlags;
	c_bjWFlags = _c_bjWFlags;
	z_lepWFlags = _z_lepWFlags;
	c_lepWFlags = _c_lepWFlags;
	jet1_HiggsWFlags = _jet1_HiggsWFlags;
	jet2_HiggsWFlags = _jet2_HiggsWFlags;

	in_mpx[0] = _in_mpx[0];
	in_mpx[1] = _in_mpx[1];
	in_mpy[0] = _in_mpy[0];
	in_mpy[1] = _in_mpy[1];
	in_mpz[0] = _in_mpz[0];
	in_mpz[1] = _in_mpz[1];

	MissPx = _MissPx;
	MissPy = _MissPy;

	t_mass[0] = _t_mass[0];
	t_mass[1] = _t_mass[1];
	w_mass[0] = _w_mass[0];
	w_mass[1] = _w_mass[1];

	hasSolution = 0;
}

// Constructor
DilepInput::DilepInput (const DilepInput &other) {
	in_mpx[0] = other.getInMpx(0);
	in_mpx[1] = other.getInMpx(1);
	in_mpy[0] = other.getInMpy(0);
	in_mpy[1] = other.getInMpy(1);
	in_mpz[0] = other.getInMpz(0);
	in_mpz[1] = other.getInMpz(1);

	t_mass[0] = other.getTmass(0);
	t_mass[1] = other.getTmass(1);
	w_mass[0] = other.getWmass(0);
	w_mass[1] = other.getWmass(1);

	MissPx = other.getMissPx();
	MissPy = other.getMissPy();

	z_lep = other.getZlep();
	c_lep = other.getClep();
	z_bj = other.getZbj();
	c_bj = other.getCbj();

	z_bjWFlags = other.getZbjW();
	c_bjWFlags = other.getCbjW();
	z_lepWFlags = other.getZlepW();
	c_lepWFlags = other.getClepW();

	jet1_HiggsWFlags = other.getJet1HiggsW();
	jet2_HiggsWFlags = other.getJet2HiggsW();

	z_bl = other.getZbl();
	c_bl = other.getCbl();

	hasSolution = other.getHasSol();

	result = other.getResult();
}

// Print the private variables of the class for debug purposes
void DilepInput::print (char* filename) {
	ofstream file (filename, fstream::app);

	file << "in_mpx " << in_mpx[0] << " " << in_mpx[1] << endl;
	file << "in_mpy " << in_mpy[0] << " " << in_mpy[1] << endl;
	file << "in_mpz " << in_mpz[0] << " " << in_mpz[1] << endl;

	file << "MissPx " << MissPx << endl;
	file << "MissPy " << MissPy << endl;

	file << "t_mass " << t_mass[0] << " " << t_mass[1] << endl;
	file << "w_mass " << w_mass[0] << " " << w_mass[1] << endl;

	file << "z_lep " << z_lep.Px() << " " << z_lep.Py() << " " << z_lep.Pz() << " " << z_lep.E() << endl;
	file << "c_lep " << c_lep.Px() << " " << c_lep.Py() << " " << c_lep.Pz() << " " << c_lep.E() << endl;
	file << "z_bj " << z_bj.Px() << " " << z_bj.Py() << " " << z_bj.Pz() << " " << z_bj.E() << endl;
	file << "c_bj " << c_bj.Px() << " " << c_bj.Py() << " " << c_bj.Pz() << " " << c_bj.E() << endl;

	file << "z_bjWFlags " << z_bjWFlags.Px() << " " << z_bjWFlags.Py() << " " << z_bjWFlags.Pz() << " " << z_bjWFlags.E() << endl;
	file << "c_bjWFlags " << c_bjWFlags.Px() << " " << c_bjWFlags.Py() << " " << c_bjWFlags.Pz() << " " << c_bjWFlags.E() << endl;
	file << "z_lepWFlags " << z_lepWFlags.Px() << " " << z_lepWFlags.Py() << " " << z_lepWFlags.Pz() << " " << z_lepWFlags.E() << endl;
	file << "c_lepWFlags " << c_lepWFlags.Px() << " " << c_lepWFlags.Py() << " " << c_lepWFlags.Pz() << " " << c_lepWFlags.E() << endl;

	file << "jet1_HiggsWFlags " << jet1_HiggsWFlags.Px() << " " << jet1_HiggsWFlags.Py() << " " << jet1_HiggsWFlags.Pz() << " " << jet1_HiggsWFlags.E() << endl;
	file << "jet2_HiggsWFlags " << jet2_HiggsWFlags.Px() << " " << jet2_HiggsWFlags.Py() << " " << jet2_HiggsWFlags.Pz() << " " << jet2_HiggsWFlags.E() << endl;

	file.close();
}

// Apply variance to the inputs of the dilep function a given number of times
// resulting a vector of varied inputs
vector<DilepInput> applyVariance (vector<DilepInput> &vdi, float res, int amount, int seed) {
	vector<DilepInput> inputs;

	for (int i = 0; i < vdi.size(); ++i) {
		DilepInput di = vdi[i];

		for (int j = 0; j < amount; ++j) {
			DilepInput aux (di);
			aux.applyVariance (res, seed + i * j);

			inputs.push_back (aux);		}
	}

	return inputs;
}

// Apply variance to the inputs of the dilep function, given a resolution and a seed,
// following a Gaussian distribution
void DilepInput::applyVariance (float res, int seed) {
	// Resolution values
	double Sx_e = res, Sy_e = res, Sz_e = res, St_e = res, Se_e = res;  // electrons
	double Sx_m = res, Sy_m = res, Sz_m = res, St_m = res, Se_m = res;  // muons
	double Sx_j = res, Sy_j = res, Sz_j = res, St_j = res, Se_j = res;  // jets

	// new four-vectors	
	double n_Px, n_Py, n_Pz, n_Pt, n_E;	
	double delPx, delPy;

	// Initialize the random number generator
	TRandom3 random (seed);

	// Vary!

	// _______________________________
	// _______z_lep___________________
	// _______________________________
	if (  abs(  z_lepWFlags.isb  )  ==  11  ){ //___electrons____
		n_Px = z_lepWFlags.Px() * ( 1. + random.Gaus( 0., Sx_e ) );
		n_Py = z_lepWFlags.Py() * ( 1. + random.Gaus( 0., Sy_e ) );
		n_Pz = z_lepWFlags.Pz() * ( 1. + random.Gaus( 0., Sz_e ) );
		n_Pt = z_lepWFlags.Pt() * ( 1. + random.Gaus( 0., St_e ) );
		n_E  = z_lepWFlags.E()  * ( 1. + random.Gaus( 0., Se_e ) );
	} else if (  abs(z_lepWFlags.isb) == 13 ){ //_____muons______
		n_Px = z_lepWFlags.Px() * ( 1. + random.Gaus( 0., Sx_m ) );
		n_Py = z_lepWFlags.Py() * ( 1. + random.Gaus( 0., Sy_m ) );
		n_Pz = z_lepWFlags.Pz() * ( 1. + random.Gaus( 0., Sz_m ) );
		n_Pt = z_lepWFlags.Pt() * ( 1. + random.Gaus( 0., St_m ) );
		n_E  = z_lepWFlags.E()  * ( 1. + random.Gaus( 0., Se_m ) );
	}
	// Recalculate z_lep
	n_E = sqrt ( n_Px*n_Px + n_Py*n_Py + n_Pz*n_Pz + z_lepWFlags.M()*z_lepWFlags.M() );
	z_lep.SetPx( n_Px );	// Change Px 				
	z_lep.SetPy( n_Py ); 	// Change Py 	
	z_lep.SetPz( n_Pz ); 	// Change Pz 
	z_lep.SetE(  n_E  ); 	// Change E 
	// Propagate to MissPx and MissPy
	delPx = z_lepWFlags.Px() - n_Px; 
	delPy = z_lepWFlags.Py() - n_Py;			
	in_mpx[0] = MissPx + delPx; in_mpx[1] = MissPx + delPx; // initialize miss(Px,Py) neutrino 1
	in_mpy[0] = MissPy + delPy; in_mpy[1] = MissPy + delPy; // initialize miss(Px,Py) neutrino 2
	in_mpz[0] = 0.            ; in_mpz[1] = 0.;		// initialize neutrinos Pz to zero

	// _______________________________
	// _______c_lep___________________
	// _______________________________
	if (  abs(  c_lepWFlags.isb  )  ==  11  ){ //___electrons____
		n_Px = c_lepWFlags.Px() * ( 1. + random.Gaus( 0., Sx_e ) );
		n_Py = c_lepWFlags.Py() * ( 1. + random.Gaus( 0., Sy_e ) );
		n_Pz = c_lepWFlags.Pz() * ( 1. + random.Gaus( 0., Sz_e ) );
		n_Pt = c_lepWFlags.Pt() * ( 1. + random.Gaus( 0., St_e ) );
		n_E  = c_lepWFlags.E()  * ( 1. + random.Gaus( 0., Se_e ) );
	} else if (  abs(c_lepWFlags.isb) == 13 ){ //_____muons______
		n_Px = c_lepWFlags.Px() * ( 1. + random.Gaus( 0., Sx_m ) );
		n_Py = c_lepWFlags.Py() * ( 1. + random.Gaus( 0., Sy_m ) );
		n_Pz = c_lepWFlags.Pz() * ( 1. + random.Gaus( 0., Sz_m ) );
		n_Pt = c_lepWFlags.Pt() * ( 1. + random.Gaus( 0., St_m ) );
		n_E  = c_lepWFlags.E()  * ( 1. + random.Gaus( 0., Se_m ) );
	}
	// Recalculate c_lep
	n_E = sqrt ( n_Px*n_Px + n_Py*n_Py + n_Pz*n_Pz + c_lepWFlags.M()*c_lepWFlags.M() );
	c_lep.SetPx( n_Px );	// Change Px 				
	c_lep.SetPy( n_Py ); 	// Change Py 	
	c_lep.SetPz( n_Pz ); 	// Change Pz 
	c_lep.SetE(  n_E  ); 	// Change E 
	// Propagate to MissPx and MissPy
	delPx = c_lepWFlags.Px() - n_Px; 
	delPy = c_lepWFlags.Py() - n_Py;			
	in_mpx[0] += delPx; in_mpx[1] += delPx; // correct miss(Px,Py) neutrino 1
	in_mpy[0] += delPy; in_mpy[1] += delPy; // correct miss(Px,Py) neutrino 2
	in_mpz[0] += 0.   ; in_mpz[1] += 0.;	// initialize neutrinos Pz to zero

	// _______________________________
	// _______z_bj____________________
	// _______________________________
	n_Px = z_bjWFlags.Px() * ( 1. + random.Gaus( 0., Sx_j ) );
	n_Py = z_bjWFlags.Py() * ( 1. + random.Gaus( 0., Sy_j ) );
	n_Pz = z_bjWFlags.Pz() * ( 1. + random.Gaus( 0., Sz_j ) );
	n_Pt = z_bjWFlags.Pt() * ( 1. + random.Gaus( 0., St_j ) );
	n_E  = z_bjWFlags.E()  * ( 1. + random.Gaus( 0., Se_j ) );
	// Recalculate z_bj
	n_E = sqrt ( n_Px*n_Px + n_Py*n_Py + n_Pz*n_Pz + z_bjWFlags.M()*z_bjWFlags.M() );
	z_bj.SetPx( n_Px );	// Change Px 				
	z_bj.SetPy( n_Py ); 	// Change Py 	
	z_bj.SetPz( n_Pz ); 	// Change Pz 
	z_bj.SetE(  n_E  ); 	// Change E 
	// Propagate to MissPx and MissPy
	delPx = z_bjWFlags.Px() - n_Px; 
	delPy = z_bjWFlags.Py() - n_Py;			
	in_mpx[0] += delPx; in_mpx[1] += delPx; // correct miss(Px,Py) neutrino 1
	in_mpy[0] += delPy; in_mpy[1] += delPy; // correct miss(Px,Py) neutrino 2
	in_mpz[0] += 0.   ; in_mpz[1] += 0.;	// initialize neutrinos Pz to zero


	// _______________________________
	// _______c_bj____________________
	// _______________________________
	n_Px = c_bjWFlags.Px() * ( 1. + random.Gaus( 0., Sx_j ) );
	n_Py = c_bjWFlags.Py() * ( 1. + random.Gaus( 0., Sy_j ) );
	n_Pz = c_bjWFlags.Pz() * ( 1. + random.Gaus( 0., Sz_j ) );
	n_Pt = c_bjWFlags.Pt() * ( 1. + random.Gaus( 0., St_j ) );
	n_E  = c_bjWFlags.E()  * ( 1. + random.Gaus( 0., Se_j ) );
	// Recalculate c_bj
	n_E = sqrt ( n_Px*n_Px + n_Py*n_Py + n_Pz*n_Pz + c_bjWFlags.M()*c_bjWFlags.M() );
	c_bj.SetPx( n_Px );	// Change Px 				
	c_bj.SetPy( n_Py ); 	// Change Py 	
	c_bj.SetPz( n_Pz ); 	// Change Pz 
	c_bj.SetE(  n_E  ); 	// Change E 
	// Propagate to MissPx and MissPy
	delPx = c_bjWFlags.Px() - n_Px; 
	delPy = c_bjWFlags.Py() - n_Py;			
	in_mpx[0] += delPx; in_mpx[1] += delPx; // correct miss(Px,Py) neutrino 1
	in_mpy[0] += delPy; in_mpy[1] += delPy; // correct miss(Px,Py) neutrino 2
	in_mpz[0] += 0.   ; in_mpz[1] += 0.;	// initialize neutrinos Pz to zero


	// ---------------------------------------
	// Define TLorentzVectors for (b,l) system
	// ---------------------------------------
	z_bl = z_bj + z_lep;
	c_bl = c_bj + c_lep;
}

// Getters

TLorentzVector DilepInput::getZlep (void) const {
	return z_lep;
}

TLorentzVector DilepInput::getClep (void) const {
	return c_lep;
}

TLorentzVector DilepInput::getZbj (void) const {
	return z_bj;
}

TLorentzVector DilepInput::getCbj (void) const {
	return c_bj;
}

TLorentzVector DilepInput::getZbl (void) const {
	return z_bl;
}

TLorentzVector DilepInput::getCbl (void) const {
	return c_bl;
}

TLorentzVectorWFlags DilepInput::getZlepW (void) const {
	return z_lepWFlags;
}

TLorentzVectorWFlags DilepInput::getClepW (void) const {
	return c_lepWFlags;
}

TLorentzVectorWFlags DilepInput::getZbjW (void) const {
	return z_bjWFlags;
}

TLorentzVectorWFlags DilepInput::getCbjW (void) const {
	return c_bjWFlags;
}

TLorentzVectorWFlags DilepInput::getJet1HiggsW (void) const {
	return jet1_HiggsWFlags;
}

TLorentzVectorWFlags DilepInput::getJet2HiggsW (void) const {
	return jet2_HiggsWFlags;
}

double DilepInput::getMissPx (void) const {
	return MissPx;
}

double DilepInput::getMissPy (void) const {
	return MissPy;
}

double DilepInput::getInMpx (int x) const {
	return in_mpx[x];
}

double DilepInput::getInMpy (int x) const {
	return in_mpy[x];
}

double DilepInput::getInMpz (int x) const {
	return in_mpz[x];
}

double DilepInput::getTmass (int x) const {
	return t_mass[x];
}

double DilepInput::getWmass (int x) const {
	return w_mass[x];
}

int DilepInput::getHasSol (void) const {
	return hasSolution;
}

vector<myvector> DilepInput::getResult (void) const {
	return result;
}

// Setters

void DilepInput::setHasSol (int x) {
	hasSolution = x;
}

void DilepInput::setResult (vector<myvector> *x) {
	result = *x;
}
