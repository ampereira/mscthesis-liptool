#include <vector>
#include <iostream>

using std::cout;
using std::endl;


#ifndef myvector_h
#define myvector_h


//==============================================
class myvector
//==============================================
{
	public:

	myvector(double px, double py, double pz, double npz){apx = px; apy = py; apz = pz; bpz = npz;};
  virtual ~myvector() {} ;

	inline double aPz() const {return bpz;};
	inline double Px() const {return apx;};
	inline double Py() const {return apy;};
	inline double Pz() const {return apz;};

	private:
		double bpz;
		double apx;
		double apy;
		double apz;

};

#endif

