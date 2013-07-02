#include <cmath>


#define STRIDE2(a,i) a[tid * 2 + i]
#define STRIDE5(a,i) a[tid * 5 + i]


inline void __attribute__((target(mic))) Csqrt(double _ar, double _ai, double _my[]);
inline void __attribute__((target(mic))) cubic(double a[], double rr[], double ri[]);
inline void __attribute__((target(mic))) toz(double k[], double l[], double g[]);
inline void __attribute__((target(mic))) my_qu( double my_in[], double my_val[]);
inline double __attribute__((target(mic))) calcMass(double x, double y, double z, double e);

inline void __attribute__((target(mic))) calc_dilep(double t_mass[], double w_mass[], 
											double in_mpx[], double in_mpy[], double _lep_a[], 
											double _lep_b[], double _bl_a[], double _bl_b[], 
											double nc[], int a[], int size);