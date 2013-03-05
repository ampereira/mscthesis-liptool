#ifndef UTILITIES
#define UTILITIES

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/time.h>

#define APP_TIME_RESOLUTION 1000000.0	// time measuring resolution (us)
#define DEFAULT_DILEP 1					// default number of dilep iterations
#define DEFAULT_THREADS 1				// default number of threads


namespace Dilep {

	// timer functions
	long long int startTimer (void);
	long long int stopTimer (long long int);
	void printTimer (void);
}

namespace ttH {

	void defineDilepIterations (void);
	void defineNumThreads (void);
	long long int startTimer (void);
	long long int stopTimer (long long int);

	namespace KinFit {
		long long int startTimer (void);
		long long int stopTimer (long long int);
		void printTimer (void);
	}
}

#endif