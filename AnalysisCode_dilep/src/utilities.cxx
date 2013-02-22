// Utilities for time measurement and other stuffs
#include "utilities.h"

using namespace std;


int dilep_iterations;
static long long int totaltime = 0;
static long long int num_measurements = 0;

namespace ttH {

	// Defines the number of dilep iterations per event
	void defineDilepIterations (void) {
		char *num = getenv("DILEP_ITER");

		if (num != NULL) {
			dilep_iterations = atoi(num);
			cout << "\033[0;32mRuning analysis for " << dilep_iterations;
			cout << " dilep iterations per event\033[0m" << endl << endl;
		} else {
			dilep_iterations = DEFAULT_DILEP;
			cout << "\033[0;31mNumber of dilep iterations not defined!" << endl;
			cout << "Running analysis for the default value of ";
			cout << dilep_iterations << " dilep iterations\033[0m" << endl << endl;
		}
	}

	// Time measurement functions
	long long int startTimer (void) {
		char *flag = getenv("MEASURE_APP");
		long long int time;
		timeval t;

		if (flag != NULL) {
			cout << "\033[0;31m Measuring time with microssecond precision! \033[0m" << endl << endl;

			gettimeofday(&t, NULL);
			time = t.tv_sec * APP_TIME_RESOLUTION + t.tv_usec;

			return time;
		} else
			return -1;
	}

	long long int stopTimer (long long int init) {
		timeval t;
		long long int end;
		char *buff = NULL;
		char *flag = NULL;

		gettimeofday(&t, NULL);
		flag = getenv("MEASURE_APP");

		if (flag != NULL) {
			ofstream file;
			string filename = "time_";

			stringstream ss;
			ss << dilep_iterations;

			filename.append(ss.str());
			filename.append(".txt");

			end = t.tv_sec * APP_TIME_RESOLUTION + t.tv_usec;
			end -= init;

			file.precision(15);
			file.open(filename.c_str(), fstream::app);
			file << end << endl;
			file.close();
		} else
			return -1;
	}


	namespace KinFit {
		int measurements = 0;
		int totaltimeKinFit = 0;

		long long int startTimer (void) {
			timeval t;
			long long int time;

			gettimeofday(&t, NULL);
			time = t.tv_sec * APP_TIME_RESOLUTION + t.tv_usec;

			return time;
		}

		long long int stopTimer (long long int init) {
			timeval t;
			long long int end;

			gettimeofday(&t, NULL);

			end = t.tv_sec * APP_TIME_RESOLUTION + t.tv_usec;
			end -= init;

			totaltimeKinFit += end;
			++measurements;
		}

		void printTimer (void) {
			ofstream file;
			string filename = "ttDilepKinFit_time_";

			stringstream ss;
			ss << dilep_iterations;

			filename.append(ss.str());
			filename.append(".txt");

			file.precision(15);
			file.open(filename.c_str(), fstream::app);
			file << "Total time: " << totaltimeKinFit << " usec" << endl;
			file << "Average time: " << totaltimeKinFit / measurements << " usec" << endl;
			file.close();
		}
	}
}

namespace Dilep {

	void printTimer (void) {
		ofstream file;
		string filename = "dileptime_";

		stringstream ss;
		ss << dilep_iterations;

		filename.append(ss.str());
		filename.append(".txt");

		file.precision(15);
		file.open(filename.c_str(), fstream::app);
		file << "Total time: " << totaltime << " usec" << endl;
		file << "Average time: " << totaltime / num_measurements << " usec" << endl;
		file.close();
	}

	long long int startTimer (void) {
		long long int time;
		timeval t;

		gettimeofday(&t, NULL);
		time = t.tv_sec * APP_TIME_RESOLUTION + t.tv_usec;

		return time;
	}
	
	long long int stopTimer (long long int init) {
		timeval t;
		long long int end;

		gettimeofday(&t, NULL);

		end = t.tv_sec * APP_TIME_RESOLUTION + t.tv_usec;
		end -= init;

		totaltime += end;
		++num_measurements;
		
		return end;
	}
}
