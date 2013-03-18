// Utilities for time measurement and other stuffs
#include "utilities.h"

using namespace std;


int dilep_iterations;
int num_threads;
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
	// Defines the number of threads
	void defineNumThreads (void) {
		char *num = getenv("NUM_THREADS");

		if (num != NULL) {
			num_threads = atoi(num);
			cout << "\033[0;32mRuning analysis with " << num_threads;
			cout << " threads\033[0m" << endl << endl;
		} else {
			num_threads = DEFAULT_THREADS;
			cout << "\033[0;31mNumber of threads not defined!" << endl;
			cout << "Running analysis for the default value of ";
			cout << num_threads << " threads\033[0m" << endl << endl;
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

		#ifdef OMP
		
		ttDKF_Best_Sol reduce (ttDKF_Best_Sol list[]) {
			unsigned size = omp_get_num_threads();
			float tdp = log2f(size);
			unsigned depth = (tdp > (int) tdp) ? tdp + 1 : tdp;
			unsigned tid = omp_get_thread_num();
			ttDKF_Best_Sol *list2;

			memcpy(list2, list, sizeof(ttDKF_Best_Sol) * size);

			/*#pragma omp critical
			{
				ofstream of ("sols.txt", fstream::app);
				of << tid << " - " << list[tid].getProb() << endl;
				of.close();
			}*/

			#pragma omp barrier
			
			// Cycle through all levels of the reduction tree
			for (unsigned i = 0; i < depth; ++i) {
				// First level of the tree is a special scenario
				if (i == 0) {
					// Checks if there is any thread to the right
					if ((tid % 2) == 0 && (tid + 1) < size)
						if (list[tid].getProb() < list[tid + 1].getProb())
							list[tid] = list[tid + 1];
				} else {
					unsigned stride = 2 * i;

					if ((tid % stride) == 0 && (tid + stride) < size)
						if (list[tid].getProb() < list[tid + stride].getProb())
							list[tid] = list[tid + stride];
				}
				#pragma omp barrier
			}
			/*#pragma omp master
			{
				ofstream of ("sols.txt", fstream::app);
				of << "Best - " << list[0].getProb() << endl;
				of.close();
			}
			#pragma omp barrier
			exit(0);*/
			
			return list[0];
		}
		
		#endif
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
