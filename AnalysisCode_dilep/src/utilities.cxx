// Utilities for time measurement and other stuffs
#include "utilities.h"
#define SEED 123456.0

using namespace std;


int dilep_iterations;
int num_threads;
static TRandom3 *t_rnd;
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

		t_rnd = new TRandom3 [num_threads];

		for (int i = 0; i < num_threads; ++i)
			t_rnd[i].SetSeed(SEED + i);
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
			std::queue<unsigned> on_hold;


			#pragma omp barrier
			
			// Cycle through all levels of the reduction tree
			for (unsigned i = 0; i < depth; ++i) {
				// First level of the tree is a special case
				if (i == 0) {
					// Checks if there is any thread to the right
					if ((tid % 2) == 0) {
						// If the tid has the partner
						if (tid + 1 < size) {
							if (list[tid].getProb() < list[tid + 1].getProb())
								list[tid] = list[tid + 1];
						} else {
							// If there is no partner the thread is put on hold
							on_hold.push(tid);
						}
					}
				} else {
					unsigned stride = pow(2, i + 1);

					if ((tid % stride) == 0) {
						if (tid + stride / 2 < size) {
							if (list[tid].getProb() < list[tid + stride / 2].getProb())
								list[tid] = list[tid + stride / 2];
						} else {
							// If there is any thread on hold it is paired with the current
							if (!on_hold.empty()) {
								unsigned remain = on_hold.front();
								if (list[tid].getProb() < list[remain].getProb())
									list[tid] = list[remain];
								on_hold.pop();
							} else {
								on_hold.push(tid);
							}
						}
					}
				}
				// To ensure memory consistency
				#pragma omp barrier
			}
			
			/*#pragma omp master
			{
				ttDKF_Best_Sol best = list2[0];
				for (int i = 1; i < size; ++i) {
					if (list2[i].getProb() < best.getProb())
							best = list2[i];
				}

				if (best.getProb() != list[0].getProb()) {
					ofstream of ("sols.txt", fstream::app);
					of << "Best - " << best.getProb() << endl;
					of << "Sol - " << list[0].getProb() << endl;
					of.close();
				}
			}*/
			
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
