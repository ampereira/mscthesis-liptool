#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <omp.h>
#include <unistd.h>
#include <sstream>
#include <sys/time.h>
#include <fstream>
//#include <hwloc.h>
#include <sched.h>
#include <stdio.h>
#include <linux/unistd.h>
#include <sys/syscall.h>
#include <unistd.h>
#include "app.h"

unsigned iterations;


using namespace std;

// Time measurement functions
long long int startTimer (void) {
	long long int time;
	timeval t;

	
	gettimeofday(&t, NULL);
	time = t.tv_sec * 1000000 + t.tv_usec;

	return time;
}

long long int stopTimer (long long int init) {
	timeval t;
	long long int end;

	gettimeofday(&t, NULL);

	ofstream file;
	string filename = "time_scheduler_";

	stringstream ss;
	ss << iterations;

	filename.append(ss.str());
	filename.append(".txt");

	end = t.tv_sec * 1000000 + t.tv_usec;
	end -= init;

	file.precision(15);
	file.open(filename.c_str(), fstream::app);
	file << end << endl;
	file.close();

	return end;
}


void setup (unsigned its, unsigned threads) {
	stringstream a1, a2;
	a1 << its;
	a2 << threads;

	setenv("DILEP_ITER", a1.str().c_str(), true);
	setenv("NUM_THREADS", a2.str().c_str(), true);
	//setenv("OMP_PROC_BIND", "true", true);
}

// No software multithreading for now
// Optimized for 2 x 8core systems with MT
string* setInputs (unsigned processes, unsigned threads) {
	string inputs[processes];

	switch (processes) {
		case 1: break;
		case 2: switch (threads) {
				// 1 CPU
				case 1: inputs[0] = "\"0\"";
						inputs[1] = "\"1\"";
						break;
				case 2: inputs[0] = "\"0 1\"";
						inputs[1] = "\"2 3\"";
						break;
				case 4: inputs[0] = "\"0 1 2 3\"";
						inputs[1] = "\"4 5 6 7\"";
						break;
				// 2 CPUs
				case 8: inputs[0] = "\"0 1 2 3 4 5 6 7\"";      
						inputs[1] = "\"8 9 10 11 12 13 14 15\"";
						break;
				case 16:inputs[0] = "\"0 1 2 3 4 5 6 7 16 17 18 19 20 21 22 23\"";
						inputs[1] = "\"8 9 10 11 12 13 14 15 24 25 26 27 28 29 30 31\"";
						break;
			}
		case 4: switch (threads) {
				// 1 CPU
				case 1: inputs[0] = "\"0\"";
						inputs[1] = "\"1\"";
						inputs[2] = "\"2\"";
						inputs[3] = "\"3\"";
						break;
				case 2: inputs[0] = "\"0 1\"";
						inputs[1] = "\"2 3\"";
						inputs[2] = "\"4 5\"";
						inputs[3] = "\"6 7\"";
						break;
				// 2 CPUs
				case 4: inputs[0] = "\"0 1 2 3\"";
						inputs[1] = "\"4 5 6 7\"";
						inputs[2] = "\"8 9 10 11\"";
						inputs[3] = "\"12 13 14 15\"";
						break;
				case 8: inputs[0] = "\"0 1 2 3 16 17 18 19\"";
						inputs[1] = "\"4 5 6 7 20 21 22 23\"";
						inputs[2] = "\"8 9 10 11 24 25 26 27\"";
						inputs[3] = "\"12 13 14 15 28 29 30 31\"";
						break;
			}
		case 8: switch (threads) {
				// 1 CPU
				case 1: inputs[0] = "\"0\"";
						inputs[1] = "\"1\"";
						inputs[2] = "\"2\"";
						inputs[3] = "\"3\"";
						inputs[4] = "\"4\"";
						inputs[5] = "\"5\"";
						inputs[6] = "\"6\"";
						inputs[7] = "\"7\"";
						break;
				// 2 CPUs
				case 2: inputs[0] = "\"0 1\"";
						inputs[1] = "\"2 3\"";
						inputs[2] = "\"4 5\"";
						inputs[3] = "\"6 7\"";
						inputs[4] = "\"8 9\"";
						inputs[5] = "\"10 11\"";
						inputs[6] = "\"12 13\"";
						inputs[7] = "\"14 15\"";
						break;
				case 4: inputs[0] = "\"0 1 16 17\"";
						inputs[1] = "\"2 3 18 19\"";
						inputs[2] = "\"4 5 20 21\"";
						inputs[3] = "\"6 7 22 23\"";
						inputs[4] = "\"8 9 24 25\"";
						inputs[5] = "\"10 11 26 27\"";
						inputs[6] = "\"12 13 28 29\"";
						inputs[7] = "\"14 15 30 31\"";
						break;
			}
		case 16:switch (threads) {
				// 2 CPUs
				case 1: inputs[0] = "\"0\"";
						inputs[1] = "\"1\"";
						inputs[2] = "\"2\"";
						inputs[3] = "\"3\"";
						inputs[4] = "\"4\"";
						inputs[5] = "\"5\"";
						inputs[6] = "\"6\"";
						inputs[7] = "\"7\"";
						inputs[8] = "\"8\"";
						inputs[9] = "\"9\"";
						inputs[10] = "\"10\"";
						inputs[11] = "\"11\"";
						inputs[12] = "\"12\"";
						inputs[13] = "\"13\"";
						inputs[14] = "\"14\"";
						inputs[15] = "\"15\"";
						break;
				case 2: inputs[0] = "\"0 16\"";
						inputs[1] = "\"1 17\"";
						inputs[2] = "\"2 18\"";
						inputs[3] = "\"3 19\"";
						inputs[4] = "\"4 20\"";
						inputs[5] = "\"5 21\"";
						inputs[6] = "\"6 22\"";
						inputs[7] = "\"7 23\"";
						inputs[8] = "\"8 24\"";
						inputs[9] = "\"9 25\"";
						inputs[10] = "\"10 26\"";
						inputs[11] = "\"11 27\"";
						inputs[12] = "\"12 28\"";
						inputs[13] = "\"13 29\"";
						inputs[14] = "\"14 30\"";
						inputs[15] = "\"15 31\"";
						break;
			}
	}

	return inputs;
}

int main (int argc, char **argv) {
	// inputs
	unsigned num_threads;
	unsigned num_parallel_apps;
	unsigned num_total_runs = 2;	// number of files
	string app ("run_aff.sh");
	string *inputs;

	if (argc < 4) {
		cout << "Not enough arguments" << endl;
		exit(-1);
	}
	stringstream s1, s2, s3;
	s1 << argv[1];
	s2 << argv[2];
	s3 << argv[3];
	s1 >> iterations;
	s2 >> num_threads;
	s3 >> num_parallel_apps;

	setup(iterations, num_threads);
	inputs = setInputs(num_parallel_apps, num_threads);

	cout << endl << "-= Execution specifications =-" << endl;
	cout << "\tFiles: " << num_total_runs << endl;
	cout << "\tProcesses: " << num_parallel_apps << endl;
	cout << "\tThreads: " << num_threads << endl;
	cout << "\tIterations: " << iterations << endl << endl;

	long long int t = startTimer ();
	vector<App> applications;

	// build apps vector
	for (unsigned i = 0; i < num_total_runs; ++i) {
		App a (app, inputs[i]);
		applications.push_back(a);
	}

	omp_set_num_threads(num_parallel_apps);

	#pragma omp parallel
	{
		#pragma omp for schedule(dynamic)
		for (unsigned i = 0; i < applications.size(); ++i) {
			applications[i].run();
		}
	}

	stopTimer (t);

	return 0;

}
