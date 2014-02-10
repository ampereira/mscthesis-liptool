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
vector<string> setInputs (unsigned processes, unsigned threads) {
	vector<string> inputs;

	switch (processes) {
		case 1: break;
		case 2: switch (threads) {
				// 1 CPU
				case 1: inputs.push_back("\"0\"");
						inputs.push_back("\"1\"");
						break;
				case 2: inputs.push_back("\"0 1\"");
						inputs.push_back("\"2 3\"");
						break;
				case 4: inputs.push_back("\"0 1 2 3\"");
						inputs.push_back("\"4 5 6 7\"");
						break;
				// 2 CPUs
				case 8: inputs.push_back("\"0 1 2 3 4 5 6 7\"");      
						inputs.push_back("\"8 9 10 11 12 13 14 15\"");
						break;
				case 16:inputs.push_back("\"0 1 2 3 4 5 6 7 16 17 18 19 20 21 22 23\"");
						inputs.push_back("\"8 9 10 11 12 13 14 15 24 25 26 27 28 29 30 31\"");
						break;
				}
				break;
		case 4: switch (threads) {
				// 1 CPU
				case 1: inputs.push_back("\"0\"");
						inputs.push_back("\"1\"");
						inputs.push_back("\"2\"");
						inputs.push_back("\"3\"");
						break;
				case 2: inputs.push_back("\"0 1\"");
						inputs.push_back("\"2 3\"");
						inputs.push_back("\"4 5\"");
						inputs.push_back("\"6 7\"");
						break;
				// 2 CPUs
				case 4: inputs.push_back("\"0 1 2 3\"");
						inputs.push_back("\"4 5 6 7\"");
						inputs.push_back("\"8 9 10 11\"");
						inputs.push_back("\"12 13 14 15\"");
						break;
				case 8: inputs.push_back("\"0 1 2 3 16 17 18 19\"");
						inputs.push_back("\"4 5 6 7 20 21 22 23\"");
						inputs.push_back("\"8 9 10 11 24 25 26 27\"");
						inputs.push_back("\"12 13 14 15 28 29 30 31\"");
						break;
				}
				break;
		case 8: switch (threads) {
				// 1 CPU
				case 1: inputs.push_back("\"0\"");
						inputs.push_back("\"1\"");
						inputs.push_back("\"2\"");
						inputs.push_back("\"3\"");
						inputs.push_back("\"4\"");
						inputs.push_back("\"5\"");
						inputs.push_back("\"6\"");
						inputs.push_back("\"7\"");
						break;
				// 2 CPUs
				case 2: inputs.push_back("\"0 1\"");
						inputs.push_back("\"2 3\"");
						inputs.push_back("\"4 5\"");
						inputs.push_back("\"6 7\"");
						inputs.push_back("\"8 9\"");
						inputs.push_back("\"10 11\"");
						inputs.push_back("\"12 13\"");
						inputs.push_back("\"14 15\"");
						break;
				case 4: inputs.push_back("\"0 1 16 17\"");
						inputs.push_back("\"2 3 18 19\"");
						inputs.push_back("\"4 5 20 21\"");
						inputs.push_back("\"6 7 22 23\"");
						inputs.push_back("\"8 9 24 25\"");
						inputs.push_back("\"10 11 26 27\"");
						inputs.push_back("\"12 13 28 29\"");
						inputs.push_back("\"14 15 30 31\"");
						break;
				}
				break;
		case 16:switch (threads) {
				// 2 CPUs
				case 1: inputs.push_back("\"0\"");
						inputs.push_back("\"1\"");
						inputs.push_back("\"2\"");
						inputs.push_back("\"3\"");
						inputs.push_back("\"4\"");
						inputs.push_back("\"5\"");
						inputs.push_back("\"6\"");
						inputs.push_back("\"7\"");
						inputs.push_back("\"8\"");
						inputs.push_back("\"9\"");
						inputs.push_back("\"10\"");
						inputs.push_back("\"11\"");
						inputs.push_back("\"12\"");
						inputs.push_back("\"13\"");
						inputs.push_back("\"14\"");
						inputs.push_back("\"15\"");
						break;
				case 2: inputs.push_back("\"0 16\"");
						inputs.push_back("\"1 17\"");
						inputs.push_back("\"2 18\"");
						inputs.push_back("\"3 19\"");
						inputs.push_back("\"4 20\"");
						inputs.push_back("\"5 21\"");
						inputs.push_back("\"6 22\"");
						inputs.push_back("\"7 23\"");
						inputs.push_back("\"8 24\"");
						inputs.push_back("\"9 25\"");
						inputs.push_back("\"10 26\"");
						inputs.push_back("\"11 27\"");
						inputs.push_back("\"12 28\"");
						inputs.push_back("\"13 29\"");
						inputs.push_back("\"14 30\"");
						inputs.push_back("\"15 31\"");
						break;
				} break;
		case 32:inputs.push_back("\"0\"");
				inputs.push_back("\"1\"");
				inputs.push_back("\"2\"");
				inputs.push_back("\"3\"");
				inputs.push_back("\"4\"");
				inputs.push_back("\"5\"");
				inputs.push_back("\"6\"");
				inputs.push_back("\"7\"");
				inputs.push_back("\"8\"");
				inputs.push_back("\"9\"");
				inputs.push_back("\"10\"");
				inputs.push_back("\"11\"");
				inputs.push_back("\"12\"");
				inputs.push_back("\"13\"");
				inputs.push_back("\"14\"");
				inputs.push_back("\"15\"");
				inputs.push_back("\"16\"");
				inputs.push_back("\"17\"");
				inputs.push_back("\"18\"");
				inputs.push_back("\"19\"");
				inputs.push_back("\"20\"");
				inputs.push_back("\"21\"");
				inputs.push_back("\"22\"");
				inputs.push_back("\"23\"");
				inputs.push_back("\"24\"");
				inputs.push_back("\"25\"");
				inputs.push_back("\"26\"");
				inputs.push_back("\"27\"");
				inputs.push_back("\"28\"");
				inputs.push_back("\"29\"");
				inputs.push_back("\"30\"");
				inputs.push_back("\"31\"");
				break;
	}

	return inputs;
}

int main (int argc, char **argv) {
	// inputs
	unsigned num_threads;
	unsigned num_parallel_apps;
	unsigned num_total_runs = 32;	// number of files
	string app ("run_aff.sh");
	vector<string> inputs;

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


	num_total_runs = num_parallel_apps;


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
