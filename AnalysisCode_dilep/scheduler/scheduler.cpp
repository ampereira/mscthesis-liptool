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

// necessario alterar para receber os argumentos como input
pid_t gettid( void )
{
    return syscall( __NR_gettid );
}

int main (int argc, char **argv) {
	// inputs
	unsigned num_threads;
	unsigned num_parallel_apps;
	unsigned num_total_runs = 2;	// number of files
	string app ("run_aff.sh");
	string inputs[num_parallel_apps];

	string s1 ("\"0 1\"");
	string s2 ("\"2 3\"");

	// Change according to the 
	inputs[0] = s1;
	inputs[1] = s2;

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

	cout << "Execution specifications" << endl;
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
