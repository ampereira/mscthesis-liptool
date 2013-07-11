#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <omp.h>
#include <unistd.h>
#include <sstream>

#include "app.h"

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
	char *buff = NULL;
	char *flag = NULL;

	gettimeofday(&t, NULL);

	ofstream file;
	string filename = "time_scheduler";

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
}

void setup (unsigned its, unsigned threads) {
	stringstream a1, a2;
	a1 << its;
	a2 << threads;

	setenv("DILEP_ITER", a1.str().c_str(), true);
	setenv("NUM_THREADS", a2.str().c_str(), true);
}

// necessario alterar para receber os argumentos como input

int main (int argc, char **argv) {
	// inputs
	unsigned iterations;
	unsigned num_threads;
	unsigned num_parallel_apps = 2;
	unsigned num_total_runs = 10;
	string app ("ttH_dilep_omp");
	string inputs ("--OutputFileName=ttH125_dilepbb_em --SetSystematicsFileName=../../RefSys/Ref.txt --Sample=901 --User=\"CutTriggerEleMuo=1\" --User=\"lepSample=23\"");

	if (argc < 3) {
		cout << "Not enough arguments" << endl;
		exit(-1);
	}
	stringstream s1, s2;
	s1 << argv[1];
	s2 << argv[2];
	s1 >> iterations;
	s2 >> num_threads;
	
	setup(iterations, num_threads);

	long long int t = startTimer ();
	vector<App> applications;

	// build apps vector
	for (int i = 0; i < num_total_runs; ++i) {
		App a (app, inputs);
		applications.push_back(a);
	}

	omp_set_num_threads(num_parallel_apps);

	#pragma omp parallel
	{
		#pragma omp for schedule(dynamic)
		for (int i = 0; i < applications.size(); ++i) {
			applications[i].run();
		}
	}

	stopTimer (t);

	return 0;

}
