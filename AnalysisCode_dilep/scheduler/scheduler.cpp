#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <omp.h>
#include <unistd.h>

#include "app.h"

using namespace std;

void setup (unsigned its, unsigned threads) {
	string c1 = "export DILEP_ITER=" + its;
	string c2 = "export NUM_THREADS=" + threads;

	system(c1.c_str());
	system(c2.c_str());
}

// necessario alterar para receber os argumentos como input

int main (void) {
	// inputs
	unsigned iterations = 1;
	unsigned num_threads = 4;
	unsigned num_parallel_apps = 2;
	unsigned num_total_runs = 10;
	string app ("ttH_dilep_omp");
	string inputs ("--OutputFileName=ttH125_dilepbb_em --SetSystematicsFileName=../../RefSys/Ref.txt --Sample=901 --User=\"CutTriggerEleMuo=1\" --User=\"lepSample=23\"");


	vector<App> applications;

	// build apps vector
	for (int i = 0; i < num_total_runs; ++i) {
		App a (app, inputs);
		applications.push_back(a);
	}

	setup();

	omp_set_num_threads(num_parallel_apps);

	#pragma omp parallel
	{
		#pragma omp for schedule(dynamic)
		for (int i = 0; i < applications.size(); ++i) {
			applications[i].run();
		}
	}

}