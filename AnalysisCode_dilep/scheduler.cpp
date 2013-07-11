#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <omp.h>

#include "app.h"

using namespace std;

// necessario alterar para receber os argumentos como input

int main (void) {
	// inputs
	unsigned iterations = 1;
	unsigned num_threads = 4;
	unsigned num_parallel_apps = 2;
	unsigned num_total_runs = 10;
	string app ("ttH_dilep_omp");
	string inputs ("--OutputFileName=ttH125_dilepbb_em --SetSystematicsFileName=../../RefSys/Ref.txt --Sample=901 --User=\"CutTriggerEleMuo=1\" --User=\"lepSample=23\"");


	vector<App> applications (num_total_runs);

	// build apps vector
	for (int i = 0; i < num_total_runs; ++i) {
		App a (app, inputs);
		a.run();
		applications.push_back(a);
	}

}