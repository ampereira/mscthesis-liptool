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
#define _GNU_SOURCE             /* See feature_test_macros(7) */
       #include <sched.h>

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
	char *buff = NULL;
	char *flag = NULL;

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
	unsigned num_total_runs = 16;	// number of files
	string app ("ttH_dilep_omp");
	string inputs ("--OutputFileName=ttH125_dilepbb_em --SetSystematicsFileName=../../RefSys/Ref.txt --Sample=901 --User=\"CutTriggerEleMuo=1\" --User=\"lepSample=23\"");

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
			cpu_set_t set;
			int proc_num = omp_get_thread_num();
			CPU_ZERO( &set );
			CPU_SET( proc_num, &set );
			printf("proc_num=(%d)\n",proc_num) ;
			if (sched_setaffinity( gettid(), sizeof( cpu_set_t ), &set ))
			{
			    perror( "sched_setaffinity" );
			  //  return NULL;
			}
			//applications[i].run();
		}
	}

	stopTimer (t);

	return 0;

}
