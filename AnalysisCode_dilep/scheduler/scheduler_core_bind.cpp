#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <omp.h>
#include <unistd.h>
#include <sstream>
#include <sys/time.h>
#include <fstream>
#include <pthread.h>
#include <sys/syscall.h>
#include "app.h"

unsigned iterations;
unsigned data_size;
unsigned index;
unsigned num_threads;
unsigned num_parallel_apps;

// Data structure holding each application information
vector<App> applications;
// Index of the data structure being processed
unsigned counter = 0;
// Mutex for handling concurrency on the counter variable
pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
// Threads to run
pthread_t *threads;
int *thread_ids;
int t_id = 0;


using namespace std;

// Time measurement functions
long long int startTimer (void) {
	long long int time;
	timeval t;

	
	gettimeofday(&t, NULL);
	time = t.tv_sec * 1000000 + t.tv_usec;

	return time;
}

void stopTimer (long long int init) {
	timeval t;
	long long int end;

	gettimeofday(&t, NULL);

	ofstream file;
	string filename = "time_scheduler";

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

// Initial application setup
void setup (unsigned its, unsigned threads) {
	stringstream a1, a2;
	a1 << its;
	a2 << threads;

	setenv("DILEP_ITER", a1.str().c_str(), true);
	setenv("NUM_THREADS", a2.str().c_str(), true);
}

// Reads the input configuration file
void readInputs (int argc, char **argv) {
	string app ("ttH_dilep_omp");
	string inputs ("--OutputFileName=ttH125_dilepbb_em --SetSystematicsFileName=../../RefSys/Ref.txt --Sample=901 --User=\"CutTriggerEleMuo=1\" --User=\"lepSample=23\"");

	data_size = 2;	// number of files to process

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

	// Builds the applications vector
	for (unsigned i = 0; i < data_size; ++i) {
		App a (app, inputs);
		applications.push_back(a);
	}
}

// What each thread will execute
void* worker (void *ptr) {
	int ret;
	pthread_t self = pthread_self();

	// get a thread id, fifo
	pthread_mutex_lock(&mutex);

	int id = t_id;

	pthread_mutex_unlock(&mutex);

	cpu_set_t cpuset;
	CPU_ZERO(&cpuset);
	CPU_SET(id, &cpuset);
	pthread_setaffinity_np(self, sizeof(cpu_set_t), &cpuset);

	while (true) {
		pthread_mutex_lock(&mutex);
		index = counter++;

		cout << "Thread with id " << id << " working on index " << index << endl;
		cout << "Core: " << pthread_getaffinity_np(self, sizeof(cpu_set_t), &cpuset) << endl << endl;

		pthread_mutex_unlock(&mutex);

		// checks if the dataset is processed
		if (index >= data_size)
			pthread_exit(&ret);

		//applications[index].run();
	}
}

// Setups and starts the threads
void setupWorkers (void) {
	threads = new pthread_t [num_parallel_apps];
	thread_ids = new int [num_parallel_apps];

	for (unsigned i = 0; i < num_parallel_apps; ++i) {
		pthread_create(&threads[i], NULL, worker, &cores[i]);
	}
}

// Waits for thread completion
void destroyWorkers (void) {
	
	for (unsigned i = 0; i < num_parallel_apps; ++i) {
		pthread_join(threads[i], NULL);
	}
}


int main (int argc, char **argv) {

	long long int t = startTimer ();
	

	readInputs(argc, argv);

	setupWorkers();

	destroyWorkers();


	stopTimer (t);

	return 0;
}
