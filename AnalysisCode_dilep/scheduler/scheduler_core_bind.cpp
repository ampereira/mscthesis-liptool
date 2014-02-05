#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <omp.h>
#include <unistd.h>
#include <sstream>
#include <sys/time.h>
#include <fstream>
#include <hwloc.h>

#include "app.h"

unsigned iterations;
unsigned data_size;
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
hwloc_topology_t topo;


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

	data_size = 16;

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
	unsigned index;
	int ret;

	// Prints additional control information
	#ifdef VERBOSE

	pthread_t self = pthread_self();
	hwloc_bitmap_t cpuset = hwloc_bitmap_alloc();
	char *s;

	hwloc_get_cpubind(topo, cpuset, HWLOC_CPUBIND_THREAD);
	hwloc_bitmap_asprintf(&s, cpuset);

	pthread_mutex_lock(&mutex);
	cout << "Thread with id " << self << " bound on " << s << " core" << endl;
	pthread_mutex_unlock(&mutex);
	free(s);

	#endif

	while (true) {
		pthread_mutex_lock(&mutex);
		index = counter++;
		#ifdef VERBOSE
		cout << "Thread with id " << self << " working on index " << index << endl;
		#endif
		pthread_mutex_unlock(&mutex);

		// checks if the dataset is processed
		if (index >= data_size)
			pthread_exit(&ret);

		applications[index].run();
	}
}

// Setups and starts the threads
void setupWorkers (void) {
	threads = new pthread_t [num_parallel_apps];
	thread_ids = new int [num_parallel_apps];
	unsigned nbcores, cores[num_parallel_apps];
	hwloc_obj_t core_binds[num_parallel_apps];

	int err;
	
	err = hwloc_topology_init(&topo);
	assert(!err);
	err = hwloc_topology_load(topo);
	assert(!err);

	switch (num_parallel_apps) {
		case 2 : cores[0] = 0; cores[1] = 8; break;
		case 4 : cores[0] = 0; cores[1] = 4; cores[2] = 8; cores[3] = 12; break;
		case 8 : cores[0] = 0; cores[1] = 2; cores[2] = 4; cores[3] = 6;
				 cores[4] = 8; cores[5] = 10; cores[6] = 12; cores[7] = 14; break;
	}

	nbcores = hwloc_get_nbobjs_by_type(topo, HWLOC_OBJ_CORE);

	for (unsigned i = 0; i < num_parallel_apps; ++i) {
		thread_ids[i] = pthread_create(&threads[i], NULL, worker, NULL);
		core_binds[i] = hwloc_get_obj_by_type(topo, HWLOC_OBJ_CORE, cores[i]);
		//hwloc_set_thread_cpubind(topo, threads[i], core_binds[i]->cpuset, HWLOC_CPUBIND_THREAD);
	}
}

// Waits for thread completion
void destroyWorkers (void) {
	
	for (unsigned i = 0; i < num_parallel_apps; ++i) {
		pthread_join(threads[i], NULL);
	}

	hwloc_topology_destroy(topo);
}


int main (int argc, char **argv) {

	long long int t = startTimer ();
	

	readInputs(argc, argv);

	setupWorkers();

	destroyWorkers();


	stopTimer (t);

	return 0;
}
