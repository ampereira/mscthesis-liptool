#!/bin/sh

CXX = g++
CXXFLAGS = -Wall -Wextra

default: scheduler_bind

debug: CXXFLAGS += -g -DVERBOSE
debug: scheduler_bind

optim: CXXFLAGS += -O3 -mavx
optim: scheduler_bind

app.o: app.cpp app.h
	$(CXX) $(CXXFLAGS) -c app.cpp

scheduler: scheduler.cpp app.o
	$(CXX) $(CXXFLAGS) -fopenmp scheduler.cpp app.o -o scheduler

scheduler_bind: scheduler.cpp app.o
	$(CXX) $(CXXFLAGS) -pthread scheduler_core_bind.cpp app.o -o scheduler

clean:
	rm -rf *.o scheduler
