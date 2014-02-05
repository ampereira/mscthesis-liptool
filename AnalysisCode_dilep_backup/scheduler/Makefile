#!/bin/sh

CXX = g++
CXXFLAGS = -O3 -fopenmp

default: scheduler

app.o: app.cpp app.h
	$(CXX) $(CXXFLAGS) -c app.cpp

scheduler: scheduler.cpp app.o
	$(CXX) $(CXXFLAGS) scheduler.cpp app.o -o scheduler
