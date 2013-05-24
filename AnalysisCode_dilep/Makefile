################################################################################
# Make file for LipMiniAnalysis
################################################################################

SHELL = /bin/sh

DEFINES = -Dextname

NVCC = nvcc
CXX  = g++
LD   = g++

ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
#ROOTGLIBS  = $(shell root-config --glibs) -lMinuit -lHtml -lEG -lPhysics -lTreePlayer
ROOTGLIBS  = $(shell root-config --glibs) -lMinuit -lEG -lPhysics -lTreePlayer

NVCCFLAGS  = -arch=sm_20 -DCUDA 
CXXFLAGS   = -O3 $(ROOTCFLAGS)
#CXXFLAGS   = -g -DCUDA $(ROOTCFLAGS)
# specific time measurement
#CXXFLAGS  += -DMEASURE_KINFIT #-DMEASURE_DILEP

LIBS       = $(ROOTLIBS)
GLIBS      = $(ROOTGLIBS)

INCLUDES = -I$(incdir) -I$(ROOTSYS)/include -I$(ROOTCOREDIR)/include

################################################################################
# analysis code
################################################################################

LipMiniAnalysis = ../LipMiniAnalysis
COMMONANALYSISCODE = src/define_samples_simulation_MarkOwen.cxx src/define_samples_data_MarkOwen.cxx src/UserCommandLineOptions.cxx
LIBDIRS = -L$(LipMiniAnalysis) -lLipMiniAnalysis $(LIBS) $(GLIBS)

################################################################################
# Rules
################################################################################

seq: CXXFLAGS += -DSEQ
seq: ttH_dilep_seq

omp: CXXFLAGS += -DOMP -fopenmp
omp: ttH_dilep_omp

cuda: NVCCFLAGS += $(filter-out -pthread,$(CXXFLAGS))
cuda: LIBDIRS2 = $(filter-out -pthread,$(LIBDIRS))
cuda: LIBDIRSCUDA = $(filter-out -rdynamic,$(LIBDIRS2))
cuda: CXXFLAGS += -DCUDA
cuda: ttH_dilep_cuda

papi: CXXFLAGS += -DPAPI
papi: ttH_dilep_papi

sse: CXXFLAGS += -DSSE -mfpmath=sse -msse4.2
sse: ttH_dilep_sse

dilep_time: CXXFLAGS += -DMEASURE_DILEP
dilep_time: seq


default: seq
all: seq omp cuda sse

build/dilep_input.o: src/dilep_input.cxx src/dilep_input.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c src/dilep_input.cxx -o build/dilep_input.o

build/utilities.o: src/utilities.cxx src/utilities.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c src/utilities.cxx -o build/utilities.o

# Sequential application
build/neut_seq.o: src/seq/neut.cxx src/myvector.h src/seq/neut.h
	$(CXX) -Wall -Wextra $(CXXFLAGS) $(INCLUDES) -c src/seq/neut.cxx -o build/neut_seq.o

build/neut_cuda.o: src/cuda/neut.cxx src/myvector.h src/cuda/neut.cuh
	$(NVCC) $(NVCCFLAGS) $(INCLUDES) -c src/cuda/neut.cu -o build/neut_cuda.o
#	$(CXX) $(CXXFLAGS) $(INCLUDES) -c src/cuda/neut.cxx -o build/neut_cuda.o

ttH_dilep_seq: src/ttH_dilep_old.cxx src/ttH_dilep.h $(COMMONANALYSISCODE) build/dilep_input.o build/utilities.o build/neut_seq.o $(LipMiniAnalysis)/libLipMiniAnalysis.a
	$(CXX) $(CXXFLAGS) -o bin/ttH_dilep_seq $(INCLUDES) src/ttH_dilep_old.cxx build/neut_seq.o build/dilep_input.o build/utilities.o $(LIBDIRS)

# OpenMP application
#build/neut_omp.o: src/omp/neut.cxx src/myvector.h src/omp/neut.h
#	$(CXX) -Wall -Wextra $(CXXFLAGS) $(INCLUDES) -c src/omp/neut.cxx -o build/neut_omp.o

build/ttDKF_Best_Sol.o: src/ttDKF_Best_Sol.cxx src/ttDKF_Best_Sol.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c src/ttDKF_Best_Sol.cxx -o build/ttDKF_Best_Sol.o

ttH_dilep_omp: src/ttH_dilep_omp.cxx src/ttH_dilep.h $(COMMONANALYSISCODE) build/dilep_input.o build/utilities.o build/neut_seq.o build/ttDKF_Best_Sol.o $(LipMiniAnalysis)/libLipMiniAnalysis.a
	$(CXX) $(CXXFLAGS) -o bin/ttH_dilep_omp $(INCLUDES) src/ttH_dilep_omp.cxx build/neut_seq.o build/dilep_input.o build/utilities.o build/ttDKF_Best_Sol.o $(LIBDIRS)

ttH_dilep_cuda: src/ttH_dilep_cuda.cxx src/ttH_dilep.h $(COMMONANALYSISCODE) build/dilep_input.o build/utilities.o build/neut_cuda.o build/ttDKF_Best_Sol.o $(LipMiniAnalysis)/libLipMiniAnalysis.a
	$(NVCC) $(NVCCFLAGS) -o bin/ttH_dilep_cuda $(INCLUDES) src/ttH_dilep_cuda.cu build/neut_cuda.o build/dilep_input.o build/utilities.o build/ttDKF_Best_Sol.o $(LIBDIRSCUDA)
#	$(CXX) $(CXXFLAGS) -o bin/ttH_dilep_cuda $(INCLUDES) src/ttH_dilep_cuda.cxx build/neut_cuda.o build/dilep_input.o build/utilities.o build/ttDKF_Best_Sol.o $(LIBDIRSCUDA)

clean:
	rm -rf build/* bin/ttH_dilep*

