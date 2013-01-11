################################################################################
# Make file for LipMiniAnalysis
################################################################################

SHELL = /bin/sh

DEFINES = -Dextname

CXX        = g++
LD         = g++

ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
#ROOTGLIBS  = $(shell root-config --glibs) -lMinuit -lHtml -lEG -lPhysics -lTreePlayer
ROOTGLIBS  = $(shell root-config --glibs) -lMinuit -lEG -lPhysics -lTreePlayer

CXXFLAGS   = -O3 $(ROOTCFLAGS)
#CXXFLAGS   = -g $(ROOTCFLAGS)

LIBS       = $(ROOTLIBS)
GLIBS      = $(ROOTGLIBS)

INCLUDES = -I$(incdir) -I$(ROOTSYS)/include -I$(ROOTCOREDIR)/include

################################################################################
# analysis code
################################################################################

LipMiniAnalysis = ../LipMiniAnalysis
COMMONANALYSISCODE = define_samples_simulation_MarkOwen.cxx define_samples_data_MarkOwen.cxx UserCommandLineOptions.cxx

################################################################################
# Rules
################################################################################

all : neut.o ttH_dilep

neut.o: neut.cxx myvector.h neut.h
	$(CXX) $(CXXFLAGS) -I$(INCLUDES) -c neut.cxx

ttH_dilep : ttH_dilep.cxx ttH_dilep.h $(COMMONANALYSISCODE) neut.o $(LipMiniAnalysis)/libLipMiniAnalysis.a
	$(CXX) $(CXXFLAGS) -o ttH_dilep -I$(INCLUDES) ttH_dilep.cxx neut.o \
	-L$(LipMiniAnalysis) -lLipMiniAnalysis \
	$(LIBS) $(GLIBS) -lMinuit -lPhysics

clean:
	rm -f neut.o ttH_dilep.o ttH_dilep
