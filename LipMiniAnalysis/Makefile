################################################################################
# Make file for LipMiniAnalysis
################################################################################

SHELL = /bin/sh

DEFINES = -Dextname

CXX        = g++
LD         = g++

ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS  = $(shell root-config --glibs) -lMinuit -lHtml -lEG -lPhysics -lTreePlayer

CXXFLAGS   = -Wall -Wextra -O3 $(ROOTCFLAGS)
#CXXFLAGS   = -g $(ROOTCFLAGS)

LIBS       = $(ROOTLIBS)
GLIBS      = $(ROOTGLIBS)

INCLUDES = -I$(incdir) -I$(ROOTSYS)/include -I$(ROOTCOREDIR)/include

################################################################################
# Rules
################################################################################

all :	libLipMiniAnalysis.a

TLorentzVectorWFlags.o: TLorentzVectorWFlags.cxx TLorentzVectorWFlags.h
	$(CXX) $(CXXFLAGS) -I$(INCLUDES) -c TLorentzVectorWFlags.cxx

Ntu.o: Ntu.cxx Ntu.h TLorentzVectorWFlags.o
	$(CXX) $(CXXFLAGS) -I$(INCLUDES) -c Ntu.cxx

MiniTTHReader.o: Ntu.o MiniTTHReader.cxx MiniTTHReader.h
	$(CXX) $(CXXFLAGS) -I$(INCLUDES) -c MiniTTHReader.cxx

LipMiniAnalysis.o: LipMiniAnalysis.cxx LipMiniAnalysis.h TLorentzVectorWFlags.o Ntu.o MiniTTHReader.o
	$(CXX) $(CXXFLAGS) -DMEASURE_READFILE -I$(INCLUDES) -c LipMiniAnalysis.cxx

libLipMiniAnalysis.a: LipMiniAnalysis.o
	ar -r libLipMiniAnalysis.a \
	LipMiniAnalysis.o \
	MiniTTHReader.o \
	Ntu.o \
	TLorentzVectorWFlags.o

clean:
	rm -f *.o libLipMiniAnalysis.a 

# DO NOT DELETE
