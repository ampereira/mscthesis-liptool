#!/bin/bash
#PBS -l walltime=12:00:00
#PBS -j oe
#PBS -N ttH_dilep_seq
#PBS -l nodes=1:ppn=24:tesla
#PBS -M lazeroth89@gmail.com
#PBS -m ae
#PBS -V

cd ATLAS/analises/ttH_new/LipTool/AnalysisCode_dilep/test_seq

./run_sandy.sh

