#!/bin/bash
#PBS -l walltime=5:00:00
#PBS -j oe
#PBS -N ttH_dilep_omp1
#PBS -l nodes=1:ppn=24
#PBS -M lazeroth89@gmail.com
#PBS -m ae
#PBS -V

cd ~/ATLAS/analises/ttH/AnalysisCode_dilep/test_omp1

export NUM_THREADS=1

ruby kbest.rb 1
ruby kbest.rb 10
ruby kbest.rb 100
ruby kbest.rb 1000

