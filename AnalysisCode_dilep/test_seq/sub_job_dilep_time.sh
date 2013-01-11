#!/bin/bash
#PBS -l walltime=5:00:00
#PBS -j oe
#PBS -N ttH_dilep_seq
#PBS -l nodes=1:ppn=24
#PBS -M lazeroth89@gmail.com
#PBS -m ae
#PBS -V

cd ~/ATLAS/analises/ttH/AnalysisCode_dilep/test_seq

./dilep_time.sh

