#!/bin/bash
#PBS -l walltime=15:00:00
#PBS -j oe
#PBS -N ttH_dilep_call
#PBS -l nodes=1:ppn=1:hex
#PBS -M lazeroth89@gmail.com
#PBS -m ae
#PBS -V

cd ~/ATLAS/analises/ttH_new/LipTool/AnalysisCode_dilep/test_seq

export DILEP_ITER=128

#valgrind --tool=callgrind ./ttH_dilep_old_seq --OutputFileName=ttH125_dilepbb_em --SetSystematicsFileName=../../RefSys/Ref.txt --Sample=901 --User="CutTriggerEleMuo=1" --User="lepSample=23"

valgrind --tool=callgrind ./ttH_dilep_seq --OutputFileName=ttH125_dilepbb_em --SetSystematicsFileName=../../RefSys/Ref.txt --Sample=901 --User="CutTriggerEleMuo=1" --User="lepSample=23"
