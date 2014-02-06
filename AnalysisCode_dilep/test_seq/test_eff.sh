#!/bin/sh

#export MEASURE_APP=1
export NUM_THREADS=1

for i in 1 2 4 8 16 32 64 128 256 512
do
	unset DILEP_ITER
	export DILEP_ITER=$i

	for j in {1..3}
	do
		./ttH_dilep_omp --OutputFileName=ttH125_dilepbb_em --SetSystematicsFileName=../../RefSys/Ref.txt --Sample=901 --User="CutTriggerEleMuo=1" --User="lepSample=23"
	done
done

