#!/bin/sh

# Receives number of dilep iterations as input

unset DILEP_ITER
unset MEASURE_APP

export DILEP_ITER=$1
export MEASURE_APP=1

for i in {1..4}
do
	./run.sh
done

