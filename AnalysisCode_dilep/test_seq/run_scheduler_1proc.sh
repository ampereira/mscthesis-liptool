#!/bin/sh

# iterates through the number of threads
for j in 1 2 4 8 16 32 64
do
	# iterates through the dilep variations
	for i in 1 2 4 8 16 32 64 128 256 512 1024
	do
		ruby kbest_scheduler.rb $i $j 1
	done
done
