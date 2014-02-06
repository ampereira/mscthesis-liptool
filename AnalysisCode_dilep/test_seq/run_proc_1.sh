#!/bin/sh

touch log.txt

for i in 1
do
	ruby kbest.rb 1 $i
	ruby kbest.rb 2 $i
	ruby kbest.rb 4 $i
	ruby kbest.rb 8 $i
	ruby kbest.rb 16 $i
	ruby kbest.rb 32 $i
	ruby kbest.rb 64 $i
	ruby kbest.rb 128 $i
	ruby kbest.rb 256 $i
	ruby kbest.rb 512 $i
#ruby kbest.rb 1024 $i

	echo "acabou remessa $i" >> log.txt
done

