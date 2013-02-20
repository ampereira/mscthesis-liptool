currfile=ttH_dilep_seq
binfile=../bin/$currfile

if [ $binfile -nt $currfile ]
then
	cp $binfile .
	echo "Using newer version of $currfile"
fi

time ./$currfile --OutputFileName=ttH125_dilepbb_em --SetSystematicsFileName=../../RefSys/Ref.txt --Sample=901 --User="CutTriggerEleMuo=1" --User="lepSample=23" 
