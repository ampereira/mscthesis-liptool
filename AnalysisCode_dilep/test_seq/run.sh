
case $1 in
	"omp")
		currfile=ttH_dilep_omp
		;;
	"mic")
		currfile=ttH_dilep_mic
		;;
	"cuda")
		currfile=ttH_dilep_cuda
		;;
	*)
		currfile=ttH_dilep_seq
		;;
esac

binfile=../bin/$currfile

#if [ $binfile -nt $currfile ]
#then
#	cp $binfile .
#	echo "Using newer version of $currfile"
#fi

export GOMP_CPU_AFFINITY=$2

echo "Running for $GOMP_CPU_AFFINITY affinity"

./$currfile --OutputFileName=ttH125_dilepbb_em --SetSystematicsFileName=../../RefSys/Ref.txt --Sample=901 --User="CutTriggerEleMuo=1" --User="lepSample=23" 
