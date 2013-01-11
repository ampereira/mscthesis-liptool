#!/bin/bash 
################################################################################
# Script used to run all data (both streams Egamma,Muon)
# Command:    antonioonofre> run_ALL_MY_DATA_ttH.sh 
################################################################################

echo
echo ======== Running Code Locally ===========
# 	__m+m channel____________ 
. TTbarLF_Alpgen_dilep_mm.sh
# 	__e+e channel____________
. TTbarLF_Alpgen_dilep_ee.sh
# 	__e+m channel____________ 
. TTbarLF_Alpgen_dilep_em.sh

