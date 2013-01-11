#!/bin/bash 
################################################################################
# Script used to run all data (both streams Egamma,Muon)
# Command:    antonioonofre> run_ALL_MY_ttH_Pythia.sh 
################################################################################

echo
echo ======== Running Code Locally ===========
# 	__m+m channel____________ 
. ttH_dilep_mm.sh
# 	__e+e channel____________
. ttH_dilep_ee.sh
# 	__e+m channel____________ 
. ttH_dilep_em.sh

