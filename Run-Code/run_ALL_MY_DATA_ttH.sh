#!/bin/bash 
################################################################################
# Script used to run all data (both streams Egamma,Muon)
# Command:    antonioonofre> run_ALL_MY_DATA_ttH.sh 
################################################################################

echo
echo ======== Run Code Locally ===========
# 	__m+m channel____________ 
. run_data_mm.sh
# 	__e+e channel____________
. run_data_ee.sh
# 	__e+m channel____________ 
. run_data_em_muon.sh
. run_data_em_egamma.sh

