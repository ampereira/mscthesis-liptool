#!/bin/bash 
################################################################################


# -----------------------------------------------------------
#  Copy codes for cross-check
# -----------------------------------------------------------
mkdir ../Output/LipMiniAnalysis
mkdir ../Output/AnalysisCode_dilep
cp -r ../LipMiniAnalysis/*    ../Output/LipMiniAnalysis/
cp -r ../AnalysisCode_dilep/* ../Output/AnalysisCode_dilep/
#
echo
echo
date
echo detected host: `hostname`,  os = `uname`
echo
#
# -----------------------------------------------------------
#  Compile if necessary
# -----------------------------------------------------------
#echo compiling...
#cd ../LipMiniAnalysis
#make 
#cd ../AnalysisCode_dilep
#make
#
echo
echo running analysis...
echo
echo

#
################################################################################
cd ../AnalysisCode_dilep
time ttH_dilep --OutputFileName=ttHdilep_data_Egamma_em_disc --SetSystematicsFileName=../RefSys/Ref.txt --isData=1 --User="alldata_ele=1" --User="CutTriggerEle=1" --User="lepSample=23"
cd ../Run-Code
################################################################################
