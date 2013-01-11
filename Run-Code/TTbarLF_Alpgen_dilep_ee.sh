#!/bin/bash 
################################################################################


# -----------------------------------------------------------
#  Copy codes for cross-check
# -----------------------------------------------------------
#mkdir ../Output/LipMiniAnalysis
#mkdir ../Output/AnalysisCode_dilep
#cp -r ../LipMiniAnalysis/*    ../Output/LipMiniAnalysis/
#cp -r ../AnalysisCode_dilep/* ../Output/AnalysisCode_dilep/
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

#============
# MONTE CARLO
#============
# ttH dileptonic
cd ../AnalysisCode_dilep
time ./ttH_dilep --OutputFileName=TTbarLF_Alpgen_dilepbb_ee --SetSystematicsFileName=../RefSys/Ref.txt --Sample=801 --User="CutTriggerEle=1" --User="lepSample=21" 
cd ../Run-Code

