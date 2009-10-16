#!/bin/bash

TagLLum=STARTUP31X_V7
TagHLum=MC_31X_V8
CMSEvn=$CMSSW_BASE/src

cd $CMSEvn
cvs co Configuration/GenProduction
cp ~/PYCard*.py $CMSEvn/Configuration/GenProduction/python/
cd $CMSEvn/Configuration
scramv1 b
cd $CMSEvn/Configuration/GenProduction/python/
for i in `ls PYCard*.py`
do
  ROOTFile_Name=${i%.py}
  ROOTFile_Name=/tmp/zhangjin/${ROOTFile_Name#PYCard}.root
  echo $ROOTFile_Name
  cd $CMSEvn
    if [ $2 = "LowL" ]||[ $2 = "lowl" ]; then
	echo "High Level Trigger 8E29 (under Low Lumiosity)"
	cmsDriver.py Configuration/GenProduction/python/$i -s GEN:ProductionFilterSequence,SIM,DIGI,L1,DIGI2RAW,HLT --conditions FrontierConditions_GlobalTag,$TagLLum::All --fileout $ROOTFile_Name --number $1 --mc --no_exec --datatier 'GEN-SIM-RAW' --eventcontent RAWSIM --processName HLT8E29
    else
	echo "High Level Trigger 1E31 (under High Lumiosity)"
	cmsDriver.py Configuration/GenProduction/python/$i -s GEN:ProductionFilterSequence,SIM,DIGI,L1,DIGI2RAW,HLT:1E31 --fileout $ROOTFile_Name --conditions FrontierConditions_GlobalTag,$TagHLum::All --datatier 'GEN-SIM-RAW'  --eventcontent RAWSIM  -n $1 --no_exec 
    fi
done
rename PYCard Produce PYCard*.py
mv Produce*.py ~/
rm $CMSEvn/Configuration/GenProduction/python/PYCard*.py
