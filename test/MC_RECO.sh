#!/bin/bash

CMSEvn=~/CMSSW_3_2_4/src/
cd $CMSEvn
eval `scramv1 runtime -sh`
cp ~/RECO_MC.py /tmp/zhangjin/
cd /tmp/zhangjin/
cmsRun -p RECO_MC.py
rfcp *.root /castor/cern.ch/user/z/zhangjin/
rm *.root
rm *.py
