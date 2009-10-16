#!/bin/bash

cd ~/CMSSW_3_2_6/src/
eval `scramv1 runtime -sh`
cd ~/
cmsRun CSCSkim_trial_cfg.py
cd /tmp/zhangjin
rfcp *.root /castor/cern.ch/user/z/zhangjin/
rm *.root
