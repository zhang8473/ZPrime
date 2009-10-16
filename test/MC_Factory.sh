#!/bin/bash

cd ~/CMSSW_3_2_6/src/
eval `scramv1 runtime -sh`
cd ~/
mv Produce*.py /tmp/zhangjin/
cd /tmp/zhangjin/
for i in `ls Produce*.py`; do
	echo Running $i
	cmsRun -p $i
        if [ $1 == "local" ]; then 
	    cp *.root ~/
	else
	    rfcp *.root /castor/cern.ch/user/z/zhangjin/
        fi
	rm *.root
done
