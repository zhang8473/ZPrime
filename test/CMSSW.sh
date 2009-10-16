#!/bin/bash
if test -e $1; then 
    CMSSWDIR="CMSSW_"3_2_7
else 
    CMSSWDIR="CMSSW_"$1
fi
export SCRAM_ARCH=slc4_ia32_gcc345
GROUP_DIR=/afs/cern.ch/group/zh
CMS_SYS=amd64_linux26
source /afs/cern.ch/cms/sw/cmsset_default.sh
source $GROUP_DIR/group_env.sh
export CVSROOT=:gserver:cmscvs.cern.ch:/cvs_server/repositories/CMSSW
cd ~
if ! test -d $CMSSWDIR; then
        echo -e "\e[42mBuilding New $CMSSWDIR Directory \e[30m"
	scramv1 project CMSSW $CMSSWDIR
fi
cd $CMSSWDIR/src/
eval `scramv1 runtime -sh`
HLTTriggerPath=/afs/cern.ch/cms/sw/slc4_ia32_gcc345/cms/cmssw/$CMSSWDIR/src/HLTrigger/Configuration/python/
AFS_HOME=/afs/cern.ch/user/z/zhangjin/
export PS1='\e[42m[\u]\[\e[0;30m\] \w>'