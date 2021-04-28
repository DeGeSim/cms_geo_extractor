#!/bin/bash
echo "foobar42 testjob"
echo $HOSTNAME $USER 



source ~/.profile
sendmsg.py "starting full geoextraction"
source ~/utils/bashFunctionCollection.sh 


export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh 

export CMSVER=CMSSW_11_3_0_pre5
export SCRAM_ARCH=slc7_amd64_gcc900

cd /nfs/dust/cms/user/mscham/CMSSW_11_3_0_pre5/src
eval `scramv1 runtime -sh`
scram b clean
scram b
logandrun cmsRun EDAnalyzers/GeoExtractor/python/ConfFile_cfg.py

sendmsg.py "full geoextraction done"
