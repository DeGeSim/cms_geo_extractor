Extracts the geometry of cells in the HGCal from CMSSW.

```
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh 

export CMSVER=CMSSW_11_3_0_pre5
export SCRAM_ARCH=slc7_amd64_gcc900

alias cmsenv='eval `scramv1 runtime -sh`'


cmsrel $CMSVER
mv $CMSVER geo_ex

cd geo_ex/src/
mkdir output
cmsenv
scram b ProjectRename
scram b clean

git clone git@github.com:mova/geo_extractor.git .

scram b && cmsRun EDAnalyzers/GeoExtractor/python/ConfFile_cfg.py
```