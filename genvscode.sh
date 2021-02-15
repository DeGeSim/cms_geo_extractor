#!/bin/bash
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh 

export CMSVER=CMSSW_11_2_0_pre10
export SCRAM_ARCH=slc7_amd64_gcc820




echo '"C_Cpp.default.compilerPath": "'$(scram tool tag gcc-cxxcompiler CXX)'",'

echo '"C_Cpp.default.includePath": ['
## compiler libs
echo '"'$(echo $(dirname /cvmfs/cms.cern.ch/slc7_amd64_gcc820/external/gcc/8.2.0-bcolbf/bin/c++ | xargs dirname)/include/c++/*/)'",'
env | grep -E '(CMSSW_FWLITE_INCLUDE_PATH|SHERPA_INCLUDE_PATH|ROOT_INCLUDE_PATH)' | grep -v '_SCRAMRTDEL' |  sd '.*=(.*)' '$1' | tr ':' '\n' | grep -v "/external/" | sd '(.*)' '"$1",' | sed '$d'

## CLHEP...
echo '"'${CLHEP_PARAM_PATH}/include'",' 


## external libs
for extlib in tbb CLHEP fmt boost; do 
env | grep -E '(CMSSW_FWLITE_INCLUDE_PATH|SHERPA_INCLUDE_PATH|ROOT_INCLUDE_PATH)' | grep -v '_SCRAMRTDEL' |  sd '.*=(.*)' '$1' | tr ':' '\n' | grep $extlib | sd '(.*)' '"$1",' | sed '$d'
done



echo ']'

# fd --type f '.*\.h' /cvmfs/cms.cern.ch/$SCRAM_ARCH/cms/cmssw/$CMSVER/src/ | grep "interface" | xargs dirname  | sort -u | sd '(.*)' '"$1",' | sed '$d'
# env | grep -E '(CMSSW_FWLITE_INCLUDE_PATH|SHERPA_INCLUDE_PATH|ROOT_INCLUDE_PATH)' | grep -v '_SCRAMRTDEL' |  sd '.*=(.*)' '$1' | tr ':' '\n' | grep -v "/external/" | sd '(.*)' '"$1",' | sed '$d'