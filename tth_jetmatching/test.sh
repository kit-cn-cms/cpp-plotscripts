#!/bin/bash

#export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
#source $VO_CMS_SW_DIR/cmsset_default.sh
#cd /afs/desy.de/user/h/hmildner/CMSSW_7_2_3/src
#eval `scram runtime -sh`
#cd -

export FILENAMES=" ../../data/tthbb_full.root "
export OUTFILENAME="test.root"
./a.out
