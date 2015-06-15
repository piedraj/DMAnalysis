1. Get ROOT
====

    ssh -Y gridui.ifca.es -o ServerAliveInterval=240
    cd /gpfs/csic_projects/cms/$USER
    mkdir b2g
    cd b2g

    source /cvmfs/cms.cern.ch/cmsset_default.sh
    export SCRAM_ARCH=slc6_amd64_gcc491
    cmsrel CMSSW_7_4_4_ROOT5
    cd CMSSW_7_4_4_ROOT5/src
    cmsenv

    scram b -j 8
     cd DMAnalysis/DMAnalysisTreeMaker/test/
     cmsRun ../python/ConfFile_cfg.py
