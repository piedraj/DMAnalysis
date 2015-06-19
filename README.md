1. Get CMSSW
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


2. Get the material
====

    git clone https://github.com/piedraj/DMAnalysis


3. Run (local)
====

    scram b -j 8
    cd DMAnalysis/DMAnalysisTreeMaker/test/
    cmsRun ConfFile_cfg.py


4. Run with CRAB
====	

    source /cvmfs/cms.cern.ch/crab3/crab.sh
    cd crab/
    python multicrab.py samples_b2g.py
    crab status crab_projects_ferrero/crab_<sample name at "samples_b2g.py">
