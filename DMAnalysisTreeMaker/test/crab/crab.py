from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferOutputs = True
config.General.requestName = 'DM_ferrero04'

config.section_('JobType')
config.JobType.pluginName  = 'Analysis'
config.JobType.psetName    = '../ConfFile_cfg.py'
config.JobType.outputFiles = ['DM.root']
config.JobType.allowUndistributedCMSSW = True

config.section_('Data')
config.Data.inputDBS = 'phys03'
config.Data.inputDataset = '/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/decosa-TT_Phys14DR-PU20bx25_PHYS14_25_V1-v1_EDMNtuple_18May2015_v2-4fb54d969ee15e14416bf4390fbb44b0/USER'
config.Data.splitting    = 'FileBased'
config.Data.unitsPerJob  = 20
config.Data.outLFNDirBase = '/store/user/jgarciaf'

config.section_('Site')
config.Site.storageSite = 'T2_ES_IFCA'
