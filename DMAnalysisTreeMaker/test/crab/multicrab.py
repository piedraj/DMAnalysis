import os
from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferLogs = True
config.General.workArea     = 'crab_data'  # Make sure you set this parameter

config.section_('JobType')
config.JobType.pluginName  = 'Analysis'
config.JobType.psetName    = '../ConfFile_cfg.py'
config.JobType.outputFiles = ['DM.root']
config.JobType.allowUndistributedCMSSW = True

config.section_('Data')    
config.Data.inputDBS      = 'phys03'
config.Data.splitting     = 'FileBased'
config.Data.unitsPerJob   = 20
config.Data.outLFNDirBase = '/store/user/jgarciaf'

config.section_('Site')
config.Site.storageSite = 'T2_ES_IFCA'


import sys

if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand

    def submit(config):
        print " to do: ",config
        res = crabCommand('submit', config = config)

    ######### From now on this is what users should modify. It is the a-la-CRAB2 configuration part.
   
    print sys.argv
    if len(sys.argv) <= 1 :
       print "no arguments?"
       print "Usage: python multicrab.py test.py"
       exit()
       

    samples = {}
    SamplesFile = sys.argv[1]
    print " SamplesFile = ", SamplesFile
    
    if os.path.exists(SamplesFile):
       handle = open(SamplesFile,'r')
       exec(handle)
       handle.close()
                
    # samples to be analysed
                   
    for key, value in samples.iteritems():
        print key, ' -> ', value
        
        config.General.requestName = key
        config.Data.inputDataset = value[0]
#       config.JobType.pyCfgParams = list(pyCfgParams)
#       config.JobType.pyCfgParams.extend(value[1])
        submit(config)

