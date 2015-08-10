import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))   # correr sobre los 100 primeros sucesos

process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 10000   # report poc pantalla cada x sucesos

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("DM.root"),
                                   closeFileFast = cms.untracked.bool(True))

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring('file:/gpfs/csic_projects/cms/jgarciaf/B2G/B2GEDMNtuple_10.root')   # muestra de juguete
### fileNames = cms.untracked.vstring('file:/afs/cern.ch/user/p/piedra/public/forB2G/B2GEDMNtuple_1.root')   # muestra de juguete
    )

process.demo = cms.EDAnalyzer('DMAnalysisTreeMaker',
			      readGen = cms.untracked.bool(False)   # data vs mc
)

process.demo = cms.EDAnalyzer('DMAnalysisTreeMaker',
			      weight = cms.untracked.bool(False)   # negative weights vs normal weights
)

process.p = cms.Path(process.demo)
