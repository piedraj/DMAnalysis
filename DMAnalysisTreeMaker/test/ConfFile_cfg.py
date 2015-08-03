import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))   # correr sobre los 1000 primeros sucesos

process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000   # report poc pantalla cada x sucesos

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("DM.root"),
                                   closeFileFast = cms.untracked.bool(True))

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring('file:/gpfs/csic_projects/cms/jgarciaf/B2G/B2GEDMNtuple_10.root')   # muestra de juguete
    )

process.demo = cms.EDAnalyzer('DMAnalysisTreeMaker',
			      readGen = cms.untracked.bool(False)   # data vs mc
)

process.p = cms.Path(process.demo)
