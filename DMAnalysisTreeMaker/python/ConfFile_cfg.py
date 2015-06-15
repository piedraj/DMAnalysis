import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("test.root"),
                                   closeFileFast = cms.untracked.bool(True))

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring('file:/gpfs/csic_projects/cms/jgarciaf/B2G/B2GEDMNtuple_1.root')
    )

process.demo = cms.EDAnalyzer('DMAnalysisTreeMaker')

process.p = cms.Path(process.demo)
