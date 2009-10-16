import FWCore.ParameterSet.Config as cms

process = cms.Process("MuonSelector")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load('Configuration/StandardSequences/Services_cff')
process.load('Configuration/StandardSequences/GeometryIdeal_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('HLTrigger/Configuration/HLT_1E31_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'MC_31X_V8::All'

process.source = cms.Source("PoolSource",
                           # fileNames = cms.untracked.vstring()
    fileNames = cms.untracked.vstring('file:~/RECO_StZp_mumu_HighLum.root')
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.MuonAnalyzer = cms.EDAnalyzer('Analyzer',
                triggerEventTag = cms.untracked.InputTag("TriggerResults","","HLT"),
                hltFilterNames = cms.VInputTag(cms.InputTag("hltRPCMuonNormaL1Filtered0","","HLT"))
					)

process.hltTrigReport = cms.EDAnalyzer( "HLTrigReport",
HLTriggerResults = cms.InputTag( 'TriggerResults','','HLT')
)

process.hltL1GtTrigReport = cms.EDAnalyzer( "L1GtTrigReport",
    UseL1GlobalTriggerRecord = cms.bool( False ),
    L1GtRecordInputTag = cms.InputTag( "hltGtDigis" ) ## gtDigis
)

process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('MuonSelector')
process.MessageLogger.categories.append('TriggerSummaryProducerAOD')
#process.MessageLogger.categories.append('L1GtTrigReport')
#process.MessageLogger.categories.append('HLTrigReport')

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('Muons.root')
)
