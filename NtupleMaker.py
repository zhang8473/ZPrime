import FWCore.ParameterSet.Config as cms
process = cms.Process("MuonSelector")

#process.load("Configuration.StandardSequences.MagneticField_cff")
#process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")
#process.TrackAssociatorByHits.Cut_RecoToSim = cms.double(0.5)
#process.load("SimTracker.TrackAssociation.TrackAssociatorByChi2_cfi")
#process.TrackAssociatorByChi2ESProducer.chi2cut = cms.double(1e9)

#process.options = cms.untracked.PSet(fileMode = cms.untracked.string('NOMERGE'))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring("rfio:/castor/cern.ch/user/z/zhangjin/CMSSW_3_8_4/StuZp_MuMu_M600_Epsilon06_7TeV_n4000_START38V10_RECO.root")
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.MuonAnalyzer = cms.EDFilter('DiMuonSelector',
                                    FileName = cms.string("StuZpMuMu.root"),
                                    CrossSection = cms.double(1.4010E+003),
                                    MuonPtSelection = cms.untracked.double(15.),
                                    MuonHighestPtSelection = cms.untracked.double(40.),
                                    TriggerResultsTag = cms.untracked.InputTag('TriggerResults','','HLT'),
                                    triggerEventTag = cms.untracked.InputTag('hltTriggerSummaryAOD','','HLT'),
                                    hltFilterNames = cms.VInputTag(cms.InputTag('hltSingleMu9L3Filtered9','','HLT')),
                                    StandardMuonCuts = cms.untracked.vstring("GlobalMuonPromptTight","TMLastStationLoose","TMLastStationTight","TMLastStationAngLoose","TMLastStationAngTight"),
                                    maxChamberDist = cms.untracked.double(-3.),
                                    maxChamberDistPull = cms.untracked.double(-3.),
                                    minTrackHits = cms.untracked.uint32(3)
                                    )

#HLT
#process.load("FWCore.MessageService.MessageLogger_cfi")
#process.hltTrigReport = cms.EDAnalyzer( "HLTrigReport",HLTriggerResults = cms.InputTag( 'TriggerResults','','HLT'))
#process.MessageLogger.categories.append('HLTrigReport')
#process.MessageLogger.categories.append('MuonSelector')
process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",moduleMemorySummary = cms.untracked.bool(True))
process.p = cms.Path(process.MuonAnalyzer)
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 500 #quench the message logger (optional)
