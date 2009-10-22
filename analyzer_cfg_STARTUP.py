import FWCore.ParameterSet.Config as cms

process = cms.Process("MuonSelector")

process.source = cms.Source("PoolSource",
#                            fileNames = cms.untracked.vstring()
    fileNames = cms.untracked.vstring('file:~/RECO_StZp_all_LowLum.root')
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.MuonAnalyzer = cms.EDAnalyzer('Analyzer',
                PrimaryVertices = cms.untracked.vstring("offlinePrimaryVertices","offlinePrimaryVerticesWithBS"),
                TriggerResultsTag = cms.untracked.InputTag("TriggerResults::HLT8E29"),
                triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD","","HLT8E29"),
                hltFilterNames = cms.VInputTag(cms.InputTag("hltL2Mu11L2Filtered11","","HLT8E29"),
                                               cms.InputTag("hltSingleMuIsoL3IsoFiltered3","","HLT8E29")
                                               )
					)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('Muons.root')
)
#HLT
#process.load("FWCore.MessageService.MessageLogger_cfi")
process.hltTrigReport = cms.EDAnalyzer( "HLTrigReport",
                                        HLTriggerResults = cms.InputTag( 'TriggerResults','','HLT')
                                       )
#process.MessageLogger.categories.append('HLTrigReport')
#process.MessageLogger.categories.append('MuonSelector')

process.p = cms.Path(process.MuonAnalyzer + process.hltTrigReport)
