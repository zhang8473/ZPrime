import FWCore.ParameterSet.Config as cms

process = cms.Process("MuonSelector")

process.source = cms.Source("PoolSource",
#                            fileNames = cms.untracked.vstring()
    fileNames = cms.untracked.vstring('file:~/Test_CMSSW333_DYmumu_Mcut500_7TeV_n10_MC31XV9_RECO.root')
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.MuonAnalyzer = cms.EDAnalyzer('Analyzer',
                FileName = cms.string("Muons.root"),
                PrimaryVertices = cms.untracked.vstring("offlinePrimaryVertices","offlinePrimaryVerticesWithBS"),
                TriggerResultsTag = cms.untracked.InputTag("TriggerResults::HLT"),
                triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD","","HLT"),
                hltFilterNames = cms.VInputTag(cms.InputTag("hltSingleMu15L3PreFiltered15","","HLT"),
                                               cms.InputTag("hltSingleMuIsoL3IsoFiltered9","","HLT"),
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
