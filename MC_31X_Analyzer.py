import FWCore.ParameterSet.Config as cms

process = cms.Process("MuonSelector")

process.source = cms.Source("PoolSource",
#                            fileNames = cms.untracked.vstring()
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
    fileNames = cms.untracked.vstring("rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut200-MC_31X_V3/PYTHIA6_DYmumu_Mcut200_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_1.root",
                                      "rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut200-MC_31X_V3/PYTHIA6_DYmumu_Mcut200_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_2.root",
                                      "rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut200-MC_31X_V3/PYTHIA6_DYmumu_Mcut200_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_3.root",
                                      "rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut200-MC_31X_V3/PYTHIA6_DYmumu_Mcut200_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_4.root",
                                      "rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut200-MC_31X_V3/PYTHIA6_DYmumu_Mcut200_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_5.root",
                                      "rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut200-MC_31X_V3/PYTHIA6_DYmumu_Mcut200_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_6.root",
                                      "rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut200-MC_31X_V3/PYTHIA6_DYmumu_Mcut200_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_7.root",
                                      "rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut200-MC_31X_V3/PYTHIA6_DYmumu_Mcut200_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_8.root",
                                      "rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut200-MC_31X_V3/PYTHIA6_DYmumu_Mcut200_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_9.root",
                                      "rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut200-MC_31X_V3/PYTHIA6_DYmumu_Mcut200_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_10.root",
                                      "rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut200-MC_31X_V3/PYTHIA6_DYmumu_Mcut200_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_11.root",
                                      "rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut200-MC_31X_V3/PYTHIA6_DYmumu_Mcut200_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_12.root",
                                      "rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut200-MC_31X_V3/PYTHIA6_DYmumu_Mcut200_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_13.root",
                                      "rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut200-MC_31X_V3/PYTHIA6_DYmumu_Mcut200_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_14.root",
                                      "rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut200-MC_31X_V3/PYTHIA6_DYmumu_Mcut200_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_15.root",
                                      "rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut200-MC_31X_V3/PYTHIA6_DYmumu_Mcut200_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_16.root",
                                      "rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut200-MC_31X_V3/PYTHIA6_DYmumu_Mcut200_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_17.root",
                                      "rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut200-MC_31X_V3/PYTHIA6_DYmumu_Mcut200_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_18.root",
                                      "rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut200-MC_31X_V3/PYTHIA6_DYmumu_Mcut200_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_19.root",
                                      "rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut200-MC_31X_V3/PYTHIA6_DYmumu_Mcut200_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_20.root",
                                      "rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut200-MC_31X_V3/PYTHIA6_DYmumu_Mcut200_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_21.root",
                                      "rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut200-MC_31X_V3/PYTHIA6_DYmumu_Mcut200_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_22.root",
                                      "rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut200-MC_31X_V3/PYTHIA6_DYmumu_Mcut200_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_23.root",
                                      "rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut200-MC_31X_V3/PYTHIA6_DYmumu_Mcut200_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_24.root",
                                      "rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut200-MC_31X_V3/PYTHIA6_DYmumu_Mcut200_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_25.root")
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.MuonAnalyzer = cms.EDAnalyzer('Analyzer',
                PrimaryVertices = cms.untracked.vstring("offlinePrimaryVertices","offlinePrimaryVerticesWithBS"),
                TriggerResultsTag = cms.untracked.InputTag("TriggerResults::HLT"),
                triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD","","HLT"),
                hltFilterNames = cms.VInputTag(cms.InputTag("hltSingleMu15L3PreFiltered15","","HLT"),
                                               cms.InputTag("hltSingleMuIsoL3IsoFiltered9","","HLT"),
                                               )
					)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('~/scratch0/MyNtuple_DYmumu_PYTHIA6_Mcut200_10TeV_HLT1E31_MC31X_V3.root')
)
#HLT
#process.load("FWCore.MessageService.MessageLogger_cfi")
process.hltTrigReport = cms.EDAnalyzer( "HLTrigReport",
                                        HLTriggerResults = cms.InputTag( 'TriggerResults','','HLT')
                                       )
#process.MessageLogger.categories.append('HLTrigReport')
#process.MessageLogger.categories.append('MuonSelector')

process.p = cms.Path(process.MuonAnalyzer + process.hltTrigReport)
