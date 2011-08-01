import FWCore.ParameterSet.Config as cms

import sys
if 'castor' in sys.argv[2]:
    input_filename="rfio:"+sys.argv[2]
else:
    input_filename="file:"+sys.argv[2]
output_filename=input_filename.split('/')
output_filename=output_filename[-1]
output_filename=output_filename.split('_')
output_filename="ntuple_"+output_filename[0]+'_'+output_filename[1]+'_'+output_filename[2]+'.root'
output_filename="/tmp/zhangjin/"+output_filename

print input_filename,"-> Ntuple:",output_filename
        
process = cms.Process("MuonSelector")
HLTProc = "HLT"

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.GeometryDB_cff")
#process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")
#process.TrackAssociatorByHits.Cut_RecoToSim = cms.double(0.5)
#process.load("SimTracker.TrackAssociation.TrackAssociatorByChi2_cfi")
#process.TrackAssociatorByChi2ESProducer.chi2cut = cms.double(1e9)
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.PyReleaseValidation.autoCond import autoCond
process.GlobalTag.globaltag = autoCond['mc']

process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(input_filename)
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.MuonAnalyzer = cms.EDAnalyzer('MuESelector',
                                    FileName = cms.string(output_filename),
                                    MuonPtCut = cms.untracked.double(15.),
                                    DiMuonInvarMassCut = cms.untracked.double(40.),
                                    TriggerResultsTag = cms.InputTag('TriggerResults','',HLTProc),
                                    triggerEventTag = cms.untracked.InputTag('hltTriggerSummaryAOD','',HLTProc),
                                    HLTObj_HLTNames = cms.vstring("HLT_Mu24_v100"),
                                    hltFilterNames = cms.VInputTag(cms.InputTag('hltSingleMu24L3Filtered24','',HLTProc)),
                                    StandardMuonCuts = cms.untracked.vstring("GlobalMuonPromptTight","TMLastStationLoose","TMLastStationTight","TMLastStationAngLoose","TMLastStationAngTight"),
                                    maxChamberDist = cms.untracked.double(-3.),
                                    maxChamberDistPull = cms.untracked.double(-3.),
                                    minTrackHits = cms.untracked.uint32(6),
                                    )
########### Event cleaning and  Trigger selection ###########
# Select events based on the HLTtriggers....singleJet and BTag triggers
# Use the instructions provided at:
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/TriggerResultsFilter
# This eases the trigger selection for different HLT menus and also takes care of wildcard and trigger versioning
#############################################################
process.MuonAnalyzer.VertexFilterPSet = cms.PSet ( vertexCollection = cms.InputTag("offlinePrimaryVertices"),
                                                   minimumNDOF = cms.uint32(4),
                                                   maxAbsZ = cms.double(24),
                                                   maxd0 = cms.double(2.0)
                                                 )

process.MuonAnalyzer.noscrapingPSet = cms.PSet( numtrack = cms.untracked.uint32(10),
                                thresh = cms.untracked.double(0.2)
                                )

process.MuonAnalyzer.HLTFilterPSet = cms.PSet ( triggerConditions = cms.vstring(
                                                #SingleMuStream Triggers 
                                                "HLT_IsoMu*",
                                                "HLT_L1SingleMu*",
                                                "HLT_L1DoubleMu*",
                                                "HLT_L2Mu*",
                                                "HLT_Mu??_v*",
                                                "HLT_Mu?_v*",
                                                #DoubleMuStream Triggers
                                                "HLT_DoubleMu*",
                                                "HLT_L1DoubleMu*",
                                                "HLT_L2DoubleMu*",
                                                "HLT_Mu*_Jet*_v*",
                                                "HLT_TripleMu*"
                                                ),
                                              hltResults = cms.InputTag('TriggerResults','',HLTProc),
                                              l1tResults = cms.InputTag( "" ),
                                              l1tIgnoreMask = cms.bool(False),
                                              l1techIgnorePrescales   = cms.bool(False),
                                              daqPartitions           = cms.uint32(0x01),
                                              throw = cms.bool( False ) #set to false to deal with missing triggers while running over different trigger menus
                                            )


#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",moduleMemorySummary = cms.untracked.bool(True))
process.p = cms.Path(process.MuonAnalyzer)

########### Messager Logger ###########
process.load('FWCore.MessageService.MessageLogger_cfi')
#process.MessageLogger.categories.extend( cms.untracked.vstring('DataLost','DataFormat','RefNULL','InvalidTag','CodeWrong') )
process.MessageLogger.cerr.FwkReport.reportEvery = 500 #quench the message logger (optional)
process.MessageLogger.cerr_stats.threshold = cms.untracked.string('WARNING')
#process.MessageLogger.destinations.append('MSGOUT')
#process.MessageLogger.MSGOUT=cms.untracked.PSet(threshold = cms.untracked.string('WARNING'))
#process.MessageLogger.statistics.append('MSG_stats'),
#process.MessageLogger.MSG_stats = cms.untracked.PSet(
#    optionalPSet = cms.untracked.bool(True),
#    threshold = cms.untracked.string('WARNING'),
#    output = cms.untracked.string('MSGOUT')
#)
#process.hltTrigReport = cms.EDAnalyzer( "HLTrigReport",HLTriggerResults = cms.InputTag( 'TriggerResults','',HLTProc))
#process.MessageLogger.categories.append('HLTrigReport')


process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

#grep Events *STDOUT| awk '{print $5}' to print the total events number processed
