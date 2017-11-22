import FWCore.ParameterSet.Config as cms

process = cms.Process("TestTriggerReport")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['TriggerReport.txt']
process.MessageLogger.categories.append('HLTrigReport')
process.MessageLogger.categories.append('L1GtTrigReport')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace this file with the one you actually want to use
    fileNames = cms.untracked.vstring(
        'root://eospublic.cern.ch//eos/opendata/cms/Run2011A/DoubleMu/AOD/12Oct2013-v1/20000/045CCED6-033F-E311-9E93-003048678F74.root'
    )
)

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:hltonline_8E33v2', '')

process.load("HLTrigger.HLTanalyzers.hlTrigReport_cfi")
process.hlTrigReport.HLTriggerResults = cms.InputTag( 'TriggerResults','','HLT') 

process.load("L1Trigger.GlobalTriggerAnalyzer.l1GtTrigReport_cfi")
process.l1GtTrigReport.L1GtRecordInputTag = cms.InputTag( 'gtDigis' )

process.HLTAnalyzerEndpath = cms.EndPath( process.l1GtTrigReport + process.hlTrigReport )
