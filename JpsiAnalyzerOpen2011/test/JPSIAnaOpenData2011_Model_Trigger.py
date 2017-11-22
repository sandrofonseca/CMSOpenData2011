# For Data Trigger
import FWCore.ParameterSet.Config as cms
import FWCore.PythonUtilities.LumiList as LumiList

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    'root://eospublic.cern.ch//eos/opendata/cms/Run2011A/DoubleMu/AOD/12Oct2013-v1/20000/045CCED6-033F-E311-9E93-003048678F74.root'       
    )
)

goodJSON = 'Cert_160404-180252_7TeV_ReRecoNov08_Collisions11_JSON.txt'
myLumis = LumiList.LumiList(filename = goodJSON).getCMSSWString().split(',')
process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange()
process.source.lumisToProcess.extend(myLumis)

process.TFileService = cms.Service ('TFileService',
    fileName = cms.string (
    'histojpsi.root'
    )
)
 
process.demo = cms.EDAnalyzer('JpsiAnalyzerOpen2011',
		#trigger_on = cms.bool(True),
        verbose = cms.bool(False),
	triggerflag = cms.bool(False), # False = Data and True = MC		      
	# Trigger
	TriggerResultsTag = cms.untracked.InputTag("TriggerResults", "", "HLT"),
    TriggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD", "", "HLT"),
    PathName = cms.untracked.string("HLT_DoubleMu7_v"),
    # RECO Labels
    primaryVertexProducer = cms.InputTag("offlinePrimaryVertices"),
	recoMuonsLabel = cms.InputTag("muons"),   
	# RECO Configs
    minMuPt = cms.double(2.0),# in GeV
    maxMuEta = cms.double(2.4),
    minMuonLeadPt = cms.double(20.0),# in GeV
	minMuonTrailPt = cms.double(4.0), # in GeV
	minJPsiMass = cms.double(2.95),# in GeV
	maxJPsiMass = cms.double(3.25)# in GeV 
)


process.p = cms.Path(process.demo)
