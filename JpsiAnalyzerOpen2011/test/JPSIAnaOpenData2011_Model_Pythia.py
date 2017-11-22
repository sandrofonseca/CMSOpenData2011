import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    'root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2011/Summer11LegDR/JPsiToMuMu_2MuPEtaFilter_7TeV-pythia6-evtgen-v2/AODSIM/PU_S13_START53_LV6-v1/00000/0062051B-8C33-E411-953B-00215AEDFC8E.root'


    )
)

process.TFileService = cms.Service ('TFileService',
    fileName = cms.string (
    'histojpsiMC.root'
    )
)

process.demo = cms.EDAnalyzer('JpsiAnalyzerOpen2011',
        verbose = cms.bool(False),
	triggerflag = cms.bool(True), # False = Data and True = MC
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
