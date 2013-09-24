import FWCore.ParameterSet.Config as cms

process = cms.Process("Ntuple")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.MessageLogger = cms.Service(
    "MessageLogger",
    destinations = cms.untracked.vstring('cerr', 'cout', 'message'),
    categories = cms.untracked.vstring('myBu'),
    cerr = cms.untracked.PSet(threshold = cms.untracked.string('WARNING')),
    cout = cms.untracked.PSet(
        threshold = cms.untracked.string('INFO'),
        INFO = cms.untracked.PSet(limit = cms.untracked.int32(0)), 
        myBu = cms.untracked.PSet(limit = cms.untracked.int32(-1)), 
    ), 
     message = cms.untracked.PSet(
         threshold = cms.untracked.string('INFO'),
         INFO = cms.untracked.PSet(limit = cms.untracked.int32(0)), 
         myBu = cms.untracked.PSet(limit = cms.untracked.int32(-1)), 
     )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.load("PhysicsTools.PatAlgos.patSequences_cff")

from PhysicsTools.PatAlgos.tools.coreTools import *
removeMCMatching(process, ['All'], outputInProcess = False)

process.ntuple = cms.EDAnalyzer(
    'BToKMuMu',
    OutputFileName = cms.string("BuToKMuMu.root"),
    BuildBuToKMuMu = cms.untracked.bool(True), 

    # labels
    TriggerResultsLabel = cms.InputTag("TriggerResults","", 'HLT'),
    BeamSpotLabel = cms.InputTag('offlineBeamSpot'),
    VertexLabel = cms.InputTag('offlinePrimaryVertices'),
    MuonLabel = cms.InputTag('cleanPatMuonsTriggerMatch'),
    TrackLabel = cms.InputTag('cleanPatTrackCands'), 
    TriggerNames = cms.vstring([]),
    LastFilterNames = cms.vstring([]),

    # HLT-trigger cuts 
    MuonMinPt = cms.untracked.double(3.0), # [GeV]
    MuonMaxEta = cms.untracked.double(2.2),  
    MuonMaxDcaBs = cms.untracked.double(2.0), # [cm]
    
)

# Remove not used from PAT 
process.patDefaultSequence.remove(process.patJetCorrFactors)
process.patDefaultSequence.remove(process.patJetCharge)
process.patDefaultSequence.remove(process.patJetPartonMatch)
process.patDefaultSequence.remove(process.patJetGenJetMatch)
process.patDefaultSequence.remove(process.patJetPartons)
process.patDefaultSequence.remove(process.patJetPartonAssociation)
process.patDefaultSequence.remove(process.patJetFlavourAssociation)
process.patDefaultSequence.remove(process.patJets)
process.patDefaultSequence.remove(process.metJESCorAK5CaloJet)
process.patDefaultSequence.remove(process.metJESCorAK5CaloJetMuons)
process.patDefaultSequence.remove(process.patMETs)
process.patDefaultSequence.remove(process.selectedPatJets)
process.patDefaultSequence.remove(process.cleanPatJets)
process.patDefaultSequence.remove(process.countPatJets)

process.p = cms.Path(process.patDefaultSequence * process.ntuple)

