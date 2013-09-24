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

    # particle properties 
    MuonMass = cms.untracked.double(0.10565837), 
    MuonMassErr = cms.untracked.double(0.10565837*1e-6), 
    KaonMass = cms.untracked.double(0.493677), 
    KaonMassErr = cms.untracked.double(1.6e-5),
    BuMass = cms.untracked.double(5.27925),

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
 
    MuMuMinPt = cms.untracked.double(6.9),      # [GeV/c]
    MuMuMinInvMass = cms.untracked.double(1.0), # [GeV/c2]
    MuMuMaxInvMass = cms.untracked.double(4.8), # [GeV/c2]

    MuMuMinVtxCl = cms.untracked.double(0.05), 
    MuMuMinLxySigmaBs = cms.untracked.double(3.0), 
    MuMuMaxDca = cms.untracked.double(0.5), # [cm]
    MuMuMinCosAlphaBs = cms.untracked.double(0.9),

    # pre-selection cuts 
    TrkMinPt = cms.untracked.double(0.6), #[GeV/c]
    TrkMaxDcaSigBs = cms.untracked.double(1.2), # hadron DCA/sigma w/respect to BS 
    TrkMaxR = cms.untracked.double(110.0), # [cm]
    TrkMaxZ = cms.untracked.double(280.0), # [cm]
   
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

