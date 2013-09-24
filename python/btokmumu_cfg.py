import FWCore.ParameterSet.Config as cms
from btokmumu_cfi import process 


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring('file:/afs/cern.ch/user/x/xshi/work/cms/afb/dat/aod/data/Run2011/2011May10ReReco_numEvent1000.root' )
    )

process.GlobalTag.globaltag = cms.string('GR_R_42_V25::All')


# do trigger matching for muons
triggerProcessName = 'HLT'

process.cleanMuonTriggerMatchHLT0 = cms.EDProducer(
    # match by DeltaR only (best match by DeltaR)
    'PATTriggerMatcherDRLessByR',                         
    src                   = cms.InputTag('cleanPatMuons'),
    # default producer label as defined in
    # PhysicsTools/PatAlgos/python/triggerLayer1/triggerProducer_cfi.py
    matched               = cms.InputTag('patTrigger'),
    matchedCuts           = cms.string('path("HLT_Dimuon6p5_LowMass_Displaced_v*",0,0)'),
    maxDeltaR             = cms.double(0.1),
    # only one match per trigger object
    resolveAmbiguities    = cms.bool(True),
    # take best match found per reco object (by DeltaR here, see above)       
    resolveByMatchQuality = cms.bool(False))

process.cleanMuonTriggerMatchHLT1 = cms.EDProducer(
    'PATTriggerMatcherDRLessByR',
    src                   = cms.InputTag('cleanPatMuons'),
    matched               = cms.InputTag('patTrigger'),
    matchedCuts           = cms.string('path("HLT_Dimuon7_LowMass_Displaced_v*")'),
    maxDeltaR             = cms.double(0.1),
    resolveAmbiguities    = cms.bool(True),
    resolveByMatchQuality = cms.bool(False))
process.cleanMuonTriggerMatchHLT2 = cms.EDProducer(
    'PATTriggerMatcherDRLessByR',
    src                   = cms.InputTag('cleanPatMuons'),
    matched               = cms.InputTag('patTrigger'),
    matchedCuts           = cms.string('path("HLT_DoubleMu4_LowMass_Displaced_v*")'),
    maxDeltaR             = cms.double(0.1),
    resolveAmbiguities    = cms.bool(True),
    resolveByMatchQuality = cms.bool(False))
process.cleanMuonTriggerMatchHLT3 = cms.EDProducer(
    'PATTriggerMatcherDRLessByR',
    src                   = cms.InputTag('cleanPatMuons'),
    matched               = cms.InputTag('patTrigger'),
    matchedCuts           = cms.string('path("HLT_DoubleMu4p5_LowMass_Displaced_v*")'),
    maxDeltaR             = cms.double(0.1),
    resolveAmbiguities    = cms.bool(True),
    resolveByMatchQuality = cms.bool(False))
process.cleanMuonTriggerMatchHLT4 = cms.EDProducer(
    'PATTriggerMatcherDRLessByR',
    src                   = cms.InputTag('cleanPatMuons'),
    matched               = cms.InputTag('patTrigger'),
    matchedCuts           = cms.string('path("HLT_DoubleMu5_LowMass_Displaced_v*")'),
    maxDeltaR             = cms.double(0.1),
    resolveAmbiguities    = cms.bool(True),
    resolveByMatchQuality = cms.bool(False))

from PhysicsTools.PatAlgos.tools.trigTools import *
# switchOnTriggerMatchEmbedding(process, triggerMatchers = ['cleanMuonTriggerMatchHLT0'],
#                               hltProcess = triggerProcessName, outputModule = '')

switchOnTriggerMatchEmbedding(process, triggerMatchers = ['cleanMuonTriggerMatchHLT1','cleanMuonTriggerMatchHLT2','cleanMuonTriggerMatchHLT3','cleanMuonTriggerMatchHLT4'], hltProcess = triggerProcessName, outputModule = '')

g_TriggerNames_LastFilterNames = [
    ('HLT_Dimuon6p5_LowMass_Displaced',
     'hltDisplacedmumuFilterLowMass'), #5E32

    ('HLT_Dimuon7_LowMass_Displaced',
      'hltDisplacedmumuFilterLowMass'), #1E33, 1.4E33
      
    ('HLT_DoubleMu4_LowMass_Displaced',
     'hltDisplacedmumuFilterLowMass'), #2E33

    ('HLT_DoubleMu4p5_LowMass_Displaced',
     'hltDisplacedmumuFilterDoubleMu4p5LowMass'), #3E33, 5E33
     
    ('HLT_DoubleMu5_LowMass_Displaced',
     'hltDisplacedmumuFilterDoubleMu5LowMass'), #3E33, 5E33
    ]

g_TriggerNames = [i[0] for i in g_TriggerNames_LastFilterNames]
g_LastFilterNames = [i[1] for i in g_TriggerNames_LastFilterNames]

process.ntuple.TriggerNames = cms.vstring(g_TriggerNames)
process.ntuple.LastFilterNames = cms.vstring(g_LastFilterNames)

