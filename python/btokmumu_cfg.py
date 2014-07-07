import FWCore.ParameterSet.Config as cms
from btokmumu_2012_cfi import process 


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
                               'file:/afs/cern.ch/work/n/nsahoo/BPH-ANALYSIS/afb/SE/aod/data/2012_test_data_Run2012B_Jan2013_1000.root' )
#                                'file:/afs/cern.ch/work/n/nsahoo/BPH-ANALYSIS/afb/SE/aod/data/2012_test_data_Run2012B_15000.root')
                            
    )

process.GlobalTag.globaltag = cms.string('START53_V7G::All')


# do trigger matching for muons
triggerProcessName = 'HLT'


process.cleanMuonTriggerMatchHLT0 = cms.EDProducer(
        # match by DeltaR only (best match by DeltaR)
        'PATTriggerMatcherDRLessByR',
            src                   = cms.InputTag('cleanPatMuons'),
            # default producer label as defined in
            # PhysicsTools/PatAlgos/python/triggerLayer1/triggerProducer_cfi.py
            matched               = cms.InputTag('patTrigger'),
            matchedCuts           = cms.string('path("HLT_DoubleMu3p5_LowMass_Displaced*",0,0)'),
            maxDeltaR             = cms.double(0.1),
            # only one match per trigger object
            resolveAmbiguities    = cms.bool(True),
            # take best match found per reco object (by DeltaR here, see above)
            resolveByMatchQuality = cms.bool(False))


from PhysicsTools.PatAlgos.tools.trigTools import *
switchOnTriggerMatchEmbedding(process, triggerMatchers = ['cleanMuonTriggerMatchHLT0'],
                               hltProcess = triggerProcessName, outputModule = '')

g_TriggerNames_LastFilterNames = [
        ('HLT_DoubleMu3p5_LowMass_Displaced',  'hltDisplacedmumuFilterDoubleMu3p5LowMass')
        ]




g_TriggerNames = [i[0] for i in g_TriggerNames_LastFilterNames]
g_LastFilterNames = [i[1] for i in g_TriggerNames_LastFilterNames]

process.ntuple.TriggerNames = cms.vstring(g_TriggerNames)
process.ntuple.LastFilterNames = cms.vstring(g_LastFilterNames)

