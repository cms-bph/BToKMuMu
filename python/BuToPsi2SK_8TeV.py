import FWCore.ParameterSet.Config as cms
from btokmumu_2012_cfi import process


#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
                               #'file:/afs/cern.ch/user/x/xshi/work/cms/afb/dat/aod/mc/BuToKstarMuMu_7TeV_PYTHIA6.root',
# 'file:/afs/cern.ch/user/x/xshi/work/cms/afb/dat/aod/mc/BuToKstarMuMu_7TeV_PYTHIA6_1_1_96m.root',
# 'file:/afs/cern.ch/user/x/xshi/work/cms/afb/dat/aod/mc/BuToKstarMuMu_7TeV_PYTHIA6_2_1_dfI.root',
# 'file:/afs/cern.ch/user/x/xshi/work/cms/afb/dat/aod/mc/BuToKstarMuMu_7TeV_PYTHIA6_3_1_ud5.root',
# 'file:/afs/cern.ch/user/x/xshi/work/cms/afb/dat/aod/mc/BuToKstarMuMu_7TeV_PYTHIA6_4_1_neu.root',
# 'file:/afs/cern.ch/user/x/xshi/work/cms/afb/dat/aod/mc/BuToKstarMuMu_7TeV_PYTHIA6_5_1_8Kc.root',
#                                'file:/afs/cern.ch/user/x/xshi/work/cms/afb/dat/aod/mc/BuToKstarJPsi_7TeV_PYTHIA6_100.root',
#  'file:/afs/cern.ch/user/w/wguo/work/CMSSW_5_3_11/src/wguo_xin/BToKMuMu/python/2012testMC_numEvent5000.root',
#  'file:/afs/cern.ch/work/n/nsahoo/BPH-ANALYSIS/afb/SE/aod/mc/2012_test_BToKMuMu_500.root',
#  'file:/afs/cern.ch/work/n/nsahoo/BPH-ANALYSIS/afb/SE/aod/mc/2012_test_BToKMuMu_1000.root',
#'root://xrootd.unl.edu//store/mc/Summer12_DR53X/Bu2MuMuK_TuneZ2star_8TeV_Pythia6/AODSIM/PU_RD2_START53_V19F-v1/00000/000771A0-9408-E411-94C0-002590DB92A8.root',
'root://xrootd.unl.edu//store/mc/Summer12_DR53X/BuToPsiK_KFilter_TuneZ2star_SVS_8TeV-pythia6-evtgen/AODSIM/PU_RD2_START53_V19F-v1/00000/001E2794-0711-E411-B5F0-02163E008EBC.root',
#    '/store/mc/Summer12_DR53X/BuToMuMuK_EtaPtFilter_8TeV-pythia6-evtgen/AODSIM/PU_S10_START53_V7C-v1/20000/9C4AFE3E-BA65-E211-821F-AC162DACC3E8.root',
    
                            )
                        )
#process.GlobalTag.globaltag = cms.string('GR_R_42_V25::All')r
#process.GlobalTag.globaltag = cms.string('START53_V23::All')
#process.GlobalTag.globaltag = cms.string('START53_V7G::All')
process.GlobalTag.globaltag = cms.string('START53_V19F::All')
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
process.ntuple.IsMonteCarlo = cms.untracked.bool(True)
process.ntuple.KeepGENOnly  = cms.untracked.bool(False)
