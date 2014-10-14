import FWCore.ParameterSet.Config as cms
from btokmumu_2012_cfi import process

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
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
  'file:/afs/cern.ch/user/g/gechen/gechen/work/BToKMuMu/Ana/RootFiles/MC_GENOnly/BToKMuMu_MC_OnlyGEN_8TeV_9.root',
#    '/store/mc/Summer12_DR53X/BuToMuMuK_EtaPtFilter_8TeV-pythia6-evtgen/AODSIM/PU_S10_START53_V7C-v1/20000/9C4AFE3E-BA65-E211-821F-AC162DACC3E8.root',
    
                            )
                        )
#process.GlobalTag.globaltag = cms.string('GR_R_42_V25::All')r
#process.GlobalTag.globaltag = cms.string('START53_V23::All')
process.GlobalTag.globaltag = cms.string('START53_V7G::All')

process.ntuple.IsMonteCarlo = cms.untracked.bool(True)
process.ntuple.KeepGENOnly  = cms.untracked.bool(True)
process.p = cms.Path(process.ntuple)
