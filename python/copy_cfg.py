import FWCore.ParameterSet.Config as cms

# Give the process a name
process = cms.Process("PickEvent")

# Tell the process which files to use as the sourdce
process.source = cms.Source ("PoolSource",
          fileNames = cms.untracked.vstring ("/store/data/Run2011A/MuOnia/AOD/May10ReReco-v1/0000/002443C7-6B7E-E011-9B01-0018F3D09644.root")
                             
)

# tell the process to only run over 500 events (-1 would mean run over
#  everything
process.maxEvents = cms.untracked.PSet(
            input = cms.untracked.int32 (500)

)

# Tell the process what filename to use to save the output
process.Out = cms.OutputModule("PoolOutputModule",

         fileName = cms.untracked.string ("/afs/cern.ch/work/n/nsahoo/BPH-ANALYSIS/afb/SE/aod/data/2011_test_data_Run2011A_May10ReReco_500.root")             

)

# make sure everything is hooked up
process.end = cms.EndPath(process.Out)
