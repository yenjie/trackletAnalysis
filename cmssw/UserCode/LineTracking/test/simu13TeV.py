import FWCore.ParameterSet.Config as cms

process = cms.Process("LineTracking")

process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.StandardSequences.Digi_cff")
process.load("RecoLocalTracker.Configuration.RecoLocalTracker_cff")

###############################################################################
# Source
process.source = cms.Source("PoolSource",
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
    skipEvents = cms.untracked.uint32(0),
    fileNames = cms.untracked.vstring(
        'file:///tmp/sikler/gen_sim_raw.root'
    ),
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

###############################################################################
# Paths
process.lreco = cms.Path(process.RawToDigi
                       * process.siStripZeroSuppression
                       * process.siStripClusters
                       * process.siStripMatchedRecHits)

process.load("UserCode.LineTracking.LineTracking_cfi")
process.greco = cms.Path(process.lineTracking)

process.out = cms.OutputModule("PoolOutputModule",
                fileName = cms.untracked.string("lineTrackingSimu.root")
        )
process.outp = cms.EndPath(process.out)

###############################################################################
# Global tag
process.GlobalTag.globaltag = 'MCRUN2_74_V8::All'

###############################################################################
# Schedule
process.schedule = cms.Schedule(process.lreco,
                                process.greco,
                                process.outp)

from SLHCUpgradeSimulations.Configuration.postLS1Customs import customisePostLS1
process = customisePostLS1(process)

