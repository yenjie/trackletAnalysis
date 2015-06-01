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
       'file:///tmp/sikler/1C871E00-3301-E511-95EE-02163E01453D.root'
    ),
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

###############################################################################
# Paths
process.load("UserCode.LineTracking.LineTracking_cfi")

process.greco = cms.Path(process.siStripMatchedRecHits
                       * process.lineTracking)

process.out = cms.OutputModule("PoolOutputModule",
                fileName = cms.untracked.string("lineTrackingData.root")
        )
process.outp = cms.EndPath(process.out)

###############################################################################
# Global tag
process.GlobalTag.globaltag = 'GR_P_V54::All'

###############################################################################
# Schedule
process.schedule = cms.Schedule(process.greco,
                                process.outp)

from SLHCUpgradeSimulations.Configuration.postLS1Customs import customisePostLS1
process = customisePostLS1(process)

