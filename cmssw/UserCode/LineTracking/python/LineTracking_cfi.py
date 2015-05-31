import FWCore.ParameterSet.Config as cms

lineTracking = cms.EDProducer("LineTrackingProducer",
  maxClusters         = cms.int32(6000),
  maxClusterWidthDiff = cms.int32(2),
  nRounds             = cms.int32(2),
  maxFirstHitRadius   = cms.double(31.),
  maxSharedDets       = cms.int32(3),
  maxAbsoluteZ0       = cms.double(20.),
  dMaxForVertexing    = cms.double(50.)
)
