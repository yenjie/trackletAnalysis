import FWCore.ParameterSet.Config as cms

process = cms.Process("Pixel")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("RecoVertex.BeamSpotProducer.BeamSpot_cfi")

# Timing service
#process.Timing = cms.Service("Timing")

# MC Globaltag for 2016 dN/deta analysis
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_Express_v15', '')
process.GlobalTag.snapshotTime = cms.string("9999-12-31 23:59:59.000")

process.GlobalTag.toGet.extend([
    cms.PSet(record = cms.string("HeavyIonRcd"),
        connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
        ## 5.02 TeV Centrality Tables
        #tag = cms.string("CentralityTable_HFtowers200_HydjetDrum5_v740x01_mc"),
        #label = cms.untracked.string("HFtowersHydjetDrum5")
        ## 2.76 TeV Centrality Tables for data
        tag = cms.string("CentralityTable_HFtowers200_Glauber2010A_eff99_run1v750x01_offline"),
        label = cms.untracked.string("HFtowers")
    ),
])

process.pixelVertexFromClusters = cms.EDProducer('PixelVertexProducerClusters')

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
                            fileNames = cms.untracked.vstring(
'___C___'
    )
)

process.ana = cms.EDAnalyzer('PixelHitAnalyzer',
                             vertexSrc = cms.vstring('offlinePrimaryVerticesWithBS'),
                             trackSrc = cms.untracked.InputTag('generalTracks'),
                             doTracking = cms.untracked.bool(False),
                             doMC = cms.untracked.bool(False),
                             )

process.load('HeavyIonsAnalysis.EventAnalysis.hltanalysis_cff')

process.TFileService = cms.Service('TFileService',
                                   fileName = cms.string('PixelTree-___D___.root'))

process.analyze = cms.Path(process.hltanalysis *
                           process.siPixelRecHits *
                           process.ana)
