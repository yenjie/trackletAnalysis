import FWCore.ParameterSet.Config as cms

process = cms.Process("myRECO")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.MagneticField_0T_cff")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("RecoVertex.BeamSpotProducer.BeamSpot_cfi")

# Timing service
process.Timing = cms.Service("Timing") 

# MC Globaltag for 2015 dN/deta analysis
process.GlobalTag.globaltag = 'MCRUN2_74_V8::All'

process.pixelVertexFromClusters = cms.EDProducer('PixelVertexProducerClusters')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#    'file:RelVal_MinBias_13TeV_28478DD9-99A9-E411-891C-0025905B861C.root'
#'file:step3_RAW2DIGI_L1Reco_RECO.root'
'file:step3_EPOS_100k.root'

    #'RelValMinBias_314_STARTUP31X_V2-v1-Reco.root'
    )
)

# Centrality
process.load("RecoHI.HiCentralityAlgos.pACentrality_cfi") 
process.pACentrality.producePixelTracks = cms.bool(False)
#process.pACentrality.produceTracks = cms.bool(False)


# Centrality Binning
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi") 
process.centralityBin.Centrality = cms.InputTag("pACentrality")
process.centralityBin.centralityVariable = cms.string("HFtowers")
process.centralityBin.nonDefaultGlauberModel = cms.string("HydjetDrum5")

# Add the HeavyIon Record: it is for PbPb cent binning, so we shoud not
# trust the centrality bin and only use the variables from the centrality
# provider

process.GlobalTag.toGet.extend([
   cms.PSet(record = cms.string("HeavyIonRcd"),
   tag = cms.string("CentralityTable_HFtowers200_HydjetDrum5_v740x01_mc"),
   connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_PHYSICSTOOLS"),
   label = cms.untracked.string("HFtowersHydjetDrum5")
   ),
])

process.ana = cms.EDAnalyzer('PixelHitAnalyzer',
                             vertexSrc = cms.vstring('pixelVertexFromClusters'),
                             trackSrc = cms.untracked.InputTag('generalTracks'),
                             doTracking = cms.untracked.bool(False),
                             doCentrality = cms.untracked.bool(True)
                             )

process.anaStrip = cms.EDAnalyzer('StripHitAnalyzer',
                             vertexSrc = cms.vstring('pixelVertexFromClusters'),
                             trackSrc = cms.untracked.InputTag('generalTracks'),
                             doTracking = cms.untracked.bool(False),
                             doCentrality = cms.untracked.bool(True),
     RecHitCollections = cms.VInputTag(
     # cms.InputTag('siStripMatchedRecHits','rphiRecHit'), 
                                             cms.InputTag('siStripMatchedRecHits','matchedRecHit')
                             )
                             )

#process.SiStripRecHitsAnalyzer = cms.EDAnalyzer('SiStripRecHitsAnalyzer',
#     RecHitCollections = cms.VInputTag( cms.InputTag('siStripMatchedRecHits','rphiRecHit'), 
#                                             cms.InputTag('siStripMatchedRecHits','stereoRecHit')
#                             )
#)                                                                                       
#process.load("HeavyIonsAnalysis.EventAnalysis.hievtanalyzer_mc_cfi")
process.load("HLTrigger.HLTanalyzers.HLTBitAnalyser_cfi")

process.hltbitanalysis.UseTFileService = cms.untracked.bool(True)
process.hltanalysis = process.hltbitanalysis.clone(
    l1GtReadoutRecord    = cms.InputTag("gtDigis"),
    l1GctHFBitCounts     = cms.InputTag("gctDigis"),
    l1GctHFRingSums      = cms.InputTag("gctDigis"),
    l1extramu            = cms.string('l1extraParticles'),
    l1extramc            = cms.string('l1extraParticles'),
    hltresults           = cms.InputTag("TriggerResults","","HLT"),
)

process.TFileService = cms.Service('TFileService',
                                   fileName = cms.string('PixelTree-EPOS.root')
                                   )


process.analyze = cms.Path(
       process.siPixelRecHits*
       process.siStripMatchedRecHits*
      process.pixelVertexFromClusters*
#       process.hiSelectedVertex*
       process.pACentrality*
       process.centralityBin*
       process.hltanalysis*
#       process.hiEvtAnalyzer*
       process.ana*
       process.anaStrip
#       process.SiStripRecHitsAnalyzer
      )
     
