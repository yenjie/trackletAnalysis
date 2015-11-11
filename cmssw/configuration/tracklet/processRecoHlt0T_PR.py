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
process.GlobalTag.globaltag = 'GR_P_V54::All'

process.pixelVertexFromClusters = cms.EDProducer('PixelVertexProducerClusters')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#    'file:RelVal_MinBias_13TeV_28478DD9-99A9-E411-891C-0025905B861C.root'
#'file:step3_RAW2DIGI_L1Reco_RECO.root'
'file:promptReco/0880DA14-3201-E511-8164-02163E014458.root',
'file:promptReco/1C871E00-3301-E511-95EE-02163E01453D.root',
'file:promptReco/240A958B-3201-E511-B276-02163E014686.root',
'file:promptReco/26D30546-3201-E511-8338-02163E01432C.root',
'file:promptReco/2EE8A2AB-3A01-E511-9F37-02163E012324.root',
'file:promptReco/3252E700-3301-E511-A790-02163E014614.root',
'file:promptReco/3647A68F-3201-E511-9400-02163E012A0D.root',
'file:promptReco/38DE00BB-3101-E511-8365-02163E0141CE.root',
'file:promptReco/42C54FCF-3201-E511-B436-02163E014458.root',
'file:promptReco/4623CD4A-3201-E511-8CD4-02163E012462.root',
'file:promptReco/48E81A34-3401-E511-8B6C-02163E0141F9.root',
'file:promptReco/4AFBD421-3B01-E511-ABEC-02163E011DDC.root',
'file:promptReco/52A184AE-3701-E511-92A0-02163E0141D2.root',
'file:promptReco/54027BDB-3301-E511-AF9D-02163E013476.root',
'file:promptReco/549E845E-3301-E511-B6B3-02163E014668.root',
'file:promptReco/5662C54A-3101-E511-93CD-02163E013476.root',
'file:promptReco/569C7E47-3201-E511-8ACC-02163E0138BB.root',
'file:promptReco/5E4E1C4F-3201-E511-82C3-02163E014613.root',
'file:promptReco/66828D38-3201-E511-A2C8-02163E0122E1.root',

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
                             doMC = cms.untracked.bool(False),
                             doCentrality = cms.untracked.bool(True)
                             )

process.anaStrip = cms.EDAnalyzer('StripHitAnalyzer',
                             vertexSrc = cms.vstring('pixelVertexFromClusters'),
                             trackSrc = cms.untracked.InputTag('generalTracks'),
                             doTracking = cms.untracked.bool(False),
                             doMC = cms.untracked.bool(False),
                             doCentrality = cms.untracked.bool(True),
     RecHitCollections = cms.VInputTag( 
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
                                   fileName = cms.string('PixelTree_PR.root')
                                   )


import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltHighPtPA = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltHighPtPA.HLTPaths = ['HLT_ZeroBias_part0_v1',
'HLT_ZeroBias_part1_v1',
'HLT_ZeroBias_part2_v1',
'HLT_ZeroBias_part3_v1',
'HLT_ZeroBias_part4_v1',
'HLT_ZeroBias_part5_v1',
'HLT_ZeroBias_part6_v1',
'HLT_ZeroBias_part7_v1'
                                          ]
process.hltHighPtPA.throw = False
process.hltHighPtPA.andOr = True

                                                                                                                                                  #

process.analyze = cms.Path(
       process.hltHighPtPA*
       process.siPixelRecHits*
       process.siStripMatchedRecHits*
      process.pixelVertexFromClusters*
#       process.hiSelectedVertex*
       process.pACentrality*
       process.centralityBin*
       process.hltanalysis*
#       process.hiEvtAnalyzer*
#       process.ana*
       process.anaStrip
#       process.SiStripRecHitsAnalyzer
      )
     

