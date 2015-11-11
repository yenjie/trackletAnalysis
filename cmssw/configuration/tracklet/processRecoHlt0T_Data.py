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
process.GlobalTag.globaltag = 'GR_E_V47::All'

process.pixelVertexFromClusters = cms.EDProducer('PixelVertexProducerClusters')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#    'file:RelVal_MinBias_13TeV_28478DD9-99A9-E411-891C-0025905B861C.root'
#'file:step3_RAW2DIGI_L1Reco_RECO.root'
#'file:Run246908/CE86B688-D509-E511-8C80-02163E014177.root'
#'/store/express/Run2015A/ExpressPhysics/FEVT/Express-v1/000/247/324/00000/A46E2596-100D-E511-AE7A-02163E01295D.root'
'/store/express/Run2015A/ExpressPhysics/FEVT/Express-v1/000/247/324/00000/5803DFA2-100D-E511-9714-02163E012925.root',
'/store/express/Run2015A/ExpressPhysics/FEVT/Express-v1/000/247/324/00000/A4A5C5D7-100D-E511-A5FD-02163E0138F3.root',
'/store/express/Run2015A/ExpressPhysics/FEVT/Express-v1/000/247/324/00000/90D5FF79-130D-E511-AC72-02163E0138F6.root',
'/store/express/Run2015A/ExpressPhysics/FEVT/Express-v1/000/247/324/00000/E2008FAC-100D-E511-8EA0-02163E013729.root',
'/store/express/Run2015A/ExpressPhysics/FEVT/Express-v1/000/247/324/00000/DA43257E-130D-E511-826E-02163E013952.root',
'/store/express/Run2015A/ExpressPhysics/FEVT/Express-v1/000/247/324/00000/3E1776FA-130D-E511-9A8D-02163E0143D5.root',
'/store/express/Run2015A/ExpressPhysics/FEVT/Express-v1/000/247/324/00000/5E4F37F0-130D-E511-B10F-02163E013952.root',
'/store/express/Run2015A/ExpressPhysics/FEVT/Express-v1/000/247/324/00000/EC4552F1-130D-E511-891C-02163E0138F6.root',
'/store/express/Run2015A/ExpressPhysics/FEVT/Express-v1/000/247/324/00000/EE0E81F7-130D-E511-8E65-02163E012951.root',
'/store/express/Run2015A/ExpressPhysics/FEVT/Express-v1/000/247/324/00000/3082CF1B-140D-E511-AE5E-02163E013949.root',
'/store/express/Run2015A/ExpressPhysics/FEVT/Express-v1/000/247/324/00000/7EDAC114-140D-E511-B586-02163E01475C.root',
'/store/express/Run2015A/ExpressPhysics/FEVT/Express-v1/000/247/324/00000/72ADDA02-140D-E511-BE79-02163E011D3A.root',
'/store/express/Run2015A/ExpressPhysics/FEVT/Express-v1/000/247/324/00000/266D2B11-140D-E511-9E79-02163E014569.root',
'/store/express/Run2015A/ExpressPhysics/FEVT/Express-v1/000/247/324/00000/F24FB4D7-100D-E511-9DB4-02163E013684.root',
'/store/express/Run2015A/ExpressPhysics/FEVT/Express-v1/000/247/324/00000/0CCF7303-140D-E511-A982-02163E011BB0.root',
'/store/express/Run2015A/ExpressPhysics/FEVT/Express-v1/000/247/324/00000/A252D90E-110D-E511-A517-02163E014613.root',
'/store/express/Run2015A/ExpressPhysics/FEVT/Express-v1/000/247/324/00000/1E44709C-100D-E511-9087-02163E011E08.root',
'/store/express/Run2015A/ExpressPhysics/FEVT/Express-v1/000/247/324/00000/F6D82CE7-140D-E511-84E4-02163E013671.root',
'/store/express/Run2015A/ExpressPhysics/FEVT/Express-v1/000/247/324/00000/C477480A-140D-E511-B421-02163E0146F6.root',
'/store/express/Run2015A/ExpressPhysics/FEVT/Express-v1/000/247/324/00000/FE6B8B97-100D-E511-B2E4-02163E012A7A.root',
'/store/express/Run2015A/ExpressPhysics/FEVT/Express-v1/000/247/324/00000/4019E36A-140D-E511-9436-02163E01438B.root',
'/store/express/Run2015A/ExpressPhysics/FEVT/Express-v1/000/247/324/00000/86FFEAFA-130D-E511-A524-02163E014204.root',
'/store/express/Run2015A/ExpressPhysics/FEVT/Express-v1/000/247/324/00000/9EF61809-140D-E511-B0D7-02163E012355.root'

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
                                   fileName = cms.string('PixelTree_dataTest.root')
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
#       process.hltHighPtPA*
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
     

