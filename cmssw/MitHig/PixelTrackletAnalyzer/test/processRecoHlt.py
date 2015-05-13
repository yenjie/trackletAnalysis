import FWCore.ParameterSet.Config as cms

process = cms.Process("myRECO")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("RecoVertex.BeamSpotProducer.BeamSpot_cfi")

# Timing service
process.Timing = cms.Service("Timing") 

# MC Globaltag for 2015 dN/deta analysis
process.GlobalTag.globaltag = 'MCRUN2_74_V6B::All'

process.pixelVertexFromClusters = cms.EDProducer('PixelVertexProducerClusters')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#    'file:RelVal_MinBias_13TeV_28478DD9-99A9-E411-891C-0025905B861C.root'
'/store/relval/CMSSW_7_4_0_pre6/RelValMinBias_13/GEN-SIM-RECO/MCRUN2_74_V1-v1/00000/28478DD9-99A9-E411-891C-0025905B861C.root',
'/store/relval/CMSSW_7_4_0_pre6/RelValMinBias_13/GEN-SIM-RECO/MCRUN2_74_V1-v1/00000/30A73DF2-C8A8-E411-9D7A-003048FFCB9E.root',
'/store/relval/CMSSW_7_4_0_pre6/RelValMinBias_13/GEN-SIM-RECO/MCRUN2_74_V1-v1/00000/42CC9E00-CDA8-E411-95F4-002618943832.root',
'/store/relval/CMSSW_7_4_0_pre6/RelValMinBias_13/GEN-SIM-RECO/MCRUN2_74_V1-v1/00000/A0259BF5-CEA8-E411-9AA1-002618943959.root',
'/store/relval/CMSSW_7_4_0_pre6/RelValMinBias_13/GEN-SIM-RECO/MCRUN2_74_V1-v1/00000/C889930B-98A9-E411-820E-0025905A60A6.root',
'/store/relval/CMSSW_7_4_0_pre6/RelValMinBias_13/GEN-SIM-RECO/MCRUN2_74_V1-v1/00000/CC0FAFEF-C8A8-E411-B985-002618943983.root',
'/store/relval/CMSSW_7_4_0_pre6/RelValMinBias_13/GEN-SIM-RECO/MCRUN2_74_V1-v1/00000/D6EF821B-CBA8-E411-9FDD-003048FFCBA8.root',
'/store/relval/CMSSW_7_4_0_pre6/RelValMinBias_13/GEN-SIM-RECO/MCRUN2_74_V1-v1/00000/EE770C1A-CBA8-E411-AFDE-0025905964CC.root'


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
                             vertexSrc = cms.vstring('offlinePrimaryVerticesWithBS'),
                             trackSrc = cms.InputTag('generalTracks'),
                             doTracking = cms.untracked.bool(False),
                             doCentrality = cms.untracked.bool(True)
                             )

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
                                   fileName = cms.string('PixelTree.root')
                                   )


process.analyze = cms.Path(
       process.siPixelRecHits*
#       process.hiSelectedVertex*
       process.pACentrality*
       process.centralityBin*
       process.hltanalysis*
#       process.hiEvtAnalyzer*
       process.ana)


