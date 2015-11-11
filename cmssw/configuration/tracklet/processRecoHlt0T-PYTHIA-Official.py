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
process.GlobalTag.globaltag = 'MCRUN2_740TV0::All'

process.pixelVertexFromClusters = cms.EDProducer('PixelVertexProducerClusters')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
 '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/00000/0C201EA4-FC0F-E511-B628-D48564594FB4.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/00000/12AB9BD3-5910-E511-AE61-0025905C2CA6.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/00000/162F095B-FB0F-E511-A283-A0369F30FFD2.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/00000/1A3DC3AC-5A10-E511-BDCF-008CFA007CE0.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/00000/1AB7DCD8-5910-E511-B176-0002C94D4656.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/00000/3ED338D3-5910-E511-9937-B083FED04276.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/00000/5087E1F0-5910-E511-8F9B-003048344C68.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/00000/521B7EE2-F90F-E511-BF9F-002590200B4C.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/00000/528DA5EA-FD0F-E511-B54C-002590147CA2.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/00000/52D3E8A8-FC0F-E511-A0E9-002590200A28.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/00000/5C5C3D18-FB0F-E511-B2D8-001E67398E12.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/00000/645342F1-5910-E511-AD83-0025901ABD30.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/00000/6AE90826-FA0F-E511-B22A-00266CFFC13C.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/00000/6CBD389E-5310-E511-B727-00266CFFBF94.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/00000/784AB821-FE0F-E511-B591-001E673983A4.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/00000/84B1B311-FF0F-E511-B003-002590200B4C.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/00000/8A6A9FC0-FA0F-E511-A227-AC162DACC3F8.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/00000/8EA88560-FB0F-E511-A9BF-00266CFEFC5C.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/00000/901AACCC-F90F-E511-8602-0002C90B7F5A.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/00000/948FF1E0-F90F-E511-988E-00266CFE8A04.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/00000/98C3FABC-5910-E511-BF76-0025B31E330A.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/00000/9CEC6FE8-5910-E511-9E5E-842B2B760921.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/00000/A03BAA4D-5310-E511-8FB1-0002C94CD2CE.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/00000/A09A4BAA-FB0F-E511-9362-0002C90C5A4E.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/00000/A6E99916-FB0F-E511-9D33-7845C4FC3A40.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/00000/A8CA4CF9-5910-E511-BA61-20CF305B053E.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/00000/B0237DF7-FC0F-E511-A533-F45214CEF24A.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/00000/B07391CD-F90F-E511-A26F-3417EBE64561.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/00000/BEB03759-FA0F-E511-AB47-002590DB91C6.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/00000/C4F95FDE-F80F-E511-A620-0002C90B39A4.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/00000/E440FB01-0010-E511-8C07-001E673983A4.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/00000/E4652DD8-5910-E511-9846-AC162DAB0B08.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/00000/F694957E-5A10-E511-B2A6-047D7BD6DF26.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/10000/02A242A4-3C10-E511-B138-B8CA3A709648.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/10000/04710838-3A10-E511-81C3-B8CA3A708F98.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/10000/066C726F-3C10-E511-BE1C-0CC47A0AD456.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/10000/0C09EA2D-3B10-E511-ACB0-F45214C748C4.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/10000/0E938593-3A10-E511-A22B-B8CA3A70B608.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/10000/1039806D-3B10-E511-816B-0025901AC3C6.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/10000/10A912B9-3810-E511-ACBD-A0040420FE80.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/10000/1AE03FB4-3810-E511-B44F-B8CA3A708F98.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/10000/1E894088-3E10-E511-9A40-0025901AC3CE.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/10000/2000422E-4110-E511-800B-A0369F301924.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/10000/28690499-A210-E511-970C-0025901AEDA0.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/10000/2E33FDC1-3810-E511-85D7-0025B3E05C32.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/10000/36CB7636-3710-E511-88F7-A0040420FE80.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/10000/3C37A7CE-3810-E511-A0EB-002590200A68.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/10000/48DFEB2B-3810-E511-A2BB-A0040420FE80.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/10000/4A2A3435-A310-E511-8389-A0369F3102F6.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/10000/4CB552D6-3810-E511-8B99-001E673982E6.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/10000/5A4C7325-3810-E511-8287-0002C92DB44E.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/10000/5A8433A4-3910-E511-9A43-001E673982E6.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/10000/6010EE75-3910-E511-B309-001E67396A1D.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/10000/76CE99C4-3810-E511-A97E-0002C90F7FDE.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/10000/7A4BB122-3810-E511-993D-002590A3C970.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/10000/82D40A29-3810-E511-8538-002590A3711C.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/10000/84934CF2-5210-E511-B8AE-002590200A1C.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/10000/98108C23-3810-E511-9231-F45214C748CE.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/10000/9A2B675C-3910-E511-882C-002590D9D9F0.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/10000/9C0968BD-3810-E511-BD48-001E672CC1E7.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/10000/AC3CCEAF-3D10-E511-BB1E-00304867FDDF.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/10000/AE2A66C8-3810-E511-8FEC-F45214C748C0.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/10000/C8D273A2-3910-E511-BAED-0002C90F8088.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/10000/D439FA43-3A10-E511-923B-001E67396E3C.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/10000/DE0D81D6-A110-E511-A79C-001E67396C9D.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/10000/E6FD8413-3C10-E511-8C9E-0002C90F8088.root',
       '/store/mc/RunIISpring15DR74/MinBias_TuneMonash13_13TeV-pythia8/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/10000/FED517B7-3810-E511-95D2-B8CA3A70B608.root'

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
##                             doTrackingParticle = cms.untracked.bool(True),
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
                                   fileName = cms.string('PixelTree-PYTHIA8-official.root')
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
     
