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
'file:data/D66AA394-A0FF-E411-A948-02163E011A20.root',
'file:data/D89C55C6-97FF-E411-A722-02163E0135C3.root',
'file:data/D8A612CB-9EFF-E411-B25D-02163E012927.root',
'file:data/D8DFFBCD-97FF-E411-B957-02163E01257B.root',
'file:data/DA81AD6A-9BFF-E411-98AD-02163E01383E.root',
'file:data/DAECA1C8-97FF-E411-9DF1-02163E0146A5.root',
'file:data/DC2BB5BA-97FF-E411-AA84-02163E014548.root',
'file:data/DE0CD14D-9CFF-E411-AEA5-02163E014281.root',
'file:data/E052A280-9BFF-E411-B3A6-02163E012A04.root',
'file:data/E26D70EC-99FF-E411-B3C0-02163E013476.root',
'file:data/E2F5B081-9BFF-E411-B483-02163E01268A.root',
'file:data/E4158982-A0FF-E411-937C-02163E014206.root',
'file:data/E4344C95-9CFF-E411-BF6E-02163E0133F2.root',
'file:data/E4FF7EE5-94FF-E411-810D-02163E012AC1.root',
'file:data/E6309482-A0FF-E411-8C81-02163E013809.root',
'file:data/EC85B1D6-97FF-E411-8EC6-02163E011DA4.root',
'file:data/EE21B3E9-99FF-E411-A02F-02163E0138F8.root',
'file:data/F01B307E-9BFF-E411-8DAA-02163E0138F8.root',
'file:data/F66662E3-94FF-E411-90F9-02163E0142B3.root',
'file:data/F6EE427E-9BFF-E411-AF01-02163E0140EE.root',
'file:data/0266BF83-A0FF-E411-B7DD-02163E0143AD.root',
'file:data/04C79FC9-9EFF-E411-B5D9-02163E011A20.root',
'file:data/06B4300D-9AFF-E411-9303-02163E014295.root',
'file:data/0C7358DF-94FF-E411-8EB8-02163E012324.root',
'file:data/0ECA3BCD-97FF-E411-A8C7-02163E01432C.root',
'file:data/10789C86-A0FF-E411-82F4-02163E0120B1.root',
'file:data/10A14B1A-94FF-E411-AFCF-02163E013809.root',
'file:data/10AC2C86-9BFF-E411-AC8A-02163E012927.root',
'file:data/10CF00D2-99FF-E411-A052-02163E014281.root',
'file:data/10D940C9-9EFF-E411-932F-02163E01375C.root',
'file:data/10FD326B-9BFF-E411-82F5-02163E012011.root',
'file:data/12BC4169-9BFF-E411-9946-02163E01432C.root',
'file:data/12DD0A11-96FF-E411-89EC-02163E012927.root',
'file:data/1857869C-A0FF-E411-8A58-02163E0121DA.root',
'file:data/1891A2AC-99FF-E411-814C-02163E01257B.root',
'file:data/18A08FCE-9EFF-E411-8027-02163E0133E3.root',
'file:data/18ED8886-A0FF-E411-8E54-02163E012927.root',
'file:data/18F36A85-A0FF-E411-A4B0-02163E0136B3.root',
'file:data/1A2813BE-97FF-E411-90E1-02163E01369F.root',
'file:data/1A7513C9-9EFF-E411-9076-02163E011D70.root',
'file:data/1AC5DE83-A0FF-E411-91A4-02163E011DE2.root',
'file:data/1E9B4386-A0FF-E411-8C08-02163E0133D9.root',
'file:data/1EAD2677-9CFF-E411-AE63-02163E011DE2.root',
'file:data/2015E3CB-9EFF-E411-ACE2-02163E0144DC.root',
'file:data/20D1C684-A0FF-E411-9870-02163E01382B.root',
'file:data/24559B8C-9BFF-E411-B96D-02163E01262E.root',
'file:data/24959BE9-99FF-E411-BAA9-02163E01432C.root',
'file:data/24E5766C-9BFF-E411-9A84-02163E0133D9.root',
'file:data/285AD4CB-97FF-E411-867A-02163E014652.root',
'file:data/2ADA8977-9BFF-E411-96DD-02163E014295.root',
'file:data/2ADAFD2F-97FF-E411-AB60-02163E0142E6.root',
'file:data/2C0CE7CB-97FF-E411-B69A-02163E0142B3.root',
'file:data/2E8C82CC-97FF-E411-B0C0-02163E013476.root',
'file:data/30A76C76-9BFF-E411-9350-02163E014295.root',
'file:data/30C3D0B8-97FF-E411-8E52-02163E01289B.root',
'file:data/32F2C17E-9BFF-E411-A0FF-02163E0133F2.root',
'file:data/346204C7-9EFF-E411-8C01-02163E01384F.root',
'file:data/3498EECB-97FF-E411-94DD-02163E01383E.root',
'file:data/3834DE17-94FF-E411-AEC8-02163E0138F8.root',
'file:data/3864932C-95FF-E411-AB83-02163E0143AD.root',
'file:data/38FCE789-A0FF-E411-B123-02163E012A04.root',
'file:data/3A76A0C7-97FF-E411-B83D-02163E0142DF.root',
'file:data/3C0E6AEE-99FF-E411-96CB-02163E012011.root',
'file:data/3C11AA0A-9AFF-E411-A2D7-02163E01391B.root',
'file:data/3C337869-9BFF-E411-BFE8-02163E011D70.root',
'file:data/3C8DA6A8-9CFF-E411-ABB9-02163E014274.root',
'file:data/3EE178D5-97FF-E411-BEDF-02163E013809.root',
'file:data/40ACAE7E-9BFF-E411-9805-02163E011DE2.root',
'file:data/40FF2A69-9BFF-E411-ABC1-02163E0142B3.root',
'file:data/423F95AA-99FF-E411-85A2-02163E01201A.root',
'file:data/44011ADC-97FF-E411-8A2A-02163E012011.root',
'file:data/44E8492C-95FF-E411-B19F-02163E014652.root',
'file:data/44F094B3-97FF-E411-8FC0-02163E012324.root',
'file:data/44F65BB9-97FF-E411-B1A4-02163E014668.root',
'file:data/4673FD4D-9CFF-E411-A61E-02163E01466B.root',
'file:data/467F857E-9BFF-E411-9F16-02163E014556.root',
'file:data/46B08FE5-94FF-E411-A4D1-02163E013476.root',
'file:data/481A37D4-9CFF-E411-B38F-02163E014218.root',
'file:data/4842382B-95FF-E411-A32F-02163E0120B1.root',
'file:data/48FD8C1C-94FF-E411-940E-02163E012A5E.root',
'file:data/4AA980C5-9EFF-E411-B98E-02163E01384F.root',
'file:data/4ADB4C84-A0FF-E411-BBDB-02163E0144DC.root',
'file:data/4C0120C7-9EFF-E411-8108-02163E01391B.root',
'file:data/4C1F5417-A2FF-E411-978D-02163E014613.root',
'file:data/4C726EEE-99FF-E411-B0A8-02163E011D70.root',
'file:data/4CDADFE7-99FF-E411-B13C-02163E01289B.root',
'file:data/4CDED692-A8FF-E411-8511-02163E01201A.root',
'file:data/4CE77CED-99FF-E411-A136-02163E014218.root',
'file:data/52D1A4EE-99FF-E411-B7AF-02163E0133E3.root',
'file:data/52D57FEC-99FF-E411-8AF2-02163E011A20.root',
'file:data/54579120-94FF-E411-9CA8-02163E014274.root',
'file:data/54953C9C-9DFF-E411-B43D-02163E0133E8.root',
'file:data/585429A8-A0FF-E411-863C-02163E014652.root',
'file:data/5A03E59F-9BFF-E411-8D1C-02163E01201A.root',
'file:data/5A096497-9CFF-E411-9894-02163E0133E3.root',
'file:data/5A761CC2-97FF-E411-B8AF-02163E01457D.root',
'file:data/5AFF9DAD-9CFF-E411-9D3D-02163E0137F0.root',
'file:data/5CBE07B6-96FF-E411-BF2A-02163E014150.root',
'file:data/5CF0C0CE-97FF-E411-BBD7-02163E01201A.root',
'file:data/5E1D9931-97FF-E411-9418-02163E014239.root',
'file:data/5E692FB4-A0FF-E411-B7AC-02163E01466B.root',
'file:data/5EB2D9E5-94FF-E411-AFA3-02163E0133D9.root',
'file:data/605F7EEC-99FF-E411-B155-02163E014218.root',
'file:data/606FC16A-9BFF-E411-BADB-02163E01383E.root',
'file:data/625C00D2-97FF-E411-93A4-02163E013662.root',
'file:data/629047F2-97FF-E411-84FF-02163E0133E3.root',
'file:data/62AE8ABE-97FF-E411-8943-02163E014146.root',
'file:data/6412B87C-9BFF-E411-97AA-02163E012387.root',
'file:data/64C7EC96-9CFF-E411-A9AF-02163E013809.root',
'file:data/66072120-98FF-E411-BDE2-02163E014295.root',
'file:data/68CCF618-94FF-E411-A2CC-02163E0137FE.root',
'file:data/6A5A89E6-99FF-E411-A039-02163E012927.root',
'file:data/6AFA4997-9CFF-E411-95AD-02163E012A04.root',
'file:data/6CB723CC-97FF-E411-BFBD-02163E0143AD.root',
'file:data/6E109A83-A0FF-E411-A0B6-02163E01384F.root',
'file:data/702C1CC9-97FF-E411-9E36-02163E014175.root',
'file:data/70A463BC-9EFF-E411-8921-02163E01425D.root',
'file:data/72FFD185-A0FF-E411-90F4-02163E01201A.root',
'file:data/767174C1-97FF-E411-9594-02163E01251E.root',
'file:data/76CFAD85-A0FF-E411-BB40-02163E012324.root',
'file:data/7805853B-95FF-E411-8DDE-02163E01424A.root',
'file:data/78B624D2-9EFF-E411-9B6C-02163E01382B.root',
'file:data/7ADB89EB-99FF-E411-BB60-02163E01201A.root',
'file:data/7AEFA569-9BFF-E411-A045-02163E01432C.root',
'file:data/7E0F3070-9BFF-E411-9A92-02163E0142B3.root',
'file:data/7EC84DB1-97FF-E411-BCBB-02163E0146A9.root',
'file:data/80406453-A6FF-E411-956C-02163E0145AE.root',
'file:data/823C23CB-9EFF-E411-BD52-02163E013809.root',
'file:data/828DBEC7-92FF-E411-9D5E-02163E0137FE.root',
'file:data/82BB4690-A0FF-E411-8A8F-02163E012B06.root',
'file:data/845FBD87-9BFF-E411-AABE-02163E014652.root',
'file:data/86ACB38D-A0FF-E411-9618-02163E012B02.root',
'file:data/8A18AC7F-9BFF-E411-B119-02163E012824.root',
'file:data/8C029A78-9BFF-E411-B486-02163E0120B1.root',
'file:data/8C6D47EE-99FF-E411-868A-02163E0133E3.root',
'file:data/8C768BC9-9EFF-E411-839E-02163E0133D9.root',
'file:data/8C99E62C-95FF-E411-ABAB-02163E01384F.root',
'file:data/8E19A4E8-94FF-E411-B0BE-02163E01257B.root',
'file:data/8E71B7C8-9EFF-E411-8F3D-02163E0142B3.root',
'file:data/9000CE0C-9AFF-E411-9F14-02163E0133D9.root',
'file:data/9006FCEC-99FF-E411-9336-02163E0121DA.root',
'file:data/90BF7680-9BFF-E411-9E06-02163E01391B.root',
'file:data/9420827E-A0FF-E411-80BD-02163E014565.root',
'file:data/94C6DF7F-9BFF-E411-B867-02163E01289B.root',
'file:data/962788EE-99FF-E411-A1D7-02163E011D70.root',
'file:data/982691F3-99FF-E411-A7D3-02163E012011.root',
'file:data/98FB8A9E-9DFF-E411-A71A-02163E012011.root',
'file:data/9A4151D6-A2FF-E411-8F3E-02163E01432C.root',
'file:data/9A444B77-9BFF-E411-979B-02163E0142B3.root',
'file:data/9E4CDC0C-96FF-E411-9378-02163E0133C3.root',
'file:data/A02D44F5-99FF-E411-8990-02163E0120B1.root',
'file:data/A0A56AA5-A0FF-E411-BEDF-02163E0133E8.root',
'file:data/A0EB82E4-94FF-E411-B2A4-02163E0137FE.root',
'file:data/A2D7C8CB-97FF-E411-90E4-02163E014652.root',
'file:data/A2FEA36B-9BFF-E411-BDF6-02163E01443D.root',
'file:data/A41A4BE4-94FF-E411-887B-02163E014556.root',
'file:data/A47111C7-9EFF-E411-ABE0-02163E01383E.root',
'file:data/A4F5F57F-9BFF-E411-A60E-02163E014175.root',
'file:data/A80D6DC1-9EFF-E411-AB78-02163E0121DA.root',
'file:data/AA4D66CE-97FF-E411-AE27-02163E012A5E.root',
'file:data/AAA5ABC0-9EFF-E411-9A02-02163E0133E8.root',
'file:data/B04A819C-9CFF-E411-AB7B-02163E012011.root',
'file:data/B0DE3C6C-9BFF-E411-9A7A-02163E0133D9.root',
'file:data/B238B06A-9BFF-E411-BFFA-02163E01416C.root',
'file:data/B2E200C7-9EFF-E411-ABD9-02163E0143AD.root',
'file:data/B496EE64-A8FF-E411-855C-02163E014613.root',
'file:data/B6B5E0BF-9EFF-E411-AA40-02163E011DE2.root',
'file:data/B6DE42E7-A5FF-E411-922B-02163E012BD8.root',
'file:data/B82ECAE7-94FF-E411-BEC4-02163E01432C.root',
'file:data/B84715ED-99FF-E411-ABCA-02163E012028.root',
'file:data/B875E6CC-97FF-E411-997B-02163E01384F.root',
'file:data/B8C3E70C-9AFF-E411-8A28-02163E0133D9.root',
'file:data/B8DBD3CE-97FF-E411-84CA-02163E014556.root',
'file:data/B8E0CAC6-9EFF-E411-9ED7-02163E014556.root',
'file:data/BE48B779-9BFF-E411-B631-02163E013476.root',
'file:data/BEE1AE7F-9BFF-E411-BC24-02163E01382B.root',
'file:data/C0E371E2-9EFF-E411-80CE-02163E012AC1.root',
'file:data/C464E07E-9BFF-E411-B54A-02163E01375C.root',
'file:data/C48506C3-97FF-E411-ADE5-02163E0143AD.root',
'file:data/C4C54AC1-97FF-E411-BACD-02163E013476.root',
'file:data/C4EB9989-A0FF-E411-849C-02163E014556.root',
'file:data/C825206B-9BFF-E411-8F2C-02163E012011.root',
'file:data/C82D6F82-9BFF-E411-92B3-02163E0133E3.root',
'file:data/C8E667CF-97FF-E411-A7FE-02163E012927.root',
'file:data/CA6BE9EC-99FF-E411-81D5-02163E011A20.root',
'file:data/CCC5B5CB-A1FF-E411-9C10-02163E012927.root',
'file:data/CCFB7E82-A0FF-E411-AAAD-02163E01432C.root',
'file:data/D0314A87-A0FF-E411-8F91-02163E012A5E.root',
'file:data/D0E1412C-A4FF-E411-BC61-02163E0133E3.root',
'file:data/D2311CAB-99FF-E411-9433-02163E01469E.root',
'file:data/D271DB86-A0FF-E411-A956-02163E012AC1.root',
'file:data/D2966785-A0FF-E411-8929-02163E0133E3.root',
'file:data/D2A0E9E6-94FF-E411-9946-02163E014274.root',
'file:data/D2BB1882-A0FF-E411-B4AF-02163E0133F2.root',
'file:data/D41BD39A-9DFF-E411-89DB-02163E0146C6.root',
'file:data/F8BD4FED-99FF-E411-BEBA-02163E013476.root',
'file:data/F8C73A96-9CFF-E411-AA32-02163E01384F.root',
'file:data/FC673D87-A0FF-E411-891B-02163E012795.root',
'file:data/FE0BB7B7-9BFF-E411-A15E-02163E0133E8.root',
'file:data/FE101979-9BFF-E411-9518-02163E0120B1.root',

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
#       process.ana*
       process.anaStrip
#       process.SiStripRecHitsAnalyzer
      )
     

