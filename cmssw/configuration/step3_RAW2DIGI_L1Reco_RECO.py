# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step3 -n 100 --conditions MCRUN2_74_V8 -s RAW2DIGI,L1Reco,RECO --customise SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1 --eventcontent FEVTSIM --magField 0T --beamspot NominalCollision2015 --geometry Extended2015 --no_exec --filein file:step2_DIGI_L1_DIGI2RAW_HLT_RAW2DIGI_L1Reco.root
import FWCore.ParameterSet.Config as cms

process = cms.Process('RECO')
process.Timing = cms.Service("Timing")
# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2015Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_0T_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

process.load('MagneticField.ParametrizedEngine.ParabolicParametrizedField_cfi')
process.ParametrizedMagneticFieldProducer.label = "ParabolicMf"


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:step2_DIGI_L1_DIGI2RAW_HLT_RAW2DIGI_L1Reco.root'),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step3 nevts:100'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.FEVTSIMoutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string(''),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    fileName = cms.untracked.string('step3_RAW2DIGI_L1Reco_RECO.root'),
    outputCommands = process.FEVTSIMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'MCRUN2_74_V8', '')

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
#process.reconstruction_step = cms.Path(process.reconstruction)
#process.reconstruction_step = cms.Path(process.localreco*process.offlineBeamSpot*process.MeasurementTrackerEventPreSplitting*process.siPixelClusterShapeCachePreSplitting*process.standalonemuontracking*process.trackingGlobalReco*process.vertexreco*process.hcalGlobalRecoSequence*process.ecalClusters*process.caloTowersRec)
process.reconstruction_step = cms.Path(process.localreco*process.offlineBeamSpot*process.MeasurementTrackerEventPreSplitting*process.siPixelClusterShapeCachePreSplitting*process.standalonemuontracking*process.trackingGlobalReco*process.vertexreco*process.hcalGlobalRecoSequence*process.caloTowersRec)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTSIMoutput_step = cms.EndPath(process.FEVTSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.endjob_step,process.FEVTSIMoutput_step)

# customisation of the process.

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.postLS1Customs
from SLHCUpgradeSimulations.Configuration.postLS1Customs import customisePostLS1 

#call to customisation function customisePostLS1 imported from SLHCUpgradeSimulations.Configuration.postLS1Customs
process = customisePostLS1(process)

# End of customisation functions

