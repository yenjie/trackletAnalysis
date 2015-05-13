from WMCore.Configuration import Configuration

config = Configuration()

config.section_("General")

config.General.requestName   = 'Ana_PYTHIA8_CUETP8M1'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName  = 'Analysis'
config.JobType.psetName    = 'processReco.py'

config.section_("Data")
config.Data.inputDataset = '/MinBias_TuneCUETP8M1_13TeV-pythia8/azsigmon-MCRUN2_74_V7-v1-GEN-SIM-RAW-RECO-d15047e93fede7d97e3b31895963835c/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.outLFN = '/store/group/phys_heavyions/yjlee/ppMC2015'

config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'

