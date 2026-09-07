import CRABClient
from CRABAPI.RawCommand import crabCommand
from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'mmg_Run3Summer22MiniAODv3_InclusiveDileptonMinBias_28JUL27'

config.JobType.pluginName = 'Analysis'
# Name of the CMSSW configuration file
config.JobType.psetName = 'muMuGammaTree_MC.py'

config.Data.inputDataset = "/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/Run3Summer22MiniAODv3-Pilot_124X_mcRun3_2022_realistic_v12-v5/MINIAODSIM"
config.Data.inputDBS = 'global'
#config.Data.outputDatasetTag = config.General.requestName
config.Data.inputDBS = 'global'
config.Data.splitting = "FileBased"
config.Data.unitsPerJob = 5

# Where the output files will be transmitted to
config.Site.storageSite = 'T3_CH_CERNBOX'


crabCommand('submit', config = config)
