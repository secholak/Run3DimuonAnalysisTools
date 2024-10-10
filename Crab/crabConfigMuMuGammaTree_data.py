import CRABClient
from CRABClient.UserUtilities import config
config = config()

config.JobType.pluginName = 'Analysis'
# Name of the CMSSW configuration file
config.JobType.psetName = 'muMuGammaTree_data.py'

#config.Data.inputDataset = '/ParkingDoubleMuonLowMass1/Run2022F-PromptReco-v1/MINIAOD'
#config.Data.inputDataset = '/ScoutingPFMonitor/Run2022F-v1/RAW'

# These values only make sense for processing data
#    Select input data based on a lumi mask
# 2022
config.Data.lumiMask = 'Cert_Collisions2022_355100_362760_Golden.json'
# 2023
# config.Data.lumiMask = 'Cert_Collisions2023_366442_370790_Golden.json'

# 2024 Golden is not finalised yet, 
# most fo the eras are there:/eos/user/c/cmsdqm/www/CAF/certification/Collisions24/*_era*_Golden*

# Where the output files will be transmitted to
config.Site.storageSite = 'T2_US_MIT'

if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand

    datasets = [



        # no Golden file yet for 2024
        # 2024B, 0.130/fb
        # '/ParkingDoubleMuonLowMass0/Run2024B-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass1/Run2024B-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass2/Run2024B-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass3/Run2024B-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass4/Run2024B-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass5/Run2024B-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass6/Run2024B-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass7/Run2024B-PromptReco-v1/MINIAOD',

        # 2024C, 7.238/fb
        # '/ParkingDoubleMuonLowMass0/Run2024C-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass1/Run2024C-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass2/Run2024C-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass3/Run2024C-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass4/Run2024C-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass5/Run2024C-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass6/Run2024C-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass7/Run2024C-PromptReco-v1/MINIAOD',

        ## 2024D, 7.957/fb
        # '/ParkingDoubleMuonLowMass0/Run2024D-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass1/Run2024D-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass2/Run2024D-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass3/Run2024D-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass4/Run2024D-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass5/Run2024D-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass6/Run2024D-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass7/Run2024D-PromptReco-v1/MINIAOD',

        ## 2024E, 11.319/fb
        # '/ParkingDoubleMuonLowMass0/Run2024E-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass1/Run2024E-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass2/Run2024E-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass3/Run2024E-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass4/Run2024E-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass5/Run2024E-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass6/Run2024E-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass7/Run2024E-PromptReco-v1/MINIAOD',

        # '/ParkingDoubleMuonLowMass0/Run2024E-PromptReco-v2/MINIAOD',
        # '/ParkingDoubleMuonLowMass1/Run2024E-PromptReco-v2/MINIAOD',
        # '/ParkingDoubleMuonLowMass2/Run2024E-PromptReco-v2/MINIAOD',
        # '/ParkingDoubleMuonLowMass3/Run2024E-PromptReco-v2/MINIAOD',
        # '/ParkingDoubleMuonLowMass4/Run2024E-PromptReco-v2/MINIAOD',
        # '/ParkingDoubleMuonLowMass5/Run2024E-PromptReco-v2/MINIAOD',
        # '/ParkingDoubleMuonLowMass6/Run2024E-PromptReco-v2/MINIAOD',
        # '/ParkingDoubleMuonLowMass7/Run2024E-PromptReco-v2/MINIAOD',

        ## 2024F, 25.790/fb
        # '/ParkingDoubleMuonLowMass0/Run2024F-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass1/Run2024F-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass2/Run2024F-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass3/Run2024F-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass4/Run2024F-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass5/Run2024F-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass6/Run2024F-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass7/Run2024F-PromptReco-v1/MINIAOD',


        ## 2024G, 5.477/fb
        # '/ParkingDoubleMuonLowMass0/Run2024G-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass1/Run2024G-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass2/Run2024G-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass3/Run2024G-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass4/Run2024G-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass5/Run2024G-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass6/Run2024G-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass7/Run2024G-PromptReco-v1/MINIAOD',





        # ######### 2023

        # 23B
        # '/ParkingDoubleMuonLowMass0/Run2023B-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass1/Run2023B-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass2/Run2023B-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass3/Run2023B-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass4/Run2023B-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass5/Run2023B-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass6/Run2023B-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass7/Run2023B-PromptReco-v1/MINIAOD',

        # 23C
        # '/ParkingDoubleMuonLowMass0/Run2023C-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass1/Run2023C-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass2/Run2023C-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass3/Run2023C-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass4/Run2023C-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass5/Run2023C-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass6/Run2023C-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass7/Run2023C-PromptReco-v1/MINIAOD',

        # '/ParkingDoubleMuonLowMass0/Run2023C-PromptReco-v2/MINIAOD',
        # '/ParkingDoubleMuonLowMass1/Run2023C-PromptReco-v2/MINIAOD',
        # '/ParkingDoubleMuonLowMass2/Run2023C-PromptReco-v2/MINIAOD',
        # '/ParkingDoubleMuonLowMass3/Run2023C-PromptReco-v2/MINIAOD',
        # '/ParkingDoubleMuonLowMass4/Run2023C-PromptReco-v2/MINIAOD',
        # '/ParkingDoubleMuonLowMass5/Run2023C-PromptReco-v2/MINIAOD',
        # '/ParkingDoubleMuonLowMass6/Run2023C-PromptReco-v2/MINIAOD',
        # '/ParkingDoubleMuonLowMass7/Run2023C-PromptReco-v2/MINIAOD',

        # '/ParkingDoubleMuonLowMass0/Run2023C-PromptReco-v3/MINIAOD',
        # '/ParkingDoubleMuonLowMass1/Run2023C-PromptReco-v3/MINIAOD',
        # '/ParkingDoubleMuonLowMass2/Run2023C-PromptReco-v3/MINIAOD',
        # '/ParkingDoubleMuonLowMass3/Run2023C-PromptReco-v3/MINIAOD',
        # '/ParkingDoubleMuonLowMass4/Run2023C-PromptReco-v3/MINIAOD',
        # '/ParkingDoubleMuonLowMass5/Run2023C-PromptReco-v3/MINIAOD',
        # '/ParkingDoubleMuonLowMass6/Run2023C-PromptReco-v3/MINIAOD',
        # '/ParkingDoubleMuonLowMass7/Run2023C-PromptReco-v3/MINIAOD',

        # '/ParkingDoubleMuonLowMass0/Run2023C-PromptReco-v4/MINIAOD',
        # '/ParkingDoubleMuonLowMass1/Run2023C-PromptReco-v4/MINIAOD',
        # '/ParkingDoubleMuonLowMass2/Run2023C-PromptReco-v4/MINIAOD',
        # '/ParkingDoubleMuonLowMass3/Run2023C-PromptReco-v4/MINIAOD',
        # '/ParkingDoubleMuonLowMass4/Run2023C-PromptReco-v4/MINIAOD',
        # '/ParkingDoubleMuonLowMass5/Run2023C-PromptReco-v4/MINIAOD',
        # '/ParkingDoubleMuonLowMass6/Run2023C-PromptReco-v4/MINIAOD',
        # '/ParkingDoubleMuonLowMass7/Run2023C-PromptReco-v4/MINIAOD',

        # 23D
        # '/ParkingDoubleMuonLowMass0/Run2023D-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass1/Run2023D-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass2/Run2023D-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass3/Run2023D-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass4/Run2023D-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass5/Run2023D-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass6/Run2023D-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass7/Run2023D-PromptReco-v1/MINIAOD',

        # '/ParkingDoubleMuonLowMass0/Run2023D-PromptReco-v2/MINIAOD',
        # '/ParkingDoubleMuonLowMass1/Run2023D-PromptReco-v2/MINIAOD',
        # '/ParkingDoubleMuonLowMass2/Run2023D-PromptReco-v2/MINIAOD',
        # '/ParkingDoubleMuonLowMass3/Run2023D-PromptReco-v2/MINIAOD',
        # '/ParkingDoubleMuonLowMass4/Run2023D-PromptReco-v2/MINIAOD',
        # '/ParkingDoubleMuonLowMass5/Run2023D-PromptReco-v2/MINIAOD',
        # '/ParkingDoubleMuonLowMass6/Run2023D-PromptReco-v2/MINIAOD',
        # '/ParkingDoubleMuonLowMass7/Run2023D-PromptReco-v2/MINIAOD',


        # ##### ##### ##### #####  2022

        ## 2022C, 
        '/ParkingDoubleMuonLowMass0/Run2022C-10Dec2022-v2/MINIAOD',
        '/ParkingDoubleMuonLowMass1/Run2022C-10Dec2022-v3/MINIAOD',
        '/ParkingDoubleMuonLowMass2/Run2022C-10Dec2022-v2/MINIAOD',
        '/ParkingDoubleMuonLowMass3/Run2022C-10Dec2022-v2/MINIAOD',
        '/ParkingDoubleMuonLowMass4/Run2022C-10Dec2022-v2/MINIAOD',
        '/ParkingDoubleMuonLowMass5/Run2022C-10Dec2022-v2/MINIAOD',
        '/ParkingDoubleMuonLowMass6/Run2022C-10Dec2022-v2/MINIAOD',
        '/ParkingDoubleMuonLowMass7/Run2022C-10Dec2022-v2/MINIAOD',

        # 2022D
        # '/ParkingDoubleMuonLowMass0/Run2022D-10Dec2022-v2/MINIAOD',
        # '/ParkingDoubleMuonLowMass1/Run2022D-10Dec2022-v3/MINIAOD',
        # '/ParkingDoubleMuonLowMass2/Run2022D-10Dec2022-v2/MINIAOD',
        # '/ParkingDoubleMuonLowMass3/Run2022D-10Dec2022-v2/MINIAOD',
        # '/ParkingDoubleMuonLowMass4/Run2022D-10Dec2022-v2/MINIAOD',
        # '/ParkingDoubleMuonLowMass5/Run2022D-10Dec2022-v2/MINIAOD',
        # '/ParkingDoubleMuonLowMass6/Run2022D-10Dec2022-v2/MINIAOD',
        # '/ParkingDoubleMuonLowMass7/Run2022D-10Dec2022-v2/MINIAOD',
        
        # 2022E
        # '/ParkingDoubleMuonLowMass0/Run2022E-10Dec2022-v2/MINIAOD',
        # '/ParkingDoubleMuonLowMass1/Run2022E-10Dec2022-v2/MINIAOD',
        # '/ParkingDoubleMuonLowMass2/Run2022E-10Dec2022-v2/MINIAOD',
        # '/ParkingDoubleMuonLowMass3/Run2022E-10Dec2022-v2/MINIAOD',
        # '/ParkingDoubleMuonLowMass4/Run2022E-10Dec2022-v2/MINIAOD',
        # '/ParkingDoubleMuonLowMass5/Run2022E-10Dec2022-v2/MINIAOD',
        # '/ParkingDoubleMuonLowMass6/Run2022E-10Dec2022-v2/MINIAOD',
        # '/ParkingDoubleMuonLowMass7/Run2022E-10Dec2022-v2/MINIAOD',

        # 2022F
        # '/ParkingDoubleMuonLowMass0/Run2022F-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass1/Run2022F-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass2/Run2022F-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass3/Run2022F-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass4/Run2022F-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass5/Run2022F-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass6/Run2022F-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass7/Run2022F-PromptReco-v1/MINIAOD',
        
        # 2022G     
        # '/ParkingDoubleMuonLowMass0/Run2022G-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass1/Run2022G-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass2/Run2022G-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass3/Run2022G-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass4/Run2022G-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass5/Run2022G-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass6/Run2022G-PromptReco-v1/MINIAOD',
        # '/ParkingDoubleMuonLowMass7/Run2022G-PromptReco-v1/MINIAOD',

    ]
    
    for dataset in datasets:
        config.Data.inputDataset = dataset
        config.General.requestName = "muMuGamma_10Oct24_" + dataset.split('/')[1] +  dataset.split('/')[2]
        crabCommand('submit', config = config)
