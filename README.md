# Run3DimuonAnalysisTools
### Repository for run3 analysis tools

#### Setup
Setup a CMSSW working area and clone the `Run3DimuonAnalysisTools` repo:
```
cmsrel CMSSW_12_4_2
cd CMSSW_12_4_2/src
cmsenv
git cms-init
git clone git@github.com:jafriesen/Run3DimuonAnalysisTools.git
scram b -j 8
```

#### Run ntuplizer on ParkingDoubleMuonLowMass0 dataset
Run on one file from the Run2022F/ParkingDoubleMuonLowMass0 dataset:
```
cmsenv
voms-proxy-init --voms cms
cmsRun muMuGammaTree.py
```

#### Script settings

x-check if `L1Info` list has changed if adding a new dataset, for ex 2025.


#### Run on all ParkingDoubleMuonLowMass datasets
Start parallel CRAB jobs on Run2022F/ParkingDoubleMuonLowMass0-7 datasets:
(Modify datasets and job name in crabConfigMuMuGammaTree.py as necessary)
```
cd Crab
cmsenv
source /cvmfs/cms.cern.ch/common/crab-setup.sh
voms-proxy-init --voms cms
python3 crabConfigMuMuGammaTree.py
```

#### Collect job outputs

Check if the jobs are finished with `CRAB/get_jobs_status.py`.
Then collect a list of the output files to be used in the histograms production with `CRAB/get_output_list.py`.
