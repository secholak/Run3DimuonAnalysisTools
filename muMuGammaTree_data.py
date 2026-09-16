import FWCore.ParameterSet.Config as cms

process = cms.Process("SKIM")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)
process.MessageLogger.cerr.FwkSummary.reportEvery = 1000
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    '/store/data/Run2022F/ParkingDoubleMuonLowMass0/MINIAOD/PromptReco-v1/000/360/390/00000/05607e2d-f8e5-41a7-8392-bb5d44c16a6d.root',
    '/store/data/Run2025B/ParkingDoubleMuonLowMass0/MINIAOD/PromptReco-v1/000/392/112/00000/f4e5f4b8-8472-4b8a-9f9e-7505ecf746ff.root', #16% HLt effcy?????
    '/store/data/Run2025D/ParkingDoubleMuonLowMass1/MINIAOD/PromptReco-v1/000/395/948/00000/0deb3098-9c2f-49c5-a0ec-1c225008e891.root',
    '/store/data/Run2025G/ParkingDoubleMuonLowMass0/MINIAOD/PromptReco-v1/000/398/827/00000/2ad0ae6b-9918-48e1-8e4c-22f085708317.root',
    '/store/data/Run2024C/ParkingDoubleMuonLowMass0/MINIAOD/PromptReco-v1/000/380/235/00000/62589bf4-050e-4741-a0cb-a654da2d6998.root',
    '/store/data/Run2024D/ParkingDoubleMuonLowMass0/MINIAOD/PromptReco-v1/000/380/945/00000/e7578f67-7e90-44a9-9fb3-00b692731564.root',    
    # '/store/data/Run2024D/BTagMu/MINIAOD/PromptReco-v1/000/380/306/00000/2b1e5777-67f5-443e-b50e-1a5547eb09ab.root',
 )
)

#process.load("Run3ScoutingAnalysisTools.ScoutingFilter.ScoutingFilter_cff")

#process.load("EventFilter.L1TRawToDigi.gtStage2Digis_cfi")
#process.gtStage2Digis.InputLabel = cms.InputTag( "hltFEDSelectorL1" )

process.TFileService = cms.Service("TFileService", 
    fileName = cms.string("mmgTree.root")
)

#process.ScoutingFilterPath = cms.Path(process.scoutingFilter)

process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.Geometry.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '124X_dataRun3_Prompt_v4', '')

# 2022: run 362616
# L1Info = [
#     "L1_DoubleMu3er2p0_SQ_OS_dR_Max1p4", 
#     "L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6", 
#     "L1_DoubleMu0er1p4_OQ_OS_dEta_Max1p6", 
#     "L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p5", 
#     "L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4", 
#     "L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4", 
#     "L1_DoubleMu4p5_SQ_OS_dR_Max1p2", 
#     "L1_DoubleMu4_SQ_OS_dR_Max1p2"
#     ]


# # 2023: run 370725
# L1Info =  [
#     'L1_DoubleMu0er1p4_SQ_OS_dEta_Max1p2', 
#     'L1_DoubleMu4er2p0_SQ_OS_dR_Max1p6', 
#     'L1_DoubleMu5_SQ_OS_dR_Max1p6', 
#     'L1_DoubleMu3er2p0_SQ_OS_dR_Max1p6', 
#     'L1_DoubleMu0er1p5_SQ_OS_dEta_Max1p2', 
#     'L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6', 
#     'L1_DoubleMu0er1p4_OQ_OS_dEta_Max1p6', 
#     'L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p5', 
#     'L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4', 
#     'L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4', 
#     'L1_DoubleMu4p5_SQ_OS_dR_Max1p2', 
#     'L1_DoubleMu4_SQ_OS_dR_Max1p2'
#     ]

# # 2024: run 383779
L1Info =  [
    "L1_DoubleMu0er1p4_SQ_OS_dEta_Max1p2", 
    "L1_DoubleMu4er2p0_SQ_OS_dR_Max1p6", 
    "L1_DoubleMu5_SQ_OS_dR_Max1p6", 
    "L1_DoubleMu3er2p0_SQ_OS_dR_Max1p6", 
    "L1_DoubleMu0er1p5_SQ_OS_dEta_Max1p2", 
    "L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6", 
    "L1_DoubleMu0er1p4_OQ_OS_dEta_Max1p6", 
    "L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p5", 
    "L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4", 
    "L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4", 
    "L1_DoubleMu4p5_SQ_OS_dR_Max1p2", 
    "L1_DoubleMu4_SQ_OS_dR_Max1p2"
]


# 2025 HLT_DoubleMu4_3_LowMass_v17:
# L1Info = [
#     'L1_DoubleMu0er1p4_SQ_OS_dEta_Max1p2',
#     'L1_DoubleMu4er2p0_SQ_OS_dR_Max1p6',
#     'L1_DoubleMu5_SQ_OS_dR_Max1p6',
#     'L1_DoubleMu3er2p0_SQ_OS_dR_Max1p6',
#     'L1_DoubleMu0er1p5_SQ_OS_dEta_Max1p2',
#     'L1_DoubleMu0er1p4_OQ_OS_dEta_Max1p6',
#     'L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p5',
#     'L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4',
#     'L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4',
#     'L1_DoubleMu4p5_SQ_OS_dR_Max1p2',
#     'L1_DoubleMu4_SQ_OS_dR_Max1p2',

# ]

# this is 25 HLT_DoubleMu4_LowMass_Displaced_v17, seems the same as the prompt one
# ['L1_DoubleMu0er1p4_SQ_OS_dEta_Max1p2',
#  'L1_DoubleMu4er2p0_SQ_OS_dR_Max1p6',
#  'L1_DoubleMu5_SQ_OS_dR_Max1p6',
#  'L1_DoubleMu3er2p0_SQ_OS_dR_Max1p6',
#  'L1_DoubleMu0er1p5_SQ_OS_dEta_Max1p2',
#  'L1_DoubleMu0er1p4_OQ_OS_dEta_Max1p6',
#  'L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p5',
#  'L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4',
#  'L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4',
#  'L1_DoubleMu4p5_SQ_OS_dR_Max1p2',
#  'L1_DoubleMu4_SQ_OS_dR_Max1p2']


process.tree = cms.EDAnalyzer('MuMuGammaTreeMaker',
                                      triggerresults   = cms.InputTag("TriggerResults", "", "HLT"),
                                      ReadPrescalesFromFile = cms.bool( False ),
                                      AlgInputTag       = cms.InputTag("gtStage2Digis"),
                                      l1tAlgBlkInputTag = cms.InputTag("gtStage2Digis"),
                                      l1tExtBlkInputTag = cms.InputTag("gtStage2Digis"),
                                      doL1              = cms.bool( True ),
                                      l1Seeds           = cms.vstring(L1Info),
                                      muons             = cms.InputTag("slimmedMuons"),
                                      electrons         = cms.InputTag("slimmedElectrons"),
                                      photons           = cms.InputTag("slimmedPhotons"),
                                      pfcands           = cms.InputTag("packedPFCandidates"),
                                      prunedGenParticles  = cms.InputTag("prunedGenParticles"),
                                      packedGenParticles  = cms.InputTag("packedGenParticles"),
                                      doGEN             = cms.bool( False ),
                                      doSoftMuonMVA     = cms.bool( True ),
                                      softMvaRun3Model  = cms.string("RecoMuon/MuonIdentification/data/Run2022-20231030-1731-Event0"),
                                      primaryVertices   = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                      displacedVertices = cms.InputTag("slimmedSecondaryVertices"),
                                  )

#process.p = cms.Path(process.gtStage2Digis+process.scoutingTree)
process.p = cms.Path(process.tree)
