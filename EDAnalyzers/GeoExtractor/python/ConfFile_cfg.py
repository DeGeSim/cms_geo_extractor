import os
import FWCore.ParameterSet.Config as cms
import yaml

processName = "Demo"

from Configuration.Eras.Era_Phase2C22I13M9_cff import Phase2C22I13M9

process = cms.Process('Demo',Phase2C22I13M9)

# import of standard configurations
process.load("Configuration.StandardSequences.Services_cff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.RawToDigi_cff")
#process.load('Geometry.HGCalCommonData.testGeometryV14_cff')
process.load('Configuration.Geometry.GeometryExtended2026D113Reco_cff')
process.load("Configuration.StandardSequences.L1Reco_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.RecoSim_cff")
process.load("Configuration.StandardSequences.EndOfProcess_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

from Configuration.AlCa.GlobalTag import GlobalTag

process.GlobalTag = GlobalTag(process.GlobalTag, "auto:phase2_realistic_T38", "")



############################## Parse arguments ##############################

process.load("FWCore.MessageLogger.MessageLogger_cfi")
# To adjust the loglevel, change it in the GeoExtractor.cc and recompile.

process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(False))

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(1))

process.source = cms.Source("EmptySource")
process.TFileService = cms.Service(
    "TFileService", fileName=cms.string("DetIdLUT.root")
)

with open("settings.yaml", "r") as f:
    settingsD = yaml.load(f,Loader=yaml.SafeLoader)

for key in settingsD:
    print("setting "+str(key)+" <- "+str(settingsD[key]))
    if type(settingsD[key]) is int:
        settingsD[key] = cms.int32(settingsD[key])
    if type(settingsD[key]) is float:
        settingsD[key] = cms.double(settingsD[key])
    if type(settingsD[key]) is bool:
        settingsD[key] = cms.bool(settingsD[key])
    if type(settingsD[key]) is str:
        settingsD[key] = cms.string(settingsD[key])


process.analyzer = cms.EDAnalyzer("GeoExtractor", **settingsD)

process.p = cms.Path(process.analyzer)

process.schedule = cms.Schedule()
process.schedule.insert(0, process.p)
