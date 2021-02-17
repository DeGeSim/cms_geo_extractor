import os
import FWCore.ParameterSet.Config as cms

processName = "Demo"

from Configuration.StandardSequences.Eras import eras

process = cms.Process(processName, eras.Phase2C9)

# import of standard configurations
process.load("Configuration.StandardSequences.Services_cff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.StandardSequences.L1Reco_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.RecoSim_cff")
process.load("Configuration.StandardSequences.EndOfProcess_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

from Configuration.AlCa.GlobalTag import GlobalTag

process.GlobalTag = GlobalTag(process.GlobalTag, "auto:phase2_realistic_T15", "")

process.load("Configuration.Geometry.GeometryExtended2026D49Reco_cff")
process.load("Configuration.Geometry.GeometryExtended2026D49_cff")


############################## Parse arguments ##############################

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = "INFO"
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.categories.append("Demo")
process.MessageLogger.cerr.INFO = cms.untracked.PSet(limit=cms.untracked.int32(-1))

process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(True))
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(1))


# inputFileNames = cms.untracked.vstring(["file:step1.root"])
# process.source = cms.Source(
#     "PoolSource",
#     fileNames=inputFileNames
# )
process.source = cms.Source("EmptySource")
process.TFileService = cms.Service("TFileService", fileName=cms.string("DetIdLUT.root"))

process.analyzer =cms.EDAnalyzer("GeoExtractor")
process.p = cms.Path(process.analyzer)

process.schedule = cms.Schedule()
process.schedule.insert(0, process.p)

print "\n"
print "*" * 50
print "process.schedule:", process.schedule
print "*" * 50
print "\n"


process.MessageLogger = cms.Service(
    "MessageLogger",
    destinations=cms.untracked.vstring(
        "cerr",
    ),
    cerr=cms.untracked.PSet(
        # threshold  = cms.untracked.string("ERROR"),
        DEBUG=cms.untracked.PSet(limit=cms.untracked.int32(0)),
        WARNING=cms.untracked.PSet(limit=cms.untracked.int32(0)),
        ERROR=cms.untracked.PSet(limit=cms.untracked.int32(0)),
    ),
)