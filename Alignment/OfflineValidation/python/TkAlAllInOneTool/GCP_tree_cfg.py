import FWCore.ParameterSet.Config as cms
#import FWCore.PythonUtilities.LumiList as LumiList

from FWCore.ParameterSet.VarParsing import VarParsing

import json
import os

##Define process
process = cms.Process("validation")

##Argument parsing
options = VarParsing()
options.register("config", "", VarParsing.multiplicity.singleton, VarParsing.varType.string , "AllInOne config")

options.parseArguments()

###Set validation mode
#valiMode = "StandAlone"

##Read in AllInOne config in JSON format
with open(options.config, "r") as configFile:
    config = json.load(configFile)

#Global tag
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag,config["alignment"]["comp"]["globaltag"])


process.load("Configuration.Geometry.GeometryRecoDB_cff")

process.load("CondCore.CondDB.CondDB_cfi")

process.MessageLogger = cms.Service("MessageLogger",
    destinations = cms.untracked.vstring('detailedInfo', 
        'cout')
)

process.source = cms.Source("EmptySource",
    firstRun=cms.untracked.uint32(276237)
    )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)

process.siStripQualityESProducer.ListOfRecordToMerge=cms.VPSet(
	cms.PSet(record = cms.string('SiStripDetCablingRcd'),
            tag = cms.string('')), 
        cms.PSet(record = cms.string('RunInfoRcd'),
            tag = cms.string('')), 
        cms.PSet(record = cms.string('SiStripBadChannelRcd'),
            tag = cms.string('')), 
        cms.PSet(record = cms.string('SiStripBadFiberRcd'),
            tag = cms.string('')), 
        cms.PSet(record = cms.string('SiStripBadModuleRcd'),
            tag = cms.string('')), 
        cms.PSet(record = cms.string('SiStripBadStripRcd'),
            tag = cms.string(''))
)

process.load("DQM.SiStripCommon.TkHistoMap_cff")

  # configuration of the Tracker Geometry Comparison Tool
  # Tracker Geometry Comparison
process.load("Alignment.OfflineValidation.TrackerGeometryCompare_cfi")
  # the input "IDEAL" is special indicating to use the ideal geometry of the release

process.TrackerGeometryCompare.inputROOTFile1 = config['input_ref']+'/Ntuples.root'
process.TrackerGeometryCompare.inputROOTFile2 = config['input_comp']+'/Ntuples.root'
process.TrackerGeometryCompare.moduleList = 'emptyModuleList.txt'
process.TrackerGeometryCompare.outputFile = "GCPtree.root"

process.load("CommonTools.UtilAlgos.TFileService_cfi")  
process.TFileService.fileName = cms.string("TkSurfDeform.root") 

process.TrackerGeometryCompare.levels = [ config['validation']['GCP']['levels'] ]

  ##FIXME!!!!!!!!!
  ##replace TrackerGeometryCompare.writeToDB = False
  ##removed: dbOutputService

process.p = cms.Path(process.TrackerGeometryCompare)

