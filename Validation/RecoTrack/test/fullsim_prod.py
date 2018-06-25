# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step3 --runUnscheduled --conditions auto:run2_mc -s RAW2DIGI,L1Reco,RECO,VALIDATION:tracksValidationTrackingOnly --datatier AODSIM,DQMIO -n 100 --era Run2_2016 --eventcontent AODSIM,DQM --python_filename fullsim.py --fileout fullsim.root --no_exec
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('RECO',eras.Run2_2016)

# import of standard configurations
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.Validation_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)

)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:step2_DIGI_L1_DIGI2RAW_HLT.root'),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step3 nevts:1000'),
#
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.AODSIMoutput = cms.OutputModule("PoolOutputModule",
#
#
#
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(4),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('AODSIM'),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(31457280),
    fileName = cms.untracked.string('fullsim.root'),
    outputCommands = process.AODSIMEventContent.outputCommands
)

process.DQMoutput = cms.OutputModule("DQMRootOutputModule",
#
#
#
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('DQMIO'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('fullsim_inDQM.root'),
    outputCommands = process.DQMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
process.mix.playback = True
process.mix.digitizers = cms.PSet()
for a in process.aliases: delattr(process, a)
process.RandomNumberGeneratorService.restoreStateLabel=cms.untracked.string("randomEngineStateProducer")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

process.rechitanalysis = cms.EDAnalyzer("TracksProducer",
        track_label = cms.InputTag("generalTracks"),
        verbose = cms.untracked.int32(1),
        outfile = cms.string('FullSim_Analyzer.root'),
        UseAssociators = cms.bool(True),
        associators = cms.untracked.VInputTag("quickTrackAssociatorByHits"),
        label = cms.VInputTag("generalTracks"),
        parametersDefiner = cms.string('LhcParametersDefinerForTP'),
        ignoremissingtrackcollection = cms.untracked.bool(False),
        label_tp_effic = cms.InputTag("mix","MergedTrackTruth"),
        label_tp_fake = cms.InputTag("mix","MergedTrackTruth"),
        label_tp_effic_refvector = cms.bool(False),
        label_tp_fake_refvector = cms.bool(False),
        beamSpot = cms.InputTag("offlineBeamSpot"),
        ptMinTP = cms.double(0.9),
        ptMaxTP = cms.double(1e100),
        minRapidityTP = cms.double(-2.5),
        maxRapidityTP = cms.double(2.5),
        tipTP = cms.double(3.5),
        lipTP = cms.double(30.),
        minHitTP = cms.int32(0),
        signalOnlyTP = cms.bool(True),
        intimeOnlyTP = cms.bool(False),
        chargedOnlyTP = cms.bool(True),
        stableOnlyTP = cms.bool(False),
        pdgIdTP = cms.vint32(),
)

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)



process.validation_step = cms.EndPath(process.tracksValidationTrackingOnly)
process.AODSIMoutput_step = cms.EndPath(process.AODSIMoutput)
process.DQMoutput_step = cms.EndPath(process.DQMoutput)
process.analysis_step = cms.EndPath(process.rechitanalysis)

process.source.fileNames = [
 'filetoreplace.root',
]

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.validation_step,process.AODSIMoutput_step,process.DQMoutput_step,process.analysis_step)

from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# customisation of the process.

# Automatic addition of the customisation function from SimGeneral.MixingModule.fullMixCustomize_cff
from SimGeneral.MixingModule.fullMixCustomize_cff import setCrossingFrameOn 

#call to customisation function setCrossingFrameOn imported from SimGeneral.MixingModule.fullMixCustomize_cff
process = setCrossingFrameOn(process)

# End of customisation functions
#do not add changes to your config after this point (unless you know what you are doing)
from FWCore.ParameterSet.Utilities import convertToUnscheduled
process=convertToUnscheduled(process)


# Customisation from command line

#Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
