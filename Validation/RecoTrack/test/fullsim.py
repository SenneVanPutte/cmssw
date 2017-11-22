# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step3 --conditions auto:run2_mc -n 10 --era Run2_25ns --magField 38T_PostLS1 -s RAW2DIGI,RECO,VALIDATION:tracksValidationTrackingOnly --datatier AODSIM,DQMIO --eventcontent AODSIM,DQM --filein das:/RelValTTbar_13/CMSSW_9_3_0_pre1-92X_mcRun2_asymptotic_v2-v1/GEN-SIM-DIGI-RAW-HLTDEBUG --python_filename fullsim.py --fileout fullsim.root --no_exec
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('RECO',eras.Run2_25ns)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.Validation_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('/store/relval/CMSSW_9_3_0_pre1/RelValTTbar_13/GEN-SIM-DIGI-RAW-HLTDEBUG/92X_mcRun2_asymptotic_v2-v1/00000/24956E50-7C63-E711-BA0C-003048FFD71C.root', 
        '/store/relval/CMSSW_9_3_0_pre1/RelValTTbar_13/GEN-SIM-DIGI-RAW-HLTDEBUG/92X_mcRun2_asymptotic_v2-v1/00000/488F6E50-7C63-E711-B8F5-003048FFD71C.root', 
        '/store/relval/CMSSW_9_3_0_pre1/RelValTTbar_13/GEN-SIM-DIGI-RAW-HLTDEBUG/92X_mcRun2_asymptotic_v2-v1/00000/589EEE98-7B63-E711-AF06-0CC47A4D767A.root', 
        '/store/relval/CMSSW_9_3_0_pre1/RelValTTbar_13/GEN-SIM-DIGI-RAW-HLTDEBUG/92X_mcRun2_asymptotic_v2-v1/00000/72E48785-7C63-E711-9547-0CC47A78A4B0.root', 
        '/store/relval/CMSSW_9_3_0_pre1/RelValTTbar_13/GEN-SIM-DIGI-RAW-HLTDEBUG/92X_mcRun2_asymptotic_v2-v1/00000/F026299D-7B63-E711-9530-0025905B859E.root'),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step3 nevts:10'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.AODSIMoutput = cms.OutputModule("PoolOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(4),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('AODSIM'),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
    fileName = cms.untracked.string('fullsim.root'),
    outputCommands = process.AODSIMEventContent.outputCommands
)

process.DQMoutput = cms.OutputModule("DQMRootOutputModule",
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

process.rechitanalysis = cms.EDAnalyzer("TrackDensityValidator",
        track_label = cms.InputTag("generalTracks"),
        verbose = cms.untracked.int32(5),
        outfile = cms.string('FullSim_Analyzer.root'),
)

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.reconstruction_step = cms.Path(process.reconstruction)
process.validation_step = cms.EndPath(process.tracksValidationTrackingOnly)
process.AODSIMoutput_step = cms.EndPath(process.AODSIMoutput)
process.DQMoutput_step = cms.EndPath(process.DQMoutput)
process.analysis_step = cms.EndPath(process.rechitanalysis)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.reconstruction_step,process.analysis_step,process.validation_step,process.AODSIMoutput_step,process.DQMoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# customisation of the process.

# Automatic addition of the customisation function from SimGeneral.MixingModule.fullMixCustomize_cff
from SimGeneral.MixingModule.fullMixCustomize_cff import setCrossingFrameOn 

#call to customisation function setCrossingFrameOn imported from SimGeneral.MixingModule.fullMixCustomize_cff
process = setCrossingFrameOn(process)

# End of customisation functions

# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
