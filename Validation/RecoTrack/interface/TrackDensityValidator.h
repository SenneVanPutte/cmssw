#ifndef TrackDensityValidator_h
#define TrackDensityValidator_h


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "Validation/RecoTrack/interface/MultiTrackValidatorBase.h"
#include "Validation/RecoTrack/interface/MTVHistoProducerAlgoForTracker.h"
#include "SimDataFormats/Associations/interface/VertexToTrackingVertexAssociator.h"
#include "CommonTools/RecoAlgos/interface/RecoTrackSelectorBase.h"
#include "SimTracker/TrackAssociation/interface/ParametersDefinerForTP.h"
#include "CommonTools/Utils/interface/DynArray.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include "MagneticField/Engine/interface/MagneticField.h" 
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h" 

#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"

#include "FWCore/ServiceRegistry/interface/Service.h"

#include "SimTracker/Common/interface/TrackingParticleSelector.h"
#include "CommonTools/RecoAlgos/interface/CosmicTrackingParticleSelector.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include <DQMServices/Core/interface/DQMStore.h>

#include "FWCore/Framework/interface/ConsumesCollector.h"

class TrackDensityValidator : public edm::EDAnalyzer {
   public:
      explicit TrackDensityValidator(const edm::ParameterSet&);
      ~TrackDensityValidator();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

      typedef std::vector<PSimHit> PSimHitCollection;


   protected:
      //std::vector<edm::InputTag> associators;
      edm::EDGetTokenT<TrackingParticleCollection> label_tp_effic;
      edm::EDGetTokenT<TrackingParticleCollection> label_tp_fake;
      edm::EDGetTokenT<TrackingParticleRefVector> label_tp_effic_refvector;
      edm::EDGetTokenT<TrackingParticleRefVector> label_tp_fake_refvector;
      //edm::EDGetTokenT<TrackingVertexCollection> label_tv;
      //edm::EDGetTokenT<std::vector<PileupSummaryInfo> > label_pileupinfo;

      //std::vector<edm::EDGetTokenT<std::vector<PSimHit> > > simHitTokens_;

      //std::vector<edm::InputTag> label;
      std::vector<edm::EDGetTokenT<edm::View<reco::Track> > > labelToken;


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      edm::ParameterSet iPset;
      edm::EDGetTokenT<reco::TrackCollection> track_label;
      int verbose;
      std::string fname;
      bool UseAssociators;
      std::vector<edm::InputTag> associators;
      std::vector<edm::EDGetTokenT<reco::TrackToTrackingParticleAssociator>> associatorTokens;
      std::vector<edm::EDGetTokenT<reco::SimToRecoCollection>> associatormapStRs;
      std::vector<edm::EDGetTokenT<reco::RecoToSimCollection>> associatormapRtSs;
      std::vector<edm::InputTag> label;
      std::string parametersDefiner;
      edm::EDGetTokenT<edm::View<reco::Track> > trackCollectionHandle;
      const bool ignoremissingtkcollection_;
      TFile * fout;
      TH1D * TrackDensity_higher;
      TH1D * TrackDensity_mean;
      TH2D * hPhiVsEta;
      TH1D * TrackDensity_tp_higher;
      TH1D * TrackDensity_tp_mean;
      TH2D * hPhiVsEta_tp;

};

#endif
