#ifndef TracksProducer_h
#define TracksProducer_h


#include "CommonTools/RecoAlgos/interface/CosmicTrackingParticleSelector.h"
#include "CommonTools/RecoAlgos/interface/RecoTrackSelectorBase.h"
#include "CommonTools/Utils/interface/DynArray.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "DQMServices/Core/interface/DQMStore.h"

#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include "MagneticField/Engine/interface/MagneticField.h" 
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h" 

#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"
#include "SimDataFormats/Associations/interface/VertexToTrackingVertexAssociator.h"
#include "SimTracker/Common/interface/TrackingParticleSelector.h"
#include "SimTracker/TrackAssociation/interface/ParametersDefinerForTP.h"

#include "Validation/RecoTrack/interface/MTVHistoProducerAlgoForTracker.h"
//#include "Validation/RecoTrack/interface/MultiTrackValidatorBase.h"

#include "TTree.h"
#include "TBranch.h"

const Int_t nMaxTrack = 400;

class TracksProducer : public edm::EDAnalyzer {
   public:
      explicit TracksProducer(const edm::ParameterSet&);
      ~TracksProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

      typedef std::vector<PSimHit> PSimHitCollection;

   protected:
      edm::ParameterSet iPset;
      edm::EDGetTokenT<reco::TrackCollection> track_label;
      int verbose;
      std::string fname;
      const bool ignoremissingtkcollection_;
      const bool useAssociators_;
      std::vector<edm::InputTag> associators;
      edm::EDGetTokenT<TrackingParticleCollection> label_tp_effic;
      edm::EDGetTokenT<TrackingParticleCollection> label_tp_fake;
      edm::EDGetTokenT<TrackingParticleRefVector> label_tp_effic_refvector;
      edm::EDGetTokenT<TrackingParticleRefVector> label_tp_fake_refvector;
      std::vector<edm::InputTag> label;
      std::string parametersDefiner;
      std::vector<edm::EDGetTokenT<edm::View<reco::Track> > > labelToken;
      edm::EDGetTokenT<reco::BeamSpot> bsSrc;


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      
     void tpParametersAndSelection(const TrackingParticleRefVector& tPCeff,
                                   const ParametersDefinerForTP& parametersDefinerTP,
                                   const edm::Event& event, const edm::EventSetup& setup,
                                   const reco::BeamSpot& bs,
                                   std::vector<std::tuple<TrackingParticle::Vector, TrackingParticle::Point> >& momVert_tPCeff,
                                   std::vector<size_t>& selected_tPCeff) const;
      //edm::ParameterSet iPset;
      //edm::EDGetTokenT<reco::TrackCollection> track_label;
      //int verbose;
      //std::string fname;
      //bool UseAssociators;
      //std::vector<edm::InputTag> associators;
      std::vector<edm::EDGetTokenT<reco::TrackToTrackingParticleAssociator>> associatorTokens;
      std::vector<edm::EDGetTokenT<reco::SimToRecoCollection>> associatormapStRs;
      std::vector<edm::EDGetTokenT<reco::RecoToSimCollection>> associatormapRtSs;
      //std::vector<edm::InputTag> label;
      //edm::EDGetTokenT<std::string> parametersDefiner;
      TrackingParticleSelector tpSelector;
      edm::EDGetTokenT<edm::View<reco::Track> > trackCollectionHandle;
      //const bool ignoremissingtkcollection_;
      TFile * fout;

      TTree * tree;
      
      int eventy;
      int nrt;
      int nst;

      //std::vector<float> rt_pt;
      //std::vector<float> rt_eta;
      //std::vector<float> rt_phi;
      
      //std::vector<float> st_pt;
      //std::vector<float> st_eta;
      //std::vector<float> st_phi;
      
      Float_t rt_pt[nMaxTrack];
      Float_t rt_eta[nMaxTrack];
      Float_t rt_phi[nMaxTrack];
      
      Float_t st_pt[nMaxTrack];
      Float_t st_eta[nMaxTrack];
      Float_t st_phi[nMaxTrack];
      
      TBranch * Bevent;
      TBranch * BnRT;
      TBranch * BnST;

      TBranch * BRT_pt;
      TBranch * BRT_eta;
      TBranch * BRT_phi;
      
      TBranch * BST_pt;
      TBranch * BST_eta;
      TBranch * BST_phi;
};

#endif
