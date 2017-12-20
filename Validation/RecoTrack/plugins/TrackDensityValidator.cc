// system include files
#include <memory>
#include <cmath>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/EDConsumerBase.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include "CommonTools/Utils/interface/associationMapFilterValues.h"

//RecHits
#include "DataFormats/TrackerRecHit2D/interface/FastTrackerRecHit.h"
#include "DataFormats/Common/interface/OwnVector.h" 
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/EncodedEventId/interface/EncodedEventId.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"
#include "SimTracker/TrackAssociation/plugins/ParametersDefinerForTPESProducer.h"
#include "SimTracker/TrackAssociation/plugins/CosmicParametersDefinerForTPESProducer.h"
#include "SimTracker/TrackAssociation/interface/TrackingParticleIP.h"

// PSimHits
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"

#include "Validation/RecoTrack/interface/TrackDensityValidator.h"

//ROOT
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TString.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TStyle.h"

//C++
#include <vector>
#include <iostream>
#include <sstream>
#include <string>
using namespace std;
//
// class declaration
//
// look in the header (interface/TrackDensityValidator.h)
//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
TrackDensityValidator::TrackDensityValidator(const edm::ParameterSet& iConfig): iPset(iConfig), 
	track_label(consumes<reco::TrackCollection>(iPset.getParameter<edm::InputTag>("track_label"))),
	verbose(iConfig.getUntrackedParameter<int>("verbose",5)),
	fname( iPset.getParameter<string>("outfile")),
        UseAssociators(iPset.getParameter< bool >("UseAssociators")),
        associators(iPset.getUntrackedParameter< std::vector<edm::InputTag> >("associators")),
        label(iPset.getParameter< std::vector<edm::InputTag> >("label")),
        parametersDefiner(iPset.getParameter<std::string>("parametersDefiner")),
        ignoremissingtkcollection_(iPset.getUntrackedParameter<bool>("ignoremissingtrackcollection",false))
{
  const edm::InputTag& label_tp_effic_tag = iPset.getParameter< edm::InputTag >("label_tp_effic");
  const edm::InputTag& label_tp_fake_tag = iPset.getParameter< edm::InputTag >("label_tp_fake");

  if(iPset.getParameter<bool>("label_tp_effic_refvector")) {
    label_tp_effic_refvector = consumes<TrackingParticleRefVector>(label_tp_effic_tag);
  }
  else {
    label_tp_effic = consumes<TrackingParticleCollection>(label_tp_effic_tag);
  }
  if(iPset.getParameter<bool>("label_tp_fake_refvector")) {
    label_tp_fake_refvector = consumes<TrackingParticleRefVector>(label_tp_fake_tag);
  }
  else {
    label_tp_fake = consumes<TrackingParticleCollection>(label_tp_fake_tag);
  }
  for (auto& itag : label) {
    labelToken.push_back(consumes<edm::View<reco::Track> >(itag));
  }

  //usesResource("TFileService");
  if(UseAssociators) {
    for (auto const& src: associators) {
      associatorTokens.push_back(consumes<reco::TrackToTrackingParticleAssociator>(src));
    }
  } else {
    for (auto const& src: associators) {
      associatormapStRs.push_back(consumes<reco::SimToRecoCollection>(src));
      associatormapRtSs.push_back(consumes<reco::RecoToSimCollection>(src));
    }
  }

}


TrackDensityValidator::~TrackDensityValidator()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


// ------------ method called for each iEvent  ------------
void
TrackDensityValidator::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    using namespace std;
    TH1::SetDefaultSumw2(kTRUE);

    if (verbose > 4)
      cout<<">>>>>>>>>>>>>>>>>>> Begining of the event loop"<<endl;

    Handle<reco::TrackCollection> track_handle;
    iEvent.getByToken(track_label, track_handle);

    int iTrack = 0;
    reco::TrackCollection::const_iterator itTrk = track_handle->begin();
    reco::TrackCollection::const_iterator trkEnd = track_handle->end();
    if (verbose > 3)
	cout<<">>>>>>>>>>>>>>>>>>> Begining of the track loop"<<endl;
    TH2D * PhiVsEta = new TH2D("Track_scatter","Track Scatter", 10,-3.15,3.15, 10,-2.4,2.4);
    TH1D * hTrackDensity = new TH1D("hTrackDensity","hTrackDensity", 100,0,100);
    for(; itTrk != trkEnd; itTrk++){
      iTrack++;
      auto && p = itTrk->momentum();
      float pt = sqrt(p.perp2());
      float phi = p.phi();
      float eta = p.eta();
      PhiVsEta->Fill(phi,eta);
      hPhiVsEta->Fill(phi,eta);
    }
    if(verbose > 3)
      cout<<"                        END OF TRACK LOOP!"<<endl;
    if (verbose > 4)
      cout<< iTrack <<" TRACKS !!"<<endl;
    double higherdensity = 0;
    //double ndensity = 0;
    //double sumdensity = 0;
    for (int iBinx = 0; iBinx <= 10; iBinx++){
      for (int iBiny = 0; iBiny <= 10; iBiny++){
        double density = PhiVsEta->GetBinContent(iBinx,iBiny);
        if (density== 0) continue;
        hTrackDensity->Fill(density);
        //sumdensity = sumdensity + density;
        //ndensity++;
        //cout<<"PhiBin = "<< iBinx << " , EtaBin = " << iBiny << ", Density = " << density <<endl;
        if (higherdensity < density) higherdensity = density;
      }
    }
    double meandensity = hTrackDensity->GetMean();
    if(verbose > 3)
      cout<<"higherdensity = " << higherdensity <<"; meandensity = " << meandensity<<endl;
    TrackDensity_higher->Fill(higherdensity);
    TrackDensity_mean->Fill(meandensity);

    ////////////////////////////////////////////////////
    // RECO SIM implementation from tracking particles
    ///////////////////////////////////////////////////
    TH2D * PhiVsEta_tp = new TH2D("PhiVsEta_tp","PhiVsEta_tp", 10,-3.15,3.15, 10,-2.4,2.4);
    TH1D * hTrackDensity_tp = new TH1D("hTrackDensity_tp","hTrackDensity_tp", 100,0,100);

    edm::ESHandle<ParametersDefinerForTP> parametersDefinerTPHandle;
    iSetup.get<TrackAssociatorRecord>().get(parametersDefiner,parametersDefinerTPHandle);
    auto parametersDefinerTP = parametersDefinerTPHandle->clone();

    int iTrack_tp = 0; 
    int w=0; //counter counting the number of sets of histograms
    TrackingParticleRefVector tmpTPeff;
    TrackingParticleRefVector tmpTPfake;
    const TrackingParticleRefVector *tmpTPeffPtr = nullptr;
    const TrackingParticleRefVector *tmpTPfakePtr = nullptr;

    edm::Handle<TrackingParticleCollection>  TPCollectionHeff;
    edm::Handle<TrackingParticleRefVector>  TPCollectionHeffRefVector;

    const bool tp_effic_refvector = label_tp_effic.isUninitialized();
    if(!tp_effic_refvector) {
      iEvent.getByToken(label_tp_effic, TPCollectionHeff);
      for(size_t i=0, size=TPCollectionHeff->size(); i<size; ++i) {
        tmpTPeff.push_back(TrackingParticleRef(TPCollectionHeff, i));
      }
      tmpTPeffPtr = &tmpTPeff;
    }
    else {
      iEvent.getByToken(label_tp_effic_refvector, TPCollectionHeffRefVector);
      tmpTPeffPtr = TPCollectionHeffRefVector.product();
    }
    if(!label_tp_fake.isUninitialized()) {
      edm::Handle<TrackingParticleCollection> TPCollectionHfake ;
      iEvent.getByToken(label_tp_fake,TPCollectionHfake);
      for(size_t i=0, size=TPCollectionHfake->size(); i<size; ++i) {
      tmpTPfake.push_back(TrackingParticleRef(TPCollectionHfake, i));
      }
      tmpTPfakePtr = &tmpTPfake;
    }
    else {
      edm::Handle<TrackingParticleRefVector> TPCollectionHfakeRefVector;
      iEvent.getByToken(label_tp_fake_refvector, TPCollectionHfakeRefVector);
      tmpTPfakePtr = TPCollectionHfakeRefVector.product();
    }

    TrackingParticleRefVector const & tPCeff = *tmpTPeffPtr;
    TrackingParticleRefVector const & tPCfake = *tmpTPfakePtr;

    std::vector<size_t> selected_tPCeff;
    std::vector<std::tuple<TrackingParticle::Vector, TrackingParticle::Point>> momVert_tPCeff;
    //tpParametersAndSelection(tPCeff, *parametersDefinerTP, event, setup, bs, momVert_tPCeff, selected_tPCeff);

    for (unsigned int ww=0;ww<associators.size();ww++){
      // run value filtering of recoToSim map already here as it depends only on the association, not track collection
      reco::RecoToSimCollection const * recSimCollP=nullptr;
      reco::RecoToSimCollection recSimCollL;
      if(!UseAssociators) {
        Handle<reco::SimToRecoCollection > simtorecoCollectionH;
        iEvent.getByToken(associatormapStRs[ww], simtorecoCollectionH);

        Handle<reco::RecoToSimCollection > recotosimCollectionH;
        iEvent.getByToken(associatormapRtSs[ww],recotosimCollectionH);
        recSimCollP = recotosimCollectionH.product();

        // We need to filter the associations of the fake-TrackingParticle
        // collection only from RecoToSim collection, otherwise the
        // RecoToSim histograms get false entries
        recSimCollL = associationMapFilterValues(*recSimCollP, tPCfake);
        recSimCollP = &recSimCollL;
      }

      for (unsigned int www=0;www<label.size();www++, w++){   
        edm::Handle<View<reco::Track> >  trackCollectionHandle;
        if(!iEvent.getByToken(labelToken[www], trackCollectionHandle)&&ignoremissingtkcollection_)continue;
        const edm::View<reco::Track>& trackCollection = *trackCollectionHandle;
        reco::SimToRecoCollection const * simRecCollP=nullptr;
        reco::SimToRecoCollection simRecCollL;
        //associate tracks
        if(UseAssociators){
          edm::Handle<reco::TrackToTrackingParticleAssociator> theAssociator;
          iEvent.getByToken(associatorTokens[ww], theAssociator);
  
          edm::RefToBaseVector<reco::Track> trackRefs;
          for(edm::View<reco::Track>::size_type i=0; i<trackCollection.size(); ++i) {
            trackRefs.push_back(trackCollection.refAt(i));
          }
  
          recSimCollL = std::move(theAssociator->associateRecoToSim(trackRefs, tPCfake));
          recSimCollP = &recSimCollL;
          // It is necessary to do the association wrt. fake TPs.
          simRecCollL = std::move(theAssociator->associateSimToReco(trackRefs, tPCfake));
          simRecCollP = &simRecCollL;
        }
        reco::RecoToSimCollection const & recSimColl = *recSimCollP;
        reco::SimToRecoCollection const & simRecColl = *simRecCollP;

        size_t j=0;
        for(auto const& tpr: tPCeff) {
          const TrackingParticle& tp = *tpr;
          // TODO: do we want to fill these from all TPs that include IT
          // and OOT (as below), or limit to IT+OOT TPs passing tpSelector
          // (as it was before)? The latter would require another instance
          // of tpSelector with intimeOnly=False.
          selected_tPCeff.push_back(j);
          TrackingParticle::Vector momentum = parametersDefinerTP.momentum(iEvent,iSetup,tpr);
          TrackingParticle::Point vertex = parametersDefinerTP.vertex(iEvent,iSetup,tpr);
          momVert_tPCeff.emplace_back(momentum, vertex);
          ++j;
        }
   
        // ########################################################
        // fill simulation histograms (LOOP OVER TRACKINGPARTICLES)
        // ########################################################
        //loop over already-selected TPs for tracking efficiency
        if (verbose > 3)
          cout<<">>>>>>>>>>>>>>>>>>> Begining of the track loop"<<endl;
        for(size_t i=0; i<selected_tPCeff.size(); ++i) {
          cout << "# selected tp" << i << endl;
          size_t iTP = selected_tPCeff[i];
          const TrackingParticleRef& tpr = tPCeff[iTP];
          const TrackingParticle& tp = *tpr;
          auto const& momVert = momVert_tPCeff[i];
          TrackingParticle::Vector momentumTP;
          float pt = sqrt(momentumTP.perp2());
          float phi = momentumTP.phi();
          float eta = momentumTP.eta();
          if(simRecColl.find(tpr) != simRecColl.end()){
            auto const & rt = simRecColl[tpr];
            if (rt.size()!=0) {
              iTrack_tp++; //This counter counts the number of simTracks that have a recoTrack associated
              PhiVsEta_tp->Fill(phi,eta);
              hPhiVsEta_tp->Fill(phi,eta);
            }
          }
        }
        if(verbose > 3)
          cout<<"                        END OF TP TRACK LOOP!"<<endl;
        if (verbose > 4)
          cout<< iTrack_tp <<" TP TRACKS !!"<<endl;
        double higherdensity_tp = 0;
        for (int iBinx = 0; iBinx <= 10; iBinx++){
          for (int iBiny = 0; iBiny <= 10; iBiny++){
            double density = PhiVsEta_tp->GetBinContent(iBinx,iBiny);
            if (density== 0) continue;
            hTrackDensity_tp->Fill(density);
            if (higherdensity_tp < density) higherdensity_tp = density;
          }
        }
        double meandensity_tp = hTrackDensity_tp->GetMean();
        if(verbose > 3)
          cout<<"tp higherdensity = " << higherdensity_tp <<"; tp meandensity = " << meandensity<<endl;
        TrackDensity_tp_higher->Fill(higherdensity_tp);
        TrackDensity_tp_mean->Fill(meandensity_tp);
      }
    }    

}

// ------------ method called once each job just before starting iEvent loop  ------------
void 
TrackDensityValidator::beginJob()
{
   TrackDensity_higher = new TH1D("TrackDensity_higher","TrackDensity_higher", 100,0,100);
   TrackDensity_mean = new TH1D("TrackDensity_mean","TrackDensity_mean", 100,0,10);
   hPhiVsEta = new TH2D("hPhiVsEta","hPhiVsEta", 10,-3.15,3.15, 10,-2.4,2.4);
   TrackDensity_tp_higher = new TH1D("TrackDensity_tp_higher","TrackDensity_tp_higher", 100,0,100);
   TrackDensity_tp_mean = new TH1D("TrackDensity_tp_mean","TrackDensity_tp_mean", 100,0,10);
   hPhiVsEta_tp = new TH2D("hPhiVsEta_tp","hPhiVsEta_tp", 10,-3.15,3.15, 10,-2.4,2.4);
}

// ------------ method called once each job just after ending the iEvent loop  ------------
void 
TrackDensityValidator::endJob() 
{
   fout = new TFile(fname.c_str(),"recreate");
   fout->cd();
   TrackDensity_higher->Write();
   TrackDensity_mean->Write();
   hPhiVsEta->Write();
   TrackDensity_tp_higher->Write();
   TrackDensity_tp_mean->Write();
   hPhiVsEta_tp->Write();
   fout->Close();
}

// ------------ method called when starting to processes a run  ------------
/*
void 
TrackDensityValidator::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
TrackDensityValidator::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
TrackDensityValidator::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
TrackDensityValidator::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TrackDensityValidator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
