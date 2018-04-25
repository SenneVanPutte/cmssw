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
        ignoremissingtkcollection_(iPset.getUntrackedParameter<bool>("ignoremissingtrackcollection",false)),
        useAssociators_(iPset.getParameter< bool >("UseAssociators")),
        associators(iPset.getUntrackedParameter< std::vector<edm::InputTag> >("associators")),
        label(iPset.getParameter< std::vector<edm::InputTag> >("label")),
        parametersDefiner(iPset.getParameter<std::string>("parametersDefiner"))
        //parametersDefiner(consumes<edm::ESHandle<ParametersDefinerForTP>>(iPset.getParameter<edm::InputTag>("parametersDefiner")))
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

  edm::InputTag beamSpotTag = iPset.getParameter<edm::InputTag>("beamSpot");
  bsSrc = consumes<reco::BeamSpot>(beamSpotTag);
  
  tpSelector = TrackingParticleSelector(iPset.getParameter<double>("ptMinTP"),
                                        iPset.getParameter<double>("ptMaxTP"),
                                        iPset.getParameter<double>("minRapidityTP"),
                                        iPset.getParameter<double>("maxRapidityTP"),
                                        iPset.getParameter<double>("tipTP"),
                                        iPset.getParameter<double>("lipTP"),
                                        iPset.getParameter<int>("minHitTP"),
                                        iPset.getParameter<bool>("signalOnlyTP"),
                                        iPset.getParameter<bool>("intimeOnlyTP"),
                                        iPset.getParameter<bool>("chargedOnlyTP"),
                                        iPset.getParameter<bool>("stableOnlyTP"),
                                        iPset.getParameter<std::vector<int> >("pdgIdTP"));

  //usesResource("TFileService");
  if(useAssociators_) {
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

void TrackDensityValidator::tpParametersAndSelection(const TrackingParticleRefVector& tPCeff,
                                                   const ParametersDefinerForTP& parametersDefinerTP,
                                                   const edm::Event& iEvent, const edm::EventSetup& iSetup,
                                                   const reco::BeamSpot& bs,
                                                   std::vector<std::tuple<TrackingParticle::Vector, TrackingParticle::Point> >& momVert_tPCeff,
                                                   std::vector<size_t>& selected_tPCeff) const {
  selected_tPCeff.reserve(tPCeff.size());
  momVert_tPCeff.reserve(tPCeff.size());
  size_t j=0;
  for(auto const& tpr: tPCeff) {
    const TrackingParticle& tp = *tpr;
    // TODO: do we want to fill these from all TPs that include IT
    // and OOT (as below), or limit to IT+OOT TPs passing tpSelector
    // (as it was before)? The latter would require another instance
    // of tpSelector with intimeOnly=False.
    if(tpSelector(tp)) {
      selected_tPCeff.push_back(j);
      TrackingParticle::Vector momentum = parametersDefinerTP.momentum(iEvent,iSetup,tpr);
      TrackingParticle::Point vertex = parametersDefinerTP.vertex(iEvent,iSetup,tpr);
      momVert_tPCeff.emplace_back(momentum, vertex);
    }
    ++j;
  }
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
    TH2D * PtVsEta = new TH2D("PtVsEta","PtVsEta", 30,0,15, 10,-2.4,2.4);
    TH1D * Pt = new TH1D("Pt","Pt", 100,0,100);
    TH1D * Eta = new TH1D("Eta","Eta", 100,-2.4,2.4);
    TH1D * hTrackDensity = new TH1D("hTrackDensity","hTrackDensity", 100,0,100);
    TH1D * hTrackDensity_pt = new TH1D("hTrackDensity_pt","hTrackDensity_pt", 100,0,100);
    for(; itTrk != trkEnd; itTrk++){
      iTrack++;
      auto && p = itTrk->momentum();
      float pt = sqrt(p.perp2());
      float phi = p.phi();
      float eta = p.eta();
      Pt->Fill(pt);
      Eta->Fill(eta);
      PhiVsEta->Fill(phi,eta);
      hPhiVsEta->Fill(phi,eta);
      PtVsEta->Fill(pt,eta);
      hPtVsEta->Fill(pt,eta);
    }
    if(verbose > 3)
      cout<<"                        END OF TRACK LOOP!"<<endl;
    if (verbose > 4)
      cout<< iTrack <<" TRACKS !!"<<endl;
    double higherdensity = 0;
    double higherdensity_pt = 0;
    //double pt_higherdensity = 0;
    //double eta_higherdensity = 0;
    for (int iBinx = 0; iBinx <= 10; iBinx++){
      for (int iBiny = 0; iBiny <= 10; iBiny++){
        double density = PhiVsEta->GetBinContent(iBinx,iBiny);
        if (density== 0) continue;
        hTrackDensity->Fill(density);
        if (higherdensity < density) higherdensity = density;
      }
    }
    for (int iBinx = 0; iBinx <= 20; iBinx++){
      for (int iBiny = 0; iBiny <= 10; iBiny++){
        double density = PtVsEta->GetBinContent(iBinx,iBiny);
        if (density== 0) continue;
        hTrackDensity_pt->Fill(density);
        if (higherdensity_pt < density){ 
          higherdensity_pt = density;
          //pt_higherdensity = iBinx*2.5+1.25;
          //eta_higherdensity = iBiny*0.48-1.92;
        }
      }
    }
    double pt_mean = Pt->GetMean();
    double eta_mean = Eta->GetMean();
    double pt_stdev = Pt->GetStdDev();
    double meandensity = hTrackDensity->GetMean();
    double meandensity_pt = hTrackDensity_pt->GetMean();
    if(verbose > 3)
      cout << "higherdensity = " << higherdensity <<"; meandensity = " << meandensity << endl;
    if(verbose > 3)
      cout << "higherdensity_pt = " << higherdensity_pt <<"; meandensity_pt = " << meandensity_pt << endl;
    if(verbose > 3)
      cout<<"pt  = " << pt_mean <<" +/- " << pt_stdev <<endl;
    TrackDensity_higher->Fill(higherdensity);
    TrackDensity_mean->Fill(meandensity);
    hPt_mean_hd ->Fill(pt_mean,higherdensity);
    hPt_stdev_hd->Fill(pt_stdev,higherdensity);
    hPt_mean_md ->Fill(pt_mean,meandensity);
    hPt_stdev_md->Fill(pt_stdev,meandensity);
    TrackDensity_pt_higher->Fill(higherdensity_pt);
    TrackDensity_pt_mean->Fill(meandensity_pt);
    hPt_mean_hdpt->Fill(pt_mean,higherdensity_pt);
    hEta_mean_hdpt->Fill(eta_mean,higherdensity_pt);
    hPt_mean_mdpt->Fill(pt_mean,meandensity_pt);
    hEta_mean_mdpt->Fill(eta_mean,meandensity_pt);


    ////////////////////////////////////////////////////
    // RECO SIM implementation /from tracking particles
    ///////////////////////////////////////////////////
    TH2D * PhiVsEta_tp = new TH2D("PhiVsEta_tp","PhiVsEta_tp", 10,-3.15,3.15, 10,-2.4,2.4);
    TH2D * PtVsEta_tp = new TH2D("PtVsEta_tp","PtVsEta_tp", 30,0,15, 10,-2.4,2.4);
    TH1D * Pt_tp = new TH1D("Pt_tp","Pt_tp", 100,0,100);
    TH1D * Eta_tp = new TH1D("Eta_tp","Eta_tp", 100,-2.4,2.4);
    TH1D * hTrackDensity_tp = new TH1D("hTrackDensity_tp","hTrackDensity_tp", 100,0,100);
    TH1D * hTrackDensity_tp_pt = new TH1D("hTrackDensity_tp_pt","hTrackDensity_tp_pt", 100,0,100);

    edm::ESHandle<ParametersDefinerForTP> parametersDefinerTPHandle;
    iSetup.get<TrackAssociatorRecord>().get(parametersDefiner,parametersDefinerTPHandle);
    auto parametersDefinerTP = parametersDefinerTPHandle->clone();
    //Since we modify the object, we must clone it

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

    edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
    iEvent.getByToken(bsSrc,recoBeamSpotHandle);
    reco::BeamSpot const & bs = *recoBeamSpotHandle;

    std::vector<size_t> selected_tPCeff;
    std::vector<std::tuple<TrackingParticle::Vector, TrackingParticle::Point>> momVert_tPCeff;
    tpParametersAndSelection(tPCeff, *parametersDefinerTP, iEvent, iSetup, bs, momVert_tPCeff, selected_tPCeff);

    for (unsigned int ww=0;ww<associators.size();ww++){
      // run value filtering of recoToSim map already here as it depends only on the association, not track collection
      reco::RecoToSimCollection const * recSimCollP=nullptr;
      reco::RecoToSimCollection recSimCollL;
      if(!useAssociators_) {
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
        if(useAssociators_){
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
          if(tpSelector(tp)){
            selected_tPCeff.push_back(j);
            TrackingParticle::Vector momentum;
            momentum = parametersDefinerTP->momentum(iEvent,iSetup,tpr);
            TrackingParticle::Point vertex;
            vertex = parametersDefinerTP->vertex(iEvent,iSetup,tpr);
            momVert_tPCeff.emplace_back(momentum, vertex);
          }
          ++j;
        }
   
        // ########################################################
        // fill simulation histograms (LOOP OVER TRACKINGPARTICLES)
        // ########################################################
        //loop over already-selected TPs for tracking efficiency
        if (verbose > 3)
          cout<<">>>>>>>>>>>>>>>>>>> Begining of the track loop"<<endl;
        for(size_t i=0; i<selected_tPCeff.size(); ++i) {
          size_t iTP = selected_tPCeff[i];
          const TrackingParticleRef& tpr = tPCeff[iTP];
          const TrackingParticle& tp = *tpr;
          auto const& momVert = momVert_tPCeff[i];
          const TrackingParticle::Vector& momentumTP = std::get<TrackingParticle::Vector>(momVert);
          float pt = sqrt(momentumTP.perp2());
          float phi = momentumTP.phi();
          float eta = momentumTP.eta();
          if(simRecColl.find(tpr) != simRecColl.end()){
            auto const & rt = simRecColl[tpr];
            if (rt.size()!=0) {
              iTrack_tp++; //This counter counts the number of simTracks that have a recoTrack associated
              PhiVsEta_tp->Fill(phi,eta);
              PtVsEta_tp->Fill(pt,eta);
              Pt_tp->Fill(pt);
              Eta_tp->Fill(eta);
              hPhiVsEta_tp->Fill(phi,eta);
              hPtVsEta_tp->Fill(pt,eta);
            }
          }
        }
        if(verbose > 3)
          cout<<"                        END OF TP TRACK LOOP!"<<endl;
        if (verbose > 4)
          cout<< iTrack_tp <<" TP TRACKS !!"<<endl;
        double higherdensity_tp = 0;
        double higherdensity_tp_pt = 0;
        for (int iBinx = 0; iBinx <= 10; iBinx++){
          for (int iBiny = 0; iBiny <= 10; iBiny++){
            double density = PhiVsEta_tp->GetBinContent(iBinx,iBiny);
            if (density== 0) continue;
            hTrackDensity_tp->Fill(density);
            if (higherdensity_tp < density) higherdensity_tp = density;
          }
        }
        for (int iBinx = 0; iBinx <= 20; iBinx++){
          for (int iBiny = 0; iBiny <= 10; iBiny++){
            double density = PtVsEta_tp->GetBinContent(iBinx,iBiny);
            if (density== 0) continue;
            hTrackDensity_tp_pt->Fill(density);
            if (higherdensity_tp_pt < density) higherdensity_tp_pt = density;
          }
        }
        double pt_mean_tp = Pt_tp->GetMean();
        double eta_mean_tp = Eta_tp->GetMean();
        double pt_stdev_tp = Pt_tp->GetStdDev();
        double meandensity_tp = hTrackDensity_tp->GetMean();
        double meandensity_tp_pt = hTrackDensity_tp_pt->GetMean();
        if(verbose > 3)
          cout << "tp higherdensity = " << higherdensity_tp << "; tp meandensity = " << meandensity_tp << endl;
         if(verbose > 3)
          cout << "tp higherdensity_pt = " << higherdensity_tp_pt << "; tp meandensity_pt = " << meandensity_tp_pt << endl;
        if(verbose > 3)
          cout<<"tp pt  = " << pt_mean_tp <<" +/- " << pt_stdev_tp <<endl;
        TrackDensity_tp_higher->Fill(higherdensity_tp);
        TrackDensity_tp_mean->Fill(meandensity_tp);
        hPt_mean_tp_hd ->Fill(pt_mean_tp,higherdensity_tp);
        hPt_stdev_tp_hd->Fill(pt_stdev_tp,higherdensity_tp); 
        hPt_mean_tp_md ->Fill(pt_mean_tp,meandensity_tp);
        hPt_stdev_tp_md->Fill(pt_stdev_tp,meandensity_tp);
        TrackDensity_tp_pt_higher->Fill(higherdensity_tp_pt);
        TrackDensity_tp_pt_mean->Fill(meandensity_tp_pt);
        hPt_mean_tp_hdpt->Fill(pt_mean_tp,higherdensity_tp_pt);
        hEta_mean_tp_hdpt->Fill(eta_mean_tp,higherdensity_tp_pt);
        hPt_mean_tp_mdpt->Fill(pt_mean_tp,meandensity_tp_pt);
        hEta_mean_tp_mdpt->Fill(eta_mean_tp,meandensity_tp_pt);
      }
    }

}

// ------------ method called once each job just before starting iEvent loop  ------------
void 
TrackDensityValidator::beginJob()
{
   hPt_mean_hd  = new TH2D("hPt_mean_hd","hPt_mean_hd", 20,0,20,15,0,30);
   hPt_stdev_hd = new TH2D("hPt_stdev_hd","hPt_stdev_hd", 20,0,20,25,0,50);
   hPt_mean_md  = new TH2D("hPt_mean_md","hPt_mean_md", 20,0,20,15,0,15);
   hPt_stdev_md = new TH2D("hPt_stdev_md","hPt_stdev_md", 20,0,20,15,0,15);
   TrackDensity_higher = new TH1D("TrackDensity_higher","TrackDensity_higher", 100,0,40);
   TrackDensity_mean = new TH1D("TrackDensity_mean","TrackDensity_mean", 100,0,40);
   hPhiVsEta = new TH2D("hPhiVsEta","hPhiVsEta", 10,-3.15,3.15, 10,-2.4,2.4);
   hPt_mean_tp_hd  = new TH2D("hPt_mean_tp_hd","hPt_mean_tp_hd", 20,0,20,15,0,30);
   hPt_stdev_tp_hd = new TH2D("hPt_stdev_tp_hd","hPt_stdev_tp_hd", 20,0,20,25,0,50);
   hPt_mean_tp_md  = new TH2D("hPt_mean_tp_md","hPt_mean_tp_md", 20,0,20,15,0,15);
   hPt_stdev_tp_md = new TH2D("hPt_stdev_tp_md","hPt_stdev_tp_md", 20,0,20,15,0,15);
   TrackDensity_tp_higher = new TH1D("TrackDensity_tp_higher","TrackDensity_tp_higher", 100,0,40);
   TrackDensity_tp_mean = new TH1D("TrackDensity_tp_mean","TrackDensity_tp_mean", 100,0,40);
   hPhiVsEta_tp = new TH2D("hPhiVsEta_tp","hPhiVsEta_tp", 10,-3.15,3.15, 10,-2.4,2.4);
/// PT and Eta study
   hPt_mean_hdpt  = new TH2D("hPt_mean_hdpt","hPt_mean_hdpt", 20,0,20,15,0,30);
   hEta_mean_hdpt  = new TH2D("hEta_mean_hdpt","hEta_mean_hdpt", 10,-2.4,2.4,15,0,30);
   hPt_mean_mdpt  = new TH2D("hPt_mean_mdpt","hPt_mean_mdpt", 20,0,20,15,0,15);
   hEta_mean_mdpt  = new TH2D("hEta_mean_mdpt","hEta_mean_mdpt", 10,-2.4,2.4,15,0,15);
   TrackDensity_pt_higher = new TH1D("TrackDensity_pt_higher","TrackDensity_pt_higher", 100,0,40); 
   TrackDensity_pt_mean = new TH1D("TrackDensity_pt_mean","TrackDensity_pt_mean", 100,0,40);
   hPtVsEta = new TH2D("hPtVsEta","hPtVsEta", 30,0,15, 10,-2.4,2.4);
   hPt_mean_tp_hdpt  = new TH2D("hPt_mean_tp_hdpt","hPt_mean_tp_hdpt", 20,0,20,15,0,30);
   hEta_mean_tp_hdpt = new TH2D("hEta_mean_tp_hdpt","hEta_mean_tp_hdpt", 10,-2.4,2.4,15,0,30);
   hPt_mean_tp_mdpt  = new TH2D("hPt_mean_tp_mdpt","hPt_mean_tp_mdpt", 20,0,20,15,0,15);
   hEta_mean_tp_mdpt = new TH2D("hEta_mean_tp_mdpt","hEta_mean_tp_mdpt", 10,-2.4,2.4,15,0,15);
   TrackDensity_tp_pt_higher = new TH1D("TrackDensity_tp_pt_higher","TrackDensity_tp_pt_higher", 100,0,40);
   TrackDensity_tp_pt_mean = new TH1D("TrackDensity_tp_pt_mean","TrackDensity_tp_pt_mean", 100,0,40);
   hPtVsEta_tp = new TH2D("hPtVsEta_tp","hPtVsEta_tp", 30,0,15, 10,-2.4,2.4); 
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
   hPt_mean_hd->Write();
   hPt_stdev_hd->Write();
   hPt_mean_md->Write();
   hPt_stdev_md->Write();
   TrackDensity_tp_higher->Write();
   TrackDensity_tp_mean->Write();
   hPhiVsEta_tp->Write();
   hPt_mean_tp_hd->Write();
   hPt_stdev_tp_hd->Write();
   hPt_mean_tp_md->Write(); 
   hPt_stdev_tp_md->Write();
   TrackDensity_pt_higher->Write();
   TrackDensity_pt_mean->Write();
   hPtVsEta->Write();
   hPt_mean_hdpt->Write();
   hEta_mean_hdpt->Write();
   hPt_mean_mdpt->Write();
   hEta_mean_mdpt->Write();
   TrackDensity_tp_pt_higher->Write();
   TrackDensity_tp_pt_mean->Write();
   hPtVsEta_tp->Write();
   hPt_mean_tp_hdpt->Write();
   hEta_mean_tp_hdpt->Write();
   hPt_mean_tp_mdpt->Write();
   hEta_mean_tp_mdpt->Write();
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
