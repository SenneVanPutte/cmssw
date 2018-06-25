// system include files
#include <memory>
#include <cmath>
#include <stdexcept>

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

//TFile service
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

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

#include "Validation/RecoTrack/interface/TracksProducer.h"

//ROOT
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
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
// look in the header (interface/TracksProducer.h)
//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
TracksProducer::TracksProducer(const edm::ParameterSet& iConfig): iPset(iConfig), 
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


TracksProducer::~TracksProducer()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

void TracksProducer::tpParametersAndSelection(const TrackingParticleRefVector& tPCeff,
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
TracksProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    using namespace std;
    TH1::SetDefaultSumw2(kTRUE);
   
    //TD_high.clear();
    //TD_mean.clear();
    //TPtD_high.clear();
    //TPtD_mean.clear();
    //Pt_mean.clear();
    //Pt_stdev.clear();
    //Eta_mean.clear();
       
    //TD_high_tp.clear();
    //TD_mean_tp.clear();
    //TPtD_high_tp.clear();
    //TPtD_mean_tp.clear();
    //Pt_mean_tp.clear();
    //Pt_stdev_tp.clear();
    //Eta_mean_tp.clear();

    if (verbose > 4)
      cout<<">>>>>>>>>>>>>>>>>>> Begining of the event loop"<<endl;

    Handle<reco::TrackCollection> track_handle;
    iEvent.getByToken(track_label, track_handle);

    eventy = iEvent.id().event();
    cout<<"====Event: " << eventy << endl; 

    nrt = 0;
    reco::TrackCollection::const_iterator itTrk = track_handle->begin();
    reco::TrackCollection::const_iterator trkEnd = track_handle->end();
    if (verbose > 3) cout<<">>>>>>>>>>>>>>>>>>> Begining of the track loop"<<endl;
    for(; itTrk != trkEnd; itTrk++){
      auto && p = itTrk->momentum();
      rt_pt[nrt] = sqrt(p.perp2());
      rt_phi[nrt] = p.phi();
      rt_eta[nrt] = p.eta();
      if (p.phi() > 3.15) { 
          std::cout << "pt: " << sqrt(p.perp2()) << ", eta: " <<  p.eta() << ", phi: " << p.phi() << std::endl;
          std::cout << *rt_eta << std::endl;
          throw std::invalid_argument( "phi to high" );
      }
      //std::cout << "rt_pt = " << sqrt(p.perp2()) << " rt_phi = " << p.phi() << " rt_eta = " << p.eta() << std::endl;
      nrt++;
    }
    if(verbose > 3) cout<<"                        END OF TRACK LOOP!"<<endl;
    if (verbose > 4) cout<< nrt <<" TRACKS !!"<<endl;
    //std::cout << "_____" << std::endl;
    ////////////////////////////////////////////////////
    // RECO SIM implementation /from tracking particles
    ///////////////////////////////////////////////////
    edm::ESHandle<ParametersDefinerForTP> parametersDefinerTPHandle;
    iSetup.get<TrackAssociatorRecord>().get(parametersDefiner,parametersDefinerTPHandle);
    auto parametersDefinerTP = parametersDefinerTPHandle->clone();
    //Since we modify the object, we must clone it

    nst = 0; 
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
        if (verbose > 3) cout<<">>>>>>>>>>>>>>>>>>> Begining of the track loop"<<endl;
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
              st_pt[nst] = pt;
              st_eta[nst] = eta;
              st_phi[nst] = phi;
              //std::cout << "st_pt = " << pt << " st_phi = " << phi << " st_eta = " << eta << std::endl;
              nst++; //This counter counts the number of simTracks that have a recoTrack associated
            }
          }
        }
        if(verbose > 3) cout<<"                        END OF TP TRACK LOOP!"<<endl;
        if (verbose > 4) cout<< nst <<" TP TRACKS !!"<<endl;
        tree->Fill();
      }
    }

}

// ------------ method called once each job just before starting iEvent loop  ------------
void 
TracksProducer::beginJob()
{

   fout = new TFile(fname.c_str(),"recreate");

   tree = new TTree("TrackDensity", "TrackDensity");

   Bevent = tree->Branch("event", &eventy, "eventy/I");
   BnRT = tree->Branch("nRT", &nrt, "nRT/I");
   BnST = tree->Branch("nST", &nst, "nST/I");

   BRT_pt = tree->Branch("RT_pt", &rt_pt, "RT_pt[nRT]/F");
   BRT_eta = tree->Branch("RT_eta", &rt_eta, "RT_eta[nRT]/F");
   BRT_phi = tree->Branch("RT_phi", &rt_phi, "RT_phi[nRT]/F");
   
   BST_pt = tree->Branch("ST_pt", &st_pt, "ST_pt[nST]/F");
   BST_eta = tree->Branch("ST_eta", &st_eta, "ST_eta[nST]/F");
   BST_phi = tree->Branch("ST_phi", &st_phi, "ST_phi[nST]/F");

}

// ------------ method called once each job just after ending the iEvent loop  ------------
void 
TracksProducer::endJob() 
{
   fout->cd();
   tree->Write();
   fout->Close();
}

// ------------ method called when starting to processes a run  ------------
/*
void 
TracksProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
TracksProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
TracksProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
TracksProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TracksProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
