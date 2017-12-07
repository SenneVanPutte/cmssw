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

//RecHits
#include "DataFormats/TrackerRecHit2D/interface/FastTrackerRecHit.h"
#include "DataFormats/Common/interface/OwnVector.h" 
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

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
	SIM_track_label(consumes<SimTrack>(iPset.getParameter<edm::InputTag>("SIM_track_label"))),
	verbose(iConfig.getUntrackedParameter<int>("verbose",5)),
	fname( iPset.getParameter<string>("outfile"))
{
  //usesResource("TFileService");
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

    Handle<SimTrack> SIM_track_handle;
    iEvent.getByToken(SIM_track_label, SIM_track_handle);

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

    //-----------  SIM TRACKS!!!
    int iSIMTrack = 0;
    SimTrack::const_iterator SIM_itTrk =  SIM_track_handle->begin();
    SimTrack::const_iterator SIM_trkEnd = SIM_track_handle->end();
    if (verbose > 3)
        cout<<">>>>>>>>>>>>>>>>>>> Begining of the SIM track loop"<<endl;
    TH2D * SIMPhiVsEta = new TH2D("Track_scatter","Track Scatter", 10,-3.15,3.15, 10,-2.4,2.4);
    TH1D * hSIMTrackDensity = new TH1D("hTrackDensity","hTrackDensity", 100,0,100);
    for(; SIM_itTrk != SIM_trkEnd; SIM_itTrk++){
      iSIMTrack++;
      auto && p = itTrk->momentum();
      float pt = sqrt(p.perp2());
      float phi = p.phi();
      float eta = p.eta();
      SIMPhiVsEta->Fill(phi,eta);
      hSIMPhiVsEta->Fill(phi,eta);
    }
    if(verbose > 3)
      cout<<"                        END OF SIM TRACK LOOP!"<<endl;
    if (verbose > 4)
      cout<< iTrack <<" TRACKS !!"<<endl;
    double SIMhigherdensity = 0;
    //double ndensity = 0;
    //double sumdensity = 0;
    for (int iBinx = 0; iBinx <= 10; iBinx++){
      for (int iBiny = 0; iBiny <= 10; iBiny++){
        double density = SIMPhiVsEta->GetBinContent(iBinx,iBiny);
        if (density== 0) continue;
        hSIMTrackDensity->Fill(density);
        //sumdensity = sumdensity + density;
        //ndensity++;
        //cout<<"PhiBin = "<< iBinx << " , EtaBin = " << iBiny << ", Density = " << density <<endl;
        if (SIMhigherdensity < density) SIMhigherdensity = density;
      }
    }
    double SIMmeandensity = hSIMTrackDensity->GetMean();
    if(verbose > 3)
      cout<<"SIMhigherdensity = " << SIMhigherdensity <<"; SIMmeandensity = " << meandensity<<endl;
    SIMTrackDensity_higher->Fill(SIMhigherdensity);
    SIMTrackDensity_mean->Fill(SIMmeandensity);

}

// ------------ method called once each job just before starting iEvent loop  ------------
void 
TrackDensityValidator::beginJob()
{
   TrackDensity_higher = new TH1D("TrackDensity_higher","TrackDensity_higher", 100,0,100);
   TrackDensity_mean = new TH1D("TrackDensity_mean","TrackDensity_mean", 100,0,10);
   hPhiVsEta = new TH2D("hPhiVsEta","hPhiVsEta", 10,-3.15,3.15, 10,-2.4,2.4);
   SIMTrackDensity_higher = new TH1D("SIMTrackDensity_higher","SIMTrackDensity_higher", 100,0,100);
   SIMTrackDensity_mean = new TH1D("SIMTrackDensity_mean","SIMTrackDensity_mean", 100,0,10);
   hSIMPhiVsEta = new TH2D("hSIMPhiVsEta","hSIMPhiVsEta", 10,-3.15,3.15, 10,-2.4,2.4);
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
   SIMTrackDensity_higher->Write();
   SIMTrackDensity_mean->Write();
   hSIMPhiVsEta->Write();
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
