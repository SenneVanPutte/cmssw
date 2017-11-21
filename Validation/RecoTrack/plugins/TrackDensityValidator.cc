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

class TrackDensityValidator : public edm::EDAnalyzer {
   public:
      explicit TrackDensityValidator(const edm::ParameterSet&);
      ~TrackDensityValidator();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      
      typedef std::vector<PSimHit> PSimHitCollection;


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      edm::ParameterSet iPset;
      edm::EDGetTokenT<reco::TrackCollection> track_label;
      int verbose;
      std::string fname;
      TFile * fout;
      TH1D * TrackDensity;
};

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

   if (verbose > 4)
    cout<<">>>>>>>>>>>>>>>>>>> Begining of the event loop"<<endl;

    Handle<reco::TrackCollection> track_handle;
    iEvent.getByToken(track_label, track_handle);

    int iTrack = 0;
    reco::TrackCollection::const_iterator itTrk = track_handle->begin();
    reco::TrackCollection::const_iterator trkEnd = track_handle->end();
    if (verbose > 3)
	cout<<">>>>>>>>>>>>>>>>>>> Begining of the track loop"<<endl;
    TH2D * PtVsEta = new TH2D("Track_scatter","Track Scatter", 20,0,200, 20,-5,5.);
    for(; itTrk != trkEnd; itTrk++){
      iTrack++;
      auto && p = itTrk->momentum();
      float pt = sqrt(p.perp2());
      float eta = p.eta();
      PtVsEta->Fill(pt,eta);
    }
    if(verbose > 3)
      cout<<"                        END OF TRACK LOOP!"<<endl;
    if (verbose > 4)
      cout<< iTrack <<" TRACKS !!"<<endl;
    double higherdensity = 0;
    for ( int iBinx = 0, iBiny = 0; iBinx && iBiny < 10; iBinx++, iBiny++ ){
      if (higherdensity > PtVsEta->GetBinContent(iBinx,iBiny)) higherdensity = PtVsEta->GetBinContent(iBinx,iBiny);
    }
    TrackDensity->Fill(higherdensity);
}

// ------------ method called once each job just before starting iEvent loop  ------------
void 
TrackDensityValidator::beginJob()
{
   TrackDensity = new TH1D("TrackDensity","TrackDensity", 100,0,1000);
}

// ------------ method called once each job just after ending the iEvent loop  ------------
void 
TrackDensityValidator::endJob() 
{
   fout = new TFile(fname.c_str(),"recreate");
   fout->cd();
   TrackDensity->Write();
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
