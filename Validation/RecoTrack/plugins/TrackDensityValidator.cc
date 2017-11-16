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
//#include "DataFormats/TrackerRecHit2D/interface/SiTrackerGSRecHit2D.h"
//#include "DataFormats/TrackerRecHit2D/interface/SiTrackerGSRecHit2DCollection.h"
//#include "DataFormats/TrackerRecHit2D/interface/SiTrackerGSMatchedRecHit2D.h"
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
class AuxilaryPlots{
public:
      AuxilaryPlots(TString name_,  bool isFS = true, int nBinYW = 125, int nBinZ = 200, int nBinX = 200, int nBinY = 200): name(name_){
	DX = new TH1D("DX", "#Delta x",nBinX,-0.02,0.02);
	DY = new TH1D("DYWide", "#Delta y",nBinYW,-10,10);
	DYLR = new TH1D("DY", "#Delta y",nBinY,-0.02,0.02);
	DZ = new TH1D("DZ", "#Delta z",nBinZ,-0.00002,0.00002);
	const double nBin = 125.;
	if(isFS){
		recErrXX = new TH1D("recXX", "xx",nBin,0,0.02);
		recErrXY = new TH1D("recXY", "xy",nBin,-0.01,0.01);
		recErrYY = new TH1D("recYY", "yy",nBin,0,0.02);
	} else {
		double xxL = 0, xxH = 0, xyL = 0, xyH = 0;
		if(name.Contains("TIB")){
			xxL = 0.;
			xxH = 0.02;
			xyL = -0.01;
			xyH = 0.01;
		} else if(name.Contains("TID")){
			xxL = 0.;
			xxH = 0.5;
			xyL = 0;
			xyH = 3;
		} else if(name.Contains("TOB")){
                        xxL = 0.;
                        xxH = 0.2;
                        xyL = -0.01;
                        xyH = 0.01;
                } else if(name.Contains("TEC")){
                        xxL = 0.;
                        xxH = 0.5;
                        xyL = 0;
                        xyH = 3;
                } else {
                        xxL = 0.;
                        xxH = 0.02;
                        xyL = -0.01;
                        xyH = 0.01;
		}
		recErrXX = new TH1D("recXX", "xx",nBin,xxL,xxH);
		recErrXY = new TH1D("recXY", "xy",nBin,xyL,xyH);
		recErrYY = new TH1D("recYY", "yy",100,0,1);
	}
      }
      void Fill(double dx, double dy, double dz){
	DX->Fill(dx);
        DY->Fill(dy);
        DYLR->Fill(dy);
        DZ->Fill(dz);
      }
      void FillRecHitErr(double xx, double xy, double yy){
	recErrXX->Fill(sqrt(xx));
        recErrXY->Fill(sqrt(xy));
        recErrYY->Fill(sqrt(yy));
      }
      void Write(TDirectory * f){
	f->cd();
	TDirectory * dir = f->mkdir(name); 
	dir->cd();
	TCanvas C(name, name,65,106,650,634);
	gStyle->SetOptFit(1);
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	C.Range(17.40968,-8.734181,66.40232,63.92405);
	C.SetFillColor(0);
	C.SetBorderMode(0);
	C.SetBorderSize(2);
	C.SetLeftMargin(0.1517028);
	C.SetRightMargin(0.03560372);
	C.SetTopMargin(0.05400697);
	C.SetBottomMargin(0.1202091);
	C.SetFrameFillStyle(0);
	C.SetFrameBorderMode(0);
	C.SetFrameFillStyle(0);
	C.SetFrameBorderMode(0);
	C.Divide(4,2);
	C.cd(1);
	DX->Draw();
	C.cd(2);
	DYLR->Draw();
	C.cd(3);
	DZ->Draw();
	C.cd(4);
	recErrXX->Draw();
	C.cd(5);
	recErrXY->Draw();
	C.cd(6);
	recErrYY->Draw();
	C.cd(7);
	DY->Draw();
	C.Write();
        DX->Write();
        DY->Write();
        DYLR->Write();
        DZ->Write();
        //PtVsEta->Write();
        dir->mkdir("recHitErr")->cd();
        recErrXX->Write();
        recErrXY->Write();
        recErrYY->Write();
	dir->cd();
	f->cd();
      }
      TString name;
private:
      TH1D * DX;
      TH1D * DY;
      TH1D * DZ;
      TH1D * DYLR;
      TH1D * recErrXX;
      TH1D * recErrXY;
      TH1D * recErrYY;
};
class SubDetLayersHistograms{
public:
	SubDetLayersHistograms(TString name_, unsigned int n, bool isFS,int nBinYW = 125, int nBinZ = 2000, int nBinX = 200, int nBinY = 200): name(name_),nLayer(n){
		stringstream s;
		for(unsigned int i = 0; i < nLayer; i++){
			s.str("");
			s << name <<"_"<<i+1;
			layerPlots.push_back(new AuxilaryPlots(s.str().c_str(), isFS, nBinYW, nBinZ, nBinX, nBinY));		
		}
	}
	~SubDetLayersHistograms(){}
	void Fill(int ID, double dx, double dy, double dz){
		if(ID > (int)layerPlots.size()){
			cout<<"INVALID LAYER FOR "<<name<<endl;	
			return;
		}
		layerPlots[ID-1]->Fill(dx,dy,dz);
	}
	void FillRecHitErr(int ID, double xx, double xy, double yy){
	        if(ID > (int)layerPlots.size()){
                        cout<<"INVALID LAYER FOR "<<name<<endl;
                        return;
                }
                layerPlots[ID-1]->FillRecHitErr(xx,xy,yy);
        }

	void Write(TDirectory * f){
		f->cd();
		TDirectory * dir = f->mkdir(name);
		for(unsigned int i = 0; i < nLayer; i++){
			layerPlots[i]->Write(dir);
		}
	}
private:
	TString name;
	unsigned int nLayer;
	std::vector<AuxilaryPlots*> layerPlots;
};

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
      void FillAuxilaryPlots(DetId, unsigned int, const TrackerTopology*, double, double, double, double, double, double);

      // ----------member data ---------------------------
      edm::ParameterSet iPset;
      //edm::InputTag track_label;
      edm::EDGetTokenT<reco::TrackCollection> track_label;
      std::vector<edm::EDGetTokenT<PSimHitCollection> > simhit_label;
      std::vector<TString> subDetNames;
      //std::vector<edm::InputTag> simhit_label;
      int verbose;
      bool isFS;
      std::string fname;
      TFile * fout;
      SubDetLayersHistograms * PixleBarrel;
      SubDetLayersHistograms * PixleFwd;
      SubDetLayersHistograms * TIB;
      SubDetLayersHistograms * TID;
      SubDetLayersHistograms * TOB;
      SubDetLayersHistograms * TEC;
      TH1D * NSimRecSameDetId;
      TH1D * TrackDensity;
      TH1D * Rfs;
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
//	track_label( iPset.getParameter<edm::InputTag>("track_label")),
	track_label(consumes<reco::TrackCollection>(iPset.getParameter<edm::InputTag>("track_label"))),
	verbose(iConfig.getUntrackedParameter<int>("verbose",5)),
	isFS( iPset.getParameter<bool>("isFastSimOnly")),
	fname( iPset.getParameter<string>("outfile"))

{
	for (edm::InputTag const& tag : iPset.getParameter<std::vector<edm::InputTag> >("simhit_label")){
		simhit_label.push_back(consumes<PSimHitCollection>(tag));
		subDetNames.push_back(tag.instance());
	}
}


TrackDensityValidator::~TrackDensityValidator()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

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
//    iEvent.getByLabel(track_label, track_handle);
   

    edm::ESHandle<TrackerTopology> tTopoHand;
    iSetup.get<TrackerTopologyRcd>().get(tTopoHand);
    const TrackerTopology *tTopo=tTopoHand.product();

    TH2D * PtVsEta;
    PtVsEta = new TH2D("Track_scatter","Track Scatter", 10,0,200, 10,0,3.);

    int iTrack = 0;
    reco::TrackCollection::const_iterator itTrk = track_handle->begin();
    reco::TrackCollection::const_iterator trkEnd = track_handle->end();
    if (verbose > 4)
	cout<<">>>>>>>>>>>>>>>>>>> Begining of the track loop"<<endl;
    for(; itTrk != trkEnd; itTrk++){
      trackingRecHit_iterator itRecHit = itTrk->recHitsBegin();
      trackingRecHit_iterator rhEnd = itTrk->recHitsEnd();
      int iRechit = 0;
      if (verbose > 4)
	cout<<">>>>>>>>>>>>>>>>>>> Begining of the recHit loop"<<endl;
      for(; itRecHit != rhEnd; itRecHit++){
	if(!(*itRecHit)->isValid()) continue;
	if(verbose > 3)
          std::cout<<"Track "<< iTrack<< " and RecHit "<<iRechit<<": "<<(*itRecHit)->localPosition().x() <<"\t" << (*itRecHit)->localPosition().y() <<"\t"<< (*itRecHit)->localPosition().z() <<std::endl;
	TrackingRecHit * tempHit = (*itRecHit)->clone();
	float xRec = tempHit->localPosition().x();
	float yRec = tempHit->localPosition().y();
	float zRec = tempHit->localPosition().z();
	float pt   = tempHit->pt();
	float eta  = tempHit->eta();
        PtVsEta->Fill(pt,eta);
	float xxRecErr = tempHit->localPositionError().xx();
	float xyRecErr = tempHit->localPositionError().xy();
	float yyRecErr = tempHit->localPositionError().yy();
	if(verbose > 3)
	  std::cout<<"Track "<< iTrack<< " and RecHit "<<iRechit<<": "<<xRec <<"\t" << yRec <<"\t"<< zRec <<std::endl;
	if(verbose > 3)
	  std::cout<<"Looking for detId ";
	DetId detId = tempHit->geographicalId().rawId();
	delete tempHit;
	double minDR = 99999999;
	const PSimHit* simHit = NULL;
        int nSimRecSameDetId = 0;
        std::stringstream labelName;
	TString SelectedLabel = "";
	if (verbose > 4)
	    cout<<">>>>>>>>>>>>>>>>>>> Begining of the simHit collections loop"<<endl;
	for(unsigned int iSimLabel = 0; iSimLabel < simhit_label.size(); iSimLabel++){
	  labelName.str("");
	  labelName << subDetNames[iSimLabel];
          edm::Handle<PSimHitCollection> simhit_handle; 
	  //iEvent.getByLabel(simhit_label[iSimLabel], simhit_handle);  
	  iEvent.getByToken(simhit_label[iSimLabel], simhit_handle);
	  PSimHitCollection::const_iterator itSimHit = simhit_handle->begin();
	  PSimHitCollection::const_iterator itSimHitEnd = simhit_handle->end();
	  if(verbose > 3)
	    std::cout<<simhit_handle->size()<<endl;
	  if(verbose > 0)
		cout<<">>> Subdetector name: "<<labelName.str()<<endl;
	  if (verbose > 4)
	    cout<<">>>>>>>>>>>>>>>>>>> Begining of the simHit loop named "<<labelName.str()<<endl;
	  for (; itSimHit!= itSimHitEnd; itSimHit++) {
	    if (detId == itSimHit->detUnitId()){
	      nSimRecSameDetId++;
	      float xSim = itSimHit->localPosition().x();
	      float ySim = itSimHit->localPosition().y();
	      float zSim = itSimHit->localPosition().z();
	      if(verbose > 0)
		cout<<"SimHit Coordinates: "<<xSim <<"\t"<<ySim<<"\t"<<zSim<<endl;
	      double DR = sqrt(pow((xRec-xSim),2)+pow((yRec-ySim),2)+pow((zRec-zSim),2));
	      if(verbose > 0)
		cout<< "DR is "<<DR<<endl;
	      if(DR < minDR){
		if(verbose > 0)
		  cout<< "keep new simHit :-) "<<endl;
		minDR = DR;
		SelectedLabel = labelName.str().c_str();
		simHit = &(*itSimHit);
	      }
	    }
	  }
	  if(verbose > 3)
	    cout<<"================================= END OF SIMHIT LOOP NAMED "<<labelName.str()<<" ==============================="<<endl;
	}
	if(verbose > 3)
            cout<<"================================= END OF SIMHIT LOOP COLLECTION ==============================="<<endl;    
	if(verbose > 0)
	  std::cout<<"After simHit loop, the closest simHit is: "<<simHit<<endl;
	if(isFS){
	  //const SiTrackerGSMatchedRecHit2D * rechit = (const SiTrackerGSMatchedRecHit2D*) (*itRecHit);
	  const FastTrackerRecHit * rechit = (const FastTrackerRecHit*) (*itRecHit);
	  if(verbose > 0){
	    cout<<"From SiTrackerGSMatched:\n";
	    if(rechit->isMatched()) {
	      cout<<"\tThe RecHit is Matched with SimHit ======= "<<endl;
	      cout<<"position: "<<rechit->localPosition().x()<<"\t"<<rechit->localPosition().y()<<"\t"<<rechit->localPosition().z()<<endl;
	    } else  cout<<"\tThe RecHit is NOT Matched with SimHit"<<endl;
	  }
	  if(rechit->isMatched()){
	    Rfs->Fill(sqrt(pow(rechit->localPosition().x(),2)+pow(rechit->localPosition().y(),2)+pow(rechit->localPosition().z(),2)));
	  }
	}
	if(simHit == NULL)
	  continue;
	if(verbose > 0)
	  std::cout<<"simHit coordinates are: ";
	float xSim = simHit->localPosition().x();
	float ySim = simHit->localPosition().y();
	float zSim = simHit->localPosition().z();
	if(verbose > 0)
	  std::cout<<xSim<<"\t"<<ySim<<"\t"<<zSim<<endl;
	double dx = xRec-xSim;
	double dy = yRec-ySim;
	double dz = zRec-zSim;
	if(verbose > 0)
	  std::cout<<dx<<"\t"<<dy<<"\t"<<dz<<endl;
	//Layer definitions and ID
	unsigned int subDetId = DetId(simHit->detUnitId()).subdetId();
	detId = simHit->detUnitId();
        FillAuxilaryPlots(detId, subDetId, tTopo,dx,dy,dz,xxRecErr,xyRecErr,yyRecErr);
	NSimRecSameDetId->Fill(nSimRecSameDetId);
	iRechit++;
      }
      if(verbose > 3)
	cout<<"================================= END OF RECHITS ==============================="<<endl;
    }
    double higherdensity = 0;
    for ( int iBinx = 0, iBiny = 0; iBinx && iBiny < 10; iBinx++, iBiny++ ){
      if (higherdensity > PtVsEta->GetBinContent(iBinx,iBiny)) higherdensity = PtVsEta->GetBinContent(iBinx,iBiny);
    }
    //Fill tracker density
    TrackDensity->Fill(higherdensity);
    if(verbose > 3)
      cout<<"================================= END OF TRACKS ==============================="<<endl;
  }


// ------------ method called once each job just before starting iEvent loop  ------------
void 
TrackDensityValidator::beginJob()
{
   PixleBarrel = new SubDetLayersHistograms("PixelBarrel",3, isFS);
   PixleFwd = new SubDetLayersHistograms("PixelFwd",2, isFS);
   TIB = new SubDetLayersHistograms("TIB", 4, isFS);
   TID = new SubDetLayersHistograms("TID",3, isFS);
   TOB = new SubDetLayersHistograms("TOB",6, isFS);
   TEC = new SubDetLayersHistograms("TEC",7, isFS);
   NSimRecSameDetId = new TH1D("NSimRecSameDetId","N (SimHit,RecHit)_{same detid}", 11, -5.5, 5.5);
   TrackDensity = new TH1D("TrackDensity","TrackDensity", 100,0,1000);
   Rfs = new TH1D("Rfs","Matched SiTrack2D position (R)", 2000, -10, 10);
}

// ------------ method called once each job just after ending the iEvent loop  ------------
void 
TrackDensityValidator::endJob() 
{
   fout = new TFile(fname.c_str(),"recreate");
   fout->cd();
   PixleBarrel->Write(fout);
   PixleFwd->Write(fout);
   TIB->Write(fout);
   TID->Write(fout);
   TOB->Write(fout);
   TEC->Write(fout);
   fout->cd();  
   NSimRecSameDetId->Write();
   TrackDensity->Write();
   Rfs->Write();
   fout->Close();
}

void
TrackDensityValidator::FillAuxilaryPlots(DetId detId, unsigned int subDetId, const TrackerTopology * tt, double dx, double dy, double dz, double recXX, double recXY, double recYY){
   using namespace std;
   if(subDetId > 6)
	return;
   switch(subDetId){
	case 1:
	    {
                if(verbose > 3 )
			cout << "\tPixel Barrel Layer " << tt->pxbLayer(detId)<<endl;
		PixleBarrel->Fill((int)tt->pxbLayer(detId),dx,dy,dz);
		PixleBarrel->FillRecHitErr((int)tt->pxbLayer(detId),recXX,recXY,recYY);
		break;
	    }
        case 2:
	    {
                if(verbose > 3 )
                	cout << "\tPixel Forward Disk " << tt->pxfDisk(detId)<<endl;
		PixleFwd->Fill((int)tt->pxfDisk(detId),dx,dy,dz);
		PixleFwd->FillRecHitErr((int)tt->pxfDisk(detId),recXX,recXY,recYY);
		break;
	    }
	case 3:
	    {
                if(verbose > 3 )
                	cout << "\tTIB Layer " << tt->tibLayer(detId)<<endl;
		TIB->Fill((int)tt->tibLayer(detId),dx,dy,dz);
		TIB->FillRecHitErr((int)tt->tibLayer(detId),recXX,recXY,recYY);
		break;
	    }
        case 4:
	    {
                if(verbose > 3 )
                	cout << "\tTID " << tt->tidRing(detId)<<endl;
		TID->Fill((int)tt->tidRing(detId), dx,dy,dz);
		TID->FillRecHitErr((int)tt->tidRing(detId),recXX,recXY,recYY);
		break;
	    }
	case 5:
	    {
                if(verbose > 3 )
             		cout << "\tTOB Layer " << tt->tobLayer(detId)<<endl;
		TOB->Fill((int)tt->tobLayer(detId), dx,dy,dz);
		TOB->FillRecHitErr((int)tt->tobLayer(detId),recXX,recXY,recYY);
		break;
	    }
        case 6:
	    {
                if(verbose > 3 )
                	cout << "\tTEC " << tt->tecRing(detId)<<endl;
		TEC->Fill((int)tt->tecRing(detId), dx,dy,dz);
		TEC->FillRecHitErr((int)tt->tecRing(detId),recXX,recXY,recYY);
		break;
	    }
	default:
	    {
   	        cout<< "\t INVALID SUBDETECTOR =========================="<<endl;
	        return;
		break;
	    }
    }    

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

//define this as a plug-in
DEFINE_FWK_MODULE(TrackDensityValidator);
