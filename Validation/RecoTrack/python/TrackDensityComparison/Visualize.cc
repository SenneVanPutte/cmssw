TFile *fComp = new TFile("compare.root");

TIter next(fComp->GetListOfKeys());
TKey *key;

while ((key = (TKey*)next())) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TH1")) continue;
    TH1 *h = (TH1*)key->ReadObj();
   }
//TObjArray obj = fComp->GetListOfKeys()->Print();

TH2D *RT_Phi = (TH2D*)fComp->Get("RT_Phi");
TH2D *ST_Phi = (TH2D*)fComp->Get("ST_Phi");

TH2D *RT_Eta = (TH2D*)fComp->Get("RT_Eta");
TH2D *ST_Eta = (TH2D*)fComp->Get("ST_Eta");

TH2D *RT_Pt = (TH2D*)fComp->Get("RT_Pt");
TH2D *ST_Pt = (TH2D*)fComp->Get("ST_Pt");

TH1D *RT_All_Phi_Eta = (TH1D*)fComp->Get("RT_All_Phi_Eta");
TH1D *ST_All_Phi_Eta = (TH1D*)fComp->Get("ST_All_Phi_Eta");

double cw = 900;
double ch = 600;

void test(){
    fComp->GetListOfKeys()->Print();
}

void heat_maps(){
    TCanvas *c_rt_phi = new TCanvas("c_rt_phi", "c_rt_phi", cw, ch);
    RT_Phi->Draw("COLZ");

    TCanvas *c_st_phi = new TCanvas("c_st_phi", "c_st_phi", cw, ch);
    ST_Phi->Draw("COLZ");
    
    TCanvas *c_rt_eta = new TCanvas("c_rt_eta", "c_rt_eta", cw, ch);
    RT_Eta->Draw("COLZ");

    TCanvas *c_st_eta = new TCanvas("c_st_eta", "c_st_eta", cw, ch);
    ST_Eta->Draw("COLZ");

    TCanvas *c_rt_pt = new TCanvas("c_rt_pt", "c_rt_pt", cw, ch);
    RT_Pt->Draw("COLZ");

    TCanvas *c_st_pt = new TCanvas("c_st_pt", "c_st_pt", cw, ch);
    ST_Pt->Draw("COLZ");
}

void logs(){
    TCanvas *c_rt_all_phi_eta = new TCanvas("c_rt_all_phi_eta", "c_rt_all_phi_eta", cw, ch);
    c_rt_all_phi_eta->SetLogy();
    RT_All_Phi_Eta->Draw();    
    
    TCanvas *c_st_all_phi_eta = new TCanvas("c_st_all_phi_eta", "c_st_all_phi_eta", cw, ch);
    c_st_all_phi_eta->SetLogy();
    ST_All_Phi_Eta->Draw();    
}

