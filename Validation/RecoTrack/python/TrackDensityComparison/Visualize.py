import ROOT

fComp = ROOT.TFile.Open("compare_100.root")
fFull = ROOT.TFile.Open("FullSIM_merged.root")
fFast = ROOT.TFile.Open("FastSIM_merged.root")

tFull = fFull.Get("TrackDensity")
tFast = fFast.Get("TrackDensity")

cw = 900
ch = 600

keys = fComp.GetListOfKeys()
th2 = []
th1 = []
for key in keys:
    if 'TH2' in key.GetClassName():
        th2.append(key.GetName())
    if 'TH1' in key.GetClassName():
        th1.append(key.GetName())
canv_l = []

def nRT_nST():
    ROOT.gStyle.SetOptStat(0)
    c_RT = ROOT.TCanvas('c_RT', 'c_RT', cw, ch)
    c_ST = ROOT.TCanvas('c_ST', 'c_ST', cw, ch)

    full_proxy = ROOT.TH1D('ST', 'ST', 600,0,300)  
    fast_proxy = ROOT.TH1D('RT', 'RT', 600,0,300)  
    #for i in range(100):
    #    full_proxy.Fill(200)
    #    fast_proxy.Fill(200)
    full_proxy.SetLineColor(1)
    fast_proxy.SetLineColor(2)

    legend = ROOT.TLegend(0.62, 0.8, 0.9, 0.9)
    #legend = ROOT.TLegend(0.52, 0.7, 0.9, 0.9)
    #legend.SetHeader('Title', 'C')
    legend.AddEntry(full_proxy, 'Full SIM', 'L')
    legend.AddEntry(fast_proxy, 'Fast SIM', 'L')


    tFull.SetLineColor(1)
    tFast.SetLineColor(2)
 
    c_RT.cd() 
    tFull.Draw('nRT')
    tFast.Draw('nRT', '', 'same')
    full_proxy.Draw('same')
    fast_proxy.Draw('same')
    legend.Draw()
    
    c_ST.cd() 
    tFull.Draw('nST')#>>hsqrt(300,0,300)
    tFast.Draw('nST', '', 'same')
    full_proxy.Draw('same')
    fast_proxy.Draw('same')
    legend.Draw()
    print(type(legend))

def heat_map():
    for hist_n in th2:
        canv_n = 'c_' + hist_n
        TH2 = fComp.Get(hist_n)
        c = ROOT.TCanvas(canv_n, canv_n, cw, ch)
        canv_l.append(c)
        TH2.Draw('COLZ')


if __name__ == '__main__':
    nRT_nST()
    heat_map()
    print('Press enter to exit.')
    text = raw_input()

'''
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
'''

'''
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
'''
