import ROOT
import sys
import time
import math
import numpy as np

#---------------------------Useful-Func---------------------------
PI = 3.1415927

def truc_phidist(dist):
    if dist > PI: return truc_phidist(dist - 2*PI)
    elif dist < -PI: return truc_phidist(dist + 2*PI)
    else: return dist

#---------------------------Start-Prosess---------------------------

ifn_full = 'FullSIM_merged.root'
ifn_fast = 'FastSIM_merged.root'
ofn = 'compare_100.root'

print("Opening: " + ifn_full + " and " + ifn_fast)

if_full = ROOT.TFile.Open(ifn_full)
if_fast = ROOT.TFile.Open(ifn_fast)

full_tree = if_full.Get('TrackDensity')
fast_tree = if_fast.Get('TrackDensity')

#otree = ROOT.TTree("TrackDensity", "full vs fast")

print("Tree's loaded")

n_full = full_tree.GetEntries()
n_fast = fast_tree.GetEntries()

print('Fullsim tree has: ' + str(n_full) + ' events')
print('Fastsim tree has: ' + str(n_fast) + ' events')

print('Initializing density histograms')
n_phi = 20
n_eta = 20
# Density histograms on the fly
fast_RT_phiVeta_TH2 = ROOT.TH2D("Fast_RT_phiVeta","Fast_RT_phiVeta", n_phi,-PI,PI, n_eta,-2.4,2.4) 
fast_ST_phiVeta_TH2 = ROOT.TH2D("Fast_ST_phiVeta","Fast_ST_phiVeta", n_phi,-PI,PI, n_eta,-2.4,2.4) 
fast_RT_ptVeta_TH2 = ROOT.TH2D("Fast_RT_ptVeta","Fast_RT_ptVeta", 30,0,15, n_eta,-2.4,2.4) 
fast_ST_ptVeta_TH2 = ROOT.TH2D("Fast_ST_ptVeta","Fast_ST_ptVeta", 30,0,15, n_eta,-2.4,2.4) 

fast_RT_pt_TH1 = ROOT.TH1D("Fast_RT_pt", "Fast_RT_pt", 100,0,100)
fast_ST_pt_TH1 = ROOT.TH1D("Fast_ST_pt", "Fast_ST_pt", 100,0,100)
fast_RT_eta_TH1 = ROOT.TH1D("Fast_RT_eta", "Fast_RT_eta", 100,-2.4,2.4)
fast_ST_eta_TH1 = ROOT.TH1D("Fast_ST_eta", "Fast_ST_eta", 100,-2.4,2.4)

full_RT_phiVeta_TH2 = ROOT.TH2D("Full_RT_phiVeta","Full_RT_phiVeta", n_phi,-PI,PI, n_eta,-2.4,2.4) 
full_ST_phiVeta_TH2 = ROOT.TH2D("Full_ST_phiVeta","Full_ST_phiVeta", n_phi,-PI,PI, n_eta,-2.4,2.4) 
full_RT_ptVeta_TH2 = ROOT.TH2D("Full_RT_ptVeta","Full_RT_ptVeta", 30,0,15, n_eta,-2.4,2.4) 
full_ST_ptVeta_TH2 = ROOT.TH2D("Full_ST_ptVeta","Full_ST_ptVeta", 30,0,15, n_eta,-2.4,2.4) 

full_RT_pt_TH1 = ROOT.TH1D("Full_RT_pt", "Full_RT_pt", 100,0,100)
full_ST_pt_TH1 = ROOT.TH1D("Full_ST_pt", "Full_ST_pt", 100,0,100)
full_RT_eta_TH1 = ROOT.TH1D("Full_RT_eta", "Full_RT_eta", 100,-2.4,2.4)
full_ST_eta_TH1 = ROOT.TH1D("Full_ST_eta", "Full_ST_eta", 100,-2.4,2.4)

temp_full_RT_TD_TH1 = ROOT.TH1D("Temp_Full_RT_TD", "Temp_Full_RT_TD", 100,0,100)
temp_full_ST_TD_TH1 = ROOT.TH1D("Temp_Full_ST_TD", "Temp_Full_ST_TD", 100,0,100)
temp_fast_RT_TD_TH1 = ROOT.TH1D("Temp_Fast_RT_TD", "Temp_Fast_RT_TD", 100,0,100)
temp_fast_ST_TD_TH1 = ROOT.TH1D("Temp_Fast_ST_TD", "Temp_Fast_ST_TD", 100,0,100)

temp_full_RT_TptD_TH1 = ROOT.TH1D("Temp_Full_RT_TptD", "Temp_Full_RT_TptD", 100,0,100)
temp_full_ST_TptD_TH1 = ROOT.TH1D("Temp_Full_ST_TptD", "Temp_Full_ST_TptD", 100,0,100)
temp_fast_RT_TptD_TH1 = ROOT.TH1D("Temp_Fast_RT_TptD", "Temp_Fast_RT_TptD", 100,0,100)
temp_fast_ST_TptD_TH1 = ROOT.TH1D("Temp_Fast_ST_TptD", "Temp_Fast_ST_TptD", 100,0,100)

# Density histograms for keeping
fast_RT_TD_TH1 = ROOT.TH1D("Fast_Bin_RT_TD", "Fast_RT_TD", 20,0,20)
fast_ST_TD_TH1 = ROOT.TH1D("Fast_Bin_ST_TD", "Fast_ST_TD", 20,0,20)
fast_RT_TptD_TH1 = ROOT.TH1D("Fast_Bin_RT_TptD", "Fast_RT_TptD", 20,0,20)
fast_ST_TptD_TH1 = ROOT.TH1D("Fast_Bin_ST_TptD", "Fast_ST_TptD", 20,0,20)

full_RT_TD_TH1 = ROOT.TH1D("Full_Bin_RT_TD", "Full_RT_TD", 20,0,20)
full_ST_TD_TH1 = ROOT.TH1D("Full_Bin_ST_TD", "Full_ST_TD", 20,0,20)
full_RT_TptD_TH1 = ROOT.TH1D("Full_Bin_RT_TptD", "Full_RT_TptD", 20,0,20)
full_ST_TptD_TH1 = ROOT.TH1D("Full_Bin_ST_TptD", "Full_ST_TptD", 20,0,20)

full_RT_meanTD_TH1 = ROOT.TH1D("Full_RT_meanTD", "Full_RT_meanTD", 20,0,20)
full_ST_meanTD_TH1 = ROOT.TH1D("Full_ST_meanTD", "Full_ST_meanTD", 20,0,20)
fast_RT_meanTD_TH1 = ROOT.TH1D("Fast_RT_meanTD", "Fast_RT_meanTD", 20,0,20)
fast_ST_meanTD_TH1 = ROOT.TH1D("Fast_ST_meanTD", "Fast_ST_meanTD", 20,0,20)

full_RT_meanTptD_TH1 = ROOT.TH1D("Full_RT_meanTptD", "Full_RT_meanTptD", 20,0,20)
full_ST_meanTptD_TH1 = ROOT.TH1D("Full_ST_meanTptD", "Full_ST_meanTptD", 20,0,20)
fast_RT_meanTptD_TH1 = ROOT.TH1D("Fast_RT_meanTptD", "Fast_RT_meanTptD", 20,0,20)
fast_ST_meanTptD_TH1 = ROOT.TH1D("Fast_ST_meanTptD", "Fast_ST_meanTptD", 20,0,20)

full_RT_maxTD_TH1 = ROOT.TH1D("Full_RT_maxTD", "Full_RT_maxTD", 50,0,50)
full_ST_maxTD_TH1 = ROOT.TH1D("Full_ST_maxTD", "Full_ST_maxTD", 50,0,50)
fast_RT_maxTD_TH1 = ROOT.TH1D("Fast_RT_maxTD", "Fast_RT_maxTD", 50,0,50)
fast_ST_maxTD_TH1 = ROOT.TH1D("Fast_ST_maxTD", "Fast_ST_maxTD", 50,0,50)

full_RT_maxTptD_TH1 = ROOT.TH1D("Full_RT_maxTptD", "Full_RT_maxTptD", 50,0,50)
full_ST_maxTptD_TH1 = ROOT.TH1D("Full_ST_maxTptD", "Full_ST_maxTptD", 50,0,50)
fast_RT_maxTptD_TH1 = ROOT.TH1D("Fast_RT_maxTptD", "Fast_RT_maxTptD", 50,0,50)
fast_ST_maxTptD_TH1 = ROOT.TH1D("Fast_ST_maxTptD", "Fast_ST_maxTptD", 50,0,50)

fast_RT_meanPtVmaxTD_TH2 = ROOT.TH2D("Fast_RT_meanPtVmaxTD","Fast_RT_meanPtVmaxTD", 20,0,20, 15,0,30)
fast_ST_meanPtVmaxTD_TH2 = ROOT.TH2D("Fast_ST_meanPtVmaxTD","Fast_ST_meanPtVmaxTD", 20,0,20, 15,0,30)
full_RT_meanPtVmaxTD_TH2 = ROOT.TH2D("Full_RT_meanPtVmaxTD","Full_RT_meanPtVmaxTD", 20,0,20, 15,0,30)
full_ST_meanPtVmaxTD_TH2 = ROOT.TH2D("Full_ST_meanPtVmaxTD","Full_ST_meanPtVmaxTD", 20,0,20, 15,0,30)

fast_RT_meanPtVmaxTptD_TH2 = ROOT.TH2D("Fast_RT_meanPtVmaxTptD","Fast_RT_meanPtVmaxTptD", 20,0,20, 15,0,30)
fast_ST_meanPtVmaxTptD_TH2 = ROOT.TH2D("Fast_ST_meanPtVmaxTptD","Fast_ST_meanPtVmaxTptD", 20,0,20, 15,0,30)
full_RT_meanPtVmaxTptD_TH2 = ROOT.TH2D("Full_RT_meanPtVmaxTptD","Full_RT_meanPtVmaxTptD", 20,0,20, 15,0,30)
full_ST_meanPtVmaxTptD_TH2 = ROOT.TH2D("Full_ST_meanPtVmaxTptD","Full_ST_meanPtVmaxTptD", 20,0,20, 15,0,30)

fast_RT_meanPtVmeanTD_TH2 = ROOT.TH2D("Fast_RT_meanPtVmeanTD","Fast_RT_meanPtVmeanTD", 20,0,20, 15,0,15)
fast_ST_meanPtVmeanTD_TH2 = ROOT.TH2D("Fast_ST_meanPtVmeanTD","Fast_ST_meanPtVmeanTD", 20,0,20, 15,0,15)
full_RT_meanPtVmeanTD_TH2 = ROOT.TH2D("Full_RT_meanPtVmeanTD","Full_RT_meanPtVmeanTD", 20,0,20, 15,0,15)
full_ST_meanPtVmeanTD_TH2 = ROOT.TH2D("Full_ST_meanPtVmeanTD","Full_ST_meanPtVmeanTD", 20,0,20, 15,0,15)

fast_RT_meanPtVmeanTptD_TH2 = ROOT.TH2D("Fast_RT_meanPtVmeanTptD","Fast_RT_meanPtVmeanTptD", 20,0,20, 15,0,15)
fast_ST_meanPtVmeanTptD_TH2 = ROOT.TH2D("Fast_ST_meanPtVmeanTptD","Fast_ST_meanPtVmeanTptD", 20,0,20, 15,0,15)
full_RT_meanPtVmeanTptD_TH2 = ROOT.TH2D("Full_RT_meanPtVmeanTptD","Full_RT_meanPtVmeanTptD", 20,0,20, 15,0,15)
full_ST_meanPtVmeanTptD_TH2 = ROOT.TH2D("Full_ST_meanPtVmeanTptD","Full_ST_meanPtVmeanTptD", 20,0,20, 15,0,15)

fast_RT_stdevPtVmaxTD_TH2 = ROOT.TH2D("Fast_RT_stdevPtVmaxTD","Fast_RT_stdevPtVmaxTD", 20,0,20, 25,0,25)
fast_ST_stdevPtVmaxTD_TH2 = ROOT.TH2D("Fast_ST_stdevPtVmaxTD","Fast_ST_stdevPtVmaxTD", 20,0,20, 25,0,25)
full_RT_stdevPtVmaxTD_TH2 = ROOT.TH2D("Full_RT_stdevPtVmaxTD","Full_RT_stdevPtVmaxTD", 20,0,20, 25,0,25)
full_ST_stdevPtVmaxTD_TH2 = ROOT.TH2D("Full_ST_stdevPtVmaxTD","Full_ST_stdevPtVmaxTD", 20,0,20, 25,0,25)

fast_RT_stdevPtVmeanTD_TH2 = ROOT.TH2D("Fast_RT_stdevPtVmeanTD","Fast_RT_stdevPtVmeanTD", 20,0,20, 15,0,15)
fast_ST_stdevPtVmeanTD_TH2 = ROOT.TH2D("Fast_ST_stdevPtVmeanTD","Fast_ST_stdevPtVmeanTD", 20,0,20, 15,0,15)
full_RT_stdevPtVmeanTD_TH2 = ROOT.TH2D("Full_RT_stdevPtVmeanTD","Full_RT_stdevPtVmeanTD", 20,0,20, 15,0,15)
full_ST_stdevPtVmeanTD_TH2 = ROOT.TH2D("Full_ST_stdevPtVmeanTD","Full_ST_stdevPtVmeanTD", 20,0,20, 15,0,15)

fast_RT_meanEtaVmaxTptD_TH2 = ROOT.TH2D("Fast_RT_meanEtaVmaxTptD","Fast_RT_meanEtaVmaxTptD", 20,0,20, 15,0,30)
fast_ST_meanEtaVmaxTptD_TH2 = ROOT.TH2D("Fast_ST_meanEtaVmaxTptD","Fast_ST_meanEtaVmaxTptD", 20,0,20, 15,0,30)
full_RT_meanEtaVmaxTptD_TH2 = ROOT.TH2D("Full_RT_meanEtaVmaxTptD","Full_RT_meanEtaVmaxTptD", 20,0,20, 15,0,30)
full_ST_meanEtaVmaxTptD_TH2 = ROOT.TH2D("Full_ST_meanEtaVmaxTptD","Full_ST_meanEtaVmaxTptD", 20,0,20, 15,0,30)

fast_RT_meanEtaVmeanTptD_TH2 = ROOT.TH2D("Fast_RT_meanEtaVmeanTptD","Fast_RT_meanEtaVmeanTptD", 20,0,20, 15,0,15)
fast_ST_meanEtaVmeanTptD_TH2 = ROOT.TH2D("Fast_ST_meanEtaVmeanTptD","Fast_ST_meanEtaVmeanTptD", 20,0,20, 15,0,15)
full_RT_meanEtaVmeanTptD_TH2 = ROOT.TH2D("Full_RT_meanEtaVmeanTptD","Full_RT_meanEtaVmeanTptD", 20,0,20, 15,0,15)
full_ST_meanEtaVmeanTptD_TH2 = ROOT.TH2D("Full_ST_meanEtaVmeanTptD","Full_ST_meanEtaVmeanTptD", 20,0,20, 15,0,15)

fast_RT_integrated_phiVeta_TH2 = ROOT.TH2D("Fast_RT_integrated_phiVeta","Fast_RT_integrated_phiVeta", n_phi,-PI,PI, n_eta,-2.4,2.4) 
fast_ST_integrated_phiVeta_TH2 = ROOT.TH2D("Fast_ST_integrated_phiVeta","Fast_ST_integrated_phiVeta", n_phi,-PI,PI, n_eta,-2.4,2.4) 
full_RT_integrated_phiVeta_TH2 = ROOT.TH2D("Full_RT_integrated_phiVeta","Full_RT_integrated_phiVeta", n_phi,-PI,PI, n_eta,-2.4,2.4) 
full_ST_integrated_phiVeta_TH2 = ROOT.TH2D("Full_ST_integrated_phiVeta","Full_ST_integrated_phiVeta", n_phi,-PI,PI, n_eta,-2.4,2.4) 

fast_RT_integrated_ptVeta_TH2 = ROOT.TH2D("Fast_RT_integrated_ptVeta","Fast_RT_integrated_ptVeta", 30,0,15, n_eta,-2.4,2.4) 
fast_ST_integrated_ptVeta_TH2 = ROOT.TH2D("Fast_ST_integrated_ptVeta","Fast_ST_integrated_ptVeta", 30,0,15, n_eta,-2.4,2.4) 
full_RT_integrated_ptVeta_TH2 = ROOT.TH2D("Full_RT_integrated_ptVeta","Full_RT_integrated_ptVeta", 30,0,15, n_eta,-2.4,2.4) 
full_ST_integrated_ptVeta_TH2 = ROOT.TH2D("Full_ST_integrated_ptVeta","Full_ST_integrated_ptVeta", 30,0,15, n_eta,-2.4,2.4) 

print('Initializing comparison histograms')
# Comparison histograms (fastTD - fullTD or distance between)
nB_phiVeta_phi = fast_RT_phiVeta_TH2.GetNbinsX() 
nB_phiVeta_eta = fast_RT_phiVeta_TH2.GetNbinsY() 
phi_step = (PI + PI)/(nB_phiVeta_phi + 0.0)
eta_step = (2.4 + 2.4)/(nB_phiVeta_eta + 0.0)

nB_ptVeta_pt = fast_RT_ptVeta_TH2.GetNbinsX() 
nB_ptVeta_eta = fast_RT_ptVeta_TH2.GetNbinsY() 
pt_step = (15.0/(nB_ptVeta_pt + 0.0))
 
diff_up = 30
diff_dn = -30
nB_diff = int(diff_up - diff_dn)

dist_up = 10
dist_dn = 0
nB_dist = 40

dist_phi_up = PI
dist_phi_dn = 0
nB_dist_phi = int(n_phi/2)

ffDiff_RT_TDphi = {}
#for iPhi in range(1,nB_phiVeta_phi+1):
#    phiname = 'RT_Phi_' + str(iPhi)
#    ffDiff_RT_TDphi[phiname] = ROOT.TH1D(phiname, phiname, nB_diff,diff_dn,diff_up) 
    #for iEta in range(1,nB_phiVeta_eta+1):
        #name = 'Phi_' + str(iPhi) + '_Eta_' + str(iEta)
        #ffDiff_RT_TDphi[name] = ROOT.TH1D(name, name, nB_diff,diff_dn,diff_up)
ffDiff_RT_TDphi['RT_All_Phi_Eta'] = ROOT.TH1D('RT_All_Phi_Eta', 'RT_All_Phi_Eta', nB_diff,diff_dn,diff_up) 
ffDiff_RT_TDphi_TH2 = ROOT.TH2D('RT_Phi', 'RT_Phi', nB_diff,diff_dn,diff_up, nB_phiVeta_phi,-PI,PI) 
ffDiff_RT_TDeta_TH2 = ROOT.TH2D('RT_Eta', 'RT_Eta', nB_diff,diff_dn,diff_up, nB_phiVeta_eta,-2.4,2.4) 

ffDiff_RT_TDpt = {}
#for iPt in range(1,nB_ptVeta_pt+1):
#    ptname = 'RT_Pt_' + str(iPt)
#    ffDiff_RT_TDpt[ptname] = ROOT.TH1D(ptname, ptname, nB_diff,diff_dn,diff_up) 
    #for iEta in range(1,nB_ptVeta_eta+1):
        #name = 'Pt_' + str(iPt) + '_Eta_' + str(iEta)
        #ffDiff_RT_TDpt[name] = ROOT.TH1D(name, name, nB_diff,diff_dn,diff_up)
ffDiff_RT_TDpt['RT_All_Pt_Eta'] = ROOT.TH1D('RT_All_Pt_Eta', 'RT_All_Pt_Eta', nB_diff,diff_dn,diff_up) 
ffDiff_RT_TDpt_TH2 = ROOT.TH2D('RT_Pt', 'RT_Pt', nB_diff,diff_dn,diff_up, nB_ptVeta_pt,0,15) 

ffDiff_ST_TDphi = {}
#for iPhi in range(1,nB_phiVeta_phi+1):
#    phiname = 'ST_Phi_' + str(iPhi)
#    ffDiff_ST_TDphi[phiname] = ROOT.TH1D(phiname, phiname, nB_diff,diff_dn,diff_up) 
    #for iEta in range(1,nB_phiVeta_eta+1):
        #name = 'Phi_' + str(iPhi) + '_Eta_' + str(iEta)
        #ffDiff_ST_TDphi[name] = ROOT.TH1D(name, name, nB_diff,diff_dn,diff_up)
ffDiff_ST_TDphi['ST_All_Phi_Eta'] = ROOT.TH1D('ST_All_Phi_Eta', 'ST_All_Phi_Eta', nB_diff,diff_dn,diff_up) 
ffDiff_ST_TDphi_TH2 = ROOT.TH2D('ST_Phi', 'ST_Phi', nB_diff,diff_dn,diff_up, nB_phiVeta_phi,-PI,PI) 
ffDiff_ST_TDeta_TH2 = ROOT.TH2D('ST_Eta', 'ST_Eta', nB_diff,diff_dn,diff_up, nB_phiVeta_eta,-2.4,2.4) 

ffDiff_ST_TDpt = {}
#for iPt in range(1,nB_ptVeta_pt+1):
#    ptname = 'ST_Pt_' + str(iPt)
#    ffDiff_ST_TDpt[ptname] = ROOT.TH1D(ptname, ptname, nB_diff,diff_dn,diff_up) 
    #for iEta in range(1,nB_ptVeta_eta+1):
        #name = 'Pt_' + str(iPt) + '_Eta_' + str(iEta)
        #ffDiff_ST_TDpt[name] = ROOT.TH1D(name, name, nB_diff,diff_dn,diff_up)
ffDiff_ST_TDpt['ST_All_Pt_Eta'] = ROOT.TH1D('ST_All_Pt_Eta', 'ST_All_Pt_Eta', nB_diff,diff_dn,diff_up) 
ffDiff_ST_TDpt_TH2 = ROOT.TH2D('ST_Pt', 'ST_Pt', nB_diff,diff_dn,diff_up, nB_ptVeta_pt,0,15) 

ffDiff_RT_TDmean_TH1 = ROOT.TH1D('RT_mean_TD_diff', 'RT_mean_TD_diff', 20,-5,5)
ffDiff_ST_TDmean_TH1 = ROOT.TH1D('ST_mean_TD_diff', 'ST_mean_TD_diff', 20,-5,5)
ffDiff_RT_TptDmean_TH1 = ROOT.TH1D('RT_mean_TptD_diff', 'RT_mean_TptD_diff', 20,-5,5)
ffDiff_ST_TptDmean_TH1 = ROOT.TH1D('ST_mean_TptD_diff', 'ST_mean_TptD_diff', 20,-5,5)

ffDiff_RT_TDmax_TH1 = ROOT.TH1D('RT_max_TD_diff', 'RT_max_TD_diff', nB_diff,diff_dn,diff_up)
ffDiff_ST_TDmax_TH1 = ROOT.TH1D('ST_max_TD_diff', 'ST_max_TD_diff', nB_diff,diff_dn,diff_up)
ffDiff_RT_TptDmax_TH1 = ROOT.TH1D('RT_max_TptD_diff', 'RT_max_TptD_diff', nB_diff,diff_dn,diff_up)
ffDiff_ST_TptDmax_TH1 = ROOT.TH1D('ST_max_TptD_diff', 'ST_max_TptD_diff', nB_diff,diff_dn,diff_up)

ffDist_RT_TDmax_TH1 = ROOT.TH1D('RT_max_TD_dist', 'RT_max_TD_dist', nB_dist,dist_dn,dist_up)
ffDist_ST_TDmax_TH1 = ROOT.TH1D('ST_max_TD_dist', 'ST_max_TD_dist', nB_dist,dist_dn,dist_up)
ffDistPhi_RT_TDmax_TH1 = ROOT.TH1D('RT_max_TD_distPhi', 'RT_max_TD_distPhi', nB_dist_phi,dist_phi_dn,dist_phi_up)
ffDistPhi_ST_TDmax_TH1 = ROOT.TH1D('ST_max_TD_distPhi', 'ST_max_TD_distPhi', nB_dist_phi,dist_phi_dn,dist_phi_up)
ffDistEta_RT_TDmax_TH1 = ROOT.TH1D('RT_max_TD_distEta', 'RT_max_TD_distEta', nB_dist,dist_dn,dist_up)
ffDistEta_ST_TDmax_TH1 = ROOT.TH1D('ST_max_TD_distEta', 'ST_max_TD_distEta', nB_dist,dist_dn,dist_up)

ffBDistPhi_RT_TDmax_TH1 = ROOT.TH1D('RT_max_TD_BdistPhi', 'RT_max_TD_BdistPhi', (n_phi/2)+1,0,(n_phi/2)+1)
ffBDistPhi_ST_TDmax_TH1 = ROOT.TH1D('ST_max_TD_BdistPhi', 'ST_max_TD_BdistPhi', (n_phi/2)+1,0,(n_phi/2)+1)
ffBDistEta_RT_TDmax_TH1 = ROOT.TH1D('RT_max_TD_BdistEta', 'RT_max_TD_BdistEta', n_eta,0,n_eta)
ffBDistEta_ST_TDmax_TH1 = ROOT.TH1D('ST_max_TD_BdistEta', 'ST_max_TD_BdistEta', n_eta,0,n_eta)

ffDiff_fastloc_RT_TDmax_TH1 = ROOT.TH1D('RT_max_fastloc_TD_diff', 'RT_max_TD_diff', nB_diff,diff_dn,diff_up)
ffDiff_fastloc_ST_TDmax_TH1 = ROOT.TH1D('ST_max_fastloc_TD_diff', 'ST_max_TD_diff', nB_diff,diff_dn,diff_up)
ffDiff_fastloc_RT_TptDmax_TH1 = ROOT.TH1D('RT_max_fastloc_TptD_diff', 'RT_max_TptD_diff', nB_diff,diff_dn,diff_up)
ffDiff_fastloc_ST_TptDmax_TH1 = ROOT.TH1D('ST_max_fastloc_TptD_diff', 'ST_max_TptD_diff', nB_diff,diff_dn,diff_up)

ffDiff_FastMaxMinFullMean_RT_TDmax_TH1 = ROOT.TH1D('RT_max_FastMaxMinFullMean_TD_diff', 'RT_max_TD_diff', nB_diff,diff_dn,diff_up)
ffDiff_FastMaxMinFullMean_ST_TDmax_TH1 = ROOT.TH1D('ST_max_FastMaxMinFullMean_TD_diff', 'ST_max_TD_diff', nB_diff,diff_dn,diff_up)
ffDiff_FastMaxMinFullMean_RT_TptDmax_TH1 = ROOT.TH1D('RT_max_FastMaxMinFullMean_TptD_diff', 'RT_max_TptD_diff', nB_diff,diff_dn,diff_up)
ffDiff_FastMaxMinFullMean_ST_TptDmax_TH1 = ROOT.TH1D('ST_max_FastMaxMinFullMean_TptD_diff', 'ST_max_TptD_diff', nB_diff,diff_dn,diff_up)

ffDiff_FullMaxMinFastMean_RT_TDmax_TH1 = ROOT.TH1D('RT_max_FullMaxMinFastMean_TD_diff', 'RT_max_TD_diff', nB_diff,diff_dn,diff_up)
ffDiff_FullMaxMinFastMean_ST_TDmax_TH1 = ROOT.TH1D('ST_max_FullMaxMinFastMean_TD_diff', 'ST_max_TD_diff', nB_diff,diff_dn,diff_up)
ffDiff_FullMaxMinFastMean_RT_TptDmax_TH1 = ROOT.TH1D('RT_max_FullMaxMinFastMean_TptD_diff', 'RT_max_TptD_diff', nB_diff,diff_dn,diff_up)
ffDiff_FullMaxMinFastMean_ST_TptDmax_TH1 = ROOT.TH1D('ST_max_FullMaxMinFastMean_TptD_diff', 'ST_max_TptD_diff', nB_diff,diff_dn,diff_up)

ffDiff_fullloc_RT_TDmax_TH1 = ROOT.TH1D('RT_max_fullloc_TD_diff', 'RT_max_TD_diff', nB_diff,diff_dn,diff_up)
ffDiff_fullloc_ST_TDmax_TH1 = ROOT.TH1D('ST_max_fullloc_TD_diff', 'ST_max_TD_diff', nB_diff,diff_dn,diff_up)
ffDiff_fullloc_RT_TptDmax_TH1 = ROOT.TH1D('RT_max_fullloc_TptD_diff', 'RT_max_TptD_diff', nB_diff,diff_dn,diff_up)
ffDiff_fullloc_ST_TptDmax_TH1 = ROOT.TH1D('ST_max_fullloc_TptD_diff', 'ST_max_TptD_diff', nB_diff,diff_dn,diff_up)

ffDiff_nRT_TH1 = ROOT.TH1D('nRT_diff', 'nRT_diff', 400,-200,200)
ffDiff_nST_TH1 = ROOT.TH1D('nST_diff', 'nST_diff', 400,-200,200)

#---------------------------Start-Loop---------------------------
print('Loop starting: ')
nEntries = 5000 # n_full
nupd = 40
upd = (nEntries + 0.0)/nupd
npr = min(50,upd)
start = time.time()
for i in range(nEntries):
    full_tree.GetEntry(i)
    fast_tree.GetEntry(i)
    
    if abs(full_tree.event - fast_tree.event) != 0:
        raise IOError('Tree(s) improperly ordered in event: event FullSIM = ' + str(full_tree.event) + ', event FastSIM = ' + str(fast_tree.event) + ', current #event = ' + str(i + 1))

    #print('Number of initial entries ' + str(full_RT_pt_TH1.GetEntries()))

    # Fill density histograms
    fast_RT_pt = [] 
    fast_ST_pt = [] 
    full_RT_pt = [] 
    full_ST_pt = [] 
    
    fast_RT_eta = [] 
    fast_ST_eta = [] 
    full_RT_eta = [] 
    full_ST_eta = [] 
    
    for ifuRT in range(full_tree.nRT):
        full_RT_pt.append(full_tree.RT_pt[ifuRT]) 
        full_RT_eta.append(full_tree.RT_eta[ifuRT]) 
        full_RT_pt_TH1.Fill(full_tree.RT_pt[ifuRT])
        full_RT_eta_TH1.Fill(full_tree.RT_eta[ifuRT])
        full_RT_phiVeta_TH2.Fill(full_tree.RT_phi[ifuRT], full_tree.RT_eta[ifuRT])
        full_RT_ptVeta_TH2.Fill(full_tree.RT_pt[ifuRT], full_tree.RT_eta[ifuRT])
        full_RT_integrated_phiVeta_TH2.Fill(full_tree.RT_phi[ifuRT], full_tree.RT_eta[ifuRT])
        full_RT_integrated_ptVeta_TH2.Fill(full_tree.RT_pt[ifuRT], full_tree.RT_eta[ifuRT])
    for ifuST in range(full_tree.nST):
        full_ST_pt.append(full_tree.ST_pt[ifuST]) 
        full_ST_eta.append(full_tree.ST_eta[ifuST]) 
        full_ST_pt_TH1.Fill(full_tree.ST_pt[ifuST])
        full_ST_eta_TH1.Fill(full_tree.ST_eta[ifuST])
        full_ST_phiVeta_TH2.Fill(full_tree.ST_phi[ifuST], full_tree.ST_eta[ifuST])
        full_ST_ptVeta_TH2.Fill(full_tree.ST_pt[ifuST], full_tree.ST_eta[ifuST])
        full_ST_integrated_phiVeta_TH2.Fill(full_tree.ST_phi[ifuST], full_tree.ST_eta[ifuST])
        full_ST_integrated_ptVeta_TH2.Fill(full_tree.ST_pt[ifuST], full_tree.ST_eta[ifuST])

    for ifaRT in range(fast_tree.nRT):
        fast_RT_pt.append(fast_tree.RT_pt[ifaRT]) 
        fast_RT_eta.append(fast_tree.RT_eta[ifaRT]) 
        fast_RT_pt_TH1.Fill(fast_tree.RT_pt[ifaRT])
        fast_RT_eta_TH1.Fill(fast_tree.RT_eta[ifaRT])
        fast_RT_phiVeta_TH2.Fill(fast_tree.RT_phi[ifaRT], fast_tree.RT_eta[ifaRT])
        fast_RT_ptVeta_TH2.Fill(fast_tree.RT_pt[ifaRT], fast_tree.RT_eta[ifaRT])
        fast_RT_integrated_phiVeta_TH2.Fill(fast_tree.RT_phi[ifaRT], fast_tree.RT_eta[ifaRT])
        fast_RT_integrated_ptVeta_TH2.Fill(fast_tree.RT_pt[ifaRT], fast_tree.RT_eta[ifaRT])
    for ifaST in range(fast_tree.nST):
        fast_ST_pt.append(fast_tree.ST_pt[ifaST]) 
        fast_ST_eta.append(fast_tree.ST_eta[ifaST]) 
        fast_ST_pt_TH1.Fill(fast_tree.ST_pt[ifaST])
        fast_ST_eta_TH1.Fill(fast_tree.ST_eta[ifaST])
        fast_ST_phiVeta_TH2.Fill(fast_tree.ST_phi[ifaST], fast_tree.ST_eta[ifaST])
        fast_ST_ptVeta_TH2.Fill(fast_tree.ST_pt[ifaST], fast_tree.ST_eta[ifaST])
        fast_ST_integrated_phiVeta_TH2.Fill(fast_tree.ST_phi[ifaST], fast_tree.ST_eta[ifaST])
        fast_ST_integrated_ptVeta_TH2.Fill(fast_tree.ST_pt[ifaST], fast_tree.ST_eta[ifaST])

    #print('Number of final entries ' + str(full_RT_pt_TH1.GetEntries()) + ' == ' + str(full_tree.nRT) + '?')
    # Start Comparing
    full_RT_TD_max = 0
    fast_RT_TD_max = 0
    full_RT_TD_max_idxphi = -1
    fast_RT_TD_max_idxphi = -1
    full_RT_TD_max_idxeta = -1
    fast_RT_TD_max_idxeta = -1
    for iPhi in range(1,nB_phiVeta_phi+1):
        for iEta in range(1,nB_phiVeta_eta+1):
            fast_dens = fast_RT_phiVeta_TH2.GetBinContent(iPhi,iEta)
            full_dens = full_RT_phiVeta_TH2.GetBinContent(iPhi,iEta)
            if fast_dens == 0 and full_dens == 0: continue
            if fast_RT_TD_max < fast_dens:
                fast_RT_TD_max = fast_dens
                fast_RT_TD_max_idxphi = iPhi
                fast_RT_TD_max_idxeta = iEta
            if full_RT_TD_max < full_dens:
                full_RT_TD_max = full_dens
                full_RT_TD_max_idxphi = iPhi
                full_RT_TD_max_idxeta = iEta
            diff = fast_RT_phiVeta_TH2.GetBinContent(iPhi,iEta) - full_RT_phiVeta_TH2.GetBinContent(iPhi,iEta)
            ##ffDiff_RT_TDphi['Phi_' + str(iPhi) + '_Eta_' + str(iEta)].Fill(diff)
            #ffDiff_RT_TDphi['RT_Phi_' + str(iPhi)].Fill(diff)
            ffDiff_RT_TDphi_TH2.Fill(diff, -PI + (iPhi-0.5)*phi_step)
            ffDiff_RT_TDeta_TH2.Fill(diff, -2.4 + (iEta-0.5)*eta_step)
            ffDiff_RT_TDphi['RT_All_Phi_Eta'].Fill(diff)
            if fast_dens != 0: 
                fast_RT_TD_TH1.Fill(fast_dens)
                temp_fast_RT_TD_TH1.Fill(fast_dens)
            if full_dens != 0: 
                full_RT_TD_TH1.Fill(full_dens)
                temp_full_RT_TD_TH1.Fill(full_dens)

    full_RT_TptD_max = 0
    fast_RT_TptD_max = 0
    full_RT_TptD_max_idxpt = -1
    fast_RT_TptD_max_idxpt = -1
    full_RT_TptD_max_idxeta = -1
    fast_RT_TptD_max_idxeta = -1
    for iPt in range(1,nB_ptVeta_pt+1):
        for iEta in range(1,nB_ptVeta_eta+1):
            fast_dens = fast_RT_ptVeta_TH2.GetBinContent(iPt,iEta)
            full_dens = full_RT_ptVeta_TH2.GetBinContent(iPt,iEta)
            if fast_dens == 0 and full_dens == 0: continue
            if fast_RT_TptD_max < fast_dens:
                fast_RT_TptD_max = fast_dens
                fast_RT_TptD_max_idxpt = iPt
                fast_RT_TptD_max_idxeta = iEta
            if full_RT_TptD_max < full_dens:
                full_RT_TptD_max = full_dens
                full_RT_TptD_max_idxpt = iPt
                full_RT_TptD_max_idxeta = iEta
            diff = fast_RT_ptVeta_TH2.GetBinContent(iPt,iEta) - full_RT_ptVeta_TH2.GetBinContent(iPt,iEta)
            ##ffDiff_RT_TDpt['Pt_' + str(iPt) + '_Eta_' + str(iEta)].Fill(diff)
            #ffDiff_RT_TDpt['RT_Pt_' + str(iPt)].Fill(diff)
            ffDiff_RT_TDpt_TH2.Fill(diff, 0 + (iPt-0.5)*pt_step)
            ffDiff_RT_TDpt['RT_All_Pt_Eta'].Fill(diff)
            if fast_dens != 0: 
                fast_RT_TptD_TH1.Fill(fast_dens)
                temp_fast_RT_TptD_TH1.Fill(fast_dens)
            if full_dens != 0: 
                full_RT_TptD_TH1.Fill(full_dens)
                temp_full_RT_TptD_TH1.Fill(full_dens)

    full_ST_TD_max = 0
    fast_ST_TD_max = 0
    full_ST_TD_max_idxphi = -1
    fast_ST_TD_max_idxphi = -1
    full_ST_TD_max_idxeta = -1
    fast_ST_TD_max_idxeta = -1
    for iPhi in range(1,nB_phiVeta_phi+1):
        for iEta in range(1,nB_phiVeta_eta+1):
            fast_dens = fast_ST_phiVeta_TH2.GetBinContent(iPhi,iEta)
            full_dens = full_ST_phiVeta_TH2.GetBinContent(iPhi,iEta)
            if fast_dens == 0 and full_dens == 0: continue
            if fast_ST_TD_max < fast_dens:
                fast_ST_TD_max = fast_dens
                fast_ST_TD_max_idxphi = iPhi
                fast_ST_TD_max_idxeta = iEta
            if full_ST_TD_max < full_dens:
                full_ST_TD_max = full_dens
                full_ST_TD_max_idxphi = iPhi
                full_ST_TD_max_idxeta = iEta
            diff = fast_ST_phiVeta_TH2.GetBinContent(iPhi,iEta) - full_ST_phiVeta_TH2.GetBinContent(iPhi,iEta)
            ##ffDiff_ST_TDphi['Phi_' + str(iPhi) + '_Eta_' + str(iEta)].Fill(diff)
            #ffDiff_ST_TDphi['ST_Phi_' + str(iPhi)].Fill(diff)
            ffDiff_ST_TDphi_TH2.Fill(diff, -PI + (iPhi-0.5)*phi_step)
            ffDiff_ST_TDeta_TH2.Fill(diff, -2.4 + (iEta-0.5)*eta_step)
            ffDiff_ST_TDphi['ST_All_Phi_Eta'].Fill(diff)
            if fast_dens != 0:
                fast_ST_TD_TH1.Fill(fast_dens)
                temp_fast_ST_TD_TH1.Fill(fast_dens)
            if full_dens != 0: 
                full_ST_TD_TH1.Fill(full_dens)
                temp_full_ST_TD_TH1.Fill(full_dens)

    full_ST_TptD_max = 0
    fast_ST_TptD_max = 0
    full_ST_TptD_max_idxpt = -1
    fast_ST_TptD_max_idxpt = -1
    full_ST_TptD_max_idxeta = -1
    fast_ST_TptD_max_idxeta = -1
    for iPt in range(1,nB_ptVeta_pt+1):
        for iEta in range(1,nB_ptVeta_eta+1):
            fast_dens = fast_ST_ptVeta_TH2.GetBinContent(iPhi,iEta)
            full_dens = full_ST_ptVeta_TH2.GetBinContent(iPhi,iEta)
            if fast_dens == 0 and full_dens == 0: continue
            if fast_ST_TptD_max < fast_dens:
                fast_ST_TptD_max = fast_dens
                fast_ST_TptD_max_idxpt = iPt
                fast_ST_TptD_max_idxeta = iEta
            if full_ST_TptD_max < full_dens:
                full_ST_TptD_max = full_dens
                full_ST_TptD_max_idxpt = iPt
                full_ST_TptD_max_idxeta = iEta
            diff = fast_ST_ptVeta_TH2.GetBinContent(iPt,iEta) - full_ST_ptVeta_TH2.GetBinContent(iPt,iEta)
            ##ffDiff_ST_TDpt['Pt_' + str(iPt) + '_Eta_' + str(iEta)].Fill(diff)
            #ffDiff_ST_TDpt['ST_Pt_' + str(iPt)].Fill(diff)
            ffDiff_ST_TDpt_TH2.Fill(diff, 0 + (iPt-0.5)*pt_step)
            ffDiff_ST_TDpt['ST_All_Pt_Eta'].Fill(diff)
            if fast_dens != 0: 
                fast_ST_TptD_TH1.Fill(fast_dens)
                temp_fast_ST_TptD_TH1.Fill(fast_dens)
            if full_dens != 0: 
                full_ST_TptD_TH1.Fill(full_dens)
                temp_full_ST_TptD_TH1.Fill(full_dens)
    

    ffDiff_nRT_TH1.Fill(fast_tree.nRT - full_tree.nRT)
    ffDiff_nST_TH1.Fill(fast_tree.nST - full_tree.nST)

    fast_RT_Pt_mean = np.mean(fast_RT_pt)
    full_RT_Pt_mean = np.mean(full_RT_pt)
    fast_ST_Pt_mean = np.mean(fast_ST_pt)
    full_ST_Pt_mean = np.mean(full_ST_pt)
    
    fast_RT_Pt_stdev = np.std(fast_RT_pt)
    full_RT_Pt_stdev = np.std(full_RT_pt)
    fast_ST_Pt_stdev = np.std(fast_ST_pt)
    full_ST_Pt_stdev = np.std(full_ST_pt)
    
    fast_RT_Eta_mean = np.mean(fast_RT_eta)
    full_RT_Eta_mean = np.mean(full_RT_eta)
    fast_ST_Eta_mean = np.mean(fast_ST_eta)
    full_ST_Eta_mean = np.mean(full_ST_eta)
    
    fast_RT_TD_mean = temp_fast_RT_TD_TH1.GetMean()
    full_RT_TD_mean = temp_full_RT_TD_TH1.GetMean()
    fast_ST_TD_mean = temp_fast_ST_TD_TH1.GetMean()
    full_ST_TD_mean = temp_full_ST_TD_TH1.GetMean()

    full_RT_meanTD_TH1.Fill(full_RT_TD_mean)
    full_ST_meanTD_TH1.Fill(full_ST_TD_mean)
    fast_RT_meanTD_TH1.Fill(fast_RT_TD_mean)
    fast_ST_meanTD_TH1.Fill(fast_ST_TD_mean)
    
    full_RT_maxTD_TH1.Fill(full_RT_TD_max)
    full_ST_maxTD_TH1.Fill(full_ST_TD_max)
    fast_RT_maxTD_TH1.Fill(fast_RT_TD_max)
    fast_ST_maxTD_TH1.Fill(fast_ST_TD_max)

    ffDiff_RT_TDmean_TH1.Fill(fast_RT_TD_mean - full_RT_TD_mean)
    ffDiff_ST_TDmean_TH1.Fill(fast_ST_TD_mean - full_ST_TD_mean)
    
    fast_RT_TptD_mean = temp_fast_RT_TptD_TH1.GetMean()
    full_RT_TptD_mean = temp_full_RT_TptD_TH1.GetMean()
    fast_ST_TptD_mean = temp_fast_ST_TptD_TH1.GetMean()
    full_ST_TptD_mean = temp_full_ST_TptD_TH1.GetMean()

    full_RT_meanTptD_TH1.Fill(full_RT_TptD_mean)
    full_ST_meanTptD_TH1.Fill(full_ST_TptD_mean)
    fast_RT_meanTptD_TH1.Fill(fast_RT_TptD_mean)
    fast_ST_meanTptD_TH1.Fill(fast_ST_TptD_mean)
    
    full_RT_maxTptD_TH1.Fill(full_RT_TptD_max)
    full_ST_maxTptD_TH1.Fill(full_ST_TptD_max)
    fast_RT_maxTptD_TH1.Fill(fast_RT_TptD_max)
    fast_ST_maxTptD_TH1.Fill(fast_ST_TptD_max)

    ffDiff_RT_TptDmean_TH1.Fill(fast_RT_TptD_mean - full_RT_TptD_mean)
    ffDiff_ST_TptDmean_TH1.Fill(fast_ST_TptD_mean - full_ST_TptD_mean)

    ffDiff_RT_TDmax_TH1.Fill(fast_RT_TD_max - full_RT_TD_max)
    ffDiff_ST_TDmax_TH1.Fill(fast_ST_TD_max - full_ST_TD_max)
    ffDiff_RT_TptDmax_TH1.Fill(fast_RT_TptD_max - full_RT_TptD_max)
    ffDiff_ST_TptDmax_TH1.Fill(fast_ST_TptD_max - full_ST_TptD_max)
    
    ffDist_RT_TDmax_TH1.Fill(math.sqrt((truc_phidist((fast_RT_TD_max_idxphi - full_RT_TD_max_idxphi)*phi_step))**2 + ((fast_RT_TD_max_idxeta - full_RT_TD_max_idxeta)*eta_step)**2))
    ffDist_ST_TDmax_TH1.Fill(math.sqrt((truc_phidist((fast_ST_TD_max_idxphi - full_ST_TD_max_idxphi)*phi_step))**2 + ((fast_ST_TD_max_idxeta - full_ST_TD_max_idxeta)*eta_step)**2))
    ffDistPhi_RT_TDmax_TH1.Fill(abs(truc_phidist((fast_RT_TD_max_idxphi - full_RT_TD_max_idxphi)*phi_step)))
    ffDistPhi_ST_TDmax_TH1.Fill(abs(truc_phidist((fast_ST_TD_max_idxphi - full_ST_TD_max_idxphi)*phi_step)))
    ffDistEta_RT_TDmax_TH1.Fill(abs((fast_RT_TD_max_idxeta - full_RT_TD_max_idxeta)*eta_step))
    ffDistEta_ST_TDmax_TH1.Fill(abs((fast_ST_TD_max_idxeta - full_ST_TD_max_idxeta)*eta_step))
    
    ffBDistPhi_RT_TDmax_TH1.Fill(min(abs(fast_RT_TD_max_idxphi - full_RT_TD_max_idxphi), n_phi - abs(fast_RT_TD_max_idxphi - full_RT_TD_max_idxphi)))
    ffBDistPhi_ST_TDmax_TH1.Fill(min(abs(fast_ST_TD_max_idxphi - full_ST_TD_max_idxphi), n_phi - abs(fast_ST_TD_max_idxphi - full_ST_TD_max_idxphi)))
    ffBDistEta_RT_TDmax_TH1.Fill(abs(fast_RT_TD_max_idxeta - full_RT_TD_max_idxeta))
    ffBDistEta_ST_TDmax_TH1.Fill(abs(fast_ST_TD_max_idxeta - full_ST_TD_max_idxeta))
    
    ffDiff_fastloc_RT_TDmax_TH1.Fill(fast_RT_TD_max - full_RT_phiVeta_TH2.GetBinContent(fast_ST_TD_max_idxphi, fast_ST_TD_max_idxeta))
    ffDiff_fastloc_ST_TDmax_TH1.Fill(fast_ST_TD_max - full_ST_phiVeta_TH2.GetBinContent(fast_ST_TD_max_idxphi, fast_ST_TD_max_idxeta))
    ffDiff_fastloc_RT_TptDmax_TH1.Fill(fast_RT_TptD_max - full_RT_ptVeta_TH2.GetBinContent(fast_ST_TptD_max_idxpt, fast_ST_TptD_max_idxeta))
    ffDiff_fastloc_ST_TptDmax_TH1.Fill(fast_ST_TptD_max - full_ST_ptVeta_TH2.GetBinContent(fast_ST_TptD_max_idxpt, fast_ST_TptD_max_idxeta))
    
    ffDiff_FastMaxMinFullMean_RT_TDmax_TH1.Fill(fast_RT_TD_max - full_RT_TD_mean)
    ffDiff_FastMaxMinFullMean_ST_TDmax_TH1.Fill(fast_ST_TD_max - full_ST_TD_mean)
    ffDiff_FastMaxMinFullMean_RT_TptDmax_TH1.Fill(fast_RT_TptD_max - full_RT_TptD_mean)
    ffDiff_FastMaxMinFullMean_ST_TptDmax_TH1.Fill(fast_ST_TptD_max - full_ST_TptD_mean)
    
    ffDiff_FullMaxMinFastMean_RT_TDmax_TH1.Fill(fast_RT_TD_mean - full_RT_TD_max)
    ffDiff_FullMaxMinFastMean_ST_TDmax_TH1.Fill(fast_ST_TD_mean - full_ST_TD_max)
    ffDiff_FullMaxMinFastMean_RT_TptDmax_TH1.Fill(fast_RT_TptD_mean - full_RT_TptD_max)
    ffDiff_FullMaxMinFastMean_ST_TptDmax_TH1.Fill(fast_ST_TptD_mean - full_ST_TptD_max)
 
    ffDiff_fullloc_RT_TDmax_TH1.Fill(fast_RT_phiVeta_TH2.GetBinContent(full_ST_TD_max_idxphi, full_ST_TD_max_idxeta) - full_RT_TD_max)
    ffDiff_fullloc_ST_TDmax_TH1.Fill(fast_ST_phiVeta_TH2.GetBinContent(full_ST_TD_max_idxphi, full_ST_TD_max_idxeta) - full_ST_TD_max)
    ffDiff_fullloc_RT_TptDmax_TH1.Fill(fast_RT_ptVeta_TH2.GetBinContent(full_ST_TptD_max_idxpt, full_ST_TptD_max_idxeta) - full_RT_TptD_max)
    ffDiff_fullloc_ST_TptDmax_TH1.Fill(fast_ST_ptVeta_TH2.GetBinContent(full_ST_TptD_max_idxpt, full_ST_TptD_max_idxeta) - full_ST_TptD_max)
    
    fast_RT_meanPtVmaxTD_TH2.Fill(fast_RT_Pt_mean, fast_RT_TD_max)
    fast_ST_meanPtVmaxTD_TH2.Fill(fast_ST_Pt_mean, fast_ST_TD_max)
    full_RT_meanPtVmaxTD_TH2.Fill(full_RT_Pt_mean, full_RT_TD_max)
    full_ST_meanPtVmaxTD_TH2.Fill(full_ST_Pt_mean, full_ST_TD_max)
    
    fast_RT_stdevPtVmaxTD_TH2.Fill(fast_RT_Pt_stdev, fast_RT_TD_max)
    fast_ST_stdevPtVmaxTD_TH2.Fill(fast_ST_Pt_stdev, fast_ST_TD_max)
    full_RT_stdevPtVmaxTD_TH2.Fill(full_RT_Pt_stdev, full_RT_TD_max)
    full_ST_stdevPtVmaxTD_TH2.Fill(full_ST_Pt_stdev, full_ST_TD_max)
    
    fast_RT_meanPtVmeanTD_TH2.Fill(fast_RT_Pt_mean, fast_RT_TD_mean)
    fast_ST_meanPtVmeanTD_TH2.Fill(fast_ST_Pt_mean, fast_ST_TD_mean)
    full_RT_meanPtVmeanTD_TH2.Fill(full_RT_Pt_mean, full_RT_TD_mean)
    full_ST_meanPtVmeanTD_TH2.Fill(full_ST_Pt_mean, full_ST_TD_mean)
    
    fast_RT_meanPtVmaxTptD_TH2.Fill(fast_RT_Pt_mean, fast_RT_TptD_max)
    fast_ST_meanPtVmaxTptD_TH2.Fill(fast_ST_Pt_mean, fast_ST_TptD_max)
    full_RT_meanPtVmaxTptD_TH2.Fill(full_RT_Pt_mean, full_RT_TptD_max)
    full_ST_meanPtVmaxTptD_TH2.Fill(full_ST_Pt_mean, full_ST_TptD_max)
    
    fast_RT_meanEtaVmaxTptD_TH2.Fill(fast_RT_Eta_mean, fast_RT_TptD_max)
    fast_ST_meanEtaVmaxTptD_TH2.Fill(fast_ST_Eta_mean, fast_ST_TptD_max)
    full_RT_meanEtaVmaxTptD_TH2.Fill(full_RT_Eta_mean, full_RT_TptD_max)
    full_ST_meanEtaVmaxTptD_TH2.Fill(full_ST_Eta_mean, full_ST_TptD_max)
    
    fast_RT_meanEtaVmeanTptD_TH2.Fill(fast_RT_Eta_mean, fast_RT_TptD_mean)
    fast_ST_meanEtaVmeanTptD_TH2.Fill(fast_ST_Eta_mean, fast_ST_TptD_mean)
    full_RT_meanEtaVmeanTptD_TH2.Fill(full_RT_Eta_mean, full_RT_TptD_mean)
    full_ST_meanEtaVmeanTptD_TH2.Fill(full_ST_Eta_mean, full_ST_TptD_mean)
    
    # Reset density histograms
    fast_RT_phiVeta_TH2.Reset()
    fast_RT_ptVeta_TH2.Reset()
    fast_RT_pt_TH1.Reset()
    fast_RT_eta_TH1.Reset()
    temp_fast_RT_TD_TH1.Reset()
    temp_fast_RT_TptD_TH1.Reset()

    fast_ST_phiVeta_TH2.Reset()
    fast_ST_ptVeta_TH2.Reset()
    fast_ST_pt_TH1.Reset()
    fast_ST_eta_TH1.Reset()
    temp_fast_ST_TD_TH1.Reset()
    temp_fast_ST_TptD_TH1.Reset()

    full_RT_phiVeta_TH2.Reset()
    full_RT_ptVeta_TH2.Reset()
    full_RT_pt_TH1.Reset()
    full_RT_eta_TH1.Reset()
    temp_full_RT_TD_TH1.Reset()
    temp_full_RT_TptD_TH1.Reset()

    full_ST_phiVeta_TH2.Reset()
    full_ST_ptVeta_TH2.Reset()
    full_ST_pt_TH1.Reset()
    full_ST_eta_TH1.Reset()
    temp_full_ST_TD_TH1.Reset()
    temp_full_ST_TptD_TH1.Reset()

    if i%npr < 1:
       now = time.time()
       evt_str = 'event (%6i/%6i)' % (i, nEntries)
       prog_bar = '[' + ((i*nupd)//nEntries)*'=' + (nupd -((i*nupd)//nEntries))*' ' + ']'
       prc_str = '%4d%%' % ((i + 1.0)*100/nEntries)
       time_str = 'time =  %4ds' % (now - start)
       eta_str = 'eta = %4ds' % (((now - start)/(i + 1))*(nEntries - i - 1))
       sys.stdout.write(' --Progress: ' + evt_str + ', ' + prog_bar + ', ' + prc_str + ', ' + time_str + ', ' + eta_str + ' \r')
       sys.stdout.flush()
now = time.time()
print(' --Done:     event (%6i/%6i), [' % (nEntries, nEntries) + nupd*'=' + '],  %4d%%, time = %4ds, eta = %4ds \r' % (100., now - start, 0.))
print('Loop done')
#---------------------------Finish-Loop---------------------------

print('Writing histograms to: ' + ofn)
ofile = ROOT.TFile(ofn, 'recreate')
for name in ffDiff_RT_TDphi:
    ffDiff_RT_TDphi[name].Write()

for name in ffDiff_RT_TDpt:
    ffDiff_RT_TDpt[name].Write()

for name in ffDiff_ST_TDphi:
    ffDiff_ST_TDphi[name].Write()

for name in ffDiff_ST_TDpt:
    ffDiff_ST_TDpt[name].Write()

ffDiff_RT_TDphi_TH2.Write()
ffDiff_ST_TDphi_TH2.Write()
ffDiff_RT_TDeta_TH2.Write()
ffDiff_ST_TDeta_TH2.Write()
ffDiff_RT_TDpt_TH2.Write()
ffDiff_ST_TDpt_TH2.Write()

ffDiff_RT_TDmean_TH1.Write()
ffDiff_ST_TDmean_TH1.Write()
ffDiff_RT_TptDmean_TH1.Write()
ffDiff_ST_TptDmean_TH1.Write()

ffDiff_RT_TDmax_TH1.Write()
ffDiff_ST_TDmax_TH1.Write()
ffDiff_RT_TptDmax_TH1.Write()
ffDiff_ST_TptDmax_TH1.Write()

ffDist_RT_TDmax_TH1.Write()
ffDist_ST_TDmax_TH1.Write()
ffDistPhi_RT_TDmax_TH1.Write()
ffDistPhi_ST_TDmax_TH1.Write()
ffDistEta_RT_TDmax_TH1.Write()
ffDistEta_ST_TDmax_TH1.Write()

ffBDistPhi_RT_TDmax_TH1.Write()
ffBDistPhi_ST_TDmax_TH1.Write()
ffBDistEta_RT_TDmax_TH1.Write()
ffBDistEta_ST_TDmax_TH1.Write()

ffDiff_fastloc_RT_TDmax_TH1.Write()
ffDiff_fastloc_ST_TDmax_TH1.Write()
ffDiff_fastloc_RT_TptDmax_TH1.Write()
ffDiff_fastloc_ST_TptDmax_TH1.Write()

ffDiff_FastMaxMinFullMean_RT_TDmax_TH1.Write()
ffDiff_FastMaxMinFullMean_ST_TDmax_TH1.Write()
ffDiff_FastMaxMinFullMean_RT_TptDmax_TH1.Write()
ffDiff_FastMaxMinFullMean_ST_TptDmax_TH1.Write()

ffDiff_FullMaxMinFastMean_RT_TDmax_TH1.Write()
ffDiff_FullMaxMinFastMean_ST_TDmax_TH1.Write()
ffDiff_FullMaxMinFastMean_RT_TptDmax_TH1.Write()
ffDiff_FullMaxMinFastMean_ST_TptDmax_TH1.Write()

ffDiff_fullloc_RT_TDmax_TH1.Write()
ffDiff_fullloc_ST_TDmax_TH1.Write()
ffDiff_fullloc_RT_TptDmax_TH1.Write()
ffDiff_fullloc_ST_TptDmax_TH1.Write()

fast_RT_TD_TH1.Write()
fast_ST_TD_TH1.Write()
full_RT_TD_TH1.Write()
full_ST_TD_TH1.Write()

fast_RT_TptD_TH1.Write()
fast_ST_TptD_TH1.Write()
full_RT_TptD_TH1.Write()
full_ST_TptD_TH1.Write()

fast_RT_meanTD_TH1.Write()
fast_ST_meanTD_TH1.Write()
full_RT_meanTD_TH1.Write()
full_ST_meanTD_TH1.Write()

fast_RT_meanTptD_TH1.Write()
fast_ST_meanTptD_TH1.Write()
full_RT_meanTptD_TH1.Write()
full_ST_meanTptD_TH1.Write()

fast_RT_maxTD_TH1.Write()
fast_ST_maxTD_TH1.Write()
full_RT_maxTD_TH1.Write()
full_ST_maxTD_TH1.Write()

fast_RT_maxTptD_TH1.Write()
fast_ST_maxTptD_TH1.Write()
full_RT_maxTptD_TH1.Write()
full_ST_maxTptD_TH1.Write()

full_RT_integrated_phiVeta_TH2.Write()
full_ST_integrated_phiVeta_TH2.Write()
fast_RT_integrated_phiVeta_TH2.Write()
fast_ST_integrated_phiVeta_TH2.Write()

fast_RT_integrated_ptVeta_TH2.Write()
fast_ST_integrated_ptVeta_TH2.Write()
full_RT_integrated_ptVeta_TH2.Write()
full_ST_integrated_ptVeta_TH2.Write()
    
fast_RT_meanPtVmaxTD_TH2.Write()
fast_ST_meanPtVmaxTD_TH2.Write()
full_RT_meanPtVmaxTD_TH2.Write()
full_ST_meanPtVmaxTD_TH2.Write()

fast_RT_stdevPtVmaxTD_TH2.Write()
fast_ST_stdevPtVmaxTD_TH2.Write()
full_RT_stdevPtVmaxTD_TH2.Write()
full_ST_stdevPtVmaxTD_TH2.Write()

fast_RT_meanPtVmeanTD_TH2.Write()
fast_ST_meanPtVmeanTD_TH2.Write()
full_RT_meanPtVmeanTD_TH2.Write()
full_ST_meanPtVmeanTD_TH2.Write()

fast_RT_meanPtVmaxTptD_TH2.Write()
fast_ST_meanPtVmaxTptD_TH2.Write()
full_RT_meanPtVmaxTptD_TH2.Write()
full_ST_meanPtVmaxTptD_TH2.Write()

fast_RT_meanEtaVmaxTptD_TH2.Write()
fast_ST_meanEtaVmaxTptD_TH2.Write()
full_RT_meanEtaVmaxTptD_TH2.Write()
full_ST_meanEtaVmaxTptD_TH2.Write()

fast_RT_meanEtaVmeanTptD_TH2.Write()
fast_ST_meanEtaVmeanTptD_TH2.Write()
full_RT_meanEtaVmeanTptD_TH2.Write()
full_ST_meanEtaVmeanTptD_TH2.Write()

ffDiff_nRT_TH1.Write()
ffDiff_nST_TH1.Write()

print('Closing files')
if_full.Close()
if_fast.Close()
ofile.Close()


