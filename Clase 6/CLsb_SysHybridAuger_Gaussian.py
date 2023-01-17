#import libraries
import numpy as np
import matplotlib.pyplot as plt
import os,sys,time
import pickle
import ROOT
from ROOT import TH1D, TRandom3, TCanvas, TF1, gStyle, TFitResultPtr, kRed, kBlue, kGreen, kGray, kCyan, kMagenta, TGraphErrors, gPad, TLine, kDashed, double, gROOT, TMath, TPaveLabel, gMinuit, TTree, TH1F, TLegend, TFile, TArrow, TVirtualFitter, TMatrixDSym, TGraph, kBlack, TLatex
from decimal import *
import array


nObs=0
sMax=2.8
bObs=3.5
tau=10.0
nObsEff=100
nTotalEff=100
alpha=0.05

rnd = TRandom3(0);

# eObs=double(nObsEff)/nTotalEff; # Esperanza de la Binomial
eObs = 1.0;

# eErr=np.sqrt(1/nObsEff-1/nTotalEff); # Relative error (sigma/mean) for Binomial given nObsEff and nTotalEff
eErr = 0.1; # Relative error (sigma/mean) for Binomial given nObsEff and nTotalEff

# bErr=np.sqrt(bObs/tau)/(bObs/tau); # Relative error (sigma/mean) for Poisson given bObs and tau
bErr = 0.0;

nToys = 300000;  # number of toys
nscan_points= 100;

g1 = TGraph(nscan_points);
for i in range(0,nscan_points):

    signal=-1+sMax/(nscan_points-1) * i;
    nexp_sb = 0;  # counter for N<=NObs S+B events 

    for iToy in range(0,nToys):
        efi=0.0
   
        # get Efficiency from truncated Gaussian
        while (efi<=0): efi=rnd.Gaus(eObs,eErr*eObs); # Sigma for Binomial is above

        b=0.0;
        # get Background from truncated Gaussian 
        while (b<=0): b=rnd.Gaus(bObs,bObs*bErr); # Sigma for Poisson is sqrt(mean)

        # Generate Nr of events for S+B experiments 
        n_sb = rnd.Poisson(efi*(signal+b));

        # Count N<=nObs events for S+B experiments 
        if(n_sb<=nObs): nexp_sb += 1;
    

    CLsb = double(nexp_sb) / nToys;
    g1.SetPoint(i,CLsb,signal);

print("For Nobs {0:.2}, bkg={1:.2}+/-{2:.2}%, eficiencia={3:.2}+/-{4:.2}\n\nSignal UL = {5:.4} \n".format(Decimal(nObs),Decimal(bObs),Decimal(bErr*100),Decimal(eObs),Decimal(eErr*100),Decimal(g1.Eval(alpha))))


# # Create TCanvas & TFrame, and a TLine at alpha.

# canvas = TCanvas("canvas","canvas",700,500);
# gPad.DrawFrame(0.001,-1,1,sMax-1,"UL with systematics: hybrid;1-CL;Signal UL");
# gPad.SetGridy();
# gPad.SetLogx();
# line = TLine(alpha,-1,alpha,sMax-1);
# line.SetLineColor(kRed);
# line.SetLineWidth(2);
# line.Draw();
# g1.Draw("L");
  