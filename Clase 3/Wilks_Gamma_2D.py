#import libraries
import numpy as np
import matplotlib.pyplot as plt
import os,sys,time
import pickle
import ROOT
from ROOT import TH1D, TRandom3, TCanvas, TF1, gStyle, TFitResultPtr, kRed, kBlue, kGreen, TGraphErrors, gPad, TLine, kDashed, double, gROOT, TMath, TPaveLabel, gMinuit, TTree
from decimal import *
import array


#Test Wilks para Gamma:  LLR = -2*[LL(alfa,beta)-LL(alfaHat,betaHat)] es chisq_2 ?
#No hay formula cerrada para alfaHat, hay que obtener (alfaHat,betaHat) numericamente

n=10
beta=0.5
alfa=3.0
nToys = 10000

#Set style defaults
gStyle.SetOptStat(10);      # only histo name and Nentries
gStyle.SetOptFit(1111);     # display fit results
gStyle.SetPadTickX(1);      # right ticks also
gStyle.SetPadTickY(1);      # upper ticks also
gStyle.SetFuncWidth(3);     # thicker function lines
gStyle.SetHistLineWidth(2); # thickrr histo lines

r = TRandom3(0);          # Initialize random generator
x=array.array('d',np.zeros(n))    # Will store the n random variables

# Log Likelihood LL(alfa,beta), x se pasa a la lambda
LL=(lambda alfa,beta: np.sum([np.log(ROOT.Math.gamma_pdf(i,alfa,1/beta)) for i in x])) 

fGamma = TF1("fGamma","ROOT::Math::gamma_pdf(x,[0],1./[1])",0.001,4*alfa/beta);
fGamma.SetParameter(0,alfa);
fGamma.SetParameter(1,beta);

# Statistics and histos to fill for each experiment
hbetaHat = TH1D("hbetaHat","Estimador de beta",200,0,4*beta);
halfaHat = TH1D("halfaHat","Estimador de alfa",200,0,4*alfa);
hLL0     = TH1D("hLL0","-2log[L(x|alfa,beta)] y -2log[L(x|alfaHat,betaHat)] (azul/rojo)",100,0,n*10);
hLL1     = TH1D("hLL1","-2log[L(x|alfa,beta)] y -2log[L(x|alfaHat,betaHat)] (azul/rojo)",100,0,n*10);
hLLR     = TH1D("hLLR","LLR: -2*log[L(alfa,beta)/L(alfaHat,betaHat)]",200,0,10);
hPvalue  = TH1D("hPvalue","Distribucion of p-value",100,0,1);

c1 = TCanvas("c1");

#####################################
##     Loop over experiments       ##
#####################################

for iexp in range(0,nToys):
    # Generate n gamma random numbers
    fGamma.SetParameter(0,alfa);
    fGamma.SetParameter(1,beta);
    for i in range(0,n): x[i] = fGamma.GetRandom()
    
    #Save them on TTree
    rndm = array.array('d', [0.])
    tree = TTree("tree","tree")
    tree.Branch("rndm",rndm,"rndm/D")
    for i in range(0,n): rndm[0]=x[i]; tree.Fill();

    # Fit Tree to fGamma
    tree.UnbinnedFit("fGamma","rndm","","Q");       # E:Minos Q:quiet
    alfaHat = fGamma.GetParameter(0);
    betaHat = fGamma.GetParameter(1);

    # Wilks: distribucion de -2log[LL(alfa,beta)/LL(alfaHat,betaHat)] es chi2(ndf=2) ?

    LL0 = LL(alfa,beta);
    LL1 = LL(alfaHat,betaHat);
    LLR = (-2*LL0)-(-2*LL1);
    pvalue = ROOT.Math.chisquared_cdf(LLR,2);
    
    hLLR.Fill(LLR);
    hLL0.Fill(-2*LL0);
    hLL1.Fill(-2*LL1);
    hPvalue.Fill(pvalue);
    halfaHat.Fill(alfaHat);
    hbetaHat.Fill(betaHat);
                                                                                   

halfaHat.SetLineColor(kBlue);
if (n>5): halfaHat.Fit("gaus");
halfaHat.Draw();
gPad.Update();
gPad.WaitPrimitive();

hbetaHat.SetLineColor(kBlue);
if (n>5): hbetaHat.Fit("gaus");
hbetaHat.Draw();
gPad.Update();
gPad.WaitPrimitive();

hLL0.SetLineColor(kBlue);
hLL0.Draw("E");
hLL1.SetLineColor(kRed);
hLL1.Draw("E same");
gPad.Update();
gPad.WaitPrimitive();

gPad.SetLogy(0);
hLLR.Scale(1./nToys);
hLLR.Draw("E");
gPad.Update();
gPad.WaitPrimitive();

fchi2 = TF1("chi2","(10./200)*ROOT::Math::chisquared_pdf(x,2)");
fchi2.SetRange(0,10);
fchi2.Draw("same");
gPad.Update();
gPad.WaitPrimitive();

gPad.SetLogy(0);
hPvalue.SetMinimum(0);
hPvalue.Draw("E");
hPvalue.Fit("pol0");
gPad.Update();



errorup=array.array('d', [0.])
errordn=array.array('d', [0.])
errparabolic=array.array('d', [0.])
corr=array.array('d', [0.])
gMinuit.mnerrs(0,errorup,errordn,errparabolic,corr);

print("TTree.UnbinnedFit: {0:.4} {1:.4} \n +/- {2:.4} [ {3:.3}, + {4:.3}] {5:.3}".format(Decimal(alfaHat),Decimal(betaHat),Decimal(fGamma.GetParError(0)),Decimal(errordn[0]),Decimal(errorup[0]),Decimal(errparabolic[0])))
