#import libraries
import numpy as np
import matplotlib.pyplot as plt
import os,sys,time
import pickle
import ROOT
from ROOT import TH1D, TRandom3, TCanvas, TF1, gStyle, TFitResultPtr, kRed, kBlue, kGreen, kGray, kCyan, kMagenta, TGraphErrors, gPad, TLine, kDashed, double, gROOT, TMath, TPaveLabel, gMinuit, TTree, TH1F, TLegend, TFile, TArrow
from decimal import *
import array

# Ajuste no bineado a n numeros aleatorios con distribucion exponencial usando la verosimilitud
# El estimador es siempre el de maxima verosimilitud y se comparan dos formas de estimar la incerteza 
# dejando que un paquete de ROOT use la verosimilitud y haciendolo a mano.
#
#   (1) TTree::UnbinnedFit Likelihood  
#   (2) LL(ci) = LLmax-1/2*n^2 (Likelihood)

#include "TMinuit.h"  // needed for TTree::UnbinnedFit	

n=5
nexperiments=1000
CL=0.9544 # 0.683 0.9544 0.9973
tau=4.0; #exponential distribution parameter


## Cosmetic stuff ##
gStyle.SetOptStat(11);      # only histo name and Nentries
gStyle.SetOptFit(1111);     # display fit results
gStyle.SetPadTickX(1);      # right ticks also
gStyle.SetPadTickY(1);      # upper ticks also
gStyle.SetFuncWidth(3);     # thicker function lines
gStyle.SetHistLineWidth(2); # thickrr histo lines
#####################

r=TRandom3(0)    # Initialize random generator
x=array.array('d',np.zeros(n))    # Will store the n random variables

# Unbinned Log Likelihood LL(tau): 
# x[n] is fixed and passed as reference.
# p[0] is the independent variable. There are no parameters.

LL=(lambda p,aux: np.sum([np.log((1./p[0]) * np.exp(-x_i/p[0])) for x_i in x]))


# TF1 LL: ultimos tres argumentos xmin,xmax,npar
# Aca se define LL(tau), entonces: xmin=taumin, xmax=taumax, npar=0
fLL = TF1("LL",LL,0.001,8*tau,1);
fLL.SetTitle("Log-Likelihood LL(tau|x) for exponential case;tau;LL");

# Funcion exponencial. Usada por TTree::UnbinnedFit
fExpo = TF1("fExpo", "(1/[0])*exp(-x/[0])",0.001,8*tau);

# Función "parabólica" (log(gaussiana)) 
fPar = TF1("fPar", "-(1/2)*((x-[0])/[1])*((x-[0])/[1])",0.001,8*tau);

#Transform input CL => number of sigmas => Delta(logL)
nsigmas = ROOT.Math.gaussian_quantile((1+CL)/2,1)
delta = nsigmas*nsigmas/2;

#/////////////////////////////////////////////
#//                                         //
#//   Comienza el bucle sobre experimentos  //
#//                                         //
#/////////////////////////////////////////////

coverage1 = 0.0
coverage2 = 0.0

for iexp in range(0,nexperiments):

    # Generate n exponential randon numbers
    for i in range(0,n): x[i] = r.Exp(tau);

    # (0) TTreeUnbinnedFit  

    tree = TTree("tree","tree")
    rndm = array.array('d', [0.])
    tree.Branch("rndm",rndm,"rndm/D")

    for i in range(0,n): rndm[0]=x[i]; tree.Fill();

    fExpo.SetParameter(0,4);
    ROOT.Math.MinimizerOptions.SetDefaultErrorDef(delta); # Esta linea es la que hace que ande para mas de un sigma
    tree.UnbinnedFit("fExpo","rndm","","EQ");       # E:Minos Q:quiet

       
    # Ahora le pedimos que calcule los errores por el metodo de la verosimilitud (errorup,errordn)
    # y tambien le pedimos el error usando la aproximacion parabolica (errparabolic)
    errorup=array.array('d', [0.])
    errordn=array.array('d', [0.])
    errparabolic=array.array('d', [0.])
    corr=array.array('d', [0.])
    gMinuit.mnerrs(0,errorup,errordn,errparabolic,corr) 
              
    # (1) Intervalo por el metodo de la verosimilitud
    tau1= fExpo.GetParameter(0); # Tau estimado por MLM;
    ci1_dn = tau1+errordn[0] # Limite inferior del intervalo
    ci1_up = tau1+errorup[0] # Limite superior del intervalo

       	
    # (2) LL(ci) = LLmax-delta (Likelihood) ///////////////////////////
    tau2   = fLL.GetMaximumX();
    LL_max = fLL.Eval(tau2);
    ci2_dn = fLL.GetX(LL_max-delta,0.0001,tau2);
    ci2_up = fLL.GetX(LL_max-delta,tau2,5*tau2);

    print("Resultado utilizando la verosimilitud:   {0:.3} [ {1:.3}, +{2:.3}]\n".format(Decimal(tau1),Decimal(ci1_dn),Decimal(ci1_up)))
    print("Metodo de la verosimilitud hecho a mano:   {0:.3} [ {1:.3}, +{2:.3}]\n".format(Decimal(tau2),Decimal(ci2_dn),Decimal(ci2_up)))
                            
    if (tau > ci1_dn and tau < ci1_up): coverage1+=1;  # Contamos cuantas veces el intervalo incluye el valor real
    if (tau > ci2_dn and tau < ci2_up): coverage2+=1;                                                                                

print("\n Cobertura del intervalo de confianza:\n")                                                                     

print("Resultado utilizando la verosimilitud: {0:.4} %".format(Decimal(coverage1/nexperiments*100)))
print("Metodo de la verosimilitud hecho a mano: {0:.4} % \n".format(Decimal(coverage2/nexperiments*100)))
