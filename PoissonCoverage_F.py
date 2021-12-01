#import libraries
import numpy as np
import matplotlib.pyplot as plt
import os,sys,time
import pickle
import ROOT
from ROOT import TH1D, TRandom3, TCanvas, TF1, gStyle, TFitResultPtr, kRed, kBlue, kGreen, kGray, kCyan, kMagenta, TGraphErrors, gPad, TLine, kDashed, double, gROOT, TMath, TPaveLabel, gMinuit, TTree, TH1F, TLegend, TFile, TArrow, TVirtualFitter, TMatrixDSym, TGraph, kBlack, TLatex
from decimal import *
import array
from scipy.stats import poisson
from scipy.integrate import quad
# Este script calcula el intervalo frecuentista para el mu de una variable aleatoria de Poisson.

mu_min = 0.0001
mu_max = 6.5
CL = 0.6827
nscan_points= 100
CoverageMin = 0.0

# Transformamos el CL en una numero de sigmas de la Gaussiana N(0, 1)
# esto es de utilidad para el calculo de intervalo usando la varianza y la verosimilitud (ambos pendientes de implementar)

# double nsigmas = ROOT::Math::gaussian_quantile((1+CL)/2,1);
# double delta = nsigmas*nsigmas/2; 

# (1) Intervalo usando la varianza:  nobs +- sqrt(nobs) /////////////////////
# A completar!
def IsInside1(nobs,mu):
    sigma = np.sqrt(nobs)
    xmin = nobs - sigma
    xmax = nobs + sigma
    return (xmin <= mu and mu <= xmax)

# (2) Intervalo usando la verosimilitud: LL(ci_limits) = LLmax-1/2 ///////////
# A completar!
    
def LL(nobs,mu):
    if mu>0:
        return -mu+nobs*np.log(mu)
    else:
        return 0
     

def LLvec(nobs,mus):
    if len(mus)>0:
        return -mus+nobs*np.log(mus)
    else:
        return np.zeros((len(mus)))

#def fLL(x,par):
    #if par[0] > 0:
    #return -par[0]+x*np.log(par[0])
    #else:
        #return 0

   

def IsInside2(nobs,mu):
    ########## python ########################################################
    # LLmax = LL(nobs,nobs)
    # dom = np.linspace(0.001, mu+10*np.sqrt(mu),1000)
    # img = LLvec(nobs,dom); y = LLmax - 1/2
    # arr1 = np.where(img >= y)
    # xmin = dom[arr1[0][0]]; xmax = dom[arr1[0][-1]]
    # print('python puro',xmin, xmax)
    
    ######### pyROOT #########################################################
    delta = 1/2
    rootLL = TF1("rootLL", "-x+[n]*log(x)",0.0001,mu+10*np.sqrt(mu),1);
    rootLL.SetParameter(0,nobs)
    muMLE   = rootLL.GetMaximumX();
    # print('nobs',muMLE,nobs)
    LL_max = rootLL.Eval(muMLE);
    if nobs==0:
        # LL = -x; LLmax = 0 porque x>0; LLmax - 1/2 = -1/2
        xmin = 0 #     porque mu>0
        xmax = delta # porque LLmax-1/2 = = -1/2 = -x
    else:
        xmin = rootLL.GetX(LL_max-delta,0.0001,muMLE);
        xmax = rootLL.GetX(LL_max-delta,muMLE,5*muMLE);
    # print('pyROOT', xmin, xmax)
    
    return (xmin <= mu and mu <= xmax)


# (3) Intervalo frecuentista (analitico) ///////////////////////////////////
def IsInside3(nobs,mu):
    alpha = 1- CL;
    xmin = 0.5*ROOT.Math.chisquared_quantile( alpha/2, 2*nobs) if (nobs>0) else 0
    xmax = 0.5*ROOT.Math.chisquared_quantile_c( alpha/2, 2*(nobs+1));
    return (xmin <= mu and mu <= xmax)

def poisson_bayes(x, nobs):
    return poisson.pmf(nobs, x)

# (4) Intervalo bayesiano
def IsInside4(nobs,mu):
    alpha = (1 - CL)/2;
    
    if nobs == 0: # Calculo exactamente el caso nobs = 0
        xmin = 0
        xmax = -np.log(1-CL)
        return (xmin <= mu and mu <= xmax)
    
    ########## python ## IN PROCESS ############################################
    #x = np.linspace(mu_min, mu_max, nscan_points)
    #N = quad(poisson_bayes, 0, np.inf, args=(nobs,))[0]
    

    
    #cdf_bayes = 0
    #i = 0
    #while cdf_bayes <= alpha and i<nscan_points:
    #    print('while 1',i)
    #    m = x[i]
    #    cdf_bayes = quad(poisson_bayes, 0, m, args=(nobs,))[0]/N
    #    i += 1
    #if i == len(x):i -= 1
    #xmin = x[i]
    #
    #cdf_bayes = 0
    #while cdf_bayes <= alpha and i<nscan_points:
    #    print('while 2',i)
    #    m = x[i]
    #    cdf_bayes = quad(poisson_bayes, m, np.inf, args=(nobs,))[0]/N
    #    i += 1
    #if i == len(x):i -= 1
    #xmax = x[i]
    
    ########## pyROOT #########################################################
    fP = TF1("fP", "[0]*TMath::PoissonI([1],x)",0.0001,mu+10*np.sqrt(mu),1);
    fP.SetParameter(0,1) # lo fijo en 1 para integrar y luego normalizo
    fP.SetParameter(1,nobs) 
    N = fP.Integral(0,np.inf)
    fP.SetParameter(0,1/N) # normalizo
    
    xq=array.array('d',[alpha, 1-alpha]) # array con las prob cuantiles
    yq=array.array('d',[0.,0.]) # array donde guardar los cuantiles
    
    fP.GetQuantiles(2,yq,xq)
    print(mu, yq)
    xmin, xmax = yq
    
    return (xmin <= mu and mu <= xmax)
    
###########################################################################
# Cobertura vs mu a ser guardadas en cuatro objetos de la clase TGraphs
g1 = TGraph(nscan_points);
g2 = TGraph(nscan_points);
g3 = TGraph(nscan_points);
g4 = TGraph(nscan_points);

# Bucle sobre mu, desde mu_min a mu_max.
# Para cada mu calcula la cobertura del intervalo y lo guarda en el TGraphs.
for i in range(0,nscan_points):
    mu = mu_min + (mu_max-mu_min)/(nscan_points-1) * i;
    #print(mu)
   
    Nmax = int(mu+10*np.sqrt(mu))+1 # esperanza mas 10 veces sigma

    # Inicializar las variables que faltan
    probInside1 = 0; probInside2 = 0; probInside3 = 0; probInside4 = 0

    # Barro desde n observado hasta la esperanza mas 10 veces sigma    
    for nobs in range(0,Nmax):
        # print(mu)
        inside1 = IsInside1(nobs,mu); inside2 = IsInside2(nobs,mu); 
        inside3 = IsInside3(nobs,mu); inside4 = IsInside4(nobs,mu);
        prob = ROOT.Math.poisson_pdf(nobs,mu);
        # Completar con lo correspondiente a los otros dos intervalos

    # Si nobs esta dentro del intervalo, agrega la probabilidd de ese caso particular, i.e. Poisson(nobs,mu).
        if (inside1): probInside1 += prob; if (inside2): probInside2 += prob;
        if (inside3): probInside3 += prob; if (inside4): probInside4 += prob;

    Offset = 0.003; # Vertical offset between plots to avoid overlap   
    g1.SetPoint(i,mu,probInside1); g2.SetPoint(i,mu,probInside2);
    g3.SetPoint(i,mu,probInside3); g4.SetPoint(i,mu,probInside4);
   

########################################################################
# Create TCanvas & TFrame, and a TLine at CL for reference 
gStyle.SetOptStat(0);

canvas = TCanvas("canvas","canvas",900,700);
canvas.Divide(2,2, 0.01, 0.01);

gPad.DrawFrame(mu_min,CoverageMin,mu_max,1,"Comparacion de los Intervalos de Confianza para Poisson;\\text{Parametro de Poisson }\\mu;Cobertura");

canvas.cd(1);

l = TLine(mu_min,CL,mu_max,CL);
l.SetLineStyle(kDashed);
l.Draw(); # Draw a reference line at y-axis = CL

# The 4 TGraphs can be plotted together (I added a little offset between them),
# but first better look at each of them one at a time. Uncomment as necessary:
 
# Naive Poisson interval (Red) 
g1.SetLineWidth(2);
g1.SetLineColor(kRed);
g1.Draw("L");
gPad.Update();
gPad.WaitPrimitive();

canvas.cd(2);

l = TLine(mu_min,CL,mu_max,CL);
l.SetLineStyle(kDashed);
l.Draw(); # Draw a reference line at y-axis = CL

# LogLikelihood interval (Blue line) 
g2.SetLineWidth(2);
g2.SetLineColor(kBlue);
g2.Draw("L");
gPad.Update();
gPad.WaitPrimitive();
  
canvas.cd(3);

l = TLine(mu_min,CL,mu_max,CL);
l.SetLineStyle(kDashed);
l.Draw(); # Draw a reference line at y-axis = CL

# Analytic frequentist interval (Black line) 
g3.SetLineWidth(2);
g3.SetLineColor(kBlack);
g3.Draw("L");
gPad.Update();
gPad.WaitPrimitive();

canvas.cd(4);

l = TLine(mu_min,CL,mu_max,CL);
l.SetLineStyle(kDashed);
l.Draw(); # Draw a reference line at y-axis = CL

# Bayesian interval (Green line)
g4.SetLineWidth(2);
g4.SetLineColor(kGreen);
g4.Draw("L");
gPad.Update();
gPad.WaitPrimitive();

legen = TLegend();
legen.AddEntry(g1,"Aprox gaussiana");
legen.AddEntry(g2,"Log Likelihood");
legen.AddEntry(g3,"Frecuentista");
legen.AddEntry(g4,"Bayesiano");
legen.Draw("same");
gPad.Update();

nombre = input("Presiona una tecla para terminar")  # Asi en python3
# nombre = raw_input("Presiona una tecla para terminar...")  # Asi en python2
