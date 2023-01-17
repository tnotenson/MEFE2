#import libraries
import numpy as np
import matplotlib.pyplot as plt
import os,sys,time
import pickle
import ROOT
from ROOT import TH1D, TRandom3, TCanvas, TF1, gStyle, TFitResultPtr, kRed, kBlue, kGreen, TGraphErrors, gPad, TLine, kDashed, double, gROOT, TMath, TPaveLabel
import array

# Test Wilks para Gamma, con alfa fijo  LLR = -2*[LL(beta)-LL(betaHat)] es chisq_1 ?
# la formula para beta_hat analitica es betahat = alfa/x_mean

n=5
beta=0.5
alpha=3.0
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

# x[n] is fixed and passed as global argument.
# val[0] is beta, the independent variable, There are no parameters.
LL=(lambda val,par: np.sum([np.log(ROOT.Math.gamma_pdf(i,alpha,1/val[0])) for i in x])) 

# TF1 LL: last 3 arguments xmin,xmax,npar
# Here xmin=betamin, xmax=betamax, npar=0
fLL = TF1("LL",LL,0.001,8*beta,0);
fLL.SetTitle("Log-Likelihood LL(x|#alpha,#beta) for gamma pdf;beta;LL");

# Gamma function, to generate random numbers
fGamma = TF1("fGamma","ROOT::Math::gamma_pdf(x,[0],1./[1])",0.001,4*alpha/beta);
fGamma.SetParameter(0,alpha);
fGamma.SetParameter(1,beta);

# Statistics and histos to fill for each experiment
hbetahat = TH1D("hbeta1"    ,"Distribucion de beta estimator",200,0,4*beta);
hLLR =     TH1D("hLLR"      ,"Distribucion de LLR",200,0,10);
hLLbeta =  TH1D("hLLbeta"   ,"Distribucion de -2logL(beta) y -2logL(betahat) (azul/rojo)",100,0,n*10);
hLLbetahat = TH1D("hLLbetahat","Distribucion de -2logL(beta) y -2logL(betahat) (azul/rojo)",100,0,n*10)
hPvalue = TH1D("hPvalue"   ,"Distribucion de p-value",100,0,1);

c1 = TCanvas("c1");
c1.Draw();

#####################################
##     Loop over experiments       ##
#####################################

for iexp in range(0,nToys):
    # Generate n gamma(alfa,beta) random numbers
    
    for i in range(0,n): x[i] = fGamma.GetRandom()
    # Wilks: distribution of -2log[LL(beta)/LL(betahat)] is chi2(ndf=1) ?
    x_sum = np.sum([i for i in x])
    betahat = alpha * n/x_sum;
    LLbeta = -2*fLL.Eval(beta);
    LLbetahat = -2*fLL.Eval(betahat);
    LLR = LLbeta - LLbetahat;
    pvalue = 1 - ROOT.Math.chisquared_cdf(LLR,1);
    hLLR.Fill(LLR);
    hLLbeta.Fill(LLbeta);
    hLLbetahat.Fill(LLbetahat);
    hPvalue.Fill(pvalue);
    hbetahat.Fill(betahat);
                                                                      

hbetahat.SetLineColor(kBlue);
if (n>5): hbetahat.Fit("gaus");
hbetahat.Draw();
gPad.Update();
gPad.WaitPrimitive();

hLLbeta.Draw("E");

hLLbetahat.SetLineColor(kRed);
hLLbetahat.Draw("E same");
gPad.Update();
gPad.WaitPrimitive();

gPad.SetLogy(1);
hLLR.Scale(1./nToys);
hLLR.Draw("E");
gPad.Update();
gPad.WaitPrimitive();

fchi2 = TF1("chi2","(10./200)*ROOT::Math::chisquared_pdf(x,1)");
fchi2.SetRange(0,10);
fchi2.Draw("same");
gPad.Update();
gPad.WaitPrimitive();

gPad.SetLogy(0);
hPvalue.SetMinimum(0);
hPvalue.Draw("E");
hPvalue.Fit("pol0");
gPad.Update();
gPad.WaitPrimitive();



# # La distribucion gamma es:
# # f(x|alfa,beta) = [beta^alfa/Gamma(alfa)] x^(alfa-1) e^(-beta*x)

# # donde la funcion Gamma(x) en root se calcula como ROOT::Math::tgamma(x)

# # Tomemos como alfa=3 y beta=0.5

# # Sea X una VA con distribucion f(x|alfa,beta), y defino una nueva variable Y=-2log[f(X|alfa,beta)],
# # la distribucion de Y es extranha, no baja de 4 y diverge en 4 viniendo de la derecha, hLLbeta

# # Esto se entiende analizando al funcion -2log[f(X|alfa,beta):

# f = TF1("f","-2*([a]*log([b])-log(ROOT::Math::tgamma([a]))+([a]-1)*log(x)-[b]*x)",0.1,10)
# f.SetParameters(3,0.5);
# f.Draw();

# # donde se ve que a bajo x domina ([a]-1)*log(x) y a alto x domina -[b]*x



# Transform input CL => number of sigmas => Delta(logL)
# nsigmas = ROOT.Math.gaussian_quantile((1+CL)/2,1);
# delta = nsigmas*nsigmas/2;  
#
# To Draw logL for each experiment
# fLLDraw = TF1("log L(beta|xi)",LL,0.001,8*beta,0);
# fLLDraw.SetTitle(";Parameter #beta;log L(#beta|x)");

