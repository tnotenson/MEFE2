#import libraries
import numpy as np
import matplotlib.pyplot as plt
import os,sys,time
import pickle
import array
from ROOT import gStyle, TRandom3, TH1D, TF1, kRed, kBlue, kGreen, TGraphErrors, TCanvas, gPad
from decimal import *

# Parabola_Fit_Linear(0.05,0.05)
# Subir hasta a=0.10
 
npoints=11
sigma=0.1
a=0.05
xx = array.array('d',[-5,-4,-3,-2,-1,0,+1,+2,+3,+4,+5])
ex = array.array('d',np.zeros(len(xx))) # sin errores en X
ey = array.array('d',np.zeros(len(xx))) # errores y
yyL = array.array('d',np.zeros(len(xx)))
yyP = array.array('d',np.zeros(len(xx)))

b=0.3

r=TRandom3(0)

for i in range(0,npoints):
    x=xx[i]
    ey[i]=sigma
    yyP[i] = r.Gaus(b*x+a*x*x,ey[i]) # data de H1: parabola

grL = TGraphErrors(npoints,xx,yyL,ex,ey);
grP = TGraphErrors(npoints,xx,yyP,ex,ey);
grQ = TGraphErrors(npoints,xx,yyP,ex,ey);

gStyle.SetOptStat(0);
gStyle.SetOptFit(11);

c1 = TCanvas("c1","Lineal vc Parabola",700,1000);
c1.Divide(1,2);
c1.cd(1);

grP.SetMarkerColor(1);
grP.SetMarkerStyle(20);
grP.SetTitle(";X;Y");
grP.Draw("apm");
gPad.SetGridx();
gPad.SetGridy();

xmin = grL.GetXaxis().GetXmin();
xmax = grL.GetXaxis().GetXmax();
fLinear = TF1("fLinear","        [b]*x+[c]",xmin,xmax);
fParabo = TF1("fParabo","[A]*x*x+[B]*x+[C]",xmin,xmax);
fLinear.SetParameters(b,1);
fParabo.SetParameters(0.03,b,1);
fLinear.SetLineColor(kGreen+2);
fParabo.SetLineColor(kBlue);
grP.Fit(fLinear);

c1.cd(2);
grQ.Fit(fParabo);
grQ.SetMarkerColor(1);
grQ.SetMarkerStyle(20);
grQ.SetTitle(";X;Y");
grQ.Draw("apm");
gPad.SetGridx();
gPad.SetGridy();