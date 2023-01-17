#import libraries
import numpy as np
import matplotlib.pyplot as plt
import os,sys,time
import pickle
import ROOT
from ROOT import TH1D, TRandom3, TCanvas, TF1, gStyle, TFitResultPtr, kRed, kBlue, kGreen, TGraphErrors, gPad, TLine, kDashed, double, gROOT, TMath, TPaveLabel, gMinuit, TTree, TH1F, TLegend, TFile
from decimal import *
import array

# ToH: hipotesis simples

# H0: los datos vienen de una dada recta.

# H1: los datos vienen de una dada parabola. 

#

# STD es un test de chi cuadrado de los datos contra cada hipotesis

a=0.03
Nexp=20000

# "a" es el coeficiente cuadratico de la parabola
npoints=11
xx = array.array('d',[-5,-4,-3,-2,-1,0,+1,+2,+3,+4,+5])
ex = array.array('d',np.zeros(len(xx))) # sin errores en X
ey = array.array('d',np.zeros(len(xx))) # errores y
yyL = array.array('d',np.zeros(len(xx)))
yyP = array.array('d',np.zeros(len(xx))) 

b=0.3;

r = TRandom3(0);
for i in range(0,npoints):
   x = xx[i];
   ey[i] = 0.5;
   yyL[i] = r.Gaus (      b*x,ey[i]);  # data de H0: lineal
   yyP[i] = r.Gaus (a*x*x+b*x,ey[i]);  # data de H1: parabola

grL = TGraphErrors(npoints,xx,yyL,ex,ey);
# grP = TGraphErrors(npoints,xx,yyP,ex,ey);

gStyle.SetOptStat(0);
c1 = TCanvas("c1","Lineal vc Parabola",700,1000);
c1.Divide(1,2);
c1.cd(1);

grL.SetMarkerColor(1);
grL.SetMarkerStyle(20);
grL.SetTitle(";X;Y");
grL.Draw("apm");

gPad.SetGrid();

xmin = grL.GetXaxis().GetXmin();
xmax = grL.GetXaxis().GetXmax();
fLinear = TF1("fLinear","[0]*x"        ,xmin,xmax);
fParabo = TF1("fParabo","[0]*x+[1]*x*x",xmin,xmax);
fLinear.SetParameter(0,b);
fParabo.SetParameter(0,b);
fParabo.SetParameter(1,a);
fLinear.SetLineColor(kGreen+2);
fParabo.SetLineColor(kBlue);
fLinear.DrawCopy("same");
fParabo.DrawCopy("same");

c1.Print("Simple_Hipotesis.pdf["); # open file
c1.Print("Simple_Hipotesis.pdf"); # print canvas
gPad.Update();
gPad.WaitPrimitive();

########################################################################
########################################################################
########################################################################

# Lo anterior sirvio para mostrar una realizacion particular de los 11 los puntos y dos modelos
# Ahora se inicia el bucle de Nexp experimentos
# en cada uno de ellos se generan 11 puntos alrededor de la recta 
# y otros 11 alrededor de la parabola
 
hTestSTD_H0 = TH1F("hTestSTD_H0","",200,0,50); # histograma donde voy a guardar los chi2 asumiendo H0 verdadera
hTestSTD_H1 = TH1F("hTestSTD_H1","",200,0,50); # histograma donde voy a guardar los chi2 asumiendo H1 verdadera

for iexp in range(0,Nexp):
    for i in range(0,npoints):
        x = xx[i];
        ey[i] = 0.5;
        yyL[i] = r.Gaus (      b*x,ey[i]);  # data de H0: lineal
        yyP[i] = r.Gaus (a*x*x+b*x,ey[i]);  # data de H1: parabola

    grL2 = TGraphErrors(npoints,xx,yyL,ex,ey);
    grP2 = TGraphErrors(npoints,xx,yyP,ex,ey);

    chi2_H0_Linear = grL2.Chisquare(fLinear);  # Calculo chi2 de los puntos generados a partir de la recta respecto de una recta
    chi2_H1_Linear = grP2.Chisquare(fLinear);  # Calculo chi2 de los puntos generados a partir de la parabola respecto de una recta

    hTestSTD_H0.Fill(chi2_H0_Linear); # voy llenando un histograma con los chi2 que resultan cuando H0 es cierta
    hTestSTD_H1.Fill(chi2_H1_Linear); # voy llenando un histograma con los chi2 que resultan cuando H1 es cierta


hTestSTD_H0.SetLineWidth(2);
hTestSTD_H0.SetLineColor(kBlue);
hTestSTD_H1.SetLineColor(2);
hTestSTD_H1.SetLineWidth(2);
   
###############################
# Calcular Quantiles para STD #
###############################

xq=array.array('d',[0.95]) # Obtener el quantil 0.95
yq=array.array('d',[0.]) # array donde guardarlo

###############################
# Calcular Quantiles para STD #
###############################

h = hTestSTD_H0;
# Calculo el error tipo 1
h.GetQuantiles(1,yq,xq); # Busco el cuantil en la distribucion obtenida bajo H0 verdadera
QuantileSTD = yq[0];
# Integral a izquierda de la distribucion cuando H0 es verdadera
STDerrT1 = h.Integral(h.FindBin(0),h.FindBin(QuantileSTD),""); 

h = hTestSTD_H1;
# Calculo el error tipo 2
# La integral a izquierda de la distribucion cuando H1 es verdadera (usando el mismo cuantil que antes!)
STDerrT2 = h.Integral(h.FindBin(   0),h.FindBin(QuantileSTD),""); 

print("\n STD   Corte: +{0:.3}  CL:{1:.3}   Error T2:{2:.4} \n".format(Decimal(QuantileSTD),Decimal(STDerrT1/Nexp),Decimal(STDerrT2/Nexp)))

#### Dibujar STD ####

XaxisSTD = hTestSTD_H0.GetXaxis();
QuantilBinSTD = XaxisSTD.FindBin(QuantileSTD);
QuantilModSTD = XaxisSTD.GetBinUpEdge(QuantilBinSTD);

#### Crear histogramas de los errores Tipo1 y Tipo2

RejRegionSTD_H0 = hTestSTD_H0.Clone("RejRegionH0");
RejRegionSTD_H0.GetXaxis().SetRangeUser(QuantilModSTD,50);
RejRegionSTD_H0.SetFillStyle(3335);
RejRegionSTD_H0.SetFillColor(4);

RejRegionSTD_H1 = hTestSTD_H1.Clone("RejRegionH1");
RejRegionSTD_H1.GetXaxis().SetRangeUser(0,QuantilModSTD);
RejRegionSTD_H1.SetFillStyle(3353);
RejRegionSTD_H1.SetFillColor(2);

# Dibujar los histogramas de STD

c1.cd(2);
gPad.SetGrid();

hTestSTD_H0.Draw();
c1.Print("Simple_Hipotesis.pdf"); # print canvas
hTestSTD_H1.Draw("same");
c1.Print("Simple_Hipotesis.pdf"); # print canvas
RejRegionSTD_H0.Draw("same");
c1.Print("Simple_Hipotesis.pdf"); # print canvas
RejRegionSTD_H1.Draw("same");
c1.Print("Simple_Hipotesis.pdf"); # print canvas

# Dibujar la leyenda

gStyle.SetLegendBorderSize(0);
gStyle.SetLegendFillColor(0);
gStyle.SetLegendFont(42);
gStyle.SetLegendTextSize(0.05);

legend1 = TLegend(0.7,0.75,0.9,0.89);
legend1.AddEntry(RejRegionSTD_H0,"Test STD|H0","f");
legend1.AddEntry(RejRegionSTD_H1,"Test STD|H1","f");
legend1.Draw();
c1.Print("Simple_Hipotesis.pdf"); # 'print canvas

# Guarda en el pdf y lo cierra.
c1.Print("Simple_Hipotesis.pdf"); # print canvas
c1.Print("Simple_Hipotesis.pdf]"); # close file

# Guarda histogramas para hacer las curvas ROC
f = TFile("Histos.root","recreate");
hTestSTD_H0.Write();
hTestSTD_H1.Write();

#nombre = input("Presiona una tecla para terminar")  # Asi en python3
nombre = raw_input("Presiona una tecla para terminar...")  # Asi en python2

f.Close();
