#import libraries
import numpy as np
import matplotlib.pyplot as plt
import os,sys,time
import pickle
import ROOT
from ROOT import TH1D, TRandom3, TCanvas, TF1, gStyle, TFitResultPtr, kRed, kBlue, kGreen, kGray, kCyan, kMagenta, TGraphErrors, gPad, TLine, kDashed, double, gROOT, TMath, TPaveLabel, gMinuit, TTree, TH1F, TLegend, TFile, TArrow, TVirtualFitter, TMatrixDSym, RooEllipse
from decimal import *
import array

# Conjunto de pares de datos con sus errores ###########################

n = 59;

x = array.array('d',[15,30,45,60,75,90,105,120,135,150,165,180,195,210,225,240,255,270,285,300,315,330,345,360,375,390,405,420,435,450,465,480,495,510,525,540,555,570,585,600,615,630,645,660,675,690,705,720,735,750,765,780,795,810,825,840,855,870,885])
y = array.array('d',[775,479,380,302,185,157,137,119,110,89,74,61,66,68,48,54,51,46,55,29,28,37,49,26,35,29,31,24,25,35,24,30,26,28,21,18,20,27,17,17,14,17,24,11,22,17,12,10,13,16, 9, 9,14,21,17,13,12,18,10])
ex = array.array('d',[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]);
ey = array.array('d',[27.8,21.9,19.5,17.4,13.6,12.5,11.7,10.9,10.5,9.4,8.6,7.8,8.1,8.2,6.9,7.3,7.1,6.8,7.4,5.4,5.3,6.1,7.0,5.1,5.9,5.4,5.6,4.9,5.0,5.9,4.9,5.5,5.1,5.3,4.6,4.2,4.5,5.2,4.1,4.1,3.7,4.1,4.9,3.3,4.7,4.1,3.5,3.2,3.6,4,3,3.0,3.7,4.6,4.1,3.6,3.5,4.2,3.2]);
	    	    
################################################
# Set default Minuit in order to use gMinuit
TVirtualFitter.SetDefaultFitter("Minuit");

c1 = TCanvas("c1","Double Exponential Decay",200,10,1500,1000);
c1.Divide(1,2)
c1.cd(1)

gr = TGraphErrors(n,x,y,ex,ey);
   
## Cosmetic stuff ##
gr.SetLineColor(4);
gr.SetLineWidth(2);
gr.SetMarkerColor(4);
gr.SetMarkerSize(0.5);
gr.SetMarkerStyle(20);
gr.SetTitle("Double Exponential Decay;time (sec);Counts");
gr.Draw("AP");

gPad.SetLogy();
gPad.SetGridx();
gPad.SetGridy();

################################################

gPad.Print("Out_Exp_2D.pdf["); # open file
gPad.Print("Out_Exp_2D.pdf"); 

################################################

## Define a double exponential function to be used for fitting

double_exp=(lambda x,p: p[0]+p[1]*np.exp(-x[0]/p[3])+p[2]*np.exp(-x[0]/p[4]))

f1 = TF1("double_exp",double_exp,0,900,5);
# Valores iniciales de los parametros
f1.SetParameters(10,1000,100,40,200); 

# Calculo del mejor conjunto de parametros
gPad.WaitPrimitive()
r = gr.Fit(f1,"SE"); #TFitResultPtr 
sigma = r.GetCovarianceMatrix()
r.Print("V")
gPad.WaitPrimitive()
gPad.Print("Out_Exp_2D.pdf")

# Valores de los parametros obtenidos del ajuste
p0_hat=f1.GetParameter(0); 
p1_hat=f1.GetParameter(1);
p2_hat=f1.GetParameter(2);
p3_hat=f1.GetParameter(3);
p4_hat=f1.GetParameter(4);

# Fijo en sus valores optimos a los parametros que no juegan en las elipses
f1.FixParameter(0, p0_hat); # parametro [0] ahora esta fijo en p0_hat
f1.FixParameter(1, p1_hat); # parametro [1] ahora esta fijo en p1_hat
f1.FixParameter(2, p2_hat); # parametro [2] ahora esta fijo en p2_hat

gr.Fit(f1); # Ajusta de nuevo

################################################

c1.cd(2)

gMinuit.SetErrorDef(1); # 1-sigma
# Hace el contorno de 3 y 4 usando 80 puntos
gr1b = gMinuit.Contour(80,3,4); #TGraph

gMinuit.SetErrorDef(4); # 2-sigma  
# Porque ponemos 4 para dos sigmas aca?
gr2b = gMinuit.Contour(80,3,4); #TGraph

gr2b.SetTitle("68.3% and 95.4% CL regions for parameters tau1 & tau2;tau1;tau2");
gr2b.SetFillColor(42);
gr2b.Draw("Alf"); # Dibuja "elipse" para 2-sigma

gr1b.SetFillColor(38);
gr1b.Draw("Lf"); # Dibuja "elipse" para 1-sigma

gPad.SetLogy(0);
gPad.SetGridx()
gPad.SetGridy()

axis = gr2b.GetXaxis()
axis.SetLimits(28.0,40.0)

gr2b.GetHistogram().SetMinimum(120.); 
gr2b.GetHistogram().SetMaximum(320.); 

gPad.Update();
gPad.Print("Out_Exp_2D.pdf"); 

print("\nPresionar enter para ver las elipses")
sys.stdin.readline()

#### Dibuja las elipses ####

sigma11=np.sqrt(sigma(3,3)); 
sigma22=np.sqrt(sigma(4,4));
cor12=sigma(3,4)/np.sqrt(sigma(3,3)*sigma(4,4));

ell1 = RooEllipse("Parab CI", p3_hat, p4_hat, 1*sigma11, 1*sigma22, cor12, 100);
ell2 = RooEllipse("Parab CI", p3_hat, p4_hat, 2*sigma11, 2*sigma22, cor12, 100);

# Cosmetica  ////////////////////////////
ell1.SetLineColor(2);
ell2.SetLineColor(2);
ell1.SetLineWidth(2);
ell2.SetLineWidth(2);
ell1.SetLineStyle(1);
ell2.SetLineStyle(1);
ell1.SetFillColor(6);
# ///////////////////////////////////////

ell1.Draw("lsame");
ell2.Draw("lsame");

print("\nNo se parecen ... algo va mal, no? modifica minimamente el codigo para arreglarlo! \n \n")
#nombre = input("Presiona una tecla para terminar")  # Asi en python3
nombre = raw_input("Presiona una tecla para terminar...")  # Asi en python2

gPad.Print("Out_Exp_2D.pdf")
gPad.Print("Out_Exp_2D.pdf]");  # cierra el archivo
# Este script crea un pdf de dos paginas
