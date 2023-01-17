#import libraries
import numpy as np
import matplotlib.pyplot as plt
import os,sys,time
import pickle
import ROOT
from ROOT import TH1D, TRandom3, TCanvas, TF1, gStyle, TFitResultPtr, kRed, kBlue, kGreen, kGray, kCyan, kMagenta, TGraphErrors, gPad, TLine, kDashed, double, gROOT, TMath, TPaveLabel, gMinuit, TTree, TH1F, TLegend, TFile, TArrow, TVirtualFitter, TMatrixDSym, TGraph, kBlack, TLatex,TF2
from decimal import *
import array

# Este Scripts compara la performance de Bootstrap contra Toys Montecarlo
# Ayuda a visualizar que tamanio de muestra y numero de replicas son suficientes
# para obtener una estimacion razonable de la incerteza en el promedio de una muestra Poissoniana

sample_size=20;
replicas=200;

mu=50.0;
sigma = np.sqrt(50.0);
sigmam = sigma/np.sqrt(sample_size);
emin=mu-3*sigma;
emax=mu+3*sigma;
nsigmas=10.0;
	
bines = int(emax-emin);

population=[];  
sample=[];  
bootstrap_sample=[]; 
	
X_Random_variable = TRandom3(0); 

###############################################################################################
###############################################################################################

gROOT.SetStyle("Plain");
canvas = TCanvas("canvas","canvas",2000,1000);
canvas.Divide(2,2);

###############################################################################################
###############################################################################################

canvas.cd(1);
canvas.cd(1).SetGridx();
canvas.cd(1).SetGridy();

pdf = TF1("pdf","TMath::Poisson(x,50)",emin,emax);
pdf.SetLineColor(kRed);
pdf.Draw();


############################################################
#### Aca empiezo los Toys
	
histo_toy  =  TH1D("histo_toy","", bines, mu-nsigmas*sigmam, mu+nsigmas*sigmam);
histo_toy.Sumw2();

canvas.cd(2);
canvas.cd(2).SetGridx();
canvas.cd(2).SetGridy();

################################################################
# Hago los toys   

toys_mean=0.0;
for i in range(0,replicas+1):
    t_sample=0.0
    t_sample_suma=0.0

    for j in range(1,sample_size+1):
        t_sample = X_Random_variable.Poisson(mu);
        t_sample_suma+= t_sample;

    mean=double(t_sample_suma)/sample_size;
    t_sample_suma=0;
    statistic=mean;
    toys_mean+=mean;
    histo_toy.Fill(statistic);

toys_mean=toys_mean/replicas;
print("El promedio de los toys es:   "+str(toys_mean))

histo_toy.SetTitle("Toys del promedio");
histo_toy.GetXaxis().SetTitle("promedio de la replica");
histo_toy.GetYaxis().SetTitle("#");
histo_toy.SetLineWidth(2);

histo_toy.SetLineColor(kRed);
histo_toy.Draw("E0 HIST");

###############################################################
# Aca empiezo el bootstrap

###############################################################
###############################################################

histo_sample  =  TH1D("histo_sample","", bines, emin, emax);
histo_sample.Sumw2();

canvas.cd(3);
canvas.cd(3).SetGridx();
canvas.cd(3).SetGridy();


###############################################################
# Tomo la muestra 

x_sample_mean=0.0;
for i in range(1,sample_size+1):
	x_sample=0.0;
	#x_sample=X_Random_variable.Uniform(emin,emax);
	x_sample=X_Random_variable.Poisson(mu);
	x_sample_mean+=x_sample;
	sample.append(x_sample);
	histo_sample.Fill(x_sample);

x_sample_mean=double(x_sample_mean)/sample_size;
print("El promedio de la muestra es: "+str(x_sample_mean))

histo_sample.SetTitle("Muestra");
histo_sample.GetXaxis().SetTitle("realizacion");
histo_sample.GetYaxis().SetTitle("#");
histo_sample.SetLineWidth(2);

histo_sample.SetLineColor(kBlue);
histo_sample.Draw("E0 HIST");

###############################################################
###############################################################
	
histo_bootstrap  =  TH1D("histo_bootstrap","", bines, mu-nsigmas*sigmam, mu+nsigmas*sigmam);
histo_bootstrap.Sumw2();

canvas.cd(4);
canvas.cd(4).SetGridx();
canvas.cd(4).SetGridy();

###############################################################
# Hago las replicas, estas son con reposicion

bootstrap_mean=0.0
for i in range(1,replicas+1):
    b_sample=0.0;
    b_sample_suma=0.0;
    for j in range(1,sample_size+1):
        indice=int(X_Random_variable.Uniform(1,sample_size))
        b_sample=sample[indice];
        b_sample_suma+=b_sample;

    mean=b_sample_suma/sample_size;
    b_sample_suma=0.0;
    statistic=mean;
    bootstrap_mean+=mean;
    # print("Estadistico: "+str(statistic))
    bootstrap_sample.append(statistic);
    histo_bootstrap.Fill(statistic);

bootstrap_mean=bootstrap_mean/replicas;
print("El promedio de las replicas es:   "+str(bootstrap_mean))
bias=abs(bootstrap_mean-x_sample_mean);
print("por lo tanto el bias es:          "+str(bias))

histo_bootstrap.SetTitle("Bootstrap para el promedio");
histo_bootstrap.GetXaxis().SetTitle("promedio de la replica");
histo_bootstrap.GetYaxis().SetTitle("#");
histo_bootstrap.SetLineWidth(2);

histo_bootstrap.SetLineColor(kBlue);
histo_bootstrap.Draw("E0 HIST");

#nombre = input("Presiona una tecla para terminar")  # Asi en python3
nombre = raw_input("Presiona una tecla para terminar...")  # Asi en python2