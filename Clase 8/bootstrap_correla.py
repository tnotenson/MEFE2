#import libraries
import numpy as np
import matplotlib.pyplot as plt
import os,sys,time
import pickle
import ROOT
from ROOT import TH1D, TRandom3, TCanvas, TF1, gStyle, TFitResultPtr, kRed, kBlue, kGreen, kGray, kCyan, kMagenta, TGraphErrors, gPad, TLine, kDashed, double, gROOT, TMath, TPaveLabel, gMinuit, TTree, TH1F, TLegend, TFile, TArrow, TVirtualFitter, TMatrixDSym, TGraph, kBlack, TLatex,TF2
from decimal import *
import array

# Bootstrap para la correlacion

def correlacion(xdata, ydata):

    xmean = np.mean(xdata)
    #print("x mean: "+str(xmean))
        
    ymean = np.mean(ydata)
    #print("x mean: "+str(ymean))
        
    numerador=0.0;
    for i in range(0,len(xdata)): numerador+=(xdata[i]-xmean)*(ydata[i]-ymean)
        
    sumax=0
    sumay=0
    for i in range(0,len(xdata)): sumax+=np.power(xdata[i]-xmean,2)
    for i in range(0,len(xdata)): sumay+=np.power(ydata[i]-ymean,2)
    denominador=np.power(sumax*sumay,0.5);
	
    corr_xy = numerador/denominador;
    ucorr_xy = (1-np.power(corr_xy,2))/np.power(len(xdata)-3,0.5); # Bajo la aproximacion de que x, y es normal bivariada

    return corr_xy, ucorr_xy

########################################################################
########################################################################

xpoblacion=[];  
ypoblacion=[];  
xsample=[];  
ysample=[];  
bootstrap_sample=[]; 
    
X_Random_variable= TRandom3 (); 

gROOT.SetStyle("Plain");
canvas = TCanvas("canvas","canvas",2000,500);
canvas.Divide(3,1);

########################################################################
########################################################################

# Cargamos los datos de la poblacion
tiempo=[45, 40, 35, 20, 60, 30, 60, 40, 35, 45, 35, 15, 45, 20, 30, 60, 25, 15, 45]; # size = population_size
distancia=[10, 9, 6, 3.1, 12.7, 10, 9, 7.6, 5.8, 16, 4.9, 6, 6.2, 5, 5, 10, 7.5, 2, 5]; # size = population_size

# Elegimos los indices de la poblacion que determinan la muestra
indice=[1, 5, 6, 7, 9, 12, 16]; #size = sample_size

population_size=len(tiempo);
sample_size=len(indice);
replicas=200;

canvas.cd(1);
canvas.cd(1).SetGridx();
canvas.cd(1).SetGridy();

for i in range(0,population_size):
    xpoblacion.append(tiempo[i]);
    ypoblacion.append(distancia[i]);


gr_population = TGraph(len(xpoblacion), array.array('d',xpoblacion), array.array('d',ypoblacion));
gr_population.SetTitle("Poblacion");
gr_population.GetXaxis().SetTitle("Tiempo (minutos)");
gr_population.GetYaxis().SetTitle("Distancia (km)"); 
gr_population.SetMarkerStyle(21);

gr_population.Draw("AP");

corr_xy=0.0;
ucorr_xy=0.0;
corr_xy, ucorr_xy = correlacion(xpoblacion, ypoblacion);
print("Correlacion de la poblacion: "+str(corr_xy)+"  +/- "+str(ucorr_xy))

########################################################################
########################################################################

canvas.cd(2);
canvas.cd(2).SetGridx();
canvas.cd(2).SetGridy();


########################################################################
# Tomo la muestra sin reposicion

for i in range(0,sample_size):
    xsample.append(xpoblacion[indice[i]]);
    ysample.append(ypoblacion[indice[i]]);

gr_sample= TGraph(len(xsample), array.array('d',xsample), array.array('d',ysample))
gr_sample.SetTitle("Muestra");
gr_sample.GetXaxis().SetTitle("Tiempo (minutos)");
gr_sample.GetYaxis().SetTitle("Distancia (km)"); 
gr_sample.SetMarkerStyle(21);

gr_sample.Draw("AP");

corr_xy_sample=0.0;
ucorr_xy_sample=0.0;
corr_xy_sample, ucorr_xy_sample = correlacion(xsample, ysample);
print("Correlacion de la muestra: "+str(corr_xy_sample)+"  +/- "+str(ucorr_xy_sample))

########################################################################
########################################################################
    
histo_bootstrap  =  TH1D("histo_bootstrap","", 1000, 0, 1);
histo_bootstrap.Sumw2();

canvas.cd(3);
canvas.cd(3).SetGridx();
canvas.cd(3).SetGridy();

b_xsample=[];  
b_ysample=[]; 
corr_xy_replica=0.0;
ucorr_xy_replica=0.0;

for i in range(0,replicas):
    b_sample=0.0;
    b_sample_suma=0.0;
    for j in range(0,sample_size):
        indice=int(X_Random_variable.Uniform(0,sample_size))
        b_xsample.append(xsample[indice]);
        b_ysample.append(ysample[indice]);	
    

    corr_xy_replica, ucorr_xy_replica = correlacion(b_xsample, b_ysample);
    
    statistic=corr_xy_replica;
    bootstrap_sample.append(statistic);
    histo_bootstrap.Fill(statistic);


mean = np.mean(bootstrap_sample)

var=0.0;
for i in range(0,replicas):
    var+=np.power(bootstrap_sample[i]-mean,2)/(replicas-1);


ucorr_xy_bootstrap=np.power(var,0.5);

histo_bootstrap.SetTitle("Histograma del rho de las replicas");
histo_bootstrap.GetXaxis().SetTitle("correlacion de la replica");
histo_bootstrap.GetYaxis().SetTitle("#");
histo_bootstrap.SetLineWidth(2);

histo_bootstrap.Draw("HIST");

ucorr_xy_bootstrap2=histo_bootstrap.GetStdDev();
print("Incerteza en la correlacion obtenida por bootstrap: "+str(ucorr_xy_bootstrap))
print("Incerteza en la correlacion obtenida por bootstrap del histograma: "+str(ucorr_xy_bootstrap2))