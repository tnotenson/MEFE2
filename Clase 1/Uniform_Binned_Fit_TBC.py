#import libraries
import numpy as np
import matplotlib.pyplot as plt
import os,sys,time
import pickle
from ROOT import TH1D, TRandom3, TCanvas, TF1, gStyle, TFitResultPtr, kRed, kBlue, kGreen, TGraphErrors, gPad, TLine, kDashed, double, gROOT, TMath

# Ajuste a una funcion constante
# Muestra que Neyman y Pearson estan sesgados con una cobertura pobre
# El ajuste usando la verosimilitud no esta sesgado pero aparentemente tiene 100 por ciento de cobertura.
# Preguntas:
#    1. Por que Neyman subestima y Pearson sobre estima?
#    2. Por que el error en el ajuste con verosimilitud parece estar mal?
#    3. Cambiar el codigo para que el error y la cobertura luzcan correctos.  

# Fill histogram with N random numbers, uniform in (0,100).
N=1000
h1 = TH1D("h1","Constant Distribution",100,0,100)
r=TRandom3(0) #seed=0  ->  different numbers every time

for i in range(0,N): h1.Fill(r.Uniform(0,100))

# If object "canvas" exists, retrieve it in *canvas.
# If not (ie, canvas==0), create it, (for instance in first loop)
# By reusing the *canvas object, focus remains on ROOT terminal.

canvas = gROOT.FindObject("canvas");
if (canvas==None): canvas = TCanvas("canvas","canvas",20,20,1400,1200);
canvas.Clear();

# Plotear el histograma llenado al azar con distribucion uniforme 

canvas.Divide(2,3);
canvas.cd(1);
gStyle.SetOptStat(0); # o no incluye estadistica y 10 solo el numero de eventos
h1.Draw("e");

# Define una funcion constante
gStyle.SetOptFit(0);
f1 = TF1("f1","[0]");

# Choose Minuit2 and setup TFitResultPtr pointers for each fit
#ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2"); 
#TFitResultPtr resultN, resultP, resultL;

# Ajusta el histograma con la funcion f1 usando Neyman, Pearson y la verosimilitud.
# (1) Neyman chi2 (default fit, using observed errors)
f1.SetLineColor(kRed);
resultN = h1.Fit(f1,"S"); 
# (2) P: Use Pearson chi2 (using expected errors instead of observed errors)
f1.SetLineColor(kBlue);
resultP = h1.Fit(f1,"S P+ ");
# (3) L: Use Loglikelihood method (default is chisquare method)
f1.SetLineColor(kGreen+2);
resultL = h1.Fit(f1,"S L+");

# cout<<endl<<"resultN"<<endl; resultN->Print();
# cout<<endl<<"resultP"<<endl; resultP->Print();
# cout<<endl<<"resultL"<<endl; resultL->Print();

# Resume los tres resultados en un TGraph

g1 = TGraphErrors(1);
g2 = TGraphErrors(1);
g3 = TGraphErrors(1);

g1.SetMarkerStyle(20);
g1.SetMarkerColor(kRed);
g1.SetLineColor(kRed);
g1.SetLineWidth(2);
g1.SetPoint(0,1,resultN.Parameter(0));
g1.SetPointError(0,0,resultN.ParError(0));

g2.SetMarkerStyle(20);
g2.SetMarkerColor(kBlue);
g2.SetLineColor(kBlue);
g2.SetLineWidth(2);
g2.SetPoint(0,2,resultP.Parameter(0));
g2.SetPointError(0,0,resultP.ParError(0));

g3.SetMarkerStyle(20);
g3.SetMarkerColor(kGreen+2);
g3.SetLineColor(kGreen+2);
g3.SetLineWidth(2);
g3.SetPoint(0,3,resultL.Parameter(0));
g3.SetPointError(0,0,resultL.ParError(0));

canvas.cd(3);
trueValue = double(N)/h1.GetNbinsX(); # true result of fit
gPad.DrawFrame(0.5,trueValue*0.85,3.5,trueValue*1.1,"Resultado del fit;;Fit Bias");
#https://root.cern.ch/root/htmldoc/guides/users-guide/Histograms.html
g1.Draw("P0 "); 
g2.Draw("P0 ");
g3.Draw("P0 ");

line = TLine(0.5,trueValue,3.5,trueValue);
line.SetLineStyle(kDashed);
line.SetLineColor(1);
line.Draw();
gPad.Update();

# cout << "Neyman     fit bias = " << resultN->Parameter(0)-trueValue << endl;
# cout << "Pearson    fit bias = " << resultP->Parameter(0)-trueValue << endl;
# cout << "Likelihood fit bias = " << resultL->Parameter(0)-trueValue << endl;

#print("(1) <tau_hat>: {0:.4} +- {1:.4}".format(resultN.Parameter(0),resultN.ParError(0)))
#print("(2) <tau_hat>: {0:.4} +- {1:.4}".format(resultP.Parameter(0),resultP.ParError(0)))
#print("(3) <tau_hat>: {0:.4} +- {1:.4}".format(resultL.Parameter(0),resultL.ParError(0)))

# sys.exit(0) # Descomentar para ver plots intermediarios

########################################################################

nToys = 10000;
Uniform = TF1("Uniform","[a]");

coverage1 = 0.0
coverage2 = 0.0
coverage3 = 0.0

# double rangeUp = N/100 * 1.4;
# double rangeDn = N/100 * 0.6;
rangeUp = N/100 + 6* np.sqrt(N)/100;
rangeDn = N/100 - 6* np.sqrt(N)/100;

htau1 = TH1D("htau1","Distribucion de ajustes a la Unforme",100,rangeDn,rangeUp);
htau2 = TH1D("htau2","Distribucion de ajustes a la Unforme",100,rangeDn,rangeUp);
htau3 = TH1D("htau3","Distribucion de ajustes a la Unforme",100,rangeDn,rangeUp);

pvalueN = TH1D("pvalueN","Neyman p-value"       ,100,0,1);
pvalueP = TH1D("pvalueP","Pearson p-value"      ,100,0,1);
pvalueL = TH1D("pvalueL","Baker-Cousins p-value",100,0,1);

for i in range(0,nToys):
    h1.Reset()
    Nentries = N
    #Nentries = r.Poisson(N)
    
    for i in range(0,Nentries): h1.Fill(r.Uniform(0,100))
    Uniform.SetParameter(0,1);

    resultN = h1.Fit(Uniform,"SQN");    # option Q avoids too much printout
    resultP = h1.Fit(Uniform,"SQN P");  # option N avoids adding function histogram
    resultL = h1.Fit(Uniform,"SQN L"); 
    
    const_hat_1 = resultN.Parameter(0); #  +1.11;
    const_err_1 = resultN.ParError(0);
    # Agregar lo necesario para llenar los histogramas con los pvalores

    const_hat_2 = resultP.Parameter(0); # -0.479;
    const_err_2 = resultP.ParError(0);
    # Agregar lo necesario para llenar los histogramas con los pvalores
    
    const_hat_3 = resultL.Parameter(0);    
    const_err_3 = resultL.ParError(0);
    # Agregar lo necesario para llenar los histogramas con los pvalores
    
    print("N {0:.4} {1:.4}".format(const_hat_1,const_err_1))
    print("P {0:.4} {1:.4}".format(const_hat_2,const_err_2))
    print("L {0:.4} {1:.4}\n".format(const_hat_3,const_err_3))

    htau1.Fill(const_hat_1);
    htau2.Fill(const_hat_2);
    htau3.Fill(const_hat_3);

# Agregar acá los contadores coverage1, coverage2 y coverage3
# con la condición que corresponda ...

print("(1) <const_hat>: {0:.4} ".format(htau1.GetMean()))
print("(2) <const_hat>: {0:.4} ".format(htau2.GetMean()))
print("(3) <const_hat>: {0:.4} ".format(htau3.GetMean()))

print("Coverage for the 68% Confidence Intervals:\n")

print("(1) {0:.3} %".format(float(coverage1)/nToys*100))
print("(2) {0:.3} %".format(float(coverage2)/nToys*100))
print("(3) {0:.3} %".format(float(coverage3)/nToys*100))

Rhtau1 = htau1.Fit("gaus","S0");
Rhtau2 = htau2.Fit("gaus","S0");
Rhtau3 = htau3.Fit("gaus","S0");

#Rhtau1.Print();
#Rhtau2.Print();
#Rhtau3.Print();

canvas.cd(5);
htau3.SetLineColor(kGreen+2);
htau3.GetFunction("gaus").SetLineColor(kGreen+2);
htau3.Draw("");
htau2.SetLineColor(kBlue);
htau2.GetFunction("gaus").SetLineColor(kBlue);
htau2.Draw("same");
htau1.SetLineColor(kRed);
htau1.GetFunction("gaus").SetLineColor(kRed);
htau1.Draw("same");

canvas.cd(2);
pvalueN.Draw();
canvas.cd(4);
pvalueP.Draw();
canvas.cd(6);
pvalueL.Draw();

#nombre = input("Presiona una tecla para terminar")  # Asi en python3
nombre = raw_input("Presiona una tecla para terminar...")  # Asi en python2

# Resultados esperados:

#(1) <const_hat>: 8.890
#(2) <const_hat>: 10.479
#(3) <const_hat>: 9.996

# Ayuda:

# chi2_obs = 2.* resultL.MinFcnValue();
# numeros de grados de libertad: resultL.Ndf()
# Integral de la chi2 de Ndf grados de libertad
# desde chi2_obs hasta infinito: TMath.Prob(chi2_obs,ndfL));