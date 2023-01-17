#import libraries
import numpy as np
import matplotlib.pyplot as plt
import os,sys,time
import pickle
from ROOT import TH1D, TRandom3, TCanvas, TF1, gStyle, TFitResultPtr, kRed, kBlue, kGreen, TGraphErrors, gPad, TLine, kDashed, double

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


# Plot histogram
canvas = TCanvas("canvas","canvas",10,20,700,800);
canvas.Divide(1,2)
canvas.cd(1)
gStyle.SetOptStat(10) # Only number of entries
h1.Draw("e")

# Define a constant function

f1 = TF1("f1","[0]")

#Choose Minuit2 and setup TFitResultPtr pointers for each fit

#ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2"); 
result_neyman = TFitResultPtr()
result_pearson = TFitResultPtr()
result_likelihood = TFitResultPtr()

# Fit histo with f1 using Neyman, Pearson and Likelihood
# (1) Neyman chi2 (default fit, using observed errors)
f1.SetLineColor(kRed);
result_neyman = h1.Fit(f1,"S"); 
#(2) P: Use Pearson chi2 (using expected errors instead of observed errors)
f1.SetLineColor(kBlue);
result_pearson = h1.Fit(f1,"S P+ ");
#(3) L: Use Loglikelihood method (default is chisquare method)
f1.SetLineColor(kGreen+1);
result_likelihood = h1.Fit(f1,"S L+");

# Summarize the three results in a TGraph

g1=TGraphErrors(1)
g2=TGraphErrors(1)
g3=TGraphErrors(1)

g1.SetMarkerColor(kRed)
g1.SetLineColor(kRed)
g1.SetPoint(0,1,result_neyman.Value(0))
g1.SetMarkerStyle(20)
g1.SetPointError(0,0,result_neyman.Error(0))

g2.SetMarkerStyle(20)
g2.SetMarkerColor(kBlue)
g2.SetLineColor(kBlue)
g2.SetPoint(0,2,result_pearson.Value(0))
g2.SetPointError(0,0,result_pearson.Error(0))

g3.SetMarkerStyle(20)
g3.SetMarkerColor(kGreen+1)
g3.SetLineColor(kGreen+1)
g3.SetPoint(0,3,result_likelihood.Value(0))
g3.SetPointError(0,0,result_likelihood.Error(0))

canvas.cd(2);
# double trueValue = double(N)/h1->GetNbinsX();  // true result of fit
trueValue = double(N)/h1.GetNbinsX();  # true result of fit
gPad.DrawFrame(0.5,trueValue*0.85,3.5,trueValue*1.1,"Summary of fit results;;Fit Bias");
g1.Draw("EP0"); #https://root.cern.ch/root/htmldoc/guides/users-guide/Histograms.html
g2.Draw("EP0");
g3.Draw("EP0");

line = TLine(0.5,trueValue,3.5,trueValue);
line.SetLineStyle(kDashed);
line.SetLineColor(1);
line.Draw();
gPad.Update();

print("Neyman     fit bias = "+str(result_neyman.Value(0)-trueValue))
print("Pearson    fit bias = "+str(result_pearson.Value(0)-trueValue))
print("Likelihood fit bias = "+str(result_likelihood.Value(0)-trueValue))