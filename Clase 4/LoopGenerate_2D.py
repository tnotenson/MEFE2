#import libraries
import numpy as np
import matplotlib.pyplot as plt
import os,sys,time
import pickle
import ROOT
from ROOT import TH1D, TRandom3, TCanvas, TF1, gStyle, TFitResultPtr, kRed, kBlue, kGreen, kGray, kCyan, kMagenta, TGraphErrors, gPad, TLine, kDashed, double, gROOT, TMath, TPaveLabel, gMinuit, TTree, TH1F, TLegend, TFile, TArrow, TVirtualFitter, TMatrixDSym, RooEllipse, kWarning, TMarker
from decimal import *
import array

Ntoys = 10

# Datos a Fitear
Hdata  = TH1D("Hdata" ,"Dos fuentes radioactivas",60,0,900);
Hdata.SetTitle("Dos exponenciales;Tiempo (seg);Cuentas / 15 seg");

# Funcion con la que fitear
fTrue = TF1("fTrue","[0]+[1]*exp(-x/[3])+[2]*exp(-x/[4])",0,900);

PrintCanvases = True;   # "false" para no imprimir un archivo pdf

c1 = TCanvas("c1","Dos Fuentes Radioactivas",700,700);
if (PrintCanvases): c1.Print("Out_Generate_2D.pdf["); # open output file

stop_at_plot = True;  # Hacer pausa despues de cada plot

MinuitStatus=0;           # signal of problems
gErrorIgnoreLevel = kWarning; # Dont print "Current canvas added to pdf file"

# Estudio de coberturas: (p3,p4,p34) (CR,LR) (1s,2s,3s)
ntot=0;
p3CR1s=p3LR1s=p3CR2s=p3LR2s=p3CR3s=p3LR3s=0;
p4CR1s=p4LR1s=p4CR2s=p4LR2s=p4CR3s=p4LR3s=0;

InElipse1=InElipse2=InElipse3=0;
InMinos1 =InMinos2 =InMinos3 =0;

###### LOOP SOBRE PSEUDO-EXPERIMENTOS ######

for toy in range(0,Ntoys):

    fTrue.SetParameters(10,1000,100,40,200); 

    r = TRandom3(0); 
    Nevents = r.Poisson(4000);
    for i in range(0,Nevents): Hdata.Fill(fTrue.GetRandom())

    gStyle.SetOptStat(0);
    #Hdata.Draw("E");
    c1.SetLogy();
    c1.SetGrid();

    #Hdata.Fit(fTrue,"N"); 
    Ptr = Hdata.Fit(fTrue,"SN");
    cov = Ptr.GetCovarianceMatrix();

    # Fitted values of lifetimes
    p3_hat=fTrue.GetParameter(3);
    p4_hat=fTrue.GetParameter(4);

    # Donde guarda los errores del fit 
    errorup1        = errorup2        = errorup3        = array.array('d', [0.])
    errordn1        = errordn2        = errordn3        = array.array('d', [0.])
    errparabolic1   = errparabolic2   = errparabolic3   = array.array('d', [0.])
    corr1           = corr2           = corr3           = array.array('d', [0.])
    
    # Ajuste y obtencion de errores para 1,2,3 sigmas

    gMinuit.SetErrorDef(2.296); # 1-sigma     2.296
    gr1b = gMinuit.Contour(80,3,4); 
    gMinuit.mnerrs(4,errorup1,errordn1,errparabolic1,corr1);
    # 
    # 
    # gMinuit.SetErrorDef(6.180);   # 2-sigma    6.180
    # gr2b = gMinuit.Contour(80,3,4); 
    # gMinuit.mnerrs(4,errorup2,errordn2,errparabolic2,corr2);
# 
    # gMinuit.SetErrorDef(11.829);   # 3-sigma    11.829
    # gr3b = gMinuit.Contour(80,3,4); 
    # gMinuit.mnerrs(4,errorup3,errordn3,errparabolic3,corr3);
    

    
    #Codigo de arriba para 2 y 3 sigmas remplazado, pues hay que saltear casos 
    #que el fit no converge con este mensaje de error:
    #    Warning in <TM::Contour>: Cannot find more than 4 points, no TGraph returned
    #    Warning in <TMinuit::Contour>: Returning a TGraph with 4 points only
    

    # Ajuste y obtencion de errores para 2,3 sigmas, salteando si fit no OK

    gMinuit.SetErrorDef(6.180); # 2-sigma  6.180
    gr2b = gMinuit.Contour(80,3,4); 
    MinuitStatus = gMinuit.GetStatus();
    if (MinuitStatus!=0): print("MinuitStatus"+str(MinuitStatus)+"\n")
    if (MinuitStatus!=0): continue;
    if (gr2b.GetN() < 40 ): continue;
    gMinuit.mnerrs(4,errorup2,errordn2,errparabolic2,corr2);
    
    print("Tau4 95.4 CL: {0:.3}+/-{1:.3} [{2:.3},{3:.3}] CI: ({4:.3}, {5:.3})\n".format(Decimal(p4_hat),Decimal(errparabolic2[0]),Decimal(errordn2[0]),Decimal(errorup2[0]),Decimal(p4_hat+errordn2[0]),Decimal(p4_hat+errorup2[0])))

    gMinuit.SetErrorDef(11.829); # 3-sigma 11.829
    gr3b = gMinuit.Contour(80,3,4); 
    MinuitStatus = gMinuit.GetStatus();
    if (MinuitStatus!=0): print("MinuitStatus"+str(MinuitStatus)+"\n")
    if (MinuitStatus!=0): continue;
    if (gr3b.GetN() < 40 ): continue;
    gMinuit.mnerrs(4,errorup3,errordn3,errparabolic3,corr3);
    print("Tau4 95.4 CL: {0:.3}+/-{1:.3} [{2:.3},{3:.3}] CI: ({4:.3}, {5:.3})\n\n".format(Decimal(p4_hat),Decimal(errparabolic3[0]),Decimal(errordn3[0]),Decimal(errorup3[0]),Decimal(p4_hat+errordn3[0]),Decimal(p4_hat+errorup3[0])))

    ### DIBUJA CONTORNOS DE MINOS ###
    
    c1.SetLogy(0);
    c1.DrawFrame(30,0,50,800,"68.3% 95.4% 99.7% CL regions for tau3 & tau4;tau3;tau4");

    gr3b.SetFillColorAlpha(590, 0.5);
    #gr3b->SetFillColorAlpha(45, 0.5);
    gr3b.Draw("lf"); # Draw "contour" for 3-sigma
    gr2b.SetFillColorAlpha(406, 0.5);
    #gr2b->SetFillColorAlpha(42, 0.5);
    gr2b.Draw("Lf");  # Draw "contour" for 2-sigma
    gr1b.SetFillColorAlpha(390, 0.5);
    gr1b.Draw("Lf");  # Draw "contour" for 1-sigma

    ### DIBUJA ELIPSES ###

    cov11=np.sqrt(cov(3,3));
    cov22=np.sqrt(cov(4,4));
    cor12=cov(3,4)/np.sqrt(cov(3,3)*cov(4,4));
    
    ell1 = RooEllipse("Parab CI", p3_hat, p4_hat, 1*cov11, 1*cov22, cor12, 100);
    ell2 = RooEllipse("Parab CI", p3_hat, p4_hat, 2*cov11, 2*cov22, cor12, 100);
    ell3 = RooEllipse("Parab CI", p3_hat, p4_hat, 3*cov11, 3*cov22, cor12, 100);
    
    ell1.SetLineColor(2);
    ell2.SetLineColor(2);
    ell3.SetLineColor(2);
    ell1.SetLineWidth(2);
    ell2.SetLineWidth(2);
    ell3.SetLineWidth(2);
    ell1.SetLineStyle(1);
    ell2.SetLineStyle(1);
    ell3.SetLineStyle(1);
    ell1.SetFillColor(6);
    ell1.Draw("lsame");
    ell2.Draw("lsame");
    ell3.Draw("lsame");

    ### DIBUJA EL PUNTO FITEADO Y EL VERDADRO ###
    
    FittedResult = TMarker(p3_hat,p4_hat,20);
    FittedResult.SetMarkerColor(kRed);
    FittedResult.SetMarkerSize(1.);
    FittedResult.Draw();

    TrueValue = TMarker(40,200,20);
    TrueValue.SetMarkerColor(kBlue);
    TrueValue.SetMarkerStyle(71);
    TrueValue.SetMarkerSize(2);
    TrueValue.Draw();

    if (PrintCanvases): c1.Print("Out_Generate_2D.pdf"); 
    c1.Update();

    ntot+=1;

    # Check if inside Minos Contours and inside Ellipses
    # (uncomment next 6 lines to print which cases are inside)

    # if (TMath.IsInside(40., 200.,  80, gr1b.GetX(), gr1b.GetY())) print("Minos1")
    # if (TMath.IsInside(40., 200.,  80, gr2b.GetX(), gr2b.GetY())) print("Minos2")
    # if (TMath.IsInside(40., 200.,  80, gr3b.GetX(), gr3b.GetY())) print("Minos3")
    # if (TMath.IsInside(40., 200., 101, ell1.GetX(), ell1.GetY())) print("Elipse1")
    # if (TMath.IsInside(40., 200., 101, ell2.GetX(), ell2.GetY())) print("Elipse2")
    # if (TMath.IsInside(40., 200., 101, ell3.GetX(), ell3.GetY())) print("Elipse3")

    if (TMath.IsInside(40., 200.,  80, gr1b.GetX(), gr1b.GetY())): InMinos1+=1;
    if (TMath.IsInside(40., 200.,  80, gr2b.GetX(), gr2b.GetY())): InMinos2+=1;
    if (TMath.IsInside(40., 200.,  80, gr3b.GetX(), gr3b.GetY())): InMinos3+=1;
    if (TMath.IsInside(40., 200., 101, ell1.GetX(), ell1.GetY())): InElipse1+=1;
    if (TMath.IsInside(40., 200., 101, ell2.GetX(), ell2.GetY())): InElipse2+=1;
    if (TMath.IsInside(40., 200., 101, ell3.GetX(), ell3.GetY())): InElipse3+=1;

    if( p4_hat-errparabolic1[0] < 200 and p4_hat+errparabolic1[0] > 200): p4CR1s+=1;
    if( p4_hat+errordn1[0]      < 200 and p4_hat+errorup1[0]      > 200): p4LR1s+=1;
    if( p4_hat-errparabolic2[0] < 200 and p4_hat+errparabolic2[0] > 200): p4CR2s+=1;
    if( p4_hat+errordn2[0]      < 200 and p4_hat+errorup2[0]      > 200): p4LR2s+=1;
    if( p4_hat-errparabolic3[0] < 200 and p4_hat+errparabolic3[0] > 200): p4CR3s+=1;
    if( p4_hat+errordn3[0]      < 200 and p4_hat+errorup3[0]      > 200): p4LR3s+=1;

    if (sys.stdin.read(1)=='.'): stop_at_plot = False # Hit "." to stop drawing

    Hdata.Reset();
    del gr1b;
    del gr2b;
    del gr3b;
    del FittedResult;
    del TrueValue;
# Fin de los Toys

# Cobertura p4 CR (Crame-Rao) y LR (Likelihood-Ratio) error 1-2-3 sigmas

CLp4CR1s = 100.*p4CR1s/ntot;
CLp4LR1s = 100.*p4LR1s/ntot;
CLp4CR2s = 100.*p4CR2s/ntot;
CLp4LR2s = 100.*p4LR2s/ntot;
CLp4CR3s = 100.*p4CR3s/ntot;
CLp4LR3s = 100.*p4LR3s/ntot;

# Cobertura Minos y Elipse,  1-2-3 sigmas
CL_InMinos1  = 100.*InMinos1 /ntot;
CL_InMinos2  = 100.*InMinos2 /ntot;
CL_InMinos3  = 100.*InMinos3 /ntot;
CL_InElipse1 = 100.*InElipse1/ntot;
CL_InElipse2 = 100.*InElipse2/ntot;
CL_InElipse3 = 100.*InElipse3/ntot;

### PRINT SUMMARY STATISTICS ###

print("Ntot: "+str(ntot)+"\n")

print("Parametro 4: Contador CR 1s,2s,3s:  ")
print(str(p4CR1s)+"  "+str(p4CR2s)+"  "+str(p4CR3s))
print("Parametro 4: Contador LR 1s,2s,3s:  ")
print(str(p4LR1s)+"  "+str(p4LR2s)+"  "+str(p4LR3s)+"\n")
print("Parametro 4: Coverage CR 1s,2s,3s:  ")
print(str(CLp4CR1s)+"  "+str(CLp4CR2s)+"  "+str(CLp4CR3s))
print("Parametro 4: Coverage LR 1s,2s,3s:  ")
print(str(CLp4LR1s)+"  "+str(CLp4LR2s)+"  "+str(CLp4LR3s)+"\n")
print("Contador MINOS 2D para 1s,2s,3s:   ")
print(str(InMinos1)+" "+str(InMinos2)+" "+str(InMinos3))
print("Coverage MINOS 2D para 1s,2s,3s:   ")
print(str(CL_InMinos1)+"  "+str(CL_InMinos2)+"  "+str(CL_InMinos3)+"\n")
print("Contador Elipse 2D para 1s,2s,3s:  ")
print(str(InElipse1)+" "+str(InElipse2)+" "+str(InElipse3))
print("Coverage Elipse 2D para 1s,2s,3s:  ")
print(str(CL_InElipse1)+"  "+str(CL_InElipse2)+"  "+str(CL_InElipse3)+"\n")

if (PrintCanvases): c1.Print("Out_Generate_2D.pdf]");  # close file