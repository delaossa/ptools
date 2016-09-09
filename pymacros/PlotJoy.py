#!/usr/bin/env python

from ROOT import PData, PGlobals, gStyle, TObject, TFile, TCanvas, TH1F

pData = PData('rake-v10kA.G.SR2.RI.3D')

pData.LoadFileNames(20)

PGlobals.Initialize()

gStyle.SetCanvasColor(1)
gStyle.SetFrameFillColor(1)
gStyle.SetFrameFillColor(1)
gStyle.SetHistLineColor(0)
gStyle.SetHistLineWidth(1)
gStyle.SetPadLeftMargin(0.20)
  
ofile = TFile("joy.root","RECREATE")
NbinsX1 = pData.GetX1N()

hDen1D = []
#TH1F **hDen1D = new TH1F*[NbinsX1]

TOP = 100
shift = 0.2
MAX = -999
MAX2 = 2*(NbinsX1 * shift) 
  
C = TCanvas("C","",800,1280)
  
j =  (pData.GetX3iMin() + pData.GetX3iMax()) / 2

for i in range(pData.GetX1iMin(),pData.GetX1iMax()) :
    
    hName = "hDen1D_%i" % i
    if pData.Is3D() :
        hDen1D.append(pData.GetH1SliceX3D(pData.GetChargeFileName(0).c_str(),"charge",i,i,j,j,"avg"))

    hDen1D[i].SetName(hName)

    for k in range (1,hDen1D[i].GetNbinsX()+1) :

        if(hDen1D[i].GetBinContent(k)<TOP) :
	    hDen1D[i].SetBinContent(k,hDen1D[i].GetBinContent(k) + i*shift)
        else :
	    hDen1D[i].SetBinContent(k,TOP + i*shift)

        if(MAX<hDen1D[i].GetMaximum()) :
            MAX = hDen1D[i].GetMaximum()
    
#        hDen1D[i].Write(hName,TObject.kOverwrite)
    
  
hFrame = hDen1D[0].Clone("hFrame")
hFrame.Reset()
hFrame.GetYaxis().SetRangeUser(0.0,MAX2)

hFrame.Draw("AC")
for i in range(pData.GetX1iMax()-1,pData.GetX1iMin(), -NbinsX1/50):
    hDen1D[i].SetFillColor(1)
    hDen1D[i].Draw("FC same")
  
PGlobals.imgconv(C,"joy","pdf")
