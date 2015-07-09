#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>

#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraphErrors.h>
#include <TMarker.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TMath.h>
#include <TPaveText.h>
#include <TExec.h>
#include <TGaxis.h>
#include <TPaletteAxis.h>

#include "PData.hh"
#include "PlasmaGlob.hh"
#include "PPalette.hh"

void PlotRakeEvolutions(const TString &sim, const TString &options="png") { 
  
#ifdef __CINT__  
  gSystem->Load("libplasma.so");
#endif

  PlasmaGlob::Initialize();
  
  // Palettes!
  gROOT->Macro("PlasmaPalettes.C");

  TString opt = options;
  
  // More makeup            
  Float_t margins[4] = {0.15,0.15,0.20,0.10};
  gStyle->SetPadLeftMargin(margins[0]);  // Margin left axis  
  gStyle->SetPadRightMargin(margins[2]);
  gStyle->SetPadTopMargin(margins[3]);  // Margin left axis  
  gStyle->SetPadBottomMargin(margins[1]);

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);


  if(opt.Contains("grid")) {
    gStyle->SetPadGridX(1);
    gStyle->SetPadGridY(1);
  }

  const Int_t Nspaces = 4;
  TString phaname[Nspaces] = {"p1x1","p2x2","p3x3","x2x1"};
  TGraph *gXmean[Nspaces]; 
  TGraph *gYmean[Nspaces]; 
  TGraph *gXrms[Nspaces]; 
  TGraph *gYrms[Nspaces]; 
  TGraph *gEmit[Nspaces];
  TGraph *gCharge = NULL;
  // Special graph with an Energy spread band:
  TGraphErrors *gEneRms = NULL;
 
  Float_t maxEmit[Nspaces] = { -999., -999., -999., -999.};
  Float_t minEmit[Nspaces] = { 999., 999.,999., 999.};
  Float_t maxXmean[Nspaces] = { -999., -999., -999., -999.};
  Float_t minXmean[Nspaces] = { 999., 999., 999., 999.};
  Float_t maxXrms[Nspaces] = { -999., -999., -999., -999.};
  Float_t minXrms[Nspaces] = { 999., 999., 999., 999.};
  Float_t maxYmean[Nspaces] = { -999., -999., -999., -999.};
  Float_t minYmean[Nspaces] = { 999., 999., 999.};
  Float_t maxYrms[Nspaces] = { -999., -999., -999., -999.};
  Float_t minYrms[Nspaces] = { 999., 999., 999., 999.};
  Float_t maxCharge =-999.;
  Float_t minCharge = 999.;
  
  // Resolution:
  Int_t sizex = 600;
  Int_t sizey = 800;
  if(opt.Contains("hres")) {
    Int_t sizex = 1024;
    Int_t sizey = 768;    
  }
  

  for(Int_t i=0;i<Nspaces;i++) {
    TString filename;
    filename = Form("./%s/Plots/EmittanceEvolution/Evolutions-%s-%s.root",sim.Data(),sim.Data(),phaname[i].Data());
    
    TFile  *ifile = (TFile*) gROOT->GetListOfFiles()->FindObject(filename.Data());
    if (!ifile) ifile = new TFile(filename,"READ");
    
    
    gEmit[i] = (TGraph*) ifile->Get("gEmitvsTime");  
    gXmean[i] = (TGraph*) ifile->Get("gXmeanvsTime");
    gXrms[i] = (TGraph*) ifile->Get("gXrmsvsTime");
    gYmean[i] = (TGraph*) ifile->Get("gYmeanvsTime");
    gYrms[i] = (TGraph*) ifile->Get("gYrmsvsTime");

    // Energy spread
    if(i==0) {
      gCharge = (TGraph*) ifile->Get("gChargevsTime");  
      Int_t Npoints = gCharge->GetN();
      Double_t *yCharge = gCharge->GetY();
      for(Int_t j=0;j<Npoints;j++) {
	if(yCharge[j]>maxCharge)
	  maxCharge = yCharge[j];
	if(yCharge[j]<minCharge)
	  minCharge = yCharge[j];
      }
      Double_t *xValues = gYmean[i]->GetX();
      Double_t *yMean = gYmean[i]->GetY();
      Double_t *yRms = gYrms[i]->GetY();
      Npoints =  gYmean[i]->GetN();

      gEneRms = new TGraphErrors(Npoints,xValues,yMean,0,yRms);
      // cout << "NPOints = " << gEneRms->GetN() << endl;
      // for(Int_t j=0;j<gEneRms->GetN();j++) {
      // 	gEneRms->SetPointError(j,0.0,yRms[j]);
      // 	cout << "eoooo" << endl;
      // }
    }
    // Calculate the max and min of every set of graphs:
  
    Int_t Npoints = gEmit[i]->GetN();
    Double_t *yEmit = gEmit[i]->GetY();
    for(Int_t j=0;j<Npoints;j++) {
      if(yEmit[j]>maxEmit[i])
	maxEmit[i] = yEmit[j];
      if(yEmit[j]<minEmit[i])
	minEmit[i] = yEmit[j];
    }

    Npoints = gXmean[i]->GetN();
    Double_t *yXmean = gXmean[i]->GetY();
    for(Int_t j=0;j<Npoints;j++) {
      if(yXmean[j]>maxXmean[i])
	maxXmean[i] = yXmean[j];
      if(yXmean[j]<minXmean[i])
	minXmean[i] = yXmean[j];
    }

    Npoints = gXrms[i]->GetN();
    Double_t *yXrms = gXrms[i]->GetY();
    for(Int_t j=0;j<Npoints;j++) {
      if(yXrms[j]>maxXrms[i])
	maxXrms[i] = yXrms[j];
      if(yXrms[j]<minXrms[i])
	minXrms[i] = yXrms[j];
    }

    Npoints = gYmean[i]->GetN();
    Double_t *yYmean = gYmean[i]->GetY();
    for(Int_t j=0;j<Npoints;j++) {
      if(yYmean[j]>maxYmean[i])
	maxYmean[i] = yYmean[j];
      if(yYmean[j]<minYmean[i])
	minYmean[i] = yYmean[j];
    }

    Npoints = gYrms[i]->GetN();
    Double_t *yYrms = gYrms[i]->GetY();
    for(Int_t j=0;j<Npoints;j++) {
      if(yYrms[j]>maxYrms[i])
	maxYrms[i] = yYrms[j];
      if(yYrms[j]<minYrms[i])
	minYrms[i] = yYrms[j];
    }
    
  }

  // --------------------

  // Canvas setup
  TCanvas *C1 = new TCanvas("C1","Evolution of Emittance",sizex,sizey);
  C1->cd();
  
  // C1->Divide(1,Nspaces);
  // PlasmaGlob::CanvasPartition(C1,pad,Nspaces);
  
  const Int_t Npads = 3;

  TH1F *hFrame[Npads];
  TPad **pad = new TPad*[Npads];
  
  // Setup Pad layout:
  Double_t lMargin = 0.15;
  Double_t rMargin = 0.15;
  Double_t bMargin = 0.08;
  Double_t tMargin = 0.08;
  Double_t vSpacing = 0.01; 
  Double_t hStep = (1.-lMargin-rMargin);
  Double_t vStep = (1.- bMargin - tMargin - (Npads-1) * vSpacing) / Npads;
  
  Float_t *vposd = new Float_t[Npads];
  Float_t *vposu = new Float_t[Npads];
  Float_t *vmard = new Float_t[Npads];
  Float_t *vmaru = new Float_t[Npads];
  Float_t *vfactor = new Float_t[Npads];
  Float_t *hposl = new Float_t[Npads];
  Float_t *hposr = new Float_t[Npads];
  Float_t *hmarl = new Float_t[Npads];
  Float_t *hmarr = new Float_t[Npads];
  Float_t *hfactor = new Float_t[Npads];
  
  for(Int_t i=0;i<Npads;i++) {
    
    hposl[i] = 0.0;
    hposr[i] = 1.0;
    hmarl[i] = lMargin;
    hmarr[i] = rMargin;
    
    if(i==0) {
      vposd[i] = 0.0;
      vposu[i] = bMargin + vStep;
      vfactor[i] = vposu[i]-vposd[i];  
      vmard[i] = bMargin / vfactor[i];
      vmaru[i] = 0.0;
    } else if(i == Npads-1) {
      vposd[i] = vposu[i-1] + vSpacing;
      vposu[i] = vposd[i] + vStep + tMargin;
      vfactor[i] = vposu[i]-vposd[i];   
      vmard[i] = 0.0;
      vmaru[i] = tMargin / (vposu[i]-vposd[i]);
    } else {
      vposd[i] = vposu[i-1] + vSpacing;
      vposu[i] = vposd[i] + vStep; 
      vfactor[i] = vposu[i]-vposd[i];
      vmard[i] = 0.0;
      vmaru[i] = 0.0;
    } 
    hfactor[i] = hposl[i]-hposr[i];
    
    C1->cd();
    
    char name[16];
    sprintf(name,"pad_%i",i);
    cout << endl << Form("%s :   %4.2f  %4.2f  %4.2f  %4.2f",name,hposl[i],vposd[i],hposr[i],vposu[i]) << endl;
    pad[i] = new TPad(name,"",hposl[i],vposd[i],hposr[i],vposu[i]);
    
    pad[i]->SetLeftMargin(hmarl[i]);
    pad[i]->SetRightMargin(hmarr[i]);  
    pad[i]->SetBottomMargin(vmard[i]);
    pad[i]->SetTopMargin(vmaru[i]);
    pad[i]->Draw();
    
    pad[i]->cd();
  
    sprintf(name,"hFrame_%i",i);  
    hFrame[i] = new TH1F(name,"",100,0,20);


    hFrame[i]->GetYaxis()->SetLabelSize(0.02/vfactor[i]);
    hFrame[i]->GetYaxis()->SetLabelOffset(0.01/vfactor[i]);

    hFrame[i]->GetYaxis()->SetTitleSize(0.02/vfactor[i]);
    hFrame[i]->GetYaxis()->SetTitleOffset(1.0-vfactor[i]);///vfactor[i]);

    hFrame[i]->GetXaxis()->SetLabelSize(0.07);
    hFrame[i]->GetXaxis()->SetLabelOffset(0.01/vfactor[i]);

    hFrame[i]->GetXaxis()->SetTitleSize(0.07);
    hFrame[i]->GetXaxis()->SetTitleOffset(1.1);///vfactor[i]);
    
    hFrame[i]->Draw("axis");    

    if(i==2) {
      Float_t yMin = 0.0001; // minYmean[0] - (maxYmean[0]-minYmean[0])*0.1;
      Float_t yMax = maxYmean[0] + (maxYmean[0]-minYmean[0])*0.1;
      hFrame[i]->GetYaxis()->SetRangeUser(yMin,yMax);
      hFrame[i]->GetXaxis()->SetTitle("propagation length [mm]");
      hFrame[i]->GetYaxis()->SetTitle("p_{z} [GeV/c]");
      //    PlasmaGlob::SetH1LabelSize(hFrame[i]);
      
      gEneRms->SetFillColor(kGray);
      gEneRms->Draw("3");
      gYmean[0]->Draw("L");  

      // Charge on right axis:
      Float_t yrMin = 0.00001; // minCharge - (maxCharge-minCharge)*0.1;
      Float_t yrMax = maxCharge + (maxCharge-minCharge)*0.1;
      Float_t slope = (yMax-yMin)/(yrMax-yrMin);
      
      Double_t *x = gCharge->GetX();
      Double_t *y = gCharge->GetY();
      for(Int_t j=0;j<gCharge->GetN();j++) {
      	gCharge->SetPoint(j,x[j],slope*(y[j]-yrMin)+yMin);
      }

      //hFrame[i]->GetYaxis()->SetRangeUser(yrMin,yrMax);
      gCharge->SetLineStyle(2);
      gCharge->SetLineWidth(3);
      gCharge->SetLineColor(kAzure-8);      
      gCharge->Draw("L");

      pad[i]->Update();

    } else if(i==1) {
      Float_t yMin,yMax;
      if(minXrms[1]<minXrms[0])
	yMin = minXrms[1];
      else
	yMin = minXrms[0];
      if(maxXrms[1]>maxXrms[0])
	yMax = maxXrms[1];
      else
	yMax = maxXrms[0];

      Float_t yDist = yMax - yMin;
      yMin -= 0.1*yDist; 
      yMax += 0.1*yDist;
      hFrame[i]->GetYaxis()->SetRangeUser(yMin,yMax);
      hFrame[i]->GetXaxis()->SetTitle("propagation length [mm]");
      hFrame[i]->GetYaxis()->SetTitle("rms size [#mum]");
      //    PlasmaGlob::SetH1LabelSize(hFrame[i]);

      gXrms[0]->SetLineStyle(1);
      gXrms[0]->SetLineWidth(3);
      gXrms[0]->SetLineColor(kOrange+10);      
      gXrms[0]->Draw("L");
      gXrms[1]->SetLineStyle(2);
      gXrms[1]->SetLineWidth(3);
      gXrms[1]->SetLineColor(kAzure-8);
      gXrms[1]->Draw("L");
 
      pad[i]->Update();
    } else if(i==0) {
      Float_t yMin,yMax;
      if(minEmit[1]<minEmit[0])
	yMin = minEmit[1];
      else
	yMin = minEmit[0];
      if(maxEmit[1]>maxEmit[0])
	yMax = maxEmit[1];
      else
	yMax = maxEmit[0];

      Float_t yDist = yMax - yMin;
      yMin -= 0.1*yDist; 
      yMax += 0.1*yDist;
      hFrame[i]->GetYaxis()->SetRangeUser(yMin,yMax);
      hFrame[i]->GetXaxis()->SetTitle("propagation length [mm]");
      hFrame[i]->GetYaxis()->SetTitle("trans. emittance [(MeV/c) #mum]");
      // PlasmaGlob::SetH1LabelSize(hFrame[i]);
      
      gEmit[0]->SetLineStyle(1);
      gEmit[0]->SetLineWidth(3);
      gEmit[0]->SetLineColor(kOrange+10);
      gEmit[0]->Draw("L");
      gEmit[1]->SetLineStyle(2);
      gEmit[1]->SetLineWidth(3);
      gEmit[1]->SetLineColor(kAzure-8);
      gEmit[1]->Draw("L");
 
      pad[i]->Update();

    }
  }
    
  C1->cd();
  
  // Print to a file
  // Output file
  TString fOutName = Form("./%s/Plots/RakeEvolution/RakeEvolution-%s",sim.Data(),sim.Data());
  PlasmaGlob::imgconv(C1,fOutName,opt);
  
}
