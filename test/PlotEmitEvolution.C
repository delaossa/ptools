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

void PlotEmitEvolution(const TString &sim, const TString &phaname = "p1x1", const TString &options="png") { 
  
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

  // Resolution:
  Int_t sizex = 800;
  Int_t sizey = 600;
  if(opt.Contains("hres")) {
    Int_t sizex = 1024;
    Int_t sizey = 768;    
  }
  
  TString filename;
  filename = Form("./%s/Plots/EmittanceEvolution/Evolutions-%s-%s.root",sim.Data(),sim.Data(),phaname.Data());
  
  TFile  *ifile = (TFile*) gROOT->GetListOfFiles()->FindObject(filename.Data());
  if (!ifile) ifile = new TFile(filename,"READ");


  Float_t maxEmit = -999.;
  Float_t minEmit = 999.;
  Float_t maxTime = -999.;
  Float_t minTime = 999.;
  TGraph *gEmit = NULL;
  gEmit = (TGraph*) ifile->Get("gEmitvsTime");
       
  if( !gEmit ) return;
      
  // Calculate the max and min of every set of graphs:
  
  Int_t Npoints = gEmit->GetN();
  Double_t *yEmit = gEmit->GetY();
  Double_t *xTime = gEmit->GetX();
  
  for(Int_t j=0;j<Npoints;j++) {
    if(yEmit[j]>maxEmit)
      maxEmit = yEmit[j];
    if(yEmit[j]<minEmit)
      maxEmit = yEmit[j];
    if(xTime[j]>maxTime)
      maxTime = xTime[j];
    if(xTime[j]<minTime)
      maxTime = xTime[j];
  }

  // --------------------

  // Canvas setup
  TCanvas *C1 = new TCanvas("C1","Evolution of Emittance",sizex,sizey);
  
  C1->cd(0);
  gPad->SetFrameLineWidth(2);  
  gPad->SetRightMargin(0.10);
  
  Float_t marginY = (maxEmit - minEmit)/10;
  Float_t E0   = minEmit-marginY;
  Float_t E1   = maxEmit+marginY;
  Float_t marginX = (maxTime - minTime)/10;
  Float_t T0   = minTime-marginX;
  Float_t T1   = maxTime+marginX;
  TH1F *hFrame = new TH1F("hFrame","",100,0,30);
  hFrame->GetYaxis()->SetRangeUser(1.5,2.5);
  hFrame->GetXaxis()->SetTitle("propagation length [mm]");
  hFrame->GetYaxis()->SetTitle("Emittance [MeV/c #mum]");
  PlasmaGlob::SetH1LabelSize(hFrame);
 
  hFrame->Draw("axis");

  gPad->Update();
  gEmit->Draw("C");
   
  C1->cd(0);
  
  // Print to a file
  // Output file
  TString fOutName = Form("./%s/Plots/EmittanceEvolution/Emittance-%s-%s",sim.Data(),sim.Data(),phaname.Data());
  PlasmaGlob::imgconv(C1,fOutName,opt);
  
}
