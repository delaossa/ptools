#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>

#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TCut.h>
#include <TChain.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TPaletteAxis.h>
#include <TExec.h>

#include "PData.hh"
#include "PlasmaGlob.hh"
#include "PPalette.hh"

using namespace std;

void PlotBunch( const TString &sim, Int_t time, Int_t index = 0, const TString &options="") {
  
#ifdef __CINT__  
  gSystem->Load("libplasma.so");
#endif

  PlasmaGlob::Initialize();

  TString opt = options;
 
  // Palettes!
  gROOT->Macro("PlasmaPalettes.C");

  if(opt.Contains("grid")) {
    gStyle->SetPadGridX(1);
    gStyle->SetPadGridY(1);
  }
  gStyle->SetLabelFont(42,"xyz");
  gStyle->SetTextFont(62);
 

  // Load PData
  PData *pData = PData::Get(sim.Data());
  pData->LoadFileNames(time);
  if(!pData->IsInit()) return;


  Int_t Nspecies = pData->NSpecies();
  if(index>Nspecies-1) {
    return;
  }
  if(!pData->GetRawFileName(index)) {
    return;    
  }


  Bool_t CYL = kFALSE;
  if(sim.Contains("cyl")) CYL = kTRUE; 
    
  Bool_t ThreeD = kFALSE;
  if(sim.Contains("3D")) ThreeD = kTRUE; 

  // Some plasma constants
  Double_t n0 = pData->GetPlasmaDensity();
  Double_t kp = pData->GetPlasmaK();
  Double_t skindepth = 1.;
  if(kp!=0.0) skindepth = 1/kp;
  Double_t E0 = pData->GetPlasmaE0();

  // Time in OU
  Float_t Time = pData->GetRealTime();
  // z start of the plasma in normalized units.
  Float_t zStartPlasma = pData->GetPlasmaStart()*kp;
  // z start of the beam in normalized units.
  Float_t zStartBeam = pData->GetBeamStart()*kp;
  
  opt += "comovcenter";

  // Centering time and z position:
  Double_t shiftz = pData->Shift(opt);
  TString sshiftz = Form("(x1-%f)",shiftz);

  if(opt.Contains("center")) {
    Time -= zStartPlasma;
    if(opt.Contains("comov"))      // Centers on the head of the beam.
      Time += zStartBeam;
  } 

  // Bining, intervals, labels, etc.
  Int_t xNbin = 200;
  Int_t pNbin = 200;

  // Spatial resolution
  Float_t dx1 = pData->GetDX(0);
  Float_t dx2 = pData->GetDX(1);
  Float_t dx3 = pData->GetDX(2);
  
  // Spatial coordinates intervals:
  Float_t x1Min = -7.8;
  Float_t x1Max = -7.0;
  Float_t x2Min = -0.5;
  Float_t x2Max =  0.5;
  Float_t x3Min = -0.5;
  Float_t x3Max =  0.5;

  // Momentum coordinates intervals:
  Float_t p1Min =  4050.01;
  Float_t p1Max =  6199.99;
  Float_t p2Min = -15.0;
  Float_t p2Max =  15.0;
  Float_t p3Min = -15.0;
  Float_t p3Max =  15.0;

  char hName[8];
  sprintf(hName,"hP1X1");
  TH2F *hP1X1 = (TH2F*) gROOT->FindObject(hName);
  if(hP1X1) delete hP1X1;
  hP1X1 = new TH2F(hName,"",xNbin,x1Min,x1Max,pNbin,p1Min,p1Max);

  sprintf(hName,"hP2X2");
  TH2F *hP2X2 =  (TH2F*) gROOT->FindObject(hName);
  if(hP2X2) delete hP2X2;
  hP2X2 = new TH2F(hName,"",xNbin,x2Min,x2Max,pNbin,p2Min,p2Max);

  // Sliced quantities:
  // --------------------------------------------------------------------------
 
  UInt_t SNbin = 40;
  Float_t x1BinMin = -7.54;
  Float_t x1BinMax = -7.14;
  if(sim.Contains("DDR.") || sim.Contains(".DR.")) {
    SNbin = 10;
    x1BinMin = -4.7;
    x1BinMax = -3.8;
  } else if(sim.Contains("flash")) {
    SNbin = 20;
    x1BinMin = -4.15;
    x1BinMax = -4.02;
  } else if(sim.Contains("facet_v23kA.G.RI.e")) {
    SNbin = 30;
    x1BinMin = -7.60;
    x1BinMax = -7.35;
  } else if(sim.Contains("facet_v23kA.G.RI")) {
    SNbin = 40;
    x1BinMin = -7.54;
    x1BinMax = -7.14;
  } else if(sim.Contains("facet_v23kA.G")) {
    SNbin = 40;
    x1BinMin = -7.55;
    x1BinMax = -7.15;
  } 
  
  // Set the binning
  Float_t *sBinLim = new Float_t[SNbin+1];
  sBinLim[0] = x1BinMin;
  sBinLim[SNbin] = x1BinMax;
  Float_t slbinSize = (sBinLim[SNbin] - sBinLim[0])/SNbin;
  for(UInt_t i=1;i<SNbin;i++) {
    sBinLim[i] = sBinLim[i-1] + slbinSize;
  }
  
  TH1F **hP1sl = new TH1F*[SNbin];
  TH2F **hP2X2sl = new TH2F*[SNbin];
  for(UInt_t k=0;k<SNbin;k++) {
    sprintf(hName,"hP2X2sl_%2i",k);
    hP2X2sl[k] = (TH2F*) gROOT->FindObject(hName);
    if(hP2X2sl[k]) delete hP2X2sl[k];
    hP2X2sl[k] = new TH2F(hName,"",xNbin,x2Min,x2Max,pNbin,p2Min,p2Max);

    hP2X2sl[k]->GetXaxis()->SetTitle("k_{p} y");
    hP2X2sl[k]->GetYaxis()->SetTitle("p_{y}/mc");


    sprintf(hName,"hP1sl_%2i",k);
    hP1sl[k] = (TH1F*) gROOT->FindObject(hName);
    if(hP1sl[k]) delete hP1sl[k];
    hP1sl[k] = new TH1F(hName,"",pNbin,p1Min,p1Max);

  }


  //Float_t **var = NULL;
  const UInt_t Nvar = 8;
  Float_t *var[Nvar];  
  UInt_t Np = pData->GetRawArray(pData->GetRawFileName(index)->c_str(),var);  
  
  // Filling histos
  cout << endl;
  cout << " 2. Filling histograms from file : " << pData->GetRawFileName(index)->c_str() << endl;
  for(UInt_t i=0;i<Np;i++) {
    var[5][i] = var[5][i] - shiftz;
    
    if(var[5][i]<x1Min || var[5][i]>x1Max ) continue; 
    if(var[6][i]<x2Min || var[6][i]>x2Max ) continue; 
    if(var[7][i]<x3Min || var[7][i]>x3Max ) continue; 
    if(var[1][i]<p1Min || var[1][i]>p1Max ) continue; 
    if(var[2][i]<p2Min || var[2][i]>p2Max ) continue; 
    if(var[3][i]<p3Min || var[3][i]>p3Max ) continue; 
    
    hP1X1->Fill(var[5][i],var[1][i],TMath::Abs(var[4][i]));
    
  }

  
  cout << "\n7. Plotting... " << endl;

  // Canvas setup
  // Create the canvas and the pads before the Frame loop
  // Resolution:
  Int_t sizex = 800;
  Int_t sizey = 600;
  if(opt.Contains("hres")) {
    Int_t sizex = 1600;
    Int_t sizey = 1200;    
  }
  
  TCanvas *C = new TCanvas("C","Bunch properties",sizex,sizey);
  C->cd();

  // Set palette:
  PPalette * pPalette = (PPalette*) gROOT->FindObject("electron");
  pPalette->cd();

  hP1X1->Draw("colz");
  
  // Print to a file
  // Output file
  TString fOutName = Form("./%s/Plots/Bunch/Bunch",sim.Data());
  fOutName += Form("-%s_%i",sim.Data(),time);

  PlasmaGlob::imgconv(C,fOutName,opt);
  // ---------------------------------------------------------
  
  
  
}
