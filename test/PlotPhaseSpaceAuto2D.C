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

void PlotPhaseSpaceAuto2D( const TString &sim, Int_t time, Int_t index = 0,TString dcom = "p1:x1", const TString &options="") {
  
#ifdef __CINT__  
  gSystem->Load("libplasma.so");
#endif

  PlasmaGlob::Initialize();

  TString opt = options;
 
  // Palettes!
  gROOT->Macro("PlasmaPalettes.C");
  
  if(opt.Contains("gridx")) {
    gStyle->SetPadGridX(1);
  } else if(opt.Contains("gridy")) {
    gStyle->SetPadGridY(1);
  }
  
  gStyle->SetTitleAlign(22);
  gStyle->SetPadRightMargin(0.17);   // Margin for palettes in 2D histos

  gStyle->SetLabelFont(42,"xyz");
  gStyle->SetTitleFont(42,"xyz");

  // Load PData
  PData *pData = PData::Get(sim.Data());
  pData->LoadFileNames(time);
  if(!pData->IsInit()) return;

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
  
  //opt += "comovcenter";

  // Centering time and z position:
  Double_t shiftz = pData->Shift(opt);
  TString sshiftz = Form("(x1-%f)",shiftz);
  
  //cout << Form("string = %s",sshiftz.Data()) << endl; 
  TString sx1 = "x1";
  dcom.ReplaceAll(sx1,sshiftz);
  //cout << Form("%s replaced by %s result = %s",sx1.Data(),sshiftz.Data(),dcom.Data()) << endl;
  // ----

  // Weighting by the macroparticle charge:
  char cutString[128];
  sprintf(cutString,"TMath::Abs(q)"); 
  
  TCut Cut = cutString;
  

  // Bining, intervals, labels, etc.
  Int_t xNbin = 200;
  Int_t yNbin = 200;
  
  // Pointer to data TTree
  TTree *tree = pData->GetTreeRaw(pData->GetRawFileName(index)->c_str(),opt);
  
  // Get phasespace histos
  TH2F *hPha2D = NULL;
  char hName[24]; 
  sprintf(hName,"hPha2D");
  char dCommand[128];
  sprintf(dCommand,"%s>>%s(%i,%i)",dcom.Data(),hName,xNbin,yNbin);

  tree->Draw(dCommand,Cut,"goff");
  hPha2D = (TH2F*) gROOT->FindObject(hName);    
  hPha2D->GetXaxis()->CenterTitle();
  hPha2D->GetYaxis()->CenterTitle();
  hPha2D->GetZaxis()->CenterTitle();

  // hPha2D->GetXaxis()->SetTitle("#zeta [c/#omega_{p}]");
  // hPha2D->GetYaxis()->SetTitle("y [c/#omega_{p}]");
  // hPha2D->GetZaxis()->SetTitle("dn/d#zetady");
  
  Float_t xMean = hPha2D->GetMean(1);
  Float_t xRms  = hPha2D->GetRMS(1);
  Float_t yMean = hPha2D->GetMean(2);
  Float_t yRms  = hPha2D->GetRMS(2);

  cout << Form(" xMean = %f  xRms = %f   yMean = %f  yRms = %f",xMean,xRms,yMean,yRms) << endl;

  // Plotting
  // -----------------------------------------------
    
  // Canvas setup
  TCanvas *C;
  if(opt.Contains("hres") && !opt.Contains("pdf")) // high resolution for plain graphics output.
    C = new TCanvas("C","Phasespace",1200,750);
  else
    C = new TCanvas("C","Phasespace",960,600);

    
  // Text objects
  TPaveText *textTime = new TPaveText(0.55,0.85,0.82,0.9,"NDC");
  PlasmaGlob::SetPaveTextStyle(textTime,32); 
  char ctext[128];
  if(opt.Contains("units") && pData->GetPlasmaDensity()) 
    sprintf(ctext,"z = %5.1f mm", Time * skindepth / PUnits::mm);
  else
    sprintf(ctext,"t = %5.1f #omega_{p}^{-1}",Time);
  textTime->AddText(ctext);
 

  // Actual Plotting!
  // ------------------------------------------------------------

  // Output file
  TString fOutName = Form("./%s/Plots/PhaseSpaceAuto2D/PhaseSpaceAuto2D",sim.Data());
  fOutName += Form("-%s_%i",sim.Data(),time);

  C->cd();
 
  gPad->SetFrameLineWidth(2);  

  if(opt.Contains("logz")) {
    gPad->SetLogz(1);
  } else {
    gPad->SetLogz(0);
  }

  // Set palette:
  PPalette * pPalette = (PPalette*) gROOT->FindObject("electron");
  pPalette->cd();

  hPha2D->Draw("colz");
  
  textTime->Draw();
 
  gPad->RedrawAxis(); 

  C->cd();

  // Print to a file
  PlasmaGlob::imgconv(C,fOutName,opt);
  // ---------------------------------------------------------

}
