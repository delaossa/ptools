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
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TExec.h>

#include "PData.hh"
#include "PlasmaGlob.hh"
#include "PPalette.hh"

void PlotPlasmaDensity( const TString &sim, UInt_t time, Int_t zoom=2, const TString &opt="") {
  
  PlasmaGlob::Initialize();

  // Palettes!
  gROOT->Macro("PlasmaPalettes.C");

  // More makeup
  gStyle->SetPadRightMargin(0.12); // Margin for palettes in 2D histos
  
  // Load PData
  PData *pData = PData::Get(sim.Data());
  pData->LoadFileNames(time);
  if(!pData->IsInit()) return;

  Bool_t ThreeD = kFALSE;
  if(sim.Contains("3D")) ThreeD = kTRUE; 

  // Get charge density histos
  Int_t index = 0;
  if(opt.Contains("beam")) index = 1;

  TH2F *hDen2D = NULL;
  if(!ThreeD)
    hDen2D = pData->GetCharge(index);
  else {
    if(opt.Contains("sum"))
      hDen2D = pData->GetChargeR(index,"sum");
    else
      hDen2D = pData->GetChargeR(index);
  }
  char hName[24];
  sprintf(hName,"hDen_%i",0);
  hDen2D->SetName(hName);
  hDen2D->GetXaxis()->CenterTitle();
  hDen2D->GetYaxis()->CenterTitle();
  hDen2D->GetXaxis()->SetTitle("z [c/#omega_{p}]");
  hDen2D->GetYaxis()->SetTitle("y [c/#omega_{p}]");

  // Tunning the Histograms
  // ---------------------
  
  Float_t Time = hDen2D->GetXaxis()->GetXmin();
  
  // Set the range of the histogram for maximum constrast
  Float_t density = 1;
  if(opt.Contains("units") && pData->GetPlasmaDensity())
    density = 1e-15 * 1e-6 * pData->GetPlasmaDensity();
  
  if(!opt.Contains("beam")) {
    Float_t Max  = 1.1 * hDen2D->GetMaximum();
    Float_t Base = density;
    Float_t Min  = 2.* Base - Max;
    if(Max >= 2. * Base) {
      Min = 0;
    } else if(Max<1.0 * Base) {
      Max = 1.1 * Base;
      Min = 0.;
    }
    
    hDen2D->GetZaxis()->SetRangeUser(Min,Max);  
  }

  // Zoom
  Float_t range    = (hDen2D->GetYaxis()->GetXmax() - hDen2D->GetYaxis()->GetXmin())/zoom;
  Float_t midPoint = (hDen2D->GetYaxis()->GetXmax() + hDen2D->GetYaxis()->GetXmin())/2.;
  
  hDen2D->GetYaxis()->SetRangeUser(midPoint-range/2,midPoint+range/2);  

  // Plotting
  // -----------------------------------------------

  // Canvas setup
  TCanvas *C = new TCanvas("C","Charge density and Electric field",850,500);

  // Text objects
  TPaveText *text1 = new TPaveText(0.6,0.94,0.88,1.0,"NDC");
  PlasmaGlob::SetPaveStyle(text1);
  text1->AddText("Charge density [a.u.]");
  
  TPaveText *textTime = new TPaveText(0.7,0.82,0.85,0.88,"NDC");
  PlasmaGlob::SetPaveStyle(textTime); 
  char ctext[128];
  sprintf(ctext,"T = %5.1f [1/#omega_{p}]",Time);
  textTime->AddText(ctext);
  
  
  // Actual Plotting!
  // ------------------------------------------------------------
  
  // Output file
  TString fOutName = Form("./%s/Plots/PlasmaDensity/PlasmaDensity",sim.Data());
  fOutName += Form("-%s_%i",sim.Data(),time);
  
  C->cd(0);
  
  C->cd(1); // <--- Top Plot
  gPad->SetFrameLineWidth(3);  
  
  PPalette * electronPalette = (PPalette*) gROOT->FindObject("electron");
  PPalette * plasmaPalette = (PPalette*) gROOT->FindObject("plasma");
  if(index == 0)
    plasmaPalette->cd();
  else
    electronPalette->cd();
  
  hDen2D->Draw("colz");
  
  text1->Draw();
  textTime->Draw();

  C->cd();

  // Print to a file
  PlasmaGlob::imgconv(C,fOutName,opt);
  // ---------------------------------------------------------
}
