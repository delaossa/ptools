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
#include <TPaletteAxis.h>
#include <TExec.h>

#include "PData.hh"
#include "PDataHiP.hh"
#include "PGlobals.hh"
#include "PPalette.hh"

void PlotCharge3D( const TString &sim, Int_t time, Float_t zoom=2, const TString &options="") {
    
#ifdef __CINT__  
  gSystem->Load("libptools.so");
#endif
  
  PData *pData = PData::Get(sim.Data());
  if(pData->isHiPACE()) {
    delete pData; pData = NULL;
    pData = PDataHiP::Get(sim.Data());
  }
  
  pData->LoadFileNames(time);
  if(!pData->IsInit()) return;

  if(!pData->Is3D()) {
    cout << " It is not a 3D simulation: exiting..." << endl;
    return;
  }
  
  // Refresh and Style
  PGlobals::Initialize();

  // Coloured palettes
  gROOT->Macro("PPalettes.C");

  TString opt = options;
 
  // Style
  gStyle->SetCanvasPreferGL(1);
  if(opt.Contains("grid")) {
    gStyle->SetPadGridX(1);
    gStyle->SetPadGridY(1);
  }
  gStyle->SetNumberContours(64);

  gStyle->SetTitleFont(42);
  gStyle->SetStatFont(42);
  gStyle->SetTextFont(42);
  gStyle->SetTitleFont(42,"xyz");
  gStyle->SetLabelFont(42,"xyz");

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
  
  if(opt.Contains("center")) {
    Time -= zStartPlasma;
    if(opt.Contains("comov"))      // Centers on the head of the beam.
      Time += zStartBeam;
  }
  Float_t shiftz = pData->Shift(opt);
  
  // Zoom box
  // ----------------------------------------------------------------
  Double_t xRange = (pData->GetXMax(1) - pData->GetXMin(1))/zoom;
  Double_t xMid   = (pData->GetXMax(1) + pData->GetXMin(1))/2.;
  Double_t xMin = xMid - xRange/2.0;
  Double_t xMax = xMid + xRange/2.0;

  pData->SetX2Min(xMin);
  pData->SetX2Max(xMax);

  Double_t yRange = (pData->GetXMax(2) - pData->GetXMin(2))/zoom;
  Double_t yMid   = (pData->GetXMax(2) + pData->GetXMin(2))/2.;
  Double_t yMin = yMid - yRange/2.0;
  Double_t yMax = yMid + yRange/2.0;

  pData->SetX3Min(yMin);
  pData->SetX3Max(yMax);

  Double_t zMin = pData->GetX1Min();
  Double_t zMax = pData->GetX1Max();
  // ----------------------------------------------------------------
  

  // Get charge density histos
  Int_t Nspecies = pData->NSpecies();

  TH3F **hDen3D = new TH3F*[Nspecies];
  for(Int_t i=0;i<Nspecies;i++) {
    hDen3D[i] = NULL;
 
    //    if(i!=index) continue;

    if(!pData->GetChargeFileName(i)) 
      continue;

    cout << Form("Loaded histo from species %i: %s", i, pData->GetSpeciesName(i).c_str()) << endl;
    
    char hName[24];
    sprintf(hName,"hDen3D_%i",i);
    hDen3D[i] = (TH3F*) gROOT->FindObject(hName);
    if(hDen3D[i]) delete hDen3D[i];

    hDen3D[i] = pData->GetH3(pData->GetChargeFileName(i)->c_str(),"charge");
    
  
    hDen3D[i]->SetName(hName);
    hDen3D[i]->GetXaxis()->CenterTitle();
    hDen3D[i]->GetYaxis()->CenterTitle();
    hDen3D[i]->GetZaxis()->CenterTitle();
    hDen3D[i]->GetXaxis()->SetTitle("k_{p} z");
    hDen3D[i]->GetYaxis()->SetTitle("k_{p} x");
    hDen3D[i]->GetZaxis()->SetTitle("k_{p} y");
    
    // -----------------------------------------------------------------
  }


  // Save output to file
    
  TString filename = Form("./%s/Plots/Charge3D/Charge3D-%s_%i.root",sim.Data(),sim.Data(),time);
  // cout << Form("\n Saving snapshot objects to %s",filename.Data()) << endl;
  
  TString f = filename;
  TString dir2 = f.Remove( f.Last( '/' ), f.Length() - f.Last( '/' ) );
  TString dir1 = f.Remove( f.Last( '/' ), f.Length() - f.Last( '/' ) );
      
  gSystem->mkdir( dir1 );
  gSystem->mkdir( dir2 );
      
  TFile *ofile = new TFile(filename,"RECREATE");
  char hName[24];
  for(Int_t i=0;i<Nspecies;i++) {
    if(hDen3D[i]) {
      sprintf(hName,"hDen3D_%i",i);
      hDen3D[i]->Write(hName,TObject::kOverwrite);
    }
  }

  
  // // Dynamic plasma palette
  // PPalette * plasmaPalette = (PPalette*) gROOT->FindObject("plasma");
  // if(!plasmaPalette) {
  //   plasmaPalette = new PPalette("plasma");
  // }
  // plasmaPalette->SetPalette("electron0");


  // // Plotting
  // // -----------------------------------------------

  // // Canvas setup
  // TCanvas *C = new TCanvas("C","3D Charge density",800,500);
  
  // // Output file
  // TString fOutName = Form("./%s/Plots/Charge3D/Charge3D-%s",sim.Data(),pData->GetSpeciesName(index).c_str());
  // fOutName += Form("-%s_%i",sim.Data(),time);
  
  // C->cd();

  // // Setup Pad layout:
  // Float_t lMargin = 0.15;
  // Float_t rMargin = 0.18;
  // Float_t bMargin = 0.20;
  // Float_t tMargin = 0.06;
  // gPad->SetLeftMargin(lMargin);
  // gPad->SetRightMargin(rMargin);
  // gPad->SetBottomMargin(bMargin);
  // gPad->SetTopMargin(tMargin);

  // hDen3D[index]->Draw("glcol");

  // C->cd();

  // // Print to a file
  // PGlobals::imgconv(C,fOutName,opt);
  // // ---------------------------------------------------------

}
