#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>

#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>

#include <TH1D.h>
#include <TVirtualFFT.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TMath.h>

#include "PData.hh"
#include "PDataHiP.hh"
#include "PGlobals.hh"
#include "PPalette.hh"

void PlotVectorPotential( const TString &sim, Int_t timestep, const TString &options="") {
  
#ifdef __CINT__  
  gSystem->Load("libptools.so");
#endif

  PData *pData = PData::Get(sim.Data());
  if(pData->isHiPACE()) {
    delete pData; pData = NULL;
    pData = PDataHiP::Get(sim.Data());
  }
  
  pData->LoadFileNames(timestep);
  if(!pData->IsInit()) return;

  PGlobals::Initialize();

  TString opt = options;

  // Open snapshot file and get the histograms
  TString filename;
  filename = Form("./%s/Plots/Snapshots/Snapshot-%s_%i.root",sim.Data(),sim.Data(),timestep);
  
  TFile  *ifile = (TFile*) gROOT->GetListOfFiles()->FindObject(filename.Data());
  if (!ifile) ifile = new TFile(filename,"READ");

  ifile->cd();


  // Time in OU
  Double_t Time = pData->GetRealTime();

  char hName[36];
  sprintf(hName,"hDen2D_0"); 
  TH2F *hDen2D = (TH2F*) ifile->Get(hName);
  if(hDen2D) {
    if(hDen2D->GetMaximum() < 1E-4) {
      delete hDen2D;
      hDen2D = NULL;
    }
  }
  
  // Get the sliced 1D histograms
  Int_t NBinsX = hDen2D->GetXaxis()->GetNbins();
  TH1D **hDen1D = new TH1D*[NBinsX]; 
  TH1 **hDen1Dt = new TH1*[NBinsX];
  TVirtualFFT::SetTransform(0);
  for(Int_t i=0; i<hDen2D->GetYaxis()->GetNbins(); i++) {
    hDen1D[i] = hDen2D->ProjectionY("",i,i); 
    sprintf(hName,"hDen1D_x_%i",i); 
    hDen1D[i]->SetName(hName);
    hDen1Dt[i] = 0;
    hDen1Dt[i] = hDen1D[i]->FFT(hDen1Dt[i],"MAG");
    
  }
  
}
