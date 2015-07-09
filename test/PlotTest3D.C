#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>

#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TH3.h>
#include <TExec.h>

#include "PData.hh"
#include "PlasmaGlob.hh"
#include "PPalette.hh"

void PlotTest3D(const TString &sim, Int_t time, Float_t zoom=2, Int_t Nbins=2, const TString &options="") {
  
#ifdef __CINT__  
  gSystem->Load("libplasma.so");
#endif

  PlasmaGlob::Initialize();

  // Palettes!
  gROOT->Macro("PlasmaPalettes.C");

  TH3F *h3 = new TH3F("h3","",100,-1.,1.,100,-1.,1.,100,-1.,1.);
  h3->FillRandom("gaus", 100000);

  TCanvas *C = new TCanvas("C","",1024,576);

  TExec *exDriver   = new TExec("exDriver","redelectronPalette->cd();");  

  C->cd();

  h3->Draw("iso");
}
