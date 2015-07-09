
#ifndef __CINT__
#include "PFunctions.hh"
#endif

#include "PData.hh"
#include "PlasmaGlob.hh"

void testZOOM(const TString &sim, Int_t time, const TString &opt=""){

#ifdef __CINT__
  gSystem->Load("libplasma.so");
#endif
  PlasmaGlob::Initialize();  

  // Palettes!
  gROOT->Macro("PlasmaPalettes.C");
  
  // Init Units table
  PUnits::UnitsTable::Get();
  
  // Load PData
  PData *pData = PData::Get(sim.Data());
  pData->LoadFileNames(time);
  if(!pData->IsInit()) return;

  pData->SetX1Min(130);
  pData->SetX1Max(135);
  pData->SetX2Min(-5);
  pData->SetX2Max(5);
  
  opt += "avg";
  TH2F *h2 = pData->GetH2SliceZYplus(pData->GetChargeFileName(0)->c_str(),"charge",-1,1,opt);
  h2->Print();
  
  TCanvas *C = new TCanvas("C","TEST",800,600);

  gPad->SetLogz(1);
  h2->Draw("colz");
}
