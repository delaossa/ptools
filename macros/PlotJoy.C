#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>

#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TH1F.h>

#include "PData.hh"
#include "PDataHiP.hh"
#include "PGlobals.hh"


void PlotJoy( const TString &sim, Int_t time, Int_t index = 0, Float_t zoom=2, const TString &opt="") {
    
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
  
  // Refresh and Style
  PGlobals::Initialize();

  gStyle->SetCanvasColor(1);
  gStyle->SetFrameFillColor(1);
  gStyle->SetFrameFillColor(1);
  gStyle->SetHistLineColor(0);
  gStyle->SetHistLineWidth(1);
  gStyle->SetPadLeftMargin(0.20);
  
  TFile *ofile = new TFile("joy.root","RECREATE");
  Int_t NbinsX1 = pData->GetX1N();
  TH1F **hDen1D = new TH1F*[NbinsX1];

  Double_t TOP = 100;
  Double_t shift = 0.2;
  Double_t MAX = -999;
  Double_t MAX2 = 2*(NbinsX1 * shift) ;
  
  TCanvas *C = new TCanvas("C","",800,1280);
  
  Int_t j =  (pData->GetX3iMin() + pData->GetX3iMax()) / 2;
  for(Int_t i=pData->GetX1iMin();i<pData->GetX1iMax();i++) {
    
    char hName[24];
    sprintf(hName,"hDen1D_%i",i);
    hDen1D[i] = (TH1F*) gROOT->FindObject(hName);
    if(hDen1D[i]) delete hDen1D[i];
    
    // 1D histograms
    if(pData->Is3D()) {
      hDen1D[i] = pData->GetH1SliceX3D(pData->GetChargeFileName(0)->c_str(),"charge",i,i,j,j,opt+"avg");
    } else {
      return;
    }
    hDen1D[i]->SetName(hName);

    for(Int_t k=1;k<=hDen1D[i]->GetNbinsX();k++) {
      if(hDen1D[i]->GetBinContent(k)<TOP)
	hDen1D[i]->SetBinContent(k,hDen1D[i]->GetBinContent(k) + i*shift);
      else
	hDen1D[i]->SetBinContent(k,TOP + i*shift);

    }

    if(MAX<hDen1D[i]->GetMaximum())
      MAX = hDen1D[i]->GetMaximum();
    
    hDen1D[i]->Write(hName,TObject::kOverwrite);
    
  }
  
  TH1F *hFrame = (TH1F*) hDen1D[0]->Clone("hFrame");
  hFrame->Reset();
  hFrame->GetYaxis()->SetRangeUser(0.0,MAX2);

  hFrame->Draw("AC");
  for(Int_t i=pData->GetX1iMax()-1;i>=pData->GetX1iMin();i -= NbinsX1/50) {
    hDen1D[i]->SetFillColor(1);
    hDen1D[i]->Draw("FC same");
  }
   
  //  hDen1D[NbinsX1/2]->Draw();
  PGlobals::imgconv(C,"joy","pdf");
}
