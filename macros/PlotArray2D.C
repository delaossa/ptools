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
#include "PGlobals.hh"
#include "PPalette.hh"

void PlotArray2D( const TString &sim, Int_t time, Float_t zoom=1, Int_t Nbins=1, const TString &options="") {
    
#ifdef __CINT__  
  gSystem->Load("libptools.so");
#endif
  
  // Load PData
  PData *pData = PData::Get(sim.Data());
  pData->LoadFileNames(time);
  if(!pData->IsInit()) return;

  // Refresh and Style
  PGlobals::Initialize();

  // Coloured palettes
  gROOT->Macro("PPalettes.C");

  TString opt = options;
 
  // More makeup
  if(opt.Contains("grid")) {
    gStyle->SetPadGridX(1);
    gStyle->SetPadGridY(1);
  }

  // Some plasma constants
  Double_t n0 = pData->GetPlasmaDensity();
  Double_t kp = pData->GetPlasmaK();
  Double_t skindepth = 1.;
  if(kp!=0.0) skindepth = 1/kp;
  Double_t E0 = pData->GetPlasmaE0();

  // Some beam properties:
  Float_t Ebeam = pData->GetBeamEnergy() * PUnits::MeV;
  Float_t gamma = Ebeam / PConst::ElectronMassE;
  Float_t vbeam = TMath::Sqrt(1 - 1/(gamma*gamma));
  // cout << Form(" - Bunch gamma      = %8.4f", gamma ) << endl;
  // cout << Form(" - Bunch velocity   = %8.4f c", vbeam ) << endl;
  Float_t nb = pData->GetBeamDensity();
  
  // Time in OU
  Float_t Time = pData->GetRealTime();
  // z start of the plasma in normalized units.
  Float_t zStartPlasma = pData->GetPlasmaStart()*kp;
  // z start of the beam in normalized units.
  Float_t zStartBeam = pData->GetBeamStart()*kp;
  // z start of the neutral in normalized units.
  Float_t zStartNeutral = pData->GetNeutralStart()*kp;
  // z end of the neutral in normalized units.
  Float_t zEndNeutral = pData->GetNeutralEnd()*kp;
  
  if(opt.Contains("center")) {
    Time -= zStartPlasma;
    if(opt.Contains("comov"))      // Centers on the head of the beam.
      Time += zStartBeam;
  }
  Float_t shiftz = pData->Shift(opt);
  //  cout << "Shift = " << shiftz << endl;
  
  // Calculate the "axis range" in number of bins. 
  // If Nbins==0 a RMS width is taken.
  Double_t rms0 = pData->GetBeamRmsY() * kp;
  if(pData->IsCyl())  rms0  = pData->GetBeamRmsR() * kp;
  if(Nbins==0) { 
    Nbins =  TMath::Nint(rms0 / pData->GetDX(1));
  }
  // -------------------------------------------------------------------------------

  // Zoom window:
  
  // Take just the middle bin in 3rd dimension
  Float_t x3Min = (pData->GetXMax(2) + pData->GetXMin(2))/2.;
  Float_t x3Max = x3Min + pData->GetDX(2);
  pData->SetX3Min(x3Min);
  pData->SetX3Max(x3Max);
    
  // Perform a center zoom in 2nd dimension.
  Float_t x2Range    = (pData->GetXMax(1) - pData->GetXMin(1))/zoom;
  Float_t midPoint  = (pData->GetXMax(1) + pData->GetXMin(1))/2.;
  Float_t x2Min = midPoint - x2Range/2.0;
  Float_t x2Max = midPoint + x2Range/2.0;
  pData->SetX2Min(x2Min);
  pData->SetX2Max(x2Max);
  
  // keep the whole 1st dimension
  // Float_t x1Min = pData->GetXMin(0);
  // Float_t x1Max = pData->GetXMax(0);
  // pData->SetX1Min(x1Min);
  // pData->SetX1Max(x1Max);
  
  cout << Form(" Plotting range:  %.2f < x1 < %.2f ,  %.2f < x2 < %.2f ,  %.2f < x2 < %.2f",
	       pData->GetX1Min(),pData->GetX1Max(),pData->GetX2Min(),pData->GetX2Max(),pData->GetX3Min(),pData->GetX3Max()) << endl;
  cout << Form("                  %3i [%3i,%3i] , %3i [%3i,%3i] , %3i [%3i,%3i]",
	       pData->GetX1N(),pData->GetX1iMin(),pData->GetX1iMax(),pData->GetX2N(),pData->GetX2iMin(),pData->GetX2iMax(),pData->GetX3N(),pData->GetX3iMin(),pData->GetX3iMax()) << endl;
  
  // ----------------------------------------------------------------------------------
  
  
  // Get array
  UInt_t  dim[3];
  Float_t *data = NULL;
  if(pData->GetChargeFileName(0)) {
    data = pData->Get3Darray(pData->GetChargeFileName(0)->c_str(),"charge",dim);
  }
  
  cout << Form("%i %i %i", dim[0],dim[1],dim[2]) << endl;
  
  // if(data) {
  //   for(UInt_t i=0;i<dim[0];i++)
  //     for(UInt_t j=0;j<dim[1];j++)
  // 	for(UInt_t k=0;k<dim[2];k++)
  // 	  cout << Form("%i, %i, %i : %.2f",i,j,k,data[i * dim[2] * dim[1] + j * dim[2] + k]) << endl; }
  
  TH2F *h2D = new TH2F("h2D","",dim[2],pData->GetX1Min(),pData->GetX1Max(),dim[1],pData->GetX2Min(),pData->GetX2Max());
    
  if(data) {
    for(UInt_t j=0;j<dim[1];j++)
      for(UInt_t k=0;k<dim[2];k++) {
	h2D->SetBinContent(k+1,j+1,-data[ (dim[0]/2) * dim[2] * dim[1] + j * dim[2] + k] );
      }
  }
  
  h2D->Draw("colz");
  
}
