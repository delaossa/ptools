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
#include <TRandom3.h>

#include "PData.hh"
#include "PGlobals.hh"
#include "PPalette.hh"

void PlotTracks( const TString &sim, Int_t index = 0, const TString &options="") {
  
#ifdef __CINT__  
  gSystem->Load("libptools.so");
#endif

  PGlobals::Initialize();

  TString opt = options;
 
  if(opt.Contains("grid")) {
    gStyle->SetPadGridX(1);
    gStyle->SetPadGridY(1);
  }

  // Load PData
  PData *pData = PData::Get(sim.Data());
  pData->LoadFileNames(0);
  pData->PrintData();
  string *trackfile = pData->GetTrackFileName(index);
  //  cout << "\n EEEEEEEEEEEEEEE \n" << endl;
  // if(trackfile != NULL) {
  //   cout << Form("\n Track file name : %s \n",trackfile->c_str()) << endl;
  // } else {
  //   return 0;
  // }
  
  Int_t NTracks = 1;
  if(trackfile != NULL) 
    TTree **tracktree = pData->GetTrackTree(trackfile->c_str(),NTracks);
  
  TTree *rawtree = pData->GetRawTree(pData->GetRawFileName(index)->c_str());
  
  
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
  
  // opt += "comovcenter";

  // Centering time and z position:
  Double_t shiftz = 0.0;
  if(opt.Contains("center")) {
    Time -= zStartPlasma;
    if(opt.Contains("comov")) {     // Centers on the head of the beam.
      Time += zStartBeam;
      shiftz += zStartBeam;
    } else {
      shiftz += zStartPlasma;
    }
  } 
  if(opt.Contains("comov")) {
    Double_t v = pData->GetBeamVelocity();    
    if(v==0) v = 1.0; // If equals to 0 (default), then set to c.
    shiftz += v * pData->GetRealTime();
  }   
  
  // Pointer to data TTree
  // TTree *tree = pData->GetTreeRaw(pData->GetRawFileName(index)->c_str(),opt);
  
}
