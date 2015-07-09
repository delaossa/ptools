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
#include "PGlobals.hh"
#include "PPalette.hh"
#include "H5Cpp.h"

using namespace std;
using namespace H5;


int main(int argc,char *argv[]) {
  if(argc<=2) {
    printf("\n Usage: %s <simulation name> <-t(time)>\n",argv[0]);
    printf("      <-i(initial time)> <-f(final time)> <-s(time step)>\n");
    printf("      <--png> <--pdf> <--eps> <--units> \n");
    printf("      <--file> <--loop>\n");
    return 0;
  }

  // General options
  TString   sim = "";
  Int_t    time = 0;
  Int_t  iStart = -1;
  Int_t    iEnd = -1;
  Int_t   iStep = 1;
  TString   opt = "";
  
  // Options for Spectrum
  Float_t Pmin =  99999.;
  Float_t Pmax = -99999.;

  // Interfacing command line:
  for(int l=1;l<argc;l++){
    TString arg = argv[l];

    if(arg.Contains("--pdf")) {
      opt += "pdf";
    } else if(arg.Contains("--eps")) {
      opt += "eps";
    } else if(arg.Contains("--png")) {
      opt += "png";
    } else if(arg.Contains("--tiff")) {
      opt += "tiff";
    } else if(arg.Contains("--comov")){
      opt += "comov";
    } else if(arg.Contains("--units")){
      opt += "units";
    } else if(arg.Contains("--center")){
      opt += "center";
    } else if(arg.Contains("--grid")){
      opt += "grid"; 
    } else if(arg.Contains("--logz")){
      opt += "logz"; 
    } else if(arg.Contains("--autop")){
      opt += "autop"; 
    } else if(arg.Contains("--auto")){
      opt += "auto"; 
    } else if(arg.Contains("--loop")){
      opt += "loop"; 
    } else if(arg.Contains("--file")){
      opt += "file"; 
    } else if(arg.Contains("--notext")){
      opt += "notext"; 
    } else if(arg.Contains("--noinfo")){
      opt += "noinfo"; 
    } else if(arg.Contains("-t")) {
      char ss[2];
      sscanf(arg,"%2s%i",ss,&time);
    } else if(arg.Contains("-i")) {
      char ss[2];
      sscanf(arg,"%2s%i",ss,&iStart);
    } else if(arg.Contains("-f")) {
      char ss[2];
      sscanf(arg,"%2s%i",ss,&iEnd);
    } else if(arg.Contains("-s")) {
      char ss[2];
      sscanf(arg,"%2s%i",ss,&iStep);
    } else if(arg.Contains("-pmin")) {
      char ss[5];
      sscanf(arg,"%5s%f",ss,&Pmin);
    } else if(arg.Contains("-pmax")) {
      char ss[5];
      sscanf(arg,"%5s%f",ss,&Pmax);
    } else if( !(arg.Contains("pitz") || arg.Contains("flash") || arg.Contains("regae") || arg.Contains("Gauss") || arg.Contains("pwfa") || arg.Contains("facet") || arg.Contains("FACET") || arg.Contains("rake")) ) {
      cout << Form("\t Invalid argument (%i): exiting...\n",l) << endl;
      return 0;
    } else {
      sim = arg;
    }
  }
  

  PGlobals::Initialize();

  // Palettes!
  gROOT->Macro("PPalettes.C");

  if(opt.Contains("grid")) {
    gStyle->SetPadGridX(1);
    gStyle->SetPadGridY(1);
  }

  // Load PData
  PData *pData = PData::Get(sim.Data());

  if(iStart<0) iStart = time;
  if(iEnd<=iStart) iEnd = iStart;
  
  // Some plasma constants
  Float_t n0 = pData->GetPlasmaDensity();
  Float_t kp = pData->GetPlasmaK();
  Float_t skindepth = 1.;
  if(kp!=0.0) skindepth = 1/kp;
  // Float_t E0 = pData->GetPlasmaE0();

  // z start of the plasma in normalized units.
  Float_t zStartPlasma = pData->GetPlasmaStart()*kp;
  // z start of the beam in normalized units.
  Float_t zStartBeam = pData->GetBeamStart()*kp;

  opt += "comovcenter";
  
  // Time looper
  for(Int_t i=iStart; i<iEnd+1; i+=iStep) {

    time = i;
    pData->LoadFileNames(time);    
    
    if(!pData->IsInit()) continue;
    
    // Time in OU
    Float_t Time = pData->GetRealTime();

    // Centering time and z position:
    Float_t shiftz = pData->Shift(opt);

    if(opt.Contains("center")) {
      Time -= zStartPlasma;
      if(opt.Contains("comov"))      // Centers on the head of the beam.
	Time += zStartBeam;
    }

    // ----------------------------------------------------------------------------------
  
    Int_t Nspecies = pData->NSpecies();
    Float_t **ene = new Float_t*[Nspecies];
    Float_t **x1 = new Float_t*[Nspecies];
    Float_t **q = new Float_t*[Nspecies];
    UInt_t *NP = new UInt_t[Nspecies];
    Float_t eneMin = 99999;
    Float_t eneMax = -99999;
    for(Int_t i=0;i<Nspecies;i++) {
      ene[i] = NULL;
      x1[i] = NULL;
      NP[i] = 0;
    
      if(!pData->GetRawFileName(i)) continue;
      NP[i] = pData->GetRawSingleArray(pData->GetRawFileName(i)->c_str(),&x1[i],"x1");
      NP[i] = pData->GetRawSingleArray(pData->GetRawFileName(i)->c_str(),&ene[i],"ene");
      NP[i] = pData->GetRawSingleArray(pData->GetRawFileName(i)->c_str(),&q[i],"q");
    
      // for(UInt_t j=0;j<1000;j++)
      // 	cout << Form(" %6f  %6f  %6f",x1[i][j],ene[i][j],q[i][j]) << endl;
    
      for(UInt_t j=0;j<NP[i];j++) {
      	if(ene[i][j]<eneMin) eneMin = ene[i][j];
      	if(ene[i][j]>eneMax) eneMax = ene[i][j];
      }

      
    }    

    Float_t x1Min = pData->GetX1Min() - shiftz;
    Float_t x1Max = pData->GetX1Max() - shiftz;
    Int_t Nx1Bin  = pData->GetX1N();

    Float_t eMin = eneMin - (eneMax-eneMin)*0.1;
    Float_t eMax = eneMax + (eneMax-eneMin)*0.1;
    Int_t NeBin  = Nx1Bin;

    char hName[24];
    TH2F **hEneVsZ  = new TH2F*[Nspecies];
    for(Int_t i=0;i<Nspecies;i++) {
    
      hEneVsZ[i] = NULL;
      if(ene[i]==NULL) continue;
    
      sprintf(hName,"EneVsZ_%i",i); 
      hEneVsZ[i] = new TH2F(hName,"",Nx1Bin,x1Min,x1Max,NeBin,eMin,eMax);

      for(UInt_t j=0;j<NP[i];j++) {
	hEneVsZ[i]->Fill(x1[i][j]-shiftz,ene[i][j],-q[i][j]);
      }
      
    }

    TCanvas *C = new TCanvas("C");
    C->SetFillStyle(4000);

    gPad->SetLogz(1);

    TH2F *hEneVsZjoint = (TH2F*) hEneVsZ[1]->Clone("hEneVsZjoint");
    hEneVsZjoint->Add(hEneVsZ[2]);
  
    TExec *exElec   = new TExec("exElec","redelectronPalette->cd();");
  
    exElec->Draw();
  
    hEneVsZjoint->Draw("colz");
  
    // Print to a file
    // Output file
    TString fOutName = Form("./%s/Plots/BeamsEnergy/BeamsEnergy-%s_%i",pData->GetPath().c_str(),pData->GetName(),time);
  
    PGlobals::imgconv(C,fOutName,opt);
    // ---------------------------------------------------------
  

  }
}
