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


  cout << Form(" Getting histos..." );

  // Electron density (plasma)
  char hName[36];
  sprintf(hName,"hDen2D_0"); 
  TH2F *hDen2D = (TH2F*) ifile->Get(hName);
  
  // Electron density species 2 (if any)
  sprintf(hName,"hDen2D_1"); 
  TH2F *hDen2D_1 = (TH2F*) ifile->Get(hName);
  if(hDen2D_1)
    hDen2D->Add(hDen2D_1,1);

  // Electron density species 3 (if any)
  sprintf(hName,"hDen2D_2"); 
  TH2F *hDen2D_2 = (TH2F*) ifile->Get(hName);
  if(hDen2D_2)
    hDen2D->Add(hDen2D_2,1);

  
  cout << Form(" done!" ) << endl;

  // Get the sliced 1D histograms 
  Int_t NBinsZ = hDen2D->GetXaxis()->GetNbins();
  Float_t zmin = hDen2D->GetXaxis()->GetBinLowEdge(1);
  Float_t zmax = hDen2D->GetXaxis()->GetBinUpEdge(NBinsZ);
  Float_t zrange = zmax-zmin;

  Int_t NBinsX = hDen2D->GetYaxis()->GetNbins();
  Float_t xmin = hDen2D->GetYaxis()->GetBinLowEdge(1);
  Float_t xmax = hDen2D->GetYaxis()->GetBinUpEdge(NBinsX);
  Float_t xrange = xmax-xmin;

  // cout << Form(" Creating 2D histos..." );
  // cout << Form(" done!" ) << endl;

  cout << Form(" Allocating array with real data..." );

  Int_t dims[2] = {NBinsZ,NBinsX};
  Double_t *data = new Double_t[NBinsZ*NBinsX];
  
  // TVirtualFFT::SetTransform(0);
  
  cout << Form(" done!" ) << endl;

  cout << Form(" Filling data aray..." );
  for(Int_t i=0; i<NBinsZ; i++) {
    
    for(Int_t j=0; j<NBinsX; j++) {

      Int_t index =  i * NBinsX + j;
      
      data[index]  = hDen2D->GetBinContent(i+1,j+1);

      // substract ion background
      data[index] -= 1;
    }
  }
  
  cout << Form("   done!" ) << endl;
  
  cout << Form(" Fourier transform ..." );

  // TH2F *hFFTREb = new TH2F("hFFTREb","",NBinsZ,kzmin,kzmax,NBinsX,kxmin,kxmax);
  // TH2F *hFFTIMb = new TH2F("hFFTIMb","",NBinsZ,kzmin,kzmax,NBinsX,kxmin,kxmax);
  // TH2F *hFFTREb = new TH2F("hFFTREb","",NBinsZ,0.,NBinsZ,NBinsX,0.,NBinsX);
  // TH2F *hFFTIMb = new TH2F("hFFTIMb","",NBinsZ,0.,NBinsZ,NBinsX,0.,NBinsX);
  
  // TH2F *hFFTRE = 0;
  // hFFTRE = (TH2F*) TH1::TransformHisto(fft, hFFTRE, "RE");
  // hFFTRE->SetName("hFFTRE");
  // TH2F *hFFTIM = 0;
  // hFFTIM = (TH2F*) TH1::TransformHisto(fft, hFFTIM, "IM");
  // hFFTIM->SetName("hFFTIM");

  TVirtualFFT *fft = TVirtualFFT::FFT(2, dims, "R2C ES K");
  fft->SetPoints(data);
  fft->Transform();  
  cout << Form("   done!" ) << endl;
  
   cout << Form(" Solving equation for potential in fourier space" )  << endl;
  // \phi(kz,kx) = 1/(kz^2 + kx^2) 
  
  fft->GetPoints(data);    
  for(Int_t i=0; i<dims[0]; i++) {
    for(Int_t j=0; j<(dims[1]/2+1); j++) {

      Int_t index =  2 * (i * (dims[1]/2+1)  + j);
      
      Float_t kz;
      if(i<dims[0]/2)
	kz = TMath::Pi() * (i+0.5) / zrange;
      else
	kz = -TMath::Pi() * (dims[0]-i+0.5) / zrange;     
      Float_t kx = TMath::TwoPi() * (j+0.5) / xrange;
      
      Float_t k2 = kx*kx + kz*kz;
      
      data[index] /= k2;
      data[index+1] /= k2;
      
    }
    
  }
  cout << Form("   done!" ) << endl;
  
  cout << Form(" Inverse Fourier transform ..." );
  // backward transform:
  TVirtualFFT *fft_back = TVirtualFFT::FFT(2, dims, "C2R ES K");
  fft_back->SetPoints(data);
  fft_back->Transform();
  
  cout << Form(" done!" ) << endl;

  Double_t *re_back = fft_back->GetPointsReal();
  
  TH2F *hPhi2D = new TH2F("hPhi2D","",NBinsZ,zmin,zmax,NBinsX,xmin,xmax);
  TH2F *hE2D_1_ifft = new TH2F("hE2D_1_ifft","",NBinsZ,zmin,zmax,NBinsX,xmin,xmax);
  Double_t dx = hPhi2D->GetYaxis()->GetBinWidth(1);
  for(Int_t i=0; i<NBinsZ; i++) {
    for(Int_t j=0; j<NBinsX; j++) {
      Int_t index =  i * NBinsX + j;
      
      hPhi2D->SetBinContent(i+1,j+1,re_back[index]);
      
      Double_t der = 0.;
      if(j>2 && j<=NBinsX-2) {
	Int_t indjp1 =  i * NBinsX + (j+1);
	Int_t indjp2 =  i * NBinsX + (j+2);
	Int_t indjm1 =  i * NBinsX + (j-1);
	Int_t indjm2 =  i * NBinsX + (j-2);
	
	der =  ( 4.0 / 3.0 * (re_back[indjp1] - re_back[indjm1]) / (2.0 * dx)
	       - 1.0 / 3.0 * (re_back[indjp2] - re_back[indjm2]) / (4.0 * dx) );
	
      }
      
      hE2D_1_ifft->SetBinContent(i+1,j+1,der);
    }
  }
  hPhi2D->Scale(1.0/(NBinsZ*NBinsX));    
  hE2D_1_ifft->Scale(1.0/(NBinsZ*NBinsX));


  TH2F *hE2D_0_ifft = new TH2F("hE2D_0_ifft","",NBinsZ,zmin,zmax,NBinsX,xmin,xmax);
  Double_t dz = hPhi2D->GetXaxis()->GetBinWidth(1);
  for(Int_t j=0; j<NBinsX; j++) {
    for(Int_t i=0; i<NBinsZ; i++) {
      Double_t der = 0.;
      if(i>2 && j<=NBinsZ-2) {
	Int_t indip1 =  (i+1) * NBinsX + j;
	Int_t indip2 =  (i+2) * NBinsX + j;
	Int_t indim1 =  (i-1) * NBinsX + j;
	Int_t indim2 =  (i-2) * NBinsX + j;
	
	der =  ( 4.0 / 3.0 * (re_back[indip1] - re_back[indim1]) / (2.0 * dz)
		 - 1.0 / 3.0 * (re_back[indip2] - re_back[indim2]) / (4.0 * dz) );
	
      }
      
      hE2D_0_ifft->SetBinContent(i+1,j+1,der);
    }
  }
  hE2D_0_ifft->Scale(1.0/(NBinsZ*NBinsX));
  
  
  // TH2F *hPhi2D = 0;
  // hPhi2D = (TH2F*) TH1::TransformHisto(fft_back, hPhi2D, "RE");
  // hPhi2D->SetName("hPhi2D");
    
    
}
