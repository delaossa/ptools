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

  // Electron density species 4 (if any)
  sprintf(hName,"hDen2D_3"); 
  TH2F *hDen2D_3 = (TH2F*) ifile->Get(hName);
  if(hDen2D_3)
    hDen2D->Add(hDen2D_3,1);



  
  
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

  Float_t kzmin = 0.;
  Float_t kzmax = (TMath::TwoPi() / zrange);
  
  Float_t kxmin = 0.;
  Float_t kxmax = (TMath::TwoPi() / xrange);


  // cout << Form(" Creating 2D histos..." );
  // cout << Form(" done!" ) << endl;

  cout << Form(" Allocating data array ..." );

  Int_t dims[2] = {NBinsZ,NBinsX};
  //  Double_t *data = new Double_t[NBinsZ*NBinsX];
  Double_t *data = new Double_t[NBinsZ*2*(NBinsX/2+1)];
  // Extra padding!
  // See http://www.fftw.org/fftw3_doc/Multi_002dDimensional-DFTs-of-Real-Data.html#Multi_002dDimensional-DFTs-of-Real-Data 
  
  cout << Form(" done!" ) << endl;

  cout << Form(" Filling data aray ..." );
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
  TVirtualFFT *fft = TVirtualFFT::FFT(2, dims, "R2C ES K");
  fft->SetPoints(data);
  fft->Transform();  
  cout << Form("   done!" ) << endl;

  TH2F *hFFTRE = 0;
  hFFTRE = (TH2F*) TH1::TransformHisto(fft, hFFTRE, "RE");
  hFFTRE->SetName("hFFTRE");
  hFFTRE->Scale(1.0/TMath::Sqrt(NBinsZ*NBinsX));
  hFFTRE->SetBins(NBinsZ,kzmin,kzmax,NBinsX,kxmin,kxmax);
  
  TH2F *hFFTIM = 0;
  hFFTIM = (TH2F*) TH1::TransformHisto(fft, hFFTIM, "IM");
  hFFTIM->SetName("hFFTIM");
  hFFTIM->SetBins(NBinsZ,kzmin,kzmax,NBinsX,kxmin,kxmax);
  hFFTIM->GetZaxis()->SetRangeUser(-1,1);
  
  cout << Form(" Solving equation for potential in fourier space" )  << endl;
  // Equation: \phi(kz,kx) = \rho(kz,kx)/(kz^2 + kx^2) 
  
  fft->GetPoints(data);    
  for(Int_t i=0; i<dims[0]; i++) {
    for(Int_t j=0; j<(dims[1]/2+1); j++) {

      Int_t index =  2 * (i * (dims[1]/2+1)  + j);
      //      cout << Form(" i = %4i  j = %4i  index = %10i",i,j,index) << endl;
      
      Float_t kz;
      if(i<dims[0]/2)
	kz = kzmax * (i);
      else
	kz =  kzmax * (i-dims[0]);
      Float_t kx = (kxmax) * (j) ;
      
      Float_t k2 = kx*kx + kz*kz;
      if(k2==0) {
	data[index] = 0.0;
	data[index+1] = 0.0; 
      } else {
	data[index] /= k2;
	data[index+1] /= k2;      
      }
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
  
  cout << Form(" Take the gradient of the potential ..." );
  TH2F *hPhi2D = new TH2F("hPhi2D","",NBinsZ,zmin,zmax,NBinsX,xmin,xmax);
  TH2F *hE2D_1_ifft = new TH2F("hE2D_1_ifft","",NBinsZ,zmin,zmax,NBinsX,xmin,xmax);
  Double_t dx = hPhi2D->GetYaxis()->GetBinWidth(1);
  for(Int_t i=0; i<NBinsZ; i++) {
    for(Int_t j=0; j<NBinsX; j++) {
      Int_t index =  i * NBinsX + j;
      
      hPhi2D->SetBinContent(i+1,j+1,re_back[index]);
      
      Float_t der = 0.;
      if(j>2 && j<NBinsX-2) {
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
  TH1D *hE1D_1_ifft = hE2D_1_ifft->ProjectionX("hE1D_1_ifft",hE2D_1_ifft->GetNbinsY()/2,hE2D_1_ifft->GetNbinsY()/2);


  TH2F *hE2D_0_ifft = new TH2F("hE2D_0_ifft","",NBinsZ,zmin,zmax,NBinsX,xmin,xmax);
  Double_t dz = hPhi2D->GetXaxis()->GetBinWidth(1);
  for(Int_t j=0; j<NBinsX; j++) {
    for(Int_t i=0; i<NBinsZ; i++) {
      Double_t der = 0.;
      if(i>2 && i<NBinsZ-2) {
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
  TH1D *hE1D_0_ifft = hE2D_0_ifft->ProjectionX("hE1D_0_ifft",hE2D_0_ifft->GetNbinsY()/2,hE2D_0_ifft->GetNbinsY()/2);

  cout << Form("   done!" ) << endl;
    
  // TH2F *hPhi2D = 0;
  // hPhi2D = (TH2F*) TH1::TransformHisto(fft_back, hPhi2D, "RE");
  // hPhi2D->SetName("hPhi2D");

  
  // Extract radiation from the total electric field
  // by substracting the electrocstatic fields
  // Get Ex from simulation
  sprintf(hName,"hE2D_1"); 
  TH2F *hE2D_1 = (TH2F*) ifile->Get(hName);
  sprintf(hName,"hE1D_1"); 
  TH1F *hE1D_1 = (TH1F*) ifile->Get(hName);
  
  TH2F *hE2D_1_rad = new TH2F("hE2D_1_rad","",NBinsZ,zmin,zmax,NBinsX,xmin,xmax);
  for(Int_t i=0; i<NBinsZ; i++) {
    for(Int_t j=0; j<NBinsX; j++) {
      Float_t content = hE2D_1->GetBinContent(i+1,j+1) - hE2D_1_ifft->GetBinContent(i+1,j+1) ;
      hE2D_1_rad->SetBinContent(i+1,j+1,fabs(content));     
    }
  }
  
  TH1D *hE1D_1_rad = hE2D_1_rad->ProjectionX("hE1D_1_rad",hE2D_1_rad->GetNbinsY()/2,hE2D_1_rad->GetNbinsY()/2);
  
  //  Float_t omega0 = 417.368;
  Float_t omega0 = 20.87;

  // Vector potential
  TH2F *ha2D = (TH2F*) hE2D_1_rad->Clone("ha2D");
  ha2D->Scale(1.0/omega0);

  TH1F *ha1D = (TH1F*) hE1D_1_rad->Clone("ha1D");
  ha1D->Scale(1.0/omega0);
  
  ha2D->GetZaxis()->SetRangeUser(0.0101,ha2D->GetMaximum());
  fieldTPalette->cd();
  ha2D->Draw("colz");
}
