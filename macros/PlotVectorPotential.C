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


Double_t *Derivative(Double_t a[], Int_t n)
{
    Double_t firstDer(0.0);
    Double_t t(0.001);
    Double_t *dx = new Double_t[n];
    
    for (Int_t i=0;i<n;i++) {
      if(i<2 || i>=n-2)  dx[i] = 0.0;
      dx[i] = (   4.0 / 3.0 * (a[i + 1] - a[i - 1]) / (2.0 * t)
		  - 1.0 / 3.0 * (a[i + 2] - a[i - 2]) / (4.0 * t) );
    }
    return dx;
}

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

  sprintf(hName,"hDen2D_1"); 
  TH2F *hDen2D_1 = (TH2F*) ifile->Get(hName);
  if(hDen2D_1)
    hDen2D->Add(hDen2D_1,1);

  sprintf(hName,"hDen2D_2"); 
  TH2F *hDen2D_2 = (TH2F*) ifile->Get(hName);
  if(hDen2D_2)
    hDen2D->Add(hDen2D_2,1);

  sprintf(hName,"hE2D_1"); 
  TH2F *hEx2D = (TH2F*) ifile->Get(hName);
  TH2F *hEx2Dsta = (TH2F*) hEx2D->Clone("hEx2Dsta");
  hEx2Dsta->Reset();
  
  // Get the sliced 1D histograms 
  Int_t NBinsZ = hDen2D->GetXaxis()->GetNbins();
  Int_t NBinsX = hDen2D->GetYaxis()->GetNbins();
  Double_t xmin = hDen2D->GetYaxis()->GetBinLowEdge(1);
  Double_t xmax = hDen2D->GetYaxis()->GetBinUpEdge(NBinsX);
  Double_t xrange = xmax-xmin;


  TH1D **hDen1D = new TH1D*[NBinsZ]; 
  TH1  **hDen1Dt = new TH1*[NBinsZ];
  TH1  **hPhi1D  = new TH1*[NBinsZ];
  TH1  **hEx1D   = new TH1*[NBinsZ];
  TH1  **hDen1Db = new TH1*[NBinsZ];
  TH1D *hDen1Dion = 0;

  Double_t kmin, kmax;
  
  //  TVirtualFFT::SetTransform(0);
  for(Int_t i=NBinsZ; i>0; i--) {
    TVirtualFFT::SetTransform(0);
    // cout << Form(" Bin = %i",i) << endl;

    // Get 1D histograms for each longitudinal slice
    hDen1D[i-1] = hDen2D->ProjectionY("",i,i); 
    sprintf(hName,"hDen1D_x_%i",i); 
    hDen1D[i-1]->SetName(hName);

    // Get ion distribution
    if(i==NBinsZ) {
      hDen1Dion = (TH1D*) hDen1D[i-1]->Clone("hDen1Dion");      
    }

    // and substract it to the electron density 
    hDen1D[i-1]->Add(hDen1Dion,-1);
    
    // Get Transform
    hDen1Dt[i-1] = 0;
    hDen1Dt[i-1] = hDen1D[i-1]->FFT(hDen1Dt[i-1],"Re");
    sprintf(hName,"hDen1Dt_x_%i",i); 
    hDen1Dt[i-1]->SetName(hName);

    TVirtualFFT *fft = TVirtualFFT::GetCurrentTransform();
    Int_t n = NBinsX;

    // Set the frequency range
    if(i==NBinsZ) {
      kmin = 0.0;
      kmax = TMath::TwoPi() * hDen1Dt[i-1]->GetNbinsX()/xrange;
    }    
    hDen1Dt[i-1]->SetBins(NBinsX,kmin,kmax);
    hDen1Dt[i-1]->Scale(TMath::Sqrt(1.0/NBinsX));

    Double_t *re_full = new Double_t[n];
    Double_t *im_full = new Double_t[n];
    fft->GetPointsComplex(re_full,im_full);

    
    for(Int_t j=0;j<n;j++) {
      //Double_t k = (j*(n-1.0)/n - (n/2.0)) * (TMath::TwoPi()/xrange);
      //Double_t k = (j*(n-1.0)/n - (n/2.0)) * (TMath::Pi()/xrange);
      Double_t k = hDen1Dt[i-1]->GetBinCenter(j+1);
      //Double_t k = ((j+1)/xrange);
      Double_t k2 = k * k;
      // if(i==NBinsZ)
      //  	cout << Form(" k = %.4e  k2 = %.4e",k,k2) << endl;
      re_full[j] /= k2;
      im_full[j] /= k2;
    }
    
    //Now let's make a backward transform:
    TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &n, "C2R EX K");
    fft_back->SetPointsComplex(re_full,im_full);
    fft_back->Transform();

    hPhi1D[i-1] = 0;
    hPhi1D[i-1] = TH1::TransformHisto(fft_back,hPhi1D[i-1],"RE");
    sprintf(hName,"hPhi1D_x_%i",i); 
    hPhi1D[i-1]->SetName(hName);
    hPhi1D[i-1]->SetBins(NBinsX,xmin,xmax);
    hPhi1D[i-1]->Scale(1.0/NBinsX);    

    // Derivative
    sprintf(hName,"hEx1D_x_%i",i); 
    hEx1D[i-1] = new TH1D(hName,"",NBinsX,xmin,xmax);

    Double_t dt = hPhi1D[i-1]->GetBinWidth(1);
    for(Int_t j=0;j<hPhi1D[i-1]->GetNbinsX();j++) {
      Double_t der = 0.;
      if(j>2 && j<=n-2) {
	
	der =  (   4.0 / 3.0 * (hPhi1D[i-1]->GetBinContent(j+1) - hPhi1D[i-1]->GetBinContent(j-1))/ (2.0 * dt)
		   
		   - 1.0 / 3.0 * (hPhi1D[i-1]->GetBinContent(j+2) - hPhi1D[i-1]->GetBinContent(j-2)) / (4.0 * dt) );
	
      }
      
      hEx1D[i-1]->SetBinContent(j+1,der);
      hEx2Dsta->SetBinContent(i,j+1,der);
    }
    
  }
  
}
