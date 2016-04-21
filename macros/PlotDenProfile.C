#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>

#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TMath.h>
#include <TPaveText.h>
#include <TBox.h>

#include "PData.hh"
#include "PGlobals.hh"
#include "PPalette.hh"


Double_t DensityGauss(Double_t *x, Double_t *par) {
  Double_t f;

  Double_t n0 = par[0];
  Double_t z0 = par[1];
  Double_t sigma0 = par[2];
  Double_t n1 = par[3];
  Double_t z1 = par[4];
  Double_t sigma1 = par[5];
  Double_t n2 = par[6];
  Double_t z2 = par[7];
  Double_t sigma2 = par[8];
  
  if( x[0] < z0) f = n0 * TMath::Gaus(x[0],z0,sigma0);
  else if( x[0] < z1) f = n0;
  else if( x[0] >= z1 && x[0]<z2) f = (n0-n1) * TMath::Gaus(x[0],z1,sigma1) + n1;
  else if( x[0] >= z2) f = (n1-n2) *  TMath::Gaus(x[0],z2,sigma2) + n2;
  
  return f;
}

Double_t DensityGaussSum(Double_t *x, Double_t *par) {
  Double_t f = DensityGauss(x,par) +  DensityGauss(x,&par[9]);

  
  return f;
}


void PlotDenProfile(const TString &sim, const Float_t zmin0 = 0, const Float_t zmax0 = 100, const TString &options="pdf") { 
  
#ifdef __CINT__  
  gSystem->Load("libptools.so");
#endif

  PGlobals::Initialize();
  

  TString opt = options;

  // More makeup            
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  gStyle->SetNumberContours(255);

  if(opt.Contains("grid")) {
    gStyle->SetPadGridX(1);
    gStyle->SetPadGridY(1);
  }


  Int_t myBlue = TColor::GetColor((Float_t) 0.16, (Float_t) 0.83, (Float_t) 0.5);
  Int_t myNaranja = TColor::GetColor((Float_t) 0.992157, (Float_t) 0.411765, (Float_t) 0.027451);

  // Load first simulation data (for instance)
  PData *pData = PData::Get(sim.Data());
  Double_t np  = pData->GetPlasmaDensity() / (1E17/PUnits::cm3);
  Double_t skd = pData->GetPlasmaSkinDepth() / PUnits::mm;

  Double_t shiftz = 60.0;
  Double_t  zmin = zmin0 - shiftz; 
  Double_t  zmax = zmax0 - shiftz; 

  if(opt.Contains("units")) {
    zmin *= skd;
    zmax *= skd;
  }
  
  cout << Form("  Plasma density = %.2f E17 cm^{-3}", np) << endl;
  
  Double_t n0 = 0.1;
  Double_t z0 = 34.85 - shiftz;
  Double_t sigma0 = 0.1;
  Double_t n1 = 0.0;
  Double_t z1 = 47.42 - shiftz;
  Double_t sigma1 = 0.1;
  Double_t n2 = 0;
  Double_t z2 = 1000 - shiftz;
  Double_t sigma2 = 13.31;

  if(opt.Contains("units")) {
    n0 *= np;
    z0 *= skd;
    sigma0 *= skd;
    n1 *= np;
    z1 *= skd;
    sigma1 *= skd;
    n2 *= np;
    z2 *= skd;
    sigma2 *= skd;
  }

  TF1 *fDenHe = new TF1("fDenHe",DensityGauss,zmin,zmax,9);
  fDenHe->SetParameters(n0,z0,sigma0,n1,z1,sigma1,n2,z2,sigma2);
  fDenHe->SetNpx(10000);

  Double_t n0p = 10.0;
  Double_t z0p = 60 - shiftz;
  Double_t sigma0p = 1.88;
  Double_t n1p = 1.0;
  Double_t z1p = 60 - shiftz;
  Double_t sigma1p = 2.5;
  Double_t n2p = 1.0;
  Double_t z2p = 10000 - shiftz;
  Double_t sigma2p = 0;

  if(opt.Contains("units")) {    
    n0p *= np;
    z0p *= skd;
    sigma0p *= skd;
    n1p *= np;
    z1p *= skd;
    sigma1p *= skd;
    n2p *= np;
    z2p *= skd;
    sigma2p *= skd;
  }

  TF1 *fDenH = new TF1("fDenH",DensityGauss,zmin,zmax,9);
  fDenH->SetParameters(n0p,z0p,sigma0p,n1p,z1p,sigma1p,n2p,z2p,sigma2p);
  fDenH->SetNpx(1000);
  
  TF1 *fDenSum = new TF1("fDenSum",DensityGaussSum,zmin,zmax,18);

  Double_t pars[18] = {n0,z0,sigma0,n1,z1,sigma1,n2,z2,sigma2,
		      n0p,z0p,sigma0p,n1p,z1p,sigma1p,n2p,z2p,sigma2p};
  fDenSum->SetParameters(pars);
  fDenSum->SetNpx(1000);
  
  
  // Canvas setup
  // Create the canvas and the pads before the Frame loop
  // Resolution:
  Int_t sizex = 600;
  Int_t sizey = 300;
  char cName[32];
  sprintf(cName,"C");     
  TCanvas *C = (TCanvas*) gROOT->FindObject(cName);
  if(C==NULL) C = new TCanvas("C","",sizex,sizey);
  C->SetFillStyle(4000);
  C->cd();
  C->Clear();

  Float_t lMargin = 0.12;
  Float_t rMargin = 0.04;
  Float_t bMargin = 0.20;
  Float_t tMargin = 0.04;

  gPad->SetLeftMargin(lMargin);
  gPad->SetRightMargin(rMargin);  
  gPad->SetBottomMargin(bMargin);
  gPad->SetTopMargin(tMargin);

  Int_t fonttype = 43;
  Int_t fontsize = 18;
  Int_t tfontsize = 20;
  Int_t txsize = tfontsize;
  Int_t lxsize = fontsize;
  Int_t tysize = tfontsize;
  Int_t lysize = fontsize;
  Float_t txoffset = 1.1;
  Float_t lxoffset = 0.02;
  Float_t tyoffset = 0.7;
  Float_t lyoffset = 0.01;
  Float_t tylength = 0.01;
  Float_t txlength = 0.02;

  char name[16];
  sprintf(name,"hFrame");
  TH1F *hFrame = new TH1F(name,"",10,zmin,zmax);

  hFrame->GetXaxis()->SetLabelFont(fonttype);
  hFrame->GetXaxis()->SetLabelSize(lxsize);
  hFrame->GetXaxis()->SetLabelOffset(lxoffset);
  hFrame->GetXaxis()->SetTitleFont(fonttype);
  hFrame->GetXaxis()->SetTitleSize(txsize);
  hFrame->GetXaxis()->SetTitleOffset(txoffset);   
  hFrame->GetXaxis()->SetTickLength(txlength);      
  
  hFrame->GetYaxis()->SetLabelFont(fonttype);
  hFrame->GetYaxis()->SetLabelSize(lysize);
  hFrame->GetYaxis()->SetLabelOffset(lyoffset);
  hFrame->GetYaxis()->SetTitleFont(fonttype);
  hFrame->GetYaxis()->SetTitleSize(tysize);
  hFrame->GetYaxis()->SetTitleOffset(tyoffset);
  hFrame->GetYaxis()->SetTickLength(tylength);

  hFrame->GetXaxis()->CenterTitle();
  hFrame->GetYaxis()->CenterTitle();

  Float_t maxHe = fDenHe->GetMaximum();
  Float_t maxH = fDenH->GetMaximum();

  Float_t maxDen = maxHe ;
  if( maxH>maxHe) maxDen = maxH;
  
  hFrame->GetYaxis()->SetRangeUser(0.0, 1.2 * maxDen);  
  
  // hFrame->GetXaxis()->SetTitle("n_{He} [10^{15} e/cm^{3}]");
  hFrame->GetXaxis()->SetTitle("k_{p} z");
  hFrame->GetYaxis()->SetTitle("n/n_{0}");

  if(opt.Contains("units")) {
    hFrame->GetXaxis()->SetTitle("z [mm]");
    hFrame->GetYaxis()->SetTitle("n [10^{17} cm^{-3}]");
  }
  
  hFrame->Draw("AXIS");

  fDenSum->SetLineColor(kGray+2);
  fDenSum->SetLineWidth(1);
  //  fDenSum->Draw("same C");

  fDenHe->SetLineColor(myNaranja);
  fDenHe->SetLineWidth(2);
  // fDenHe->Draw("same C");

  fDenH->SetLineColor(myBlue);
  fDenH->SetLineWidth(2);
  fDenH->Draw("same C");

  gPad->Update();

  TBox *lFrame = new TBox(gPad->GetUxmin(), gPad->GetUymin(),
  			  gPad->GetUxmax(), gPad->GetUymax());
  lFrame->SetFillStyle(0);
  lFrame->SetLineColor(kGray+3);
  lFrame->SetLineWidth(2);
  lFrame->Draw();
  
  gPad->RedrawAxis("g");   

  // Print to file --------------------------------------
  
  C->cd();
  
  // Print to a file
  // Output file
  TString fOutName = Form("./%s/Plots/DenProf/DenProf-%s",sim.Data(),sim.Data());
  PGlobals::imgconv(C,fOutName,opt);
  // ---------------------------------------------------------

  
}
  
