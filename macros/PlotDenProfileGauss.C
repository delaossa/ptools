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
  Double_t f = 0;

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

Double_t DensityGaussTap(Double_t *x, Double_t *par) {
  Double_t f = 0;

  Double_t n0 = par[0];
  Double_t z0 = par[1];
  Double_t sigma0 = par[2];
  Double_t n1 = par[3];
  Double_t z1 = par[4];
  Double_t sigma1 = par[5];
  Double_t nt = par[6];
  Double_t zt = par[7];
  Double_t sigmat = par[8];
  
  if( x[0] < z0) f = n0 * TMath::Gaus(x[0],z0,sigma0);
  else if( x[0] < z1) f = n0;
  else if( x[0] >= z1 && x[0]<zt) {
    f = (n0-n1) * TMath::Gaus(x[0],z1,sigma1) + n1
      + (nt-n1) *  TMath::Gaus(x[0],zt,sigmat);
  } else if( x[0] >= zt) f = nt;
  
  return f;
}


Double_t DensityGaussSum(Double_t *x, Double_t *par) {
  Double_t f = DensityGauss(x,par) +  DensityGauss(x,&par[9]);

  
  return f;
}


void PlotDenProfileGauss(const TString &sim, const Float_t zmin0 = 0, const Float_t zmax0 = 100, const TString &options="pdf") { 
  
#ifdef __CINT__  
  gSystem->Load("libptools.so");
#endif

  PGlobals::Initialize();
  

  TString opt = options;

  // More makeup            
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetLineWidth(3);  
  
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
  
  Double_t n0 = 10.0;
  Double_t z0 = 60.0 - shiftz;
  Double_t sigma0 = 7.5;
  Double_t n1 = 1.0;
  Double_t z1 = z0;
  Double_t sigma1 = 2.5;
  Double_t n2 = 0;
  Double_t z2 = 10000 - shiftz;
  Double_t sigma2 = 0.0;
  Double_t nt = 1.5;
  Double_t zt = 100 - shiftz;
  Double_t sigmat = 5.0;

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
    nt *= np;
    zt *= skd;
    sigmat *= skd;
  }

  TF1 *fDenProf = new TF1("fDenProf",DensityGauss,zmin,zmax,9);
  fDenProf->SetParameters(n0,z0,sigma0,n1,z1,sigma1,n2,z2,sigma2);
  fDenProf->SetNpx(10000);

  TF1 *fDenProf2 = new TF1("fDenProf2",DensityGauss,zmin,zmax,9);
  fDenProf2->SetParameters(n0,z0,sigma0,n1,z1,2*sigma1,n2,z2,sigma2);
  fDenProf2->SetNpx(10000);

  TF1 *fDenProf3 = new TF1("fDenProf3",DensityGauss,zmin,zmax,9);
  fDenProf3->SetParameters(n0,z0,sigma0,n1,z1,3*sigma1,n2,z2,sigma2);
  fDenProf3->SetNpx(10000);

  TF1 *fDenProf4 = new TF1("fDenProf4",DensityGauss,zmin,zmax,9);
  fDenProf4->SetParameters(n0,z0,sigma0,n1,z1,4*sigma1,n2,z2,sigma2);
  fDenProf4->SetNpx(10000);

  TF1 *fDenProf5 = new TF1("fDenProf5",DensityGaussTap,zmin,zmax,9);
    fDenProf5->SetParameters(n0,z0,sigma0,n1,z1,sigma1,nt,zt,sigmat);
  fDenProf5->SetNpx(10000);

  const Int_t Nsim = 4;
  Int_t color[Nsim];
  PPalette *pal = new PPalette("pal");
  pal->SetPalette("oli");
  for(Int_t i=0;i<Nsim;i++) {
    
    Int_t index =  i * pal->GetNColors() / Nsim;
    color[i] = pal->GetColorIndex(index);

    // switch(i) {
    // case 0 :
    //   color[i] = kGray+2;
    //   break;
    // case 1 :
    //   color[i] = kRed-7;
    //   break;
    // case 2 :
    //   color[i] = kGray+3;
    //   break;
    // case 3 :
    //   color[i] = kRed;
    //   break;
    // }
  }
  
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
  Float_t tylength = 0.015;
  Float_t txlength = 0.03;

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

  Float_t maxDen = fDenProf->GetMaximum();
  hFrame->GetYaxis()->SetRangeUser(0.0, 1.2 * maxDen);  
  
  // hFrame->GetXaxis()->SetTitle("n_{He} [10^{15} e/cm^{3}]");
  hFrame->GetXaxis()->SetTitle("k_{p}^{0} z");
  hFrame->GetYaxis()->SetTitle("n/n_{0}");

  if(opt.Contains("units")) {
    hFrame->GetXaxis()->SetTitle("z [mm]");
    hFrame->GetYaxis()->SetTitle("n [10^{17} cm^{-3}]");
  }
  
  hFrame->Draw("AXIS");

  // Lines and guides
  Float_t xMin = hFrame->GetXaxis()->GetXmin();
  Float_t xMax = hFrame->GetXaxis()->GetXmax();
  Float_t yMin = hFrame->GetYaxis()->GetXmin();
  Float_t yMax = 1.2 * maxDen;

  TLine *linezero = new TLine(0.0,yMin,0.0,yMax);
  linezero->SetLineColor(kGray);
  linezero->SetLineWidth(2);
  linezero->SetLineStyle(3);
  linezero->Draw();

  TLine *linemax = new TLine(xMin,maxDen,xMax,maxDen);
  linemax->SetLineColor(kGray);
  linemax->SetLineWidth(2);
  linemax->SetLineStyle(3);
  linemax->Draw();

  TLine *lineplateau = new TLine(xMin,n1,xMax,n1);
  lineplateau->SetLineColor(kGray);
  lineplateau->SetLineWidth(2);
  lineplateau->SetLineStyle(3);
  lineplateau->Draw();

  // Functions


  //  fDenProf4->SetLineColor(kGray);
  fDenProf4->SetLineColor(color[3]);
  fDenProf4->SetLineWidth(3);
  fDenProf4->Draw("same C");

  //  fDenProf3->SetLineColor(kGray+1);
  fDenProf3->SetLineColor(color[2]);
  fDenProf3->SetLineWidth(3);
  fDenProf3->Draw("same C");

  //  fDenProf2->SetLineColor(kGray+2);
  fDenProf2->SetLineColor(color[1]);
  fDenProf2->SetLineWidth(3);
  fDenProf2->Draw("same C");

  //  fDenProf->SetLineColor(kBlack);
  fDenProf->SetLineColor(color[0]);
  fDenProf->SetLineWidth(3);
  fDenProf->Draw("same C");

  //  fDenProf5->SetLineColor(kRed);
  fDenProf5->SetLineColor(color[4]);
  fDenProf5->SetLineWidth(3);
  // fDenProf5->SetLineStyle(3);
  // fDenProf5->Draw("same C");
  
  hFrame->Draw("AXIS same");

  gPad->Update();

  TBox *lFrame = new TBox(gPad->GetUxmin(), gPad->GetUymin(),
  			  gPad->GetUxmax(), gPad->GetUymax());
  lFrame->SetFillStyle(0);
  lFrame->SetLineColor(kBlack);
  lFrame->SetLineWidth(2);
  lFrame->Draw();

  TLegend *Leg = new TLegend(0.65,0.35,0.85,0.80);
  PGlobals::SetPaveStyle(Leg);
  Leg->SetTextAlign(12);
  Leg->SetTextColor(kGray+3);
  Leg->SetTextFont(43);
  Leg->SetTextSize(14);
  Leg->SetLineColor(kGray+1);
  Leg->SetBorderSize(0);
  Leg->SetLineWidth(1);
  Leg->SetFillColor(0);
  Leg->SetFillStyle(1001);
  Leg->SetFillStyle(0); // Hollow

  Leg->AddEntry(fDenProf,"(a): k_{p}^{0}#sigma_{l} = 2.5","L");	
  Leg->AddEntry(fDenProf2,"(b): k_{p}^{0}#sigma_{l} = 5.0","L");	
  Leg->AddEntry(fDenProf3,"(c): k_{p}^{0}#sigma_{l} = 7.5","L");	
  Leg->AddEntry(fDenProf4,"(d): k_{p}^{0}#sigma_{l} = 10","L");	
  //  Leg->AddEntry(fDenProf5,"(e): k_{p}^{0}#sigma_{l} = 2.5 [tap]","L");	
  Leg->Draw();
  
  gPad->RedrawAxis("g");   

  // Print to file --------------------------------------
  
  C->cd();
  
  // Print to a file
  // Output file
  TString fOutName = Form("./DenProf/DenProf-%s",sim.Data());
  PGlobals::imgconv(C,fOutName,opt);
  // ---------------------------------------------------------

  
}
  
