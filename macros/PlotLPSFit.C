#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>

#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TLegend.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TProfile.h>
#include <TPaletteAxis.h>

#include "PData.hh"
#include "PDataHiP.hh"
#include "PGlobals.hh"
#include "PPalette.hh"

void PlotLPSFit(const TString &sim, Int_t timestep, Int_t index = 1) {
  
  gStyle->SetPadRightMargin(0.20);

  TCanvas *C = new TCanvas("C","Longitudinal phasespace",800,600);
  C->SetFillStyle(4000);
  C->cd();
  C->Clear();
  
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

  // Some plasma constants
  Double_t n0 = pData->GetPlasmaDensity();
  Double_t kp = pData->GetPlasmaK() * TMath::Sqrt(0.1);
  Double_t skd= (1/kp)/PUnits::mm;
  
  cout << Form(" - Skindepth = %.2f mm",skd) << endl;
  
  // Open snapshot file and get the histograms
  TString filename = Form("./%s/Plots/Bunch/%s/Bunch-%s-%s_%i.root",
			  sim.Data(),pData->GetRawSpeciesName(index).c_str(),
			  pData->GetRawSpeciesName(index).c_str(),
			  sim.Data(),
			  timestep);
  TFile *file = new TFile(filename,"OPEN");
  

  char hName[16];
  sprintf(hName,"hP1X1");
  TH2F *hP1X1 = (TH2F*) file->Get(hName);
  sprintf(hName,"hP1");
  TH1F *hP1 = (TH1F*) file->Get(hName);
  hP1->ResetStats();
  
  // Set palette:
  PPalette * pPalette = (PPalette*) gROOT->FindObject("electron0");
  if(!pPalette)
    pPalette = new PPalette("electron0");
  pPalette->SetPalette("electron0");

  hP1X1->GetXaxis()->CenterTitle();
  hP1X1->GetYaxis()->CenterTitle();
  hP1X1->Draw("colz");
  
  gPad->Update();
  TPaletteAxis *palette = (TPaletteAxis*)hP1X1->GetListOfFunctions()->FindObject("palette");
  if(palette) {
    Double_t y1 = gPad->GetBottomMargin();
    Double_t y2 = 1 - gPad->GetTopMargin();
    Double_t x1 = 1 - gPad->GetRightMargin();
    
    Double_t x1b = x1 + 0.01;
    Double_t x2b = x1 + 0.035;
    Double_t y1b = y1 + 0.03;
    Double_t y2b = y2 - 0.03;
    palette->SetX1NDC(x1b);
    palette->SetY1NDC(y1b);
    palette->SetX2NDC(x2b);
    palette->SetY2NDC(y2b);
    palette->SetBorderSize(2);
    palette->SetLineColor(1);
    
    TPave *pFrame = new TPave(x1b,y1b,x2b,y2b,1,"NDCL");
    pFrame->SetFillStyle(0);
    pFrame->SetLineColor(kBlack);
    pFrame->SetShadowColor(0);
    pFrame->Draw();
  }
  
  TString pname = hName;
  pname += "_pfx";

  TProfile *hP1X1prof = hP1X1->ProfileX("_pfx",1,-1,"s");
  //  hP1X1prof->Draw("same");

  Int_t Npoints = hP1X1prof->GetNbinsX();
  Int_t Np = 0;
  Double_t *x = new Double_t[Npoints];
  Double_t *y = new Double_t[Npoints];
  Double_t *yup = new Double_t[Npoints];
  Double_t *ydo = new Double_t[Npoints];
  Double_t *ex = new Double_t[Npoints];
  Double_t *ey = new Double_t[Npoints];

  Double_t xmin = -9;
  Double_t xmax = -0.5;
  
  TGraphErrors *gP1X1 = (TGraphErrors*) gROOT->FindObject("gP1X1");
  if(gP1X1) delete gP1X1;
  
  for(Int_t j=0;j<Npoints;j++) {

    if(hP1X1prof->GetBinCenter(j)<xmin || hP1X1prof->GetBinCenter(j)>xmax) {
      continue;
    }
    
    x[Np] = hP1X1prof->GetBinCenter(j);
    y[Np] = hP1X1prof->GetBinContent(j);
    ex[Np] = 0;
    ey[Np] = hP1X1prof->GetBinError(j);
    yup[Np] = y[Np]+ey[Np];
    ydo[Np] = y[Np]-ey[Np];
    
    Np++;
  }
  
  gP1X1 = new TGraphErrors(Np,x,y,ex,ey);
  gP1X1->SetName("gP1X1");

  TGraph *gP1X1avg = new TGraph(Np,x,y); 
  TGraph *gP1X1up = new TGraph(Np,x,yup); 
  TGraph *gP1X1do = new TGraph(Np,x,ydo); 

  gP1X1up->GetXaxis()->SetRangeUser(xmin,xmax);
  gP1X1do->GetXaxis()->SetRangeUser(xmin,xmax);

  //gP1X1->Draw("3");
  gP1X1avg->SetMarkerStyle(20);
  gP1X1avg->SetMarkerSize(0.4);
  gP1X1avg->SetMarkerColor(kWhite);
  gP1X1avg->SetLineStyle(1);
  gP1X1avg->SetLineWidth(1);
  gP1X1avg->SetLineColor(kWhite);  
  gP1X1avg->Draw("C");
  gP1X1up->SetMarkerStyle(20);
  gP1X1up->SetMarkerSize(0.2);
  gP1X1up->SetMarkerColor(kGray+2);
  gP1X1up->SetLineStyle(3);
  gP1X1up->SetLineWidth(2);
  gP1X1up->SetLineColor(kGray);
  gP1X1up->Draw("C");
  gP1X1do->SetMarkerStyle(20);
  gP1X1do->SetMarkerSize(0.2);
  gP1X1do->SetMarkerColor(kGray+2);
  gP1X1do->SetLineStyle(3);
  gP1X1do->SetLineWidth(2);
  gP1X1do->SetLineColor(kGray);
  gP1X1do->Draw("C");

  // gP1X1up->Draw("L");
  // gP1X1do->Draw("L");
  // gPad->Update();
  
  char fitName[64];
  sprintf(fitName,"fitFunction_%i",timestep);
  TF1 *fitFunction = (TF1*) gROOT->FindObject(fitName);
  if(!fitFunction) 
    fitFunction = new TF1(fitName,"[0]+[4]*x+[1]*sin([2]*(x-[3]))",xmin,xmax);

  Double_t p1max = -99;
  for(Int_t j=0;j<Np;j++) {
    if(y[j]>p1max) p1max = y[j];
  }
  Double_t p1mean = hP1->GetMean();
  TLine lp1mean(hP1X1->GetXaxis()->GetXmin(),p1mean,hP1X1->GetXaxis()->GetXmax(),p1mean);
  lp1mean.SetLineColor(kGray+2);
  lp1mean.SetLineStyle(2);
  lp1mean.Draw();

  cout << Form(" p1mean = %.2f",p1mean) << endl;

  
  fitFunction->SetParameters(p1mean,1.0,1./skd,0,0);
  fitFunction->SetParLimits(1,0.1,1.5);
  
  Int_t res = gP1X1->Fit(fitName,"R0");

  // Retrieve parameters of the fit
  Double_t chi2 = fitFunction->GetChisquare();
  Double_t NDF = fitFunction->GetNDF();
  const Double_t *par = fitFunction->GetParameters();
  const Double_t *parErr = fitFunction->GetParErrors();
  
  // cout << " x1 = " << x1 << "  x2 = " << x2 << "  k = " << k << endl;
  if(res==0) {
    cout << "Fit results: " << endl;
    cout << Form("chi2/NDF = %5.2f/%5.2f : A = %5.2f  B = %5.2f  C = %5.2f  D  = %5.2f",chi2,NDF,par[0],par[1],par[2],par[3]) << endl;
  }

  
  fitFunction->SetLineStyle(2);
  fitFunction->SetLineWidth(2);
  fitFunction->SetLineColor(kRed);
  fitFunction->Draw("same");

  Double_t wavelength = TMath::TwoPi()/par[2];
  Double_t waveerr    = parErr[2] * TMath::TwoPi()/(par[2]*par[2]);

  TLegend *leg = new TLegend(0.5,0.75,0.78,0.92);
  PGlobals::SetPaveStyle(leg);
  leg->SetTextAlign(12);
  leg->SetTextColor(kGray+3);
  leg->SetTextFont(43);
  leg->SetTextSize(16);
  leg->SetLineColor(1);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(1001);
  leg->SetFillStyle(0); // Hollow
  leg->SetMargin(0.1);

  // leg->AddEntry(gP1X1avg,"profile", "l");
  leg->AddEntry(fitFunction,"fit: p0 + p4 #zeta + p1 sin(p2 (#zeta-p3))", "l");
  leg->AddEntry((TObject*)0,Form("Amplitude  = %.2f #pm %.2f MeV/c",par[1],parErr[1]), "");
  leg->AddEntry((TObject*)0,Form("Wavelength = %.2f #pm %.2f mm",wavelength,waveerr), "");
  leg->AddEntry((TObject*)0,Form("#chi^{2}/Ndf = %.2f",chi2/NDF), "");
  leg->Draw();
  
  TString fOutName = Form("./%s/Plots/Bunch/%s/Bunch-%s-%s-P1X1fit_%i",
			  sim.Data(),pData->GetRawSpeciesName(index).c_str(),
			  pData->GetRawSpeciesName(index).c_str(),
			  sim.Data(),
			  timestep);
  PGlobals::imgconv(C,fOutName,"pdf");
  
}
