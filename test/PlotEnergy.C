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
#include <TExec.h>

#include "PData.hh"
#include "PlasmaGlob.hh"
#include "PPalette.hh"

#ifndef __CINT__
#include "PFunctions.hh"
#endif

void PlotEnergy( const TString &sim, Int_t time, Float_t Emin, Float_t Emax, Int_t zoom=2, const TString &options="") {
  
#ifdef __CINT__  
  gSystem->Load("libplasma.so");
#endif

  PlasmaGlob::Initialize();

  // Palettes!
  gROOT->Macro("PlasmaPalettes.C");


  TString opt = options;

  // More makeup
  gStyle->SetPadTopMargin(0.07); 
  gStyle->SetPadLeftMargin(0.10);  // Margin left axis 
  gStyle->SetPadRightMargin(0.14); // Margin for palettes in 2D histos
  gStyle->SetTitleSize(0.06, "x");
  gStyle->SetTitleSize(0.06, "y");
  gStyle->SetTitleSize(0.06, "z");
  gStyle->SetTitleOffset(1.2, "x");
  gStyle->SetTitleOffset(0.7, "y");
  gStyle->SetTitleOffset(0.6, "z");
  if(opt.Contains("grid")) {
    gStyle->SetPadGridX(1);
    gStyle->SetPadGridY(1);
  }
  
  // Load PData
  PData *pData = PData::Get(sim.Data());
  pData->LoadFileNames(time);
  if(!pData->IsInit()) return;
  
  Bool_t CYL = kFALSE;
  if(sim.Contains("cyl"))  { CYL = kTRUE; opt += "cyl"; } 
  
  Bool_t ThreeD = kFALSE;
  if(sim.Contains("3D")) ThreeD = kTRUE; 
  
  // Extract histogram with the plasma density
  TH2F *hDen2D = NULL;
  if(pData->GetChargeFileName(1)) {

    if(!ThreeD)
      hDen2D = pData->GetCharge(1);
    else
      hDen2D = pData->GetCharge2DSliceZY(1);

  } else {
    cout << " Plasma density file : " << pData->GetChargeFileName(0) 
	 << " not present!! -> EXIT. " << endl;
  }
  
  // Some plasma constants
  Double_t n0 = pData->GetPlasmaDensity();
  Double_t kp = pData->GetPlasmaK();
  Double_t skindepth = 1/kp;
  Double_t E0 = pData->GetPlasmaE0();

  // Some beam properties:
  Double_t Ebeam = pData->GetBeamEnergy();
  Double_t gamma = Ebeam / PConst::ElectronMassE;
  Double_t vbeam = TMath::Sqrt(1 - 1/(gamma*gamma));

  // Time in OU
  Double_t Time = pData->GetRealTime();
  
  // z start of the plasma in normalized units.
  Double_t zStartPlasma = pData->GetPlasmaStart()*kp;
  
  // z start of the beam in normalized units.
  Double_t zStartBeam = pData->GetBeamStart()*kp;
   
  // For center and comov options
  Float_t shiftz = 0;
  if(opt.Contains("center") || opt.Contains("comov")) {
    Time -= zStartPlasma - zStartBeam;
    if(opt.Contains("center")) {
      shiftz += zStartPlasma;
    } if(opt.Contains("comov")) {
      shiftz += pData->GetBeamVelocity() * Time;
    }   
  }
  
  // Binning for 2D hitograms:
  // We get this values from the 2D density histogram.
  Int_t x1Nbin, x2Nbin; x1Nbin = x2Nbin = 0;
  Double_t x1Range, x1Mid, x1Min, x1Max;
  x1Range = x1Mid = x1Min = x1Max = 0.0;
  Double_t x2Range, x2Mid, x2Min, x2Max;
  x2Range = x2Mid = x2Min = x2Max = 0.0;

  if(hDen2D) {
    x1Nbin    = hDen2D->GetNbinsX();
    x1Range   = (hDen2D->GetXaxis()->GetXmax() - hDen2D->GetXaxis()->GetXmin())/zoom;
    x1Mid     = (hDen2D->GetXaxis()->GetXmax() + hDen2D->GetXaxis()->GetXmin())/2. - shiftz;
    x1Min     = hDen2D->GetXaxis()->GetXmin() - shiftz;
    x1Max     = hDen2D->GetXaxis()->GetXmax() - shiftz;
    
    x2Nbin    = hDen2D->GetNbinsY();      
    x2Range   = (hDen2D->GetYaxis()->GetXmax() - hDen2D->GetYaxis()->GetXmin())/zoom;
    x2Mid     = (hDen2D->GetYaxis()->GetXmax() + hDen2D->GetYaxis()->GetXmin())/2.;
    x2Min     = x2Mid - x2Range/2;
    x2Max     = x2Mid + x2Range/2;
    if(CYL) {
      x2Min     = 0.0;
      x2Max     = x2Range;
      x2Mid     = (x2Min+x2Max)/2.0;
    }
  }
  
  // Binning for Energy histograms:
  Int_t eNbin = 200;
  Double_t eMin = Emin;
  Double_t eMax = Emax;

  // Define pointers to Histograms
  Int_t Nspecies = pData->NSpecies();
  TH1F **hEnergy = new TH1F*[Nspecies];    // Energy spectrum
  TH2F **hEvsX1 = new TH2F*[Nspecies];     // Energy vs x1
  TH2F **hEvsX2 = new TH2F*[Nspecies];     // Energy vs x2
  for(Int_t i=0;i<Nspecies;i++) {
    hEnergy[i] = NULL;
    hEvsX1[i] = hEvsX2[i] = NULL;

    if(!pData->GetRawFileName(i))
      continue;

    char hName[24];
    sprintf(hName,"hEnergy_%i",i);
    hEnergy[i] = (TH1F*) gROOT->FindObject(hName);
    if(hEnergy[i]) { delete hEnergy[i]; hEnergy[i] = NULL; }
    hEnergy[i] = new TH1F(hName,"",eNbin,eMin,eMax);
    
    pData->GetH1Raw(pData->GetRawFileName(i)->c_str(),"energy",hEnergy[i],opt);
   
    hEnergy[i]->GetXaxis()->CenterTitle();
    hEnergy[i]->GetYaxis()->CenterTitle();
    hEnergy[i]->GetXaxis()->SetTitle("Energy [MeV]");
    hEnergy[i]->GetYaxis()->SetTitle("Charge [a.u.]");
    PlasmaGlob::SetH1Style(hEnergy[i],i);
    
    sprintf(hName,"hEvsX1_%i",i);
    hEvsX1[i] = (TH2F*) gROOT->FindObject(hName);
    if(hEvsX1[i]) { delete hEvsX1[i]; hEvsX1[i] = NULL; }
    hEvsX1[i] = new TH2F(hName,"",x1Nbin,x1Min,x1Max,eNbin,eMin,eMax);
    pData->GetH2Raw(pData->GetRawFileName(i)->c_str(),"x1","energy",hEvsX1[i],opt);
        
    hEvsX1[i]->GetXaxis()->CenterTitle();
    hEvsX1[i]->GetYaxis()->CenterTitle();
    hEvsX1[i]->GetZaxis()->CenterTitle();
    hEvsX1[i]->GetXaxis()->SetTitle("z [c/#omega_{p}]");
    hEvsX1[i]->GetYaxis()->SetTitle("Energy [MeV]");
    hEvsX1[i]->GetZaxis()->SetTitle("Charge [a.u.]");
    PlasmaGlob::SetH1Style(hEvsX1[i],i);
  
    sprintf(hName,"hEvsX2_%i",i);
    hEvsX2[i] = (TH2F*) gROOT->FindObject(hName);
    if(hEvsX2[i]) { delete hEvsX2[i]; hEvsX2[i] = NULL; }
    hEvsX2[i] = new TH2F(hName,"",x2Nbin,x2Min,x2Max,eNbin,eMin,eMax);
    pData->GetH2Raw(pData->GetRawFileName(i)->c_str(),"x2","energy",hEvsX2[i],opt);
        
    hEvsX2[i]->GetXaxis()->CenterTitle();
    hEvsX2[i]->GetYaxis()->CenterTitle();
    hEvsX2[i]->GetZaxis()->CenterTitle();
    if(CYL)
      hEvsX2[i]->GetXaxis()->SetTitle("r [c/#omega_{p}]");
    else
      hEvsX2[i]->GetXaxis()->SetTitle("y [c/#omega_{p}]");
    hEvsX2[i]->GetYaxis()->SetTitle("Energy [MeV]");
    hEvsX2[i]->GetZaxis()->SetTitle("Charge [a.u.]");
    PlasmaGlob::SetH1Style(hEvsX2[i],i);
    
  }
  
  // Vertical graphs: Displayed at a side of 2D histograms
  TGraph *gEnergyX1 = NULL;
  TGraph *gEnergyX2 = NULL;
  if(hEnergy[1]) {
    Int_t Nbin   = hEnergy[1]->GetNbinsX();
    Float_t *y   = new Float_t[Nbin];
    Float_t *x   = new Float_t[Nbin];

    // This is for the right side:
    // Float_t xMax = x1Min + (x1Max-x1Min) * 0.9;
    // Float_t xMin = x1Max;
    // And this for left:
    Float_t xMin = x1Min;
    Float_t xMax = x1Min + (x1Max-x1Min) * 0.2;
    Float_t Emax = hEnergy[1]->GetMaximum();
    for(Int_t i=0; i<Nbin; i++) {
      y[i] = hEnergy[1]->GetBinCenter(i+1);
      x[i] = ((xMax-xMin)/Emax)*hEnergy[1]->GetBinContent(i+1) + xMin;
    }
    gEnergyX1 = new TGraph(Nbin,x,y);
    gEnergyX1->SetLineColor(PlasmaGlob::elecLine);
    gEnergyX1->SetLineWidth(2);
    gEnergyX1->SetFillStyle(1001);
    gEnergyX1->SetFillColor(PlasmaGlob::elecFill);

    if(CYL) {
      // Right side:
      xMax = x2Min + (x2Max-x2Min) * 0.8;
      xMin = x2Max;
    } else {
      // Left side:
      xMin = x2Min;
      xMax = x2Min +  (x2Max-x2Min) * 0.2;    
    }
    
    for(Int_t i=0; i<Nbin; i++) {
      y[i] = hEnergy[1]->GetBinCenter(i+1);
      x[i] = ((xMax-xMin)/Emax)*hEnergy[1]->GetBinContent(i+1) + xMin;
    }
    gEnergyX2 = new TGraph(Nbin,x,y);
    gEnergyX2->SetLineColor(PlasmaGlob::elecLine);
    gEnergyX2->SetLineWidth(2);
    gEnergyX2->SetFillStyle(1001);
    gEnergyX2->SetFillColor(PlasmaGlob::elecFill);
  }
  
  // Plotting
  // -----------------------------------------------


  // Canvas setup
  TCanvas *C = new TCanvas("C","2D Energy distribution",750,1000);
  C->Divide(1,3);

  // Palettes setup
  TExec *exElec = new TExec("exElec","electronPalette->cd();");
  
  // Text objects
  TPaveText *text1 = new TPaveText(0.6,0.93,0.86,1.0,"NDC");
  PlasmaGlob::SetPaveTextStyle(text1,32);
  text1->AddText("Energy distribution");
  
  TPaveText *textTime = new TPaveText(0.7,0.82,0.85,0.88,"NDC");
  PlasmaGlob::SetPaveTextStyle(textTime); 
  char ctext[128];
  if(opt.Contains("units") && pData->GetPlasmaDensity()) 
    sprintf(ctext,"Z = %5.1f mm", 1e3 * pData->GetPlasmaSkinDepth() * Time);
  else
    sprintf(ctext,"T = %5.1f #omega_{p}^{-1}",Time);
  textTime->AddText(ctext);
  
  // Actual Plotting!
  // ------------------------------------------------------------

  // Output file
  TString fOutName = Form("./%s/Plots/Energy/Energy",pData->GetPath().c_str());
  fOutName += Form("-%s_%i",pData->GetName(),time);

  C->cd(1); // <--- Top Plot

  hEnergy[1]->GetYaxis()->SetNdivisions(505);
 
  if(hEnergy[1])
    hEnergy[1]->Draw("L");
  if(hEnergy[2] && Nspecies>=3)
    hEnergy[2]->Draw("same L");
  
  text1->Draw();
  textTime->Draw();

  gPad->RedrawAxis(); 
  gPad->RedrawAxis("G"); 


  TLine *lineEne0 = new TLine(Ebeam/PUnits::MeV,0.0,
			      Ebeam/PUnits::MeV,hEnergy[1]->GetMaximum());
  lineEne0->SetLineColor(kGray+1);
  lineEne0->SetLineStyle(2);
  lineEne0->Draw();

  C->cd(2); // <--- Mid Plot

  TH2F *hFrame2 = (TH2F*) hEvsX1[1]->Clone("hFrame2");
  hFrame2->Reset();
  hFrame2->GetYaxis()->SetNdivisions(505);
  hFrame2->Draw("col");

  exElec->Draw();
  if(hEvsX1[1]) 
    hEvsX1[1]->Draw("col same z");
  if(hEvsX1[2] && Nspecies>=3)
    hEvsX1[2]->Draw("col same");  
  
  gEnergyX1->Draw("F");
  gEnergyX1->Draw("L");
  
  textTime->Draw();
  gPad->RedrawAxis(); 
  gPad->RedrawAxis("G"); 

  TLine *lineEne0X = new TLine(hFrame2->GetXaxis()->GetXmin(),Ebeam/PUnits::MeV,
			       hFrame2->GetXaxis()->GetXmax(),Ebeam/PUnits::MeV);
  lineEne0X->SetLineColor(kGray+1);
  lineEne0X->SetLineStyle(2);
  lineEne0X->Draw();

  C->cd(3); // <--- Bottom Plot

  TH2F *hFrame3 = (TH2F*) hEvsX2[1]->Clone("hFrame3");
  hFrame3->Reset();
  hFrame3->GetYaxis()->SetNdivisions(505);
  hFrame3->Draw("col");

  exElec->Draw();
  if(hEvsX2[1]) 
    hEvsX2[1]->Draw("col same z");
  if(hEvsX2[2] && Nspecies>=3)
    hEvsX2[2]->Draw("col same");  
  
  gEnergyX2->Draw("F");
  gEnergyX2->Draw("L");

  textTime->Draw();
  gPad->RedrawAxis(); 
  gPad->RedrawAxis("G"); 

  TLine *lineEne0Y = new TLine(hFrame3->GetXaxis()->GetXmin(),Ebeam/PUnits::MeV,
			       hFrame3->GetXaxis()->GetXmax(),Ebeam/PUnits::MeV);
  lineEne0Y->SetLineColor(kGray+1);
  lineEne0Y->SetLineStyle(2);
  lineEne0Y->Draw();

  C->Update();

  C->cd();
  
  // Print to a file
  PlasmaGlob::imgconv(C,fOutName,opt);

  // ---------------------------------------------------------
}
