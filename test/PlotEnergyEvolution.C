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
#include <TMarker.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TMath.h>
#include <TPaveText.h>
#include <TExec.h>
#include <TGaxis.h>
#include <TF1.h>

#include "PData.hh"
#include "PlasmaGlob.hh"


// CONSTANTS
// Speed of light
const Float_t c = 299792458;              // m/s
// electron charge
const Float_t e = 1.602176565e-19;        // (C)oulombs
// Vaccum permitivity
const Float_t eps0 = 8.854187817620e-12;  // F/m
// Electron mass
const Float_t me = 9.10938291e-31;        // Kg  

void PlotEnergyEvolution(const TString &opt="png") { 
  
  PlasmaGlob::Initialize();
    
  // More makeup
  // gStyle->SetPadTopMargin(0.02);   // Margin Top
  // gStyle->SetPadLeftMargin(0.10);  // Margin left axis 
  gStyle->SetPadRightMargin(0.05); // Margin for palettes in 2D histos
  gStyle->SetTitleSize(0.06, "x");
  gStyle->SetTitleSize(0.06, "y");
  // gStyle->SetTitleOffset(1., "x");
  // gStyle->SetTitleOffset(0.6, "y");
  gStyle->SetPadGridY(0);
  gStyle->SetPadGridX(1);
  
  // Get Files from different PITZ simulations
  const Int_t Nsim = 8;
  TString simNames[Nsim] = {"pitz_FT01","pitz_FT02","pitz_FT04","pitz_FT06","pitz_FT08","pitz_FT","pitz_FT2","pitz_FT4"};//,"pitz_FT6","pitz_FT8"};
  Int_t simColors[Nsim] = {kGray,kGray+1,kGray+2,kOrange-2,kOrange+7,kRed,kAzure-4,kAzure-5};//,kAzure-9,kGreen+1};
  Float_t den[Nsim] = {0.1e15,0.2e15,0.4e15,0.6e15,0.8e15,1.0e15,2.0e15,4.0e15};//,6.0e15,8.0e15};
  TString filename[Nsim];
  TFile  *ifile[Nsim];
  TGraphErrors *gEmean[Nsim];
  TGraphErrors *gEamp[Nsim];
  TGraphErrors *gChi2[Nsim];
  TGraphErrors *gLambda[Nsim];
  TGraphErrors *gPhase[Nsim];
  TGraphErrors *gAsym[Nsim];
  TF1 *fFit[Nsim];

  Float_t maxEmean = -999.;
  Float_t minEmean = 999.;
  Float_t maxEamp = -999.;
  Float_t minEamp = 999.;
  Float_t maxLambda = -999.;
  Float_t minLambda = 999.;
  Float_t maxPhase = -999.;
  Float_t minPhase = 999.;
  Float_t maxChi2 = -999.;
  Float_t minChi2 = 999.;  
  Float_t maxAsym = -999.;
  Float_t minAsym = 999.;
  Float_t maxFitVal = -999.;
  Float_t minFitVal = 999.;

  for(Int_t i=0;i<Nsim;i++) {
    filename[i] = Form("./%s/Plots/EnergyModulation/EnergyModulation-%s.root",simNames[i].Data(),simNames[i].Data());
    ifile[i] = (TFile*) gROOT->GetListOfFiles()->FindObject(filename[i].Data());
    if (!ifile[i]) ifile[i] = new TFile(filename[i],"READ");
    
    // Get the needed Graphs for every simulation on the list 
    gEmean[i] = (TGraphErrors*) ifile[i]->Get("gEmean");
    gEamp[i] = (TGraphErrors*) ifile[i]->Get("gEamp");
    gChi2[i] = (TGraphErrors*) ifile[i]->Get("gChi2");
    gLambda[i] = (TGraphErrors*) ifile[i]->Get("gLambda");
    gPhase[i] = (TGraphErrors*) ifile[i]->Get("gPhase");
    gAsym[i] = (TGraphErrors*) ifile[i]->Get("gAsym");

    char fname[64];
    sprintf(fname,"fitFunction_%i",800);
    fFit[i] = (TF1*) ifile[i]->Get(fname);
    
    // Calculate extremes for the set of Graphs
    Double_t *yEmean = gEmean[i]->GetY();
    Double_t *yEamp = gEamp[i]->GetY();
    Double_t *yChi2 = gChi2[i]->GetY();
    Double_t *yLambda = gLambda[i]->GetY();
    Double_t *yPhase = gPhase[i]->GetY();
    Double_t *yAsym = gAsym[i]->GetY();
    
    // Change lambda to normalized coordinates:
    Double_t *x    = gLambda[i]->GetX();
    Double_t *yErr = gLambda[i]->GetEY();
    Double_t skinDepth = 1e3 * c / TMath::Sqrt((1e6 * den[i] * e * e)/(eps0 * me));
    // Change the Graph to plasma units:
    for(Int_t j=0;j<gLambda[i]->GetN();j++) {
      Float_t value = yLambda[j]/(skinDepth*TMath::TwoPi());
      Float_t error = yErr[j]/(skinDepth*TMath::TwoPi());
      gLambda[i]->SetPoint(j,x[j],value);
      gLambda[i]->SetPointError(j,0.,error);
    }
    yLambda = gLambda[i]->GetY();
    

    for(Int_t j=0;j<gEmean[i]->GetN();j++) {
      if(yEmean[j]>maxEmean)
	maxEmean = yEmean[j];
      if(yEmean[j]<minEmean)
	minEmean = yEmean[j];
      
      if(yEamp[j]>maxEamp)
	maxEamp = yEamp[j];
      if(yEamp[j]<minEamp)
	minEamp = yEamp[j];

      if(yLambda[j]>maxLambda)
	maxLambda = yLambda[j];
      if(yLambda[j]<minLambda)
	minLambda = yLambda[j];

      if(yPhase[j]>maxPhase)
       	maxPhase = yPhase[j];
      if(yPhase[j]<minPhase)
       	minPhase = yPhase[j];
      
      if(yChi2[j]>maxChi2)
	maxChi2 = yChi2[j];
      if(yChi2[j]<minChi2)
	minChi2 = yChi2[j];

      if(yAsym[j]>maxAsym)
       	maxAsym = yAsym[j];
      if(yAsym[j]<minAsym)
       	minAsym = yAsym[j];
    }

    if(maxChi2>50.) maxChi2=50.;
    
    if(fFit[i]->GetMaximum()>maxFitVal)
      maxFitVal = fFit[i]->GetMaximum();
    if(fFit[i]->GetMinimum()>minFitVal)
      minFitVal = fFit[i]->GetMinimum();
    
    // Graph's attributes
    gEmean[i]->SetLineColor(simColors[i]);
    gEmean[i]->SetLineWidth(2);
    gEmean[i]->SetFillColor(simColors[i]);
    gEmean[i]->SetMarkerColor(simColors[i]);
    gEmean[i]->SetMarkerSize(0.2);
    gEmean[i]->GetXaxis()->CenterTitle();
    gEmean[i]->GetYaxis()->CenterTitle();
    
    gEamp[i]->SetLineColor(simColors[i]);
    gEamp[i]->SetLineWidth(2);
    gEamp[i]->SetFillColor(simColors[i]);
    gEamp[i]->SetMarkerColor(simColors[i]);
    gEamp[i]->SetMarkerSize(0.2);
    gEamp[i]->GetXaxis()->CenterTitle();
    gEamp[i]->GetYaxis()->CenterTitle();

    gChi2[i]->SetLineColor(simColors[i]);
    gChi2[i]->SetMarkerColor(simColors[i]);
    gChi2[i]->SetLineWidth(2);
    gChi2[i]->SetMarkerSize(0.2);
    gChi2[i]->GetXaxis()->CenterTitle();
    gChi2[i]->GetYaxis()->CenterTitle();

    gLambda[i]->SetLineColor(simColors[i]);
    gLambda[i]->SetLineWidth(2);
    gLambda[i]->SetFillColor(simColors[i]);   
    gLambda[i]->SetMarkerColor(simColors[i]);
    gLambda[i]->SetMarkerSize(0.2);
    gLambda[i]->GetXaxis()->CenterTitle();
    gLambda[i]->GetYaxis()->CenterTitle();

    gPhase[i]->SetLineColor(simColors[i]);
    gPhase[i]->SetLineWidth(2);
    gPhase[i]->SetFillColor(simColors[i]);
    gPhase[i]->SetMarkerColor(simColors[i]);
    gPhase[i]->SetMarkerSize(0.2);
    gPhase[i]->GetXaxis()->CenterTitle();
    gPhase[i]->GetYaxis()->CenterTitle();

    gAsym[i]->SetLineColor(simColors[i]);
    gAsym[i]->SetLineWidth(2);
    gAsym[i]->SetFillColor(simColors[i]);
    gAsym[i]->SetMarkerColor(simColors[i]);
    gAsym[i]->SetMarkerSize(0.2);
    gAsym[i]->GetXaxis()->CenterTitle();
    gAsym[i]->GetYaxis()->CenterTitle();

    fFit[i]->SetLineColor(simColors[i]);
    fFit[i]->SetLineWidth(2);
    fFit[i]->SetFillColor(simColors[i]);
    fFit[i]->SetMarkerColor(simColors[i]);
    fFit[i]->SetMarkerSize(0.2);
    fFit[i]->GetXaxis()->CenterTitle();
    fFit[i]->GetYaxis()->CenterTitle();
   
  }
  
  // Output file
  TString fOutName = Form("./pitz_All/Plots/EnergyEvolution");
 
  // Canvas setup
  TCanvas *C = new TCanvas("C","Evolution of Energy modulation",1500,1000);
  C->Divide(2,3);

  TLegend *Leg = new TLegend(0.17,0.20,0.3,0.60);
  PlasmaGlob::SetPaveStyle(Leg);
  Leg->SetTextAlign(22);
  Leg->SetTextColor(kGray+3);
  Leg->SetLineColor(1);
  Leg->SetBorderSize(1);
  Leg->SetFillColor(0);
  Leg->SetFillStyle(1001);
  //Leg-> SetNColumns(2);
  for(Int_t i=0;i<Nsim;i++) {
    Leg->AddEntry(gChi2[i],Form("%.1E/cc",den[i]),"L");
  }
  Leg->SetTextColor(kGray+3);

  TPaveText *text1 = new TPaveText(0.7,0.85,0.92,0.94,"NDC");
  PlasmaGlob::SetPaveTextStyle(text1);
  text1->AddText("Z #approx 84 mm");

  TPaveText *textF = new TPaveText(0.6,0.20,0.92,0.30,"NDC");
  PlasmaGlob::SetPaveTextStyle(textF);
  textF->AddText("f_{E}(#zeta) = A + B * sin(C * (#zeta-D))");
  
  // Actual Plotting!
  // ------------------------------------------------------------
  
  // More makeup
  C->cd(1);
  gPad->SetGridy(0);
  gPad->SetGridx(1);
  gPad->SetFrameLineWidth(1);  
  
  Float_t margin = (maxFitVal-minFitVal)/10.;
  TH1F *hFrameFit = new TH1F("hFrameFit","",10,-6.,0.);
  hFrameFit->GetYaxis()->SetRangeUser(24.,26.);
  hFrameFit->GetXaxis()->CenterTitle();
  hFrameFit->GetYaxis()->CenterTitle();
  hFrameFit->GetXaxis()->SetTitle("#zeta [mm]");
  hFrameFit->GetYaxis()->SetTitle("f_{E}(#zeta) [MeV]");
  
  hFrameFit->Draw("axis");
  
  for(Int_t i=0;i<Nsim;i++) {
    
    fFit[i]->Draw("AC same");

    // The spline creates a smooth curve passing through all points in the graph
    // TSpline3 *grSpline = new TSpline3("grSpline",gEmean[i]);
    // grSpline->Draw("Csame");
  }
  
  gPad->RedrawAxis("G");
  text1->Draw();
  textF->Draw();

  // More makeup
  C->cd(2);
  gPad->SetGridy(0);
  gPad->SetGridx(1);
  gPad->SetFrameLineWidth(1);  
  
  margin = (maxEmean - minEmean)/10;
  TH1F *hFrame = new TH1F("hFrame","",10,20.,170);
  hFrame->GetYaxis()->SetRangeUser(minEmean-margin,maxEmean+margin);
  hFrame->GetXaxis()->CenterTitle();
  hFrame->GetYaxis()->CenterTitle();
  hFrame->GetXaxis()->SetTitle("Z [mm]");
  hFrame->GetYaxis()->SetTitle("A [MeV]");
  
  hFrame->Draw("axis");

  for(Int_t i=0;i<Nsim;i++) {

    gEmean[i]->Draw("3");

    // The spline creates a smooth curve passing through all points in the graph
    // TSpline3 *grSpline = new TSpline3("grSpline",gEmean[i]);
    // grSpline->Draw("Csame");
  }
  
  gPad->RedrawAxis("G");
  Leg->Draw();

  C->cd(4);
  gPad->SetGridy(0);
  gPad->SetGridx(1);
  gPad->SetFrameLineWidth(1);  

  margin = (maxEamp - minEamp)/10;
  TH1F *hFrame2 = new TH1F("hFrame2","",10,20.,170);
  hFrame2->GetXaxis()->CenterTitle();
  hFrame2->GetYaxis()->CenterTitle();
  hFrame2->GetYaxis()->SetRangeUser(minEamp-margin,maxEamp+margin);
  hFrame2->GetXaxis()->SetTitle("Z [mm]");
  hFrame2->GetYaxis()->SetTitle("B [MeV]");
  hFrame2->Draw("axis");
  
  for(Int_t i=0;i<Nsim;i++) {
    gEamp[i]->Draw("CX");
    gEamp[i]->Draw("3");
  }

  gPad->RedrawAxis("G");
  
  C->cd(3);
  gPad->SetGridy(0);
  gPad->SetGridx(1);
  gPad->SetFrameLineWidth(1);
  
  margin = (maxChi2 - minChi2)/10;
  TH1F *hFrame3 = new TH1F("hFrame3","",10,20.,170);
  hFrame3->GetXaxis()->CenterTitle();
  hFrame3->GetYaxis()->CenterTitle();
  hFrame3->GetYaxis()->SetRangeUser(minChi2-margin,maxChi2+margin);
  hFrame3->GetXaxis()->SetTitle("Z [mm]");
  hFrame3->GetYaxis()->SetTitle("#chi^{2}/ndf");
  hFrame3->Draw("axis");

  for(Int_t i=0;i<Nsim;i++) {
    gChi2[i]->Draw("C");
  }
  gPad->RedrawAxis("G");

 
  C->cd(6);
  gPad->SetGridy(0);
  gPad->SetGridx(1);
  gPad->SetFrameLineWidth(1);
  
  margin = (maxLambda - minLambda)/10;
  TH1F *hFrame4 = new TH1F("hFrame4","",10,20.,170);
  hFrame4->GetXaxis()->CenterTitle();
  hFrame4->GetYaxis()->CenterTitle();
  hFrame4->GetYaxis()->SetRangeUser(minLambda-margin,maxLambda+margin);
  hFrame4->GetXaxis()->SetTitle("Z [mm]");
  hFrame4->GetYaxis()->SetTitle("C [c/#omega_{p}]");
  hFrame4->Draw("axis");

  for(Int_t i=0;i<Nsim;i++) {
    gLambda[i]->Draw("CX");
    gLambda[i]->Draw("3");
  }
  gPad->RedrawAxis("G");

 
  C->cd(5);
  gPad->SetGridy(0);
  gPad->SetGridx(1);
  gPad->SetFrameLineWidth(1);
  
  margin = (maxAsym - minAsym)/10;
  TH1F *hFrame5 = new TH1F("hFrame5","",10,20.,170);
  hFrame5->GetXaxis()->CenterTitle();
  hFrame5->GetYaxis()->CenterTitle();
  hFrame5->GetYaxis()->SetRangeUser(minAsym-margin,maxAsym+margin);
  hFrame5->GetXaxis()->SetTitle("Z [mm]");
  hFrame5->GetYaxis()->SetTitle("Asymmetry");
  hFrame5->Draw("axis");

  for(Int_t i=0;i<Nsim;i++) {
    gAsym[i]->Draw("CX");
    gAsym[i]->Draw("3");
  }
  gPad->RedrawAxis("G");

 
  // C->cd(6);
  // gPad->SetGridy(0);
  // gPad->SetGridx(1);
  // gPad->SetFrameLineWidth(1);
  
  // margin = (maxPhase - minPhase)/10;
  // TH1F *hFrame6 = new TH1F("hFrame6","",10,20.,170);
  // hFrame6->GetXaxis()->CenterTitle();
  // hFrame6->GetYaxis()->CenterTitle();
  // hFrame6->GetYaxis()->SetRangeUser(minPhase-margin,maxPhase+margin);
  // hFrame6->GetXaxis()->SetTitle("Z [mm]");
  // hFrame6->GetYaxis()->SetTitle("#psi [mm]");
  // hFrame6->Draw("axis");

  // for(Int_t i=0;i<Nsim;i++) {
  //   gPhase[i]->Draw("CX");
  //   gPhase[i]->Draw("3");
  // }
  // gPad->RedrawAxis("G");

  // Print to a file
  PlasmaGlob::imgconv(C,fOutName,opt);
  // ---------------------------------------------------------
  
  for(Int_t i=0;i<Nsim;i++)
    ifile[i]->Close();
  
}
