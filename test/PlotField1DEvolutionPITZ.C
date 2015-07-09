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

#include "PData.hh"
#include "PlasmaGlob.hh"

void PlotField1DEvolution(const TString &opt="png") { 
  
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
  gStyle->SetPadGridX(0);
  
  // Get Files from different PITZ simulations
  const Int_t Nsim = 8;
  TString simNames[Nsim] = {"pitz_FT01","pitz_FT02","pitz_FT04","pitz_FT06","pitz_FT08","pitz_FTHD","pitz_FT2","pitz_FT4"}; // ,"pitz_FT3D"};//,"pitz_FT6","pitz_FT8"};
  Int_t simColors[Nsim] = {kGray,kGray+1,kGray+2,kOrange-2,kOrange+7,kRed,kAzure-4,kAzure-5}; // ,kGreen};//,kAzure-9,kGreen+1};
  Float_t n0[Nsim] = {0.1e15,0.2e15,0.4e15,0.6e15,0.8e15,1.0e15,2.0e15,4.0e15}; // ,1.0e15};//,6.0e15,8.0e15};

  TString filename[Nsim];
  TFile  *ifile[Nsim];
  TGraph *gE0MaxPos[Nsim];
  // TGraph *gE1MaxPos[Nsim];
  TGraph *gE0MaxValue[Nsim];
  // TGraph *gE1MaxValue[Nsim];
  TGraph *gE0Vel[Nsim];
  TGraph *gVel[Nsim];
  TGraph *gVelT[Nsim];
  TGraph *gPhaseT[Nsim];
  
  Float_t maxEPhase = -999.;
  Float_t minEPhase = 999.;
  Float_t maxField = -999.;
  Float_t minField = 999.;
  Float_t maxEmean = -999.;
  Float_t minEmean = 999.;  
  Float_t maxVel = -999.;
  Float_t minVel = 999.;  
  Float_t maxVelT = -999.;
  Float_t minVelT = 999.;  
  Float_t maxPhaseT = -999.;
  Float_t minPhaseT = 999.;  

  // beam parameters
  Double_t nb = 1.05e13 / PUnits::cm3;
  Double_t r0 = 34.259 * PUnits::um;
  Double_t E = 25.0 * PUnits::MeV;
  Double_t gamma = E / PConst::ElectronMassE;  
  Double_t beta  = TMath::Sqrt(1. - (1./(gamma*gamma)));
  
  cout << " beam density = " << nb * PUnits::cm3 << " e/cm3 " << endl;
  cout << " beta beam = " << beta << endl;

  for(Int_t i=0;i<Nsim;i++) {
    filename[i] = Form("./%s/Plots/Field1D/Field1D-%s.root",simNames[i].Data(),simNames[i].Data());
    ifile[i] = (TFile*) gROOT->GetListOfFiles()->FindObject(filename[i].Data());
    if (!ifile[i]) ifile[i] = new TFile(filename[i],"READ");
    
    // Get the needed Graphs for every simulation on the list 
    gE0MaxPos[i] = (TGraph*) ifile[i]->Get("gEMaxPos_0");
    // gE1MaxPos[i] = (TGraph*) ifile[i]->Get("gEMaxPos_1");
    gE0MaxValue[i] = (TGraph*) ifile[i]->Get("gEMaxValue_0");
    // gE1MaxValue[i] = (TGraph*) ifile[i]->Get("gEMaxValue_1");
    
    // Make relative phases
    Double_t *x = gE0MaxPos[i]->GetX();
    Double_t *y = gE0MaxPos[i]->GetY();
    Double_t skinDepth = PFunc::PlasmaSkindepth(n0[i]/PUnits::cm3) / PUnits::mm;
   
    // normalize variables and parameters
    Double_t n0n = 1.0;
    Double_t nbn = nb / ( n0[i] / PUnits::cm3 );
    Double_t r0n = (r0 / PUnits::mm) / skinDepth;

    // Create graphs using theoretical expressions
    char sName[64];    
    Int_t Np = gE0MaxPos[i]->GetN();

    gVel[i] = new TGraph(Np-1);
    sprintf(sName,"gVel_%i",i);
    gVel[i]->SetName(sName);

    gVelT[i] = new TGraph(Np);
    sprintf(sName,"gVelT_%i",i);
    gVelT[i]->SetName(sName);

    gPhaseT[i] = new TGraph(Np);
    sprintf(sName,"gPhaseT_%i",i);
    gPhaseT[i]->SetName(sName);
    
    // Double_t skinDepth = 1.;
    // Change the Graph to plasma units:
    for(Int_t j=0;j<Np;j++) {
      
      Double_t vel = 0.;
      if(j+1<Np) {
	vel = beta + (y[j+1] - y[j]) / (x[j+1] -x[j]) ;
	gVel[i]->SetPoint(j,x[j],vel);
      }
      
      // phase velocity of the corresponding point
      // cout << Form("  z = %5.2f  zg = %5.2f  gamma = %6.4f  nb = %6.4f  n0 = %5.2f",x[j]/skinDepth,y[j]/skinDepth,gamma,nbn,n0n);
      // Double_t vph = PFunc::PhaseVelocity2(x[j]/skinDepth,y[j]/skinDepth,gamma,nbn,n0n);
      Double_t vph = PFunc::PhaseVelocity(x[j]/skinDepth,y[j]/skinDepth,r0n,gamma,nbn,n0n);
      gVelT[i]->SetPoint(j,x[j],vph);
      // cout << " velocity = " << vph << endl;
      if(j==0) 
	gPhaseT[i]->SetPoint(j,x[j],0.);
      else {
	Double_t *v  = gVelT[i]->GetY();
	Double_t *ph = gPhaseT[i]->GetY();
	Float_t dz = (x[j] - x[j-1])/skinDepth;
	// Adding velocities:
	// Float_t vprime = (v[j-1] - beta) / (1 - v[j-1]*beta);
	Float_t dphase = ph[j-1] + (v[j-1] - beta) * dz ;
	gPhaseT[i]->SetPoint(j,x[j],dphase);
      }
      
    }
    
    Double_t *xE0MaxPos = gE0MaxPos[i]->GetX();
    Double_t *yE0MaxPos = gE0MaxPos[i]->GetY();
    // Double_t *xE1MaxPos = gE1MaxPos[i]->GetX();
    // Double_t *yE1MaxPos = gE1MaxPos[i]->GetY();
    
    // Shift to relative phases
    Double_t yE0MaxPos0 = yE0MaxPos[0];
    //  Double_t yE1MaxPos0 = yE1MaxPos[0];
    for(Int_t j=0;j<gE0MaxPos[i]->GetN();j++) {
      gE0MaxPos[i]->SetPoint(j,xE0MaxPos[j],(yE0MaxPos[j]-yE0MaxPos0)/skinDepth);
      // gE1MaxPos[i]->SetPoint(j,xE1MaxPos[j],(yE1MaxPos[j]-yE1MaxPos0)/skinDepth);
    }

    // Calculate extremes for the set of Graphs
    yE0MaxPos = gE0MaxPos[i]->GetY();
    // yE1MaxPos = gE1MaxPos[i]->GetY();
    Double_t *yE0MaxValue = gE0MaxValue[i]->GetY();
    // Double_t *yEmeanMaxValue = gE1MaxValue[i]->GetY();
    Double_t *yVel = gVel[i]->GetY();
    Double_t *yVelT = gVelT[i]->GetY();
    Double_t *yPhaseT = gPhaseT[i]->GetY();
    
    for(Int_t j=0;j<gE0MaxPos[i]->GetN();j++) {
      if(yE0MaxPos[j]>maxEPhase)
	maxEPhase = yE0MaxPos[j];
      if(yE0MaxPos[j]<minEPhase)
	minEPhase = yE0MaxPos[j];
      
      if(yE0MaxValue[j]>maxField)
	maxField = yE0MaxValue[j];
      if(yE0MaxValue[j]<minField)
	minField = yE0MaxValue[j];
      
      // if(yEmeanMaxValue[j]>maxEmean)
      // 	maxEmean = yEmeanMaxValue[j];
      // if(yEmeanMaxValue[j]<minEmean)
      // 	minEmean = yEmeanMaxValue[j];
      
      if(yVel[j]>maxVel)
	maxVel = yVel[j];
      if(yVel[j]<minVel)
	minVel = yVel[j];

      if(yVelT[j]>maxVelT)
	maxVelT = yVelT[j];
      if(yVelT[j]<minVelT)
	minVelT = yVelT[j];

      if(yPhaseT[j]>maxPhaseT)
	maxPhaseT = yPhaseT[j];
      if(yPhaseT[j]<minPhaseT)
	minPhaseT = yPhaseT[j];
    }

    // Graph's attributes
    gE0MaxPos[i]->SetLineColor(simColors[i]);
    gE0MaxPos[i]->SetLineWidth(2);
    gE0MaxPos[i]->SetMarkerColor(simColors[i]);
    gE0MaxPos[i]->SetMarkerSize(0.2);
    gE0MaxPos[i]->GetXaxis()->SetTitle("Z [mm]");
    gE0MaxPos[i]->GetXaxis()->CenterTitle();
    gE0MaxPos[i]->GetYaxis()->SetTitle("#Delta#zeta [mm]");
    gE0MaxPos[i]->GetYaxis()->CenterTitle();
    
    // gE1MaxPos[i]->SetLineColor(simColors[i]);
    // gE1MaxPos[i]->SetLineWidth(2);
    // gE1MaxPos[i]->SetMarkerColor(simColors[i]);
    // gE1MaxPos[i]->SetMarkerSize(0.2);
    // gE1MaxPos[i]->GetXaxis()->CenterTitle();
    // gE1MaxPos[i]->GetYaxis()->CenterTitle();

    gE0MaxValue[i]->SetLineColor(simColors[i]);
    gE0MaxValue[i]->SetLineWidth(2);
    gE0MaxValue[i]->SetMarkerColor(simColors[i]);
    gE0MaxValue[i]->SetMarkerSize(0.2);
    gE0MaxValue[i]->GetXaxis()->SetTitle("Z [mm]");
    gE0MaxValue[i]->GetXaxis()->CenterTitle();
    gE0MaxValue[i]->GetYaxis()->SetTitle("E_{min} [GeV/m]");
    gE0MaxValue[i]->GetYaxis()->CenterTitle();

    // gE1MaxValue[i]->SetLineColor(simColors[i]);
    // gE1MaxValue[i]->SetLineWidth(2);
    // gE1MaxValue[i]->SetMarkerColor(simColors[i]);
    // gE1MaxValue[i]->SetMarkerSize(0.2);
    // gE1MaxValue[i]->GetXaxis()->SetTitle("Z [mm]");
    // gE1MaxValue[i]->GetXaxis()->CenterTitle();
    // gE1MaxValue[i]->GetYaxis()->SetTitle("E_{min} [MeV]");
    // gE1MaxValue[i]->GetYaxis()->CenterTitle();

    gVel[i]->SetLineColor(simColors[i]);
    gVel[i]->SetLineWidth(2);
    gVel[i]->SetMarkerColor(simColors[i]);
    gVel[i]->SetMarkerSize(0.2);
    gVel[i]->GetXaxis()->SetTitle("Z [mm]");
    gVel[i]->GetXaxis()->CenterTitle();
    gVel[i]->GetYaxis()->SetTitle("v_{ph} [c]");
    gVel[i]->GetYaxis()->CenterTitle();

    gVelT[i]->SetLineColor(simColors[i]);
    gVelT[i]->SetLineWidth(1);
    gVelT[i]->SetLineStyle(7);
    gVelT[i]->SetMarkerColor(simColors[i]);
    gVelT[i]->SetMarkerSize(0.2);
    gVelT[i]->GetXaxis()->SetTitle("Z [mm]");
    gVelT[i]->GetXaxis()->CenterTitle();
    gVelT[i]->GetYaxis()->SetTitle("v_{ph} [c]");
    gVelT[i]->GetYaxis()->CenterTitle();

    gPhaseT[i]->SetLineColor(simColors[i]);
    gPhaseT[i]->SetLineWidth(1);
    gPhaseT[i]->SetLineStyle(7);
    gPhaseT[i]->SetMarkerColor(simColors[i]);
    gPhaseT[i]->SetMarkerSize(0.2);
    gPhaseT[i]->GetXaxis()->SetTitle("Z [mm]");
    gPhaseT[i]->GetXaxis()->CenterTitle();
    gPhaseT[i]->GetYaxis()->SetTitle("v_{ph} [c]");
    gPhaseT[i]->GetYaxis()->CenterTitle();
  }
  
  // Output file
  TString fOutName = Form("./pitz_All/Plots/FieldEvolution");
 
  // Canvas setup
  TCanvas *C = new TCanvas("C","Evolution of Electric fields",850,1000);
  C->Divide(1,2);

  TLegend *Leg = new TLegend(0.17,0.20,0.35,0.60);
  PlasmaGlob::SetPaveStyle(Leg);
  Leg->SetTextAlign(22);
  Leg->SetTextColor(kGray+3);
  Leg->SetLineColor(kGray+1);
  Leg->SetBorderSize(1);
  Leg->SetFillColor(0);
  Leg->SetFillStyle(1001);
  //Leg-> SetNColumns(2);
  for(Int_t i=0;i<Nsim;i++) {
    if(!simNames[i].Contains("3D"))
      Leg->AddEntry(gE0MaxPos[i],Form("%2.1f #times 10^{15} cm^{-3}",n0[i]/1e15),"L");
    else
      Leg->AddEntry(gE0MaxPos[i],Form("%2.1f #times 10^{15} cm^{-3} 3D",n0[i]/1e15),"L");  
  }
  Leg->SetTextColor(kGray+3);


  // Actual Plotting!
  // ------------------------------------------------------------
  
  // More makeup

  C->cd(1);
  gPad->SetGridy(0);
  gPad->SetGridx(0);
  gPad->SetFrameLineWidth(2);  
  
  Float_t margin = (maxEPhase - minEPhase)/10;
  TH1F *hFrame = new TH1F("hFrame","",10,0,170);
  Float_t minPhase = minEPhase;
  // if(minPhaseT<minPhase) minPhase = minPhaseT;
  Float_t maxPhase = maxEPhase;
  // if(maxPhaseT<maxPhase) maxPhase = maxPhaseT;

  hFrame->GetYaxis()->SetRangeUser(minPhase-margin,maxPhase+margin);
  hFrame->GetXaxis()->CenterTitle();
  hFrame->GetYaxis()->CenterTitle();
  hFrame->GetXaxis()->SetTitle("Z [mm]");
  hFrame->GetYaxis()->SetTitle("#Delta#zeta [c/#omega_{p}]");
  
  hFrame->Draw("axis");

  // De-phasing lines:
  TLine *lpi0 = new TLine(hFrame->GetXaxis()->GetXmin(),0.,hFrame->GetXaxis()->GetXmax(),0.);
  lpi0->SetLineColor(kGray+1);
  lpi0->SetLineStyle(2);
  lpi0->Draw();
  TLine *lpi = new TLine(hFrame->GetXaxis()->GetXmin(),-TMath::Pi(),hFrame->GetXaxis()->GetXmax(),-TMath::Pi());
  lpi->SetLineColor(kGray+1);
  lpi->SetLineStyle(2);
  lpi->Draw();
  TLine *lpi2 = new TLine(hFrame->GetXaxis()->GetXmin(),-TMath::TwoPi(),hFrame->GetXaxis()->GetXmax(),-TMath::TwoPi());
  lpi2->SetLineColor(kGray+1);
  lpi2->SetLineStyle(2);
  lpi2->Draw();
  for(Int_t i=0;i<Nsim;i++) {

    gE0MaxPos[i]->Draw("L");

    // The spline creates a smooth curve passing through all points in the graph
    // TSpline3 *grSpline = new TSpline3("grSpline",gE0MaxPos[i]);
    // grSpline->Draw("Csame");

    //  gPhaseT[i]->Draw("L");
  
  }
  
  gPad->RedrawAxis("G");

  C->cd(2);
  gPad->SetGridy(0);
  gPad->SetGridx(0);
  gPad->SetFrameLineWidth(2);  

  margin = (maxField - minField)/10;
  TH1F *hFrame2 = new TH1F("hFrame2","",10,0,170);
  hFrame2->GetXaxis()->CenterTitle();
  hFrame2->GetYaxis()->CenterTitle();
  hFrame2->GetYaxis()->SetRangeUser(minField-margin,maxField+margin);
  hFrame2->GetXaxis()->SetTitle("Z [mm]");
  hFrame2->GetYaxis()->SetTitle("E_{z}^{min} [GV/m]");
  hFrame2->GetYaxis()->SetNdivisions(505);
  hFrame2->Draw("axis");

  TLine *line0 = new TLine(hFrame2->GetXaxis()->GetXmin(),0.,hFrame2->GetXaxis()->GetXmax(),0.);
  line0->SetLineColor(kGray+1);
  line0->SetLineStyle(2);
  line0->Draw();
  
  for(Int_t i=0;i<Nsim;i++) {
    gE0MaxValue[i]->Draw("L");
  }
  gPad->RedrawAxis("G");
  
  // Leg->Draw();

  // C->cd(3);
  // gPad->SetGridy(0);
  // gPad->SetGridx(1);
  // gPad->SetFrameLineWidth(2);
  
  // Float_t margin = (maxVelT - minVelT)/10;
  // TH1F *hFrame3 = new TH1F("hFrame3","",10,0,170);
  // hFrame3->GetXaxis()->CenterTitle();
  // hFrame3->GetYaxis()->CenterTitle();
  // hFrame3->GetYaxis()->SetRangeUser(minVelT-margin,maxVelT+margin);
  // hFrame3->GetXaxis()->SetTitle("Z [mm]");
  // hFrame3->GetYaxis()->SetTitle("v_{ph} [c]");
  // hFrame3->Draw("axis");
  
  // for(Int_t i=0;i<Nsim;i++) {
  //   gVel[i]->Draw("L");
  //   gVelT[i]->Draw("L");
  // }
  // gPad->RedrawAxis("G");


  C->cd(0);
  // Print to a file
  PlasmaGlob::imgconv(C,fOutName,opt);
  // ---------------------------------------------------------
  
  for(Int_t i=0;i<Nsim;i++)
    ifile[i]->Close();
  
}
