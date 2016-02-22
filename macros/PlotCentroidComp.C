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
#include <TF1.h>
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

using namespace std;

Double_t Density(Double_t *x, Double_t *par) {
  Double_t f;
  Double_t skd = par[0];
  Double_t length = (391.4 - 15.1) * skd;
  Double_t slope = 1.0/length;
  if(x[0] < -length) f = 0;
  else if(x[0] < 0.0) f = slope * (x[0]+length);
  else f = 1.0;
  return f;
}


void PlotCentroidComp(const TString &options="") {

#ifdef __CINT__  
  gSystem->Load("libptools.so");
#endif

  PGlobals::Initialize();

  TString opt = options;

  const Int_t Nsim = 4;
  char sName[Nsim][36] = {"flash-500pC_S2E_R56-max.3D", 
			  "flash-500pC_S2E_R56-max.Be2mm.3D", 
			  "flash-500pC_S2E_R56-max.Be4mm.3D", 
			  "flash-500pC_S2E_R56-max.Be6mm.3D"};

  char lName[Nsim][36] = {"No spoiler", 
			  "Spoiler (Be 2mm) ", 
			  "Spoiler (Be 4mm) ", 
			  "Spoiler (Be 6mm)"};

  
  // Load first simulation data (for instance)
  PData *pData = PData::Get(sName[0]);
  //if(!pData->IsInit()) return;
  
  Float_t skd = pData->GetPlasmaSkinDepth() / PUnits::mm;
  cout << Form("\n skd = %.6f mm",skd) << endl;
  
  TFile *sFile[Nsim];
  TGraph *gX[Nsim];
  TGraph *gXrms[Nsim];
  TGraph *gXup[Nsim];
  TGraph *gXdo[Nsim];
  TGraph *gY[Nsim];
  TGraph *gYrms[Nsim];
  TGraph *gYup[Nsim];
  TGraph *gYdo[Nsim];

  Int_t N[Nsim];
  Double_t *T[Nsim];
  Double_t *X[Nsim];
  Double_t *Xrms[Nsim];
  Double_t *Xup[Nsim]; 
  Double_t *Xdo[Nsim];
  Double_t *Y[Nsim];
  Double_t *Yrms[Nsim];
  Double_t *Yup[Nsim]; 
  Double_t *Ydo[Nsim];

  
  Int_t color[Nsim];

  PPalette *pal = new PPalette("pal");
  pal->SetPalette("oli");

  Double_t Max = -99999;
  Double_t Min = 99999;
  
  for(Int_t i=0;i<Nsim;i++) {

    Int_t index =  i * pal->GetNColors() / Nsim;
    color[i] = pal->GetColorIndex(index);

    TString filename = Form("./%s/Plots/Bunch/beam-driver/Bunch-Evolution-%s.root",sName[i],sName[i]);
    sFile[i] = new TFile(filename,"READ");
    
    // Load histos and graphs
    gX[i]    = (TGraph*) sFile[i]->Get("gXmeanvsTime"); 
    gXrms[i] = (TGraph*) sFile[i]->Get("gXrmsvsTime");

    T[i] =  gX[i]->GetX();
    X[i] =  gX[i]->GetY();
    Xrms[i] = gXrms[i]->GetY();
    N[i] = gX[i]->GetN();

    Xup[i] = new Double_t[N[i]];
    Xdo[i] = new Double_t[N[i]];
  
    for(Int_t j=0;j<N[i];j++) {
      
      Xup[i][j] = X[i][j] + Xrms[i][j];
      Xdo[i][j] = X[i][j] - Xrms[i][j];

      if( Xup[i][j] > Max) Max =  Xup[i][j];
      if( Xdo[i][j] < Min) Min =  Xdo[i][j];
      
    }

    gXup[i] = new TGraph(N[i],T[i],Xup[i]);
    gXdo[i] = new TGraph(N[i],T[i],Xdo[i]);

    gXup[i]->SetLineWidth(1);
    gXup[i]->SetLineColor(color[i]);
    gXup[i]->SetLineStyle(2);
    
    gXdo[i]->SetLineWidth(1);
    gXdo[i]->SetLineColor(color[i]);
    gXdo[i]->SetLineStyle(2);

    gX[i]->SetLineWidth(3);
    gX[i]->SetLineColor(color[i]);

    gY[i]    = (TGraph*) sFile[i]->Get("gYmeanvsTime"); 
    gYrms[i] = (TGraph*) sFile[i]->Get("gYrmsvsTime");

    T[i] =  gY[i]->GetX();
    Y[i] =  gY[i]->GetY();
    Yrms[i] = gYrms[i]->GetY();
    N[i] = gY[i]->GetN();

    Yup[i] = new Double_t[N[i]];
    Ydo[i] = new Double_t[N[i]];
  
    for(Int_t j=0;j<N[i];j++) {
      
      Yup[i][j] = Y[i][j] + Yrms[i][j];
      Ydo[i][j] = Y[i][j] - Yrms[i][j];

      if( Yup[i][j] > Max) Max =  Yup[i][j];
      if( Ydo[i][j] < Min) Min =  Ydo[i][j];
      
    }

    gYup[i] = new TGraph(N[i],T[i],Yup[i]);
    gYdo[i] = new TGraph(N[i],T[i],Ydo[i]);

    gYup[i]->SetLineWidth(1);
    gYup[i]->SetLineColor(color[i]);
    gYup[i]->SetLineStyle(2);
    
    gYdo[i]->SetLineWidth(1);
    gYdo[i]->SetLineColor(color[i]);
    gYdo[i]->SetLineStyle(2);

    gY[i]->SetLineWidth(3);
    gY[i]->SetLineColor(color[i]);

  }
  
  // Canvas setup
  // Create the canvas and the pads before the Frame loop
  // Resolution:
  Int_t sizex = 600;
  Int_t sizey = 600;
  char cName[32];
  sprintf(cName,"C");     
  TCanvas *C = (TCanvas*) gROOT->FindObject(cName);
  if(C==NULL) C = new TCanvas("C","",sizex,sizey);
  C->SetFillStyle(4000);
  C->cd();
  C->Clear();

  TLine *lineZero = new TLine(T[0][0],0.0,T[0][N[0]-1],0.0);
  lineZero->SetLineColor(kGray+2);
  lineZero->SetLineStyle(2);
 
  // Setup Pad layout: 
  const Int_t NPad = 3;
  TPad *pad[NPad];
  TH1F *hFrame[NPad];

  Float_t lMargin = 0.12;
  Float_t rMargin = 0.04;
  Float_t bMargin = 0.10;
  Float_t tMargin = 0.04;
  Float_t factor = 0.2;    
  PGlobals::CanvasAsymPartition(C,NPad,lMargin,rMargin,bMargin,tMargin,factor,0.01);

  // Define the frames for plotting
  Int_t fonttype = 43;
  Int_t fontsize = 18;
  Int_t tfontsize = 20;
  Float_t txoffset = 2.2;
  Float_t lxoffset = 0.02;
  Float_t tyoffset = 1.4;
  Float_t lyoffset = 0.01;
  Float_t tylength = 0.015;
  Float_t txlength = 0.02;
  for(Int_t i=0;i<NPad;i++) {
    char name[16];
    sprintf(name,"pad_%i",i);
    pad[i] = (TPad*) gROOT->FindObject(name);
    pad[i]->SetFrameLineWidth(2);  
    pad[i]->SetTickx(1);
    pad[i]->SetTicky(1);

    sprintf(name,"hFrame_%i",i);

    hFrame[i] = new TH1F(name,"",10,T[0][0],T[0][N[0]-1]);
    
    Float_t xFactor = pad[0]->GetAbsWNDC()/pad[i]->GetAbsWNDC();
    Float_t yFactor = pad[0]->GetAbsHNDC()/pad[i]->GetAbsHNDC();

    // Format for y axis
    hFrame[i]->GetYaxis()->SetTitleFont(fonttype);
    hFrame[i]->GetYaxis()->SetTitleSize(tfontsize);
    hFrame[i]->GetYaxis()->SetTitleOffset(tyoffset);
    hFrame[i]->GetYaxis()->SetLabelFont(fonttype);
    hFrame[i]->GetYaxis()->SetLabelSize(fontsize);
    hFrame[i]->GetYaxis()->SetLabelOffset(lyoffset);
    hFrame[i]->GetYaxis()->SetTickLength(xFactor*tylength/yFactor);
    hFrame[i]->GetYaxis()->CenterTitle();
    hFrame[i]->GetYaxis()->SetNdivisions(505);

    // Format for x axis
    hFrame[i]->GetXaxis()->SetTitleFont(fonttype);
    hFrame[i]->GetXaxis()->SetTitleSize(tfontsize+2);
    hFrame[i]->GetXaxis()->SetTitleOffset(txoffset);
    hFrame[i]->GetXaxis()->SetLabelFont(fonttype);
    hFrame[i]->GetXaxis()->SetLabelSize(fontsize+2);
    hFrame[i]->GetXaxis()->SetLabelOffset(lxoffset);
    hFrame[i]->GetXaxis()->CenterTitle();
    hFrame[i]->GetXaxis()->SetTickLength(yFactor*txlength/xFactor);      
  }

  Double_t mfactor = 0.1;
  if(Max>-Min) Min = -Max;
  if(-Min>Max) Max = -Min;

  Double_t xrange = T[0][N[0]-1] - T[0][0];
  Double_t yrange = Max - Min;
  TLegend *Leg = new TLegend(T[0][0] + 0.02 * xrange, Max - 0.3*yrange , T[0][0] + 0.4 * xrange, Max - 0.01 * yrange,"","tl");
  PGlobals::SetPaveStyle(Leg);
  Leg->SetTextAlign(12);
  Leg->SetTextColor(kGray+3);
  Leg->SetTextFont(42);
  Leg->SetLineColor(1);
  Leg->SetBorderSize(0);
  Leg->SetFillColor(0);
  Leg->SetFillStyle(1001);
  Leg->SetFillStyle(0); // Hollow  

  C->cd(0);
  pad[0]->Draw();
  pad[0]->cd();

  hFrame[0]->GetYaxis()->SetRangeUser(Min-(Max-Min)*mfactor,Max+(Max-Min)*mfactor);  
  
  hFrame[0]->GetXaxis()->SetTitle("z [mm]");
  hFrame[0]->GetYaxis()->SetTitle("y [#mum]");

  hFrame[0]->Draw("AXIS");

  lineZero->Draw();

  for(Int_t i=0;i<Nsim;i++) {
    gY[i]->Draw("L");
    gYup[i]->Draw("L");
    gYdo[i]->Draw("L");
  }
  
  C->cd(0);
  pad[1]->Draw();
  pad[1]->cd();

  hFrame[1]->GetYaxis()->SetRangeUser(Min-(Max-Min)*mfactor,Max+(Max-Min)*mfactor);  
  
  hFrame[1]->GetXaxis()->SetTitle("z [mm]");
  hFrame[1]->GetYaxis()->SetTitle("x [#mum]");

  hFrame[1]->Draw("AXIS");

  lineZero->Draw();

  for(Int_t i=0;i<Nsim;i++) {
    gX[i]->Draw("L");
    gXup[i]->Draw("L");
    gXdo[i]->Draw("L");

    Leg->AddEntry(gX[i],lName[i],"L");
  }
  Leg->Draw();
  
  C->cd(0);
  pad[2]->Draw();
  pad[2]->cd();

  TF1 *fDen = new TF1("fDen",Density,T[0][0],T[0][N[0]-1],1);
  fDen->SetParameter(0,skd);
  hFrame[2]->GetYaxis()->SetRangeUser(0.01,1.4);  
  
  // hFrame[2]->GetXaxis()->SetTitle("n_{He} [10^{15} e/cm^{3}]");
  hFrame[2]->GetXaxis()->SetTitle("");
  hFrame[2]->GetYaxis()->SetTitle("n/n_{0}");
  
  hFrame[2]->Draw("AXIS");

  lineZero->Draw();

  fDen->SetLineColor(kGray+2);
  fDen->SetLineWidth(2);
  fDen->Draw("same C");
  
  // Print to file --------------------------------------
  
  C->cd();
  
  // Print to a file
  // Output file
  TString fOutName = Form("./CentroidComp/CentroidComp");
  //fOutName += Form("_%i",time);
  
  PGlobals::imgconv(C,fOutName,opt);
  // ---------------------------------------------------------
  
  
}


