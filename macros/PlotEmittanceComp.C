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

void PlotEmittanceComp(const TString &options="") {
  
#ifdef __CINT__  
  gSystem->Load("libptools.so");
#endif
  
  PGlobals::Initialize();

  TString opt = options;

  const Int_t Nsim = 4;
  char sName[Nsim][56] = {"flash_v2.5kA.G.ZH.DDR.10.FinalSet.3D",
  			  "flash_v2.5kA.G.ZH.DDR.20.FinalSet.3D",
 			  "flash_v2.5kA.G.ZH.DDR.30.3D",
			  "flash_v2.5kA.G.ZH.DDR.40.3D"};

  char lName[Nsim][56] = {"(a): k_{p}^{0}#sigma_{l} = 2.5",
   			  "(b): k_{p}^{0}#sigma_{l} = 5.0",
   			  "(c): k_{p}^{0}#sigma_{l} = 7.5",
			  "(d): k_{p}^{0}#sigma_{l} = 10.0"};
  
  
  // Load first simulation data (for instance)
  PData *pData = PData::Get(sName[0]);
  //if(!pData->IsInit()) return;
  
  Float_t skd = pData->GetPlasmaSkinDepth() / PUnits::mm;
  cout << Form("\n skd = %.6f mm",skd) << endl;
  
  TFile *sFile[Nsim];
  TGraph *gEmitx[Nsim];
  TGraph *gEmity[Nsim];

  Int_t N[Nsim];
  Double_t *T[Nsim];
  Double_t *Ex[Nsim];
  Double_t *Ey[Nsim];
  
  Int_t color[Nsim];

  PPalette *pal = new PPalette("pal");
  pal->SetPalette("oli");

  Double_t eMax = -99999;
  Double_t eMin = 99999;
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
    
    TString filename = Form("./%s/Plots/Bunch/plasma/Bunch-Evolution-%s.root",sName[i],sName[i]);
    sFile[i] = new TFile(filename,"READ");
    
    // Load histos and graphs
    gEmitx[i] = (TGraph*) sFile[i]->Get("gEmitxcorevsTime");
    Ex[i] = gEmitx[i]->GetY();

    gEmitx[i]->SetLineWidth(3);
    gEmitx[i]->SetLineColor(color[i]);

    gEmity[i] = (TGraph*) sFile[i]->Get("gEmitycorevsTime");
    Ey[i] = gEmity[i]->GetY();

    gEmity[i]->SetLineWidth(3);
    gEmity[i]->SetLineColor(color[i]);
    gEmity[i]->SetLineStyle(3);
    
    T[i]  = gEmitx[i]->GetX();
    N[i]  = gEmitx[i]->GetN();
    
    for(Int_t j=0;j<N[i];j++) {
      
      if( Ex[i][j] > eMax) eMax =  Ex[i][j];
      if( Ex[i][j] < eMin) eMin =  Ex[i][j];

      if( Ey[i][j] > eMax) eMax =  Ey[i][j];
      if( Ey[i][j] < eMin) eMin =  Ey[i][j];

    }

    
  }
  
  // Canvas setup
  // Create the canvas and the pads before the Frame loop
  // Resolution:
  Int_t sizex = 800;
  Int_t sizey = 500;

  char cName[32];
  sprintf(cName,"C");     
  TCanvas *C = (TCanvas*) gROOT->FindObject(cName);
  if(C==NULL) C = new TCanvas("C","",sizex,sizey);
  C->SetFillStyle(4000);
  C->cd();
  C->Clear();

  Float_t zmin = 999.0;
  Float_t zmax = -999.0;
  for(Int_t i=0;i<Nsim;i++) {

    if(T[i][0]<zmin) zmin = T[i][0];
    if(T[i][N[i]-1]>zmax) zmax = T[i][N[i]-1];
  }

  zmin = 5;
  zmax = 75;
   
  // Setup Pad layout: 
  Int_t NPad = 1;
  TPad **pad = new TPad*[NPad];
  TH1F **hFrame = new TH1F*[NPad];
  
  Float_t lMargin = 0.12;
  Float_t rMargin = 0.04;
  Float_t bMargin = 0.15;
  Float_t tMargin = 0.04;
  Float_t factor = 0.2;    

  // Define the frames for plotting
  Int_t fonttype = 43;
  Int_t fontsize = 24;
  Int_t tfontsize = 26;
  Float_t txoffset = 1.2;
  Float_t lxoffset = 0.02;
  Float_t tyoffset = 0.8;
  Float_t lyoffset = 0.01;
  Float_t tylength = 0.015;
  Float_t txlength = 0.02;

  PGlobals::CanvasAsymPartition(C,NPad,lMargin,rMargin,bMargin,tMargin);
  for(Int_t i=0;i<NPad;i++) {
    char name[16];
    sprintf(name,"pad_%i",i);
    pad[i] = (TPad*) gROOT->FindObject(name);
    pad[i]->SetFrameLineWidth(2);  
    pad[i]->SetTickx(1);
    pad[i]->SetTicky(1);

    sprintf(name,"hFrame_%i",i);
    hFrame[i] = new TH1F(name,"",10,zmin,zmax);
    
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
  Double_t emin =  eMin-(eMax-eMin)*mfactor;
  Double_t emax =  eMax+(eMax-eMin)*mfactor;
  emin = 0;
  emax = 1.99;
  Double_t erange = emax - emin;
  Double_t zrange = zmax - zmin;

  TLegend *Leg = new TLegend(zmin + 0.02 * zrange, emax - 0.3*erange , zmin + 0.4 * zrange, emax - 0.05 * erange,"","tl");
  PGlobals::SetPaveStyle(Leg);
  Leg->SetTextAlign(12);
  Leg->SetTextColor(kGray+3);
  Leg->SetTextFont(43);
  Leg->SetTextSize(16);
  Leg->SetLineColor(1);
  Leg->SetBorderSize(0);
  Leg->SetFillColor(0);
  Leg->SetFillStyle(1001);
  Leg->SetFillStyle(0); // Hollow  

  Int_t ipad = 0;

  C->cd(0);
  pad[ipad]->Draw();
  pad[ipad]->cd();

  hFrame[ipad]->GetYaxis()->SetRangeUser(emin,emax);
  
  hFrame[ipad]->GetXaxis()->SetTitle("z [mm]");
  hFrame[ipad]->GetYaxis()->SetTitle("#varepsilon_{n} [#mum]");

  hFrame[ipad]->Draw("AXIS");
  
  TLine *lineZero = new TLine(zmin,0.0,zmax,0.0);
  lineZero->SetLineColor(kGray+2);
  lineZero->SetLineStyle(2);
  // lineZero->Draw();

  for(Int_t i=0;i<Nsim;i++) {
    gEmitx[i]->Draw("L");
    gEmity[i]->Draw("L");
    Leg->AddEntry(gEmitx[i],lName[i],"L");
  }
  
  Leg->Draw();

  gPad->Update();
  TBox *lFrame = new TBox(gPad->GetUxmin(), gPad->GetUymin(),
			  gPad->GetUxmax(), gPad->GetUymax());
  lFrame->SetFillStyle(0);
  lFrame->SetLineColor(kBlack);
  lFrame->SetLineWidth(2);
  lFrame->Draw();
  
  gPad->RedrawAxis("g");

  
  ipad++;
  
  C->cd(0);
  C->cd();
  
  // Print to file -------------------------------------------

  TString fOutName = Form("./EmittanceCompDDR/EmittanceCompDDR");
  PGlobals::imgconv(C,fOutName,opt);

  // ---------------------------------------------------------
  
  
}


