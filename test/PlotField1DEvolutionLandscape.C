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
#include <TPaletteAxis.h>

#include "PData.hh"
#include "PlasmaGlob.hh"
#include "PPalette.hh"

void PlotField1DEvolutionLandscape(const TString &sim, Int_t Nosc = 1, const Int_t Nfields=1, const TString &options="png") { 
  
#ifdef __CINT__  
  gSystem->Load("libplasma.so");
#endif

  PlasmaGlob::Initialize();
  
  // Palettes!
  gROOT->Macro("PlasmaPalettes.C");
  
  // Load PData
  PData *pData = PData::Get(sim.Data());
 
  TString opt = options;

  // Some plasma constants
  Float_t n0 = pData->GetPlasmaDensity();
  Float_t skindepth = pData->GetPlasmaSkinDepth();
  Float_t E0 = pData->GetPlasmaE0();
 
  // More makeup            
  Float_t margins[4] = {0.15,0.25,0.05,0.05};
  gStyle->SetPadLeftMargin(margins[0]);  // Margin left axis  
  gStyle->SetPadBottomMargin(margins[1]);
  gStyle->SetPadRightMargin(margins[2]);
  gStyle->SetPadTopMargin(margins[3]);  // Margin left axis  
 
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);


  if(opt.Contains("grid")) {
    gStyle->SetPadGridX(1);
    gStyle->SetPadGridY(1);
  }

  // Resolution:
  Int_t sizex = 1024;
  Int_t sizey = 320;
   
  Float_t maxEextr  = -999.;
  Float_t minEextr  = 999.;
 
  // In every oscillation there are 2 crossings
  Int_t Ncross = 2 * Nosc;

  TGraph *gEextr[Nfields][30];

  TString filename;
  filename = Form("./%s/Plots/Field1DEvolutions/Field1DEvolutions-%s.root",sim.Data(),sim.Data());
  
  TFile  *ifile = (TFile*) gROOT->GetListOfFiles()->FindObject(filename.Data());
  if (!ifile) ifile = new TFile(filename,"READ");

  for(Int_t i=0;i<Nfields;i++) {
  
    char hName[24];
    for(Int_t ic=0;ic<Ncross;ic++) {
      char gName[24];
      
      sprintf(gName,"gEextr_%i_%i",i,ic); 
      gEextr[i][ic] = (TGraph*) ifile->Get(gName);

      if( !gEextr[i][ic] ) continue;

      // if(i==0 && ic==0) continue;
      // if(i==1 && ic==1) continue;
      
      // Calculate the max and min of every set of graphs:
      
      Int_t Npoints = gEextr[i][ic]->GetN();
      Double_t *yEextr  = gEextr[i][ic]->GetY();
      Double_t *xEextr   = gEextr[i][ic]->GetX();

      if(opt.Contains("units")) {
	for(Int_t ip=0;ip<Npoints;ip++) {
	  yEextr[ip]  *= E0 / (PUnits::GV/PUnits::m);
	  xEextr[ip]     *= skindepth / PUnits::um;
	}
      }

      // only takes into account the minimums of the accelerating field:
      // if(ic%2-1) continue;
     
      for(Int_t j=0;j<Npoints;j++) {
	if(yEextr[j]>maxEextr)
	  maxEextr = yEextr[j];
	if(yEextr[j]<minEextr)
	  minEextr = yEextr[j];
      }
      
    }

  }

  // Manual coloring:
  const Int_t NCOLORS = 5;
  Int_t colors[NCOLORS] = {kOrange+10,kMagenta+2,kBlue,kYellow+2,kCyan+2};
  for(Int_t i=0;i<Nfields;i++) { 
    for(Int_t ic=0;ic<Ncross;ic++) {
      Int_t index = ic/2;
      if(index>=NCOLORS) index = NCOLORS-1;
      gEextr[i][ic]->SetLineColor(colors[index]);
      gEextr[i][ic]->SetLineWidth(2);
    }
  }

  // Canvas setup
  TCanvas *C1 = new TCanvas("C1","Evolution of Electric fields",sizex,sizey);
  C1->SetFillStyle(4000);
  C1->SetFrameFillStyle(4000);

  TLegend *Leg = new TLegend(margins[0]+0.05, margins[1] + 0.05, margins[0] + 0.15, margins[1]+0.25);
  PlasmaGlob::SetPaveStyle(Leg);
  Leg->SetTextAlign(22);
  Leg->SetTextColor(kGray+3);
  Leg->SetLineColor(kGray+1);
  Leg->SetBorderSize(0);
  Leg->SetFillColor(0);
  Leg->SetFillStyle(4000);
  Leg->SetTextFont(42);
  for(Int_t i=0;i<Nfields;i++) {
    char name[10];
    if(i==0) sprintf(name,"E_{z}"); 
    if(i==1) sprintf(name,"E_{y}"); 
    Leg->AddEntry(gEextr[i][0],name,"L");
  }
  Leg->SetTextColor(kGray+3);
  
  
  C1->cd(0);
  gPad->SetFrameLineWidth(1);  
  // gPad->SetRightMargin(0.10);
  
  Float_t margin = (maxEextr - minEextr)/10;
  Float_t Eextr0   = minEextr-margin;
  Float_t Eextr1   = maxEextr+margin;

  Int_t Npoints = gEextr[0][0]->GetN();
  Double_t *xEextr   = gEextr[0][0]->GetX();
  
  TH1F *hFrame2 = new TH1F("hFrame2","",100,0.0,xEextr[Npoints-1]);
  
  hFrame2->GetYaxis()->SetRangeUser(Eextr0,Eextr1);
 
  if(opt.Contains("units")) {
    hFrame2->GetXaxis()->SetTitle("z [#mum]");
    hFrame2->GetYaxis()->SetTitle("E [GV/m]");   
  } else {
    hFrame2->GetXaxis()->SetTitle("k_{p}z");
    hFrame2->GetYaxis()->SetTitle("E/E_{0}");
  }
  PlasmaGlob::SetH1LabelSize(hFrame2);

  hFrame2->GetYaxis()->SetNdivisions(503);
 
  hFrame2->GetYaxis()->SetTitleSize(0.08);
  hFrame2->GetYaxis()->SetTitleOffset(0.8);

  hFrame2->GetYaxis()->SetLabelSize(0.08);
  hFrame2->GetYaxis()->SetLabelOffset(0.01);
  hFrame2->GetYaxis()->SetTickLength(0.02);

  hFrame2->GetYaxis()->SetTitleFont(42);

  hFrame2->GetXaxis()->SetTitleSize(0.08);
  hFrame2->GetXaxis()->SetTitleOffset(1.2);

  hFrame2->GetXaxis()->SetLabelSize(0.085);
  hFrame2->GetXaxis()->SetLabelOffset(0.02);
  hFrame2->GetXaxis()->SetTickLength(0.04);

  hFrame2->GetXaxis()->SetTitleFont(42);

  hFrame2->Draw("axis");

  TLine *lineYzero = new TLine(0.0,0.0,xEextr[Npoints-1],0.0);
  lineYzero->SetLineColor(kGray+2);
  lineYzero->SetLineStyle(2);
  lineYzero->Draw();

  Float_t HeIon  = 92.75;  // GV/m
  Float_t HeIon2 = 234.96;  // GV/m
  if(!opt.Contains("units")) {
    HeIon  /= ( E0 / (PUnits::GV/PUnits::m));
    HeIon2 /= ( E0 / (PUnits::GV/PUnits::m));
  }
  
  TLine *lineHe = new TLine(0.0,-HeIon,xEextr[Npoints-1],-HeIon);
  lineHe->SetLineColor(kOrange+10);
  lineHe->SetLineStyle(2);
  lineHe->Draw();

  TLine *lineHe2 = new TLine(0.0,-HeIon2,xEextr[Npoints-1],-HeIon2);
  lineHe2->SetLineColor(kOrange+3);
  lineHe2->SetLineStyle(2);
  //lineHe2->Draw();
  
  
  for(Int_t i=0;i<Nfields;i++) {
    for(Int_t ic=0;ic<Ncross;ic++) {

      if(i>0) gEextr[i][ic]->SetLineColor(kGray+2);
      
      if(i==0 && ic==0) continue;
      if(i==1 && ic==1) continue;
     
      if( !gEextr[i][ic] ) continue;
      
      // Remove crazy points
      Double_t *yEextr  = gEextr[i][ic]->GetY();
      for(Int_t k=0; k<gEextr[i][ic]->GetN(); k++) {
	if(fabs(yEextr[k])>1000) {
	  gEextr[i][ic]->RemovePoint(k);
	  k--;
	}
      }
      
      gEextr[i][ic]->Draw("L");
    }
    Leg->Draw();
  } 


  C1->cd(0);
  
  // Print to a file
  // Output file
  TString fOutName = Form("./%s/Plots/Field1DEvolutions/Field1DEvolutions-Landscape-%s",sim.Data(),sim.Data());
  PlasmaGlob::imgconv(C1,fOutName,opt);
  
  // ---------------------------------------------------------

  ifile->Close();
  cout << endl;
}
  
