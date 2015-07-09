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

void PlotField1DEvolution(const TString &sim, Int_t Nosc = 1, const Int_t Nfields=1, const TString &options="png") { 
  
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
  Float_t margins[4] = {0.15,0.15,0.20,0.10};
  gStyle->SetPadLeftMargin(margins[0]);  // Margin left axis  
  gStyle->SetPadRightMargin(margins[2]);
  gStyle->SetPadTopMargin(margins[3]);  // Margin left axis  
  gStyle->SetPadBottomMargin(margins[1]);

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);


  if(opt.Contains("grid")) {
    gStyle->SetPadGridX(1);
    gStyle->SetPadGridY(1);
  }

  // Resolution:
  Int_t sizex = 800;
  Int_t sizey = 600;
  if(opt.Contains("hres")) {
    Int_t sizex = 1024;
    Int_t sizey = 768;    
  }
  
  Float_t maxzEcross = -999.;
  Float_t minzEcross = 999.;
  Float_t maxEextr  = -999.;
  Float_t minEextr  = 999.;
  Float_t maxzEextr  = -999.;
  Float_t minzEextr  = 999.;
  Float_t maxEdephas  = -999.;
  Float_t minEdephas  = 999.;
 
  // In every oscillation there are 2 crossings
  Int_t Ncross = 2 * Nosc;

  TGraph *gzEcross[Nfields][30];
  TGraph *gEextr[Nfields][30];
  TGraph *gzEextr[Nfields][30];
  TGraph *gEdephas[Nfields][30];
  TH2F   *hEvsTime[Nfields];
  TGraph *gTRatio[Nfields];

  TString filename;
  filename = Form("./%s/Plots/Field1DEvolutions/Field1DEvolutions-%s.root",sim.Data(),sim.Data());
  
  TFile  *ifile = (TFile*) gROOT->GetListOfFiles()->FindObject(filename.Data());
  if (!ifile) ifile = new TFile(filename,"READ");

  for(Int_t i=0;i<Nfields;i++) {
  
    char hName[24];
    sprintf(hName,"hEvsTime_%i",i); 
    hEvsTime[i] = (TH2F*) ifile->Get(hName);
    if(!hEvsTime[i]) {
      cout<< "NO MAIN HISTO: canceling..." <<endl;
      return;
    }
    
    if(opt.Contains("units")) {
      Int_t NbinsX = hEvsTime[i]->GetNbinsX();
      Float_t xMin = skindepth * hEvsTime[i]->GetXaxis()->GetXmin() / PUnits::um;
      Float_t xMax = skindepth * hEvsTime[i]->GetXaxis()->GetXmax() / PUnits::um;
      Int_t NbinsY = hEvsTime[i]->GetNbinsY();
      Float_t ymin = skindepth * hEvsTime[i]->GetYaxis()->GetXmin() / PUnits::um;
      Float_t ymax = skindepth * hEvsTime[i]->GetYaxis()->GetXmax() / PUnits::um;
      hEvsTime[i]->SetBins(NbinsX,xMin,xMax,NbinsY,ymin,ymax);
      for(Int_t j=0;j<=hEvsTime[i]->GetNbinsX();j++) {
      	for(Int_t k=0;k<=hEvsTime[i]->GetNbinsY();k++) {
      	  hEvsTime[i]->SetBinContent(j,k, hEvsTime[i]->GetBinContent(j,k) * ( E0 / (PUnits::GV/PUnits::m) ) );
      	}
      }

      if(opt.Contains("units")) {
	if(i==0) {
	  hEvsTime[i]->GetXaxis()->SetTitle("z [#mum]");
	  hEvsTime[i]->GetYaxis()->SetTitle("E_{z} [GV/m]");   
	} else if (i==1) {
	  hEvsTime[i]->GetXaxis()->SetTitle("z [#mum]");
	  hEvsTime[i]->GetYaxis()->SetTitle("E_{y} [GV/m]");   	  
	}
      }
      
    }
    
    Int_t counter = 0;
    for(Int_t ic=0;ic<Ncross;ic++) {
      char gName[24];
      
      sprintf(gName,"gzEcross_%i_%i",i,ic); 
      gzEcross[i][ic] = (TGraph*) ifile->Get(gName);
            
      sprintf(gName,"gEextr_%i_%i",i,ic); 
      gEextr[i][ic] = (TGraph*) ifile->Get(gName);

      sprintf(gName,"gzEextr_%i_%i",i,ic); 
      gzEextr[i][ic] = (TGraph*) ifile->Get(gName);
      
      if( !gzEcross[i][ic] || !gEextr[i][ic] || !gzEextr[i][ic] ) continue;
      
      // Calculate the max and min of every set of graphs:
      
      Int_t Npoints = gzEcross[i][ic]->GetN();
      Double_t *yzEcross = gzEcross[i][ic]->GetY();
      Double_t *yEextr  = gEextr[i][ic]->GetY();
      Double_t *yzEextr  = gzEextr[i][ic]->GetY();
     
      Double_t *xzEcross = gzEcross[i][ic]->GetX();
      Double_t *xEextr   = gEextr[i][ic]->GetX();
      Double_t *xzEextr  = gzEextr[i][ic]->GetX();

      if(opt.Contains("units")) {
	for(Int_t ip=0;ip<Npoints;ip++) {
	  yzEcross[ip] *= skindepth / PUnits::um;
	  yzEextr[ip]  *= skindepth / PUnits::um;

	  yEextr[ip]  *= E0 / (PUnits::GV/PUnits::m);
	  
	  xzEcross[ip]   *= skindepth / PUnits::um;
	  xEextr[ip]     *= skindepth / PUnits::um;
	  xzEextr[ip]    *= skindepth / PUnits::um;
	  
	}
	
      }

      Double_t *yEdephas = new Double_t[Npoints];
      for(Int_t j=0;j<Npoints;j++) {
	yEdephas[j] = yzEcross[j] - yzEcross[0];
      }
      gEdephas[i][ic] = new TGraph(Npoints,xEextr,yEdephas);
      sprintf(gName,"gEdephas_%i_%i",i,ic); 
      gEdephas[i][ic]->SetName(gName);
           
      // Transformer ratio vs time:
      // Take the ratio of the minimum over the maximum of the first oscillation.
      Double_t *TR = new Double_t[Npoints];
      if(ic==1) {
	Double_t *yExtrPrev = gEextr[i][ic-1]->GetY();
	for(Int_t j=0;j<Npoints;j++) {
	  TR[j] = TMath::Abs(yEextr[j]/yExtrPrev[j]);	
	}
	
	gTRatio[i] = new TGraph(Npoints,xEextr,TR);
	gTRatio[i]->SetName("TransRatio");
      }
      
      // Graph's attributes
      gzEcross[i][ic]->SetLineStyle(2);
      gzEcross[i][ic]->SetLineColor(kGray+2);
      gzEcross[i][ic]->SetMarkerColor(kGray+2);
      gzEcross[i][ic]->SetLineWidth(1);
      gzEcross[i][ic]->SetMarkerSize(0.2);
      
      // only takes into account the minimums of the accelerating field:
      // if(ic%2-1) continue;
      counter++;

      for(Int_t j=0;j<Npoints;j++) {
	if(yzEcross[j]>maxzEcross)
	  maxzEcross = yzEcross[j];
	if(yzEcross[j]<minzEcross)
	  minzEcross = yzEcross[j];
	
	if(yEextr[j]>maxEextr)
	  maxEextr = yEextr[j];
	if(yEextr[j]<minEextr)
	  minEextr = yEextr[j];

	if(yzEextr[j]>maxzEextr)
	  maxzEextr = yzEextr[j];
	if(yzEextr[j]<minzEextr)
	  minzEextr = yzEextr[j];
	
	if(i>0) continue;
	if(yEdephas[j]>maxEdephas)
	  maxEdephas = yEdephas[j];
	if(yEdephas[j]<minEdephas)
	  minEdephas = yEdephas[j];
	
      }
      
    }

  }

  // Set the color of the different evolutions according to a palette
  // UInt_t np = 50;
  // PPalette * colorPalette = (PPalette*) gROOT->FindObject("colorPalette");
  // if(!colorPalette) {
  //   const UInt_t Number = 3;
  //   Double_t Red[Number] = { 1.00, 0.00, 0.00};
  //   Double_t Green[Number]  = { 0.00, 1.00, 0.00};
  //   Double_t Blue[Number]   = { 1.00, 0.00, 1.00};
  //   Double_t Length[Number] = { 0.00, 0.50, 1.00 };
  //   colorPalette = new PPalette("colorPalette");
  //   colorPalette->CreateGradientColorTable(Number,Length,Red,Green,Blue,np);
  // }
  
  // for(Int_t i=0;i<Nfields;i++) { 
  //   for(Int_t ic=0;ic<Ncross;ic++) {
  //     Float_t step = (np/Nosc);
  //     Int_t icolor = TMath::Nint( ((ic+1)/2) * step - 1 );
  //     gEextr[i][ic]->SetLineColor(colorPalette->GetColor(icolor));
  //     gEextr[i][ic]->SetLineWidth(2);
  //     gEdephas[i][ic]->SetLineColor(colorPalette->GetColor(icolor));
  //     gEdephas[i][ic]->SetLineWidth(2);
  //   }
  // }
  // -------------------------------------------------------

  // Manual coloring:
  const Int_t NCOLORS = 5;
  Int_t colors[NCOLORS] = {kMagenta+2,kRed,kBlue,kYellow+2,kCyan+2};
  for(Int_t i=0;i<Nfields;i++) { 
    for(Int_t ic=0;ic<Ncross;ic++) {
      Int_t index = ic/2;
      if(index>=NCOLORS) index = NCOLORS-1;
      gzEcross[i][ic]->SetLineColor(colors[index]);
      gzEcross[i][ic]->SetLineWidth(2);
      gEextr[i][ic]->SetLineColor(colors[index]);
      gEextr[i][ic]->SetLineWidth(2);
      gzEextr[i][ic]->SetLineColor(colors[index]);
      gzEextr[i][ic]->SetLineWidth(2);
      gzEextr[i][ic]->SetLineStyle(1);
      gEdephas[i][ic]->SetLineColor(colors[index]);
      gEdephas[i][ic]->SetLineWidth(2);
    }
  }

  // palettes for drawing
  PPalette * rbowPalette = (PPalette*) gROOT->FindObject("rbow");
  PPalette * rbowgrayPalette = (PPalette*) gROOT->FindObject("rbowgray");
  PPalette * rbow2whitePalette = (PPalette*) gROOT->FindObject("rbow2white");
  PPalette * electronPalette = (PPalette*) gROOT->FindObject("electron0");
  PPalette * barsaPalette = (PPalette*) gROOT->FindObject("barsa");
  
  // Canvas setup
  TCanvas *C = new TCanvas("C","Evolution of Electric fields",sizex,sizey);
  
  
  for(Int_t i=0;i<Nfields;i++) {

    C->cd(0);
    C->Clear();
    gPad->SetFrameLineWidth(2);  

    PlasmaGlob::SetH1LabelSize(hEvsTime[i]);
 
    // Change the range of z axis for the fields to be symmetric.
    Float_t Emax = hEvsTime[i]->GetMaximum();
    Float_t Emin = hEvsTime[i]->GetMinimum();
    if(Emax > TMath::Abs(Emin))
      Emin = -Emax;
    else
      Emax = -Emin;
    hEvsTime[i]->GetZaxis()->SetRangeUser(Emin,Emax); 
    
    rbow2whitePalette->cd();
    
    hEvsTime[i]->Draw("colz");
    
    gPad->Update();
    TPaletteAxis *palette = (TPaletteAxis*)hEvsTime[i]->GetListOfFunctions()->FindObject("palette");
    palette->SetX1NDC(1-margins[2]+0.02);
    palette->SetX2NDC(1-margins[2]+0.05);
    palette->SetY1NDC(margins[1]);
    palette->SetY2NDC(1-margins[3]);
    C->Modified();
    C->Update();
    
    for(Int_t ic=0;ic<Ncross;ic++) {
      if( gzEcross[i][ic] ) {
  	gzEcross[i][ic]->SetLineColor(kGray+2);
	gzEcross[i][ic]->Draw("L");
      }
      
      if( gzEextr[i][ic] ) {
  	gzEextr[i][ic]->Draw("L");
      }
    }
      
    gPad->RedrawAxis("G");
    
    // Print to a file
    // Output file
    TString fOutName = Form("./%s/Plots/Field1DEvolutions/Field1DEvolutions-E%1i-%s",sim.Data(),i,sim.Data());
    PlasmaGlob::imgconv(C,fOutName,opt);
    // ---------------------------------------------------------

  }

  // Canvas setup
  TCanvas *C1 = new TCanvas("C1","Evolution of Electric fields",sizex,sizey);

  TLegend *Leg = new TLegend(0.20,0.20,0.35,0.40);
  PlasmaGlob::SetPaveStyle(Leg);
  Leg->SetTextAlign(22);
  Leg->SetTextColor(kGray+3);
  Leg->SetLineColor(kGray+1);
  Leg->SetBorderSize(0);
  Leg->SetFillColor(0);
  Leg->SetFillStyle(1001);
  for(Int_t ic=0;ic<Ncross;ic++) {
    // Only minimums:
    if(ic%2-1) continue;
    Leg->AddEntry(gEextr[0][ic],Form("Minimum %2i",ic/2+1),"L");
  }
  Leg->SetTextColor(kGray+3);
  
  
  C1->cd(0);
  gPad->SetFrameLineWidth(2);  
  gPad->SetRightMargin(0.10);
  
  Float_t margin = (maxEextr - minEextr)/10;
  Float_t Eextr0   = minEextr-margin;
  Float_t Eextr1   = maxEextr+margin;
  TH1F *hFrame2 = new TH1F("hFrame2","",100,hEvsTime[0]->GetXaxis()->GetXmin(),hEvsTime[0]->GetXaxis()->GetXmax());
  
  hFrame2->GetYaxis()->SetRangeUser(Eextr0,Eextr1);
 
  if(opt.Contains("units")) {
    hFrame2->GetXaxis()->SetTitle("z [#mum]");
    hFrame2->GetYaxis()->SetTitle("E_{z} [GV/m]");   
  } else {
    hFrame2->GetXaxis()->SetTitle("k_{p}z");
    hFrame2->GetYaxis()->SetTitle("E_{z}/E_{0}");
  }
  PlasmaGlob::SetH1LabelSize(hFrame2);
 
  hFrame2->Draw("axis");

  gPad->Update();
  
  for(Int_t i=0;i<Nfields;i++) {
    for(Int_t ic=0;ic<Ncross;ic++) {

      if(i>0) gEextr[i][ic]->SetLineColor(kGray+2);
      
      // Only minimums:
      //  if(ic%2-1) continue;
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
    // Leg->Draw();
  }

  C1->cd(0);
  
  // Print to a file
  // Output file
  TString fOutName = Form("./%s/Plots/Field1DEvolutions/Field1DEvolutions-Emax-%s",sim.Data(),sim.Data());
  PlasmaGlob::imgconv(C1,fOutName,opt);
  
  
  // Canvas setup
  TCanvas *C2 = new TCanvas("C2","Evolution of the dephasing",sizex,sizey);
  C2->cd(0);
  gPad->SetFrameLineWidth(2);  
  gPad->SetRightMargin(0.10);
 
  margin = (maxEdephas - minEdephas)/10;
  Float_t D0   = minEdephas-margin;
  Float_t D1   = maxEdephas+margin;
  TH1F *hFrame3 = new TH1F("hFrame3","",100,hEvsTime[0]->GetXaxis()->GetXmin(),hEvsTime[0]->GetXaxis()->GetXmax());
  
  hFrame3->GetYaxis()->SetRangeUser(D0,D1);

  if(opt.Contains("units")) {
    hFrame3->GetXaxis()->SetTitle("z [#mum]");
    hFrame3->GetYaxis()->SetTitle("#Delta#zeta [#mum]");   
  } else {
    hFrame3->GetXaxis()->SetTitle("k_{p}z");
    hFrame3->GetYaxis()->SetTitle("#Delta#zeta");
  }
  PlasmaGlob::SetH1LabelSize(hFrame3);
  hFrame3->Draw("axis");
  
  gPad->Update();
     
  for(Int_t i=0;i<Nfields;i++) {
    
    // Only E_z
    if(i>0) continue;
   
    Int_t counter = 0;
    for(Int_t ic=0;ic<Ncross;ic++) {
      
      // Only minimums:
      if(ic%2-1) continue;
      if( !gEdephas[i][ic] ) continue;
      
      gEdephas[i][ic]->Draw("L");
    }
    
    // Leg->Draw();
  }
  
  gPad->RedrawAxis("G");
  
  C2->cd(0);
  
  // Print to a file
  // Output file
  fOutName = Form("./%s/Plots/Field1DEvolutions/Field1DEvolutions-dephasing-%s",sim.Data(),sim.Data());
  PlasmaGlob::imgconv(C2,fOutName,opt);
  // ---------------------------------------------------------

  ifile->Close();
  cout << endl;
}
  
