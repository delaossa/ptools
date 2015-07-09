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

void PlotEvolutionsTimon(const TString &sim, Int_t Nosc = 1, const Int_t Nfields=1, const TString &options="png") { 
  
#ifdef __CINT__  
  gSystem->Load("libplasma.so");
#endif

  PlasmaGlob::Initialize();
  
  // Palettes!
  gROOT->Macro("PlasmaPalettes.C");

  TString opt = options;
  
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
  
  Float_t maxEcross = -999.;
  Float_t minEcross = 999.;
  Float_t maxEextr  = -999.;
  Float_t minEextr  = 999.;
  Float_t maxEdephas  = -999.;
  Float_t minEdephas  = 999.;
 
  // In every oscillation there are 2 crossings
  Int_t Ncross = 2 * Nosc;

  TGraph *gEcross[Nfields][30];
  TGraph *gEextr[Nfields][30];
  TGraph *gEdephas[Nfields][30];
  TH2F   *hEvsTime[Nfields];
  TGraph *gTRatio[Nfields];

  TString filename;
  filename = Form("./%s/Plots/Evolutions/Evolutions-%s.root",sim.Data(),sim.Data());
  
  TFile  *ifile = (TFile*) gROOT->GetListOfFiles()->FindObject(filename.Data());
  if (!ifile) ifile = new TFile(filename,"READ");

  for(Int_t i=0;i<Nfields;i++) {
  
    char hName[24];
    sprintf(hName,"hEvsTime_%i",i); 
    hEvsTime[i] = (TH2F*) ifile->Get(hName);
    if(!hEvsTime[i]) {
      cout<< "NO MAIN HISTO: cancelling..." <<endl;
      return;
    }


    Int_t counter = 0;
    for(Int_t ic=0;ic<Ncross;ic++) {
      char gName[24];
      
      sprintf(gName,"gEcross_%i_%i",i,ic); 
      gEcross[i][ic] = (TGraph*) ifile->Get(gName);
            
      sprintf(gName,"gEextr_%i_%i",i,ic); 
      gEextr[i][ic] = (TGraph*) ifile->Get(gName);
      
      if( !gEcross[i][ic] || !gEextr[i][ic] ) continue;
      
      // Calculate the max and min of every set of graphs:
      
      Int_t Npoints = gEcross[i][ic]->GetN();
      Double_t *yEcross = gEcross[i][ic]->GetY();
      Double_t *yEextr  = gEextr[i][ic]->GetY();
      Double_t *xEextr  = gEextr[i][ic]->GetX();

      Double_t *yEdephas = new Double_t[Npoints];
      for(Int_t j=0;j<Npoints;j++) {
	yEdephas[j] = yEcross[j] - yEcross[0];
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
      if(ic%2-1) gEcross[i][ic]->SetLineStyle(2);
      else gEcross[i][ic]->SetLineStyle(1);
      gEcross[i][ic]->SetLineColor(kGray+2);
      gEcross[i][ic]->SetMarkerColor(kGray+2);
      gEcross[i][ic]->SetLineWidth(1);
      gEcross[i][ic]->SetMarkerSize(0.2);
      
      // only takes into account the minimums of the accelerating field:
      if(ic%2-1) continue;
      counter++;

      for(Int_t j=0;j<Npoints;j++) {
	if(yEcross[j]>maxEcross)
	  maxEcross = yEcross[j];
	if(yEcross[j]<minEcross)
	  minEcross = yEcross[j];
	
	if(yEextr[j]>maxEextr)
	  maxEextr = yEextr[j];
	if(yEextr[j]<minEextr)
	  minEextr = yEextr[j];
	
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
      gEcross[i][ic]->SetLineColor(colors[index]);
      gEextr[i][ic]->SetLineColor(colors[index]);
      gEextr[i][ic]->SetLineWidth(1);
      gEdephas[i][ic]->SetLineColor(colors[index]);
      gEdephas[i][ic]->SetLineWidth(1);
    }
  }

  // palettes for drawing
  PPalette * rbowPalette = (PPalette*) gROOT->FindObject("rbow");
  PPalette * rbowgrayPalette = (PPalette*) gROOT->FindObject("rbowgray");
  PPalette * rbowwhitePalette = (PPalette*) gROOT->FindObject("rbowwhite");
  PPalette * electronPalette = (PPalette*) gROOT->FindObject("electron0");
  PPalette * barsaPalette = (PPalette*) gROOT->FindObject("barsa");
  
  // Canvas setup
  TCanvas *C = new TCanvas("C","Evolution of Electric fields",sizex,sizey);
  
  C->cd(0);
  gPad->SetFrameLineWidth(2);  
  
  PlasmaGlob::SetH1LabelSize(hEvsTime[0]);
 
  // Change the range of z axis for the fields to be symmetric.
  Float_t Emax = hEvsTime[0]->GetMaximum();
  Float_t Emin = hEvsTime[0]->GetMinimum();
  if(Emax > TMath::Abs(Emin))
    Emin = -Emax;
  else
    Emax = -Emin;
  hEvsTime[0]->GetZaxis()->SetRangeUser(Emin,Emax); 
 
  rbowwhitePalette->cd();
  
  hEvsTime[0]->Draw("colz");

  gPad->Update();
  TPaletteAxis *palette = (TPaletteAxis*)hEvsTime[0]->GetListOfFunctions()->FindObject("palette");
  palette->SetX1NDC(1-margins[2]+0.02);
  palette->SetX2NDC(1-margins[2]+0.05);
  palette->SetY1NDC(margins[1]);
  palette->SetY2NDC(1-margins[3]);
  C->Modified();
  C->Update();

  for(Int_t i=0;i<Nfields;i++) {
    // if(i>0) continue;
    for(Int_t ic=0;ic<Ncross;ic++) {
      if( gEcross[i][ic] )
  	gEcross[i][ic]->Draw("C");
    }
  }
  
  gPad->RedrawAxis("G");

  // Print to a file
  // Output file
  TString fOutName = Form("./%s/Plots/Evolutions/Evolutions-A-%s",sim.Data(),sim.Data());
  PlasmaGlob::imgconv(C,fOutName,opt);
  // ---------------------------------------------------------

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
  Float_t E0   = minEextr-margin;
  Float_t E1   = maxEextr+margin;
  TH1F *hFrame2 = new TH1F("hFrame2","",100,hEvsTime[0]->GetXaxis()->GetXmin(),hEvsTime[0]->GetXaxis()->GetXmax());
  
  hFrame2->GetYaxis()->SetRangeUser(E0,E1);
  hFrame2->GetXaxis()->SetTitle("k_{p}z");
  hFrame2->GetYaxis()->SetTitle("E_{z}/E_{0}");
  PlasmaGlob::SetH1LabelSize(hFrame2);
 
  hFrame2->Draw("axis");

  gPad->Update();
  
  for(Int_t i=0;i<Nfields;i++) {
    for(Int_t ic=0;ic<Ncross;ic++) {
      
      // Only minimums:
      if(ic%2-1) continue;
      if( !gEextr[i][ic] ) continue;
      
      // Remove crazy points
      Double_t *yEextr  = gEextr[i][ic]->GetY();
      for(Int_t k=0; k<gEextr[i][ic]->GetN(); k++) {
	if(fabs(yEextr[k])>1000) {
	  gEextr[i][ic]->RemovePoint(k);
	  k--;
	}
      }
      
      gEextr[i][ic]->Draw("C");
    }
    // Leg->Draw();
  }

  C1->cd(0);
  
  // Print to a file
  // Output file
  fOutName = Form("./%s/Plots/Evolutions/Evolutions-B-%s",sim.Data(),sim.Data());
  PlasmaGlob::imgconv(C1,fOutName,opt);
  
  
  // ----------------------------------------------------------- Bottom plot -----
  
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
  hFrame3->GetXaxis()->SetTitle("k_{p}z");
  hFrame3->GetYaxis()->SetTitle("k_{p}#Delta#zeta");
  PlasmaGlob::SetH1LabelSize(hFrame3);
  hFrame3->Draw("axis");
  
  gPad->Update();
     
  for(Int_t i=0;i<Nfields;i++) {
    // if(i>0) continue;
    Int_t counter = 0;
    for(Int_t ic=0;ic<Ncross;ic++) {
      
      // Only minimums:
      if(ic%2-1) continue;
      if( !gEdephas[i][ic] ) continue;
	
      gEdephas[i][ic]->Draw("C");
    }
    
    // Leg->Draw();
  }
  
  gPad->RedrawAxis("G");
  
  C2->cd(0);
  
  // Print to a file
  // Output file
  fOutName = Form("./%s/Plots/Evolutions/Evolutions-C-%s",sim.Data(),sim.Data());
  PlasmaGlob::imgconv(C2,fOutName,opt);
  // ---------------------------------------------------------


  TH2F *hGvsTime = NULL;
  char hName[24];
  sprintf(hName,"hGvsTime"); 
  hGvsTime = (TH2F*) ifile->Get(hName);
  
  if(hGvsTime) {
    // Canvas setup
    TCanvas *C3 = new TCanvas("C3","Evolution of the dephasing",sizex,sizey);
    C3->cd(0);
    gPad->SetFrameLineWidth(2);  

    Float_t minGamma =  43.07 - 0.5;
    Float_t maxGamma =  43.07 + 0.5;
  
    hGvsTime->GetZaxis()->SetRangeUser(minGamma,maxGamma); 
    PlasmaGlob::SetH1LabelSize(hGvsTime);
  
    // barsaPalette->cd();
    rbowwhitePalette->cd();
    
    hGvsTime->Draw("colz");

    gPad->Update();
  
    palette = (TPaletteAxis*)hGvsTime->GetListOfFunctions()->FindObject("palette");
    palette->SetX1NDC(1-margins[2]+0.02);
    palette->SetX2NDC(1-margins[2]+0.05);
    palette->SetY1NDC(margins[1]);
    palette->SetY2NDC(1-margins[3]);
    C3->Modified();  
    gPad->RedrawAxis("G");
  
    C3->cd(0);
  
    // Print to a file
    // Output file
    fOutName = Form("./%s/Plots/Evolutions/Evolutions-D-%s",sim.Data(),sim.Data());
    PlasmaGlob::imgconv(C3,fOutName,opt);
    // ---------------------------------------------------------
  }

  // Canvas setup
  TCanvas *C4 = new TCanvas("C4","Evolution of the transverse RMS",sizex,sizey);
  C4->cd(0);
  gPad->SetFrameLineWidth(2);  

  TH2F *hRmsvsTime = NULL;
  sprintf(hName,"hRmsvsTime"); 
  hRmsvsTime = (TH2F*) ifile->Get(hName);
  
  Float_t Rmsmax = hRmsvsTime->GetMaximum();
  
  electronPalette->cd();

  hRmsvsTime->GetZaxis()->SetRangeUser(0,Rmsmax);   
  PlasmaGlob::SetH1LabelSize(hRmsvsTime);
  hRmsvsTime->Draw("colz");

  gPad->Update();
  
  palette = (TPaletteAxis*)hRmsvsTime->GetListOfFunctions()->FindObject("palette");
  palette->SetX1NDC(1-margins[2]+0.02);
  palette->SetX2NDC(1-margins[2]+0.05);
  palette->SetY1NDC(margins[1]);
  palette->SetY2NDC(1-margins[3]);
  C4->Modified();  
  gPad->RedrawAxis("G");
  
  C4->cd(0);
  
  // Print to a file
  // Output file
  fOutName = Form("./%s/Plots/Evolutions/Evolutions-E-%s",sim.Data(),sim.Data());
  PlasmaGlob::imgconv(C4,fOutName,opt);
  // ---------------------------------------------------------

  // Canvas setup
  TCanvas *C5 = new TCanvas("C5","Evolution of the on-axis density",sizex,sizey);
  C5->cd(0);
  gPad->SetFrameLineWidth(2);  

  TH2F *hDen1DvsTime = NULL;
  sprintf(hName,"hDen1DvsTime"); 
  hDen1DvsTime = (TH2F*) ifile->Get(hName);
  
  Float_t Denmax = hDen1DvsTime->GetMaximum();
  
  electronPalette->cd();

  hDen1DvsTime->GetZaxis()->SetRangeUser(0,Denmax);   
  PlasmaGlob::SetH1LabelSize(hDen1DvsTime); 
  hDen1DvsTime->Draw("colz");

  gPad->Update();
  
  palette = (TPaletteAxis*)hDen1DvsTime->GetListOfFunctions()->FindObject("palette");
  palette->SetX1NDC(1-margins[2]+0.02);
  palette->SetX2NDC(1-margins[2]+0.05);
  palette->SetY1NDC(margins[1]);
  palette->SetY2NDC(1-margins[3]);
  C5->Modified();  
  gPad->RedrawAxis("G");
  
  C5->cd(0);
  
  // Print to a file
  // Output file
  fOutName = Form("./%s/Plots/Evolutions/Evolutions-F-%s",sim.Data(),sim.Data());
  PlasmaGlob::imgconv(C5,fOutName,opt);
  // ---------------------------------------------------------

  
  ifile->Close();
  cout << endl;
}
  
