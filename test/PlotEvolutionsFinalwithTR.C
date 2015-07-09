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

void PlotEvolutionsFinalwithTR(const TString &sim, Int_t Nosc = 2, const Int_t Nfields=2, const TString &options="png") { 
  
#ifdef __CINT__  
  gSystem->Load("libplasma.so");
#endif

  PlasmaGlob::Initialize();
  
  // Palettes!
  gROOT->Macro("PlasmaPalettes.C");

  TString opt = options;
  
  // More makeup
  gStyle->SetPadLeftMargin(0.10);  // Margin left axis  
  gStyle->SetPadRightMargin(0.12);
  gStyle->SetTitleSize(0.06, "xy");
  gStyle->SetTitleSize(0.05, "z");
  gStyle->SetLabelSize(0.05, "xy");
  gStyle->SetLabelSize(0.04, "z");
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetTitleOffset(0.8,"y");
  gStyle->SetTitleOffset(0.5,"z");

  if(opt.Contains("grid")) {
    gStyle->SetPadGridX(1);
    gStyle->SetPadGridY(1);
  }
  
  Float_t maxEcross = -999.;
  Float_t minEcross = 999.;
  Float_t maxEextr  = -999.;
  Float_t minEextr  = 999.;

  // In every oscillation there are 2 crossings
  Int_t Ncross = 2 * Nosc;
 
  TGraph *gEcross[Nfields][20];
  TGraph *gEextr[Nfields][20];
  TH2F   *hEvsTime[Nfields];
  TGraph **gE1D[Nfields];
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
      cout << "NO MAIN HISTO: " << hName << " cancelling..." <<endl;
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

      // Transformer ratio vs time:
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
	
      }
      
    }

    Int_t NT = hEvsTime[i]->GetNbinsX();
    gE1D[i] = new TGraph*[NT];
    for(Int_t it = 0; it<NT; it++) {

      Int_t NZ = hEvsTime[i]->GetNbinsY();
      Float_t *x = new Float_t[NZ];
      Float_t *y = new Float_t[NZ];
      for(Int_t iz = 0; iz<NZ; iz++) {
	x[iz] = hEvsTime[i]->GetYaxis()->GetBinCenter(iz);
	y[iz] = hEvsTime[i]->GetBinContent(it,iz);
      }
      
      gE1D[i][it] = new TGraph(NZ,x,y);
    }
  }


  // Manual coloring:
  const Int_t NCOLORS = 5;
  Int_t colors[NCOLORS] = {kGray+3,kGray+2,kGray+1,kCyan+2};
  for(Int_t i=0;i<Nfields;i++) { 
    for(Int_t ic=0;ic<Ncross;ic++) {
      Int_t index = ic/2;
      if(index>=NCOLORS) index = NCOLORS-1;
      gEcross[i][ic]->SetLineColor(colors[index]);
      gEextr[i][ic]->SetLineColor(colors[index]);
    }
  }
  
  // Output file
  TString fOutName = Form("./%s/Plots/Evolutions/EvolutionsFinal-%s",sim.Data(),sim.Data());
  
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
  Leg->AddEntry(gEcross[0][0],Form("Accelerating field"),"L");
  Leg->AddEntry(gEcross[1][0],Form("Focusing field"),"L"); 
  Leg->SetTextColor(kGray+3);

  
  // Actual Plotting!
  // ------------------------------------------------------------
  
  // More makeup

  C->cd(1);
  gPad->SetFrameLineWidth(2);  
  
  // Float_t margin = (maxEcross - minEcross)/10;
  // TH1F *hFrame = new TH1F("hFrame","",10,500,3000);
  // hFrame->GetYaxis()->SetRangeUser(minEcross-margin,maxEcross+margin);
  // hFrame->GetXaxis()->CenterTitle();
  // hFrame->GetYaxis()->CenterTitle();
  // hFrame->GetXaxis()->SetTitle("time [#omega_{p}^{-1}]");
  // hFrame->GetYaxis()->SetTitle("#zeta [c/#omega_{p}]");
  
  // Change the range of z axis for the fields to be symmetric.
  Float_t Emax = hEvsTime[0]->GetMaximum();
  Float_t Emin = hEvsTime[0]->GetMinimum();
  if(Emax > TMath::Abs(Emin))
    Emin = -Emax;
  else
    Emax = -Emin;
  hEvsTime[0]->GetZaxis()->SetRangeUser(Emin,Emax); 
  PPalette * rbowPalette = (PPalette*) gROOT->FindObject("rbow2");
  rbowPalette->cd();

  hEvsTime[0]->GetZaxis()->SetTitleOffset(0.6);  
  hEvsTime[0]->Draw("colz");

  gPad->Update();
  TPaletteAxis *palette = (TPaletteAxis*)hEvsTime[0]->GetListOfFunctions()->FindObject("palette");
  palette->SetX1NDC(0.9);
  palette->SetX2NDC(0.92);
  palette->SetY1NDC(0.2);
  palette->SetY2NDC(0.9);
  C->Modified();
  C->Update();

  // De-phasing lines:
  // TLine *lpi0 = new TLine(hFrame->GetXaxis()->GetXmin(),0.,hFrame->GetXaxis()->GetXmax(),0.);
  // lpi0->SetLineColor(kGray+1);
  // lpi0->SetLineStyle(2);
  // lpi0->Draw();
  // TLine *lpi = new TLine(hFrame->GetXaxis()->GetXmin(),-TMath::Pi(),hFrame->GetXaxis()->GetXmax(),-TMath::Pi());
  // lpi->SetLineColor(kGray+1);
  // lpi->SetLineStyle(2);
  // lpi->Draw();
  // TLine *lpi2 = new TLine(hFrame->GetXaxis()->GetXmin(),-TMath::TwoPi(),hFrame->GetXaxis()->GetXmax(),-TMath::TwoPi());
  // lpi2->SetLineColor(kGray+1);
  // lpi2->SetLineStyle(2);
  // lpi2->Draw();
  // for(Int_t i=0;i<Nsim;i++) {
  
  for(Int_t i=0;i<Nfields;i++) {
    // if(i>0) continue;
    for(Int_t ic=0;ic<Ncross;ic++) {
      if( gEcross[i][ic] )
	gEcross[i][ic]->Draw("L");
    }
  }
  
  gPad->RedrawAxis("G");

  // ----------------------------------------------------------- Bottom plot -----
  C->cd(2);
  gPad->SetFrameLineWidth(2);  
  
  Float_t margin = (maxEextr - minEextr)/10;
  Float_t E0   = minEextr-margin;
  Float_t E1   = maxEextr+margin;
  TH1F *hFrame2 = new TH1F("hFrame2","",100,hEvsTime[0]->GetXaxis()->GetXmin(),hEvsTime[0]->GetXaxis()->GetXmax());
  
  hFrame2->GetYaxis()->SetRangeUser(E0,E1);
  hFrame2->GetXaxis()->CenterTitle();
  hFrame2->GetYaxis()->CenterTitle();
  hFrame2->GetXaxis()->SetTitle("#omega_{p}t");
  hFrame2->GetYaxis()->SetTitle("E_{z}/E_{0}");
  
  hFrame2->Draw("axis");

  gPad->Update();
  
  // Texts for the Field evolution lines
  TPaveText *textF[20];
   
  for(Int_t i=0;i<Nfields;i++) {
    // if(i>0) continue;
    Int_t counter = 0;
    for(Int_t ic=Ncross-1;ic>=0;ic--) {
      
      if(ic%2-1) continue;
      if( !gEextr[i][ic] ) continue;
      
      gEextr[i][ic]->Draw("L");
      counter++;
      
      Double_t slope = 0.8/(E1-E0);
      Double_t t,E;	
      gEextr[i][ic]->GetPoint(gEextr[i][ic]->GetN()-1,t,E);
      textF[ic] = new TPaveText(0.7,slope*(E-E0)+0.05,0.85,slope*(E-E0)+0.15,"NDC");
      //cout << "E value: " << E << "  y value: " << slope*(E-E0) + 0.15 << endl;
      PlasmaGlob::SetPaveTextStyle(textF[ic],32);
      textF[ic]->AddText(Form("Minimum number %2i",counter));
      textF[ic]->SetTextColor(PlasmaGlob::fieldLine);
      //  textF[ic]->Draw();
      
    }
    
    Double_t rightmin = 0.0;
    Double_t rightmax = 12.0;
    Double_t slope = (gPad->GetUymax() - gPad->GetUymin())/(rightmax-rightmin); 
    
    Double_t *xTR = gTRatio[i]->GetX();
    Double_t *yTR = gTRatio[i]->GetY();
    for(Int_t j=0;j<gTRatio[i]->GetN();j++) {
      gTRatio[i]->SetPoint(j,xTR[j],yTR[j]*slope + gPad->GetUymin());
    }
    
    gTRatio[i]->SetLineWidth(2);
    gTRatio[i]->SetLineColor(kRed);
    gTRatio[i]->Draw("L");

    //draw an axis on the right side
    TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),gPad->GetUxmax(),
			      gPad->GetUymax(),rightmin,rightmax,505,"+L");
    
    axis->SetLineWidth(1);
    axis->SetLineColor(gTRatio[i]->GetLineColor());
    axis->SetLabelColor(gTRatio[i]->GetLineColor());
    axis->SetLabelSize(0.05);
    axis->SetTitleSize(0.06);
    axis->SetTitleOffset(0.7);
    axis->SetTitle("TR");
    axis->CenterTitle();
    axis->SetTitleColor(gTRatio[i]->GetLineColor());
    
    axis->Draw();
  }
  
  gPad->RedrawAxis("G");
  
  C->cd(0);
  // Print to a file
  PlasmaGlob::imgconv(C,fOutName,opt);
  // ---------------------------------------------------------
  
  ifile->Close();
  cout << endl;
}
  
