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
#include "PGlobals.hh"
#include "PPalette.hh"


Float_t ntop0 = 15.63;
Float_t sigmal0 = 68.49;

Double_t DenGauss(Double_t z, Double_t sigmal, Double_t ntop) {
  
  Double_t nnorm =  1.0 + (ntop - 1.0) * TMath::Exp(-(z*z)/(2*sigmal*sigmal)) ;
  return nnorm;

}

Double_t DerDenGauss(Double_t z, Double_t sigmal, Double_t ntop) {
  
  Double_t dernnorm =  -(z/(sigmal*sigmal)) * (ntop - 1) * TMath::Exp(-(z*z)/(2*sigmal*sigmal));
  return dernnorm;
  
}

Double_t betaph(Double_t z, Double_t phase, Double_t sigmal, Double_t ntop) {

  Double_t beta = 1.0 / ( 1.0 + (phase/2) * TMath::Power(DenGauss(z,sigmal,ntop),-3/2.) * DerDenGauss(z,sigmal,ntop)  ) ;
  return beta;
    
}

Double_t betapsi(Double_t psi) {
  Double_t beta = - ( 2*psi + psi*psi ) / (2*psi + psi*psi + 2) ;
  return beta;
}


void PlotDDR(const TString &sim, const TString &options="png") { 
  
#ifdef __CINT__  
  gSystem->Load("libptools.so");
#endif

  PGlobals::Initialize();
  
  // Palettes!
  gROOT->Macro("PPalettes.C");

  TString opt = options;

  // Colors
  Int_t myBlue = TColor::GetColor((Float_t) 0.16, (Float_t) 0.83, (Float_t) 0.5);
  Int_t myNaranja = TColor::GetColor((Float_t) 0.992157, (Float_t) 0.411765, (Float_t) 0.027451);

  // More makeup            
  Float_t margins[4] = {0.15,0.15,0.20,0.10};
  gStyle->SetPadLeftMargin(margins[0]);  // Margin left axis  
  gStyle->SetPadRightMargin(margins[2]);
  gStyle->SetPadTopMargin(margins[3]);  // Margin left axis  
  gStyle->SetPadBottomMargin(margins[1]);

  gStyle->SetPadTickX(0);
  gStyle->SetPadTickY(0);

  gStyle->SetNumberContours(255);

  if(opt.Contains("grid")) {
    gStyle->SetPadGridX(1);
    gStyle->SetPadGridY(1);
  }

  // Load first simulation data (for instance)
  PData *pData = PData::Get(sim.Data());
  Double_t E0 = pData->GetPlasmaE0();
  
  Float_t maxEcross = -999.;
  Float_t minEcross = 999.;
  Float_t maxEextr  = -999.;
  Float_t minEextr  = 999.;
  Float_t maxEdephas  = -999.;
  Float_t minEdephas  = 999.;
  Float_t maxVextr = -999.;
  Float_t minVextr = 999.;

  const Int_t Nfields = 1;
  TH2F   *hEvsTime[Nfields];
  Int_t NCross[Nfields];
  TGraph **gEcross[Nfields]; 
  TGraph **gEextr[Nfields]; 
  TGraph **gEdephas[Nfields]; 
  TGraph *gTRatio;

  TH2F   *hVvsTime;
  TGraph **gVcross = NULL;
  TGraph **gVextr = NULL;
  TGraph **gVextr_alt = NULL;
  TGraph **gVextr_avg = NULL;
  TGraph *gVminvsDen = NULL;
  TGraph *gVmaxvsDen = NULL;
    
  char hName[24];
  char gName[24];
  
  TString filename;
  filename = Form("./%s/Plots/Evolutions/Evolutions-%s.root",sim.Data(),sim.Data());
  
  TFile  *ifile = (TFile*) gROOT->GetListOfFiles()->FindObject(filename.Data());
  if (!ifile) ifile = new TFile(filename,"READ");
  
  for(Int_t i=0;i<Nfields;i++) {
    
    sprintf(hName,"hEvsTime_%i",i); 
    hEvsTime[i] = (TH2F*) ifile->Get(hName);
    if(!hEvsTime[i]) continue;
    
    cout << Form("ANALYZING FIELD %i ...",i) << endl;

    Int_t NTBins = hEvsTime[i]->GetNbinsX();
    for(Int_t it=NTBins;it>0;it--) {
      
      // 1D field at certain timestep "it".
      TH1F *hE1D = (TH1F*) hEvsTime[i]->ProjectionY("_py",it,it);

      Float_t x1max = hE1D->GetXaxis()->GetXmax();

      Int_t MAXCROSS = 2;
      Float_t *Cross = new Float_t[MAXCROSS];
      Float_t *Extr = new Float_t[MAXCROSS];
      memset(Cross,0,sizeof(Float_t)*MAXCROSS);
      memset(Extr,0,sizeof(Float_t)*MAXCROSS);

      Int_t auxNcross = PGlobals::HCrossings(hE1D,Cross,Extr,MAXCROSS,0.,x1max/2);
      // cout << Form("  -> Number of crossings for histogram \"%s\" : %i ",hE1D->GetName(),auxNcross) << endl;
      // for(Int_t ic=0;ic<auxNcross;ic++) {
      // 	cout << Form(" %2i:  cross = %6.4f  extreme = %6.4f", ic, Cross[ic], Extr[ic]) << endl; 
      // }
      
      if(it==NTBins) {
	NCross[i] = auxNcross;
	
	gEcross[i] = new TGraph*[NCross[i]];
	gEextr[i] = new TGraph*[NCross[i]];

	for(Int_t ic = 0;ic<NCross[i];ic++) {
	  gEcross[i][ic] = new TGraph(NTBins);
	  sprintf(gName,"gEcross_%i_%i",i,ic); 
	  gEcross[i][ic]->SetName(gName);
	  
	  gEextr[i][ic] = new TGraph(NTBins);
	  sprintf(gName,"gEextr_%i_%i",i,ic); 
	  gEextr[i][ic]->SetName(gName);
	}

      }

      Float_t time = hEvsTime[i]->GetXaxis()->GetBinCenter(it);
      // cout << Form("Time step %i (%.2f): %i crossings",it,time,NCross[i]) << endl;
     
      for(Int_t ic=0;ic<NCross[i];ic++) {
	// cout << Form("  - Adding %i crossing: cross = %6.4f extreme = %6.4f",ic,Cross[ic],Extr[ic]) << endl;
	
	gEcross[i][ic]->SetPoint(it-1,time,Cross[ic]);
	gEextr[i][ic]->SetPoint(it-1,time,Extr[ic]);
      }
      
    }

    
    // Calculate the max and min of every crossing.
    // Also calculates dephasing:
    gEdephas[i] = new TGraph*[NCross[i]];
    for(Int_t ic = 0;ic<NCross[i];ic++) {

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
      
      for(Int_t j=0;j<Npoints;j++) {
	if(yEcross[j]>maxEcross)
	  maxEcross = yEcross[j];
	if(yEcross[j]<minEcross)
	  minEcross = yEcross[j];
	
	if(yEextr[j]>maxEextr)
	  maxEextr = yEextr[j];
	if(yEextr[j]<minEextr)
	  minEextr = yEextr[j];
      
	// Only takes into account the minimums of the accelerating field:
	if(ic%2 || i!=0) continue;
	if(yEdephas[j]>maxEdephas)
	  maxEdephas = yEdephas[j];
	if(yEdephas[j]<minEdephas)
	  minEdephas = yEdephas[j];
	
      }
    
    }
  }

  // Transformer ratio vs time:
  // Take the ratio of the minimum over the maximum of the first oscillation of the Ez.
  if(gEextr[0][1] && gEextr[0][0]) {
    Int_t Npoints = gEextr[0][1]->GetN();
    Double_t *TR = new Double_t[Npoints];
    Double_t *yExtrPrev = gEextr[0][0]->GetY();
    Double_t *yExtr = gEextr[0][1]->GetY();
    Double_t *xExtr = gEextr[0][1]->GetX();
    for(Int_t j=0;j<Npoints;j++) {
      TR[j] = TMath::Abs(yExtr[j]/yExtrPrev[j]);	
    }
    
    gTRatio = new TGraph(Npoints,xExtr,TR);
    sprintf(gName,"gTRatio"); 
    gTRatio->SetName(gName);
  }
  

  sprintf(hName,"hVvsTime"); 
  hVvsTime = (TH2F*) ifile->Get(hName);
  Int_t NVCross = 0;
  cout << Form("ANALYZING POTENTIAL") << endl;

  Int_t NTBins = hVvsTime->GetNbinsX();
  for(Int_t it=NTBins;it>0;it--) {
      
    // 1D field at certain timestep "it".
    TH1F *hV1D = (TH1F*) hVvsTime->ProjectionY("_py",it,it);
      
    Int_t MAXCROSS = 2;
    Float_t *Cross = new Float_t[MAXCROSS];
    Float_t *Extr = new Float_t[MAXCROSS];
    memset(Cross,0,sizeof(Float_t)*MAXCROSS);
    memset(Extr,0,sizeof(Float_t)*MAXCROSS);
    
    Int_t auxNcross = PGlobals::HCrossings(hV1D,Cross,Extr,MAXCROSS,0.,0.);
    // cout << Form("  -> Number of crossings for histogram \"%s\" : %i ",hV1D->GetName(),auxNcross) << endl;
    // for(Int_t ic=0;ic<auxNcross;ic++) {
    // 	cout << Form(" %2i:  cross = %6.4f  extreme = %6.4f", ic, Cross[ic], Extr[ic]) << endl; 
    // }
    
    if(it==NTBins) {
      NVCross = auxNcross;
	
      gVcross = new TGraph*[NVCross];
      gVextr = new TGraph*[NVCross];
      gVextr_alt = new TGraph*[NVCross];
      gVextr_avg = new TGraph*[NVCross];
      
      for(Int_t ic = 0;ic<NVCross;ic++) {
	gVcross[ic] = new TGraph(NTBins);
	sprintf(gName,"gVcross_%i",ic); 
	gVcross[ic]->SetName(gName);
	  
	gVextr[ic] = new TGraph(NTBins);
	sprintf(gName,"gVextr_%i",ic); 
	gVextr[ic]->SetName(gName);

	gVextr_alt[ic] = new TGraph(NTBins);
	sprintf(gName,"gVextr_alt_%i",ic); 
	gVextr_alt[ic]->SetName(gName);

	gVextr_avg[ic] = new TGraph(NTBins);
	sprintf(gName,"gVextr_avg_%i",ic); 
	gVextr_avg[ic]->SetName(gName);
      }

    }

    Float_t time = hVvsTime->GetXaxis()->GetBinCenter(it);
    // cout << Form("Time step %i (%.2f): %i crossings",it,time,NVCross) << endl;
     
    for(Int_t ic=0;ic<NVCross;ic++) {
      // cout << Form("  - Adding %i crossing: cross = %6.4f extreme = %6.4f",ic,Cross[ic],Extr[ic]) << endl;
	
      gVcross[ic]->SetPoint(it-1,time,Cross[ic]);
      gVextr[ic]->SetPoint(it-1,time,Extr[ic]);

      if(Extr[ic]>maxVextr)
	maxVextr = Extr[ic];
      if(Extr[ic]<minVextr)
	minVextr = Extr[ic];

      if(gEcross[0][ic]) {
	Double_t x,y;
	gEcross[0][ic]->GetPoint(it-1,x,y);      
	Double_t value = hV1D->GetBinContent(hV1D->FindBin(y));
	gVextr_alt[ic]->SetPoint(it-1,time,value);
	gVextr_avg[ic]->SetPoint(it-1,time,(value+Extr[ic])/2);
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
  // --------------------------------------------------------------------------

  // Manual coloring:
  const Int_t NCOLORS = 5;
  //  Int_t colors[NCOLORS] = {kMagenta+2,kRed,kBlue,kYellow+2,kCyan+2};
  Int_t colors[NCOLORS] = {kGray+3,kGray+2,kGray+1,kGray};
  for(Int_t i=0;i<Nfields;i++) { 
    for(Int_t ic=0;ic<NCross[i];ic++) {

      if( !gEcross[i][ic] || !gEextr[i][ic] ) continue;
            
      Int_t index = ic/2;
      if(index>=NCOLORS) index = NCOLORS-1;
      gEcross[i][ic]->SetLineColor(colors[index]);
      gEextr[i][ic]->SetLineColor(colors[index]);
      gEextr[i][ic]->SetLineWidth(1);
      gEdephas[i][ic]->SetLineColor(colors[index]);
      gEdephas[i][ic]->SetLineWidth(1);
    }
  }

  for(Int_t ic = 0;ic<NVCross;ic++) {
    // Graph's attributes
    Int_t index = ic/2;
    if(index>=NCOLORS) index = NCOLORS-1;
    gVcross[ic]->SetLineColor(colors[index]);
    gVcross[ic]->SetLineStyle(3);
    gVextr[ic]->SetLineStyle(1);
    gVextr[ic]->SetLineWidth(2);
    gVextr_alt[ic]->SetLineStyle(1);
    gVextr_alt[ic]->SetLineColor(2);
    gVextr_alt[ic]->SetLineWidth(2);

    gVextr_avg[ic]->SetLineStyle(1);
    gVextr_avg[ic]->SetLineWidth(3);

    if(ic==0)
      gVextr_avg[ic]->SetLineColor(myNaranja);
    else if(ic==1)
      gVextr_avg[ic]->SetLineColor(myBlue);
    
  }
  

  // palettes for drawing
  PPalette *fieldPalette = (PPalette*) gROOT->FindObject("field");
  if(!fieldPalette) {
    fieldPalette = new PPalette("field");
  }
  fieldPalette->SetPalette("rbow0");
    
  // Canvas setup
  Int_t sizex = 1024;
  Int_t sizey = 640;
  TCanvas *C = new TCanvas("C","Evolution of Ez (on-axis)",sizex,sizey);
  
  C->cd(0);
  gPad->SetFrameLineWidth(2);  
  gPad->SetTickx(1);
  gPad->SetTicky(1);
  // gPad->SetFillStyle(4000);
  // gPad->SetFrameFillStyle(4000);

  PGlobals::SetH1LabelSize(hEvsTime[0]);
 
  // Change the range of z axis for the fields to be symmetric.
  Float_t Emax = hEvsTime[0]->GetMaximum();
  Float_t Emin = hEvsTime[0]->GetMinimum();
  if(Emax > TMath::Abs(Emin))
    Emin = -Emax;
  else
    Emax = -Emin;
  hEvsTime[0]->GetZaxis()->SetRangeUser(Emin,Emax); 
 
  hEvsTime[0]->GetYaxis()->SetLabelFont(43);
  hEvsTime[0]->GetYaxis()->SetLabelSize(32);
  //  hEvsTime[0]->GetYaxis()->SetLabelOffset(0.8);
  hEvsTime[0]->GetYaxis()->SetTitleFont(43);
  hEvsTime[0]->GetYaxis()->SetTitleSize(36);
  hEvsTime[0]->GetYaxis()->SetTitleOffset(0.9);

  hEvsTime[0]->GetXaxis()->SetLabelFont(43);
  hEvsTime[0]->GetXaxis()->SetLabelSize(32);
  //  hEvsTime[0]->GetXaxis()->SetLabelOffset(0.8);
  hEvsTime[0]->GetXaxis()->SetTitleFont(43);
  hEvsTime[0]->GetXaxis()->SetTitleSize(36);
  hEvsTime[0]->GetXaxis()->SetTitleOffset(0.95);


  hEvsTime[0]->GetZaxis()->SetLabelFont(43);
  hEvsTime[0]->GetZaxis()->SetLabelSize(32);
  // hEvsTime[0]->GetZaxis()->SetLabelOffset(0.8);
  hEvsTime[0]->GetZaxis()->SetTitleFont(43);
  hEvsTime[0]->GetZaxis()->SetTitleSize(36);
  hEvsTime[0]->GetZaxis()->SetTitleOffset(0.9);

  fieldPalette->cd();
  
  hEvsTime[0]->Draw("colz");

  gPad->Update();
  TPaletteAxis *palette = (TPaletteAxis*)hEvsTime[0]->GetListOfFunctions()->FindObject("palette");
  palette->SetX1NDC(1-margins[2]+0.02);
  palette->SetX2NDC(1-margins[2]+0.05);
  palette->SetY1NDC(margins[1]);
  palette->SetY2NDC(1-margins[3]);
  C->Modified();
  C->Update();

  // Plot the crossings
  if(!opt.Contains("nocross")) {
    for(Int_t ic=0;ic<NCross[0];ic++) {
      if( gEcross[0][ic] )
	gEcross[0][ic]->Draw("L");
    }
  }
  
  
  gPad->RedrawAxis("G");

  // Print to a file
  // Output file
  TString fOutName = Form("./%s/Plots/DDR/DDR-Ez-%s",sim.Data(),sim.Data());
  PGlobals::imgconv(C,fOutName,opt);
  // ---------------------------------------------------------


  Float_t zMin = hEvsTime[0]->GetXaxis()->GetXmin();
  Float_t zMax = hEvsTime[0]->GetXaxis()->GetXmax();
  
  if(hEvsTime[0]) {
    
    // Canvas setup
    TCanvas *C1 = new TCanvas("C1","Evolution of the Ez extremes and transformer ratio",sizex,sizey);
    gPad->SetTickx(0);
    gPad->SetTicky(0);


    TLegend *Leg = new TLegend(0.20,0.20,0.35,0.40);
    PGlobals::SetPaveStyle(Leg);
    Leg->SetTextAlign(22);
    Leg->SetTextColor(kGray+3);
    Leg->SetLineColor(kGray+1);
    Leg->SetBorderSize(0);
    Leg->SetFillColor(0);
    Leg->SetFillStyle(1001);
    for(Int_t ic=0;ic<NCross[0];ic++) {
      // Only minimums:
      // if(ic%2-1) continue;
      Leg->AddEntry(gEextr[0][ic],Form("Minimum %2i",ic/2+1),"L");
    }
    Leg->SetTextColor(kGray+3);
    
    
    C1->cd(0);
  
    Float_t margin = (maxEextr - minEextr)/10;
    Float_t Ez0   = minEextr - margin;
    Float_t Ez1   = maxEextr + margin;
    TH1F *hFrame2 = new TH1F("hFrame2","",100,zMin,zMax);
    // TH2F *hFrame2 = (TH2F*) hEvsTime[0]->Clone("hFrame2");
  
    hFrame2->GetYaxis()->SetRangeUser(Ez0,Ez1);
    hFrame2->GetXaxis()->SetTitle(hEvsTime[0]->GetXaxis()->GetTitle());
    hFrame2->GetYaxis()->SetTitle("E_{z}/E_{0}");
    PGlobals::SetH1LabelSize(hFrame2);
    
    hFrame2->Draw("axis");

    gPad->Update();

    if(!opt.Contains("nocross")) {

      for(Int_t ic=0;ic<NCross[0];ic++) {
      
	// Only minimums:
	// if(ic%2-1) continue;
	if( !gEextr[0][ic] ) continue;

	// Remove crazy points
	Double_t *yEextr  = gEextr[0][ic]->GetY();

	for(Int_t k=0; k < gEextr[0][ic]->GetN(); k++) {
	  if(fabs(yEextr[k])>1000) {
	    gEextr[0][ic]->RemovePoint(k);
	    k--;
	  }
	}
      
	gEextr[0][ic]->Draw("L");
      
      }
    }


    // Transformer ratio
    Double_t rightmin = 0.0;
    Double_t rightmax = 12.0;
    Double_t slope = (gPad->GetUymax() - gPad->GetUymin())/(rightmax-rightmin); 
  
    Double_t *xTR = gTRatio->GetX();
    Double_t *yTR = gTRatio->GetY();
    for(Int_t j=0;j<gTRatio->GetN();j++) {
      gTRatio->SetPoint(j,xTR[j],yTR[j]*slope + gPad->GetUymin());
    }
  
    gTRatio->SetLineWidth(2);
    gTRatio->SetLineColor(kRed);
    gTRatio->Draw("L");
  
    //draw an axis on the right side
    TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),gPad->GetUxmax(),
			      gPad->GetUymax(),rightmin,rightmax,505,"+L");
  
    axis->SetLineWidth(1);
    axis->SetLineColor(gTRatio->GetLineColor());
    axis->SetLabelColor(gTRatio->GetLineColor());
    axis->SetTitle("R_{max}");
    axis->CenterTitle();
    axis->SetTitleColor(gTRatio->GetLineColor());
  
    axis->Draw();


    C1->cd(0);
  
    // Print to a file
    // Output file
    fOutName = Form("./%s/Plots/DDR/DDR-EzMax-%s",sim.Data(),sim.Data());
    PGlobals::imgconv(C1,fOutName,opt);
  }

  if(hVvsTime) {
    // Canvas setup
    TCanvas *CV = new TCanvas("CV","Evolution of the wakefield potential",sizex,sizey);
  
    CV->cd(0);
    gPad->SetFrameLineWidth(2);  
    gPad->SetTickx(1);
    gPad->SetTicky(1);

    PGlobals::SetH1LabelSize(hVvsTime);

    // Change the range of z axis for the fields to be symmetric.
    Float_t Vmax = hVvsTime->GetMaximum();
    Float_t Vmin = hVvsTime->GetMinimum();
    if(Vmax > TMath::Abs(Vmin))
      Vmin = -Vmax;
    else
      Vmax = -Vmin;
    hVvsTime->GetZaxis()->SetRangeUser(Vmin,Vmax); 
 
    hVvsTime->GetYaxis()->SetLabelFont(43);
    hVvsTime->GetYaxis()->SetLabelSize(32);
    //  hVvsTime->GetYaxis()->SetLabelOffset(0.8);
    hVvsTime->GetYaxis()->SetTitleFont(43);
    hVvsTime->GetYaxis()->SetTitleSize(36);
    hVvsTime->GetYaxis()->SetTitleOffset(0.9);

    hVvsTime->GetXaxis()->SetLabelFont(43);
    hVvsTime->GetXaxis()->SetLabelSize(32);
    //  hVvsTime->GetXaxis()->SetLabelOffset(0.8);
    hVvsTime->GetXaxis()->SetTitleFont(43);
    hVvsTime->GetXaxis()->SetTitleSize(36);
    hVvsTime->GetXaxis()->SetTitleOffset(0.95);


    hVvsTime->GetZaxis()->SetLabelFont(43);
    hVvsTime->GetZaxis()->SetLabelSize(32);
    // hVvsTime->GetZaxis()->SetLabelOffset(0.8);
    hVvsTime->GetZaxis()->SetTitleFont(43);
    hVvsTime->GetZaxis()->SetTitleSize(36);
    hVvsTime->GetZaxis()->SetTitleOffset(0.9);

    fieldPalette->cd();
  
    // Float_t Vmax = hVvsTime->GetMaximum();
    // Float_t Vmin = hVvsTime->GetMinimum();
    // hVvsTime->GetZaxis()->SetRangeUser(Vmin,Vmax); 
    // if(Vmax<0.1) Vmax = 0.1;
  
    // // Dynamic potential palette
    // const Int_t potPNRGBs = 6;
    // const Int_t potPNCont = 64;
    // zeroPos = -Vmin/(Vmax-Vmin);
  
    // Double_t potPStops[potPNRGBs] = { 0.00, zeroPos-3.0/potPNCont,zeroPos-1.0/potPNCont, zeroPos, zeroPos+3.0/potPNCont, 1.00 };
    // Double_t potPRed[potPNRGBs]   = { 0.518, 0.965, 0.90,0.90, 0.498, 0.106 };
    // Double_t potPGreen[potPNRGBs] = { 0.078, 0.925, 0.90,0.90, 0.718, 0.078 };
    // Double_t potPBlue[potPNRGBs]  = { 0.106, 0.353, 0.90,0.90, 0.780, 0.518 };
   
    // PPalette * potentialPalette = (PPalette*) gROOT->FindObject("rbowinv");
    // potentialPalette->CreateGradientColorTable(potPNRGBs, potPStops, 
    // 					     potPRed, potPGreen, potPBlue, potPNCont);
  
    hVvsTime->Draw("colz");

    gPad->Update();
    palette = (TPaletteAxis*)hVvsTime->GetListOfFunctions()->FindObject("palette");
    palette->SetX1NDC(1-margins[2]+0.02);
    palette->SetX2NDC(1-margins[2]+0.05);
    palette->SetY1NDC(margins[1]);
    palette->SetY2NDC(1-margins[3]);
    CV->Modified();
    CV->Update();

    // Plot the crossings
    if(!opt.Contains("nocross")) {
    
      for(Int_t ic=0;ic<NVCross;ic++) {
	if( gVcross[ic] ) {	  
	  gVcross[ic]->Draw("L");
	}
      }
    
      for(Int_t ic=0;ic<NCross[0];ic++) {
	if( gEcross[0][ic] ) {
	  if(ic==0) {
	    gEcross[0][ic]->SetLineColor(myNaranja);
	    gEcross[0][ic]->SetLineWidth(3);
	  } else {
	    gEcross[0][ic]->SetLineColor(myBlue);
	    gEcross[0][ic]->SetLineWidth(3);
	  }

	  //	  gEcross[0][ic]->SetLineColor(kWhite);
	  gEcross[0][ic]->Draw("L");
	  
	}
	
      }
      
    }
  
  
    gPad->RedrawAxis("G");

    // Print to a file
    // Output file
    fOutName = Form("./%s/Plots/DDR/DDR-Psi-%s",sim.Data(),sim.Data());
    PGlobals::imgconv(CV,fOutName,opt);
    // ---------------------------------------------------------

    
    // Canvas setup
    TCanvas *CV1 = new TCanvas("CV1","Evolution of potential extremes",sizex,sizey);
    gPad->SetTickx(0);
    gPad->SetTicky(0);


    TLegend *Leg = new TLegend(0.20,0.20,0.35,0.40);
    PGlobals::SetPaveStyle(Leg);
    Leg->SetTextAlign(22);
    Leg->SetTextColor(kGray+3);
    Leg->SetLineColor(kGray+1);
    Leg->SetBorderSize(0);
    Leg->SetFillColor(0);
    Leg->SetFillStyle(1001);
    for(Int_t ic=0;ic<NVCross;ic++) {
      // Only minimums:
      // if(ic%2-1) continue;
      Leg->AddEntry(gVextr[ic],Form("Minimum %2i",ic/2+1),"L");
    }
    Leg->SetTextColor(kGray+3);
  
  
    CV1->cd(0);
  
    Float_t margin = (maxVextr - minVextr)/10;
    Float_t V0   = minVextr - margin;
    Float_t V1   = maxVextr + margin;
    TH1F *hFrameV = new TH1F("hFrameV","",100,hVvsTime->GetXaxis()->GetXmin(),hVvsTime->GetXaxis()->GetXmax());
    // TH2F *hFrameV = (TH2F*) hVvsTime->Clone("hFrameV");
  
    hFrameV->GetYaxis()->SetRangeUser(V0,V1);
    hFrameV->GetXaxis()->SetTitle(hVvsTime->GetXaxis()->GetTitle());
    hFrameV->GetYaxis()->SetTitle("#psi");
    PGlobals::SetH1LabelSize(hFrameV);

    hFrameV->GetYaxis()->SetLabelFont(43);
    hFrameV->GetYaxis()->SetLabelSize(32);
    //  hFrameV->GetYaxis()->SetLabelOffset(0.8);
    hFrameV->GetYaxis()->SetTitleFont(43);
    hFrameV->GetYaxis()->SetTitleSize(36);
    hFrameV->GetYaxis()->SetTitleOffset(0.9);

    hFrameV->GetXaxis()->SetLabelFont(43);
    hFrameV->GetXaxis()->SetLabelSize(32);
    //  hFrameV->GetXaxis()->SetLabelOffset(0.8);
    hFrameV->GetXaxis()->SetTitleFont(43);
    hFrameV->GetXaxis()->SetTitleSize(36);
    hFrameV->GetXaxis()->SetTitleOffset(0.95);


    hFrameV->GetZaxis()->SetLabelFont(43);
    hFrameV->GetZaxis()->SetLabelSize(32);
    // hFrameV->GetZaxis()->SetLabelOffset(0.8);
    hFrameV->GetZaxis()->SetTitleFont(43);
    hFrameV->GetZaxis()->SetTitleSize(36);
    hFrameV->GetZaxis()->SetTitleOffset(0.9);

    
    hFrameV->Draw("axis");

    gPad->Update();

    if(!opt.Contains("nocross")) {

      for(Int_t ic=0;ic<NVCross;ic++) {
      
	// Only minimums:
	// if(ic%2-1) continue;
	if( !gVextr[ic] ) continue;

	// Remove crazy points
	Double_t *yVextr  = gVextr[ic]->GetY();

	for(Int_t k=0; k < gVextr[ic]->GetN(); k++) {
	  if(fabs(yVextr[k])>1000) {
	    gVextr[ic]->RemovePoint(k);
	    k--;
	  }
	}
      
	// gVextr[ic]->Draw("L");
	// gVextr_alt[ic]->Draw("L");	
	gVextr_avg[ic]->Draw("C");
      
      }
    }


    CV1->cd(0);
  
    // Print to a file
    // Output file
    fOutName = Form("./%s/Plots/DDR/DDR-PsiMax-%s",sim.Data(),sim.Data());
    PGlobals::imgconv(CV1,fOutName,opt);
  


    // Special plot
    
    // Canvas setup
    TCanvas *CVDen = new TCanvas("CVDen","Evolution of potential extremes",sizex,sizey);
    gPad->SetTickx(0);
    gPad->SetTicky(0);
  
    CVDen->cd(0);

    gPad->SetLogx(1);

    
    TH1F *hFrame3 = new TH1F("hFrame3","",100,1.0,16.);
  
    hFrame3->GetYaxis()->SetRangeUser(V0,V1);
    hFrame3->GetXaxis()->SetTitle("n/n_{0}");
    hFrame3->GetYaxis()->SetTitle("#psi");
    PGlobals::SetH1LabelSize(hFrame3);

    hFrame3->GetYaxis()->SetLabelFont(43);
    hFrame3->GetYaxis()->SetLabelSize(32);
    //  hFrame3->GetYaxis()->SetLabelOffset(0.8);
    hFrame3->GetYaxis()->SetTitleFont(43);
    hFrame3->GetYaxis()->SetTitleSize(36);
    hFrame3->GetYaxis()->SetTitleOffset(0.9);

    hFrame3->GetXaxis()->SetLabelFont(43);
    hFrame3->GetXaxis()->SetLabelSize(32);
    //  hFrame3->GetXaxis()->SetLabelOffset(0.8);
    hFrame3->GetXaxis()->SetTitleFont(43);
    hFrame3->GetXaxis()->SetTitleSize(36);
    hFrame3->GetXaxis()->SetTitleOffset(0.95);


    hFrame3->GetZaxis()->SetLabelFont(43);
    hFrame3->GetZaxis()->SetLabelSize(32);
    // hFrame3->GetZaxis()->SetLabelOffset(0.8);
    hFrame3->GetZaxis()->SetTitleFont(43);
    hFrame3->GetZaxis()->SetTitleSize(36);
    hFrame3->GetZaxis()->SetTitleOffset(0.9);

    hFrame3->Draw("axis");

    gPad->Update();

    Int_t Np = gVextr_avg[1]->GetN();
    Double_t *yVmin  = gVextr_avg[1]->GetY();
    Double_t *yVmax  = gVextr_avg[0]->GetY();
    Double_t *xVextr  = gVextr_avg[1]->GetX();
    Double_t *den = new Double_t[Np];
    
    for(Int_t ip=0;ip<Np;ip++) {
      den[ip] = 1 + (ntop0-1) * TMath::Exp(-(xVextr[ip]*xVextr[ip])/(2*sigmal0*sigmal0));
    }
    
    gVminvsDen = new TGraph(Np,den,yVmin);
    gVminvsDen->SetLineStyle(1);
    gVminvsDen->SetLineColor(myBlue);
    gVminvsDen->SetLineWidth(3);

    gVminvsDen->Draw("C");

    gVmaxvsDen = new TGraph(Np,den,yVmax);
    gVmaxvsDen->SetLineStyle(1);
    gVmaxvsDen->SetLineColor(myNaranja);
    gVmaxvsDen->SetLineWidth(3);

    gVmaxvsDen->Draw("C");

    CVDen->cd(0);
  
    // Print to a file
    // Output file
    fOutName = Form("./%s/Plots/DDR/DDR-PsivsDen-%s",sim.Data(),sim.Data());
    PGlobals::imgconv(CVDen,fOutName,opt);
    
  }
  
  // Canvas setup
  TCanvas *CDZ = new TCanvas("CDZ","Evolution of the dephasing",sizex,sizey);
  CDZ->cd(0);
 
  gPad->SetTickx(0);
  gPad->SetTicky(0);

  Float_t margin = (maxEdephas - minEdephas)/10;
  Float_t D0   = minEdephas-margin;
  Float_t D1   = maxEdephas+margin;
  TH1F *hFrame4 = new TH1F("hFrame4","",100,zMin,zMax);
  
  hFrame4->GetYaxis()->SetRangeUser(D0,D1);
  hFrame4->GetXaxis()->SetTitle(hEvsTime[0]->GetXaxis()->GetTitle());
  hFrame4->GetYaxis()->SetTitle("#Delta#zeta");

  PGlobals::SetH1LabelSize(hFrame4);
  hFrame4->Draw("axis");
  
  gPad->Update();
     
  for(Int_t i=0;i<Nfields;i++) {
    if(i>0) continue;
    
    for(Int_t ic=0;ic<NCross[i];ic++) {
      
      // Only minimums:
      if(ic%2) continue;
      if( !gEdephas[i][ic] ) continue;
      
      gEdephas[i][ic]->Draw("L");
    }
    
    // Leg->Draw();
  }
  
  gPad->RedrawAxis("G");
  
  CDZ->cd(0);
  
  // Print to a file
  // Output file
  fOutName = Form("./%s/Plots/DDR/DDR-EzDephas-%s",sim.Data(),sim.Data());
  PGlobals::imgconv(CDZ,fOutName,opt);
  // ---------------------------------------------------------


  {
    Int_t sizex = 383;
    Int_t sizey = 403;
    
    // Canvas setup
    TCanvas *CDDR = new TCanvas("CDDR","DDR vs phase velocity",sizex,sizey);
    CDDR->cd(0);
 
    gPad->SetTickx(0);
    gPad->SetTicky(0);
    gPad->SetFillStyle(4000);
    gPad->SetFrameFillStyle(4000);

    Int_t Np = 100;
    Double_t *zarray = new Double_t[Np];
    Double_t *denarray = new Double_t[Np];
    Double_t *betapharray = new Double_t[Np];
    Double_t ntop   = 10;
    Double_t sigmal = 2.0;
    Double_t phase = -TMath::TwoPi();

    Float_t dMin = -0.5;
    Float_t dMax = ntop+0.5;
    
    for(Int_t i=0; i<Np; i++) {
      zarray[i] = (i + 0.5) * 6*sigmal/Np;
      denarray[i] = DenGauss(zarray[i],sigmal,ntop);
      betapharray[i] = betaph(zarray[i],phase,sigmal,ntop);
    }

    TGraph *denvsz = new TGraph(Np,zarray,denarray);
    denvsz->SetLineColor(myNaranja);
    denvsz->SetLineWidth(3);
      
    TH1F *hFrameD = new TH1F("hFrameD","",100,0.0,6*sigmal);
    hFrameD->GetYaxis()->SetRangeUser(dMin,dMax);
    hFrameD->GetXaxis()->SetTitle("k_{p}^{0} z");
    hFrameD->GetYaxis()->SetTitle("n/n_{0}");

    PGlobals::SetH1LabelSize(hFrameD);
    hFrameD->GetYaxis()->SetNdivisions(6);
    hFrameD->GetYaxis()->SetAxisColor(myNaranja);
    hFrameD->GetYaxis()->SetLabelColor(myNaranja);
    hFrameD->GetYaxis()->SetTitleColor(myNaranja);

    hFrameD->Draw("axis");

    denvsz->Draw("C");

    gPad->Update();

    // Draw velocities
    Double_t rightmin = -0.05;
    Double_t rightmax = 1.05;
    Double_t slope = (gPad->GetUymax() - gPad->GetUymin())/(rightmax-rightmin); 
    for(Int_t i=0; i<Np; i++) {
      betapharray[i] = (betapharray[i]-rightmin) * slope +  gPad->GetUymin();
    }
    
    TGraph *betaphvsz = new TGraph(Np,zarray,betapharray); 
    //   betaphvsz->SetLineColor(myBlue);
    betaphvsz->SetLineColor(kGray+2);
    betaphvsz->SetLineWidth(3);
    betaphvsz->SetLineStyle(2);
    betaphvsz->Draw("C");


    // Draw beta from simulation 
    Double_t *betaarray = new Double_t[Np];
    Double_t *den0array = new Double_t[Np];
    Double_t nfactor = (ntop0 - 1 + TMath::Exp(1./2.)) / (ntop - 1 + TMath::Exp(1./2.));
    //    Double_t nfactor = 4;

    Int_t Npp = gVminvsDen->GetN();
    Double_t *psival = gVminvsDen->GetY();
    Double_t *nval = gVminvsDen->GetX();
    // cout << Form(" Density factor = %f ",nfactor) << endl;

    // for(Int_t j=0; j<Npp; j++)
    //   cout << Form(" nval = %f",nval[j]) << endl;
    
    for(Int_t i=0; i<Np; i++) {
      den0array[i] = denarray[i] * nfactor;

      Bool_t found = kFALSE;
      for(Int_t j=1; j<Npp; j++) {
	if( den0array[i]<nval[j-1] && den0array[i]>nval[j]) {
	  // interpolate
	  Double_t psi_int = (psival[j]-psival[j-1])/(nval[j]-nval[j-1]) * (den0array[i] - nval[j-1]) + psival[j-1];

	  betaarray[i] = betapsi(psi_int);
	  betaarray[i] = (betaarray[i]-rightmin) * slope +  gPad->GetUymin();

	  //	  cout << Form(" Point %i: %f --> %.4f (%.4f) %.4f (%.4f) -> %.4f ",i,den0array[i],nval[j-1],psival[j-1],nval[j],psival[j],psi_int) << endl;
	  
	  
	  found = kTRUE;
      	}

	if(found) break;
      }

      if(found) continue;
      
      //  cout << Form("  %f ", den0array[i]) << endl;
    }
    
    TGraph *betavsz = new TGraph(Np,zarray,betaarray); 
    betavsz->SetLineColor(myBlue);
    betavsz->SetLineWidth(3);
    betavsz->SetLineStyle(1);
    betavsz->Draw("C");
    
    //draw an axis on the right side
    TGaxis *axisb = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),gPad->GetUxmax(),
			      gPad->GetUymax(),rightmin,rightmax,206,"+L");
  
    axisb->SetLineWidth(1);
    axisb->SetLineColor(betaphvsz->GetLineColor());
    axisb->SetLabelColor(betaphvsz->GetLineColor());
    axisb->SetTitle("#beta");
    axisb->CenterTitle();
    axisb->SetTitleColor(betaphvsz->GetLineColor());
    axisb->SetLabelFont(43);
    axisb->SetLabelColor(kGray+2);
    axisb->SetLabelSize(14);
    axisb->SetTitleOffset(1.6);
    axisb->SetLabelOffset(0.02);
    axisb->SetTickSize(0.02);
    axisb->Draw();

    TBox *lFrame = new TBox(gPad->GetUxmin(), gPad->GetUymin(),
			    gPad->GetUxmax(), gPad->GetUymax());
    lFrame->SetFillStyle(0);
    lFrame->SetLineColor(PGlobals::frameColor);
    lFrame->SetLineWidth(PGlobals::frameWidth);
    lFrame->Draw();
    
    TLine laxis(gPad->GetUxmin(), gPad->GetUymin(),
		gPad->GetUxmin(), gPad->GetUymax());
    laxis.SetLineColor(myNaranja);
    laxis.SetLineWidth(3);
    laxis.Draw();
    
    TLine raxis(gPad->GetUxmax(), gPad->GetUymin(),
		gPad->GetUxmax(), gPad->GetUymax());
    raxis.SetLineColor(kGray+2);
    raxis.SetLineWidth(3);
    raxis.Draw();

    
    gPad->Update();
     
  
  
    gPad->RedrawAxis("G");
  
    CDDR->cd(0);
  
    // Print to a file
    // Output file
    fOutName = Form("./%s/Plots/DDR/DDR-DenVsVel-%s",sim.Data(),sim.Data());
    PGlobals::imgconv(CDDR,fOutName,opt);
    // ---------------------------------------------------------
  }
  

  
  ifile->Close();
  cout << endl;
}
  
