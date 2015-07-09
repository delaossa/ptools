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
#include <TProfile.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TExec.h>
#include <TGaxis.h>

#include "PData.hh"
#include "PlasmaGlob.hh"
#include "PPalette.hh"

void PlotPitzPanel( const TString &sim, Int_t time, Float_t Emin, Float_t Emax, Int_t zoom=2, Int_t NYbins=2, const TString &opt="") {
  
  PlasmaGlob::Initialize();

  // Palettes!
  gROOT->Macro("PlasmaPalettes.C");

  // More makeup
  gStyle->SetPadTopMargin(0.07); 

  // Load PData
  PData *pData = PData::Get(sim.Data());
  pData->LoadFileNames(time);
  if(!pData->IsInit()) return;

  Bool_t ThreeD = kFALSE;
  if(sim.Contains("3D")) ThreeD = kTRUE; 

  // Binning for Spectrum
  Int_t x1Nbin, x2Nbin; x1Nbin = x2Nbin = 0;
  Double_t x1Range, x1Mid, x1Min, x1Max;
  x1Range = x1Mid = x1Min = x1Max = 0.0;
  Double_t x2Range, x2Mid, x2Min, x2Max;
  x2Range = x2Mid = x2Min = x2Max = 0.0;
  
  // Get charge density histos
  Int_t Nspecies = pData->NSpecies();
  TH2F **hDen2D = new TH2F*[Nspecies];
  TH1F **hDen1D = new TH1F*[Nspecies];
  TH1F **hEnergy = new TH1F*[Nspecies];    // Energy spectrum
  TH2F **hEvsX1 = new TH2F*[Nspecies];     // Energy vs x1
  TH2F **hEvsX2 = new TH2F*[Nspecies];     // Energy vs x2
  for(Int_t i=0;i<Nspecies;i++) {
    hDen2D[i] = NULL;
    hDen1D[i] = NULL;
    hEnergy[i] = NULL;
    hEvsX1[i] = hEvsX2[i] = NULL;

    if(!pData->GetChargeFileName(i)) 
      continue;
    
    if(!ThreeD)
      hDen2D[i] = pData->GetCharge(i);
    else
      hDen2D[i] = pData->GetCharge2DSliceZY(i);
    
    char hName[24];
    sprintf(hName,"hDen2D_%i",i);
    hDen2D[i]->SetName(hName);
    hDen2D[i]->GetXaxis()->CenterTitle();
    hDen2D[i]->GetYaxis()->CenterTitle();
    hDen2D[i]->GetXaxis()->SetTitle("z [c/#omega_{p}]");
    hDen2D[i]->GetYaxis()->SetTitle("y [c/#omega_{p}]");

    if(i==1) {
      x1Nbin    = hDen2D[i]->GetNbinsX();
      x2Nbin    = hDen2D[i]->GetNbinsY();
      
      x2Range   = (hDen2D[i]->GetYaxis()->GetXmax() - hDen2D[i]->GetYaxis()->GetXmin())/zoom;
      x2Mid     = (hDen2D[i]->GetYaxis()->GetXmax() + hDen2D[i]->GetYaxis()->GetXmin())/2.;
      x2Min     = x2Mid - x2Range/2;
      x2Max     = x2Mid + x2Range/2;
      
      x1Min     = hDen2D[i]->GetXaxis()->GetXmin();
      x1Max     = hDen2D[i]->GetXaxis()->GetXmax();
    }
    
    { // Center the y coordinate
      Float_t midY = 0.5 * (hDen2D[i]->GetYaxis()->GetXmin()+hDen2D[i]->GetYaxis()->GetXmax());
      Int_t NbinsX = hDen2D[i]->GetNbinsX();
      Float_t xMin = hDen2D[i]->GetXaxis()->GetXmin();
      Float_t xMax = hDen2D[i]->GetXaxis()->GetXmax();
      Int_t NbinsY = hDen2D[i]->GetNbinsY();
      Float_t yMin = hDen2D[i]->GetYaxis()->GetXmin() - midY;
      Float_t yMax = hDen2D[i]->GetYaxis()->GetXmax() - midY;
      hDen2D[i]->SetBins(NbinsX,xMin,xMax,NbinsY,yMin,yMax);
    }

    // Change to co-moving coordinate
    if(opt.Contains("comov")) {
      Int_t NbinsX = hDen2D[i]->GetNbinsX();
      Float_t xMin = hDen2D[i]->GetXaxis()->GetXmin()-hDen2D[i]->GetXaxis()->GetXmax();
      Float_t xMax = hDen2D[i]->GetXaxis()->GetXmax()-hDen2D[i]->GetXaxis()->GetXmax();
      Int_t NbinsY = hDen2D[i]->GetNbinsY();
      Float_t yMin = hDen2D[i]->GetYaxis()->GetXmin();
      Float_t yMax = hDen2D[i]->GetYaxis()->GetXmax();
      hDen2D[i]->SetBins(NbinsX,xMin,xMax,NbinsY,yMin,yMax);
      hDen2D[i]->GetXaxis()->SetTitle("#zeta [c/#omega_{p}]");
    }
    
    // Chaning to user units: 
    if(opt.Contains("units") && pData->GetPlasmaDensity()) {
      Int_t NbinsX = hDen2D[i]->GetNbinsX();
      Float_t xMin = 1e3 * pData->GetPlasmaSkinDepth() * hDen2D[i]->GetXaxis()->GetXmin();
      Float_t xMax = 1e3 * pData->GetPlasmaSkinDepth() * hDen2D[i]->GetXaxis()->GetXmax();
      Int_t NbinsY = hDen2D[i]->GetNbinsY();
      Float_t yMin = 1e3 * pData->GetPlasmaSkinDepth() * hDen2D[i]->GetYaxis()->GetXmin();
      Float_t yMax = 1e3 * pData->GetPlasmaSkinDepth() * hDen2D[i]->GetYaxis()->GetXmax();
      hDen2D[i]->SetBins(NbinsX,xMin,xMax,NbinsY,yMin,yMax);
      
      for(Int_t j=0;j<hDen2D[i]->GetNbinsX();j++) {
	for(Int_t k=0;k<hDen2D[i]->GetNbinsY();k++) {
	  hDen2D[i]->SetBinContent(j,k,1e-15 * 1e-6 * pData->GetPlasmaDensity() * hDen2D[i]->GetBinContent(j,k));
	}
      }
      
      hDen2D[i]->GetYaxis()->SetTitle("y [mm]");      
      if(opt.Contains("comov"))
	hDen2D[i]->GetXaxis()->SetTitle("#zeta [mm]");
      else
	hDen2D[i]->GetXaxis()->SetTitle("z [mm]");
      
    }

    // 1D - on axis:
    Int_t FirstyBin = hDen2D[i]->GetNbinsY()/2 - NYbins;
    Int_t LastyBin =  hDen2D[i]->GetNbinsY()/2 + (NYbins-1);
    hDen1D[i] = (TH1F*) hDen2D[i]->ProjectionX("hDen1D",FirstyBin,LastyBin);
    hDen1D[i]->Scale(1.0/(LastyBin-FirstyBin+1));
    
    sprintf(hName,"hDen1D_%i",i);
    hDen1D[i]->SetName(hName);
    hDen1D[i]->GetXaxis()->SetTitle("z [c/#omega_{p}]");
    hDen1D[i]->GetYaxis()->SetTitle("dQ [ou]");
    
    if(opt.Contains("comov")) {
      hDen1D[i]->GetXaxis()->SetTitle("#zeta [c/#omega_{p}]");
    }
    
    // Chaning to user units: 
    if(opt.Contains("units") && pData->GetPlasmaDensity()) {
      hDen1D[i]->GetYaxis()->SetTitle("dQ [10^{15}/cc]");
      
      if(opt.Contains("comov"))
	hDen1D[i]->GetXaxis()->SetTitle("#zeta [mm]");
      else
	hDen1D[i]->GetXaxis()->SetTitle("z [mm]");
      
    }
    
    if(!pData->GetRawFileName(i))
      continue;
    
    sprintf(hName,"hEnergy_%i",i);
    hEnergy[i] = (TH1F*) gROOT->FindObject(hName);
    if(hEnergy[i]) { delete hEnergy[i]; hEnergy[i] = NULL; }
    hEnergy[i] = new TH1F(hName,"",200,Emin,Emax);
    
    pData->GetH1Raw(pData->GetRawFileName(i)->c_str(),"ene",hEnergy[i]);
    hEnergy[i]->GetXaxis()->CenterTitle();
    hEnergy[i]->GetYaxis()->CenterTitle();
    hEnergy[i]->GetXaxis()->SetTitle("Ene [MeV]");
    hEnergy[i]->GetYaxis()->SetTitle("dQ/dE [a.u.]");
    PlasmaGlob::SetH1Style(hEnergy[i],i);
    
    sprintf(hName,"hEvsX1_%i",i);
    hEvsX1[i] = (TH2F*) gROOT->FindObject(hName);
    if(hEvsX1[i]) { delete hEvsX1[i]; hEvsX1[i] = NULL; }
    hEvsX1[i] = new TH2F(hName,"",x1Nbin,x1Min,x1Max,200,Emin,Emax);
    
    pData->GetH2Raw(pData->GetRawFileName(i)->c_str(),"x1","ene",hEvsX1[i]);
   
    hEvsX1[i]->GetXaxis()->CenterTitle();
    hEvsX1[i]->GetYaxis()->CenterTitle();
    hEvsX1[i]->GetXaxis()->SetTitle("z [c/#omega_{p}]");
    hEvsX1[i]->GetYaxis()->SetTitle("Ene [MeV]");
    PlasmaGlob::SetH1Style(hEvsX1[i],i);
        
    // Change to co-moving coordinate
    if(opt.Contains("comov")) {
      Int_t NbinsX = hEvsX1[i]->GetNbinsX();
      Float_t xMin = hEvsX1[i]->GetXaxis()->GetXmin()-hEvsX1[i]->GetXaxis()->GetXmax();
      Float_t xMax = hEvsX1[i]->GetXaxis()->GetXmax()-hEvsX1[i]->GetXaxis()->GetXmax();
      Int_t NbinsY = hEvsX1[i]->GetNbinsY();
      Float_t yMin = hEvsX1[i]->GetYaxis()->GetXmin();
      Float_t yMax = hEvsX1[i]->GetYaxis()->GetXmax();
      hEvsX1[i]->SetBins(NbinsX,xMin,xMax,NbinsY,yMin,yMax);
      hEvsX1[i]->GetXaxis()->SetTitle("#zeta [c/#omega_{p}]");
    }
    
    // Chaning to user units: 
    if(opt.Contains("units") && pData->GetPlasmaDensity()) {
      Int_t NbinsX = hEvsX1[i]->GetNbinsX();
      Float_t xMin = 1e3 * pData->GetPlasmaSkinDepth() * hEvsX1[i]->GetXaxis()->GetXmin();
      Float_t xMax = 1e3 * pData->GetPlasmaSkinDepth() * hEvsX1[i]->GetXaxis()->GetXmax();
      Int_t NbinsY = hEvsX1[i]->GetNbinsY();
      Float_t yMin = hEvsX1[i]->GetYaxis()->GetXmin();
      Float_t yMax = hEvsX1[i]->GetYaxis()->GetXmax();
      hEvsX1[i]->SetBins(NbinsX,xMin,xMax,NbinsY,yMin,yMax);
      
      for(Int_t j=0;j<hEvsX1[i]->GetNbinsX();j++) {
	for(Int_t k=0;k<hEvsX1[i]->GetNbinsY();k++) {
	  hEvsX1[i]->SetBinContent(j,k,1e-15 * 1e-6 * pData->GetPlasmaDensity() * hEvsX1[i]->GetBinContent(j,k));
	}
      }
      
      if(opt.Contains("comov"))
	hEvsX1[i]->GetXaxis()->SetTitle("#zeta [mm]");
      else
	hEvsX1[i]->GetXaxis()->SetTitle("z [mm]");
      
    }
    
    sprintf(hName,"hEvsX2_%i",i);
    hEvsX2[i] = (TH2F*) gROOT->FindObject(hName);
    if(hEvsX2[i]) { delete hEvsX2[i]; hEvsX2[i] = NULL; }
    hEvsX2[i] = new TH2F(hName,"",x2Nbin,x2Min,x2Max,200,Emin,Emax);
    
    pData->GetH2Raw(pData->GetRawFileName(i)->c_str(),"x2","ene",hEvsX2[i]);
    
    hEvsX2[i]->GetXaxis()->CenterTitle();
    hEvsX2[i]->GetYaxis()->CenterTitle();
    hEvsX2[i]->GetXaxis()->SetTitle("y [c/#omega_{p}]");
    hEvsX2[i]->GetYaxis()->SetTitle("Ene [MeV]");
    PlasmaGlob::SetH1Style(hEvsX2[i],i);
    
    { // Center the y coordinate
      Float_t midY = 0.5 * (hEvsX2[i]->GetXaxis()->GetXmin()+hEvsX2[i]->GetXaxis()->GetXmax());
      Int_t NbinsX = hEvsX2[i]->GetNbinsX();
      Float_t xMin = hEvsX2[i]->GetXaxis()->GetXmin() - midY;
      Float_t xMax = hEvsX2[i]->GetXaxis()->GetXmax() - midY;
      Int_t NbinsY = hEvsX2[i]->GetNbinsY();
      Float_t yMin = hEvsX2[i]->GetYaxis()->GetXmin();
      Float_t yMax = hEvsX2[i]->GetYaxis()->GetXmax();
      hEvsX2[i]->SetBins(NbinsX,xMin,xMax,NbinsY,yMin,yMax);
    }

    // Chaning to user units: 
    if(opt.Contains("units") && pData->GetPlasmaDensity()) {
      Int_t NbinsX = hEvsX2[i]->GetNbinsX();
      Float_t xMin = 1e3 * pData->GetPlasmaSkinDepth() * hEvsX2[i]->GetXaxis()->GetXmin();
      Float_t xMax = 1e3 * pData->GetPlasmaSkinDepth() * hEvsX2[i]->GetXaxis()->GetXmax();
      Int_t NbinsY = hEvsX2[i]->GetNbinsY();
      Float_t yMin = hEvsX2[i]->GetYaxis()->GetXmin();
      Float_t yMax = hEvsX2[i]->GetYaxis()->GetXmax();
      hEvsX2[i]->SetBins(NbinsX,xMin,xMax,NbinsY,yMin,yMax);
      
      for(Int_t j=0;j<hEvsX2[i]->GetNbinsX();j++) {
	for(Int_t k=0;k<hEvsX2[i]->GetNbinsY();k++) {
	  hEvsX2[i]->SetBinContent(j,k,1e-15 * 1e-6 * pData->GetPlasmaDensity() * hEvsX2[i]->GetBinContent(j,k));
	}
      }

      hEvsX2[i]->GetXaxis()->SetTitle("y [mm]");      
      
    }
    
  }
  
  TProfile *hDen2Dprof = NULL;
  TH1F *hRms = NULL;
  if(hDen2D[1]) {
    TString pname = hDen2D[1]->GetName();
    pname += "_pfx";

    hDen2Dprof =  (TProfile*) gROOT->FindObject(pname.Data());
    if(hDen2Dprof) delete hDen2Dprof;
    hDen2Dprof = hDen2D[1]->ProfileX("_pfx",1,-1,"s");

    hRms = (TH1F*) gROOT->FindObject("hRms");
    if(hRms) delete hRms;

    hRms = new TH1F("hRms","",hDen2D[1]->GetNbinsX(),hDen2D[1]->GetXaxis()->GetXmin(),
		    hDen2D[1]->GetXaxis()->GetXmax());

    for(Int_t j=0;j<hRms->GetNbinsX();j++) {
      if(opt.Contains("units") && pData->GetBeamRmsR())
	hRms->SetBinContent(j,hDen2Dprof->GetBinError(j)/(pData->GetBeamRmsR() * 1e-3) ); 
      else if(pData->GetBeamRmsZ())
	hRms->SetBinContent(j,hDen2Dprof->GetBinError(j) / pData->GetBeamRmsR() 
			    * pData->GetPlasmaSkinDepth() * 1e6 ); 
      
    }
    
    hRms->GetXaxis()->SetTitle("z [c/#omega_{p}]");
    if(opt.Contains("comov"))
      hRms->GetXaxis()->SetTitle("#zeta [c/#omega_{p}]");
    
    if(opt.Contains("units") && pData->GetBeamRmsR()) {
      hRms->GetXaxis()->SetTitle("z [mm]");
      if(opt.Contains("comov"))
	hRms->GetXaxis()->SetTitle("#zeta [mm]");
    }

    hRms->GetYaxis()->SetTitle("r_{b}/r_{0}");
    
  }


  // Get electric fields
  const Int_t Nfields = 2;
  TH2F **hE2D = new TH2F*[Nfields];
  TH1F **hE1D = new TH1F*[Nfields];
  for(Int_t i=0;i<Nfields;i++) {
    hE2D[i] = NULL;

    if(!pData->GetEfieldFileName(i))
      continue;
    
    if(!ThreeD)
      hE2D[i] = pData->GetEField(i);
    else
	hE2D[i] = pData->GetEField2DSliceZY(i);
    
    char hName[24];
    sprintf(hName,"hE2D_%i",i);
    hE2D[i]->SetName(hName);   
    hE2D[i]->GetXaxis()->CenterTitle();
    hE2D[i]->GetYaxis()->CenterTitle();
    hE2D[i]->GetXaxis()->SetTitle("z [c/#omega_{p}]");
    hE2D[i]->GetYaxis()->SetTitle("y [c/#omega_{p}]");
    
    { // Center the y coordinate
      Float_t midY = 0.5 * (hE2D[i]->GetYaxis()->GetXmin()+hE2D[i]->GetYaxis()->GetXmax());
      Int_t NbinsX = hE2D[i]->GetNbinsX();
      Float_t xMin = hE2D[i]->GetXaxis()->GetXmin();
      Float_t xMax = hE2D[i]->GetXaxis()->GetXmax();
      Int_t NbinsY = hE2D[i]->GetNbinsY();
      Float_t yMin = hE2D[i]->GetYaxis()->GetXmin() - midY;
      Float_t yMax = hE2D[i]->GetYaxis()->GetXmax() - midY;
      hE2D[i]->SetBins(NbinsX,xMin,xMax,NbinsY,yMin,yMax);
    }

    // Change to co-moving coordinate
    if(opt.Contains("comov")) {
      Int_t NbinsX = hE2D[i]->GetNbinsX();
      Float_t xMin = hE2D[i]->GetXaxis()->GetXmin()-hE2D[i]->GetXaxis()->GetXmax();
      Float_t xMax = hE2D[i]->GetXaxis()->GetXmax()-hE2D[i]->GetXaxis()->GetXmax();
      Int_t NbinsY = hE2D[i]->GetNbinsY();
      Float_t yMin = hE2D[i]->GetYaxis()->GetXmin();
      Float_t yMax = hE2D[i]->GetYaxis()->GetXmax();
      hE2D[i]->SetBins(NbinsX,xMin,xMax,NbinsY,yMin,yMax);
      hE2D[i]->GetXaxis()->SetTitle("#zeta [c/#omega_{p}]");
    }
  
    // Chaning to user units: 
    if(opt.Contains("units") && pData->GetPlasmaDensity()) {
      Int_t NbinsX = hE2D[i]->GetNbinsX();
      Float_t xMin = 1e3 * pData->GetPlasmaSkinDepth() * hE2D[i]->GetXaxis()->GetXmin();
      Float_t xMax = 1e3 * pData->GetPlasmaSkinDepth() * hE2D[i]->GetXaxis()->GetXmax();
      Int_t NbinsY = hE2D[i]->GetNbinsY();
      Float_t yMin = 1e3 * pData->GetPlasmaSkinDepth() * hE2D[i]->GetYaxis()->GetXmin();
      Float_t yMax = 1e3 * pData->GetPlasmaSkinDepth() * hE2D[i]->GetYaxis()->GetXmax();
      hE2D[i]->SetBins(NbinsX,xMin,xMax,NbinsY,yMin,yMax);
      
      for(Int_t j=0;j<hE2D[i]->GetNbinsX();j++) {
	for(Int_t k=0;k<hE2D[i]->GetNbinsY();k++) {
	  hE2D[i]->SetBinContent(j,k,1e-9 * pData->GetPlasmaE0() * hE2D[i]->GetBinContent(j,k));
	}
      }

      hE2D[i]->GetYaxis()->SetTitle("y [mm]");	
      if(opt.Contains("comov"))
	hE2D[i]->GetXaxis()->SetTitle("#zeta [mm]");
      else
	hE2D[i]->GetXaxis()->SetTitle("z [mm]");
	
    }
    
    // 1D - on axis:
    Int_t FirstyBin = hE2D[i]->GetNbinsY()/2 - NYbins;
    Int_t LastyBin =  hE2D[i]->GetNbinsY()/2 + (NYbins-1);
    hE1D[i] = (TH1F*) hE2D[i]->ProjectionX("hE1D",FirstyBin,LastyBin);
    hE1D[i]->Scale(1.0/(LastyBin-FirstyBin+1));
    
    sprintf(hName,"hE1D_%i",i);
    hE1D[i]->SetName(hName);
    hE1D[i]->GetXaxis()->SetTitle("z [c/#omega_{p}]");
    hE1D[i]->GetYaxis()->SetTitle("E [E_{0}]");
    
    if(opt.Contains("comov")) {
      hE1D[i]->GetXaxis()->SetTitle("#zeta [c/#omega_{p}]");
    }
    
    // Chaning to user units: 
    if(opt.Contains("units") && pData->GetPlasmaDensity()) {
      hE1D[i]->GetYaxis()->SetTitle("E [GeV/m]");
      
      if(opt.Contains("comov"))
	hE1D[i]->GetXaxis()->SetTitle("#zeta [mm]");
      else
	hE1D[i]->GetXaxis()->SetTitle("z [mm]");
      
    }
    
  }
  
  // Tunning the Histograms
  // ---------------------
  
  Float_t Time = pData->GetRealTime();
   
  // Set the range of the histogram for maximum constrast
  Float_t density = 1;
  if(opt.Contains("units") && pData->GetPlasmaDensity())
    density = 1e-15 * 1e-6 * pData->GetPlasmaDensity();
  
  Float_t Max  = 1.1 * hDen2D[0]->GetMaximum();
  Float_t Base = density;
  Float_t Min  = 2.* Base - Max;
  if(Max >= 2. * Base) {
    Max = 2. * Base;
    Min = 2. * Base  - Max;
  } else if(Max<1.0 * Base) {
    Max = 1.1 * Base;
    Min = 0.;
  }
  hDen2D[0]->GetZaxis()->SetRangeUser(Min,Max);  

  Max  = 1.1 * hDen1D[0]->GetMaximum();
  Min  = 2.* Base - Max;
  if(Max >= 2. * Base) {
    Max = 2. * Base;
    Min = 2. * Base  - Max;
  } else if(Max<1.0 * Base) {
    Max = 1.1 * Base;
    Min = 0.;
  }
  hDen1D[0]->GetYaxis()->SetRangeUser(Min,Max);  

  // Zoom
  Float_t range    = (hDen2D[0]->GetYaxis()->GetXmax() - hDen2D[0]->GetYaxis()->GetXmin())/zoom;
  Float_t midPoint = (hDen2D[0]->GetYaxis()->GetXmax() + hDen2D[0]->GetYaxis()->GetXmin())/2.;
  
  hDen2D[0]->GetYaxis()->SetRangeUser(midPoint-range/2,midPoint+range/2);  
  hE2D[0]->GetYaxis()->SetRangeUser(midPoint-range/2,midPoint+range/2);  
  hE2D[1]->GetYaxis()->SetRangeUser(midPoint-range/2,midPoint+range/2);  
  
  
  // Redefine the ranges according to the newly created histograms
  for(Int_t i=0;i<Nspecies;i++) {
    if(hEvsX2[i]) {
      x2Range   = hEvsX2[i]->GetXaxis()->GetXmax() - hEvsX2[i]->GetXaxis()->GetXmin();
      x2Mid     = (hEvsX2[i]->GetXaxis()->GetXmax() + hEvsX2[i]->GetXaxis()->GetXmin())/2.;
      x2Min     = x2Mid - x2Range/2;
      x2Max     = x2Mid + x2Range/2;
    }
    
    if(hEvsX1[i]) {
      x1Min     = hEvsX1[i]->GetXaxis()->GetXmin();
      x1Max     = hEvsX1[i]->GetXaxis()->GetXmax();
    }
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


    // Right side:
    // xMax = x2Min + (x2Max-x2Min) * 0.9;
    // xMin = x2Max;
    // Left side:
    xMin = x2Min;
    xMax = x2Min +  (x2Max-x2Min) * 0.2;    
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
  TCanvas *C = new TCanvas("C","2D Charge density and Electric field",1500,1000);
  C->Divide(2,3);

  // Palettes setup
  TExec *exPlasma = new TExec("exPlasma","plasmaPalette->cd();"); 
  TExec *exElec   = new TExec("exElec","electronPalette->cd();");
  TExec *exField  = new TExec("exField","rbowPalette->cd();");

  // Text objects
  TPaveText *text1 = new TPaveText(0.6,0.94,0.88,1.0,"NDC");
  PlasmaGlob::SetPaveTextStyle(text1);
  if(opt.Contains("units") && pData->GetPlasmaDensity()) 
    text1->AddText("dQ/dydz [10^{15}/cc]");
  else
    text1->AddText("dQ/dydz [a.u.]");

  TPaveText *text2 = new TPaveText(0.6,0.94,0.88,1.0,"NDC");
  PlasmaGlob::SetPaveTextStyle(text2);  
  if(opt.Contains("units") && pData->GetPlasmaDensity()) 
    text2->AddText("Longitudinal E field [GeV/m]");
  else
    text2->AddText("Longitudinal E field [E_{0}]");

  TPaveText *text3 = new TPaveText(0.6,0.94,0.88,1.0,"NDC");
  PlasmaGlob::SetPaveTextStyle(text3);  
  if(opt.Contains("units") && pData->GetPlasmaDensity()) 
    text3->AddText("Transverse E field [GeV/m]");
  else
    text3->AddText("Transverse E field [E_{0}]");

  
  TPaveText *text4 = new TPaveText(0.6,0.94,0.88,1.0,"NDC");
  PlasmaGlob::SetPaveTextStyle(text4);  
  if(opt.Contains("units") && pData->GetPlasmaDensity())
    if(opt.Contains("comov"))
      text4->AddText("dQ/dEd#zeta [10^{15}/cc]");
    else
      text4->AddText("dQ/dEdz [10^{15}/cc]");
  else
    if(opt.Contains("comov"))
      text4->AddText("dQ/dEd#zeta [a.u.]");
    else
      text4->AddText("dQ/dEdz [a.u.]");

  TPaveText *text5 = new TPaveText(0.6,0.94,0.88,1.0,"NDC");
  PlasmaGlob::SetPaveTextStyle(text5); 
  if(opt.Contains("units") && pData->GetPlasmaDensity()) 
    text5->AddText("dQ/dEdy [10^{15}/cc]");
  else
    text5->AddText("dQ/dEdy [a.u.]");
 
  TPaveText *textTime = new TPaveText(0.7,0.82,0.85,0.88,"NDC");
  PlasmaGlob::SetPaveTextStyle(textTime); 
  char ctext[128];
  if(opt.Contains("units") && pData->GetPlasmaDensity()) 
    sprintf(ctext,"L = %5.1f mm", 1e3 * pData->GetPlasmaSkinDepth() * Time);
  else
    sprintf(ctext,"T = %5.1f 1/#omega_{p}",Time);
  textTime->AddText(ctext);

  TPaveText *textDen = new TPaveText(0.10,0.94,0.35,1.00,"NDC");
  PlasmaGlob::SetPaveTextStyle(textDen); 
  textDen->SetTextColor(kRed+1);
  sprintf(ctext,"#rho = %5.2f x 10^{15} / cc", 1e-6 * 1e-15 * pData->GetPlasmaDensity());
  textDen->AddText(ctext);
  

  
  // Actual Plotting!
  // ------------------------------------------------------------

  // Output file
  TString fOutName = Form("./%s/Plots/PitzPanel/PitzPanel",sim.Data());
  fOutName += Form("-%s_%i",sim.Data(),time);

  C->cd(1); // <--- Top Plot

  gPad->SetFrameLineWidth(2);  

  
  TH2F *hFrame = (TH2F*) gROOT->FindObject("hFrame1");
  if(hFrame) delete hFrame;
  hFrame = (TH2F*) hDen2D[0]->Clone("hFrame1");
  hFrame->Reset();
  hFrame->Draw("col");
  
  exPlasma->Draw();
  hDen2D[0]->Draw("col same");
  
  if(hDen2D[1]) {
    exElec->Draw();
    hDen2D[1]->Draw("colz same");
  }
  
  text1->Draw();
  textTime->Draw();
  textDen->Draw();

  gPad->RedrawAxis(); 

  C->cd(3); // <--- Mid Plot
  gPad->SetFrameLineWidth(2);  
  
  hFrame = (TH2F*) gROOT->FindObject("hFrame2");
  if(hFrame) delete hFrame;
  hFrame = (TH2F*) hE2D[0]->Clone("hFrame2");
  hFrame->Reset();

  hFrame->Draw("col");
  exField->Draw();
  hE2D[0]->Draw("colz same");

  text2->Draw();
  textTime->Draw();

  gPad->RedrawAxis(); 

  C->cd(5); // <--- Mid Plot
  gPad->SetFrameLineWidth(2); 
   
  hFrame = (TH2F*) gROOT->FindObject("hFrame8");
  if(hFrame) delete hFrame;
  hFrame = (TH2F*) hE2D[1]->Clone("hFrame8");
  hFrame->Reset();

  hFrame->Draw("col");
  exField->Draw();
  hE2D[1]->Draw("colz same");

  text3->Draw();
  textTime->Draw();

  gPad->RedrawAxis(); 

  C->cd(2);
  hFrame = (TH2F*) gROOT->FindObject("hFrame3");
  if(hFrame) delete hFrame;
  hFrame = (TH2F*) hEvsX1[1]->Clone("hFrame3");
  hFrame->Reset();
  hFrame->Draw("col");

  exElec->Draw();
  if(hEvsX1[1]) 
    hEvsX1[1]->Draw("col same z");
  if(hEvsX1[2] && Nspecies>=3)
    hEvsX1[2]->Draw("col same");  
  
  gEnergyX1->Draw("F");
  gEnergyX1->Draw("L");
  
  text4->Draw();
  textTime->Draw();
  gPad->RedrawAxis(); 
  gPad->RedrawAxis("G"); 

  C->cd(4);
  hFrame = (TH2F*) gROOT->FindObject("hFrame4");
  if(hFrame) delete hFrame;
  hFrame = (TH2F*) hEvsX2[1]->Clone("hFrame4");
  hFrame->Reset();
  hFrame->Draw("col");
  
  exElec->Draw();
  if(hEvsX2[1]) 
    hEvsX2[1]->Draw("col same z");
  if(hEvsX2[2] && Nspecies>=3)
    hEvsX2[2]->Draw("col same");  
  
  gEnergyX2->Draw("F");
  gEnergyX2->Draw("L");
  
  text5->Draw();
  textTime->Draw();
  gPad->RedrawAxis(); 
  gPad->RedrawAxis("G"); 
  
  C->cd(6); // <--- Bottom Plot
  gPad->SetGridy(0);
  gPad->SetGridx(1);
  gPad->SetFrameLineWidth(1);  

  Int_t beamC   = kAzure-5;
  Int_t fieldC  = kGray+1; // PlasmaGlob::fieldLine;
  Int_t fieldCb = kGray+1; // PlasmaGlob::fieldLine;
  
  hE1D[0]->SetLineWidth(1);
  hE1D[0]->GetYaxis()->CenterTitle();
  hE1D[0]->GetXaxis()->CenterTitle();
  hE1D[0]->SetLineStyle(1);
  hE1D[0]->SetLineColor(fieldC);
  hE1D[0]->SetMarkerStyle(20);

  hE1D[1]->GetYaxis()->CenterTitle();
  hE1D[1]->GetXaxis()->CenterTitle();
  hE1D[1]->SetLineStyle(3);
  hE1D[1]->SetLineColor(fieldCb);
  hE1D[1]->SetMarkerStyle(24);

  Float_t factor = 1.5;
  Float_t minimum = factor * hE1D[0]->GetMinimum();
  Float_t maximum = factor * hE1D[0]->GetMaximum();
    
  if(hE1D[1]->GetMaximum() > hE1D[0]->GetMaximum()) {
    maximum = factor * hE1D[1]->GetMaximum();
  } 

  if(hE1D[1]->GetMinimum() < hE1D[0]->GetMinimum()) {
    minimum = factor * hE1D[1]->GetMinimum();
  } 

  if( maximum >= TMath::Abs(minimum)) minimum = -maximum;
  else maximum = - minimum;
  
  hE1D[0]->GetYaxis()->SetRangeUser(minimum,maximum);
  hE1D[0]->Draw("C");
  hE1D[1]->Draw("C same");
  
  C->Update();
  
  TLine *line0 = new TLine(hE1D[0]->GetXaxis()->GetXmin(),
			   (gPad->GetUymin()+gPad->GetUymax())/2.,
			   hE1D[0]->GetXaxis()->GetXmax(),
			   (gPad->GetUymin()+gPad->GetUymax())/2.);
  line0->SetLineColor(kGray+1);
  line0->SetLineStyle(2);
  line0->Draw();
  
  Float_t rightmin = 0;
  Float_t rightmax = 2.5 * hDen1D[1]->GetMaximum();
  Float_t slope = (gPad->GetUymax() - gPad->GetUymin())/(rightmax-rightmin);

  for(Int_t i=0;i<hDen1D[1]->GetNbinsX();i++) {
    hDen1D[1]->SetBinContent(i+1,hDen1D[1]->GetBinContent(i+1)*slope + gPad->GetUymin());
  }

  hDen1D[1]->SetLineWidth(1);
  hDen1D[1]->SetLineColor(beamC);
  hDen1D[1]->Draw("same C");
  
  //draw an axis on the right side
  TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),gPad->GetUxmax(),
			    gPad->GetUymax(),0,rightmax,505,"+L");
  
  axis->SetLineWidth(1);
  axis->SetLineColor(beamC);
  axis->SetLabelColor(beamC);
  if(opt.Contains("units") && pData->GetPlasmaDensity()) 
    axis->SetTitle("dQ [10^{15}/cc]");
  else
    axis->SetTitle("dQ [ou]");
  axis->CenterTitle();
  axis->SetTitleColor(beamC);
  axis->SetTitleOffset(-0.5);
  
  axis->Draw();

  // Bunch RMS
  // Trick for nicer plotting
  Float_t bin1 = -1;
  Float_t bin2 = -1;
  for(Int_t i=0;i<hRms->GetNbinsX();i++) 
    if(bin1==-1 && hRms->GetBinContent(i+1)>0) bin1 = hRms->GetBinCenter(i+1);  
  for(Int_t i=hRms->GetNbinsX()-1;i>=0;i--) 
    if(bin2==-1 && hRms->GetBinContent(i+1)>0) bin2 = hRms->GetBinCenter(i+1);
  hRms->GetXaxis()->SetRangeUser(bin1,bin2);

  rightmin = 0;
  rightmax = 2.5 * hRms->GetMaximum();
  slope = (gPad->GetUymax() - gPad->GetUymin())/(rightmax-rightmin);
  
  for(Int_t i=0;i<hRms->GetNbinsX();i++) {
    hRms->SetBinContent(i+1,(hRms->GetBinContent(i+1)-rightmin)*slope + gPad->GetUymin());
  }
  
  hRms->SetLineStyle(3);
  hRms->SetLineColor(beamC);
  hRms->Draw("same C");

  //draw an axis on the right side
  Float_t rightmargin = 0.06;
  Float_t ux = gPad->PixeltoX(gPad->UtoPixel(1-rightmargin));
  TGaxis *axisRMS = new TGaxis(ux,gPad->GetUymin(),
			       ux,gPad->GetUymax(),rightmin,rightmax,505,"+L");
  
  axisRMS->SetLineWidth(1);
  axisRMS->SetLineColor(beamC);
  axisRMS->SetLabelColor(beamC);
  axisRMS->SetTitle("r_{b}/r_{0}");
  axisRMS->CenterTitle();
  axisRMS->SetTitleColor(beamC);
  axisRMS->SetTitleOffset(0.22);
  axisRMS->Draw();

  C->cd();

  // Print to a file
  PlasmaGlob::imgconv(C,fOutName,opt);
  // ---------------------------------------------------------

}
