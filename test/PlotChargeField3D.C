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
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TPaletteAxis.h>
#include <TExec.h>

#include "PData.hh"
#include "PlasmaGlob.hh"
#include "PPalette.hh"

void PlotChargeField3D( const TString &sim, Int_t time, Int_t zoom=2, Int_t Nbins=1, const TString &options="") {
  
  PlasmaGlob::Initialize();

  // Palettes!
  gROOT->Macro("PlasmaPalettes.C");

  // Init Units table
  PUnits::UnitsTable::Get();

  // Load PData
  PData *pData = PData::Get(sim.Data());
  pData->LoadFileNames(time);
  if(!pData->IsInit()) return;

  Bool_t ThreeD = kFALSE;
  if(sim.Contains("3D")) ThreeD = kTRUE; 
  else {
    cout << " It is not a 3D simulation: exiting..." << endl;
    return;
  }

  TString opt = options;
  
  // More makeup
  gStyle->SetPadGridY(0);
  if(opt.Contains("grid")) {
    gStyle->SetPadGridX(1);
    gStyle->SetPadGridY(1);
  }
  
  // Some plasma constants
  Double_t n0 = pData->GetPlasmaDensity();
  Double_t kp = pData->GetPlasmaK();
  Double_t skindepth = 1/kp;
  Double_t E0 = pData->GetPlasmaE0();

  // Some beam properties:
  Float_t Ebeam = pData->GetBeamEnergy() * PUnits::MeV;
  Float_t gamma = Ebeam / PConst::ElectronMassE;
  Float_t vbeam = TMath::Sqrt(1 - 1/(gamma*gamma));
  // cout << Form(" - Bunch gamma      = %8.4f", gamma ) << endl;
  // cout << Form(" - Bunch velocity   = %8.4f c", vbeam ) << endl;
  
  // Time in OU
  Float_t Time = pData->GetRealTime();

  // z start of the plasma in normalized units.
  Float_t zStartPlasma = pData->GetPlasmaStart()*kp;

  // z start of the beam in normalized units.
  Float_t zStartBeam = pData->GetBeamStart()*kp;
  Time -= zStartPlasma - zStartBeam;
  
  // Get charge density histos
  Int_t Nspecies = pData->NSpecies();
  TH2F **hDen2D = new TH2F*[Nspecies];
  TH2F **hDen2DXY = new TH2F*[Nspecies];
  for(Int_t i=0;i<Nspecies;i++) {
    hDen2D[i] = NULL;
    hDen2DXY[i] = NULL;

    if(!pData->GetChargeFileName(i)) 
      continue;
    
    char hName[24];
    sprintf(hName,"hDen_%i",i);
    hDen2D[i] = (TH2F*) gROOT->FindObject(hName);
    if(hDen2D[i]) delete hDen2D[i];
   
    hDen2D[i] = pData->GetCharge2DSliceZY(i,-1,Nbins);
  
    hDen2D[i]->SetName(hName);
    hDen2D[i]->GetXaxis()->CenterTitle();
    hDen2D[i]->GetYaxis()->CenterTitle();
    hDen2D[i]->GetZaxis()->CenterTitle();
    hDen2D[i]->GetXaxis()->SetTitle("z [c/#omega_{p}]");
    hDen2D[i]->GetYaxis()->SetTitle("y [c/#omega_{p}]");
    if(i==0)
      hDen2D[i]->GetZaxis()->SetTitle("#LTn_{e}#GT [n_{0}]");
    else
      hDen2D[i]->GetZaxis()->SetTitle("#LTn_{b}#GT [n_{0}]");

    sprintf(hName,"hDenXY_%i",i);
    hDen2DXY[i] = (TH2F*) gROOT->FindObject(hName);
    if(hDen2DXY[i]) delete hDen2D[i];

    hDen2DXY[i] = pData->GetCharge2DSliceXY(i,-1,Nbins);
    
    hDen2DXY[i]->SetName(hName);
    hDen2DXY[i]->GetXaxis()->CenterTitle();
    hDen2DXY[i]->GetYaxis()->CenterTitle();
    hDen2DXY[i]->GetZaxis()->CenterTitle();
    hDen2DXY[i]->GetXaxis()->SetTitle("x [c/#omega_{p}]");
    hDen2DXY[i]->GetYaxis()->SetTitle("y [c/#omega_{p}]");
    if(i==0)
      hDen2DXY[i]->GetZaxis()->SetTitle("#LTn_{e}#GT [n_{0}]");
    else
      hDen2DXY[i]->GetZaxis()->SetTitle("#LTn_{b}#GT [n_{0}]");
    

    { // Center the y and z coordinate
      Int_t NbinsX = hDen2D[i]->GetNbinsX();
      Float_t xMin = hDen2D[i]->GetXaxis()->GetXmin() - zStartPlasma;
      Float_t xMax = hDen2D[i]->GetXaxis()->GetXmax() - zStartPlasma;
      Int_t NbinsY = hDen2D[i]->GetNbinsY();
      Float_t yMin = hDen2D[i]->GetYaxis()->GetXmin();
      Float_t yMax = hDen2D[i]->GetYaxis()->GetXmax();

      if(!opt.Contains("cyl")) {
	Float_t midY = 0.5 * (hDen2D[i]->GetYaxis()->GetXmin()+hDen2D[i]->GetYaxis()->GetXmax());
	yMin -= midY;
	yMax -= midY;
      }

      hDen2D[i]->SetBins(NbinsX,xMin,xMax,NbinsY,yMin,yMax);

      NbinsX = hDen2DXY[i]->GetNbinsX();
      xMin = hDen2DXY[i]->GetXaxis()->GetXmin();
      xMax = hDen2DXY[i]->GetXaxis()->GetXmax();
      NbinsY = hDen2DXY[i]->GetNbinsY();
      yMin = hDen2DXY[i]->GetYaxis()->GetXmin();
      yMax = hDen2DXY[i]->GetYaxis()->GetXmax();

      if(!opt.Contains("cyl")) {
	Float_t midY = 0.5 * (hDen2DXY[i]->GetYaxis()->GetXmin()+hDen2DXY[i]->GetYaxis()->GetXmax());
	yMin -= midY;
	yMax -= midY;

	Float_t midX = 0.5 * (hDen2DXY[i]->GetXaxis()->GetXmin()+hDen2DXY[i]->GetXaxis()->GetXmax());
	xMin -= midX;
	xMax -= midX;
      }
      
      hDen2DXY[i]->SetBins(NbinsX,xMin,xMax,NbinsY,yMin,yMax);
    }


    // Change to co-moving coordinate
    if(opt.Contains("comov")) {
      Int_t NbinsX = hDen2D[i]->GetNbinsX();
      Float_t xMin = hDen2D[i]->GetXaxis()->GetXmin() - vbeam*Time;
      Float_t xMax = hDen2D[i]->GetXaxis()->GetXmax() - vbeam*Time;
      Int_t NbinsY = hDen2D[i]->GetNbinsY();
      Float_t yMin = hDen2D[i]->GetYaxis()->GetXmin();
      Float_t yMax = hDen2D[i]->GetYaxis()->GetXmax();
      hDen2D[i]->SetBins(NbinsX,xMin,xMax,NbinsY,yMin,yMax);
      hDen2D[i]->GetXaxis()->SetTitle("#zeta [c/#omega_{p}]");
    }
    
    // Chaning to user units: 
    if(opt.Contains("units") && pData->GetPlasmaDensity()) {
      Int_t NbinsZ = hDen2D[i]->GetNbinsX();
      Float_t zMin = 1e3 * pData->GetPlasmaSkinDepth() * hDen2D[i]->GetXaxis()->GetXmin();
      Float_t zMax = 1e3 * pData->GetPlasmaSkinDepth() * hDen2D[i]->GetXaxis()->GetXmax();
      Int_t NbinsY = hDen2D[i]->GetNbinsY();
      Float_t yMin = 1e3 * pData->GetPlasmaSkinDepth() * hDen2D[i]->GetYaxis()->GetXmin();
      Float_t yMax = 1e3 * pData->GetPlasmaSkinDepth() * hDen2D[i]->GetYaxis()->GetXmax();
      Int_t NbinsX = hDen2DXY[i]->GetNbinsX();
      Float_t xMin = 1e3 * pData->GetPlasmaSkinDepth() * hDen2DXY[i]->GetXaxis()->GetXmin();
      Float_t xMax = 1e3 * pData->GetPlasmaSkinDepth() * hDen2DXY[i]->GetXaxis()->GetXmax();

      hDen2D[i]->SetBins(NbinsZ,zMin,zMax,NbinsY,yMin,yMax);
      hDen2DXY[i]->SetBins(NbinsX,xMin,xMax,NbinsY,yMin,yMax);
      
      for(Int_t j=0;j<hDen2D[i]->GetNbinsX();j++) {
	for(Int_t k=0;k<hDen2D[i]->GetNbinsY();k++) {
	  hDen2D[i]->SetBinContent(j,k,1e-15 * 1e-6 * pData->GetPlasmaDensity() * hDen2D[i]->GetBinContent(j,k));
	}
      }

      for(Int_t j=0;j<hDen2DXY[i]->GetNbinsX();j++) {
	for(Int_t k=0;k<hDen2DXY[i]->GetNbinsY();k++) {
	  hDen2DXY[i]->SetBinContent(j,k,1e-15 * 1e-6 * pData->GetPlasmaDensity() * hDen2DXY[i]->GetBinContent(j,k));
	}
      }
      
      hDen2D[i]->GetYaxis()->SetTitle("y [mm]");      
      if(opt.Contains("comov"))
	hDen2D[i]->GetXaxis()->SetTitle("#zeta [mm]");
      else
	hDen2D[i]->GetXaxis()->SetTitle("z [mm]");
     
      if(i==0)
	hDen2D[i]->GetZaxis()->SetTitle("#LTn_{e}#GT [10^{15}/cc]");
      else
	hDen2D[i]->GetZaxis()->SetTitle("#LTn_{b}#GT [10^{15}/cc]"); 

      hDen2DXY[i]->GetXaxis()->SetTitle("x [mm]"); 
      hDen2DXY[i]->GetYaxis()->SetTitle("y [mm]");      
      
      if(i==0)
	hDen2DXY[i]->GetZaxis()->SetTitle("#LTn_{e}#GT [10^{15}/cc]");
      else
	hDen2DXY[i]->GetZaxis()->SetTitle("#LTn_{b}#GT [10^{15}/cc]"); 
      
    }
  }
  
  // Get electric fields
  const Int_t Nfields = 1;
  TH2F **hE2D = new TH2F*[Nfields];
  TH2F **hE2DXY = new TH2F*[Nfields];
  for(Int_t i=0;i<Nfields;i++) {
    hE2D[i] = NULL;
    hE2DXY[i] = NULL;

    if(!pData->GetEfieldFileName(i))
      continue;

    char hName[24];
    sprintf(hName,"hE2D_%i",i);
    hE2D[i] = (TH2F*) gROOT->FindObject(hName);
    if(hE2D[i]) delete hE2D[i];

    hE2D[i] = pData->GetEField2DSliceZY(i,-1,Nbins);
    
    hE2D[i]->SetName(hName);   
    hE2D[i]->GetXaxis()->CenterTitle();
    hE2D[i]->GetYaxis()->CenterTitle();
    hE2D[i]->GetZaxis()->CenterTitle();
    hE2D[i]->GetXaxis()->SetTitle("z [c/#omega_{p}]");
    hE2D[i]->GetYaxis()->SetTitle("y [c/#omega_{p}]");
    if(i==0)
      hE2D[i]->GetZaxis()->SetTitle("E_{z} [E_{0}]");
    else if(i==1)
      hE2D[i]->GetZaxis()->SetTitle("E_{y} [E_{0}]");
    else if(i==2)
      hE2D[i]->GetZaxis()->SetTitle("E_{x} [E_{0}]");
    
    sprintf(hName,"hEXY_%i",i);
    hE2DXY[i] = (TH2F*) gROOT->FindObject(hName);
    if(hE2DXY[i]) delete hE2D[i];

    hE2DXY[i] = pData->GetEField2DSliceXY(i,-1,Nbins);
    
    hE2DXY[i]->SetName(hName);
    hE2DXY[i]->GetXaxis()->CenterTitle();
    hE2DXY[i]->GetYaxis()->CenterTitle();
    hE2DXY[i]->GetZaxis()->CenterTitle();
    hE2DXY[i]->GetXaxis()->SetTitle("x [c/#omega_{p}]");
    hE2DXY[i]->GetYaxis()->SetTitle("y [c/#omega_{p}]");
    if(i==0)
      hE2DXY[i]->GetZaxis()->SetTitle("E_{z} [E_{0}]");
    else if(i==1)
      hE2DXY[i]->GetZaxis()->SetTitle("E_{y} [E_{0}]");
    else if(i==2)
      hE2DXY[i]->GetZaxis()->SetTitle("E_{x} [E_{0}]");
  

    { // Center the y and z coordinate
      Int_t NbinsX = hE2D[i]->GetNbinsX();
      Float_t xMin = hE2D[i]->GetXaxis()->GetXmin() - zStartPlasma;
      Float_t xMax = hE2D[i]->GetXaxis()->GetXmax() - zStartPlasma;
      Int_t NbinsY = hE2D[i]->GetNbinsY();
      Float_t yMin = hE2D[i]->GetYaxis()->GetXmin();
      Float_t yMax = hE2D[i]->GetYaxis()->GetXmax();

      if(!opt.Contains("cyl")) {
	Float_t midY = 0.5 * (hE2D[i]->GetYaxis()->GetXmin()+hE2D[i]->GetYaxis()->GetXmax());
	yMin -= midY;
	yMax -= midY;
      }

      hE2D[i]->SetBins(NbinsX,xMin,xMax,NbinsY,yMin,yMax);

      NbinsX = hE2DXY[i]->GetNbinsX();
      xMin = hE2DXY[i]->GetXaxis()->GetXmin();
      xMax = hE2DXY[i]->GetXaxis()->GetXmax();
      NbinsY = hE2DXY[i]->GetNbinsY();
      yMin = hE2DXY[i]->GetYaxis()->GetXmin();
      yMax = hE2DXY[i]->GetYaxis()->GetXmax();

      if(!opt.Contains("cyl")) {
	Float_t midY = 0.5 * (hE2DXY[i]->GetYaxis()->GetXmin()+hE2DXY[i]->GetYaxis()->GetXmax());
	yMin -= midY;
	yMax -= midY;

	Float_t midX = 0.5 * (hE2DXY[i]->GetXaxis()->GetXmin()+hE2DXY[i]->GetXaxis()->GetXmax());
	xMin -= midX;
	xMax -= midX;
      }
      
      hE2DXY[i]->SetBins(NbinsX,xMin,xMax,NbinsY,yMin,yMax);
    }

    // Change to co-moving coordinate
    if(opt.Contains("comov")) {
      Int_t NbinsX = hE2D[i]->GetNbinsX();
      Float_t xMin = hE2D[i]->GetXaxis()->GetXmin() - vbeam*Time;
      Float_t xMax = hE2D[i]->GetXaxis()->GetXmax() - vbeam*Time;
      Int_t NbinsY = hE2D[i]->GetNbinsY();
      Float_t yMin = hE2D[i]->GetYaxis()->GetXmin();
      Float_t yMax = hE2D[i]->GetYaxis()->GetXmax();
      hE2D[i]->SetBins(NbinsX,xMin,xMax,NbinsY,yMin,yMax);
      hE2D[i]->GetXaxis()->SetTitle("#zeta [c/#omega_{p}]");
    }
  
    // Chaning to user units: 
    if(opt.Contains("units") && pData->GetPlasmaDensity()) {
      Int_t NbinsZ = hE2D[i]->GetNbinsX();
      Float_t zMin = 1e3 * pData->GetPlasmaSkinDepth() * hE2D[i]->GetXaxis()->GetXmin();
      Float_t zMax = 1e3 * pData->GetPlasmaSkinDepth() * hE2D[i]->GetXaxis()->GetXmax();
      Int_t NbinsY = hE2D[i]->GetNbinsY();
      Float_t yMin = 1e3 * pData->GetPlasmaSkinDepth() * hE2D[i]->GetYaxis()->GetXmin();
      Float_t yMax = 1e3 * pData->GetPlasmaSkinDepth() * hE2D[i]->GetYaxis()->GetXmax();
      Int_t NbinsX = hE2DXY[i]->GetNbinsX();
      Float_t xMin = 1e3 * pData->GetPlasmaSkinDepth() * hE2DXY[i]->GetXaxis()->GetXmin();
      Float_t xMax = 1e3 * pData->GetPlasmaSkinDepth() * hE2DXY[i]->GetXaxis()->GetXmax();

      hE2D[i]->SetBins(NbinsZ,zMin,zMax,NbinsY,yMin,yMax);
      hE2DXY[i]->SetBins(NbinsX,xMin,xMax,NbinsY,yMin,yMax);
      
      
      for(Int_t j=0;j<hE2D[i]->GetNbinsX();j++) {
	for(Int_t k=0;k<hE2D[i]->GetNbinsY();k++) {
	  hE2D[i]->SetBinContent(j,k,1e-9 * pData->GetPlasmaE0() * hE2D[i]->GetBinContent(j,k));
	}
      }

      for(Int_t j=0;j<hE2DXY[i]->GetNbinsX();j++) {
	for(Int_t k=0;k<hE2DXY[i]->GetNbinsY();k++) {
	  hE2DXY[i]->SetBinContent(j,k,1e-15 * 1e-6 * pData->GetPlasmaDensity() * hE2DXY[i]->GetBinContent(j,k));
	}
      }

      hE2D[i]->GetYaxis()->SetTitle("y [mm]");	
      if(opt.Contains("comov"))
	hE2D[i]->GetXaxis()->SetTitle("#zeta [mm]");
      else
	hE2D[i]->GetXaxis()->SetTitle("z [mm]");

      if(i==0)
	hE2D[i]->GetZaxis()->SetTitle("E_{z} [GeV/m]");
      else if(i==1)
	hE2D[i]->GetZaxis()->SetTitle("E_{y} [GeV/m}]");
      else if(i==2)
	hE2D[i]->GetZaxis()->SetTitle("E_{x} [GeV/m}]");

      hE2DXY[i]->GetYaxis()->SetTitle("y [mm]");	
      hE2DXY[i]->GetXaxis()->SetTitle("x [mm]");

      if(i==0)
	hE2DXY[i]->GetZaxis()->SetTitle("E_{z} [GeV/m]");
      else if(i==1)
	hE2DXY[i]->GetZaxis()->SetTitle("E_{y} [GeV/m}]");
      else if(i==2)
	hE2DXY[i]->GetZaxis()->SetTitle("E_{x} [GeV/m}]");
	
    }
    
  }
  
  // Tunning the Histograms
  // ---------------------
  
  // Set the range of the histogram for maximum constrast
  Float_t density = 1;
  if(opt.Contains("units") && pData->GetPlasmaDensity())
    density = 1e-15 * 1e-6 * pData->GetPlasmaDensity();
  
  Float_t Max  = 1.1 * hDen2D[0]->GetMaximum();
  Float_t Base = density;
  Float_t Min  = 2.* Base - Max;
  if(Max >= 2. * Base) {
    Min = 0;
  } else if(Max<1.0 * Base) {
    Max = 1.1 * Base;
    Min = 0.;
  }
  
  hDen2D[0]->GetZaxis()->SetRangeUser(Min,Max);  
  hDen2DXY[0]->GetZaxis()->SetRangeUser(Min,Max);  

  // Zoom
  Float_t range    = (hDen2D[0]->GetYaxis()->GetXmax() - hDen2D[0]->GetYaxis()->GetXmin())/zoom;
  Float_t midPoint = (hDen2D[0]->GetYaxis()->GetXmax() + hDen2D[0]->GetYaxis()->GetXmin())/2.;
  
  hDen2D[0]->GetYaxis()->SetRangeUser(midPoint-range/2,midPoint+range/2);  
  hE2D[0]->GetYaxis()->SetRangeUser(midPoint-range/2,midPoint+range/2);  

  range    = (hDen2DXY[0]->GetYaxis()->GetXmax() - hDen2DXY[0]->GetYaxis()->GetXmin())/zoom;
  midPoint = (hDen2DXY[0]->GetYaxis()->GetXmax() + hDen2DXY[0]->GetYaxis()->GetXmin())/2.;
  hDen2DXY[0]->GetYaxis()->SetRangeUser(midPoint-range/2,midPoint+range/2);  
  hE2DXY[0]->GetYaxis()->SetRangeUser(midPoint-range/2,midPoint+range/2);  

  range    = (hDen2DXY[0]->GetXaxis()->GetXmax() - hDen2DXY[0]->GetXaxis()->GetXmin())/zoom;
  midPoint = (hDen2DXY[0]->GetXaxis()->GetXmax() + hDen2DXY[0]->GetXaxis()->GetXmin())/2.;
  hDen2DXY[0]->GetXaxis()->SetRangeUser(midPoint-range/2,midPoint+range/2);  
  hE2DXY[0]->GetXaxis()->SetRangeUser(midPoint-range/2,midPoint+range/2);  

  Float_t MinE = 1.1*hE2D[0]->GetMinimum();
  Float_t MaxE = 1.1*hE2D[0]->GetMaximum();
  hE2D[0]->GetZaxis()->SetRangeUser(MinE,MaxE);   
  hE2DXY[0]->GetZaxis()->SetRangeUser(MinE,MaxE);  
  
  // Plotting
  // -----------------------------------------------

  // Canvas setup
  TCanvas *C = new TCanvas("C","3D Charge density and Electric field",1000,700);

  // Palettes setup
  TExec *exPlasma;
  TExec *exElec;
  if( (Max >= 2. * Base) || (Max<1.0 * Base) ) {
    exPlasma= new TExec("exPlasma","plasmaHPalette->cd();"); 
    exElec   = new TExec("exElec","electronPalette->cd();");
  } else {
    exPlasma= new TExec("exPlasma","plasmaPalette->cd();");
    exElec   = new TExec("exElec","electronPalette->cd();");
  }
  TExec *exField  = new TExec("exField","rbowPalette->cd();");

  // Text objects
  TPaveText *textTime = new TPaveText(0.78,0.85,0.95,0.92,"NDC");
  PlasmaGlob::SetPaveTextStyle(textTime); 
  char ctext[128];
  if(opt.Contains("units") && pData->GetPlasmaDensity()) 
    sprintf(ctext,"Z = %5.1f mm", 1e3 * pData->GetPlasmaSkinDepth() * Time);
  else
    sprintf(ctext,"T = %5.1f 1/#omega_{p}",Time);
  textTime->AddText(ctext);
 
  TPaveText *textDen = new TPaveText(0.15,0.85,0.40,0.92,"NDC");
  PlasmaGlob::SetPaveTextStyle(textDen,12); 
  textDen->SetTextColor(kOrange+10);
  if(opt.Contains("units") && pData->GetPlasmaDensity())
    sprintf(ctext,"n_{0} = %5.2f x 10^{15} / cc", 1e-6 * 1e-15 * pData->GetPlasmaDensity());
  else if(pData->GetBeamDensity() && pData->GetPlasmaDensity())
    sprintf(ctext,"n_{b}/n_{0} = %5.2f", pData->GetBeamDensity()/pData->GetPlasmaDensity());
  
  textDen->AddText(ctext);
  
  // Actual Plotting!
  // ------------------------------------------------------------

  // Output file
  TString fOutName = Form("./%s/Plots/ChargeField3D/ChargeField3D",sim.Data());
  fOutName += Form("-%s_%i",sim.Data(),time);

  // Setup Pad layout:
  Double_t lMargin = 0.08;
  Double_t rMargin = 0.30;
  Double_t bMargin = 0.10;
  Double_t tMargin = 0.02;
  Double_t vSpacing = 0.01; 
  Double_t hSpacing = 0.01; 
  Double_t hStepL = (1.-lMargin-rMargin-hSpacing) * 0.90;
  Double_t hStepR = (1.-lMargin-rMargin-hSpacing) * 0.10;
  Double_t vStep = (1.-bMargin-tMargin-vSpacing)/2.;
  
  TPad *pad[4];

  // top plots
  pad[0] = new TPad("padt", "padt",  0.00, bMargin + vStep + vSpacing,
		    lMargin + hStepL, 1.00);
  pad[0]->SetLeftMargin(1./(lMargin+hStepL)*lMargin);
  pad[0]->SetRightMargin(0.0);  
  pad[0]->SetBottomMargin(0.0);                                   
  pad[0]->SetTopMargin(1./(tMargin+vStep)*tMargin);
  pad[0]->Draw();

  pad[1] = new TPad("padtr", "padtr",lMargin + hStepL + hSpacing, bMargin + vStep + vSpacing,
		    1.00, 1.00);
  pad[1]->SetLeftMargin(0.);
  pad[1]->SetRightMargin(1./(1 - lMargin + hStepL + hSpacing) * rMargin);  
  pad[1]->SetBottomMargin(0.0);                                   
  pad[1]->SetTopMargin(1./(tMargin+vStep)*tMargin);
  pad[1]->Draw();

  // bottom plots
  pad[2] = new TPad("padb", "padb", 0.00, 0.00,
		    lMargin + hStepL, bMargin + vStep);
  pad[2]->SetLeftMargin(1./(lMargin+hStepL)*lMargin);
  pad[2]->SetRightMargin(0.0);
  pad[2]->SetBottomMargin(1./(bMargin+vStep)*bMargin);
  pad[2]->SetTopMargin(0.);
  pad[2]->Draw();  

  pad[3] = new TPad("padbr", "padbr",lMargin + hStepL + hSpacing, 0.00,
		    1.0, bMargin + vStep);
  pad[3]->SetLeftMargin(0.0);
  pad[3]->SetRightMargin(1./(1-lMargin + hStepL + hSpacing) * rMargin);
  pad[3]->SetBottomMargin(1./(bMargin+vStep)*bMargin);
  pad[3]->SetTopMargin(0.);
  pad[3]->Draw();  

  Float_t yFactor = pad[0]->GetAbsHNDC()/pad[1]->GetAbsHNDC();
  Float_t xFactor = pad[0]->GetAbsWNDC()/pad[1]->GetAbsWNDC();


  pad[0]->cd(); // <--- Top Plot
  pad[0]->SetFrameLineWidth(3);  

  TH2F *hFrame = (TH2F*) gROOT->FindObject("hFrame0");
  if(hFrame) delete hFrame;
  hFrame = (TH2F*) hDen2D[0]->Clone("hFrame0");
  hFrame->Reset();
  hFrame->GetYaxis()->SetNdivisions(505);
  hFrame->GetYaxis()->SetLabelSize(0.075);
  hFrame->GetYaxis()->SetTitleSize(0.075);
  hFrame->GetYaxis()->SetTitleOffset(0.5);
  hFrame->GetXaxis()->SetLabelOffset(0.10);
  hFrame->Draw("col");
  
  
  hDen2D[0]->GetZaxis()->SetLabelSize(0.04);  
  hDen2D[0]->GetZaxis()->SetTitleSize(0.04);
  //hDen2D[0]->GetZaxis()->SetTitleOffset(0.06);
  hDen2D[0]->GetZaxis()->SetNdivisions(505);

  exPlasma->Draw();
  hDen2D[0]->Draw("col same");
  
  if(hDen2D[1]) {
    hDen2D[1]->GetZaxis()->SetLabelSize(0.04);  
    hDen2D[1]->GetZaxis()->SetTitleSize(0.04);
    hDen2D[1]->GetZaxis()->SetNdivisions(505);

    exElec->Draw();
    hDen2D[1]->Draw("col same");
  }
  
  textTime->Draw();
  textDen->Draw();
  pad[0]->RedrawAxis(); 


  pad[1]->cd(); // <--- Top Plot
  pad[1]->SetFrameLineWidth(3);  

  hFrame = (TH2F*) gROOT->FindObject("hFrame1");
  if(hFrame) delete hFrame;
  hFrame = (TH2F*) hDen2DXY[0]->Clone("hFrame1");
  hFrame->Reset();
  hFrame->GetYaxis()->SetNdivisions(505);
  hFrame->GetYaxis()->SetTitleOffset(-99);
  hFrame->GetYaxis()->SetLabelOffset(-99);
  hFrame->GetXaxis()->SetLabelOffset(0.10);
  hFrame->Draw("col");
  
  hDen2DXY[0]->GetZaxis()->SetLabelSize(0.05);  
  hDen2DXY[0]->GetZaxis()->SetTitleSize(0.05);
  hDen2DXY[0]->GetZaxis()->SetTitleOffset(1.3);
  hDen2DXY[0]->GetZaxis()->SetNdivisions(505);

  exPlasma->Draw();
  hDen2DXY[0]->Draw("colz same");
  
  if(hDen2DXY[1]) {
    hDen2DXY[1]->GetZaxis()->SetLabelSize(0.05);  
    hDen2DXY[1]->GetZaxis()->SetTitleSize(0.05);
    hDen2DXY[1]->GetZaxis()->SetTitleOffset(1.3);
    hDen2DXY[1]->GetZaxis()->SetNdivisions(505);

    exElec->Draw();
    hDen2DXY[1]->Draw("colz same");
  }
   
  pad[1]->Update();

  TPaletteAxis *palette = (TPaletteAxis*)hDen2DXY[1]->GetListOfFunctions()->FindObject("palette");
  
  Float_t y1 = pad[1]->GetBottomMargin();
  Float_t y2 = 1 - pad[1]->GetTopMargin();
  palette->SetY2NDC(y2);
  palette->SetY1NDC(0.5*(y1+y2));
  //palette->SetTitleOffset(0.8);
  palette->SetBorderSize(2);
  palette->SetLineColor(1);
  
  TPaletteAxis *palette2 = (TPaletteAxis*)hDen2DXY[0]->GetListOfFunctions()->FindObject("palette");
  palette2->SetY2NDC(0.5*(y1+y2));
  palette2->SetY1NDC(y1);

  pad[1]->RedrawAxis(); 
 
  pad[2]->cd(); // <--- Mid Plot
  pad[2]->SetFrameLineWidth(3);  
  
  hFrame = (TH2F*) gROOT->FindObject("hFrame2");
  if(hFrame) delete hFrame;
  hFrame = (TH2F*) hE2D[0]->Clone("hFrame2");
  hFrame->Reset();
  hFrame->GetYaxis()->SetNdivisions(505);
  hFrame->GetYaxis()->SetLabelSize(yFactor*0.075);
  hFrame->GetYaxis()->SetTitleSize(yFactor*0.075);
  hFrame->GetYaxis()->SetTitleOffset(0.5/yFactor);
  
  hFrame->GetXaxis()->SetTitleSize(0.07);
  hFrame->GetXaxis()->SetLabelSize(0.06);
  hFrame->GetXaxis()->SetTitleOffset(1.2);
 
  hFrame->Draw("col");

  hE2D[0]->GetZaxis()->SetLabelSize(yFactor*0.06);
  hE2D[0]->GetZaxis()->SetTitleSize(yFactor*0.06);
  hE2D[0]->GetZaxis()->SetTitleOffset(yFactor*1.4);
  hE2D[0]->GetZaxis()->SetNdivisions(510);

  exField->Draw();
  hE2D[0]->Draw("col same");

  pad[2]->RedrawAxis(); 

  pad[3]->cd(); // <--- Mid Plot
  pad[3]->SetFrameLineWidth(3);  
  
  hFrame = (TH2F*) gROOT->FindObject("hFrame3");
  if(hFrame) delete hFrame;
  hFrame = (TH2F*) hE2DXY[0]->Clone("hFrame3");
  hFrame->Reset();
  hFrame->GetYaxis()->SetNdivisions(505);
  hFrame->GetYaxis()->SetTitleOffset(-99);
  hFrame->GetYaxis()->SetLabelOffset(-99);
  
  hFrame->GetXaxis()->SetTitleSize(0.07);
  hFrame->GetXaxis()->SetLabelSize(0.06);
  hFrame->GetXaxis()->SetTitleOffset(1.2);
 
  hFrame->Draw("col");

  hE2DXY[0]->GetZaxis()->SetLabelSize(yFactor*0.06);
  hE2DXY[0]->GetZaxis()->SetTitleSize(yFactor*0.06);
  hE2DXY[0]->GetZaxis()->SetTitleOffset(yFactor*1.4);
  hE2DXY[0]->GetZaxis()->SetNdivisions(510);

  exField->Draw();
  hE2DXY[0]->Draw("colz same");

  pad[3]->RedrawAxis(); 

  C->cd();

  // Print to a file
  PlasmaGlob::imgconv(C,fOutName,opt);
  // ---------------------------------------------------------

}
