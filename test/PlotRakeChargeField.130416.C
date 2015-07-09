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

void PlotRakeChargeField( const TString &sim, Int_t time, Int_t zoom=2, Int_t Nbins=2, const TString &options="") {
  
#ifdef __CINT__  
  gSystem->Load("libplasma.so");
#endif

  PlasmaGlob::Initialize();

  // Palettes!
  gROOT->Macro("PlasmaPalettes.C");

  // Init Units table
  PUnits::UnitsTable::Get();
  
  // Load PData
  PData *pData = PData::Get(sim.Data());
  pData->LoadFileNames(time);
  if(!pData->IsInit()) return;

  TString opt = options;
 
  // More makeup
  gStyle->SetPadGridY(0);
  if(opt.Contains("grid")) {
    gStyle->SetPadGridX(1);
    gStyle->SetPadGridY(1);
  }
  gStyle->SetTitleFont(42);
 
  

  Bool_t CYL = kFALSE;
  if(sim.Contains("cyl")) { CYL = kTRUE; opt += "cyl"; } 
  
  Bool_t ThreeD = kFALSE;
  if(sim.Contains("3D")) ThreeD = kTRUE; 
  
  // Some plasma constants
  Double_t n0 = pData->GetPlasmaDensity();
  Double_t kp = pData->GetPlasmaK();
  Double_t skindepth = 1.;
  if(kp!=0.0) skindepth = 1/kp;
  Double_t E0 = pData->GetPlasmaE0();
  
  // Some beam properties:
  Double_t Ebeam = pData->GetBeamEnergy();
  Double_t gamma = pData->GetBeamGamma();
  Double_t vbeam = pData->GetBeamVelocity();

  cout << Form(" - Bunch gamma      = %8.4f", gamma ) << endl;
  cout << Form(" - Bunch velocity   = %8.4f c", vbeam ) << endl;
  
  // Other parameters
  Float_t trapPotential = 1.0 - (1.0/gamma);
  cout << Form(" - Trap. potential  = %8.4f mc2/e",trapPotential) << endl;
  
  // Time in OU
  Float_t Time = pData->GetRealTime();
  // z start of the plasma in normalized units.
  Float_t zStartPlasma = pData->GetPlasmaStart()*kp;
  // z start of the beam in normalized units.
  Float_t zStartBeam = pData->GetBeamStart()*kp;
  
  if(opt.Contains("center")) {
    Time -= zStartPlasma;
    if(opt.Contains("comov"))      // Centers on the head of the beam.
      Time += zStartBeam;
  }
  
  
  // Get charge density histos
  Int_t Nspecies = pData->NSpecies();
  TH2F **hDen2D = new TH2F*[Nspecies];
  for(Int_t i=0;i<Nspecies;i++) {
    hDen2D[i] = NULL;
    
    if(!pData->GetChargeFileName(i)) 
      continue;
    
    char hName[24];
    sprintf(hName,"hDen_%i",i);
    hDen2D[i] = (TH2F*) gROOT->FindObject(hName);
    if(hDen2D[i]) delete hDen2D[i];
    
    if(!ThreeD)
      hDen2D[i] = pData->GetCharge(i,opt);
    else
      hDen2D[i] = pData->GetCharge2DSliceZY(i,-1,Nbins,opt+"avg");
    
    hDen2D[i]->SetName(hName);
    hDen2D[i]->GetXaxis()->CenterTitle();
    hDen2D[i]->GetYaxis()->CenterTitle();
    hDen2D[i]->GetZaxis()->CenterTitle();
    
    if(opt.Contains("comov"))
      hDen2D[i]->GetXaxis()->SetTitle("k_{p} #zeta");
    else
      hDen2D[i]->GetXaxis()->SetTitle("k_{p} z");
    
    if(CYL) 
      hDen2D[i]->GetYaxis()->SetTitle("k_{p} r");
    else
      hDen2D[i]->GetYaxis()->SetTitle("k_{p} y");
    
    hDen2D[i]->GetZaxis()->SetTitle("n [n_{0}]");
  }
  
  
  // Slice width limits.
  Int_t FirstyBin = 0;
  Int_t LastyBin = 0;
  if(!CYL) {
    FirstyBin = hDen2D[0]->GetNbinsY()/2 + 1 - Nbins;
    LastyBin =  hDen2D[0]->GetNbinsY()/2 + Nbins;
  } else {
    FirstyBin = 1; 
    LastyBin  = Nbins;
  }
  
  // Get charge density on-axis
  TH1F **hDen1D = new TH1F*[Nspecies];
  for(Int_t i=0;i<Nspecies;i++) {
    hDen1D[i] = NULL;
    
    if(!pData->GetEfieldFileName(i))
      continue;
    
    char hName[24];
    sprintf(hName,"hDen1D_%i",i);
    hDen1D[i] = (TH1F*) gROOT->FindObject(hName);
    if(hDen1D[i]) delete hDen1D[i];
      
    // 1D histograms
    if(ThreeD) {
      hDen1D[i] = pData->GetH1SliceZ3D(pData->GetChargeFileName(i)->c_str(),"charge",-1,Nbins,-1,Nbins,opt+"avg");
    } else if(CYL) { // Cylindrical: The first bin with r>0 is actually the number 1 (not the 0).
      hDen1D[i] = pData->GetH1SliceZ(pData->GetChargeFileName(i)->c_str(),"charge",1,Nbins,opt+"avg");
    } else { // 2D cartesian
      hDen1D[i] = pData->GetH1SliceZ(pData->GetChargeFileName(i)->c_str(),"charge",-1,Nbins,opt+"avg");
    }
    hDen1D[i]->SetName(hName); 
   
    // if(hDen1D[i]) delete hDen1D[i];
    // hDen1D[i] = (TH1F*) hE2D[i]->ProjectionX(hName,FirstyBin,LastyBin);
    // hDen1D[i]->Scale(1.0/(LastyBin-FirstyBin+1));
    
    if(opt.Contains("comov"))
      hDen1D[i]->GetXaxis()->SetTitle("#zeta [c/#omega_{p}]");
    else
      hDen1D[i]->GetXaxis()->SetTitle("z [c/#omega_{p}]");

    if(i==0)
      hDen1D[i]->GetYaxis()->SetTitle("n/n_{0}");
    else if(i==1)
      hDen1D[i]->GetYaxis()->SetTitle("n_{b}/n_{0}");
    else   
      hDen1D[i]->GetYaxis()->SetTitle("n_{i}/n_{0}");
  }


  // Get electric fields 2D
  const Int_t Nfields = 1;
  TH2F **hE2D = new TH2F*[Nfields];
  for(Int_t i=0;i<Nfields;i++) {
    hE2D[i] = NULL;
    
    if(!pData->GetEfieldFileName(i))
      continue;
    
    char hName[24];
    sprintf(hName,"hE2D_%i",i);
    hE2D[i] = (TH2F*) gROOT->FindObject(hName);
    if(hE2D[i]) delete hE2D[i];
    
    if(!ThreeD)
      hE2D[i] = pData->GetEField(i,opt);
    else
      hE2D[i] = pData->GetEField2DSliceZY(i,-1,Nbins,opt+"avg");
    
    hE2D[i]->SetName(hName);   
    hE2D[i]->GetXaxis()->CenterTitle();
    hE2D[i]->GetYaxis()->CenterTitle();
    hE2D[i]->GetZaxis()->CenterTitle();
    if(opt.Contains("comov"))
      hE2D[i]->GetXaxis()->SetTitle("k_{p} #zeta");
    else
      hE2D[i]->GetXaxis()->SetTitle("k_{p} z");
    
    if(CYL) 
      hE2D[i]->GetYaxis()->SetTitle("k_{p} r");
    else
      hE2D[i]->GetYaxis()->SetTitle("k_{p} y");
    
    if(i==0)
      hE2D[i]->GetZaxis()->SetTitle("E_{z}/E_{0}");
    else if(i==1)
      hE2D[i]->GetZaxis()->SetTitle("E_{y}/E_{0}");
    else if(i==2)
      hE2D[i]->GetZaxis()->SetTitle("E_{x}/E_{0}");
    
  }
  
  // Get electric fields
  TH1F **hE1D = new TH1F*[Nfields];
  TH1F **hV1D = new TH1F*[Nfields];
  for(Int_t i=0;i<Nfields;i++) {
    hE1D[i] = NULL;
    
    if(!pData->GetEfieldFileName(i))
      continue;
    
    char hName[24];
    sprintf(hName,"hE1D_%i",i);
    hE1D[i] = (TH1F*) gROOT->FindObject(hName);
    if(hE1D[i]) delete hE1D[i];
      
    // 1D histograms
    TString opth1 = opt;
    opth1 += "avg";
    
    char nam[3]; sprintf(nam,"e%i",i+1);
    if(ThreeD) {
      
      if(i==0) 
	hE1D[i] = pData->GetH1SliceZ3D(pData->GetEfieldFileName(i)->c_str(),nam,-1,Nbins,-1,Nbins,opt+"avg");
      else  
	hE1D[i] = pData->GetH1SliceZ3D(pData->GetEfieldFileName(i)->c_str(),nam,-Nbins,Nbins,-Nbins,Nbins,opt+"avg");
   
    } else if(CYL) { // Cylindrical: The first bin with r>0 is actually the number 1 (not the 0).
    
	hE1D[i] = pData->GetH1SliceZ(pData->GetEfieldFileName(i)->c_str(),nam,1,Nbins,opt+"avg");
    
    } else { // 2D cartesian
    
      if(i==0) 
	hE1D[i] = pData->GetH1SliceZ(pData->GetEfieldFileName(i)->c_str(),nam,-1,Nbins,opt+"avg");
      else 
	hE1D[i] = pData->GetH1SliceZ(pData->GetEfieldFileName(i)->c_str(),nam,-Nbins,Nbins,opt+"avg");
    
    }
    
    // Alternative
    // if(hE1D[i]) delete hE1D[i];
    // hE1D[i] = (TH1F*) hE2D[i]->ProjectionX(hName,FirstyBin,LastyBin);
    // hE1D[i]->Scale(1.0/(LastyBin-FirstyBin+1));
    
    Int_t   Nbins = hE1D[i]->GetNbinsX();
    Float_t dz = hE1D[i]->GetBinWidth(Nbins);
    sprintf(hName,"hV1D_%i",i);
    hV1D[i] = (TH1F*) hE1D[i]->Clone(hName);
    hV1D[i]->Reset();
    Double_t integral = 0.0;
    for(Int_t k=Nbins;k>0;k--) {
      integral += - hE1D[i]->GetBinContent(k) * dz;
      hV1D[i]->SetBinContent(k,integral);
    }
  
    hE1D[i]->SetName(hName);
    if(opt.Contains("comov"))
      hE1D[i]->GetXaxis()->SetTitle("#zeta [c/#omega_{p}]");
    else
      hE1D[i]->GetXaxis()->SetTitle("z [c/#omega_{p}]");
    
    if(i==0)
      hE1D[i]->GetYaxis()->SetTitle("E_{z} [E_{0}]");
    else if(i==1)
      hE1D[i]->GetYaxis()->SetTitle("E_{y} [E_{0}]");
    else if(i==2)
      hE1D[i]->GetYaxis()->SetTitle("E_{x} [E_{0}]");
  }

  // Tunning the Histograms
  // ---------------------
  
  // Chaning to user units: 
  // --------------------------
  
  if(opt.Contains("units") && n0) {
    
    for(Int_t i=0;i<Nspecies;i++) {
      Int_t NbinsX = hDen2D[i]->GetNbinsX();
      Float_t xMin = skindepth * hDen2D[i]->GetXaxis()->GetXmin() / PUnits::um;
      Float_t xMax = skindepth * hDen2D[i]->GetXaxis()->GetXmax() / PUnits::um;
      Int_t NbinsY = hDen2D[i]->GetNbinsY();
      Float_t yMin = skindepth * hDen2D[i]->GetYaxis()->GetXmin() / PUnits::um;
      Float_t yMax = skindepth * hDen2D[i]->GetYaxis()->GetXmax() / PUnits::um;
      hDen2D[i]->SetBins(NbinsX,xMin,xMax,NbinsY,yMin,yMax);
      // for(Int_t j=0;j<hDen2D[i]->GetNbinsX();j++) {
      // 	for(Int_t k=0;k<hDen2D[i]->GetNbinsY();k++) {
      // 	  hDen2D[i]->SetBinContent(j,k, hDen2D[i]->GetBinContent(j,k) * n0 / (1e17/PUnits::cm3) );
      // 	}
      // }

      if(CYL)
	hDen2D[i]->GetYaxis()->SetTitle("r [#mum]");      
      else
	hDen2D[i]->GetYaxis()->SetTitle("y [#mum]");      

      if(opt.Contains("comov"))
	hDen2D[i]->GetXaxis()->SetTitle("#zeta [#mum]");
      else
	hDen2D[i]->GetXaxis()->SetTitle("z [#mum]");
      
      // if(i==0)
      // 	hDen2D[i]->GetZaxis()->SetTitle("n_{e} [10^{17}/cm^{3}]");
      // else if(i==1)
      // 	hDen2D[i]->GetZaxis()->SetTitle("n_{b} [10^{17}/cm^{3}]"); 
      // else
      // 	hDen2D[i]->GetZaxis()->SetTitle("n_{i} [10^{17}/cm^{3}]"); 
    }

    for(Int_t i=0;i<Nspecies;i++) {
      Int_t NbinsX = hDen1D[i]->GetNbinsX();
      Float_t xMin = skindepth * hDen1D[i]->GetXaxis()->GetXmin() / PUnits::um;
      Float_t xMax = skindepth * hDen1D[i]->GetXaxis()->GetXmax() / PUnits::um;
      hDen1D[i]->SetBins(NbinsX,xMin,xMax);
      // for(Int_t j=0;j<hDen1D[i]->GetNbinsX();j++) {
      // 	hDen1D[i]->SetBinContent(j, hDen1D[i]->GetBinContent(j) * n0 / (1e17/PUnits::cm3) );
      // }
      

      if(opt.Contains("comov"))
    	hDen1D[i]->GetXaxis()->SetTitle("#zeta [#mum]");
      else
    	hDen1D[i]->GetXaxis()->SetTitle("z [#mum]");
      
      if(i==0)
    	hDen1D[i]->GetYaxis()->SetTitle("n_{e} [10^{17}/cm^{3}]");
      else if(i==1)
    	hDen1D[i]->GetYaxis()->SetTitle("n_{b} [10^{17}/cm^{3}]"); 
      else
 	hDen1D[i]->GetYaxis()->SetTitle("n_{i} [10^{17}/cm^{3}]"); 	
    }
    
    
    for(Int_t i=0;i<Nfields;i++) {
      Int_t NbinsX = hE2D[i]->GetNbinsX();
      Float_t xMin = skindepth * hE2D[i]->GetXaxis()->GetXmin() / PUnits::um;
      Float_t xMax = skindepth * hE2D[i]->GetXaxis()->GetXmax() / PUnits::um;
      Int_t NbinsY = hE2D[i]->GetNbinsY();
      Float_t yMin = skindepth * hE2D[i]->GetYaxis()->GetXmin() / PUnits::um;
      Float_t yMax = skindepth * hE2D[i]->GetYaxis()->GetXmax() / PUnits::um;
      hE2D[i]->SetBins(NbinsX,xMin,xMax,NbinsY,yMin,yMax);
      for(Int_t j=0;j<hE2D[i]->GetNbinsX();j++) {
	for(Int_t k=0;k<hE2D[i]->GetNbinsY();k++) {
	  hE2D[i]->SetBinContent(j,k, hE2D[i]->GetBinContent(j,k) * ( E0 / (PUnits::GV/PUnits::m) ) );
	}
      }
      
      if(CYL)
	hE2D[i]->GetYaxis()->SetTitle("r [#mum]");      
      else
	hE2D[i]->GetYaxis()->SetTitle("y [#mum]");      
 
      if(opt.Contains("comov"))
	hE2D[i]->GetXaxis()->SetTitle("#zeta [#mum]");
      else
	hE2D[i]->GetXaxis()->SetTitle("z [#mum]");
      
      if(i==0)
	hE2D[i]->GetZaxis()->SetTitle("E_{z} [GV/m]");
      else if(i==1)
	hE2D[i]->GetZaxis()->SetTitle("E_{y} [GV/m]");
      else if(i==2)
	hE2D[i]->GetZaxis()->SetTitle("E_{x} [GV/m]");
    }

    for(Int_t i=0;i<Nfields;i++) {
      Int_t NbinsX = hE1D[i]->GetNbinsX();
      Float_t xMin = skindepth * hE1D[i]->GetXaxis()->GetXmin() / PUnits::um;
      Float_t xMax = skindepth * hE1D[i]->GetXaxis()->GetXmax() / PUnits::um;

      hE1D[i]->SetBins(NbinsX,xMin,xMax);
      hV1D[i]->SetBins(NbinsX,xMin,xMax);
      
      for(Int_t j=0;j<NbinsX;j++) {
	hE1D[i]->SetBinContent(j, hE1D[i]->GetBinContent(j) * ( E0 / (PUnits::GV) ) );
	hV1D[i]->SetBinContent(j, hV1D[i]->GetBinContent(j) * ( E0 / (PUnits::GV) ) );
      }
      
      if(opt.Contains("comov"))
	hE1D[i]->GetXaxis()->SetTitle("#zeta [mm]");
      else
	hE1D[i]->GetXaxis()->SetTitle("z [mm]");
      
      if(i==0)
	hE1D[i]->GetYaxis()->SetTitle("E_{z} [GV/m]");
      else if(i==1)
	hE1D[i]->GetYaxis()->SetTitle("E_{y} [GV/m]");
      else if(i==2)
	hE1D[i]->GetYaxis()->SetTitle("E_{x} [GV/m]");
    }

    
  }

  
  
  // Zoom
  Float_t range    = (hDen2D[0]->GetYaxis()->GetXmax() - hDen2D[0]->GetYaxis()->GetXmin())/zoom;
  Float_t midPoint = (hDen2D[0]->GetYaxis()->GetXmax() + hDen2D[0]->GetYaxis()->GetXmin())/2.;
  
  Float_t yMin,yMax;
  if(!CYL) {
    yMin = midPoint-range/2;
    yMax = midPoint+range/2;
  } else {
    yMin = 0.;
    yMax = range;
  }
  hDen2D[0]->GetYaxis()->SetRangeUser(yMin,yMax);
  hE2D[0]->GetYaxis()->SetRangeUser(yMin,yMax);  
  
  // ------------- z Zoom --------------------------------- Plasma palette -----------
  // Set the range of the plasma charge density histogram for maximum constrast 
  // using a dynamic palette wich adjust the nominal value to a certain color.
  
  Float_t density = 1;
  // if(opt.Contains("units") && n0)
  //   density = n0 / (1e17/PUnits::cm3);
  Float_t Base = density;

  Float_t Max  = hDen2D[1]->GetMaximum();
  Float_t Min  = (0.1001) * density;
  
  hDen2D[0]->GetZaxis()->SetRangeUser(Min,Max); 

  if(hDen2D[1]) {
  
    Float_t Minb = (0.1001) * density;
    Float_t Maxb = hDen2D[1]->GetMaximum();
  
    // hDen2D[1]->GetZaxis()->SetRangeUser(Minb,Maxb);
    hDen2D[1]->GetZaxis()->SetRangeUser(Min,Max);

  } 

  if(Nspecies>=3) {
    if(hDen2D[2]) {
      
      Float_t Minb = (0.1001) * density;
      Float_t Maxb = hDen2D[2]->GetMaximum();
      
      // hDen2D[2]->GetZaxis()->SetRangeUser(Minb,Maxb);
      hDen2D[2]->GetZaxis()->SetRangeUser(Min,Max);
      
    } 
  }
  
  // Dynamic plasma palette
  const Int_t plasmaDNRGBs = 3;
  const Int_t plasmaDNCont = 128;
  Double_t basePos = 0.5;
  if(Max!=Min) {
    if(opt.Contains("logz")) {
      Float_t a = 1.0/(TMath::Log10(Max)-TMath::Log10(Min));
      Float_t b = TMath::Log10(Min);
      basePos = a*(TMath::Log10(Base) - b);
      //  cout << Form("Min = %f Base = %f Max = %f",a*(TMath::Log10(Min) - b),basePos,a*(TMath::Log10(Max) - b)) << endl;     
    } else {
      basePos = (1.0/(Max-Min))*(Base - Min);
    }
  }
  
  Double_t plasmaDStops[plasmaDNRGBs] = { 0.00, basePos, 1.00 };
  Double_t plasmaDRed[plasmaDNRGBs]   = { 0.99, 0.90, 0.00 };
  Double_t plasmaDGreen[plasmaDNRGBs] = { 0.99, 0.90, 0.00 };
  Double_t plasmaDBlue[plasmaDNRGBs]  = { 0.99, 0.90, 0.00 };

  PPalette * plasmaPalette = (PPalette*) gROOT->FindObject("plasma");
  plasmaPalette->CreateGradientColorTable(plasmaDNRGBs, plasmaDStops, 
					   plasmaDRed, plasmaDGreen, plasmaDBlue, plasmaDNCont);
  
  // Change the range of z axis for the fields to be symmetric.
  Float_t Emax = hE2D[0]->GetMaximum();
  Float_t Emin = hE2D[0]->GetMinimum();
  if(Emax > TMath::Abs(Emin))
    Emin = -Emax;
  else
    Emax = -Emin;
  hE2D[0]->GetZaxis()->SetRangeUser(Emin,Emax); 


  if(opt.Contains("units")) {
    trapPotential *=  ( E0 / (PUnits::GV) ); 
  }  
  
  Float_t Vmax = hV1D[0]->GetMaximum();
  Float_t Vmin = hV1D[0]->GetMinimum();
  
  { // Shift potential
    Int_t   Nbins = hV1D[0]->GetNbinsX();
    for(Int_t k=Nbins;k>0;k--) {
      Float_t bincontent = hV1D[0]->GetBinContent(k);
      hV1D[0]->SetBinContent(k,bincontent - Vmax + trapPotential);// - Vmin - trapPotential);
    }
  }
  
  
  
  // "Axis range" in Osiris units:
  Double_t ymin = hDen2D[1]->GetYaxis()->GetBinLowEdge(FirstyBin);
  Double_t ymax = hDen2D[1]->GetYaxis()->GetBinUpEdge(LastyBin);
  TLine *lineYup = new TLine(hDen2D[1]->GetXaxis()->GetXmin(),ymax,
			     hDen2D[1]->GetXaxis()->GetXmax(),ymax);
  lineYup->SetLineColor(kGray+1);
  lineYup->SetLineStyle(2);
 
  TLine *lineYdown = new TLine(hDen2D[1]->GetXaxis()->GetXmin(),ymin,
			     hDen2D[1]->GetXaxis()->GetXmax(),ymin);
  lineYdown->SetLineColor(kGray+1);
  lineYdown->SetLineStyle(2);


  // Plotting
  // -----------------------------------------------

  // Canvas setup
  TCanvas *C;
  if(opt.Contains("hres") && !opt.Contains("pdf")) // high resolution for plain grahics output.
    C = new TCanvas("C","2D Charge density and Electric field",1126,1000);
  else
    C = new TCanvas("C","2D Charge density and Electric field",750,750);

  // Palettes setup
  TExec *exPlasma = new TExec("exPlasma","plasmaPalette->cd();");
  TExec *exElec   = new TExec("exElec","electronPalette->cd();");
  TExec *exHot    = new TExec("exHot","hotPalette->cd();");
  TExec *exField  = new TExec("exField","rbow2Palette->cd();");
  
  // Text objects
  TPaveText *textTime = new TPaveText(0.68,0.85,0.83,0.92,"NDC");
  PlasmaGlob::SetPaveTextStyle(textTime); 
  char ctext[128];
  if(opt.Contains("units") && pData->GetPlasmaDensity()) 
    sprintf(ctext,"z = %5.0f #mum", pData->GetPlasmaSkinDepth() * Time / PUnits::um);
  else
    sprintf(ctext,"t = %5.1f #omega_{p}^{-1}",Time);
  textTime->AddText(ctext);
 
  TPaveText *textDen = new TPaveText(0.13,0.85,0.38,0.92,"NDC");
  PlasmaGlob::SetPaveTextStyle(textDen,12); 
  textDen->SetTextColor(kOrange+10);
  if(opt.Contains("units") && pData->GetPlasmaDensity())
    sprintf(ctext,"n_{0} = %5.2f x 10^{17} / cm^{3}", pData->GetPlasmaDensity() / (1e17/PUnits::cm3));
  else if(pData->GetBeamDensity() && pData->GetPlasmaDensity())
    sprintf(ctext,"n_{b}/n_{0} = %5.2f", pData->GetBeamDensity()/pData->GetPlasmaDensity());
  textDen->AddText(ctext);

  TPaveText *textWav = new TPaveText(0.13,0.05,0.38,0.12,"NDC");
  PlasmaGlob::SetPaveTextStyle(textWav,12); 
  textWav->SetTextColor(kGray+2);
  sprintf(ctext,"#lambda_{p} = %5.3f #mum", pData->GetPlasmaWaveLength()/PUnits::um);
  textWav->AddText(ctext);

  
  // Actual Plotting!
  // ------------------------------------------------------------

  // Output file
  TString fOutName = Form("./%s/Plots/RakeChargeField/RakeChargeField",pData->GetPath().c_str());
  fOutName += Form("-%s_%i",pData->GetName(),time);

  // Setup Pad layout:
  Double_t lMargin = 0.10;
  Double_t rMargin = 0.14;
  Double_t bMargin = 0.15;
  Double_t tMargin = 0.02;
  Double_t vSpacing = 0.01; 
  Double_t hStep = (1.-lMargin-rMargin);
  Double_t vStep1 = (1.-bMargin-tMargin)*(4./10.);
  Double_t vStep2 = (1.-bMargin-tMargin)*(6./10.);
 
  TPad *pad[2];
  TString sLabels[] = {"(a)","(b)"};
  // Text objects
  TPaveText **textLabel = new TPaveText*[2];

  // top plots
  pad[0] = new TPad("padt", "padt",  0.00,   bMargin + vStep1 + vSpacing,
		    lMargin+hStep+rMargin,   1.00);
  pad[0]->SetLeftMargin(1./(lMargin+hStep)*lMargin);
  pad[0]->SetRightMargin(1./(rMargin+hStep)*rMargin);  
  pad[0]->SetBottomMargin(0.0);                                   
  pad[0]->SetTopMargin(1./(tMargin+vStep2)*tMargin);
  pad[0]->Draw();

  // bottom plots
  pad[1] = new TPad("padb", "padb",  0.00,   0.00,
		    lMargin+hStep+rMargin,   bMargin + vStep1);
  pad[1]->SetLeftMargin(1./(lMargin+hStep)*lMargin);
  pad[1]->SetRightMargin((1./(rMargin+hStep)*rMargin));
  pad[1]->SetBottomMargin(1./(bMargin+vStep1)*bMargin);
  pad[1]->SetTopMargin(0.);
  pad[1]->Draw();  

  Float_t yFactor = pad[0]->GetAbsHNDC()/pad[1]->GetAbsHNDC();
  Float_t xFactor = pad[0]->GetAbsWNDC()/pad[1]->GetAbsWNDC();


  pad[0]->cd(); // <---------------------------------------------- Top Plot ---------
  if(opt.Contains("logz")) {
    pad[0]->SetLogz(1);
  } else {
    pad[0]->SetLogz(0);
  }
  pad[0]->SetFrameLineWidth(3);  


  // Re-range
  hDen2D[0]->GetYaxis()->SetRangeUser(yMin - range/2, yMax);
  
  TH2F *hFrame = (TH2F*) gROOT->FindObject("hFrame1");
  if(hFrame) delete hFrame;
  hFrame = (TH2F*) hDen2D[0]->Clone("hFrame1");
  hFrame->Reset();

  hFrame->SetLabelFont(42,"xyz");
  hFrame->SetTitleFont(42,"xyz");

  //  hFrame->GetYaxis()->SetNdivisions(505);
  hFrame->GetYaxis()->SetLabelSize(0.075);
  hFrame->GetYaxis()->SetTitleSize(0.075);
  hFrame->GetYaxis()->SetTitleOffset(0.7);
  hFrame->GetXaxis()->SetLabelOffset(0.10);
  // Frame asymmetry:
  hFrame->Draw("col");
  
  if(opt.Contains("units") && pData->GetPlasmaDensity()) {
    hFrame->GetYaxis()->SetNdivisions(510);
  }

  // hDen2D[0]->GetZaxis()->SetNdivisions(505);


  // Injected electrons if any
  if(Nspecies>=3) {
    if(hDen2D[2]) {
      exHot->Draw();
      hDen2D[2]->Draw("colz same");
    }
  }

  // Plasma
  hDen2D[0]->GetZaxis()->SetTitleFont(42);
  exPlasma->Draw();
  hDen2D[0]->Draw("colz same");
  
  // Beam driver.
  if(hDen2D[1]) {
    //    hDen2D[1]->GetZaxis()->SetNdivisions(505);

    exElec->Draw();
    hDen2D[1]->Draw("colz same");
  }

  // lineYdown->Draw();
  // lineYup->Draw();
   
  // Palettes re-arrangement
  pad[0]->Update();
  Float_t y1 = pad[0]->GetBottomMargin();
  Float_t y2 = 1 - pad[0]->GetTopMargin();
  Float_t x1 = pad[0]->GetLeftMargin();
  Float_t x2 = 1 - pad[0]->GetRightMargin();
 
  TPaletteAxis *palette = (TPaletteAxis*)hDen2D[2]->GetListOfFunctions()->FindObject("palette");
  if(palette) {
    palette->SetY2NDC(y2 - 0.00);
    palette->SetY1NDC(0.66*(y1+y2) + 0.00);
    palette->SetX1NDC(x2 + 0.005);
    palette->SetX2NDC(x2 + 0.03);
    //  palette->SetTitleFont(42);
    //  palette->SetTitleOffset(0.85);
    palette->SetTitleOffset(999.9);
    palette->SetTitleSize(0.075);
    palette->SetLabelFont(42);
    palette->SetLabelSize(0.075);
    palette->SetLabelOffset(0.001);
    palette->SetBorderSize(2);
    palette->SetLineColor(1);
  }
 
  palette = (TPaletteAxis*)hDen2D[0]->GetListOfFunctions()->FindObject("palette");
  if(palette) {
    palette->SetY2NDC(0.66*(y1+y2) - 0.00);
    palette->SetY1NDC(0.33*(y1+y2) + 0.00);
    palette->SetX1NDC(x2 + 0.005);
    palette->SetX2NDC(x2 + 0.03);
    // palette->SetTitleFont(42);
    palette->SetTitleOffset(0.80);
    palette->SetTitleSize(0.075);
    palette->SetLabelFont(42);
    palette->SetLabelSize(0.075);
    palette->SetLabelOffset(0.001);
    palette->SetBorderSize(2);
    palette->SetLineColor(1);
  } 
  
  palette = (TPaletteAxis*)hDen2D[1]->GetListOfFunctions()->FindObject("palette");
  if(palette) {
    palette->SetY2NDC(0.33*(y1+y2) - 0.00);
    palette->SetY1NDC(y1 + 0.00);
    palette->SetX1NDC(x2 + 0.005);
    palette->SetX2NDC(x2 + 0.03);
    //palette->SetTitleFont(42);
    //palette->SetTitleOffset(0.85);
    palette->SetTitleOffset(999.9);
    palette->SetTitleSize(0.075);
    palette->SetLabelFont(42);
    palette->SetLabelSize(0.075);
    palette->SetLabelOffset(0.001);
    palette->SetBorderSize(2);
    palette->SetLineColor(1);
  }
  
  // 1D charge density plots:
  Float_t yaxismin  =  pad[0]->GetUymin();
  Float_t yaxismax  =  pad[0]->GetUymin() + 0.33*(pad[0]->GetUymax() - pad[0]->GetUymin()) - 0.00;
  Float_t denmin = (0.1001) * density;
  Float_t denmax = hDen1D[1]->GetMaximum();
  if(opt.Contains("logz")) {
    denmin = TMath::Log10(denmin);
    denmax = TMath::Log10(denmax);
  }
  for(Int_t i=0;i<Nspecies;i++) {
    if(!hDen1D[i]) continue;
    
    Float_t slope = (yaxismax - yaxismin)/(denmax - denmin);
    
    for(Int_t j=0;j<hDen1D[i]->GetNbinsX();j++) {
      Float_t content = hDen1D[i]->GetBinContent(j+1);
      if(opt.Contains("logz")) content = TMath::Log10(content); 
      
      if(content<denmin) 
	hDen1D[i]->SetBinContent(j+1,yaxismin);
      else 
	hDen1D[i]->SetBinContent(j+1,(content - denmin) * slope + yaxismin);
      
      
    }
    
  }


  // hDen1D[0]->SetLineWidth(2);
  // hDen1D[0]->SetLineColor(kGray+1);
  // // // PlasmaGlob::SetH1Style(hDen1D[0],1);
  // hDen1D[0]->Draw("same C");  
  
  hDen1D[1]->SetLineWidth(2);
  hDen1D[1]->SetLineColor(PlasmaGlob::elecLine);
  // PlasmaGlob::SetH1Style(hDen1D[1],1);
  hDen1D[1]->Draw("same C");
  
  if(Nspecies>=3) {
    if(hDen1D[2]) {
      hDen1D[2]->SetLineWidth(2);
      //      hDen1D[2]->SetLineColor(kGray+2);
      //      hDen1D[2]->SetLineColor(kMagenta-5);
      hDen1D[2]->SetLineColor(kOrange+8);
      hDen1D[2]->Draw("same C");
    }
  }

  TLine *line1 = new TLine(hFrame->GetXaxis()->GetXmin(),yaxismax,
			   hFrame->GetXaxis()->GetXmax(),yaxismax);
  line1->SetLineColor(kGray+2);
  line1->SetLineStyle(2);
  line1->Draw();

  textTime->Draw();
  textDen->Draw();
  if(opt.Contains("units"))
    textWav->Draw();

  textLabel[0] = new TPaveText(x2-0.10,y1+0.04,x2-0.02,y1+(0.12),"NDC"); 
  PlasmaGlob::SetPaveTextStyle(textLabel[0],32); 
  textLabel[0]->SetTextFont(42);
  textLabel[0]->AddText(sLabels[0]);
  textLabel[0]->Draw();
  

  
  pad[0]->RedrawAxis(); 

 
  pad[1]->cd(); // <--- Bottom Plot
  pad[1]->SetFrameLineWidth(3);  
  
  hFrame = (TH2F*) gROOT->FindObject("hFrame2");
  if(hFrame) delete hFrame;
  hFrame = (TH2F*) hE2D[0]->Clone("hFrame2");
  hFrame->Reset();
  
  // hFrame->GetXaxis()->SetNdivisions(505);
  // cout << Form(" Axis range : %8.4f - %8.4f",hFrame->GetXaxis()->GetXmin(),hFrame->GetXaxis()->GetXmax()) << endl;

  
  hFrame->GetYaxis()->SetLabelSize(yFactor*0.075);
  hFrame->GetYaxis()->SetTitleSize(yFactor*0.075);
  hFrame->GetYaxis()->SetTitleOffset(0.7/yFactor);
  
  hFrame->GetXaxis()->SetTitleSize(0.10);
  hFrame->GetXaxis()->SetLabelSize(0.08);
  hFrame->GetXaxis()->SetTitleOffset(1.1);

  hFrame->SetLabelFont(42,"xyz");
  hFrame->SetTitleFont(42,"xyz");

  hFrame->Draw("col");

  //  hE2D[0]->GetZaxis()->SetNdivisions(505);
  hE2D[0]->GetZaxis()->SetTitleFont(42);
  exField->Draw();
  hE2D[0]->Draw("colz same");

  TLine *line0 = new TLine(hE1D[0]->GetXaxis()->GetXmin(),
			   0.0,
			   hE1D[0]->GetXaxis()->GetXmax(),
			   0.0);
  line0->SetLineColor(kGray+3);
  line0->SetLineStyle(2);
  line0->Draw();

  // Fit the V1D in the E2D pad:
  // Float_t rightmin = Vmin;
  // Float_t rightmax = Vmax;
  Float_t rightmin = Emin;
  Float_t rightmax = Emax;
  Float_t slope = yMax/rightmax;
  
  for(Int_t j=0;j<hV1D[0]->GetNbinsX();j++) {
    hV1D[0]->SetBinContent(j+1,hV1D[0]->GetBinContent(j+1)*slope);
  }
  hV1D[0]->SetLineStyle(1);
  hV1D[0]->SetLineWidth(2);
  hV1D[0]->SetLineColor(PlasmaGlob::elecLine);

  hV1D[0]->Draw("sameC");
  
  // Line for trapping potential
  TLine *lineTrap = new TLine(hV1D[0]->GetXaxis()->GetXmin(),
			      (trapPotential)*slope,
			      hV1D[0]->GetXaxis()->GetXmax(),
			      (trapPotential)*slope);
  lineTrap->SetLineColor(PlasmaGlob::elecLine);
  lineTrap->SetLineStyle(2);
  lineTrap->SetLineWidth(1);
  lineTrap->Draw();
  
  // Fit the E1D in the E2D pad:
  rightmin = Emin;
  rightmax = Emax;
  slope = yMax/rightmax;
  
  for(Int_t j=0;j<hE1D[0]->GetNbinsX();j++) {
    hE1D[0]->SetBinContent(j+1,hE1D[0]->GetBinContent(j+1)*slope);
  }
  hE1D[0]->SetLineStyle(1);
  hE1D[0]->SetLineWidth(2);
  hE1D[0]->SetLineColor(kOrange+10);

  hE1D[0]->Draw("sameC");

  Float_t HeIon = 36.0; // GV/m
  if(!opt.Contains("units"))
    HeIon /= ( E0 / (PUnits::GV/PUnits::m) ); 
  
  TLine *lineHe = new TLine(hE1D[0]->GetXaxis()->GetXmin(),
			    HeIon*slope,
			    hE1D[0]->GetXaxis()->GetXmax(),
			    HeIon*slope);
  lineHe->SetLineColor(kGray+3);
  lineHe->SetLineStyle(2);
  // lineHe->Draw();


  TLine *lineHe2 = new TLine(hE1D[0]->GetXaxis()->GetXmin(),
			   -HeIon*slope,
			   hE1D[0]->GetXaxis()->GetXmax(),
			   -HeIon*slope);
  lineHe2->SetLineColor(kOrange+10);
  lineHe2->SetLineStyle(2);
  lineHe2->Draw();
  
  //lineYdown->Draw();
  //lineYup->Draw();

  pad[1]->Update();
  TPaletteAxis *palette3 = (TPaletteAxis*)hE2D[0]->GetListOfFunctions()->FindObject("palette");
  
  y1 = pad[1]->GetBottomMargin();
  y2 = 1 - pad[1]->GetTopMargin();
  x1 = pad[1]->GetLeftMargin();
  x2 = 1 - pad[1]->GetRightMargin();
  palette3->SetY2NDC(y2 - 0.00);
  palette3->SetY1NDC(y1 + 0.00);
  palette3->SetX1NDC(x2 + 0.005);
  palette3->SetX2NDC(x2 + 0.03);
  // palette3->SetTitleFont(42);
  palette3->SetTitleSize(yFactor*0.075);
  palette3->SetTitleOffset(0.80/yFactor);
  palette3->SetLabelSize(yFactor*0.075);
  palette3->SetLabelFont(42);
  palette3->SetLabelOffset(0.01/yFactor);
  palette3->SetBorderSize(2);
  palette3->SetLineColor(1);
   
  pad[1]->RedrawAxis(); 

  // hE1D[0]->GetYaxis()->SetRangeUser(Emin,Emax);;
  // hE1D[0]->Draw();


  textLabel[1] = new TPaveText(x2-0.10,y1+(0.04),x2-0.02,y1+(0.12),"NDC"); 
  PlasmaGlob::SetPaveTextStyle(textLabel[1],32); 
  textLabel[1]->SetTextFont(42);
  textLabel[1]->AddText(sLabels[1]);
  textLabel[1]->Draw();
  


  C->cd();

  // Print to a file
  PlasmaGlob::imgconv(C,fOutName,opt);
  // ---------------------------------------------------------


  // delete objects:
  delete pad[1];
  delete pad[0];
  delete C;

  delete textWav;
  delete textDen;
  delete textTime;
  
  delete exPlasma;
  delete exElec;
  delete exField;
}
