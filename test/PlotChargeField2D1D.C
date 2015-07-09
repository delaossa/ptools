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

void PlotChargeField2D1D( const TString &sim, Int_t time, Int_t zoom=2, Int_t Nbins=2, const TString &options="") {
  
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
    
    if(i==0)
      hDen2D[i]->GetZaxis()->SetTitle("n_{e}/n_{0}");
    else
      hDen2D[i]->GetZaxis()->SetTitle("n_{b}/n_{0}");
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
    
    // if(hE1D[i]) delete hE1D[i];
    // hE1D[i] = (TH1F*) hE2D[i]->ProjectionX(hName,FirstyBin,LastyBin);
    // hE1D[i]->Scale(1.0/(LastyBin-FirstyBin+1));
    
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
      for(Int_t j=0;j<hDen2D[i]->GetNbinsX();j++) {
	for(Int_t k=0;k<hDen2D[i]->GetNbinsY();k++) {
	  hDen2D[i]->SetBinContent(j,k, hDen2D[i]->GetBinContent(j,k) * n0 / (1e17/PUnits::cm3) );
	}
      }

      if(CYL)
	hDen2D[i]->GetYaxis()->SetTitle("r [#mum]");      
      else
	hDen2D[i]->GetYaxis()->SetTitle("y [#mum]");      

      if(opt.Contains("comov"))
	hDen2D[i]->GetXaxis()->SetTitle("#zeta [#mum]");
      else
	hDen2D[i]->GetXaxis()->SetTitle("z [#mum]");
      
      if(i==0)
	hDen2D[i]->GetZaxis()->SetTitle("n_{e} [10^{17}/cm^{3}]");
      else if(i==1)
	hDen2D[i]->GetZaxis()->SetTitle("n_{b} [10^{17}/cm^{3}]"); 
      else
	hDen2D[i]->GetZaxis()->SetTitle("n_{i} [10^{17}/cm^{3}]"); 
    }

    for(Int_t i=0;i<Nspecies;i++) {
      Int_t NbinsX = hDen1D[i]->GetNbinsX();
      Float_t xMin = skindepth * hDen1D[i]->GetXaxis()->GetXmin() / PUnits::um;
      Float_t xMax = skindepth * hDen1D[i]->GetXaxis()->GetXmax() / PUnits::um;
      hDen1D[i]->SetBins(NbinsX,xMin,xMax);
      for(Int_t j=0;j<hDen1D[i]->GetNbinsX();j++) {
    	hDen1D[i]->SetBinContent(j, hDen1D[i]->GetBinContent(j) * n0 / (1e17/PUnits::cm3) );
      }
      

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
      
      for(Int_t j=0;j<NbinsX;j++) {
	hE1D[i]->SetBinContent(j, hE1D[i]->GetBinContent(j) * ( E0 / (PUnits::GV/PUnits::m) ) );
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
  if(opt.Contains("units") && n0)
    density = n0 / (1e17/PUnits::cm3);

  Float_t Max  = hDen2D[0]->GetMaximum();
  Float_t Min  = hDen2D[0]->GetMinimum();

  Float_t Base = density;
  if(Max>Base && Max<2*Base) {
    Min  = 2.* Base - Max;
  } else if(Max <= Base)  {
    Max = 1.01 * Base;
    Min = 0;
  } else if(Max >= 2.0*Base) {
    Min = 0;
  }
  
  
  // if(opt.Contains("fixz")) {
  //   if(sim.Contains("pitz")){
  //     Min = 0.75;
  //     Max = 1.25;
      
  //     Minb = 0.0001;
  //     Maxb = 0.065;
  //   } else if(sim.Contains("flash")) {
  //     Min = 0.001;
  //     Max = 10.0;
      
  //     Minb = 0.01;
  //     Maxb = 6.00;
  //   }
  // }
  
  hDen2D[0]->GetZaxis()->SetRangeUser(Min,Max); 

  if(hDen2D[1]) {
  
    Float_t Minb = 0.005 * hDen2D[1]->GetMaximum();
    Float_t Maxb = hDen2D[1]->GetMaximum();
  
    hDen2D[1]->GetZaxis()->SetRangeUser(Minb,Maxb);

  } 

  if(Nspecies>=3) {
    if(hDen2D[2]) {
      
      Float_t Minb = 0.01 * hDen2D[2]->GetMaximum();
      Float_t Maxb = hDen2D[2]->GetMaximum();
      
      hDen2D[2]->GetZaxis()->SetRangeUser(Minb,Maxb);
      
    } 
  }

  // Dynamic plasma palette
  const Int_t plasmaDNRGBs = 3;
  const Int_t plasmaDNCont = 64;
  Double_t basePos = 0.5;
  if(Max!=Min) {
    if(opt.Contains("logz")) {
      // cout << Form(" log Min = %f  log Base = %f  log Max = %f",TMath::Log10(Min+1),TMath::Log10(Base+1),TMath::Log10(Max+1)) << endl;
      basePos = (1.0/(TMath::Log10(Max)-TMath::Log10(Min+1)))*(TMath::Log10(Base+1) - TMath::Log10(Min+1));
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
    C = new TCanvas("C","2D Charge density and Electric field",750,666);

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
    sprintf(ctext,"Z = %5.0f #mum", pData->GetPlasmaSkinDepth() * Time / PUnits::um);
  else
    sprintf(ctext,"T = %5.1f #omega_{p}^{-1}",Time);
  textTime->AddText(ctext);
 
  TPaveText *textDen = new TPaveText(0.13,0.85,0.38,0.92,"NDC");
  PlasmaGlob::SetPaveTextStyle(textDen,12); 
  textDen->SetTextColor(kOrange+10);
  if(opt.Contains("units") && pData->GetPlasmaDensity())
    sprintf(ctext,"n_{0} = %5.2f x 10^{17} / cm^{3}", pData->GetPlasmaDensity() / (1e17/PUnits::cm3));
  else if(pData->GetBeamDensity() && pData->GetPlasmaDensity())
    sprintf(ctext,"n_{b}/n_{0} = %5.2f", pData->GetBeamDensity()/pData->GetPlasmaDensity());
  textDen->AddText(ctext);

  TPaveText *textWav = new TPaveText(0.13,0.75,0.38,0.82,"NDC");
  PlasmaGlob::SetPaveTextStyle(textWav,12); 
  textWav->SetTextColor(kGray+2);
  sprintf(ctext,"#lambda_{p} = %5.3f #mum", pData->GetPlasmaWaveLength()/PUnits::um);
  textWav->AddText(ctext);

  
  // Actual Plotting!
  // ------------------------------------------------------------

  // Output file
  TString fOutName = Form("./%s/Plots/ChargeField2D1D/ChargeField2D1D",pData->GetPath().c_str());
  fOutName += Form("-%s_%i",pData->GetName(),time);

  // Setup Pad layout:
  Double_t lMargin = 0.08;
  Double_t rMargin = 0.14;
  Double_t bMargin = 0.10;
  Double_t tMargin = 0.02;
  Double_t vSpacing = 0.01; 
  Double_t hStep = (1.-lMargin-rMargin);
  Double_t vStep = (1.-bMargin-tMargin)/2.;
 
  TPad *pad[2];

  // top plots
  pad[0] = new TPad("padt", "padt",  0.00,   bMargin + vStep + vSpacing,
		    lMargin+hStep+rMargin,   1.00);
  pad[0]->SetLeftMargin(1./(lMargin+hStep)*lMargin);
  pad[0]->SetRightMargin(1./(rMargin+hStep)*rMargin);  
  pad[0]->SetBottomMargin(0.0);                                   
  pad[0]->SetTopMargin(1./(tMargin+vStep)*tMargin);
  pad[0]->Draw();

  // bottom plots
  pad[1] = new TPad("padb", "padb",  0.00,   0.00,
		    lMargin+hStep+rMargin,   bMargin+vStep);
  pad[1]->SetLeftMargin(1./(lMargin+hStep)*lMargin);
  pad[1]->SetRightMargin((1./(rMargin+hStep)*rMargin));
  pad[1]->SetBottomMargin(1./(bMargin+vStep)*bMargin);
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


  TH2F *hFrame = (TH2F*) gROOT->FindObject("hFrame1");
  if(hFrame) delete hFrame;
  hFrame = (TH2F*) hDen2D[0]->Clone("hFrame1");
  hFrame->Reset();
  //  hFrame->GetYaxis()->SetNdivisions(505);
  hFrame->GetYaxis()->SetLabelSize(0.075);
  hFrame->GetYaxis()->SetTitleSize(0.075);
  hFrame->GetYaxis()->SetTitleOffset(0.5);
  hFrame->GetXaxis()->SetLabelOffset(0.10);
  hFrame->Draw("col");
   
  hDen2D[0]->GetZaxis()->SetLabelSize(0.06);  
  hDen2D[0]->GetZaxis()->SetTitleSize(0.06);
  hDen2D[0]->GetZaxis()->SetTitleOffset(0.06);
  //  hDen2D[0]->GetZaxis()->SetNdivisions(505);


  if(Nspecies>=3) {
    if(hDen2D[2]) {
      exHot->Draw();
      hDen2D[2]->Draw("colz same");
    }
  }


  exPlasma->Draw();
  hDen2D[0]->Draw("col same");
  
  if(hDen2D[1]) {
    hDen2D[1]->GetZaxis()->SetLabelSize(0.06);  
    hDen2D[1]->GetZaxis()->SetTitleSize(0.06);
    hDen2D[1]->GetZaxis()->SetTitleOffset(0.06);
    //    hDen2D[1]->GetZaxis()->SetNdivisions(505);

    exElec->Draw();
    hDen2D[1]->Draw("colz same");
  }

  // lineYdown->Draw();
  // lineYup->Draw();
   
  pad[0]->Update();
  TPaletteAxis *palette = (TPaletteAxis*)hDen2D[2]->GetListOfFunctions()->FindObject("palette");
  
  Float_t y1 = pad[0]->GetBottomMargin();
  Float_t y2 = 1 - pad[0]->GetTopMargin();
  Float_t x1 = 1 - pad[0]->GetRightMargin();
  palette->SetY2NDC(y2 - 0.01);
  palette->SetY1NDC(0.5*(y1+y2) + 0.01);
  palette->SetX1NDC(x1 + 0.005);
  palette->SetX2NDC(x1 + 0.03);
  palette->SetTitleOffset(0.7);
  palette->SetTitleSize(0.07);
  palette->SetBorderSize(2);
  palette->SetLineColor(1);
  
  TPaletteAxis *palette2 = (TPaletteAxis*)hDen2D[1]->GetListOfFunctions()->FindObject("palette");
  palette2->SetY2NDC(0.5*(y1+y2) - 0.01);
  palette2->SetY1NDC(y1 + 0.01);
  palette2->SetX1NDC(x1 + 0.005);
  palette2->SetX2NDC(x1 + 0.03);
  palette2->SetTitleOffset(0.7);
  palette2->SetTitleSize(0.07);
  palette2->SetBorderSize(2);
  palette2->SetLineColor(1);

  // 1D charge density plots:
  Float_t sfactor = 3.0;
  for(Int_t i=0;i<Nspecies;i++) {
    if(!hDen1D[i]) continue;
    
    Float_t rightmin = 0;
    Float_t rightmax = sfactor * hDen1D[i]->GetMaximum();
    Float_t slope = (pad[0]->GetUymax() - pad[0]->GetUymin())/(rightmax-rightmin);
    
    for(Int_t j=0;j<hDen1D[i]->GetNbinsX();j++) {
      hDen1D[i]->SetBinContent(j+1,hDen1D[i]->GetBinContent(j+1)*slope + pad[0]->GetUymin());
    }        
  }
  
  
  hDen1D[1]->SetLineWidth(2);
  hDen1D[1]->SetLineColor(PlasmaGlob::elecLine);
  // PlasmaGlob::SetH1Style(hDen1D[1],1);
  hDen1D[1]->Draw("same C");
  
  if(Nspecies>=3) {
    if(hDen1D[2]) {
      hDen1D[2]->SetLineWidth(2);
      //      hDen1D[2]->SetLineColor(kGray+2);
      //      hDen1D[2]->SetLineColor(kMagenta-5);
      hDen1D[2]->SetLineColor(kOrange+10);
      hDen1D[2]->Draw("same C");
    }
  }

  TLine *line1 = new TLine(hFrame->GetXaxis()->GetXmin(),pad[0]->GetUymin()/sfactor,
			   hFrame->GetXaxis()->GetXmax(),pad[0]->GetUymin()/sfactor);
  line1->SetLineColor(kGray+2);
  line1->SetLineStyle(2);
  line1->Draw();

  textTime->Draw();
  textDen->Draw();
  if(opt.Contains("units"))
    textWav->Draw();
  
  pad[0]->RedrawAxis(); 

 
  pad[1]->cd(); // <--- Mid Plot
  pad[1]->SetFrameLineWidth(3);  
  
  hFrame = (TH2F*) gROOT->FindObject("hFrame2");
  if(hFrame) delete hFrame;
  hFrame = (TH2F*) hE2D[0]->Clone("hFrame2");
  hFrame->Reset();
  
  // hFrame->GetXaxis()->SetNdivisions(505);
  // cout << Form(" Axis range : %8.4f - %8.4f",hFrame->GetXaxis()->GetXmin(),hFrame->GetXaxis()->GetXmax()) << endl;
  
  hFrame->GetYaxis()->SetLabelSize(yFactor*0.075);
  hFrame->GetYaxis()->SetTitleSize(yFactor*0.075);
  hFrame->GetYaxis()->SetTitleOffset(0.5/yFactor);
  
  hFrame->GetXaxis()->SetTitleSize(0.07);
  hFrame->GetXaxis()->SetLabelSize(0.06);
  hFrame->GetXaxis()->SetTitleOffset(1.1);
 
  hFrame->Draw("col");

  hE2D[0]->GetZaxis()->SetLabelSize(yFactor*0.06);
  hE2D[0]->GetZaxis()->SetTitleSize(yFactor*0.06);
  hE2D[0]->GetZaxis()->SetTitleOffset(yFactor*1.);
  //  hE2D[0]->GetZaxis()->SetNdivisions(505);

  exField->Draw();
  hE2D[0]->Draw("colz same");


  // Fit the E1D in the E2D pad:
  Float_t rightmin = Emin;
  Float_t rightmax = Emax;
  Float_t slope = yMax/rightmax;
  
  for(Int_t j=0;j<hE1D[0]->GetNbinsX();j++) {
    hE1D[0]->SetBinContent(j+1,hE1D[0]->GetBinContent(j+1)*slope);
  }
  hE1D[0]->SetLineStyle(1);
  hE1D[0]->SetLineWidth(2);
  hE1D[0]->SetLineColor(kOrange+10);

  hE1D[0]->Draw("sameC");
  pad[1]->Update();
  
  TLine *line0 = new TLine(hE1D[0]->GetXaxis()->GetXmin(),
			   (pad[1]->GetUymin()+pad[1]->GetUymax())/2.,
			   hE1D[0]->GetXaxis()->GetXmax(),
			   (pad[1]->GetUymin()+pad[1]->GetUymax())/2.);
  line0->SetLineColor(kGray+1);
  line0->SetLineStyle(2);
  line0->Draw();

  Float_t HeIon = 52.5; // GV/m
  if(!opt.Contains("units"))
    HeIon /= ( E0 / (PUnits::GV/PUnits::m) ); 

  TLine *lineHe = new TLine(hE1D[0]->GetXaxis()->GetXmin(),
			   HeIon*slope,
			   hE1D[0]->GetXaxis()->GetXmax(),
			   HeIon*slope);
  lineHe->SetLineColor(kGray+3);
  lineHe->SetLineStyle(2);
  lineHe->Draw();


  Float_t HeIon2 = -52.5; // GV/m
  if(!opt.Contains("units"))
    HeIon2 /=  ( E0 / (PUnits::GV/PUnits::m) ); 

  TLine *lineHe2 = new TLine(hE1D[0]->GetXaxis()->GetXmin(),
			   HeIon2*slope,
			   hE1D[0]->GetXaxis()->GetXmax(),
			   HeIon2*slope);
  lineHe2->SetLineColor(kGray+3);
  lineHe2->SetLineStyle(2);
  lineHe2->Draw();


  //lineYdown->Draw();
  //lineYup->Draw();

  pad[1]->Update();
  TPaletteAxis *palette3 = (TPaletteAxis*)hE2D[0]->GetListOfFunctions()->FindObject("palette");
  
  y1 = pad[1]->GetBottomMargin();
  y2 = 1 - pad[1]->GetTopMargin();
  x1 = 1 - pad[1]->GetRightMargin();
  palette3->SetY2NDC(y2 - 0.01);
  palette3->SetY1NDC(y1 + 0.01);
  palette3->SetX1NDC(x1 + 0.005);
  palette3->SetX2NDC(x1 + 0.03);
  palette3->SetTitleOffset(0.7);
  palette3->SetTitleSize(0.07);
  palette3->SetBorderSize(2);
  palette3->SetLineColor(1);
   
  pad[1]->RedrawAxis(); 

  // hE1D[0]->GetYaxis()->SetRangeUser(Emin,Emax);;
  // hE1D[0]->Draw();

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
