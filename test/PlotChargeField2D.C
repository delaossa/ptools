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

void PlotChargeField2D( const TString &sim, Int_t time, Int_t zoom=2, Int_t Nbins=2, const TString &options="") {
  
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
  if(opt.Contains("grid")) {
    gStyle->SetPadGridX(1);
    gStyle->SetPadGridY(1);
  }
 
  gStyle->SetTitleFont(42);
  gStyle->SetStatFont(42);
  gStyle->SetTextFont(42);
  gStyle->SetTitleFont(42,"xyz");
  gStyle->SetLabelFont(42,"xyz");

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
  Float_t nb = pData->GetBeamDensity();
  
  // Time in OU
  Float_t Time = pData->GetRealTime();
  // z start of the plasma in normalized units.
  Float_t zStartPlasma = pData->GetPlasmaStart()*kp;
  // z start of the beam in normalized units.
  Float_t zStartBeam = pData->GetBeamStart()*kp;
  // z start of the neutral in normalized units.
  Float_t zStartNeutral = pData->GetNeutralStart()*kp;
  // z end of the neutral in normalized units.
  Float_t zEndNeutral = pData->GetNeutralEnd()*kp;
  
  if(opt.Contains("center")) {
    Time -= zStartPlasma;
    if(opt.Contains("comov"))      // Centers on the head of the beam.
      Time += zStartBeam;
  }
  Float_t shiftz = pData->Shift(opt);
  //  cout << "Shift = " << shiftz << endl;
  
  // Calculate the "axis range" in number of bins. If Nbins==0 a RMS width is taken.
  Double_t rms0 = pData->GetBeamRmsY() * kp;
  if(pData->IsCyl())  rms0  = pData->GetBeamRmsR() * kp;
  
  Int_t FirstyBin = 0;
  Int_t LastyBin = 0;
  if(Nbins==0) { 
    Nbins =  TMath::Nint(rms0 / pData->GetDX(1));
  }
  
  // Slice width limits.
  if(!pData->IsCyl()) {
    FirstyBin = pData->GetNX(1)/2 + 1 - Nbins;
    LastyBin =  pData->GetNX(1)/2 + Nbins;
  } else {
    FirstyBin = 1; 
    LastyBin  = Nbins;
  }
  // -------------------------------------------------------------------------------

  // Get charge density histos
  Int_t Nspecies = pData->NSpecies();
  TH2F **hDen2D = new TH2F*[Nspecies];
  TH1F **hDen1D = new TH1F*[Nspecies];
  for(Int_t i=0;i<Nspecies;i++) {

    hDen2D[i] = NULL;
    hDen1D[i] = NULL;
    
    if(!pData->GetChargeFileName(i)) 
      continue;

    cout << Form(" Getting charge density of specie: ") << i << endl;

    
    char hName[24];
    sprintf(hName,"hDen2D_%i",i);
    hDen2D[i] = (TH2F*) gROOT->FindObject(hName);
    if(hDen2D[i]) delete hDen2D[i];

    if(!pData->Is3D())
      hDen2D[i] = pData->GetCharge(i,opt);
    else
      hDen2D[i] = pData->GetCharge2DSliceZY(i,-1,Nbins,opt+"avg");

    //  cout << Form(" Charge density of specie '%s' loaded into histogram '%s'",pData->GetSpeciesName(i).c_str(),hName) << endl;
    
    hDen2D[i]->SetName(hName);
    hDen2D[i]->GetXaxis()->CenterTitle();
    hDen2D[i]->GetYaxis()->CenterTitle();
    hDen2D[i]->GetZaxis()->CenterTitle();
    
    if(opt.Contains("comov"))
      hDen2D[i]->GetXaxis()->SetTitle("#zeta [c/#omega_{p}]");
    else
      hDen2D[i]->GetXaxis()->SetTitle("z [c/#omega_{p}]");
    
    if(pData->IsCyl()) 
      hDen2D[i]->GetYaxis()->SetTitle("r [c/#omega_{p}]");
    else
      hDen2D[i]->GetYaxis()->SetTitle("y [c/#omega_{p}]");

    hDen2D[i]->GetZaxis()->SetTitle("n [n_{0}]");
    
    if(!opt.Contains("1dline")) continue;
    
    sprintf(hName,"hDen1D_%i",i);
    hDen1D[i] = (TH1F*) gROOT->FindObject(hName);
    if(hDen1D[i]) delete hDen1D[i];
    
    if(pData->Is3D()) {
      hDen1D[i] = pData->GetH1SliceZ3D(pData->GetChargeFileName(i)->c_str(),"charge",-1,Nbins,-1,Nbins,opt+"avg");
    } else if(pData->IsCyl()) { // Cylindrical: The first bin with r>0 is actually the number 1 (not the 0).
      hDen1D[i] = pData->GetH1SliceZ(pData->GetChargeFileName(i)->c_str(),"charge",1,Nbins,opt+"avg");
    } else { // 2D cartesian
      hDen1D[i] = pData->GetH1SliceZ(pData->GetChargeFileName(i)->c_str(),"charge",-1,Nbins,opt+"avg");
    }
    
    //  cout << Form(" Charge density of specie '%s' loaded into histogram '%s'",pData->GetSpeciesName(i).c_str(),hName) << endl;
    
    hDen1D[i]->SetName(hName);
    hDen1D[i]->GetXaxis()->CenterTitle();
    hDen1D[i]->GetYaxis()->CenterTitle();
       
    if(opt.Contains("comov"))
      hDen1D[i]->GetXaxis()->SetTitle("#zeta [c/#omega_{p}]");
    else
      hDen1D[i]->GetXaxis()->SetTitle("z [c/#omega_{p}]");
    
    hDen1D[i]->GetYaxis()->SetTitle("n [n_{0}]");


  }
  
  // Get electric fields
  const Int_t Nfields = 1;
  TH2F **hE2D = new TH2F*[Nfields];
  TH1F **hE1D = new TH1F*[Nfields];
  for(Int_t i=0;i<Nfields;i++) {
    hE2D[i] = NULL;
    hE1D[i] = NULL;

    if(!pData->GetEfieldFileName(i))
      continue;

    cout << Form(" Getting electric field number ") << i+1 << endl;


    char hName[24];
    sprintf(hName,"hE2D_%i",i);
    hE2D[i] = (TH2F*) gROOT->FindObject(hName);
    if(hE2D[i]) delete hE2D[i];

    if(!pData->Is3D())
      hE2D[i] = pData->GetEField(i,opt);
    else
      hE2D[i] = pData->GetEField2DSliceZY(i,-1,Nbins,opt+"avg");
    
   
    hE2D[i]->SetName(hName);   
    hE2D[i]->GetXaxis()->CenterTitle();
    hE2D[i]->GetYaxis()->CenterTitle();
    hE2D[i]->GetZaxis()->CenterTitle();
    if(opt.Contains("comov"))
      hE2D[i]->GetXaxis()->SetTitle("#zeta [c/#omega_{p}]");
    else
      hE2D[i]->GetXaxis()->SetTitle("z [c/#omega_{p}]");
    
    if(pData->IsCyl()) 
      hE2D[i]->GetYaxis()->SetTitle("r [c/#omega_{p}]");
    else
      hE2D[i]->GetYaxis()->SetTitle("y [c/#omega_{p}]");
    
    if(i==0)
      hE2D[i]->GetZaxis()->SetTitle("E_{z} [E_{0}]");
    else if(i==1)
      hE2D[i]->GetZaxis()->SetTitle("E_{y} [E_{0}]");
    else if(i==2)
      hE2D[i]->GetZaxis()->SetTitle("E_{x} [E_{0}]");

    if(!opt.Contains("1dline")) continue;
     
    sprintf(hName,"hE1D_%i",i);
    hE1D[i] = (TH1F*) gROOT->FindObject(hName);
    if(hE1D[i]) delete hE1D[i];
    
    // 1D histograms
    char nam[3]; sprintf(nam,"e%i",i+1);
    if(pData->Is3D()) {
      
      if(i==0) 
	hE1D[i] = pData->GetH1SliceZ3D(pData->GetEfieldFileName(i)->c_str(),nam,-1,Nbins,-1,Nbins,opt+"avg");
      else  
	hE1D[i] = pData->GetH1SliceZ3D(pData->GetEfieldFileName(i)->c_str(),nam,-Nbins,Nbins,-Nbins,Nbins,opt+"avg");
      
    } else if(pData->IsCyl()) { // Cylindrical: The first bin with r>0 is actually the number 1 (not the 0).
      
      hE1D[i] = pData->GetH1SliceZ(pData->GetEfieldFileName(i)->c_str(),nam,1,Nbins,opt+"avg");
      
    } else { // 2D cartesian
      
      if(i==0) 
	hE1D[i] = pData->GetH1SliceZ(pData->GetEfieldFileName(i)->c_str(),nam,-1,Nbins,opt+"avg");
      else 
	hE1D[i] = pData->GetH1SliceZ(pData->GetEfieldFileName(i)->c_str(),nam,-Nbins,Nbins,opt+"avg");    
      
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
  
  // Now, combine the electric field components into the total |E|
  // and calculate ionization probability for He:
  // Outter Helium electron
  Double_t Eion0 = 24.59 * PUnits::eV;
  Double_t Z     = 1;
  Double_t l     = 0;
  Double_t m     = 0;
  
  TH1F *hIonProb1D = NULL;
  if(opt.Contains("ionprob")) {
    hIonProb1D = (TH1F*) hE1D[0]->Clone("hIonProb1D");
    hIonProb1D->Reset();
    Int_t NbinsX = hE1D[0]->GetNbinsX();
    for(Int_t j=1;j<=NbinsX;j++) {
      Double_t E1 = hE1D[0]->GetBinContent(j);
      Double_t E2 = hE1D[1]->GetBinContent(j);
      Double_t E3 = hE1D[2]->GetBinContent(j);
      Double_t E  = TMath::Sqrt(E1*E1+E2*E2+E3*E3);
      
      E *= E0;
      
      Double_t IonProb = (PFunc::ADK(E,Eion0,Z,l,m)/PUnits::atomictime)*PUnits::femtosecond;
      hIonProb1D->SetBinContent(j,IonProb);
    }
    hIonProb1D->GetYaxis()->SetTitle("W_{ADK} [fs^{-1}]");
  }
  
  // Tunning the Histograms
  // ---------------------
  
  // Chaning to user units: 
  // --------------------------
  
  if(opt.Contains("units") && n0) {
    
    for(Int_t i=0;i<Nspecies;i++) {

      if(!hDen2D[i]) continue;

      Int_t NbinsX = hDen2D[i]->GetNbinsX();
      Float_t xMin = skindepth * hDen2D[i]->GetXaxis()->GetXmin() / PUnits::um;
      Float_t xMax = skindepth * hDen2D[i]->GetXaxis()->GetXmax() / PUnits::um;
      Int_t NbinsY = hDen2D[i]->GetNbinsY();
      Float_t yMin = skindepth * hDen2D[i]->GetYaxis()->GetXmin() / PUnits::um;
      Float_t yMax = skindepth * hDen2D[i]->GetYaxis()->GetXmax() / PUnits::um;
      hDen2D[i]->SetBins(NbinsX,xMin,xMax,NbinsY,yMin,yMax);
      // for(Int_t j=0;j<hDen2D[i]->GetNbinsX();j++) {
      // 	for(Int_t k=0;k<hDen2D[i]->GetNbinsY();k++) {
      // 	  hDen2D[i]->SetBinContent(j,k, hDen2D[i]->GetBinContent(j,k) * n0 / (1e15/PUnits::cm3) );
      // 	}
      // }

      if(pData->IsCyl())
	hDen2D[i]->GetYaxis()->SetTitle("r [#mum]");      
      else
	hDen2D[i]->GetYaxis()->SetTitle("y [#mum]");      

      if(opt.Contains("comov"))
	hDen2D[i]->GetXaxis()->SetTitle("#zeta [#mum]");
      else
	hDen2D[i]->GetXaxis()->SetTitle("z [#mum]");
      
      //      hDen2D[i]->GetZaxis()->SetTitle("n [10^{15}/cm^{3}]");
      
      if(!hDen1D[i]) continue;
      
      hDen1D[i]->SetBins(NbinsX,xMin,xMax);
      
      if(opt.Contains("comov"))
    	hDen1D[i]->GetXaxis()->SetTitle("#zeta [#mum]");
      else
    	hDen1D[i]->GetXaxis()->SetTitle("z [#mum]");
      
    }

    for(Int_t i=0;i<Nfields;i++) {
      Int_t NbinsX = hE2D[i]->GetNbinsX();
      Float_t xMin = skindepth * hE2D[i]->GetXaxis()->GetXmin() / PUnits::um;
      Float_t xMax = skindepth * hE2D[i]->GetXaxis()->GetXmax() / PUnits::um;
      Int_t NbinsY = hE2D[i]->GetNbinsY();
      Float_t yMin = skindepth * hE2D[i]->GetYaxis()->GetXmin() / PUnits::um;
      Float_t yMax = skindepth * hE2D[i]->GetYaxis()->GetXmax() / PUnits::um;
      hE2D[i]->SetBins(NbinsX,xMin,xMax,NbinsY,yMin,yMax);
      if(hE1D[i]) 
	hE1D[i]->SetBins(NbinsX,xMin,xMax);
	

      for(Int_t j=0;j<hE2D[i]->GetNbinsX();j++) {	
	for(Int_t k=0;k<hE2D[i]->GetNbinsY();k++) {
	  hE2D[i]->SetBinContent(j,k, hE2D[i]->GetBinContent(j,k) * ( E0 / (PUnits::GV/PUnits::m) ) );
	}
	
	if(!hE1D[i]) continue;
	hE1D[i]->SetBinContent(j, hE1D[i]->GetBinContent(j) * ( E0 / (PUnits::GV) ) );
	
      }
      
      if(pData->IsCyl())
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


      if(!hE1D[i]) continue;
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
    if(hIonProb1D) {
      Int_t NbinsX = hIonProb1D->GetNbinsX();
      Float_t xMin = skindepth * hIonProb1D->GetXaxis()->GetXmin() / PUnits::um;
      Float_t xMax = skindepth * hIonProb1D->GetXaxis()->GetXmax() / PUnits::um;
      hIonProb1D->SetBins(NbinsX,xMin,xMax);
    }
  }
  
  // --------------------------------------------------- Vertical Zoom ------------

  Float_t range    = (hDen2D[0]->GetYaxis()->GetXmax() - hDen2D[0]->GetYaxis()->GetXmin())/zoom;
  Float_t midPoint = (hDen2D[0]->GetYaxis()->GetXmax() + hDen2D[0]->GetYaxis()->GetXmin())/2.;
  Double_t ymin = midPoint-range/2;
  Double_t ymax = midPoint+range/2;
  if(pData->IsCyl()) {
    ymin = hDen2D[0]->GetYaxis()->GetXmin();
    ymax = range;
  }

  for(Int_t i=0;i<Nspecies;i++) {
    if(!hDen2D[i]) continue;
    hDen2D[i]->GetYaxis()->SetRangeUser(ymin,ymax);
  }
  
  for(Int_t i=0;i<Nfields;i++) {
    if(!hE2D[i]) continue;
    hE2D[i]->GetYaxis()->SetRangeUser(ymin,ymax);
  }
    
  // ------------- z Zoom --------------------------------- Plasma palette -----------
  // Set the range of the plasma charge density histogram for maximum constrast 
  // using a dynamic palette wich adjust the nominal value to a certain color.
 

  Float_t density = 1;
  Float_t Base  = density;
  
  Float_t *Max = new Float_t[Nspecies];
  Float_t *Min = new Float_t[Nspecies];
  
  for(Int_t i=0;i<Nspecies;i++) {
    if(!hDen2D[i]) continue;
    
    Max[i] = hDen2D[i]->GetMaximum();
    Min[i] = 1.01E-1 * Base;  

    if(i==2) Min[i] = 1.01E-3 * Base;
    if(sim.Contains("pitz")) {
      if(i==0) {
	if(Max[i]<2) {
	  Min[i] = Base - (Max[i] - Base);
	}
      } else if(i==1) {
	Min[i] = 0.501E-3 * Base;
      }
    }
    
    hDen2D[i]->GetZaxis()->SetRangeUser(Min[i],Max[i]);
  }

  // Dynamic plasma palette
  const Int_t plasmaDNRGBs = 3;
  const Int_t plasmaDNCont = 64;
  Double_t basePos = 0.5;
  if(Max[0]!=Min[0]) {
    if(opt.Contains("logz")) {
      Float_t a = 1.0/(TMath::Log10(Max[0])-TMath::Log10(Min[0]));
      Float_t b = TMath::Log10(Min[0]);
      basePos = a*(TMath::Log10(Base) - b);
      
    } else {
      basePos = (1.0/(Max[0]-Min[0]))*(Base - Min[0]);
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
  for(Int_t i=0;i<Nfields;i++) {
    if(!hE2D[i]) continue;
    Float_t Emax = hE2D[i]->GetMaximum();
    Float_t Emin = hE2D[i]->GetMinimum();
    if(Emax > TMath::Abs(Emin))
      Emin = -Emax;
    else
      Emax = -Emin;
    hE2D[i]->GetZaxis()->SetRangeUser(Emin,Emax); 
  }
  
  // "Axis range" in Osiris units:
  Double_t ylow  = hDen2D[1]->GetYaxis()->GetBinLowEdge(FirstyBin);
  Double_t yup = hDen2D[1]->GetYaxis()->GetBinUpEdge(LastyBin);
  Double_t xmin = hDen2D[1]->GetXaxis()->GetXmin();
  Double_t xmax = hDen2D[1]->GetXaxis()->GetXmax();

  TLine *lineYzero = new TLine(xmin,0.0,xmax,0.0);
  lineYzero->SetLineColor(kGray+2);
  lineYzero->SetLineStyle(2);

  TLine *lineYup = new TLine(xmin,yup,xmax,yup);
  lineYup->SetLineColor(kGray+1);
  lineYup->SetLineStyle(2);
 
  TLine *lineYdown = new TLine(xmin,ylow,xmax,ylow);
  lineYdown->SetLineColor(kGray+1);
  lineYdown->SetLineStyle(2);

  zStartPlasma -= shiftz; 
  zStartNeutral -= shiftz; 
  zEndNeutral -= shiftz; 
  
  if(opt.Contains("units")) {
    zStartPlasma *= skindepth / PUnits::um;
    zStartNeutral *= skindepth / PUnits::um;
    zEndNeutral *= skindepth / PUnits::um;
  }

  //  cout << "Start plasma = " << zStartPlasma << endl;
  TLine *lineStartPlasma = new TLine(zStartPlasma,ymin,zStartPlasma,ymax);
  lineStartPlasma->SetLineColor(kRed);
  lineStartPlasma->SetLineStyle(1);
  lineStartPlasma->SetLineWidth(2);

  //  cout << "Start plasma = " << zStartNeutral << endl;
  TLine *lineStartNeutral = new TLine(zStartNeutral,ymin,zStartNeutral,ymax);
  lineStartNeutral->SetLineColor(kGray+1);
  lineStartNeutral->SetLineStyle(1);
  lineStartNeutral->SetLineWidth(2);

  //  cout << "End plasma = " << zEndNeutral << endl;
  TLine *lineEndNeutral = new TLine(zEndNeutral,ymin,zEndNeutral,ymax);
  lineEndNeutral->SetLineColor(kGray+1);
  lineEndNeutral->SetLineStyle(2);
  lineEndNeutral->SetLineWidth(2);
  
  // Plotting
  // -----------------------------------------------

  // Canvas setup
  TCanvas *C;
  if(opt.Contains("hres") && !opt.Contains("pdf")) // high resolution for plain grahics output.
    C = new TCanvas("C","2D Charge density, Accelerating and focusing fields",1126,1000);
  else
    C = new TCanvas("C","2D Charge density, Accelerating and focusing fields",750,666);
 

  // Palettes setup
  TExec *exPlasma = new TExec("exPlasma","plasmaPalette->cd();");
  TExec *exHot    = new TExec("exHot","hotPalette->cd();");
  TExec *exElec   = new TExec("exElec","electronPalette->cd();");
  TExec *exField  = new TExec("exField","rbow2Palette->cd();");
    
  
  // Actual Plotting!
  // ------------------------------------------------------------

  // Output file
  TString fOutName = Form("./%s/Plots/ChargeField2D/ChargeField2D",pData->GetPath().c_str());
  fOutName += Form("-%s_%i",pData->GetName(),time);
 
  // Setup Pad layout:
  Float_t lMargin = 0.15;
  Float_t rMargin = 0.18;
  Float_t bMargin = 0.12;
  Float_t tMargin = 0.04;
  PlasmaGlob::CanvasPartition(C,2,lMargin,rMargin,bMargin,tMargin);

  TPad *pad[2];

  
  // Draw!
  C->cd(0);
  char pname[16];
  sprintf(pname,"pad_%i",1);
  pad[0] = (TPad*) gROOT->FindObject(pname);
  pad[0]->Draw();
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

  hFrame->GetXaxis()->SetLabelOffset(999);

  hFrame->GetYaxis()->SetTitleSize(0.075);
  hFrame->GetYaxis()->SetTitleOffset(0.8);
  hFrame->GetYaxis()->SetLabelSize(0.065);
  hFrame->GetYaxis()->SetLabelOffset(0.02);
  hFrame->GetYaxis()->SetTickLength(0.02);
 
  hFrame->GetZaxis()->SetTitleSize(0.06);                        
  hFrame->GetZaxis()->SetTitleOffset(0.45);
  hFrame->GetZaxis()->SetLabelSize(0.06);  
  hFrame->GetZaxis()->SetTickLength(0.02);
  //  hFrame->GetZaxis()->SetNdivisions(505);

  hFrame->Draw("col");            

  if(Nspecies>=3) {
    if(hDen2D[2]) {
      exHot->Draw();
      hDen2D[2]->Draw("colz same");
    }
  }

  exPlasma->Draw();
  hDen2D[0]->Draw("colz same");
  
  if(hDen2D[1]) {
    exElec->Draw();
    hDen2D[1]->Draw("colz same");
  }
  
  if(opt.Contains("1dline")) {
    lineYzero->Draw();
    lineYdown->Draw();
    lineYup->Draw();
  }

  if(zStartPlasma>xmin && zStartPlasma<xmax)
    lineStartPlasma->Draw();
  if(zStartNeutral>xmin && zStartNeutral<xmax)
    lineStartNeutral->Draw();
  if(zEndNeutral>xmin && zEndNeutral<xmax)
    lineEndNeutral->Draw();
  
  pad[0]->Update();
  TPaletteAxis *palette = NULL;
  for(Int_t i=0;i<Nspecies;i++) {
    
    if(!hDen2D[i]) continue;
    palette = (TPaletteAxis*) hDen2D[i]->GetListOfFunctions()->FindObject("palette");
    if(!palette) continue;
    
    Float_t y1 = gPad->GetBottomMargin();
    Float_t y2 = 1 - gPad->GetTopMargin();
    Float_t x1 = gPad->GetLeftMargin();
    Float_t x2 = 1 - gPad->GetRightMargin();
    palette->SetY2NDC( (Nspecies-i)*(y2-y1)/Nspecies + y1);
    palette->SetY1NDC( (Nspecies-i-1)*(y2-y1)/Nspecies + y1);
    palette->SetX1NDC(x2 + 0.005);
    palette->SetX2NDC(x2 + 0.03);
    palette->SetTitleOffset(0.8);
    palette->SetTitleSize(0.07);
    palette->SetLabelSize(0.065);
    palette->SetBorderSize(2);
    palette->SetLineColor(1);
  
  }


  // 1D charge density plots:
  Float_t yaxismin  =  pad[0]->GetUymin();
  Float_t yaxismax  =  pad[0]->GetUymin() + (1.0/Nspecies)*(pad[0]->GetUymax() - pad[0]->GetUymin()) - 0.00;
  
  //  Float_t denmin = (0.1001) * density;
  Float_t denmin = Min[1];
  Float_t denmax = Max[1];
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

    hDen1D[i]->SetLineWidth(2);
    if(i==1) {
      hDen1D[i]->SetLineColor(PlasmaGlob::elecLine);
      hDen1D[i]->Draw("same C");
    } else if (i==2) {
      hDen1D[i]->SetLineColor(kOrange+8);
      hDen1D[i]->Draw("same C");
    }
  }

  // Text objects
  Float_t y1 = gPad->GetBottomMargin();
  Float_t y2 = 1 - gPad->GetTopMargin();
  Float_t x1 = gPad->GetLeftMargin();
  Float_t x2 = 1 - gPad->GetRightMargin();
  
  TPaveText *textTime = new TPaveText(x2-0.17,y2-0.12,x2-0.02,y2-0.02,"NDC");
  PlasmaGlob::SetPaveTextStyle(textTime,32); 
  char ctext[128];
  if(opt.Contains("units") && n0) 
    sprintf(ctext,"z = %5.1f #mum", Time * skindepth / PUnits::um);
  else
    sprintf(ctext,"t = %5.1f #omega_{p}^{-1}",Time);
  textTime->AddText(ctext);
  
  TPaveText *textDen = new TPaveText(x1+0.02,y2-0.12,x1+0.22,y2-0.02,"NDC");
  PlasmaGlob::SetPaveTextStyle(textDen,12); 
  textDen->SetTextColor(kOrange+10);
  if(opt.Contains("units") && n0)
    sprintf(ctext,"n_{0} = %5.2f x 10^{17} / cm^{3}", 1e-17 * n0 * PUnits::cm3);
  else if(pData->GetBeamDensity() && n0)
    sprintf(ctext,"n_{b}/n_{0} = %5.2f", pData->GetBeamDensity()/n0);
  textDen->AddText(ctext);

    
  textTime->Draw();
  textDen->Draw();
    
  pad[0]->RedrawAxis(); 

  C->cd(0);
  sprintf(pname,"pad_%i",0);
  pad[1] = (TPad*) gROOT->FindObject(pname);
  pad[1]->Draw();
  pad[1]->cd(); // <---------------------------------------------------------- bottom Plot
  pad[1]->SetFrameLineWidth(3);  
    
  TH2F *hFrame2 = (TH2F*) gROOT->FindObject("hFrame2");
  if(hFrame2) delete hFrame2;
  hFrame2 = (TH2F*) hE2D[0]->Clone("hFrame2");
  hFrame2->Reset();

  Float_t yFactor = pad[0]->GetAbsHNDC()/pad[1]->GetAbsHNDC(); 

  hFrame2->GetXaxis()->SetTitleSize(0.07);
  hFrame2->GetXaxis()->SetLabelSize(0.06);
  hFrame2->GetXaxis()->SetLabelOffset(0.03);
 

  hFrame2->GetYaxis()->SetTitleSize(0.075*yFactor);
  hFrame2->GetYaxis()->SetTitleOffset(0.8/yFactor);
  hFrame2->GetYaxis()->SetLabelSize(0.0650*yFactor);
  hFrame2->GetYaxis()->SetLabelOffset(0.02/yFactor);
  hFrame2->GetYaxis()->SetTickLength(0.02/yFactor);
  
  hE2D[0]->GetZaxis()->SetTitleSize(0.06*yFactor);                        
  hE2D[0]->GetZaxis()->SetTitleOffset(0.45/yFactor);
  hE2D[0]->GetZaxis()->SetLabelSize(0.06*yFactor);  
  hE2D[0]->GetZaxis()->SetTickLength(0.02/yFactor);
  // hFrame2->GetZaxis()->SetNdivisions(505);
  
  hFrame2->Draw("col");

  exField->Draw();
  hE2D[0]->Draw("colz same");

  if(opt.Contains("1dline")) {
    lineYzero->Draw();
    lineYdown->Draw();
    lineYup->Draw();
  }
  if(zStartPlasma>xmin && zStartPlasma<xmax)
    lineStartPlasma->Draw();
  if(zStartNeutral>xmin && zStartNeutral<xmax)
    lineStartNeutral->Draw();
  if(zEndNeutral>xmin && zEndNeutral<xmax)
    lineEndNeutral->Draw();
  
  pad[1]->Update();
  
  palette = (TPaletteAxis*) hE2D[0]->GetListOfFunctions()->FindObject("palette");
  if(palette) {
    Float_t y1 = gPad->GetBottomMargin();
    Float_t y2 = 1 - gPad->GetTopMargin();
    Float_t x1 = gPad->GetLeftMargin();
    Float_t x2 = 1 - gPad->GetRightMargin();
    palette->SetY2NDC(y2 - 0.01);
    palette->SetY1NDC(y1 + 0.01);
    palette->SetX1NDC(x2 + 0.005);
    palette->SetX2NDC(x2 + 0.03);
    palette->SetTitleSize(0.07*yFactor);
    palette->SetTitleOffset(0.8/yFactor);
    palette->SetLabelSize(0.065*yFactor);
    palette->SetBorderSize(2);
    palette->SetLineColor(1);
  }
  
  // 1D electric field plots:
  if(hE1D[0]) {
    // yaxismin  =  pad[1]->GetUymin() + 0.16666*(pad[1]->GetUymax() - pad[1]->GetUymin());
    // yaxismax  =  pad[1]->GetUymax() - 0.16666*(pad[1]->GetUymax() - pad[1]->GetUymin());
    Float_t yaxismin  =  pad[1]->GetUymin();
    Float_t yaxismax  =  pad[1]->GetUymax();
    
    Float_t emin = hE1D[0]->GetMinimum();
    Float_t emax = hE1D[0]->GetMaximum();
    
    //Float_t slope = (yaxismax - yaxismin)/(emax - emin);
    Float_t slope = yaxismin/emin;
    
    for(Int_t j=0;j<hE1D[0]->GetNbinsX();j++) {
      Float_t content = hE1D[0]->GetBinContent(j+1);
      hE1D[0]->SetBinContent(j+1,(content - emin) * slope + yaxismin);
  }    
    
    hE1D[0]->SetLineWidth(2);
    hE1D[0]->SetLineColor(PlasmaGlob::elecLine);
    // hE1D[0]->SetLineColor(kGray+2);
    hE1D[0]->Draw("same C");
  }

  // Plot ionization probability:
  if(hIonProb1D) {
    Float_t ionmin = 0;
    Float_t ionmax = hIonProb1D->GetMaximum();
    
    Float_t slope = yaxismax/ionmax;
    
    for(Int_t j=0;j<hIonProb1D->GetNbinsX();j++) {
      Float_t content = hIonProb1D->GetBinContent(j+1);
      hIonProb1D->SetBinContent(j+1,content * slope);
    }    
    
    hIonProb1D->SetLineWidth(2);
    hIonProb1D->SetLineColor(kGray+2);
    hIonProb1D->Draw("same C");
  }
  
  pad[1]->RedrawAxis(); 

  C->cd();

  // Print to a file
  PlasmaGlob::imgconv(C,fOutName,opt);
  // ---------------------------------------------------------
 
}
