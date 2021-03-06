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
#include <TClonesArray.h>

#include "PData.hh"
#include "PlasmaGlob.hh"
#include "PPalette.hh"

void PlotChargeFieldFocus2D( const TString &sim, Int_t time, Int_t zoom=2, Int_t Nbins=2, const TString &options="") {
  
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
  if(opt.Contains("gridx")) {
    gStyle->SetPadGridX(1);
  }
  if(opt.Contains("gridy")) {
    gStyle->SetPadGridY(1);
  }
 
  
  // Some plasma constants
  Double_t n0 = pData->GetPlasmaDensity();
  Double_t omegap = pData->GetPlasmaFrequency();
  Double_t timedepth = 1.;
  if(omegap!=0.0) timedepth = 1/omegap;
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
  cout << endl;

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
    if(rms0>0.0)
      Nbins =  TMath::Nint(rms0 / pData->GetDX(1));
    else
      Nbins = 1;
  }
  
  // Slice width limits.
  if(!pData->IsCyl()) {
    FirstyBin = pData->GetNX(1)/2 + 1 - Nbins;
    LastyBin =  pData->GetNX(1)/2 + Nbins;
  } else {
    FirstyBin = 1; 
    LastyBin  = Nbins;
  }
  // ----------------------------------------------------------------------------------
  
  // Get charge density histos
  Int_t Nspecies = pData->NSpecies();
  TH2F **hDen2D = new TH2F*[Nspecies];
  // Get charge density on-axis
  TH1F **hDen1D = new TH1F*[Nspecies];
  // And electric current (integrated)
  TH1F **hCur1D = new TH1F*[Nspecies];
  for(Int_t i=0;i<Nspecies;i++) {
     
    hDen2D[i] = NULL;
    
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
    
    hDen2D[i]->SetName(hName);
    hDen2D[i]->GetXaxis()->CenterTitle();
    hDen2D[i]->GetYaxis()->CenterTitle();
    hDen2D[i]->GetZaxis()->CenterTitle();
    
    if(opt.Contains("comov"))
      hDen2D[i]->GetXaxis()->SetTitle("k_{p} #zeta");
    else
      hDen2D[i]->GetXaxis()->SetTitle("k_{p} z");
    
    if(pData->IsCyl()) 
      hDen2D[i]->GetYaxis()->SetTitle("k_{p} r");
    else
      hDen2D[i]->GetYaxis()->SetTitle("k_{p} y");
    
    hDen2D[i]->GetZaxis()->SetTitle("n [n_{0}]");
    
    
    hDen1D[i] = NULL;
    hCur1D[i] = NULL;
    
    if(!pData->GetEfieldFileName(i))
      continue;
    
    sprintf(hName,"hDen1D_%i",i);
    hDen1D[i] = (TH1F*) gROOT->FindObject(hName);
    if(hDen1D[i]) delete hDen1D[i];
      
    // 1D histograms
    if(pData->Is3D()) {
      hDen1D[i] = pData->GetH1SliceZ3D(pData->GetChargeFileName(i)->c_str(),"charge",-1,Nbins,-1,Nbins,opt+"avg");
    } else if(pData->IsCyl()) { // Cylindrical: The first bin with r>0 is actually the number 1 (not the 0).
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

    // Get the current:
    if(i==0) continue;

    sprintf(hName,"hCur1D_%i",i);
    hCur1D[i] = (TH1F*) gROOT->FindObject(hName);
    if(hCur1D[i]) delete hCur1D[i];
      
    if(opt.Contains("curr")) {
      // To get the current is needed to read in a wider transverse range which includes all the charge.
      Int_t NbinsT = 100;
      if(pData->Is3D()) {
	hCur1D[i] = pData->GetH1SliceZ3D(pData->GetChargeFileName(i)->c_str(),"charge",-1,NbinsT,-1,NbinsT,opt+"int");
      } else if(pData->IsCyl()) { // Cylindrical: The first bin with r>0 is actually the number 1 (not the 0).
	hCur1D[i] = pData->GetH1SliceZ(pData->GetChargeFileName(i)->c_str(),"charge",1,NbinsT,opt+"int");
      } else { // 2D cartesian
	hCur1D[i] = pData->GetH1SliceZ(pData->GetChargeFileName(i)->c_str(),"charge",-1,NbinsT,opt+"int");
      }
      hCur1D[i]->SetName(hName); 
      
      if(opt.Contains("comov")) {
	hCur1D[i]->GetXaxis()->SetTitle("#zeta [c/#omega_{p}]");
	hCur1D[i]->GetYaxis()->SetTitle("dn/d#zeta [(n_{0}/k_{p}^{3}) (#omega_{p}/c)]");  
      } else {
	hCur1D[i]->GetXaxis()->SetTitle("z [c/#omega_{p}]");
	hCur1D[i]->GetYaxis()->SetTitle("dn/dz [(n_{0}/k_{p}^{3}) (#omega_{p}/c)]");  
      }
      
      Int_t NB = hCur1D[i]->GetNbinsX();
      Float_t dx = (hCur1D[i]->GetBinLowEdge(1)-hCur1D[i]->GetBinLowEdge(NB+1))/NB;
      
      // hCur1D[i]->Scale(dx);
      Float_t Charge = hCur1D[i]->Integral() * dx; 
      
      cout << Form(" Integrated charge of specie %3i = %8.4f n0 * kp^-3",i,Charge) << endl;
    }
  }

  // Get electric fields 2D
  const Int_t Nfields = 2;
  TH2F **hE2D = new TH2F*[Nfields];
  TH1F **hE1D = new TH1F*[Nfields];
  TH1F *hV1D = NULL;
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
      hE2D[i]->GetXaxis()->SetTitle("k_{p} #zeta");
    else
      hE2D[i]->GetXaxis()->SetTitle("k_{p} z");
    
    if(pData->IsCyl()) 
      hE2D[i]->GetYaxis()->SetTitle("k_{p} r");
    else
      hE2D[i]->GetYaxis()->SetTitle("k_{p} y");
    
    if(i==0)
      hE2D[i]->GetZaxis()->SetTitle("E_{z}/E_{0}");
    else if(i==1)
      hE2D[i]->GetZaxis()->SetTitle("E_{y}/E_{0}");
    else if(i==2)
      hE2D[i]->GetZaxis()->SetTitle("E_{x}/E_{0}");
    
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
    
    // Alternative
    // if(hE1D[i]) delete hE1D[i];
    // hE1D[i] = (TH1F*) hE2D[i]->ProjectionX(hName,FirstyBin,LastyBin);
    // hE1D[i]->Scale(1.0/(LastyBin-FirstyBin+1));
    
    if(i==0) {
      Int_t   NbinsX = hE2D[i]->GetNbinsX();
      Int_t   NbinsY = hE2D[i]->GetNbinsY();

      Float_t dx = pData->GetDX(0);
     
      sprintf(hName,"hV1D");
      hV1D = (TH1F*) hE1D[i]->Clone(hName);
      hV1D->Reset();

      Double_t integral = 0.0;
      for(Int_t k=NbinsX;k>0;k--) {
	integral += hE1D[i]->GetBinContent(k) * dx;
	hV1D->SetBinContent(k,integral);
      }

    }
  
  }

  // Get magnetic fields 2D
  const Int_t NBfields = 3;
  TH2F **hB2D = new TH2F*[NBfields];
  for(Int_t i=0;i<NBfields;i++) {
    hB2D[i] = NULL;
    if(i<2) continue;  // Just get the third component.

    if(!pData->GetBfieldFileName(i))
      continue;

    cout << Form(" Getting magnetic field number ") << i+1 << endl;
    
    char hName[24];
    sprintf(hName,"hB2D_%i",i);
    hB2D[i] = (TH2F*) gROOT->FindObject(hName);
    if(hB2D[i]) delete hB2D[i];
    
    if(!pData->Is3D())
      hB2D[i] = pData->GetBField(i,opt);
    else
      hB2D[i] = pData->GetBField2DSliceZY(i,-1,Nbins,opt+"avg");
    
    hB2D[i]->SetName(hName);   
    hB2D[i]->GetXaxis()->CenterTitle();
    hB2D[i]->GetYaxis()->CenterTitle();
    hB2D[i]->GetZaxis()->CenterTitle();
    if(opt.Contains("comov"))
      hB2D[i]->GetXaxis()->SetTitle("k_{p} #zeta");
    else
      hB2D[i]->GetXaxis()->SetTitle("k_{p} z");
    
    if(pData->IsCyl()) 
      hB2D[i]->GetYaxis()->SetTitle("k_{p} r");
    else
      hB2D[i]->GetYaxis()->SetTitle("k_{p} y");
    
    if(i==0)
      hB2D[i]->GetZaxis()->SetTitle("B_{z}/E_{0}");
    else if(i==1)
      hB2D[i]->GetZaxis()->SetTitle("B_{y}/E_{0}");
    else if(i==2)
      hB2D[i]->GetZaxis()->SetTitle("B_{x}/E_{0}");
      
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
      Float_t ymin = skindepth * hDen2D[i]->GetYaxis()->GetXmin() / PUnits::um;
      Float_t ymax = skindepth * hDen2D[i]->GetYaxis()->GetXmax() / PUnits::um;
      hDen2D[i]->SetBins(NbinsX,xMin,xMax,NbinsY,ymin,ymax);
      // for(Int_t j=0;j<hDen2D[i]->GetNbinsX();j++) {
      // 	for(Int_t k=0;k<hDen2D[i]->GetNbinsY();k++) {
      // 	  hDen2D[i]->SetBinContent(j,k, hDen2D[i]->GetBinContent(j,k) * n0 / (1e17/PUnits::cm3) );
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
      
      // if(i==0)
      // 	hDen2D[i]->GetZaxis()->SetTitle("n_{e} [10^{17}/cm^{3}]");
      // else if(i==1)
      // 	hDen2D[i]->GetZaxis()->SetTitle("n_{b} [10^{17}/cm^{3}]"); 
      // else
      // 	hDen2D[i]->GetZaxis()->SetTitle("n_{i} [10^{17}/cm^{3}]"); 
      
      hDen1D[i]->SetBins(NbinsX,xMin,xMax);
      // for(Int_t j=0;j<hDen1D[i]->GetNbinsX();j++) {
      // 	hDen1D[i]->SetBinContent(j, hDen1D[i]->GetBinContent(j) * n0 / (1e17/PUnits::cm3) );
      // }
      
      
      if(opt.Contains("comov"))
    	hDen1D[i]->GetXaxis()->SetTitle("#zeta [#mum]");
      else
    	hDen1D[i]->GetXaxis()->SetTitle("z [#mum]");
      
      if(hCur1D[i]) {

	hCur1D[i]->SetBins(NbinsX,xMin,xMax);
	Double_t binSize = (xMax - xMin)/NbinsX;  // bin size in um.
	
	Double_t dV = skindepth * skindepth * skindepth;
	Double_t lightspeed =  PConst::c_light / (PUnits::um/PUnits::femtosecond);
	
	hCur1D[i]->Scale(TMath::Abs(n0 * dV * (PConst::ElectronCharge/PUnits::picocoulomb) * (kp * PConst::c_light * PUnits::femtosecond)));
	
	hCur1D[i]->GetYaxis()->SetTitle("I[kA]");
	hCur1D[i]->GetYaxis()->SetTitle("");
	if(opt.Contains("comov"))
	  hCur1D[i]->GetXaxis()->SetTitle("#zeta [#mum]");
	else
	  hCur1D[i]->GetXaxis()->SetTitle("z [#mum]");
	
	
	Float_t Charge = hCur1D[i]->Integral() * (binSize / lightspeed); 
	cout << Form(" Integrated charge of specie %3i = %8f pC",i,Charge) << endl;
      }
    }

    
    for(Int_t i=0;i<Nfields;i++) {
      Int_t NbinsX = hE2D[i]->GetNbinsX();
      Float_t xMin = skindepth * hE2D[i]->GetXaxis()->GetXmin() / PUnits::um;
      Float_t xMax = skindepth * hE2D[i]->GetXaxis()->GetXmax() / PUnits::um;
      Int_t NbinsY = hE2D[i]->GetNbinsY();
      Float_t ymin = skindepth * hE2D[i]->GetYaxis()->GetXmin() / PUnits::um;
      Float_t ymax = skindepth * hE2D[i]->GetYaxis()->GetXmax() / PUnits::um;
      hE2D[i]->SetBins(NbinsX,xMin,xMax,NbinsY,ymin,ymax);
      hE1D[i]->SetBins(NbinsX,xMin,xMax);
            
      for(Int_t j=0;j<hE2D[i]->GetNbinsX();j++) {
	for(Int_t k=0;k<hE2D[i]->GetNbinsY();k++) {
	  hE2D[i]->SetBinContent(j,k, hE2D[i]->GetBinContent(j,k) * ( E0 / (PUnits::GV/PUnits::m) ) );
	}
	hE1D[i]->SetBinContent(j, hE1D[i]->GetBinContent(j) * ( E0 / (PUnits::GV/PUnits::m) ) );
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
      
      
      if(i==0) {    
	hV1D->SetBins(NbinsX,xMin,xMax);

	for(Int_t j=0;j<NbinsX;j++) {
	  hV1D->SetBinContent(j, hV1D->GetBinContent(j) * ( E0 * skindepth / (PUnits::MV) ) );
	}
	
   
	if(opt.Contains("comov")) {
	  hV1D->GetXaxis()->SetTitle("#zeta [#mum]");
	} else {
	  hV1D->GetXaxis()->SetTitle("z [#mum]");
	}
	 
	hV1D->GetYaxis()->SetTitle("#Psi [MV]");
      }
    }

    for(Int_t i=0;i<NBfields;i++) {
      if(!hB2D[i]) continue;
      
      Int_t NbinsX = hB2D[i]->GetNbinsX();
      Float_t xMin = skindepth * hB2D[i]->GetXaxis()->GetXmin() / PUnits::um;
      Float_t xMax = skindepth * hB2D[i]->GetXaxis()->GetXmax() / PUnits::um;
      Int_t NbinsY = hB2D[i]->GetNbinsY();
      Float_t ymin = skindepth * hB2D[i]->GetYaxis()->GetXmin() / PUnits::um;
      Float_t ymax = skindepth * hB2D[i]->GetYaxis()->GetXmax() / PUnits::um;
      hB2D[i]->SetBins(NbinsX,xMin,xMax,NbinsY,ymin,ymax);
      
      for(Int_t j=0;j<hB2D[i]->GetNbinsX();j++) {
	for(Int_t k=0;k<hB2D[i]->GetNbinsY();k++) {
	  hB2D[i]->SetBinContent(j,k, hB2D[i]->GetBinContent(j,k) * ( E0 / (PUnits::GV/PUnits::m) ) );
	}

      }
      
      if(pData->IsCyl())
	hB2D[i]->GetYaxis()->SetTitle("r [#mum]");      
      else
	hB2D[i]->GetYaxis()->SetTitle("y [#mum]");      
      
      if(opt.Contains("comov"))
	hB2D[i]->GetXaxis()->SetTitle("#zeta [#mum]");
      else
	hB2D[i]->GetXaxis()->SetTitle("z [#mum]");
      
      if(i==0)
	hB2D[i]->GetZaxis()->SetTitle("B_{z} [GV/m]");
      else if(i==1)
	hB2D[i]->GetZaxis()->SetTitle("B_{y} [GV/m]");
      else if(i==2)
	hB2D[i]->GetZaxis()->SetTitle("B_{x} [GV/m]");
      
      
    }


  }

  // Focusing field
  TH2F *hFocus2D = (TH2F*) hE2D[1]->Clone("hFocus2D");
  if(hB2D[2]) {
    hFocus2D->Add(hB2D[2],-1);
    if(opt.Contains("units")) {
      hFocus2D->GetZaxis()->SetTitle("E_{y}-B_{x} [GV/m]");
    } else {
      hFocus2D->GetZaxis()->SetTitle("E_{y}-B_{x} [E_{0}]");
    }
  }
  

  // --------------------------------------------------- Vertical Zoom ------------
  
  Float_t yRange    = (hDen2D[0]->GetYaxis()->GetXmax() - hDen2D[0]->GetYaxis()->GetXmin())/zoom;
  Float_t midPoint = (hDen2D[0]->GetYaxis()->GetXmax() + hDen2D[0]->GetYaxis()->GetXmin())/2.;
  Float_t yMin = midPoint-yRange/2;
  Float_t yMax = midPoint+yRange/2;
  if(pData->IsCyl()) {
    yMin = pData->GetXMin(1);
    yMax = yRange;
  }

  for(Int_t i=0;i<Nspecies;i++) {
    if(!hDen2D[i]) continue;
    hDen2D[i]->GetYaxis()->SetRangeUser(yMin,yMax);
  }
  
  for(Int_t i=0;i<Nfields;i++) {
    if(!hE2D[i]) continue;
    hE2D[i]->GetYaxis()->SetRangeUser(yMin,yMax);
  }
  
  hFocus2D->GetYaxis()->SetRangeUser(yMin,yMax);


  Float_t xMin = hDen2D[0]->GetXaxis()->GetXmin();
  Float_t xMax = hDen2D[0]->GetXaxis()->GetXmax();
  Float_t xRange = xMax - xMin;
  
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
    if(i==1) Min[i] = 1.01E-1 * Base;
    if(i==2) Min[i] = 1.01E-4 * Base;
    hDen2D[i]->GetZaxis()->SetRangeUser(Min[i],Max[i]);
  }
  
  // Dynamic plasma palette
  const Int_t plasmaDNRGBs = 3;
  const Int_t plasmaDNCont = 64;
  Float_t basePos = 0.5;
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
  Float_t *Emax = new Float_t[Nfields];
  Float_t *Emin = new Float_t[Nfields];
  for(Int_t i=0;i<Nfields;i++) {
    Emax[i] = hE2D[i]->GetMaximum();
    Emin[i] = hE2D[i]->GetMinimum();
    if(Emax[i] > TMath::Abs(Emin[i]))
      Emin[i] = -Emax[i];
    else
      Emax[i] = -Emin[i];
    hE2D[i]->GetZaxis()->SetRangeUser(Emin[i],Emax[i]); 
  }
  
  Float_t Fmax = hFocus2D->GetMaximum();
  Float_t Fmin = hFocus2D->GetMinimum();
  if(Fmax > TMath::Abs(Fmin))
    Fmin = -Fmax;
  else
    Fmax = -Fmin;
  hFocus2D->GetZaxis()->SetRangeUser(Fmin,Fmax); 
  
  // Potential 
  if(opt.Contains("units")) {
    trapPotential *=  ( E0 * skindepth / (PUnits::MV) ); 
  }  
  
  Float_t Vmin = hV1D->GetMinimum();    
  { // Shift potential
    Int_t   NbinsX = hV1D->GetNbinsX(); 
    for(Int_t j=0;j<NbinsX;j++) {
      hV1D->SetBinContent(j, hV1D->GetBinContent(j) - Vmin -trapPotential);
    }
  }
  
  Vmin = hV1D->GetMinimum();   
  Float_t Vmax = hV1D->GetMaximum();    

  // "Axis range" in Osiris units:
  Double_t ylow  = hDen2D[0]->GetYaxis()->GetBinLowEdge(FirstyBin);
  Double_t yup = hDen2D[0]->GetYaxis()->GetBinUpEdge(LastyBin);
  Double_t xmin = hDen2D[0]->GetXaxis()->GetXmin();
  Double_t xmax = hDen2D[0]->GetXaxis()->GetXmax();

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
  TLine *lineStartPlasma = new TLine(zStartPlasma,yMin,zStartPlasma,yMax);
  lineStartPlasma->SetLineColor(kGray+2);
  lineStartPlasma->SetLineStyle(2);
  lineStartPlasma->SetLineWidth(3);

  //  cout << "Start plasma = " << zStartNeutral << endl;
  TLine *lineStartNeutral = new TLine(zStartNeutral,yMin,zStartNeutral,yMax);
  lineStartNeutral->SetLineColor(kGray+1);
  lineStartNeutral->SetLineStyle(2);
  lineStartNeutral->SetLineWidth(3);

  //  cout << "End plasma = " << zEndNeutral << endl;
  TLine *lineEndNeutral = new TLine(zEndNeutral,yMin,zEndNeutral,yMax);
  lineEndNeutral->SetLineColor(kGray+1);
  lineEndNeutral->SetLineStyle(2);
  lineEndNeutral->SetLineWidth(3);


  // Plotting
  // -----------------------------------------------

  // Canvas setup
  TCanvas *C = new TCanvas("C","2D Charge density and Electric field",858,1000);
 
  // Palettes setup
  TExec *exPlasma = new TExec("exPlasma","plasmaPalette->cd();");
  TExec *exElec   = new TExec("exElec","redelectronPalette->cd();");
  TExec *exHot    = new TExec("exHot","hotPalette->cd();");
  TExec *exField  = new TExec("exField","rbow2Palette->cd();");
  TExec *exFieldT = new TExec("exFieldT","redPalette->cd();");
  TExec *exIonP   = new TExec("exIonP","redPalette->cd();");
  TExec *exPot    = new TExec("exPot","rbow2invPalette->cd();");
     
  // Actual Plotting!
  // ------------------------------------------------------------

  // Output file
  TString fOutName = Form("./%s/Plots/ChargeFieldFocus2D/ChargeFieldFocus2D",pData->GetPath().c_str());
  fOutName += Form("-%s_%i",pData->GetName(),time);

  // Setup Pad layout:
  Float_t lMargin = 0.15;
  Float_t rMargin = 0.18;
  Float_t bMargin = 0.10;
  Float_t tMargin = 0.04;
  Float_t factor = 1.0;    
  PlasmaGlob::CanvasAsymPartition(C,3,lMargin,rMargin,bMargin,tMargin,factor);
  
  TPad *pad[3];
  TString sLabels[] = {"(a)","(b)","(c)"};
  // Text objects
  TPaveText **textLabel = new TPaveText*[3];

  C->cd(0);
  char pname[16];
  sprintf(pname,"pad_%i",2);
  pad[0] = (TPad*) gROOT->FindObject(pname);
  pad[0]->Draw();
  pad[0]->cd(); // <---------------------------------------------- Top Plot ---------
  if(opt.Contains("logz")) {
    pad[0]->SetLogz(1);
  } else {
    pad[0]->SetLogz(0);
  }
  pad[0]->SetFrameLineWidth(3);  
  pad[0]->SetTickx(1);


  // Re-range:
  for(Int_t i=0;i<Nspecies;i++) {
    if(!hDen2D[i]) continue;
    hDen2D[i]->GetYaxis()->SetRangeUser(yMin -(factor-1)*yRange, yMax);
  }

  
  TH2F *hFrame = (TH2F*) gROOT->FindObject("hFrame1");
  if(hFrame) delete hFrame;
  hFrame = (TH2F*) hDen2D[0]->Clone("hFrame1");
  hFrame->Reset();

  hFrame->SetLabelFont(42,"xyz");
  hFrame->SetTitleFont(42,"xyz");

  hFrame->GetYaxis()->SetNdivisions(505);
  hFrame->GetYaxis()->SetLabelSize(0.085);
  hFrame->GetYaxis()->SetTitleSize(0.09);
  hFrame->GetYaxis()->SetTitleOffset(0.7);
  hFrame->GetYaxis()->SetTickLength(0.02);

  hFrame->GetXaxis()->SetLabelOffset(999.);
  hFrame->GetXaxis()->SetTitleOffset(999.);
  hFrame->GetXaxis()->SetTickLength(0.04);

  // Frame asymmetry:
  hFrame->Draw("col");
  
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

  if(opt.Contains("1dline")) {
    lineYzero->Draw();
    lineYdown->Draw();
    lineYup->Draw();
  }
  
  if(opt.Contains("sline")) {
    if(zStartPlasma>xmin && zStartPlasma<xmax)
      lineStartPlasma->Draw();
    if(zStartNeutral>xmin && zStartNeutral<xmax)
      lineStartNeutral->Draw();
    if(zEndNeutral>xmin && zEndNeutral<xmax)
      lineEndNeutral->Draw();
  }

  // lineYdown->Draw();
  // lineYup->Draw();
   
  // Palettes re-arrangement
  pad[0]->Update();
  Float_t y1 = pad[0]->GetBottomMargin();
  Float_t y2 = 1 - pad[0]->GetTopMargin();
  Float_t x1 = pad[0]->GetLeftMargin();
  Float_t x2 = 1 - pad[0]->GetRightMargin();
  
  TPaletteAxis *palette = NULL;
  if(Nspecies>=3) {
    if(hDen2D[2]) {
      palette = (TPaletteAxis*)hDen2D[2]->GetListOfFunctions()->FindObject("palette");
    }
  }
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
    if(opt.Contains("logz")) 
      palette->SetLabelOffset(0.001);
    else
      palette->SetLabelOffset(0.01);
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
    if(opt.Contains("logz")) 
      palette->SetLabelOffset(0.001);
    else
      palette->SetLabelOffset(0.01);
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
    if(opt.Contains("logz")) 
      palette->SetLabelOffset(0.001);
    else
      palette->SetLabelOffset(0.01);
    palette->SetBorderSize(2);
    palette->SetLineColor(1);
  }

  
  // 1D charge density plots:
  Float_t yaxismin  =  pad[0]->GetUymin();
  Float_t yaxismax  =  pad[0]->GetUymin() + 0.33*(pad[0]->GetUymax() - pad[0]->GetUymin()) - 0.00;
  Float_t denmin = Min[1];
  Float_t denmax = Max[1];
  if(opt.Contains("logz")) {
    denmin = TMath::Log10(denmin);
    denmax = TMath::Log10(denmax);
  }

  Float_t curmin = 0.0;
  Float_t curmax = 0.0;
  if(opt.Contains("curr")) {
    curmin = 0.0;
    curmax = hCur1D[1]->GetMaximum();
    
    cout << Form(" Maximum driver  current = %6.2f kA ", curmax) << endl ;
    if(Nspecies>=3)
      if(hCur1D[2])
	cout << Form(" Maximum witness current = %6.2f kA ", hCur1D[2]->GetMaximum()) << endl ;

    // Round for better plotting
    curmax = 0.1*TMath::Nint(curmax*10);
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

    if(hCur1D[i]) {
      slope = (yaxismax - yaxismin)/(curmax - curmin);
      
      for(Int_t j=0;j<hCur1D[i]->GetNbinsX();j++) {
	Float_t content = hCur1D[i]->GetBinContent(j+1);
	
	if(content<curmin) 
	  hCur1D[i]->SetBinContent(j+1,yaxismin);
	else 
	  hCur1D[i]->SetBinContent(j+1,(content - curmin) * slope + yaxismin);	
      }
      
    }
    
  }
  
  // Plasma on-axis density:
  // hDen1D[0]->SetLineWidth(2);
  // hDen1D[0]->SetLineColor(kGray+1);
  // // // PlasmaGlob::SetH1Style(hDen1D[0],1);
  // hDen1D[0]->Draw("same C");  
  
  
  if(opt.Contains("curr")) {
    hCur1D[1]->SetLineWidth(2);
    hCur1D[1]->SetLineColor(PlasmaGlob::elecLine);
    hCur1D[1]->Draw("same C");
  } else {
    hDen1D[1]->SetLineWidth(2);
    hDen1D[1]->SetLineColor(PlasmaGlob::elecLine);
    //    hDen1D[1]->Draw("same C");
  }
  
  if(Nspecies>=3) {
    if(hDen1D[2]) {
      if(opt.Contains("curr")) {
	hCur1D[2]->SetLineWidth(2);
	hCur1D[2]->SetLineColor(kOrange+8);
	hCur1D[2]->Draw("same C");
      } else {
	hDen1D[2]->SetLineWidth(2);
	hDen1D[2]->SetLineColor(kOrange+8);
	//   hDen1D[2]->Draw("same C");
      }
    }
  }
  
  // Current axis
  TGaxis *axis = NULL;
  if(opt.Contains("curr")) {
    axis = new TGaxis(xMax-xRange/6.0,yMin - (factor-1)*yRange,
		      xMax-xRange/6.0,yaxismax,
		      0.001,curmax,503,"+LS");
    
    axis->SetLineWidth(1);
    axis->SetLineColor(kGray+3);//PlasmaGlob::elecLine);
    axis->SetLabelColor(kGray+3);//PlasmaGlob::elecLine);
    axis->SetLabelSize(0.06);
    axis->SetLabelOffset(0.01);
    axis->SetLabelFont(42);
    axis->SetTitleColor(kGray+3);//PlasmaGlob::elecLine);
    axis->SetTitleSize(0.06);
    axis->SetTitleOffset(0.6);
    axis->SetTitleFont(42);
    axis->SetTickSize(0.03);
    axis->SetTitle("I [kA]");
    axis->CenterTitle();
    axis->SetNdivisions(505);
    
    axis->Draw();
  }

  
  TPaveText *textTime = new TPaveText(xMax - 0.3*xRange, yMax-0.15*yRange, xMax-0.1, yMax-0.05*yRange);
  //x2-0.17,y2-0.12,x2-0.02,y2-0.02,"NDC");
  PlasmaGlob::SetPaveTextStyle(textTime,32); 
  char ctext[128];
  if(opt.Contains("units") && n0) 
    sprintf(ctext,"z = %5.1f #mum", Time * skindepth / PUnits::um);
  else
    sprintf(ctext,"t = %5.1f #omega_{p}^{-1}",Time);
  textTime->SetTextFont(42);
  textTime->AddText(ctext);

  textTime->Draw();
  // textDen->Draw();
  // if(opt.Contains("units"))
  //   textWav->Draw();

  textLabel[0] = new TPaveText(xMin + 0.02*xRange, yMax-0.2*yRange, xMin+0.30*xRange, yMax-0.05*yRange);
  PlasmaGlob::SetPaveTextStyle(textLabel[0],12); 
  textLabel[0]->SetTextFont(42);
  textLabel[0]->AddText(sLabels[0]);
  textLabel[0]->Draw();
  
  
  pad[0]->RedrawAxis(); 

 
  C->cd(0);
  sprintf(pname,"pad_%i",1);
  pad[1] = (TPad*) gROOT->FindObject(pname);
  pad[1]->Draw();
  pad[1]->cd(); // <---------------------------------------------------------- Middle Plot
  // if(opt.Contains("logz")) {
  //   pad[1]->SetLogz(1);
  // } else {
  //   pad[1]->SetLogz(0);
  // }
  pad[1]->SetFrameLineWidth(3);  
  pad[1]->SetTickx(1);
   
  hFrame = (TH2F*) gROOT->FindObject("hFrame2");
  if(hFrame) delete hFrame;
  hFrame = (TH2F*) hE2D[0]->Clone("hFrame2");
  hFrame->Reset();
  
  Float_t yFactor = pad[0]->GetAbsHNDC()/pad[1]->GetAbsHNDC();

  hFrame->GetYaxis()->SetLabelSize(yFactor*0.085);
  hFrame->GetYaxis()->SetTitleSize(yFactor*0.09);
  hFrame->GetYaxis()->SetTitleOffset(0.7/yFactor);
  hFrame->GetYaxis()->SetTickLength(0.02/yFactor);

  hFrame->GetXaxis()->SetLabelOffset(999.);
  hFrame->GetXaxis()->SetTitleOffset(999.);
  hFrame->GetXaxis()->SetTickLength(0.04*yFactor);
  
  hFrame->SetLabelFont(42,"xyz");
  hFrame->SetTitleFont(42,"xyz");

  hFrame->Draw("col");

  //  hE2D[0]->GetZaxis()->SetNdivisions(505);
  hE2D[0]->GetZaxis()->SetTitleFont(42);
  hE2D[0]->GetZaxis()->SetTickLength(0.02);
  
  exField->Draw();
  hE2D[0]->Draw("colz same");

  lineYzero->Draw();
  if(opt.Contains("1dline")) {
    lineYzero->Draw();
    lineYdown->Draw();
    lineYup->Draw();
  }

  if(opt.Contains("sline")) {
    if(zStartPlasma>xmin && zStartPlasma<xmax)
      lineStartPlasma->Draw();
    if(zStartNeutral>xmin && zStartNeutral<xmax)
      lineStartNeutral->Draw();
    if(zEndNeutral>xmin && zEndNeutral<xmax)
      lineEndNeutral->Draw();
  }
  
  // Fit the V1D in the E2D pad:
  // Float_t rightmin = Vmin;
  // Float_t rightmax = Vmax;
  Float_t rightmin = Vmin;
  Float_t rightmax = Vmax;
  Float_t slope = (yMax/2)/rightmax;
  
  for(Int_t j=0;j<hV1D->GetNbinsX();j++) {
    hV1D->SetBinContent(j+1,hV1D->GetBinContent(j+1)*slope);
  }
  hV1D->SetLineStyle(1);
  hV1D->SetLineWidth(2);
  hV1D->SetLineColor(PlasmaGlob::elecLine);

  if(opt.Contains("pot"))
    hV1D->Draw("sameC");
  
  // Line for trapping potential
  TLine *lineTrap = new TLine(hV1D->GetXaxis()->GetXmin(),
			      -(trapPotential)*slope,
			      hV1D->GetXaxis()->GetXmax(),
			      -(trapPotential)*slope);
  lineTrap->SetLineColor(PlasmaGlob::elecLine);
  lineTrap->SetLineStyle(2);
  lineTrap->SetLineWidth(1);
  lineTrap->Draw();
  
  // Fit the E1D in the E2D pad:
  rightmin = Emin[0];
  rightmax = Emax[0];
  slope = yMax/rightmax;
  
  for(Int_t j=0;j<hE1D[0]->GetNbinsX();j++) {
    hE1D[0]->SetBinContent(j+1,hE1D[0]->GetBinContent(j+1)*slope);
  }
  hE1D[0]->SetLineStyle(1);
  hE1D[0]->SetLineWidth(2);
  hE1D[0]->SetLineColor(kOrange+10);

  hE1D[0]->Draw("sameC");
  
  //lineYdown->Draw();
  //lineYup->Draw();

  pad[1]->Update();
  
  y1 = pad[1]->GetBottomMargin();
  y2 = 1 - pad[1]->GetTopMargin();
  x1 = pad[1]->GetLeftMargin();
  x2 = 1 - pad[1]->GetRightMargin();
  
  palette = (TPaletteAxis*) hE2D[0]->GetListOfFunctions()->FindObject("palette");  
  if(palette) {
    palette->SetY2NDC(y2 - 0.00);
    palette->SetY1NDC(y1 + 0.00);
    palette->SetX1NDC(x2 + 0.005);
    palette->SetX2NDC(x2 + 0.03);
    palette->SetLabelSize(0.075*yFactor);
    palette->SetLabelFont(42);
    palette->SetLabelOffset(0.01/yFactor);
    palette->SetTitleSize(0.075*yFactor);
    palette->SetTitleOffset(0.80/yFactor);
    palette->SetBorderSize(2);
    palette->SetLineColor(1);
  }
   
  pad[1]->RedrawAxis(); 

  
  if(opt.Contains("1dline")) {
    lineYzero->Draw();
  }

  if(opt.Contains("sline")) {
    if(zStartPlasma>xmin && zStartPlasma<xmax)
      lineStartPlasma->Draw();
    if(zStartNeutral>xmin && zStartNeutral<xmax)
      lineStartNeutral->Draw();
    if(zEndNeutral>xmin && zEndNeutral<xmax)
      lineEndNeutral->Draw();
  }
 

  textLabel[1] = new TPaveText(xMin + 0.02*xRange, yMax-0.2*yRange, xMin+0.30*xRange, yMax-0.05*yRange);
  PlasmaGlob::SetPaveTextStyle(textLabel[1],12); 
  textLabel[1]->SetTextFont(42);
  textLabel[1]->AddText(sLabels[1]);
  textLabel[1]->Draw();


  C->cd(0);
  sprintf(pname,"pad_%i",0);
  pad[2] = (TPad*) gROOT->FindObject(pname);
  pad[2]->Draw();
  pad[2]->cd(); // <--------------------------------------------------------- Bottom Plot
  pad[2]->SetFrameLineWidth(3);  
  pad[2]->SetTickx(1);

  hFrame = (TH2F*) gROOT->FindObject("hFrame3");
  if(hFrame) delete hFrame;
  hFrame = (TH2F*) hE2D[0]->Clone("hFrame3");
  hFrame->Reset();

  yFactor = pad[0]->GetAbsHNDC()/pad[2]->GetAbsHNDC();
  
  
  hFrame->GetYaxis()->SetLabelSize(yFactor*0.085);
  hFrame->GetYaxis()->SetTitleSize(yFactor*0.09);
  hFrame->GetYaxis()->SetTitleOffset(0.7/yFactor);
  hFrame->GetYaxis()->SetTickLength(0.02/yFactor);

  hFrame->GetXaxis()->SetTitleSize(0.10);
  hFrame->GetXaxis()->SetLabelSize(0.08);
  hFrame->GetXaxis()->SetLabelOffset(0.02);
  hFrame->GetXaxis()->SetTitleOffset(1.0);
  hFrame->GetXaxis()->SetTickLength(0.04*yFactor);
  
  hFrame->SetLabelFont(42,"xyz");
  hFrame->SetTitleFont(42,"xyz");
  
  hFrame->Draw("col");

  //  hFocus2D->GetZaxis()->SetNdivisions(505);
  hFocus2D->GetZaxis()->SetTitleFont(42);
  hFocus2D->GetZaxis()->SetTickLength(0.02/yFactor);

  exField->Draw();
  hFocus2D->Draw("colz same");


  if(opt.Contains("1dline")) {
    lineYzero->Draw();
    lineYdown->Draw();
    lineYup->Draw();
  }

  if(opt.Contains("sline")) {
    if(zStartPlasma>xmin && zStartPlasma<xmax)
      lineStartPlasma->Draw();
    if(zStartNeutral>xmin && zStartNeutral<xmax)
      lineStartNeutral->Draw();
    if(zEndNeutral>xmin && zEndNeutral<xmax)
      lineEndNeutral->Draw();
  }
  
    
  pad[2]->Update();

  y1 = pad[2]->GetBottomMargin();
  y2 = 1 - pad[2]->GetTopMargin();
  x1 = pad[2]->GetLeftMargin();
  x2 = 1 - pad[2]->GetRightMargin();
 
  palette = (TPaletteAxis*)hFocus2D->GetListOfFunctions()->FindObject("palette");
  if(palette) {
    palette->SetY2NDC(y2 - 0.00);
    palette->SetY1NDC(y1 + 0.00);
    palette->SetX1NDC(x2 + 0.005);
    palette->SetX2NDC(x2 + 0.03);
    palette->SetLabelSize(yFactor*0.075);
    palette->SetLabelFont(42);
    palette->SetLabelOffset(0.01/yFactor);
    palette->SetTitleSize(yFactor*0.075);
    palette->SetTitleOffset(0.80/yFactor);
    palette->SetBorderSize(2);
    palette->SetLineColor(1);
  }



  pad[2]->RedrawAxis(); 

  textLabel[2] = new TPaveText(xMin + 0.02*xRange, yMax-0.2*yRange, xMin+0.30*xRange, yMax-0.05*yRange);
  PlasmaGlob::SetPaveTextStyle(textLabel[2],12); 
  textLabel[2]->SetTextFont(42);
  textLabel[2]->AddText(sLabels[2]);
  textLabel[2]->Draw();

  C->cd();

  // Print to a file
  PlasmaGlob::imgconv(C,fOutName,opt);
  // ---------------------------------------------------------

  PlasmaGlob::DestroyCanvases();
}
  
