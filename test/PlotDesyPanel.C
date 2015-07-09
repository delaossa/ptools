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
#include <TPad.h>
#include <TPaletteAxis.h>

#ifndef __CINT__
#include "PFunctions.hh"
#endif

#include "PData.hh"
#include "PlasmaGlob.hh"
#include "PPalette.hh"

void PlotDesyPanel( const TString &sim, Int_t time, Double_t Emin, Double_t Emax, Int_t zoom=2, Int_t Nbins=2, const TString &options="") {

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
  // gStyle->SetPadLeftMargin(0.12);  // Margin left axis 
  // gStyle->SetPadRightMargin(0.12); // Margin for palettes in 2D histos
  // gStyle->SetTitleSize(0.06, "x");
  // gStyle->SetTitleSize(0.06, "y");
  // gStyle->SetTitleOffset(1.2, "x");
  // gStyle->SetTitleOffset(0.6,"y");
  gStyle->SetPadGridY(0);
  gStyle->SetPadGridX(0);
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
  Double_t skindepth = 1/kp;
  Double_t E0 = pData->GetPlasmaE0();
  
  // Some beam properties:
  Float_t nb = pData->GetBeamDensity();
  Double_t Ebeam = pData->GetBeamEnergy();
  Double_t gamma = pData->GetBeamGamma();
  Double_t vbeam = pData->GetBeamVelocity();
  Double_t kbeta = PFunc::BeamBetatronWavenumber(gamma,n0);
  Double_t rms0  = pData->GetBeamRmsY() * kp;
  if(CYL)  rms0  = pData->GetBeamRmsR() * kp;
  
  cout << Form(" - Bunch gamma      = %8.4f", gamma ) << endl;
  cout << Form(" - Bunch velocity   = %8.4f c", vbeam ) << endl;
  cout << Form(" - Bunch betatron k = %8.4f mm-1", kbeta * PUnits::mm) << endl;
  cout << endl;
  
  // Time in OU
  Float_t Time = pData->GetRealTime();
  // z start of the plasma in normalized units.
  Float_t zStartPlasma = pData->GetPlasmaStart() * kp;
  // z start of the beam in normalized units.
  Float_t zStartBeam = pData->GetBeamStart() * kp;
  
  if(opt.Contains("center")) {
    Time -= zStartPlasma;
    if(opt.Contains("comov"))      // Centers on the head of the beam.
      Time += zStartBeam;
  }
  Float_t shiftz = pData->Shift(opt);
  //  cout << "Shift = " << shiftz << endl;
  
  // Binning for 2D histograms:
  Int_t   x1Nbin = pData->GetNX(0);
  Float_t x1Min  = pData->GetXMin(0)-shiftz; 
  Float_t x1Max  = pData->GetXMax(0)-shiftz; 

  Int_t   x2Nbin = pData->GetNX(1);
  Float_t x2Min  = pData->GetXMin(1); 
  Float_t x2Max  = pData->GetXMax(1); 

  // Calculate the "axis range" in number of bins. If Nbins==0 a RMS width is taken.
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
  
  // axis pos
  Double_t axisPos = 0.0;
  
  // Binning for Energy histograms:
  Int_t eNbin = 100;
  Double_t eMin = Emin;
  Double_t eMax = Emax;
  
  // Get charge density histos
  Int_t Nspecies = pData->NSpecies();
  TH2F **hDen2D = new TH2F*[Nspecies];
  TH1F **hDen1D = new TH1F*[Nspecies];
  TH1F **hDen1Db = new TH1F*[Nspecies];
  TH1F **hEnergy = new TH1F*[Nspecies];    // Energy spectrum
  TH2F **hEvsZ = new TH2F*[Nspecies];      // Energy vs x1
  TH2F **hEvsY = new TH2F*[Nspecies];      // Energy vs x2
  // RMS (vs z) of the beam's charge distribution: 
  TProfile **hDen2Dprof = new TProfile*[Nspecies];
  TH1F **hRms = new TH1F*[Nspecies];
  TH1F **hRms2 = new TH1F*[Nspecies];
  TH1F **hRmsNorm = new TH1F*[Nspecies];     // normalized in respect with the initital value.
  // others
  Double_t *Q = new Double_t[Nspecies];
  Double_t *Qraw = new Double_t[Nspecies];
  TGraph **gEnergyX1 = new TGraph*[Nspecies];
  TGraph **gEnergyX2 = new TGraph*[Nspecies];
  
  for(Int_t i=0;i<Nspecies;i++) {
    hDen2D[i] = NULL;
    hDen1D[i] = NULL;
    hDen1Db[i] = NULL;
    hEnergy[i] = NULL;
    hEvsZ[i] = hEvsY[i] = NULL;
    hDen2Dprof[i] = NULL;
    hRms[i] = hRmsNorm[i] = hRms2[i] = NULL;
    Q[i] = Qraw[i] = 0.0;
    
    if(!pData->GetChargeFileName(i)) 
      continue;
    
    char hName[24];
    sprintf(hName,"hDen2D_%i",i);
    hDen2D[i] = (TH2F*) gROOT->FindObject(hName);
    if(hDen2D[i]) delete hDen2D[i];

    if(!ThreeD)
      hDen2D[i] = pData->GetCharge(i,opt);
    else
      hDen2D[i] = pData->GetCharge2DSliceZY(i,-1,1,opt+"avg");
     
    hDen2D[i]->SetName(hName);
    hDen2D[i]->GetXaxis()->CenterTitle();
    hDen2D[i]->GetYaxis()->CenterTitle();
    hDen2D[i]->GetZaxis()->CenterTitle();
    if(opt.Contains("comov"))
      hDen2D[i]->GetXaxis()->SetTitle("#zeta [c/#omega_{p}]");
    else
      hDen2D[i]->GetXaxis()->SetTitle("z [c/#omega_{p}]");
    hDen2D[i]->GetYaxis()->SetTitle("y [c/#omega_{p}]");
    if(CYL) hDen2D[i]->GetYaxis()->SetTitle("r [c/#omega_{p}]");
    
    if(i==0)
      hDen2D[i]->GetZaxis()->SetTitle("n [n_{0}]");
    else if(i==1)
      hDen2D[i]->GetZaxis()->SetTitle("n_{b} [n_{0}]");
    else
      hDen2D[i]->GetZaxis()->SetTitle("n_{w} [n_{0}]");
      
    // 1D histograms
    
    // b. mode:
    sprintf(hName,"hDen1Db_%i",i);
    hDen1Db[i] = (TH1F*) gROOT->FindObject(hName);
    if(hDen1Db[i]) delete hDen1Db[i];
    hDen1Db[i] = (TH1F*) hDen2D[i]->ProjectionX(hName,FirstyBin,LastyBin);
    hDen1Db[i]->Scale(1.0/(LastyBin-FirstyBin+1));
    
    // a. mode:   
    TString opth1 = opt;
    opth1 += "avg";
       
    sprintf(hName,"hDen1D_%i",i);
    hDen1D[i] = (TH1F*) gROOT->FindObject(hName);
    if(hDen1D[i]) delete hDen1D[i];

    if(ThreeD) {
      hDen1D[i] = pData->GetH1SliceZ3D(pData->GetChargeFileName(i)->c_str(),"charge",-1,Nbins,-1,Nbins,opth1.Data());
    } else if(CYL) { // Cylindrical: The first bin with r>0 is actually the number 1 (not the 0).
      hDen1D[i] = pData->GetH1SliceZ(pData->GetChargeFileName(i)->c_str(),"charge",1,Nbins,opth1.Data());
    } else { // 2D cartesian
      hDen1D[i] = pData->GetH1SliceZ(pData->GetChargeFileName(i)->c_str(),"charge",-1,Nbins,opth1.Data());
    }
    hDen1D[i]->SetName(hName); 
   
    if(opt.Contains("comov"))
      hDen1D[i]->GetXaxis()->SetTitle("#zeta [c/#omega_{p}]");
    else
      hDen1D[i]->GetXaxis()->SetTitle("z [c/#omega_{p}]");
    if(i==0)
      hDen1D[i]->GetYaxis()->SetTitle("n [n_{0}]");
    else
      hDen1D[i]->GetYaxis()->SetTitle("n_{b} [n_{0}]");
        
    // Histograms from RAW data.
    // --------------------------

    TString sName = pData->GetSpeciesName(i).c_str();
    if(!pData->GetRawFileName(i)) // || sName.Contains("plasma") )
      continue;
    
    sprintf(hName,"hEnergy_%i",i);
    hEnergy[i] = (TH1F*) gROOT->FindObject(hName);
    if(hEnergy[i]) { delete hEnergy[i]; hEnergy[i] = NULL; }
    hEnergy[i] = new TH1F(hName,"",eNbin,eMin,eMax);
    pData->GetH1Raw(pData->GetRawFileName(i)->c_str(),"energy",hEnergy[i],opt);
    
    if(hEnergy[i]->GetEntries()==0) {
      continue;
    }
    
    hEnergy[i]->GetXaxis()->CenterTitle();
    hEnergy[i]->GetYaxis()->CenterTitle();
    hEnergy[i]->GetXaxis()->SetTitle("Energy [GeV]");

    if(i==0)
      hEnergy[i]->GetYaxis()->SetTitle("n [n_{0}]");
    else
      hEnergy[i]->GetYaxis()->SetTitle("n_{b} [n_{0}]");
    PlasmaGlob::SetH1Style(hEnergy[i],i);
    
    sprintf(hName,"hEvsZ_%i",i);
    hEvsZ[i] = (TH2F*) gROOT->FindObject(hName);
    if(hEvsZ[i]) { delete hEvsZ[i]; hEvsZ[i] = NULL; }
    hEvsZ[i] = new TH2F(hName,"",x1Nbin,x1Min,x1Max,eNbin,eMin,eMax);
    pData->GetH2Raw(pData->GetRawFileName(i)->c_str(),"x1","energy",hEvsZ[i],opt);
    
    hEvsZ[i]->GetXaxis()->CenterTitle();
    hEvsZ[i]->GetYaxis()->CenterTitle();
    hEvsZ[i]->GetZaxis()->CenterTitle();
    hEvsZ[i]->GetYaxis()->SetTitle("Energy [GeV]");
    if(opt.Contains("comov")) {
      hEvsZ[i]->GetXaxis()->SetTitle("#zeta [c/#omega_{p}]");
      hEvsZ[i]->GetZaxis()->SetTitle("dN/d#zetadE [a.u.]");
    }  else {
      hEvsZ[i]->GetXaxis()->SetTitle("z [c/#omega_{p}]");
      hEvsZ[i]->GetZaxis()->SetTitle("dN/dzdE [a.u.]");
    }
    
    sprintf(hName,"hEvsY_%i",i);
    hEvsY[i] = (TH2F*) gROOT->FindObject(hName);
    if(hEvsY[i]) { delete hEvsY[i]; hEvsY[i] = NULL; }
    hEvsY[i] = new TH2F(hName,"",x2Nbin,x2Min,x2Max,eNbin,eMin,eMax);
    // Variable "ene" in RAW data stands for kinetic energy, given in gamma units.
    // To get the total energy do: E = m * ( "ene" + 1);
    // This is automatically done by PData, when the string "energy" is given (instead of "ene").
    pData->GetH2Raw(pData->GetRawFileName(i)->c_str(),"x2","energy",hEvsY[i],opt);
    
    hEvsY[i]->GetXaxis()->CenterTitle();
    hEvsY[i]->GetYaxis()->CenterTitle();
    hEvsY[i]->GetXaxis()->SetTitle("y [c/#omega_{p}]");
    hEvsY[i]->GetYaxis()->SetTitle("Energy [GeV]");
    PlasmaGlob::SetH1Style(hEvsY[i],i);
   
    if(hDen2D[i]) {
      TString pname = hDen2D[i]->GetName();
      pname += "_pfx";
      
      hDen2Dprof[i] =  (TProfile*) gROOT->FindObject(pname.Data());
      if(hDen2Dprof[i]) { delete hDen2Dprof[i]; hDen2Dprof[i] = NULL; }
      hDen2Dprof[i] = hDen2D[i]->ProfileX("_pfx",1,-1,"s");
      
      sprintf(hName,"hRms_%i",i);
      hRms[i] = (TH1F*) gROOT->FindObject(hName);
      if(hRms[i]) delete hRms[i];
      
      hRms[i] = new TH1F(hName,"",x1Nbin,x1Min,x1Max);
      
      sprintf(hName,"hRmsNorm_%i",i);
      hRmsNorm[i] = (TH1F*) gROOT->FindObject(hName);
      if(hRmsNorm[i]) delete hRmsNorm[i];
      
      hRmsNorm[i] = new TH1F(hName,"",x1Nbin,x1Min,x1Max);
      
      if(CYL) axisPos = 0.0;
      
      for(Int_t j=0;j<hRms[i]->GetNbinsX();j++) {
	Double_t rms = 0;
	Double_t total = 0;
	for(Int_t k=1;k<=x2Nbin;k++) {
	  Double_t value  = hDen2D[i]->GetBinContent(j,k);
	  Double_t radius = hDen2D[i]->GetYaxis()->GetBinCenter(k) - axisPos;
	  if(CYL) {
	    rms += radius*radius*radius*value;
	    total += radius*value;
	  } else {
	    rms += radius*radius*value;
	    total += value;
	  }
	  // cout << Form(" (%i,%i) -> radius = %7.4f ,  density = %7.4f",j,k,radius,value) << endl;
	}
	
	rms /= total;
	rms = sqrt(rms);
	
	hRmsNorm[i]->SetBinContent(j,rms/rms0); 
	hRms[i]->SetBinContent(j,rms); 
	
      }
      
      hRms[i]->GetXaxis()->SetTitle("z [c/#omega_{p}]");
      if(opt.Contains("comov"))
	hRms[i]->GetXaxis()->SetTitle("#zeta [c/#omega_{p}]");
      
      hRms[i]->GetYaxis()->SetTitle("r_{b} [c/#omega_{p}] ");
      
      hRmsNorm[i]->GetXaxis()->SetTitle("z [c/#omega_{p}]");
      if(opt.Contains("comov"))
	hRmsNorm[i]->GetXaxis()->SetTitle("#zeta [c/#omega_{p}]");
      
      hRmsNorm[i]->GetYaxis()->SetTitle("r_{b}/r_{0}");
      
    }
    
  
    // INTEGRATED Beam's Charge:
    for(Int_t j=1;j<=x1Nbin;j++) {
      for(Int_t k=1;k<=x2Nbin;k++) {
	Double_t value  = hDen2D[i]->GetBinContent(j,k);
	if(CYL) {
	  Double_t radius = hDen2D[i]->GetYaxis()->GetBinCenter(k);
	  Q[i] += radius * value;
	  // cout << Form(" (%i,%i) -> radius = %7.4f , value = %7.4f",j,k,radius,value) << endl;
	} else {
	  Q[i] += value;
	}
      }    
    }
    Double_t xbinsize = hDen2D[i]->GetXaxis()->GetBinWidth(1);
    Double_t ybinsize = hDen2D[i]->GetYaxis()->GetBinWidth(1); 
    Q[i] *= xbinsize * ybinsize;
    
    if(!CYL && !ThreeD) {
      Q[i] *= TMath::Sqrt(2*TMath::Pi()) * rms0; 
    } else if(CYL) {
      Q[i] *= 2*TMath::Pi();
    }
    
    if(opt.Contains("units")) {
      Double_t dV = skindepth * skindepth * skindepth;
      Q[i] *= n0 * dV;
      Q[i] *= (PConst::ElectronCharge/PUnits::picocoulomb); 
      cout << Form(" Integrated charge of specie %3i = %8i pC",i, TMath::Nint(Q[i])) << endl;
    } else {
      cout << Form(" Integrated charge of specie %3i = %8.4f n0 * kp^-3",i,Q[i]) << endl;
    }
    // ------------------------------------------------------------------------------------
    
    // Get integrated charge from RAW data:
    for(Int_t j=0; j<eNbin; j++) {
      Qraw[i] += hEnergy[i]->GetBinContent(j+1);
    }
    
    Qraw[i] *= xbinsize * ybinsize/2;
    if(!CYL && !ThreeD) {
      Qraw[i] *= TMath::Sqrt(2*TMath::Pi()) * rms0; 
    } else if(CYL) {
      Qraw[i] *= 2*TMath::Pi();
    } else {
      Qraw[i] *= ybinsize;
    }
    
    if(opt.Contains("units")) {    
      Double_t dV = skindepth * skindepth * skindepth;
      Qraw[i] *= n0 * dV;
      Qraw[i] *= (PConst::ElectronCharge/PUnits::picocoulomb);
    }
    
    cout << Form(" Integrated charge (RAW) of specie %3i = %8.4f a.u. ",i, Qraw[i]) << endl;
    // ------------------------------------------------------------------------------------
    
  }
  

  // cjewklcdms
  
  // Get electric fields
  const Int_t Nfields = 2;
  TH1F **hE1D = new TH1F*[Nfields];
  for(Int_t i=0;i<Nfields;i++) {
    hE1D[i] = NULL;

    if(!pData->GetEfieldFileName(i))
      continue;
      
    // 1D histograms
    TString opth1 = opt;
    opth1 += "avg";
    
    char nam[3]; sprintf(nam,"e%i",i+1);
    if(ThreeD) {
    
      if(i==0) 
	hE1D[i] = pData->GetH1SliceZ3D(pData->GetEfieldFileName(i)->c_str(),nam,-1,Nbins,-1,Nbins,opth1.Data());
      else  
	hE1D[i] = pData->GetH1SliceZ3D(pData->GetEfieldFileName(i)->c_str(),nam,-Nbins,Nbins,-Nbins,Nbins,opth1.Data());
   
    } else if(CYL) { // Cylindrical: The firt bin with r>0 is actually the number 1 (not the 0).
    
	hE1D[i] = pData->GetH1SliceZ(pData->GetEfieldFileName(i)->c_str(),nam,1,Nbins,opth1.Data());
    
    } else { // 2D cartesian
    
      if(i==0) 
	hE1D[i] = pData->GetH1SliceZ(pData->GetEfieldFileName(i)->c_str(),nam,-1,Nbins,opth1.Data());
      else 
	hE1D[i] = pData->GetH1SliceZ(pData->GetEfieldFileName(i)->c_str(),nam,-Nbins,Nbins,opth1.Data());
    
    }
    
    char hName[24];
    sprintf(hName,"hE1D_%i",i);
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
    
    axisPos *= skindepth / PUnits::um;
    x1Min *= skindepth / PUnits::um;
    x1Max *= skindepth / PUnits::um;
    x2Min *= skindepth / PUnits::um;
    x2Max *= skindepth / PUnits::um;
    
    for(Int_t i=0;i<Nspecies;i++) {
      hDen2D[i]->SetBins(x1Nbin,x1Min,x1Max,x2Nbin,x2Min,x2Max);
      // for(Int_t j=0;j<=x1Nbin;j++) {
      // 	for(Int_t k=0;k<=x2Nbin;k++) {
      // 	  hDen2D[i]->SetBinContent(j,k, hDen2D[i]->GetBinContent(j,k) * n0 / (1e17/PUnits::cm3) );
      // 	}
      // }
      
      hDen2D[i]->GetYaxis()->SetTitle("y [#mum]");  
      if(CYL) hDen2D[i]->GetYaxis()->SetTitle("r [#mum]");
      
      if(opt.Contains("comov"))
	hDen2D[i]->GetXaxis()->SetTitle("#zeta [#mum]");
      else
	hDen2D[i]->GetXaxis()->SetTitle("z [#mum]");
      
      // if(i==0)
      // 	hDen2D[i]->GetZaxis()->SetTitle("n [10^{17}/cm^{3}]");
      // else
      // 	hDen2D[i]->GetZaxis()->SetTitle("n_{b} [10^{17}/cm^{3}]"); 
      
      hDen1D[i]->SetBins(x1Nbin,x1Min,x1Max);
      // for(Int_t j=0;j<x1Nbin;j++) {
      // 	Double_t bincontent = (hDen1D[i]->GetBinContent(j) * n0 / (1e17/PUnits::cm3));
      // 	hDen1D[i]->SetBinContent(j,bincontent);
      // }
      
      // hDen1D[i]->GetYaxis()->SetTitle("n [10^{17}/cm^{3}]");
      
      if(opt.Contains("comov"))
	hDen1D[i]->GetXaxis()->SetTitle("#zeta [#mum]");
      else
	hDen1D[i]->GetXaxis()->SetTitle("z [#mum]");


      hDen1Db[i]->SetBins(x1Nbin,x1Min,x1Max);
      // for(Int_t j=0;j<x1Nbin;j++) {
      // 	Double_t bincontent = (hDen1Db[i]->GetBinContent(j) * n0 / (1e17/PUnits::cm3));
      // 	hDen1Db[i]->SetBinContent(j,bincontent);
      // }
      
      // hDen1Db[i]->GetYaxis()->SetTitle("n [10^{17}/cm^{3}]");
      
      if(opt.Contains("comov"))
	hDen1Db[i]->GetXaxis()->SetTitle("#zeta [#mum]");
      else
	hDen1Db[i]->GetXaxis()->SetTitle("z [#mum]");
      
      if(hEvsZ[i]) {
	hEvsZ[i]->SetBins(x1Nbin,x1Min,x1Max,eNbin,eMin/1000,eMax/1000);
	for(Int_t j=0;j<x1Nbin;j++) {
	  for(Int_t k=0;k<eNbin;k++) {
	    hEvsZ[i]->SetBinContent(j,k,hEvsZ[i]->GetBinContent(j,k) * n0 / (1e17/PUnits::cm3));
	  }
	}
	
	if(opt.Contains("comov")) {
	  hEvsZ[i]->GetXaxis()->SetTitle("#zeta [#mum]");
	  hEvsZ[i]->GetZaxis()->SetTitle("dN/d#zetadE [a.u.]");
	}  else {
	  hEvsZ[i]->GetXaxis()->SetTitle("z [#mum]");
	  hEvsZ[i]->GetZaxis()->SetTitle("dN/dzdE [a.u.]");
	}
      }
	
      if(hEvsY[i]) {
	hEvsY[i]->SetBins(x2Nbin,x2Min,x2Max,eNbin,eMin,eMax);	
	for(Int_t j=0;j<x2Nbin;j++) {
	  for(Int_t k=0;k<eNbin;k++) {
	    hEvsY[i]->SetBinContent(j,k,hEvsY[i]->GetBinContent(j,k) * n0 / (1e17/PUnits::cm3)) ;
	  }
	}
	
	hEvsY[i]->GetXaxis()->SetTitle("y [#mum]");  
	if(CYL) hEvsY[i]->GetXaxis()->SetTitle("r [#mum]");  
      }
        
      if(hRms[i]) {
	hRms[i]->SetBins(x1Nbin,x1Min,x1Max);
	for(Int_t j=0;j<x1Nbin;j++) {
	  Double_t bincontent = (hRms[i]->GetBinContent(j) * skindepth / PUnits::um);
	  hRms[i]->SetBinContent(j,bincontent);
	}
	hRms[i]->GetYaxis()->SetTitle(" r_{b} _{rms} [#mum]");         
	if(opt.Contains("comov"))
	  hRms[i]->GetXaxis()->SetTitle("#zeta [#mum]");
	else
	  hRms[i]->GetXaxis()->SetTitle("z [#mum]");
      }
      
      if(hRmsNorm[i]) {
	hRmsNorm[i]->SetBins(x1Nbin,x1Min,x1Max);    
	if(opt.Contains("comov"))
	  hRmsNorm[i]->GetXaxis()->SetTitle("#zeta [#mum]");
	else
	  hRmsNorm[i]->GetXaxis()->SetTitle("z [#mum]");
      }
    }
    
    for(Int_t i=0;i<Nfields;i++) {
      hE1D[i]->SetBins(x1Nbin,x1Min,x1Max);
      
      for(Int_t j=0;j<x1Nbin;j++) {
	hE1D[i]->SetBinContent(j, hE1D[i]->GetBinContent(j) * ( E0 / (PUnits::GV/PUnits::m) ) );
      }
      
      if(opt.Contains("comov"))
	hE1D[i]->GetXaxis()->SetTitle("#zeta [#mum]");
      else
	hE1D[i]->GetXaxis()->SetTitle("z [#mum]");
      
      if(i==0)
	hE1D[i]->GetYaxis()->SetTitle("E_{z} [GV/m]");
      else if(i==1)
	hE1D[i]->GetYaxis()->SetTitle("E_{y} [GV/m]");
      else if(i==2)
	hE1D[i]->GetYaxis()->SetTitle("E_{x} [GV/m]");
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
    
  for(Int_t i=0;i<Nspecies;i++) {
     
    // Trick for nicer plotting
    // It uses hDen1D[i] to setup the valid range
    Double_t bin1 = -1;
    Double_t bin2 = -1;
    if(hDen1D[i] && hRms[i]) {
      
      for(Int_t j=0;j<hRms[i]->GetNbinsX();j++) 
	if(bin1==-1 && hDen1D[i]->GetBinContent(j+1)>0.0005) {
	  bin1 = hRms[i]->GetBinCenter(j+1); break;
	} 
      
      for(Int_t j=hRms[i]->GetNbinsX()-1;j>=0;j--) 
	if(bin2==-1 && hDen1D[i]->GetBinContent(j+1)>0.0005) {
	  bin2 = hRms[i]->GetBinCenter(j+1); break;
	}
      
      hRms[i]->GetXaxis()->SetRangeUser(bin1,bin2);
      hRmsNorm[i]->GetXaxis()->SetRangeUser(bin1,bin2);
      
      // Clone and reverse hRms
      char hName[24];
      sprintf(hName,"hRms2_%i",i);
      hRms2[i] = (TH1F*) hRms[i]->Clone(hName);
      for(Int_t j=0;j<hRms[i]->GetNbinsX();j++) {
	Double_t rms = hRms[i]->GetBinContent(j);
	hRms2[i]->SetBinContent(j,axisPos - rms);
	hRms[i]->SetBinContent(j,axisPos + rms);
      }
    }
    
    // Vertical graphs: Displayed at a side of 2D histograms
   
    if(hEnergy[i]) {
      Double_t *y   = new Double_t[eNbin];
      Double_t *x   = new Double_t[eNbin];
      
      // This is for the right side:
      // Double_t xMax = x1Min + (x1Max-x1Min) * 0.9;
      // Double_t xMin = x1Max;
      // And this for left:
      Double_t xMin = x1Min;
      Double_t xMax = x1Min + (x1Max-x1Min) * 0.2;
      Double_t EneMax = hEnergy[i]->GetMaximum();
      for(Int_t j=0; j<eNbin; j++) {
	y[j] = hEnergy[i]->GetBinCenter(j+1);
	x[j] = ((xMax-xMin)/EneMax)*hEnergy[i]->GetBinContent(j+1) + xMin;
      }
      gEnergyX1[i] = new TGraph(eNbin,x,y);
      gEnergyX1[i]->SetLineColor(PlasmaGlob::elecLine);
      gEnergyX1[i]->SetLineWidth(2);
      gEnergyX1[i]->SetFillStyle(1001);
      gEnergyX1[i]->SetFillColor(PlasmaGlob::elecFill);
      
      // Right side:
      // xMax = x2Min + (x2Max-x2Min) * 0.9;
      // xMin = x2Max;
      // Left side:
      xMin = x2Min;
      xMax = x2Min +  (x2Max-x2Min) * 0.2;    
      for(Int_t j=0; j<eNbin; j++) {
	y[j] = hEnergy[i]->GetBinCenter(j+1);
	x[j] = ((xMax-xMin)/EneMax)*hEnergy[i]->GetBinContent(j+1) + xMin;
      }
      gEnergyX2[i] = new TGraph(eNbin,x,y);
      gEnergyX2[i]->SetLineColor(PlasmaGlob::elecLine);
      gEnergyX2[i]->SetLineWidth(2);
      gEnergyX2[i]->SetFillStyle(1001);
      gEnergyX2[i]->SetFillColor(PlasmaGlob::elecFill);
      
    }
  }

  // Plotting
  // -----------------------------------------------

  // Canvas setup
  TCanvas *C;
  if(opt.Contains("hres") && !opt.Contains("pdf")) // high resolution for plain grahics output.
    C = new TCanvas("C","2D Charge density, Electric field and Energy",1500,2000);
  else
    C = new TCanvas("C","2D Charge density, Electric field and Energy",750,1000);
  C->Divide(1,3);

  // Palettes setup
  TExec *exPlasma = new TExec("exPlasma","plasmaPalette->cd();"); 
  TExec *exHot    = new TExec("exHot","hotPalette->cd();");
  TExec *exElec   = new TExec("exElec","electronPalette->cd();");
  TExec *exField  = new TExec("exField","rbowgrayPalette->cd();");

  // Text objects
  TPaveText *textTime = new TPaveText(0.7,0.85,0.85,0.92,"NDC");
  PlasmaGlob::SetPaveTextStyle(textTime,32); 
  textTime->SetTextColor(kGray+3);
  char ctext[128];
  if(opt.Contains("units") && n0) 
    sprintf(ctext,"z = %5.1f mm", Time * skindepth / PUnits::mm);
  else
    sprintf(ctext,"t = %5.1f #omega_{p}^{-1}",Time);
  textTime->AddText(ctext);

  TPaveText *textDen = new TPaveText(0.13,0.85,0.38,0.92,"NDC");
  PlasmaGlob::SetPaveTextStyle(textDen,12); 
  textDen->SetTextColor(kOrange+10);
  if(opt.Contains("units") && n0)
    sprintf(ctext,"n_{0} = %5.2f x 10^{17} / cm^{3}",1e-17 * n0 * PUnits::cm3);
  else if(pData->GetBeamDensity() && n0)
    sprintf(ctext,"n_{b}/n_{0} = %5.2f", pData->GetBeamDensity()/n0);
  textDen->AddText(ctext);

  TPaveText *textCharge = new TPaveText(0.6,0.05,0.85,0.20,"NDC");
  PlasmaGlob::SetPaveTextStyle(textCharge,32); 
  textCharge->SetTextColor(kGray+1);
  for(Int_t i=1;i<Nspecies;i++) {
    if(opt.Contains("units")) 
      sprintf(ctext,"Q_{%1i} = %5i pC",i,TMath::Nint(Q[i]));
    else
      sprintf(ctext,"Q_{%1i} = %8.4f (c/#omega_{p})^{3} #times n_{0}",i,Q[i]);
    
    textCharge->AddText(ctext);
  }

  TPaveText *textWav = new TPaveText(0.13,0.75,0.38,0.82,"NDC");
  PlasmaGlob::SetPaveTextStyle(textWav,12); 
  textWav->SetTextColor(kGray+2);
  sprintf(ctext,"#lambda_{p} = %5.3f um", pData->GetPlasmaWaveLength()/PUnits::um);
  textWav->AddText(ctext);

  TPaveText *textLabel[3];
  for(Int_t i=0;i<3;i++) {
 
    if(i==0) {
      textLabel[i] = new TPaveText(0.14,0.08,0.19,0.17,"NDC");
      textLabel[i]->AddText("(a)");
    } else if(i==1) {
      textLabel[i] = new TPaveText(0.14,0.08,0.19,0.17,"NDC");
      textLabel[i]->AddText("(b)");
    } else if(i==2) {
      textLabel[i] = new TPaveText(0.14,0.23,0.19,0.30,"NDC");
      textLabel[i]->AddText("(c)");
    }

    
    PlasmaGlob::SetPaveTextStyle(textLabel[i]); 
    textLabel[i]->SetTextColor(kGray+3);
  }
  
  
  // Actual Plotting!
  // ------------------------------------------------------------

  // Output file
  TString fOutName = Form("./%s/Plots/DesyPanel/DesyPanel",pData->GetPath().c_str());
  fOutName += Form("-%s_%i",pData->GetName(),time);

  // Setup Pad layout:
  Double_t lMargin = 0.10;
  Double_t rMargin = 0.12;
  Double_t bMargin = 0.07;
  Double_t tMargin = 0.02;
  Double_t vSpacing = 0.01; 
  Double_t hStep = (1.-lMargin-rMargin);
  Double_t vStep = (1.-bMargin-tMargin)/3.;
 
  TPad *pad[3];

  // top plots
  pad[0] = new TPad("padt", "padt",0.00, bMargin + 2.*vStep + vSpacing,
		    lMargin+hStep+rMargin, 1.00);
  pad[0]->SetLeftMargin(1./(lMargin+hStep)*lMargin);
  pad[0]->SetRightMargin(1./(rMargin+hStep)*rMargin);  
  pad[0]->SetBottomMargin(0.0);                                   
  pad[0]->SetTopMargin(1./(tMargin+vStep)*tMargin);
  pad[0]->Draw();

  // middle plots
  pad[1] = new TPad("padm", "padm",0.00, bMargin + vStep + vSpacing,
		    lMargin + hStep + rMargin, bMargin + 2.*vStep );
  pad[1]->SetLeftMargin(1./(lMargin+hStep)*lMargin);
  pad[1]->SetRightMargin((1./(rMargin+hStep)*rMargin));
  pad[1]->SetBottomMargin(0.0);                                   
  pad[1]->SetTopMargin(0.);
  pad[1]->Draw();          

  // bottom plots
  pad[2] = new TPad("padb", "padb",0.00, 0.,lMargin+hStep+rMargin,bMargin+vStep);
  pad[2]->SetLeftMargin(1./(lMargin+hStep)*lMargin);
  pad[2]->SetRightMargin((1./(rMargin+hStep)*rMargin));
  pad[2]->SetBottomMargin(1./(bMargin+vStep)*bMargin);
  pad[2]->SetTopMargin(0.);
  pad[2]->Draw();       


  pad[0]->cd(); // <--------------------------------------------- Top Plot -----------
  if(opt.Contains("logz")) {
    pad[0]->SetLogz(1);
  } else {
    pad[0]->SetLogz(0);
  }
  pad[0]->SetFrameLineWidth(3);  

  TH2F *hFrame = (TH2F*) gROOT->FindObject("hFrame");
  if(hFrame) delete hFrame;
  hFrame = (TH2F*) hDen2D[0]->Clone("hFrame");
  hFrame->Reset();

  hFrame->GetXaxis()->SetLabelOffset(0.10);
 
  hFrame->GetYaxis()->SetTitleSize(0.075);
  hFrame->GetYaxis()->SetTitleOffset(0.65);
  hFrame->GetYaxis()->SetLabelSize(0.06);
  //hFrame->GetYaxis()->SetLabelOffset(0.05);

  hFrame->GetZaxis()->SetLabelSize(0.06);  
  hFrame->GetZaxis()->SetTitleSize(0.06);                        
  hFrame->GetZaxis()->SetTitleOffset(0.06);
  //  hFrame->GetZaxis()->SetNdivisions(505);

  hFrame->Draw("col");
  
  exPlasma->Draw();
   hDen2D[0]->Draw("colz same");
  
  if(hDen2D[1]) {
    exElec->Draw();
    hDen2D[1]->Draw("colz same");

    hRms[1]->SetLineStyle(2);
    hRms[1]->SetLineColor(kGray+2);
    //hRms[1]->Draw("same C");
    
    hRms2[1]->SetLineStyle(2);
    hRms2[1]->SetLineColor(kGray+2);
    //hRms2[1]->Draw("same C");
  }

  if(Nspecies>=3) {
    if(hDen2D[2]) {
      exHot->Draw();
      hDen2D[2]->Draw("col same");
    }
  }
  
  TLine *lineY0 = new TLine(hDen2D[1]->GetXaxis()->GetXmin(),0.,
			    hDen2D[1]->GetXaxis()->GetXmax(),0.);
  lineY0->SetLineColor(kGray+1);
  lineY0->SetLineStyle(2);
  //lineY0->Draw();

  // "Axis range" in Osiris units:
  Double_t ymin = hDen2D[1]->GetYaxis()->GetBinLowEdge(FirstyBin);
  Double_t ymax = hDen2D[1]->GetYaxis()->GetBinUpEdge(LastyBin);
  TLine *lineYup = new TLine(hDen2D[1]->GetXaxis()->GetXmin(),ymax,
			     hDen2D[1]->GetXaxis()->GetXmax(),ymax);
  lineYup->SetLineColor(kGray+1);
  lineYup->SetLineStyle(2);
  //lineYup->Draw();  
  TLine *lineYdown = new TLine(hDen2D[1]->GetXaxis()->GetXmin(),ymin,
			     hDen2D[1]->GetXaxis()->GetXmax(),ymin);
  lineYdown->SetLineColor(kGray+1);
  lineYdown->SetLineStyle(2);
  //lineYdown->Draw();  

  pad[0]->Update();
  TPaletteAxis *palette = (TPaletteAxis*)hDen2D[1]->GetListOfFunctions()->FindObject("palette");
  Float_t y1 = pad[0]->GetBottomMargin();
  Float_t y2 = 1 - pad[0]->GetTopMargin();
  Float_t x1 = 1 - pad[0]->GetRightMargin();
  palette->SetY2NDC(y2 - 0.01);
  palette->SetY1NDC(0.5*(y1+y2) + 0.01);
  palette->SetX1NDC(x1 + 0.005);
  palette->SetX2NDC(x1 + 0.03);
  palette->SetTitleOffset(0.65);
  palette->SetTitleSize(0.07);
  palette->SetBorderSize(2);
  palette->SetLineColor(1);
  
  TPaletteAxis *palette2 = (TPaletteAxis*)hDen2D[0]->GetListOfFunctions()->FindObject("palette");
  palette2->SetY2NDC(0.5*(y1+y2) - 0.01);
  palette2->SetY1NDC(y1 + 0.01);
  palette2->SetX1NDC(x1 + 0.005);
  palette2->SetX2NDC(x1 + 0.03);
  palette2->SetTitleOffset(0.65);
  palette2->SetTitleSize(0.07);
  palette2->SetBorderSize(2);
  palette2->SetLineColor(1);
 
  //text1->Draw();
  textTime->Draw();
  textDen->Draw();
  if(opt.Contains("units"))
    textWav->Draw();
  textLabel[0]->Draw();
  //textCharge->Draw();
  
  pad[0]->RedrawAxis(); 

  pad[1]->cd(); // <----------------------------------------- Mid Plot ---------------
  pad[1]->SetFrameLineWidth(2);   

  Int_t beamC   = kAzure-5;
  Int_t fieldC  = kOrange+10; // PlasmaGlob::fieldLine;
  Int_t fieldCb = kGray+1; // PlasmaGlob::fieldLine;
  
  hE1D[0]->SetLineWidth(3);
  hE1D[0]->GetYaxis()->CenterTitle();
  hE1D[0]->GetYaxis()->SetAxisColor(fieldC);
  hE1D[0]->GetYaxis()->SetLabelColor(fieldC);
  hE1D[0]->GetYaxis()->SetTitleColor(fieldC);
  hE1D[0]->GetXaxis()->CenterTitle();

  hE1D[0]->GetXaxis()->SetLabelOffset(0.10);
  
  Float_t factor = 1.15;
  hE1D[0]->GetYaxis()->SetTitleSize(0.075*factor);
  hE1D[0]->GetYaxis()->SetTitleOffset(0.65/factor);
  hE1D[0]->GetYaxis()->SetLabelSize(0.06*factor);
  // hE1D[0]->GetYaxis()->SetLabelOffset(0.05/factor);
 
  hE1D[0]->GetYaxis()->SetNdivisions(505);

  hE1D[0]->SetLineStyle(1);
  hE1D[0]->SetLineColor(fieldC);
  hE1D[0]->SetMarkerStyle(20);

  hE1D[1]->GetYaxis()->CenterTitle();
  hE1D[1]->GetXaxis()->CenterTitle();
  hE1D[1]->SetLineStyle(1);
  hE1D[1]->SetLineColor(fieldC);
  hE1D[1]->SetMarkerStyle(24);

  if(opt.Contains("smooth")) {
    hE1D[0]->Smooth(10);
    hE1D[1]->Smooth(10);
  }

  Double_t factor2 = 1.7;
  Double_t minimum = factor2 * hE1D[0]->GetMinimum();
  Double_t maximum = factor2 * hE1D[0]->GetMaximum();
    
  // if(hE1D[1]->GetMaximum() > hE1D[0]->GetMaximum()) {
  //   maximum = factor2 * hE1D[1]->GetMaximum();
  // } 

  // if(hE1D[1]->GetMinimum() < hE1D[0]->GetMinimum()) {
  //   minimum = factor2 * hE1D[1]->GetMinimum();
  // } 

  if( maximum >= TMath::Abs(minimum)) minimum = -maximum;
  else maximum = - minimum;
  
  hE1D[0]->GetYaxis()->SetRangeUser(minimum,maximum);
  
  if(opt.Contains("units")) 
    hE1D[0]->GetYaxis()->SetTitle("E [GV/m]");
  else
    hE1D[0]->GetYaxis()->SetTitle("E [E_{0}]");
  
  hE1D[0]->Draw("C");
  //hE1D[1]->Draw("C same");
  
  pad[1]->Update();
  
  TLine *line0 = new TLine(hE1D[0]->GetXaxis()->GetXmin(),
			   (pad[1]->GetUymin()+pad[1]->GetUymax())/2.,
			   hE1D[0]->GetXaxis()->GetXmax(),
			   (pad[1]->GetUymin()+pad[1]->GetUymax())/2.);
  line0->SetLineColor(kGray+1);
  line0->SetLineStyle(2);
  line0->Draw();
  
  TGaxis **axis = new TGaxis*[Nspecies];
  for(Int_t i=0;i<Nspecies;i++) {
    if(!hDen1D[i]) continue;
    
    Float_t rightmin = 0;
    Float_t rightmax = 2.5 * hDen1D[i]->GetMaximum();
    Float_t slope = (pad[1]->GetUymax() - pad[1]->GetUymin())/(rightmax-rightmin);
    
    for(Int_t j=0;j<hDen1D[i]->GetNbinsX();j++) {
      hDen1D[i]->SetBinContent(j+1,hDen1D[i]->GetBinContent(j+1)*slope + pad[1]->GetUymin());
      hDen1Db[i]->SetBinContent(j+1,hDen1Db[i]->GetBinContent(j+1)*slope + pad[1]->GetUymin());
    }
    
    axis[i] = new TGaxis(pad[1]->GetUxmax(),pad[1]->GetUymin(),pad[1]->GetUxmax(),
			 pad[1]->GetUymax(),0.00001,rightmax,505,"+L");
    
  }
  
  
  hDen1D[1]->SetLineWidth(3);
  hDen1D[1]->SetLineColor(PlasmaGlob::elecLine);
  // PlasmaGlob::SetH1Style(hDen1D[1],1);
  hDen1D[1]->Draw("same C");
  
  hDen1Db[1]->SetLineColor(kPink);
  //hDen1Db[1]->Draw("same C");

  // if(hDen1D[0]) {
  //   hDen1D[0]->SetLineWidth(2);
  //   hDen1D[0]->SetLineColor(kGray+2);
  //   hDen1D[0]->Draw("same C");
  // }
  
  if(Nspecies>=3) {
    if(hDen1D[2]) {
      hDen1D[2]->SetLineWidth(2);
      hDen1D[2]->SetLineColor(kGray+2);
      hDen1D[2]->Draw("same C");
    }
  }
  
  
  axis[1]->SetLineWidth(1);
  axis[1]->SetLineColor(PlasmaGlob::elecLine);
  axis[1]->SetLabelColor(PlasmaGlob::elecLine);
  axis[1]->SetLabelSize(0.06);
  axis[1]->SetTitleSize(0.07);
  axis[1]->SetTitleOffset(0.8);
  // if(opt.Contains("units") && n0) 
  //   axis[1]->SetTitle("n_{b} [10^{17}/cm^{3}]");
  // else
    axis[1]->SetTitle("n_{b} [n_{0}]");
  axis[1]->CenterTitle();
  axis[1]->SetTitleColor(PlasmaGlob::elecLine);
  
  axis[1]->Draw();
  
  // Bunch RMS
  Float_t rightmin = 0;
  Float_t rightmax = 2.5 * hRmsNorm[1]->GetMaximum();
  Float_t slope = (pad[1]->GetUymax() - pad[1]->GetUymin())/(rightmax-rightmin);
  
  for(Int_t j=0;j<hRmsNorm[1]->GetNbinsX();j++) {
    hRmsNorm[1]->SetBinContent(j+1,(hRmsNorm[1]->GetBinContent(j+1)-rightmin)*slope + pad[1]->GetUymin());
    
    if(Nspecies>=3) {
      if(hRmsNorm[2]) {
	hRmsNorm[2]->SetBinContent(j+1,(hRmsNorm[2]->GetBinContent(j+1)-rightmin)*slope + pad[1]->GetUymin());
      }
    }
  }
  
  hRmsNorm[1]->SetLineStyle(1);
  hRmsNorm[1]->SetLineColor(PlasmaGlob::elecLine);
  hRmsNorm[1]->Draw("same C");

  if(Nspecies>=3) {
    if(hRmsNorm[2]) {
      hRmsNorm[2]->SetLineStyle(1);
      hRmsNorm[2]->SetLineColor(kGray+2);
      // hRmsNorm[2]->Draw("same C");
    }
  }

  TLine *line1 = new TLine(hFrame->GetXaxis()->GetXmin(), (1.0-rightmin)*slope + pad[1]->GetUymin(),
			   hFrame->GetXaxis()->GetXmax(), (1.0-rightmin)*slope + pad[1]->GetUymin());
  line1->SetLineColor(PlasmaGlob::elecLine);
  line1->SetLineStyle(2);
  line1->Draw();

  //draw an axis on the right side
  Double_t rightmargin = 0.06;
  Double_t ux = pad[1]->PixeltoX(pad[1]->UtoPixel(1-rightmargin));
  TGaxis *axisRMS = new TGaxis(ux,pad[1]->GetUymin(),
			       ux,pad[1]->GetUymax(),rightmin,rightmax,505,"+L");
  
  axisRMS->SetLineWidth(1);
  axisRMS->SetLineColor(PlasmaGlob::elecLine);
  axisRMS->SetLabelColor(PlasmaGlob::elecLine);
  axisRMS->SetTitle("r_{b}/r_{0}");
  axisRMS->CenterTitle();
  axisRMS->SetTitleColor(PlasmaGlob::elecLine);
  axisRMS->SetTitleOffset(0.22);
  // axisRMS->Draw();

  textLabel[1]->Draw();

  // --------------------------------------------------------- Bottom plot! --------------
 
  pad[2]->cd();
  pad[2]->SetFrameLineWidth(2); 

  TH2F *hFrame3 = (TH2F*) gROOT->FindObject("hFrame3");
  if(hFrame3) delete hFrame3;
  hFrame3 = (TH2F*) hEvsZ[1]->Clone("hFrame3");
  hFrame3->Reset();

  factor = 0.90;
  hFrame3->GetXaxis()->SetTitleSize(0.07*factor);
  hFrame3->GetXaxis()->SetLabelSize(0.07*factor);

  hFrame3->GetYaxis()->SetTitleSize(0.070*factor);
  hFrame3->GetYaxis()->SetLabelSize(0.06*factor);
  hFrame3->GetYaxis()->SetTitleOffset(0.75/factor);

  hEvsZ[1]->GetZaxis()->SetTitleSize(0.05*factor);                        
  hEvsZ[1]->GetZaxis()->SetTitleOffset(0.7/factor);
  hEvsZ[1]->GetZaxis()->SetLabelSize(0.05*factor);  
  hEvsZ[1]->GetZaxis()->SetTickLength(0.02/factor);
  // hEvsZ[1]->GetZaxis()->SetNdivisions(505);

  hFrame3->Draw();

  exElec->Draw();
  hEvsZ[1]->Draw("col same z");
    
  //gEnergyX1[1]->Draw("F");
  gEnergyX1[1]->Draw("L");

  if(hEvsZ[0]) {
    exElec->Draw();
    hEvsZ[0]->Draw("col same");

    if(gEnergyX1[0]) {
      
      gEnergyX1[0]->SetLineColor(kGray+2);
      
      //gEnergyX1[0]->Draw("F");
      gEnergyX1[0]->Draw("L");
      
    }
  }
 
  if(Nspecies>=3) {
    if(hEvsZ[2]) {
      exElec->Draw();
      hEvsZ[2]->Draw("col same");
    }

    if(gEnergyX1[2]) {
      
      gEnergyX1[2]->SetLineColor(kGray+2);
      
      //gEnergyX1[2]->Draw("F");
      gEnergyX1[2]->Draw("L");
      
    }
    
  }
  
  TLine *lineEne0 = new TLine(hFrame3->GetXaxis()->GetXmin(),Ebeam/PUnits::GeV,
			    hFrame3->GetXaxis()->GetXmax(),Ebeam/PUnits::GeV);
  lineEne0->SetLineColor(kGray+1);
  lineEne0->SetLineStyle(2);
  lineEne0->Draw();
  
  // text4->Draw();
  // textTime->Draw();
  textLabel[2]->Draw();

  pad[2]->Update();
  TPaletteAxis *palette3 = (TPaletteAxis*)hEvsZ[1]->GetListOfFunctions()->FindObject("palette");
  y1 = pad[2]->GetBottomMargin();
  y2 = 1 - pad[2]->GetTopMargin();
  x1 = 1 - pad[2]->GetRightMargin();
  palette3->SetY2NDC(y2 - 0.01);
  palette3->SetY1NDC(y1 + 0.01);
  palette3->SetX1NDC(x1 + 0.005);
  palette3->SetX2NDC(x1 + 0.03);
  palette3->SetTitleOffset(0.65);
  palette3->SetTitleSize(0.07);
  palette3->SetBorderSize(2);
  palette3->SetLineColor(1);

  pad[2]->RedrawAxis(); 
  pad[2]->RedrawAxis("G"); 

  C->cd();

  // Print to a file
  PlasmaGlob::imgconv(C,fOutName,opt);
  // ---------------------------------------------------------

}
