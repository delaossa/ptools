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
#include "PGlobals.hh"
#include "PPalette.hh"

void PlotRake5Panel2D( const TString &sim, Int_t time, Float_t zoom=1, Int_t Nbins=1, const TString &options="") {
  
#ifdef __CINT__  
  gSystem->Load("libptools.so");
#endif

  PGlobals::Initialize();

  // Palettes!
  gROOT->Macro("PPalettes.C");

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
  Float_t n0 = pData->GetPlasmaDensity();
  Float_t omegap = pData->GetPlasmaFrequency();
  Float_t timedepth = pData->GetPlasmaTimeDepth();
  Float_t kp = pData->GetPlasmaK();
  Float_t skindepth = pData->GetPlasmaSkinDepth();
  Float_t E0 = pData->GetPlasmaE0();
  
  // Some beam properties:
  Float_t Ebeam = pData->GetBeamEnergy();
  Float_t gamma = pData->GetBeamGamma();
  Float_t vbeam = pData->GetBeamVelocity();

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
  // cout << "Shift = " << shiftz << endl;
  

  // Calculate the "axis range" in number of bins. If Nbins==0 a RMS width is taken.
  Float_t rms0 = pData->GetBeamRmsY() * kp;
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
  const Int_t Nfields = 3;
  TH2F **hE2D = new TH2F*[Nfields];
  TH1F **hE1D = new TH1F*[Nfields];
  TH2F *hV2D = NULL;
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
     
      sprintf(hName,"hV2D");
      hV2D = (TH2F*) hE2D[i]->Clone(hName);
      hV2D->Reset();

      sprintf(hName,"hV1D");
      hV1D = (TH1F*) hE1D[i]->Clone(hName);
      hV1D->Reset();

      for(Int_t j=NbinsY;j>0;j--) {
	Float_t integral = 0.0;
	for(Int_t k=NbinsX;k>0;k--) {
	  integral += hE2D[i]->GetBinContent(k,j) * dx;
	  hV2D->SetBinContent(k,j,integral);
	}
      }

      Float_t integral = 0.0;
      for(Int_t k=NbinsX;k>0;k--) {
	integral += hE1D[i]->GetBinContent(k) * dx;
	hV1D->SetBinContent(k,integral);
      }

    }
  
  }

  // Now, combine the electric field components into the total |E|
  // and calculate ionization probability for He:
  // Outter Helium electron
  Float_t Eion0 = 24.59 * PUnits::eV;
  Float_t Z     = 1;

  if(opt.Contains("He2")) { // inner He electron
    Eion0 = 54.4 * PUnits::eV;
    Z = 2;
  } else if (opt.Contains("H")) { // inner He electron
    Eion0 = 13.6 * PUnits::eV;
    Z = 1;
  } 
  
  Float_t z10 = 99.0;
  Float_t z100 = 99.0;
    
    
  TH2F *hETotal2D = (TH2F*) hE2D[0]->Clone("hETotal2D");
  hETotal2D->Reset();
  TH2F *hIonRate2D = (TH2F*) hE2D[0]->Clone("hIonRate2D");
  hIonRate2D->Reset();
  TH2F *hIonProb2D = (TH2F*) hE2D[0]->Clone("hIonProb2D");
  hIonProb2D->Reset();
  {
    Int_t NbinsX = hE2D[0]->GetNbinsX();
    Int_t NbinsY = hE2D[0]->GetNbinsY();

    Float_t dx = pData->GetDX(0);
    // cout << "dz = " << dx << endl;
    
    for(Int_t j=0;j<NbinsY;j++) {
    
      Float_t integral = 0.0;

      for(Int_t k=NbinsX;k>0;k--) {
	Float_t E1 = hE2D[0]->GetBinContent(k,j);
	Float_t E2 = hE2D[1]->GetBinContent(k,j);
	Float_t E3 = hE2D[2]->GetBinContent(k,j);
	Float_t E  = TMath::Sqrt(E1*E1+E2*E2+E3*E3);
    
	hETotal2D->SetBinContent(k,j,E);
    
	E *= E0;

	if(E<10) continue;

	// Float_t IonRate = (PFunc::ADK(E,Eion0,Z,l,m)/PUnits::atomictime)*PUnits::femtosecond;
	Float_t IonRate = PFunc::ADK_ENG(E,Eion0,Z) * PUnits::femtosecond;
	hIonRate2D->SetBinContent(k,j,IonRate);

	if(integral==-99) continue;

	integral += IonRate * (dx * timedepth / PUnits::femtosecond);	
	if(integral>1) {
	  hIonProb2D->SetBinContent(k,j,100.0);	  
	  integral = -99;
	} else {
	  hIonProb2D->SetBinContent(k,j,integral*100);
	}
	
	if(j==NbinsY/2) {  // on axis
	  if(integral>0.1 && z10>0) {
	    z10 = hIonRate2D->GetXaxis()->GetBinCenter(k);
	  }
	  
	  if(integral>1.0 && z100>0) {
	    z100 = hIonRate2D->GetXaxis()->GetBinCenter(k);
	  }
	  
	  
	}
	
      }
    }
  }
  
  if(opt.Contains("units") && n0 ) {
    z10 *= skindepth / PUnits::um; 
    z100 *= skindepth / PUnits::um; 
  }
  cout << Form("  -->  10%% ionization found at z = %.4f",z10) << endl;
  cout << Form("  --> 100%% ionization found at z = %.4f",z100) << endl;
  
  hETotal2D->GetZaxis()->SetTitle("E [E_{0}]");
  hIonRate2D->GetZaxis()->SetTitle("W_{ADK} [fs^{-1}]");
  hIonProb2D->GetZaxis()->SetTitle("Prob_{ADK} [%]");

  
  TH1F *hETotal1D = (TH1F*) hE1D[0]->Clone("hETotal1D");
  hETotal1D->Reset();
  TH1F *hIonRate1D = (TH1F*) hE1D[0]->Clone("hIonRate1D");
  hIonRate1D->Reset();
  TH1F *hIonProb1D = (TH1F*) hE1D[0]->Clone("hIonProb1D");
  hIonProb1D->Reset();
  {

    Float_t integral = 0.0;
    Float_t dx = pData->GetDX(0);
    
    Int_t NbinsX = hE1D[0]->GetNbinsX();
    for(Int_t j=NbinsX;j>=0;j--) {
      Float_t E1 = hE1D[0]->GetBinContent(j);
      Float_t E2 = hE1D[1]->GetBinContent(j);
      Float_t E3 = hE1D[2]->GetBinContent(j);
      Float_t E  = TMath::Sqrt(E1*E1+E2*E2+E3*E3);
    
      hETotal1D->SetBinContent(j,E);
    
      E *= E0;
    
      if(E<10) continue;

      // Float_t IonRate = (PFunc::ADK(E,Eion0,Z,l,m)/PUnits::atomictime)*PUnits::femtosecond;
      Float_t IonRate = PFunc::ADK_ENG(E,Eion0,Z) * PUnits::femtosecond;
      hIonRate1D->SetBinContent(j,IonRate);

      if(integral==-99) continue;
	 
      integral += IonRate * (dx * timedepth / PUnits::femtosecond);	
      if(integral>1) {
	hIonProb1D->SetBinContent(j,100.0);	  
	integral = -99;
      } else {
	hIonProb1D->SetBinContent(j,integral*100);
      } 
      
    }
  } 
  hETotal1D->GetYaxis()->SetTitle("E [E_{0}]");
  hIonRate1D->GetYaxis()->SetTitle("W_{ADK} [fs^{-1}]");
  hIonProb1D->GetYaxis()->SetTitle("Prob_{ADK} [%]");


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
	Float_t binSize = (xMax - xMin)/NbinsX;  // bin size in um.
	
	Float_t dV = skindepth * skindepth * skindepth;
	Float_t lightspeed =  PConst::c_light / (PUnits::um/PUnits::femtosecond);
	
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
            
      for(Int_t j=0;j<=hE2D[i]->GetNbinsX();j++) {
	for(Int_t k=0;k<=hE2D[i]->GetNbinsY();k++) {
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
	hV2D->SetBins(NbinsX,xMin,xMax,NbinsY,ymin,ymax);
	hETotal2D->SetBins(NbinsX,xMin,xMax,NbinsY,ymin,ymax);
	hIonRate2D->SetBins(NbinsX,xMin,xMax,NbinsY,ymin,ymax);
	hIonProb2D->SetBins(NbinsX,xMin,xMax,NbinsY,ymin,ymax);
	hV1D->SetBins(NbinsX,xMin,xMax);
	hETotal1D->SetBins(NbinsX,xMin,xMax);
	hIonRate1D->SetBins(NbinsX,xMin,xMax);
	hIonProb1D->SetBins(NbinsX,xMin,xMax);
	for(Int_t j=0;j<=NbinsX;j++) {
	  for(Int_t k=0;k<=NbinsY;k++) {
	    hV2D->SetBinContent(j,k, hV2D->GetBinContent(j,k) * E0 * skindepth / (PUnits::MV));
	    hETotal2D->SetBinContent(j,k, hETotal2D->GetBinContent(j,k) * ( E0 / (PUnits::GV/PUnits::m) ) );
	  }
	  hV1D->SetBinContent(j, hV1D->GetBinContent(j) * ( E0 * skindepth / (PUnits::MV) ) );
	  hETotal1D->SetBinContent(j, hETotal1D->GetBinContent(j) * ( E0 / (PUnits::GV/PUnits::m) ) );
	}
	
	if(pData->IsCyl()) {
	  hV2D->GetYaxis()->SetTitle("r [#mum]");      
	  hETotal2D->GetYaxis()->SetTitle("r [#mum]");      
	} else {
	  hV2D->GetYaxis()->SetTitle("y [#mum]");      
	  hETotal2D->GetYaxis()->SetTitle("y [#mum]");      
	
	}
   
	if(opt.Contains("comov")) {
	  hV2D->GetXaxis()->SetTitle("#zeta [#mum]");
	  hV1D->GetXaxis()->SetTitle("#zeta [#mum]");
	  hETotal2D->GetXaxis()->SetTitle("#zeta [#mum]");
	  hETotal1D->GetXaxis()->SetTitle("#zeta [#mum]");
	} else {
	  hV2D->GetXaxis()->SetTitle("z [#mum]");
	  hV2D->GetXaxis()->SetTitle("z [#mum]");
	  hETotal2D->GetXaxis()->SetTitle("z [#mum]");
	  hETotal1D->GetXaxis()->SetTitle("z [#mum]");
	}
	 
	hV2D->GetZaxis()->SetTitle("#Psi-#Psi_{t} [MV]");
	hV1D->GetYaxis()->SetTitle("#Psi-#Psi_{t} [MV]");
	hETotal2D->GetZaxis()->SetTitle("E [GV/m]");
	hETotal1D->GetYaxis()->SetTitle("E [GV/m]");
      }
    }
  }

  // More calculations
  // get the hE1D[0] crossings
  const Int_t NcrossMax = 100;
  Float_t Ecross[NcrossMax] = {0.0};
  Float_t Eextr[NcrossMax]  = {0.0};
  Int_t Ncross = 0;
  
  for(Int_t ip=hE1D[0]->GetNbinsX()-10;ip>1;ip--) {
    
    Float_t Z2 = hE1D[0]->GetBinCenter(ip-1);
    Float_t E1 = hE1D[0]->GetBinContent(ip);
    Float_t E2 = hE1D[0]->GetBinContent(ip-1);
    Float_t Z1 = hE1D[0]->GetBinCenter(ip);
    
    // cout << Form("Z1 = %6.4f  Z2 = %6.4f   E1 = %6.4f   E2 = %6.4f", Z1, Z2, E1, E2) << endl; 
    
    if(E1*E2 >= 0) { // No change of sign means we are in a side of the zero axis.
      if(fabs(E2)>fabs(Eextr[Ncross])) {
	Eextr[Ncross] = E2;
      } 
    }
    
    if(E1*E2 < 0) { // change of sign means a crossing!
      
      // The next crossing has to be far enough from the previous one:
      Float_t zcross =  -E1 * ( (Z2-Z1)/(E2-E1) ) + Z1;
      if(Ncross>0 && fabs(Ecross[Ncross-1]-zcross)<1 ) continue;	
      // cout << " CROSS! " << endl;
      
      // add the point
      Ecross[Ncross] = zcross;
      Ncross++;
    }
  }
  
  cout << "  -> Number of crossings for field " << 0 << " : " << Ncross << endl;
  for(Int_t ic=0;ic<Ncross;ic++) {
    cout << Form(" %2i:  zeta = %6.4f  E = %6.4f", ic, Ecross[ic], Eextr[ic]) << endl; 
  }  
  
  Float_t R = 1.0;
  if(Ncross>=2) {
    R = fabs(Eextr[1]/Eextr[0]);
  }
  cout << Form("  -> Transformer ratio = %.2f ",R) << endl;

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

  hETotal2D->GetYaxis()->SetRangeUser(yMin,yMax);

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
    if(i==0) {
      if(Max[i]<1) {
	Max[i] = 1.1*Base;
      } else if(Max[i]<2) {
	Min[i] = 2*Base - Max[i];
      } else {
	Max[i] = 0.8*hDen2D[i]->GetMaximum(); // enhance plasma contrast.
      }
    }
        
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
    // hE2D[i]->GetZaxis()->SetRangeUser(Emin[i],Emax[i]); 
    hE2D[i]->GetZaxis()->SetRangeUser(Emin[0],Emax[0]); 
  }

  Float_t ETmin = 0.001;  
  Float_t ETmax = hETotal2D->GetMaximum();
  hETotal2D->GetZaxis()->SetRangeUser(ETmin,ETmax);
  
  // Potential 
  if(opt.Contains("units") && n0) {
    trapPotential *=  ( E0 * skindepth / (PUnits::MV) ); 
  }  
  
  Float_t Vmin = hV1D->GetMinimum();    
  { // Shift potential
    Int_t   NbinsX = hV2D->GetNbinsX(); 
    Int_t   NbinsY = hV2D->GetNbinsY(); 
    for(Int_t j=0;j<NbinsX;j++) {
      for(Int_t k=0;k<NbinsY;k++) {
  	hV2D->SetBinContent(j,k, hV2D->GetBinContent(j,k) - Vmin -trapPotential);
      }
      hV1D->SetBinContent(j, hV1D->GetBinContent(j) - Vmin -trapPotential);
    }
  }

  Vmin = hV1D->GetMinimum();   
  Float_t Vmax = hV1D->GetMaximum();    
  if(Vmax<0.1) Vmax = 0.1;

  // Dynamic potential palette
  const Int_t potPNRGBs = 6;
  const Int_t potPNCont = 64;
  Float_t zeroPos = -Vmin/(Vmax-Vmin);

  Double_t potPStops[potPNRGBs] = { 0.00, zeroPos-3.0/potPNCont,zeroPos-1.0/potPNCont, zeroPos, zeroPos+3.0/potPNCont, 1.00 };
  Double_t potPRed[potPNRGBs]   = { 0.518, 0.965, 0.90,0.90, 0.498, 0.106 };
  Double_t potPGreen[potPNRGBs] = { 0.078, 0.925, 0.90,0.90, 0.718, 0.078 };
  Double_t potPBlue[potPNRGBs]  = { 0.106, 0.353, 0.90,0.90, 0.780, 0.518 };
   
  PPalette * potentialPalette = (PPalette*) gROOT->FindObject("rbowinv");
  potentialPalette->CreateGradientColorTable(potPNRGBs, potPStops, 
					     potPRed, potPGreen, potPBlue, potPNCont);
  
  // Extract contours
  TCanvas* c = new TCanvas("c","Contour List",0,0,600,600);
  c->cd();
  
  // Potential
  TH2F *hV2Dc = (TH2F*) hV2D->Clone("hV2Dc");
  const Int_t Ncontours = 25;
  Double_t contours[Ncontours];
  for(Int_t i=0; i<Ncontours; i++) {
    contours[i] = i*(trapPotential/5.0) - trapPotential; 
  }
  hV2Dc->SetContour(Ncontours, contours);
  hV2Dc->Draw("cont list");
  
  c->Update();
  TObjArray *contsV2D = (TObjArray*) gROOT->GetListOfSpecials()->FindObject("contours");
  TClonesArray graphsV2D("TGraph",Ncontours);
  TClonesArray graphsV2D_main("TGraph",Ncontours);
  {
    Int_t ncontours = contsV2D->GetSize();
    TList* clist = NULL;
    Int_t nGraphs = 0;
    Int_t nGraphs_main = 0;
    TGraph *gr = NULL;
    for(Int_t i = 0; i < ncontours; i++){
      if(i==0) continue;
      
      clist = (TList*) contsV2D->At(i);
    
      for(Int_t j = 0 ; j < clist->GetSize(); j++) {
	gr = (TGraph*) clist->At(j);
	if(!gr) continue;
      
	gr->SetLineWidth(1);
	gr->SetLineColor(kGray+1);
      
	if( !((i-1)%5) && (contours[i]>=0)) {
	  gr->SetLineWidth(2);
	  gr->SetLineColor(kGray+2);
	} 
      	new(graphsV2D[nGraphs]) TGraph(*gr) ;
	nGraphs++;
      
	if(i-1==5) {
	  TGraph *grm = new(graphsV2D_main[nGraphs_main]) TGraph(*gr);
	  nGraphs_main++;
	  grm->SetLineWidth(1);
	  grm->SetLineStyle(2);
	  grm->SetLineColor(kGray+3);	  
	}

	if(i-1==10) {
	  TGraph *grm = new(graphsV2D_main[nGraphs_main]) TGraph(*gr);
	  nGraphs_main++;
	  grm->SetLineWidth(1);
	  grm->SetLineStyle(2);
	  grm->SetLineColor(kGray+3);	  
	}

      }
    }
  }
  
  // Ion probability
  // hIonRate2D->GetZaxis()->SetRangeUser(0.00501,80);
  
  // TH2F *hIonRate2Dc = (TH2F*) hIonRate2D->Clone("hIonRate2Dc");
  // const Int_t NcontI = 4;
  // Double_t contI[NcontI] = {0.01,0.1,1.0,10.0};
  // //const Int_t NcontI = 1;
  // //Float_t contI[NcontI] = {0.1};
  // hIonRate2Dc->SetContour(NcontI, contI);
  // hIonRate2Dc->Draw("cont list");
  
  Float_t IPmin = 0.001;
  Float_t IPmax = 99.999;
  hIonProb2D->GetZaxis()->SetRangeUser(IPmin,IPmax);
  
  TH2F *hIonProb2Dc = (TH2F*) hIonProb2D->Clone("hIonProb2Dc");
  const Int_t NcontI = 2;
  Double_t contI[NcontI] = {10.0,100.0};
  //const Int_t NcontI = 1;
  //Float_t contI[NcontI] = {0.1};
  hIonProb2Dc->SetContour(NcontI, contI);
  hIonProb2Dc->Draw("cont list");
  
  c->Update();
  TObjArray *contsI2D = (TObjArray*) gROOT->GetListOfSpecials()->FindObject("contours");
  TClonesArray graphsI2D("TGraph",NcontI);
  TClonesArray graphsI2D_main("TGraph",NcontI);
  {
    Int_t ncontours = contsI2D->GetSize();
    TList* clist = NULL;
    Int_t nGraphs = 0;
    Int_t nGraphs_main = 0;
    TGraph *gr = NULL;
    for(Int_t i = 0; i < ncontours; i++){
      clist = (TList*) contsI2D->At(i);
      
      for(Int_t j = 0 ; j < clist->GetSize(); j++) {
	gr = (TGraph*) clist->At(j);
	if(!gr) continue;
	
	if( i==0 ) {
	  TGraph *grm = new(graphsI2D_main[nGraphs_main]) TGraph(*gr);
	  nGraphs_main++;
	  grm->SetLineWidth(1);
	  grm->SetLineStyle(2);
	  grm->SetLineColor(kGray+1);

	  TGraph *gr2 = new(graphsI2D[nGraphs]) TGraph(*gr) ;
	  nGraphs++;
	  if(j==0) {
	    gr2->SetLineWidth(2);
	    gr2->SetLineStyle(2);
	    //	    gr2->SetLineColor(PGlobals::elecLine);
	    gr2->SetLineColor(kGray+2);
	  } else {
	    gr2->SetLineWidth(2);
	    gr2->SetLineStyle(1);
	    gr2->SetLineColor(kGray+2);
	  }
	  
	} 
		
      }
    }
  }
    
  // "Axis range" in Osiris units:
  Float_t ylow  = hDen2D[0]->GetYaxis()->GetBinLowEdge(FirstyBin);
  Float_t yup = hDen2D[0]->GetYaxis()->GetBinUpEdge(LastyBin);
  Float_t xmin = hDen2D[0]->GetXaxis()->GetXmin();
  Float_t xmax = hDen2D[0]->GetXaxis()->GetXmax();

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
  
  if(opt.Contains("units") && n0) {
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
  lineStartNeutral->SetLineColor(kGray+3);
  lineStartNeutral->SetLineStyle(3);
  lineStartNeutral->SetLineWidth(2);

  //  cout << "End plasma = " << zEndNeutral << endl;
  TLine *lineEndNeutral = new TLine(zEndNeutral,yMin,zEndNeutral,yMax);
  lineEndNeutral->SetLineColor(kGray+3);
  lineEndNeutral->SetLineStyle(3);
  lineEndNeutral->SetLineWidth(2);


  // Plotting
  // -----------------------------------------------

  // Canvas setup
  TCanvas *C = new TCanvas("C","Rake 5 panel",1050,1680);
  

  // Palettes setup
  TExec *exPlasma = new TExec("exPlasma","plasmaPalette->cd();");
  TExec *exElec   = new TExec("exElec","redelectronPalette->cd();");
  TExec *exHot    = new TExec("exHot","hotPalette->cd();");
  TExec *exField  = new TExec("exField","rbowPalette->cd();");
  TExec *exFieldT = new TExec("exFieldT","redPalette->cd();");
  TExec *exIonP   = new TExec("exIonP","grayPalette->cd();");
  TExec *exPot    = new TExec("exPot","rbowinvPalette->cd();");
     
  // Actual Plotting!
  // ------------------------------------------------------------

  // Output file
  TString fOutName = Form("./%s/Plots/Rake5Panel2D/Rake5Panel2D",pData->GetPath().c_str());
  fOutName += Form("-%s_%i",pData->GetName(),time);

  // Setup Pad layout:
  const Int_t NPad = 5;
  TPad *pad[NPad];
  TH2F *hFrame[NPad];
  TString sLabels[] = {"(e)","(d)","(c)","(b)","(a)"};
  // Text objects
  TPaveText **textLabel = new TPaveText*[NPad];
 
  Float_t lMargin = 0.15;
  Float_t rMargin = 0.18;
  Float_t bMargin = 0.10;
  Float_t tMargin = 0.04;
  Float_t vSpacing = 0.008;
  PGlobals::CanvasPartition(C,NPad,lMargin,rMargin,bMargin,tMargin,vSpacing);
  
  // Define the frames for plotting
  Int_t fonttype = 43;
  Int_t fontsize = 40;
  Int_t tfontsize = 44;
  Float_t txoffset = 3.5;
  Float_t lxoffset = 0.02;
  Float_t tyoffset = 2.5;
  Float_t lyoffset = 0.01;
  Float_t tylength = 0.02;
  Float_t txlength = 0.04;
  for(Int_t i=0;i<NPad;i++) {
    char name[16];
    sprintf(name,"pad_%i",i);
    pad[i] = (TPad*) gROOT->FindObject(name);
    pad[i]->SetFrameLineWidth(3);  
    pad[i]->SetTickx(1);
    pad[i]->SetTicky(0);

    sprintf(name,"hFrame_%i",i);
    hFrame[i] = (TH2F*) gROOT->FindObject(name);
    if(hFrame[i]) delete hFrame[i];
    hFrame[i] = (TH2F*) hDen2D[0]->Clone(name);
    hFrame[i]->Reset();
    
    Float_t xFactor = pad[0]->GetAbsWNDC()/pad[i]->GetAbsWNDC();
    Float_t yFactor = pad[0]->GetAbsHNDC()/pad[i]->GetAbsHNDC();

    // Format for y axis
    hFrame[i]->GetYaxis()->SetTitleFont(fonttype);
    hFrame[i]->GetYaxis()->SetTitleSize(tfontsize);
    hFrame[i]->GetYaxis()->SetTitleOffset(tyoffset);
    hFrame[i]->GetYaxis()->SetLabelFont(fonttype);
    hFrame[i]->GetYaxis()->SetLabelSize(fontsize);
    hFrame[i]->GetYaxis()->SetLabelOffset(lyoffset);

    hFrame[i]->GetYaxis()->SetTickLength(xFactor*tylength/yFactor);

    // Format for x axis
    hFrame[i]->GetXaxis()->SetTitleFont(fonttype);
    hFrame[i]->GetXaxis()->SetTitleSize(tfontsize+8);
    hFrame[i]->GetXaxis()->SetTitleOffset(txoffset);
    hFrame[i]->GetXaxis()->SetLabelFont(fonttype);
    hFrame[i]->GetXaxis()->SetLabelSize(fontsize+4);
    hFrame[i]->GetXaxis()->SetLabelOffset(lxoffset);
    
    hFrame[i]->GetXaxis()->SetTickLength(yFactor*txlength/xFactor);      

    // Labels for the frames

  }

  C->cd(0);
  pad[4]->Draw();
  pad[4]->cd(); // <---------------------------------------------- 1st panel
  if(opt.Contains("logz")) {
    pad[4]->SetLogz(1);
  } else {
    pad[4]->SetLogz(0);
  }

  hFrame[4]->Draw("col");

  // Injected electrons if any
  if(Nspecies>=3) {
    if(hDen2D[2]) {
      exHot->Draw();
      hDen2D[2]->GetZaxis()->SetNdivisions(505);
      hDen2D[2]->Draw("colz same");
    }
  }
  
  // Plasma
  hDen2D[0]->GetZaxis()->SetNdivisions(505);
  hDen2D[0]->GetZaxis()->SetTitleFont(fonttype);
  exPlasma->Draw();
  hDen2D[0]->Draw("colz same");
  
  // Beam driver.
  if(hDen2D[1]) {
    hDen2D[1]->GetZaxis()->SetNdivisions(505);
    exElec->Draw();
    hDen2D[1]->Draw("colz same");
  }


  // ADK MAIN contour  
  for(Int_t i=0;i<graphsI2D_main.GetEntriesFast();i++) {
    TGraph *gr = (TGraph*) graphsI2D_main.At(i);
    if(!gr) continue;
    gr->Draw("C");
  }  

  // PSI MAIN contour  
  for(Int_t i=0;i<graphsV2D_main.GetEntriesFast();i++) {
    TGraph *gr = (TGraph*) graphsV2D_main.At(i);
    if(!gr) continue;
    gr->Draw("C");
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
  pad[4]->Update();
  Float_t y1 = pad[4]->GetBottomMargin();
  Float_t y2 = 1 - pad[4]->GetTopMargin();
  Float_t x1 = pad[4]->GetLeftMargin();
  Float_t x2 = 1 - pad[4]->GetRightMargin();
  Float_t gap = 0.005;

  TPaletteAxis *palette = NULL;
  if(Nspecies>=3) {
    if(hDen2D[2]) {
      palette = (TPaletteAxis*)hDen2D[2]->GetListOfFunctions()->FindObject("palette");
    }
  }
 
    if(palette) {
    palette->SetY2NDC(y2 - gap);
    palette->SetY1NDC(0.66*(y1+y2) + gap);
    palette->SetX1NDC(x2 + 0.005);
    palette->SetX2NDC(x2 + 0.03);
    //  palette->SetTitleOffset(0.85);
    palette->SetTitleOffset(999.9);
    palette->SetTitleSize(tfontsize);
    palette->SetLabelFont(fonttype);
    palette->SetLabelSize(fontsize);
    if(opt.Contains("logz")) 
      palette->SetLabelOffset(0);
    else
      palette->SetLabelOffset(lyoffset);
    palette->SetBorderSize(2);
    palette->SetLineColor(1);
  }
 
  palette = (TPaletteAxis*)hDen2D[0]->GetListOfFunctions()->FindObject("palette");
  if(palette) {
    palette->SetY2NDC(0.66*(y1+y2) - gap);
    palette->SetY1NDC(0.33*(y1+y2) + gap);
    palette->SetX1NDC(x2 + 0.005);
    palette->SetX2NDC(x2 + 0.03);
    palette->SetTitleOffset(tyoffset);
    // palette->SetTitleOffset(999.9);
    palette->SetTitleSize(tfontsize);
    palette->SetLabelFont(fonttype);
    palette->SetLabelSize(fontsize);
    if(opt.Contains("logz")) 
      palette->SetLabelOffset(0);
    else
      palette->SetLabelOffset(lyoffset);
    palette->SetBorderSize(2);
    palette->SetLineColor(1);
  } 
  
  palette = (TPaletteAxis*)hDen2D[1]->GetListOfFunctions()->FindObject("palette");
  if(palette) {
    palette->SetY2NDC(0.33*(y1+y2) - gap);
    palette->SetY1NDC(y1 + gap);
    palette->SetX1NDC(x2 + 0.005);
    palette->SetX2NDC(x2 + 0.03);
    //palette->SetTitleOffset(tyoffset);
    palette->SetTitleOffset(999.9);
    palette->SetTitleSize(tfontsize);
    palette->SetLabelFont(fonttype);
    palette->SetLabelSize(fontsize);
    if(opt.Contains("logz")) 
      palette->SetLabelOffset(0);
    else
      palette->SetLabelOffset(lyoffset);
    palette->SetBorderSize(2);
    palette->SetLineColor(1);
  }

  // 1D charge density plots:
  Float_t yaxismin  =  pad[4]->GetUymin();
  Float_t yaxismax  =  pad[4]->GetUymin() + 0.33*(pad[4]->GetUymax() - pad[4]->GetUymin()) - 0.00;
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
  // // // PGlobals::SetH1Style(hDen1D[0],1);
  // hDen1D[0]->Draw("same C");  
  
  
  if(opt.Contains("curr")) {
    hCur1D[1]->SetLineWidth(2);
    hCur1D[1]->SetLineColor(PGlobals::elecLine);
    hCur1D[1]->Draw("same C");
  } else {
    hDen1D[1]->SetLineWidth(2);
    hDen1D[1]->SetLineColor(PGlobals::elecLine);
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
    axis = new TGaxis(xMax-xRange/6.0,yMin,
		      xMax-xRange/6.0,yaxismax,
		      0.001,curmax,503,"+LS");
    
    axis->SetLineWidth(1);
    axis->SetLineColor(kGray+3);//PGlobals::elecLine);
    axis->SetLabelColor(kGray+3);//PGlobals::elecLine);
    axis->SetLabelSize(0.06);
    axis->SetLabelOffset(0.01);
    axis->SetLabelFont(42);
    axis->SetTitleColor(kGray+3);//PGlobals::elecLine);
    axis->SetTitleSize(0.06);
    axis->SetTitleOffset(0.6);
    axis->SetTitleFont(42);
    axis->SetTickSize(0.03);
    axis->SetTitle("I [kA]");
    axis->CenterTitle();
    axis->SetNdivisions(505);
    
    axis->Draw();
  }

  
  TPaveText *textTime = new TPaveText(xMin+0.05*xRange, yMax-0.15*yRange, xMin+0.4*xRange, yMax-0.05*yRange);
  //x2-0.17,y2-0.12,x2-0.02,y2-0.02,"NDC");
  PGlobals::SetPaveTextStyle(textTime,12); 
  char ctext[128];
  if(opt.Contains("units") && n0) 
    sprintf(ctext,"z = %5.1f #mum", Time * skindepth / PUnits::um);
  else
    sprintf(ctext,"t = %5.1f #omega_{p}^{-1}",Time);
  textTime->SetTextFont(42);
  textTime->AddText(ctext);
  if(!opt.Contains("notime"))
    textTime->Draw();
  // textDen->Draw();
  // if(opt.Contains("units"))
  //   textWav->Draw();
  
  textLabel[4] = new TPaveText(xMax - 0.07 * xRange, yMax - 0.2 * yRange, xMax - 0.01 * xRange, yMax - 0.05 * yRange);
  PGlobals::SetPaveTextStyle(textLabel[4],12); 
  textLabel[4]->SetTextFont(42);
  textLabel[4]->AddText(sLabels[4]);
  textLabel[4]->Draw();
  
  pad[4]->RedrawAxis(); 


  C->cd(0);
  pad[3]->Draw();
  pad[3]->cd(); // <---------------------------------------------------------- 2nd panel
  // if(opt.Contains("logz")) {
  //   pad[3]->SetLogz(1);
  // } else {
  //   pad[3]->SetLogz(0);
  // }

  hFrame[3]->Draw("col");

  //  hE2D[0]->GetZaxis()->SetNdivisions(505);
  hE2D[0]->GetZaxis()->SetTitleFont(fonttype);
  Float_t xFactor = pad[0]->GetAbsWNDC()/pad[3]->GetAbsWNDC();
  Float_t yFactor = pad[0]->GetAbsHNDC()/pad[3]->GetAbsHNDC();
  hE2D[0]->GetZaxis()->SetTickLength(xFactor*tylength/yFactor);
  
  exField->Draw();
  hE2D[0]->Draw("colz same");

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
   
  // Fit the E1D in the E2D pad:
  TLine *lineV = NULL;
  {
    Float_t rightmin = Emin[0];
    Float_t rightmax = Emax[0];
    Float_t slope = (yMax-yMin)/(rightmax-rightmin);
    
    for(Int_t j=0;j<hE1D[0]->GetNbinsX();j++) {
      hE1D[0]->SetBinContent(j+1,slope*(hE1D[0]->GetBinContent(j+1)-rightmin)+yMin);
    }
    
    lineV = new TLine(xMin,slope*(-rightmin)+yMin,
		      xMax,slope*(-rightmin)+yMin);
  }

  lineV->SetLineColor(kOrange+10);
  lineV->SetLineStyle(2);
  lineV->Draw();


  hE1D[0]->SetLineStyle(1);
  hE1D[0]->SetLineWidth(2);
  hE1D[0]->SetLineColor(kOrange+10);

  hE1D[0]->Draw("sameL");
  
  pad[3]->Update();
  
  y1 = pad[3]->GetBottomMargin();
  y2 = 1 - pad[3]->GetTopMargin();
  x1 = pad[3]->GetLeftMargin();
  x2 = 1 - pad[3]->GetRightMargin();
  
  palette = (TPaletteAxis*) hE2D[0]->GetListOfFunctions()->FindObject("palette");  
  if(palette) {
    palette->SetY2NDC(y2 - gap);
    palette->SetY1NDC(y1 + gap);
    palette->SetX1NDC(x2 + 0.005);
    palette->SetX2NDC(x2 + 0.03);
    palette->SetTitleOffset(tyoffset);
    palette->SetTitleSize(tfontsize);
    palette->SetLabelFont(fonttype);
    palette->SetLabelSize(fontsize);
    palette->SetLabelOffset(lyoffset);
    palette->SetBorderSize(2);
    palette->SetLineColor(1);
  }
   
  pad[3]->RedrawAxis(); 
  
  
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
 

  textLabel[3] = new TPaveText(xMax - 0.07 * xRange, yMax - 0.2 * yRange, xMax - 0.01 * xRange, yMax - 0.05 * yRange);
  PGlobals::SetPaveTextStyle(textLabel[3],12); 
  textLabel[3]->SetTextFont(42);
  textLabel[3]->AddText(sLabels[3]);
  textLabel[3]->Draw();

 
  C->cd(0);
  pad[2]->Draw();
  pad[2]->cd(); // <---------------------------------------------------------- 3rd panel
  // if(opt.Contains("logz")) {
  //   pad[2]->SetLogz(1);
  // } else {
  //   pad[2]->SetLogz(0);
  // }

  hFrame[2]->Draw("col");

  //  hETotal2D->GetZaxis()->SetNdivisions(505);
  hETotal2D->GetZaxis()->SetTitleFont(fonttype);
  xFactor = pad[0]->GetAbsWNDC()/pad[2]->GetAbsWNDC();
  yFactor = pad[0]->GetAbsHNDC()/pad[2]->GetAbsHNDC();
  hETotal2D->GetZaxis()->SetTickLength(xFactor*tylength/yFactor);

  exFieldT->Draw();
  
  hETotal2D->Draw("colz same");

  // ADK contours  
  // for(Int_t i=0;i<graphsI2D.GetEntriesFast();i++) {
  //   TGraph *gr = (TGraph*) graphsI2D.At(i);
  //   if(!gr) continue;
  //   gr->Draw("C");
  // }

  // Fit the E1D in the E2D pad:
  Float_t HeIon  = 92.75;  // GV/m
  Float_t HeIon2 = 234.96;  // GV/m
  Float_t HIon   = 33.8;  // GV/m 
  if(!opt.Contains("units") || !n0 ) {
    HeIon  /= ( E0 / (PUnits::GV/PUnits::m));
    HeIon2 /= ( E0 / (PUnits::GV/PUnits::m));
    HIon   /= ( E0 / (PUnits::GV/PUnits::m));
  }
  TLine *lineET = NULL;
  TLine *lineET2 = NULL;
  TLine *lineETH = NULL;
  {
    Float_t rightmin = ETmin;
    Float_t rightmax = ETmax;
    Float_t slope = (yMax-yMin)/(rightmax-rightmin);
    
    for(Int_t j=0;j<hETotal1D->GetNbinsX();j++) {
      hETotal1D->SetBinContent(j+1,slope*(hETotal1D->GetBinContent(j+1)-rightmin)+yMin);
    }
    
    lineET  = new TLine(xMin,slope*(HeIon-rightmin)+yMin,
			xMax,slope*(HeIon-rightmin)+yMin);
    lineET2 = new TLine(xMin,slope*(HeIon2-rightmin)+yMin,
			xMax,slope*(HeIon2-rightmin)+yMin);

    lineETH = new TLine(xMin,slope*(HIon-rightmin)+yMin,
			xMax,slope*(HIon-rightmin)+yMin);
  

    lineET->SetLineColor(kOrange+2);
    lineET->SetLineStyle(2);
    lineET->Draw();


    lineET2->SetLineColor(kOrange+3);
    lineET2->SetLineStyle(2);
    lineET2->Draw();

    lineETH->SetLineColor(kAzure-4);
    lineETH->SetLineStyle(2);
    lineETH->Draw();

  }
  
  hETotal1D->SetLineStyle(1);
  hETotal1D->SetLineWidth(2);
  hETotal1D->SetLineColor(kOrange+10);

  hETotal1D->Draw("sameL");
    
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
  
  palette = (TPaletteAxis*) hETotal2D->GetListOfFunctions()->FindObject("palette");  
  if(palette) {
    palette->SetY2NDC(y2 - gap);
    palette->SetY1NDC(y1 + gap);
    palette->SetX1NDC(x2 + 0.005);
    palette->SetX2NDC(x2 + 0.03);
    palette->SetTitleOffset(tyoffset);
    palette->SetTitleSize(tfontsize);
    palette->SetLabelFont(fonttype);
    palette->SetLabelSize(fontsize);
    palette->SetLabelOffset(lyoffset);
    palette->SetBorderSize(2);
    palette->SetLineColor(1);
  }
   
  pad[2]->RedrawAxis(); 

  textLabel[2] = new TPaveText(xMax - 0.07 * xRange, yMax - 0.2 * yRange, xMax - 0.01 * xRange, yMax - 0.05 * yRange);
  PGlobals::SetPaveTextStyle(textLabel[2],12); 
  textLabel[2]->SetTextFont(42);
  textLabel[2]->AddText(sLabels[2]);
  textLabel[2]->Draw();



  C->cd(0);
  pad[1]->Draw();
  pad[1]->cd(); // <---------------------------------------------------------- 4th panel
  // if(opt.Contains("logz")) {
  //   pad[1]->SetLogz(1);
  // } else {
  //   pad[1]->SetLogz(0);
  // }

  hFrame[1]->Draw("col");

  hIonProb2D->GetZaxis()->SetNdivisions(505);
  hIonProb2D->GetZaxis()->SetTitleFont(fonttype);
  xFactor = pad[0]->GetAbsWNDC()/pad[1]->GetAbsWNDC();
  yFactor = pad[0]->GetAbsHNDC()/pad[1]->GetAbsHNDC();
  hIonProb2D->GetZaxis()->SetTickLength(xFactor*tylength/yFactor);
  
  exIonP->Draw();
  hIonProb2D->Draw("colz same");

  // ADK contours  
  for(Int_t i=0;i<graphsI2D.GetEntriesFast();i++) {
    TGraph *gr = (TGraph*) graphsI2D.At(i);
    if(!gr) continue;
    gr->Draw("C");
  }

  // Ion Prob contour plot
  // TH2F *hIonProb2Dc = (TH2F*) hIonProb2D->Clone("hIonProb2Dc");
  // {
  //   const Int_t Ncontours = 11;
  //   Double_t contours[Ncontours];
  //   for(Int_t i=0; i<Ncontours; i++) {
  //     contours[i] = i*(Ncontours-1);
  //   }
  //   hIonProb2Dc->SetContour(Ncontours, contours);
  // }
  // hIonProb2Dc->Draw("cont2 z same");

  {
    Float_t rightmin = IPmin;
    Float_t rightmax = IPmax;
    Float_t slope = (yMax-yMin)/(rightmax-rightmin);
    
    for(Int_t j=0;j<hIonProb1D->GetNbinsX();j++) {
      hIonProb1D->SetBinContent(j+1,slope*(hIonProb1D->GetBinContent(j+1)-rightmin)+yMin);
    }
   }
  
  hIonProb1D->SetLineStyle(1);
  hIonProb1D->SetLineWidth(2);
  hIonProb1D->SetLineColor(kOrange+10);

  hIonProb1D->Draw("sameL");
    
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
  
  pad[1]->Update();
  
  y1 = pad[1]->GetBottomMargin();
  y2 = 1 - pad[1]->GetTopMargin();
  x1 = pad[1]->GetLeftMargin();
  x2 = 1 - pad[1]->GetRightMargin();
  
  palette = (TPaletteAxis*) hIonProb2D->GetListOfFunctions()->FindObject("palette");  
  if(palette) {
    palette->SetY2NDC(y2 - gap);
    palette->SetY1NDC(y1 + gap);
    palette->SetX1NDC(x2 + 0.005);
    palette->SetX2NDC(x2 + 0.03);
    palette->SetTitleOffset(tyoffset);
    palette->SetTitleSize(tfontsize);
    palette->SetLabelFont(fonttype);
    palette->SetLabelSize(fontsize);
    palette->SetLabelOffset(lyoffset);
    palette->SetBorderSize(2);
    palette->SetLineColor(1);
  }
   
  pad[1]->RedrawAxis(); 

  textLabel[1] = new TPaveText(xMax - 0.07 * xRange, yMax - 0.2 * yRange, xMax - 0.01 * xRange, yMax - 0.05 * yRange);
  PGlobals::SetPaveTextStyle(textLabel[1],12); 
  textLabel[1]->SetTextFont(42);
  textLabel[1]->AddText(sLabels[1]);
  textLabel[1]->Draw();


  C->cd(0);
  pad[0]->Draw();
  pad[0]->cd(); // <--------------------------------------------------------- 5th panel
 
  hFrame[0]->Draw("col");

  hV2D->GetZaxis()->SetTitleFont(fonttype);
  xFactor = pad[0]->GetAbsWNDC()/pad[0]->GetAbsWNDC();
  yFactor = pad[0]->GetAbsHNDC()/pad[0]->GetAbsHNDC();
  hV2D->GetZaxis()->SetTickLength(xFactor*tylength/yFactor);
 
  exPot->Draw();
  
  hV2D->GetZaxis()->SetRangeUser(Vmin,Vmax);
  hV2D->Draw("col z same");
  
  for(Int_t i=0;i<graphsV2D.GetEntriesFast();i++) {
    TGraph *gr = (TGraph*) graphsV2D.At(i);
    if(!gr) continue;
    gr->Draw("C");
  }
  
  // PSI MAIN contour  
  // for(Int_t i=0;i<graphsV2D_main.GetEntriesFast();i++) {
  //   TGraph *gr = (TGraph*) graphsV2D_main.At(i);
  //   if(!gr) continue;
  //   gr->Draw("C");
  // }

  // Fit the V1D in the V2D pad:
  TLine *lineT = NULL;
  {
    Float_t rightmin = Vmin;
    Float_t rightmax = Vmax;
    Float_t slope = (yMax-yMin)/(rightmax-rightmin);
    
    for(Int_t j=0;j<hV1D->GetNbinsX();j++) {
      hV1D->SetBinContent(j+1,slope*(hV1D->GetBinContent(j+1)-rightmin)+yMin);
    }

    lineT = new TLine(xMin,slope*(-rightmin)+yMin,
		      xMax,slope*(-rightmin)+yMin);
    lineT->SetLineColor(kOrange+10);
    lineT->SetLineStyle(2);
    lineT->Draw();

  }

  hV1D->SetLineStyle(1);
  hV1D->SetLineWidth(2);
  //  hV1D->SetLineColor(PGlobals::elecLine);
  hV1D->SetLineColor(kOrange+10);

  hV1D->Draw("sameL");

  // ADK MAIN contour  
  // for(Int_t i=0;i<graphsI2D_main.GetEntriesFast();i++) {
  //   TGraph *gr = (TGraph*) graphsI2D_main.At(i);
  //   if(!gr) continue;
  //   gr->Draw("C");
  // }

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
  
    
  pad[0]->Update();

  y1 = pad[0]->GetBottomMargin();
  y2 = 1 - pad[0]->GetTopMargin();
  x1 = pad[0]->GetLeftMargin();
  x2 = 1 - pad[0]->GetRightMargin();
 
  palette = (TPaletteAxis*)hV2D->GetListOfFunctions()->FindObject("palette");
  if(palette) {
    palette->SetY2NDC(y2 - gap);
    palette->SetY1NDC(y1 + gap);
    palette->SetX1NDC(x2 + 0.005);
    palette->SetX2NDC(x2 + 0.03);
    palette->SetTitleOffset(tyoffset);
    // palette->SetTitleOffset(999.9);
    palette->SetTitleSize(tfontsize);
    palette->SetLabelFont(fonttype);
    palette->SetLabelSize(fontsize);
    palette->SetLabelOffset(lyoffset);
    palette->SetBorderSize(2);
    palette->SetLineColor(1);
  }



  pad[0]->RedrawAxis(); 

  textLabel[0] = new TPaveText(xMax - 0.07 * xRange, yMax - 0.2 * yRange, xMax - 0.01 * xRange, yMax - 0.05 * yRange);
  PGlobals::SetPaveTextStyle(textLabel[0],12); 
  textLabel[0]->SetTextFont(42);
  textLabel[0]->AddText(sLabels[0]);
  textLabel[0]->Draw();

  C->cd();

  // Print to a file
  PGlobals::imgconv(C,fOutName,opt);
  // ---------------------------------------------------------

  PGlobals::DestroyCanvases();
}
  
