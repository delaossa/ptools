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

void PlotRakeFancy2D( const TString &sim, Int_t time, Float_t zoom=2, Int_t Nbins=2, const TString &options="") {
  
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
      Double_t dx = (hCur1D[i]->GetBinLowEdge(1)-hCur1D[i]->GetBinLowEdge(NB+1))/NB;
      
      // hCur1D[i]->Scale(dx);
      Float_t Charge = hCur1D[i]->Integral() * dx; 
      
      cout << Form(" Integrated charge of specie %3i = %8.4f n0 * kp^-3",i,Charge) << endl;
    }
  }
  

  // Get electric fields 2D
  const Int_t Nfields = 1;
  TH2F **hE2D = new TH2F*[Nfields];
  TH1F **hE1D = new TH1F*[Nfields];
  TH1F **hV1D = new TH1F*[Nfields];
  for(Int_t i=0;i<Nfields;i++) {
    hE2D[i] = NULL;
    hE1D[i] = NULL;
    hV1D[i] = NULL;

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
    
    Int_t   NbinsX = hE1D[i]->GetNbinsX();
    Float_t dz = hE1D[i]->GetBinWidth(Nbins);
    sprintf(hName,"hV1D_%i",i);
    hV1D[i] = (TH1F*) hE1D[i]->Clone(hName);
    hV1D[i]->Reset();
    Double_t integral = 0.0;
    for(Int_t k=NbinsX;k>0;k--) {
      integral += - hE1D[i]->GetBinContent(k) * dz;
      hV1D[i]->SetBinContent(k,integral);
    }
  
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
      hV1D[i]->SetBins(NbinsX,xMin,xMax);
      
      for(Int_t j=0;j<hE2D[i]->GetNbinsX();j++) {
	for(Int_t k=0;k<hE2D[i]->GetNbinsY();k++) {
	  hE2D[i]->SetBinContent(j,k, hE2D[i]->GetBinContent(j,k) * ( E0 / (PUnits::GV/PUnits::m) ) );
	}
	hE1D[i]->SetBinContent(j, hE1D[i]->GetBinContent(j) * ( E0 / (PUnits::GV) ) );
	hV1D[i]->SetBinContent(j, hV1D[i]->GetBinContent(j) * ( E0 / (PUnits::GV) ) );
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
      
    }
  }


 // --------------------------------------------------- Vertical Zoom ------------
  
  Float_t yRange    = (hDen2D[0]->GetYaxis()->GetXmax() - hDen2D[0]->GetYaxis()->GetXmin())/zoom;
  Float_t midPoint = (hDen2D[0]->GetYaxis()->GetXmax() + hDen2D[0]->GetYaxis()->GetXmin())/2.;
  Float_t yMin = midPoint-yRange/2.0;
  Float_t yMax = midPoint+yRange/2.0;
  if(pData->IsCyl()) {
    yMin = pData->GetXMin(1);
    yMax = yRange;
  }

  //  cout << Form(" %f  %f   %f  %f",hDen2D[0]->GetYaxis()->GetXmin(),hDen2D[0]->GetYaxis()->GetXmax(),yMin,yMax) << endl;

  for(Int_t i=0;i<Nspecies;i++) {
    if(!hDen2D[i]) continue;
    hDen2D[i]->GetYaxis()->SetRangeUser(yMin,yMax);
  }
  
  for(Int_t i=0;i<Nfields;i++) {
    if(!hE2D[i]) continue;
    hE2D[i]->GetYaxis()->SetRangeUser(yMin,yMax);
  }

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
      if(Max[i]<2) {
	Min[i] = 2*Base - Max[i];
      } else {
	Max[i] = 0.1*hDen2D[i]->GetMaximum(); // enhance plasma contrast.
      }
    }
        
    if(i==1) Min[i] = 3.01E-1 * Base;
    if(i==2) {
      if(time<30)
	Min[i] = 6.01E-3 * Base;
      else
	Min[i] = 1.01E-1 * Base;	
    }    
    hDen2D[i]->GetZaxis()->SetRangeUser(Min[i],Max[i]);
  }
  
  // Dynamic plasma palette
  const Int_t plasmaDNRGBs = 3;
  const Int_t plasmaDNCont = 128;
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
  Double_t plasmaDRed[plasmaDNRGBs]   = { 0.04, 0.09, 1.00 };
  Double_t plasmaDGreen[plasmaDNRGBs] = { 0.04, 0.17, 1.00 };
  Double_t plasmaDBlue[plasmaDNRGBs]  = { 0.04, 0.32, 1.00 };
   
  PPalette * plasmaPalette = (PPalette*) gROOT->FindObject("plasma");
  plasmaPalette->CreateGradientColorTable(plasmaDNRGBs, plasmaDStops, 
					  plasmaDRed, plasmaDGreen, plasmaDBlue, plasmaDNCont);


  const Int_t redelectronNRGBs = 4;
  const Int_t redelectronNCont = 128;
  Double_t redelectronStops[redelectronNRGBs] = { 0.30, 0.40, 0.60, 1.00};
  Double_t redelectronRed[redelectronNRGBs] =   { 0.09, 0.39, 0.70, 1.00};
  Double_t redelectronGreen[redelectronNRGBs] = { 0.17, 0.05, 0.20, 1.00};
  Double_t redelectronBlue[redelectronNRGBs] =  { 0.32, 0.33, 0.30, 0.20};
  // const Int_t redelectronNRGBs = 3;
  // const Int_t redelectronNCont = 128;
  // Double_t redelectronStops[redelectronNRGBs] = { 0.50, 0.60, 1.00};
  // Double_t redelectronRed[redelectronNRGBs] =   { 0.39, 0.70, 1.00};
  // Double_t redelectronGreen[redelectronNRGBs] = { 0.05, 0.20, 1.00};
  // Double_t redelectronBlue[redelectronNRGBs] =  { 0.33, 0.30, 0.20};

  PPalette * redelectronPalette = (PPalette*) gROOT->FindObject("redelectron");
  redelectronPalette->CreateGradientColorTable(redelectronNRGBs, redelectronStops, 
					       redelectronRed, redelectronGreen, redelectronBlue, redelectronNCont);
  

  const Int_t hotNRGBs = 3;
  const Int_t hotNCont = 128;
  Double_t hotStops[hotNRGBs] =  { 0.25, 0.5, 1.00 };
  Double_t hotRed[hotNRGBs] =   { 1.00, 1.000, 1.000 };
  Double_t hotGreen[hotNRGBs] = { 0.15, 0.984, 1.000 };
  Double_t hotBlue[hotNRGBs] =  { 0.00, 0.000, 1.000 };
  // Double_t hotRed[hotNRGBs] =   { 0.39, 0.70, 1.00 };
  // Double_t hotGreen[hotNRGBs] = { 0.05, 0.20, 1.00 };
  // Double_t hotBlue[hotNRGBs] =  { 0.33, 0.30, 0.20 };
  PPalette * hotPalette = (PPalette*) gROOT->FindObject("hot");
  hotPalette->CreateGradientColorTable(hotNRGBs, hotStops,hotRed, hotGreen, hotBlue, hotNCont);

  
    
  // Change the range of z axis for the fields to be symmetric.
  Float_t *Emax = new Float_t[Nfields];
  Float_t *Emin = new Float_t[Nfields];
  for(Int_t i=0;i<Nfields;i++) {
    Emax[i] = hE1D[i]->GetMaximum();
    Emin[i] = hE1D[i]->GetMinimum();
    if(Emax[i] > TMath::Abs(Emin[i]))
      Emin[i] = -Emax[i];
    else
      Emax[i] = -Emin[i];
    hE1D[i]->GetZaxis()->SetRangeUser(Emin[i],Emax[i]); 
  }
  
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
  Double_t ylow  = hDen2D[1]->GetYaxis()->GetBinLowEdge(FirstyBin);
  Double_t yup = hDen2D[1]->GetYaxis()->GetBinUpEdge(LastyBin);

  zStartPlasma -= shiftz; 
  zStartNeutral -= shiftz; 
  zEndNeutral -= shiftz; 
  
  if(opt.Contains("units")) {
    zStartPlasma *= skindepth / PUnits::um;
    zStartNeutral *= skindepth / PUnits::um;
    zEndNeutral *= skindepth / PUnits::um;
  }

  // Plotting
  // -----------------------------------------------

  // Canvas setup
  TCanvas *C;

  if(opt.Contains("pdf"))
    C = new TCanvas("C","",1024,576);
  else
    C = new TCanvas("C","",1680,945);
 
  // Palettes setup
  TExec *exPlasma = new TExec("exPlasma","plasmaPalette->cd();");
  TExec *exDriver   = new TExec("exDriver","redelectronPalette->cd();");
  TExec *exWitness    = new TExec("exWitness","hotPalette->cd();");
     
  // Actual Plotting!
  // ------------------------------------------------------------

  // Output file
  TString fOutName = Form("./%s/Plots/RakeFancy2D/RakeFancy2D",pData->GetPath().c_str());
  fOutName += Form("-%s_%i",pData->GetName(),time);

  // Setup Pad layout:
  const Int_t NPad = 1;
  TPad *pad[NPad];
  TH2F *hFrame[NPad];

  // Text objects
  TPaveText **textLabel = new TPaveText*[NPad];
  Float_t lMargin = 0.0;
  Float_t rMargin = 0.0;
  Float_t bMargin = 0.0;
  Float_t tMargin = 0.0;
  Float_t factor = 1.0;    
  PlasmaGlob::CanvasAsymPartition(C,NPad,lMargin,rMargin,bMargin,tMargin,factor);
  
  // Re-range:
  for(Int_t i=0;i<Nspecies;i++) {
    if(!hDen2D[i]) continue;
    hDen2D[i]->GetYaxis()->SetRangeUser(yMin -(factor-1) * yRange, yMax);
  }


  // Define the frames for plotting
  Int_t fonttype = 43;
  Int_t fontsize = 24;
  Int_t tfontsize = 28;
  Float_t txoffset = 2.0;
  Float_t lxoffset = 0.02;
  Float_t tyoffset = 1.3;
  Float_t lyoffset = 0.01;
  Float_t tylength = 0.0;
  Float_t txlength = 0.0;
  for(Int_t i=0;i<NPad;i++) {
    char name[16];
    sprintf(name,"pad_%i",i);
    pad[i] = (TPad*) gROOT->FindObject(name);
    pad[i]->SetFrameLineWidth(2);  
    // pad[i]->SetTickx(1);
    // pad[i]->SetTicky(1);

    sprintf(name,"hFrame_%i",i);
    hFrame[i] = (TH2F*) gROOT->FindObject(name);
    if(hFrame[i]) delete hFrame[i];
    if(i==NPad-1)
      hFrame[i] = (TH2F*) hDen2D[0]->Clone(name);
    else
      hFrame[i] = (TH2F*) hE2D[0]->Clone(name);
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
    hFrame[i]->GetXaxis()->SetTitleSize(tfontsize+2);
    hFrame[i]->GetXaxis()->SetTitleOffset(txoffset);
    hFrame[i]->GetXaxis()->SetLabelFont(fonttype);
    hFrame[i]->GetXaxis()->SetLabelSize(fontsize+2);
    hFrame[i]->GetXaxis()->SetLabelOffset(lxoffset);
    
    hFrame[i]->GetXaxis()->SetTickLength(yFactor*txlength/xFactor);      
  }


  C->cd(0);
  pad[0]->Draw();
  pad[0]->cd(); // <---------------------------------------------- Top Plot ---------
  
  Int_t backcolor = TColor::GetColor(10,10,10);
  pad[0]->SetFillColor(backcolor);

  if(opt.Contains("logz")) {
    pad[0]->SetLogz(1);
  } else {
    pad[0]->SetLogz(0);
  }

  hFrame[0]->Draw("AH");
    
  // Plasma
  hDen2D[0]->GetZaxis()->SetTitleFont(fonttype);
  exPlasma->Draw();
  hDen2D[0]->Draw("col same");
  
  // Beam driver.
  if(hDen2D[1]) {
    //    hDen2D[1]->GetZaxis()->SetNdivisions(505);
    exDriver->Draw();
    hDen2D[1]->Draw("col same");
  }

  // Injected electrons if any
  if(Nspecies>=3) {
    if(hDen2D[2]) {
      exWitness->Draw();
      hDen2D[2]->Draw("col same");
    }
  }

  // Fit the E1D in the pad:
  // Float_t y1 = yMin + yRange/10.;
  // Float_t y2 = y1 + yRange/2. - 2.*yRange/10.;
  Float_t y1 = yMin + yRange/80.;
  Float_t y2 = y1 + yRange/2.8 - 2*yRange/80.;
  Float_t slope = (y2-y1)/(Emax[0]-Emin[0]);
  
  Int_t lineColor = TColor::GetColor(196,30,78);
  //Int_t lineColor = TColor::GetColor(236,78,35);

  TLine *lineEzero = new TLine(xMin,(0.0-Emin[0])*slope + y1,xMax,(0.0-Emin[0])*slope + y1);
  lineEzero->SetLineColor(lineColor);
  lineEzero->SetLineStyle(2);
  if(opt.Contains("pdf"))
    lineEzero->SetLineWidth(1);
  else
    lineEzero->SetLineWidth(2);
  lineEzero->Draw();
  

  for(Int_t j=0;j<hE1D[0]->GetNbinsX();j++) {
    hE1D[0]->SetBinContent(j+1,(hE1D[0]->GetBinContent(j+1)-Emin[0])*slope + y1);
  }
  hE1D[0]->SetLineStyle(1);
  if(opt.Contains("pdf"))
    hE1D[0]->SetLineWidth(3);
  else
    hE1D[0]->SetLineWidth(4);
    
  hE1D[0]->SetLineColor(lineColor);

  hE1D[0]->Draw("sameC");



  pad[0]->RedrawAxis(); 
 
  C->cd(0);
  C->cd();

  // Print to a file
  PlasmaGlob::imgconv(C,fOutName,opt);
  // ---------------------------------------------------------

  PlasmaGlob::DestroyCanvases();
}
  
