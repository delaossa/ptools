#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>

#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TCut.h>
#include <TChain.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TPaletteAxis.h>
#include <TExec.h>
#include <TVirtualFFT.h>


#include "PData.hh"
#include "PDataHiP.hh"
#include "PGlobals.hh"
#include "PFunctions.hh"
//#include "PPalette.hh"
#include "H5Cpp.h"

using namespace std;
using namespace H5;


int main(int argc,char *argv[]) {
  if(argc<=2) {
    printf("\n Usage: %s <simulation name> <-t(time)>\n",argv[0]);
    printf("      <-i(initial time)> <-f(final time)> <-s(time step)>\n");
    printf("      <--center> <--comov> <--units>\n");
    printf("      <-z(zoom factor)> <-non(# on-axis)> <-nof(# off-axis)> <-avg(# bins)>\n");
    printf("      <-msk(mask)> <-opt(plot options)>\n");
    printf("      \n");
    return 0;
  }

  // General options
  TString   sim = "";
  Int_t    time = 0;
  Int_t  iStart = -1;
  Int_t    iEnd = -1;
  Int_t   iStep = 1;
  TString   opt = "";
  
  // Options for Density 2D
  Float_t  zoom = 1;

  // Options for Density 1D and 2D
  Int_t   NonBin = 1;
  Int_t   NofBin = 10;

  // Number of bins for averaging in z direction
  Int_t   Navg = 1;
  
  // Options for plotting mode
  UInt_t    mask = 0x3;  // First two bins on:  000011

  // Interfacing command line:
  for(int l=1;l<argc;l++){
    TString arg = argv[l];

    if(l==1) {
      sim = arg;
      continue;
    }
    
    if(arg.Contains("--pdf")) {
      opt += "pdf";
    } else if(arg.Contains("--eps")) {
      opt += "eps";
    } else if(arg.Contains("--png")) {
      opt += "png";
    } else if(arg.Contains("-msk")) {
      char ss[4];
      sscanf(arg,"%4s%i",ss,&mask);
    } else if(arg.Contains("--comov")){
      opt += "comov";
    } else if(arg.Contains("--units")){
      opt += "units";
    } else if(arg.Contains("--center")){
      opt += "center";
    } else if(arg.Contains("--curr")){
      opt += "curr";	
    } else if(arg.Contains("--alt1D")){
      opt += "alt1D";	
    } else if(arg.Contains("-t")) {
      char ss[2];
      sscanf(arg,"%2s%i",ss,&time);
    } else if(arg.Contains("-i")) {
      char ss[2];
      sscanf(arg,"%2s%i",ss,&iStart);
    } else if(arg.Contains("-f")) {
      char ss[2];
      sscanf(arg,"%2s%i",ss,&iEnd);
    } else if(arg.Contains("-s")) {
      char ss[2];
      sscanf(arg,"%2s%i",ss,&iStep);
    } else if(arg.Contains("-z")) {
      char ss[2];
      sscanf(arg,"%2s%f",ss,&zoom);
    } else if(arg.Contains("-non")) {
      char ss[2];
      sscanf(arg,"%4s%i",ss,&NonBin);
    } else if(arg.Contains("-nof")) {
      char ss[2];
      sscanf(arg,"%4s%i",ss,&NofBin);
    } else if(arg.Contains("-avg")){
      char ss[4];
      sscanf(arg,"%4s%i",ss,&Navg);
    } else if(arg.Contains("-opt")){
      char ss[4], sopt[64];
      sscanf(arg,"%4s%s",ss,sopt);
      opt += sopt;
    } else {
      cout << Form("\t Invalid argument (%i): exiting...\n",l) << endl;
      return 0;
    }
  }
  
  
  PGlobals::Initialize();
  
  // Load PData
  PData *pData = PData::Get(sim.Data());
  if(pData->isHiPACE()) {
    delete pData; pData = NULL;
    pData = PDataHiP::Get(sim.Data());
    opt += "comov";
  }
  
  pData->SetNavg(Navg);

  if(iStart<0) iStart = time;
  if(iEnd<=iStart) iEnd = iStart;
  
  // Some plasma constants
  Double_t n0 = pData->GetPlasmaDensity();
  //  Double_t omegap = pData->GetPlasmaFrequency();
  //  Double_t timedepth = pData->GetPlasmaTimeDepth();
  Double_t period = pData->GetPlasmaPeriod();
  Double_t kp = pData->GetPlasmaK();
  Double_t skindepth = pData->GetPlasmaSkinDepth();
  Double_t wavelength = pData->GetPlasmaWaveLength();
  Double_t E0 = pData->GetPlasmaE0();

  // Some beam properties:
  // Double_t Ebeam = pData->GetBeamEnergy();
  // Double_t gamma = pData->GetBeamGamma();
  // Double_t vbeam = pData->GetBeamVelocity();
  
  // Other parameters
  // Double_t trapPotential = 1.0;
  // Double_t trapPotential = 1.0 - (1.0/gamma);

  // Time looper
  for(Int_t i=iStart; i<iEnd+1; i+=iStep) {

    time = i;

    cout << Form("\n Processing time step %i",time) << endl;
    cout << Form(" -------------------------\n") << endl;

    pData->LoadFileNames(time);
    
    if(!pData->IsInit()) continue;
    if(time==iStart) pData->PrintData();
    
    // Time in OU
    Double_t Time = pData->GetRealTime();
    cout << Form(" Real time = %.2f  ",Time);
    Time += pData->ShiftT(opt);
    cout << Form(" Shifted time = %.2f  ",Time) << endl;
    
    // Off-axis definition    
    Double_t rms0 = pData->GetBeamRmsX() * kp;
    if(pData->IsCyl()) rms0 = pData->GetBeamRmsR() * kp;
    
    // Calculate the "axis range" in number of bins. 
    if(NonBin==0) {  // If NonBin==0 a RMS width is taken.
      
      if(rms0>0.0)
	NonBin =  TMath::Nint(rms0 / pData->GetDX(1));
      else
	NonBin = 1;
    
    } else if(NonBin<0) { // If negative, the number is taken in units of 10/kp
      NonBin = -TMath::Nint(float(NonBin*0.1) / pData->GetDX(1));
    }
  
    // Off-axis bin
    if(NofBin==0) {  // If NofBin==0 a RMS width is taken.
      if(pData->IsCyl()) rms0 = pData->GetBeamRmsR() * kp;
      
      if(rms0>0.0)
	NofBin =  TMath::Nint(rms0 / pData->GetDX(1));
      else
	NofBin = 1;
    
    } else if(NofBin<0) { // If negative, the number is taken in units of (kp^-1)/10
      NofBin = -TMath::Nint(float(NofBin*0.1) / pData->GetDX(1));
    }

    Double_t xon  = NonBin *  pData->GetDX(1); // onaxis interval in norm. units
    Double_t xoff = NofBin *  pData->GetDX(1); // offaxis distance in norm. units

    // Range for current calculation
    Int_t NonBinT = pData->GetNX(1)/4;
    if(rms0>0) {    // This adapts the size of the data chunk to the initial rms size of the beam
      NonBinT = TMath::Nint(3.5 * rms0 / pData->GetDX(1));
      // cout << Form(" Getting current within an on-axis range of NonBinT = %i",NonBinT) << endl;
    }
    Double_t xint = NonBinT * pData->GetDX(1);

    
    // Slice width limits.
    Int_t FirstxBin = 0;
    Int_t LastxBin = 0;
    if(!pData->IsCyl()) {
      FirstxBin = pData->GetNX(1)/2 + 1 - NonBin;
      LastxBin =  pData->GetNX(1)/2 + NonBin;
    } else {
      FirstxBin = 1; 
      LastxBin  = NonBin;
    }
    Int_t FirstOffxBin = FirstxBin + NofBin;
    Int_t LastOffxBin = LastxBin + NofBin;

    // Zoom window:  
    Double_t xRange = (pData->GetXMax(1) - pData->GetXMin(1))/zoom;
    Double_t xMid   = (pData->GetXMax(1) + pData->GetXMin(1))/2.;
    Double_t xMin = xMid - xRange/2.0;
    Double_t xMax = xMid + xRange/2.0;
    if(pData->IsCyl()) {
      xMin = pData->GetXMin(1);
      xMax = xRange;
    }
    pData->SetX2Min(xMin);
    pData->SetX2Max(xMax);
  
    Double_t zMin = pData->GetX1Min();
    Double_t zMax = pData->GetX1Max();
    // Double_t zRange = zMax - zMin;  
    // pData->SetX1Min(zMin);
    // pData->SetX1Max(zMax);

    cout << Form(" --- doSnapshot ---\n") << endl;
  
    cout << Form(" Plotting range:  %.2f < x1 < %.2f ,  %.2f < x2 < %.2f",zMin,zMax,xMin,xMax) << endl;

    // ----------------------------------------------------------------------------------
  

    cout << Form("\n Reading simulation data: ") << endl; 

    Float_t maxCur = -999.0;
    
    // Get charge density histos
    Int_t Nspecies = pData->NSpecies();
    TH2F **hDen2D = new TH2F*[Nspecies];
    // Get charge density on-axis
    TH1F **hDen1D = new TH1F*[Nspecies];
    // And electric current (integrated)
    TH1F **hCur1D = new TH1F*[Nspecies];

    for(Int_t i=0;i<Nspecies;i++) {
     
      hDen2D[i] = NULL;
      hDen1D[i] = NULL;
      hCur1D[i] = NULL;
      
      if(!pData->GetChargeFileName(i)) 
	continue;

      cout << Form(" -> Getting charge density of specie %i  (%s)", i, pData->GetSpeciesName(i).c_str())  << endl;
      
    
      char hName[24];
      sprintf(hName,"hDen2D_%i",i);
      hDen2D[i] = (TH2F*) gROOT->FindObject(hName);
      if(hDen2D[i]) delete hDen2D[i];

      if(pData->Is3D()) {
	hDen2D[i] = pData->GetCharge2DSliceZX(i,-1,NonBin,opt+"avg");
	//hDen2D[i] = pData->GetH2SliceZX(pData->GetChargeFileName(i)->c_str(),"charge",-1,NonBin,opt+"avg");
      } else {
	hDen2D[i] = pData->GetCharge(i,opt);
      }

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
	hDen2D[i]->GetYaxis()->SetTitle("k_{p} x");
    
      if(i==0)
	hDen2D[i]->GetZaxis()->SetTitle("n_{p}/n_{0}");
      else 
	hDen2D[i]->GetZaxis()->SetTitle(Form("n_{%i}/n_{0}",i));
      
      sprintf(hName,"hDen1D_%i",i);
      hDen1D[i] = (TH1F*) gROOT->FindObject(hName);
      if(hDen1D[i]) delete hDen1D[i];
      
      // 1D histograms
      if(pData->Is3D()) {
	hDen1D[i] = pData->GetH1SliceZ3D(pData->GetChargeFileName(i)->c_str(),"charge",-1,NonBin,-1,NonBin,opt+"avg");
      } else if(pData->IsCyl()) { // Cylindrical: The first bin with r>0 is actually the number 1 (not the 0).
	hDen1D[i] = pData->GetH1SliceZ(pData->GetChargeFileName(i)->c_str(),"charge",1,NonBin,opt+"avg");
      } else { // 2D cartesian
	hDen1D[i] = pData->GetH1SliceZ(pData->GetChargeFileName(i)->c_str(),"charge",-1,NonBin,opt+"avg");
      }
      hDen1D[i]->SetName(hName); 
   
      // if(hDen1D[i]) delete hDen1D[i];
      // hDen1D[i] = (TH1F*) hE2D[i]->ProjectionX(hName,FirstxBin,LastxBin);
      // hDen1D[i]->Scale(1.0/(LastxBin-FirstxBin+1));
    
      if(opt.Contains("comov"))
	hDen1D[i]->GetXaxis()->SetTitle("k_{p} #zeta");
      else
	hDen1D[i]->GetXaxis()->SetTitle("k_{p} z");

      if(i==0)
	hDen1D[i]->GetYaxis()->SetTitle("n/n_{0}");
      else if(i==1)
	hDen1D[i]->GetYaxis()->SetTitle("n_{b}/n_{0}");
      else   
	hDen1D[i]->GetYaxis()->SetTitle("n_{i}/n_{0}");

      // Get the current:
      if(i==0 && pData->GetSpeciesName(i).find("beam") == string::npos) continue;
      
      sprintf(hName,"hCur1D_%i",i);
      hCur1D[i] = (TH1F*) gROOT->FindObject(hName);
      if(hCur1D[i]) delete hCur1D[i];
      
      if(opt.Contains("curr")) {

	// To get the current is needed to read in a wider transverse range which includes all the charge.	
	if(pData->Is3D()) {
	  hCur1D[i] = pData->GetH1SliceZ3D(pData->GetChargeFileName(i)->c_str(),"charge",-1,NonBinT,-1,NonBinT,opt+"int");
	} else if(pData->IsCyl()) { // Cylindrical: The first bin with r>0 is actually the number 1 (not the 0).
	  hCur1D[i] = pData->GetH1SliceZ(pData->GetChargeFileName(i)->c_str(),"charge",1,NonBinT,opt+"int");
	} else { // 2D cartesian
	  hCur1D[i] = pData->GetH1SliceZ(pData->GetChargeFileName(i)->c_str(),"charge",-1,NonBinT,opt+"int");
	}
	hCur1D[i]->SetName(hName); 
      
	// Normalized current:
	Double_t dV = skindepth * skindepth * skindepth;
	if(!pData->Is3D()) {
	  hCur1D[i]->Scale(TMath::Abs( (5./3.) * n0 * dV * PConst::ElectronCharge * kp * PConst::c_light) / PConst::I0);
	} else {
	  hCur1D[i]->Scale(TMath::Abs(n0 * dV * PConst::ElectronCharge * kp * PConst::c_light) / PConst::I0);
	}

	if(hCur1D[i]->GetMaximum()>maxCur) {
	  maxCur = hCur1D[i]->GetMaximum();
	}
	
	hCur1D[i]->GetYaxis()->SetTitle("#Lambda_{b}");  
	if(opt.Contains("comov")) {
	  hCur1D[i]->GetXaxis()->SetTitle("k_{p} #zeta");
	} else {
	  hCur1D[i]->GetXaxis()->SetTitle("k_{p} z");
	}
	
      }

    }

    // Current density
    TH2F **hJz2D = new TH2F*[Nspecies];
    for(Int_t i=0;i<Nspecies;i++) {
      
      hJz2D[i] = NULL;
      if(!pData->GetCurrentFileName(i)) 
	continue;
      
      cout << Form(" -> Getting current density of specie %i  (%s)", i, pData->GetSpeciesName(i).c_str())  << endl;
      
      char hName[24];
      sprintf(hName,"hJz2D_%i",i);
      hJz2D[i] = (TH2F*) gROOT->FindObject(hName);
      if(hJz2D[i]) delete hJz2D[i];

      if(pData->Is3D()) {
	//hJz2D[i] = pData->GetCharge2DSliceZX(i,-1,NonBin,opt+"avg");
	hJz2D[i] = pData->GetH2SliceZX(pData->GetCurrentFileName(i)->c_str(),"j1",-1,NonBin,opt+"avg");
      } else {
	hJz2D[i] = pData->GetCharge(i,opt);
      }

      hJz2D[i]->SetName(hName);
      hJz2D[i]->GetXaxis()->CenterTitle();
      hJz2D[i]->GetYaxis()->CenterTitle();
      hJz2D[i]->GetZaxis()->CenterTitle();
    
      if(opt.Contains("comov"))
	hJz2D[i]->GetXaxis()->SetTitle("k_{p} #zeta");
      else
	hJz2D[i]->GetXaxis()->SetTitle("k_{p} z");
    
      if(pData->IsCyl()) 
	hJz2D[i]->GetYaxis()->SetTitle("k_{p} r");
      else
	hJz2D[i]->GetYaxis()->SetTitle("k_{p} x");
    
      hJz2D[i]->GetZaxis()->SetTitle(Form("j_{z,%i}",i));
      
      
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

      cout << Form(" -> Getting electric field number ") << i+1 << endl;
    
      char hName[24];
      sprintf(hName,"hE2D_%i",i);
      hE2D[i] = (TH2F*) gROOT->FindObject(hName);
      if(hE2D[i]) delete hE2D[i];
  
      
      char enam[3]; sprintf(enam,"e%i",i+1);
      if(!pData->Is3D())
	hE2D[i] = pData->GetEField(i,opt);
      else
	hE2D[i] = pData->GetEField2DSliceZX(i,-1,NonBin,opt+"avg");
      //hE2D[i] = pData->GetH2SliceZX(pData->GetEfieldFileName(i)->c_str(),enam,-1,NonBin,opt+"avg");

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
	hE2D[i]->GetYaxis()->SetTitle("k_{p} x");
    
      if(i==0)
	hE2D[i]->GetZaxis()->SetTitle("E_{z}/E_{0}");
      else if(i==1)
	hE2D[i]->GetZaxis()->SetTitle("E_{x}/E_{0}");
      else if(i==2)
	hE2D[i]->GetZaxis()->SetTitle("E_{y}/E_{0}");
    
      sprintf(hName,"hE1D_%i",i);
      hE1D[i] = (TH1F*) gROOT->FindObject(hName);
      if(hE1D[i]) delete hE1D[i];
      
      // 1D histograms
      char nam[3];
      if(pData->isHiPACE()){
	if(i==0)
	  sprintf(nam,"Ez");
	else if(i==1)
	  sprintf(nam,"Ex");
	else if(i==2)
	  sprintf(nam,"Ey");
      } else
	sprintf(nam,"e%i",i+1);
      
      if(pData->Is3D()) {
      
	if(!opt.Contains("alt1D")){
	  if(i==0) 
	    hE1D[i] = pData->GetH1SliceZ3D(pData->GetEfieldFileName(i)->c_str(),nam,-1,NonBin,-1,NonBin,opt+"avg");
	  else // In case of transverse fields, the 1D line is taken off-axis 
	    hE1D[i] = pData->GetH1SliceZ3D(pData->GetEfieldFileName(i)->c_str(),nam,-NofBin,NonBin,-1,NonBin,opt+"avg");
	
	} else { // Alternative
	  cout << " Alternative! " << endl;
	  if(i==0) {
	    hE1D[i] = (TH1F*) hE2D[i]->ProjectionX(hName,FirstxBin,LastxBin);
	    hE1D[i]->Scale(1.0/(LastxBin-FirstxBin+1));
	  } else { // In case of transverse fields, the 1D line is taken off-axis 
	    hE1D[i] = (TH1F*) hE2D[i]->ProjectionX(hName,FirstOffxBin,LastOffxBin);
	    hE1D[i]->Scale(1.0/(LastOffxBin-FirstOffxBin+1));
	  }
	}
	
      } else if(pData->IsCyl()) { // Cylindrical: The first bin with r>0 is actually the number 1 (not the 0).
      
	hE1D[i] = pData->GetH1SliceZ(pData->GetEfieldFileName(i)->c_str(),nam,1,NonBin,opt+"avg");
      
      } else { // 2D cartesian
      
	if(i==0) 
	  hE1D[i] = pData->GetH1SliceZ(pData->GetEfieldFileName(i)->c_str(),nam,-1,NonBin,opt+"avg");
	else 
	  hE1D[i] = pData->GetH1SliceZ(pData->GetEfieldFileName(i)->c_str(),nam,-NofBin,NonBin,opt+"avg");    
      }
    
      hE1D[i]->SetName(hName);
      if(opt.Contains("comov"))
	hE1D[i]->GetXaxis()->SetTitle("#zeta [c/#omega_{p}]");
      else
	hE1D[i]->GetXaxis()->SetTitle("z [c/#omega_{p}]");
    
      if(i==0)
	hE1D[i]->GetYaxis()->SetTitle("E_{z}/E_{0}");
      else if(i==1)
	hE1D[i]->GetYaxis()->SetTitle("E_{x}/E_{0}");
      else if(i==2)
	hE1D[i]->GetYaxis()->SetTitle("E_{y}/E_{0}");

      // Alternative
      // if(hE1D[i]) delete hE1D[i];
      // hE1D[i] = (TH1F*) hE2D[i]->ProjectionX(hName,FirstxBin,LastxBin);
      // hE1D[i]->Scale(1.0/(LastxBin-FirstxBin+1));
     
    }

    // Get magnetic fields 2D
    const Int_t NBfields = 3;
    TH2F **hB2D = new TH2F*[NBfields];
    TH1F **hB1D = new TH1F*[NBfields];
    for(Int_t i=0;i<NBfields;i++) {
      hB2D[i] = NULL;
      hB1D[i] = NULL;
      if(i<2) continue;  // Just get the third component.

      if(!pData->GetBfieldFileName(i)) 
	continue;

      cout << Form(" -> Getting magnetic field number ") << i+1 << endl;
    
      char hName[24];
      sprintf(hName,"hB2D_%i",i);
      hB2D[i] = (TH2F*) gROOT->FindObject(hName);
      if(hB2D[i]) delete hB2D[i];
    
      char bnam[3]; sprintf(bnam,"b%i",i+1);
      if(!pData->Is3D())
	hB2D[i] = pData->GetBField(i,opt);
      else
	hB2D[i] = pData->GetBField2DSliceZX(i,-1,NonBin,opt+"avg");
      //hB2D[i] = pData->GetH2SliceZX(pData->GetBfieldFileName(i)->c_str(),bnam,-1,NonBin,opt+"avg");

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
	hB2D[i]->GetYaxis()->SetTitle("k_{p} x");
    
      if(i==0)
	hB2D[i]->GetZaxis()->SetTitle("B_{z}/E_{0}");
      else if(i==1)
	hB2D[i]->GetZaxis()->SetTitle("B_{x}/E_{0}");
      else if(i==2)
	hB2D[i]->GetZaxis()->SetTitle("B_{y}/E_{0}");


      
      // 1D histograms
      char nam[3];
      if(pData->isHiPACE()){
	if(i==0)
	  sprintf(nam,"Bz");
	else if(i==1)
	  sprintf(nam,"Bx");
	else if(i==2)
	  sprintf(nam,"By");
      } else
	sprintf(nam,"b%i",i+1);
      
      if(pData->Is3D()) {
	
	if(!opt.Contains("alt1D")){
	  if(i==0) 
	    hB1D[i] = pData->GetH1SliceZ3D(pData->GetBfieldFileName(i)->c_str(),nam,-1,NonBin,-1,NonBin,opt+"avg");
	  else // In case of transverse fields, the 1D line is taken offaxis
	    hB1D[i] = pData->GetH1SliceZ3D(pData->GetBfieldFileName(i)->c_str(),nam,-NofBin,NonBin,-1,NonBin,opt+"avg");
	  
	} else {  // Alternative
	  cout << " Alternative! " << endl;
	  if(i==0) {
	    hB1D[i] = (TH1F*) hB2D[i]->ProjectionX(hName,FirstxBin,LastxBin);
	    hB1D[i]->Scale(1.0/(LastxBin-FirstxBin+1));
	  } else { // In case of transverse fields, the 1D line is taken off-axis 
	    hB1D[i] = (TH1F*) hB2D[i]->ProjectionX(hName,FirstOffxBin,LastOffxBin);
	    hB1D[i]->Scale(1.0/(LastOffxBin-FirstOffxBin+1));
	  }
	}


      } else if(pData->IsCyl()) { // Cylindrical: The first bin with r>0 is actually the number 1 (not the 0).
      
	hB1D[i] = pData->GetH1SliceZ(pData->GetBfieldFileName(i)->c_str(),nam,1,NonBin,opt+"avg");
      
      } else { // 2D cartesian
      
	if(i==0) 
	  hB1D[i] = pData->GetH1SliceZ(pData->GetBfieldFileName(i)->c_str(),nam,-1,NonBin,opt+"avg");
	else 
	  hB1D[i] = pData->GetH1SliceZ(pData->GetBfieldFileName(i)->c_str(),nam,-NonBin,NonBin,opt+"avg");    
      }
    
      hB1D[i]->SetName(hName);
      if(opt.Contains("comov"))
	hB1D[i]->GetXaxis()->SetTitle("#zeta [c/#omega_{p}]");
      else
	hB1D[i]->GetXaxis()->SetTitle("z [c/#omega_{p}]");
    
      if(i==0)
	hB1D[i]->GetYaxis()->SetTitle("B_{z}/E_{0}");
      else if(i==1)
	hB1D[i]->GetYaxis()->SetTitle("B_{x}/E_{0}");
      else if(i==2)
	hB1D[i]->GetYaxis()->SetTitle("B_{y}/E_{0}");
      
      
    }


    // Get vector potential of the radiation
    // (substract the electro-estatic fields from the total)
    if(hDen2D[0]) {

      char hName[24];
      sprintf(hName,"hRho2D");
      TH2F *hRho2D = (TH2F*) hDen2D[0]->Clone(hName);

      for(Int_t i=1;i<Nspecies;i++) {
	if(hDen2D[i]) hRho2D->Add(hDen2D[i],1);
      }

      Int_t NBinsZ = hRho2D->GetXaxis()->GetNbins();
      Float_t zmin = hRho2D->GetXaxis()->GetBinLowEdge(1);
      Float_t zmax = hRho2D->GetXaxis()->GetBinUpEdge(NBinsZ);
      Float_t zrange = zmax-zmin;
      
      Int_t NBinsX = hRho2D->GetYaxis()->GetNbins();
      Float_t xmin = hRho2D->GetYaxis()->GetBinLowEdge(1);
      Float_t xmax = hRho2D->GetYaxis()->GetBinUpEdge(NBinsX);
      Float_t xrange = xmax-xmin;
      
      Float_t kzmin = 0.;
      Float_t kzmax = (TMath::TwoPi() / zrange);
      
      Float_t kxmin = 0.;
      Float_t kxmax = (TMath::TwoPi() / xrange);
      
      Int_t dims[2] = {NBinsZ,NBinsX};
      //  Double_t *data = new Double_t[NBinsZ*NBinsX];
      Double_t *data = new Double_t[NBinsZ*2*(NBinsX/2+1)];
      // Extra padding!
      // See http://www.fftw.org/fftw3_doc/Multi_002dDimensional-DFTs-of-Real-Data.html#Multi_002dDimensional-DFTs-of-Real-Data 
      
      for(Int_t i=0; i<NBinsZ; i++) {
	
	for(Int_t j=0; j<NBinsX; j++) {
	  
	  Int_t index =  i * NBinsX + j;
	  
	  data[index]  = hRho2D->GetBinContent(i+1,j+1);
	  
	  // substract ion background
	  data[index] -= 1;
	}
      }
  
      // cout << Form(" Fourier transform ..." );
      TVirtualFFT *fft = TVirtualFFT::FFT(2, dims, "R2C ES K");
      fft->SetPoints(data);
      fft->Transform();  
      fft->GetPoints(data);    

      // cout << Form(" Solving equation for potential in fourier space" )  << endl;
      // Equation: \phi(kz,kx) = \rho(kz,kx)/(kz^2 + kx^2) 
      for(Int_t i=0; i<dims[0]; i++) {
	for(Int_t j=0; j<(dims[1]/2+1); j++) {

	  Int_t index =  2 * (i * (dims[1]/2+1)  + j);
	  //      cout << Form(" i = %4i  j = %4i  index = %10i",i,j,index) << endl;
      
	  Float_t kz;
	  if(i<dims[0]/2)
	    kz = kzmax * (i);
	  else
	    kz =  kzmax * (i-dims[0]);
	  Float_t kx = (kxmax) * (j) ;
      
	  Float_t k2 = kx*kx + kz*kz;
	  if(k2==0) {
	    data[index] = 0.0;
	    data[index+1] = 0.0; 
	  } else {
	    data[index] /= k2;
	    data[index+1] /= k2;      
	  }
	}
      }

      //  cout << Form(" Inverse Fourier transform ..." );
      // backward transform:
      TVirtualFFT *fft_back = TVirtualFFT::FFT(2, dims, "C2R ES K");
      fft_back->SetPoints(data);
      fft_back->Transform();  
      
      Double_t *re_back = fft_back->GetPointsReal();
  
      //  cout << Form(" Take the gradient of the potential ..." );
      TH2F *hE2D_1_ifft = new TH2F("hE2D_1_ifft","",NBinsZ,zmin,zmax,NBinsX,xmin,xmax);
      Double_t dx = hE2D_1_ifft->GetYaxis()->GetBinWidth(1);
      for(Int_t i=0; i<NBinsZ; i++) {
	for(Int_t j=0; j<NBinsX; j++) {
	  Int_t index =  i * NBinsX + j;
	  
	  Float_t der = 0.;
	  if(j>2 && j<NBinsX-2) {
	    Int_t indjp1 =  i * NBinsX + (j+1);
	    Int_t indjp2 =  i * NBinsX + (j+2);
	    Int_t indjm1 =  i * NBinsX + (j-1);
	    Int_t indjm2 =  i * NBinsX + (j-2);
	    
	    der =  ( 4.0 / 3.0 * (re_back[indjp1] - re_back[indjm1]) / (2.0 * dx)
		     - 1.0 / 3.0 * (re_back[indjp2] - re_back[indjm2]) / (4.0 * dx) );
	    
	  }
	  
	  hE2D_1_ifft->SetBinContent(i+1,j+1,der);
	}
      }
      
      hE2D_1_ifft->Scale(1.0/(NBinsZ*NBinsX));
      TH1D *hE1D_1_ifft = hE2D_1_ifft->ProjectionX("hE1D_1_ifft",hE2D_1_ifft->GetNbinsY()/2,hE2D_1_ifft->GetNbinsY()/2);
      
    }
    
    // Chaning to user units: 
    // --------------------------
  
    Double_t denUnit, spaUnit, propUnit, timUnit, eUnit, curUnit;
    string denSUnit, spaSUnit, propSUnit, timSUnit, eSUnit, curSUnit;
    
    if(opt.Contains("units") && n0) {

      cout << endl << Form("Changing to SI units:") << endl;
      
      // Get the best units for each quantity
      PUnits::BestUnit bdenSUnit(n0,"PartDensity");
      bdenSUnit.GetBestUnits(denUnit,denSUnit);
      cout << Form(" n0 = %.2f %s", n0/denUnit, denSUnit.c_str()) << endl;
      // cout << bdenSUnit << endl;
  
      PUnits::BestUnit bspaSUnit(wavelength,"Length");
      bspaSUnit.GetBestUnits(spaUnit,spaSUnit);
      cout << Form(" L  = %.2f %s", wavelength/spaUnit, spaSUnit.c_str()) << endl;

      PUnits::BestUnit btimSUnit(period,"Time");
      btimSUnit.GetBestUnits(timUnit,timSUnit);
      cout << Form(" T  = %.2f %s", period/timUnit, timSUnit.c_str()) << endl;

      propUnit = PUnits::mm;
      propSUnit = "mm";

      
      Float_t Emax =  hE1D[0]->GetMaximum();
      if(Emax < 10E-3) Emax = 10E-3;
      PUnits::BestUnit beSUnit(E0 * Emax,"Efield");
      beSUnit.GetBestUnits(eUnit,eSUnit);
      cout << Form(" E0 = %.2f %s", E0 * Emax/eUnit, eSUnit.c_str()) << endl;

      for(Int_t i=0;i<Nspecies;i++) {

	if(!hDen2D[i]) continue;
    
	Int_t NbinsX = hDen2D[i]->GetNbinsX();
	Double_t zMin = skindepth * hDen2D[i]->GetXaxis()->GetXmin() / spaUnit;
	Double_t zMax = skindepth * hDen2D[i]->GetXaxis()->GetXmax() / spaUnit;
	Int_t NbinsY = hDen2D[i]->GetNbinsY();
	Double_t ymin = skindepth * hDen2D[i]->GetYaxis()->GetXmin() / spaUnit;
	Double_t ymax = skindepth * hDen2D[i]->GetYaxis()->GetXmax() / spaUnit;
	hDen2D[i]->SetBins(NbinsX,zMin,zMax,NbinsY,ymin,ymax);
	// for(Int_t j=0;j<hDen2D[i]->GetNbinsX();j++) {
	// 	for(Int_t k=0;k<hDen2D[i]->GetNbinsY();k++) {
	// 	  hDen2D[i]->SetBinContent(j,k, hDen2D[i]->GetBinContent(j,k) * n0 / (1e17/PUnits::cm3) );
	// 	}
	// }
	if(pData->IsCyl())
	  hDen2D[i]->GetYaxis()->SetTitle(Form("r [%s]",spaSUnit.c_str()));      
	else
	  hDen2D[i]->GetYaxis()->SetTitle(Form("x [%s]",spaSUnit.c_str()));      
      
	if(opt.Contains("comov"))
	  hDen2D[i]->GetXaxis()->SetTitle(Form("#zeta [%s]",spaSUnit.c_str()));
	else
	  hDen2D[i]->GetXaxis()->SetTitle(Form("z [%s]",spaSUnit.c_str()));
      
	// if(i==0)
	// 	hDen2D[i]->GetZaxis()->SetTitle("n_{e} [10^{17}/cm^{3}]");
	// else if(i==1)
	// 	hDen2D[i]->GetZaxis()->SetTitle("n_{b} [10^{17}/cm^{3}]"); 
	// else
	// 	hDen2D[i]->GetZaxis()->SetTitle("n_{i} [10^{17}/cm^{3}]"); 

	if(hJz2D[i]) {
	  hJz2D[i]->SetBins(NbinsX,zMin,zMax,NbinsY,ymin,ymax);
	  
	  if(pData->IsCyl())
	    hJz2D[i]->GetYaxis()->SetTitle(Form("r [%s]",spaSUnit.c_str()));      
	  else
	    hJz2D[i]->GetYaxis()->SetTitle(Form("x [%s]",spaSUnit.c_str()));      
	  
	  if(opt.Contains("comov"))
	    hJz2D[i]->GetXaxis()->SetTitle(Form("#zeta [%s]",spaSUnit.c_str()));
	  else
	    hJz2D[i]->GetXaxis()->SetTitle(Form("z [%s]",spaSUnit.c_str()));
	}
	
	
	hDen1D[i]->SetBins(NbinsX,zMin,zMax);
	// for(Int_t j=0;j<hDen1D[i]->GetNbinsX();j++) {
	// 	hDen1D[i]->SetBinContent(j, hDen1D[i]->GetBinContent(j) * n0 / (1e17/PUnits::cm3) );
	// }
      
      
	if(opt.Contains("comov"))
	  hDen1D[i]->GetXaxis()->SetTitle(Form("#zeta [%s]",spaSUnit.c_str()));
	else
	  hDen1D[i]->GetXaxis()->SetTitle(Form("z [%s]",spaSUnit.c_str()));
      
	if(hCur1D[i]) {

	  hCur1D[i]->SetBins(NbinsX,zMin,zMax);
	
	  Double_t binSize = ((zMax - zMin)/NbinsX) * spaUnit;
	  Double_t Charge = hCur1D[i]->Integral() * binSize * PConst::I0 / PConst::c_light;
	  cout << Form(" Integrated charge of specie %3i = %8f pC",i,Charge/PUnits::picocoulomb) << endl;
	  
	  // PUnits::BestUnit bcurSUnit(hCur1D[i]->GetMaximum() * PConst::I0,"Current");
	  // bcurSUnit.GetBestUnits(curUnit,curSUnit);
	  PUnits::BestUnit bcurSUnit(maxCur * PConst::I0,"Current");
	  bcurSUnit.GetBestUnits(curUnit,curSUnit);
	  
	  hCur1D[i]->Scale(PConst::I0 / curUnit);
	  hCur1D[i]->GetYaxis()->SetTitle(Form("I [%s]",curSUnit.c_str()));
	  // if(i==0)
	  //   hCur1D[i]->GetYaxis()->SetTitle("I_{p} [kA]");
	  // else if(i==1)
	  //   hCur1D[i]->GetYaxis()->SetTitle("I_{b} [kA]");
	  // else if(i==2) 
	  //   hCur1D[i]->GetYaxis()->SetTitle("I_{i} [kA]");
	  // else if(i==3) 
	  //   hCur1D[i]->GetYaxis()->SetTitle("I_{x} [kA]");
	  
	  if(opt.Contains("comov"))
	    hCur1D[i]->GetXaxis()->SetTitle(Form("#zeta [%s]",spaSUnit.c_str()));
	  else
	    hCur1D[i]->GetXaxis()->SetTitle(Form("z [%s]",spaSUnit.c_str()));
	
	}
      }
      
      
      for(Int_t i=0;i<Nfields;i++) {
	if(!hE2D[i]) continue;

	Int_t NbinsX = hE2D[i]->GetNbinsX();
	Double_t zMin = skindepth * hE2D[i]->GetXaxis()->GetXmin() / spaUnit;
	Double_t zMax = skindepth * hE2D[i]->GetXaxis()->GetXmax() / spaUnit;
	Int_t NbinsY = hE2D[i]->GetNbinsY();
	Double_t ymin = skindepth * hE2D[i]->GetYaxis()->GetXmin() / spaUnit;
	Double_t ymax = skindepth * hE2D[i]->GetYaxis()->GetXmax() / spaUnit;
	hE2D[i]->SetBins(NbinsX,zMin,zMax,NbinsY,ymin,ymax);
	hE1D[i]->SetBins(NbinsX,zMin,zMax);
            
	for(Int_t j=0;j<=hE2D[i]->GetNbinsX();j++) {
	  for(Int_t k=0;k<=hE2D[i]->GetNbinsY();k++) {
	    hE2D[i]->SetBinContent(j,k, hE2D[i]->GetBinContent(j,k) * E0 / eUnit );
	  }
	  hE1D[i]->SetBinContent(j, hE1D[i]->GetBinContent(j) * E0 / eUnit );
	}
      
	if(pData->IsCyl())
	  hE2D[i]->GetYaxis()->SetTitle(Form("r [%s]",spaSUnit.c_str()));
	else
	  hE2D[i]->GetYaxis()->SetTitle(Form("x [%s]",spaSUnit.c_str()));
      
	if(opt.Contains("comov"))
	  hE2D[i]->GetXaxis()->SetTitle(Form("#zeta [%s]",spaSUnit.c_str()));
	else
	  hE2D[i]->GetXaxis()->SetTitle(Form("z [%s]",spaSUnit.c_str()));
      
	if(i==0)
	  hE2D[i]->GetZaxis()->SetTitle(Form("E_{z} [%s]",eSUnit.c_str()));
	else if(i==1)
	  hE2D[i]->GetZaxis()->SetTitle(Form("E_{x} [%s]",eSUnit.c_str()));
	else if(i==2)
	  hE2D[i]->GetZaxis()->SetTitle(Form("E_{y} [%s]",eSUnit.c_str()));
      
      
	if(opt.Contains("comov"))
	  hE1D[i]->GetXaxis()->SetTitle(Form("#zeta [%s]",spaSUnit.c_str()));
	else
	  hE1D[i]->GetXaxis()->SetTitle(Form("z [%s]",spaSUnit.c_str()));
      
	if(i==0)
	  hE1D[i]->GetYaxis()->SetTitle(Form("E_{z} [%s]",eSUnit.c_str()));
	else if(i==1)
	  hE1D[i]->GetYaxis()->SetTitle(Form("E_{x} [%s]",eSUnit.c_str()));
	else if(i==2)
	  hE1D[i]->GetYaxis()->SetTitle(Form("E_{y} [%s]",eSUnit.c_str()));
	
      }

      for(Int_t i=0;i<NBfields;i++) {
	if(!hB2D[i]) continue;

	Int_t NbinsX = hB2D[i]->GetNbinsX();
	Double_t zMin = skindepth * hB2D[i]->GetXaxis()->GetXmin() / spaUnit;
	Double_t zMax = skindepth * hB2D[i]->GetXaxis()->GetXmax() / spaUnit;
	Int_t NbinsY = hB2D[i]->GetNbinsY();
	Double_t ymin = skindepth * hB2D[i]->GetYaxis()->GetXmin() / spaUnit;
	Double_t ymax = skindepth * hB2D[i]->GetYaxis()->GetXmax() / spaUnit;
	hB2D[i]->SetBins(NbinsX,zMin,zMax,NbinsY,ymin,ymax);
	hB1D[i]->SetBins(NbinsX,zMin,zMax);
      
	for(Int_t j=0;j<=hB2D[i]->GetNbinsX();j++) {
	  for(Int_t k=0;k<=hB2D[i]->GetNbinsY();k++) {
	    hB2D[i]->SetBinContent(j,k, hB2D[i]->GetBinContent(j,k) * E0  / eUnit );
	  }
	  hB1D[i]->SetBinContent(j, hB1D[i]->GetBinContent(j) * E0  / eUnit );
	  
	}
	
	if(pData->IsCyl())
	  hB2D[i]->GetYaxis()->SetTitle(Form("r [%s]",spaSUnit.c_str()));
	else
	  hB2D[i]->GetYaxis()->SetTitle(Form("x [%s]",spaSUnit.c_str()));
      
	if(opt.Contains("comov"))
	  hB2D[i]->GetXaxis()->SetTitle(Form("#zeta [%s]",spaSUnit.c_str()));
	else
	  hB2D[i]->GetXaxis()->SetTitle(Form("z [%s]",spaSUnit.c_str()));
      
	if(i==0)
	  hB2D[i]->GetZaxis()->SetTitle(Form("B_{z} [%s]",eSUnit.c_str()));
	else if(i==1)
	  hB2D[i]->GetZaxis()->SetTitle(Form("B_{x} [%s]",eSUnit.c_str()));
	else if(i==2)
	  hB2D[i]->GetZaxis()->SetTitle(Form("B_{y} [%s]",eSUnit.c_str()));
      
      
	if(opt.Contains("comov"))
	  hB1D[i]->GetXaxis()->SetTitle(Form("#zeta [%s]",spaSUnit.c_str()));
	else
	  hB1D[i]->GetXaxis()->SetTitle(Form("z [%s]",spaSUnit.c_str()));
      
	if(i==0)
	  hB1D[i]->GetYaxis()->SetTitle(Form("B_{z} [%s]",eSUnit.c_str()));
	else if(i==1)
	  hB1D[i]->GetYaxis()->SetTitle(Form("B_{x} [%s]",eSUnit.c_str()));
	else if(i==2)
	  hB1D[i]->GetYaxis()->SetTitle(Form("B_{y} [%s]",eSUnit.c_str()));
	
      }
    }

    // Focusing field
    TH2F *hFocus2D = NULL;
    if(hE2D[1]) {
      hFocus2D = (TH2F*) hE2D[1]->Clone("hFocus2D");
      
      if(hB2D[2]) {
	hFocus2D->Add(hB2D[2],-1);
	if(opt.Contains("units")) {
	  hFocus2D->GetZaxis()->SetTitle(Form("E_{x}-cB_{y} [%s]",eSUnit.c_str()));
	} else {
	  hFocus2D->GetZaxis()->SetTitle("(E_{x}-cB_{y})/E_{0}");
	}
      } else {
	delete hFocus2D;
	hFocus2D = NULL;
      }
      
    } else if(pData->GetWfieldFileName(0)) {

      cout << Form(" -> Getting transverse wakefield ") << endl;
    
      char hName[24];
      sprintf(hName,"hFocus2D");
      hFocus2D = (TH2F*) gROOT->FindObject(hName);
      if(hFocus2D) delete hFocus2D;
    
      if(!pData->Is3D())
	hFocus2D = pData->GetBField(0,opt);
      else
	hFocus2D = pData->GetWField2DSliceZX(0,-1,NonBin,opt+"avg");

      hFocus2D->SetName(hName);   
      hFocus2D->GetXaxis()->CenterTitle();
      hFocus2D->GetYaxis()->CenterTitle();
      hFocus2D->GetZaxis()->CenterTitle();

      if(opt.Contains("units") && n0) {
	Int_t NbinsX = hFocus2D->GetNbinsX();
	Double_t zMin = skindepth * hFocus2D->GetXaxis()->GetXmin() / spaUnit;
	Double_t zMax = skindepth * hFocus2D->GetXaxis()->GetXmax() / spaUnit;
	Int_t NbinsY = hFocus2D->GetNbinsY();
	Double_t ymin = skindepth * hFocus2D->GetYaxis()->GetXmin() / spaUnit;
	Double_t ymax = skindepth * hFocus2D->GetYaxis()->GetXmax() / spaUnit;
	hFocus2D->SetBins(NbinsX,zMin,zMax,NbinsY,ymin,ymax);
	
	for(Int_t j=0;j<=hFocus2D->GetNbinsX();j++) {
	  for(Int_t k=0;k<=hFocus2D->GetNbinsY();k++) {
	    hFocus2D->SetBinContent(j,k, hFocus2D->GetBinContent(j,k) * E0 / eUnit );
	  }
	}
      
	if(pData->IsCyl())
	  hFocus2D->GetYaxis()->SetTitle(Form("r [%s]",spaSUnit.c_str()));
	else
	  hFocus2D->GetYaxis()->SetTitle(Form("x [%s]",spaSUnit.c_str()));
      
	if(opt.Contains("comov"))
	  hFocus2D->GetXaxis()->SetTitle(Form("#zeta [%s]",spaSUnit.c_str()));
	else
	  hFocus2D->GetXaxis()->SetTitle(Form("z [%s]",spaSUnit.c_str()));
      
	hFocus2D->GetZaxis()->SetTitle(Form("E_{x}-cB_{y} [%s]",eSUnit.c_str()));

      } else {
	if(opt.Contains("comov"))
	  hFocus2D->GetXaxis()->SetTitle("k_{p} #zeta");
	else
	  hFocus2D->GetXaxis()->SetTitle("k_{p} z");
	
	if(pData->IsCyl()) 
	  hFocus2D->GetYaxis()->SetTitle("k_{p} r");
	else
	  hFocus2D->GetYaxis()->SetTitle("k_{p} x");
	
	hFocus2D->GetZaxis()->SetTitle("(E_{x}-cB_{y})/E_{0}");
      }
    }
    
    TH1F *hFocus1D = NULL;
    if(hE1D[1]) {
      hFocus1D = (TH1F*) hE1D[1]->Clone("hFocus1D");

      if(hB1D[2]) {
	hFocus1D->Add(hB1D[2],-1);
	if(opt.Contains("units")) {
	  hFocus1D->GetZaxis()->SetTitle(Form("E_{x}-cB_{y} [%s]",eSUnit.c_str()));
	} else {
	  hFocus1D->GetZaxis()->SetTitle("(E_{x}-cB_{y})/E_{0}");
	}
      } else {
	delete hFocus1D;
	hFocus1D = NULL;
      }

    } else if(pData->GetWfieldFileName(0)) {
      
      // 1D histograms
      if(pData->Is3D()) {
	// In case of transverse fields, the 1D line is taken offaxis
	hFocus1D = pData->GetH1SliceZ3D(pData->GetWfieldFileName(0)->c_str(),"ExmBy",-NofBin,NonBin,-1,NonBin,opt+"avg");
      } else if(pData->IsCyl()) { // Cylindrical: The first bin with r>0 is actually the number 1 (not the 0).
	hFocus1D = pData->GetH1SliceZ(pData->GetWfieldFileName(0)->c_str(),"ExmBy",1,NonBin,opt+"avg");
      } else { // 2D cartesian
	hFocus1D = pData->GetH1SliceZ(pData->GetWfieldFileName(0)->c_str(),"ExmBy",-NonBin,NonBin,opt+"avg");    
      }
      hFocus1D->SetName("hFocus1D");

      if(opt.Contains("units") && n0) {
	Int_t NbinsX = hFocus1D->GetNbinsX();
	Double_t zMin = skindepth * hFocus1D->GetXaxis()->GetXmin() / spaUnit;
	Double_t zMax = skindepth * hFocus1D->GetXaxis()->GetXmax() / spaUnit;
	hFocus1D->SetBins(NbinsX,zMin,zMax);
	
	for(Int_t j=0;j<=hFocus1D->GetNbinsX();j++) {
	  hFocus1D->SetBinContent(j,hFocus1D->GetBinContent(j) * E0 / eUnit );
	}
      
	if(opt.Contains("comov"))
	  hFocus1D->GetXaxis()->SetTitle(Form("#zeta [%s]",spaSUnit.c_str()));
	else
	  hFocus1D->GetXaxis()->SetTitle(Form("z [%s]",spaSUnit.c_str()));
	
	hFocus1D->GetYaxis()->SetTitle(Form("E_{x}-cB_{y} [%s]",eSUnit.c_str()));
	
      } else {
	if(opt.Contains("comov"))
	  hFocus1D->GetXaxis()->SetTitle("#zeta [c/#omega_{p}]");
	else
	  hFocus1D->GetXaxis()->SetTitle("z [c/#omega_{p}]");
	
	hFocus1D->GetYaxis()->SetTitle("(E_{x}-cB_{y})/E_{0}");
      }
    }
      
    
    // RMS (vs z) of the beam's charge distribution: 
    TProfile *hDen2Dprof = NULL;
    TH1F *hRms = NULL;
    Double_t axisPos = xMid;
    for(Int_t i=0;i<Nspecies;i++) {
      
      if(hDen2D[i] && hDen1D[i]) {
	if(pData->GetSpeciesName(i).find("beam") == string::npos)
	  continue;
	
	TString pname = hDen2D[i]->GetName();
	pname += "_pfx";
      
	hDen2Dprof =  (TProfile*) gROOT->FindObject(pname.Data());
	if(hDen2Dprof) { delete hDen2Dprof; hDen2Dprof = NULL; }
	hDen2Dprof = hDen2D[i]->ProfileX("_pfx",1,-1,"s");
      
	char hName[24];
	sprintf(hName,"hRms_%i",1);
	hRms = (TH1F*) gROOT->FindObject(hName);
	if(hRms) delete hRms;
      
	hRms = (TH1F*) hDen1D[i]->Clone(hName);
	hRms->Reset();
      
	if(pData->IsCyl()) axisPos = 0.0;
      
	for(Int_t j=0;j<hRms->GetNbinsX();j++) {
	  Double_t rms = 0;
	  Double_t total = 0;
	  for(Int_t k=1;k<=hDen2D[i]->GetNbinsY();k++) {
	    Double_t value  = hDen2D[i]->GetBinContent(j,k);
	    Double_t radius = hDen2D[i]->GetYaxis()->GetBinCenter(k) - axisPos;
	    if(pData->IsCyl()) {
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
	
	  hRms->SetBinContent(j,rms); 
	
	}
      
	if(opt.Contains("units")) {
	  hRms->GetYaxis()->SetTitle(Form("#Deltay_{rms} [%s]",spaSUnit.c_str()));
	} else {
	  hRms->GetYaxis()->SetTitle("k_{p}#Deltay_{rms}");
	}
      }
    }


    // Now, combine the electric field components into the total |E|
    TH2F *hETotal2D = NULL;
    // if(pData->Is3D() && hE2D[0] && hE2D[1] && hE2D[2])
    //   hETotal2D = (TH2F*) hE2D[0]->Clone("hETotal2D");
    // else if(!pData->Is3D() && hE2D[0] && hE2D[1])
    //   hETotal2D = (TH2F*) hE2D[0]->Clone("hETotal2D");

    if(hE2D[0] && hE2D[1])
      hETotal2D = (TH2F*) hE2D[0]->Clone("hETotal2D");

    
    if(hETotal2D) {
      hETotal2D->Reset();
    
      Int_t NbinsX = hE2D[0]->GetNbinsX();
      Int_t NbinsY = hE2D[0]->GetNbinsY();    
      for(Int_t j=1;j<=NbinsY;j++) {     
	for(Int_t k=NbinsX;k>0;k--) {
	  Double_t E1 = 0;
	  if(hE2D[0])
	    E1 = hE2D[0]->GetBinContent(k,j);
	  Double_t E2 = 0;
	  if(hE2D[1])
	    E2 = hE2D[1]->GetBinContent(k,j);
	  Double_t E3 = 0;
	  if(hE2D[2])
	    E3 = hE2D[2]->GetBinContent(k,j);
	  
	  Double_t E  = TMath::Sqrt(E1*E1+E2*E2+E3*E3);
	  hETotal2D->SetBinContent(k,j,E);
	}
      }
         
      if(opt.Contains("units")) {
	hETotal2D->GetZaxis()->SetTitle(Form("E [%s]",eSUnit.c_str()));
      } else {
	hETotal2D->GetZaxis()->SetTitle("E/E_{0}");
      }
    }
  
    TH1F *hETotal1D = NULL;
    if(hE1D[0] && hE1D[1])
      hETotal1D = (TH1F*) hE1D[0]->Clone("hETotal1D");

    if(hETotal1D) {

      hETotal1D->Reset();
      
      Int_t NbinsX = hE1D[0]->GetNbinsX();
      for(Int_t j=NbinsX;j>=0;j--) {
	Double_t E1 = 0;
	if(hE1D[0])
	  E1 = hE1D[0]->GetBinContent(j);
	Double_t E2 = 0;
	if(hE1D[1])
	  E2 = hE1D[1]->GetBinContent(j);
	Double_t E3 = 0;
	if(hE1D[2])
	  E3 = hE1D[2]->GetBinContent(j);
	
	Double_t E  = TMath::Sqrt(E1*E1+E2*E2+E3*E3);
	
	hETotal1D->SetBinContent(j,E);	
      }
      
      if(opt.Contains("units")) {
	hETotal1D->GetYaxis()->SetTitle(Form("E [%s]",eSUnit.c_str()));
      } else {
	hETotal1D->GetYaxis()->SetTitle("E/E_{0}");
      }
    }
    
    // Potential
    if(hE2D[0] && hE1D[0]) {
      Int_t   NbinsX = hE2D[0]->GetNbinsX();
      Int_t   NbinsY = hE2D[0]->GetNbinsY();
      
      // Double_t dx = pData->GetDX(0);
      Double_t dx = hE2D[0]->GetXaxis()->GetBinWidth(0);
      //      if(opt.Contains("units")) dx *= skindepth / PUnits::m;
      
      char hName[24];
      sprintf(hName,"hV2D");
      hV2D = (TH2F*) hE2D[0]->Clone(hName);
      hV2D->Reset();
      
      sprintf(hName,"hV1D");
      hV1D = (TH1F*) hE1D[0]->Clone(hName);
      hV1D->Reset();
      
      for(Int_t j=NbinsY;j>0;j--) {
	Double_t integral = 0.0;
	for(Int_t k=NbinsX;k>0;k--) {
	  Double_t value = hE2D[0]->GetBinContent(k,j);
	  if(opt.Contains("units")) { 
	    integral += (value * dx) * eUnit * spaUnit / (PUnits::MV);
	  } else {
	    integral += (value * dx);
	  }
	  hV2D->SetBinContent(k,j,integral);
	}
      }
      
      Double_t integral = 0.0;
      for(Int_t k=NbinsX;k>0;k--) {
	Double_t value = hE1D[0]->GetBinContent(k);
	if(opt.Contains("units")) { 
	  integral += (value * dx) * eUnit * spaUnit / (PUnits::MV);
	} else {
	  integral += (value * dx);
	}
	hV1D->SetBinContent(k,integral);
      }
    }

    // Shift potential
    // Double_t trapPot = trapPotential;
    // if(opt.Contains("units") && n0) {
    //   trapPot *=  ( E0 * skindepth / (PUnits::MV) ); 
    // }  
  
    // Double_t Vmin = hV1D->GetMinimum();    
    // { // Shift potential
    //   Int_t   NbinsX = hV2D->GetNbinsX(); 
    //   Int_t   NbinsY = hV2D->GetNbinsY(); 
    //   for(Int_t j=1;j<=NbinsX;j++) {
    // 	for(Int_t k=1;k<=NbinsY;k++) {
    // 	  hV2D->SetBinContent(j,k, hV2D->GetBinContent(j,k) - Vmin -trapPot);
    // 	}
    // 	hV1D->SetBinContent(j, hV1D->GetBinContent(j) - Vmin -trapPot);
    //   }
    // }

    if(opt.Contains("units")) {
      hV2D->GetZaxis()->SetTitle("#Psi [MV]");
      hV1D->GetYaxis()->SetTitle("#Psi [MV]");
    } else {
      hV2D->GetZaxis()->SetTitle("#psi");
      hV1D->GetYaxis()->SetTitle("#psi");
    }


    // Change the range of z axis for the fields to be symmetric.
    Double_t *Emax = new Double_t[Nfields];
    Double_t *Emin = new Double_t[Nfields];
    for(Int_t i=0;i<Nfields;i++) {

      if(!hE2D[i]) continue;
      
      Emax[i] = hE2D[i]->GetMaximum();
      Emin[i] = hE2D[i]->GetMinimum();
      if(Emax[i] > TMath::Abs(Emin[i]))
	Emin[i] = -Emax[i];
      else
	Emax[i] = -Emin[i];
      hE2D[i]->GetZaxis()->SetRangeUser(Emin[i],Emax[i]); 
      //hE2D[i]->GetZaxis()->SetRangeUser(Emin[0],Emax[0]); 
    }

    // Double_t ETmin = 0.001;  
    // Double_t ETmax = hETotal2D->GetMaximum();
    // hETotal2D->GetZaxis()->SetRangeUser(ETmin,ETmax);

    // Save output to file
    
    TString filename = Form("./%s/Plots/Snapshots/Snapshot-%s_%i.root",sim.Data(),sim.Data(),time);
    cout << Form("\n Saving snapshot objects to %s",filename.Data()) << endl;

    TString f = filename;
    TString dir2 = f.Remove( f.Last( '/' ), f.Length() - f.Last( '/' ) );
    TString dir1 = f.Remove( f.Last( '/' ), f.Length() - f.Last( '/' ) );
      
    gSystem->mkdir( dir1 );
    gSystem->mkdir( dir2 );
      
    TFile *ofile = new TFile(filename,"RECREATE");

    char hName[24];
    for(Int_t i=0;i<Nspecies;i++) {
      if(hDen2D[i]) {
	sprintf(hName,"hDen2D_%i",i);
	hDen2D[i]->Write(hName,TObject::kOverwrite);
      }

      if(hDen1D[i]) {
	sprintf(hName,"hDen1D_%i",i);
	hDen1D[i]->Write(hName,TObject::kOverwrite);
      }

      if(hCur1D[i]) {
	sprintf(hName,"hCur1D_%i",i);
	hCur1D[i]->Write(hName,TObject::kOverwrite);
      }

      if(hJz2D[i]) {
	sprintf(hName,"hJz2D_%i",i);
	hJz2D[i]->Write(hName,TObject::kOverwrite);
      }
	
    }

      

    for(Int_t i=0;i<Nfields;i++) {

      if(hE2D[i]) {
	sprintf(hName,"hE2D_%i",i);
	hE2D[i]->Write(hName,TObject::kOverwrite);
      }

      if(hE1D[i]) {
	sprintf(hName,"hE1D_%i",i);
	hE1D[i]->Write(hName,TObject::kOverwrite);
      }

    }

    for(Int_t i=0;i<NBfields;i++) {

      if(hB2D[i]) {
	sprintf(hName,"hB2D_%i",i);
	hB2D[i]->Write(hName,TObject::kOverwrite);
      }

      if(hB1D[i]) {
	sprintf(hName,"hB1D_%i",i);
	hB1D[i]->Write(hName,TObject::kOverwrite);
      }
	

    }

    if(hFocus2D) {
      sprintf(hName,"hFocus2D");
      hFocus2D->Write(hName,TObject::kOverwrite);
    }

    if(hFocus1D) {
      sprintf(hName,"hFocus1D");
      hFocus1D->Write(hName,TObject::kOverwrite);
    }

    if(hV2D) {
      sprintf(hName,"hV2D");
      hV2D->Write(hName,TObject::kOverwrite);
    }
     
    if(hETotal2D) {
      sprintf(hName,"hETotal2D");
      hETotal2D->Write(hName,TObject::kOverwrite);
    }

    if(hV1D) {
      sprintf(hName,"hV1D");
      hV1D->Write(hName,TObject::kOverwrite);
    }
 
    if(hETotal1D) {
      sprintf(hName,"hETotal1D");
      hETotal1D->Write(hName,TObject::kOverwrite);
    }

    // Save TTree with options information for further processing.
    TTree *infotree = new TTree("infoTree","options tree from doSnapshot");
    infotree->Branch("xon",&xon,"xon/D");
    infotree->Branch("xoff",&xoff,"xoff/D");
    infotree->Branch("xint",&xint,"xint/D");
    infotree->Branch("denUnit",&denUnit,"denUnit/D");
    infotree->Branch("denSUnit","string",&denSUnit);
    infotree->Branch("spaUnit",&spaUnit,"spaUnit/D");
    infotree->Branch("spaSUnit","string",&spaSUnit);
    infotree->Branch("propUnit",&propUnit,"propUnit/D");
    infotree->Branch("propSUnit","string",&propSUnit);
    infotree->Branch("timUnit",&timUnit,"timUnit/D");
    infotree->Branch("timSUnit","string",&timSUnit);
    infotree->Branch("eUnit",&eUnit,"eUnit/D");
    infotree->Branch("eSUnit","string",&eSUnit);
    infotree->Branch("curUnit",&curUnit,"curUnit/D");
    infotree->Branch("curSUnit","string",&curSUnit);

    TString *optpoint = &opt;
    infotree->Branch("options","TString",&optpoint);

    infotree->Fill();
    if(infotree) {
      infotree->Write("infotree",TObject::kOverwrite);
    }
    
    
    // Plot the snapshot!
    if(opt.Contains("pdf") || opt.Contains("eps") || opt.Contains("png")) {
      
      TString command = Form(".x rootlogon.C");
      gROOT->ProcessLine(command);

      command = Form(".x PlotSnapshot.C(\"%s\", %i , %i, \"%s\")",
		     sim.Data(),time,mask,opt.Data());
      gROOT->ProcessLine(command);
      
    }
    
    ofile->Close();
        
  }



}
