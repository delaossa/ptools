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

#include "PData.hh"
#include "PGlobals.hh"
#include "PFunctions.hh"
#include "PPalette.hh"
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

  // z start of the plasma in normalized units.
  Double_t zStartPlasma = pData->GetPlasmaStart()*kp;
  // z start of the beam in normalized units.
  Double_t zStartBeam = pData->GetBeamStart()*kp;
 
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

    if(opt.Contains("center")) {
      Time -= zStartPlasma;
      if(opt.Contains("comov"))      // Centers on the head of the beam.
	Time += zStartBeam;
    } 
    
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
    Int_t FirstyBin = 0;
    Int_t LastyBin = 0;
    if(!pData->IsCyl()) {
      FirstyBin = pData->GetNX(1)/2 + 1 - NonBin;
      LastyBin =  pData->GetNX(1)/2 + NonBin;
    } else {
      FirstyBin = 1; 
      LastyBin  = NonBin;
    }
    Int_t FirstOffyBin = FirstyBin + NofBin;
    Int_t LastOffyBin = LastyBin + NofBin;

    // Zoom window:
  
    Double_t yRange = (pData->GetXMax(1) - pData->GetXMin(1))/zoom;
    Double_t yMid   = (pData->GetXMax(1) + pData->GetXMin(1))/2.;
    //Double_t yMin = pData->GetXMin(1);
    //Double_t yMax = pData->GetXMax(1);
    Double_t yMin = yMid - yRange/2.0;
    Double_t yMax = yMid + yRange/2.0;
    if(pData->IsCyl()) {
      yMin = pData->GetXMin(1);
      yMax = yRange;
    }
    pData->SetX2Min(yMin);
    pData->SetX2Max(yMax);
  
    Double_t zMin = pData->GetX1Min();
    Double_t zMax = pData->GetX1Max();
    // Double_t zRange = zMax - zMin;  
    // pData->SetX1Min(zMin);
    // pData->SetX1Max(zMax);

    cout << Form(" --- doSnapshot ---\n") << endl;
  
    cout << Form(" Plotting range:  %.2f < x1 < %.2f ,  %.2f < x2 < %.2f",zMin,zMax,yMin,yMax) << endl;

    // ----------------------------------------------------------------------------------
  

    cout << Form("\n Reading simulation data: ") << endl; 

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
      // hDen1D[i] = (TH1F*) hE2D[i]->ProjectionX(hName,FirstyBin,LastyBin);
      // hDen1D[i]->Scale(1.0/(LastyBin-FirstyBin+1));
    
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
	hCur1D[i]->Scale(TMath::Abs(n0 * dV * PConst::ElectronCharge * kp * PConst::c_light) / PConst::I0);
	
	hCur1D[i]->GetYaxis()->SetTitle("#Lambda_{b}");  
	if(opt.Contains("comov")) {
	  hCur1D[i]->GetXaxis()->SetTitle("k_{p} #zeta");
	} else {
	  hCur1D[i]->GetXaxis()->SetTitle("k_{p} z");
	}
     	
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
      char nam[3]; sprintf(nam,"e%i",i+1);
      if(pData->Is3D()) {
      
	if(!opt.Contains("alt1D")){
	  if(i==0) 
	    hE1D[i] = pData->GetH1SliceZ3D(pData->GetEfieldFileName(i)->c_str(),nam,-1,NonBin,-1,NonBin,opt+"avg");
	  else // In case of transverse fields, the 1D line is taken off-axis 
	    hE1D[i] = pData->GetH1SliceZ3D(pData->GetEfieldFileName(i)->c_str(),nam,-NofBin,NonBin,-1,NonBin,opt+"avg");
	
	} else { // Alternative
	  cout << " Alternative! " << endl;
	  if(i==0) {
	    hE1D[i] = (TH1F*) hE2D[i]->ProjectionX(hName,FirstyBin,LastyBin);
	    hE1D[i]->Scale(1.0/(LastyBin-FirstyBin+1));
	  } else { // In case of transverse fields, the 1D line is taken off-axis 
	    hE1D[i] = (TH1F*) hE2D[i]->ProjectionX(hName,FirstOffyBin,LastOffyBin);
	    hE1D[i]->Scale(1.0/(LastOffyBin-FirstOffyBin+1));
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
      // hE1D[i] = (TH1F*) hE2D[i]->ProjectionX(hName,FirstyBin,LastyBin);
      // hE1D[i]->Scale(1.0/(LastyBin-FirstyBin+1));
     
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
      if(pData->Is3D()) {
	
	if(!opt.Contains("alt1D")){
	  if(i==0) 
	    hB1D[i] = pData->GetH1SliceZ3D(pData->GetBfieldFileName(i)->c_str(),bnam,-1,NonBin,-1,NonBin,opt+"avg");
	  else // In case of transverse fields, the 1D line is taken offaxis
	    hB1D[i] = pData->GetH1SliceZ3D(pData->GetBfieldFileName(i)->c_str(),bnam,-NofBin,NonBin,-1,NonBin,opt+"avg");
	  
	} else {  // Alternative
	  cout << " Alternative! " << endl;
	  if(i==0) {
	    hB1D[i] = (TH1F*) hB2D[i]->ProjectionX(hName,FirstyBin,LastyBin);
	    hB1D[i]->Scale(1.0/(LastyBin-FirstyBin+1));
	  } else { // In case of transverse fields, the 1D line is taken off-axis 
	    hB1D[i] = (TH1F*) hB2D[i]->ProjectionX(hName,FirstOffyBin,LastOffyBin);
	    hB1D[i]->Scale(1.0/(LastOffyBin-FirstOffyBin+1));
	  }
	}


      } else if(pData->IsCyl()) { // Cylindrical: The first bin with r>0 is actually the number 1 (not the 0).
      
	hB1D[i] = pData->GetH1SliceZ(pData->GetBfieldFileName(i)->c_str(),bnam,1,NonBin,opt+"avg");
      
      } else { // 2D cartesian
      
	if(i==0) 
	  hB1D[i] = pData->GetH1SliceZ(pData->GetBfieldFileName(i)->c_str(),bnam,-1,NonBin,opt+"avg");
	else 
	  hB1D[i] = pData->GetH1SliceZ(pData->GetBfieldFileName(i)->c_str(),bnam,-NonBin,NonBin,opt+"avg");    
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

    // Chaning to user units: 
    // --------------------------
  
    Double_t denUnit, spaUnit, timUnit, eUnit, curUnit;
    string denSUnit, spaSUnit, timSUnit, eSUnit, curSUnit;
    
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
	  
	  PUnits::BestUnit bcurSUnit(hCur1D[i]->GetMaximum() * PConst::I0,"Current");
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
    TH2F *hFocus2D = (TH2F*) hE2D[1]->Clone("hFocus2D");
    if(hB2D[2]) {
      hFocus2D->Add(hB2D[2],-1);
      if(opt.Contains("units")) {
	hFocus2D->GetZaxis()->SetTitle(Form("E_{x}-cB_{y} [%s]",eSUnit.c_str()));
      } else {
	hFocus2D->GetZaxis()->SetTitle("(E_{x}-cB_{y})/E_{0}");
      }
    }

    TH1F *hFocus1D = (TH1F*) hE1D[1]->Clone("hFocus1D");
    if(hB1D[2]) {
      hFocus1D->Add(hB1D[2],-1);
      if(opt.Contains("units")) {
	hFocus1D->GetZaxis()->SetTitle(Form("E_{x}-cB_{y} [%s]",eSUnit.c_str()));
      } else {
	hFocus1D->GetZaxis()->SetTitle("(E_{x}-cB_{y})/E_{0}");
      }
    }
    
    // RMS (vs z) of the beam's charge distribution: 
    TProfile *hDen2Dprof = NULL;
    TH1F *hRms = NULL;
    Double_t axisPos = yMid;
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
    TH2F *hETotal2D = (TH2F*) hE2D[0]->Clone("hETotal2D");
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
  
    TH1F *hETotal1D = (TH1F*) hE1D[0]->Clone("hETotal1D");
    hETotal1D->Reset();

    NbinsX = hE1D[0]->GetNbinsX();
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
