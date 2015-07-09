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
#include "PPalette.hh"
#include "H5Cpp.h"

using namespace std;
using namespace H5;


int main(int argc,char *argv[]) {
  if(argc<=2) {
    printf("\n Usage: %s <simulation name> <-t(time)>\n",argv[0]);
    printf("      <-i(initial time)> <-f(final time)> <-s(time step)>\n");
    printf("      <-z(zoom factor)> <-non(# on-axis)> <-nof(# off-axis)>");
    printf("      <--center> <--comov>\n");
    printf("      <--units>\n");
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
  
  // Interfacing command line:
  for(int l=1;l<argc;l++){
    TString arg = argv[l];

    if(arg.Contains("--comov")){
      opt += "comov";
    } else if(arg.Contains("--units")){
      opt += "units";
    } else if(arg.Contains("--center")){
      opt += "center";
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
    } else if( !(arg.Contains("pitz") || arg.Contains("flash") || arg.Contains("regae") || arg.Contains("Gauss") || arg.Contains("pwfa") || arg.Contains("vacuum") || arg.Contains("seed") || arg.Contains("facet") || arg.Contains("FACET") || arg.Contains("rake")) ) {
      cout << Form("\t Invalid argument (%i): exiting...\n",l) << endl;
      return 0;
    } else {
      sim = arg;
    }
  }
  
  // This program always has to use the comov option.
  opt += "comov";
  

  PGlobals::Initialize();
  
  // Palettes!
  gROOT->Macro("PPalettes.C");

  if(opt.Contains("grid")) {
    gStyle->SetPadGridX(1);
    gStyle->SetPadGridY(1);
  }

  // Load PData
  PData *pData = PData::Get(sim.Data());

  if(iStart<0) iStart = time;
  if(iEnd<=iStart) iEnd = iStart;
  
  // Some plasma constants
  Float_t n0 = pData->GetPlasmaDensity();
  Float_t omegap = pData->GetPlasmaFrequency();
  Float_t timedepth = pData->GetPlasmaTimeDepth();
  Float_t kp = pData->GetPlasmaK();
  Float_t skindepth = pData->GetPlasmaSkinDepth();
  Float_t E0 = pData->GetPlasmaE0();

  // Some beam properties:
  Float_t gamma = pData->GetBeamGamma();
  
  // Other parameters
  Float_t trapPotential = 1.0 - (1.0/gamma);
  cout << Form(" - Trap. potential  = %8.4f mc2/e",trapPotential) << endl;
  cout << endl;

  // z start of the plasma in normalized units.
  Float_t zStartPlasma = pData->GetPlasmaStart()*kp;
  // z start of the beam in normalized units.
  Float_t zStartBeam = pData->GetBeamStart()*kp;
 
  // Time looper
  for(Int_t i=iStart; i<iEnd+1; i+=iStep) {

    time = i;

    cout << Form("\n Processing time step %i",time) << endl;
    cout << Form(" -------------------------\n") << endl;

    pData->LoadFileNames(time);    
    //if(time==iStart) pData->PrintData();
    
    if(!pData->IsInit()) continue;

    // Time in OU
    Float_t Time = pData->GetRealTime();

    if(opt.Contains("center")) {
      Time -= zStartPlasma;
      if(opt.Contains("comov"))      // Centers on the head of the beam.
	Time += zStartBeam;
    } 

    // Centering time and z position:
    // Float_t shiftz = pData->Shift(opt);
    Float_t rms0 = pData->GetBeamRmsY() * kp;

    // Calculate the "axis range" in number of bins. 
    if(NonBin==0) {  // If NonBin==0 a RMS width is taken.
      if(pData->IsCyl()) rms0 = pData->GetBeamRmsR() * kp;
      
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

    // Zoom window:  
    Float_t yRange = (pData->GetXMax(1) - pData->GetXMin(1))/zoom;
    Float_t yMid   = (pData->GetXMax(1) + pData->GetXMin(1))/2.;
    //Float_t yMin = pData->GetXMin(1);
    //Float_t yMax = pData->GetXMax(1);
    Float_t yMin = yMid - yRange/2.0;
    Float_t yMax = yMid + yRange/2.0;
    if(pData->IsCyl()) {
      yMin = pData->GetXMin(1);
      yMax = yRange;
    }
    pData->SetX2Min(yMin);
    pData->SetX2Max(yMax);
  
    Float_t zMin = pData->GetX1Min();
    Float_t zMax = pData->GetX1Max();
    // pData->SetX1Min(zMin);
    // pData->SetX1Max(zMax);

    cout << Form(" Plotting range:  %.2f < x1 < %.2f ,  %.2f < x2 < %.2f",zMin,zMax,yMin,yMax) << endl;

    // ----------------------------------------------------------------------------------
  
  
    // Get charge density histos
    Int_t Nspecies = pData->NSpecies();
    TH2F **hDen2D = new TH2F*[Nspecies];
    // Get charge density on-axis
    TH1F **hDen1D = new TH1F*[Nspecies];
    for(Int_t i=0;i<Nspecies;i++) {
     
      hDen2D[i] = NULL;
      hDen1D[i] = NULL;
            
      if(!pData->GetChargeFileName(i)) 
	continue;

      cout << Form(" Getting charge density of specie: ") << i << endl;
    
      char hName[24];
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

      // Get the 2D den for the beam only:
      if(i==0) continue;

      sprintf(hName,"hDen2D_%i",i);
      hDen2D[i] = (TH2F*) gROOT->FindObject(hName);
      if(hDen2D[i]) delete hDen2D[i];
      
      if(!pData->Is3D())
	hDen2D[i] = pData->GetCharge(i,opt);
      else {
      	//hDen2D[i] = pData->GetCharge2DSliceZY(i,-1,NonBin,opt+"avg");
	hDen2D[i] = pData->GetH2SliceZY(pData->GetChargeFileName(i)->c_str(),"charge",-1,NonBin,opt+"avg");
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
      
      hDen2D[i]->GetZaxis()->SetTitle("n/n_{0}");
      
    }
  
  
    // Get electric fields 1D
    const Int_t Nfields = 3;
    TH1F **hE1D = new TH1F*[Nfields];
    TH1F *hV1D = NULL;
    for(Int_t i=0;i<Nfields;i++) {
      hE1D[i] = NULL;

      if(!pData->GetEfieldFileName(i))
	continue;

      cout << Form(" Getting electric field number ") << i+1 << endl;
    
      char hName[24];
      sprintf(hName,"hE1D_%i",i);
      hE1D[i] = (TH1F*) gROOT->FindObject(hName);
      if(hE1D[i]) delete hE1D[i];
      
      // 1D histograms
      char nam[3]; sprintf(nam,"e%i",i+1);
      if(pData->Is3D()) {
      
	if(i==0) 
	  hE1D[i] = pData->GetH1SliceZ3D(pData->GetEfieldFileName(i)->c_str(),nam,-1,NonBin,-1,NonBin,opt+"avg");
	else // In case of transverse fields, the 1D line is taken 1 k_p^-1 aside of the axis
	  hE1D[i] = pData->GetH1SliceZ3D(pData->GetEfieldFileName(i)->c_str(),nam,-NofBin,NonBin,-1,NonBin,opt+"avg");
      
      } else if(pData->IsCyl()) { // Cylindrical: The first bin with r>0 is actually the number 1 (not the 0).
      
	hE1D[i] = pData->GetH1SliceZ(pData->GetEfieldFileName(i)->c_str(),nam,1,NonBin,opt+"avg");
      
      } else { // 2D cartesian
      
	if(i==0) 
	  hE1D[i] = pData->GetH1SliceZ(pData->GetEfieldFileName(i)->c_str(),nam,-1,NonBin,opt+"avg");
	else 
	  hE1D[i] = pData->GetH1SliceZ(pData->GetEfieldFileName(i)->c_str(),nam,-NonBin,NonBin,opt+"avg");    
      }
    
      hE1D[i]->SetName(hName);
      if(opt.Contains("comov"))
	hE1D[i]->GetXaxis()->SetTitle("k_{p} #zeta");
      else
	hE1D[i]->GetXaxis()->SetTitle("k_{p} z");
    
      if(i==0)
	hE1D[i]->GetYaxis()->SetTitle("E_{z}/E_{0}");
      else if(i==1)
	hE1D[i]->GetYaxis()->SetTitle("E_{x}/E_{0}");
      else if(i==2)
	hE1D[i]->GetYaxis()->SetTitle("E_{y}/E_{0}");
    
    }

    // Get magnetic fields 1D
    const Int_t NBfields = 3;
    TH1F **hB1D = new TH1F*[NBfields];
    for(Int_t i=0;i<NBfields;i++) {
      hB1D[i] = NULL;
      if(i<2) continue;  // Just get the third component.
      
      if(!pData->GetBfieldFileName(i))
	continue;
      
      cout << Form(" Getting magnetic field number ") << i+1 << endl;
      
      char hName[24];
      sprintf(hName,"hB1D_%i",i);
      hB1D[i] = (TH1F*) gROOT->FindObject(hName);
      if(hB1D[i]) delete hB1D[i];
      
      char bnam[3]; sprintf(bnam,"b%i",i+1);
      // 1D histograms
      if(pData->Is3D()) {
	if(i==0) 
	  hB1D[i] = pData->GetH1SliceZ3D(pData->GetBfieldFileName(i)->c_str(),bnam,-1,NonBin,-1,NonBin,opt+"avg");
	else // In case of transverse fields, the 1D line is taken offaxis
	  hB1D[i] = pData->GetH1SliceZ3D(pData->GetBfieldFileName(i)->c_str(),bnam,-NofBin,NonBin,-1,NonBin,opt+"avg");
	
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
	hB1D[i]->GetXaxis()->SetTitle("k_{p} #zeta");
      else
	hB1D[i]->GetXaxis()->SetTitle("k_{p} z");
    
      if(i==0)
	hB1D[i]->GetYaxis()->SetTitle("B_{z}/E_{0}");
      else if(i==1)
	hB1D[i]->GetYaxis()->SetTitle("B_{x}/E_{0}");
      else if(i==2)
	hB1D[i]->GetYaxis()->SetTitle("B_{y}/E_{0}");
      
      
    }

    
    // Chaning to user units: 
    // --------------------------
  
    if(opt.Contains("units") && n0) {
    
      for(Int_t i=0;i<Nspecies;i++) {
	if(!hDen1D[i]) continue;
    
	Int_t NbinsX = hDen1D[i]->GetNbinsX();
	Float_t zMin = skindepth * hDen1D[i]->GetXaxis()->GetXmin() / PUnits::um;
	Float_t zMax = skindepth * hDen1D[i]->GetXaxis()->GetXmax() / PUnits::um;
	hDen1D[i]->SetBins(NbinsX,zMin,zMax);      
	if(opt.Contains("comov"))
	  hDen1D[i]->GetXaxis()->SetTitle("#zeta [#mum]");
	else
	  hDen1D[i]->GetXaxis()->SetTitle("z [#mum]");
	 
	if(!hDen2D[i]) continue;
	Int_t NbinsY = hDen2D[i]->GetNbinsY();
	Float_t ymin = skindepth * hDen2D[i]->GetYaxis()->GetXmin() / PUnits::um;
	Float_t ymax = skindepth * hDen2D[i]->GetYaxis()->GetXmax() / PUnits::um;
	hDen2D[i]->SetBins(NbinsX,zMin,zMax,NbinsY,ymin,ymax);

	if(pData->IsCyl())
	  hDen2D[i]->GetYaxis()->SetTitle("r [#mum]");      
	else
	  hDen2D[i]->GetYaxis()->SetTitle("x [#mum]");      
      
	if(opt.Contains("comov"))
	  hDen2D[i]->GetXaxis()->SetTitle("#zeta [#mum]");
	else
	  hDen2D[i]->GetXaxis()->SetTitle("z [#mum]");
	 
      }
       
       
      for(Int_t i=0;i<Nfields;i++) {
	if(!hE1D[i]) continue;
	 
	Int_t NbinsX = hE1D[i]->GetNbinsX();
	Float_t zMin = skindepth * hE1D[i]->GetXaxis()->GetXmin() / PUnits::um;
	Float_t zMax = skindepth * hE1D[i]->GetXaxis()->GetXmax() / PUnits::um;
	hE1D[i]->SetBins(NbinsX,zMin,zMax);
            
	for(Int_t j=0;j<=hE1D[i]->GetNbinsX();j++) {
	  hE1D[i]->SetBinContent(j, hE1D[i]->GetBinContent(j) * ( E0 / (PUnits::GV/PUnits::m) ) );
	}
      
	if(opt.Contains("comov"))
	  hE1D[i]->GetXaxis()->SetTitle("#zeta [mm]");
	else
	  hE1D[i]->GetXaxis()->SetTitle("z [mm]");
      
	if(i==0)
	  hE1D[i]->GetYaxis()->SetTitle("E_{z} [GV/m]");
	else if(i==1)
	  hE1D[i]->GetYaxis()->SetTitle("E_{x} [GV/m]");
	else if(i==2)
	  hE1D[i]->GetYaxis()->SetTitle("E_{y} [GV/m]");
	
      }

      for(Int_t i=0;i<NBfields;i++) {
	if(!hB1D[i]) continue;

	Int_t NbinsX = hB1D[i]->GetNbinsX();
	Float_t zMin = skindepth * hB1D[i]->GetXaxis()->GetXmin() / PUnits::um;
	Float_t zMax = skindepth * hB1D[i]->GetXaxis()->GetXmax() / PUnits::um;
	hB1D[i]->SetBins(NbinsX,zMin,zMax);
      
	for(Int_t j=0;j<=hB1D[i]->GetNbinsX();j++) {
	  hB1D[i]->SetBinContent(j, hB1D[i]->GetBinContent(j) * ( E0 / (PUnits::GV/PUnits::m) ) );
	}
	
	if(opt.Contains("comov"))
	  hB1D[i]->GetXaxis()->SetTitle("#zeta [mm]");
	else
	  hB1D[i]->GetXaxis()->SetTitle("z [mm]");
	
	if(i==0)
	  hB1D[i]->GetYaxis()->SetTitle("B_{z} [GV/m]");
	else if(i==1)
	  hB1D[i]->GetYaxis()->SetTitle("B_{x} [GV/m]");
	else if(i==2)
	  hB1D[i]->GetYaxis()->SetTitle("B_{y} [GV/m]");
	
      }
    }
    
    // Focusing field
    TH1F *hFocus1D = (TH1F*) hE1D[1]->Clone("hFocus1D");
    if(hB1D[2]) {
      hFocus1D->Add(hB1D[2],-1);
      if(opt.Contains("units")) {
	hFocus1D->GetZaxis()->SetTitle("E_{x}-B_{y} [GV/m]");
      } else {
	hFocus1D->GetZaxis()->SetTitle("E_{x}-B_{y}/E_{0}");
      }
    }
     

    // RMS (vs z) of the beam's charge distribution: 
    TProfile *hDen2Dprof = NULL;
    TH1F *hRms = NULL;
    Double_t axisPos = yMid;
    if(hDen2D[1] && hDen1D[1]) {
      TString pname = hDen2D[1]->GetName();
      pname += "_pfx";
      
      hDen2Dprof =  (TProfile*) gROOT->FindObject(pname.Data());
      if(hDen2Dprof) { delete hDen2Dprof; hDen2Dprof = NULL; }
      hDen2Dprof = hDen2D[1]->ProfileX("_pfx",1,-1,"s");
      
      char hName[24];
      sprintf(hName,"hRms_%i",1);
      hRms = (TH1F*) gROOT->FindObject(hName);
      if(hRms) delete hRms;
      
      hRms = (TH1F*) hDen1D[1]->Clone(hName);
      hRms->Reset();
      
      if(pData->IsCyl()) axisPos = 0.0;
      
      for(Int_t j=0;j<hRms->GetNbinsX();j++) {
	Double_t rms = 0;
	Double_t total = 0;
	for(Int_t k=1;k<=hDen2D[1]->GetNbinsY();k++) {
	  Double_t value  = hDen2D[1]->GetBinContent(j,k);
	  Double_t radius = hDen2D[1]->GetYaxis()->GetBinCenter(k) - axisPos;
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
	hRms->GetYaxis()->SetTitle("#sigma_{x} [#mum]");
      } else {
	hRms->GetYaxis()->SetTitle("k_{p} #sigma_{x}");
      }
    }

    
    // Now, combine the electric field components into the total |E|
    // and calculate ionization probability for the different atomic species:
    const Int_t NAtoms = 3;
    char atNames[NAtoms][4] = {"H","He","He2"};
    char atAxNames[NAtoms][8] = {"H","He","He^{+}"};
    Float_t Eion0[NAtoms] = {13.6 * PUnits::eV, 24.59 * PUnits::eV, 54.4 * PUnits::eV};
    Float_t Z[NAtoms]     = { 1.0,  1.0,  2.0};

    TH1F *hETotal1D = (TH1F*) hE1D[0]->Clone("hETotal1D");
    hETotal1D->Reset();

    TH1F *hIonRate1D[NAtoms] = {NULL};
    TH1F *hIonProb1D[NAtoms] = {NULL};
    for(Int_t iat=0; iat<NAtoms; iat++) {
      char hName[24];
      sprintf(hName,"hIonRate1D_%s",atNames[iat]);
      hIonRate1D[iat] = (TH1F*) hE1D[0]->Clone(hName);
      hIonRate1D[iat]->Reset();
      sprintf(hName,"hIonProb1D_%s",atNames[iat]);
      hIonProb1D[iat] = (TH1F*) hE1D[0]->Clone(hName);
      hIonProb1D[iat]->Reset();
    }

    {
      Float_t integral[NAtoms] = {0.0};
      Float_t dx = pData->GetDX(0);
      if(opt.Contains("units")) dx *=  timedepth / PUnits::femtosecond;
      
      Int_t NbinsX = hE1D[0]->GetNbinsX();
      for(Int_t j=NbinsX;j>=0;j--) {
	Float_t E1 = hE1D[0]->GetBinContent(j);
	Float_t E2 = hE1D[1]->GetBinContent(j);
	Float_t E3 = hE1D[2]->GetBinContent(j);
	Float_t E  = TMath::Sqrt(E1*E1+E2*E2+E3*E3);
    
	hETotal1D->SetBinContent(j,E);
	
	if(opt.Contains("units")) E *= PUnits::GV/PUnits::m;
	else E *= E0;
	if(E<10) continue; // This is a filter to not to crash the ADK function.

	for(Int_t iat=0; iat<NAtoms; iat++) {
	  Float_t IonRate = PFunc::ADK_ENG(E,Eion0[iat],Z[iat]);
	  
	  if(opt.Contains("units"))
	    IonRate /= 1.0/PUnits::femtosecond;
	  else
	    IonRate /= omegap;
	  
	  hIonRate1D[iat]->SetBinContent(j,IonRate);
	  
	  if(integral[iat]==-99) continue;
	  
	  integral[iat] += IonRate * dx;	
	  if(integral[iat]>=1) integral[iat] = 1;
	  hIonProb1D[iat]->SetBinContent(j,(1-integral[iat])*IonRate);
	  
	}
      }
    } 
    
     if(opt.Contains("units")) {
      hETotal1D->GetYaxis()->SetTitle("E_{T} [GV/m]");
      for(Int_t iat=0; iat<NAtoms; iat++) {
	char axName[24];
	sprintf(axName,"W_{%s} [fs^{-1}]",atAxNames[iat]);
	hIonRate1D[iat]->GetZaxis()->SetTitle(axName);
	sprintf(axName,"#Gamma_{%s} [fs^{-1}]",atAxNames[iat]);
	hIonProb1D[iat]->GetZaxis()->SetTitle(axName);
      }
    } else {
      hETotal1D->GetYaxis()->SetTitle("E_{}/E_{0}");
      for(Int_t iat=0; iat<NAtoms; iat++) {
	char axName[24];
	sprintf(axName,"W_{%s}/#omega_{p}",atAxNames[iat]);
	hIonRate1D[iat]->GetZaxis()->SetTitle(axName);
	sprintf(axName,"#Gamma_{%s}/#omega_{p}",atAxNames[iat]);
	hIonProb1D[iat]->GetZaxis()->SetTitle(axName);
      }
    }
 
    // Potential
    if(hE1D[0]) {
      Int_t   NbinsX = hE1D[0]->GetNbinsX();
      Float_t dx = pData->GetDX(0);
      // if(opt.Contains("units")) dx *= skindepth / PUnits::m;
      
      char hName[24];
      sprintf(hName,"hV1D");
      hV1D = (TH1F*) hE1D[0]->Clone(hName);
      hV1D->Reset();
      
      Float_t integral = 0.0;
      for(Int_t k=NbinsX;k>0;k--) {
	Float_t value = hE1D[0]->GetBinContent(k);
	if(opt.Contains("units")) { 
	  value /= E0 / (PUnits::GV/PUnits::m);
	  integral += (value * dx) * E0 * skindepth / (PUnits::MV);
	} else {
	  integral += (value * dx);
	}
	hV1D->SetBinContent(k,integral);
      }
    }
    
    Float_t trapPot = trapPotential;
    if(opt.Contains("units") && n0) {
      trapPot *=  ( E0 * skindepth / (PUnits::MV) ); 
    }  
  
    Float_t Vmin = hV1D->GetMinimum();    
    { // Shift potential
      Int_t   NbinsX = hV1D->GetNbinsX(); 
      for(Int_t j=0;j<NbinsX;j++) {
	hV1D->SetBinContent(j, hV1D->GetBinContent(j) - Vmin -trapPot);
      }
    }
    
    Vmin = hV1D->GetMinimum();   
    Float_t Vmax = hV1D->GetMaximum();    
    if(Vmax<0.1) Vmax = 0.1;
    // --

    if(opt.Contains("units")) {
      hV1D->GetYaxis()->SetTitle("#Psi-#Psi_{T} [MV]");
    } else {
      hV1D->GetYaxis()->SetTitle("#psi-#psi_{T}");
    }

    
    {
      
      // OUTPUT ROOT FILE WITH THE Evolution PLOTS:
      TString filename = Form("./%s/Plots/Evolutions/Evolutions-%s.root",sim.Data(),sim.Data());
      TFile * loopfile = (TFile*) gROOT->GetListOfFiles()->FindObject(filename);
      // if doesn't exist the directory should be created
      if (!loopfile) {
	TString f = filename;
	TString dir2 = f.Remove( f.Last( '/' ), f.Length() - f.Last( '/' ) );
	TString dir1 = f.Remove( f.Last( '/' ), f.Length() - f.Last( '/' ) );
	gSystem->mkdir( dir1 );
	gSystem->mkdir( dir2 );
	loopfile = new TFile(filename,"UPDATE");
      }  

      cout << " Saving time looping objects to " << filename << endl;

    
      if(opt.Contains("units")) Time *= skindepth / PUnits::mm;
    
      Int_t nBins   = 1;
      Float_t edge0 = Time-0.5;
      Float_t edge1 = Time+0.5;
      char hName[24];

      TH2F **hEvsTime = new TH2F*[Nfields]; 
      for(Int_t i=0;i<Nfields;i++) {
	sprintf(hName,"hEvsTime_%i",i);
	TH2F *hEvsTimeOld = (TH2F*) loopfile->Get(hName);
	if(hEvsTimeOld) {
	  nBins = hEvsTimeOld->GetNbinsX()+1;
	  Float_t binwidth =  (Time - hEvsTimeOld->GetXaxis()->GetBinCenter(1))/(nBins-1);
	  edge0 = hEvsTimeOld->GetXaxis()->GetBinCenter(1) - binwidth/2.;
	  edge1 = Time + binwidth/2.;
	}
	hEvsTime[i] = new TH2F("temp","",nBins,edge0,edge1,
			       hE1D[i]->GetNbinsX(),
			       hE1D[i]->GetBinLowEdge(1),
			       hE1D[i]->GetBinLowEdge(hE1D[i]->GetNbinsX()+1));
      
	for(Int_t ix=1;ix<hEvsTime[i]->GetNbinsX();ix++) {
	  for(Int_t iy=1;iy<hEvsTime[i]->GetNbinsY();iy++) {
	    hEvsTime[i]->SetBinContent(ix,iy,hEvsTimeOld->GetBinContent(ix,iy));
	  }
	}  
	delete hEvsTimeOld;
      
	// Fill last bin with the newest values.
	for(Int_t iy=1;iy<=hE1D[i]->GetNbinsX();iy++) {
	  hEvsTime[i]->SetBinContent(nBins,iy,hE1D[i]->GetBinContent(iy));
	}   
      
	if(opt.Contains("units")) {
	  if(i==0) 
	    hEvsTime[i]->GetZaxis()->SetTitle("E_{z} [GV/m]");
	  else if(i==1)
	    hEvsTime[i]->GetZaxis()->SetTitle("E_{x} [GV/m]");
	  else if(i==2)
	    hEvsTime[i]->GetZaxis()->SetTitle("E_{y} [GV/m]");
	
	  hEvsTime[i]->GetYaxis()->SetTitle("#zeta [#mum]");
	  hEvsTime[i]->GetXaxis()->SetTitle("z [mm]");
	} else {
	  if(i==0) 
	    hEvsTime[i]->GetZaxis()->SetTitle("E_{z}/E_{0}");
	  else if(i==1)
	    hEvsTime[i]->GetZaxis()->SetTitle("E_{x}/E_{0}");
	  else if(i==2)
	    hEvsTime[i]->GetZaxis()->SetTitle("E_{y}/E_{0}");
	
	  hEvsTime[i]->GetYaxis()->SetTitle("k_{p} #zeta");
	  hEvsTime[i]->GetXaxis()->SetTitle("k_{p} z");
	}

	hEvsTime[i]->GetZaxis()->CenterTitle();
	hEvsTime[i]->GetYaxis()->CenterTitle();
	hEvsTime[i]->GetXaxis()->CenterTitle();
	hEvsTime[i]->SetName(hName);
      
	// Change the range of z axis for the fields to be symmetric.
	Float_t Emax = hEvsTime[i]->GetMaximum();
	Float_t Emin = hEvsTime[i]->GetMinimum();
	if(Emax > TMath::Abs(Emin))
	  Emin = -Emax;
	else
	  Emax = -Emin;
	hEvsTime[i]->GetZaxis()->SetRangeUser(Emin,Emax); 
      
	hEvsTime[i]->Write(hName,TObject::kOverwrite);
      }

      // Total field
      TH2F *hETotalvsTime = NULL;
      sprintf(hName,"hETotalvsTime");
      TH2F *hETotalvsTimeOld = (TH2F*) loopfile->Get(hName);
      if(hETotalvsTimeOld) {
	nBins = hETotalvsTimeOld->GetNbinsX()+1;
	Float_t binwidth =  (Time - hETotalvsTimeOld->GetXaxis()->GetBinCenter(1))/(nBins-1);
	edge0 = hETotalvsTimeOld->GetXaxis()->GetBinCenter(1) - binwidth/2.;
	edge1 = Time + binwidth/2.;
      }
      hETotalvsTime = new TH2F("temp","",nBins,edge0,edge1,
			       hETotal1D->GetNbinsX(),
			       hETotal1D->GetBinLowEdge(1),
			       hETotal1D->GetBinLowEdge(hETotal1D->GetNbinsX()+1));
      
      for(Int_t ix=1;ix<hETotalvsTime->GetNbinsX();ix++) {
	for(Int_t iy=1;iy<hETotalvsTime->GetNbinsY();iy++) {
	  hETotalvsTime->SetBinContent(ix,iy,hETotalvsTimeOld->GetBinContent(ix,iy));
	}
      }  
      delete hETotalvsTimeOld;
      
      // Fill last bin with the newest values.
      for(Int_t iy=1;iy<=hETotal1D->GetNbinsX();iy++) {
	hETotalvsTime->SetBinContent(nBins,iy,hETotal1D->GetBinContent(iy));
      }   
      
      if(opt.Contains("units")) {
	hETotalvsTime->GetZaxis()->SetTitle("E [GV/m]");
	hETotalvsTime->GetYaxis()->SetTitle("#zeta [#mum]");
	hETotalvsTime->GetXaxis()->SetTitle("z [mm]");
      } else {
	hETotalvsTime->GetZaxis()->SetTitle("E/E_{0}");
	hETotalvsTime->GetYaxis()->SetTitle("k_{p} #zeta");
	hETotalvsTime->GetXaxis()->SetTitle("k_{p} z");
      }

      hETotalvsTime->GetZaxis()->CenterTitle();
      hETotalvsTime->GetYaxis()->CenterTitle();
      hETotalvsTime->GetXaxis()->CenterTitle();
      hETotalvsTime->SetName(hName);

      hETotalvsTime->Write(hName,TObject::kOverwrite);

      // Ion probability
      TH2F *hIonProbvsTime[NAtoms] = {NULL};
      for(Int_t iat=0; iat<NAtoms; iat++) {
	sprintf(hName,"hIonProbvsTime_%s",atNames[iat]);
	TH2F *hIonProbvsTimeOld = (TH2F*) loopfile->Get(hName);
	if(hIonProbvsTimeOld) {
	  nBins = hIonProbvsTimeOld->GetNbinsX()+1;
	  Float_t binwidth =  (Time - hIonProbvsTimeOld->GetXaxis()->GetBinCenter(1))/(nBins-1);
	  edge0 = hIonProbvsTimeOld->GetXaxis()->GetBinCenter(1) - binwidth/2.;
	  edge1 = Time + binwidth/2.;
	}
	hIonProbvsTime[iat] = new TH2F("temp","",nBins,edge0,edge1,
				       hIonProb1D[iat]->GetNbinsX(),
				       hIonProb1D[iat]->GetBinLowEdge(1),
				       hIonProb1D[iat]->GetBinLowEdge(hIonProb1D[iat]->GetNbinsX()+1));
	
	for(Int_t ix=1;ix<hIonProbvsTime[iat]->GetNbinsX();ix++) {
	  for(Int_t iy=1;iy<hIonProbvsTime[iat]->GetNbinsY();iy++) {
	    hIonProbvsTime[iat]->SetBinContent(ix,iy,hIonProbvsTimeOld->GetBinContent(ix,iy));
	  }
	}  
	delete hIonProbvsTimeOld;
      
	// Fill last bin with the newest values.
	for(Int_t iy=1;iy<=hIonProb1D[iat]->GetNbinsX();iy++) {
	  hIonProbvsTime[iat]->SetBinContent(nBins,iy,hIonProb1D[iat]->GetBinContent(iy));
	}   
	
	if(opt.Contains("units")) {
	  char axName[24];
	  sprintf(axName,"#Gamma_{%s} [fs^{-1}]",atAxNames[iat]);
	  hIonProbvsTime[iat]->GetZaxis()->SetTitle(axName);
	  hIonProbvsTime[iat]->GetYaxis()->SetTitle("#zeta [#mum]");
	  hIonProbvsTime[iat]->GetXaxis()->SetTitle("z [mm]");
	} else {
	  char axName[24];
	  sprintf(axName,"#Gamma_{%s}/#omega_{p}",atAxNames[iat]);
	  hIonProbvsTime[iat]->GetZaxis()->SetTitle(axName);
	  hIonProbvsTime[iat]->GetYaxis()->SetTitle("k_{p} #zeta");
	  hIonProbvsTime[iat]->GetXaxis()->SetTitle("k_{p} z");
	}
	
	hIonProbvsTime[iat]->GetZaxis()->CenterTitle();
	hIonProbvsTime[iat]->GetYaxis()->CenterTitle();
	hIonProbvsTime[iat]->GetXaxis()->CenterTitle();
	hIonProbvsTime[iat]->SetName(hName);

	hIonProbvsTime[iat]->Write(hName,TObject::kOverwrite);
      }


      // Potential
      TH2F *hVvsTime = NULL;
      sprintf(hName,"hVvsTime");
      TH2F *hVvsTimeOld = (TH2F*) loopfile->Get(hName);
      if(hVvsTimeOld) {
	nBins = hVvsTimeOld->GetNbinsX()+1;
	Float_t binwidth =  (Time - hVvsTimeOld->GetXaxis()->GetBinCenter(1))/(nBins-1);
	edge0 = hVvsTimeOld->GetXaxis()->GetBinCenter(1) - binwidth/2.;
	edge1 = Time + binwidth/2.;
      }
      hVvsTime = new TH2F("temp","",nBins,edge0,edge1,
			  hV1D->GetNbinsX(),
			  hV1D->GetBinLowEdge(1),
			  hV1D->GetBinLowEdge(hV1D->GetNbinsX()+1));
      
      for(Int_t ix=1;ix<hVvsTime->GetNbinsX();ix++) {
	for(Int_t iy=1;iy<hVvsTime->GetNbinsY();iy++) {
	  hVvsTime->SetBinContent(ix,iy,hVvsTimeOld->GetBinContent(ix,iy));
	}
      }  
      delete hVvsTimeOld;
      
      // Fill last bin with the newest values.
      for(Int_t iy=1;iy<=hV1D->GetNbinsX();iy++) {
	hVvsTime->SetBinContent(nBins,iy,hV1D->GetBinContent(iy));
      }   
      
      if(opt.Contains("units")) {
	hVvsTime->GetZaxis()->SetTitle("#Psi-#Psi_{T} [MV]");
	hVvsTime->GetYaxis()->SetTitle("#zeta [#mum]");
	hVvsTime->GetXaxis()->SetTitle("z [mm]");
      } else {
	hVvsTime->GetZaxis()->SetTitle("#psi-#psi_{T}");
	hVvsTime->GetYaxis()->SetTitle("k_{p} #zeta");
	hVvsTime->GetXaxis()->SetTitle("k_{p} z");
      }

      hVvsTime->GetZaxis()->CenterTitle();
      hVvsTime->GetYaxis()->CenterTitle();
      hVvsTime->GetXaxis()->CenterTitle();
      hVvsTime->SetName(hName);

      hVvsTime->Write(hName,TObject::kOverwrite);


      // Focusing field
      TH2F *hFocusvsTime = NULL;
      sprintf(hName,"hFocusvsTime");
      TH2F *hFocusvsTimeOld = (TH2F*) loopfile->Get(hName);
      if(hFocusvsTimeOld) {
	nBins = hFocusvsTimeOld->GetNbinsX()+1;
	Float_t binwidth =  (Time - hFocusvsTimeOld->GetXaxis()->GetBinCenter(1))/(nBins-1);
	edge0 = hFocusvsTimeOld->GetXaxis()->GetBinCenter(1) - binwidth/2.;
	edge1 = Time + binwidth/2.;
      }
      hFocusvsTime = new TH2F("temp","",nBins,edge0,edge1,
			      hFocus1D->GetNbinsX(),
			      hFocus1D->GetBinLowEdge(1),
			      hFocus1D->GetBinLowEdge(hFocus1D->GetNbinsX()+1));
      
      for(Int_t ix=1;ix<hFocusvsTime->GetNbinsX();ix++) {
	for(Int_t iy=1;iy<hFocusvsTime->GetNbinsY();iy++) {
	  hFocusvsTime->SetBinContent(ix,iy,hFocusvsTimeOld->GetBinContent(ix,iy));
	}
      }  
      delete hFocusvsTimeOld;
      
      // Fill last bin with the newest values.
      for(Int_t iy=1;iy<=hFocus1D->GetNbinsX();iy++) {
	hFocusvsTime->SetBinContent(nBins,iy,hFocus1D->GetBinContent(iy));
      }   
      
      if(opt.Contains("units")) {
	hFocusvsTime->GetZaxis()->SetTitle("E_{x}-B_{y} [GV/m]");
	hFocusvsTime->GetYaxis()->SetTitle("#zeta [#mum]");
	hFocusvsTime->GetXaxis()->SetTitle("z [mm]");
      } else {
	hFocusvsTime->GetZaxis()->SetTitle("E_{x}-B_{y}/E_{0}");
	hFocusvsTime->GetYaxis()->SetTitle("k_{p} #zeta");
	hFocusvsTime->GetXaxis()->SetTitle("k_{p} z");
      }

      hFocusvsTime->GetZaxis()->CenterTitle();
      hFocusvsTime->GetYaxis()->CenterTitle();
      hFocusvsTime->GetXaxis()->CenterTitle();
      hFocusvsTime->SetName(hName);

      hFocusvsTime->Write(hName,TObject::kOverwrite);

      
      // On-axis beam density vs \zeta vs time! _________________________________
      TH2F **hDen1DvsTime = new TH2F*[Nspecies];
      for(Int_t i=0;i<Nspecies;i++) {
	
	if(hDen1D[i]) {

	  char hName[24];
	  sprintf(hName,"hDenvsTime_%i",i);
	  TH2F *hDen1DvsTimeOld = (TH2F*) loopfile->Get(hName);
	  
	  Int_t nBins   = 1;
	  Float_t edge0 = Time-0.5;
	  Float_t edge1 = Time+0.5;
	  if(hDen1DvsTimeOld!=NULL) {
	    nBins = hDen1DvsTimeOld->GetNbinsX()+1;
	    Float_t binwidth =  (Time - hDen1DvsTimeOld->GetXaxis()->GetBinCenter(1))/(nBins-1);
	    edge0 = hDen1DvsTimeOld->GetXaxis()->GetBinCenter(1) - binwidth/2.;
	    edge1 = Time + binwidth/2.;
	  }
	  hDen1DvsTime[i] = new TH2F("temp","",nBins,edge0,edge1,
				     hDen1D[i]->GetNbinsX(),
				     hDen1D[i]->GetBinLowEdge(1),
				     hDen1D[i]->GetBinLowEdge(hDen1D[i]->GetNbinsX()+1));
	  
	  for(Int_t ix=1;ix<hDen1DvsTime[i]->GetNbinsX();ix++) {
	    for(Int_t iy=1;iy<hDen1DvsTime[i]->GetNbinsY();iy++) {
	      hDen1DvsTime[i]->SetBinContent(ix,iy,hDen1DvsTimeOld->GetBinContent(ix,iy));
	    }
	  }  
	  delete hDen1DvsTimeOld;
	  
	  // Fill last bin with the newest values.
	  for(Int_t iy=1;iy<=hDen1D[i]->GetNbinsX();iy++) {
	    hDen1DvsTime[i]->SetBinContent(nBins,iy,hDen1D[i]->GetBinContent(iy));
	  }   
	  
	  hDen1DvsTime[i]->GetZaxis()->SetTitle("n_{b}/n_{0}");
	  if(opt.Contains("units")) {
	    hDen1DvsTime[i]->GetYaxis()->SetTitle("#zeta [#mum]");
	    hDen1DvsTime[i]->GetXaxis()->SetTitle("z [mm]"); 
	  } else {
	    hDen1DvsTime[i]->GetYaxis()->SetTitle("k_{p} #zeta");
	    hDen1DvsTime[i]->GetXaxis()->SetTitle("k_{p} z"); 
	  }
	  hDen1DvsTime[i]->GetZaxis()->CenterTitle();
	  hDen1DvsTime[i]->GetYaxis()->CenterTitle();
	  hDen1DvsTime[i]->GetXaxis()->CenterTitle();
	  hDen1DvsTime[i]->SetName(hName);
	
	  // Change the range of z axis 
	  Float_t Denmax = hDen1DvsTime[i]->GetMaximum();
	  hDen1DvsTime[i]->GetZaxis()->SetRangeUser(0,Denmax); 
	  hDen1DvsTime[i]->Write(hName,TObject::kOverwrite);
	}
      }

      // Transverse charge RMS vs \zeta vs time! _________________________________
      TH2F *hRmsvsTime = NULL; 
      if(hRms) {
	char hName[24];
	sprintf(hName,"hRmsvsTime_1");
	TH2F *hRmsvsTimeOld = (TH2F*) loopfile->Get(hName);

	Int_t nBins   = 1;
	Float_t edge0 = Time-0.5;
	Float_t edge1 = Time+0.5;
	if(hRmsvsTimeOld!=NULL) {
	  nBins = hRmsvsTimeOld->GetNbinsX()+1;
	  Float_t binwidth =  (Time - hRmsvsTimeOld->GetXaxis()->GetBinCenter(1))/(nBins-1);
	  edge0 = hRmsvsTimeOld->GetXaxis()->GetBinCenter(1) - binwidth/2.;
	  edge1 = Time + binwidth/2.;
	}
	hRmsvsTime = new TH2F("temp","",nBins,edge0,edge1,
			      hRms->GetNbinsX(),
			      hRms->GetBinLowEdge(1),
			      hRms->GetBinLowEdge(hRms->GetNbinsX()+1));
	
	for(Int_t ix=1;ix<hRmsvsTime->GetNbinsX();ix++) {
	  for(Int_t iy=1;iy<hRmsvsTime->GetNbinsY();iy++) {
	    hRmsvsTime->SetBinContent(ix,iy,hRmsvsTimeOld->GetBinContent(ix,iy));
	  }
	}  
	delete hRmsvsTimeOld;
	
	// Fill last bin with the newest values.
	for(Int_t iy=1;iy<=hRms->GetNbinsX();iy++) {
	  hRmsvsTime->SetBinContent(nBins,iy,hRms->GetBinContent(iy));
	}   
	
	if(opt.Contains("units")) {
	  hRmsvsTime->GetYaxis()->SetTitle("#zeta [#mum]");
	  hRmsvsTime->GetXaxis()->SetTitle("z [mm]"); 
	  hRmsvsTime->GetZaxis()->SetTitle("#sigma_{x} [#mum]");
	} else {
	  hRmsvsTime->GetYaxis()->SetTitle("k_{p} #zeta");
	  hRmsvsTime->GetXaxis()->SetTitle("k_{p} z"); 
	  hRmsvsTime->GetZaxis()->SetTitle("k_{p} #sigma_{x}");
	}

	hRmsvsTime->GetZaxis()->CenterTitle();
	hRmsvsTime->GetYaxis()->CenterTitle();
	hRmsvsTime->GetXaxis()->CenterTitle();
	hRmsvsTime->SetName(hName);

	// Change the range of z axis
	Float_t Rmsmax = hRmsvsTime->GetMaximum();
	hRmsvsTime->GetZaxis()->SetRangeUser(0,Rmsmax); 
	hRmsvsTime->Write(hName,TObject::kOverwrite);
      }
      
      //--
      
      
      loopfile->Close();
    }
  
    cout << Form("\n End. (%i) -----------\n",time) << endl;

  }
  

}




