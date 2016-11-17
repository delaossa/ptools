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
  Float_t period = pData->GetPlasmaPeriod();
  Float_t timedepth = pData->GetPlasmaTimeDepth();
  Float_t kp = pData->GetPlasmaK();
  Float_t skindepth = pData->GetPlasmaSkinDepth();
  Float_t wavelength = pData->GetPlasmaWaveLength();
  Float_t E0 = pData->GetPlasmaE0();

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

    cout << Form(" --- doEvolutions ---\n") << endl;
  
    cout << Form(" Plotting range:  %.2f < x1 < %.2f ,  %.2f < x2 < %.2f",zMin,zMax,yMin,yMax) << endl;

    // ----------------------------------------------------------------------------------
      
    // Get charge density histos
    Int_t Nspecies = pData->NSpecies();
    TH2F **hDen2D = new TH2F*[Nspecies];
    // Get charge density on-axis
    TH1F **hDen1D = new TH1F*[Nspecies];
    // And electric current (integrated)
    TH1F **hCur1D = new TH1F*[Nspecies];
    // Get Jz histos
    TH1F **hJz1D = new TH1F*[Nspecies];
    for(Int_t i=0;i<Nspecies;i++) {

      hDen2D[i] = NULL;
      hDen1D[i] = NULL;
      hCur1D[i] = NULL;
      hJz1D[i] = NULL;
      
      if(!pData->GetChargeFileName(i)) 
	continue;
      
      cout << Form(" -> Getting charge density of specie %i  (%s)", i, pData->GetSpeciesName(i).c_str())  << endl;

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

      // on-axis current density
      if(pData->GetCurrentFileName(i)) {
	sprintf(hName,"hJz1D_%i",i);
	hJz1D[i] = (TH1F*) gROOT->FindObject(hName);
	if(hJz1D[i]) delete hJz1D[i];
	
	// 1D histograms
	if(pData->Is3D()) {
	  hJz1D[i] = pData->GetH1SliceZ3D(pData->GetCurrentFileName(i)->c_str(),"j1",-1,NonBin,-1,NonBin,opt+"avg");
	} else if(pData->IsCyl()) { // Cylindrical: The first bin with r>0 is actually the number 1 (not the 0).
	  hJz1D[i] = pData->GetH1SliceZ(pData->GetCurrentFileName(i)->c_str(),"j1",1,NonBin,opt+"avg");
	} else { // 2D cartesian
	  hJz1D[i] = pData->GetH1SliceZ(pData->GetCurrentFileName(i)->c_str(),"j1",-1,NonBin,opt+"avg");
	}
	hJz1D[i]->SetName(hName); 
	
	// if(hJz1D[i]) delete hJz1D[i];
	// hJz1D[i] = (TH1F*) hE2D[i]->ProjectionX(hName,FirstxBin,LastxBin);
	// hJz1D[i]->Scale(1.0/(LastxBin-FirstxBin+1));
	
	if(opt.Contains("comov"))
	  hJz1D[i]->GetXaxis()->SetTitle("k_{p} #zeta");
	else
	  hJz1D[i]->GetXaxis()->SetTitle("k_{p} z");
	
	hJz1D[i]->GetYaxis()->SetTitle(Form("j_{z,%i}",i));
      }

      
      // Only for the beams here on
      if(i==0 && pData->GetSpeciesName(i).find("beam") == string::npos) continue;

      // Get the current
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

      cout << Form(" -> Getting electric field number ") << i+1 << endl;
    
      char hName[24];
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
	
        if(i==0) 
	  hE1D[i] = pData->GetH1SliceZ3D(pData->GetEfieldFileName(i)->c_str(),nam,-1,NonBin,-1,NonBin,opt+"avg");
	else // In case of transverse fields, the 1D line is taken off-axis 
	  hE1D[i] = pData->GetH1SliceZ3D(pData->GetEfieldFileName(i)->c_str(),nam,-NofBin,NonBin,-1,NonBin,opt+"avg");
	
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

    // -> Get magnetic fields 1D
    const Int_t NBfields = 3;
    TH1F **hB1D = new TH1F*[NBfields];
    for(Int_t i=0;i<NBfields;i++) {
      hB1D[i] = NULL;
      if(i<2) continue;  // Just get the third component.
      
      if(!pData->GetBfieldFileName(i))
	continue;
      
      cout << Form(" -> Getting magnetic field number ") << i+1 << endl;
      
      char hName[24];
      sprintf(hName,"hB1D_%i",i);
      hB1D[i] = (TH1F*) gROOT->FindObject(hName);
      if(hB1D[i]) delete hB1D[i];
      
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
       
      // 1D histograms
      if(pData->Is3D()) {
	if(i==0) 
	  hB1D[i] = pData->GetH1SliceZ3D(pData->GetBfieldFileName(i)->c_str(),nam,-1,NonBin,-1,NonBin,opt+"avg");
	else // In case of transverse fields, the 1D line is taken offaxis
	  hB1D[i] = pData->GetH1SliceZ3D(pData->GetBfieldFileName(i)->c_str(),nam,-NofBin,NonBin,-1,NonBin,opt+"avg");
	
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
      Emax = 1.0;
      PUnits::BestUnit beSUnit(E0 * Emax,"Efield");
      beSUnit.GetBestUnits(eUnit,eSUnit);
      cout << Form(" E0 = %.2f %s", E0 * Emax/eUnit, eSUnit.c_str()) << endl;

      for(Int_t i=0;i<Nspecies;i++) {

	if(!hDen1D[i]) continue;
    
	Int_t NbinsX = hDen1D[i]->GetNbinsX();
	Double_t zMin = skindepth * hDen1D[i]->GetXaxis()->GetXmin() / spaUnit;
	Double_t zMax = skindepth * hDen1D[i]->GetXaxis()->GetXmax() / spaUnit;
	hDen1D[i]->SetBins(NbinsX,zMin,zMax);
      
	if(opt.Contains("comov"))
	  hDen1D[i]->GetXaxis()->SetTitle(Form("#zeta [%s]",spaSUnit.c_str()));
	else
	  hDen1D[i]->GetXaxis()->SetTitle(Form("z [%s]",spaSUnit.c_str()));

	if(hJz1D[i]) {
	  hJz1D[i]->SetBins(NbinsX,zMin,zMax);
	  if(opt.Contains("comov"))
	    hJz1D[i]->GetXaxis()->SetTitle(Form("#zeta [%s]",spaSUnit.c_str()));
	  else
	    hJz1D[i]->GetXaxis()->SetTitle(Form("z [%s]",spaSUnit.c_str()));
	}	
      
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
	if(!hE1D[i]) continue;

	Int_t NbinsX = hE1D[i]->GetNbinsX();
	Double_t zMin = skindepth * hE1D[i]->GetXaxis()->GetXmin() / spaUnit;
	Double_t zMax = skindepth * hE1D[i]->GetXaxis()->GetXmax() / spaUnit;
	hE1D[i]->SetBins(NbinsX,zMin,zMax);
            
	for(Int_t j=0;j<=hE1D[i]->GetNbinsX();j++) {
	  hE1D[i]->SetBinContent(j, hE1D[i]->GetBinContent(j) * E0 / eUnit );
	}
      
      
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
	if(!hB1D[i]) continue;

	Int_t NbinsX = hB1D[i]->GetNbinsX();
	Double_t zMin = skindepth * hB1D[i]->GetXaxis()->GetXmin() / spaUnit;
	Double_t zMax = skindepth * hB1D[i]->GetXaxis()->GetXmax() / spaUnit;
	hB1D[i]->SetBins(NbinsX,zMin,zMax);
      
	for(Int_t j=0;j<=hB1D[i]->GetNbinsX();j++) {
	  hB1D[i]->SetBinContent(j, hB1D[i]->GetBinContent(j) * E0  / eUnit );
	}
	
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
    TH1F *hFocus1D = (TH1F*) hE1D[1]->Clone("hFocus1D");
    if(hB1D[2]) {
      hFocus1D->Add(hB1D[2],-1);
      if(opt.Contains("units")) {
	hFocus1D->GetZaxis()->SetTitle(Form("E_{x}-cB_{y} [%s]",eSUnit.c_str()));
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
	hRms->GetYaxis()->SetTitle(Form("#Deltay_{rms} [%s]",spaSUnit.c_str()));
      } else {
	hRms->GetYaxis()->SetTitle("k_{p}#Deltay_{rms}");
      }
      	
    }

    // Now, combine the electric field components into the total |E|
    TH1F *hETotal1D = NULL;
    if(hE1D[0] && hE1D[1] && hE1D[2])
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
    if(hE1D[0]) {
      Int_t   NbinsX = hE1D[0]->GetNbinsX();
      Float_t dx = pData->GetDX(0);
      
      char hName[24];
      sprintf(hName,"hV1D");
      hV1D = (TH1F*) hE1D[0]->Clone(hName);
      hV1D->Reset();
      
      Float_t integral = 0.0;
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
    
    if(opt.Contains("units")) {
      hV1D->GetYaxis()->SetTitle("#Psi [MV]");
    } else {
      hV1D->GetYaxis()->SetTitle("#psi");
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

      Double_t zUnits = PUnits::mm;
      string zSUnit = "mm";
      if(opt.Contains("units")) Time *= skindepth / zUnits;

      Int_t nBins   = 1;
      Float_t edge0 = Time-0.5;
      Float_t edge1 = Time+0.5;
      char hName[24];

      TH2F **hEvsTime = new TH2F*[Nfields]; 
      for(Int_t i=0;i<Nfields;i++) {
	if(!hE1D[i]) continue;
	
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
	    hEvsTime[i]->GetZaxis()->SetTitle(Form("E_{z} [%s]",eSUnit.c_str()));
	  else if(i==1)
	    hEvsTime[i]->GetZaxis()->SetTitle(Form("E_{x} [%s]",eSUnit.c_str()));
	  else if(i==2)
	    hEvsTime[i]->GetZaxis()->SetTitle(Form("E_{y} [%s]",eSUnit.c_str()));
	  
	  hEvsTime[i]->GetYaxis()->SetTitle(Form("#zeta [%s]",spaSUnit.c_str()));
	  hEvsTime[i]->GetXaxis()->SetTitle(Form("z [%s]",zSUnit.c_str()));
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
      if(hETotal1D) {
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
	  hETotalvsTime->GetZaxis()->SetTitle(Form("E [%s]",eSUnit.c_str()));
	  hETotalvsTime->GetYaxis()->SetTitle(Form("#zeta [%s]",spaSUnit.c_str()));
	  hETotalvsTime->GetXaxis()->SetTitle(Form("z [%s]",zSUnit.c_str()));
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
      }
      
      // Potential
      if(hV1D) {
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
	  hVvsTime->GetZaxis()->SetTitle("#Psi [MV]");
	  hVvsTime->GetYaxis()->SetTitle(Form("#zeta [%s]",spaSUnit.c_str()));
	  hVvsTime->GetXaxis()->SetTitle(Form("z [%s]",zSUnit.c_str()));
	} else {
	  hVvsTime->GetZaxis()->SetTitle("#psi");
	  hVvsTime->GetYaxis()->SetTitle("k_{p} #zeta");
	  hVvsTime->GetXaxis()->SetTitle("k_{p} z");
	}

	hVvsTime->GetZaxis()->CenterTitle();
	hVvsTime->GetYaxis()->CenterTitle();
	hVvsTime->GetXaxis()->CenterTitle();
	hVvsTime->SetName(hName);

	hVvsTime->Write(hName,TObject::kOverwrite);
      }

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
	hFocusvsTime->GetZaxis()->SetTitle(Form("E_{x}-B_{y} [%s]",eSUnit.c_str()));
	hFocusvsTime->GetYaxis()->SetTitle(Form("#zeta [%s]",spaSUnit.c_str()));
	hFocusvsTime->GetXaxis()->SetTitle(Form("z [%s]",zSUnit.c_str()));
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

      
      // On-axis density vs \zeta vs time! _________________________________
      TH2F **hDen1DvsTime = new TH2F*[Nspecies];
      TH2F **hJz1DvsTime = new TH2F*[Nspecies];
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
	    hDen1DvsTime[i]->GetYaxis()->SetTitle(Form("#zeta [%s]",spaSUnit.c_str()));
	    hDen1DvsTime[i]->GetXaxis()->SetTitle(Form("z [%s]",zSUnit.c_str()));
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

	if(hJz1D[i]) {

	  char hName[24];
	  sprintf(hName,"hJzvsTime_%i",i);
	  TH2F *hJz1DvsTimeOld = (TH2F*) loopfile->Get(hName);
	  
	  Int_t nBins   = 1;
	  Float_t edge0 = Time-0.5;
	  Float_t edge1 = Time+0.5;
	  if(hJz1DvsTimeOld!=NULL) {
	    nBins = hJz1DvsTimeOld->GetNbinsX()+1;
	    Float_t binwidth =  (Time - hJz1DvsTimeOld->GetXaxis()->GetBinCenter(1))/(nBins-1);
	    edge0 = hJz1DvsTimeOld->GetXaxis()->GetBinCenter(1) - binwidth/2.;
	    edge1 = Time + binwidth/2.;
	  }
	  hJz1DvsTime[i] = new TH2F("temp","",nBins,edge0,edge1,
				    hJz1D[i]->GetNbinsX(),
				    hJz1D[i]->GetBinLowEdge(1),
				    hJz1D[i]->GetBinLowEdge(hJz1D[i]->GetNbinsX()+1));
	  
	  for(Int_t ix=1;ix<hJz1DvsTime[i]->GetNbinsX();ix++) {
	    for(Int_t iy=1;iy<hJz1DvsTime[i]->GetNbinsY();iy++) {
	      hJz1DvsTime[i]->SetBinContent(ix,iy,hJz1DvsTimeOld->GetBinContent(ix,iy));
	    }
	  }  
	  delete hJz1DvsTimeOld;
	  
	  // Fill last bin with the newest values.
	  for(Int_t iy=1;iy<=hJz1D[i]->GetNbinsX();iy++) {
	    hJz1DvsTime[i]->SetBinContent(nBins,iy,hJz1D[i]->GetBinContent(iy));
	  }   
	  
	  hJz1DvsTime[i]->GetZaxis()->SetTitle(Form("j_{z,%i}",i));
	  if(opt.Contains("units")) {
	    hJz1DvsTime[i]->GetYaxis()->SetTitle(Form("#zeta [%s]",spaSUnit.c_str()));
	    hJz1DvsTime[i]->GetXaxis()->SetTitle(Form("z [%s]",zSUnit.c_str()));
	  } else {
	    hJz1DvsTime[i]->GetYaxis()->SetTitle("k_{p} #zeta");
	    hJz1DvsTime[i]->GetXaxis()->SetTitle("k_{p} z"); 
	  }
	  hJz1DvsTime[i]->GetZaxis()->CenterTitle();
	  hJz1DvsTime[i]->GetYaxis()->CenterTitle();
	  hJz1DvsTime[i]->GetXaxis()->CenterTitle();
	  hJz1DvsTime[i]->SetName(hName);
	
	  // Change the range of z axis 
	  Float_t Jzmax = hJz1DvsTime[i]->GetMaximum();
	  Float_t Jzmin = hJz1DvsTime[i]->GetMinimum();

	  // if(Jzmax >  TMath::Abs(Jzmin))
	  //   Jzmin = -Jzmax;
	  // else
	  //   Jzmax = -Jzmin;
	  
	  hJz1DvsTime[i]->GetZaxis()->SetRangeUser(Jzmin,Jzmax); 
	  hJz1DvsTime[i]->Write(hName,TObject::kOverwrite);
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
	  hRmsvsTime->GetYaxis()->SetTitle(Form("#zeta [%s]",spaSUnit.c_str()));
	  hRmsvsTime->GetXaxis()->SetTitle(Form("z [%s]",zSUnit.c_str()));
	  hRmsvsTime->GetZaxis()->SetTitle(Form("#sigma_{x} [%s]",spaSUnit.c_str()));
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




