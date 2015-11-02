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
#include <TEllipse.h>

#include "PData.hh"
#include "PDataHiP.hh"
#include "PGlobals.hh"
#include "PPalette.hh"
#include "H5Cpp.h"

using namespace std;
using namespace H5;

void FindLimits(TH1F *h,Double_t &xmin, Double_t &xmax, Double_t factor = 0.11) {
  Double_t maxValue = h->GetBinContent(h->GetMaximumBin());
  Int_t   lBin = -1;
  for(Int_t i=1;i<=h->GetNbinsX();i++) {
    Double_t binValue = h->GetBinContent(i);
    if(binValue>maxValue*factor) {
      lBin = i;
      xmin = h->GetBinCenter(i);
      break;
    }
  }
  
  Int_t rBin = -1;	
  for(Int_t i=h->GetNbinsX();i>0;i--) {
    Double_t binValue = h->GetBinContent(i);
    if(binValue>maxValue*factor) {
      rBin = i;
      xmax = h->GetBinCenter(i);
      break;
    }
  } 
}

int main(int argc,char *argv[]) {
  if(argc<=2) {
    printf("\n Usage: %s <simulation name> <-t(time)>\n",argv[0]);
    printf("      <-index(species index)>\n");
    printf("      <-i(initial time)> <-f(final time)> <-s(time step)>\n");
    printf("      <--png> <--pdf> <--eps> <--slcs>\n");
    printf("      <--center> <--comov>\n");
    printf("      <--units> <--logz>\n");
    printf("      <--file> <--loop>\n");
    return 0;
  }

  // General options
  TString   sim = "";
  Int_t    time = 0;
  Int_t  iStart = -1;
  Int_t    iEnd = -1;
  Int_t   iStep = 1;
  TString   opt = "";
  Int_t   index = 1;

  // Option for raw fraction correction
  Float_t rawf = 1;
  
  // Options for Spectrum
  Float_t Pmin =  99999.;
  Float_t Pmax = -99999.;

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
    } else if(arg.Contains("--tiff")) {
      opt += "tiff";
    } else if(arg.Contains("--slcs")) {
      opt += "slcs";
    } else if(arg.Contains("--comov")){
      opt += "comov";
    } else if(arg.Contains("--units")){
      opt += "units";
    } else if(arg.Contains("--center")){
      opt += "center";
    } else if(arg.Contains("--best")){
      opt += "best";
    } else if(arg.Contains("--grid")){
      opt += "grid"; 
    } else if(arg.Contains("--logz")){
      opt += "logz"; 
    } else if(arg.Contains("--autop")){
      opt += "autop"; 
    } else if(arg.Contains("--auto")){
      opt += "auto"; 
    } else if(arg.Contains("--loop")){
      opt += "loop"; 
    } else if(arg.Contains("--file")){
      opt += "file"; 
    } else if(arg.Contains("--trans")){
      opt += "trans"; 
    } else if(arg.Contains("--notext")){
      opt += "notext"; 
    } else if(arg.Contains("--noinfo")){
      opt += "noinfo"; 
    } else if(arg.Contains("--bw")){
      opt += "bw"; 
    } else if(arg.Contains("--nospec")){
      opt += "nospec"; 
    } else if(arg.Contains("-index")) {
      char ss[6];
      sscanf(arg,"%6s%i",ss,&index);
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
    } else if(arg.Contains("-rawf")) {
      char ss[5];
      sscanf(arg,"%5s%f",ss,&rawf);
    } else if(arg.Contains("-pmin")) {
      char ss[5];
      sscanf(arg,"%5s%f",ss,&Pmin);
    } else if(arg.Contains("-pmax")) {
      char ss[5];
      sscanf(arg,"%5s%f",ss,&Pmax);
    } else {
      cout << Form("\t Invalid argument (%i): exiting...\n",l) << endl;
      return 0;
    }
  }
  
  
  PGlobals::Initialize();

  // Palettes!
  gROOT->Macro("PPalettes.C");

  if(opt.Contains("grid")) {
    gStyle->SetPadGridX(1);
    gStyle->SetPadGridY(1);
  }

  // Load PData
  PData *pData = PData::Get(sim.Data());
  if(pData->isHiPACE()) {
    delete pData; pData = NULL;
    pData = PDataHiP::Get(sim.Data());
    opt += "comov";
  }
  
  if(iStart<0) iStart = time;
  if(iEnd<=iStart) iEnd = iStart;
  
  // Some plasma constants
  Double_t n0 = pData->GetPlasmaDensity();
  Double_t kp = pData->GetPlasmaK();
  Double_t skindepth = pData->GetPlasmaSkinDepth();
  Double_t wavelength = pData->GetPlasmaWaveLength();

  // z start of the plasma in normalized units.
  Double_t zStartPlasma = pData->GetPlasmaStart()*kp;
  // z start of the beam in normalized units.
  Double_t zStartBeam = pData->GetBeamStart()*kp;

  //  opt += "comovcenter";

  // Units
  Double_t chargeUnit,  denUnit,  spaUnit,  tspaUnit,  eneUnit,  teneUnit,  curUnit,   emitUnit,   ermsUnit,  betaUnit;
  string  chargeSUnit, denSUnit, spaSUnit, tspaSUnit, eneSUnit, teneSUnit, curSUnit,  emitSUnit,  ermsSUnit,  betaSUnit;

  if(opt.Contains("units") && n0) {

    //  cout << endl << Form("Changing to SI units:") << endl;
      
    // Get the best units for each quantity
    PUnits::BestUnit bdenSUnit(n0,"PartDensity");
    bdenSUnit.GetBestUnits(denUnit,denSUnit);
    //    cout << Form(" n0 = %.2f %s", n0/denUnit, denSUnit.c_str()) << endl;

    PUnits::BestUnit bspaSUnit(wavelength,"Length");
    bspaSUnit.GetBestUnits(spaUnit,spaSUnit);
    //    cout << Form(" L  = %.2f %s", wavelength/spaUnit, spaSUnit.c_str()) << endl;
  }

  // Current units
  curUnit = PUnits::kA;
  curSUnit = "kA";
  
  // Charge units
  chargeUnit = PUnits::picocoulomb;
  chargeSUnit = "pC";
  
  // Longitudinal momentum units
  eneUnit = PUnits::GeV;
  eneSUnit = "GeV";

  // Transverse momentum units
  teneUnit = PUnits::MeV;
  teneSUnit = "MeV";

  // Transverse length units
  tspaUnit = PUnits::um;
  tspaSUnit = "#mum";
  
  // Emittance units
  emitUnit = PUnits::um;
  emitSUnit = "#mum";
  // emitUnit = 10 * PUnits::um;
  // emitSUnit = "10 #mum";

  // Relative energy spread units
  ermsUnit = PUnits::perCent; 
  ermsSUnit = "%";

  
  // Time looper
  for(Int_t i=iStart; i<iEnd+1; i+=iStep) {

    time = i;
    pData->LoadFileNames(time);    
    //if(time==iStart) pData->PrintData();
    
    if(!pData->IsInit()) continue;

    Int_t Nspecies = pData->NSpecies();
    if(index>Nspecies-1) {
      return 0;
    }
    if(!pData->GetRawFileName(index)) {
      return 0;    
    }

    // Time in OU
    Double_t Time = pData->GetRealTime();
    Time += pData->ShiftT(opt);
    
    // Centering time and z position:
    Double_t shiftz = pData->Shift(opt);

    // BOX limits
    Double_t X1MIN = pData->GetXMin(0);
    Double_t X1MAX = pData->GetXMax(0);
    Double_t X2MIN = pData->GetXMin(1);
    Double_t X2MAX = pData->GetXMax(1);
    Double_t X3MIN = -999;
    Double_t X3MAX = -999;
    if(pData->Is3D()) {
      X3MIN = pData->GetXMin(2);
      X3MAX = pData->GetXMax(2);
    }
    
    // Spatial resolution
    Double_t dx1, dx2, dx3;
    
    dx1 = pData->GetDX(0);
    dx2 = pData->GetDX(1);    
    if(pData->Is3D())
      dx3 = pData->GetDX(2);
    else
      dx3 = 1.0;
    
    // Bining, intervals, labels, etc.
    Int_t x1Nbin = 200;
    Int_t p1Nbin = 200;
    Int_t x2Nbin = 200;
    Int_t p2Nbin = 200;
    Int_t x3Nbin = 200;
    Int_t p3Nbin = 200;
    
    // Slices
    Int_t SNbin = 100;
    Double_t x1BinMin = -4.5;
    Double_t x1BinMax = -4.0;

    // Spatial coordinates intervals:
    // Double_t x1Min = -7.8;
    Double_t x1Min = -7.75;
    Double_t x1Max = -7.0;
    Double_t x2Min = -0.5;
    Double_t x2Max =  0.5;
    Double_t x3Min = -0.5;
    Double_t x3Max =  0.5;

    // Momentum coordinates intervals:
    Double_t p1Min =  Pmin;
    Double_t p1Max =  Pmax;
    Double_t p2Min = -15.0;
    Double_t p2Max =  15.0;
    Double_t p3Min = -15.0;
    Double_t p3Max =  15.0;

    // Specific initializations:
    // The ranges in x1 are defined in the comoving frame, centered at the driver's position.
    // We setup a dummy shift to make them consitent with the plotting options.
    if(sim.Contains("flash")) {
    
      x1Min = -5.5;
      x1Max = -3.2; 

      SNbin = 50;
      x1BinMin = -4.80;
      x1BinMax = -3.70;
    
      if(sim.Contains("v8.0kA.C.DDR")) {

	x1Min = -5.69;
	x1Max = -2.81; 

	SNbin = 50;
	x1BinMin = -5.1;
	x1BinMax = -3.2;
      }

     if(sim.Contains("v7.5kA.T.BI")) {

       x1Min = -4.3;
       x1Max = -3.9;
       
       SNbin = 20;
       x1BinMin = -4.15;
       x1BinMax = -4.02;
       
     }

     if(sim.Contains("v7.5kA.G.BI")) {

	x1Min = -6.5;
	x1Max = -5.0; 

	SNbin = 20;
	x1BinMin = -5.9;
	x1BinMax = -5.3;
      }
    
     if(sim.Contains("v2.5kA.G.JG.DDR")) {

       x1Nbin = 300;
       p1Nbin = 200;
       x2Nbin = 200;
       p2Nbin = 200;
       
       x1Min = -5.2;
       x1Max = -3.2; 
       
       SNbin = 100;
       x1BinMin = -4.70;
       x1BinMax = -3.60;
      
     }
     
     if(sim.Contains("v2.5kA.G.ZH.DDR")) {
       
       x1Nbin = 300;
       p1Nbin = 200;
       x2Nbin = 200;
       p2Nbin = 200;
       
       x1Min = -5.4;
       x1Max = -2.5; 
       
       SNbin = 100;
       x1BinMin = -4.65;
       x1BinMax = -3.40;
      
     }
     

    } else if(sim.Contains("facet")) {
      
      if(sim.Contains("v23kA.G.RI")) {
	x1Min = -7.75;
	x1Max = -7.0;
	
	SNbin = 40;
	x1BinMin = -7.54;
	x1BinMax = -7.14;

	if(sim.Contains("n20")) {
	  SNbin = 37;
	  x1BinMin = -7.51;
	  x1BinMax = -7.14;

	}

      }	else if (sim.Contains("v10kA.G.DDR")) {
	
	x1Min = -7.99;
	x1Max = -3.01; 
	
	SNbin = 50;
	x1BinMin = -6.40;
	x1BinMax = -4.15;
	
      } else if (sim.Contains("DDR")) {
	
	x1Min = -5.80;
	x1Max = -5.0; 
	
	SNbin = 40;
	x1BinMin = -5.65;
	x1BinMax = -5.1;
	
      }
      
    } else if (sim.Contains("rake-v10kA.G.SR2.RI.3D")) {
      
      x1Min = -6.4;
      x1Max = -5.9; 
      
      //      SNbin = 76;
      SNbin = 38;
      x1BinMin = -6.23;
      x1BinMax = -6.05;
      
    } else if (sim.Contains("FACET_5x1016_PP20140415_webBeam")) {
      
      x1Min = -6.5;
      x1Max = -4.5; 
      
      SNbin = 40;
      x1BinMin = -6.2;
      x1BinMax = -4.9;
      
    } else if (sim.Contains("BOND-24TW")) {
      
      x1Min = -8.0;
      x1Max = -4.30; 
      
      SNbin = 100;
      x1BinMin = -6.8;
      x1BinMax = -5.7;
      
    } else if (sim.Contains("pitz")) {
      x1Min = -55.0;
      x1Max = 25.0; 
      
      SNbin = 100;
      x1BinMin = -32.0;
      x1BinMax =   2.0;

      curUnit = PUnits::ampere;
      curSUnit = "A";

      eneUnit = PUnits::MeV;
      eneSUnit = "MeV";

      ermsUnit = PUnits::perMillion; 
      ermsSUnit = "0.01 %";

    }
    
    // dummy shift
    Double_t dshiftz = pData->Shift("centercomov");
    x1Min += dshiftz;
    x1Max += dshiftz;
    x1BinMin += dshiftz;
    x1BinMax += dshiftz;

    //    cout << Form(" x1Min = %f x1Max = %f x1BinMin = %f x1BinMax = %f  ", x1Min, x1Max, x1BinMin, x1BinMax) << endl;
    
    // --------------------------------------------------
    
    cout << Form("\n 1. Reading file : ") << pData->GetRawFileName(index)->c_str() << endl;

    // Double_t **var = NULL;
    UInt_t Nvar = 7;
    if(!pData->Is3D()) Nvar = 6;
    
    Float_t **var;
    var = new Float_t*[Nvar];
    UInt_t Np = pData->GetRawArray(pData->GetRawFileName(index)->c_str(),var);

    cout << Form("    Number of particles = %i",Np) << endl;
    
    // Scan histogram
    TH1F *hScanX1 = (TH1F*) gROOT->FindObject("hScanX1");
    if(hScanX1) delete hScanX1;
    hScanX1 = new TH1F("hScanX1","",pData->GetX1N(),pData->GetX1Min(),pData->GetX1Max());
    TH1F *hScanX2 = (TH1F*) gROOT->FindObject("hScanX2");
    if(hScanX2) delete hScanX2;
    hScanX2 = new TH1F("hScanX2","",pData->GetX2N(),pData->GetX2Min(),pData->GetX2Max());
    TH1F *hScanX3 = (TH1F*) gROOT->FindObject("hScanX3");
    if(hScanX3) delete hScanX3;
    if(Nvar==7)
      hScanX3 = new TH1F("hScanX3","",pData->GetX3N(),pData->GetX3Min(),pData->GetX3Max());

    TH1F *hScanP2 = (TH1F*) gROOT->FindObject("hScanP2");
    if(hScanP2) delete hScanP2;
    hScanP2 = new TH1F("hScanP2","",p2Nbin,p2Min,p2Max);
    TH1F *hScanP3 = (TH1F*) gROOT->FindObject("hScanP3");
    if(hScanP3) delete hScanP3;
    hScanP3 = new TH1F("hScanP3","",p3Nbin,p3Min,p3Max);

    // cout << Form(" BOX (N = %i):  x1Min = %f  x1Max = %f ", pData->GetX1N(), pData->GetX1Min()-shiftz, pData->GetX1Max()-shiftz) << endl;
    
    
    // -------------------------------------------------------------------------------

    // auto ranges factor
    Double_t rfactor = 0.3;
    
    if(opt.Contains("autop")) {

      cout << Form("\n Auto ranging everything but x1...") << endl;

      Double_t MinP1 = 999999;
      Double_t MaxP1 = -999999;
      Double_t MinP2 = 999999;
      Double_t MaxP2 = -999999;
      Double_t MinP3 = 999999;
      Double_t MaxP3 = -999999;

      Double_t MinX2 = 999999;
      Double_t MaxX2 = -999999;
      Double_t MinX3 = 999999;
      Double_t MaxX3 = -999999;

      for(UInt_t i=0;i<Np;i++) {

	// if(var[4][i]-shiftz<x1BinMin || var[4][i]-shiftz>x1BinMax ) continue; 
	// if(var[4][i]-shiftz<x1Min || var[4][i]-shiftz>x1Max ) continue; 
	if(var[4][i]<x1Min || var[4][i]>x1Max ) continue; 

	if(var[0][i]<MinP1) MinP1 = var[0][i];
	if(var[0][i]>MaxP1) MaxP1 = var[0][i];
	if(var[1][i]<MinP2) MinP2 = var[1][i];
	if(var[1][i]>MaxP2) MaxP2 = var[1][i];
	if(var[2][i]<MinP3) MinP3 = var[2][i];
	if(var[2][i]>MaxP3) MaxP3 = var[2][i];
	if(var[5][i]<MinX2) MinX2 = var[5][i];
	if(var[5][i]>MaxX2) MaxX2 = var[5][i];
	if(Nvar==7) {
	  if(var[6][i]<MinX3) MinX3 = var[6][i];
	  if(var[6][i]>MaxX3) MaxX3 = var[6][i];
	}

	hScanX1->Fill(var[4][i],TMath::Abs(var[3][i]));
	hScanX2->Fill(var[5][i],TMath::Abs(var[3][i]));
	if(Nvar==7)
	  hScanX3->Fill(var[6][i],TMath::Abs(var[3][i]));
	
	hScanP2->Fill(var[1][i],TMath::Abs(var[3][i]));
	hScanP3->Fill(var[2][i],TMath::Abs(var[3][i]));
      }
      
      p1Min = MinP1 - rfactor*(MaxP1-MinP1);
      p1Max = MaxP1 + rfactor*(MaxP1-MinP1);
      p2Min = MinP2 - rfactor*(MaxP2-MinP2);
      p2Max = MaxP2 + rfactor*(MaxP2-MinP2);
      p3Min = MinP3 - rfactor*(MaxP3-MinP3);
      p3Max = MaxP3 + rfactor*(MaxP3-MinP3);

      x2Min = MinX2 - rfactor*(MaxX2-MinX2);
      x2Max = MaxX2 + rfactor*(MaxX2-MinX2);
      if(Nvar==7) {
	x3Min = MinX3 - rfactor*(MaxX3-MinX3);
	x3Max = MaxX3 + rfactor*(MaxX3-MinX3);
      }

      // Set limits using Scan histograms
      Double_t peakFactor = 0.05;
      Double_t rfactor2 = 2;

      Double_t x2min = -999;
      Double_t x2max = -999;
      FindLimits(hScanX2,x2min,x2max,peakFactor);
      x2Min = x2min - rfactor2*(x2max-x2min);
      x2Max = x2max + rfactor2*(x2max-x2min);

      Double_t x3min = -999;
      Double_t x3max = -999;
      if(Nvar==7)
	FindLimits(hScanX3,x3min,x3max,peakFactor);
      x3Min = x3min - rfactor2*(x3max-x3min);
      x3Max = x3max + rfactor2*(x3max-x3min);

      Double_t p2min = -999;
      Double_t p2max = -999;
      FindLimits(hScanP2,p2min,p2max,peakFactor);
      p2Min = p2min - rfactor2*(p2max-p2min);
      p2Max = p2max + rfactor2*(p2max-p2min);

      Double_t p3min = -999;
      Double_t p3max = -999;
      FindLimits(hScanP3,p3min,p3max,peakFactor);
      p3Min = p3min - rfactor2*(p3max-p3min);
      p3Max = p3max + rfactor2*(p3max-p3min);

      
    } else if(opt.Contains("auto")) {

      cout << Form("\n Auto ranging...") << endl;
      
      Double_t MinP1 = 999999;
      Double_t MaxP1 = -999999;
      Double_t MinP2 = 999999;
      Double_t MaxP2 = -999999;
      Double_t MinP3 = 999999;
      Double_t MaxP3 = -999999;

      Double_t MinX1 = 999999;
      Double_t MaxX1 = -999999;
      Double_t MinX2 = 999999;
      Double_t MaxX2 = -999999;
      Double_t MinX3 = 999999;
      Double_t MaxX3 = -999999;

      for(UInt_t i=0;i<Np;i++) {
	if(var[0][i]<MinP1) MinP1 = var[0][i];
	if(var[0][i]>MaxP1) MaxP1 = var[0][i];
	if(var[1][i]<MinP2) MinP2 = var[1][i];
	if(var[1][i]>MaxP2) MaxP2 = var[1][i];
	if(var[2][i]<MinP3) MinP3 = var[2][i];
	if(var[2][i]>MaxP3) MaxP3 = var[2][i];
	if(var[4][i]<MinX1) MinX1 = var[4][i];
	if(var[4][i]>MaxX1) MaxX1 = var[4][i];
	if(var[5][i]<MinX2) MinX2 = var[5][i];
	if(var[5][i]>MaxX2) MaxX2 = var[5][i];
	if(Nvar==7) {
	  if(var[6][i]<MinX3) MinX3 = var[6][i];
	  if(var[6][i]>MaxX3) MaxX3 = var[6][i];
	}

	hScanX1->Fill(var[4][i],TMath::Abs(var[3][i]));
	hScanX2->Fill(var[5][i],TMath::Abs(var[3][i]));
	if(Nvar==7)
	  hScanX3->Fill(var[6][i],TMath::Abs(var[3][i]));

	hScanP2->Fill(var[1][i],TMath::Abs(var[3][i]));
	hScanP3->Fill(var[2][i],TMath::Abs(var[3][i]));
	
      }
      
      p1Min = MinP1 - rfactor*(MaxP1-MinP1);
      p1Max = MaxP1 + rfactor*(MaxP1-MinP1);
      p2Min = MinP2 - rfactor*(MaxP2-MinP2);
      p2Max = MaxP2 + rfactor*(MaxP2-MinP2);
      p3Min = MinP3 - rfactor*(MaxP3-MinP3);
      p3Max = MaxP3 + rfactor*(MaxP3-MinP3);

      x1Min = MinX1 - rfactor*(MaxX1-MinX1);
      x1Max = MaxX1 + rfactor*(MaxX1-MinX1);
      x2Min = MinX2 - rfactor*(MaxX2-MinX2);
      x2Max = MaxX2 + rfactor*(MaxX2-MinX2);

      
      if(Nvar==7) {
	x3Min = MinX3 - rfactor*(MaxX3-MinX3);
	x3Max = MaxX3 + rfactor*(MaxX3-MinX3);
      }

      // Set limits using Scan histograms
      Double_t peakFactor = 0.2;
      Double_t rfactor2 = 2;
      
      Double_t x1min = -999;
      Double_t x1max = -999;
      FindLimits(hScanX1,x1min,x1max,peakFactor);
      x1BinMin = x1min;
      x1BinMax = x1max;

      peakFactor = 0.05;
      FindLimits(hScanX1,x1min,x1max,peakFactor);
      x1Min = x1min - rfactor*(x1max-x1min);
      x1Max = x1max + rfactor*(x1max-x1min);
      
      Double_t x2min = -999;
      Double_t x2max = -999;
      FindLimits(hScanX2,x2min,x2max,peakFactor);
      x2Min = x2min - rfactor2*(x2max-x2min);
      x2Max = x2max + rfactor2*(x2max-x2min);

      Double_t x3min = -999;
      Double_t x3max = -999;
      if(Nvar==7) {
	FindLimits(hScanX3,x3min,x3max,peakFactor);
	x3Min = x3min - rfactor2*(x3max-x3min);
	x3Max = x3max + rfactor2*(x3max-x3min);
      }
      
      Double_t p2min = -999;
      Double_t p2max = -999;
      FindLimits(hScanP2,p2min,p2max,peakFactor);
      p2Min = p2min - rfactor2*(p2max-p2min);
      p2Max = p2max + rfactor2*(p2max-p2min);

      Double_t p3min = -999;
      Double_t p3max = -999;
      FindLimits(hScanP3,p3min,p3max,peakFactor);
      p3Min = p3min - rfactor2*(p3max-p3min);
      p3Max = p3max + rfactor2*(p3max-p3min);
                  
    }

    // Prevent range limits to go out of the box
    if(x1Min<X1MIN) x1Min = X1MIN;
    if(x1Max>X1MAX) x1Max = X1MAX;
    if(x2Min<X2MIN) x2Min = X2MIN;
    if(x2Max>X2MAX) x2Max = X2MAX;
    if(pData->Is3D()) {
      if(x3Min<X3MIN) x3Min = X3MIN;
      if(x3Max>X3MAX) x3Max = X3MAX;
    }

    // Adjust the binning
    x1Min = floor((x1Min-X1MIN)/dx1) * dx1 + X1MIN;  
    x1Max = floor((x1Max-X1MIN)/dx1) * dx1 + X1MIN;  
    x1Nbin = ceil ((x1Max - x1Min)/(dx1));
    // p1Nbin = x1Nbin;
      
    x2Min = floor((x2Min-X2MIN)/dx2) * dx2 + X2MIN;  
    x2Max = floor((x2Max-X2MIN)/dx2) * dx2 + X2MIN;  
    x2Nbin = ceil ((x2Max - x2Min)/(dx2));
    // p2Nbin = x2Nbin;
      
    if(pData->Is3D()) {
      x3Min = floor((x3Min-X3MIN)/dx3) * dx3 + X3MIN;  
      x3Max = floor((x3Max-X3MIN)/dx3) * dx3 + X3MIN;  
      x3Nbin = ceil ((x3Max - x3Min)/(dx3));
    }
    // p3Nbin = x3Nbin;
    // SNbin = ceil ((x1BinMax - x1BinMin)/(2.0*dx1));

    
    x1Min -= shiftz;
    x1Max -= shiftz;
    x1BinMin -= shiftz;
    x1BinMax -= shiftz;

    if(p1Min < 0.0) p1Min = 0.01;
    if(x1Nbin < 100) x1Nbin = 100;
    if(x2Nbin < 100) x2Nbin = 100;
    if(x3Nbin < 100) x3Nbin = 100;

    cout << Form(" x1 range (N = %i):  x1Min = %f  x1Max = %f  dx1 = %f", x1Nbin, x1Min, x1Max, (x1Max - x1Min)/x1Nbin) << endl;
    cout << Form(" p1 range (N = %i):  p1Min = %f  p1Max = %f  dp1 = %f", p1Nbin, p1Min, p1Max, (p1Max - p1Min)/p1Nbin) << endl;
    cout << Form(" x2 range (N = %i):  x2Min = %f  x2Max = %f  dx2 = %f", x2Nbin, x2Min, x2Max, (x2Max - x2Min)/x2Nbin) << endl;
    cout << Form(" p2 range (N = %i):  p2Min = %f  p2Max = %f  dp2 = %f", p2Nbin, p2Min, p2Max, (p2Max - p2Min)/p2Nbin) << endl;
    cout << Form(" x3 range (N = %i):  x3Min = %f  x3Max = %f  dx3 = %f", x3Nbin, x3Min, x3Max, (x3Max - x3Min)/x3Nbin) << endl;
    cout << Form(" p3 range (N = %i):  p3Min = %f  p3Max = %f  dp3 = %f", p3Nbin, p3Min, p3Max, (p3Max - p3Min)/p3Nbin) << endl;

    cout << Form("\n x1-slices (N = %i):  x1Min = %f  x1Max = %f  Dx1 = %f", SNbin, x1BinMin, x1BinMax, (x1BinMax-x1BinMin)/SNbin) << endl;

    // Histograms
    char hName[16];

    sprintf(hName,"hX1");
    TH1F *hX1 = (TH1F*) gROOT->FindObject(hName);
    if(hX1) delete hX1;
    hX1 = new TH1F(hName,"",x1Nbin,x1Min,x1Max);
    hX1->GetYaxis()->SetTitle("#Lambda");
    if(opt.Contains("comov"))
      hX1->GetXaxis()->SetTitle("k_{p}#zeta");
    else
      hX1->GetXaxis()->SetTitle("k_{p}z");
    
    sprintf(hName,"hP1");
    TH1F *hP1 = (TH1F*) gROOT->FindObject(hName);
    if(hP1) delete hP1;
    hP1 = new TH1F(hName,"",p1Nbin,p1Min,p1Max);
    hP1->GetYaxis()->SetTitle("p_{z}/mc");
    if(opt.Contains("comov"))
      hP1->GetXaxis()->SetTitle("k_{p}#zeta");
    else
      hP1->GetXaxis()->SetTitle("k_{p}z");
    
    sprintf(hName,"hP1X1");
    TH2F *hP1X1 = (TH2F*) gROOT->FindObject(hName);
    if(hP1X1) delete hP1X1;
    hP1X1 = new TH2F(hName,"",x1Nbin,x1Min,x1Max,p1Nbin,p1Min,p1Max);
    if(opt.Contains("comov"))
      hP1X1->GetXaxis()->SetTitle("k_{p}#zeta");
    else
      hP1X1->GetXaxis()->SetTitle("k_{p}z");
    hP1X1->GetYaxis()->SetTitle("p_{z}/mc");
    hP1X1->GetZaxis()->SetTitle("Charge [n_{0}dV]");
    hP1X1->GetZaxis()->CenterTitle();
  
    sprintf(hName,"hP2X2");
    TH2F *hP2X2 =  (TH2F*) gROOT->FindObject(hName);
    if(hP2X2) delete hP2X2;
    hP2X2 = new TH2F(hName,"",x2Nbin,x2Min,x2Max,p2Nbin,p2Min,p2Max);
    hP2X2->GetXaxis()->SetTitle("k_{p}x");
    hP2X2->GetYaxis()->SetTitle("p_{x}/mc");
    hP2X2->GetZaxis()->SetTitle("Charge [n_{0}dV]");

    sprintf(hName,"hP3X3");
    TH2F *hP3X3 =  (TH2F*) gROOT->FindObject(hName);
    if(hP3X3) delete hP3X3;
    if(pData->Is3D()) {
      hP3X3 = new TH2F(hName,"",x3Nbin,x3Min,x3Max,p3Nbin,p3Min,p3Max);
      hP3X3->GetXaxis()->SetTitle("k_{p}y");
      hP3X3->GetYaxis()->SetTitle("p_{y}/mc");
      hP3X3->GetZaxis()->SetTitle("Charge [n_{0}dV]");
      hP3X3->GetZaxis()->CenterTitle();
    }
    
    sprintf(hName,"hX2X1");
    TH2F *hX2X1 = (TH2F*) gROOT->FindObject(hName);
    if(hX2X1) delete hX2X1;
    hX2X1 = new TH2F(hName,"",x1Nbin,x1Min,x1Max,x2Nbin,x2Min,x2Max);
    hX2X1->GetXaxis()->SetTitle("k_{p}#zeta");
    hX2X1->GetYaxis()->SetTitle("k_{p}x");
    hX2X1->GetZaxis()->SetTitle("Charge [n_{0}dV]");

    sprintf(hName,"hX3X1");
    TH2F *hX3X1 = (TH2F*) gROOT->FindObject(hName);
    if(hX3X1) delete hX3X1;
    if(pData->Is3D()) {
      hX3X1 = new TH2F(hName,"",x1Nbin,x1Min,x1Max,x3Nbin,x3Min,x3Max);
      hX3X1->GetXaxis()->SetTitle("k_{p}#zeta");
      hX3X1->GetYaxis()->SetTitle("k_{p}y");
      hX3X1->GetZaxis()->SetTitle("Charge [n_{0}dV]");
    }
    
    sprintf(hName,"hX3X2");
    TH2F *hX3X2 = (TH2F*) gROOT->FindObject(hName);
    if(hX3X2) delete hX3X2;
    if(pData->Is3D()) {
      hX3X2 = new TH2F(hName,"",x2Nbin,x2Min,x2Max,x3Nbin,x3Min,x3Max);
      hX3X2->GetXaxis()->SetTitle("k_{p}x");
      hX3X2->GetYaxis()->SetTitle("k_{p}y");
      hX3X2->GetZaxis()->SetTitle("Charge [n_{0}dV]");
    }
    
    // Sliced quantities:
    // --------------------------------------------------------------------------
 
    // Set the binning
    Double_t *sBinLim = new Double_t[SNbin+1];
    sBinLim[0] = x1BinMin;
    sBinLim[SNbin] = x1BinMax;
    Double_t slbinSize = (sBinLim[SNbin] - sBinLim[0])/SNbin;
    for(Int_t i=1;i<SNbin;i++) {
      sBinLim[i] = sBinLim[i-1] + slbinSize;
    }

    TH1F **hP1sl = new TH1F*[SNbin];
    TH2F **hP2X2sl = new TH2F*[SNbin];
    TH2F **hP3X3sl = new TH2F*[SNbin];
    for(Int_t k=0;k<SNbin;k++) {
      sprintf(hName,"hP2X2sl_%2i",k);
      hP2X2sl[k] = (TH2F*) gROOT->FindObject(hName);
      if(hP2X2sl[k]) delete hP2X2sl[k];
      hP2X2sl[k] = new TH2F(hName,"",x2Nbin,x2Min,x2Max,p2Nbin,p2Min,p2Max);

      hP2X2sl[k]->GetXaxis()->SetTitle("k_{p}x");
      hP2X2sl[k]->GetYaxis()->SetTitle("p_{x}/mc");
      hP2X2sl[k]->GetZaxis()->SetTitle("Charge [n_{0}dV]");

      sprintf(hName,"hP3X3sl_%i",k);
      hP3X3sl[k] = (TH2F*) gROOT->FindObject(hName);
      if(hP3X3sl[k]) delete hP3X3sl[k];
      if(pData->Is3D()) {
	hP3X3sl[k] = new TH2F(hName,"",x3Nbin,x3Min,x3Max,p3Nbin,p3Min,p3Max);	
	hP3X3sl[k]->GetXaxis()->SetTitle("k_{p}y");
	hP3X3sl[k]->GetYaxis()->SetTitle("p_{y}/mc");
	hP3X3sl[k]->GetZaxis()->SetTitle("Charge [n_{0}dV]");      
      }
      
      sprintf(hName,"hP1sl_%2i",k);
      hP1sl[k] = (TH1F*) gROOT->FindObject(hName);
      if(hP1sl[k]) delete hP1sl[k];
      hP1sl[k] = new TH1F(hName,"",p1Nbin,p1Min,p1Max);

    }

    // Filling histos
    cout << Form("\n 2. Filling histograms from file ... ") ;
    
    for(UInt_t i=0;i<Np;i++) {

      var[4][i] -= shiftz;

      if(var[4][i]<x1Min || var[4][i]>x1Max ) continue; 
      if(var[5][i]<x2Min || var[5][i]>x2Max ) continue; 
      if(Nvar==7) {
	if(var[6][i]<x3Min || var[6][i]>x3Max ) continue; 
	if(var[2][i]<p3Min || var[2][i]>p3Max ) continue; 
      }
      if(var[0][i]<p1Min || var[0][i]>p1Max ) continue; 
      if(var[1][i]<p2Min || var[1][i]>p2Max ) continue; 
      
      //cout << " filling " << endl;

      hX1->Fill(var[4][i],TMath::Abs(var[3][i]));
      hP1->Fill(var[0][i],TMath::Abs(var[3][i]));

      hP1X1->Fill(var[4][i],var[0][i],TMath::Abs(var[3][i]));

      // Slices    
      if(var[4][i]<sBinLim[0] || var[4][i]>sBinLim[SNbin]) continue;
      Int_t iBin = -1;
      for(Int_t j=0; j<SNbin; j++) {

	if(var[4][i]<sBinLim[j+1]) {
	  iBin = j;
	  break;
	}
      }
      if(iBin<0) continue;

      // Projected emittance in the bunch range. (skip the tails to "match" the sliced ones)
      hP2X2->Fill(var[5][i],var[1][i],TMath::Abs(var[3][i]));
      hP2X2sl[iBin]->Fill(var[5][i],var[1][i],TMath::Abs(var[3][i]));

      if(Nvar==7){
	hP3X3->Fill(var[6][i],var[2][i],TMath::Abs(var[3][i]));
	hP3X3sl[iBin]->Fill(var[6][i],var[2][i],TMath::Abs(var[3][i]));

	hX3X1->Fill(var[4][i],var[6][i],TMath::Abs(var[3][i]));
	hX3X2->Fill(var[5][i],var[6][i],TMath::Abs(var[3][i]));      
      }

      hX2X1->Fill(var[4][i],var[5][i],TMath::Abs(var[3][i]));
      
      hP1sl[iBin]->Fill(var[0][i],TMath::Abs(var[3][i]));
      //      cout << Form("  iBin = %i   p1 = %7.3f",iBin,var[0][i]) << endl;
    
      
    }
    // cout << " done! " << endl;

    // Integrated long. emittance:
    cout << Form("\n 3. Calculating integrated quantities: ") << endl ;

    // Longitudinal phasespace (Total):
    // --

    Double_t stats[7];  // { sumw, sumw2, sumwx, sumwx2, sumwy, sumwy2, sumwxy }
    Double_t xmean,xrms2,xrms,ymean,yrms2,yrms,xyrms2,xyrms,emit2,emit;

    // P1X1 Phasespace
    hP1X1->GetStats(stats);
    
    xmean  = stats[2]/stats[0];
    xrms2  = stats[3]/stats[0] - xmean * xmean;
    xrms   = (xrms2>0.0) ? TMath::Sqrt(xrms2) : 0.0 ;
    ymean  = stats[4]/stats[0];
    yrms2  = stats[5]/stats[0] - ymean * ymean;
    yrms   = (yrms2>0.0) ? TMath::Sqrt(yrms2) : 0.0 ;
    xyrms2 = stats[6]/stats[0] - stats[2]*stats[4]/(stats[0]*stats[0]);
    xyrms  = (xyrms2>0.0) ? TMath::Sqrt(xyrms2) : 0.0 ;
    emit2 = xrms2*yrms2 - xyrms2*xyrms2;
    emit = (emit2>0.0) ? TMath::Sqrt(emit2) : 0.0 ;
    
    Double_t emitz = emit;
    Double_t zmean = xmean;
    Double_t zrms = xrms;
    Double_t pzmean = ymean;
    Double_t pzrms = yrms;
    // Double_t zpzcorr = xyrms;
    // Double_t zpzcorrrel = xyrms/pzmean;
    
    // Total relative energy spread within FWHM:
    sprintf(hName,"hP1cut");
    TH1F *hP1cut = (TH1F*) hP1->Clone(hName);
    hP1cut->Reset();

    Double_t maxValue = hP1->GetBinContent(hP1->GetMaximumBin());
    Int_t   lBin = -1;
    Double_t pzmin = -999;
    for(Int_t i=1;i<=hP1->GetNbinsX();i++) {
      Double_t binValue = hP1->GetBinContent(i);
      if(binValue>maxValue/2) {
	lBin = i;
	pzmin = hP1->GetBinCenter(i);
	break;
      }
    }

    Int_t rBin = -1;
    Double_t pzmax = -999;
    for(Int_t i=hP1->GetNbinsX();i>0;i--) {
      Double_t binValue = hP1->GetBinContent(i);
      if(binValue>maxValue/2) {
	rBin = i;
	pzmax = hP1->GetBinCenter(i);
	break;
      }
    }
    
    for(Int_t i=lBin;i<=rBin;i++) {
      Double_t binValue = hP1->GetBinContent(i);
      hP1cut->SetBinContent(i,binValue);
    }
    //  hP1cut->ResetStats();
    Double_t pzmeanFWHM = hP1cut->GetMean();
    Double_t pzrmsFWHM = hP1cut->GetRMS();
        
    // cout << Form("  zMean = %7.3f   pzMean = %7.3f",zmean,pzmean) << endl;
    // cout << Form("  zRms  = %7.3f   pzRms  = %7.3f",zrms,pzrms) << endl;
    // cout << Form("  zEmittance = %7.3f",emitz) << endl;
    
    // cout <<  Form("  Normal  : pzMean  = %7.3f   pzRms  = %7.3f",pzmean,pzrms) << endl;
    // cout <<  Form("  In FWHM : pzMean  = %7.3f   pzRms  = %7.3f -> Interval: pmin = %7.3f pmax = %7.3f  ",pzmeanFWHM,pzrmsFWHM,pzmin,pzmax) << endl;

    
    // ----------------------------------------


    // Transverse phasespace (Total)

    // P2X2 Phasespace
    hP2X2->GetStats(stats);

    xrms2  = stats[3]/stats[0] - stats[2]*stats[2]/(stats[0]*stats[0]);
    yrms2  = stats[5]/stats[0] - stats[4]*stats[4]/(stats[0]*stats[0]);
    xyrms2 = stats[6]/stats[0] - stats[2]*stats[4]/(stats[0]*stats[0]);
    emit2 = xrms2*yrms2 - xyrms2*xyrms2;
    emit = (emit2>0.0) ? TMath::Sqrt(emit2) : 0.0 ; 

    Double_t emitx = emit;
    Double_t x_mean = stats[2]/stats[0];
    Double_t x_rms = (xrms2>0.0) ? TMath::Sqrt(xrms2) : 0.0 ;
    Double_t px_mean = stats[4]/stats[0];
    Double_t px_rms = (yrms2>0.0) ? TMath::Sqrt(yrms2) : 0.0 ;

    Double_t beta = xrms2 / emit;
    Double_t gamma = yrms2 / emit;
    Double_t alpha = -xyrms2 / emit;
  
    Double_t grel = hP1->GetMean() * eneUnit/PConst::ElectronMassE;
    Double_t betax = beta * grel;

    Double_t factor =  beta*beta + 2 * beta * gamma + gamma*gamma - 4 * emit;
    if(factor<0.0) factor *= -1;
    factor = TMath::Sqrt(factor);
  
    Double_t a2 = 0.5 * (beta + gamma - factor);
    Double_t b2 = 0.5 * (beta + gamma + factor);
  
    Double_t p  = alpha / (a2-b2);
    Double_t pf = TMath::Sqrt( 1 - 4*p*p ); 
    Double_t angle = - TMath::ATan( (1+pf) / (2*p) );
    //  if (angle <0.0) angle += 2*PConst::pi;
  
    TEllipse *ellipP2X2 = new TEllipse(x_mean,px_mean,TMath::Sqrt(a2),TMath::Sqrt(b2),0.,360.,angle * 180. / PConst::pi );
    ellipP2X2->SetFillStyle(0);
    ellipP2X2->SetLineStyle(2);
    ellipP2X2->SetLineColor(2);
    ellipP2X2->SetLineWidth(1);


    // P3X3 SPACE -----------------------
    Double_t emity = 0, y_mean = 0, y_rms = 0, py_mean = 0, py_rms = 0, betay = 0;
    TEllipse *ellipP3X3 = NULL;
    if(pData->Is3D()) {
      hP3X3->GetStats(stats);

      xrms2  = stats[3]/stats[0] - stats[2]*stats[2]/(stats[0]*stats[0]);
      yrms2  = stats[5]/stats[0] - stats[4]*stats[4]/(stats[0]*stats[0]);
      xyrms2 = stats[6]/stats[0] - stats[2]*stats[4]/(stats[0]*stats[0]);
      emit2 = xrms2*yrms2 - xyrms2*xyrms2;
      emit = (emit2>0.0) ? TMath::Sqrt(emit2) : 0.0 ; 
      
      emity = emit;
      y_mean = stats[2]/stats[0];
      y_rms = (xrms2>0.0) ? TMath::Sqrt(xrms2) : 0.0 ;
      py_mean = stats[4]/stats[0];
      py_rms = (yrms2>0.0) ? TMath::Sqrt(yrms2) : 0.0 ;
      
      beta = xrms2 / emit;
      gamma = yrms2 / emit;
      alpha = -xyrms2 / emit;
      
      betay = beta * grel;

      factor =  beta*beta + 2 * beta * gamma + gamma*gamma - 4 * emit;
      if(factor<0.0) factor *= -1;
      factor = TMath::Sqrt(factor);
      
      a2 = 0.5 * (beta + gamma - factor);
      b2 = 0.5 * (beta + gamma + factor);
      
      p  = alpha / (a2-b2);
      pf = TMath::Sqrt( 1 - 4*p*p ); 
      angle = - TMath::ATan( (1+pf) / (2*p) );
      //  if (angle <0.0) angle += 2*PConst::pi;

      ellipP3X3 = new TEllipse(x_mean,px_mean,TMath::Sqrt(a2),TMath::Sqrt(b2),0.,360.,angle * 180. / PConst::pi );
      ellipP3X3->SetFillStyle(0);
      ellipP3X3->SetLineStyle(2);
      ellipP3X3->SetLineColor(2);
      ellipP3X3->SetLineWidth(1);
    }

    cout << Form("\n 4. Calculating sliced quantities.. ") << endl ;

    TGraph *gEmitx = NULL;
    TGraph *gEmity = NULL;
    TGraph *gXrms = NULL;
    TGraph *gYrms = NULL;
    TGraph *gErms = NULL;

    Double_t * zbin = new Double_t[SNbin];
    Double_t * sEmean = new Double_t[SNbin];
    Double_t * sErms = new Double_t[SNbin];

    Double_t *sx_mean = new Double_t[SNbin];
    Double_t *sx_rms = new Double_t[SNbin];
    Double_t *spx_mean = new Double_t[SNbin];
    Double_t *spx_rms = new Double_t[SNbin];
    Double_t *semitx = new Double_t[SNbin];
    Double_t *sbetax = new Double_t[SNbin];

    TEllipse **sellipP2X2 = new TEllipse*[SNbin];

    Double_t *sy_mean = NULL, *sy_rms = NULL, *spy_mean = NULL, *spy_rms = NULL, *semity = NULL, *sbetay = NULL;
    TEllipse **sellipP3X3 = NULL;

    if(pData->Is3D()) {
      sy_mean = new Double_t[SNbin];
      sy_rms = new Double_t[SNbin];
      spy_mean = new Double_t[SNbin];
      spy_rms = new Double_t[SNbin];
      semity = new Double_t[SNbin];
      sbetay = new Double_t[SNbin];
      
      sellipP3X3 = new TEllipse*[SNbin];
    }
    
    for(Int_t k=0;k<SNbin;k++) {

      Double_t sxmean = 0, sxrms2 = 0, sxrms = 0, symean = 0, syrms2 = 0, syrms = 0, sxyrms2 = 0, sxyrms = 0, semit = 0, semit2 = 0;
      
      // P2X2 slices
      sellipP2X2[k] = NULL;            
      hP2X2sl[k]->GetStats(stats);
      
      sxmean  = stats[2]/stats[0];
      sxrms2  = stats[3]/stats[0] - sxmean * sxmean;
      sxrms   = (sxrms2>0.0) ? TMath::Sqrt(sxrms2) : 0.0 ;
      symean  = stats[4]/stats[0];
      syrms2  = stats[5]/stats[0] - symean * symean;
      syrms   = (syrms2>0.0) ? TMath::Sqrt(syrms2) : 0.0 ;
      sxyrms2 = stats[6]/stats[0] - stats[2]*stats[4]/(stats[0]*stats[0]);
      sxyrms  = (sxyrms2>0.0) ? TMath::Sqrt(sxyrms2) : 0.0 ;
      semit2 = sxrms2*syrms2 - sxyrms2*sxyrms2;
      semit = (semit2>0.0) ? TMath::Sqrt(semit2) : 0.0 ;

      zbin[k] = (sBinLim[k] + sBinLim[k+1])/2.;

      sEmean[k] = hP1sl[k]->GetMean();
      sErms[k]  = hP1sl[k]->GetRMS();

      sx_mean[k] = sxmean;
      spx_mean[k] = symean;
      sx_rms[k] = sxrms;
      spx_rms[k] = syrms;
      semitx[k] = semit;
      
      // Ellipse calculation
      Double_t beta, gamma, alpha, factor, a2, b2, p, pf, angle;
      Double_t zoomEllip = 3.0;
      
      beta = sxrms2 / semit;
      gamma = syrms2 / semit ;
      alpha = - sxyrms2 / semit;

      sbetax[k] = beta * sEmean[k];
      
      factor =  beta*beta + 2 * beta * gamma + gamma*gamma - 4 * semit;
      if(factor<0.0) factor *= -1;
      factor = TMath::Sqrt(factor);

      a2 = 0.5 * (beta + gamma - factor);
      b2 = 0.5 * (beta + gamma + factor);
 
      p  = alpha / (a2-b2);
      pf = TMath::Sqrt( 1 - 4*p*p ); 
      angle = - TMath::ATan( (1+pf) / (2*p) );

      //  if (angle <0.0) angle += 2*PConst::pi;    
      sellipP2X2[k] = new TEllipse(sx_mean[k],spx_mean[k],TMath::Sqrt(zoomEllip*a2),TMath::Sqrt(zoomEllip*b2),0.,360.,angle * 180. / PConst::pi );
      

      // P3X3 slices
      if(pData->Is3D()) {
	sellipP3X3[k] = NULL;            
	hP3X3sl[k]->GetStats(stats);
	
	sxmean  = stats[2]/stats[0];
	sxrms2  = stats[3]/stats[0] - sxmean * sxmean;
	sxrms   = (sxrms2>0.0) ? TMath::Sqrt(sxrms2) : 0.0 ;
	symean  = stats[4]/stats[0];
	syrms2  = stats[5]/stats[0] - symean * symean;
	syrms   = (syrms2>0.0) ? TMath::Sqrt(syrms2) : 0.0 ;
	sxyrms2 = stats[6]/stats[0] - stats[2]*stats[4]/(stats[0]*stats[0]);
	sxyrms  = (sxyrms2>0.0) ? TMath::Sqrt(sxyrms2) : 0.0 ;
	semit2 = sxrms2*syrms2 - sxyrms2*sxyrms2;
	semit = (semit2>0.0) ? TMath::Sqrt(semit2) : 0.0 ;
	
	sy_mean[k] = sxmean;
	spy_mean[k] = symean;
	sy_rms[k] = sxrms;
	spy_rms[k] = syrms;
	semity[k] = semit;
	
	// Ellipse calculation
	beta = sxrms2 / semit;
	gamma = syrms2 / semit ;
	alpha = - sxyrms2 / semit;
	
	sbetay[k] = beta * sEmean[k];
      
	factor =  beta*beta + 2 * beta * gamma + gamma*gamma - 4 * semit;
	if(factor<0.0) factor *= -1;
	factor = TMath::Sqrt(factor);
	
	a2 = 0.5 * (beta + gamma - factor);
	b2 = 0.5 * (beta + gamma + factor);
	
	p  = alpha / (a2-b2);
	pf = TMath::Sqrt( 1 - 4*p*p ); 
	angle = - TMath::ATan( (1+pf) / (2*p) );
	
	//  if (angle <0.0) angle += 2*PConst::pi;    
	sellipP3X3[k] = new TEllipse(sy_mean[k],spy_mean[k],TMath::Sqrt(zoomEllip*a2),TMath::Sqrt(zoomEllip*b2),0.,360.,angle * 180. / PConst::pi );
      }
    }
    //    cout << " done! " << endl;
    

    // Changing to user units: 
    // --------------------------

    // Charge
    Double_t dV = dx1*dx2*dx3 * skindepth * skindepth * skindepth;
    Double_t Q0 = rawf * fabs(n0 * dV * PConst::ElectronCharge);
    Double_t Charge = hX1->Integral() * Q0;
    if(opt.Contains("best")) {
      PUnits::BestUnit bchargeSUnit(Charge,"Charge");
      bchargeSUnit.GetBestUnits(chargeUnit,chargeSUnit);
    }
    Charge /= chargeUnit;
    
    if(opt.Contains("units") && n0) {

      Time *= skindepth / spaUnit;

      x1Min *= skindepth / spaUnit;
      x1Max *= skindepth / spaUnit;
      p1Min *= PConst::ElectronMassE / eneUnit;
      p1Max *= PConst::ElectronMassE / eneUnit;
      
      x1BinMin *= skindepth / spaUnit;
      x1BinMax *= skindepth / spaUnit;
      for(Int_t k=0;k<=SNbin;k++) 
	sBinLim[k] *= skindepth / spaUnit;
      
      hP1X1->SetBins(x1Nbin,x1Min,x1Max,p1Nbin,p1Min,p1Max);
      hP1X1->Scale(Q0 / chargeUnit);
      
      // Converting electron density
      hP1->Scale(Q0 / chargeUnit);
      hP1->SetBins(p1Nbin,p1Min,p1Max);
      hP1->GetYaxis()->SetTitle(Form("p_{z} [%s/c]",eneSUnit.c_str()));

      hX1->SetBins(x1Nbin,x1Min,x1Max);
      Double_t dt = ((x1Max - x1Min)/x1Nbin) * spaUnit / PConst::c_light;
      hX1->Scale(Q0 / dt);

      if(opt.Contains("best")) {
	PUnits::BestUnit bcurSUnit(hX1->GetMaximum(),"Current");
	bcurSUnit.GetBestUnits(curUnit,curSUnit);
      }
      hX1->Scale(1/curUnit);
      
      hX1->GetYaxis()->SetTitle(Form("I[%s]",curSUnit.c_str()));
      hX1->GetYaxis()->SetTitle("");
      if(opt.Contains("comov"))
	hX1->GetXaxis()->SetTitle(Form("#zeta [%s]",spaSUnit.c_str()));
      else
	hX1->GetXaxis()->SetTitle(Form("z [%s]",spaSUnit.c_str()));
      
      
      zmean *= skindepth / spaUnit;
      zrms  *= skindepth / spaUnit;

      pzmean *= PConst::ElectronMassE / eneUnit;
      pzrms  *= PConst::ElectronMassE / eneUnit;
      pzmeanFWHM *= PConst::ElectronMassE / eneUnit;
      pzrmsFWHM *= PConst::ElectronMassE / eneUnit;
      pzmin  *= PConst::ElectronMassE / eneUnit;
      pzmax  *= PConst::ElectronMassE / eneUnit;
      emitz *= (skindepth / spaUnit);

      // Transverse phase-space
      x_mean *= skindepth / tspaUnit;
      x_rms  *= skindepth / tspaUnit;
      px_mean *= PConst::ElectronMassE / teneUnit;
      px_rms  *= PConst::ElectronMassE / teneUnit;

      x2Min *= skindepth/tspaUnit;
      x2Max *= skindepth/tspaUnit;
      p2Min *= PConst::ElectronMassE / teneUnit;
      p2Max *= PConst::ElectronMassE / teneUnit;
      
      hP2X2->SetBins(x2Nbin,x2Min,x2Max,p2Nbin,p2Min,p2Max);
      hP2X2->Scale(Q0 / chargeUnit);
      
      hP2X2->GetXaxis()->SetTitle(Form("x [%s]",spaSUnit.c_str()));
      hP2X2->GetYaxis()->SetTitle(Form("p_{x} [%s/c]",teneSUnit.c_str()));
      hP2X2->GetZaxis()->SetTitle(Form("Charge [%s]",chargeSUnit.c_str()));
      hP2X2->GetZaxis()->CenterTitle();

      y_mean *= skindepth / tspaUnit;
      y_rms  *= skindepth / tspaUnit;
      py_mean *= PConst::ElectronMassE / teneUnit;
      py_rms  *= PConst::ElectronMassE / teneUnit;      
      
      x3Min *= skindepth/tspaUnit;
      x3Max *= skindepth/tspaUnit;
      p3Min *= PConst::ElectronMassE / teneUnit;
      p3Max *= PConst::ElectronMassE / teneUnit;
      
      hX2X1->SetBins(x1Nbin,x1Min,x1Max,x2Nbin,x2Min,x2Max);
      hX2X1->Scale(Q0 / chargeUnit);
      hX2X1->GetXaxis()->SetTitle(Form("#zeta [%s]",spaSUnit.c_str()));
      hX2X1->GetYaxis()->SetTitle(Form("x [%s]",tspaSUnit.c_str()));
      hX2X1->GetZaxis()->SetTitle(Form("Charge [%s]",chargeSUnit.c_str()));
      hX2X1->GetZaxis()->CenterTitle();

      if(pData->Is3D()) {
	hP3X3->SetBins(x3Nbin,x3Min,x3Max,p3Nbin,p3Min,p3Max);
	hP3X3->Scale(Q0 / chargeUnit);
	
	hP3X3->GetXaxis()->SetTitle(Form("y [%s]",spaSUnit.c_str()));
	hP3X3->GetYaxis()->SetTitle(Form("p_{y} [%s/c]",teneSUnit.c_str()));
	hP3X3->GetZaxis()->SetTitle(Form("Charge [%s]",chargeSUnit.c_str()));
	hP3X3->GetZaxis()->CenterTitle();
	
	hX3X1->SetBins(x1Nbin,x1Min,x1Max,x3Nbin,x3Min,x3Max);
	hX3X1->Scale(Q0 / chargeUnit);
	hX3X1->GetXaxis()->SetTitle(Form("#zeta [%s]",spaSUnit.c_str()));
	hX3X1->GetYaxis()->SetTitle(Form("y [%s]",tspaSUnit.c_str()));
	hX3X1->GetZaxis()->SetTitle(Form("Charge [%s]",chargeSUnit.c_str()));
	hX3X1->GetZaxis()->CenterTitle();
	
	hX3X2->SetBins(x2Nbin,x2Min,x2Max,x3Nbin,x3Min,x3Max);
	hX3X2->Scale(Q0 / chargeUnit);
	hX3X2->GetXaxis()->SetTitle(Form("x [%s]",tspaSUnit.c_str()));
	hX3X2->GetYaxis()->SetTitle(Form("y [%s]",tspaSUnit.c_str()));
	hX3X2->GetZaxis()->SetTitle(Form("Charge [%s]",chargeSUnit.c_str()));
	hX3X2->GetZaxis()->CenterTitle();
      }
      
      emitx *= skindepth;
      if(opt.Contains("best")) {
	PUnits::BestUnit bemitSUnit(emitx,"Emittance");
	bemitSUnit.GetBestUnits(emitUnit,emitSUnit);
	// cout << bemitSUnit << endl;
      }
      emitx /= emitUnit;

      betax *= skindepth;
      if(opt.Contains("best")) {
	PUnits::BestUnit bbetaSUnit(betax,"Length");
	bbetaSUnit.GetBestUnits(betaUnit,betaSUnit);
	// cout << bbetaSUnit << endl;
      }
      betax /= betaUnit;

      if(pData->Is3D()) {
	emity *= skindepth / emitUnit;
	betay *= skindepth / betaUnit;
      }
      
      // cout << Form(" %.6f   %s",emitUnit,emitSUnit.c_str())  << endl;
      
      Double_t erelMax = -999.;
      for(Int_t k=0;k<SNbin;k++) {

	zbin[k] *= skindepth / spaUnit;
	sErms[k] = (sErms[k]/sEmean[k]) / PUnits::perCent;
	sEmean[k] *= PConst::ElectronMassE / eneUnit;

	hP2X2sl[k]->SetBins(x2Nbin,x2Min,x2Max,p2Nbin,p2Min,p2Max);
	hP2X2sl[k]->Scale(Q0 / chargeUnit);
   
	hP2X2sl[k]->GetXaxis()->SetTitle(Form("x [%s]",spaSUnit.c_str()));
	hP2X2sl[k]->GetYaxis()->SetTitle(Form("p_{x} [%s/c]",teneSUnit.c_str()));
	hP2X2sl[k]->GetZaxis()->SetTitle(Form("Charge [%s]",chargeSUnit.c_str()));

	hP2X2sl[k]->GetZaxis()->CenterTitle();

	sx_mean[k] *= skindepth / tspaUnit;
	sx_rms[k]  *= skindepth / tspaUnit;
	spx_mean[k] *= PConst::ElectronMassE / teneUnit;
	spx_rms[k] *= PConst::ElectronMassE / teneUnit;
	semitx[k] *= skindepth / emitUnit;
	sbetax[k] *= skindepth / betaUnit;
	
	if(pData->Is3D()) {
	  hP3X3sl[k]->SetBins(x3Nbin,x3Min,x3Max,p3Nbin,p3Min,p3Max);
	  hP3X3sl[k]->Scale(Q0 / chargeUnit);
	  
	  hP3X3sl[k]->GetXaxis()->SetTitle(Form("x [%s]",spaSUnit.c_str()));
	  hP3X3sl[k]->GetYaxis()->SetTitle(Form("p_{x} [%s/c]",teneSUnit.c_str()));
	  hP3X3sl[k]->GetZaxis()->SetTitle(Form("Charge [%s]",chargeSUnit.c_str()));
	  
	  hP3X3sl[k]->GetZaxis()->CenterTitle();
	  
	  sy_mean[k] *= skindepth / tspaUnit;
	  sy_rms[k]  *= skindepth / tspaUnit;
	  spy_mean[k] *= PConst::ElectronMassE / teneUnit;
	  spy_rms[k] *= PConst::ElectronMassE / teneUnit;
	  semity[k] *= skindepth / emitUnit;
	  sbetay[k] *= skindepth / betaUnit;
	}
	
	if(sErms[k]>erelMax) erelMax = sErms[k];
	
      }

      if(opt.Contains("best")) {
	PUnits::BestUnit berelSUnit(erelMax * PUnits::perCent,"Percentage");
	berelSUnit.GetBestUnits(ermsUnit,ermsSUnit);
      }
      erelMax = PUnits::perCent/ermsUnit;
      
      for(Int_t k=0;k<SNbin;k++)
       	sErms[k]  *= PUnits::perCent/ermsUnit;
      
    }
    // End of the users units module    
    
    cout << "\n  Summary _______________________________________________________ " << endl;
    if(opt.Contains("units")) {
      cout << Form("  Integrated charge (RAW) of specie %3i = %8f %s",index,Charge,chargeSUnit.c_str()) << endl;
      cout << Form("  Peak current = %6.2f %s",hX1->GetMaximum(),curSUnit.c_str()) << endl;
      cout << Form("  Total energy = %6.2f %s, rms = %3.1f %s",pzmean,eneSUnit.c_str(),(pzrms/pzmean)/ermsUnit,ermsSUnit.c_str()) << endl;
      cout << Form("  Length = %6.2f %s (rms)",zrms,spaSUnit.c_str()) << endl;
      cout << Form("  Width x = %6.2f %s (rms)",x_rms,tspaSUnit.c_str()) << endl;
      cout << Form("  Trans. emit. x = %6.2f %s",emitx,emitSUnit.c_str()) << endl;
      if(pData->Is3D()) {
	cout << Form("  Width y = %6.2f %s (rms)",y_rms,tspaSUnit.c_str()) << endl;
	cout << Form("  Trans. emit. y = %6.2f %s",emity,emitSUnit.c_str()) << endl;
      }
    }
    
    if(opt.Contains("loop")) {
      cout << Form("\n 5. Saving results to file .. ") << endl;
  
      // OUTPUT ROOT FILE WITH THE PLOTS:
      TString filename = Form("./%s/Plots/Bunch/%s/Bunch-Evolution-%s.root",sim.Data(),pData->GetSpeciesName(index).c_str(),sim.Data());
      TFile * ifile = (TFile*) gROOT->GetListOfFiles()->FindObject(filename);
      // if doesn't exist the directory should be created
      if (!ifile) {
	TString f = filename;
	TString dir3 = f.Remove( f.Last( '/' ), f.Length() - f.Last( '/' ) );
	TString dir2 = f.Remove( f.Last( '/' ), f.Length() - f.Last( '/' ) );
	TString dir1 = f.Remove( f.Last( '/' ), f.Length() - f.Last( '/' ) );
	
	gSystem->mkdir( dir1 );
	gSystem->mkdir( dir2 );
	gSystem->mkdir( dir3 );
	ifile = new TFile(filename,"UPDATE");
      }  

      Int_t nPoints = 0;
      char gName[32];
      sprintf(gName,"gEmitxvsTime");     
      TGraph *gEmitxvsTime = NULL;
      gEmitxvsTime = (TGraph*) ifile->Get(gName);
      if(gEmitxvsTime==NULL) {
	gEmitxvsTime = new TGraph();
	gEmitxvsTime->SetName(gName);
	nPoints = 0;
	// Some cosmetics at creation time:
	gEmitxvsTime->SetLineWidth(3);
	gEmitxvsTime->SetLineColor(PGlobals::fieldLine);
	gEmitxvsTime->SetMarkerStyle(20);
	gEmitxvsTime->SetMarkerSize(0.4);
	gEmitxvsTime->SetMarkerColor(PGlobals::fieldLine);	
      } else {
	nPoints = gEmitxvsTime->GetN(); 
      }  

      gEmitxvsTime->Set(nPoints+1);
      gEmitxvsTime->SetPoint(nPoints,Time,emitx);
      gEmitxvsTime->Write(gName,TObject::kOverwrite);


      sprintf(gName,"gPzmeanvsTime");     
      TGraph *gPzmeanvsTime = NULL;
      gPzmeanvsTime = (TGraph*) ifile->Get(gName);
      if(gPzmeanvsTime==NULL) {
	gPzmeanvsTime = new TGraph();
	gPzmeanvsTime->SetName(gName);
	nPoints = 0;
	// Some cosmetics at creation time:
	gPzmeanvsTime->SetLineWidth(3);
	gPzmeanvsTime->SetLineColor(PGlobals::fieldLine);
	gPzmeanvsTime->SetMarkerStyle(20);
	gPzmeanvsTime->SetMarkerSize(0.4);
	gPzmeanvsTime->SetMarkerColor(PGlobals::fieldLine);	
      } else {
	nPoints = gPzmeanvsTime->GetN(); 
      }  

      gPzmeanvsTime->Set(nPoints+1);
      gPzmeanvsTime->SetPoint(nPoints,Time,pzmean);
      gPzmeanvsTime->Write(gName,TObject::kOverwrite);


      sprintf(gName,"gPzrmsvsTime");     
      TGraph *gPzrmsvsTime = NULL;
      gPzrmsvsTime = (TGraph*) ifile->Get(gName);
      if(gPzrmsvsTime==NULL) {
	gPzrmsvsTime = new TGraph();
	gPzrmsvsTime->SetName(gName);
	nPoints = 0;
	// Some cosmetics at creation time:
	gPzrmsvsTime->SetLineWidth(3);
	gPzrmsvsTime->SetLineColor(PGlobals::fieldLine);
	gPzrmsvsTime->SetMarkerStyle(20);
	gPzrmsvsTime->SetMarkerSize(0.4);
	gPzrmsvsTime->SetMarkerColor(PGlobals::fieldLine);	
      } else {
	nPoints = gPzrmsvsTime->GetN(); 
      }  

      gPzrmsvsTime->Set(nPoints+1);
      gPzrmsvsTime->SetPoint(nPoints,Time,pzrms);
      gPzrmsvsTime->Write(gName,TObject::kOverwrite);


      sprintf(gName,"gZmeanvsTime");     
      TGraph *gZmeanvsTime = NULL;
      gZmeanvsTime = (TGraph*) ifile->Get(gName);
      if(gZmeanvsTime==NULL) {
	gZmeanvsTime = new TGraph();
	gZmeanvsTime->SetName(gName);
	nPoints = 0;
	// Some cosmetics at creation time:
	gZmeanvsTime->SetLineWidth(3);
	gZmeanvsTime->SetLineColor(PGlobals::fieldLine);
	gZmeanvsTime->SetMarkerStyle(20);
	gZmeanvsTime->SetMarkerSize(0.4);
	gZmeanvsTime->SetMarkerColor(PGlobals::fieldLine);	
      } else {
	nPoints = gZmeanvsTime->GetN(); 
      }  

      gZmeanvsTime->Set(nPoints+1);
      gZmeanvsTime->SetPoint(nPoints,Time,zmean);
      gZmeanvsTime->Write(gName,TObject::kOverwrite);

      sprintf(gName,"gZrmsvsTime");     
      TGraph *gZrmsvsTime = NULL;
      gZrmsvsTime = (TGraph*) ifile->Get(gName);
      if(gZrmsvsTime == NULL) {
	gZrmsvsTime = new TGraph();
	gZrmsvsTime->SetName(gName);
	nPoints = 0;
	// Some cosmetics at creation time:
	gZrmsvsTime->SetLineWidth(3);
	gZrmsvsTime->SetLineColor(PGlobals::fieldLine);
	gZrmsvsTime->SetMarkerStyle(20);
	gZrmsvsTime->SetMarkerSize(0.4);
	gZrmsvsTime->SetMarkerColor(PGlobals::fieldLine);	
      } else {
	nPoints = gZrmsvsTime->GetN(); 
      }  

      gZrmsvsTime->Set(nPoints+1);
      gZrmsvsTime->SetPoint(nPoints,Time,zrms);
      gZrmsvsTime->Write(gName,TObject::kOverwrite);

      sprintf(gName,"gPxmeanvsTime");     
      TGraph *gPxmeanvsTime = NULL;
      gPxmeanvsTime = (TGraph*) ifile->Get(gName);
      if(gPxmeanvsTime==NULL) {
	gPxmeanvsTime = new TGraph();
	gPxmeanvsTime->SetName(gName);
	nPoints = 0;
	// Some cosmetics at creation time:
	gPxmeanvsTime->SetLineWidth(3);
	gPxmeanvsTime->SetLineColor(PGlobals::fieldLine);
	gPxmeanvsTime->SetMarkerStyle(20);
	gPxmeanvsTime->SetMarkerSize(0.4);
	gPxmeanvsTime->SetMarkerColor(PGlobals::fieldLine);	
      } else {
	nPoints = gPxmeanvsTime->GetN(); 
      }  

      gPxmeanvsTime->Set(nPoints+1);
      gPxmeanvsTime->SetPoint(nPoints,Time,px_mean);
      gPxmeanvsTime->Write(gName,TObject::kOverwrite);

      sprintf(gName,"gPxrmsvsTime");     
      TGraph *gPxrmsvsTime = NULL;
      gPxrmsvsTime = (TGraph*) ifile->Get(gName);
      if(gPxrmsvsTime==NULL) {
	gPxrmsvsTime = new TGraph();
	gPxrmsvsTime->SetName(gName);
	nPoints = 0;
	// Some cosmetics at creation time:
	gPxrmsvsTime->SetLineWidth(3);
	gPxrmsvsTime->SetLineColor(PGlobals::fieldLine);
	gPxrmsvsTime->SetMarkerStyle(20);
	gPxrmsvsTime->SetMarkerSize(0.4);
	gPxrmsvsTime->SetMarkerColor(PGlobals::fieldLine);	
      } else {
	nPoints = gPxrmsvsTime->GetN(); 
      }  

      gPxrmsvsTime->Set(nPoints+1);
      gPxrmsvsTime->SetPoint(nPoints,Time,px_rms);
      gPxrmsvsTime->Write(gName,TObject::kOverwrite);



      sprintf(gName,"gXmeanvsTime");     
      TGraph *gXmeanvsTime = NULL;
      gXmeanvsTime = (TGraph*) ifile->Get(gName);
      if(gXmeanvsTime==NULL) {
	gXmeanvsTime = new TGraph();
	gXmeanvsTime->SetName(gName);
	nPoints = 0;
	// Some cosmetics at creation time:
	gXmeanvsTime->SetLineWidth(3);
	gXmeanvsTime->SetLineColor(PGlobals::fieldLine);
	gXmeanvsTime->SetMarkerStyle(20);
	gXmeanvsTime->SetMarkerSize(0.4);
	gXmeanvsTime->SetMarkerColor(PGlobals::fieldLine);	
      } else {
	nPoints = gXmeanvsTime->GetN(); 
      }  

      gXmeanvsTime->Set(nPoints+1);
      gXmeanvsTime->SetPoint(nPoints,Time,x_mean);
      gXmeanvsTime->Write(gName,TObject::kOverwrite);

      sprintf(gName,"gXrmsvsTime");     
      TGraph *gXrmsvsTime = NULL;
      gXrmsvsTime = (TGraph*) ifile->Get(gName);
      if(gXrmsvsTime==NULL) {
	gXrmsvsTime = new TGraph();
	gXrmsvsTime->SetName(gName);
	nPoints = 0;
	// Some cosmetics at creation time:
	gXrmsvsTime->SetLineWidth(3);
	gXrmsvsTime->SetLineColor(PGlobals::fieldLine);
	gXrmsvsTime->SetMarkerStyle(20);
	gXrmsvsTime->SetMarkerSize(0.4);
	gXrmsvsTime->SetMarkerColor(PGlobals::fieldLine);	
      } else {
	nPoints = gXrmsvsTime->GetN(); 
      }  

      gXrmsvsTime->Set(nPoints+1);
      gXrmsvsTime->SetPoint(nPoints,Time,x_rms);
      gXrmsvsTime->Write(gName,TObject::kOverwrite);


      sprintf(gName,"gChargevsTime");     
      TGraph *gChargevsTime = NULL;
      gChargevsTime = (TGraph*) ifile->Get(gName);
      if(gChargevsTime==NULL) {
	gChargevsTime = new TGraph();
	gChargevsTime->SetName(gName);
	nPoints = 0;
	// Some cosmetics at creation time:
	gChargevsTime->SetLineWidth(3);
	gChargevsTime->SetLineColor(PGlobals::fieldLine);
	gChargevsTime->SetMarkerStyle(20);
	gChargevsTime->SetMarkerSize(0.4);
	gChargevsTime->SetMarkerColor(PGlobals::fieldLine);	
      } else {
	nPoints = gChargevsTime->GetN(); 
      }  

      gChargevsTime->Set(nPoints+1);
      gChargevsTime->SetPoint(nPoints,Time,fabs(Charge));
      gChargevsTime->Write(gName,TObject::kOverwrite);

      ifile->Close();
    }

    // ------------------------------------------------------------------------------------

    // Free memory from the dynamic array of variables:
    for(UInt_t i=0;i<Nvar;i++) {
      delete var[i];
    }
    // -----------------------------------------------------------------------------------

    //    if(!opt.Contains("loop")) { 

    cout << Form("\n 6. Preparing the graphs and plotting .. ") << endl;

    // Centering in the x1 mean
    if(opt.Contains("zmean")) {
      hX1->SetBins(x1Nbin,x1Min-zmean,x1Max-zmean);
      hP1X1->SetBins(x1Nbin,x1Min-zmean,x1Max-zmean,p1Nbin,p1Min,p1Max);
      hX2X1->SetBins(x1Nbin,x1Min-zmean,x1Max-zmean,x2Nbin,x2Min,x2Max);
      if(pData->Is3D()) 
	hX3X1->SetBins(x1Nbin,x1Min-zmean,x1Max-zmean,x3Nbin,x3Min,x3Max);
   
      for(Int_t k=0;k<SNbin;k++) 
	zbin[k] -= zmean;
      
      zmean = 0.0;
      
    }
    // ------

    // Create the graph with the sliced quantities:
    gEmitx = new TGraph(SNbin,zbin,semitx);
    gXrms = new TGraph(SNbin,zbin,sx_rms);
    gErms = new TGraph(SNbin,zbin,sErms);

    if(pData->Is3D()) {
      gEmity = new TGraph(SNbin,zbin,semity);
      gYrms = new TGraph(SNbin,zbin,sy_rms);
    }

    // Profile energy for p1 vs x1:
    // TString pname = hP1X1->GetName();
    // pname += "_pfx";
    // TProfile *hP1X1prof = (TProfile*) gROOT->FindObject(pname.Data());
    // if(hP1X1prof) { delete hP1X1prof; hP1X1prof = NULL; }
    // hP1X1prof = hP1X1->ProfileX("_pfx",1,-1,"s");

    // get the errors from the profile:
    // TGraph *gErmsB = NULL;

    // if(hP1X1prof) {
    //   Int_t NP1X1Bins = hP1X1prof->GetNbinsX();
    //   Double_t *x1bins = new Double_t[NP1X1Bins];
    //   Double_t *eRms   = new Double_t[NP1X1Bins];
    //   for(Int_t i=1;i<=hP1X1prof->GetNbinsX();i++) {
    // 	x1bins[i] = hP1X1prof->GetBinCenter(i);
    // 	eRms[i] = 100 * hP1X1prof->GetBinError(i) / hP1X1prof->GetBinContent(i);
    //   }
    //   gErmsB = new TGraph(NP1X1Bins,x1bins,eRms);
    // }
    

    // Vertical Energy histogram:
    // --------------------------------------------------------------------------------   
    TGraph *gP1left = NULL;
    if(hP1) {
      Double_t *yarray   = new Double_t[p1Nbin];
      Double_t *xarray   = new Double_t[p1Nbin];

      // This is for the right side:
      // Double_t xMax = x1Min + (x1Max-x1Min) * 0.9;
      // Double_t xMin = x1Max;
      // And this for left:
      Double_t xMin = hX1->GetXaxis()->GetXmin();
      Double_t xMax = hX1->GetXaxis()->GetXmin() + 
	(hX1->GetXaxis()->GetXmax()-hX1->GetXaxis()->GetXmin()) * 0.2;
      Double_t EneMax = hP1->GetMaximum();
      // cout << Form("  EneMax = %f ", EneMax) << endl;

      for(Int_t j=0; j<p1Nbin; j++) {
	yarray[j] = hP1->GetBinCenter(j+1);
	xarray[j] = ((xMax-xMin)/EneMax)*hP1->GetBinContent(j+1) + xMin;

	// cout << Form("  x = %f  y = %f ", xarray[j],yarray[j]) << endl;
      }

      gP1left = new TGraph(p1Nbin,xarray,yarray);
      gP1left->SetLineColor(PGlobals::elecLine);
      gP1left->SetLineWidth(2);
      gP1left->SetFillStyle(1001);
      gP1left->SetFillColor(PGlobals::elecFill);

      delete yarray;
      delete xarray;

      // Ranges!!
      Double_t yMin =  999.9;
      Double_t yMax =  -999.9;

      for(Int_t k=0;k<SNbin;k++) {
	if(semitx[k]<yMin)
	  yMin = semitx[k];

	if(semitx[k]>yMax)
	  yMax = semitx[k];

	if(pData->Is3D()) {
	  if(semity[k]<yMin)
	    yMin = semity[k];
	  
	  if(semity[k]>yMax)
	    yMax = semity[k];
	}
	
	if(sErms[k]<yMin)
	  yMin = sErms[k];

	if(sErms[k]>yMax)
	  yMax = sErms[k];
      }

      for(Int_t k=1;k<=x1Nbin;k++) {
	Double_t value = hX1->GetBinContent(k);

	if(value<yMin)
	  yMin = value;

	if(value>yMax)
	  yMax = value;

      }

      yMax *= 1.1;

      // Plotting
      // -----------------------------------------------

      // Canvas setup
      // Create the canvas and the pads before the Frame loop
      // Resolution:
      Int_t sizex = 800;
      Int_t sizey = 600;

      char cName[32];
      sprintf(cName,"C");     
      TCanvas *C = (TCanvas*) gROOT->FindObject(cName);
      if(C==NULL) C = new TCanvas("C","Longitudinal phasespace",sizex,sizey);
      C->SetFillStyle(4000);
      C->cd();
      C->Clear();

      // Set palette:
      PPalette * pPalette = (PPalette*) gROOT->FindObject("electron0");
      pPalette->cd();

      // Double_t Max  = hP1X1->GetMaximum();
      // Double_t Min  = hP1X1->GetMinimum();

      // hP1X1->GetZaxis()->SetRangeUser(Min,Max); 

      // Text objects
      TPaveText *textTime =  new TPaveText(0.55,0.76,0.80,0.86,"NDC");
      PGlobals::SetPaveTextStyle(textTime,32); 
      textTime->SetTextColor(kGray+2);
      char ctext[128];
      if(opt.Contains("units") && pData->GetPlasmaDensity()) 
	// sprintf(ctext,"z = %5.1f um", Time);
	sprintf(ctext,"z = %5.1f %s", Time, spaSUnit.c_str());
      else
	sprintf(ctext,"k_{p}z = %5.1f",Time);
      textTime->AddText(ctext);
 
      TPaveText *textDen = new TPaveText(0.15,0.85,0.48,0.9,"NDC");
      PGlobals::SetPaveTextStyle(textDen,12); 
      textDen->SetTextColor(kOrange+10);
      if(opt.Contains("units") && pData->GetPlasmaDensity())
	sprintf(ctext,"n_{0} = %.2f %s", n0 / denUnit,denSUnit.c_str());
      else if(pData->GetBeamDensity() && pData->GetPlasmaDensity())
	sprintf(ctext,"n_{b}/n_{0} = %5.2f", pData->GetBeamDensity()/n0);
      textDen->AddText(ctext);

      TPaveText *textCharge = new TPaveText(0.15,0.25,0.48,0.3,"NDC");
      PGlobals::SetPaveTextStyle(textCharge,12); 
      textCharge->SetTextColor(kGray+2);
      if(opt.Contains("units") && pData->GetPlasmaDensity())
	sprintf(ctext,"Q = %5.2f %s", Charge,chargeSUnit.c_str());
      else
	sprintf(ctext,"Q = %5.2f n0#timeskp^{-3}", Charge);    
      textCharge->AddText(ctext);

      TPaveText *textMom = new TPaveText(0.55,0.035,0.80,0.14,"NDC");
      PGlobals::SetPaveTextStyle(textMom,32); 
      textMom->SetTextColor(kGray+3);
      textMom->SetTextFont(62);
      if(opt.Contains("units") && pData->GetPlasmaDensity())
	sprintf(ctext,"#LTp_{z}#GT = %5.2f %s/c", pzmeanFWHM, eneSUnit.c_str());
      else
	sprintf(ctext,"#LTp_{z}#GT = %5.2f mc", pzmeanFWHM);    
      textMom->AddText(ctext);


      TPaveText *textInfo = new TPaveText(0.55,0.35,0.80,0.75,"NDC");
      PGlobals::SetPaveTextStyle(textInfo,32); 
      textInfo->SetTextColor(kGray+2);
      textInfo->SetTextFont(42);
      sprintf(ctext,"Q = %5.2f %s",Charge,chargeSUnit.c_str());
      textInfo->AddText(ctext);
      sprintf(ctext,"#Delta#zeta = %5.2f %s",zrms,spaSUnit.c_str());
      textInfo->AddText(ctext);
      // sprintf(ctext,"#Delta#gamma/#LT#gamma#GT = %4.1f %s",(pzrmsFWHM/pzmeanFWHM)/ermsUnit,ermsSUnit.c_str());
      sprintf(ctext,"#Delta#gamma/#LT#gamma#GT = %4.1f %s",(pzrms/pzmean)/ermsUnit,ermsSUnit.c_str());
      textInfo->AddText(ctext);
      sprintf(ctext,"#varepsilon_{n,x} = %5.2f %s",emitx,emitSUnit.c_str());
      textInfo->AddText(ctext);
      if(pData->Is3D()) {
	sprintf(ctext,"#varepsilon_{n,y} = %5.2f %s",emity,emitSUnit.c_str());
	textInfo->AddText(ctext);
      }
      
      // Setup Pad layout: 
      const Int_t NPad = 2;
      TPad *pad[NPad];
      TH1F *hFrame[NPad];
      TString sLabels[] = {"(b)","(a)"};
      TPaveText **textLabel = new TPaveText*[NPad];

      Double_t lMargin = 0.15;
      Double_t rMargin = 0.18;
      Double_t bMargin = 0.15;
      Double_t tMargin = 0.04;
      Double_t factor = 1.0;    
      PGlobals::CanvasAsymPartition(C,NPad,lMargin,rMargin,bMargin,tMargin,factor,0.00);

      // Define the frames for plotting
      Int_t fonttype = 43;
      Int_t fontsize = 24;
      Int_t tfontsize = 28;
      Double_t txoffset = 2.0;
      Double_t lxoffset = 0.02;
      Double_t tyoffset = 1.3;
      Double_t lyoffset = 0.01;
      Double_t tylength = 0.02;
      Double_t txlength = 0.04;
      for(Int_t i=0;i<NPad;i++) {
	char name[16];
	sprintf(name,"pad_%i",i);
	pad[i] = (TPad*) gROOT->FindObject(name);
	pad[i]->SetFrameLineWidth(2);  
	pad[i]->SetTickx(0);
	pad[i]->SetTicky(0);
       	if(opt.Contains("trans"))
	  pad[i]->SetFillStyle(4000);
	pad[i]->SetFrameFillStyle(4000);


	sprintf(name,"hFrame_%i",i);
	hFrame[i] = (TH1F*) gROOT->FindObject(name);
	if(hFrame[i]) delete hFrame[i];
	hFrame[i] = (TH1F*) hX1->Clone(name);
	hFrame[i]->Reset();

	Double_t xFactor = pad[0]->GetAbsWNDC()/pad[i]->GetAbsWNDC();
	Double_t yFactor = pad[0]->GetAbsHNDC()/pad[i]->GetAbsHNDC();

	// Format for y axis
	hFrame[i]->GetYaxis()->SetTitleFont(fonttype);
	hFrame[i]->GetYaxis()->SetTitleSize(tfontsize);
	hFrame[i]->GetYaxis()->SetTitleOffset(tyoffset);
	hFrame[i]->GetYaxis()->SetLabelFont(fonttype);
	hFrame[i]->GetYaxis()->SetLabelSize(fontsize);
	hFrame[i]->GetYaxis()->SetLabelOffset(lyoffset);
	hFrame[i]->GetYaxis()->SetTickLength(xFactor*tylength/yFactor);
	hFrame[i]->GetYaxis()->CenterTitle();

	// Format for x axis
	hFrame[i]->GetXaxis()->SetTitleFont(fonttype);
	hFrame[i]->GetXaxis()->SetTitleSize(tfontsize+2);
	hFrame[i]->GetXaxis()->SetTitleOffset(txoffset);
	hFrame[i]->GetXaxis()->SetLabelFont(fonttype);
	hFrame[i]->GetXaxis()->SetLabelSize(fontsize+2);
	hFrame[i]->GetXaxis()->SetLabelOffset(lxoffset);
	hFrame[i]->GetXaxis()->CenterTitle();
	hFrame[i]->GetXaxis()->SetTickLength(yFactor*txlength/xFactor);      
      }

      C->cd(0);
      pad[1]->Draw();
      pad[1]->cd(); // <---------------------------------------------- Top Plot ---------

      if(opt.Contains("logz")) {
	gPad->SetLogz(1);
      } else {
	gPad->SetLogz(0);
      }

      hFrame[1]->GetYaxis()->SetRangeUser(hP1X1->GetYaxis()->GetXmin(),hP1X1->GetYaxis()->GetXmax());

      if(opt.Contains("units"))
	hFrame[1]->GetYaxis()->SetTitle("p_{z} [GeV/c]");

      hFrame[1]->Draw();

      gP1left->SetLineWidth(2);
      if(!opt.Contains("nospec")) {
	gP1left->Draw("F");
	gP1left->Draw("L");
      }
      
      TLine lZmean(zmean,hP1X1->GetYaxis()->GetXmin(),zmean,hP1X1->GetYaxis()->GetXmax());
      lZmean.SetLineColor(kGray+2);
      lZmean.SetLineStyle(2);
      lZmean.Draw();

      TLine lPmean(hP1X1->GetXaxis()->GetXmin(),pzmeanFWHM,hP1X1->GetXaxis()->GetXmax(),pzmeanFWHM);
      lPmean.SetLineColor(kGray+2);
      lPmean.SetLineStyle(2);
      lPmean.Draw();

      // lines indicating the energy interval
      // TLine pzminline(hP1X1->GetXaxis()->GetXmin(),pzmin, ((xMax-xMin)/EneMax)*(EneMax/2) + xMin,pzmin);
      TLine pzminline(hP1X1->GetXaxis()->GetXmin(),pzmin,hP1X1->GetXaxis()->GetXmax(),pzmin);
      pzminline.SetLineColor(kGray+1);
      pzminline.SetLineStyle(3);
      pzminline.Draw();
      
      TLine pzmaxline(hP1X1->GetXaxis()->GetXmin(),pzmax,hP1X1->GetXaxis()->GetXmax(),pzmax);
      pzmaxline.SetLineColor(kGray+1);
      pzmaxline.SetLineStyle(3);
      pzmaxline.Draw();

      // 2D histogram z range
      // Double_t dmax = hP1X1->GetMaximum();
      // Double_t dmin = 0.0;
      // hP1X1->GetZaxis()->SetRangeUser(dmin,dmax);

      hP1X1->GetZaxis()->SetTitleFont(fonttype);
      hP1X1->GetZaxis()->SetTickLength(0.01);
      
      hP1X1->Draw("colzsame");
      // hP1X1->SetContour(20);
      // hP1X1->Draw("contzsame");
      // hP1X1prof->SetMarkerStyle(1);
      // hP1X1prof->SetLineWidth(2);
      // hP1X1prof->Draw("zsame");

      //hP1->Draw("C");

      gPad->Update();

      TPaletteAxis *palette = (TPaletteAxis*)hP1X1->GetListOfFunctions()->FindObject("palette");
      if(palette) {
	Double_t y1 = gPad->GetBottomMargin();
	Double_t y2 = 1 - gPad->GetTopMargin();
	Double_t x1 = 1 - gPad->GetRightMargin();
	palette->SetY2NDC(y2 - 0.04);
	palette->SetY1NDC(y1 + 0.04);
	palette->SetX1NDC(x1 + 0.01);
	palette->SetX2NDC(x1 + 0.04);
	palette->SetTitleOffset(tyoffset);
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


      if(opt.Contains("notext")) {
	cout << opt;
      } else if(opt.Contains("noinfo")) {
	textTime->Draw();
	textMom->Draw();
      } else {
	textInfo->Draw();
	textTime->Draw();
	textMom->Draw();
      }


      Double_t y1 = gPad->GetBottomMargin();
      Double_t y2 = 1 - gPad->GetTopMargin();
      Double_t x1 = gPad->GetLeftMargin();
      Double_t x2 = 1 - gPad->GetRightMargin();
      Double_t yrange = y2-y1; 
      Double_t xrange = x2-x1; 

      textLabel[1] = new TPaveText(x1 + 0.02*(x2-x1), y2-0.2*(y2-y1), x1+0.30*(x2-x1), y2-0.05*(y2-y1),"NDC");
      PGlobals::SetPaveTextStyle(textLabel[1],12); 
      textLabel[1]->SetTextFont(42);
      textLabel[1]->AddText(sLabels[1]);
      // textLabel[1]->Draw();

      gPad->RedrawAxis(); 


      // Bottom plot -----------------------------------------
      C->cd(0);
      pad[0]->Draw();
      pad[0]->cd(); // <---------------------------------------------------------- Bottom Plot
      // if(opt.Contains("logz")) {
      //   pad[0]->SetLogz(1);
      // } else {
      //   pad[0]->SetLogz(0);
      // }

      TLegend *Leg;
      Leg=new TLegend(0.55,0.72,1 - 0.5*gPad->GetRightMargin() - 0.02,0.95);

      PGlobals::SetPaveStyle(Leg);
      Leg->SetTextAlign(12);
      Leg->SetTextColor(kGray+3);
      Leg->SetTextFont(42);
      Leg->SetLineColor(1);
      Leg->SetBorderSize(0);
      Leg->SetFillColor(0);
      Leg->SetFillStyle(1001);
      Leg->SetFillStyle(0); // Hollow

      char sleg[16];
      sprintf(sleg,"Current [%s]",curSUnit.c_str());
      Leg->AddEntry(hX1  ,sleg,"L");	
      //sprintf(sleg,"Energy spread [%s]",ermsUnit.c_str());
      sprintf(sleg,"E. spread [%s]",ermsSUnit.c_str());
      Leg->AddEntry(gErms,sleg,"PL");
      //      sprintf(sleg,"Emittance [%s]",emitUnit.c_str());
      sprintf(sleg,"Emitt. x [%s]",emitSUnit.c_str());
      Leg->AddEntry(gEmitx,sleg,"PL");
      if(pData->Is3D()) {
	sprintf(sleg,"Emitt. y [%s]",emitSUnit.c_str());
	Leg->AddEntry(gEmity,sleg,"PL");
      }
      //Leg->AddEntry(gXrms,"Bunch width [#mum]","PL");


      hFrame[0]->GetYaxis()->SetRangeUser(0.0,1.1*yMax);
      hFrame[0]->Draw();

      hX1->GetYaxis()->SetNdivisions(503);
      hX1->SetLineWidth(2);
      hX1->SetFillStyle(1001);
      hX1->SetFillColor(PGlobals::elecFill);
      //hX1->Smooth();
      hX1->Draw("FL same");

      TLine lZmean2(zmean,0.0,zmean,1.1*yMax);
      lZmean2.SetLineColor(kGray+2);
      lZmean2.SetLineStyle(2);
      lZmean2.Draw();

      Size_t markerSize = 1.0; 
      Width_t lineWidth  = 2.0;   

      gXrms->SetMarkerStyle(20);
      gXrms->SetLineStyle(1);
      gXrms->SetMarkerColor(kGray+1);
      gXrms->SetMarkerSize(markerSize); 
      gXrms->SetLineColor(kGray+1);
      gXrms->SetLineWidth(lineWidth);
      //gXrms->Draw("PL");

      if(pData->Is3D()) {
	gEmity->SetMarkerStyle(20);
	gEmity->SetMarkerColor(kGray+2);
	gEmity->SetMarkerSize(markerSize);
	gEmity->SetLineWidth(lineWidth);
	gEmity->SetLineColor(kGray+2);
	gEmity->Draw("PL");
      }

      gEmitx->SetMarkerStyle(20);
      gEmitx->SetMarkerColor(kGray+3);
      gEmitx->SetMarkerSize(markerSize);
      gEmitx->SetLineWidth(lineWidth);
      gEmitx->SetLineColor(kGray+3);
      gEmitx->Draw("PL");


      gErms->SetMarkerStyle(20);
      gErms->SetMarkerSize(markerSize);
      gErms->SetMarkerColor(kOrange+10);
      gErms->SetLineColor(kOrange+10);
      gErms->SetLineWidth(lineWidth);
      gErms->Draw("PL");

      if(opt.Contains("bw")) {
	gErms->SetMarkerStyle(21);
	gErms->SetMarkerSize(markerSize-0.2);
      }
   
      Leg->Draw();

      y1 = gPad->GetBottomMargin();
      y2 = 1 - gPad->GetTopMargin();
      x1 = gPad->GetLeftMargin();
      x2 = 1 - gPad->GetRightMargin();
      yrange = y2-y1; 
      xrange = x2-x1; 

      textLabel[0] = new TPaveText(x1 + 0.02*(x2-x1), y2-0.2*(y2-y1), x1+0.30*(x2-x1), y2-0.05*(y2-y1),"NDC");
      PGlobals::SetPaveTextStyle(textLabel[0],12); 
      textLabel[0]->SetTextFont(42);
      textLabel[0]->AddText(sLabels[0]);
      // textLabel[0]->Draw();

      gPad->RedrawAxis(); 

      gPad->Update();


      // Print to file --------------------------------------

      C->cd();
      
      // Output file
      TString fOutName = Form("./%s/Plots/Bunch/%s/Bunch-%s",sim.Data(),pData->GetSpeciesName(index).c_str(),sim.Data());

      TString fOutNamep1x1 = fOutName + Form("-%s_%i","p1x1",time);
      PGlobals::imgconv(C,fOutNamep1x1,opt);

      // delete [] pad;
      // delete [] hFrame;
      for(Int_t i=0;i<NPad;i++) {
	char name[16];
	sprintf(name,"pad_%i",i);
	pad[i] = (TPad*) gROOT->FindObject(name);
	if(pad[i]) delete pad[i];

	sprintf(name,"hFrame_%i",i);
	hFrame[i] = (TH1F*) gROOT->FindObject(name);
	if(hFrame[i]) delete hFrame[i];
      }
      // ---------------------------------------------------------
    
      // Transverse phasespace

      sprintf(cName,"C1");     
      TCanvas *C1 = (TCanvas*) gROOT->FindObject(cName);
      if(C1==NULL) C1 = new TCanvas("C1","Space x2-x1 and x3-x1",sizex,sizey);
      C1->cd();
      C1->Clear();

      Int_t NPad1 = 2;
      lMargin = 0.15;
      rMargin = 0.18;
      bMargin = 0.15;
      tMargin = 0.04;
      factor = 1.0;  
      txoffset = 2.0;  
      if(!pData->Is3D()) {
      	NPad1 = 1;
	txoffset = 1.2;  
      }
      PGlobals::CanvasAsymPartition(C1,NPad1,lMargin,rMargin,bMargin,tMargin,factor);
      
      for(Int_t i=0;i<NPad1;i++) {
	char name[16];
	sprintf(name,"pad_%i",i);
	pad[i] = (TPad*) gROOT->FindObject(name);
	pad[i]->SetFrameLineWidth(2);  
	pad[i]->SetTickx(1);
	pad[i]->SetTicky(1);

	sprintf(name,"hFrame_%i",i);
	hFrame[i] = (TH1F*) gROOT->FindObject(name);
	if(hFrame[i]) delete hFrame[i];
	if(i==0)
	  hFrame[i] = (TH1F*) hX2X1->Clone(name);
	else
	  hFrame[i] = (TH1F*) hX3X1->Clone(name);
      
	hFrame[i]->Reset();

	Double_t xFactor = pad[0]->GetAbsWNDC()/pad[i]->GetAbsWNDC();
	Double_t yFactor = pad[0]->GetAbsHNDC()/pad[i]->GetAbsHNDC();

	// Format for y axis
	hFrame[i]->GetYaxis()->SetTitleFont(fonttype);
	hFrame[i]->GetYaxis()->SetTitleSize(tfontsize);
	hFrame[i]->GetYaxis()->SetTitleOffset(tyoffset);
	hFrame[i]->GetYaxis()->SetLabelFont(fonttype);
	hFrame[i]->GetYaxis()->SetLabelSize(fontsize);
	hFrame[i]->GetYaxis()->SetLabelOffset(lyoffset);
	hFrame[i]->GetYaxis()->SetTickLength(xFactor*tylength/yFactor);
	hFrame[i]->GetYaxis()->CenterTitle();

	// Format for x axis
	hFrame[i]->GetXaxis()->SetTitleFont(fonttype);
	hFrame[i]->GetXaxis()->SetTitleSize(tfontsize+2);
	hFrame[i]->GetXaxis()->SetTitleOffset(txoffset);
	hFrame[i]->GetXaxis()->SetLabelFont(fonttype);
	hFrame[i]->GetXaxis()->SetLabelSize(fontsize+2);
	hFrame[i]->GetXaxis()->SetLabelOffset(lxoffset);
	hFrame[i]->GetXaxis()->CenterTitle();
	hFrame[i]->GetXaxis()->SetTickLength(yFactor*txlength/xFactor);      
      }


      C1->cd(0);
      pad[0]->Draw();
      pad[0]->cd(); 

      if(opt.Contains("logz")) {
	gPad->SetLogz(1);
      } else {
	gPad->SetLogz(0);
      }

      hFrame[0]->GetYaxis()->SetRangeUser(hX2X1->GetYaxis()->GetXmin(),hX2X1->GetYaxis()->GetXmax());
      hFrame[0]->GetYaxis()->SetTitle(hX2X1->GetYaxis()->GetTitle());
      hFrame[0]->GetXaxis()->SetRangeUser(hX2X1->GetXaxis()->GetXmin(),hX2X1->GetXaxis()->GetXmax());
      hFrame[0]->GetXaxis()->SetTitle(hX2X1->GetXaxis()->GetTitle());
      hFrame[0]->Draw("axis");

      TLine lX1mean(zmean,hFrame[0]->GetYaxis()->GetXmin(),zmean,hFrame[0]->GetYaxis()->GetXmax());
      lX1mean.SetLineColor(kGray+2);
      lX1mean.SetLineStyle(2);
      lX1mean.Draw();

      TLine lX2mean(hFrame[0]->GetXaxis()->GetXmin(),x_mean,hFrame[0]->GetXaxis()->GetXmax(),x_mean);
      lX2mean.SetLineColor(kGray+2);
      lX2mean.SetLineStyle(2);
      lX2mean.Draw();


      TH2F *hX2X1cl = (TH2F*) hX2X1->Clone("hX2X1cl");
      hX2X1cl->Draw("colz same");
      hX2X1cl->GetZaxis()->CenterTitle();
      hX2X1cl->GetZaxis()->SetTitleFont(fonttype);
      Float_t xFactor = pad[0]->GetAbsWNDC()/pad[0]->GetAbsWNDC();
      Float_t yFactor = pad[0]->GetAbsHNDC()/pad[0]->GetAbsHNDC();
      hX2X1cl->GetZaxis()->SetTickLength(xFactor*tylength/yFactor);
  
      gPad->Update();
    
      y1 = gPad->GetBottomMargin();
      y2 = 1 - gPad->GetTopMargin();
      x1 = gPad->GetLeftMargin();
      x2 = 1 - gPad->GetRightMargin();
      yrange = y2-y1; 
      xrange = x2-x1; 

      palette = (TPaletteAxis*)hX2X1cl->GetListOfFunctions()->FindObject("palette");
      if(palette) {
	Double_t y1 = gPad->GetBottomMargin();
	Double_t y2 = 1 - gPad->GetTopMargin();
	Double_t x1 = 1 - gPad->GetRightMargin();
	palette->SetY2NDC(y2 - 0.04);
	palette->SetY1NDC(y1 + 0.04);
	palette->SetX1NDC(x1 + 0.01);
	palette->SetX2NDC(x1 + 0.04);
	palette->SetTitleOffset(tyoffset);
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

      TPaveText *textInfoX2X1 = new TPaveText(x1+0.02*xrange,y2-0.30*yrange,x1+0.20*xrange,y2-0.05*yrange,"NDC");
      PGlobals::SetPaveTextStyle(textInfoX2X1,12); 
      textInfoX2X1->SetTextColor(kGray+3);
      textInfoX2X1->SetTextFont(42);

      char text[64];
      sprintf(text,"Q = %5.1f %s",Charge,chargeSUnit.c_str());
      textInfoX2X1->AddText(text);
      sprintf(text,"#Delta#zeta = %5.2f %s",zrms,spaSUnit.c_str());
      textInfoX2X1->AddText(text);
      sprintf(text,"#Deltax = %5.2f %s",x_rms,tspaSUnit.c_str());
      textInfoX2X1->AddText(text);
      // sprintf(text,"#varepsilon_{x} = %5.2f #mum",emitx);
      // textInfoX2X1->AddText(text);
      // sprintf(text,"#beta_{x} = %5.2f mm",1E-3*betax);
      // textInfoX2X1->AddText(text);
      textInfoX2X1->Draw();
    
      C1->cd(0);

      if(pData->Is3D()) {
	pad[1]->Draw();
	pad[1]->cd(); 
	
	if(opt.Contains("logz")) {
	  gPad->SetLogz(1);
	} else {
	  gPad->SetLogz(0);
	}
	
	hFrame[1]->GetYaxis()->SetRangeUser(hX3X1->GetYaxis()->GetXmin(),hX3X1->GetYaxis()->GetXmax());
	hFrame[1]->GetYaxis()->SetTitle(hX3X1->GetYaxis()->GetTitle());
	hFrame[1]->GetXaxis()->SetRangeUser(hX3X1->GetXaxis()->GetXmin(),hX3X1->GetXaxis()->GetXmax());
	hFrame[1]->GetXaxis()->SetTitle(hX3X1->GetXaxis()->GetTitle());
	hFrame[1]->Draw("axis");
	
	TLine lX1mean2(zmean,hFrame[1]->GetYaxis()->GetXmin(),zmean,hFrame[1]->GetYaxis()->GetXmax());
	lX1mean2.SetLineColor(kGray+2);
	lX1mean2.SetLineStyle(2);
	lX1mean2.Draw();
	
	TLine lX3mean(hFrame[1]->GetXaxis()->GetXmin(),y_mean,hFrame[1]->GetXaxis()->GetXmax(),y_mean);
	lX3mean.SetLineColor(kGray+2);
	lX3mean.SetLineStyle(2);
	lX3mean.Draw();
	
	
	TH2F *hX3X1cl = (TH2F*) hX3X1->Clone("hX3X1cl");
	hX3X1cl->GetZaxis()->CenterTitle();
	hX3X1cl->GetZaxis()->SetTitleFont(fonttype);
	xFactor = pad[0]->GetAbsWNDC()/pad[1]->GetAbsWNDC();
	yFactor = pad[0]->GetAbsHNDC()/pad[1]->GetAbsHNDC();
	hX3X1cl->GetZaxis()->SetTickLength(xFactor*tylength/yFactor);
	hX3X1cl->Draw("colz same");
	
	gPad->Update();
	
	y1 = gPad->GetBottomMargin();
	y2 = 1 - gPad->GetTopMargin();
	x1 = gPad->GetLeftMargin();
	x2 = 1 - gPad->GetRightMargin();
	yrange = y2-y1; 
	xrange = x2-x1; 
	
	palette = (TPaletteAxis*)hX3X1cl->GetListOfFunctions()->FindObject("palette");
	if(palette) {
	  Double_t y1 = gPad->GetBottomMargin();
	  Double_t y2 = 1 - gPad->GetTopMargin();
	  Double_t x1 = 1 - gPad->GetRightMargin();
	  palette->SetY2NDC(y2 - 0.04);
	  palette->SetY1NDC(y1 + 0.04);
	  palette->SetX1NDC(x1 + 0.01);
	  palette->SetX2NDC(x1 + 0.04);
	  palette->SetTitleOffset(tyoffset);
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

	TPaveText *textInfoX3X1 =  new TPaveText(x1+0.02*xrange,y2-0.30*yrange,x1+0.20*xrange,y2-0.05*yrange,"NDC");
	PGlobals::SetPaveTextStyle(textInfoX3X1,12); 
	textInfoX3X1->SetTextColor(kGray+3);
	textInfoX3X1->SetTextFont(42);
	
	sprintf(text,"Q = %5.1f %s",Charge,chargeSUnit.c_str());
	textInfoX3X1->AddText(text);
	sprintf(text,"#Delta#zeta = %5.2f %s",zrms,spaSUnit.c_str());
	textInfoX3X1->AddText(text);
	sprintf(text,"#Deltay = %5.2f %s",y_rms,tspaSUnit.c_str());
	textInfoX3X1->AddText(text);
	// sprintf(text,"#varepsilon_{x} = %5.2f #mum",emitx);
	// textInfoX3X1->AddText(text);
	// sprintf(text,"#beta_{x} = %5.2f mm",1E-3*betax);
	// textInfoX3X1->AddText(text);
	textInfoX3X1->Draw();
      }
      
      C1->cd();
    
      TString fOutName1 =  fOutName + Form("-%s_%i","x2x3x1",time);
      if(!pData->Is3D()) fOutName1 = fOutName + Form("-%s_%i","x2x1",time);
      
      PGlobals::imgconv(C1,fOutName1,opt);


      // Transverse phasespace
      // x3-x2

      if(pData->Is3D()) {
      
	sprintf(cName,"CX3X2");     
	TCanvas *CX3X2 = (TCanvas*) gROOT->FindObject(cName);
	if(CX3X2==NULL) CX3X2 = new TCanvas("CX3X2","Transverse space x3-x2",sizex,sizey);
	CX3X2->cd();
	CX3X2->Clear();

	const Int_t NPadX3X2 = 1;
	lMargin = 0.15;
	rMargin = 0.18;
	bMargin = 0.15;
	tMargin = 0.04;
	factor = 1.0;  
	txoffset = 1.2;  
	PGlobals::CanvasAsymPartition(CX3X2,NPadX3X2,lMargin,rMargin,bMargin,tMargin,factor);

	for(Int_t i=0;i<NPadX3X2;i++) {
	  char name[16];
	  sprintf(name,"pad_%i",i);
	  pad[i] = (TPad*) gROOT->FindObject(name);
	  pad[i]->SetFrameLineWidth(2);  
	  pad[i]->SetTickx(1);
	  pad[i]->SetTicky(1);

	  sprintf(name,"hFrame_%i",i);
	  hFrame[i] = (TH1F*) gROOT->FindObject(name);
	  if(hFrame[i]) delete hFrame[i];
	  hFrame[i] = (TH1F*) hX3X2->Clone(name);
	  hFrame[i]->Reset();

	  Double_t xFactor = pad[0]->GetAbsWNDC()/pad[i]->GetAbsWNDC();
	  Double_t yFactor = pad[0]->GetAbsHNDC()/pad[i]->GetAbsHNDC();

	  // Format for y axis
	  hFrame[i]->GetYaxis()->SetTitleFont(fonttype);
	  hFrame[i]->GetYaxis()->SetTitleSize(tfontsize);
	  hFrame[i]->GetYaxis()->SetTitleOffset(tyoffset);
	  hFrame[i]->GetYaxis()->SetLabelFont(fonttype);
	  hFrame[i]->GetYaxis()->SetLabelSize(fontsize);
	  hFrame[i]->GetYaxis()->SetLabelOffset(lyoffset);
	  hFrame[i]->GetYaxis()->SetTickLength(xFactor*tylength/yFactor);
	  hFrame[i]->GetYaxis()->CenterTitle();

	  // Format for x axis
	  hFrame[i]->GetXaxis()->SetTitleFont(fonttype);
	  hFrame[i]->GetXaxis()->SetTitleSize(tfontsize+2);
	  hFrame[i]->GetXaxis()->SetTitleOffset(txoffset);
	  hFrame[i]->GetXaxis()->SetLabelFont(fonttype);
	  hFrame[i]->GetXaxis()->SetLabelSize(fontsize+2);
	  hFrame[i]->GetXaxis()->SetLabelOffset(lxoffset);
	  hFrame[i]->GetXaxis()->CenterTitle();
	  hFrame[i]->GetXaxis()->SetTickLength(yFactor*txlength/xFactor);      
	}


	CX3X2->cd(0);
	pad[0]->Draw();
	pad[0]->cd(); 

	if(opt.Contains("logz")) {
	  gPad->SetLogz(1);
	} else {
	  gPad->SetLogz(0);
	}

	hFrame[0]->GetYaxis()->SetRangeUser(hX3X2->GetYaxis()->GetXmin(),hX3X2->GetYaxis()->GetXmax());
	hFrame[0]->GetYaxis()->SetTitle(hX3X2->GetYaxis()->GetTitle());
	hFrame[0]->GetXaxis()->SetRangeUser(hX3X2->GetXaxis()->GetXmin(),hX3X2->GetXaxis()->GetXmax());
	hFrame[0]->GetXaxis()->SetTitle(hX3X2->GetXaxis()->GetTitle());
	hFrame[0]->Draw("axis");

	TLine lX2mean2(x_mean,hFrame[0]->GetYaxis()->GetXmin(),x_mean,hFrame[0]->GetYaxis()->GetXmax());
	lX2mean2.SetLineColor(kGray+2);
	lX2mean2.SetLineStyle(2);
	lX2mean2.Draw();

	TLine lX3mean2(hFrame[0]->GetXaxis()->GetXmin(),y_mean,hFrame[0]->GetXaxis()->GetXmax(),y_mean);
	lX3mean2.SetLineColor(kGray+2);
	lX3mean2.SetLineStyle(2);
	lX3mean2.Draw();


	TH2F *hX3X2cl = (TH2F*) hX3X2->Clone("hX3X2cl");

	hX3X2cl->Draw("colz same");
	hX3X2cl->GetZaxis()->CenterTitle();
    
	y1 = gPad->GetBottomMargin();
	y2 = 1 - gPad->GetTopMargin();
	x1 = gPad->GetLeftMargin();
	x2 = 1 - gPad->GetRightMargin();
	yrange = y2-y1; 
	xrange = x2-x1; 
    
	TPaveText *textStatX3X2 =  new TPaveText(x1+0.02*xrange,y2-0.20*yrange,x1+0.30*xrange,y2-0.05*yrange,"NDC");
	PGlobals::SetPaveTextStyle(textStatX3X2,12); 
	textStatX3X2->SetTextColor(kGray+3);
	textStatX3X2->SetTextFont(42);

	sprintf(text,"Q = %5.1f %s",Charge,chargeSUnit.c_str());
	textStatX3X2->AddText(text);
	sprintf(text,"#Deltax = %5.2f %s",x_rms,tspaSUnit.c_str());
	textStatX3X2->AddText(text);
	sprintf(text,"#Deltay = %5.2f %s",y_rms,tspaSUnit.c_str());
	textStatX3X2->AddText(text);
	textStatX3X2->Draw();
    
    
	TString fOutNameX3X2 = fOutName + Form("-%s_%i","x3x2",time);
	PGlobals::imgconv(CX3X2,fOutNameX3X2,opt);
      }
      
      // p2-x2
      sprintf(cName,"C2");     
      TCanvas *C2 = (TCanvas*) gROOT->FindObject(cName);
      if(C2==NULL) C2 = new TCanvas("C2","Transverse phasespace p2-x2",sizex,sizey);
      C2->cd();
      C2->Clear();


      const Int_t NPad2 = 1;
      lMargin = 0.15;
      rMargin = 0.18;
      bMargin = 0.15;
      tMargin = 0.04;
      factor = 1.0;  
      txoffset = 1.2;  
      PGlobals::CanvasAsymPartition(C2,NPad2,lMargin,rMargin,bMargin,tMargin,factor);

      for(Int_t i=0;i<NPad2;i++) {
	char name[16];
	sprintf(name,"pad_%i",i);
	pad[i] = (TPad*) gROOT->FindObject(name);
	pad[i]->SetFrameLineWidth(2);  
	pad[i]->SetTickx(1);
	pad[i]->SetTicky(1);

	sprintf(name,"hFrame_%i",i);
	hFrame[i] = (TH1F*) gROOT->FindObject(name);
	if(hFrame[i]) delete hFrame[i];
	hFrame[i] = (TH1F*) hP2X2->Clone(name);
	hFrame[i]->Reset();

	Double_t xFactor = pad[0]->GetAbsWNDC()/pad[i]->GetAbsWNDC();
	Double_t yFactor = pad[0]->GetAbsHNDC()/pad[i]->GetAbsHNDC();

	// Format for y axis
	hFrame[i]->GetYaxis()->SetTitleFont(fonttype);
	hFrame[i]->GetYaxis()->SetTitleSize(tfontsize);
	hFrame[i]->GetYaxis()->SetTitleOffset(tyoffset);
	hFrame[i]->GetYaxis()->SetLabelFont(fonttype);
	hFrame[i]->GetYaxis()->SetLabelSize(fontsize);
	hFrame[i]->GetYaxis()->SetLabelOffset(lyoffset);
	hFrame[i]->GetYaxis()->SetTickLength(xFactor*tylength/yFactor);
	hFrame[i]->GetYaxis()->CenterTitle();

	// Format for x axis
	hFrame[i]->GetXaxis()->SetTitleFont(fonttype);
	hFrame[i]->GetXaxis()->SetTitleSize(tfontsize+2);
	hFrame[i]->GetXaxis()->SetTitleOffset(txoffset);
	hFrame[i]->GetXaxis()->SetLabelFont(fonttype);
	hFrame[i]->GetXaxis()->SetLabelSize(fontsize+2);
	hFrame[i]->GetXaxis()->SetLabelOffset(lxoffset);
	hFrame[i]->GetXaxis()->CenterTitle();
	hFrame[i]->GetXaxis()->SetTickLength(yFactor*txlength/xFactor);      
      }


      C2->cd(0);
      pad[0]->Draw();
      pad[0]->cd(); 

      if(opt.Contains("logz")) {
	gPad->SetLogz(1);
      } else {
	gPad->SetLogz(0);
      }

      hFrame[0]->GetYaxis()->SetRangeUser(hP2X2->GetYaxis()->GetXmin(),hP2X2->GetYaxis()->GetXmax());
      hFrame[0]->GetYaxis()->SetTitle(hP2X2->GetYaxis()->GetTitle());
      hFrame[0]->GetXaxis()->SetRangeUser(hP2X2->GetXaxis()->GetXmin(),hP2X2->GetXaxis()->GetXmax());
      hFrame[0]->GetXaxis()->SetTitle(hP2X2->GetXaxis()->GetTitle());
      hFrame[0]->Draw("axis");

      TLine lXmean(x_mean,hFrame[0]->GetYaxis()->GetXmin(),x_mean,hFrame[0]->GetYaxis()->GetXmax());
      lXmean.SetLineColor(kGray+2);
      lXmean.SetLineStyle(2);
      lXmean.Draw();

      TLine lPxmean(hFrame[0]->GetXaxis()->GetXmin(),px_mean,hFrame[0]->GetXaxis()->GetXmax(),px_mean);
      lPxmean.SetLineColor(kGray+2);
      lPxmean.SetLineStyle(2);
      lPxmean.Draw();


      TH2F *hP2X2cl = (TH2F*) hP2X2->Clone("hP2X2cl");
      hP2X2cl->Draw("colz same");
    
      if(opt.Contains("elli") ) 
	for(Int_t k=0;k<SNbin;k++)
	  if(sellipP2X2[k] != NULL)
	    sellipP2X2[k]->Draw();

      TPaveText *textStatInt = new TPaveText(x1+0.02,y2-0.40,x1+0.20,y2-0.05,"NDC");
      PGlobals::SetPaveTextStyle(textStatInt,12); 
      textStatInt->SetTextColor(kGray+3);
      textStatInt->SetTextFont(42);

      sprintf(text,"Q = %5.1f %s",Charge,chargeSUnit.c_str());
      textStatInt->AddText(text);
      sprintf(text,"#Deltax = %5.2f %s",x_rms,tspaSUnit.c_str());
      textStatInt->AddText(text);
      sprintf(text,"#Deltap_{x} = %5.2f %s",px_rms,teneSUnit.c_str());
      textStatInt->AddText(text);
      sprintf(text,"#varepsilon_{x} = %5.2f %s",emitx,emitSUnit.c_str());
      textStatInt->AddText(text);
      sprintf(text,"#beta_{x} = %5.2f %s",betax,betaSUnit.c_str());
      textStatInt->AddText(text);
      textStatInt->Draw();
    
    
      TString fOutName2 = fOutName + Form("-%s_%i","p2x2",time);
      PGlobals::imgconv(C2,fOutName2,opt);

      // p3-x3
      if(pData->Is3D()) {
	sprintf(cName,"C3");     
	TCanvas *C3 = (TCanvas*) gROOT->FindObject(cName);
	if(C3==NULL) C3 = new TCanvas("C3","Transverse phasespace p3-x3",sizex,sizey);
	C3->cd();
	C3->Clear();


	const Int_t NPad3 = 1;
	lMargin = 0.15;
	rMargin = 0.18;
	bMargin = 0.15;
	tMargin = 0.04;
	factor = 1.0;  
	txoffset = 1.2;  
	PGlobals::CanvasAsymPartition(C3,NPad3,lMargin,rMargin,bMargin,tMargin,factor);

	for(Int_t i=0;i<NPad3;i++) {
	  char name[16];
	  sprintf(name,"pad_%i",i);
	  pad[i] = (TPad*) gROOT->FindObject(name);
	  pad[i]->SetFrameLineWidth(2);  
	  pad[i]->SetTickx(1);
	  pad[i]->SetTicky(1);

	  sprintf(name,"hFrame_%i",i);
	  hFrame[i] = (TH1F*) gROOT->FindObject(name);
	  if(hFrame[i]) delete hFrame[i];
	  hFrame[i] = (TH1F*) hP3X3->Clone(name);
	  hFrame[i]->Reset();

	  Double_t xFactor = pad[0]->GetAbsWNDC()/pad[i]->GetAbsWNDC();
	  Double_t yFactor = pad[0]->GetAbsHNDC()/pad[i]->GetAbsHNDC();

	  // Format for y axis
	  hFrame[i]->GetYaxis()->SetTitleFont(fonttype);
	  hFrame[i]->GetYaxis()->SetTitleSize(tfontsize);
	  hFrame[i]->GetYaxis()->SetTitleOffset(tyoffset);
	  hFrame[i]->GetYaxis()->SetLabelFont(fonttype);
	  hFrame[i]->GetYaxis()->SetLabelSize(fontsize);
	  hFrame[i]->GetYaxis()->SetLabelOffset(lyoffset);
	  hFrame[i]->GetYaxis()->SetTickLength(xFactor*tylength/yFactor);
	  hFrame[i]->GetYaxis()->CenterTitle();

	  // Format for x axis
	  hFrame[i]->GetXaxis()->SetTitleFont(fonttype);
	  hFrame[i]->GetXaxis()->SetTitleSize(tfontsize+2);
	  hFrame[i]->GetXaxis()->SetTitleOffset(txoffset);
	  hFrame[i]->GetXaxis()->SetLabelFont(fonttype);
	  hFrame[i]->GetXaxis()->SetLabelSize(fontsize+2);
	  hFrame[i]->GetXaxis()->SetLabelOffset(lxoffset);
	  hFrame[i]->GetXaxis()->CenterTitle();
	  hFrame[i]->GetXaxis()->SetTickLength(yFactor*txlength/xFactor);      
	}


	C3->cd(0);
	pad[0]->Draw();
	pad[0]->cd(); 

	if(opt.Contains("logz")) {
	  gPad->SetLogz(1);
	} else {
	  gPad->SetLogz(0);
	}

	hFrame[0]->GetYaxis()->SetRangeUser(hP3X3->GetYaxis()->GetXmin(),hP3X3->GetYaxis()->GetXmax());
	hFrame[0]->GetYaxis()->SetTitle(hP3X3->GetYaxis()->GetTitle());
	hFrame[0]->GetXaxis()->SetRangeUser(hP3X3->GetXaxis()->GetXmin(),hP3X3->GetXaxis()->GetXmax());
	hFrame[0]->GetXaxis()->SetTitle(hP3X3->GetXaxis()->GetTitle());
	hFrame[0]->Draw("axis");

	TLine lYmean(y_mean,hFrame[0]->GetYaxis()->GetXmin(),y_mean,hFrame[0]->GetYaxis()->GetXmax());
	lYmean.SetLineColor(kGray+2);
	lYmean.SetLineStyle(2);
	lYmean.Draw();

	TLine lPymean(hFrame[0]->GetXaxis()->GetXmin(),py_mean,hFrame[0]->GetXaxis()->GetXmax(),py_mean);
	lPymean.SetLineColor(kGray+2);
	lPymean.SetLineStyle(2);
	lPymean.Draw();


	TH2F *hP3X3cl = (TH2F*) hP3X3->Clone("hP3X3cl");
	hP3X3cl->Draw("colz same");
        
	if(opt.Contains("elli") ) 
	  for(Int_t k=0;k<SNbin;k++)
	    if(sellipP3X3[k] != NULL)
	      sellipP3X3[k]->Draw();

	textStatInt = new TPaveText(x1+0.02,y2-0.40,x1+0.20,y2-0.05,"NDC");
	PGlobals::SetPaveTextStyle(textStatInt,12); 
	textStatInt->SetTextColor(kGray+3);
	textStatInt->SetTextFont(42);

	sprintf(text,"Q = %5.1f %s",Charge,chargeSUnit.c_str());
	textStatInt->AddText(text);
	sprintf(text,"#Deltay = %5.2f %s",y_rms,tspaSUnit.c_str());
	textStatInt->AddText(text);
	sprintf(text,"#Deltap_{y} = %5.2f %s",py_rms,teneSUnit.c_str());
	textStatInt->AddText(text);
	sprintf(text,"#varepsilon_{y} = %5.2f %s",emity,emitSUnit.c_str());
	textStatInt->AddText(text);
	sprintf(text,"#beta_{y} = %5.2f %s",betay,betaSUnit.c_str());
	textStatInt->AddText(text);
	textStatInt->Draw();
    
	TString fOutName3 = fOutName + Form("-%s_%i","p3x3",time);
	PGlobals::imgconv(C3,fOutName3,opt);
      }
      
      // SLICES
      // A pdf document containing all p2-x2 slices:
      
      if(opt.Contains("slcs")) {
      
	gStyle->cd();

	sprintf(cName,"CA4");     
	TCanvas *CA4 = (TCanvas*) gROOT->FindObject(cName);
	if(CA4==NULL) CA4 = new TCanvas("CA4","Sliced p2-x2",600,800);
	CA4->cd();
	CA4->Clear();
	
	Int_t ndiv = 4;
	CA4->Divide(1,ndiv);
	
	TString fOutName2 = Form("./%s/Plots/Bunch/%s/Bunch-%s-slp2x2_%i",sim.Data(),pData->GetSpeciesName(index).c_str(),sim.Data(),time);
	
	CA4->Print(fOutName2 + ".ps[","Portrait");
	
	hP2X2->GetXaxis()->SetLabelSize(0.08);
	hP2X2->GetXaxis()->SetTitleSize(0.08);
	hP2X2->GetXaxis()->SetTitleOffset(1.0);
	hP2X2->GetXaxis()->CenterTitle();
	
	hP2X2->GetYaxis()->SetLabelSize(0.08);
	hP2X2->GetYaxis()->SetTitleSize(0.08);
	hP2X2->GetYaxis()->SetTitleOffset(0.8);
	hP2X2->GetYaxis()->CenterTitle();
	
	hP2X2->GetZaxis()->SetLabelSize(0.08);
	hP2X2->GetZaxis()->SetTitleSize(0.08);
	hP2X2->GetZaxis()->SetTitleOffset(0.8);
	hP2X2->GetZaxis()->CenterTitle();
	
	
	CA4->cd(1);
	
	if(opt.Contains("logz")) {
	  gPad->SetLogz(1);
	} else {
	  gPad->SetLogz(0);
	}
	
	hP2X2->Draw("colz");
	
	y1 = gPad->GetBottomMargin();
	y2 = 1 - gPad->GetTopMargin();
	x1 = gPad->GetLeftMargin();
	x2 = 1 - gPad->GetRightMargin();
	
	TPaveText *textStatInt = new TPaveText(x1+0.02,y2-0.40,x1+0.20,y2-0.05,"NDC");
	PGlobals::SetPaveTextStyle(textStatInt,12); 
	textStatInt->SetTextColor(kGray+3);
	textStatInt->SetTextFont(42);
	
	char text[64];
	sprintf(text,Form("#LTx#GT_{rms} = %5.2f %s",xrms,tspaSUnit.c_str()));
	textStatInt->AddText(text);
	sprintf(text,Form("#LTp_{x}#GT_{rms} = %5.2f %s/c",yrms,teneSUnit.c_str()));
	textStatInt->AddText(text);
	sprintf(text,Form("#varepsilon_{n} = %5.2f %s",emitx,emitSUnit.c_str()));
	textStatInt->AddText(text);
	textStatInt->Draw();
	
	TPaveText **textStat = new TPaveText*[SNbin];
	Int_t pnumber = 0;
	for(Int_t k=0;k<SNbin;k++) {
	  pnumber++;
	  Int_t ic = pnumber%ndiv;
	  
	  // new page
	  if( ic==0 ) {
	    CA4->cd(0);
	    CA4->Clear();
	    CA4->Divide(1,ndiv);
	  }
	  CA4->cd(ic+1);
	  
	  hP2X2sl[k]->GetXaxis()->SetLabelSize(0.08);
	  hP2X2sl[k]->GetXaxis()->SetTitleSize(0.08);
	  hP2X2sl[k]->GetXaxis()->SetTitleOffset(1.0);
	  hP2X2sl[k]->GetXaxis()->CenterTitle();
	  
	  hP2X2sl[k]->GetYaxis()->SetLabelSize(0.08);
	  hP2X2sl[k]->GetYaxis()->SetTitleSize(0.08);
	  hP2X2sl[k]->GetYaxis()->SetTitleOffset(0.8);
	  hP2X2sl[k]->GetYaxis()->CenterTitle();
	  
	  hP2X2sl[k]->GetZaxis()->SetLabelSize(0.08);
	  hP2X2sl[k]->GetZaxis()->SetTitleSize(0.08);
	  hP2X2sl[k]->GetZaxis()->SetTitleOffset(0.8);
	  hP2X2sl[k]->GetZaxis()->CenterTitle();

	  if(opt.Contains("logz")) {
	    gPad->SetLogz(1);
	  } else {
	    gPad->SetLogz(0);
	  }
	  
	  
	  hP2X2sl[k]->Draw("colz");
	  
	  // Double_t y1 = gPad->GetBottomMargin();
	  Double_t y2 = 1 - gPad->GetTopMargin();
	  Double_t x1 = gPad->GetLeftMargin();
	  // Double_t x2 = 1 - gPad->GetRightMargin();
	  textStat[k] = new TPaveText(x1+0.02,y2-0.50,x1+0.20,y2-0.05,"NDC");
	  PGlobals::SetPaveTextStyle(textStat[k],12); 
	  textStat[k]->SetTextColor(kGray+3);
	  textStat[k]->SetTextFont(42);
	  
	  char text[64];
	  sprintf(text,Form("%5.2f %s < #zeta < %5.2f %s",sBinLim[k],spaSUnit.c_str(),sBinLim[k+1],spaSUnit.c_str()));
	  textStat[k]->AddText(text);
	  sprintf(text,Form("#LTx#GT_{rms} = %5.2f %s",sx_rms[k],tspaSUnit.c_str()));
	  textStat[k]->AddText(text);
	  sprintf(text,Form("#LTp_{x}#GT_{rms} = %5.2f %s",spx_rms[k],teneSUnit.c_str()));
	  textStat[k]->AddText(text);
	  sprintf(text,Form("#varepsilon_{n,x} = %5.2f %s",semitx[k],emitSUnit.c_str()));
	  textStat[k]->AddText(text);
	  textStat[k]->Draw();
	  
	  if(ic+1==ndiv) {
	    CA4->cd(0);
	    CA4->Print(fOutName2 + ".ps");
	  }
	}
	
	CA4->Print(fOutName2 + ".ps]");
	
	gSystem->Exec("ps2pdf " + fOutName2 + ".ps " + fOutName2 + ".pdf");
	gSystem->Exec("rm -rf " + fOutName2 + ".ps"); 
      }
      
    }
    
    if(opt.Contains("file")) {
      TString filename = Form("./%s/Plots/Bunch/%s/Bunch-%s_%i.root",sim.Data(),pData->GetSpeciesName(index).c_str(),sim.Data(),time);
      TFile *ofile = new TFile(filename,"RECREATE");

      hX1->SetLineWidth(1);
      hX1->SetLineColor(1);
      hX1->SetFillStyle(0);

      hP1->SetLineWidth(1);
      hP1->SetLineColor(1);
      hP1->SetFillStyle(0);

      hX1->Write("hX1",TObject::kOverwrite);
      hP1->Write("hP1",TObject::kOverwrite);
      hP1X1->Write("hP1X1",TObject::kOverwrite);
      hP2X2->Write("hP2X2",TObject::kOverwrite);
      hX2X1->Write("hX2X1",TObject::kOverwrite);
      
      gErms->Write("gErms",TObject::kOverwrite);
      gXrms->Write("gXrms",TObject::kOverwrite);
      gEmitx->Write("gEmitx",TObject::kOverwrite);

      if(pData->Is3D()) {
	gYrms->Write("gYrms",TObject::kOverwrite);
	gEmity->Write("gEmity",TObject::kOverwrite);

	hP3X3->Write("hP2X2",TObject::kOverwrite);
	hX3X1->Write("hX3X1",TObject::kOverwrite);
	hX3X2->Write("hX3X2",TObject::kOverwrite);
      }
	
      ofile->Close();
    }
    
    // Delete newly created vectors
    delete sBinLim;
    delete zbin;
    delete sEmean;
    delete sErms;

    delete sx_mean;
    delete sx_rms;
    delete spx_mean;
    delete spx_rms;
    delete semitx;
    delete sbetax;

    if(pData->Is3D()) {  
      delete sy_mean;
      delete sy_rms;
      delete spy_mean;
      delete spy_rms;
      delete semity;
      delete sbetay;
    }

    // end time looper
  }
}
