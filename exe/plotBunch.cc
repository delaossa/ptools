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

void InitRange(const TString &, Double_t &, Double_t &, Int_t &, Double_t &, Double_t &);

void FindLimits(TH1F *h,Double_t &xmin, Double_t &xmax, Double_t factor = 0.11) {
  Double_t maxValue = h->GetBinContent(h->GetMaximumBin());
  for(Int_t i=1;i<=h->GetNbinsX();i++) {
    Double_t binValue = h->GetBinContent(i);
    if(binValue>maxValue*factor) {
      xmin = h->GetBinCenter(i);
      break;
    }
  }
  
  for(Int_t i=h->GetNbinsX();i>0;i--) {
    Double_t binValue = h->GetBinContent(i);
    if(binValue>maxValue*factor) {
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
    printf("      <--units> <--logz> <-z(zoom factor)>\n");
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
  Int_t   indexi = 1;

  // zoom
  Float_t  zoom = -1;

  // Option for raw fraction correction
  Float_t rawf = 1;
  
  // Longitudinal range
  Float_t zmin0 =  99999.;
  Float_t zmax0 = -99999.;

  // Slices range
  Float_t zsmin0 =  99999.;
  Float_t zsmax0 = -99999.;

  // Spectrum range
  Float_t Pmin0 =  99999.;
  Float_t Pmax0 = -99999.;

  // Transverse momentum range
  Float_t Pxmax0 = 0.;
    
  // Maximum current
  Float_t Imax = -99999.;

  // Maximum transverse range 
  Float_t Xmax0 = -99999.;

  // Option for longitudinal binning
  // dzf is a fraction of the simulation binning in longitudinal direction
  Float_t dzf = 1.0;

  // Option for transverse binning
  // dxf is a fraction of the simulation binning in transverse direction
  Float_t dxf = 1.0;

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
    } else if(arg.Contains("--loop")){
      opt += "loop"; 
    } else if(arg.Contains("--file")){
      opt += "file"; 
    } else if(arg.Contains("--smooth")){
      opt += "smooth"; 
    } else if(arg.Contains("--autop")){
      opt += "autop"; 
    } else if(arg.Contains("--trans")){
      opt += "trans"; 
    } else if(arg.Contains("--notext")){
      opt += "notext"; 
    } else if(arg.Contains("--noinfo")){
      opt += "noinfo"; 
    } else if(arg.Contains("--noline")){
      opt += "noline"; 
    } else if(arg.Contains("--avge")){
      opt += "avge"; 
    } else if(arg.Contains("--ecorr")){
      opt += "ecorr"; 
    } else if(arg.Contains("--bw")){
      opt += "bw"; 
    } else if(arg.Contains("--nospec")){
      opt += "nospec"; 
    } else if(arg.Contains("--fwhm")){
      opt += "fwhm"; 
    } else if(arg.Contains("--elli")){
      opt += "elli"; 
    } else if(arg.Contains("--eztrack")){
      opt += "eztrack"; 
    } else if(arg.Contains("-index")) {
      char ss[6];
      sscanf(arg,"%6s%i",ss,&indexi);
    } else if(arg.Contains("-imax")) {
      char ss[5];
      sscanf(arg,"%5s%f",ss,&Imax);
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
    } else if(arg.Contains("-dzf")) {
      char ss[4];
      sscanf(arg,"%4s%f",ss,&dzf);
    } else if(arg.Contains("-dxf")) {
      char ss[4];
      sscanf(arg,"%4s%f",ss,&dxf);
    } else if(arg.Contains("-zmin")) {
      char ss[5];
      sscanf(arg,"%5s%f",ss,&zmin0);
    } else if(arg.Contains("-zmax")) {
      char ss[5];
      sscanf(arg,"%5s%f",ss,&zmax0);
    } else if(arg.Contains("-zsmin")) {
      char ss[6];
      sscanf(arg,"%6s%f",ss,&zsmin0);
    } else if(arg.Contains("-zsmax")) {
      char ss[6];
      sscanf(arg,"%6s%f",ss,&zsmax0);
    } else if(arg.Contains("-pmin")) {
      char ss[5];
      sscanf(arg,"%5s%f",ss,&Pmin0);
    } else if(arg.Contains("-pmax")) {
      char ss[5];
      sscanf(arg,"%5s%f",ss,&Pmax0);
    } else if(arg.Contains("-pxmax")) {
      char ss[6];
      sscanf(arg,"%6s%f",ss,&Pxmax0);
    } else if(arg.Contains("-z")) {
      char ss[2];
      sscanf(arg,"%2s%f",ss,&zoom);
    } else if(arg.Contains("-xmax")) {
      char ss[5];
      sscanf(arg,"%5s%f",ss,&Xmax0);
    } else {
      cout << Form("\t Invalid argument (%i): exiting...\n",l) << endl;
      return 0;
    }
  }

  TString simo = sim;
  
  if(opt.Contains("eztrack") && !opt.Contains("loop")) {
    cout << Form(" WARNING: 'eztrack' option only works in loop mode: add --loop flag to the program") << endl;
  }
  
  PGlobals::Initialize();

  if(opt.Contains("grid")) {
    gStyle->SetPadGridX(1);
    gStyle->SetPadGridY(1);
  }

  gStyle->SetJoinLinePS(2);

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
  Double_t chargeUnit,  denUnit,  spaUnit,  propUnit, tspaUnit,  eneUnit,  teneUnit,  curUnit,   emitUnit,   ermsUnit, ermssUnit,  betaUnit, gammaUnit, ezUnit;
  string  chargeSUnit, denSUnit, spaSUnit, propSUnit, tspaSUnit, eneSUnit, teneSUnit, curSUnit,  emitSUnit,  ermsSUnit, ermssSUnit, betaSUnit, gammaSUnit, ezSUnit;

  if(opt.Contains("units") && n0) {

    //  cout << endl << Form("Changing to SI units:") << endl;
      
    // Get the best units for each quantity
    PUnits::BestUnit bdenSUnit(n0,"PartDensity");
    bdenSUnit.GetBestUnits(denUnit,denSUnit);
    //  cout << Form(" n0 = %.2f %s", n0/denUnit, denSUnit.c_str()) << endl;

    PUnits::BestUnit bspaSUnit(wavelength,"Length");
    bspaSUnit.GetBestUnits(spaUnit,spaSUnit);
    //    cout << Form(" L  = %.2f %s", wavelength/spaUnit, spaSUnit.c_str()) << endl;
  }

  // Propagation distance units
  propUnit = PUnits::mm;
  propSUnit = "mm";

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
  // emitUnit = 0.1 * PUnits::um;
  // emitSUnit = "#times 100 nm";
  emitUnit = PUnits::um;
  emitSUnit = "#mum";
  // emitUnit = 10 * PUnits::um;
  // emitSUnit = "10 #mum";
  // emitUnit = 100 * PUnits::nm;
  // emitSUnit = "10^{2} nm";
  // emitUnit = 10 * PUnits::nm;
  // emitSUnit = "10 nm";

  // beta units
  betaUnit = PUnits::mm;
  betaSUnit = "mm";
  
  // gamma units
  gammaUnit = 1/PUnits::m;
  gammaSUnit = "m^{-1}";
  
  // Relative energy spread units
  // ermsUnit = 0.1 * PUnits::perCent; 
  // ermsSUnit = "0.1%";
  ermsUnit = PUnits::perCent; 
  ermsSUnit = "%";
  // ermsUnit = PUnits::perThousand; 
  // ermsSUnit = "0.1 %";

  // Slice elative energy spread units  
  // ermssUnit = 0.1 * PUnits::perCent; 
  // ermssSUnit = "0.1%";
  ermssUnit = PUnits::perCent; 
  ermssSUnit = "%";

  
  // Ez unit
  ezUnit = PUnits::GV/PUnits::meter;
  ezSUnit = "GV/m";
  
  if (sim.Contains("pitz")) {
    curUnit = PUnits::ampere;
    curSUnit = "A";
    
    eneUnit = PUnits::MeV;
    eneSUnit = "MeV";
    
    ermsUnit = PUnits::perMillion; 
    ermsSUnit = "0.01 %";
  }
  
  
  // Time looper
  for(Int_t i=iStart; i<iEnd+1; i+=iStep) {

    time = i;
    pData->LoadFileNames(time);    
    //if(time==iStart) pData->PrintData();
    
    if(!pData->IsInit()) continue;

    Int_t Nspecies = pData->NRawSpecies();
    Int_t index = indexi; 
    if(indexi>Nspecies-1) {
      index = Nspecies-1;
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
    
    // Binning, ranges, labels, etc.
    Int_t x1Nbin = 200;
    Int_t p1Nbin = 200;
    Int_t x2Nbin = 200;
    Int_t p2Nbin = 200;
    Int_t x3Nbin = 200;
    Int_t p3Nbin = 200;
    
    // Spatial range:
    Double_t x1Min = X1MIN;
    Double_t x1Max = X1MAX;
    Double_t x2Min = X2MIN;
    Double_t x2Max = X2MAX;
    Double_t x3Min = X3MIN;
    Double_t x3Max = X3MAX;

    x1Nbin = (X1MAX-X1MIN)/dx1;
    x2Nbin = (X2MAX-X2MIN)/dx2;
    x3Nbin = (X3MAX-X3MIN)/dx3;
    
    // Momentum range:
    Double_t p1Min =  0.;
    Double_t p1Max =  4000.;
    Double_t p2Min = -15.0;
    Double_t p2Max =  15.0;
    Double_t p3Min = -15.0;
    Double_t p3Max =  15.0;

    // z-slices
    Int_t SNbin = 100;
    Double_t x1BinMin = -4.5;
    Double_t x1BinMax = -4.0;    

    // dummy shift
    // Double_t dshiftz = pData->Shift("centercomov");
    Double_t dshiftz = pData->Shift(opt.Data());

    // Specific initializations:
    if(opt.Contains("autop")) {
      InitRange(sim,x1Min,x1Max,SNbin,x1BinMin,x1BinMax);
      
      x1Min += dshiftz;
      x1Max += dshiftz;
      x1BinMin += dshiftz;
      x1BinMax += dshiftz;
    }

    // Command line input
    Double_t Pmin = Pmin0;
    Double_t Pmax = Pmax0;
    if (Pmax0 > Pmin0) {
      
      if(opt.Contains("units")) {
	Pmin *= eneUnit / PConst::ElectronMassE;
	Pmax *= eneUnit / PConst::ElectronMassE;
      }
      
    }
    
    Double_t zmin = zmin0;
    Double_t zmax = zmax0;
    Double_t zsmin = zsmin0;
    Double_t zsmax = zsmax0;
    
    if(zmax0>zmin0) {

      if(opt.Contains("units")) {
	zmin *= spaUnit * kp;
	zmax *= spaUnit * kp;
      }
      
      x1Min = zmin + dshiftz;
      x1Max = zmax + dshiftz;
      x1BinMin = x1Min + (x1Max-x1Min)/4.0;
      x1BinMax = x1Min + 3.0*(x1Max-x1Min)/4.0;
    }

    if(zsmax0>zsmin0) {
      if(opt.Contains("units")) {
	zsmin *= spaUnit * kp;
	zsmax *= spaUnit * kp;
      }
      
      x1BinMin = zsmin + dshiftz;
      x1BinMax = zsmax + dshiftz;
    }

    // Transverse range
    Double_t Xmax = Xmax0;
    if (Xmax > 0.) {
      
      if(opt.Contains("units")) {
	Xmax *= tspaUnit * kp;
      }
      
    }

    Double_t Pxmax = Pxmax0;
    if (Pxmax > 0.) {
      
      if(opt.Contains("units")) {
	Pxmax *= teneUnit / PConst::ElectronMassE;
      }
      
    }
    
    
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
    Int_t NX1 = pData->GetX1N()*(x1Max-x1Min)/(pData->GetX1Max()-pData->GetX1Min());

    TH1F *hScanX1 = (TH1F*) gROOT->FindObject("hScanX1");
    if(hScanX1) delete hScanX1;
    hScanX1 = new TH1F("hScanX1","",NX1,x1Min,x1Max);
    
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
      
    // Set limits using extreme values

    if(Pmin<99999.)
      p1Min = Pmin;
    else
      p1Min = MinP1 - rfactor*(MaxP1-MinP1);

    if(Pmax>-99999.)
      p1Max = Pmax;
    else
      p1Max = MaxP1 + rfactor*(MaxP1-MinP1);

    // p2Min = MinP2 - rfactor*(MaxP2-MinP2);
    // p2Max = MaxP2 + rfactor*(MaxP2-MinP2);
    // p3Min = MinP3 - rfactor*(MaxP3-MinP3);
    // p3Max = MaxP3 + rfactor*(MaxP3-MinP3);

    // x1Min = MinX1 - rfactor*(MaxX1-MinX1);
    // x1Max = MaxX1 + rfactor*(MaxX1-MinX1);
    // x2Min = MinX2 - rfactor*(MaxX2-MinX2);
    // x2Max = MaxX2 + rfactor*(MaxX2-MinX2);    
    // if(Nvar==7) {
    //   x3Min = MinX3 - rfactor*(MaxX3-MinX3);
    //   x3Max = MaxX3 + rfactor*(MaxX3-MinX3);
    // }

    // Set limits using Scan histograms
    // --------------------------------
    
    Double_t peakFactor = 0.2;
    Double_t rfactor2 = 1.5;
      
    Double_t x1min = -999;
    Double_t x1max = -999;

    FindLimits(hScanX1,x1min,x1max,peakFactor);
    
    x1BinMin = x1min;
    x1BinMax = x1max;
      
    peakFactor = 0.1;
    FindLimits(hScanX1,x1min,x1max,peakFactor);
    
    x1Min = x1min - rfactor*(x1max-x1min);
    x1Max = x1max + 2*rfactor*(x1max-x1min);

    // Transverse plane
    // Zoom trasverse range
    if((zoom>0) && (Xmax<0)) {
      Double_t x2Range = (x2Max - x2Min)/zoom;
      Double_t x2Mid   = (x2Max + x2Min)/2.0;
      x2Min = x2Mid - x2Range/2.0;
      x2Max = x2Mid + x2Range/2.0;
      if(pData->Is3D()) {
	Double_t x3Range = (x3Max - x3Min)/zoom;
	Double_t x3Mid   = (x3Max + x3Min)/2.0;
	x3Min = x3Mid - x3Range/2.0;
	x3Max = x3Mid + x3Range/2.0;
      }
    } else if(Xmax>0) {
      Double_t x2Mid   = (x2Max + x2Min)/2.0;
      x2Min = x2Mid - Xmax;
      x2Max = x2Mid + Xmax;
      if(pData->Is3D()) {
	Double_t x3Mid   = (x3Max + x3Min)/2.0;
	x3Min = x3Mid - Xmax;
	x3Max = x3Mid + Xmax;
      }
    } else {
      Double_t x2min = -999;
      Double_t x2max = -999;
      
      FindLimits(hScanX2,x2min,x2max,peakFactor);
      if( x2max-x2min < 0.1) {
	x2min -= 0.1;
	x2max += 0.1;
      }
      
      x2Min = x2min - rfactor2*(x2max-x2min);
      x2Max = x2max + rfactor2*(x2max-x2min);
      
      Double_t x3min = -999;
      Double_t x3max = -999;
      if(Nvar==7) {
	FindLimits(hScanX3,x3min,x3max,peakFactor);
	x3Min = x3min - rfactor2*(x3max-x3min);
	x3Max = x3max + rfactor2*(x3max-x3min);
      }
    }
    
    if(Pxmax > 0.) {
      p2Min = -Pxmax;
      p2Max = Pxmax;

      p3Min = -Pxmax;
      p3Max = Pxmax;
    } else {
      Double_t p2min = -999;
      Double_t p2max = -999;
      FindLimits(hScanP2,p2min,p2max,peakFactor);
      p2Min = p2min - rfactor2*(p2max-p2min);
      p2Max = p2max + rfactor2*(p2max-p2min);
      if(p2Max<-p2Min) p2Max = -p2Min; 
      else p2Min = -p2Max;
      
      Double_t p3min = -999;
      Double_t p3max = -999;
      FindLimits(hScanP3,p3min,p3max,peakFactor);
      p3Min = p3min - rfactor2*(p3max-p3min);
      p3Max = p3max + rfactor2*(p3max-p3min);
      if(p3Max<-p3Min) p3Max = -p3Min; 
      else p3Min = -p3Max;
    }

    
    // Prevent range limits to go out of the box
    if(x1Min<X1MIN) x1Min = X1MIN;
    if(x1Max>X1MAX) x1Max = X1MAX;
    if(x2Min<X2MIN) x2Min = X2MIN;
    if(x2Max>X2MAX) x2Max = X2MAX;
    if(pData->Is3D()) {
      if(x3Min<X3MIN) x3Min = X3MIN;
      if(x3Max>X3MAX) x3Max = X3MAX;

      if(fabs(x3Max)<fabs(x2Max)) {
	x3Min = x2Min;
	x3Max = x2Max;	
      }
    }

    // Overrides auto z ranging if specified in command line
    if(zmax>zmin) {
      x1Min = zmin + dshiftz;
      x1Max = zmax + dshiftz;
    }

    if(zsmax>zsmin) {
      x1BinMin = zsmin + dshiftz;
      x1BinMax = zsmax + dshiftz;
    }
    
    
    // Adjust the binning to match the simulation grid
    // -----
    // cout << Form(" x1 range (N = %i) :  x1Min = %f  x1Max = %f  dx1 = %f", x1Nbin, x1Min, x1Max, (x1Max - x1Min)/x1Nbin) << endl;
    // cout << Form(" x1-slices (N = %i):  x1Min = %f  x1Max = %f  Dx1 = %f", SNbin, x1BinMin, x1BinMax, (x1BinMax-x1BinMin)/SNbin) << endl;

    Double_t ddx1 = dzf * dx1;
    x1Min = floor((x1Min-X1MIN)/ddx1) * ddx1 + X1MIN;  
    x1Max = ceil((x1Max-X1MIN)/ddx1) * ddx1 + X1MIN;
    x1Nbin = ceil ((x1Max - x1Min)/(ddx1));
    p1Nbin = x1Nbin;
    
    // slices
    x1BinMin = floor((x1BinMin-X1MIN)/ddx1) * ddx1 + X1MIN;  
    x1BinMax = ceil((x1BinMax-X1MIN)/ddx1) * ddx1 + X1MIN;
    
    SNbin  = ceil ((x1BinMax - x1BinMin)/(ddx1));

    x1Min -= shiftz;
    x1Max -= shiftz;
    x1BinMin -= shiftz;
    x1BinMax -= shiftz;
    
    // cout << Form("\n x1 range (N = %i) :  x1Min = %f  x1Max = %f  dx1 = %f", x1Nbin, x1Min, x1Max, (x1Max - x1Min)/x1Nbin) << endl;
    // cout << Form(" x1-slices (N = %i):  x1Min = %f  x1Max = %f  Dx1 = %f", SNbin, x1BinMin, x1BinMax, (x1BinMax-x1BinMin)/SNbin) << endl;

    Double_t ddx2 = dxf * dx2;      
    // Make the transverse edges to match cell boundaries
    x2Min = floor((x2Min-X2MIN)/ddx2) * ddx2 + X2MIN;  
    x2Max = ceil((x2Max-X2MIN)/ddx2) * ddx2 + X2MIN;  
    x2Nbin = ceil ((x2Max - x2Min)/(ddx2));
    // if(x2Nbin<100) x2Nbin = 100;
    p2Nbin = x2Nbin;
    
    if(pData->Is3D()) {
      Double_t ddx3 = dxf * dx3;      
      x3Min = floor((x3Min-X3MIN)/ddx3) * ddx3 + X3MIN;  
      x3Max = ceil((x3Max-X3MIN)/ddx3) * ddx3 + X3MIN;  
      x3Nbin = ceil ((x3Max - x3Min)/(ddx3));
    //   if(x3Nbin<100) x3Nbin = 100;
      p3Nbin = x3Nbin;
    } else 
      p3Nbin = x2Nbin;
    

    // other
    if(p1Min < 0.0) p1Min = 0.01;
    // if(x1Nbin < 100) x1Nbin = 100;
    // if(x2Nbin < 100) x2Nbin = 100;
    // if(x3Nbin < 100) x3Nbin = 100;

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

    hX1->SetFillStyle(1001);
    hX1->SetFillColor(PGlobals::elecFill);
    hX1->SetLineWidth(3);
    
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

    sprintf(hName,"hP1range");
    TH1F *hP1range = (TH1F*) gROOT->FindObject(hName);
    if(hP1range) delete hP1range;
    hP1range = new TH1F(hName,"",p1Nbin,p1Min,p1Max);
    hP1range->GetYaxis()->SetTitle("p_{z}/mc");
    if(opt.Contains("comov"))
      hP1range->GetXaxis()->SetTitle("k_{p}#zeta");
    else
      hP1range->GetXaxis()->SetTitle("k_{p}z");

    sprintf(hName,"hP1X1range");
    TH2F *hP1X1range = (TH2F*) gROOT->FindObject(hName);
    if(hP1X1range) delete hP1X1range;
    hP1X1range = new TH2F(hName,"",x1Nbin,x1Min,x1Max,p1Nbin,p1Min,p1Max);
    if(opt.Contains("comov"))
      hP1X1range->GetXaxis()->SetTitle("k_{p}#zeta");
    else
      hP1X1range->GetXaxis()->SetTitle("k_{p}z");
    hP1X1range->GetYaxis()->SetTitle("p_{z}/mc");
    hP1X1range->GetZaxis()->SetTitle("Charge [n_{0}dV]");
    hP1X1range->GetZaxis()->CenterTitle();

    
    sprintf(hName,"hP2X2");
    TH2F *hP2X2 =  (TH2F*) gROOT->FindObject(hName);
    if(hP2X2) delete hP2X2;
    hP2X2 = new TH2F(hName,"",x2Nbin,x2Min,x2Max,p2Nbin,p2Min,p2Max);
    hP2X2->GetXaxis()->SetTitle("k_{p}x");
    hP2X2->GetYaxis()->SetTitle("p_{x}/mc");
    hP2X2->GetZaxis()->SetTitle("Charge [n_{0}dV]");

    sprintf(hName,"hP2X2range");
    TH2F *hP2X2range = (TH2F*) gROOT->FindObject(hName);
    if(hP2X2range) delete hP2X2range;
    hP2X2range = new TH2F(hName,"",x2Nbin,x2Min,x2Max,p2Nbin,p2Min,p2Max);
    hP2X2range->GetXaxis()->SetTitle("k_{p}x");
    hP2X2range->GetYaxis()->SetTitle("p_{x}/mc");
    hP2X2range->GetZaxis()->SetTitle("Charge [n_{0}dV]");
    hP2X2range->GetZaxis()->CenterTitle();

    TH2F *hP3X3(0), *hP3X3range(0);
    if(pData->Is3D()) {
      sprintf(hName,"hP3X3");
      hP3X3 =  (TH2F*) gROOT->FindObject(hName);
      if(hP3X3) delete hP3X3;
      hP3X3 = new TH2F(hName,"",x3Nbin,x3Min,x3Max,p3Nbin,p3Min,p3Max);
      hP3X3->GetXaxis()->SetTitle("k_{p}y");
      hP3X3->GetYaxis()->SetTitle("p_{y}/mc");
      hP3X3->GetZaxis()->SetTitle("Charge [n_{0}dV]");
      hP3X3->GetZaxis()->CenterTitle();
      
      sprintf(hName,"hP3X3range");
      hP3X3range = (TH2F*) gROOT->FindObject(hName);
      if(hP3X3range) delete hP3X3range;
      hP3X3range = new TH2F(hName,"",x3Nbin,x3Min,x3Max,p3Nbin,p3Min,p3Max);
      hP3X3range->GetXaxis()->SetTitle("k_{p}y");
      hP3X3range->GetYaxis()->SetTitle("p_{y}/mc");
      hP3X3range->GetZaxis()->SetTitle("Charge [n_{0}dV]");
      hP3X3range->GetZaxis()->CenterTitle();
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
      hP2X2->Fill(var[5][i],var[1][i],TMath::Abs(var[3][i]));
      hX2X1->Fill(var[4][i],var[5][i],TMath::Abs(var[3][i]));

      if(Nvar==7){
	hP3X3->Fill(var[6][i],var[2][i],TMath::Abs(var[3][i]));
	hX3X1->Fill(var[4][i],var[6][i],TMath::Abs(var[3][i]));
	hX3X2->Fill(var[5][i],var[6][i],TMath::Abs(var[3][i]));      
      }
      
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
      // cout << Form("  iBin = %i   p1 = %7.3f",iBin,var[0][i]) << endl;

      // Fill histograms in slices range
      hP1sl[iBin]->Fill(var[0][i],TMath::Abs(var[3][i]));
      hP1range->Fill(var[0][i],TMath::Abs(var[3][i]));
      hP1X1range->Fill(var[4][i],var[0][i],TMath::Abs(var[3][i]));
      
      hP2X2sl[iBin]->Fill(var[5][i],var[1][i],TMath::Abs(var[3][i]));

      hP2X2range->Fill(var[5][i],var[1][i],TMath::Abs(var[3][i]));

      if(Nvar==7) {
	hP3X3sl[iBin]->Fill(var[6][i],var[2][i],TMath::Abs(var[3][i]));
	hP3X3range->Fill(var[6][i],var[2][i],TMath::Abs(var[3][i]));
      }
      
    }
    // cout << " done! " << endl;

    // Integrated long. emittance:
    cout << Form("\n 3. Calculating integrated quantities: ") << endl ;

    // Longitudinal phasespace (Total):
    // --

    Double_t stats[7];  // { sumw, sumw2, sumwx, sumwx2, sumwy, sumwy2, sumwxy }
    Double_t xmean,xrms2,xrms,ymean,yrms2,yrms,xyrms,emit2,emit;

    // P1X1 Phasespace
    hP1X1->GetStats(stats);
    
    xmean  = stats[2]/stats[0];
    xrms2  = stats[3]/stats[0] - xmean * xmean;
    xrms   = (xrms2>0.0) ? TMath::Sqrt(xrms2) : 0.0 ;
    ymean  = stats[4]/stats[0];
    yrms2  = stats[5]/stats[0] - ymean * ymean;
    yrms   = (yrms2>0.0) ? TMath::Sqrt(yrms2) : 0.0 ;
    xyrms = stats[6]/stats[0] - stats[2]*stats[4]/(stats[0]*stats[0]);
    emit2 = xrms2*yrms2 - xyrms*xyrms;
    emit = (emit2>0.0) ? TMath::Sqrt(emit2) : 0.0 ;
    
    Double_t emitz = emit;
    Double_t zmean = xmean;
    Double_t zrms = xrms;
    Double_t pzmean = ymean;
    Double_t pzrms = yrms;
    Double_t zpzcorr = xyrms;
    Double_t zpzcorrrel = 100 * zpzcorr/(pzmean * xrms2);
    
    // Total relative energy spread within FWHM:
    sprintf(hName,"hP1cut");
    TH1F *hP1cut = (TH1F*) hP1->Clone(hName);
    hP1cut->Reset();

    Double_t maxValue = hP1->GetBinContent(hP1->GetMaximumBin());
    Int_t   lBin = -1;
    Double_t pzmin = -999;
    for(Int_t i=1;i<=hP1->GetNbinsX();i++) {
      Double_t binValue = hP1->GetBinContent(i);
      if(binValue>maxValue/4.0) {
	lBin = i;
	pzmin = hP1->GetBinCenter(i);
	break;
      }
    }

    Int_t rBin = -1;
    Double_t pzmax = -999;
    for(Int_t i=hP1->GetNbinsX();i>0;i--) {
      Double_t binValue = hP1->GetBinContent(i);
      if(binValue>maxValue/4.0) {
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

    cout << Form(" Rel. energy spread (fwhm) = %f",pzrmsFWHM/pzmeanFWHM) << endl;
    
    // Basically the same, but in longitudinal range
    hP1X1range->GetStats(stats);

    xmean  = stats[2]/stats[0];
    xrms2  = stats[3]/stats[0] - xmean * xmean;
    xrms   = (xrms2>0.0) ? TMath::Sqrt(xrms2) : 0.0 ;
    ymean  = stats[4]/stats[0];
    yrms2  = stats[5]/stats[0] - ymean * ymean;
    yrms   = (yrms2>0.0) ? TMath::Sqrt(yrms2) : 0.0 ;
    xyrms = stats[6]/stats[0] - stats[2]*stats[4]/(stats[0]*stats[0]);
    emit2 = xrms2*yrms2 - xyrms*xyrms;
    emit = (emit2>0.0) ? TMath::Sqrt(emit2) : 0.0 ;

    // pzmean = ymean;
    // pzrms = yrms;
    // zpzcorr = xyrms;
    // zpzcorrrel = xyrms/pzmean;
    
    cout << Form("  zMean = %7.3f   pzMean = %7.3f",zmean,pzmean) << endl;
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
    xyrms = stats[6]/stats[0] - stats[2]*stats[4]/(stats[0]*stats[0]);
    emit2 = xrms2*yrms2 - xyrms*xyrms;
    emit = (emit2>0.0) ? TMath::Sqrt(emit2) : 0.0 ; 

    Double_t emitx = emit;
    Double_t x_mean = stats[2]/stats[0];
    Double_t x_rms = (xrms2>0.0) ? TMath::Sqrt(xrms2) : 0.0 ;
    Double_t px_mean = stats[4]/stats[0];
    Double_t px_rms = (yrms2>0.0) ? TMath::Sqrt(yrms2) : 0.0 ;
    Double_t xpx_rms = xyrms;

    Double_t beta = xrms2 / emit;
    Double_t gamma = yrms2 / emit;
    Double_t alpha = -xyrms / emit;
  
    //    Double_t grel = hP1->GetMean() * eneUnit/PConst::ElectronMassE;
    Double_t grel = hP1->GetMean();
    Double_t betax = beta * grel;
    Double_t gammax = gamma / grel;
    Double_t alphax = alpha;

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

    // in slices range
    hP2X2range->GetStats(stats);

    xrms2  = stats[3]/stats[0] - stats[2]*stats[2]/(stats[0]*stats[0]);
    yrms2  = stats[5]/stats[0] - stats[4]*stats[4]/(stats[0]*stats[0]);
    xyrms = stats[6]/stats[0] - stats[2]*stats[4]/(stats[0]*stats[0]);
    emit2 = xrms2*yrms2 - xyrms*xyrms;
    emit = (emit2>0.0) ? TMath::Sqrt(emit2) : 0.0 ; 

    Double_t emitx_r = emit;
    Double_t x_mean_r = stats[2]/stats[0];
    Double_t x_rms_r = (xrms2>0.0) ? TMath::Sqrt(xrms2) : 0.0 ;
    Double_t px_mean_r = stats[4]/stats[0];
    Double_t px_rms_r = (yrms2>0.0) ? TMath::Sqrt(yrms2) : 0.0 ;

    beta = xrms2 / emit;
    gamma = yrms2 / emit;
    alpha = -xyrms / emit;

    Double_t grel_r = hP1range->GetMean();
    Double_t betax_r = beta * grel_r;
    Double_t gammax_r = gamma / grel_r;
    Double_t alphax_r = alpha;

    // P3X3 SPACE -----------------------
    Double_t emity = 0, y_mean = 0, y_rms = 0, py_mean = 0, py_rms = 0, ypy_rms = 0, betay = 0;
    Double_t emity_r = 0, y_mean_r = 0, y_rms_r = 0, py_mean_r = 0, py_rms_r = 0, betay_r = 0;

    TEllipse *ellipP3X3 = NULL;
    if(pData->Is3D()) {
      hP3X3->GetStats(stats);

      xrms2  = stats[3]/stats[0] - stats[2]*stats[2]/(stats[0]*stats[0]);
      yrms2  = stats[5]/stats[0] - stats[4]*stats[4]/(stats[0]*stats[0]);
      xyrms = stats[6]/stats[0] - stats[2]*stats[4]/(stats[0]*stats[0]);
      emit2 = xrms2*yrms2 - xyrms*xyrms;
      emit = (emit2>0.0) ? TMath::Sqrt(emit2) : 0.0 ; 
      
      emity = emit;
      y_mean = stats[2]/stats[0];
      y_rms = (xrms2>0.0) ? TMath::Sqrt(xrms2) : 0.0 ;
      py_mean = stats[4]/stats[0];
      py_rms = (yrms2>0.0) ? TMath::Sqrt(yrms2) : 0.0 ;
      ypy_rms = xyrms;
      
      beta = xrms2 / emit;
      gamma = yrms2 / emit;
      alpha = -xyrms / emit;
      
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

      ellipP3X3 = new TEllipse(y_mean,py_mean,TMath::Sqrt(a2),TMath::Sqrt(b2),0.,360.,angle * 180. / PConst::pi );
      ellipP3X3->SetFillStyle(0);
      ellipP3X3->SetLineStyle(2);
      ellipP3X3->SetLineColor(2);
      ellipP3X3->SetLineWidth(1);

      // in slices range
      hP3X3range->GetStats(stats);

      xrms2  = stats[3]/stats[0] - stats[2]*stats[2]/(stats[0]*stats[0]);
      yrms2  = stats[5]/stats[0] - stats[4]*stats[4]/(stats[0]*stats[0]);
      xyrms = stats[6]/stats[0] - stats[2]*stats[4]/(stats[0]*stats[0]);
      emit2 = xrms2*yrms2 - xyrms*xyrms;
      emit = (emit2>0.0) ? TMath::Sqrt(emit2) : 0.0 ; 

      emity_r = emit;
      y_mean_r = stats[2]/stats[0];
      y_rms_r = (xrms2>0.0) ? TMath::Sqrt(xrms2) : 0.0 ;
      py_mean_r = stats[4]/stats[0];
      py_rms_r = (yrms2>0.0) ? TMath::Sqrt(yrms2) : 0.0 ;

      Double_t beta = xrms2 / emit;
      betay_r = beta * grel_r;

    }

    cout << Form("\n 4. Calculating sliced quantities.. ") << endl ;

    TGraph *gEmitx = NULL;
    TGraph *gEmity = NULL;
    TGraph *gXrms = NULL;
    TGraph *gYrms = NULL;
    TGraph *gErms = NULL;

    Double_t *zbin = new Double_t[SNbin];
    Double_t *sEmean = new Double_t[SNbin];
    Double_t *sErms = new Double_t[SNbin];

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

      Double_t sxmean = 0, sxrms2 = 0, sxrms = 0, symean = 0, syrms2 = 0, syrms = 0, sxyrms = 0, semit = 0, semit2 = 0;
      
      // P2X2 slices
      sellipP2X2[k] = NULL;            
      hP2X2sl[k]->GetStats(stats);
      
      sxmean  = stats[2]/stats[0];
      sxrms2  = stats[3]/stats[0] - sxmean * sxmean;
      sxrms   = (sxrms2>0.0) ? TMath::Sqrt(sxrms2) : 0.0 ;
      symean  = stats[4]/stats[0];
      syrms2  = stats[5]/stats[0] - symean * symean;
      syrms   = (syrms2>0.0) ? TMath::Sqrt(syrms2) : 0.0 ;
      sxyrms = stats[6]/stats[0] - stats[2]*stats[4]/(stats[0]*stats[0]);
      semit2 = sxrms2*syrms2 - sxyrms*sxyrms;
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
      alpha = - sxyrms / semit;

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
	sxyrms = stats[6]/stats[0] - stats[2]*stats[4]/(stats[0]*stats[0]);
	semit2 = sxrms2*syrms2 - sxyrms*sxyrms;
	semit = (semit2>0.0) ? TMath::Sqrt(semit2) : 0.0 ;
	
	sy_mean[k] = sxmean;
	spy_mean[k] = symean;
	sy_rms[k] = sxrms;
	spy_rms[k] = syrms;
	semity[k] = semit;
	
	// Ellipse calculation
	beta = sxrms2 / semit;
	gamma = syrms2 / semit ;
	alpha = - sxyrms / semit;
	
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

    Double_t Charge = hX1->Integral();
    if(!pData->Is3D() && pData->GetBeamRmsX()) {
      Float_t factor =  (pData->GetBeamRmsX()/skindepth) * TMath::Sqrt(2.0*TMath::Pi());
      Charge *= factor; 
    }

    Float_t dz = (x1Max-x1Min)/x1Nbin;
    hX1->Scale( (Q0/(dz*skindepth)) * PConst::c_light / PConst::I0);
    
    if(opt.Contains("units") && n0) {
      Charge *= Q0;
      
      if(opt.Contains("best")) {
	PUnits::BestUnit bchargeSUnit(Charge,"Charge");
	bchargeSUnit.GetBestUnits(chargeUnit,chargeSUnit);
      }
      
      Charge /= chargeUnit;

      Time *= skindepth / propUnit;

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
      hP1X1->GetZaxis()->SetTitle(Form("Charge [%s]",chargeSUnit.c_str()));
      hP1X1->GetYaxis()->SetTitle(Form("p_{z} [%s/c]",eneSUnit.c_str()));
      if(opt.Contains("comov"))
	hP1X1->GetXaxis()->SetTitle(Form("#zeta [%s]",spaSUnit.c_str()));
      else
	hP1X1->GetXaxis()->SetTitle(Form("z [%s]",spaSUnit.c_str()));

      // Converting electron density
      hP1->Scale(Q0 / chargeUnit);
      hP1->SetBins(p1Nbin,p1Min,p1Max);
      hP1->GetYaxis()->SetTitle(Form("p_{z} [%s/c]",eneSUnit.c_str()));

      hX1->SetBins(x1Nbin,x1Min,x1Max);
      // Double_t dt = ((x1Max - x1Min)/x1Nbin) * spaUnit / PConst::c_light;
      // hX1->Scale(Q0 / dt);

      hX1->Scale(PConst::I0);
      if(opt.Contains("best")) {
	PUnits::BestUnit bcurSUnit(hX1->GetMaximum(),"Current");
	bcurSUnit.GetBestUnits(curUnit,curSUnit);
      }
      hX1->Scale(1/curUnit);
      
      hX1->GetYaxis()->SetTitle(Form("I[%s]",curSUnit.c_str()));
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
      zpzcorr *= (PConst::ElectronMassE/eneUnit) * (skindepth / spaUnit);
      zpzcorrrel /= (skindepth / spaUnit);
      
      // Transverse phase-space
      x_mean *= skindepth / tspaUnit;
      x_rms  *= skindepth / tspaUnit;
      px_mean *= PConst::ElectronMassE / teneUnit;
      px_rms  *= PConst::ElectronMassE / teneUnit;
      xpx_rms *= (PConst::ElectronMassE / teneUnit) * (skindepth / tspaUnit);

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
      
      hP2X2range->GetXaxis()->SetTitle(Form("x [%s]",spaSUnit.c_str()));
      hP2X2range->GetYaxis()->SetTitle(Form("p_{x} [%s/c]",teneSUnit.c_str()));
      hP2X2range->GetZaxis()->SetTitle(Form("Charge [%s]",chargeSUnit.c_str()));
      hP2X2range->GetZaxis()->CenterTitle();
      
      hX2X1->SetBins(x1Nbin,x1Min,x1Max,x2Nbin,x2Min,x2Max);
      hX2X1->Scale(Q0 / chargeUnit);
      hX2X1->GetXaxis()->SetTitle(Form("#zeta [%s]",spaSUnit.c_str()));
      hX2X1->GetYaxis()->SetTitle(Form("x [%s]",tspaSUnit.c_str()));
      hX2X1->GetZaxis()->SetTitle(Form("Charge [%s]",chargeSUnit.c_str()));
      hX2X1->GetZaxis()->CenterTitle();

      if(pData->Is3D()) {
	y_mean *= skindepth / tspaUnit;
	y_rms  *= skindepth / tspaUnit;
	py_mean *= PConst::ElectronMassE / teneUnit;
	py_rms  *= PConst::ElectronMassE / teneUnit;      
	ypy_rms  *= (PConst::ElectronMassE / teneUnit) * (skindepth / tspaUnit);      
	
	x3Min *= skindepth/tspaUnit;
	x3Max *= skindepth/tspaUnit;
	p3Min *= PConst::ElectronMassE / teneUnit;
	p3Max *= PConst::ElectronMassE / teneUnit;
	
	hP3X3->SetBins(x3Nbin,x3Min,x3Max,p3Nbin,p3Min,p3Max);
	hP3X3->Scale(Q0 / chargeUnit);
	
	hP3X3->GetXaxis()->SetTitle(Form("y [%s]",spaSUnit.c_str()));
	hP3X3->GetYaxis()->SetTitle(Form("p_{y} [%s/c]",teneSUnit.c_str()));
	hP3X3->GetZaxis()->SetTitle(Form("Charge [%s]",chargeSUnit.c_str()));
	hP3X3->GetZaxis()->CenterTitle();

	hP3X3range->GetXaxis()->SetTitle(Form("y [%s]",spaSUnit.c_str()));
	hP3X3range->GetYaxis()->SetTitle(Form("p_{y} [%s/c]",teneSUnit.c_str()));
	hP3X3range->GetZaxis()->SetTitle(Form("Charge [%s]",chargeSUnit.c_str()));
	hP3X3range->GetZaxis()->CenterTitle();

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
      emitx_r *= skindepth/emitUnit;

      betax *= skindepth;
      betax_r *= skindepth;
      if(opt.Contains("best")) {
	PUnits::BestUnit bbetaSUnit(betax,"Length");
	bbetaSUnit.GetBestUnits(betaUnit,betaSUnit);
	cout << bbetaSUnit << endl;
      }
      betax /= betaUnit;
      betax_r /= betaUnit;

      gammax *= kp;
      gammax_r *= kp;

      gammax /= gammaUnit;
      gammax_r /= gammaUnit;      
      alphax *= 1.;
      alphax_r *= 1.;

      cout << Form(" Emit = %f.3 %s  Beta = %f.3 %s  gamma = %f.3 %s  alpha = %f.3",emitx,emitSUnit.c_str(),betax,betaSUnit.c_str(),gammax,gammaSUnit.c_str(),alphax) << endl;
      
      cout << Form(" Emit = %f.3 %s  Beta = %f.3 %s  gamma = %f.3 %s  alpha = %f.3",emitx_r,emitSUnit.c_str(),betax_r,betaSUnit.c_str(),gammax_r,gammaSUnit.c_str(),alphax_r) << endl;

      if(pData->Is3D()) {
	emity *= skindepth / emitUnit;
	betay *= skindepth / betaUnit;
	emity_r *= skindepth / emitUnit;
	betay_r *= skindepth / betaUnit;
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
       	sErms[k]  *= PUnits::perCent/ermssUnit;
      
    }
    // End of the users units module    

    // Sliced average quantities 
    Double_t xmeanavg = 0;
    Double_t xrmsavg = 0;
    Double_t emitxavg = 0;
    Double_t sErmsavg = 0;
    
    // Extract from selected slices:
    Double_t norm = 0;
    for(Int_t i=0;i<SNbin;i++) {
      emitxavg += hP2X2sl[i]->GetEntries() * semitx[i];
      xmeanavg += hP2X2sl[i]->GetEntries() * sx_mean[i];
      xrmsavg  += hP2X2sl[i]->GetEntries() * sx_rms[i];
      sErmsavg += hP2X2sl[i]->GetEntries() * sErms[i];
      norm += hP2X2sl[i]->GetEntries();
    }
    emitxavg /= norm;
    xmeanavg /= norm;
    xrmsavg /= norm;
    sErmsavg /= norm;
    
    Double_t ymeanavg = 0;
    Double_t yrmsavg = 0;
    Double_t emityavg = 0;
    if(pData->Is3D()) {
      norm = 0;
      for(Int_t i=0;i<SNbin;i++) {
	emityavg += hP3X3sl[i]->GetEntries() * semity[i];
	ymeanavg += hP3X3sl[i]->GetEntries() * sy_mean[i];
	yrmsavg  += hP3X3sl[i]->GetEntries() * sy_rms[i];
	norm += hP3X3sl[i]->GetEntries();
      }
      emityavg /= norm;
      ymeanavg /= norm;
      yrmsavg /= norm;
    }
    
    if(opt.Contains("units")) {
      cout << "\n  Summary _______________________________________________________ " << endl;
      cout << Form("  Integrated charge (RAW) of specie %3i = %8f %s",index,Charge,chargeSUnit.c_str()) << endl;
      cout << Form("  Peak current = %6.3f %s",hX1->GetMaximum(),curSUnit.c_str()) << endl;
      cout << Form("  Total energy = %6.3f %s, rms = %3.1f %s",pzmean,eneSUnit.c_str(),(pzrms/pzmean)/ermsUnit,ermsSUnit.c_str()) << endl;
      cout << Form("  Corr. energy spread = %.3f %%/%s",zpzcorrrel,spaSUnit.c_str()) << endl;
      cout << Form("  Minimum energy = %.3f %s",MinP1 *  PConst::ElectronMassE / eneUnit,eneSUnit.c_str()) << endl;
      cout << Form("  Length = %6.3f %s (rms)",zrms,spaSUnit.c_str()) << endl;
      cout << Form("  Width x = %6.3f %s (rms)",x_rms,tspaSUnit.c_str()) << endl;
      cout << Form("  Trans. emit. x = %6.3f %s",emitx,emitSUnit.c_str()) << endl;
      if(pData->Is3D()) {
	cout << Form("  Width y = %6.3f %s (rms)",y_rms,tspaSUnit.c_str()) << endl;
	cout << Form("  Trans. emit. y = %6.3f %s",emity,emitSUnit.c_str()) << endl;
      }
    }
    
    if(opt.Contains("loop")) {
      cout << Form("\n 5. Saving results to file .. ") << endl;
  
      // OUTPUT ROOT FILE WITH THE PLOTS:
      TString filename = Form("./%s/Plots/Bunch/%s/Bunch-Evolution-%s.root",simo.Data(),pData->GetRawSpeciesName(index).c_str(),simo.Data());
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

      
      sprintf(gName,"gPzcorrrelvsTime");     
      TGraph *gPzcorrrelvsTime = NULL;
      gPzcorrrelvsTime = (TGraph*) ifile->Get(gName);
      if(gPzcorrrelvsTime==NULL) {
	gPzcorrrelvsTime = new TGraph();
	gPzcorrrelvsTime->SetName(gName);
	nPoints = 0;
	// Some cosmetics at creation time:
	gPzcorrrelvsTime->SetLineWidth(3);
	gPzcorrrelvsTime->SetLineColor(PGlobals::fieldLine);
	gPzcorrrelvsTime->SetMarkerStyle(20);
	gPzcorrrelvsTime->SetMarkerSize(0.4);
	gPzcorrrelvsTime->SetMarkerColor(PGlobals::fieldLine);	
      } else {
	nPoints = gPzcorrrelvsTime->GetN(); 
      }  

      gPzcorrrelvsTime->Set(nPoints+1);
      gPzcorrrelvsTime->SetPoint(nPoints,Time,zpzcorrrel);
      gPzcorrrelvsTime->Write(gName,TObject::kOverwrite);

      sprintf(gName,"gErmsavgvsTime");     
      TGraph *gErmsavgvsTime = NULL;
      gErmsavgvsTime = (TGraph*) ifile->Get(gName);
      if(gErmsavgvsTime==NULL) {
	gErmsavgvsTime = new TGraph();
	gErmsavgvsTime->SetName(gName);
	nPoints = 0;
	// Some cosmetics at creation time:
	gErmsavgvsTime->SetLineWidth(3);
	gErmsavgvsTime->SetLineColor(PGlobals::fieldLine);
	gErmsavgvsTime->SetMarkerStyle(20);
	gErmsavgvsTime->SetMarkerSize(0.4);
	gErmsavgvsTime->SetMarkerColor(PGlobals::fieldLine);	
      } else {
	nPoints = gErmsavgvsTime->GetN(); 
      }  

      gErmsavgvsTime->Set(nPoints+1);
      gErmsavgvsTime->SetPoint(nPoints,Time,sErmsavg);
      gErmsavgvsTime->Write(gName,TObject::kOverwrite);


      
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

      sprintf(gName,"gxPxrmsvsTime");     
      TGraph *gxPxrmsvsTime = NULL;
      gxPxrmsvsTime = (TGraph*) ifile->Get(gName);
      if(gxPxrmsvsTime==NULL) {
	gxPxrmsvsTime = new TGraph();
	gxPxrmsvsTime->SetName(gName);
	nPoints = 0;
	// Some cosmetics at creation time:
	gxPxrmsvsTime->SetLineWidth(3);
	gxPxrmsvsTime->SetLineColor(PGlobals::fieldLine);
	gxPxrmsvsTime->SetMarkerStyle(20);
	gxPxrmsvsTime->SetMarkerSize(0.4);
	gxPxrmsvsTime->SetMarkerColor(PGlobals::fieldLine);	
      } else {
	nPoints = gxPxrmsvsTime->GetN(); 
      }  

      gxPxrmsvsTime->Set(nPoints+1);
      gxPxrmsvsTime->SetPoint(nPoints,Time,xpx_rms);
      gxPxrmsvsTime->Write(gName,TObject::kOverwrite);

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

      sprintf(gName,"gXmeanavgvsTime");     
      TGraph *gXmeanavgvsTime = NULL;
      gXmeanavgvsTime = (TGraph*) ifile->Get(gName);
      if(gXmeanavgvsTime==NULL) {
	gXmeanavgvsTime = new TGraph();
	gXmeanavgvsTime->SetName(gName);
	nPoints = 0;
	// Some cosmetics at creation time:
	gXmeanavgvsTime->SetLineWidth(3);
	gXmeanavgvsTime->SetLineColor(PGlobals::fieldLine);
	gXmeanavgvsTime->SetMarkerStyle(20);
	gXmeanavgvsTime->SetMarkerSize(0.4);
	gXmeanavgvsTime->SetMarkerColor(PGlobals::fieldLine);	
      } else {
	nPoints = gXmeanavgvsTime->GetN(); 
      }  

      gXmeanavgvsTime->Set(nPoints+1);
      gXmeanavgvsTime->SetPoint(nPoints,Time,xmeanavg);
      gXmeanavgvsTime->Write(gName,TObject::kOverwrite);

      
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

      sprintf(gName,"gXrmsavgvsTime");     
      TGraph *gXrmsavgvsTime = NULL;
      gXrmsavgvsTime = (TGraph*) ifile->Get(gName);
      if(gXrmsavgvsTime==NULL) {
	gXrmsavgvsTime = new TGraph();
	gXrmsavgvsTime->SetName(gName);
	nPoints = 0;
	// Some cosmetics at creation time:
	gXrmsavgvsTime->SetLineWidth(3);
	gXrmsavgvsTime->SetLineColor(PGlobals::fieldLine);
	gXrmsavgvsTime->SetMarkerStyle(20);
	gXrmsavgvsTime->SetMarkerSize(0.4);
	gXrmsavgvsTime->SetMarkerColor(PGlobals::fieldLine);	
      } else {
	nPoints = gXrmsavgvsTime->GetN(); 
      }  

      gXrmsavgvsTime->Set(nPoints+1);
      gXrmsavgvsTime->SetPoint(nPoints,Time,xrmsavg);
      gXrmsavgvsTime->Write(gName,TObject::kOverwrite);


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

      // Core emittance (in slices range)
      sprintf(gName,"gEmitxcorevsTime");     
      TGraph *gEmitxcorevsTime = NULL;
      gEmitxcorevsTime = (TGraph*) ifile->Get(gName);
      if(gEmitxcorevsTime==NULL) {
	gEmitxcorevsTime = new TGraph();
	gEmitxcorevsTime->SetName(gName);
	nPoints = 0;
	// Some cosmetics at creation time:
	gEmitxcorevsTime->SetLineWidth(3);
	gEmitxcorevsTime->SetLineColor(kBlue-2);
	gEmitxcorevsTime->SetMarkerStyle(20);
	gEmitxcorevsTime->SetMarkerSize(0.4);
	gEmitxcorevsTime->SetMarkerColor(kBlue-2);	
      } else {
	nPoints = gEmitxcorevsTime->GetN(); 
      }  

      gEmitxcorevsTime->Set(nPoints+1);
      gEmitxcorevsTime->SetPoint(nPoints,Time,emitx_r);
      gEmitxcorevsTime->Write(gName,TObject::kOverwrite);

      // Sliced averaged emittance
      sprintf(gName,"gEmitxavgvsTime");     
      TGraph *gEmitxavgvsTime = NULL;
      gEmitxavgvsTime = (TGraph*) ifile->Get(gName);
      if(gEmitxavgvsTime==NULL) {
	gEmitxavgvsTime = new TGraph();
	gEmitxavgvsTime->SetName(gName);
	nPoints = 0;
	// Some cosmetics at creation time:
	gEmitxavgvsTime->SetLineWidth(3);
	gEmitxavgvsTime->SetLineColor(kGray+3);
	gEmitxavgvsTime->SetMarkerStyle(20);
	gEmitxavgvsTime->SetMarkerSize(0.4);
	gEmitxavgvsTime->SetMarkerColor(kGray+3);	
      } else {
	nPoints = gEmitxavgvsTime->GetN(); 
      }  
      gEmitxavgvsTime->Set(nPoints+1);
      gEmitxavgvsTime->SetPoint(nPoints,Time,emitxavg);
      gEmitxavgvsTime->Write(gName,TObject::kOverwrite);
      // ------
     
      sprintf(gName,"gPymeanvsTime");     
      TGraph *gPymeanvsTime = NULL;
      gPymeanvsTime = (TGraph*) ifile->Get(gName);
      if(gPymeanvsTime==NULL) {
	gPymeanvsTime = new TGraph();
	gPymeanvsTime->SetName(gName);
	nPoints = 0;
	// Some cosmetics at creation time:
	gPymeanvsTime->SetLineWidth(3);
	gPymeanvsTime->SetLineColor(PGlobals::fieldLine);
	gPymeanvsTime->SetMarkerStyle(20);
	gPymeanvsTime->SetMarkerSize(0.4);
	gPymeanvsTime->SetMarkerColor(PGlobals::fieldLine);	
      } else {
	nPoints = gPymeanvsTime->GetN(); 
      }  

      gPymeanvsTime->Set(nPoints+1);
      gPymeanvsTime->SetPoint(nPoints,Time,py_mean);
      gPymeanvsTime->Write(gName,TObject::kOverwrite);

      sprintf(gName,"gPyrmsvsTime");     
      TGraph *gPyrmsvsTime = NULL;
      gPyrmsvsTime = (TGraph*) ifile->Get(gName);
      if(gPyrmsvsTime==NULL) {
	gPyrmsvsTime = new TGraph();
	gPyrmsvsTime->SetName(gName);
	nPoints = 0;
	// Some cosmetics at creation time:
	gPyrmsvsTime->SetLineWidth(3);
	gPyrmsvsTime->SetLineColor(PGlobals::fieldLine);
	gPyrmsvsTime->SetMarkerStyle(20);
	gPyrmsvsTime->SetMarkerSize(0.4);
	gPyrmsvsTime->SetMarkerColor(PGlobals::fieldLine);	
      } else {
	nPoints = gPyrmsvsTime->GetN(); 
      }  

      gPyrmsvsTime->Set(nPoints+1);
      gPyrmsvsTime->SetPoint(nPoints,Time,py_rms);
      gPyrmsvsTime->Write(gName,TObject::kOverwrite);

      sprintf(gName,"gyPyrmsvsTime");     
      TGraph *gyPyrmsvsTime = NULL;
      gyPyrmsvsTime = (TGraph*) ifile->Get(gName);
      if(gyPyrmsvsTime==NULL) {
	gyPyrmsvsTime = new TGraph();
	gyPyrmsvsTime->SetName(gName);
	nPoints = 0;
	// Some cosmetics at creation time:
	gyPyrmsvsTime->SetLineWidth(3);
	gyPyrmsvsTime->SetLineColor(PGlobals::fieldLine);
	gyPyrmsvsTime->SetMarkerStyle(20);
	gyPyrmsvsTime->SetMarkerSize(0.4);
	gyPyrmsvsTime->SetMarkerColor(PGlobals::fieldLine);	
      } else {
	nPoints = gyPyrmsvsTime->GetN(); 
      }  

      gyPyrmsvsTime->Set(nPoints+1);
      gyPyrmsvsTime->SetPoint(nPoints,Time,ypy_rms);
      gyPyrmsvsTime->Write(gName,TObject::kOverwrite);

      sprintf(gName,"gYmeanvsTime");     
      TGraph *gYmeanvsTime = NULL;
      gYmeanvsTime = (TGraph*) ifile->Get(gName);
      if(gYmeanvsTime==NULL) {
	gYmeanvsTime = new TGraph();
	gYmeanvsTime->SetName(gName);
	nPoints = 0;
	// Some cosmetics at creation time:
	gYmeanvsTime->SetLineWidth(3);
	gYmeanvsTime->SetLineColor(PGlobals::fieldLine);
	gYmeanvsTime->SetMarkerStyle(20);
	gYmeanvsTime->SetMarkerSize(0.4);
	gYmeanvsTime->SetMarkerColor(PGlobals::fieldLine);	
      } else {
	nPoints = gYmeanvsTime->GetN(); 
      }  

      gYmeanvsTime->Set(nPoints+1);
      gYmeanvsTime->SetPoint(nPoints,Time,y_mean);
      gYmeanvsTime->Write(gName,TObject::kOverwrite);

      sprintf(gName,"gYmeanavgvsTime");     
      TGraph *gYmeanavgvsTime = NULL;
      gYmeanavgvsTime = (TGraph*) ifile->Get(gName);
      if(gYmeanavgvsTime==NULL) {
	gYmeanavgvsTime = new TGraph();
	gYmeanavgvsTime->SetName(gName);
	nPoints = 0;
	// Some cosmetics at creation time:
	gYmeanavgvsTime->SetLineWidth(3);
	gYmeanavgvsTime->SetLineColor(PGlobals::fieldLine);
	gYmeanavgvsTime->SetMarkerStyle(20);
	gYmeanavgvsTime->SetMarkerSize(0.4);
	gYmeanavgvsTime->SetMarkerColor(PGlobals::fieldLine);	
      } else {
	nPoints = gYmeanavgvsTime->GetN(); 
      }  

      gYmeanavgvsTime->Set(nPoints+1);
      gYmeanavgvsTime->SetPoint(nPoints,Time,ymeanavg);
      gYmeanavgvsTime->Write(gName,TObject::kOverwrite);


      sprintf(gName,"gYrmsvsTime");     
      TGraph *gYrmsvsTime = NULL;
      gYrmsvsTime = (TGraph*) ifile->Get(gName);
      if(gYrmsvsTime==NULL) {
	gYrmsvsTime = new TGraph();
	gYrmsvsTime->SetName(gName);
	nPoints = 0;
	// Some cosmetics at creation time:
	gYrmsvsTime->SetLineWidth(3);
	gYrmsvsTime->SetLineColor(PGlobals::fieldLine);
	gYrmsvsTime->SetMarkerStyle(20);
	gYrmsvsTime->SetMarkerSize(0.4);
	gYrmsvsTime->SetMarkerColor(PGlobals::fieldLine);	
      } else {
	nPoints = gYrmsvsTime->GetN(); 
      }  

      gYrmsvsTime->Set(nPoints+1);
      gYrmsvsTime->SetPoint(nPoints,Time,y_rms);
      gYrmsvsTime->Write(gName,TObject::kOverwrite);


      sprintf(gName,"gYrmsavgvsTime");     
      TGraph *gYrmsavgvsTime = NULL;
      gYrmsavgvsTime = (TGraph*) ifile->Get(gName);
      if(gYrmsavgvsTime==NULL) {
	gYrmsavgvsTime = new TGraph();
	gYrmsavgvsTime->SetName(gName);
	nPoints = 0;
	// Some cosmetics at creation time:
	gYrmsavgvsTime->SetLineWidth(3);
	gYrmsavgvsTime->SetLineColor(PGlobals::fieldLine);
	gYrmsavgvsTime->SetMarkerStyle(20);
	gYrmsavgvsTime->SetMarkerSize(0.4);
	gYrmsavgvsTime->SetMarkerColor(PGlobals::fieldLine);	
      } else {
	nPoints = gYrmsavgvsTime->GetN(); 
      }  

      gYrmsavgvsTime->Set(nPoints+1);
      gYrmsavgvsTime->SetPoint(nPoints,Time,yrmsavg);
      gYrmsavgvsTime->Write(gName,TObject::kOverwrite);

      
      sprintf(gName,"gEmityvsTime");     
      TGraph *gEmityvsTime = NULL;
      gEmityvsTime = (TGraph*) ifile->Get(gName);
      if(gEmityvsTime==NULL) {
	gEmityvsTime = new TGraph();
	gEmityvsTime->SetName(gName);
	nPoints = 0;
	// Some cosmetics at creation time:
	gEmityvsTime->SetLineWidth(3);
	gEmityvsTime->SetLineColor(PGlobals::fieldLine);
	gEmityvsTime->SetMarkerStyle(20);
	gEmityvsTime->SetMarkerSize(0.4);
	gEmityvsTime->SetMarkerColor(PGlobals::fieldLine);	
      } else {
	nPoints = gEmityvsTime->GetN(); 
      }  

      gEmityvsTime->Set(nPoints+1);
      gEmityvsTime->SetPoint(nPoints,Time,emity);
      gEmityvsTime->Write(gName,TObject::kOverwrite);

      // Core emittance (in slices range)
      sprintf(gName,"gEmitycorevsTime");     
      TGraph *gEmitycorevsTime = NULL;
      gEmitycorevsTime = (TGraph*) ifile->Get(gName);
      if(gEmitycorevsTime==NULL) {
	gEmitycorevsTime = new TGraph();
	gEmitycorevsTime->SetName(gName);
	nPoints = 0;
	// Some cosmetics at creation time:
	gEmitycorevsTime->SetLineWidth(3);
	gEmitycorevsTime->SetLineColor(kBlue-2);
	gEmitycorevsTime->SetMarkerStyle(20);
	gEmitycorevsTime->SetMarkerSize(0.4);
	gEmitycorevsTime->SetMarkerColor(kBlue-2);	
      } else {
	nPoints = gEmitycorevsTime->GetN(); 
      }  

      gEmitycorevsTime->Set(nPoints+1);
      gEmitycorevsTime->SetPoint(nPoints,Time,emity_r);
      gEmitycorevsTime->Write(gName,TObject::kOverwrite);

      // Sliced averaged emittance
      sprintf(gName,"gEmityavgvsTime");     
      TGraph *gEmityavgvsTime = NULL;
      gEmityavgvsTime = (TGraph*) ifile->Get(gName);
      if(gEmityavgvsTime==NULL) {
	gEmityavgvsTime = new TGraph();
	gEmityavgvsTime->SetName(gName);
	nPoints = 0;
	// Some cosmetics at creation time:
	gEmityavgvsTime->SetLineWidth(3);
	gEmityavgvsTime->SetLineColor(kGray+3);
	gEmityavgvsTime->SetMarkerStyle(20);
	gEmityavgvsTime->SetMarkerSize(0.4);
	gEmityavgvsTime->SetMarkerColor(kGray+3);	
      } else {
	nPoints = gEmityavgvsTime->GetN(); 
      }  

      gEmityavgvsTime->Set(nPoints+1);
      gEmityavgvsTime->SetPoint(nPoints,Time,emityavg);
      gEmityavgvsTime->Write(gName,TObject::kOverwrite);
      // ------

      // Ez track addon (for hosing studies)
      if(opt.Contains("eztrack")) {
	// Get 3D histogram in a zoomed area to save time
	Double_t x1minr = pData->GetXMin(0);
	Double_t x1maxr = pData->GetXMax(0);
	x1maxr = (x1minr + x1maxr)/2;
	Double_t x2minr = pData->GetXMin(1);
	Double_t x2maxr = pData->GetXMax(1);
	Double_t x2width = x2maxr - x2minr;
	Double_t x2center = (x2maxr + x2minr)/2.;
	x2minr = x2center - x2width/2.;
	x2maxr = x2center + x2width/2.;
	Double_t x3minr = pData->GetXMin(2);
	Double_t x3maxr = pData->GetXMax(2);
	Double_t x3width = x3maxr - x3minr;
	Double_t x3center = (x3maxr + x3minr)/2.;
	x3minr = x3center - x3width/2.;
	x3maxr = x3center + x3width/2.;

	pData->SetX1Min(x1minr);
	pData->SetX1Max(x1maxr);
	pData->SetX2Min(x2minr);
	pData->SetX2Max(x2maxr);
	pData->SetX3Min(x3minr);
	pData->SetX3Max(x3maxr);
	
	cout << Form("\n Reading 3D Ez info for maximum tracking") << endl;
 	
	TH3F *hEz3D = pData->GetH3(pData->GetEfieldFileName(0)->c_str(),"e1",opt);
	
	Double_t ezmin = hEz3D->GetMinimum();
	Int_t im,jm,km;
	hEz3D->GetMinimumBin(im,jm,km);
	Double_t zpos = hEz3D->GetXaxis()->GetBinCenter(im);
	Double_t xpos = hEz3D->GetYaxis()->GetBinCenter(jm);
	Double_t ypos = hEz3D->GetZaxis()->GetBinCenter(km);

	// Average the values using the nearest bins
	Double_t ezminavg = 0;
	Double_t zposavg = 0;
	Double_t xposavg = 0;
	Double_t yposavg = 0;
	Int_t nbins = 4;
	for(Int_t i=im;i<im+nbins;i++)
	  for(Int_t j=jm;j<jm+nbins;j++)
	    for(Int_t k=km;k<km+nbins;k++) {
	      Double_t ez = hEz3D->GetBinContent(i,j,k);
	      ezminavg += ez;
	      zposavg += hEz3D->GetXaxis()->GetBinCenter(i) * fabs(ez);
	      xposavg += hEz3D->GetYaxis()->GetBinCenter(j) * fabs(ez);
	      yposavg += hEz3D->GetZaxis()->GetBinCenter(k) * fabs(ez);
	    }

	zposavg /= fabs(ezminavg);
	xposavg /= fabs(ezminavg);
	yposavg /= fabs(ezminavg);
	ezminavg /= nbins*nbins*nbins;
	
	if(opt.Contains("units")) {
	  Double_t E0 = pData->GetPlasmaE0();
	  ezmin *= E0 / ezUnit; 
	  ezminavg *= E0 / ezUnit; 
	  cout << Form(" - Minimum value           = %.2e %s",ezmin,ezSUnit.c_str()) << endl;
	  cout << Form(" - Minimum (average) value = %.2e %s",ezminavg,ezSUnit.c_str()) << endl;

	  zpos *= skindepth / spaUnit;
	  xpos *= skindepth / tspaUnit;
	  ypos *= skindepth / tspaUnit;

	  zposavg *= skindepth / spaUnit;
	  xposavg *= skindepth / tspaUnit;
	  yposavg *= skindepth / tspaUnit;

	  cout << Form(" - Position           = (%.4f %s, %.4f %s, %.4f %s)",zpos
		       ,spaSUnit.c_str(),xpos,tspaSUnit.c_str(),ypos,tspaSUnit.c_str()) << endl;
	  cout << Form(" - Position (average) = (%.4f %s, %.4f %s, %.4f %s)",zposavg
		       ,spaSUnit.c_str(),xposavg,tspaSUnit.c_str(),yposavg,tspaSUnit.c_str()) << endl;
	}
		
	sprintf(gName,"gEzminvsTime");     
	TGraph *gEzminvsTime = NULL;
	gEzminvsTime = (TGraph*) ifile->Get(gName);
	if(gEzminvsTime==NULL) {
	  gEzminvsTime = new TGraph();
	  gEzminvsTime->SetName(gName);
	  nPoints = 0;
	  // Some cosmetics at creation time:
	  gEzminvsTime->SetLineWidth(3);
	  gEzminvsTime->SetLineColor(PGlobals::fieldLine);
	  gEzminvsTime->SetMarkerStyle(20);
	  gEzminvsTime->SetMarkerSize(0.4);
	  gEzminvsTime->SetMarkerColor(PGlobals::fieldLine);	
	} else {
	  nPoints = gEzminvsTime->GetN(); 
	}  

	gEzminvsTime->Set(nPoints+1);
	gEzminvsTime->SetPoint(nPoints,Time,ezminavg);
	gEzminvsTime->Write(gName,TObject::kOverwrite);

	sprintf(gName,"gxEzminvsTime");     
	TGraph *gxEzminvsTime = NULL;
	gxEzminvsTime = (TGraph*) ifile->Get(gName);
	if(gxEzminvsTime==NULL) {
	  gxEzminvsTime = new TGraph();
	  gxEzminvsTime->SetName(gName);
	  nPoints = 0;
	  // Some cosmetics at creation time:
	  gxEzminvsTime->SetLineWidth(3);
	  gxEzminvsTime->SetLineColor(PGlobals::fieldLine);
	  gxEzminvsTime->SetMarkerStyle(20);
	  gxEzminvsTime->SetMarkerSize(0.4);
	  gxEzminvsTime->SetMarkerColor(PGlobals::fieldLine);	
	} else {
	  nPoints = gxEzminvsTime->GetN(); 
	}  

	gxEzminvsTime->Set(nPoints+1);
	gxEzminvsTime->SetPoint(nPoints,Time,xposavg);
	gxEzminvsTime->Write(gName,TObject::kOverwrite);

	sprintf(gName,"gyEzminvsTime");     
	TGraph *gyEzminvsTime = NULL;
	gyEzminvsTime = (TGraph*) ifile->Get(gName);
	if(gyEzminvsTime==NULL) {
	  gyEzminvsTime = new TGraph();
	  gyEzminvsTime->SetName(gName);
	  nPoints = 0;
	  // Some cosmetics at creation time:
	  gyEzminvsTime->SetLineWidth(3);
	  gyEzminvsTime->SetLineColor(PGlobals::fieldLine);
	  gyEzminvsTime->SetMarkerStyle(20);
	  gyEzminvsTime->SetMarkerSize(0.4);
	  gyEzminvsTime->SetMarkerColor(PGlobals::fieldLine);	
	} else {
	  nPoints = gyEzminvsTime->GetN(); 
	}  

	gyEzminvsTime->Set(nPoints+1);
	gyEzminvsTime->SetPoint(nPoints,Time,yposavg);
	gyEzminvsTime->Write(gName,TObject::kOverwrite);


      }

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

    // Output file
    TString fOutName = Form("./%s/Plots/Bunch/%s/Bunch-%s-%s",
			    simo.Data(),pData->GetRawSpeciesName(index).c_str(),
			    pData->GetRawSpeciesName(index).c_str(),simo.Data());
    
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
    TGraph *gP1 = NULL;
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

      if(opt.Contains("smooth"))
	hP1->Smooth(3);
         
      for(Int_t j=0; j<p1Nbin; j++) {
	yarray[j] = hP1->GetBinCenter(j+1);
	xarray[j] = ((xMax-xMin)/EneMax)*hP1->GetBinContent(j+1) + xMin;
      }

      gP1 = new TGraph(p1Nbin,xarray,yarray);
      gP1->SetLineColor(PGlobals::elecLine);
      gP1->SetLineWidth(2);
      gP1->SetFillStyle(1001);
      gP1->SetFillColor(PGlobals::elecFill);

      delete[] yarray;
      delete[] xarray;

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
      
      if(Imax > -99999.) yMax = Imax;
      
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
      if(!pPalette) {
	pPalette = new PPalette("electron0");
	pPalette->SetPalette("electron0");
	//pPalette->SetAlpha(1);
      }

      // Text objects
      TPaveText *textTime =  new TPaveText(0.55,0.76,0.84,0.9,"NDC");
      PGlobals::SetPaveTextStyle(textTime,32);
      textTime->SetTextFont(43);
      textTime->SetTextSize(22);
      textTime->SetTextColor(kGray+2);
      char ctext[256];
      if(opt.Contains("units") && pData->GetPlasmaDensity()) { 
	sprintf(ctext,"z = %5.1f %s", Time, propSUnit.c_str());
      } else {
	sprintf(ctext,"k_{p}z = %5.1f",Time);
      }
      // cout << Form("z = %5.1f %s", Time, propSUnit.c_str()) << endl;
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

      TPaveText *textMom = new TPaveText(0.55,0.07,0.84,0.17,"NDC");
      PGlobals::SetPaveTextStyle(textMom,32); 
      textMom->SetTextColor(kGray+3);
      textMom->SetTextFont(63);
      textMom->SetTextSize(22);
      if(opt.Contains("fwhm")) {
	if(opt.Contains("units") && pData->GetPlasmaDensity())
	  sprintf(ctext,"#LTp_{z}#GT = %5.2f %s/c", pzmeanFWHM, eneSUnit.c_str());
	else
	  sprintf(ctext,"#LTp_{z}#GT = %5.2f mc", pzmeanFWHM);    
      } else {
	if(opt.Contains("units") && pData->GetPlasmaDensity())
	  sprintf(ctext,"#LTp_{z}#GT = %5.2f %s/c", pzmean, eneSUnit.c_str());
	else
	  sprintf(ctext,"#LTp_{z}#GT = %5.2f mc", pzmean);    
      }
      textMom->AddText(ctext);


      TPaveText *textInfo = new TPaveText(0.55,0.25,0.84,0.75,"NDC");
      PGlobals::SetPaveTextStyle(textInfo,32); 
      textInfo->SetTextColor(kGray+2);
      textInfo->SetTextFont(43);
      textInfo->SetTextSize(22);
      if(opt.Contains("units")) 
	sprintf(ctext,"Q = %5.2f %s",Charge,chargeSUnit.c_str());
      else
	sprintf(ctext,"Q = %5.2f Q_{0}",Charge);
      
      textInfo->AddText(ctext);
      if(opt.Contains("units")) 
	sprintf(ctext,"#sigma_{#zeta} = %5.2f %s",zrms,spaSUnit.c_str());
      else
	sprintf(ctext,"k_{p} #Delta#zeta = %5.2f",zrms);
      
      textInfo->AddText(ctext);

      if(opt.Contains("ecorr")) {
	sprintf(ctext,"#deltap_{z,corr} = %4.2f %%/%s",zpzcorrrel,spaSUnit.c_str());
      } else {
	if(opt.Contains("fwhm"))
	  sprintf(ctext,"#sigma_{#gamma}/#LT#gamma#GT = %4.2f %s",(pzrmsFWHM/pzmeanFWHM)/ermsUnit,ermsSUnit.c_str());
	else
	  sprintf(ctext,"#sigma_{#gamma}/#LT#gamma#GT = %4.1f %s",(pzrms/pzmean)/ermsUnit,ermsSUnit.c_str());
      }
      textInfo->AddText(ctext);
      
      if(opt.Contains("avge")) {
	if(opt.Contains("units"))
	  sprintf(ctext,"#LT#varepsilon_{n,x}#GT = %5.2f %s",emitxavg,emitSUnit.c_str());
	else
	  sprintf(ctext,"k_{p} #LT#varepsilon_{n,x}#GT = %5.2f",emitxavg);
	
	textInfo->AddText(ctext);
	if(pData->Is3D()) {
	  if(opt.Contains("units"))
	    sprintf(ctext,"#LT#varepsilon_{n,y}#GT = %5.2f %s",emityavg,emitSUnit.c_str());
	  else
	    sprintf(ctext,"k_{p} #LT#varepsilon_{n,y}#GT = %5.2f",emityavg);
	  
	  textInfo->AddText(ctext);
	}
      } else {
	if(opt.Contains("units"))
	  sprintf(ctext,"#varepsilon_{n,x} = %5.2f %s",emitx,emitSUnit.c_str());
	else
	  sprintf(ctext,"k_{p} #varepsilon_{n,x} = %5.2f",emitx);
	
	textInfo->AddText(ctext);
	if(pData->Is3D()) {
	  if(opt.Contains("units"))
	    sprintf(ctext,"#varepsilon_{n,y} = %5.2f %s",emity,emitSUnit.c_str());
	  else
	    sprintf(ctext,"k_{p} #varepsilon_{n,y} = %5.2f",emity);
	
	  textInfo->AddText(ctext);
	}
      }
      
      // Setup Pad layout: 
      const Int_t NPad = 2;
      TPad *pad[NPad];
      TH1F *hFrame[NPad];
      TString sLabels[] = {"(b)","(a)"};
      TPaveText **textLabel = new TPaveText*[NPad];

      Double_t lMargin = 0.12;
      Double_t rMargin = 0.15;
      Double_t bMargin = 0.15;
      Double_t tMargin = 0.04;
      Double_t factor = 1.0;
      Double_t gap = 0.028;
      PGlobals::CanvasAsymPartition(C,NPad,lMargin,rMargin,bMargin,tMargin,factor,gap);

      // Define the frames for plotting
      Int_t fonttype = 43;
      Int_t fontsize = 28;
      Int_t tfontsize = 30;
      Double_t txoffset = 2.0;
      Double_t lxoffset = 0.02;
      Double_t tyoffset = 1.0;
      Double_t lyoffset = 0.01;
      Double_t tylength = 0.015;
      Double_t txlength = 0.03;
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
	hFrame[i]->GetXaxis()->SetNdivisions(510);
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
	hFrame[1]->GetYaxis()->SetTitle(Form("p_{z} [%s/c]",eneSUnit.c_str()));
      else
	hFrame[1]->GetYaxis()->SetTitle("p_{z}/mc");
      
      hFrame[1]->Draw("axis");

      hP1X1->GetZaxis()->SetTitleFont(fonttype);
      hP1X1->GetZaxis()->SetTitleOffset(tyoffset);
      hP1X1->GetZaxis()->SetTitleSize(tfontsize);
      hP1X1->GetZaxis()->SetLabelFont(fonttype);
      hP1X1->GetZaxis()->SetLabelSize(fontsize);
      if(opt.Contains("logz")) 
	hP1X1->GetZaxis()->SetLabelOffset(0);
      else
	hP1X1->GetZaxis()->SetLabelOffset(lyoffset);
            
      hP1X1->GetZaxis()->SetTickLength(0.01);      

      hP1X1->Draw("colz 0 same");
      
      gP1->SetLineWidth(2);
      if(!opt.Contains("nospec")) {
	gP1->Draw("F");
	gP1->Draw("L");
      }
      
      Float_t emean = pzmean;
      if(opt.Contains("fwhm"))
	emean = pzmeanFWHM;

      TLine lZmean(zmean,hP1X1->GetYaxis()->GetXmin(),zmean,hP1X1->GetYaxis()->GetXmax());
      lZmean.SetLineColor(kGray+2);
      lZmean.SetLineStyle(2);
     
      TLine lPmean(hP1X1->GetXaxis()->GetXmin(),emean,hP1X1->GetXaxis()->GetXmax(),emean);
      lPmean.SetLineColor(kGray+2);
      lPmean.SetLineStyle(2);
     	
      // lines indicating the energy interval
      // TLine pzminline(hP1X1->GetXaxis()->GetXmin(),pzmin, ((xMax-xMin)/EneMax)*(EneMax/2) + xMin,pzmin);
      TLine pzminline(hP1X1->GetXaxis()->GetXmin(),pzmin,hP1X1->GetXaxis()->GetXmax(),pzmin);
      pzminline.SetLineColor(kGray+1);
      pzminline.SetLineStyle(3);
      
      TLine pzmaxline(hP1X1->GetXaxis()->GetXmin(),pzmax,hP1X1->GetXaxis()->GetXmax(),pzmax);
      pzmaxline.SetLineColor(kGray+1);
      pzmaxline.SetLineStyle(3);
      
      if(!opt.Contains("noline")) {
	lZmean.Draw();
	lPmean.Draw();
	
	if(opt.Contains("fwhm")) {
	  pzmaxline.Draw();
	  pzminline.Draw();
	}
      }    
      // hP1X1prof->SetMarkerStyle(1);
      // hP1X1prof->SetLineWidth(2);
      // hP1X1prof->Draw("zsame");

      gPad->Update();
      TPaletteAxis *palette = (TPaletteAxis*)hP1X1->GetListOfFunctions()->FindObject("palette");
      if(palette) {
	Double_t y1 = gPad->GetBottomMargin();
	Double_t y2 = 1 - gPad->GetTopMargin();
	Double_t x1 = 1 - gPad->GetRightMargin();

	Double_t x1b = x1 + 0.01;
	Double_t x2b = x1 + 0.035;
	Double_t y1b = y1 + 0.03;
	Double_t y2b = y2 - 0.03;
	palette->SetX1NDC(x1b);
	palette->SetY1NDC(y1b);
	palette->SetX2NDC(x2b);
	palette->SetY2NDC(y2b);
	palette->SetBorderSize(2);
	palette->SetLineColor(1);

	TPave *pFrame = new TPave(x1b,y1b,x2b,y2b,1,"NDCL");
	pFrame->SetFillStyle(0);
	pFrame->SetLineColor(kBlack);
	pFrame->SetShadowColor(0);
	pFrame->Draw();
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

      TBox *lFrame = new TBox(gPad->GetUxmin(), gPad->GetUymin(),
			      gPad->GetUxmax(), gPad->GetUymax());
      lFrame->SetFillStyle(0);
      lFrame->SetLineColor(PGlobals::frameColor);
      lFrame->SetLineWidth(PGlobals::frameWidth);
      lFrame->Draw();
      
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
      Leg=new TLegend(0.6,0.6,1 - gPad->GetRightMargin() - 0.04 - 0.03 ,0.95);

      PGlobals::SetPaveStyle(Leg);
      Leg->SetTextAlign(12);
      Leg->SetTextColor(kGray+3);
      Leg->SetTextFont(43);
      Leg->SetTextSize(20);
      Leg->SetLineColor(1);
      Leg->SetBorderSize(0);
      Leg->SetFillColor(0);
      Leg->SetFillStyle(1001);
      Leg->SetFillStyle(0); // Hollow

      char sleg[16];
      if(opt.Contains("units")) {
	sprintf(sleg,"Current [%s]",curSUnit.c_str());
	Leg->AddEntry(hX1  ,sleg,"L");	
	//sprintf(sleg,"Energy spread [%s]",ermsUnit.c_str());
	sprintf(sleg,"E. spread [%s]",ermssSUnit.c_str());
	Leg->AddEntry(gErms,sleg,"PL");
	//      sprintf(sleg,"Emittance [%s]",emitUnit.c_str());
	sprintf(sleg,"Emitt. x [%s]",emitSUnit.c_str());
	Leg->AddEntry(gEmitx,sleg,"PL");
	if(pData->Is3D()) {
	  sprintf(sleg,"Emitt. y [%s]",emitSUnit.c_str());
	  Leg->AddEntry(gEmity,sleg,"PL");
	}
	//Leg->AddEntry(gXrms,"Bunch width [#mum]","PL");
      } else {
	sprintf(sleg,"Current");
	Leg->AddEntry(hX1  ,sleg,"L");	
	sprintf(sleg,"E. spread [%s]",ermsSUnit.c_str());
	Leg->AddEntry(gErms,sleg,"PL");
	sprintf(sleg,"k_{p} Emitt. x");
	Leg->AddEntry(gEmitx,sleg,"PL");
	if(pData->Is3D()) {
	  sprintf(sleg,"k_{p} Emitt. y");
	  Leg->AddEntry(gEmity,sleg,"PL");
	}
	
      }

      hFrame[0]->GetYaxis()->SetTitle("");
      cout << Form(" - Imax = %.2f",Imax) << endl;
      hFrame[0]->GetYaxis()->SetRangeUser(0.001,yMax);
      hFrame[0]->Draw("axis");

      if(opt.Contains("smooth"))
	hX1->Smooth(3);
      hX1->Draw("hist LF same");

      TLine lZmean2(zmean,0.0,zmean,yMax);
      lZmean2.SetLineColor(kGray+2);
      lZmean2.SetLineStyle(2);
      if(!opt.Contains("noline")) {
	lZmean2.Draw();
      }

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

      gPad->Update();

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

      TBox *lFrame2 = new TBox(gPad->GetUxmin(), gPad->GetUymin(),
			       gPad->GetUxmax(), gPad->GetUymax());
      lFrame2->SetFillStyle(0);
      lFrame2->SetLineColor(PGlobals::frameColor);
      lFrame2->SetLineWidth(PGlobals::frameWidth);
      lFrame2->Draw();

      gPad->RedrawAxis(); 

      // Print to file --------------------------------------

      C->cd();
      
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

      if(pData->Is3D()) {
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
	
	PGlobals::CanvasAsymPartition(C1,NPad1,lMargin,rMargin,bMargin,tMargin,factor,gap);

	Float_t x_min = hX2X1->GetYaxis()->GetXmin();
	Float_t x_max = hX2X1->GetYaxis()->GetXmax();
	
	Float_t y_min = hX3X1->GetYaxis()->GetXmin();
	Float_t y_max = hX3X1->GetYaxis()->GetXmax();
	
	if(x_min<y_min) y_min = x_min;
	else x_min = y_min;
	if(x_max>y_max) y_max = x_max;
	else x_max = y_max;

	for(Int_t i=0;i<NPad1;i++) {
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
	  if(i==0) {
	    hFrame[i] = new TH1F(name,"",100,x1Min,x1Max);
	    hFrame[i]->GetYaxis()->SetRangeUser(x_min,x_max);
	  } else {
	    hFrame[i] = new TH1F(name,"",100,x1Min,x1Max);
	    hFrame[i]->GetYaxis()->SetRangeUser(y_min,y_max);
	  }
	  
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

	hFrame[0]->GetYaxis()->SetRangeUser(x_min,x_max);
	hFrame[0]->GetYaxis()->SetTitle(hX2X1->GetYaxis()->GetTitle());
	hFrame[0]->GetXaxis()->SetRangeUser(hX2X1->GetXaxis()->GetXmin(),hX2X1->GetXaxis()->GetXmax());
	hFrame[0]->GetXaxis()->SetTitle(hX2X1->GetXaxis()->GetTitle());
	hFrame[0]->Draw("axis");

	TH2F *hX2X1cl = (TH2F*) hX2X1->Clone("hX2X1cl");
	Float_t xFactor = pad[0]->GetAbsWNDC()/pad[0]->GetAbsWNDC();
	Float_t yFactor = pad[0]->GetAbsHNDC()/pad[0]->GetAbsHNDC();

	hX2X1cl->GetZaxis()->CenterTitle();
	hX2X1cl->GetZaxis()->SetTitleFont(fonttype);
	hX2X1cl->GetZaxis()->SetTickLength(xFactor*tylength/yFactor);
	hX2X1cl->GetZaxis()->SetTitleOffset(tyoffset);
	hX2X1cl->GetZaxis()->SetTitleSize(tfontsize);
	hX2X1cl->GetZaxis()->SetLabelFont(fonttype);
	hX2X1cl->GetZaxis()->SetLabelSize(fontsize);
	if(opt.Contains("logz")) 
	  hX2X1cl->GetZaxis()->SetLabelOffset(0);
	else
	  hX2X1cl->GetZaxis()->SetLabelOffset(lyoffset);
	
	hX2X1cl->Draw("colz 0 same");

	TLine lX1mean(zmean,x_min,zmean,x_max);
	lX1mean.SetLineColor(kGray+2);
	lX1mean.SetLineStyle(2);
	lX1mean.Draw();

	TLine lX2mean(hFrame[0]->GetXaxis()->GetXmin(),x_mean,hFrame[0]->GetXaxis()->GetXmax(),x_mean);
	lX2mean.SetLineColor(kGray+2);
	lX2mean.SetLineStyle(2);
	lX2mean.Draw();

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

	  Double_t x1b = x1 + 0.01;
	  Double_t x2b = x1 + 0.035;
	  Double_t y1b = y1 + 0.03;
	  Double_t y2b = y2 - 0.03;
	  palette->SetX1NDC(x1b);
	  palette->SetY1NDC(y1b);
	  palette->SetX2NDC(x2b);
	  palette->SetY2NDC(y2b);
	  palette->SetBorderSize(2);
	  palette->SetLineColor(1);
	  
	  TPave *pFrame = new TPave(x1b,y1b,x2b,y2b,1,"NDCL");
	  pFrame->SetFillStyle(0);
	  pFrame->SetLineColor(kBlack);
	  pFrame->SetShadowColor(0);
	  pFrame->Draw();
	  
	}

	TPaveText *textInfoX2X1 = new TPaveText(x1+0.02*xrange,y2-0.50*yrange,
						x1+0.20*xrange,y2-0.01*yrange,"NDC");
	PGlobals::SetPaveTextStyle(textInfoX2X1,12); 
	textInfoX2X1->SetTextColor(kGray+3);
	textInfoX2X1->SetTextFont(43);

	char text[64];
	if(opt.Contains("units")) {
	  sprintf(text,"Q = %5.1f %s",Charge,chargeSUnit.c_str());
	  textInfoX2X1->AddText(text);
	  sprintf(text,"#Delta#zeta = %5.2f %s",zrms,spaSUnit.c_str());
	  textInfoX2X1->AddText(text);
	  sprintf(text,"#Deltax = %5.2f %s",x_rms,tspaSUnit.c_str());
	  textInfoX2X1->AddText(text);
	  sprintf(text,"#varepsilon_{x} = %5.2f %s",emitx,emitSUnit.c_str());
	  textInfoX2X1->AddText(text);
	  sprintf(text,"#beta_{x} = %5.2f %s",betax,betaSUnit.c_str());
	  textInfoX2X1->AddText(text);
	} else {
	  sprintf(text,"Q = %5.1f Q_{0}",Charge);
	  textInfoX2X1->AddText(text);
	  sprintf(text,"k_{p}#Delta#zeta = %5.2f",zrms);
	  textInfoX2X1->AddText(text);
	  sprintf(text,"k_{p}#Deltax = %5.2f",x_rms);
	  textInfoX2X1->AddText(text);
	  sprintf(text,"k_{p}#varepsilon_{x} = %5.2f",emitx);
	  textInfoX2X1->AddText(text);
	  sprintf(text,"k_{p}#beta_{x} = %5.2f",betax);
	  textInfoX2X1->AddText(text);
	  
	}
	  
	textInfoX2X1->Draw();

	TBox *lFrame = new TBox(gPad->GetUxmin(), gPad->GetUymin(),
				gPad->GetUxmax(), gPad->GetUymax());
	lFrame->SetFillStyle(0);
	lFrame->SetLineColor(PGlobals::frameColor);
	lFrame->SetLineWidth(PGlobals::frameWidth);
	lFrame->Draw();
	
	gPad->RedrawAxis(); 

	C1->cd(0);

	pad[1]->Draw();
	pad[1]->cd(); 
	
	if(opt.Contains("logz")) {
	  gPad->SetLogz(1);
	} else {
	  gPad->SetLogz(0);
	}
	
	hFrame[1]->GetYaxis()->SetRangeUser(y_min,y_max);
	hFrame[1]->GetYaxis()->SetTitle(hX3X1->GetYaxis()->GetTitle());
	hFrame[1]->GetXaxis()->SetRangeUser(hX3X1->GetXaxis()->GetXmin(),hX3X1->GetXaxis()->GetXmax());
	hFrame[1]->GetXaxis()->SetTitle(hX3X1->GetXaxis()->GetTitle());
	hFrame[1]->Draw("axis");
	
	TH2F *hX3X1cl = (TH2F*) hX3X1->Clone("hX3X1cl");
	xFactor = pad[0]->GetAbsWNDC()/pad[1]->GetAbsWNDC();
	yFactor = pad[0]->GetAbsHNDC()/pad[1]->GetAbsHNDC();

	hX3X1cl->GetZaxis()->CenterTitle();
	hX3X1cl->GetZaxis()->SetTitleFont(fonttype);
	hX3X1cl->GetZaxis()->SetTickLength(xFactor*tylength/yFactor);
	hX3X1cl->GetZaxis()->SetTitleOffset(tyoffset);
	hX3X1cl->GetZaxis()->SetTitleSize(tfontsize);
	hX3X1cl->GetZaxis()->SetLabelFont(fonttype);
	hX3X1cl->GetZaxis()->SetLabelSize(fontsize);
	if(opt.Contains("logz")) 
	  hX3X1cl->GetZaxis()->SetLabelOffset(0);
	else
	  hX3X1cl->GetZaxis()->SetLabelOffset(lyoffset);
	
	hX3X1cl->Draw("colz same");

	TLine lX1mean2(zmean,y_min,zmean,y_max);
	lX1mean2.SetLineColor(kGray+2);
	lX1mean2.SetLineStyle(2);
	lX1mean2.Draw();

	TLine lX3mean(hFrame[1]->GetXaxis()->GetXmin(),y_mean,hFrame[1]->GetXaxis()->GetXmax(),y_mean);
	lX3mean.SetLineColor(kGray+2);
	lX3mean.SetLineStyle(2);
	lX3mean.Draw();

	
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

	  Double_t x1b = x1 + 0.01;
	  Double_t x2b = x1 + 0.035;
	  Double_t y1b = y1 + 0.03;
	  Double_t y2b = y2 - 0.03;
	  palette->SetX1NDC(x1b);
	  palette->SetY1NDC(y1b);
	  palette->SetX2NDC(x2b);
	  palette->SetY2NDC(y2b);
	  palette->SetBorderSize(2);
	  palette->SetLineColor(1);
	  
	  TPave *pFrame = new TPave(x1b,y1b,x2b,y2b,1,"NDCL");
	  pFrame->SetFillStyle(0);
	  pFrame->SetLineColor(kBlack);
	  pFrame->SetShadowColor(0);
	  pFrame->Draw();
	}

	TPaveText *textInfoX3X1 =  new TPaveText(x1+0.02*xrange,y2-0.40*yrange,
						 x1+0.20*xrange,y2-0.05*yrange,"NDC");
	PGlobals::SetPaveTextStyle(textInfoX3X1,12); 
	textInfoX3X1->SetTextColor(kGray+3);
	textInfoX3X1->SetTextFont(42);

	if(opt.Contains("units")) {
	  sprintf(text,"Q = %5.1f %s",Charge,chargeSUnit.c_str());
	  textInfoX3X1->AddText(text);
	  sprintf(text,"#Delta#zeta = %5.2f %s",zrms,spaSUnit.c_str());
	  textInfoX3X1->AddText(text);
	  sprintf(text,"#Deltay = %5.2f %s",y_rms,tspaSUnit.c_str());
	  textInfoX3X1->AddText(text);
	  sprintf(text,"#varepsilon_{y} = %5.2f %s",emity,emitSUnit.c_str());
	  textInfoX3X1->AddText(text);
	  sprintf(text,"#beta_{y} = %5.2f %s",betay,betaSUnit.c_str());
	  textInfoX3X1->AddText(text);
	} else {
	  sprintf(text,"Q = %5.1f Q_{0}",Charge);
	  textInfoX3X1->AddText(text);
	  sprintf(text,"k_{p}#Delta#zeta = %5.2f",zrms);
	  textInfoX3X1->AddText(text);
	  sprintf(text,"k_{p}#Deltay = %5.2f",y_rms);
	  textInfoX3X1->AddText(text);
	  sprintf(text,"k_{p}#varepsilon_{y} = %5.2f",emity);
	  textInfoX3X1->AddText(text);
	  sprintf(text,"k_{p}#beta_{y} = %5.2f",betay);
	  textInfoX3X1->AddText(text);
	  
	}
	  
	textInfoX3X1->Draw();

	TBox *lFrame2 = new TBox(gPad->GetUxmin(), gPad->GetUymin(),
				 gPad->GetUxmax(), gPad->GetUymax());
	lFrame2->SetFillStyle(0);
	lFrame2->SetLineColor(PGlobals::frameColor);
	lFrame2->SetLineWidth(PGlobals::frameWidth);
	lFrame2->Draw();

	gPad->RedrawAxis(); 
	
	C1->cd();
      
	TString fOutName1 =  fOutName + Form("-%s_%i","x2x3x1",time);
	PGlobals::imgconv(C1,fOutName1,opt);
      } else {
	sprintf(cName,"C1");     
	TCanvas *C1 = (TCanvas*) gROOT->FindObject(cName);
	if(C1==NULL) C1 = new TCanvas("C1","Space x2-x1 and x3-x1",sizex,sizey);
	C1->cd();
	C1->Clear();
	
	Int_t NPad1 = 1;
	lMargin = 0.15;
	rMargin = 0.18;
	bMargin = 0.15;
	tMargin = 0.04;
	factor = 1.0;  
	txoffset = 1.2;  
      
	PGlobals::CanvasAsymPartition(C1,NPad1,lMargin,rMargin,bMargin,tMargin,factor,gap);
	
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


	Float_t x_min = hX2X1->GetYaxis()->GetXmin();
	Float_t x_max = hX2X1->GetYaxis()->GetXmax();

	TLine lX1mean3, lX2mean2;
      
	C1->cd(0);
	pad[0]->Draw();
	pad[0]->cd(); 

	if(opt.Contains("logz")) {
	  gPad->SetLogz(1);
	} else {
	  gPad->SetLogz(0);
	}

	hFrame[0]->GetYaxis()->SetRangeUser(x_min,x_max);
	hFrame[0]->GetYaxis()->SetTitle(hX2X1->GetYaxis()->GetTitle());
	hFrame[0]->GetXaxis()->SetRangeUser(hX2X1->GetXaxis()->GetXmin(),hX2X1->GetXaxis()->GetXmax());
	hFrame[0]->GetXaxis()->SetTitle(hX2X1->GetXaxis()->GetTitle());
	hFrame[0]->Draw("axis");

	lX1mean3.SetLineColor(kGray+2);
	lX1mean3.SetLineStyle(2);
	lX1mean3.DrawLine(zmean,hFrame[0]->GetYaxis()->GetXmin(),zmean,hFrame[0]->GetYaxis()->GetXmax());
      
	lX2mean2.SetLineColor(kGray+2);
	lX2mean2.SetLineStyle(2);
	lX2mean2.DrawLine(hFrame[0]->GetXaxis()->GetXmin(),x_mean,hFrame[0]->GetXaxis()->GetXmax(),x_mean);


	TH2F *hX2X1cl = (TH2F*) hX2X1->Clone("hX2X1cl");
	Float_t xFactor = pad[0]->GetAbsWNDC()/pad[0]->GetAbsWNDC();
	Float_t yFactor = pad[0]->GetAbsHNDC()/pad[0]->GetAbsHNDC();

  	hX2X1cl->GetZaxis()->CenterTitle();
	hX2X1cl->GetZaxis()->SetTitleFont(fonttype);
	hX2X1cl->GetZaxis()->SetTickLength(xFactor*tylength/yFactor);
	hX2X1cl->GetZaxis()->SetTitleOffset(tyoffset);
	hX2X1cl->GetZaxis()->SetTitleSize(tfontsize);
	hX2X1cl->GetZaxis()->SetLabelFont(fonttype);
	hX2X1cl->GetZaxis()->SetLabelSize(fontsize);
	if(opt.Contains("logz")) 
	  hX2X1cl->GetZaxis()->SetLabelOffset(0);
	else
	  hX2X1cl->GetZaxis()->SetLabelOffset(lyoffset);

	hX2X1cl->Draw("colz same");

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

	TBox *lFrame = new TBox(gPad->GetUxmin(), gPad->GetUymin(),
				gPad->GetUxmax(), gPad->GetUymax());
	lFrame->SetFillStyle(0);
	lFrame->SetLineColor(PGlobals::frameColor);
	lFrame->SetLineWidth(PGlobals::frameWidth);
	lFrame->Draw();

	gPad->RedrawAxis(); 
	
	C1->cd();

	TString fOutName1 =  fOutName + Form("-%s_%i","x2x1",time);
	PGlobals::imgconv(C1,fOutName1,opt);
      }
      
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
	PGlobals::CanvasAsymPartition(CX3X2,NPadX3X2,lMargin,rMargin,bMargin,tMargin,factor,gap);

	for(Int_t i=0;i<NPadX3X2;i++) {
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

	TH2F *hX3X2cl = (TH2F*) hX3X2->Clone("hX3X2cl");
	hX3X2cl->GetZaxis()->CenterTitle();
	hX3X2cl->GetZaxis()->SetTitleFont(fonttype);
	hX3X2cl->GetZaxis()->SetTitleOffset(tyoffset);
	hX3X2cl->GetZaxis()->SetTitleSize(tfontsize);
	hX3X2cl->GetZaxis()->SetLabelFont(fonttype);
	hX3X2cl->GetZaxis()->SetLabelSize(fontsize);
	if(opt.Contains("logz")) 
	  hX3X2cl->GetZaxis()->SetLabelOffset(0);
	else
	  hX3X2cl->GetZaxis()->SetLabelOffset(lyoffset);
	
	hX3X2cl->GetZaxis()->SetTickLength(0.01);      

	hX3X2cl->Draw("colz 0 same");

	TLine lX3mean;
	lX3mean.SetLineColor(kGray+2);
	lX3mean.SetLineStyle(2);
	lX3mean.DrawLine(hFrame[0]->GetXaxis()->GetXmin(),y_mean,hFrame[0]->GetXaxis()->GetXmax(),y_mean);
	
	TLine lX2mean;
	lX2mean.SetLineColor(kGray+2);
	lX2mean.SetLineStyle(2);
	lX2mean.DrawLine(x_mean,hFrame[0]->GetYaxis()->GetXmin(),x_mean,hFrame[0]->GetYaxis()->GetXmax());
	
	gPad->Update();
	palette = (TPaletteAxis*)hX3X2cl->GetListOfFunctions()->FindObject("palette");
	if(palette) {
	  Double_t y1 = gPad->GetBottomMargin();
	  Double_t y2 = 1 - gPad->GetTopMargin();
	  Double_t x1 = 1 - gPad->GetRightMargin();
	  
	  Double_t x1b = x1 + 0.01;
	  Double_t x2b = x1 + 0.035;
	  Double_t y1b = y1 + 0.03;
	  Double_t y2b = y2 - 0.03;
	  palette->SetX1NDC(x1b);
	  palette->SetY1NDC(y1b);
	  palette->SetX2NDC(x2b);
	  palette->SetY2NDC(y2b);
	  palette->SetBorderSize(2);
	  palette->SetLineColor(1);
	  
	  TPave *pFrame = new TPave(x1b,y1b,x2b,y2b,1,"NDCL");
	  pFrame->SetFillStyle(0);
	  pFrame->SetLineColor(kBlack);
	  pFrame->SetShadowColor(0);
	  pFrame->Draw();
	}

	y1 = gPad->GetBottomMargin();
	y2 = 1 - gPad->GetTopMargin();
	x1 = gPad->GetLeftMargin();
	x2 = 1 - gPad->GetRightMargin();
	yrange = y2-y1; 
	xrange = x2-x1; 
    
	TPaveText *textStatX3X2 =  new TPaveText(x1+0.02*xrange,y2-0.40*yrange,x1+0.30*xrange,y2-0.05*yrange,"NDC");
	PGlobals::SetPaveTextStyle(textStatX3X2,12); 
	textStatX3X2->SetTextColor(kGray+3);
	textStatX3X2->SetTextFont(42);

	char text[64];
	if(opt.Contains("units")) {
	  sprintf(text,"Q = %5.1f %s",Charge,chargeSUnit.c_str());
	  textStatX3X2->AddText(text);
	  sprintf(text,"#Deltax = %5.2f %s",x_rms,tspaSUnit.c_str());
	  textStatX3X2->AddText(text);
	  sprintf(text,"#Deltay = %5.2f %s",y_rms,tspaSUnit.c_str());
	  textStatX3X2->AddText(text);
	  sprintf(text,"#varepsilon_{x} = %5.2f %s",emitx,emitSUnit.c_str());
	  textStatX3X2->AddText(text);
	  sprintf(text,"#varepsilon_{y} = %5.2f %s",emity,emitSUnit.c_str());
	  textStatX3X2->AddText(text);
	  sprintf(text,"#beta_{x} = %5.2f %s",betax,betaSUnit.c_str());
	  textStatX3X2->AddText(text);
	  sprintf(text,"#beta_{y} = %5.2f %s",betay,betaSUnit.c_str());
	  textStatX3X2->AddText(text);
	} else {
	  sprintf(text,"Q = %5.1f %s",Charge,chargeSUnit.c_str());
	  textStatX3X2->AddText(text);
	  sprintf(text,"k_{p}#Deltax = %5.2f",x_rms);
	  textStatX3X2->AddText(text);
	  sprintf(text,"k_{p}#Deltay = %5.2f",y_rms);
	  textStatX3X2->AddText(text);
	  sprintf(text,"k_{p}#varepsilon_{x} = %5.2f",emitx);
	  textStatX3X2->AddText(text);
	  sprintf(text,"k_{p}#varepsilon_{y} = %5.2f",emity);
	  textStatX3X2->AddText(text);
	  sprintf(text,"k_{p}#beta_{x} = %5.2f",betax);
	  textStatX3X2->AddText(text);
	  sprintf(text,"k_{p}#beta_{y} = %5.2f",betay);
	  textStatX3X2->AddText(text);	  
	}
	textStatX3X2->Draw();

	TBox *lFrame = new TBox(gPad->GetUxmin(), gPad->GetUymin(),
				gPad->GetUxmax(), gPad->GetUymax());
	lFrame->SetFillStyle(0);
	lFrame->SetLineColor(PGlobals::frameColor);
	lFrame->SetLineWidth(PGlobals::frameWidth);
	lFrame->Draw();
	
	gPad->RedrawAxis(); 
	
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
	pad[i]->SetTickx(0);
	pad[i]->SetTicky(0);
	if(opt.Contains("trans"))
	  pad[i]->SetFillStyle(4000);
	pad[i]->SetFrameFillStyle(4000);
	

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

      TH2F *hP2X2cl = (TH2F*) hP2X2->Clone("hP2X2cl");

      hP2X2cl->GetZaxis()->CenterTitle();
      hP2X2cl->GetZaxis()->SetTitleFont(fonttype);
      hP2X2cl->GetZaxis()->SetTitleOffset(tyoffset);
      hP2X2cl->GetZaxis()->SetTitleSize(tfontsize);
      hP2X2cl->GetZaxis()->SetLabelFont(fonttype);
      hP2X2cl->GetZaxis()->SetLabelSize(fontsize);
      if(opt.Contains("logz")) 
	hP2X2cl->GetZaxis()->SetLabelOffset(0);
      else
	hP2X2cl->GetZaxis()->SetLabelOffset(lyoffset);
            
      hP2X2cl->GetZaxis()->SetTickLength(0.01);      
      
      hP2X2cl->Draw("colz same");

      TLine lXmean(x_mean,hFrame[0]->GetYaxis()->GetXmin(),x_mean,hFrame[0]->GetYaxis()->GetXmax());
      lXmean.SetLineColor(kGray+2);
      lXmean.SetLineStyle(2);
      lXmean.Draw();

      TLine lPxmean(hFrame[0]->GetXaxis()->GetXmin(),px_mean,hFrame[0]->GetXaxis()->GetXmax(),px_mean);
      lPxmean.SetLineColor(kGray+2);
      lPxmean.SetLineStyle(2);
      lPxmean.Draw();

    
      if(opt.Contains("elli") ) 
	for(Int_t k=0;k<SNbin;k++)
	  if(sellipP2X2[k] != NULL)
	    sellipP2X2[k]->Draw();

      gPad->Update();
      palette = (TPaletteAxis*)hP2X2cl->GetListOfFunctions()->FindObject("palette");
      if(palette) {
	Double_t y1 = gPad->GetBottomMargin();
	Double_t y2 = 1 - gPad->GetTopMargin();
	Double_t x1 = 1 - gPad->GetRightMargin();
	
	Double_t x1b = x1 + 0.01;
	Double_t x2b = x1 + 0.035;
	Double_t y1b = y1 + 0.03;
	Double_t y2b = y2 - 0.03;
	palette->SetX1NDC(x1b);
	palette->SetY1NDC(y1b);
	palette->SetX2NDC(x2b);
	palette->SetY2NDC(y2b);
	palette->SetBorderSize(2);
	palette->SetLineColor(1);

	TPave *pFrame = new TPave(x1b,y1b,x2b,y2b,1,"NDCL");
	pFrame->SetFillStyle(0);
	pFrame->SetLineColor(kBlack);
	pFrame->SetShadowColor(0);
	pFrame->Draw();
      }

      TPaveText *textStatInt = new TPaveText(x1+0.02,y2-0.38,x1+0.28,y2-0.01,"NDC");
      PGlobals::SetPaveTextStyle(textStatInt,12); 
      textStatInt->SetTextColor(kGray+3);
      textStatInt->SetTextFont(42);

      char text[64];
      if(opt.Contains("units")) {
	sprintf(text,"Q = %5.1f %s",Charge,chargeSUnit.c_str());
	textStatInt->AddText(text);
	sprintf(text,"#sigma_{x} = %5.2f %s",x_rms,tspaSUnit.c_str());
	textStatInt->AddText(text);
	sprintf(text,"#sigma_{p_{x}} = %5.2f %s/c",px_rms,teneSUnit.c_str());
	textStatInt->AddText(text);
	sprintf(text,"#varepsilon_{x} = %5.2f %s",emitx,emitSUnit.c_str());
	textStatInt->AddText(text);
	sprintf(text,"#beta_{x} = %5.2f %s",betax,betaSUnit.c_str());
	textStatInt->AddText(text);
	sprintf(text,"#gamma_{x} = %5.2f %s",gammax,gammaSUnit.c_str());
	textStatInt->AddText(text);
	sprintf(text,"#alpha_{x} = %5.2f",alphax);
	textStatInt->AddText(text);	
      } else {
	sprintf(text,"Q = %5.1f Q_{0}",Charge);
	textStatInt->AddText(text);
	sprintf(text,"k_{p}#Deltax = %5.2f",x_rms);
	textStatInt->AddText(text);
	sprintf(text,"#Deltap_{x} = %5.2f mc",px_rms);
	textStatInt->AddText(text);
	sprintf(text,"k_{p}#varepsilon_{x} = %5.2f",emitx);
	textStatInt->AddText(text);
	sprintf(text,"k_{p}#beta_{x} = %5.2f",betax);
	textStatInt->AddText(text);
      }
      
      textStatInt->Draw();

      TBox *lFrame3 = new TBox(gPad->GetUxmin(), gPad->GetUymin(),
			       gPad->GetUxmax(), gPad->GetUymax());
      lFrame3->SetFillStyle(0);
      lFrame3->SetLineColor(PGlobals::frameColor);
      lFrame3->SetLineWidth(PGlobals::frameWidth);
      lFrame3->Draw();
      
      gPad->RedrawAxis(); 
    
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
	  pad[i]->SetTickx(0);
	  pad[i]->SetTicky(0);
	  if(opt.Contains("trans"))
	    pad[i]->SetFillStyle(4000);
	  pad[i]->SetFrameFillStyle(4000);
	  
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

	TH2F *hP3X3cl = (TH2F*) hP3X3->Clone("hP3X3cl");
	hP3X3cl->GetZaxis()->CenterTitle();
	hP3X3cl->GetZaxis()->SetTitleFont(fonttype);
	hP3X3cl->GetZaxis()->SetTitleOffset(tyoffset);
	hP3X3cl->GetZaxis()->SetTitleSize(tfontsize);
	hP3X3cl->GetZaxis()->SetLabelFont(fonttype);
	hP3X3cl->GetZaxis()->SetLabelSize(fontsize);
	if(opt.Contains("logz")) 
	  hP3X3cl->GetZaxis()->SetLabelOffset(0);
	else
	  hP3X3cl->GetZaxis()->SetLabelOffset(lyoffset);
	
	hP3X3cl->GetZaxis()->SetTickLength(0.01);      
	
	hP3X3cl->Draw("colz same");
	
	TLine lYmean(y_mean,hFrame[0]->GetYaxis()->GetXmin(),y_mean,hFrame[0]->GetYaxis()->GetXmax());
	lYmean.SetLineColor(kGray+2);
	lYmean.SetLineStyle(2);
	lYmean.Draw();

	TLine lPymean(hFrame[0]->GetXaxis()->GetXmin(),py_mean,hFrame[0]->GetXaxis()->GetXmax(),py_mean);
	lPymean.SetLineColor(kGray+2);
	lPymean.SetLineStyle(2);
	lPymean.Draw();

        
	if(opt.Contains("elli") ) 
	  for(Int_t k=0;k<SNbin;k++)
	    if(sellipP3X3[k] != NULL)
	      sellipP3X3[k]->Draw();

	gPad->Update();
	palette = (TPaletteAxis*)hP3X3cl->GetListOfFunctions()->FindObject("palette");
	if(palette) {
	  Double_t y1 = gPad->GetBottomMargin();
	  Double_t y2 = 1 - gPad->GetTopMargin();
	  Double_t x1 = 1 - gPad->GetRightMargin();
	
	  Double_t x1b = x1 + 0.01;
	  Double_t x2b = x1 + 0.035;
	  Double_t y1b = y1 + 0.03;
	  Double_t y2b = y2 - 0.03;
	  palette->SetX1NDC(x1b);
	  palette->SetY1NDC(y1b);
	  palette->SetX2NDC(x2b);
	  palette->SetY2NDC(y2b);
	  palette->SetBorderSize(2);
	  palette->SetLineColor(1);

	  TPave *pFrame = new TPave(x1b,y1b,x2b,y2b,1,"NDCL");
	  pFrame->SetFillStyle(0);
	  pFrame->SetLineColor(kBlack);
	  pFrame->SetShadowColor(0);
	  pFrame->Draw();
	}

	textStatInt = new TPaveText(x1+0.02,y2-0.30,x1+0.20,y2-0.05,"NDC");
	PGlobals::SetPaveTextStyle(textStatInt,12); 
	textStatInt->SetTextColor(kGray+3);
	textStatInt->SetTextFont(42);

	if(opt.Contains("units")) {
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
	} else {
	  sprintf(text,"Q = %5.1f Q_{0}",Charge);
	  textStatInt->AddText(text);
	  sprintf(text,"k_{p}#Deltay = %5.2f",y_rms);
	  textStatInt->AddText(text);
	  sprintf(text,"#Deltap_{y} = %5.2f mc",py_rms);
	  textStatInt->AddText(text);
	  sprintf(text,"k_{p}#varepsilon_{y} = %5.2f",emity);
	  textStatInt->AddText(text);
	  sprintf(text,"k_{p}#beta_{y} = %5.2f",betay);
	  textStatInt->AddText(text);

	}

	textStatInt->Draw();

	TBox *lFrame = new TBox(gPad->GetUxmin(), gPad->GetUymin(),
				gPad->GetUxmax(), gPad->GetUymax());
	lFrame->SetFillStyle(0);
	lFrame->SetLineColor(PGlobals::frameColor);
	lFrame->SetLineWidth(PGlobals::frameWidth);
	lFrame->Draw();
	
	gPad->RedrawAxis(); 
    
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
	
	TString fOutName2 = Form("./%s/Plots/Bunch/%s/Bunch-%s-%s-slp2x2_%i",simo.Data(),pData->GetRawSpeciesName(index).c_str(),pData->GetRawSpeciesName(index).c_str(),simo.Data(),time);
	
	CA4->Print(fOutName2 + ".ps[","Portrait");

	Int_t fontsize = 18;
	Int_t tfontsize = 20;
	Double_t txoffset = 4.0;
	Double_t lxoffset = 0.02;
	Double_t tyoffset = 2.2;
	Double_t lyoffset = 0.01;
	Double_t tzoffset = 3.0;
	Double_t lzoffset = 0.01;
	Double_t tylength = 0.015;
	Double_t txlength = 0.030;

	hP2X2->GetXaxis()->SetLabelFont(fonttype);
	hP2X2->GetXaxis()->SetLabelSize(fontsize);
	hP2X2->GetXaxis()->SetTitleFont(fonttype);
	hP2X2->GetXaxis()->SetTitleSize(tfontsize);
	hP2X2->GetXaxis()->SetTitleOffset(txoffset);
	hP2X2->GetXaxis()->SetTickLength(txlength);
	hP2X2->GetXaxis()->CenterTitle();
	
	hP2X2->GetYaxis()->SetLabelFont(fonttype);
	hP2X2->GetYaxis()->SetLabelSize(fontsize);
	hP2X2->GetYaxis()->SetTitleFont(fonttype);
	hP2X2->GetYaxis()->SetTitleSize(tfontsize);
	hP2X2->GetYaxis()->SetTitleOffset(tyoffset);
	hP2X2->GetYaxis()->SetTickLength(tylength);
	hP2X2->GetYaxis()->CenterTitle();
	
	hP2X2->GetZaxis()->SetLabelFont(fonttype);
	hP2X2->GetZaxis()->SetLabelSize(fontsize);
	hP2X2->GetZaxis()->SetTitleFont(fonttype);
	hP2X2->GetZaxis()->SetTitleSize(tfontsize);
	hP2X2->GetZaxis()->SetTitleOffset(tzoffset);
	hP2X2->GetZaxis()->SetTickLength(tylength);
	hP2X2->GetZaxis()->CenterTitle();
	  
	CA4->cd(1);

	gPad->SetTopMargin(0.02);
	gPad->SetBottomMargin(0.25);
	
	if(opt.Contains("logz")) {
	  gPad->SetLogz(1);
	} else {
	  gPad->SetLogz(0);
	}
	
	hP2X2->Draw("colz");

	gPad->Update();
	TPaletteAxis *palette = NULL;
	palette = (TPaletteAxis*)hP2X2->GetListOfFunctions()->FindObject("palette");
	if(palette) {
	  Double_t y1 = gPad->GetBottomMargin();
	  Double_t y2 = 1 - gPad->GetTopMargin();
	  Double_t x1 = 1 - gPad->GetRightMargin();
	
	  Double_t x1b = x1 + 0.01;
	  Double_t x2b = x1 + 0.035;
	  Double_t y1b = y1 + 0.03;
	  Double_t y2b = y2 - 0.03;
	  palette->SetX1NDC(x1b);
	  palette->SetY1NDC(y1b);
	  palette->SetX2NDC(x2b);
	  palette->SetY2NDC(y2b);
	  palette->SetBorderSize(2);
	  palette->SetLineColor(1);

	  TPave *pFrame = new TPave(x1b,y1b,x2b,y2b,1,"NDCL");
	  pFrame->SetFillStyle(0);
	  pFrame->SetLineColor(kBlack);
	  pFrame->SetShadowColor(0);
	  pFrame->Draw();
	}

	
	
	y1 = gPad->GetBottomMargin();
	y2 = 1 - gPad->GetTopMargin();
	x1 = gPad->GetLeftMargin();
	x2 = 1 - gPad->GetRightMargin();
	
	TPaveText *textStatInt = new TPaveText(x1+0.02,y2-0.40,x1+0.20,y2-0.05,"NDC");
	PGlobals::SetPaveTextStyle(textStatInt,12); 
	textStatInt->SetTextColor(kGray+3);
	textStatInt->SetTextFont(42);
	
	char text[128];
	sprintf(text,"#LTx#GT_{rms} = %5.2f %s",xrms,tspaSUnit.c_str());
	textStatInt->AddText(text);
	sprintf(text,"#LTp_{x}#GT_{rms} = %5.2f %s/c",yrms,teneSUnit.c_str());
	textStatInt->AddText(text);
	sprintf(text,"#varepsilon_{n} = %5.2f %s",emitx,emitSUnit.c_str());
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

	  hP2X2sl[k]->GetXaxis()->SetLabelFont(fonttype);
	  hP2X2sl[k]->GetXaxis()->SetLabelSize(fontsize);
	  hP2X2sl[k]->GetXaxis()->SetTitleFont(fonttype);
	  hP2X2sl[k]->GetXaxis()->SetTitleSize(tfontsize);
	  hP2X2sl[k]->GetXaxis()->SetTitleOffset(txoffset);
	  hP2X2sl[k]->GetXaxis()->SetTickLength(txlength);
	  hP2X2sl[k]->GetXaxis()->CenterTitle();
	  
	  hP2X2sl[k]->GetYaxis()->SetLabelFont(fonttype);
	  hP2X2sl[k]->GetYaxis()->SetLabelSize(fontsize);
	  hP2X2sl[k]->GetYaxis()->SetTitleFont(fonttype);
	  hP2X2sl[k]->GetYaxis()->SetTitleSize(tfontsize);
	  hP2X2sl[k]->GetYaxis()->SetTitleOffset(tyoffset);
	  hP2X2sl[k]->GetYaxis()->SetTickLength(tylength);
	  hP2X2sl[k]->GetYaxis()->CenterTitle();

	  hP2X2sl[k]->GetZaxis()->SetLabelFont(fonttype);
	  hP2X2sl[k]->GetZaxis()->SetLabelSize(fontsize);
	  hP2X2sl[k]->GetZaxis()->SetTitleFont(fonttype);
	  hP2X2sl[k]->GetZaxis()->SetTitleSize(tfontsize);
	  hP2X2sl[k]->GetZaxis()->SetTitleOffset(tzoffset);
	  hP2X2sl[k]->GetZaxis()->SetTickLength(tylength);
	  hP2X2sl[k]->GetZaxis()->CenterTitle();


	  if(opt.Contains("logz")) {
	    gPad->SetLogz(1);
	  } else {
	    gPad->SetLogz(0);
	  }

	  gPad->SetTopMargin(0.02);
	  gPad->SetBottomMargin(0.25);
	  
	  hP2X2sl[k]->Draw("colz");

	  gPad->Update();
	  palette = (TPaletteAxis*)hP2X2sl[k]->GetListOfFunctions()->FindObject("palette");
	  if(palette) {
	    Double_t y1 = gPad->GetBottomMargin();
	    Double_t y2 = 1 - gPad->GetTopMargin();
	    Double_t x1 = 1 - gPad->GetRightMargin();
	
	    Double_t x1b = x1 + 0.01;
	    Double_t x2b = x1 + 0.035;
	    Double_t y1b = y1 + 0.03;
	    Double_t y2b = y2 - 0.03;
	    palette->SetX1NDC(x1b);
	    palette->SetY1NDC(y1b);
	    palette->SetX2NDC(x2b);
	    palette->SetY2NDC(y2b);
	    palette->SetBorderSize(2);
	    palette->SetLineColor(1);

	    TPave *pFrame = new TPave(x1b,y1b,x2b,y2b,1,"NDCL");
	    pFrame->SetFillStyle(0);
	    pFrame->SetLineColor(kBlack);
	    pFrame->SetShadowColor(0);
	    pFrame->Draw();
	  }

	  
	  // Double_t y1 = gPad->GetBottomMargin();
	  Double_t y2 = 1 - gPad->GetTopMargin();
	  Double_t x1 = gPad->GetLeftMargin();
	  // Double_t x2 = 1 - gPad->GetRightMargin();
	  textStat[k] = new TPaveText(x1+0.02,y2-0.50,x1+0.20,y2-0.05,"NDC");
	  PGlobals::SetPaveTextStyle(textStat[k],12); 
	  textStat[k]->SetTextColor(kGray+3);
	  textStat[k]->SetTextFont(42);
	  
	  char text[64];
	  sprintf(text,"%5.2f %s < #zeta < %5.2f %s",sBinLim[k],spaSUnit.c_str(),sBinLim[k+1],spaSUnit.c_str());
	  textStat[k]->AddText(text);
	  sprintf(text,"#LTx#GT_{rms} = %5.2f %s",sx_rms[k],tspaSUnit.c_str());
	  textStat[k]->AddText(text);
	  sprintf(text,"#LTp_{x}#GT_{rms} = %5.2f %s",spx_rms[k],teneSUnit.c_str());
	  textStat[k]->AddText(text);
	  sprintf(text,"#varepsilon_{n,x} = %5.2f %s",semitx[k],emitSUnit.c_str());
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
      TString filename = Form("./%s/Plots/Bunch/%s/Bunch-%s-%s_%i.root",simo.Data(),pData->GetRawSpeciesName(index).c_str(),pData->GetRawSpeciesName(index).c_str(),simo.Data(),time);
      TFile *ofile = new TFile(filename,"RECREATE");

      hX1->SetLineWidth(2);
      hX1->SetLineColor(1);
      hX1->SetFillStyle(0);

      hP1->SetLineWidth(1);
      hP1->SetLineColor(1);
      hP1->SetFillStyle(0);

      hX1->Write("hX1",TObject::kOverwrite);
      hP1->Write("hP1",TObject::kOverwrite);
      hP1X1->Write("hP1X1",TObject::kOverwrite);
      hP1range->Write("hP1range",TObject::kOverwrite);
      hP1X1range->Write("hP1X1range",TObject::kOverwrite);
      hP2X2->Write("hP2X2",TObject::kOverwrite);
      hP2X2range->Write("hP2X2range",TObject::kOverwrite);
      hX2X1->Write("hX2X1",TObject::kOverwrite);
      
      gErms->Write("gErms",TObject::kOverwrite);
      gXrms->Write("gXrms",TObject::kOverwrite);
      gEmitx->Write("gEmitx",TObject::kOverwrite);

      if(pData->Is3D()) {
	gYrms->Write("gYrms",TObject::kOverwrite);
	gEmity->Write("gEmity",TObject::kOverwrite);

	hP3X3->Write("hP3X3",TObject::kOverwrite);
	hP3X3range->Write("hP3X3range",TObject::kOverwrite);
	hX3X1->Write("hX3X1",TObject::kOverwrite);
	hX3X2->Write("hX3X2",TObject::kOverwrite);
      }
	
      ofile->Close();
    }
    
    // Delete[] newly created vectors
    delete[] sBinLim;
    delete[] zbin;
    delete[] sEmean;
    delete[] sErms;

    delete[] sx_mean;
    delete[] sx_rms;
    delete[] spx_mean;
    delete[] spx_rms;
    delete[] semitx;
    delete[] sbetax;

    if(pData->Is3D()) {  
      delete[] sy_mean;
      delete[] sy_rms;
      delete[] spy_mean;
      delete[] spy_rms;
      delete[] semity;
      delete[] sbetay;
    }

    // end time looper
  }
}

void InitRange(const TString &sim, Double_t &x1Min, Double_t &x1Max, Int_t &SNbin, Double_t &x1BinMin, Double_t &x1BinMax) {

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
      x1Min = -5.2;
      x1Max = -3.2; 
       
      SNbin = 100;
      x1BinMin = -4.70;
      x1BinMax = -3.60;
    }
     
    if(sim.Contains("v2.5kA.G.ZH.DDR")) {
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
      
  } else if (sim.Contains("BOND_betatron")) {
      
    x1Min = -9.2;
    x1Max = -8.0; 
      
    SNbin = 80;
    x1BinMin = -8.90;
    x1BinMax = -8.55;
      
  } else if (sim.Contains("pitz")) {
    x1Min = -55.0;
    x1Max = 25.0; 
      
    SNbin = 100;
    x1BinMin = -32.0;
    x1BinMax =   2.0;

  }
}



