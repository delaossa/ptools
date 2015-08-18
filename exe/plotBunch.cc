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
    printf("      <-index(species index)>\n");
    printf("      <-i(initial time)> <-f(final time)> <-s(time step)>\n");
    printf("      <--png> <--pdf> <--eps>\n");
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

    if(arg.Contains("--pdf")) {
      opt += "pdf";
    } else if(arg.Contains("--eps")) {
      opt += "eps";
    } else if(arg.Contains("--png")) {
      opt += "png";
    } else if(arg.Contains("--tiff")) {
      opt += "tiff";
    } else if(arg.Contains("--comov")){
      opt += "comov";
    } else if(arg.Contains("--units")){
      opt += "units";
    } else if(arg.Contains("--center")){
      opt += "center";
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
    } else if( !(arg.Contains("pitz") || arg.Contains("flash") || arg.Contains("regae") || arg.Contains("Gauss") || arg.Contains("pwfa") || arg.Contains("vacuum") || arg.Contains("seed") || arg.Contains("dchirp") || arg.Contains("facet") || arg.Contains("FACET") || arg.Contains("rake") || arg.Contains("LWFA") || arg.Contains("BOND") ) ) {
      cout << Form("\t Invalid argument (%i): exiting...\n",l) << endl;
      return 0;
    } else {
      sim = arg;
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

  if(iStart<0) iStart = time;
  if(iEnd<=iStart) iEnd = iStart;
  
  // Some plasma constants
  Float_t n0 = pData->GetPlasmaDensity();
  Float_t kp = pData->GetPlasmaK();
  Float_t skindepth = 1.;
  if(kp!=0.0) skindepth = 1/kp;
  // Float_t E0 = pData->GetPlasmaE0();

  // z start of the plasma in normalized units.
  Float_t zStartPlasma = pData->GetPlasmaStart()*kp;
  // z start of the beam in normalized units.
  Float_t zStartBeam = pData->GetBeamStart()*kp;

  //  opt += "comovcenter";
  
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
    Float_t Time = pData->GetRealTime();

    // Centering time and z position:
    Float_t shiftz = pData->Shift(opt);

    if(opt.Contains("center")) {
      Time -= zStartPlasma;
      if(opt.Contains("comov"))  // Centers on the head of the beam.
	Time += zStartBeam;
    } 

    // Spatial resolution
    Float_t dx1, dx2, dx3;
    
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
    Float_t x1BinMin = -4.5;
    Float_t x1BinMax = -4.0;

    // Spatial coordinates intervals:
    // Float_t x1Min = -7.8;
    Float_t x1Min = -7.75;
    Float_t x1Max = -7.0;
    Float_t x2Min = -0.5;
    Float_t x2Max =  0.5;
    Float_t x3Min = -0.5;
    Float_t x3Max =  0.5;

    // Momentum coordinates intervals:
    Float_t p1Min =  Pmin;
    Float_t p1Max =  Pmax;
    Float_t p2Min = -15.0;
    Float_t p2Max =  15.0;
    Float_t p3Min = -15.0;
    Float_t p3Max =  15.0;

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

       x1Nbin = 400;
       p1Nbin = 200;
       x2Nbin = 200;
       p2Nbin = 200;
       
       x1Min = -5.2;
       x1Max = -3.2; 
       
       SNbin = 100;
       x1BinMin = -4.70;
       x1BinMax = -3.70;
      
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
      
    }

    // dummy shift
    Float_t dshiftz = pData->Shift("centercomov");
    x1Min += dshiftz;
    x1Max += dshiftz;
    x1BinMin += dshiftz;
    x1BinMax += dshiftz;

    // --------------------------------------------------
    
    cout << Form("\n 1. Reading file : ") << pData->GetRawFileName(index)->c_str() << endl;

    // Float_t **var = NULL;
    UInt_t Nvar = 8;
    if(!pData->Is3D()) Nvar = 7;
    
    Float_t **var;
    var = new Float_t*[Nvar];
    UInt_t Np = pData->GetRawArray(pData->GetRawFileName(index)->c_str(),var);

    // Scan histogram
    TH1F *hScanX1 = (TH1F*) gROOT->FindObject("hScanX1");
    if(hScanX1) delete hScanX1;
    hScanX1 = new TH1F("hScanX1","",pData->GetX1N(),pData->GetX1Min(),pData->GetX1Max());
    // hScanX1->GetYaxis()->SetTitle("I [kA]");
    // hScanX1->GetXaxis()->SetTitle("#zeta [#mum]");

    // -------------------------------------------------------------------------------


    // auto ranges factor
    Float_t rfactor = 0.3;
    
    if(opt.Contains("autop")) {

      Float_t MinP1 = 999999;
      Float_t MaxP1 = -999999;
      Float_t MinP2 = 999999;
      Float_t MaxP2 = -999999;
      Float_t MinP3 = 999999;
      Float_t MaxP3 = -999999;

      Float_t MinX2 = 999999;
      Float_t MaxX2 = -999999;
      Float_t MinX3 = 999999;
      Float_t MaxX3 = -999999;

      for(UInt_t i=0;i<Np;i++) {

	// if(var[5][i]-shiftz<x1BinMin || var[5][i]-shiftz>x1BinMax ) continue; 
	// if(var[5][i]-shiftz<x1Min || var[5][i]-shiftz>x1Max ) continue; 
	if(var[5][i]<x1Min || var[5][i]>x1Max ) continue; 

	if(var[1][i]<MinP1) MinP1 = var[1][i];
	if(var[1][i]>MaxP1) MaxP1 = var[1][i];
	if(var[2][i]<MinP2) MinP2 = var[2][i];
	if(var[2][i]>MaxP2) MaxP2 = var[2][i];
	if(var[3][i]<MinP3) MinP3 = var[3][i];
	if(var[3][i]>MaxP3) MaxP3 = var[3][i];
	if(var[6][i]<MinX2) MinX2 = var[6][i];
	if(var[6][i]>MaxX2) MaxX2 = var[6][i];
	if(Nvar==8) {
	  if(var[7][i]<MinX3) MinX3 = var[7][i];
	  if(var[7][i]>MaxX3) MaxX3 = var[7][i];
	}
      }
      
      p1Min = MinP1 - rfactor*(MaxP1-MinP1);
      p1Max = MaxP1 + rfactor*(MaxP1-MinP1);
      p2Min = MinP2 - rfactor*(MaxP2-MinP2);
      p2Max = MaxP2 + rfactor*(MaxP2-MinP2);
      p3Min = MinP3 - rfactor*(MaxP3-MinP3);
      p3Max = MaxP3 + rfactor*(MaxP3-MinP3);

      x2Min = MinX2 - rfactor*(MaxX2-MinX2);
      x2Max = MaxX2 + rfactor*(MaxX2-MinX2);
      if(Nvar==8) {
	x3Min = MinX3 - rfactor*(MaxX3-MinX3);
	x3Max = MaxX3 + rfactor*(MaxX3-MinX3);
      }

    } else if(opt.Contains("auto")) {

      Float_t MinP1 = 999999;
      Float_t MaxP1 = -999999;
      Float_t MinP2 = 999999;
      Float_t MaxP2 = -999999;
      Float_t MinP3 = 999999;
      Float_t MaxP3 = -999999;

      Float_t MinX1 = 999999;
      Float_t MaxX1 = -999999;
      Float_t MinX2 = 999999;
      Float_t MaxX2 = -999999;
      Float_t MinX3 = 999999;
      Float_t MaxX3 = -999999;

      for(UInt_t i=0;i<Np;i++) {
	if(var[1][i]<MinP1) MinP1 = var[1][i];
	if(var[1][i]>MaxP1) MaxP1 = var[1][i];
	if(var[2][i]<MinP2) MinP2 = var[2][i];
	if(var[2][i]>MaxP2) MaxP2 = var[2][i];
	if(var[3][i]<MinP3) MinP3 = var[3][i];
	if(var[3][i]>MaxP3) MaxP3 = var[3][i];
	if(var[5][i]<MinX1) MinX1 = var[5][i];
	if(var[5][i]>MaxX1) MaxX1 = var[5][i];
	if(var[6][i]<MinX2) MinX2 = var[6][i];
	if(var[6][i]>MaxX2) MaxX2 = var[6][i];
	if(Nvar==8) {
	  if(var[7][i]<MinX3) MinX3 = var[7][i];
	  if(var[7][i]>MaxX3) MaxX3 = var[7][i];
	}

	hScanX1->Fill(var[5][i],TMath::Abs(var[4][i]));
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
      if(Nvar==8) {
	x3Min = MinX3 - rfactor*(MaxX3-MinX3);
	x3Max = MaxX3 + rfactor*(MaxX3-MinX3);
      }

      x1BinMin = MinX1;
      x1BinMax = MaxX1;

      x1Min = floor((x1Min-pData->GetXMin(0))/dx1) * dx1 + pData->GetXMin(0);  
      x1Max = floor((x1Max-pData->GetXMin(0))/dx1) * dx1 + pData->GetXMin(0);  
      x1Nbin = ceil ((x1Max - x1Min)/(dx1));
      // p1Nbin = x1Nbin;
      
      x2Min = floor((x2Min-pData->GetXMin(1))/dx2) * dx2 + pData->GetXMin(1);  
      x2Max = floor((x2Max-pData->GetXMin(1))/dx2) * dx2 + pData->GetXMin(1);  
      x2Nbin = ceil ((x2Max - x2Min)/(dx2));
      // p2Nbin = x2Nbin;

      if(pData->Is3D()) {
	x3Min = floor((x3Min-pData->GetXMin(2))/dx3) * dx3 + pData->GetXMin(2);  
	x3Max = floor((x3Max-pData->GetXMin(2))/dx3) * dx3 + pData->GetXMin(2);  
	x3Nbin = ceil ((x3Max - x3Min)/(dx3));
      }
      // p3Nbin = x3Nbin;
      
      // SNbin = ceil ((x1BinMax - x1BinMin)/(2.0*dx1));

      // Set X1 limits using hScanX1
      Float_t peakFactor = 0.2;
      Float_t x1min = -999;
      Float_t x1max = -999;
      {
	Float_t maxValue = hScanX1->GetBinContent(hScanX1->GetMaximumBin());
	Int_t   lBin = -1;
	for(Int_t i=1;i<=hScanX1->GetNbinsX();i++) {
	  Float_t binValue = hScanX1->GetBinContent(i);
	  if(binValue>maxValue*peakFactor) {
	    lBin = i;
	    x1min = hScanX1->GetBinCenter(i);
	    break;
	  }
	}
	
	Int_t rBin = -1;
	
	for(Int_t i=hScanX1->GetNbinsX();i>0;i--) {
	  Float_t binValue = hScanX1->GetBinContent(i);
	  if(binValue>maxValue*peakFactor) {
	    rBin = i;
	    x1max = hScanX1->GetBinCenter(i);
	    break;
	  }
	}
      }

      x1BinMin = x1min;
      x1BinMax = x1max;

      // x1Min = x1min - 2*rfactor*(x1max-x1min);
      // x1Max = x1max + 2*rfactor*(x1max-x1min);
    }
    

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
    hX1->GetYaxis()->SetTitle("I [kA]");
    hX1->GetXaxis()->SetTitle("#zeta [#mum]");
    
    sprintf(hName,"hP1");
    TH1F *hP1 = (TH1F*) gROOT->FindObject(hName);
    if(hP1) delete hP1;
    hP1 = new TH1F(hName,"",p1Nbin,p1Min,p1Max);
    hP1->GetYaxis()->SetTitle("p_{z} [GeV/c]");
    if(opt.Contains("comov"))
      hP1->GetXaxis()->SetTitle("#zeta [#mum]");
    else
      hP1->GetXaxis()->SetTitle("z [#mum]");
    
    sprintf(hName,"hP1X1");
    TH2F *hP1X1 = (TH2F*) gROOT->FindObject(hName);
    if(hP1X1) delete hP1X1;
    hP1X1 = new TH2F(hName,"",x1Nbin,x1Min,x1Max,p1Nbin,p1Min,p1Max);
    if(opt.Contains("comov"))
      hP1X1->GetXaxis()->SetTitle("#zeta [#mum]");
    else
      hP1X1->GetXaxis()->SetTitle("z [#mum]");
    hP1X1->GetYaxis()->SetTitle("p_{z} [GeV/c]");
    hP1X1->GetZaxis()->SetTitle("Charge [pC]");
    hP1X1->GetZaxis()->CenterTitle();
  
    sprintf(hName,"hP2X2");
    TH2F *hP2X2 =  (TH2F*) gROOT->FindObject(hName);
    if(hP2X2) delete hP2X2;
    hP2X2 = new TH2F(hName,"",x2Nbin,x2Min,x2Max,p2Nbin,p2Min,p2Max);
    hP2X2->GetXaxis()->SetTitle("x [#mum]");
    hP2X2->GetYaxis()->SetTitle("p_{x}/mc");
    hP2X2->GetZaxis()->SetTitle("Charge [pC]");
    
    // Sliced quantities:
    // --------------------------------------------------------------------------
 
    // Set the binning
    Float_t *sBinLim = new Float_t[SNbin+1];
    sBinLim[0] = x1BinMin;
    sBinLim[SNbin] = x1BinMax;
    Float_t slbinSize = (sBinLim[SNbin] - sBinLim[0])/SNbin;
    for(Int_t i=1;i<SNbin;i++) {
      sBinLim[i] = sBinLim[i-1] + slbinSize;
    }

    TH1F **hP1sl = new TH1F*[SNbin];
    TH2F **hP2X2sl = new TH2F*[SNbin];
    for(Int_t k=0;k<SNbin;k++) {
      sprintf(hName,"hP2X2sl_%2i",k);
      hP2X2sl[k] = (TH2F*) gROOT->FindObject(hName);
      if(hP2X2sl[k]) delete hP2X2sl[k];
      hP2X2sl[k] = new TH2F(hName,"",x2Nbin,x2Min,x2Max,p2Nbin,p2Min,p2Max);

      hP2X2sl[k]->GetXaxis()->SetTitle("k_{p} x");
      hP2X2sl[k]->GetYaxis()->SetTitle("p_{x}/mc");


      sprintf(hName,"hP1sl_%2i",k);
      hP1sl[k] = (TH1F*) gROOT->FindObject(hName);
      if(hP1sl[k]) delete hP1sl[k];
      hP1sl[k] = new TH1F(hName,"",p1Nbin,p1Min,p1Max);

    }

    // Filling histos
    cout << Form("\n 2. Filling histograms from file ... ") ;
    
    for(UInt_t i=0;i<Np;i++) {

      var[5][i] -= shiftz;

      if(var[5][i]<x1Min || var[5][i]>x1Max ) continue; 
      if(var[6][i]<x2Min || var[6][i]>x2Max ) continue; 
      if(Nvar==8) if(var[7][i]<x3Min || var[7][i]>x3Max ) continue; 
      if(var[1][i]<p1Min || var[1][i]>p1Max ) continue; 
      if(var[2][i]<p2Min || var[2][i]>p2Max ) continue; 
      if(var[3][i]<p3Min || var[3][i]>p3Max ) continue; 


      hX1->Fill(var[5][i],TMath::Abs(var[4][i]));
      hP1->Fill(var[1][i],TMath::Abs(var[4][i]));

      hP1X1->Fill(var[5][i],var[1][i],TMath::Abs(var[4][i]));

      // Slices    
      if(var[5][i]<sBinLim[0] || var[5][i]>sBinLim[SNbin]) continue;
      Int_t iBin = -1;
      for(Int_t j=0; j<SNbin; j++) {

	if(var[5][i]<sBinLim[j+1]) {
	  iBin = j;
	  break;
	}
      }
      if(iBin<0) continue;

      // Projected emittance in the bunch range. (skip the tails to "match" the sliced ones)
      hP2X2->Fill(var[6][i],var[2][i],TMath::Abs(var[4][i]));

      hP1sl[iBin]->Fill(var[1][i],TMath::Abs(var[4][i]));
      //      cout << Form("  iBin = %i   p1 = %7.3f",iBin,var[1][i]) << endl;
      hP2X2sl[iBin]->Fill(var[6][i],var[2][i],TMath::Abs(var[4][i]));
      
    }
    cout << " done! " << endl;

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


    // Total relative energy spread within FWHM:
    sprintf(hName,"hP1cut");
    TH1F *hP1cut = (TH1F*) hP1->Clone(hName);
    hP1cut->Reset();

    Float_t maxValue = hP1->GetBinContent(hP1->GetMaximumBin());
    Int_t   lBin = -1;
    Float_t pzmin = -999;
    for(Int_t i=1;i<=hP1->GetNbinsX();i++) {
      Float_t binValue = hP1->GetBinContent(i);
      if(binValue>maxValue/2) {
	lBin = i;
	pzmin = hP1->GetBinCenter(i);
	break;
      }
    }

    Int_t rBin = -1;
    Float_t pzmax = -999;
    for(Int_t i=hP1->GetNbinsX();i>0;i--) {
      Float_t binValue = hP1->GetBinContent(i);
      if(binValue>maxValue/2) {
	rBin = i;
	pzmax = hP1->GetBinCenter(i);
	break;
      }
    }
    
    for(Int_t i=lBin;i<=rBin;i++) {
      Float_t binValue = hP1->GetBinContent(i);
      hP1cut->SetBinContent(i,binValue);
    }
    //  hP1cut->ResetStats();
    Float_t pzmeanFWHM = hP1cut->GetMean();
    Float_t pzrmsFWHM = hP1cut->GetRMS();
        
    // cout << Form("  zMean = %7.3f   pzMean = %7.3f",zmean,pzmean) << endl;
    // cout << Form("  zRms  = %7.3f   pzRms  = %7.3f",zrms,pzrms) << endl;
    // cout << Form("  zEmittance = %7.3f",emitz) << endl;
    
    cout <<  Form("  Normal  : pzMean  = %7.3f   pzRms  = %7.3f",pzmean,pzrms) << endl;
    cout <<  Form("  In FWHM : pzMean  = %7.3f   pzRms  = %7.3f -> Interval: pmin = %7.3f pmax = %7.3f  ",pzmeanFWHM,pzrmsFWHM,pzmin,pzmax) << endl;

    
    // ----------------------------------------


    // Transverse phasespace (Total)

    // P2X2 Phasespace
    hP2X2->GetStats(stats);

    xmean  = stats[2]/stats[0];
    xrms2  = stats[3]/stats[0] - xmean * xmean;
    xrms   = (xrms2>0.0) ? TMath::Sqrt(xrms2) : 0.0 ;
    ymean  = stats[4]/stats[0];
    yrms2  = stats[5]/stats[0] - ymean * ymean;
    yrms   = (yrms2>0.0) ? TMath::Sqrt(yrms2) : 0.0 ;
    xyrms2 = stats[6]/stats[0] - stats[2]*stats[4]/(stats[0]*stats[0]);
    xyrms  = (xyrms2>0.0) ? TMath::Sqrt(xyrms2) : 0.0 ;
    emit2  = xrms2*yrms2 - xyrms2*xyrms2;
    emit   = (emit2>0.0) ? TMath::Sqrt(emit2) : 0.0 ;

    Float_t emitx = emit;
    Float_t x_mean = xmean;
    Float_t x_rms = xrms;
    Float_t pxmean = ymean;
    Float_t pxrms = yrms;


    // Charge  
    hX1->Scale(dx1*dx2*dx3);
    Float_t Charge = hX1->Integral();

    cout << Form("\n 4. Calculating sliced quantities.. ") << endl ;

    TGraph *gEmitx = NULL;
    TGraph *gXrms = NULL;
    TGraph *gErms = NULL;

    Float_t * sxmean = new Float_t[SNbin];
    Float_t * sxrms2 = new Float_t[SNbin]; 
    Float_t * sxrms = new Float_t[SNbin]; 
    Float_t * symean = new Float_t[SNbin];
    Float_t * syrms2 = new Float_t[SNbin]; 
    Float_t * syrms = new Float_t[SNbin];  
    Float_t * sxyrms2 = new Float_t[SNbin];
    Float_t * sxyrms = new Float_t[SNbin];
    Float_t * semit = new Float_t[SNbin];
    Float_t * semit2 = new Float_t[SNbin];

    Float_t * xbin = new Float_t[SNbin];

    Float_t * sEmean = new Float_t[SNbin];
    Float_t * sErms = new Float_t[SNbin];
    Float_t * sErms2 = new Float_t[SNbin];

    for(Int_t k=0;k<SNbin;k++) {
      sxmean[k] = symean[k] = sxrms2[k] = syrms2[k] = sxrms[k] = syrms[k]
	= sxyrms2[k] = sxyrms[k] = xbin[k] = semit[k] = 0.0;
      sEmean[k] = sErms[k] = sErms2[k] = 0.0;

      hP2X2sl[k]->GetStats(stats);
      
      sxmean[k]  = stats[2]/stats[0];
      sxrms2[k]  = stats[3]/stats[0] - sxmean[k] * sxmean[k];
      sxrms[k]   = (sxrms2[k]>0.0) ? TMath::Sqrt(sxrms2[k]) : 0.0 ;
      symean[k]  = stats[4]/stats[0];
      syrms2[k]  = stats[5]/stats[0] - symean[k] * symean[k];
      syrms[k]   = (syrms2[k]>0.0) ? TMath::Sqrt(syrms2[k]) : 0.0 ;
      sxyrms2[k] = stats[6]/stats[0] - stats[2]*stats[4]/(stats[0]*stats[0]);
      sxyrms[k]  = (sxyrms2[k]>0.0) ? TMath::Sqrt(sxyrms2[k]) : 0.0 ;
      semit2[k] = sxrms2[k]*syrms2[k] - sxyrms2[k]*sxyrms2[k];
      semit[k] = (semit2[k]>0.0) ? TMath::Sqrt(semit2[k]) : 0.0 ;

      xbin[k] = (sBinLim[k] + sBinLim[k+1])/2.;

      // Double_t stats1d[4];  // { sumw, sumw2, sumwx, sumwx2 }
      // hP1sl[k]->GetStats(stats1d);

      // sEmean[k] = stats[2]/stats[0];
      // sErms2[k] = stats[3]/stats[0] - sEmean[k] * sEmean[k];
      // sErms[k]  = (sErms2[k]>0.0) ? TMath::Sqrt(sErms2[k]) : 0.0 ;

      sEmean[k] = hP1sl[k]->GetMean();
      sErms[k]  = hP1sl[k]->GetRMS();
      
      // cout<< Form("\nk = %i : (x1 > %f && x1 < %f)",k,sBinLim[k],sBinLim[k+1]) << endl; 

      // cout << Form("  xMean = %7.3f   yMean = %7.3f",sxmean[k],symean[k]) << endl;
      // cout << Form("  xRms  = %7.3f   yRms  = %7.3f",sxrms[k],syrms[k]) << endl;
      // cout << Form("  Emittance = %7.3f",semit[k]) << endl;

      // cout << Form("  Emean = %7.3f   Erms = %7.3f",sEmean[k],sErms[k]) << endl;


    }
    cout << " done! " << endl;


    // Changing to user units: 
    // --------------------------

    string curUnit;
    string emitUnit;
    string ermsUnit;
    string auxUnit;
    
    Float_t  lightspeed =  PConst::c_light / (PUnits::um/PUnits::femtosecond);
    
    if(opt.Contains("units") && n0) {

      Time *= skindepth / PUnits::um;

      x1Min *= skindepth / PUnits::um;
      x1Max *= skindepth / PUnits::um;
      p1Min *= PConst::ElectronMassE / PUnits::GeV;
      p1Max *= PConst::ElectronMassE / PUnits::GeV;
      
      x1BinMin *= skindepth / PUnits::um;
      x1BinMax *= skindepth / PUnits::um;
      for(Int_t k=0;k<=SNbin;k++) 
	sBinLim[k] *= skindepth / PUnits::um;
      
      hP1X1->SetBins(x1Nbin,x1Min,x1Max,p1Nbin,p1Min,p1Max);

      // Converting electron density
      Float_t dV = skindepth * skindepth * skindepth;
      Charge *= n0 * dV * (PConst::ElectronCharge/PUnits::picocoulomb);

      for(Int_t j=0;j<hP1X1->GetNbinsX();j++) {
	for(Int_t k=0;k<hP1X1->GetNbinsY();k++) {
	  Float_t binValue =  fabs(hP1X1->GetBinContent(j,k) * dx1 * dx2 * dx3 * dV * n0 *
				   (PConst::ElectronCharge/PUnits::picocoulomb));
	  //cout << Form(" value = %f",binValue) << endl;
	  hP1X1->SetBinContent(j,k,binValue);

	}
      }

      hP1->SetBins(p1Nbin,p1Min,p1Max);
      hP1->GetYaxis()->SetTitle("p_{z} [GeV/c]");

      hX1->SetBins(x1Nbin,x1Min,x1Max);
      Float_t binSize = (x1Max - x1Min)/x1Nbin;

      hX1->Scale(TMath::Abs(n0 * dV * (PConst::ElectronCharge/PUnits::picocoulomb) * (lightspeed/binSize)));


      Float_t curscal = rawf;
      curUnit = "kA"; 

      //PUnits::UnitsTable::Get();
      stringstream stream;      
      // Float_t peakCurr = hX1->GetMaximum() * PUnits::kA;
      // stream << PUnits::BestUnit(peakCurr,"Current") << endl;
      // Float_t best;
      // stream >> best >> curUnit;
      // curscal = best/(peakCurr/PUnits::kA);      
      // cout << Form(" Peak current = %.2f %s  --> scale = %f",best,curUnit.c_str(),curscal) << endl;
      
      hX1->Scale(curscal);

      // hX1->GetYaxis()->SetTitle("I[kA]");
      hX1->GetYaxis()->SetTitle("");
      if(opt.Contains("comov"))
	hX1->GetXaxis()->SetTitle("#zeta [#mum]");
      else
	hX1->GetXaxis()->SetTitle("z [#mum]");
      
      
      zmean *= skindepth / PUnits::um;
      zrms  *= skindepth / PUnits::um;

      x_mean *= skindepth / PUnits::um;
      x_rms  *= skindepth / PUnits::um;

      xmean *= skindepth / PUnits::um;
      xrms  *= skindepth / PUnits::um;

      pzmean *= PConst::ElectronMassE / PUnits::GeV;
      pzrms  *= PConst::ElectronMassE / PUnits::GeV;
      pzmeanFWHM *= PConst::ElectronMassE / PUnits::GeV;
      pzrmsFWHM *= PConst::ElectronMassE / PUnits::GeV;
      pzmin  *= PConst::ElectronMassE / PUnits::GeV;
      pzmax  *= PConst::ElectronMassE / PUnits::GeV;

      pxmean *= PConst::ElectronMassE / PUnits::MeV;
      pxrms  *= PConst::ElectronMassE / PUnits::MeV;

      emitz *= (skindepth / PUnits::um);

      // Transverse phase-space
      x2Min *= skindepth/PUnits::um;
      x2Max *= skindepth/PUnits::um;
      p2Min *= PConst::ElectronMassE / PUnits::MeV;
      p2Max *= PConst::ElectronMassE / PUnits::MeV;
      
      hP2X2->SetBins(x2Nbin,x2Min,x2Max,p2Nbin,p2Min,p2Max);
      for(Int_t j=0;j<hP2X2->GetNbinsX();j++) {
	for(Int_t k=0;k<hP2X2->GetNbinsY();k++) {
	  Float_t binValue =  fabs(hP2X2->GetBinContent(j,k) * dx1 * dx2 * dx3 * dV * n0 *
				   (PConst::ElectronCharge/PUnits::picocoulomb));
	  //cout << Form(" value = %f",binValue) << endl;
	  hP2X2->SetBinContent(j,k,binValue);

	}
      }

      hP2X2->GetXaxis()->SetTitle("x [#mum]");
      hP2X2->GetYaxis()->SetTitle("p_{x} [MeV/c]");
      hP2X2->GetZaxis()->SetTitle("Charge [pC]");
      hP2X2->GetZaxis()->CenterTitle();

     

      //      emitx *= (skindepth / PUnits::um);
      emitx *= skindepth;

      Float_t emiscal = 1;
      emitUnit = "#mum";
      // stream << PUnits::BestUnit(emitx,"Length") << endl;
      // stream >> best >> emitUnit >> auxUnit;
      // if(!auxUnit.empty()) {
      // 	emitUnit += " ";
      // 	emitUnit += auxUnit;
      // }
      // emiscal = best/(emitx/PUnits::um);
      // // cout << Form(" Emittance = %.2f %s  --> scale = %f",best,emitUnit.c_str(),emiscal) << endl;
      emitx *= emiscal/PUnits::um;

      Float_t erelMax = -999.;
      for(Int_t k=0;k<SNbin;k++) {

	hP2X2sl[k]->SetBins(x2Nbin,x2Min,x2Max,p2Nbin,p2Min,p2Max);
	// for(Int_t j=0;j<hP2X2sl[k]->GetNbinsX();j++) {
	// 	for(Int_t l=0;l<hP2X2sl[k]->GetNbinsY();l++) {
	// 	  Float_t binValue =  fabs(hP2X2sl[k]->GetBinContent(j,l) * dx1 * dx2 * dx3 * dV * n0 *
	// 				    (PConst::ElectronCharge/PUnits::picocoulomb));
	// 	  //cout << Form(" value = %f",binValue) << endl;
	// 	  hP2X2sl[k]->SetBinContent(j,l,binValue);

	// 	}
	// }

	hP2X2sl[k]->GetZaxis()->SetTitle("Charge [a.u.]");
	hP2X2sl[k]->GetXaxis()->SetTitle("x [#mum]");
	hP2X2sl[k]->GetYaxis()->SetTitle("p_{x} [MeV/c]");

	hP2X2sl[k]->GetZaxis()->CenterTitle();

	xbin[k] *= skindepth / PUnits::um;

	sxmean[k] *= skindepth / PUnits::um;
	sxrms[k]  *= skindepth / PUnits::um;
	symean[k] *= PConst::ElectronMassE / PUnits::MeV;
	syrms[k] *= PConst::ElectronMassE / PUnits::MeV;

	semit[k] *= emiscal * (skindepth / (PUnits::um));

	sErms[k]  *= 100 / sEmean[k];
	sEmean[k] *= PConst::ElectronMassE / PUnits::GeV;

	//	cout << Form(" Energy spread = %.2f +/- %.2f  %%",sEmean[k],sErms[k]) << endl;  
	if(sErms[k]>erelMax) erelMax = sErms[k];
	
      }

      Float_t erelrmsscal = 1;
      ermsUnit = "%";
      // Float_t erelrmsscal = 1;
      // ermsUnit = "%";
      // stream << PUnits::BestUnit(erelMax,"Percentage") << endl;
      // stream >> best >> ermsUnit;
      // Float_t erelrmsscal = best/erelMax;
      // cout << Form(" Energy spread = %.2f %s  --> scale = %f",best,ermsUnit.c_str(),erelrmsscal) << endl;      
      // erelMax *= erelrmsscal/PUnits::perCent;      
      for(Int_t k=0;k<SNbin;k++)
       	sErms[k]  *= erelrmsscal;
      
    }
    // End of the users units module    

    cout << "\n  Summary _______________________________________________________ " << endl;
    if(opt.Contains("units")) {
      cout << Form("  Integrated charge (RAW) of specie %3i = %8f pC",index,Charge) << endl;
      cout << Form("  Peak current = %6.2f kA, Rms current = %6.2f kA",hX1->GetMaximum(),lightspeed*Charge/(TMath::Sqrt(2*PConst::pi)*zrms)) << endl;
      cout << Form("  Total energy = %6.2f GeV, rms = %3.1f %%",pzmean,100*pzrms/pzmean) << endl;
      cout << Form("  Trans. emit. = %6.2f  um",emitx) << endl;
      cout << Form("  Bunch length = %6.2f  um (rms), width = %6.2f  um (rms)",zrms,x_rms) << endl;
    } else {
      cout << Form("  Integrated charge (RAW) of specie %3i = %8.4f n0 * kp^-3",index,Charge) << endl;
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
      gPxmeanvsTime->SetPoint(nPoints,Time,pxmean);
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
      gPxrmsvsTime->SetPoint(nPoints,Time,pxrms);
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
      for(Int_t k=0;k<SNbin;k++) {
	xbin[k] -= zmean;
      }
      zmean = 0.0;
    }
    // ------

    // Create the graph with the sliced quantities:
    gEmitx = new TGraph(SNbin,xbin,semit);
    gXrms = new TGraph(SNbin,xbin,sxrms);
    gErms = new TGraph(SNbin,xbin,sErms);


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
    //   Float_t *x1bins = new Float_t[NP1X1Bins];
    //   Float_t *eRms   = new Float_t[NP1X1Bins];
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
      Float_t *yarray   = new Float_t[p1Nbin];
      Float_t *xarray   = new Float_t[p1Nbin];

      // This is for the right side:
      // Float_t xMax = x1Min + (x1Max-x1Min) * 0.9;
      // Float_t xMin = x1Max;
      // And this for left:
      Float_t xMin = hX1->GetXaxis()->GetXmin();
      Float_t xMax = hX1->GetXaxis()->GetXmin() + 
	(hX1->GetXaxis()->GetXmax()-hX1->GetXaxis()->GetXmin()) * 0.2;
      Float_t EneMax = hP1->GetMaximum();
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
      Float_t yMin =  999.9;
      Float_t yMax =  -999.9;

      for(Int_t k=0;k<SNbin;k++) {
	if(semit[k]<yMin)
	  yMin = semit[k];

	if(semit[k]>yMax)
	  yMax = semit[k];

	if(sErms[k]<yMin)
	  yMin = sErms[k];

	if(sErms[k]>yMax)
	  yMax = sErms[k];
      }

      for(Int_t k=1;k<=x1Nbin;k++) {
	Float_t value = hX1->GetBinContent(k);

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
      PPalette * pPalette = (PPalette*) gROOT->FindObject("electron");
      pPalette->cd();

      // Float_t Max  = hP1X1->GetMaximum();
      // Float_t Min  = hP1X1->GetMinimum();

      // hP1X1->GetZaxis()->SetRangeUser(Min,Max); 

      // Text objects
      TPaveText *textTime =  new TPaveText(0.55,0.76,0.80,0.86,"NDC");
      PGlobals::SetPaveTextStyle(textTime,32); 
      textTime->SetTextColor(kGray+2);
      char ctext[128];
      if(opt.Contains("units") && pData->GetPlasmaDensity()) 
	// sprintf(ctext,"z = %5.1f um", Time);
	sprintf(ctext,"z = %5.1f mm", Time/1000);
      else
	sprintf(ctext,"t = %5.1f #omega_{p}^{-1}",Time);
      textTime->AddText(ctext);
 
      TPaveText *textDen = new TPaveText(0.15,0.85,0.48,0.9,"NDC");
      PGlobals::SetPaveTextStyle(textDen,12); 
      textDen->SetTextColor(kOrange+10);
      if(opt.Contains("units") && pData->GetPlasmaDensity())
	sprintf(ctext,"n_{0} = %5.2f x 10^{17} / cc", n0 / (1e17/PUnits::cm3));
      else if(pData->GetBeamDensity() && pData->GetPlasmaDensity())
	sprintf(ctext,"n_{b}/n_{0} = %5.2f", pData->GetBeamDensity()/n0);
      textDen->AddText(ctext);

      TPaveText *textCharge = new TPaveText(0.15,0.25,0.48,0.3,"NDC");
      PGlobals::SetPaveTextStyle(textCharge,12); 
      textCharge->SetTextColor(kGray+2);
      if(opt.Contains("units") && pData->GetPlasmaDensity())
	sprintf(ctext,"Q = %5.2f pC", Charge);
      else
	sprintf(ctext,"Q = %5.2f n0#timeskp^{-3}", Charge);    
      textCharge->AddText(ctext);

      TPaveText *textMom = new TPaveText(0.55,0.035,0.80,0.14,"NDC");
      PGlobals::SetPaveTextStyle(textMom,32); 
      textMom->SetTextColor(kGray+3);
      textMom->SetTextFont(62);
      if(opt.Contains("units") && pData->GetPlasmaDensity())
	sprintf(ctext,"#LTp_{z}#GT = %5.2f GeV/c", pzmeanFWHM);
      else
	sprintf(ctext,"#LTp_{z}#GT = %5.2f mc", pzmeanFWHM);    
      textMom->AddText(ctext);


      TPaveText *textInfo = new TPaveText(0.55,0.4,0.80,0.75,"NDC");
      PGlobals::SetPaveTextStyle(textInfo,32); 
      textInfo->SetTextColor(kGray+2);
      textInfo->SetTextFont(42);
      sprintf(ctext,"Q = %5.2f pC",Charge);
      textInfo->AddText(ctext);
      sprintf(ctext,"#Delta#zeta = %5.2f #mum",zrms);
      textInfo->AddText(ctext);
      //     sprintf(ctext,"#LTp_{z}#GT_{rms} = %5.2f GeV/c",pzrms);
      // sprintf(ctext,"#Deltap_{z}/p_{z} = %4.1f %%",100*pzrms/pzmean);
      // sprintf(ctext,"#Deltap_{z}/p_{z} = %4.1f %%",100*pzrmsFWHM/pzmeanFWHM);
      sprintf(ctext,"#Delta#gamma/#LT#gamma#GT = %4.1f %%",100*pzrmsFWHM/pzmeanFWHM);
      textInfo->AddText(ctext);
      sprintf(ctext,"#varepsilon_{n} = %5.2f %s",emitx,emitUnit.c_str());
      //sprintf(ctext,"#varepsilon_{n} = %5.2f %s",emitx,"#mum");
      textInfo->AddText(ctext);

      // Setup Pad layout: 
      const Int_t NPad = 2;
      TPad *pad[NPad];
      TH1F *hFrame[NPad];
      TString sLabels[] = {"(b)","(a)"};
      TPaveText **textLabel = new TPaveText*[NPad];

      Float_t lMargin = 0.15;
      Float_t rMargin = 0.18;
      Float_t bMargin = 0.15;
      Float_t tMargin = 0.04;
      Float_t factor = 1.0;    
      PGlobals::CanvasAsymPartition(C,NPad,lMargin,rMargin,bMargin,tMargin,factor,0.00);

      // Define the frames for plotting
      Int_t fonttype = 43;
      Int_t fontsize = 24;
      Int_t tfontsize = 28;
      Float_t txoffset = 2.0;
      Float_t lxoffset = 0.02;
      Float_t tyoffset = 1.3;
      Float_t lyoffset = 0.01;
      Float_t tylength = 0.02;
      Float_t txlength = 0.04;
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
      gP1left->Draw("F");
      gP1left->Draw("L");

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
      // Float_t dmax = hP1X1->GetMaximum();
      // Float_t dmin = 0.0;
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
	Float_t y1 = gPad->GetBottomMargin();
	Float_t y2 = 1 - gPad->GetTopMargin();
	Float_t x1 = 1 - gPad->GetRightMargin();
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


      Float_t y1 = gPad->GetBottomMargin();
      Float_t y2 = 1 - gPad->GetTopMargin();
      Float_t x1 = gPad->GetLeftMargin();
      Float_t x2 = 1 - gPad->GetRightMargin();
      Float_t yrange = y2-y1; 
      Float_t xrange = x2-x1; 

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
      Leg=new TLegend(0.50,0.75,1 - gPad->GetRightMargin() - 0.02,0.95);

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
      sprintf(sleg,"Current [%s]",curUnit.c_str());
      Leg->AddEntry(hX1  ,sleg,"L");	
      //sprintf(sleg,"Energy spread [%s]",ermsUnit.c_str());
      sprintf(sleg,"Energy spread [%s]",ermsUnit.c_str());
      Leg->AddEntry(gErms,sleg,"PL");
      //      sprintf(sleg,"Emittance [%s]",emitUnit.c_str());
      sprintf(sleg,"Emittance [%s]",emitUnit.c_str());
      Leg->AddEntry(gEmitx,sleg,"PL");
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

      // hP2X2sl[0]->Draw("colz");
      gEmitx->SetMarkerStyle(20);
      gEmitx->SetMarkerColor(kGray+2);
      gEmitx->SetMarkerSize(markerSize);
      gEmitx->SetLineWidth(lineWidth);
      gEmitx->SetLineColor(kGray+2);
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

      // Print to a file
      // Output file
      TString fOutName = Form("./%s/Plots/Bunch/%s/Bunch",sim.Data(),pData->GetSpeciesName(index).c_str());
      fOutName += Form("-%s_%i",sim.Data(),time);

      PGlobals::imgconv(C,fOutName,opt);
      // ---------------------------------------------------------

      //    gStyle->SetOptStat(1);
      gStyle->cd();

      sprintf(cName,"CA4");     
      TCanvas *CA4 = (TCanvas*) gROOT->FindObject(cName);
      if(CA4==NULL) CA4 = new TCanvas("CA4","Sliced p2-x2",600,800);
      CA4->cd();
      CA4->Clear();

      Int_t ndiv = 4;
      CA4->Divide(1,ndiv);

      TString fOutName2 = Form("./%s/Plots/Bunch/%s/Bunch-%s-p2x2_%i",sim.Data(),pData->GetSpeciesName(index).c_str(),sim.Data(),time);

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
      sprintf(text,"#LTy#GT_{rms} = %5.2f",xrms);
      textStatInt->AddText(text);
      sprintf(text,"#LTp_{y}#GT_{rms} = %5.2f",yrms);
      textStatInt->AddText(text);
      sprintf(text,"#varepsilon_{n} = %5.2f",emitx);
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

	// Float_t y1 = gPad->GetBottomMargin();
	Float_t y2 = 1 - gPad->GetTopMargin();
	Float_t x1 = gPad->GetLeftMargin();
	// Float_t x2 = 1 - gPad->GetRightMargin();
	textStat[k] = new TPaveText(x1+0.02,y2-0.50,x1+0.20,y2-0.05,"NDC");
	PGlobals::SetPaveTextStyle(textStat[k],12); 
	textStat[k]->SetTextColor(kGray+3);
	textStat[k]->SetTextFont(42);

	char text[64];
	sprintf(text,"%5.2f #LT #zeta #LT %5.2f",sBinLim[k],sBinLim[k+1]);
	textStat[k]->AddText(text);
	sprintf(text,"#LTx#GT_{rms} = %5.2f",sxrms[k]);
	textStat[k]->AddText(text);
	sprintf(text,"#LTp_{x}#GT_{rms} = %5.2f",syrms[k]);
	textStat[k]->AddText(text);
	sprintf(text,"#varepsilon_{n} = %5.2f",semit[k]);
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
      gErms->Write("gErms",TObject::kOverwrite);
      gXrms->Write("gXrms",TObject::kOverwrite);
      gEmitx->Write("gEmitx",TObject::kOverwrite);
      ofile->Close();
    }
    
    // Delete newly created vectors
    delete sBinLim;

    delete sxmean;
    delete symean;
    delete sxrms2;  
    delete syrms2; 
    delete sxrms;  
    delete syrms;  
    delete sxyrms2;
    delete xbin;
    delete semit;
    delete sEmean;
    delete sErms;

    // delete gemit;
    // delete gXrms;
    // delete gErms;
    // delete gErmsB;

    // end time looper
  }
}
