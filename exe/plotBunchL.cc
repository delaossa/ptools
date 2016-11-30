#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <string>

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
#include <TRandom3.h>

#include "PData.hh"
#include "PGlobals.hh"
#include "PPalette.hh"
#include "H5Cpp.h"

using namespace std;
using namespace H5;

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

void DeleteOverflows(TH2 *h) {
  Int_t Nx = h->GetXaxis()->GetNbins();
  Int_t Ny = h->GetYaxis()->GetNbins();
  for(Int_t i=0;i<=Nx+1;i++) {
    for(Int_t j=0;j<=Ny+1;j++) {
      if(i == 0 || j == 0 || i == Nx+1 || j == Ny+1)
	h->SetBinContent(i,j,0.0);
    }
  }

}

int main(int argc,char *argv[]) {
  if(argc<=2) {
    printf("\n Usage: %s <input file> <--astra>\n",argv[0]);
    printf("      <--png> <--pdf> <--eps> <--center>\n");
    printf("      <--file> \n");
    return 0;
  }

  // Input file
  TString ifile = "";

  // General options
  TString opt = "";

  // Longitudinal range
  Float_t zmin0 =  99999.;
  Float_t zmax0 = -99999.;
  
  // Options for Spectrum
  Float_t Pmin =  99999.;
  Float_t Pmax = -99999.;

  // Option for emittance spoiler
  Float_t X = 1000.; // um
  TRandom3 *rndEngine = new TRandom3();

  // Interfacing command line:
  for(int l=1;l<argc;l++){
    TString arg = argv[l];

    if(arg.Contains("--pdf")) {
      opt += "pdf";
    } else if(arg.Contains("--eps")) {
      opt += "eps";
    } else if(arg.Contains("--png")) {
      opt += "png";
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
    } else if(arg.Contains("--notext")){
      opt += "notext"; 
    } else if(arg.Contains("--noinfo")){
      opt += "noinfo"; 
    } else if(arg.Contains("--bw")){
      opt += "bw"; 
    } else if(arg.Contains("--elli")){
      opt += "elli"; 
    } else if(arg.Contains("--astra")){
      opt += "astra"; 
    } else if(arg.Contains("-spoil")) {
      char ss[6];
      sscanf(arg,"%6s%f",ss,&X);
      opt += "spoil";
    } else if(arg.Contains("-zmin")) {
      char ss[5];
      sscanf(arg,"%5s%f",ss,&zmin0);
    } else if(arg.Contains("-zmax")) {
      char ss[5];
      sscanf(arg,"%5s%f",ss,&zmax0);
    } else if(arg.Contains("-pmin")) {
      char ss[5];
      sscanf(arg,"%5s%f",ss,&Pmin);
    } else if(arg.Contains("-pmax")) {
      char ss[5];
      sscanf(arg,"%5s%f",ss,&Pmax);
    } else {
      ifile = arg;
    }
  }
  

  PGlobals::Initialize();

  // Palettes!
  //  gROOT->Macro("PPalettes.C");

  if(opt.Contains("grid")) {
    gStyle->SetPadGridX(1);
    gStyle->SetPadGridY(1);
  }

  // Spatial resolution
  // Float_t dx1 = 0.04;
  // Float_t dx2 = 0.04;
  // Float_t dx3 = 0.04;

  // Bining, intervals, labels, etc.
  Int_t x1Nbin = 120;
  Int_t p1Nbin = 120;
  Int_t x2Nbin = 120;
  Int_t p2Nbin = 120;
  Int_t x3Nbin = 120;
  Int_t p3Nbin = 120;
    
  // Slices
  Int_t SNbin = 40;
  Float_t x1BinMin;
  Float_t x1BinMax;

  // Spatial coordinates intervals:
  // Float_t x1Min = -7.8;
  Float_t x1Min;
  Float_t x1Max;
  Float_t x2Min;
  Float_t x2Max;
  Float_t x3Min;
  Float_t x3Max;

  // Momentum coordinates intervals:
  Float_t p1Min;
  Float_t p1Max;
  Float_t p2Min;
  Float_t p2Max;
  Float_t p3Min;
  Float_t p3Max;

  // Spoiler
  const Float_t X0 = 35.0E4; // um (Berylium)

  // Filling histos
  cout << Form("\n 1. Reading RAW file ... ") ;

  if(opt.Contains("spoil"))
    cout << Form("\n  -> Spoiling emittance with %.1f um of Berylium foil ...\n",X);

  // --------------------------------------------------
  // READ FROM TEXT ELEGANT OR ASTRA FILE
  ifstream file(ifile.Data());
  
  Double_t Q = 0.0;
  UInt_t Np = 0;  
  const Int_t Nvar = 6;
  Double_t *var[Nvar]; 
  Double_t varMean[Nvar]; 
  Double_t varRms[Nvar]; 
  Double_t varMin[Nvar]; 
  Double_t varMax[Nvar]; 
  char varname[Nvar][8] = {{"x1"},{"x2"},{"x3"},{"p1"},{"p2"},{"p3"}};
  for(Int_t i=0;i<Nvar;i++) {
    varMean[i] = 0.0;
    varRms[i] = 0.0;
    varMin[i] = 1E20;
    varMax[i] = -1E20;
  }

  if(!opt.Contains("astra")) {
    string str; 
    Int_t irow = 0;
    while (std::getline(file, str)) {
      istringstream stream(str);

      //  cout << str << endl;
   
      // Process str    
      if(irow==0) {
	stream >> Q;
	Q *= 1E12;
      } else if(irow==1)
	stream >> Np;
      else if(irow==2) { // allocate memory for the variables
	for(Int_t i=0;i<Nvar;i++)
	  var[i] = new Double_t[Np];
      } else {
	for(Int_t i=0;i<Nvar;i++) {
	  stream >> var[i][irow-2];	
	}
 
	// transform spatial coordinates
	var[0][irow-2] *= -PConst::c_light*1E6; // um
	var[1][irow-2] *= 1E6;  // um
	var[2][irow-2] *= 1E6;  // um

	// Emittance spoiler:
	if(opt.Contains("spoil")) {
	  // \Theta_0 from the multiple scattering model:
	  // http://pdg.lbl.gov/2014/reviews/rpp2014-rev-passage-particles-matter.pdf
	
	  Float_t E0 = var[3][irow-2] * PConst::ElectronMassE/PUnits::MeV;
	  Float_t Theta0 = (13.6/E0) * TMath::Sqrt(X/X0) * (1.0 - 0.038 * TMath::Log(X/X0)); // rad
	  Float_t corr = 0.87;
	  Float_t xr1, xr2;
	  rndEngine->Rannor(xr1,xr2);
	  Float_t xplane  = (xr1 * X * Theta0 / TMath::Sqrt(12.)) + (xr2 * X * Theta0 / 2); // um
	  Float_t txplane = xr2 * Theta0; // rad
	  var[1][irow-2] += xplane;
	  var[4][irow-2] += txplane;

	  //cout << Form(" rdn1 = %6f  rdn2 = %6f  xplane = %6f  txplane = %6f",xr1,xr2,xplane,txplane) << endl;
	
	  Float_t yr1, yr2;
	  rndEngine->Rannor(yr1,yr2);
	  Float_t yplane  = yr1 * X * Theta0 / TMath::Sqrt(12.) + yr2 * X * Theta0 / 2; // um
	  Float_t typlane = yr2 * Theta0; // rad
	  var[2][irow-2] += yplane;
	  var[5][irow-2] += typlane;
	
	
	
	}
      

	var[4][irow-2] *= var[3][irow-2]; // mc
	var[5][irow-2] *= var[3][irow-2]; // mc
	var[3][irow-2] *= PConst::ElectronMassE/PUnits::GeV; // GeV

      
	for(Int_t i=0;i<Nvar;i++) {
	  varMean[i] += var[i][irow-2];
	  varRms[i]  += var[i][irow-2]*var[i][irow-2];	
	
	  if(var[i][irow-2]<varMin[i]) varMin[i] = var[i][irow-2];
	  if(var[i][irow-2]>varMax[i]) varMax[i] = var[i][irow-2];
	
	}      
      }
    
      irow++;
    }
  } else {

    cout << "   ASTRA file: " << ifile.Data() << endl;
    
    string str; 
    Int_t irow = 0;
    Double_t t,q;
    Int_t pi,ps;

    // Count the lines
    while (std::getline(file, str)) ++Np;
    for(Int_t i=0;i<Nvar;i++)
      var[i] = new Double_t[Np];

    // Rewind
    file.clear();
    file.seekg(0);
    while (std::getline(file, str)) {
      istringstream stream(str);

      //      cout << str << endl;
      
      for(Int_t i=0;i<Nvar;i++) {
	stream >> var[i][irow];	
      }
      
      stream >> t;
      stream >> q;
      q *= -1;
      stream >> pi;
      stream >> ps;
      
      if(pi!=1) continue;
      if(ps!=5) continue;
      
      // transform spatial coordinates
      Double_t aux0 = var[0][irow];
      Double_t aux1 = var[1][irow];
      var[0][irow] = var[2][irow] * 1E6;  // um;
      var[1][irow] = aux0 * 1E6;  // um;
      var[2][irow] = aux1 * 1E6;  // um;

      // transform momentum coordinates
      aux0 = var[3][irow];
      aux1 = var[4][irow];
      var[3][irow] = var[5][irow] * PUnits::eV / (PUnits::GeV);
      var[4][irow] = aux0 * PUnits::eV / PConst::ElectronMassE;
      var[5][irow] = aux1 * PUnits::eV / PConst::ElectronMassE;
      
      if(irow>0) {
	var[0][irow] += var[0][0];
	var[3][irow] += var[3][0];
      }
      
      Q += q;

      for(Int_t i=0;i<Nvar;i++) {
	varMean[i] += var[i][irow];
	varRms[i]  += var[i][irow]*var[i][irow];	
	
	if(var[i][irow]<varMin[i]) varMin[i] = var[i][irow];
	if(var[i][irow]>varMax[i]) varMax[i] = var[i][irow];
      }
      
      irow++;
    }

    Q *= 1E3; // pC
  }


  file.close();
  
  for(Int_t i=0;i<Nvar;i++) {
    varMean[i] /= Np;
    varRms[i]  /= Np;
    varRms[i] = TMath::Sqrt( varRms[i] -  varMean[i]*varMean[i] );
  }
  
  cout << Form("  %i  particles read!    Charge = %.1f pC" , Np, Q);

  cout << Form("  \n x1 = %e +/- %e \n x2 = %.4f +/- %.4f \n x3 = %.4f +/- %.4f \n p1 = %.4f +/- %.4f \n p2 = %.4f +/- %.4f \n p3 = %.4f +/- %.4f",
	       // 	       varMean[0]/PUnits::um, varRms[0]/PUnits::um, varMean[1]/PUnits::um, varRms[1]/PUnits::um, varMean[2]/PUnits::um, varRms[2]/PUnits::um, varMean[3]/PUnits::GeV, varRms[3]/PUnits::GeV, varMean[4]/PUnits::MeV, varRms[4]/PUnits::MeV, varMean[5]/PUnits::MeV, varRms[5]/PUnits::MeV) << endl;      
   	       varMean[0], varRms[0], varMean[1], varRms[1], varMean[2], varRms[2], varMean[3], varRms[3], varMean[4], varRms[4], varMean[5], varRms[5]) << endl;      
  
  
  
  // ----------------------------------------------------------------------------------------------

  // Ranges
  Float_t rangefactor = 0.2;
  
  x1BinMin = varMin[0];
  x1BinMax = varMax[0];

  x1Min = varMin[0] - rangefactor*(varMax[0]-varMin[0]);
  x1Max = varMax[0] + rangefactor*(varMax[0]-varMin[0]);

  x2Min = varMin[1] - rangefactor*(varMax[1]-varMin[1]);
  x2Max = varMax[1] + rangefactor*(varMax[1]-varMin[1]);

  x3Min = varMin[2] - rangefactor*(varMax[2]-varMin[2]);
  x3Max = varMax[2] + rangefactor*(varMax[2]-varMin[2]);

  p1Min = varMin[3] - rangefactor*(varMax[3]-varMin[3]);
  p1Max = varMax[3] + rangefactor*(varMax[3]-varMin[3]);

  p2Min = varMin[4] - rangefactor*(varMax[4]-varMin[4]);
  p2Max = varMax[4] + rangefactor*(varMax[4]-varMin[4]);

  p3Min = varMin[5] - rangefactor*(varMax[5]-varMin[5]);
  p3Max = varMax[5] + rangefactor*(varMax[5]-varMin[5]);

  if(ifile.Contains("MCELLSTART")) {
    x1Min = varMean[0] - 500;
    x1Max = varMean[0] + 300;

    x2Min = varMean[1] - 300;
    x2Max = varMean[1] + 300;

    x3Min = varMean[2] - 300;
    x3Max = varMean[2] + 300;

    p2Min = varMean[4] - 8;
    p2Max = varMean[4] + 8;

    p3Min = varMean[5] - 8;
    p3Max = varMean[5] + 8;
  }

  TH1F *hScanX1 = (TH1F*) gROOT->FindObject("hScanX1");
  if(hScanX1) delete hScanX1;
  hScanX1 = new TH1F("hScanX1","",x1Nbin,x1Min,x1Max);
  TH1F *hScanX2 = (TH1F*) gROOT->FindObject("hScanX2");
  if(hScanX2) delete hScanX2;
  hScanX2 = new TH1F("hScanX2","",x2Nbin,x2Min,x2Max);
  TH1F *hScanX3 = (TH1F*) gROOT->FindObject("hScanX3");
  if(hScanX3) delete hScanX3;
  hScanX3 = new TH1F("hScanX3","",x3Nbin,x3Min,x3Max);
  
  TH1F *hScanP1 = (TH1F*) gROOT->FindObject("hScanP1");
  if(hScanP1) delete hScanP1;
  hScanP1 = new TH1F("hScanP1","",p1Nbin,p1Min,p1Max);
  TH1F *hScanP2 = (TH1F*) gROOT->FindObject("hScanP2");
  if(hScanP2) delete hScanP2;
  hScanP2 = new TH1F("hScanP2","",p2Nbin,p2Min,p2Max);
  TH1F *hScanP3 = (TH1F*) gROOT->FindObject("hScanP3");
  if(hScanP3) delete hScanP3;
  hScanP3 = new TH1F("hScanP3","",p3Nbin,p3Min,p3Max);

  // Filling scan histos for auto ranging
  
  for(UInt_t i=0;i<Np;i++) {

    hScanX1->Fill(var[0][i],Q/Np);
    hScanP1->Fill(var[3][i],Q/Np);
    hScanX2->Fill(var[1][i],Q/Np);
    hScanP2->Fill(var[4][i],Q/Np);
    hScanX3->Fill(var[2][i],Q/Np);
    hScanP3->Fill(var[5][i],Q/Np);

  }

  // Set limits using Scan histograms (auto ranges)
  Double_t peakFactor = 0.05;
  Double_t rfactor = 0.5;

  Double_t min, max;
  FindLimits(hScanX1,min,max,peakFactor);
  x1Min = min - rfactor*(max-min);
  x1Max = max + rfactor*(max-min);

  FindLimits(hScanX2,min,max,peakFactor);
  x2Min = min - rfactor*(max-min);
  x2Max = max + rfactor*(max-min);

  FindLimits(hScanX3,min,max,peakFactor);
  x3Min = min - rfactor*(max-min);
  x3Max = max + rfactor*(max-min);

  FindLimits(hScanP1,min,max,peakFactor);
  p1Min = min - rfactor*(max-min);
  p1Max = max + 2*rfactor*(max-min);

  FindLimits(hScanP2,min,max,peakFactor);
  p2Min = min - rfactor*(max-min);
  p2Max = max + rfactor*(max-min);

  FindLimits(hScanP3,min,max,peakFactor);
  p3Min = min - rfactor*(max-min);
  p3Max = max + rfactor*(max-min);

  // Slices range
  FindLimits(hScanX1,min,max,0.1);
  x1BinMin = min;
  x1BinMax = max;

  
  cout << Form(" x1 range (N = %i):  x1Min = %f  x1Max = %f ", x1Nbin, x1Min, x1Max) << endl;
  cout << Form(" p1 range (N = %i):  p1Min = %f  p1Max = %f ", p1Nbin, p1Min, p1Max) << endl;
  cout << Form(" x2 range (N = %i):  x2Min = %f  x2Max = %f ", x2Nbin, x2Min, x2Max) << endl;
  cout << Form(" p2 range (N = %i):  p2Min = %f  p2Max = %f ", p2Nbin, p2Min, p2Max) << endl;
  cout << Form(" x3 range (N = %i):  x3Min = %f  x3Max = %f ", x3Nbin, x3Min, x3Max) << endl;
  cout << Form(" p3 range (N = %i):  p3Min = %f  p3Max = %f ", p3Nbin, p3Min, p3Max) << endl;
    
  cout <<  Form(" Number of bins = %i: x1BinMin = %f  x1BinMax = %f ", x1Nbin,x1BinMin,x1BinMax) << endl;    
    
  //  return 0;

  // Histograms
  char hName[8];

  sprintf(hName,"hX1");
  TH1F *hX1 = (TH1F*) gROOT->FindObject(hName);
  if(hX1) delete hX1;
  hX1 = new TH1F(hName,"",x1Nbin,x1Min,x1Max);

  if(ifile.Contains("MCELLSTART"))
    hX1->GetYaxis()->SetTitle("I [0.1 kA]");
  else
    hX1->GetYaxis()->SetTitle("I [kA]");
  
  hX1->GetXaxis()->SetTitle("#zeta [#mum]");

  sprintf(hName,"hP1");
  TH1F *hP1 = (TH1F*) gROOT->FindObject(hName);
  if(hP1) delete hP1;
  hP1 = new TH1F(hName,"",p1Nbin,p1Min,p1Max);
  hP1->GetYaxis()->SetTitle("p_{z} [GeV/c]");
  hP1->GetXaxis()->SetTitle("#zeta [#mum]");

  sprintf(hName,"hP1X1");
  TH2F *hP1X1 = (TH2F*) gROOT->FindObject(hName);
  if(hP1X1) delete hP1X1;
  hP1X1 = new TH2F(hName,"",x1Nbin,x1Min,x1Max,p1Nbin,p1Min,p1Max);
  hP1X1->GetXaxis()->SetTitle("#zeta [#mum]");
  hP1X1->GetYaxis()->SetTitle("p_{z} [GeV/c]");
  hP1X1->GetZaxis()->SetTitle("Charge [pC]");

  sprintf(hName,"hP2X1");
  TH2F *hP2X1 =  (TH2F*) gROOT->FindObject(hName);
  if(hP2X1) delete hP2X1;
  hP2X1 = new TH2F(hName,"",x1Nbin,x1Min,x1Max,p2Nbin,p2Min,p2Max);
  hP2X1->GetXaxis()->SetTitle("#zeta [#mum]");
  hP2X1->GetYaxis()->SetTitle("p_{x}/mc");
  hP2X1->GetZaxis()->SetTitle("Charge [pC]");
  hP2X1->GetZaxis()->CenterTitle();

  sprintf(hName,"hP3X1");
  TH2F *hP3X1 =  (TH2F*) gROOT->FindObject(hName);
  if(hP3X1) delete hP3X1;
  hP3X1 = new TH2F(hName,"",x1Nbin,x1Min,x1Max,p3Nbin,p3Min,p3Max);
  hP3X1->GetXaxis()->SetTitle("#zeta [#mum]");
  hP3X1->GetYaxis()->SetTitle("p_{y}/mc");
  hP3X1->GetZaxis()->SetTitle("Charge [pC]");
  hP3X1->GetZaxis()->CenterTitle();

  sprintf(hName,"hP2X2");
  TH2F *hP2X2 =  (TH2F*) gROOT->FindObject(hName);
  if(hP2X2) delete hP2X2;
  hP2X2 = new TH2F(hName,"",x2Nbin,x2Min,x2Max,p2Nbin,p2Min,p2Max);
  hP2X2->GetXaxis()->SetTitle("x [#mum]");
  hP2X2->GetYaxis()->SetTitle("p_{x}/mc");
  hP2X2->GetZaxis()->SetTitle("Charge [pC]");
  hP2X2->GetZaxis()->CenterTitle();

  sprintf(hName,"hP3X3");
  TH2F *hP3X3 =  (TH2F*) gROOT->FindObject(hName);
  if(hP3X3) delete hP3X3;
  hP3X3 = new TH2F(hName,"",x3Nbin,x3Min,x3Max,p3Nbin,p3Min,p3Max);
  hP3X3->GetXaxis()->SetTitle("y [#mum]");
  hP3X3->GetYaxis()->SetTitle("p_{y}/mc");
  hP3X3->GetZaxis()->SetTitle("Charge [pC]");
  hP3X3->GetZaxis()->CenterTitle();

  sprintf(hName,"hX2X1");
  TH2F *hX2X1 = (TH2F*) gROOT->FindObject(hName);
  if(hX2X1) delete hX2X1;
  hX2X1 = new TH2F(hName,"",x1Nbin,x1Min,x1Max,x2Nbin,x2Min,x2Max);
  hX2X1->GetXaxis()->SetTitle("#zeta [#mum]");
  hX2X1->GetYaxis()->SetTitle("x [#mum]");
  hX2X1->GetZaxis()->SetTitle("Charge [pC]");

  sprintf(hName,"hX3X1");
  TH2F *hX3X1 = (TH2F*) gROOT->FindObject(hName);
  if(hX3X1) delete hX3X1;
  hX3X1 = new TH2F(hName,"",x1Nbin,x1Min,x1Max,x3Nbin,x3Min,x3Max);
  hX3X1->GetXaxis()->SetTitle("#zeta [#mum]");
  hX3X1->GetYaxis()->SetTitle("y [#mum]");
  hX3X1->GetZaxis()->SetTitle("Charge [pC]");

  sprintf(hName,"hX3X2");
  TH2F *hX3X2 = (TH2F*) gROOT->FindObject(hName);
  if(hX3X2) delete hX3X2;
  hX3X2 = new TH2F(hName,"",x2Nbin,x2Min,x2Max,x3Nbin,x3Min,x3Max);
  hX3X2->GetXaxis()->SetTitle("x [#mum]");
  hX3X2->GetYaxis()->SetTitle("y [#mum]");
  hX3X2->GetZaxis()->SetTitle("Charge [pC]");


   // Filling histos
  cout << Form("\n 2. Filling general histograms ... ") ;
  

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
    sprintf(hName,"hP2X2sl_%i",k);
    hP2X2sl[k] = (TH2F*) gROOT->FindObject(hName);
    if(hP2X2sl[k]) delete hP2X2sl[k];
    hP2X2sl[k] = new TH2F(hName,"",x2Nbin,x2Min,x2Max,p2Nbin,p2Min,p2Max);

    hP2X2sl[k]->GetXaxis()->SetTitle("x [#mum]");
    hP2X2sl[k]->GetYaxis()->SetTitle("p_{x}/mc");
    hP2X2sl[k]->GetZaxis()->SetTitle("Charge [pC]");
 
    sprintf(hName,"hP3X3sl_%i",k);
    hP3X3sl[k] = (TH2F*) gROOT->FindObject(hName);
    if(hP3X3sl[k]) delete hP3X3sl[k];
    hP3X3sl[k] = new TH2F(hName,"",x3Nbin,x3Min,x3Max,p3Nbin,p3Min,p3Max);

    hP3X3sl[k]->GetXaxis()->SetTitle("y [#mum]");
    hP3X3sl[k]->GetYaxis()->SetTitle("p_{y}/mc");
    hP3X3sl[k]->GetZaxis()->SetTitle("Charge [pC]");
 
    sprintf(hName,"hP1sl_%i",k);
    hP1sl[k] = (TH1F*) gROOT->FindObject(hName);
    if(hP1sl[k]) delete hP1sl[k];
    hP1sl[k] = new TH1F(hName,"",p1Nbin,p1Min,p1Max);

  }

  // Filling histos
  cout << Form("\n 2. Filling sliced histograms ... ") ;

  
  for(UInt_t i=0;i<Np;i++) {
    
    hX1->Fill(var[0][i],Q/Np);
    hP1->Fill(var[3][i],Q/Np);
    hP1X1->Fill(var[0][i],var[3][i],Q/Np);
    
    hP2X1->Fill(var[0][i],var[4][i],Q/Np);
    hP3X1->Fill(var[0][i],var[5][i],Q/Np);

    hP2X2->Fill(var[1][i],var[4][i],Q/Np);
    hP3X3->Fill(var[2][i],var[5][i],Q/Np);

    hX2X1->Fill(var[0][i],var[1][i],Q/Np);
    hX3X1->Fill(var[0][i],var[2][i],Q/Np);
    hX3X2->Fill(var[1][i],var[2][i],Q/Np);
    
    // Slices    
    if(var[0][i]<sBinLim[0] || var[0][i]>sBinLim[SNbin]) continue;
    Int_t iBin = -1;
    for(Int_t j=0; j<SNbin; j++) {
      
      if(var[0][i]<sBinLim[j+1]) {
	iBin = j;
	break;
      }
    }
    if(iBin<0) continue;

    // Projected emittance in the bunch range. (skip the tails to "match" the sliced ones)

    hP2X2sl[iBin]->Fill(var[1][i],var[4][i],Q/Np);
    hP3X3sl[iBin]->Fill(var[2][i],var[5][i],Q/Np);
    hP1sl[iBin]->Fill(var[3][i],Q/Np);   
    
  }

  //
  DeleteOverflows(hP1X1);
  DeleteOverflows(hX2X1);
  DeleteOverflows(hX3X1);
  DeleteOverflows(hP2X1);
  DeleteOverflows(hP3X1);
  DeleteOverflows(hX3X2);
  DeleteOverflows(hP2X2);
  DeleteOverflows(hP3X3);
  
  
  // Integrated long. emittance:
  cout << Form("\n 3. Calculating integrated quantities ... ") ;

  // Longitudinal phasespace:
  // --

  Double_t xmean  = 0.0;
  Double_t ymean  = 0.0;
  Double_t x2mean = 0.0;
  Double_t y2mean = 0.0;
  Double_t xymean = 0.0;
  Double_t Ntotal = 0.0;
  for(Int_t i=1;i<=x1Nbin;i++) {
    Double_t x = hP1X1->GetXaxis()->GetBinCenter(i);
    // if(x<xmin || x>xmax) continue;
    for(Int_t j=1;j<=p1Nbin;j++) {
      Double_t y = hP1X1->GetYaxis()->GetBinCenter(j);
      // if(y<ymin || y>ymax) continue;
      Double_t value = TMath::Abs(hP1X1->GetBinContent(i,j));
      xmean += x*value;
      ymean += y*value;
      x2mean += x*x*value;
      y2mean += y*y*value;
      xymean += x*y*value;

      Ntotal += value;
    }
  }

  xmean  /= Ntotal;
  ymean  /= Ntotal;
  x2mean /= Ntotal;
  y2mean /= Ntotal;
  xymean /= Ntotal;

  

  Double_t xrms2  = x2mean - xmean*xmean;
  Double_t yrms2  = y2mean - ymean*ymean;
  Double_t xrms   = ( xrms2>0.0 ) ? TMath::Sqrt(xrms2) : 0.0 ;
  Double_t yrms   = ( yrms2>0.0 ) ? TMath::Sqrt(yrms2) : 0.0 ;
  Double_t xyrms2 = xymean - xmean*ymean;
  Double_t emit2 = xrms2*yrms2 - xyrms2*xyrms2;
  Double_t emit = (emit2>0.0) ? TMath::Sqrt(emit2) : 0.0 ;
    
  // cout << Form("  xMean = %7.3f   yMean = %7.3f",xmean,ymean) << endl;
  // cout << Form("  xRms  = %7.3f   yRms  = %7.3f",xrms,yrms) << endl;
  // cout << Form("  Emittance = %7.3f",emit) << endl;

  Double_t emitz = emit;
  Double_t zmean = xmean;
  Double_t zrms = xrms;
  Double_t pzmean = ymean;
  Double_t pzrms = yrms;
  // ----------------------------------------


  // Transverse phasespace

  Double_t stats[7];  // { sumw, sumw2, sumwx, sumwx2, sumwy, sumwy2, sumwxy }

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
  Double_t pxmean = stats[4]/stats[0];
  Double_t pxrms = (yrms2>0.0) ? TMath::Sqrt(yrms2) : 0.0 ;

  Double_t beta = xrms2 / emit;
  Double_t gamma = yrms2 / emit;
  Double_t alpha = -xyrms2 / emit;
  
  Double_t grel = hP1->GetMean() * PUnits::GeV/PConst::ElectronMassE;
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
  
  TEllipse *ellipP2X2 = new TEllipse(x_mean,pxmean,TMath::Sqrt(a2),TMath::Sqrt(b2),0.,360.,angle * 180. / PConst::pi );
  ellipP2X2->SetFillStyle(0);
  ellipP2X2->SetLineStyle(2);
  ellipP2X2->SetLineColor(2);
  ellipP2X2->SetLineWidth(1);


  // P3X3 SPACE -----------------------
  hP3X3->GetStats(stats);

  xrms2  = stats[3]/stats[0] - stats[2]*stats[2]/(stats[0]*stats[0]);
  yrms2  = stats[5]/stats[0] - stats[4]*stats[4]/(stats[0]*stats[0]);
  xyrms2 = stats[6]/stats[0] - stats[2]*stats[4]/(stats[0]*stats[0]);
  emit2 = xrms2*yrms2 - xyrms2*xyrms2;
  emit = (emit2>0.0) ? TMath::Sqrt(emit2) : 0.0 ; 

  Double_t emity = emit;
  Double_t y_mean = stats[2]/stats[0];
  Double_t y_rms = (xrms2>0.0) ? TMath::Sqrt(xrms2) : 0.0 ;
  Double_t pymean = stats[4]/stats[0];
  Double_t pyrms = (yrms2>0.0) ? TMath::Sqrt(yrms2) : 0.0 ;

  beta = xrms2 / emit;
  gamma = yrms2 / emit;
  alpha = -xyrms2 / emit;
  
  Double_t betay = beta * grel;

  factor =  beta*beta + 2 * beta * gamma + gamma*gamma - 4 * emit;
  if(factor<0.0) factor *= -1;
  factor = TMath::Sqrt(factor);
  
  a2 = 0.5 * (beta + gamma - factor);
  b2 = 0.5 * (beta + gamma + factor);
  
  p  = alpha / (a2-b2);
  pf = TMath::Sqrt( 1 - 4*p*p ); 
  angle = - TMath::ATan( (1+pf) / (2*p) );
  //  if (angle <0.0) angle += 2*PConst::pi;

  TEllipse *ellipP3X3 = new TEllipse(x_mean,pxmean,TMath::Sqrt(a2),TMath::Sqrt(b2),0.,360.,angle * 180. / PConst::pi );
  ellipP3X3->SetFillStyle(0);
  ellipP3X3->SetLineStyle(2);
  ellipP3X3->SetLineColor(2);
  ellipP3X3->SetLineWidth(1);

  // Charge  
  Double_t Charge = hX1->Integral();
  
  // cout << Form("  xMean = %7.3f   yMean = %7.3f",xmean,ymean) << endl;
  // cout << Form("  xRms  = %7.3f   yRms  = %7.3f",xrms,yrms) << endl;
  // cout << Form("  Emittance = %7.3f",emit) << endl;

  // Emittance units
  Double_t emitUnit = PUnits::um;
  string emitSUnit = "#mum";
  emitx *= emitUnit;
  emity *= emitUnit;
  //PUnits::BestUnit bemitSUnit(TMath::Sqrt(emity*emitx),"Emittance");
  PUnits::BestUnit bemitSUnit(emity,"Emittance");
  bemitSUnit.GetBestUnits(emitUnit,emitSUnit);
    
  emitx /= emitUnit;
  emity /= emitUnit;


  
  cout << Form("\n 4. Calculating sliced quantities ... ") << endl ;

  TGraph *gemitX = NULL;
  TGraph *gemitY = NULL;
  TGraph *gXrms = NULL;
  TGraph *gYrms = NULL;
  TGraph *gErms = NULL;

  TEllipse **sellipP2X2 = new TEllipse*[SNbin];

  Double_t * sxmean = new Double_t[SNbin];
  Double_t * symean = new Double_t[SNbin];
  Double_t * sx2mean = new Double_t[SNbin];
  Double_t * sy2mean = new Double_t[SNbin];
  Double_t * sxymean = new Double_t[SNbin];
  Double_t * sNtotal = new Double_t[SNbin];
  Double_t * sxrms2 = new Double_t[SNbin];  
  Double_t * syrms2 = new Double_t[SNbin]; 
  Double_t * sxrms = new Double_t[SNbin];  
  Double_t * syrms = new Double_t[SNbin];  
  Double_t * sxyrms2 = new Double_t[SNbin];

  Double_t * xbin = new Double_t[SNbin];
  Double_t * semit = new Double_t[SNbin];
  Double_t * semit2 = new Double_t[SNbin];

  Double_t * sNEtotal = new Double_t[SNbin]; 
  Double_t * sEmean = new Double_t[SNbin];

  Double_t * sErms = new Double_t[SNbin];
  Double_t * sErms2 = new Double_t[SNbin];

  Double_t * sCharge = new Double_t[SNbin];
  Double_t sChargeMax = -1E20;

  Double_t *sx_rms = new Double_t[SNbin];
  Double_t *spx_rms = new Double_t[SNbin];
  Double_t *semitx = new Double_t[SNbin];
  Double_t *sbetax = new Double_t[SNbin];
  

  for(Int_t k=0;k<SNbin;k++) {
    sxmean[k] = symean[k] = sx2mean[k] = sy2mean[k] = sxymean[k] 
      = sNtotal[k] = sxrms2[k] = syrms2[k] = sxrms[k] = syrms[k]
      = sxyrms2[k] = xbin[k] = semit[k] = 0.0;
    sNEtotal[k] = sEmean[k] = sErms[k] = 0.0;

    sellipP2X2[k] = NULL;

    xbin[k] = (sBinLim[k] + sBinLim[k+1])/2.;

    for(Int_t i=1;i<=x2Nbin;i++) {
      Double_t x = hP2X2sl[k]->GetXaxis()->GetBinCenter(i);
      // if(x<xmin || x>xmax) continue;
      for(Int_t j=1;j<=p2Nbin;j++) {
	Double_t y = hP2X2sl[k]->GetYaxis()->GetBinCenter(j);
	// if(y<ymin || y>ymax) continue;
	Double_t value = TMath::Abs(hP2X2sl[k]->GetBinContent(i,j));
	sxmean[k] += x*value;
	symean[k] += y*value;
	sx2mean[k] += x*x*value;
	sy2mean[k] += y*y*value;
	sxymean[k] += x*y*value;

	sNtotal[k] += value;
      }	
    }

    if( sNtotal[k] == 0.0 ) continue;

    sxmean[k]  /= sNtotal[k];
    symean[k]  /= sNtotal[k];
    sx2mean[k] /= sNtotal[k];
    sy2mean[k] /= sNtotal[k];
    sxymean[k] /= sNtotal[k];

    sxrms2[k]  = sx2mean[k] - sxmean[k]*sxmean[k];
    syrms2[k]  = sy2mean[k] - symean[k]*symean[k];
    sxrms[k]   = (sxrms2[k]>0.0) ? TMath::Sqrt(sxrms2[k]) : 0.0 ;
    syrms[k]   = (syrms2[k]>0.0) ? TMath::Sqrt(syrms2[k]) : 0.0 ;
    sxyrms2[k] = sxymean[k] - sxmean[k]*symean[k];
    semit2[k]  = sxrms2[k]*syrms2[k] - sxyrms2[k]*sxyrms2[k];
    semit[k]   = (semit2[k]>0.0) ? TMath::Sqrt(sxrms2[k]*syrms2[k] - sxyrms2[k]*sxyrms2[k]) : 0.0 ;
      
      
    sEmean[k]  = hP1sl[k]->GetMean();
    sErms[k]   = 1000*hP1sl[k]->GetRMS()/hP1sl[k]->GetMean();
    
    sCharge[k] = hP1sl[k]->Integral();
    if(sCharge[k]>sChargeMax) sChargeMax = sCharge[k];
    
    // cout<< Form("\nk = %i : (x1 > %f && x1 < %f)",k,sBinLim[k],sBinLim[k+1]) << endl; 

    // cout << Form("  xMean = %7.3f   yMean = %7.3f",sxmean[k],symean[k]) << endl;
    // cout << Form("  xRms  = %7.3f   yRms  = %7.3f",sxrms[k],syrms[k]) << endl;
    // cout << Form("  Emittance = %7.3f",semit[k]) << endl;

    // cout << Form("  Emean = %7.3f   Erms = %7.3f",sEmean[k],sErms[k]) << endl;
    
    Double_t beta = sxrms2[k] / semit[k];
    Double_t gamma = syrms2[k] / semit[k] ;
    Double_t alpha = - sxyrms2[k] / semit[k];

    Double_t factor =  beta*beta + 2 * beta * gamma + gamma*gamma - 4 * semit[k];
    if(factor<0.0) factor *= -1;
    factor = TMath::Sqrt(factor);

    Double_t a2 = 0.5 * (beta + gamma - factor);
    Double_t b2 = 0.5 * (beta + gamma + factor);
 
    Double_t p  = alpha / (a2-b2);
    Double_t pf = TMath::Sqrt( 1 - 4*p*p ); 
    Double_t angle = - TMath::ATan( (1+pf) / (2*p) );

    //  if (angle <0.0) angle += 2*PConst::pi;
    
    Double_t zoomEllip = 3.0;
    if(sNtotal[k] > 0.1)
      sellipP2X2[k] = new TEllipse(sxmean[k],symean[k],TMath::Sqrt(zoomEllip*a2),TMath::Sqrt(zoomEllip*b2),0.,360.,angle * 180. / PConst::pi );
    
    //    cout << Form("%i: %f %f %f %f",k+1,a2,b2,p,pf) << endl;
    // cout << Form("%i: %f %f %f %f %f",k+1,sxmean[k],symean[k],TMath::Sqrt(a2),TMath::Sqrt(b2),angle * 180. / PConst::pi) << endl; 
    sx_rms[k] = sxrms[k];
    spx_rms[k] = syrms[k];
    semitx[k] = semit[k] * PUnits::um / emitUnit;

    Double_t grel = hP1sl[k]->GetMean() * PUnits::GeV/PConst::ElectronMassE;
    sbetax[k] = beta * grel;

    DeleteOverflows(hP2X2sl[k]);
  }

  TEllipse **sellipP3X3 = new TEllipse*[SNbin];

  Double_t *sy_rms = new Double_t[SNbin];
  Double_t *spy_rms = new Double_t[SNbin];
  Double_t *semity = new Double_t[SNbin];
  Double_t *sbetay = new Double_t[SNbin];

  
  for(Int_t k=0;k<SNbin;k++) {
    sxmean[k] = symean[k] = sx2mean[k] = sy2mean[k] = sxymean[k] 
      = sNtotal[k] = sxrms2[k] = syrms2[k] = sxrms[k] = syrms[k]
      = sxyrms2[k] = xbin[k] = semit[k] = 0.0;
    sNEtotal[k] = sEmean[k] = sErms[k] = 0.0;

    sellipP3X3[k] = NULL;

    xbin[k] = (sBinLim[k] + sBinLim[k+1])/2.;

    for(Int_t i=1;i<=x3Nbin;i++) {
      Double_t x = hP3X3sl[k]->GetXaxis()->GetBinCenter(i);
      // if(x<xmin || x>xmax) continue;
      for(Int_t j=1;j<=p3Nbin;j++) {
	Double_t y = hP3X3sl[k]->GetYaxis()->GetBinCenter(j);
	// if(y<ymin || y>ymax) continue;
	Double_t value = TMath::Abs(hP3X3sl[k]->GetBinContent(i,j));
	sxmean[k] += x*value;
	symean[k] += y*value;
	sx2mean[k] += x*x*value;
	sy2mean[k] += y*y*value;
	sxymean[k] += x*y*value;

	sNtotal[k] += value;
      }	
    }

    if( sNtotal[k] == 0.0 ) continue;

    sxmean[k]  /= sNtotal[k];
    symean[k]  /= sNtotal[k];
    sx2mean[k] /= sNtotal[k];
    sy2mean[k] /= sNtotal[k];
    sxymean[k] /= sNtotal[k];

    sxrms2[k]  = sx2mean[k] - sxmean[k]*sxmean[k];
    syrms2[k]  = sy2mean[k] - symean[k]*symean[k];
    sxrms[k]   = (sxrms2[k]>0.0) ? TMath::Sqrt(sxrms2[k]) : 0.0 ;
    syrms[k]   = (syrms2[k]>0.0) ? TMath::Sqrt(syrms2[k]) : 0.0 ;
    sxyrms2[k] = sxymean[k] - sxmean[k]*symean[k];
    semit2[k]  = sxrms2[k]*syrms2[k] - sxyrms2[k]*sxyrms2[k];
    semit[k]   = (semit2[k]>0.0) ? TMath::Sqrt(sxrms2[k]*syrms2[k] - sxyrms2[k]*sxyrms2[k]) : 0.0 ;
      
      
    sEmean[k]  = hP1sl[k]->GetMean();
    sErms[k]   = 1000*hP1sl[k]->GetRMS()/hP1sl[k]->GetMean();
    
    sCharge[k] = hP1sl[k]->Integral();
    if(sCharge[k]>sChargeMax) sChargeMax = sCharge[k];
    
    // cout<< Form("\nk = %i : (x1 > %f && x1 < %f)",k,sBinLim[k],sBinLim[k+1]) << endl; 

    // cout << Form("  xMean = %7.3f   yMean = %7.3f",sxmean[k],symean[k]) << endl;
    // cout << Form("  xRms  = %7.3f   yRms  = %7.3f",sxrms[k],syrms[k]) << endl;
    // cout << Form("  Emittance = %7.3f",semit[k]) << endl;

    // cout << Form("  Emean = %7.3f   Erms = %7.3f",sEmean[k],sErms[k]) << endl;
    
    Double_t beta = sxrms2[k] / semit[k];
    Double_t gamma = syrms2[k] / semit[k] ;
    Double_t alpha = - sxyrms2[k] / semit[k];

    Double_t factor =  beta*beta + 2 * beta * gamma + gamma*gamma - 4 * semit[k];
    if(factor<0.0) factor *= -1;
    factor = TMath::Sqrt(factor);

    Double_t a2 = 0.5 * (beta + gamma - factor);
    Double_t b2 = 0.5 * (beta + gamma + factor);
 
    Double_t p  = alpha / (a2-b2);
    Double_t pf = TMath::Sqrt( 1 - 4*p*p ); 
    Double_t angle = - TMath::ATan( (1+pf) / (2*p) );

    //  if (angle <0.0) angle += 2*PConst::pi;

    Double_t zoomEllip = 3.0;
    if(sNtotal[k] > 0.1)
      sellipP3X3[k] = new TEllipse(sxmean[k],symean[k],TMath::Sqrt(zoomEllip*a2),TMath::Sqrt(zoomEllip*b2),0.,360.,angle * 180. / PConst::pi );

    sy_rms[k] = sxrms[k];
    spy_rms[k] = syrms[k];
    semity[k] = semit[k] * PUnits::um / emitUnit;

    Double_t grel = hP1sl[k]->GetMean() * PUnits::GeV/PConst::ElectronMassE;
    sbetay[k] = beta * grel;

    DeleteOverflows(hP3X3sl[k]);

  }

  
  cout << "\n done! " << endl;

  // To current.
  Float_t binSize = (x1Max - x1Min)/x1Nbin;
  Float_t lightspeed =  PConst::c_light / (PUnits::um/PUnits::femtosecond);  
  //  hX1->Scale(1000*lightspeed/binSize);
  hX1->Scale(lightspeed/binSize);


  Double_t curUnit;
  string curSUnit;
  PUnits::BestUnit bcurSUnit(hX1->GetMaximum(),"Current");
  bcurSUnit.GetBestUnits(curUnit,curSUnit);
  hX1->Scale(1/curUnit);
  hX1->GetYaxis()->SetTitle(Form("I[%s]",curSUnit.c_str()));
  
  // if(ifile.Contains("MCELLSTART")) 
  //   hX1->Scale(10);
  

  cout << "\n  Summary _______________________________________________________ " << endl;

  cout << Form("  Integrated charge (RAW) = %8f pC",Charge) << endl;
  cout << Form("  Peak current = %6.2f kA",hX1->GetMaximum()) << endl;
  cout << Form("  Total energy = %6.2f GeV, rms = %3.1f %%",pzmean,100*pzrms/pzmean) << endl;
  cout << Form("  Gamma        = %6.2f ",grel) << endl;
  cout << Form("  Trans. emit. = %6.2f  um",emitx) << endl;
  cout << Form("  Beta func.   = %6.2f  um",betax) << endl;
  cout << Form("  Bunch length = %6.2f  um (rms), width = %6.2f  um (rms)",zrms,x_rms) << endl;
  
  

  // ------------------------------------------------------------------------------------

  // Free memory from the dynamic array of variables:
  // for(UInt_t i=0;i<Nvar;i++) {
  //   delete var[i];
  // }
  // -----------------------------------------------------------------------------------
  
  cout << Form("\n 5. Preparing the graphs and plotting .. ") << endl << endl;

  // Centering in the x1 mean
  if(opt.Contains("center")) {
    hX1->SetBins(x1Nbin,x1Min-zmean,x1Max-zmean);
    hP1X1->SetBins(x1Nbin,x1Min-zmean,x1Max-zmean,p1Nbin,p1Min,p1Max);
    hX2X1->SetBins(x1Nbin,x1Min-zmean,x1Max-zmean,x2Nbin,x2Min-x_mean,x2Max-x_mean);
    hP2X1->SetBins(x1Nbin,x1Min-zmean,x1Max-zmean,p2Nbin,p2Min,p2Max);
    hP3X1->SetBins(x1Nbin,x1Min-zmean,x1Max-zmean,p3Nbin,p3Min,p3Max);
    hP2X2->SetBins(x2Nbin,x2Min-x_mean,x2Max-x_mean,p2Nbin,p2Min,p2Max);
    hP3X3->SetBins(x3Nbin,x3Min-y_mean,x3Max-y_mean,p3Nbin,p3Min,p3Max);
    hX3X1->SetBins(x1Nbin,x1Min-zmean,x1Max-zmean,x3Nbin,x3Min-y_mean,x3Max-y_mean);
    hX3X2->SetBins(x2Nbin,x2Min-x_mean,x2Max-x_mean,x3Nbin,x3Min-y_mean,x3Max-y_mean);
    
    for(Int_t k=0;k<SNbin;k++) {
      hP2X2sl[k]->SetBins(x2Nbin,x2Min-x_mean,x2Max-x_mean,p2Nbin,p2Min,p2Max);
      hP3X3sl[k]->SetBins(x3Nbin,x3Min-y_mean,x3Max-y_mean,p3Nbin,p3Min,p3Max);

      xbin[k] -= zmean;
      sBinLim[k] -= zmean;
    }
    sBinLim[SNbin] -= zmean;

    x1Min -= zmean;
    x1Max -= zmean;
    zmean = 0.0;

    if(zmax0>zmin0) {
      x1Min = zmin0;
      x1Max = zmax0;
    } else {
      if(-x1Min>x1Max) x1Min = -x1Max;
      else x1Max = -x1Min;
      
    }

    Float_t x2min = x2Min-x_mean;
    Float_t x2max = x2Max-x_mean;
    if(-x2min<x2max) x2min = -x2max;
    else x2max = -x2min;
    
    Float_t x3min = x3Min-y_mean;
    Float_t x3max = x3Max-y_mean;
    if(-x3min<x3max) x3min = -x3max;
    else x3max = -x3min;
   
    Float_t p2min = p2Min;
    Float_t p2max = p2Max;
    if(-p2min<p2max) p2min = -p2max;
    else p2max = -p2min;

    Float_t p3min = p3Min;
    Float_t p3max = p3Max;
    if(-p3min<p3max) p3min = -p3max;
    else p3max = -p3min;

    if(x2max<x3max) {
      x2max = x3max;
      x2min = x3min;
    } else {
      x3max = x2max;
      x3min = x2min;
    }

    if(p2max<p3max) {
      p2max = p3max;
      p2min = p3min;
    } else {
      p3max = p2max;
      p3min = p2min;
    }

    
    hP2X1->GetYaxis()->SetRangeUser(p2min,p2max);
    hP3X1->GetYaxis()->SetRangeUser(p3min,p3max);

    hP2X2->GetXaxis()->SetRangeUser(x2min,x2max);
    hP2X2->GetYaxis()->SetRangeUser(p2min,p2max);

    hP3X3->GetXaxis()->SetRangeUser(x3min,x3max);
    hP3X3->GetYaxis()->SetRangeUser(p3min,p3max);

    hX2X1->GetYaxis()->SetRangeUser(x2min,x2max);
    hX3X1->GetYaxis()->SetRangeUser(x3min,x3max);
    hX3X2->GetYaxis()->SetRangeUser(x3min,x3max);
    hX3X2->GetXaxis()->SetRangeUser(x2min,x2max);

    
    for(Int_t k=0;k<SNbin;k++) {
      hP2X2sl[k]->GetXaxis()->SetRangeUser(x2min,x2max);
      hP2X2sl[k]->GetYaxis()->SetRangeUser(p2min,p2max);
      
      hP3X3sl[k]->GetXaxis()->SetRangeUser(x3min,x3max);
      hP3X3sl[k]->GetYaxis()->SetRangeUser(p3min,p3max);
      
    }

    x_mean = 0;
    y_mean = 0;

    x2Min = x2min;
    x2Max = x2max;

    x3Min = x3min;
    x3Max = x3max;
    
    p2Min = p2min;
    p2Max = p2max;

    p3Min = p3min;
    p3Max = p3max;


    cout << Form(" x1 range (N = %i):  x1Min = %f  x1Max = %f ", x1Nbin, x1Min, x1Max) << endl;
    cout << Form(" p1 range (N = %i):  p1Min = %f  p1Max = %f ", p1Nbin, p1Min, p1Max) << endl;
    cout << Form(" x2 range (N = %i):  x2Min = %f  x2Max = %f ", x2Nbin, x2Min, x2Max) << endl;
    cout << Form(" p2 range (N = %i):  p2Min = %f  p2Max = %f ", p2Nbin, p2Min, p2Max) << endl;
    cout << Form(" x3 range (N = %i):  x3Min = %f  x3Max = %f ", x3Nbin, x3Min, x3Max) << endl;
    cout << Form(" p3 range (N = %i):  p3Min = %f  p3Max = %f ", p3Nbin, p3Min, p3Max) << endl;

    
  }
  // ------

  // Set palette colors
  Int_t * sPalette = new Int_t[SNbin];
  Double_t r[3]    = {0.90, 0.00, 1.00};
  Double_t g[3]    = {0.90, 0.00, 0.00};
  Double_t b[3]    = {0.90, 1.00, 0.00};
  Double_t stop[3] = {0.0, 0.3, 1.0};
  Int_t FI = TColor::CreateGradientColorTable(3, stop, r, g, b, SNbin);
  for (int i=0;i<SNbin;i++) sPalette[i] = FI+i;

  // Segmentation fault in the ellipses

  if(opt.Contains("elli") ) {
    for(Int_t k=0;k<SNbin;k++) {
      if(sellipP2X2[k] == NULL) continue;

      Int_t lineColor = sPalette[TMath::Nint( (SNbin-1) * sCharge[k]/sChargeMax)];
      sellipP2X2[k]->SetFillStyle(0);
      sellipP2X2[k]->SetLineColor(lineColor);
      sellipP2X2[k]->SetLineStyle(1);
      sellipP2X2[k]->SetLineWidth(1);
    }


    for(Int_t k=0;k<SNbin;k++) {
      if(sellipP3X3[k] == NULL) continue;

      Int_t lineColor = sPalette[TMath::Nint( (SNbin-1) * sCharge[k]/sChargeMax)];
      sellipP3X3[k]->SetLineColor(lineColor);
      sellipP3X3[k]->SetFillStyle(0);
      sellipP3X3[k]->SetLineStyle(1);
      sellipP3X3[k]->SetLineWidth(1);
    }
  }
  
  // Create the graph with the sliced quantities:
  gemitX = new TGraph(SNbin,xbin,semitx);
  gemitY = new TGraph(SNbin,xbin,semity);
  gXrms = new TGraph(SNbin,xbin,sx_rms);
  gYrms = new TGraph(SNbin,xbin,sy_rms);
  gErms = new TGraph(SNbin,xbin,sErms);

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
    Double_t xMin = x1Min;
    Double_t xMax = x1Min + (x1Max-x1Min) * 0.2;
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
    //      }

    Bool_t autoscale = kFALSE;
    // if(hX1->GetMaximum()>10) {
    // 	hX1->Scale(0.1);
    // 	autoscale = kTRUE;
    // }
      
    // Ranges!!
    Double_t yMin =  999.9;
    Double_t yMax =  -999.9;

    for(Int_t k=0;k<SNbin;k++) {
      if(semitx[k]<yMin)
	yMin = semitx[k];

      if(semitx[k]>yMax)
	yMax = semitx[k];

      if(semity[k]<yMin)
	yMin = semity[k];

      if(semity[k]>yMax)
	yMax = semity[k];

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
    C->cd();
    C->Clear();

    // Set palette:
    PPalette * pPalette = (PPalette*) gROOT->FindObject("electron0");
    if(!pPalette) {
      pPalette = new PPalette("electron0");
      pPalette->SetPalette("electron0");
      //pPalette->SetAlpha(1);
    }
       
    TString drawopt = "colz";
    TString drawopts = "colz same";
    
    // Double_t Max  = hP1X1->GetMaximum();
    // Double_t Min  = hP1X1->GetMinimum();

    // hP1X1->GetZaxis()->SetRangeUser(Min,Max); 

    // Text objects
    char ctext[128];
       
    TPaveText *textCharge = new TPaveText(0.15,0.25,0.48,0.3,"NDC");
    PGlobals::SetPaveTextStyle(textCharge,12); 
    textCharge->SetTextColor(kGray+2);
    sprintf(ctext,"Q = %5.1f pC", Charge);
    textCharge->AddText(ctext);

    TPaveText *textMom = new TPaveText(0.55,0.07,0.80,0.17,"NDC");
    PGlobals::SetPaveTextStyle(textMom,32); 
    textMom->SetTextColor(kGray+3);
    textMom->SetTextFont(62);
    sprintf(ctext,"#LTp_{z}#GT = %5.2f GeV/c", pzmean);
    textMom->AddText(ctext);


    TPaveText *textInfo = new TPaveText(0.55,0.40,0.80,0.85,"NDC");
    PGlobals::SetPaveTextStyle(textInfo,32); 
    textInfo->SetTextColor(kGray+2);
    textInfo->SetTextFont(42);
    sprintf(ctext,"Charge = %5.1f pC",Charge);
    textInfo->AddText(ctext);
    sprintf(ctext,"#Delta#zeta = %5.2f #mum",zrms);
    textInfo->AddText(ctext);
    // sprintf(ctext,"#LTp_{z}#GT_{rms} = %5.2f GeV/c",pzrms);
    sprintf(ctext,"#Delta#gamma/#gamma = %4.1f %%",100*pzrms/pzmean);
    textInfo->AddText(ctext);
    sprintf(ctext,"#varepsilon_{n,x} = %5.2f %s",emitx,emitSUnit.c_str());
    textInfo->AddText(ctext);
    sprintf(ctext,"#varepsilon_{n,y} = %5.2f %s",emity,emitSUnit.c_str());
    textInfo->AddText(ctext);
    sprintf(ctext,"#beta_{x} = %5.0f mm",1E-3*betax);
    textInfo->AddText(ctext);
    sprintf(ctext,"#beta_{y} = %5.0f mm",1E-3*betay);
    textInfo->AddText(ctext);

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
    PGlobals::CanvasAsymPartition(C,NPad,lMargin,rMargin,bMargin,tMargin,factor,0.028);

    // Define the frames for plotting
    Int_t fonttype = 43;
    Int_t fontsize = 24;
    Int_t tfontsize = 28;
    Double_t txoffset = 2.0;
    Double_t lxoffset = 0.02;
    Double_t tyoffset = 1.3;
    Double_t lyoffset = 0.01;
    Double_t tylength = 0.015;
    Double_t txlength = 0.025;
    for(Int_t i=0;i<NPad;i++) {
      char name[16];
      sprintf(name,"pad_%i",i);
      pad[i] = (TPad*) gROOT->FindObject(name);
      pad[i]->SetFrameLineWidth(2);  
      pad[i]->SetTickx(1);
      pad[i]->SetTicky(1);

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

    hFrame[1]->SetBins(x1Nbin,x1Min,x1Max);
    hFrame[1]->GetYaxis()->SetRangeUser(p1Min,p1Max);
    hFrame[1]->GetYaxis()->SetTitle("p_{z} [GeV/c]");
    hFrame[1]->Draw("axis");
    
    // 2D histogram z range
    Double_t dmax = hP1X1->GetMaximum();
    Double_t dmin = 0.0;
    hP1X1->GetZaxis()->SetRangeUser(1.1*dmin,dmax);

    hP1X1->GetZaxis()->CenterTitle();
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
    
    hP1X1->Draw(drawopts);


    TLine lZmean(zmean,p1Min,zmean,p1Max);
    lZmean.SetLineColor(kGray+2);
    lZmean.SetLineStyle(2);
    lZmean.Draw();

    TLine lPmean(x1Min,pzmean,x1Max,pzmean);
    lPmean.SetLineColor(kGray+2);
    lPmean.SetLineStyle(2);
    lPmean.Draw();

    
    gP1left->SetLineWidth(2);
    gP1left->Draw("F");
    gP1left->Draw("L");


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
      textMom->Draw();
    } else {
      textInfo->Draw();
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
    textLabel[1]->Draw();

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
    Leg=new TLegend(0.55,0.75,1 - gPad->GetRightMargin() - 0.02,0.95);

    PGlobals::SetPaveStyle(Leg);
    Leg->SetTextAlign(12);
    Leg->SetTextColor(kGray+3);
    Leg->SetTextFont(42);
    Leg->SetLineColor(1);
    Leg->SetBorderSize(0);
    Leg->SetFillColor(0);
    Leg->SetFillStyle(1001);
    Leg->SetFillStyle(0); // Hollow


    // if(autoscale) { 
    // 	Leg->AddEntry(hX1  ,"Current [10#times kA]","L");
    // } else {
    // 	Leg->AddEntry(hX1  ,"Current [kA]","L");	
    // }
    if(ifile.Contains("MCELLSTART"))
      Leg->AddEntry(hX1  ,"Current [0.1 kA]","L");
    else
      Leg->AddEntry(hX1  ,"Current [kA]","L");
    // Leg->AddEntry(gErms,"Energy spread (GeV)","PL");
    Leg->AddEntry(gErms,"Energy spread [0.1%]","PL");
    Leg->AddEntry(gemitX,Form("Emittance_{n,x} [%s]",emitSUnit.c_str()),"PL");
    Leg->AddEntry(gemitY,Form("Emittance_{n,y} [%s]",emitSUnit.c_str()),"PL");
    //Leg->AddEntry(gYrms,"Bunch width [#mum]","PL");


    hFrame[0]->SetBins(x1Nbin,x1Min,x1Max);
    hFrame[0]->GetYaxis()->SetRangeUser(0.0,1.1*yMax);
    hFrame[0]->GetXaxis()->SetTitle("#zeta [#mum]");
    hFrame[0]->Draw("axis");

    hX1->GetYaxis()->SetNdivisions(503);
    hX1->SetLineWidth(2);
    hX1->SetFillStyle(1001);
    hX1->SetFillColor(PGlobals::elecFill);
    // hX1->SetLineColor(kBlue);

    //hX1->Smooth();
    //    hX1->Draw("FL same");
    hX1->Draw("hist LF2 same");
    //hX1->Draw("C");

    TLine lZmean2(zmean,0.0,zmean,1.1*yMax);
    lZmean2.SetLineColor(kGray+2);
    lZmean2.SetLineStyle(2);
    lZmean2.Draw();

    Int_t markerSize = 1; 
    Int_t lineWidth  = 2;   

    gYrms->SetMarkerStyle(20);
    gYrms->SetLineStyle(1);
    gYrms->SetMarkerColor(kGray+1);
    gYrms->SetMarkerSize(markerSize); 
    gYrms->SetLineColor(kGray+1);
    gYrms->SetLineWidth(lineWidth);
    // gYrms->Draw("PL");

    // hP2X2sl[0]->Draw("colz");
    gemitY->SetMarkerStyle(20);
    gemitY->SetMarkerColor(kGray+2);
    gemitY->SetMarkerSize(markerSize);
    gemitY->SetLineWidth(lineWidth);
    gemitY->SetLineColor(kGray+2);
    gemitY->Draw("PL");

    gemitX->SetMarkerStyle(20);
    gemitX->SetMarkerColor(kGray+3);
    gemitX->SetMarkerSize(markerSize);
    gemitX->SetLineWidth(lineWidth);
    gemitX->SetLineColor(kGray+3);
    gemitX->Draw("PL");


    gErms->SetMarkerStyle(20);
    gErms->SetMarkerSize(markerSize);
    gErms->SetMarkerColor(kOrange+10);
    gErms->SetLineColor(kOrange+10);
    gErms->SetLineWidth(lineWidth);
    gErms->Draw("PL");
    
    if(opt.Contains("bw")) {
      gErms->SetMarkerStyle(21);
      gErms->SetMarkerSize(markerSize-0.2);
      gErms->SetLineWidth(1);
      
      gemitY->SetMarkerStyle(22);
      gemitY->SetLineWidth(1);
      gemitX->SetMarkerStyle(23);
      gemitX->SetLineWidth(1);
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
    textLabel[0]->Draw();

    lFrame = new TBox(gPad->GetUxmin(), gPad->GetUymin(),
		      gPad->GetUxmax(), gPad->GetUymax());
    lFrame->SetFillStyle(0);
    lFrame->SetLineColor(PGlobals::frameColor);
    lFrame->SetLineWidth(PGlobals::frameWidth);
    lFrame->Draw();
    
    gPad->RedrawAxis(); 

    gPad->Update();


    // Print to file --------------------------------------

    C->cd();

    // Print to a file
    // Output file
   
    TString fOutName = ifile;
    fOutName.ReplaceAll(".txt","");
    fOutName.Append("-p1x1");

    PGlobals::imgconv(C,fOutName,opt);

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

    const Int_t NPad1 = 2;
    lMargin = 0.15;
    rMargin = 0.18;
    bMargin = 0.15;
    tMargin = 0.04;
    factor = 1.0;  
    txoffset = 2.0;  
    PGlobals::CanvasAsymPartition(C1,NPad1,lMargin,rMargin,bMargin,tMargin,factor,0.028);

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

    hFrame[0]->SetBins(x1Nbin,x1Min,x1Max,x2Nbin,x2Min,x2Max);
    hFrame[0]->GetYaxis()->SetTitle(hX2X1->GetYaxis()->GetTitle());
    hFrame[0]->GetXaxis()->SetTitle(hX2X1->GetXaxis()->GetTitle());
    hFrame[0]->Draw("axis");
    
    TH2F *hX2X1cl = (TH2F*) hX2X1->Clone("hX2X1cl");
    Float_t xFactor = pad[0]->GetAbsWNDC()/pad[0]->GetAbsWNDC();
    Float_t yFactor = pad[0]->GetAbsHNDC()/pad[0]->GetAbsHNDC();
    
    hX2X1cl->GetZaxis()->CenterTitle();
    hX2X1cl->GetZaxis()->SetTitleFont(fonttype);
    hX2X1cl->GetZaxis()->SetTitleOffset(tyoffset);
    hX2X1cl->GetZaxis()->SetTitleSize(tfontsize);
    hX2X1cl->GetZaxis()->SetLabelFont(fonttype);
    hX2X1cl->GetZaxis()->SetLabelSize(fontsize);
    hX2X1cl->GetZaxis()->SetTickLength(xFactor*tylength/yFactor);
    hX2X1cl->Draw(drawopts);

    TLine lX1mean(zmean,x2Min,zmean,x2Max);
    lX1mean.SetLineColor(kGray+2);
    lX1mean.SetLineStyle(2);
    lX1mean.Draw();
    
    TLine lX2mean(x1Min,x_mean,x1Max,x_mean);
    lX2mean.SetLineColor(kGray+2);
    lX2mean.SetLineStyle(2);
    lX2mean.Draw();
    
    gPad->Update();
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

    y1 = gPad->GetBottomMargin();
    y2 = 1 - gPad->GetTopMargin();
    x1 = gPad->GetLeftMargin();
    x2 = 1 - gPad->GetRightMargin();
    yrange = y2-y1; 
    xrange = x2-x1; 

    TPaveText *textInfoX2X1 = new TPaveText(x1+0.02*xrange,y2-0.40*yrange,x1+0.20*xrange,y2-0.05*yrange,"NDC");
    PGlobals::SetPaveTextStyle(textInfoX2X1,12); 
    textInfoX2X1->SetTextColor(kGray+3);
    textInfoX2X1->SetTextFont(42);

    char text[64];
    sprintf(text,"Q = %5.1f pC",Charge);
    textInfoX2X1->AddText(text);
    sprintf(text,"#Delta#zeta = %5.2f #mum",zrms);
    textInfoX2X1->AddText(text);
    sprintf(text,"#Deltax = %5.2f #mum",x_rms);
    textInfoX2X1->AddText(text);
    sprintf(text,"#varepsilon_{x} = %5.2f #mum",emitx);
    textInfoX2X1->AddText(text);
    sprintf(text,"#beta_{x} = %5.2f mm",1E-3*betax);
    textInfoX2X1->AddText(text);
    textInfoX2X1->Draw();

    lFrame = new TBox(gPad->GetUxmin(), gPad->GetUymin(),
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

    hFrame[1]->SetBins(x1Nbin,x1Min,x1Max,x3Nbin,x3Min,x3Max);
    hFrame[1]->GetYaxis()->SetTitle(hX3X1->GetYaxis()->GetTitle());
    hFrame[1]->GetXaxis()->SetTitle(hX3X1->GetXaxis()->GetTitle());
    hFrame[1]->Draw("axis");
    
    TH2F *hX3X1cl = (TH2F*) hX3X1->Clone("hX3X1cl");
    xFactor = pad[0]->GetAbsWNDC()/pad[1]->GetAbsWNDC();
    yFactor = pad[0]->GetAbsHNDC()/pad[1]->GetAbsHNDC();
    hX3X1cl->GetZaxis()->CenterTitle();
    hX3X1cl->GetZaxis()->SetTitleFont(fonttype);
    hX3X1cl->GetZaxis()->SetTitleOffset(tyoffset);
    hX3X1cl->GetZaxis()->SetTitleSize(tfontsize);
    hX3X1cl->GetZaxis()->SetLabelFont(fonttype);
    hX3X1cl->GetZaxis()->SetLabelSize(fontsize);
    hX3X1cl->GetZaxis()->SetTickLength(xFactor*tylength/yFactor);
    hX3X1cl->Draw(drawopts);

    TLine lX1mean2(zmean,x3Min,zmean,x3Max);
    lX1mean2.SetLineColor(kGray+2);
    lX1mean2.SetLineStyle(2);
    lX1mean2.Draw();
    
    TLine lX3mean(x1Min,y_mean,x1Max,y_mean);
    lX3mean.SetLineColor(kGray+2);
    lX3mean.SetLineStyle(2);
    lX3mean.Draw();

    
    gPad->Update();
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

    y1 = gPad->GetBottomMargin();
    y2 = 1 - gPad->GetTopMargin();
    x1 = gPad->GetLeftMargin();
    x2 = 1 - gPad->GetRightMargin();
    yrange = y2-y1; 
    xrange = x2-x1; 

    
    TPaveText *textInfoX3X1 =  new TPaveText(x1+0.02*xrange,y2-0.40*yrange,x1+0.20*xrange,y2-0.05*yrange,"NDC");
    PGlobals::SetPaveTextStyle(textInfoX3X1,12); 
    textInfoX3X1->SetTextColor(kGray+3);
    textInfoX3X1->SetTextFont(42);

    sprintf(text,"Q = %5.1f pC",Charge);
    textInfoX3X1->AddText(text);
    sprintf(text,"#Delta#zeta = %5.2f #mum",zrms);
    textInfoX3X1->AddText(text);
    sprintf(text,"#Deltay = %5.2f #mum",y_rms);
    textInfoX3X1->AddText(text);
    sprintf(text,"#varepsilon_{y} = %5.2f %s",emity,emitSUnit.c_str());
    textInfoX3X1->AddText(text);
    sprintf(text,"#beta_{y} = %5.2f mm",betay/1000);
    textInfoX3X1->AddText(text);
    
    // sprintf(text,"#varepsilon_{x} = %5.2f #mum",emitx);
    // textInfoX3X1->AddText(text);
    // sprintf(text,"#beta_{x} = %5.2f mm",1E-3*betax);
    // textInfoX3X1->AddText(text);
    textInfoX3X1->Draw();

    lFrame = new TBox(gPad->GetUxmin(), gPad->GetUymin(),
		      gPad->GetUxmax(), gPad->GetUymax());
    lFrame->SetFillStyle(0);
    lFrame->SetLineColor(PGlobals::frameColor);
    lFrame->SetLineWidth(PGlobals::frameWidth);
    lFrame->Draw();

    gPad->RedrawAxis(); 

    
    C1->cd();
    
    TString fOutName1 = ifile;
    fOutName1.ReplaceAll(".txt","");
    fOutName1.Append("-x2x3x1");
    PGlobals::imgconv(C1,fOutName1,opt);


    // Transverse phasespace

    // x3-x2
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

    hFrame[0]->SetBins(x2Nbin,x2Min,x2Max,x3Nbin,x3Min,x3Max);
    hFrame[0]->GetYaxis()->SetTitle(hX3X2->GetYaxis()->GetTitle());
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

    hX3X2cl->Draw(drawopts);
        
    TLine lX2mean2(x_mean,x3Min,x_mean,x3Max);
    lX2mean2.SetLineColor(kGray+2);
    lX2mean2.SetLineStyle(2);
    lX2mean2.Draw();

    TLine lX3mean2(x2Min,y_mean,x2Max,y_mean);
    lX3mean2.SetLineColor(kGray+2);
    lX3mean2.SetLineStyle(2);
    lX3mean2.Draw();

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

    sprintf(text,"Q = %5.1f pC",Charge);
    textStatX3X2->AddText(text);
    sprintf(text,"#Deltax = %5.2f #mum",x_rms);
    textStatX3X2->AddText(text);
    sprintf(text,"#Deltay = %5.2f #mum",y_rms);
    textStatX3X2->AddText(text);
    sprintf(text,"#varepsilon_{n,x} = %5.2f %s",emitx,emitSUnit.c_str());
    textStatX3X2->AddText(text);
    sprintf(text,"#varepsilon_{n,y} = %5.2f %s",emity,emitSUnit.c_str());
    textStatX3X2->AddText(text);
    sprintf(text,"#beta_{x} = %5.2f mm",betax/1000);
    textStatX3X2->AddText(text);
    sprintf(text,"#beta_{y} = %5.2f mm",betay/1000);
    textStatX3X2->AddText(text);

    textStatX3X2->Draw();

    lFrame = new TBox(gPad->GetUxmin(), gPad->GetUymin(),
		      gPad->GetUxmax(), gPad->GetUymax());
    lFrame->SetFillStyle(0);
    lFrame->SetLineColor(PGlobals::frameColor);
    lFrame->SetLineWidth(PGlobals::frameWidth);
    lFrame->Draw();

    gPad->RedrawAxis(); 
    
    TString fOutNameX3X2 = ifile;
    fOutNameX3X2.ReplaceAll(".txt","");
    fOutNameX3X2.Append("-x3x2");
    PGlobals::imgconv(CX3X2,fOutNameX3X2,opt);


    // p2-x2
    sprintf(cName,"C2");     
    TCanvas *C2 = (TCanvas*) gROOT->FindObject(cName);
    if(C2==NULL) C2 = new TCanvas("C2","Transverse phasespace p2-x2",sizex,sizey);
    C2->cd();
    C2->Clear();

    lMargin = 0.15;
    rMargin = 0.18;
    bMargin = 0.15;
    tMargin = 0.04;
    factor = 1.0;  
    txoffset = 1.2;  

    // Format for y axis
    hP2X2->GetYaxis()->SetTitleFont(fonttype);
    hP2X2->GetYaxis()->SetTitleSize(tfontsize);
    hP2X2->GetYaxis()->SetTitleOffset(tyoffset);
    hP2X2->GetYaxis()->SetLabelFont(fonttype);
    hP2X2->GetYaxis()->SetLabelSize(fontsize);
    hP2X2->GetYaxis()->SetLabelOffset(lyoffset);
    hP2X2->GetYaxis()->SetTickLength(xFactor*tylength/yFactor);
    hP2X2->GetYaxis()->CenterTitle();
    
    // Format for x axis
    hP2X2->GetXaxis()->SetTitleFont(fonttype);
    hP2X2->GetXaxis()->SetTitleSize(tfontsize+2);
    hP2X2->GetXaxis()->SetTitleOffset(txoffset);
    hP2X2->GetXaxis()->SetLabelFont(fonttype);
    hP2X2->GetXaxis()->SetLabelSize(fontsize+2);
    hP2X2->GetXaxis()->SetLabelOffset(lxoffset);
    hP2X2->GetXaxis()->CenterTitle();
    hP2X2->GetXaxis()->SetTickLength(yFactor*txlength/xFactor);

    hP2X2->GetZaxis()->CenterTitle();
    hP2X2->GetZaxis()->SetTitleFont(fonttype);
    hP2X2->GetZaxis()->SetTitleOffset(tyoffset);
    hP2X2->GetZaxis()->SetTitleSize(tfontsize);
    hP2X2->GetZaxis()->SetLabelFont(fonttype);
    hP2X2->GetZaxis()->SetLabelSize(fontsize);
    if(opt.Contains("logz")) 
      hP2X2->GetZaxis()->SetLabelOffset(0);
    else
      hP2X2->GetZaxis()->SetLabelOffset(lyoffset);
    hP2X2->GetZaxis()->SetTickLength(0.01);      


    TH2F *hFrP2X2 = (TH2F*) hP2X2->Clone("hFrP2X2");
    hFrP2X2->Reset();
    hFrP2X2->SetBins(x2Nbin,x2Min,x2Max,p2Nbin,p2Min,p2Max);

    C2->cd(0);

    if(opt.Contains("logz")) {
      gPad->SetLogz(1);
    } else {
      gPad->SetLogz(0);
    }

    hFrP2X2->Draw("axis");
    hP2X2->Draw(drawopt + "same");

    TLine lXmean(x_mean,p2Min,x_mean,p2Max);
    lXmean.SetLineColor(kGray+2);
    lXmean.SetLineStyle(2);
    lXmean.Draw();

    TLine lPxmean(x2Min,pxmean,x2Max,pxmean);
    lPxmean.SetLineColor(kGray+2);
    lPxmean.SetLineStyle(2);
    lPxmean.Draw();

    //    hP2X2->Draw(drawopt);
    
    if(opt.Contains("elli") ) 
      for(Int_t k=0;k<SNbin;k++)
	if(sellipP2X2[k] != NULL)
	  sellipP2X2[k]->Draw();

    gPad->Update();
    palette = (TPaletteAxis*) hP2X2->GetListOfFunctions()->FindObject("palette");
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


    TPaveText *textStatInt = new TPaveText(x1+0.02,y2-0.40,x1+0.20,y2-0.05,"NDC");
    PGlobals::SetPaveTextStyle(textStatInt,12); 
    textStatInt->SetTextColor(kGray+3);
    textStatInt->SetTextFont(42);

    sprintf(text,"Q = %5.1f pC",Charge);
    textStatInt->AddText(text);
    sprintf(text,"#Deltax = %5.2f #mum",x_rms);
    textStatInt->AddText(text);
    sprintf(text,"#Deltap_{x}/mc = %5.2f",pxrms);
    textStatInt->AddText(text);
    sprintf(text,"#varepsilon_{n,x} = %5.2f %s",emitx,emitSUnit.c_str());
    textStatInt->AddText(text);
    sprintf(text,"#beta_{x} = %5.2f mm",1E-3*betax);
    textStatInt->AddText(text);
    textStatInt->Draw();

    lFrame = new TBox(gPad->GetUxmin(), gPad->GetUymin(),
		      gPad->GetUxmax(), gPad->GetUymax());
    lFrame->SetFillStyle(0);
    lFrame->SetLineColor(PGlobals::frameColor);
    lFrame->SetLineWidth(PGlobals::frameWidth);
    lFrame->Draw();
    
    gPad->RedrawAxis(); 
    
    TString fOutName2 = ifile;
    fOutName2.ReplaceAll(".txt","");
    fOutName2.Append("-p2x2");
    PGlobals::imgconv(C2,fOutName2,opt);

    // p3-x3
    sprintf(cName,"C3");     
    TCanvas *C3 = (TCanvas*) gROOT->FindObject(cName);
    if(C3==NULL) C3 = new TCanvas("C3","Transverse phasespace p3-x3",sizex,sizey);
    C3->cd();
    C3->Clear();

    // Format for y axis
    hP3X3->GetYaxis()->SetTitleFont(fonttype);
    hP3X3->GetYaxis()->SetTitleSize(tfontsize);
    hP3X3->GetYaxis()->SetTitleOffset(tyoffset);
    hP3X3->GetYaxis()->SetLabelFont(fonttype);
    hP3X3->GetYaxis()->SetLabelSize(fontsize);
    hP3X3->GetYaxis()->SetLabelOffset(lyoffset);
    hP3X3->GetYaxis()->SetTickLength(xFactor*tylength/yFactor);
    hP3X3->GetYaxis()->CenterTitle();

    // Format for x axis
    hP3X3->GetXaxis()->SetTitleFont(fonttype);
    hP3X3->GetXaxis()->SetTitleSize(tfontsize+2);
    hP3X3->GetXaxis()->SetTitleOffset(txoffset);
    hP3X3->GetXaxis()->SetLabelFont(fonttype);
    hP3X3->GetXaxis()->SetLabelSize(fontsize+2);
    hP3X3->GetXaxis()->SetLabelOffset(lxoffset);
    hP3X3->GetXaxis()->CenterTitle();
    hP3X3->GetXaxis()->SetTickLength(yFactor*txlength/xFactor);

    hP3X3->GetZaxis()->CenterTitle();
    hP3X3->GetZaxis()->SetTitleFont(fonttype);
    hP3X3->GetZaxis()->SetTitleOffset(tyoffset);
    hP3X3->GetZaxis()->SetTitleSize(tfontsize);
    hP3X3->GetZaxis()->SetLabelFont(fonttype);
    hP3X3->GetZaxis()->SetLabelSize(fontsize);
    if(opt.Contains("logz")) 
      hP3X3->GetZaxis()->SetLabelOffset(0);
    else
      hP3X3->GetZaxis()->SetLabelOffset(lyoffset);
    hP3X3->GetZaxis()->SetTickLength(0.01);      


    TH2F *hFrP3X3 = (TH2F*) hP3X3->Clone("hFrP3X3");
    hFrP3X3->Reset();
    hFrP3X3->SetBins(x3Nbin,x3Min,x3Max,p3Nbin,p3Min,p3Max);
    
    C3->cd(0);

    if(opt.Contains("logz")) {
      gPad->SetLogz(1);
    } else {
      gPad->SetLogz(0);
    }

    hFrP3X3->Draw("axis");
    hP3X3->Draw(drawopt + "same");
    
    TLine lYmean(y_mean,p3Min,y_mean,p3Max);
    lYmean.SetLineColor(kGray+2);
    lYmean.SetLineStyle(2);
    lYmean.Draw();

    TLine lPymean(x3Min,pymean,x3Max,pymean);
    lPymean.SetLineColor(kGray+2);
    lPymean.SetLineStyle(2);
    lPymean.Draw();

    //    hP3X3->Draw(drawopts);
        
    if(opt.Contains("elli") ) 
      for(Int_t k=0;k<SNbin;k++)
	if(sellipP3X3[k] != NULL)
	  sellipP3X3[k]->Draw();

    gPad->Update();
    palette = (TPaletteAxis*) hP3X3->GetListOfFunctions()->FindObject("palette");
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

    
    textStatInt = new TPaveText(x1+0.02,y2-0.40,x1+0.20,y2-0.05,"NDC");
    PGlobals::SetPaveTextStyle(textStatInt,12); 
    textStatInt->SetTextColor(kGray+3);
    textStatInt->SetTextFont(42);

    sprintf(text,"Q = %5.1f pC",Charge);
    textStatInt->AddText(text);
    sprintf(text,"#Deltay = %5.2f #mum",y_rms);
    textStatInt->AddText(text);
    sprintf(text,"#Deltap_{y}/mc = %5.2f",pyrms);
    textStatInt->AddText(text);
    sprintf(text,"#varepsilon_{n,y} = %5.2f %s",emity,emitSUnit.c_str());
    textStatInt->AddText(text);
    sprintf(text,"#beta_{y} = %5.2f mm",1E-3*betay);
    textStatInt->AddText(text);
    textStatInt->Draw();
    
    lFrame = new TBox(gPad->GetUxmin(), gPad->GetUymin(),
		      gPad->GetUxmax(), gPad->GetUymax());
    lFrame->SetFillStyle(0);
    lFrame->SetLineColor(PGlobals::frameColor);
    lFrame->SetLineWidth(PGlobals::frameWidth);
    lFrame->Draw();
    
    gPad->RedrawAxis(); 

    TString fOutName3 = ifile;
    fOutName3.ReplaceAll(".txt","");
    fOutName3.Append("-p3x3");
    PGlobals::imgconv(C3,fOutName3,opt);


    // p2-x1
    sprintf(cName,"C4");     
    TCanvas *C4 = (TCanvas*) gROOT->FindObject(cName);
    if(C4==NULL) C4 = new TCanvas("C4","p2-x1",sizex,sizey);
    C4->cd();
    C4->Clear();

    lMargin = 0.15;
    rMargin = 0.18;
    bMargin = 0.15;
    tMargin = 0.04;
    factor = 1.0;  
    txoffset = 1.2;  

    // Format for y axis
    hP2X1->GetYaxis()->SetTitleFont(fonttype);
    hP2X1->GetYaxis()->SetTitleSize(tfontsize);
    hP2X1->GetYaxis()->SetTitleOffset(tyoffset);
    hP2X1->GetYaxis()->SetLabelFont(fonttype);
    hP2X1->GetYaxis()->SetLabelSize(fontsize);
    hP2X1->GetYaxis()->SetLabelOffset(lyoffset);
    hP2X1->GetYaxis()->SetTickLength(xFactor*tylength/yFactor);
    hP2X1->GetYaxis()->CenterTitle();
    
    // Format for x axis
    hP2X1->GetXaxis()->SetTitleFont(fonttype);
    hP2X1->GetXaxis()->SetTitleSize(tfontsize+2);
    hP2X1->GetXaxis()->SetTitleOffset(txoffset);
    hP2X1->GetXaxis()->SetLabelFont(fonttype);
    hP2X1->GetXaxis()->SetLabelSize(fontsize+2);
    hP2X1->GetXaxis()->SetLabelOffset(lxoffset);
    hP2X1->GetXaxis()->CenterTitle();
    hP2X1->GetXaxis()->SetTickLength(yFactor*txlength/xFactor);

    hP2X1->GetZaxis()->CenterTitle();
    hP2X1->GetZaxis()->SetTitleFont(fonttype);
    hP2X1->GetZaxis()->SetTitleOffset(tyoffset);
    hP2X1->GetZaxis()->SetTitleSize(tfontsize);
    hP2X1->GetZaxis()->SetLabelFont(fonttype);
    hP2X1->GetZaxis()->SetLabelSize(fontsize);
    if(opt.Contains("logz")) 
      hP2X1->GetZaxis()->SetLabelOffset(0);
    else
      hP2X1->GetZaxis()->SetLabelOffset(lyoffset);
    hP2X1->GetZaxis()->SetTickLength(0.01);      


    TH2F *hFrP2X1 = (TH2F*) hP2X1->Clone("hFrP2X1");
    hFrP2X1->Reset();
    hFrP2X1->SetBins(x1Nbin,x1Min,x1Max,p2Nbin,p2Min,p2Max);


    C4->cd(0);

    if(opt.Contains("logz")) {
      gPad->SetLogz(1);
    } else {
      gPad->SetLogz(0);
    }

    hFrP2X1->Draw("axis");
    hP2X1->Draw(drawopt + "same");

    TLine lZmean3(zmean,p2Min,zmean,p2Max);
    lZmean3.SetLineColor(kGray+2);
    lZmean3.SetLineStyle(2);
    lZmean3.Draw();
    
    TLine lPxmean2(x1Min,pxmean,x1Max,pxmean);
    lPxmean2.SetLineColor(kGray+2);
    lPxmean2.SetLineStyle(2);
    lPxmean2.Draw();

    gPad->Update();
    palette = (TPaletteAxis*) hP2X1->GetListOfFunctions()->FindObject("palette");
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


    textStatInt = new TPaveText(x1+0.02,y2-0.40,x1+0.20,y2-0.05,"NDC");
    PGlobals::SetPaveTextStyle(textStatInt,12); 
    textStatInt->SetTextColor(kGray+3);
    textStatInt->SetTextFont(42);

    sprintf(text,"Q = %5.1f pC",Charge);
    textStatInt->AddText(text);
    sprintf(text,"#Delta#zeta = %5.2f #mum",zrms);
    textStatInt->AddText(text);
    sprintf(text,"#Deltap_{x}/mc = %5.2f",pxrms);
    textStatInt->AddText(text);
    sprintf(text,"#varepsilon_{n,x} = %5.2f %s",emitx,emitSUnit.c_str());
    textStatInt->AddText(text);
    sprintf(text,"#beta_{x} = %5.2f mm",1E-3*betax);
    textStatInt->AddText(text);
    textStatInt->Draw();

    lFrame = new TBox(gPad->GetUxmin(), gPad->GetUymin(),
		      gPad->GetUxmax(), gPad->GetUymax());
    lFrame->SetFillStyle(0);
    lFrame->SetLineColor(PGlobals::frameColor);
    lFrame->SetLineWidth(PGlobals::frameWidth);
    lFrame->Draw();
    
    gPad->RedrawAxis(); 

    
    TString fOutName4 = ifile;
    fOutName4.ReplaceAll(".txt","");
    fOutName4.Append("-p2x1");
    PGlobals::imgconv(C4,fOutName4,opt);
    

    // p3-x1
    sprintf(cName,"C5");     
    TCanvas *C5 = (TCanvas*) gROOT->FindObject(cName);
    if(C5==NULL) C5 = new TCanvas("C5","p2-x1",sizex,sizey);
    C5->cd();
    C5->Clear();

    lMargin = 0.15;
    rMargin = 0.18;
    bMargin = 0.15;
    tMargin = 0.04;
    factor = 1.0;  
    txoffset = 1.2;  

    // Format for y axis
    hP3X1->GetYaxis()->SetTitleFont(fonttype);
    hP3X1->GetYaxis()->SetTitleSize(tfontsize);
    hP3X1->GetYaxis()->SetTitleOffset(tyoffset);
    hP3X1->GetYaxis()->SetLabelFont(fonttype);
    hP3X1->GetYaxis()->SetLabelSize(fontsize);
    hP3X1->GetYaxis()->SetLabelOffset(lyoffset);
    hP3X1->GetYaxis()->SetTickLength(xFactor*tylength/yFactor);
    hP3X1->GetYaxis()->CenterTitle();
    
    // Format for x axis
    hP3X1->GetXaxis()->SetTitleFont(fonttype);
    hP3X1->GetXaxis()->SetTitleSize(tfontsize+2);
    hP3X1->GetXaxis()->SetTitleOffset(txoffset);
    hP3X1->GetXaxis()->SetLabelFont(fonttype);
    hP3X1->GetXaxis()->SetLabelSize(fontsize+2);
    hP3X1->GetXaxis()->SetLabelOffset(lxoffset);
    hP3X1->GetXaxis()->CenterTitle();
    hP3X1->GetXaxis()->SetTickLength(yFactor*txlength/xFactor);

    hP3X1->GetZaxis()->CenterTitle();
    hP3X1->GetZaxis()->SetTitleFont(fonttype);
    hP3X1->GetZaxis()->SetTitleOffset(tyoffset);
    hP3X1->GetZaxis()->SetTitleSize(tfontsize);
    hP3X1->GetZaxis()->SetLabelFont(fonttype);
    hP3X1->GetZaxis()->SetLabelSize(fontsize);
    if(opt.Contains("logz")) 
      hP3X1->GetZaxis()->SetLabelOffset(0);
    else
      hP3X1->GetZaxis()->SetLabelOffset(lyoffset);
    hP3X1->GetZaxis()->SetTickLength(0.01);      


    TH2F *hFrP3X1 = (TH2F*) hP3X1->Clone("hFrP3X1");
    hFrP3X1->Reset();
    hFrP3X1->SetBins(x1Nbin,x1Min,x1Max,p2Nbin,p2Min,p2Max);

    C5->cd(0);

    if(opt.Contains("logz")) {
      gPad->SetLogz(1);
    } else {
      gPad->SetLogz(0);
    }

    hFrP3X1->Draw("axis");
    hP3X1->Draw(drawopt + "same");

    TLine lZmean4(zmean,p3Min,zmean,p3Max);
    lZmean4.SetLineColor(kGray+2);
    lZmean4.SetLineStyle(2);
    lZmean4.Draw();
    
    TLine lPxmean3(x1Min,pymean,x1Max,pymean);
    lPxmean3.SetLineColor(kGray+2);
    lPxmean3.SetLineStyle(2);
    lPxmean3.Draw();

    gPad->Update();
    palette = (TPaletteAxis*) hP3X1->GetListOfFunctions()->FindObject("palette");
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


    textStatInt = new TPaveText(x1+0.02,y2-0.40,x1+0.20,y2-0.05,"NDC");
    PGlobals::SetPaveTextStyle(textStatInt,12); 
    textStatInt->SetTextColor(kGray+3);
    textStatInt->SetTextFont(42);

    sprintf(text,"Q = %5.1f pC",Charge);
    textStatInt->AddText(text);
    sprintf(text,"#Delta#zeta = %5.2f #mum",zrms);
    textStatInt->AddText(text);
    sprintf(text,"#Deltap_{y}/mc = %5.2f",pyrms);
    textStatInt->AddText(text);
    sprintf(text,"#varepsilon_{n,y} = %5.2f %s",emity,emitSUnit.c_str());
    textStatInt->AddText(text);
    sprintf(text,"#beta_{y} = %5.2f mm",1E-3*betay);
    textStatInt->AddText(text);
    textStatInt->Draw();

    lFrame = new TBox(gPad->GetUxmin(), gPad->GetUymin(),
		      gPad->GetUxmax(), gPad->GetUymax());
    lFrame->SetFillStyle(0);
    lFrame->SetLineColor(PGlobals::frameColor);
    lFrame->SetLineWidth(PGlobals::frameWidth);
    lFrame->Draw();
    
    gPad->RedrawAxis(); 

    
    TString fOutName5 = ifile;
    fOutName5.ReplaceAll(".txt","");
    fOutName5.Append("-p3x1");
    PGlobals::imgconv(C5,fOutName5,opt);

    // Slices plots
    {
      gStyle->cd();

      sprintf(cName,"CA4");     
      TCanvas *CA4 = (TCanvas*) gROOT->FindObject(cName);
      if(CA4==NULL) CA4 = new TCanvas("CA4","Sliced p2-x2",600,800);
      CA4->cd();
      CA4->Clear();

      Int_t ndiv = 4;
      CA4->Divide(1,ndiv);

      TString fOutName4 = ifile;
      fOutName4.ReplaceAll(".txt","");
      fOutName4.Append("-slp2x2");

      CA4->Print(fOutName4 + ".ps[","Portrait");

      Int_t fontsize = 18;
      Int_t tfontsize = 20;
      Double_t txoffset = 4.0;
      Double_t lxoffset = 0.02;
      Double_t tyoffset = 2.2;
      Double_t lyoffset = 0.01;
      Double_t tzoffset = 3.0;
      Double_t lzoffset = 0.01;
      Double_t tylength = 0.010;
      Double_t txlength = 0.050;

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

      hFrP2X2 = (TH2F*) hP2X2->Clone("hFrP2X2");
      hFrP2X2->Reset();
      hFrP2X2->SetBins(x2Nbin,x2Min,x2Max,p2Nbin,p2Min,p2Max);

      CA4->cd(1);
      gPad->SetTickx(1);
      gPad->SetTicky(1);
      gPad->SetTopMargin(0.02);
      gPad->SetBottomMargin(0.22);

      hFrP2X2->Draw("axis");
      hP2X2->Draw(drawopt + "same");
      
      lXmean.Draw();
      lPxmean.Draw();
      
      // hP2X2->Draw(drawopts);
      //    ellipP2X2->Draw();
      
      if(opt.Contains("elli") ) 
	for(Int_t k=0;k<SNbin;k++)
	  if(sellipP2X2[k] != NULL)
	    sellipP2X2[k]->Draw();

      gPad->Update();
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

      textStatInt = new TPaveText(x1+0.02,y2-0.40,x1+0.20,y2-0.05,"NDC");
      PGlobals::SetPaveTextStyle(textStatInt,12); 
      textStatInt->SetTextColor(kGray+3);
      textStatInt->SetTextFont(42);

      sprintf(text,"Q = %5.1f pC",Charge);
      textStatInt->AddText(text);
      sprintf(text,"#Deltax = %5.2f #mum",x_rms);
      textStatInt->AddText(text);
      sprintf(text,"#Deltap_{x}/mc = %5.2f",pxrms);
      textStatInt->AddText(text);
      sprintf(text,"#varepsilon_{n,x} = %5.2f %s",emitx,emitSUnit.c_str());
      textStatInt->AddText(text);
      sprintf(text,"#beta_{x} = %5.2f mm",1E-3*betax);
      textStatInt->AddText(text);
      textStatInt->Draw();

      lFrame = new TBox(gPad->GetUxmin(), gPad->GetUymin(),
			gPad->GetUxmax(), gPad->GetUymax());
      lFrame->SetFillStyle(0);
      lFrame->SetLineColor(PGlobals::frameColor);
      lFrame->SetLineWidth(PGlobals::frameWidth);
      lFrame->Draw();
      
      gPad->RedrawAxis(); 
      
      
      TPaveText **textStat = new TPaveText*[SNbin];
      TPaveText **textSlice = new TPaveText*[SNbin];
      TLine **slXmean = new TLine*[SNbin];
      TLine **slPxmean = new TLine*[SNbin];
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
	gPad->SetTickx(1);
	gPad->SetTicky(1);

     	gPad->SetTopMargin(0.02);
	gPad->SetBottomMargin(0.22);

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

	hFrP2X2->Draw("axis");
	hP2X2sl[k]->Draw(drawopt + "same");

	slXmean[k] = new TLine(hP2X2sl[k]->GetMean(1),p2Min,hP2X2sl[k]->GetMean(1),p2Max);
	slXmean[k]->SetLineColor(kGray+2);
	slXmean[k]->SetLineStyle(2);
	slXmean[k]->Draw();

	
	slPxmean[k] = new TLine(x2Min,hP2X2sl[k]->GetMean(2),x2Max,hP2X2sl[k]->GetMean(2));
	slPxmean[k]->SetLineColor(kGray+2);
	slPxmean[k]->SetLineStyle(2);
	slPxmean[k]->Draw();
	
	if(opt.Contains("elli") ) 
	  if(sellipP2X2[k] != NULL)
	    sellipP2X2[k]->Draw();

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

	Double_t y1 = gPad->GetBottomMargin();
	Double_t y2 = 1 - gPad->GetTopMargin();
	Double_t x1 = gPad->GetLeftMargin();
	Double_t x2 = 1 - gPad->GetRightMargin();
      
	textStat[k] = new TPaveText(x1+0.02,y2-0.40,x1+0.20,y2-0.05,"NDC");
	PGlobals::SetPaveTextStyle(textStat[k],12); 
	textStat[k]->SetTextColor(kGray+3);
	textStat[k]->SetTextFont(42);

	char text[64];
	sprintf(text,"Q = %5.1f pC",sCharge[k]);
	textStat[k]->AddText(text);
	sprintf(text,"#Deltax = %5.2f #mum",sx_rms[k]);
	textStat[k]->AddText(text);
	sprintf(text,"#Deltap_{x}/mc = %5.2f",spx_rms[k]);
	textStat[k]->AddText(text);
	sprintf(text,"#varepsilon_{n,x} = %5.2f %s",semitx[k],emitSUnit.c_str());	
	textStat[k]->AddText(text);
	sprintf(text,"#beta_{x} = %5.2f mm",1E-3*sbetax[k]);
	textStat[k]->AddText(text);
	textStat[k]->Draw();


	textSlice[k] = new TPaveText(x1+0.4,y2-0.2,x2-0.02,y2-0.0,"NDC");
	PGlobals::SetPaveTextStyle(textSlice[k],32); 
	textSlice[k]->SetTextColor(kBlack);
	textSlice[k]->SetTextFont(42);

	sprintf(text,"%5.2f #mum < #zeta < %5.2f #mum",sBinLim[k],sBinLim[k+1]);
	textSlice[k]->AddText(text);
	textSlice[k]->Draw();

	lFrame = new TBox(gPad->GetUxmin(), gPad->GetUymin(),
			  gPad->GetUxmax(), gPad->GetUymax());
	lFrame->SetFillStyle(0);
	lFrame->SetLineColor(PGlobals::frameColor);
	lFrame->SetLineWidth(PGlobals::frameWidth);
	lFrame->Draw();

	gPad->RedrawAxis(); 


	if(ic+1==ndiv) {
	  CA4->cd(0);
	  CA4->Print(fOutName4 + ".ps");
	}
      }

      CA4->Print(fOutName4 + ".ps]");

      gSystem->Exec("ps2pdf " + fOutName4 + ".ps " + fOutName4 + ".pdf");
      gSystem->Exec("rm -rf " + fOutName4 + ".ps"); 


      // ---------

      sprintf(cName,"CA5");     
      TCanvas *CA5 = (TCanvas*) gROOT->FindObject(cName);
      if(CA5==NULL) CA5 = new TCanvas("CA5","Sliced p3-x3",600,800);
      CA5->cd();
      CA5->Clear();

      ndiv = 4;
      CA5->Divide(1,ndiv);

      TString fOutName5 = ifile;
      fOutName5.ReplaceAll(".txt","");
      fOutName5.Append("-slp3x3");

      CA5->Print(fOutName5 + ".ps[","Portrait");

      hP3X3->GetXaxis()->SetLabelFont(fonttype);
      hP3X3->GetXaxis()->SetLabelSize(fontsize);
      hP3X3->GetXaxis()->SetTitleFont(fonttype);
      hP3X3->GetXaxis()->SetTitleSize(tfontsize);
      hP3X3->GetXaxis()->SetTitleOffset(txoffset);
      hP3X3->GetXaxis()->SetTickLength(txlength);
      hP3X3->GetXaxis()->CenterTitle();
	
      hP3X3->GetYaxis()->SetLabelFont(fonttype);
      hP3X3->GetYaxis()->SetLabelSize(fontsize);
      hP3X3->GetYaxis()->SetTitleFont(fonttype);
      hP3X3->GetYaxis()->SetTitleSize(tfontsize);
      hP3X3->GetYaxis()->SetTitleOffset(tyoffset);
      hP3X3->GetYaxis()->SetTickLength(tylength);
      hP3X3->GetYaxis()->CenterTitle();
	
      hP3X3->GetZaxis()->SetLabelFont(fonttype);
      hP3X3->GetZaxis()->SetLabelSize(fontsize);
      hP3X3->GetZaxis()->SetTitleFont(fonttype);
      hP3X3->GetZaxis()->SetTitleSize(tfontsize);
      hP3X3->GetZaxis()->SetTitleOffset(tzoffset);
      hP3X3->GetZaxis()->SetTickLength(tylength);
      hP3X3->GetZaxis()->CenterTitle();    

      hFrP3X3 = (TH2F*) hP3X3->Clone("hFrP3X3");
      hFrP3X3->Reset();
      hFrP3X3->SetBins(x3Nbin,x3Min,x3Max,p3Nbin,p3Min,p3Max);
      
      CA5->cd(1);
      gPad->SetTickx(1);
      gPad->SetTicky(1);
      gPad->SetTopMargin(0.02);
      gPad->SetBottomMargin(0.22);

      hFrP3X3->Draw("axis");
      hP3X3->Draw(drawopt + "same");

      lYmean.Draw();
      lPymean.Draw();
      
      //  hP3X3->Draw(drawopts);
      //    ellipP3X3->Draw();
      
      if(opt.Contains("elli") ) 
	for(Int_t k=0;k<SNbin;k++)
	  if(sellipP3X3[k] != NULL)
	    sellipP3X3[k]->Draw();

      gPad->Update();
      palette = (TPaletteAxis*)hP3X3->GetListOfFunctions()->FindObject("palette");
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

      textStatInt = new TPaveText(x1+0.02,y2-0.40,x1+0.20,y2-0.05,"NDC");
      PGlobals::SetPaveTextStyle(textStatInt,12); 
      textStatInt->SetTextColor(kGray+3);
      textStatInt->SetTextFont(42);

      sprintf(text,"Q = %5.1f pC",Charge);
      textStatInt->AddText(text);
      sprintf(text,"#Deltay = %5.2f #mum",y_rms);
      textStatInt->AddText(text);
      sprintf(text,"#Deltap_{y}/mc = %5.2f",pyrms);
      textStatInt->AddText(text);
      sprintf(text,"#varepsilon_{n,y} = %5.2f %s",emity,emitSUnit.c_str());
      textStatInt->AddText(text);
      sprintf(text,"#beta_{y} = %5.2f mm",1E-3*betay);
      textStatInt->AddText(text);
      textStatInt->Draw();

      lFrame = new TBox(gPad->GetUxmin(), gPad->GetUymin(),
			gPad->GetUxmax(), gPad->GetUymax());
      lFrame->SetFillStyle(0);
      lFrame->SetLineColor(PGlobals::frameColor);
      lFrame->SetLineWidth(PGlobals::frameWidth);
      lFrame->Draw();
      
      gPad->RedrawAxis(); 


      textStat = new TPaveText*[SNbin];
      textSlice = new TPaveText*[SNbin];
      TLine **slXmean2 = new TLine*[SNbin];
      TLine **slPxmean2 = new TLine*[SNbin];

      pnumber = 0;
      for(Int_t k=0;k<SNbin;k++) {
	pnumber++;
	Int_t ic = pnumber%ndiv;

	// new page
	if( ic==0 ) {
	  CA5->cd(0);
	  CA5->Clear();
	  CA5->Divide(1,ndiv);
	}
	CA5->cd(ic+1);
	gPad->SetTickx(1);
	gPad->SetTicky(1);
	gPad->SetTopMargin(0.02);
	gPad->SetBottomMargin(0.22);

	hP3X3sl[k]->GetXaxis()->SetLabelFont(fonttype);
	hP3X3sl[k]->GetXaxis()->SetLabelSize(fontsize);
	hP3X3sl[k]->GetXaxis()->SetTitleFont(fonttype);
	hP3X3sl[k]->GetXaxis()->SetTitleSize(tfontsize);
	hP3X3sl[k]->GetXaxis()->SetTitleOffset(txoffset);
	hP3X3sl[k]->GetXaxis()->SetTickLength(txlength);
	hP3X3sl[k]->GetXaxis()->CenterTitle();
	
	hP3X3sl[k]->GetYaxis()->SetLabelFont(fonttype);
	hP3X3sl[k]->GetYaxis()->SetLabelSize(fontsize);
	hP3X3sl[k]->GetYaxis()->SetTitleFont(fonttype);
	hP3X3sl[k]->GetYaxis()->SetTitleSize(tfontsize);
	hP3X3sl[k]->GetYaxis()->SetTitleOffset(tyoffset);
	hP3X3sl[k]->GetYaxis()->SetTickLength(tylength);
	hP3X3sl[k]->GetYaxis()->CenterTitle();
	
	hP3X3sl[k]->GetZaxis()->SetLabelFont(fonttype);
	hP3X3sl[k]->GetZaxis()->SetLabelSize(fontsize);
	hP3X3sl[k]->GetZaxis()->SetTitleFont(fonttype);
	hP3X3sl[k]->GetZaxis()->SetTitleSize(tfontsize);
	hP3X3sl[k]->GetZaxis()->SetTitleOffset(tzoffset);
	hP3X3sl[k]->GetZaxis()->SetTickLength(tylength);
	hP3X3sl[k]->GetZaxis()->CenterTitle();    

	hFrP3X3->Draw("axis");
	hP3X3sl[k]->Draw(drawopt + "same");
	
	slXmean2[k] = new TLine(hP3X3sl[k]->GetMean(1),p3Min,hP3X3sl[k]->GetMean(1),p3Max);
	slXmean2[k]->SetLineColor(kGray+2);
	slXmean2[k]->SetLineStyle(2);
	slXmean2[k]->Draw();

	
	slPxmean2[k] = new TLine(x3Min,hP3X3sl[k]->GetMean(2),x3Max,hP3X3sl[k]->GetMean(2));
	slPxmean2[k]->SetLineColor(kGray+2);
	slPxmean2[k]->SetLineStyle(2);
	slPxmean2[k]->Draw();

	
	if(opt.Contains("elli") ) 
	  if(sellipP3X3[k] != NULL)
	    sellipP3X3[k]->Draw();

	gPad->Update();
	palette = (TPaletteAxis*)hP3X3sl[k]->GetListOfFunctions()->FindObject("palette");
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

	Double_t y1 = gPad->GetBottomMargin();
	Double_t y2 = 1 - gPad->GetTopMargin();
	Double_t x1 = gPad->GetLeftMargin();
	Double_t x2 = 1 - gPad->GetRightMargin();
      
	textStat[k] = new TPaveText(x1+0.02,y2-0.40,x1+0.20,y2-0.05,"NDC");
	PGlobals::SetPaveTextStyle(textStat[k],12); 
	textStat[k]->SetTextColor(kGray+3);
	textStat[k]->SetTextFont(42);

	char text[64];
	sprintf(text,"Q = %5.1f pC",sCharge[k]);
	textStat[k]->AddText(text);
	sprintf(text,"#Deltay = %5.2f #mum",sy_rms[k]);
	textStat[k]->AddText(text);
	sprintf(text,"#Deltap_{y}/mc = %5.2f",spy_rms[k]);
	textStat[k]->AddText(text);
	sprintf(text,"#varepsilon_{n,y} = %5.2f %s",semity[k],emitSUnit.c_str());
	textStat[k]->AddText(text);
	sprintf(text,"#beta_{y} = %5.2f mm",1E-3*sbetay[k]);
	textStat[k]->AddText(text);
	textStat[k]->Draw();


	textSlice[k] = new TPaveText(x1+0.4,y2-0.2,x2-0.02,y2-0.0,"NDC");
	PGlobals::SetPaveTextStyle(textSlice[k],32); 
	textSlice[k]->SetTextColor(kBlack);
	textSlice[k]->SetTextFont(42);

	sprintf(text,"%5.2f #mum < #zeta < %5.2f #mum",sBinLim[k],sBinLim[k+1]);
	textSlice[k]->AddText(text);
	textSlice[k]->Draw();

	lFrame = new TBox(gPad->GetUxmin(), gPad->GetUymin(),
			  gPad->GetUxmax(), gPad->GetUymax());
	lFrame->SetFillStyle(0);
	lFrame->SetLineColor(PGlobals::frameColor);
	lFrame->SetLineWidth(PGlobals::frameWidth);
	lFrame->Draw();

	gPad->RedrawAxis(); 

	if(ic+1==ndiv) {
	  CA5->cd(0);
	  CA5->Print(fOutName5 + ".ps");
	}
      }

      CA5->Print(fOutName5 + ".ps]");

      gSystem->Exec("ps2pdf " + fOutName5 + ".ps " + fOutName5 + ".pdf");
      gSystem->Exec("rm -rf " + fOutName5 + ".ps"); 
    } 
  }
  
  if(opt.Contains("file")) {
    TString filename = Form("histos.root");
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
    hP3X3->Write("hP3X3",TObject::kOverwrite);
    hX2X1->Write("hX2X1",TObject::kOverwrite);
    hX3X1->Write("hX3X1",TObject::kOverwrite);
    hX3X2->Write("hX3X1",TObject::kOverwrite);

    
    gErms->Write("gErms",TObject::kOverwrite);
    gXrms->Write("gXrms",TObject::kOverwrite);
    gemitX->Write("gemitX",TObject::kOverwrite);
    gYrms->Write("gYrms",TObject::kOverwrite);
    gemitY->Write("gemitY",TObject::kOverwrite);

    for(Int_t k=0;k<SNbin;k++) {

      sprintf(hName,"hP1sl_%i",k);	
      hP1sl[k]->Write(hName,TObject::kOverwrite);

      sprintf(hName,"hP2X2sl_%i",k);	
      hP2X2sl[k]->Write(hName,TObject::kOverwrite);

      sprintf(hName,"hP3X3sl_%i",k);	
      hP3X3sl[k]->Write(hName,TObject::kOverwrite);

    }
      

    ofile->Close();
  }
    
  // Delete[] newly created vectors
  delete sBinLim;

  delete sxmean;
  delete symean;
  delete sx2mean;
  delete sy2mean;
  delete sxymean;
  delete sNtotal;
  delete sxrms2;  
  delete syrms2; 
  delete sxrms;  
  delete syrms;  
  delete sxyrms2;
  delete xbin;
  delete semit;
  delete sNEtotal; 
  delete sEmean;
  delete sErms;

  // delete gemitX;
  // delete gYrms;
  // delete gErms;
  // delete gErmsB;

}
