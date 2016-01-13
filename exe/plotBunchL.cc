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

int main(int argc,char *argv[]) {
  if(argc<=2) {
    printf("\n Usage: %s <input file> \n",argv[0]);
    printf("      <--png> <--pdf> <--eps> <--units> <--comov>\n");
    printf("      <--file> \n");
    return 0;
  }

  // Input file
  TString ifile = "";

  // General options
  TString opt = "";
  
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
    } else if(arg.Contains("--notext")){
      opt += "notext"; 
    } else if(arg.Contains("--noinfo")){
      opt += "noinfo"; 
    } else if(arg.Contains("--bw")){
      opt += "bw"; 
    } else if(arg.Contains("--zmean")){
      opt += "zmean"; 
    } else if(arg.Contains("--elli")){
      opt += "elli"; 
    } else if(arg.Contains("--astra")){
      opt += "astra"; 
    } else if(arg.Contains("-spoil")) {
      char ss[6];
      sscanf(arg,"%6s%f",ss,&X);
      opt += "spoil";
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
  gROOT->Macro("PPalettes.C");

  if(opt.Contains("grid")) {
    gStyle->SetPadGridX(1);
    gStyle->SetPadGridY(1);
  }

  // Spatial resolution
  // Float_t dx1 = 0.04;
  // Float_t dx2 = 0.04;
  // Float_t dx3 = 0.04;

  // Bining, intervals, labels, etc.
  Int_t x1Nbin = 100;
  Int_t p1Nbin = 100;
  Int_t x2Nbin = 100;
  Int_t p2Nbin = 100;
  Int_t x3Nbin = 100;
  Int_t p3Nbin = 100;
    
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
  Float_t rangefactor = 0.3;
  
  x1BinMin = varMin[0];
  x1BinMax = varMax[0];

  x1Min = varMin[0] - rangefactor*(varMax[0]-varMin[0]);
  x1Max = varMax[0] + rangefactor*(varMax[0]-varMin[0]);
  //x1Min = -72390512;
  //x1Max = -72390232;
  //x1Min = -25846100;

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
  
  for(UInt_t i=0;i<Np;i++) {

    hX1->Fill(var[0][i],Q/Np);
    hP1->Fill(var[3][i],Q/Np);

    hP1X1->Fill(var[0][i],var[3][i],Q/Np);
  }


  // Get slices interval
  Float_t iMax = hX1->GetBinContent(hX1->GetMaximumBin());
  for(Int_t i=1;i<=hX1->GetNbinsX();i++) {
    if(hX1->GetBinContent(i)>0.05*iMax) {
      x1BinMin = hX1->GetBinCenter(i);
      break;
    }
  }
  for(Int_t i=hX1->GetNbinsX();i>=0;i--) {
    if(hX1->GetBinContent(i)>0.05*iMax) {
      x1BinMax = hX1->GetBinCenter(i);
      break;
    }
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
    hP2X2->Fill(var[1][i],var[4][i],Q/Np);
    hP2X2sl[iBin]->Fill(var[1][i],var[4][i],Q/Np);

    hP3X3->Fill(var[2][i],var[5][i],Q/Np);
    hP3X3sl[iBin]->Fill(var[2][i],var[5][i],Q/Np);

    hX2X1->Fill(var[0][i],var[1][i],Q/Np);
    hX3X1->Fill(var[0][i],var[2][i],Q/Np);
    hX3X2->Fill(var[1][i],var[2][i],Q/Np);

    hP1sl[iBin]->Fill(var[3][i],Q/Np);   
    
  }
  
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
    semitx[k] = semit[k];

    Double_t grel = hP1sl[k]->GetMean() * PUnits::GeV/PConst::ElectronMassE;
    sbetax[k] = beta * grel;
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
    semity[k] = semit[k];

    Double_t grel = hP1sl[k]->GetMean() * PUnits::GeV/PConst::ElectronMassE;
    sbetay[k] = beta * grel;

    
  }

  
  cout << "\n done! " << endl;

  // To current.
  Float_t binSize = (x1Max - x1Min)/x1Nbin;
  Float_t lightspeed =  PConst::c_light / (PUnits::um/PUnits::femtosecond);  
  hX1->Scale(1000*lightspeed/binSize);


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
  if(opt.Contains("zmean")) {
    hX1->SetBins(x1Nbin,x1Min-zmean,x1Max-zmean);
    hP1X1->SetBins(x1Nbin,x1Min-zmean,x1Max-zmean,p1Nbin,p1Min,p1Max);
    hX2X1->SetBins(x1Nbin,x1Min-zmean,x1Max-zmean,x2Nbin,x2Min-x_mean,x2Max-x_mean);
    hP2X2->SetBins(x2Nbin,x2Min-x_mean,x2Max-x_mean,p2Nbin,p2Min,p2Max);
    hP3X3->SetBins(x3Nbin,x3Min-y_mean,x3Max-y_mean,p3Nbin,p3Min,p3Max);
    hX3X1->SetBins(x1Nbin,x1Min-zmean,x1Max-zmean,x3Nbin,x3Min-y_mean,x3Max-y_mean);
    hX3X2->SetBins(x2Nbin,x2Min-x_mean,x2Max-x_mean,x3Nbin,x3Min-y_mean,x3Max-y_mean);
    
    for(Int_t k=0;k<SNbin;k++) {
      xbin[k] -= zmean;
      sBinLim[k] -= zmean;
    }
    sBinLim[SNbin] -= zmean;
    zmean = 0.0;
    x_mean = 0;
    y_mean = 0;
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
    Double_t xMin = hX1->GetXaxis()->GetXmin();
    Double_t xMax = hX1->GetXaxis()->GetXmin() + (hX1->GetXaxis()->GetXmax()
						 -hX1->GetXaxis()->GetXmin()) * 0.2;
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

    TPaveText *textMom = new TPaveText(0.55,0.03,0.80,0.13,"NDC");
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
    //     sprintf(ctext,"#LTp_{z}#GT_{rms} = %5.2f GeV/c",pzrms);
    sprintf(ctext,"#Delta#gamma/#gamma = %4.1f %%",100*pzrms/pzmean);
    textInfo->AddText(ctext);
    sprintf(ctext,"#varepsilon_{x} = %5.2f #mum",emitx);
    textInfo->AddText(ctext);
    sprintf(ctext,"#varepsilon_{y} = %5.2f #mum",emity);
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
    PGlobals::CanvasAsymPartition(C,NPad,lMargin,rMargin,bMargin,tMargin,factor);

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

    hFrame[1]->GetYaxis()->SetRangeUser(hP1X1->GetYaxis()->GetXmin(),hP1X1->GetYaxis()->GetXmax());

    hFrame[1]->GetYaxis()->SetTitle("p_{z} [GeV/c]");
    
    hFrame[1]->Draw();

    gP1left->SetLineWidth(2);
    gP1left->Draw("F");
    gP1left->Draw("L");

    TLine lZmean(zmean,hP1X1->GetYaxis()->GetXmin(),zmean,hP1X1->GetYaxis()->GetXmax());
    lZmean.SetLineColor(kGray+2);
    lZmean.SetLineStyle(2);
    lZmean.Draw();

    TLine lPmean(hP1X1->GetXaxis()->GetXmin(),pzmean,hP1X1->GetXaxis()->GetXmax(),pzmean);
    lPmean.SetLineColor(kGray+2);
    lPmean.SetLineStyle(2);
    lPmean.Draw();

    // 2D histogram z range
    Double_t dmax = hP1X1->GetMaximum();
    Double_t dmin = 0.0;
    hP1X1->GetZaxis()->SetRangeUser(1.1*dmin,dmax);

    hP1X1->GetZaxis()->CenterTitle();
    hP1X1->GetZaxis()->SetTitleFont(fonttype);

    hP1X1->Draw(drawopts);
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
    Leg->AddEntry(hX1  ,"Current [kA]","L");
    // Leg->AddEntry(gErms,"Energy spread (GeV)","PL");
    Leg->AddEntry(gErms,"Energy spread [0.1%]","PL");
    Leg->AddEntry(gemitX,"Emittance_{x} [#mum]","PL");
    Leg->AddEntry(gemitY,"Emittance_{y} [#mum]","PL");
    //Leg->AddEntry(gYrms,"Bunch width [#mum]","PL");


    hFrame[0]->GetYaxis()->SetRangeUser(0.0,1.1*yMax);

    hFrame[0]->GetXaxis()->SetTitle("#zeta [#mum]");

    hFrame[0]->Draw();

    hX1->GetYaxis()->SetNdivisions(503);
    hX1->SetLineWidth(2);
    hX1->SetFillStyle(1001);
    hX1->SetFillColor(PGlobals::elecFill);
    // hX1->SetLineColor(kBlue);

    //hX1->Smooth();
    //    hX1->Draw("FL same");
    hX1->Draw("LF2 same");
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
    hX2X1cl->Draw(drawopts);
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
    sprintf(text,"Q = %5.1f pC",Charge);
    textInfoX2X1->AddText(text);
    sprintf(text,"#Delta#zeta = %5.2f #mum",zrms);
    textInfoX2X1->AddText(text);
    sprintf(text,"#Deltax = %5.2f #mum",x_rms);
    textInfoX2X1->AddText(text);
    // sprintf(text,"#varepsilon_{x} = %5.2f #mum",emitx);
    // textInfoX2X1->AddText(text);
    // sprintf(text,"#beta_{x} = %5.2f mm",1E-3*betax);
    // textInfoX2X1->AddText(text);
    textInfoX2X1->Draw();
    
    C1->cd(0);
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
    hX3X1cl->Draw(drawopts);
    
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

    sprintf(text,"Q = %5.1f pC",Charge);
    textInfoX3X1->AddText(text);
    sprintf(text,"#Delta#zeta = %5.2f #mum",zrms);
    textInfoX3X1->AddText(text);
    sprintf(text,"#Deltay = %5.2f #mum",y_rms);
    textInfoX3X1->AddText(text);
    // sprintf(text,"#varepsilon_{x} = %5.2f #mum",emitx);
    // textInfoX3X1->AddText(text);
    // sprintf(text,"#beta_{x} = %5.2f mm",1E-3*betax);
    // textInfoX3X1->AddText(text);
    textInfoX3X1->Draw();
    
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

    hX3X2cl->Draw(drawopts);
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

    sprintf(text,"Q = %5.1f pC",Charge);
    textStatX3X2->AddText(text);
    sprintf(text,"#Deltax = %5.2f #mum",x_rms);
    textStatX3X2->AddText(text);
    sprintf(text,"#Deltay = %5.2f #mum",y_rms);
    textStatX3X2->AddText(text);
    textStatX3X2->Draw();
    
    
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

    TLine lPxmean(hFrame[0]->GetXaxis()->GetXmin(),pxmean,hFrame[0]->GetXaxis()->GetXmax(),pxmean);
    lPxmean.SetLineColor(kGray+2);
    lPxmean.SetLineStyle(2);
    lPxmean.Draw();


    TH2F *hP2X2cl = (TH2F*) hP2X2->Clone("hP2X2cl");
    hP2X2cl->Draw(drawopts);
    
    if(opt.Contains("elli") ) 
      for(Int_t k=0;k<SNbin;k++)
	if(sellipP2X2[k] != NULL)
	  sellipP2X2[k]->Draw();

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
    sprintf(text,"#varepsilon_{x} = %5.2f #mum",emitx);
    textStatInt->AddText(text);
    sprintf(text,"#beta_{x} = %5.2f mm",1E-3*betax);
    textStatInt->AddText(text);
    textStatInt->Draw();
    
    
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

    TLine lPymean(hFrame[0]->GetXaxis()->GetXmin(),pymean,hFrame[0]->GetXaxis()->GetXmax(),pymean);
    lPymean.SetLineColor(kGray+2);
    lPymean.SetLineStyle(2);
    lPymean.Draw();


    TH2F *hP3X3cl = (TH2F*) hP3X3->Clone("hP3X3cl");
    hP3X3cl->Draw(drawopts);
        
    if(opt.Contains("elli") ) 
      for(Int_t k=0;k<SNbin;k++)
	if(sellipP3X3[k] != NULL)
	  sellipP3X3[k]->Draw();

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
    sprintf(text,"#varepsilon_{y} = %5.2f #mum",emity);
    textStatInt->AddText(text);
    sprintf(text,"#beta_{y} = %5.2f mm",1E-3*betay);
    textStatInt->AddText(text);
    textStatInt->Draw();
    
    
    TString fOutName3 = ifile;
    fOutName3.ReplaceAll(".txt","");
    fOutName3.Append("-p3x3");
    PGlobals::imgconv(C3,fOutName3,opt);


    // Slices plots

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
    hP2X2->Draw(drawopt);

    //    ellipP2X2->Draw();

    if(opt.Contains("elli") ) 
      for(Int_t k=0;k<SNbin;k++)
	if(sellipP2X2[k] != NULL)
	  sellipP2X2[k]->Draw();

      

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
    sprintf(text,"#varepsilon_{x} = %5.2f #mum",emitx);
    textStatInt->AddText(text);
    sprintf(text,"#beta_{x} = %5.2f mm",1E-3*betax);
    textStatInt->AddText(text);
    textStatInt->Draw();

    TPaveText **textStat = new TPaveText*[SNbin];
    TPaveText **textSlice = new TPaveText*[SNbin];
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

      hP2X2sl[k]->Draw(drawopt);

      if(opt.Contains("elli") ) 
	if(sellipP2X2[k] != NULL)
	  sellipP2X2[k]->Draw();

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
      sprintf(text,"#varepsilon_{x} = %5.2f #mum",semitx[k]);
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

    hP3X3->GetXaxis()->SetLabelSize(0.08);
    hP3X3->GetXaxis()->SetTitleSize(0.08);
    hP3X3->GetXaxis()->SetTitleOffset(1.0);
    hP3X3->GetXaxis()->CenterTitle();

    hP3X3->GetYaxis()->SetLabelSize(0.08);
    hP3X3->GetYaxis()->SetTitleSize(0.08);
    hP3X3->GetYaxis()->SetTitleOffset(0.8);
    hP3X3->GetYaxis()->CenterTitle();

    hP3X3->GetZaxis()->SetLabelSize(0.08);
    hP3X3->GetZaxis()->SetTitleSize(0.08);
    hP3X3->GetZaxis()->SetTitleOffset(0.8);
    hP3X3->GetZaxis()->CenterTitle();


    CA5->cd(1);
    hP3X3->Draw(drawopt);

    //    ellipP3X3->Draw();

    if(opt.Contains("elli") ) 
      for(Int_t k=0;k<SNbin;k++)
	if(sellipP3X3[k] != NULL)
	  sellipP3X3[k]->Draw();
      

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
    sprintf(text,"#varepsilon_{y} = %5.2f #mum",emity);
    textStatInt->AddText(text);
    sprintf(text,"#beta_{y} = %5.2f mm",1E-3*betay);
    textStatInt->AddText(text);
    textStatInt->Draw();

    textStat = new TPaveText*[SNbin];
    textSlice = new TPaveText*[SNbin];
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

      hP3X3sl[k]->GetXaxis()->SetLabelSize(0.08);
      hP3X3sl[k]->GetXaxis()->SetTitleSize(0.08);
      hP3X3sl[k]->GetXaxis()->SetTitleOffset(1.0);
      hP3X3sl[k]->GetXaxis()->CenterTitle();

      hP3X3sl[k]->GetYaxis()->SetLabelSize(0.08);
      hP3X3sl[k]->GetYaxis()->SetTitleSize(0.08);
      hP3X3sl[k]->GetYaxis()->SetTitleOffset(0.8);
      hP3X3sl[k]->GetYaxis()->CenterTitle();

      hP3X3sl[k]->GetZaxis()->SetLabelSize(0.08);
      hP3X3sl[k]->GetZaxis()->SetTitleSize(0.08);
      hP3X3sl[k]->GetZaxis()->SetTitleOffset(0.8);
      hP3X3sl[k]->GetZaxis()->CenterTitle();

      hP3X3sl[k]->Draw(drawopt);

      if(opt.Contains("elli") ) 
	if(sellipP3X3[k] != NULL)
	  sellipP3X3[k]->Draw();

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
      sprintf(text,"#varepsilon_{y} = %5.2f #mum",semity[k]);
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


      if(ic+1==ndiv) {
	CA5->cd(0);
	CA5->Print(fOutName5 + ".ps");
      }
    }

    CA5->Print(fOutName5 + ".ps]");

    gSystem->Exec("ps2pdf " + fOutName5 + ".ps " + fOutName5 + ".pdf");
    gSystem->Exec("rm -rf " + fOutName5 + ".ps"); 
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
    
  // Delete newly created vectors
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
