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
    printf("      <-it(time)> <-i(initial time)> <-f(final time)> <-s(time step)>\n");
    printf("      <--png> <--pdf> <--eps> <--units> <--comov>\n");
    printf("      <--file> <--loop>\n");
    return 0;
  }

  // General options
  TString   sim = "";
  Float_t time = 0.0;
  Int_t    itime = 0;
  Int_t  iStart = -1;
  Int_t    iEnd = -1;
  Int_t   iStep = 1;
  TString   opt = "";
  Int_t   index = 1;
  
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
    } else if(arg.Contains("-index")) {
      char ss[6];
      sscanf(arg,"%6s%i",ss,&index);
    } else if(arg.Contains("-t")) {
      char ss[2];
      sscanf(arg,"%2s%f",ss,&time);
    } else if(arg.Contains("-it")) {
      char ss[2];
      sscanf(arg,"%2s%i",ss,&itime);
    } else if(arg.Contains("-i")) {
      char ss[2];
      sscanf(arg,"%2s%i",ss,&iStart);
    } else if(arg.Contains("-f")) {
      char ss[2];
      sscanf(arg,"%2s%i",ss,&iEnd);
    } else if(arg.Contains("-s")) {
      char ss[2];
      sscanf(arg,"%2s%i",ss,&iStep);
    } else if(arg.Contains("-pmin")) {
      char ss[5];
      sscanf(arg,"%5s%f",ss,&Pmin);
    } else if(arg.Contains("-pmax")) {
      char ss[5];
      sscanf(arg,"%5s%f",ss,&Pmax);
    } else if( !(arg.Contains("pitz") || arg.Contains("flash") || arg.Contains("regae") || arg.Contains("Gauss") || arg.Contains("pwfa") || arg.Contains("facet") || arg.Contains("rake")) ) {
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
  // PData *pData = PData::Get(sim.Data());

  TString specName = "Plasma";
  if(index ==1) specName = "DriverBeam";
  if(index ==2) specName = "InjectedBeam";

  if(iStart<0) iStart = itime;
  if(iEnd<=iStart) iEnd = iStart;
  
  // Some plasma constants
  Float_t n0 = 5E16 * (1./PUnits::cm3);
  Float_t kp = PFunc::PlasmaWavenumber(n0);
  Float_t skindepth = 1.;
  if(kp!=0.0) skindepth = 1/kp;
  Float_t E0 = PFunc::PlasmaWBField(n0);
  Float_t gamma = 45009.8;
  Float_t vel = TMath::Sqrt(1. - 1./(gamma*gamma) );
  Float_t BeamMass = 0.511 * PUnits::MeV;
  
  cout << endl;
  cout << Form(" - Plasma density   = %8.4e cm3",n0 / (1./PUnits::cm3)) << endl;
  cout << Form(" - Skin depth       = %8.4f um", skindepth / PUnits::um) << endl;
  cout << Form(" - E_0              = %8.4f GV/m", E0 * PUnits::meter / PUnits::GV ) << endl;
  cout << Form(" - Electron mass    = %8.4f MeV", BeamMass/PUnits::MeV ) << endl;

  // z start of the plasma in normalized units.
  Float_t zStartPlasma = 121.;
  // z start of the beam in normalized units.
  Float_t zStartBeam = 18.;

  opt += "comovcenter";
  
  // Time looper
  for(Int_t i=iStart; i<iEnd+1; i+=iStep) {

    itime = i;
    // pData->LoadFileNames(itime);    
    // if(itime==iStart) pData->Print();
    
    // if(!pData->IsInit()) continue;

    // Int_t Nspecies = pData->NSpecies();
    // if(index>Nspecies-1) {
    //   return 0;
    // }
    // if(!pData->GetRawFileName(index)) {
    //   return 0;    
    // }

    // Time in OU
    // Float_t Time = pData->GetRealTime();
    Float_t Time = time;
    
    // Spatial resolution
    Float_t dx1 = 0.04;
    Float_t dx2 = 0.04;
    Float_t dx3 = 0.04;

    // Bining, intervals, labels, etc.
    Int_t x1Nbin = 100;
    Int_t p1Nbin = 20;
    Int_t x2Nbin = 100;
    Int_t p2Nbin = 100;
    
    // Slices
    Int_t SNbin = 40;
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
    if(sim.Contains("facet_DDR.HiPACE")) {
      
      dx1 = 15.0/1024;
      dx2 = 10.0/256;
      dx3 = 10.0/256;
      
      x1Min = 12.2;
      x1Max = 13.1; 
      
      SNbin = 20;
      x1BinMin = 12.45;
      x1BinMax = 12.90;
	
    }
        
    // --------------------------------------------------
    // READ FROM HiPACE
    TString filename;
    if(Time/100000.0>=1) 
      filename = Form("%s/DATA/%sPhaseSpace_time_%6.2f",sim.Data(),specName.Data(),Time);
    else if(Time/10000.0>=1)
      filename = Form("%s/DATA/%sPhaseSpace_time_0%6.2f",sim.Data(),specName.Data(),Time);
      
    cout << Form("\n 1. Reading file : ") << filename.Data() << endl;

    FILE * pFile;
    size_t lSize;
    double * buffer;
    size_t result;
  
    pFile = fopen (filename.Data(),"rb" );
    if (pFile==NULL) {fputs ("File error",stderr); exit (1);}
  
    // obtain file size:
    fseek (pFile , 0 , SEEK_END);
    lSize = ftell (pFile);
    rewind (pFile);

    // allocate memory to contain the whole file:
    buffer = (double*) malloc (lSize);
    if (buffer == NULL) {fputs ("Memory error",stderr); exit (2);}

    // copy the file into the buffer:
    result = fread (buffer,1,lSize,pFile);
    if (result != lSize) {fputs ("Reading error",stderr); exit (3);}
  
    /* the whole file is now loaded in the memory buffer. */
    
    // terminate
    fclose (pFile);

    cout << "Buffer read!" << endl;

    const Int_t NvarHiP = 9;
    Double_t varHiP[NvarHiP];  
    char varnameHiP[NvarHiP][8] = {{"x1"},{"x2"},{"x3"},{"p1"},{"p2"},{"p3"},{"q"},{"ipart"},{"iproc"}};    
    Int_t Np=lSize/(NvarHiP*sizeof(double));
    
    const Int_t Nvar = 8;
    Double_t var[Nvar][Np]; 
    char varname[Nvar][4] = {{"ene"},{"p1"},{"p2"},{"p3"},{"q"},{"x1"},{"x2"},{"x3"}};
    
    for(int i_part=0; i_part<Np; i_part++) {
      // cout << Form("Particle %i : ",i_part) << endl;
      for(Int_t i_var=0; i_var<NvarHiP; i_var++) {
	
	varHiP[i_var]=buffer[i_var+i_part*NvarHiP];

	Int_t i_var_os = 0;
	if(i_var<3) i_var_os = i_var+5;
	else if(i_var<=6) i_var_os = i_var-2;
	else continue;
	
	// cout << Form ("var[%i][%i] = varHiP[%i] = %e",i_part,i_var_os,i_var,varHiP[i_var]) << endl;
	
	var[i_var_os][i_part] = varHiP[i_var];
	
      }
      // cout << endl;
      var[0][i_part] = TMath::Sqrt(1 + var[1][i_part]*var[1][i_part] + var[2][i_part]*var[2][i_part] + var[3][i_part]*var[3][i_part]) - 1.0;
    }
    
    free (buffer);
    
    cout << Form("  %i  particles read! " , Np) << endl;
    
    // ----------------------------------------------------------------------------------------------
    
    // Centering time and z position:
    // Float_t shiftz = pData->Shift(opt);
    Float_t shiftz = 0;
    // if(opt.Contains("center")) {
    //   if(opt.Contains("comov"))
    // 	shiftz += zStartBeam;
    //   else
    // 	shiftz += zStartPlasma;
    // }
    
    // if(opt.Contains("comov")) 
    //   shiftz += vel * Time;
    
    
    // if(opt.Contains("center")) {
    //   Time -= zStartPlasma;
    //   if(opt.Contains("comov"))      // Centers on the head of the beam.
    // 	Time += zStartBeam;
    // } 
    

    if(opt.Contains("autop")) {

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

	if(var[5][i]-shiftz<x1BinMin || var[5][i]-shiftz>x1BinMax ) continue; 

	if(var[1][i]<MinP1) MinP1 = var[1][i];
	if(var[1][i]>MaxP1) MaxP1 = var[1][i];
	if(var[2][i]<MinP2) MinP2 = var[2][i];
	if(var[2][i]>MaxP2) MaxP2 = var[2][i];
	if(var[3][i]<MinP3) MinP3 = var[3][i];
	if(var[3][i]>MaxP3) MaxP3 = var[3][i];
	if(var[6][i]<MinX2) MinX2 = var[6][i];
	if(var[6][i]>MaxX2) MaxX2 = var[6][i];
	if(var[7][i]<MinX3) MinX3 = var[7][i];
	if(var[7][i]>MaxX3) MaxX3 = var[7][i];
      }
      
      p1Min = MinP1 - 0.3*(MaxP1-MinP1);
      p1Max = MaxP1 + 0.3*(MaxP1-MinP1);
      p2Min = MinP2 - 0.3*(MaxP2-MinP2);
      p2Max = MaxP2 + 0.3*(MaxP2-MinP2);
      p3Min = MinP3 - 0.3*(MaxP3-MinP3);
      p3Max = MaxP3 + 0.3*(MaxP3-MinP3);

      x2Min = MinX2 - 0.3*(MaxX2-MinX2);
      x2Max = MaxX2 + 0.3*(MaxX2-MinX2);
      x3Min = MinX3 - 0.3*(MaxX3-MinX3);
      x3Max = MaxX3 + 0.3*(MaxX3-MinX3);

      // p1Nbin = x1Nbin;
      // x2Nbin = ceil ((x2Max - x2Min)/(1.0*dx2));
      // p2Nbin = x2Nbin;
      
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
	if(var[5][i]-shiftz<MinX1) MinX1 = var[5][i]-shiftz;
	if(var[5][i]-shiftz>MaxX1) MaxX1 = var[5][i]-shiftz;
	if(var[6][i]<MinX2) MinX2 = var[6][i];
	if(var[6][i]>MaxX2) MaxX2 = var[6][i];
	if(var[7][i]<MinX3) MinX3 = var[7][i];
	if(var[7][i]>MaxX3) MaxX3 = var[7][i];
      }
      
      p1Min = MinP1 - 0.3*(MaxP1-MinP1);
      p1Max = MaxP1 + 0.3*(MaxP1-MinP1);
      p2Min = MinP2 - 0.3*(MaxP2-MinP2);
      p2Max = MaxP2 + 0.3*(MaxP2-MinP2);
      p3Min = MinP3 - 0.3*(MaxP3-MinP3);
      p3Max = MaxP3 + 0.3*(MaxP3-MinP3);

      x1Min = MinX1 - 0.3*(MaxX1-MinX1);
      x1Max = MaxX1 + 0.3*(MaxX1-MinX1);
      x2Min = MinX2 - 0.3*(MaxX2-MinX2);
      x2Max = MaxX2 + 0.3*(MaxX2-MinX2);
      x3Min = MinX3 - 0.3*(MaxX3-MinX3);
      x3Max = MaxX3 + 0.3*(MaxX3-MinX3);

      x1BinMin = MinX1;
      x1BinMax = MaxX1;

      x1Nbin = ceil ((x1Max - x1Min)/(1.0*dx1));
      p1Nbin = x1Nbin;
      
      x2Nbin = ceil ((x2Max - x2Min)/(1.0*dx2));
      p2Nbin = x2Nbin;
      
      SNbin = ceil ((x1BinMax - x1BinMin)/(2.0*dx1));
    } 

    if(p1Min < 0.0) p1Min = 0.01;

    cout << Form(" x1 range (N = %i):  x1Min = %f  x1Max = %f ", x1Nbin, x1Min, x1Max) << endl;
    cout << Form(" p1 range (N = %i):  p1Min = %f  p1Max = %f ", p1Nbin, p1Min, p1Max) << endl;
    cout << Form(" x2 range (N = %i):  x2Min = %f  x2Max = %f ", x2Nbin, x2Min, x2Max) << endl;
    cout << Form(" p2 range (N = %i):  p2Min = %f  p2Max = %f ", p2Nbin, p2Min, p2Max) << endl;
    

    cout <<  Form(" Number of bins = %i . resolution = %.2f", x1Nbin,(x1Max - x1Min)/x1Nbin) << endl;    
    
    // Histograms
    char hName[8];

    sprintf(hName,"hX1");
    TH1F *hX1 = (TH1F*) gROOT->FindObject(hName);
    if(hX1) delete hX1;
    hX1 = new TH1F(hName,"",x1Nbin,x1Min,x1Max);

    sprintf(hName,"hP1");
    TH1F *hP1 = (TH1F*) gROOT->FindObject(hName);
    if(hP1) delete hP1;
    hP1 = new TH1F(hName,"",p1Nbin,p1Min,p1Max);

    sprintf(hName,"hP1X1");
    TH2F *hP1X1 = (TH2F*) gROOT->FindObject(hName);
    if(hP1X1) delete hP1X1;
    hP1X1 = new TH2F(hName,"",x1Nbin,x1Min,x1Max,p1Nbin,p1Min,p1Max);

    sprintf(hName,"hP2X2");
    TH2F *hP2X2 =  (TH2F*) gROOT->FindObject(hName);
    if(hP2X2) delete hP2X2;
    hP2X2 = new TH2F(hName,"",x2Nbin,x2Min,x2Max,p2Nbin,p2Min,p2Max);

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
      sprintf(hName,"hP2X2sl_%i",k);
      hP2X2sl[k] = (TH2F*) gROOT->FindObject(hName);
      if(hP2X2sl[k]) delete hP2X2sl[k];
      hP2X2sl[k] = new TH2F(hName,"",x2Nbin,x2Min,x2Max,p2Nbin,p2Min,p2Max);

      hP2X2sl[k]->GetXaxis()->SetTitle("k_{p} y");
      hP2X2sl[k]->GetYaxis()->SetTitle("p_{y}/mc");


      sprintf(hName,"hP1sl_%i",k);
      hP1sl[k] = (TH1F*) gROOT->FindObject(hName);
      if(hP1sl[k]) delete hP1sl[k];
      hP1sl[k] = new TH1F(hName,"",p1Nbin,p1Min,p1Max);

    }

    // Filling histos
    cout << Form("\n 2. Filling histograms from file ... ") ;
    for(UInt_t i=0;i<Np;i++) {

      var[5][i] = var[5][i] - shiftz;

      if(var[5][i]<x1Min || var[5][i]>x1Max ) continue; 
      if(var[6][i]<x2Min || var[6][i]>x2Max ) continue; 
      if(var[7][i]<x3Min || var[7][i]>x3Max ) continue; 
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
      hP2X2sl[iBin]->Fill(var[6][i],var[2][i],TMath::Abs(var[4][i]));

    }
    cout << " done! " << endl;

    // Integrated long. emittance:
    cout << Form("\n 3. Calculating integrated quantities: ") << endl ;

    // Longitudinal phasespace:
    // --

    Float_t xmean  = 0.0;
    Float_t ymean  = 0.0;
    Float_t x2mean = 0.0;
    Float_t y2mean = 0.0;
    Float_t xymean = 0.0;
    Float_t Ntotal = 0.0;
    for(Int_t i=1;i<=x1Nbin;i++) {
      Float_t x = hP1X1->GetXaxis()->GetBinCenter(i);
      // if(x<xmin || x>xmax) continue;
      for(Int_t j=1;j<=p1Nbin;j++) {
	Float_t y = hP1X1->GetYaxis()->GetBinCenter(j);
	// if(y<ymin || y>ymax) continue;
	Float_t value = TMath::Abs(hP1X1->GetBinContent(i,j));
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

    Float_t xrms2  = x2mean - xmean*xmean;
    Float_t yrms2  = y2mean - ymean*ymean;
    Float_t xrms   = ( xrms2>0.0 ) ? TMath::Sqrt(xrms2) : 0.0 ;
    Float_t yrms   = ( yrms2>0.0 ) ? TMath::Sqrt(yrms2) : 0.0 ;
    Float_t xyrms2 = xymean - xmean*ymean;
    Float_t emittance2 = xrms2*yrms2 - xyrms2*xyrms2;
    Float_t emittance = (emittance2>0.0) ? TMath::Sqrt(emittance2) : 0.0 ;
    
    // cout << Form("  xMean = %7.3f   yMean = %7.3f",xmean,ymean) << endl;
    // cout << Form("  xRms  = %7.3f   yRms  = %7.3f",xrms,yrms) << endl;
    // cout << Form("  Emittance = %7.3f",emittance) << endl;

    Float_t emitz = emittance;
    Float_t zmean = xmean;
    Float_t zrms = xrms;
    Float_t pzmean = ymean;
    Float_t pzrms = yrms;
    // ----------------------------------------


    // Transverse phasespace
    xmean = 0.0;
    ymean = 0.0;
    x2mean = 0.0;
    y2mean = 0.0;
    xymean = 0.0;
    Ntotal = 0.0;
    for(Int_t i=1;i<=x2Nbin;i++) {
      Float_t x = hP2X2->GetXaxis()->GetBinCenter(i);
      // if(x<xmin || x>xmax) continue;
      for(Int_t j=1;j<=p2Nbin;j++) {
	Float_t y = hP2X2->GetYaxis()->GetBinCenter(j);
	// if(y<ymin || y>ymax) continue;
	Float_t value = TMath::Abs(hP2X2->GetBinContent(i,j));
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

    xrms2  = x2mean - xmean*xmean;
    yrms2  = y2mean - ymean*ymean;
    xrms   = ( xrms2>0.0 ) ? TMath::Sqrt(xrms2) : 0.0 ;
    yrms   = ( yrms2>0.0 ) ? TMath::Sqrt(yrms2) : 0.0 ;
    xyrms2 = xymean - xmean*ymean;
    emittance2 = xrms2*yrms2 - xyrms2*xyrms2;
    emittance = (emittance2>0.0) ? TMath::Sqrt(emittance2) : 0.0 ;

    // cout << Form("  xMean = %7.3f   yMean = %7.3f",xmean,ymean) << endl;
    // cout << Form("  xRms  = %7.3f   yRms  = %7.3f",xrms,yrms) << endl;
    // cout << Form("  Emittance = %7.3f",emittance) << endl;

    Float_t emity = emittance;
    Float_t y_mean = xmean;
    Float_t y_rms = xrms;
    Float_t pymean = ymean;
    Float_t pyrms = yrms;


    // Charge  
    hX1->Scale(dx1*dx2*dx3);
    Float_t Charge = hX1->Integral();

    cout << Form("\n 4. Calculating sliced quantities.. ") << endl ;

    TGraph *gemit = NULL;
    TGraph *gYrms = NULL;
    TGraph *gErms = NULL;

    Float_t * sxmean = new Float_t[SNbin];
    Float_t * symean = new Float_t[SNbin];
    Float_t * sx2mean = new Float_t[SNbin];
    Float_t * sy2mean = new Float_t[SNbin];
    Float_t * sxymean = new Float_t[SNbin];
    Float_t * sNtotal = new Float_t[SNbin];
    Float_t * sxrms2 = new Float_t[SNbin];  
    Float_t * syrms2 = new Float_t[SNbin]; 
    Float_t * sxrms = new Float_t[SNbin];  
    Float_t * syrms = new Float_t[SNbin];  
    Float_t * sxyrms2 = new Float_t[SNbin];

    Float_t * xbin = new Float_t[SNbin];
    Float_t * semit = new Float_t[SNbin];
    Float_t * semit2 = new Float_t[SNbin];

    Float_t * sNEtotal = new Float_t[SNbin]; 
    Float_t * sEmean = new Float_t[SNbin];

    Float_t * sErms = new Float_t[SNbin];
    Float_t * sErms2 = new Float_t[SNbin];

    for(Int_t k=0;k<SNbin;k++) {
      sxmean[k] = symean[k] = sx2mean[k] = sy2mean[k] = sxymean[k] 
	= sNtotal[k] = sxrms2[k] = syrms2[k] = sxrms[k] = syrms[k]
	= sxyrms2[k] = xbin[k] = semit[k] = 0.0;
      sNEtotal[k] = sEmean[k] = sErms[k] = 0.0;

      xbin[k] = (sBinLim[k] + sBinLim[k+1])/2.;

      for(Int_t i=1;i<=x2Nbin;i++) {
	Float_t x = hP2X2sl[k]->GetXaxis()->GetBinCenter(i);
	// if(x<xmin || x>xmax) continue;
	for(Int_t j=1;j<=p2Nbin;j++) {
	  Float_t y = hP2X2sl[k]->GetYaxis()->GetBinCenter(j);
	  // if(y<ymin || y>ymax) continue;
	  Float_t value = TMath::Abs(hP2X2sl[k]->GetBinContent(i,j));
	  sxmean[k] += x*value;
	  symean[k] += y*value;
	  sx2mean[k] += x*x*value;
	  sy2mean[k] += y*y*value;
	  sxymean[k] += x*y*value;

	  sNtotal[k] += value;
	}	
      }

      if( sNtotal[k] == 0.0 ) continue;

      // for(Int_t i=1;i<=p1Nbin;i++) {
      // 	Float_t y = hP1sl[k]->GetXaxis()->GetBinCenter(i);
      // 	Float_t value = TMath::Abs(hP1sl[k]->GetBinContent(i));
      // 	sEmean[k] += y*value;
      // 	sE2mean[k] += y*y*value;
      // 	sNEtotal[k] += value;
      // }

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
      sErms[k]   = hP1sl[k]->GetRMS();
      
      // cout<< Form("\nk = %i : (x1 > %f && x1 < %f)",k,sBinLim[k],sBinLim[k+1]) << endl; 

      // cout << Form("  xMean = %7.3f   yMean = %7.3f",sxmean[k],symean[k]) << endl;
      // cout << Form("  xRms  = %7.3f   yRms  = %7.3f",sxrms[k],syrms[k]) << endl;
      // cout << Form("  Emittance = %7.3f",semit[k]) << endl;

      // cout << Form("  Emean = %7.3f   Erms = %7.3f",sEmean[k],sErms[k]) << endl;


    }
    cout << " done! " << endl;


    // Changing to user units: 
    // --------------------------

    if(opt.Contains("units") && n0) {

      Time *= skindepth / PUnits::mm;

      x1Min *= skindepth / PUnits::um;
      x1Max *= skindepth / PUnits::um;
      p1Min *= PConst::ElectronMassE / PUnits::GeV;
      p1Max *= PConst::ElectronMassE / PUnits::GeV;

      hP1X1->SetBins(x1Nbin,x1Min,x1Max,p1Nbin,p1Min,p1Max);

      // Converting electron density
      Float_t dV = skindepth * skindepth * skindepth;
      Charge *= n0 * dV * (PConst::ElectronCharge/PUnits::picocoulomb);

      Float_t dX = (x1Max-x1Min)/x1Nbin; 
      Float_t dE = (p1Max-p1Min)/p1Nbin; 
      for(Int_t j=0;j<hP1X1->GetNbinsX();j++) {
	for(Int_t k=0;k<hP1X1->GetNbinsY();k++) {
	  Float_t binValue =  fabs(hP1X1->GetBinContent(j,k) * dx1 * dx2 * dx3 * dV * n0 *
				   (PConst::ElectronCharge/PUnits::picocoulomb));
	  //cout << Form(" value = %f",binValue) << endl;
	  hP1X1->SetBinContent(j,k,binValue);

	}
      }

      if(opt.Contains("comov"))
	hP1X1->GetXaxis()->SetTitle("#zeta [#mum]");
      else
	hP1X1->GetXaxis()->SetTitle("z [#mum]");

      hP1X1->GetYaxis()->SetTitle("p_{z} [GeV/c]");

      hP1X1->GetZaxis()->SetTitle("Charge [pC]");
      hP1X1->GetZaxis()->CenterTitle();

      hP1->SetBins(p1Nbin,p1Min,p1Max);
      hP1->GetYaxis()->SetTitle("p_{z} [GeV/c]");

      hX1->SetBins(x1Nbin,x1Min,x1Max);
      Float_t binSize = (x1Max - x1Min)/x1Nbin;

      Float_t  lightspeed =  PConst::c_light / (PUnits::um/PUnits::femtosecond);
      hX1->Scale(TMath::Abs(n0 * dV * (PConst::ElectronCharge/PUnits::picocoulomb) * (lightspeed/binSize)));

      // hX1->Scale(TMath::Abs((PUnits::um/skindepth)*(PConst::ElectronCharge/PUnits::picocoulomb)*PConst::c_light));

      // hX1->GetYaxis()->SetTitle("I[kA]");
      hX1->GetYaxis()->SetTitle("");
      if(opt.Contains("comov"))
	hX1->GetXaxis()->SetTitle("#zeta [#mum]");
      else
	hX1->GetXaxis()->SetTitle("z [#mum]");


      zmean *= skindepth / PUnits::um;
      zrms  *= skindepth / PUnits::um;

      y_mean *= skindepth / PUnits::um;
      y_rms  *= skindepth / PUnits::um;

      xmean *= skindepth / PUnits::um;
      xrms  *= skindepth / PUnits::um;

      pzmean *= PConst::ElectronMassE / PUnits::GeV;
      pzrms  *= PConst::ElectronMassE / PUnits::GeV;

      pymean *= PConst::ElectronMassE / PUnits::MeV;
      pyrms  *= PConst::ElectronMassE / PUnits::MeV;

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

      hP2X2->GetXaxis()->SetTitle("y [#mum]");
      hP2X2->GetYaxis()->SetTitle("p_{y} [MeV/c]");
      hP2X2->GetZaxis()->SetTitle("dQ/dydp_{y} [pC]");
      hP2X2->GetZaxis()->CenterTitle();

      emity *= (skindepth / PUnits::um);
 
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

	hP2X2sl[k]->GetZaxis()->SetTitle("dQ/dydp_{y} [a.u.]");
	hP2X2sl[k]->GetXaxis()->SetTitle("y [#mum]");
	hP2X2sl[k]->GetYaxis()->SetTitle("p_{y} [MeV/c]");

	hP2X2sl[k]->GetZaxis()->CenterTitle();

	xbin[k] *= skindepth / PUnits::um;

	sxmean[k] *= skindepth / PUnits::um;
	sxrms[k]  *= skindepth / PUnits::um;
	symean[k] *= PConst::ElectronMassE / PUnits::MeV;
	syrms[k] *= PConst::ElectronMassE / PUnits::MeV;

	semit[k] *= (skindepth / PUnits::um);

	sEmean[k] *= PConst::ElectronMassE / PUnits::GeV;
	sErms[k]  *= 100 * PConst::ElectronMassE / PUnits::GeV / pzmean; //sEmean[k];
	// sErms[k]  *= PConst::ElectronMassE / PUnits::GeV;

      }
    }
    // End of the users units module    

    cout << "\n  Summary _______________________________________________________ " << endl;
    if(opt.Contains("units")) {
      cout << Form("  Integrated charge (RAW) of specie %3i = %8f pC",index,Charge) << endl;
      cout << Form("  Peak current = %6.2f kA",hX1->GetMaximum()) << endl;
      cout << Form("  Total energy = %6.2f GeV, rms = %3.1f %%",pzmean,100*pzrms/pzmean) << endl;
      cout << Form("  Trans. emit. = %6.2f  um",emity) << endl;
      cout << Form("  Bunch length = %6.2f  um (rms), width = %6.2f  um (rms)",zrms,y_rms) << endl;
    } else {
      cout << Form("  Integrated charge (RAW) of specie %3i = %8.4f n0 * kp^-3",index,Charge) << endl;
    }


    if(opt.Contains("loop")) {
      cout << Form("\n 5. Saving results to file .. ") << endl;
  
      // OUTPUT ROOT FILE WITH THE PLOTS:
      TString filename = Form("./%s/Plots/Bunch/%s/Bunch-Evolution-%s.root",sim.Data(),specName.Data(),sim.Data());
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
      sprintf(gName,"gEmitvsTime");     
      TGraph *gEmitvsTime = NULL;
      gEmitvsTime = (TGraph*) ifile->Get(gName);
      if(gEmitvsTime==NULL) {
	gEmitvsTime = new TGraph();
	gEmitvsTime->SetName(gName);
	nPoints = 0;
	// Some cosmetics at creation time:
	gEmitvsTime->SetLineWidth(3);
	gEmitvsTime->SetLineColor(PGlobals::fieldLine);
	gEmitvsTime->SetMarkerStyle(20);
	gEmitvsTime->SetMarkerSize(0.4);
	gEmitvsTime->SetMarkerColor(PGlobals::fieldLine);	
      } else {
	nPoints = gEmitvsTime->GetN(); 
      }  

      gEmitvsTime->Set(nPoints+1);
      gEmitvsTime->SetPoint(nPoints,Time,emity);
      gEmitvsTime->Write(gName,TObject::kOverwrite);


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
      gPymeanvsTime->SetPoint(nPoints,Time,pymean);
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
      gPyrmsvsTime->SetPoint(nPoints,Time,pyrms);
      gPyrmsvsTime->Write(gName,TObject::kOverwrite);



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
    // for(UInt_t i=0;i<Nvar;i++) {
    //   delete var[i];
    // }
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
    gemit = new TGraph(SNbin,xbin,semit);
    gYrms = new TGraph(SNbin,xbin,sxrms);
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
      Float_t xMax = hX1->GetXaxis()->GetXmin() + (hX1->GetXaxis()->GetXmax()
						   -hX1->GetXaxis()->GetXmin()) * 0.2;
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
      //      }



      Bool_t autoscale = kFALSE;
      // if(hX1->GetMaximum()>10) {
      // 	hX1->Scale(0.1);
      // 	autoscale = kTRUE;
      // }
      
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
      if(opt.Contains("units") && n0) 
	sprintf(ctext,"z = %5.1f mm", Time);
      else
	sprintf(ctext,"t = %5.1f #omega_{p}^{-1}",Time);
      textTime->AddText(ctext);
 
      TPaveText *textDen = new TPaveText(0.15,0.85,0.48,0.9,"NDC");
      PGlobals::SetPaveTextStyle(textDen,12); 
      textDen->SetTextColor(kOrange+10);
      if(opt.Contains("units") && n0)
	sprintf(ctext,"n_{0} = %5.2f x 10^{17} / cc", n0 / (1e17/PUnits::cm3));
      textDen->AddText(ctext);
      
      TPaveText *textCharge = new TPaveText(0.15,0.25,0.48,0.3,"NDC");
      PGlobals::SetPaveTextStyle(textCharge,12); 
      textCharge->SetTextColor(kGray+2);
      if(opt.Contains("units") && n0)
	sprintf(ctext,"Q = %5.2f pC", Charge);
      else
	sprintf(ctext,"Q = %5.2f n0#timeskp^{-3}", Charge);    
      textCharge->AddText(ctext);

      TPaveText *textMom = new TPaveText(0.55,0.03,0.80,0.13,"NDC");
      PGlobals::SetPaveTextStyle(textMom,32); 
      textMom->SetTextColor(kGray+3);
      textMom->SetTextFont(62);
      if(opt.Contains("units") && n0)
	sprintf(ctext,"#LTp_{z}#GT = %5.2f GeV/c", pzmean);
      else
	sprintf(ctext,"#LTp_{z}#GT = %5.2f mc", pzmean);    
      textMom->AddText(ctext);


      TPaveText *textInfo = new TPaveText(0.55,0.50,0.80,0.75,"NDC");
      PGlobals::SetPaveTextStyle(textInfo,32); 
      textInfo->SetTextColor(kGray+2);
      textInfo->SetTextFont(42);
      sprintf(ctext,"Charge = %5.2f pC",Charge);
      textInfo->AddText(ctext);
      sprintf(ctext,"#Delta#zeta = %5.2f #mum",zrms);
      textInfo->AddText(ctext);
      //     sprintf(ctext,"#LTp_{z}#GT_{rms} = %5.2f GeV/c",pzrms);
      sprintf(ctext,"#Deltap_{z}/p_{z} = %4.1f %%",100*pzrms/pzmean);
      textInfo->AddText(ctext);
      sprintf(ctext,"#epsilon_{y} = %5.2f #mum",emity);
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
      PGlobals::CanvasAsymPartition(C,NPad,lMargin,rMargin,bMargin,tMargin,factor);

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
	pad[i]->SetTickx(1);
	pad[i]->SetTicky(1);

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

      TLine lPmean(hP1X1->GetXaxis()->GetXmin(),pzmean,hP1X1->GetXaxis()->GetXmax(),pzmean);
      lPmean.SetLineColor(kGray+2);
      lPmean.SetLineStyle(2);
      lPmean.Draw();

      // 2D histogram z range
      Float_t dmax = hP1X1->GetMaximum();
      Float_t dmin = 0.0;
      hP1X1->GetZaxis()->SetRangeUser(1.1*dmin,dmax);

      hP1X1->GetZaxis()->SetTitleFont(fonttype);

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
      Leg->AddEntry(gErms,"Energy spread [%]","PL");
      Leg->AddEntry(gemit,"Emittance [#mum]","PL");
      //Leg->AddEntry(gYrms,"Bunch width [#mum]","PL");


      hFrame[0]->GetYaxis()->SetRangeUser(0.0,1.1*yMax);
      hFrame[0]->Draw();

      hX1->GetYaxis()->SetNdivisions(503);
      hX1->SetLineWidth(2);
      hX1->SetFillStyle(1001);
      hX1->SetFillColor(PGlobals::elecFill);
      // hX1->SetLineColor(kBlue);

      //hX1->Smooth();
      hX1->Draw("FL same");
      //hX1->Draw("C");

      TLine lZmean2(zmean,0.0,zmean,1.1*yMax);
      lZmean2.SetLineColor(kGray+2);
      lZmean2.SetLineStyle(2);
      lZmean2.Draw();

      Int_t markerSize = 1.2; 
      Int_t lineWidth  = 2.0;   

      gYrms->SetMarkerStyle(20);
      gYrms->SetLineStyle(1);
      gYrms->SetMarkerColor(kGray+1);
      gYrms->SetMarkerSize(markerSize); 
      gYrms->SetLineColor(kGray+1);
      gYrms->SetLineWidth(lineWidth);
      // gYrms->Draw("PL");

      // hP2X2sl[0]->Draw("colz");
      gemit->SetMarkerStyle(20);
      gemit->SetMarkerColor(kGray+2);
      gemit->SetMarkerSize(markerSize);
      gemit->SetLineWidth(lineWidth);
      gemit->SetLineColor(kGray+2);
      gemit->Draw("PL");

      gErms->SetMarkerStyle(20);
      gErms->SetMarkerSize(markerSize);
      gErms->SetMarkerColor(kOrange+10);
      gErms->SetLineColor(kOrange+10);
      gErms->SetLineWidth(lineWidth);
      gErms->Draw("PL");


   
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
      TString fOutName = Form("./%s/Plots/Bunch/%s/Bunch",sim.Data(),specName.Data());
      fOutName += Form("-%s_%.2f",sim.Data(),time);
      
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

      TString fOutName2 = Form("./%s/Plots/Bunch/%s/Bunch-%s-p2x2_%.2f",sim.Data(),specName.Data(),sim.Data(),time);

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
      sprintf(text,"#epsilon_{y} = %5.2f",emity);
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

	hP2X2sl[k]->Draw("colz");

	Float_t y1 = gPad->GetBottomMargin();
	Float_t y2 = 1 - gPad->GetTopMargin();
	Float_t x1 = gPad->GetLeftMargin();
	Float_t x2 = 1 - gPad->GetRightMargin();
	textStat[k] = new TPaveText(x1+0.02,y2-0.40,x1+0.20,y2-0.05,"NDC");
	PGlobals::SetPaveTextStyle(textStat[k],12); 
	textStat[k]->SetTextColor(kGray+3);
	textStat[k]->SetTextFont(42);

	char text[64];
	sprintf(text,"#LTy#GT_{rms} = %5.2f",sxrms[k]);
	textStat[k]->AddText(text);
	sprintf(text,"#LTp_{y}#GT_{rms} = %5.2f",syrms[k]);
	textStat[k]->AddText(text);
	sprintf(text,"#epsilon_{y} = %5.2f",semit[k]);
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
      TString filename = Form("./%s/Plots/Bunch/%s/Bunch-%s_%.2f.root",sim.Data(),specName.Data(),sim.Data(),time);
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
      gYrms->Write("gYrms",TObject::kOverwrite);
      gemit->Write("gEmit",TObject::kOverwrite);

      for(Int_t k=0;k<SNbin;k++) {

	sprintf(hName,"hP2X2sl_%i",k);	
	hP2X2sl[k]->Write(hName,TObject::kOverwrite);

	sprintf(hName,"hP1sl_%i",k);	
	hP1sl[k]->Write(hName,TObject::kOverwrite);
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

    // delete gemit;
    // delete gYrms;
    // delete gErms;
    // delete gErmsB;

    // end time looper
  }
}
