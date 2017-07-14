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
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TPaletteAxis.h>
#include <TExec.h>
#include <TRandom3.h>

#include "PData.hh"
#include "PDataHiP.hh"
#include "PGlobals.hh"
#include "PPalette.hh"
#include "H5Cpp.h"

using namespace std;
using namespace H5;


int main(int argc,char *argv[]) {
  if(argc<=2) {
    printf("\n Usage: %s <simulation> <-t(time)>\n",argv[0]);
    printf("      <-index(species index)>\n");
    printf("      <-i(initial time)> <-f(final time)> <-s(time step)>\n");
    printf("      <--ast> <--ele> <--hip>\n");
    printf("      <--comov> <--center> <--eqw> \n");
    printf("      <--png> <--pdf> <--eps>\n");
    return 0;
  }

  // General options
  TString     sim = "";
  TString fileout = "";
  Int_t      time = 0;
  Int_t    iStart = -1;
  Int_t      iEnd = -1;
  Int_t     iStep = 1;
  TString     opt = "";
  Int_t     index = 1;
  Bool_t    noout = kTRUE;
 
  // Interfacing command line:
  for(int l=1;l<argc;l++){
    TString arg = argv[l];

    if(l==1) {
      sim = arg;
    } else if (l==2 && !noout) {
      fileout = arg;
      if(fileout.Contains(".raw")) 
	opt += "osiris";
      else if(fileout.Contains(".ast")) 
	opt += "astra";
      else if(fileout.Contains(".ele")) 
	opt += "elegant";
      else if(fileout.Contains(".hip")) 
	opt += "hipace";
      else {
	cout << Form(" No output file specified or wrong extension");
	fileout = "";
	l--;
	noout = kTRUE;
	continue;
      }
    } else if(arg.Contains("--ast")) {
      opt += "astra";
    } else if(arg.Contains("--ele")) {
      opt += "elegant";
    } else if(arg.Contains("--hip")) {
      opt += "hipace";
    } else if(arg.Contains("--pdf")) {
      opt += "plotpdf";
    } else if(arg.Contains("--eps")) {
      opt += "ploteps";
    } else if(arg.Contains("--png")) {
      opt += "plotpng";
    } else if(arg.Contains("--comov")){
      opt += "comov";
    } else if(arg.Contains("--center")){
      opt += "center";
    } else if(arg.Contains("--units")){
      opt += "units";
    } else if(arg.Contains("--zmean")){
      opt += "zmean";
    } else if(arg.Contains("--grid")){
      opt += "grid"; 
    } else if(arg.Contains("--logz")){
      opt += "logz"; 
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
    } else if(arg.Contains("--eqw")) {
      opt += "eqw";
    }

    
  }

  PGlobals::Initialize();

  // Palettes!
  // gROOT->Macro("PPalettes.C");

  if(opt.Contains("grid")) {
    gStyle->SetPadGridX(1);
    gStyle->SetPadGridY(1);
  }
  gStyle->SetPadRightMargin(0.17);   

  // Load PData
  PData *pData = PData::Get(sim.Data());
  if(pData->isHiPACE()) {
    delete pData; pData = NULL;
    pData = PDataHiP::Get(sim.Data());
    if(!opt.Contains("comov"))
      opt += "comov";
  }
  
  if(iStart<0) iStart = time;
  if(iEnd<=iStart) iEnd = iStart;
  
  // Some plasma constants
  Float_t n0 = pData->GetPlasmaDensity();
  Float_t kp = pData->GetPlasmaK();
  Float_t skindepth = 1.;
  if(kp!=0.0) skindepth = 1/kp;

  // z start of the plasma in normalized units.
  Float_t zStartPlasma = pData->GetPlasmaStart()*kp;
  // z start of the beam in normalized units.
  Float_t zStartBeam = pData->GetBeamStart()*kp;

  // Time looper
  for(Int_t i=iStart; i<=iEnd; i+=iStep) {

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
      if(opt.Contains("comov"))      // Centers on the head of the beam.
	Time += zStartBeam;
    } 

    // Spatial resolution
    Float_t dx1 = pData->GetDX(0);
    Float_t dx2 = pData->GetDX(1);
    Float_t dx3 = pData->GetDX(2);
  
    // --------------------------------------------------
    cout << Form("\n 1. Reading file : %s \n",pData->GetRawFileName(index)->c_str()) << endl;
    cout << Form("    with options = %s\n",opt.Data()) << endl;
    
    // Float_t **var = NULL;
    UInt_t Nvar = 7;
    if(!pData->Is3D()) Nvar = 6;
    Float_t **var;
    var = new Float_t*[Nvar];
    //  char varname[Nvar][4] = {{"p1"},{"p2"},{"p3"},{"q"},{"x1"},{"x2"},{"x3"}};
    UInt_t Np = pData->GetRawArray(pData->GetRawFileName(index)->c_str(),var);  
    UInt_t Npout = 0;

    cout << Form(" Number of macroparticles = %i",Np) << endl;
    
    // -------------------------------------------------------------------------------

    // First loop to analyze the data sample:
    // Get the maximums and minimums in every variable.
    
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

    Float_t MinQ = 999999;
    Float_t MaxQ = -999999;

    Float_t x1Mean = 0.0;
    Float_t qSum  = 0.0;
    
    for(UInt_t i=0;i<Np;i++) {
      if(var[0][i]<MinP1) MinP1 = var[0][i];
      if(var[0][i]>MaxP1) MaxP1 = var[0][i];
      if(var[1][i]<MinP2) MinP2 = var[1][i];
      if(var[1][i]>MaxP2) MaxP2 = var[1][i];
      if(var[2][i]<MinP3) MinP3 = var[2][i];
      if(var[2][i]>MaxP3) MaxP3 = var[2][i];
      if(fabs(var[3][i])<MinQ) MinQ = fabs(var[3][i]);
      if(fabs(var[3][i])>MaxQ) MaxQ = fabs(var[3][i]);
      if((var[4][i]-shiftz)<MinX1) MinX1 = var[4][i]-shiftz;
      if((var[4][i]-shiftz)>MaxX1) MaxX1 = var[4][i]-shiftz;
      if(var[5][i]<MinX2) MinX2 = var[5][i];
      if(var[5][i]>MaxX2) MaxX2 = var[5][i];
      if(var[6][i]<MinX3) MinX3 = var[6][i];
      if(var[6][i]>MaxX3) MaxX3 = var[6][i];
    
      x1Mean += (var[4][i]-shiftz) * fabs(var[3][i]);
      qSum  += fabs(var[3][i]);
    }
    x1Mean /= qSum;
    
    // Test histograms
    TH2F *hP1X1 = NULL;
    TH2F *hP2X2 = NULL;
    TH2F *hP1X1out = NULL;
    TH2F *hP2X2out = NULL;
    
    if(opt.Contains("plot")) {

      // Longitudinal phasespace
      Float_t x1Min = MinX1 - 0.3*(MaxX1-MinX1);
      x1Min = floor((x1Min-pData->GetXMin(0))/dx1) * dx1 + pData->GetXMin(0);
      Float_t x1Max = MaxX1 + 0.3*(MaxX1-MinX1);
      x1Max = floor((x1Max-pData->GetXMin(0))/dx1) * dx1 + pData->GetXMin(0);
      //    UInt_t x1Nbin = ceil ((x1Max - x1Min)/(dx1));
      UInt_t x1Nbin = 100;
      Float_t p1Min = MinP1 - 0.3*(MaxP1-MinP1);
      Float_t p1Max = MaxP1 + 0.3*(MaxP1-MinP1);
      UInt_t p1Nbin = x1Nbin;

      char hName[16];
      sprintf(hName,"hP1X1");
      hP1X1 = (TH2F*) gROOT->FindObject(hName);
      if(hP1X1) delete hP1X1;
      hP1X1 = new TH2F(hName,"",x1Nbin,x1Min,x1Max,p1Nbin,p1Min,p1Max);
   
      sprintf(hName,"hP1X1out");
      hP1X1out = (TH2F*) gROOT->FindObject(hName);
      if(hP1X1out) delete hP1X1out;
      hP1X1out = new TH2F(hName,"",x1Nbin,x1Min,x1Max,p1Nbin,p1Min,p1Max);
    
      cout << Form(" Test histograms created : %i\n",hP1X1 && hP1X1out);
      cout << Form(" x1 range (N = %i):  x1Min = %f  x1Max = %f ", x1Nbin, x1Min, x1Max) << endl;
      cout << Form(" p1 range (N = %i):  p1Min = %f  p1Max = %f ", p1Nbin, p1Min, p1Max) << endl;

      // Transverse phasespace
      Float_t x2Min = MinX2 - 0.3*(MaxX2-MinX2);
      x2Min = floor((x2Min-pData->GetXMin(1))/dx2) * dx2 + pData->GetXMin(1);
      Float_t x2Max = MaxX2 + 0.3*(MaxX2-MinX2);
      x2Max = floor((x2Max-pData->GetXMin(1))/dx2) * dx2 + pData->GetXMin(1);
      //    UInt_t x2Nbin = ceil ((x2Max - x2Min)/(dx2));
      UInt_t x2Nbin = 100;
      Float_t p2Min = MinP2 - 0.3*(MaxP2-MinP2);
      Float_t p2Max = MaxP2 + 0.3*(MaxP2-MinP2);
      UInt_t p2Nbin = x2Nbin;

      sprintf(hName,"hP2X2");
      hP2X2 = (TH2F*) gROOT->FindObject(hName);
      if(hP2X2) delete hP2X2;
      hP2X2 = new TH2F(hName,"",x2Nbin,x2Min,x2Max,p2Nbin,p2Min,p2Max);
   
      sprintf(hName,"hP2X2out");
      hP2X2out = (TH2F*) gROOT->FindObject(hName);
      if(hP2X2out) delete hP2X2out;
      hP2X2out = new TH2F(hName,"",x2Nbin,x2Min,x2Max,p2Nbin,p2Min,p2Max);
    
      cout << Form(" Test histograms created : %i\n",hP2X2 && hP2X2out);
      cout << Form(" x2 range (N = %i):  x2Min = %f  x2Max = %f ", x2Nbin, x2Min, x2Max) << endl;
      cout << Form(" p2 range (N = %i):  p2Min = %f  p2Max = %f ", p2Nbin, p2Min, p2Max) << endl;

    }

   
    // Output file
    if(fileout.IsNull()) {
      fileout = Form("./%s/Plots/RawConvert/RawConvert-%s_%s_%i",sim.Data(),sim.Data(),pData->GetRawSpeciesName(index).c_str(),time);
      if(opt.Contains("astra")) fileout += ".ast";
      else if(opt.Contains("elegant")) fileout += ".ele";
      else if(opt.Contains("hipace")) fileout += ".hip";
      else fileout += ".raw";
    }

    cout << Form("\n 2. Output file : %s \n",fileout.Data()) << endl;
     
    ofstream fData;
    if(opt.Contains("hipace"))
      fData.open(fileout.Data(),ios::binary | ios::out);
    else
      fData.open(fileout.Data(),ios::out);
      
    // if doesn't exist the directory should be created
    if (!fData) {
      TString f = fileout;
      TString dir2 = f.Remove( f.Last( '/' ), f.Length() - f.Last( '/' ) );
      TString dir1 = f.Remove( f.Last( '/' ), f.Length() - f.Last( '/' ) );
      gSystem->mkdir( dir1 );
      gSystem->mkdir( dir2 );
      if(opt.Contains("hipace"))
	fData.open(fileout.Data(),ios::binary | ios::out);
      else
	fData.open(fileout.Data(),ios::out);
 
    }  
    
    Double_t dV = (dx1*skindepth)*(dx2*skindepth)*(dx3*skindepth);
    Double_t Q0 = n0 * dV; // Charge normalization (from Osiris)
    Double_t QTotal = 0.00;
    Double_t QTotalOut = 0.00;

    cout << Form("\n Charge normalization: Q0 = %f pC",Q0 * PConst::ElectronCharge/PUnits::picocoulomb) << endl;

    Double_t buffer[7];
    //  cout << Form("Macroparticle buffer size %5i",(int)sizeof(buffer)) << endl;

    // Main loop for randomizing and converting
    TRandom3 *ran = new TRandom3(0);
    for(UInt_t i=0;i<Np;i++) {
   
      if( (i % 1000000) == 0 && i!=0 )
	cout << Form(" %8i macroparticles processed ",i) << endl;

      // z coordinate shift according with the options
      var[4][i] -= shiftz;
      
      QTotal += var[3][i] * Q0 * PConst::ElectronCharge/PUnits::picocoulomb;

      // Fills histogram from original distribution:
      if(hP1X1)
	hP1X1->Fill(var[4][i],var[0][i],TMath::Abs(var[3][i]));
      if(hP2X2)
	hP2X2->Fill(var[5][i],var[1][i],TMath::Abs(var[3][i]));

      if(opt.Contains("eqw")) {
	// Macro-particle charge equalizer:
	Float_t rw = ran->Uniform(MaxQ);  // Random number between 0. and "q".
	// cout << Form(" wi = %f    rw = %f   MaxQ = %f ",var[3][i],rw,MaxQ) << endl;
	if(rw>fabs(var[3][i])) continue;  
            
	var[3][i] = -MaxQ;
      }

      QTotalOut += var[3][i] * Q0 * PConst::ElectronCharge/PUnits::picocoulomb;
    
      // Fills histogram from new distribution:
      if(hP1X1out)
	hP1X1out->Fill(var[4][i],var[0][i],TMath::Abs(var[3][i]));
      if(hP2X2out)
	hP2X2out->Fill(var[5][i],var[1][i],TMath::Abs(var[3][i]));
      
      // Dumping to ASCII file:
      if(opt.Contains("astra")) {

	if(i>0) {
	  var[4][i] -= var[4][0];
	  var[0][i] -= var[0][0];
	}
	
	fData << Form("%12.6e   %12.6e   %12.6e   %12.6e   %12.6e   %12.6e   %12.6e   %12.6e   %i   %i",
		      // var[6][i] * skindepth / PUnits::um,
		      // var[5][i] * skindepth / PUnits::um,
		      // var[4][i] * skindepth / PUnits::um,
		      // var[2][i] * pData->GetBeamMass() / PUnits::MeV,
		      // var[1][i] * pData->GetBeamMass() / PUnits::MeV,
		      // var[0][i] * pData->GetBeamMass() / PUnits::MeV,
		      // -999.0,
		      // var[3][i] * Q0 * PConst::ElectronCharge/PUnits::femptocoulomb,
		      // 0.0) << endl;
		      var[6][i] * skindepth / PUnits::m,
		      var[5][i] * skindepth / PUnits::m,
		      var[4][i] * skindepth / PUnits::m,
		      var[2][i] * pData->GetBeamMass() / PUnits::eV,
		      var[1][i] * pData->GetBeamMass() / PUnits::eV,
		      var[0][i] * pData->GetBeamMass() / PUnits::eV,
		      0.0,
		      var[3][i] * Q0 * PConst::ElectronCharge/PUnits::nanocoulomb,
		      1,
		      5) << endl;
	
	
    
      } else if(opt.Contains("elegant")) {
	// fData << Form("%e   %e   %e   %e   %e   %e   %e",
	// 	      var[5][i] * skindepth / PUnits::m,
	// 	      var[1][i]/var[0][i],
	// 	      var[6][i] * skindepth / PUnits::m,
	// 	      var[2][i]/var[0][i],
	// 	      (var[4][i] * skindepth / PConst::c_light) / PUnits::second,
	// 	      var[0][i],
	// 	      var[3][i] * Q0 * PConst::ElectronCharge/PUnits::coulomb) << endl;
	fData << Form("%14e   %14e   %14e   %14e   %14e   %14e",
		      var[5][i] * skindepth / PUnits::m,
		      var[1][i]/var[0][i],
		      var[6][i] * skindepth / PUnits::m,
		      var[2][i]/var[0][i],
		      (var[4][i] * skindepth / PConst::c_light) / PUnits::second,
		      var[0][i]) << endl;
      
      } else if(opt.Contains("hipace")) {
	buffer[0] = (double) var[4][i];
	buffer[1] = (double) var[5][i];
	buffer[2] = (double) var[6][i];
	buffer[3] = (double) var[0][i];
	buffer[4] = (double) var[1][i];
	buffer[5] = (double) var[2][i];
	buffer[6] = (double) var[3][i];


	fData.write((char*)&buffer,sizeof(buffer));
	
      } else {
	// if(Npout ==0 ) {
	//   // Header
	//   for(UInt_t j=0;j<Nvar;j++) {
	//     fData << Form("%10s ",varname[j]); 
	//   }
	//   fData << endl << Form(" --------------------------------------------------------------------------------------- ") << endl;
	// }
	// for(UInt_t j=0;j<Nvar;j++) {
	//   fData << Form("%10.4f ",var[j][i]);
	// }

	// if(opt.Contains("zmean"))
	//   var[4][i] -= zMean;
	   
	fData << Form("%14e   %14e   %14e   %14e   %14e   %14e   %i",
		      var[4][i] * skindepth / PUnits::m,
		      var[5][i] * skindepth / PUnits::m,
		      var[6][i] * skindepth / PUnits::m,
		      var[0][i],
		      var[1][i],
		      var[2][i],
		      -1) << endl;
	//		      var[3][i]) << endl;
	
	//	fData << endl;
	
      }
      
      Npout++;
    }


    cout << endl;
    
    fData.close();
    

    cout << "\n End of data writing " << endl;
    
    // Plotting Section
    // -----------------------------------------------
    
    if(hP1X1 && hP1X1out) {

      cout << Form(" 3. Plotting output: \n");

      // Canvas setup
      TCanvas *C;
      C = new TCanvas("C","Phasespaces",1280,800);
    
      // Histograms
      hP1X1->GetXaxis()->CenterTitle();
      hP1X1->GetYaxis()->CenterTitle();
      hP1X1->GetZaxis()->CenterTitle();

      hP1X1->GetXaxis()->SetTitle("#zeta [k_{p}^{-1}]");
      hP1X1->GetYaxis()->SetTitle("p_{z} [mc]");
      hP1X1->GetZaxis()->SetTitle("charge");
  
      UInt_t entries = hP1X1->GetEntries();
      Float_t integral = hP1X1->Integral();
      Float_t zMean = hP1X1->GetMean(1);
      Float_t zRms  = hP1X1->GetRMS(1);
      Float_t pzMean = hP1X1->GetMean(2);
      Float_t pzRms  = hP1X1->GetRMS(2);

      cout << Form("\n Original distribution: ") << endl; 
      cout << Form(" N = %10i  Integral = %7.2f  Mean = %7.2f  Rms = %f",Np,integral,zMean,zRms) << endl;
  
      // Info for the plots
      char ctext[128];
      TPaveText *textInfo = new TPaveText(0.55,0.60,0.82,0.85,"NDC");
      PGlobals::SetPaveTextStyle(textInfo,32); 
      textInfo->SetTextColor(kGray+2);
      textInfo->SetTextFont(42);
      sprintf(ctext,"Entries = %5i",entries);
      textInfo->AddText(ctext);
      sprintf(ctext,"Integral = %5.2f",integral);
      textInfo->AddText(ctext);
      sprintf(ctext,"#LT#zeta#GT = %5.2f k_{p}^{-1}",zMean);
      textInfo->AddText(ctext);
      sprintf(ctext,"#LT#zeta#GT_{rms} = %5.2f k_{p}^{-1}",zRms);
      textInfo->AddText(ctext);
      sprintf(ctext,"#LTp_{z}#GT = %5.2f mc",pzMean);
      textInfo->AddText(ctext);
      sprintf(ctext,"#LTp_{z}#GT_{rms} = %5.2f mc",pzRms);
      textInfo->AddText(ctext);

      // Transverse phasespace histogram
      hP2X2->GetXaxis()->CenterTitle();
      hP2X2->GetYaxis()->CenterTitle();
      hP2X2->GetZaxis()->CenterTitle();

      hP2X2->GetXaxis()->SetTitle("x [k_{p}^{-1}]");
      hP2X2->GetYaxis()->SetTitle("p_{x} [mc]");
      hP2X2->GetZaxis()->SetTitle("charge");
  
      integral = hP2X2->Integral();
      entries = hP2X2->GetEntries();
      Float_t xMean = hP2X2->GetMean(1);
      Float_t xRms  = hP2X2->GetRMS(1);
      Float_t pxMean = hP2X2->GetMean(2);
      Float_t pxRms  = hP2X2->GetRMS(2);

      cout << Form("\n Original distribution: ") << endl; 
      cout << Form(" N = %10i  Integral = %7.2f  Mean = %7.2f  Rms = %f",Np,integral,xMean,xRms) << endl;
  
      // Info for the plots
      TPaveText *textInfo2 = new TPaveText(0.55,0.60,0.82,0.85,"NDC");
      PGlobals::SetPaveTextStyle(textInfo2,32); 
      textInfo2->SetTextColor(kGray+2);
      textInfo2->SetTextFont(42);
      sprintf(ctext,"Entries = %5i",entries);
      textInfo2->AddText(ctext);
      sprintf(ctext,"Integral = %5.2f",integral);
      textInfo2->AddText(ctext);
      sprintf(ctext,"#LTx#GT = %5.2f k_{p}^{-1}",xMean);
      textInfo2->AddText(ctext);
      sprintf(ctext,"#LTx#GT_{rms} = %5.2f k_{p}^{-1}",xRms);
      textInfo2->AddText(ctext);
      sprintf(ctext,"#LTp_{x}#GT = %5.2f mc",pxMean);
      textInfo2->AddText(ctext);
      sprintf(ctext,"#LTp_{x}#GT_{rms} = %5.2f mc",pxRms);
      textInfo2->AddText(ctext);


      // Output histogram
      hP1X1out->GetXaxis()->CenterTitle();
      hP1X1out->GetYaxis()->CenterTitle();
      hP1X1out->GetZaxis()->CenterTitle();

      hP1X1out->GetXaxis()->SetTitle("#zeta [k_{p}^{-1}]");
      hP1X1out->GetYaxis()->SetTitle("p_{z} [mc]");
      hP1X1out->GetZaxis()->SetTitle("charge");
  
      entries = hP1X1out->GetEntries();
      integral = hP1X1out->Integral();
      zMean = hP1X1out->GetMean(1);
      zRms  = hP1X1out->GetRMS(1);
      pzMean = hP1X1out->GetMean(2);
      pzRms  = hP1X1out->GetRMS(2);

      cout << Form("\n New distribution: ") << endl; 
      cout << Form(" N = %10i  Integral = %7.2f  Mean = %7.2f  Rms = %f",Npout,integral,zMean,zRms) << endl;
      cout << Form("\n Total charge = %10.4f pC ", QTotalOut) << endl;
  
  
      TPaveText *textInfoOut = new TPaveText(0.55,0.60,0.82,0.85,"NDC");
      PGlobals::SetPaveTextStyle(textInfoOut,32); 
      textInfoOut->SetTextColor(kGray+2);
      textInfoOut->SetTextFont(42);
      sprintf(ctext,"Entries = %5i",entries);
      textInfoOut->AddText(ctext);
      sprintf(ctext,"Integral = %5.2f",integral);
      textInfoOut->AddText(ctext);
      sprintf(ctext,"#LT#zeta#GT = %5.2f k_{p}^{-1}",zMean);
      textInfoOut->AddText(ctext);
      sprintf(ctext,"#LT#zeta#GT_{rms} = %5.2f k_{p}^{-1}",zRms);
      textInfoOut->AddText(ctext);
      sprintf(ctext,"#LTp_{z}#GT = %5.2f mc",pzMean);
      textInfoOut->AddText(ctext);
      sprintf(ctext,"#LTp_{z}#GT_{rms} = %5.2f mc",pzRms);
      textInfoOut->AddText(ctext);


      // Transverse phasespace histogram
      hP2X2out->GetXaxis()->CenterTitle();
      hP2X2out->GetYaxis()->CenterTitle();
      hP2X2out->GetZaxis()->CenterTitle();

      hP2X2out->GetXaxis()->SetTitle("x [k_{p}^{-1}]");
      hP2X2out->GetYaxis()->SetTitle("p_{x} [mc]");
      hP2X2out->GetZaxis()->SetTitle("charge");
  
      entries = hP2X2out->GetEntries();
      integral = hP2X2out->Integral();
      xMean = hP2X2out->GetMean(1);
      xRms  = hP2X2out->GetRMS(1);
      pxMean = hP2X2out->GetMean(2);
      pxRms  = hP2X2out->GetRMS(2);

      cout << Form("\n New distribution: ") << endl; 
      cout << Form(" N = %10i  Integral = %7.2f  Mean = %7.2f  Rms = %f",Npout,integral,xMean,xRms) << endl;
  
      // Info for the plots
      TPaveText *textInfo2Out = new TPaveText(0.55,0.60,0.82,0.85,"NDC");
      PGlobals::SetPaveTextStyle(textInfo2Out,32); 
      textInfo2Out->SetTextColor(kGray+2);
      textInfo2Out->SetTextFont(42);
      sprintf(ctext,"Entries = %5i",entries);
      textInfo2Out->AddText(ctext);
      sprintf(ctext,"Integral = %5.2f",integral);
      textInfo2Out->AddText(ctext);
      sprintf(ctext,"#LTx#GT = %5.2f k_{p}^{-1}",xMean);
      textInfo2Out->AddText(ctext);
      sprintf(ctext,"#LTx#GT_{rms} = %5.2f k_{p}^{-1}",xRms);
      textInfo2Out->AddText(ctext);
      sprintf(ctext,"#LTp_{x}#GT = %5.2f mc",pxMean);
      textInfo2Out->AddText(ctext);
      sprintf(ctext,"#LTp_{x}#GT_{rms} = %5.2f mc",pxRms);
      textInfo2Out->AddText(ctext);


      // Text objects
      TPaveText *textTime = new TPaveText(0.55,0.85,0.82,0.9,"NDC");
      PGlobals::SetPaveTextStyle(textTime,32); 
      if(opt.Contains("units") && pData->GetPlasmaDensity()) 
	sprintf(ctext,"z = %5.1f mm", Time * skindepth / PUnits::mm);
      else
	sprintf(ctext,"t = %5.1f #omega_{p}^{-1}",Time);
      textTime->AddText(ctext);
    

      TPaveText *textLabel = new TPaveText(0.20,0.20,0.40,0.30,"NDC");
      PGlobals::SetPaveTextStyle(textLabel,32); 
      sprintf(ctext,"Original distribution.");
      textLabel->AddText(ctext);

      TPaveText *textLabelOut = new TPaveText(0.20,0.20,0.40,0.30,"NDC");
      PGlobals::SetPaveTextStyle(textLabelOut,32); 
      textLabelOut->SetTextColor(kRed);
      sprintf(ctext,"Out distribution.");
      textLabelOut->AddText(ctext);


      // Actual Plotting!
      // ------------------------------------------------------------
    
      // Set palette:
      PPalette * pPalette = (PPalette*) gROOT->FindObject("electron0");
      if(!pPalette) {
	pPalette = new PPalette("electron0");
	pPalette->SetPalette("electron0");
	//pPalette->SetAlpha(1);
      }

      C->cd();

      C->Divide(2,2);

      C->cd(1);
    
      gPad->SetFrameLineWidth(2);  
    
      if(opt.Contains("logz")) {
	gPad->SetLogz(1);
      } else {
	gPad->SetLogz(0);
      }
        
      hP1X1->Draw("colz");

      gPad->Update();
    
      textTime->Draw();
      textInfo->Draw();
      textLabel->Draw();
    
      gPad->RedrawAxis(); 

      //---------------------------------------

      C->cd(2);
    
      gPad->SetFrameLineWidth(2);  
    
      if(opt.Contains("logz")) {
	gPad->SetLogz(1);
      } else {
	gPad->SetLogz(0);
      }
    
      hP2X2->Draw("colz");

      gPad->Update();
    
      textTime->Draw();
      textInfo2->Draw();
      textLabel->Draw();
    
      gPad->RedrawAxis(); 

      //----------------------------------------

      C->cd(3);

      gPad->SetFrameLineWidth(2);  
    
      if(opt.Contains("logz")) {
	gPad->SetLogz(1);
      } else {
	gPad->SetLogz(0);
      }
    
      hP1X1out->Draw("colz");
    
      gPad->Update();
    
      textTime->Draw();
      textInfoOut->Draw();
      textLabelOut->Draw();
    
      gPad->RedrawAxis(); 


      //----------------------------------------

      C->cd(4);

      gPad->SetFrameLineWidth(2);  
    
      if(opt.Contains("logz")) {
	gPad->SetLogz(1);
      } else {
	gPad->SetLogz(0);
      }
    
      hP2X2out->Draw("colz");
    
      gPad->Update();
    
      textTime->Draw();
      textInfo2Out->Draw();
      textLabelOut->Draw();
    
      gPad->RedrawAxis(); 
    
    
      C->cd();

      //----------------------------------------

      // Print to a file
      cout << endl;
      TString fOutName = fileout.Remove(fileout.Last('.'));
      PGlobals::imgconv(C,fOutName,opt);
      // ---------------------------------------------------------
      
      
    } 
  
    cout << endl;
  }  
}

