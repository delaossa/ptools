#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>

#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TEllipse.h>
#include <TMarker.h>
#include <TPaveText.h>
#include <TPaletteAxis.h>
#include <TExec.h>
#include <TClonesArray.h>


#include "PData.hh"
#include "PGlobals.hh"
#include "PPalette.hh"


UInt_t BitCounter(UInt_t v) {
  UInt_t c;

  for (c = 0; v; c++) {
    v &= v - 1; // clear the least significant bit set
  }

  return c;
  
}

string DecToBin(int number)
{
  if ( number == 0 ) return "0";
  if ( number == 1 ) return "1";

  if ( number % 2 == 0 )
    return DecToBin(number / 2) + "0";
  else
    return DecToBin(number / 2) + "1";
}

int BinToDec(string number)
{
  int result = 0, pow = 1;
  for ( int i = number.length() - 1; i >= 0; --i, pow <<= 1 )
    result += (number[i] - '0') * pow;

  return result;
}


void PlotSnapshot( const TString &sim, Int_t timestep, UInt_t mask = 3, const TString &options="") {
  
#ifdef __CINT__  
  gSystem->Load("libptools.so");
#endif

  string imask = DecToBin(mask);
  cout << Form("\n Plotting Snapshot with mask: %s",imask.c_str()) << endl; 
  
  PGlobals::Initialize();
  
  // Palettes!
  gROOT->Macro("PPalettes.C");

  // Init Units table
  PUnits::UnitsTable::Get();
  
  // Load PData
  PData *pData = PData::Get(sim.Data());
  pData->LoadFileNames(timestep);
  if(!pData->IsInit()) return;

  TString opt = options;
 
  // More makeup
  gStyle->SetPadGridY(0);
  if(opt.Contains("gridx")) {
    gStyle->SetPadGridX(1);
  }
  if(opt.Contains("gridy")) {
    gStyle->SetPadGridY(1);
  }
  
  // Some plasma constants
  Float_t n0 = pData->GetPlasmaDensity();
  Float_t omegap = pData->GetPlasmaFrequency();
  Float_t timedepth = pData->GetPlasmaTimeDepth();
  Float_t kp = pData->GetPlasmaK();
  Float_t skindepth = pData->GetPlasmaSkinDepth();
  Float_t E0 = pData->GetPlasmaE0();

  // Some beam properties:
  Float_t Ebeam = pData->GetBeamEnergy();
  Float_t gamma = pData->GetBeamGamma();
  Float_t vbeam = pData->GetBeamVelocity();

  // Other parameters
  Float_t trapPotential = 1.0 - (1.0/gamma);

  // Time in OU
  Float_t Time = pData->GetRealTime();
  // z start of the plasma in normalized units.
  Float_t zStartPlasma = pData->GetPlasmaStart()*kp;
  // z start of the beam in normalized units.
  Float_t zStartBeam = pData->GetBeamStart()*kp;
  // z start of the neutral in normalized units.
  Float_t zStartNeutral = pData->GetNeutralStart()*kp;
  // z end of the neutral in normalized units.
  Float_t zEndNeutral = pData->GetNeutralEnd()*kp;
  
  
  // ----------------------------------------------------------------------------------
  
  // Open snapshot file and get the histograms
  TString filename;
  filename = Form("./%s/Plots/Snapshots/Snapshot-%s_%i.root",sim.Data(),sim.Data(),timestep);
  
  TFile  *ifile = (TFile*) gROOT->GetListOfFiles()->FindObject(filename.Data());
  if (!ifile) ifile = new TFile(filename,"READ");

  ifile->cd();

  // Get options
  TTree *infotree = (TTree*) ifile->Get("infotree");
  TString *opttree = 0;
  infotree->SetBranchAddress("options",&opttree);

  Float_t xon = 0.1;
  infotree->SetBranchAddress("xon",&xon);
  Float_t xoff = 0.1;
  infotree->SetBranchAddress("xoff",&xoff);
  
  infotree->GetEntry(0);

  cout << Form(" Options from doSnapshot = %s \n", opttree->Data());

  opt += *opttree;

  cout << Form(" All options = %s \n", opt.Data());

  // Off-axis
  if(opt.Contains("units")) {
    xoff *= skindepth/PUnits::um;
  }  


  if(opt.Contains("center")) {
    Time -= zStartPlasma;
    if(opt.Contains("comov"))      // Centers on the head of the beam.
      Time += zStartBeam;
  }
  Float_t shiftz = pData->Shift(opt);
  
  // Skip one of the species
  Int_t noIndex = -1;
    
  
  
  char hName[24];

  // Get charge density histos
  Int_t Nspecies = pData->NSpecies();
  TH2F **hDen2D = new TH2F*[Nspecies];
  // Get charge density on-axis
  TH1F **hDen1D = new TH1F*[Nspecies];
  // And electric current (integrated)
  TH1F **hCur1D = new TH1F*[Nspecies];
  for(Int_t i=0;i<Nspecies;i++) {
    hCur1D[i] = NULL;
    if(i==noIndex) continue;

    sprintf(hName,"hDen2D_%i",i); 
    hDen2D[i] = (TH2F*) ifile->Get(hName);

    sprintf(hName,"hDen1D_%i",i); 
    hDen1D[i] = (TH1F*) ifile->Get(hName);

    sprintf(hName,"hCur1D_%i",i); 
    hCur1D[i] = (TH1F*) ifile->Get(hName);
    
  }
  
  
  // Get electric fields 2D
  const Int_t Nfields = 3;
  TH2F **hE2D = new TH2F*[Nfields];
  TH2F **hB2D = new TH2F*[Nfields];
  TH1F **hE1D = new TH1F*[Nfields];
  TH1F **hB1D = new TH1F*[Nfields];
  for(Int_t i=0;i<Nfields;i++) {
    hE2D[i] = hB2D[i] = NULL;
    hE1D[i] = hB1D[i] = NULL;

    sprintf(hName,"hE2D_%i",i); 
    hE2D[i] = (TH2F*) ifile->Get(hName);

    sprintf(hName,"hB2D_%i",i); 
    hB2D[i] = (TH2F*) ifile->Get(hName);

    sprintf(hName,"hE1D_%i",i); 
    hE1D[i] = (TH1F*) ifile->Get(hName);

    sprintf(hName,"hB1D_%i",i); 
    hB1D[i] = (TH1F*) ifile->Get(hName);
    
  }

  sprintf(hName,"hETotal2D");   
  TH2F *hETotal2D = (TH2F*) ifile->Get(hName);
  sprintf(hName,"hETotal1D");   
  TH1F *hETotal1D = (TH1F*) ifile->Get(hName);

  sprintf(hName,"hFocus2D"); 
  TH2F *hFocus2D = (TH2F*) ifile->Get(hName);
  sprintf(hName,"hFocus1D"); 
  TH1F *hFocus1D = (TH1F*) ifile->Get(hName);

  sprintf(hName,"hV2D"); 
  TH2F *hV2D = (TH2F*) ifile->Get(hName);
  sprintf(hName,"hV1D"); 
  TH1F *hV1D = (TH1F*) ifile->Get(hName);
  
  
  // Ionization probability rates (ADK)
  // Calculates from the total E the ionization prob. rate for a given species.
  const Int_t NAtoms = 7;
  char atNames[NAtoms][8] = {"H","He","He2","Ne","Ne2","Ne5","HIT"};
  char atAxNames[NAtoms][8] = {"H","He","He^{+}","Ne","Ne^{+}","Ne^{4+}","HIT"};
  Float_t Eion0[NAtoms] = {13.6 * PUnits::eV, 24.59 * PUnits::eV, 54.4 * PUnits::eV, 21.56 * PUnits::eV, 40.96 * PUnits::eV, 126.247 * PUnits::eV, 85.00 *  PUnits::eV};
  Float_t Z[NAtoms]     = { 1.0,  1.0,  2.0,  1.0, 2.0, 5.0, 2.0};
  // Float_t z10[NAtoms]   = {999.0, 999.0, 999.0};
  // Float_t z100[NAtoms]  = {999.0, 999.0, 999.0};

  TH2F *hIonRate2D[NAtoms] = {NULL};
  TH2F *hIonProb2D[NAtoms] = {NULL};
  TH1F *hIonRate1D[NAtoms] = {NULL};
  TH1F *hIonProb1D[NAtoms] = {NULL};

  // Ion probability: Select species index.
  UInt_t ii = 1; // Helium
  // UInt_t ii = 2; // Helium2
  // UInt_t ii = 5; // Ne4
  // UInt_t ii = 6; // Custom
  
  for(Int_t iat=0;iat<NAtoms;iat++) {

    if(iat!=ii) continue;
    
    sprintf(hName,"hIonRate2D_%s",atNames[iat]);   
    hIonRate2D[iat] = (TH2F*) hETotal2D->Clone(hName);
    hIonRate2D[iat]->Reset();
    sprintf(hName,"hIonProb2D_%s",atNames[iat]);
    hIonProb2D[iat] = (TH2F*) hETotal2D->Clone(hName);
    hIonProb2D[iat]->Reset();
    
    Int_t NbinsX = hETotal2D->GetNbinsX();
    Int_t NbinsY = hETotal2D->GetNbinsY();
    Float_t dx = hETotal2D->GetXaxis()->GetBinWidth(1);

    if(opt.Contains("units")) dx *= PUnits::um / PConst::c_light;
    else dx *= timedepth;
      
    // if(!opt.Contains("units")) dx *=  timedepth / PUnits::femtosecond; // to fs
    // else dx /= PConst::c_light / (PUnits::um/PUnits::femtosecond); // to fs
    
    for(Int_t j=1;j<=NbinsY;j++) {
      
      Float_t integral[NAtoms] = {0.0};
      
      for(Int_t k=NbinsX;k>0;k--) {

	Float_t E = hETotal2D->GetBinContent(k,j);
	if(opt.Contains("units")) E *= PUnits::GV/PUnits::m;
	else E *= E0;
	if(E<10) continue;
		
	Float_t IonRate = PFunc::ADK_ENG(E,Eion0[iat],Z[iat]);	
	hIonRate2D[iat]->SetBinContent(k,j,IonRate);
	integral[iat] += IonRate * dx;	
	if(integral[iat]>=1) integral[iat] = 1;

	if(opt.Contains("units"))
	  IonRate /= 1.0/PUnits::femtosecond;
	else
	  IonRate /= omegap;
	
	hIonProb2D[iat]->SetBinContent(k,j,(1-integral[iat])*IonRate);
	
      }
    }

    {
      sprintf(hName,"hIonRate1D_%s",atNames[iat]);
      hIonRate1D[iat] = (TH1F*) hETotal1D->Clone(hName);
      hIonRate1D[iat]->Reset();
      sprintf(hName,"hIonProb1D_%s",atNames[iat]);
      hIonProb1D[iat] = (TH1F*) hETotal1D->Clone(hName);
      hIonProb1D[iat]->Reset();

      Int_t NbinsX = hETotal1D->GetNbinsX();
      Float_t dx = hETotal1D->GetBinWidth(1);
      
      if(opt.Contains("units")) dx *= PUnits::um / PConst::c_light;
      else dx *= timedepth;

      Float_t integral[NAtoms] = {0.0};
      for(Int_t j=NbinsX;j>0;j--) {
	
	Float_t E = hETotal1D->GetBinContent(j);
	if(opt.Contains("units")) E *= PUnits::GV/PUnits::m;
	else E *= E0;
	if(E<10) continue;

	Float_t IonRate = PFunc::ADK_ENG(E,Eion0[iat],Z[iat]);
	hIonRate1D[iat]->SetBinContent(j,IonRate);
	integral[iat] += IonRate * dx;	
	if(integral[iat]>=1) integral[iat] = 1;
	
	if(opt.Contains("units"))
	  IonRate /= 1.0/PUnits::femtosecond;
	else
	  IonRate /= omegap;

	hIonProb1D[iat]->SetBinContent(j,(1-integral[iat])*IonRate);
	
      }
      
    }

    // axis labels
    if(opt.Contains("units")) {
      char axName[24];
      sprintf(axName,"W_{%s} [fs^{-1}]",atAxNames[iat]);
      hIonRate2D[iat]->GetZaxis()->SetTitle(axName);
      sprintf(axName,"#Gamma_{%s} [fs^{-1}]",atAxNames[iat]);
      hIonProb2D[iat]->GetZaxis()->SetTitle(axName);

      sprintf(axName,"W_{%s} [fs^{-1}]",atAxNames[iat]);
      hIonRate1D[iat]->GetZaxis()->SetTitle(axName);
      sprintf(axName,"#Gamma_{%s} [fs^{-1}]",atAxNames[iat]);
      hIonProb1D[iat]->GetZaxis()->SetTitle(axName);
      
    } else {
      char axName[24];
      sprintf(axName,"W_{%s} #omega_{p}",atAxNames[iat]);
      hIonRate2D[iat]->GetZaxis()->SetTitle(axName);
      sprintf(axName,"#Gamma_{%s} #omega_{p}",atAxNames[iat]);
      hIonProb2D[iat]->GetZaxis()->SetTitle(axName);

      sprintf(axName,"W_{%s} #omega_{p}",atAxNames[iat]);
      hIonRate1D[iat]->GetZaxis()->SetTitle(axName);
      sprintf(axName,"#Gamma_{%s} #omega_{p}",atAxNames[iat]);
      hIonProb1D[iat]->GetZaxis()->SetTitle(axName);
    }
    
    
    
  }
  // End of ionization histograms.

  

  // ----------- Plotting Range (ZOOM) ------------

  Float_t zoom = 1;
  Float_t yRange   = (hDen2D[0]->GetYaxis()->GetXmax() - hDen2D[0]->GetYaxis()->GetXmin())/zoom;
  Float_t midPoint = (hDen2D[0]->GetYaxis()->GetXmax() + hDen2D[0]->GetYaxis()->GetXmin())/2.;
  Float_t yMin = midPoint-yRange/2;
  Float_t yMax = midPoint+yRange/2;
  if(pData->IsCyl()) {
    yMin = pData->GetXMin(1);
    yMax = yRange;
  }
  Float_t xMin = hDen2D[0]->GetXaxis()->GetXmin();
  Float_t xMax = hDen2D[0]->GetXaxis()->GetXmax();
  Float_t xRange = xMax - xMin;
  // User range
  // xMin = -70;
  // xMax = -20;

  for(Int_t i=0;i<Nspecies;i++) {
    if(i==noIndex) continue;
    if(!hDen2D[i]) continue;
    hDen2D[i]->GetYaxis()->SetRangeUser(yMin,yMax);
    hDen2D[i]->GetXaxis()->SetRangeUser(xMin,xMax);
  }
  
  for(Int_t i=0;i<Nfields;i++) {
    if(!hE2D[i]) continue;
    hE2D[i]->GetYaxis()->SetRangeUser(yMin,yMax);
    hE2D[i]->GetXaxis()->SetRangeUser(xMin,xMax);
  }

  for(Int_t i=0;i<Nfields;i++) {
    if(!hB2D[i]) continue;
    hB2D[i]->GetYaxis()->SetRangeUser(yMin,yMax);
    hB2D[i]->GetXaxis()->SetRangeUser(xMin,xMax);
  }
  
  hETotal2D->GetYaxis()->SetRangeUser(yMin,yMax);
  hETotal2D->GetXaxis()->SetRangeUser(xMin,xMax);
      
  // ------------- z Zoom --------------- Plasma palette -----------
  // Set the range of the plasma charge density histogram for maximum constrast 
  // using a dynamic palette that adjusts the base value to a certain color.
  
  //  Float_t density = 1; 
  Float_t density = hDen1D[0]->GetBinContent(hDen1D[0]->GetNbinsX());
  Float_t Base  = 1*density;
  
  Float_t *Max = new Float_t[Nspecies];
  Float_t *Min = new Float_t[Nspecies];
  
  for(Int_t i=0;i<Nspecies;i++) {
    if(i==noIndex) continue;
    if(!hDen2D[i]) continue;
    
    Max[i] = hDen2D[i]->GetMaximum();
    // if(i==0) cout << Form("Max. plasma density = %f",Max[i]) << endl;
    Min[i] = 1.01E-1 * Base; 
    
    if(i==0) {
      if(Max[i]<Base) { 
	Max[i] = 1.1*Base;
      } else if(Max[i]<2*Base) {
	Min[i] = 2*Base - Max[i];
      } else {
	Max[i] = 0.4*hDen2D[i]->GetMaximum(); // enhance plasma contrast.
	//Max[i] = hDen2D[i]->GetMaximum(); 
	if(Max[i]<Base) Min[i] = 1.01 * Base;
      }
    } 
    
    if(i==1) {
      //  Max[i] = 1.0 * Max[i];
      Min[i] = 1.01E-1 * Base;
      //Min[i] = 5.01E-2 * Base;
      //Min[i] = 1.01E-2 * Max[i];
      // if(Base<Max[i])
      //  	Min[i] = 1.01E-1 * Base;
      // else
      //  	Min[i] = 1.01E-1 * Max[i];
    }

    if(i==2) {
      Min[i] = 1.01E-4 * Base;
      //Min[i] = 5.01E-3 * Base;
      //Min[i] = 5.01E-2 * Base;
    }
    
    hDen2D[i]->GetZaxis()->SetRangeUser(Min[i],Max[i]);  
  }
  
  
  // Dynamic plasma palette
  const Int_t plasmaDNRGBs = 3;
  const Int_t plasmaDNCont = 64;
  Float_t basePos = 0.5;
  if(Max[0]!=Min[0]) {
    if(opt.Contains("logz")) {
      Float_t a = 1.0/(TMath::Log10(Max[0])-TMath::Log10(Min[0]));
      Float_t b = TMath::Log10(Min[0]);
      basePos = a*(TMath::Log10(Base) - b);
      
    } else {
      basePos = (1.0/(Max[0]-Min[0]))*(Base - Min[0]);
    }
  }

  Double_t plasmaDStops[plasmaDNRGBs] = { 0.00, basePos, 1.00 };
  Double_t plasmaDRed[plasmaDNRGBs]   = { 0.99, 0.90, 0.00 };
  Double_t plasmaDGreen[plasmaDNRGBs] = { 0.99, 0.90, 0.00 };
  Double_t plasmaDBlue[plasmaDNRGBs]  = { 0.99, 0.90, 0.00 };
   
  PPalette * plasmaPalette = (PPalette*) gROOT->FindObject("plasma");
  plasmaPalette->CreateGradientColorTable(plasmaDNRGBs, plasmaDStops, 
					  plasmaDRed, plasmaDGreen, plasmaDBlue, plasmaDNCont,1.0);

  // Change the range of z axis for the fields to be symmetric.
  Float_t *Emax = new Float_t[Nfields];
  Float_t *Emin = new Float_t[Nfields];
  for(Int_t i=0;i<Nfields;i++) {
    if(!hE2D[i]) continue;
    
    Emax[i] = hE2D[i]->GetMaximum();
    Emin[i] = hE2D[i]->GetMinimum();
    if(Emax[i] > TMath::Abs(Emin[i]))
      Emin[i] = -Emax[i];
    else
      Emax[i] = -Emin[i];
    //    hE2D[i]->GetZaxis()->SetRangeUser(Emin[i],Emax[i]); 
    hE2D[i]->GetZaxis()->SetRangeUser(Emin[0],Emax[0]); 
  }

  Float_t ETmin = 0.001;  
  Float_t ETmax = hETotal2D->GetMaximum();
  hETotal2D->GetZaxis()->SetRangeUser(ETmin,ETmax);

  Float_t Fmin = hFocus2D->GetMinimum();
  Float_t Fmax = hFocus2D->GetMaximum();
  if(Fmax > TMath::Abs(Fmin))
    Fmin = -Fmax;
  else
    Fmax = -Fmin;
  hFocus2D->GetZaxis()->SetRangeUser(Fmin,Fmax);

  // Find the first point on-axis where Ez changes from positive to negative:
  Int_t MAXCROSS = 2;
  Float_t *EzCross = new Float_t[MAXCROSS];
  Float_t *EzExtr = new Float_t[MAXCROSS];
  memset(EzCross,0,sizeof(Float_t)*MAXCROSS);
  memset(EzExtr,0,sizeof(Float_t)*MAXCROSS);
  
  Int_t auxNcross = PGlobals::HCrossings(hE1D[0],EzCross,EzExtr,MAXCROSS,0.,0.);
  
  // for(Int_t i=0;i<auxNcross;i++) {
  //   if(opt.Contains("units"))
  //     cout << Form(" %i Ez crossing found at \\zeta = %.2f -> Ez maximum = %.2f GV/m",i,EzCross[i],EzExtr[i]) << endl;
  //   else
  //     cout << Form(" %i Ez crossing found at \\zeta = %.2f -> Ez maximum = %.2f E0",i,EzCross[i],EzExtr[i]) << endl;
      
  // }

  TH1F *hEzclone = (TH1F*) hE1D[0]->Clone("hEzclone");
  hEzclone->GetXaxis()->SetRangeUser(EzCross[0],hEzclone->GetXaxis()->GetXmax());
  Int_t binEbMax = hEzclone->GetMaximumBin();
  Float_t zEbMax = hEzclone->GetBinCenter(binEbMax);
  
  // Potential 
  if(opt.Contains("units") && n0) {
    trapPotential *=  ( E0 * skindepth / (PUnits::MV) ); 
  }  

  // Find the first point on-axis where focusing goes defocusing.
  MAXCROSS = 1;
  Float_t *FocusCross = new Float_t[MAXCROSS];
  Float_t *FocusExtr = new Float_t[MAXCROSS];
  memset(FocusCross,0,sizeof(Float_t)*MAXCROSS);
  memset(FocusExtr,0,sizeof(Float_t)*MAXCROSS);
  auxNcross = PGlobals::HCrossings(hFocus1D,FocusCross,FocusExtr,MAXCROSS,0.0,EzCross[0]);


  Int_t binVmin = hV1D->GetMinimumBin();  
  Float_t Vbase = hV1D->GetMinimum();
  
  // End of focusing
  Int_t binFcross = hFocus1D->FindBin(FocusCross[0]);
  Float_t potValueFin = hV1D->GetBinContent(binFcross);
  Float_t EmaxIon = hE1D[0]->GetBinContent(binFcross);

  if(opt.Contains("vfoc")) {
    Vbase = potValueFin;
    binVmin = binFcross;
  }
  
  Int_t binPotValueIni = -99;  
  Float_t potValueIni = -99;
  for(Int_t j=binVmin;j<hV1D->GetNbinsX();j++) {
    if(hV1D->GetBinContent(j)>=Vbase+trapPotential) {
      binPotValueIni = j;
      potValueIni = hV1D->GetBinContent(j);
      break;
    }
  }
  
  { // Shift potential
    Int_t   NbinsX = hV2D->GetNbinsX(); 
    Int_t   NbinsY = hV2D->GetNbinsY();
    //   Float_t shift = -Vbase - trapPotential;
    Float_t shift = -Vbase;
    for(Int_t j=0;j<=NbinsX;j++) {
      for(Int_t k=0;k<=NbinsY;k++) {
  	hV2D->SetBinContent(j,k, hV2D->GetBinContent(j,k) + shift);
      }
      hV1D->SetBinContent(j, hV1D->GetBinContent(j) + shift);
    }
  }

  Float_t Vzero = hV1D->GetBinContent(binPotValueIni);
  // cout << Form(" VZERO = %f",Vzero) << endl;

  Float_t Vmin = hV1D->GetMinimum();   
  Float_t Vmax = hV1D->GetMaximum();    
  if(Vmax<0.1) Vmax = 0.1;

  // Dynamic potential palette
  const Int_t potPNRGBs = 6;
  const Int_t potPNCont = 64;
  Float_t zeroPos = -(Vmin-Vzero)/(Vmax-Vmin);

  Double_t potPStops[potPNRGBs] = { 0.00, zeroPos-6.0/potPNCont,zeroPos-1.0/potPNCont, zeroPos, zeroPos+3.0/potPNCont, 1.00 };
  Double_t potPRed[potPNRGBs]   = { 0.518, 0.965, 0.90, 0.90, 0.498, 0.106 };
  Double_t potPGreen[potPNRGBs] = { 0.078, 0.925, 0.90, 0.90, 0.718, 0.078 };
  Double_t potPBlue[potPNRGBs]  = { 0.106, 0.353, 0.90, 0.90, 0.780, 0.518 };
   
  PPalette * potentialPalette = (PPalette*) gROOT->FindObject("rbowinv");
  potentialPalette->CreateGradientColorTable(potPNRGBs, potPStops, 
					     potPRed, potPGreen, potPBlue, potPNCont);
  
  // Extract contours from 2D histos
  TCanvas* c = new TCanvas("contours","Contour List",0,0,600,600);
  c->cd();
  
  // Potential 
  TH2F *hV2Dc = (TH2F*) hV2D->Clone("hV2Dc");
  const Int_t Ncontours = 25;
  Double_t contours[Ncontours];
  for(Int_t i=0; i<Ncontours; i++) {
    contours[i] = i*(trapPotential/5.0);// - trapPotential; 
  }
  hV2Dc->SetContour(Ncontours, contours);
  hV2Dc->Draw("cont list 0");
  
  c->Update();
  TObjArray *contsV2D = (TObjArray*) gROOT->GetListOfSpecials()->FindObject("contours");
  TClonesArray graphsV2D("TGraph",Ncontours);
  TClonesArray graphsV2D_main("TGraph",Ncontours);
  {
    Int_t ncontours = contsV2D->GetSize();
    TList* clist = NULL;
    Int_t nGraphs = 0;
    Int_t nGraphs_main = 0;
    TGraph *gr = NULL;
    for(Int_t i = 0; i < ncontours; i++){  
      // negative contours are returned first from low to high absolute values. 
      // e.g. -1, -2, -3, ..., 0, 1, 2, 3,...
            
      clist = (TList*) contsV2D->At(i);
    
      for(Int_t j = 0 ; j < clist->GetSize(); j++) {
	//	cout << Form(" psi contour (%i,%i)",i,j) << endl;

	gr = (TGraph*) clist->At(j);
	if(!gr) continue;
      
	gr->SetLineWidth(1);
	gr->SetLineColor(kGray+1);
      
	if( i==0 || i==5 || i==10 ) {
	  gr->SetLineWidth(2);
	  gr->SetLineColor(kGray+2);
	} 
	
	new(graphsV2D[nGraphs]) TGraph(*gr) ;
	nGraphs++;
	
	if(i==0 && j==0) {
	  TGraph *grm = new(graphsV2D_main[nGraphs_main]) TGraph(*gr);
	  nGraphs_main++;
	  grm->SetLineWidth(1);
	  grm->SetLineStyle(3);
	  grm->SetLineColor(kGray+3);	  
	}

	if(i==5 && j==0) {
	  TGraph *grm = new(graphsV2D_main[nGraphs_main]) TGraph(*gr);
	  nGraphs_main++;
	  grm->SetLineWidth(2);
	  grm->SetLineStyle(3);
	  grm->SetLineColor(kGray+3);	  
	}

      }
    }
  }
  
  Float_t IPmin = 0.00;
  Float_t IPmax = hIonProb2D[ii]->GetMaximum();
  
  IPmax = hIonProb1D[ii]->GetMaximum();
  hIonProb2D[ii]->GetZaxis()->SetRangeUser(IPmin,IPmax);
    
  TH2F *hIonProb2Dc = (TH2F*) hIonProb2D[ii]->Clone("hIonProb2Dc");
  
  const Int_t NcontI = 1;
  Double_t contI[NcontI] = {0.3*IPmax};
  // Double_t contI[NcontI] = {0.1*IPmax};
  //const Int_t NcontI = 1;
  //Float_t contI[NcontI] = {0.1};
  hIonProb2Dc->SetContour(NcontI, contI);
  hIonProb2Dc->Draw("cont list 0");
  
  c->Update();
  TObjArray *contsI2D = (TObjArray*) gROOT->GetListOfSpecials()->FindObject("contours");
  TClonesArray graphsI2D("TGraph",NcontI);
  TClonesArray graphsI2D_main("TGraph",NcontI);
  {
    Int_t ncontours = contsI2D->GetSize();
    TList* clist = NULL;
    Int_t nGraphs = 0;
    Int_t nGraphs_main = 0;
    TGraph *gr = NULL;
    for(Int_t i = 0; i < ncontours; i++){
      clist = (TList*) contsI2D->At(i);
      
      for(Int_t j = 0 ; j < clist->GetSize(); j++) {
	gr = (TGraph*) clist->At(j);
	if(!gr) continue;
	
	if( i==0 ) {
	  TGraph *grm = new(graphsI2D_main[nGraphs_main]) TGraph(*gr);
	  nGraphs_main++;
	  grm->SetLineWidth(1);
	  grm->SetLineStyle(2);
	  grm->SetLineColor(kGray+2);
	  
	  TGraph *gr2 = new(graphsI2D[nGraphs]) TGraph(*gr) ;
	  nGraphs++;
	  if(j==0) {
	    gr2->SetLineWidth(1);
	    gr2->SetLineStyle(2);
	    //	    gr2->SetLineColor(PGlobals::elecLine);
	    gr2->SetLineColor(kGray+2);
	  } else {
	    gr2->SetLineWidth(1);
	    gr2->SetLineStyle(1);
	    gr2->SetLineColor(kGray+2);
	  }
	  
	} 
	
      }
    }
  }

  delete c;

  
  // Plotting
  // ------------------------------------------------------------
  
  // Canvas setup
  UInt_t NPad = BitCounter(mask);
  if(NPad==0) NPad = 1;
  Float_t ypadsize = 250;
  Float_t ymarginsize = 200;
  if(NPad==1) ypadsize = 300;
  Float_t ysize = ypadsize * NPad + ymarginsize; 
  Float_t boom = 1.2;
  if(opt.Contains("boom"))
    ysize *= boom;
  TCanvas *C = new TCanvas("C","Snapshot",1050,ysize);
  C->SetFillStyle(4000);

  UInt_t lineColor = kOrange+10;
  //UInt_t lineColor =  TColor::GetColor(196,30,78);
  
  // Setup Pad layout:
  TPad **pad = new TPad*[NPad];
  TH2F **hFrame = new TH2F*[NPad];
  Float_t bMargin = 0.12 * (950/ysize);
  Float_t tMargin = 0.04 * (950/ysize);
  Float_t lMargin = 0.14;
  Float_t rMargin = 0.18;
  Float_t mMargin = 0.02 * (950/ysize);
  Float_t pfactor = 1.6;
  if(opt.Contains("nomar"))
    bMargin = tMargin = lMargin = rMargin = mMargin = 0.0;
  if(NPad==1)
    PGlobals::CanvasPartition(C,NPad,lMargin,rMargin,bMargin,tMargin,mMargin);
  else
    PGlobals::CanvasAsymPartition(C,NPad,lMargin,rMargin,bMargin,tMargin,pfactor,mMargin);
 
  // Define the frames for plotting
  Int_t fonttype = 43;
  Int_t fontsize = 38;
  Int_t tfontsize = 42;
  Int_t txsize = tfontsize+6;
  Int_t lxsize = fontsize+2;
  Int_t tysize = tfontsize;
  Int_t lysize = fontsize;
  Int_t tzsize = tfontsize-4;
  Int_t lzsize = fontsize-6;
  Float_t txoffset = (250/ypadsize) * 2.4 / (950/ysize);
  Float_t lxoffset = 0.015;
  Float_t tyoffset = 1.2 / (950/ysize);
  Float_t lyoffset = 0.01;
  Float_t tzoffset = 1.4 / (950/ysize);
  Float_t lzoffset = 0.01;
  Float_t tylength = 0.015;
  Float_t txlength = 0.04;
  for(Int_t i=NPad-1;i>=0;i--) {
    char name[16];
    sprintf(name,"pad_%i",i);
    pad[i] = (TPad*) gROOT->FindObject(name);
    pad[i]->SetFrameLineWidth(2);
    if(opt.Contains("tick")) {
      pad[i]->SetTickx(1);
      pad[i]->SetTicky(1);
    }
    if(opt.Contains("trans"))
      pad[i]->SetFillStyle(4000);
    pad[i]->SetFrameFillStyle(4000);

    sprintf(name,"hFrame_%i",i);
    hFrame[i] = (TH2F*) gROOT->FindObject(name);
    if(hFrame[i]) delete hFrame[i];
    hFrame[i] = (TH2F*) hDen2D[0]->Clone(name);
    hFrame[i]->Reset();
    
    Float_t xFactor = pad[NPad-1]->GetAbsWNDC()/pad[i]->GetAbsWNDC();
    Float_t yFactor = pad[NPad-1]->GetAbsHNDC()/pad[i]->GetAbsHNDC();

    // Format for y axis
    hFrame[i]->GetYaxis()->SetLabelFont(fonttype);
    hFrame[i]->GetYaxis()->SetLabelSize(lysize);
    hFrame[i]->GetYaxis()->SetLabelOffset(lyoffset);
    hFrame[i]->GetYaxis()->SetTitleFont(fonttype);
    hFrame[i]->GetYaxis()->SetTitleSize(tysize);
    hFrame[i]->GetYaxis()->SetTitleOffset(tyoffset);

    hFrame[i]->GetYaxis()->SetTickLength(xFactor*tylength/yFactor);

    // Format for x axis
    hFrame[i]->GetXaxis()->SetLabelFont(fonttype);
    hFrame[i]->GetXaxis()->SetLabelSize(lxsize);
    hFrame[i]->GetXaxis()->SetLabelOffset(lxoffset);
    hFrame[i]->GetXaxis()->SetTitleFont(fonttype);
    hFrame[i]->GetXaxis()->SetTitleSize(txsize);
    hFrame[i]->GetXaxis()->SetTitleOffset(txoffset);
    
    hFrame[i]->GetXaxis()->SetTickLength(yFactor*txlength/xFactor);      

    if(i>0) { // skip x axis labels except for the lowest one
      hFrame[i]->GetXaxis()->SetLabelSize(0.0);
      hFrame[i]->GetXaxis()->SetTitleSize(0.0);
    }

    if(opt.Contains("nomar")) {
      hFrame[i]->GetYaxis()->SetTickLength(0.0);
      hFrame[i]->GetXaxis()->SetTickLength(0.0);      
    }

    // Labels for the frames

  }

  // Text objects
  TPaveText **textLabel = new TPaveText*[NPad];
  TString pNames[] = {"Den","Ez","Ex","Ey","F","ET","I","V"};
  TString sLabels[] = {"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)"};
    
  char ctext[128];
  if(opt.Contains("units") && n0) 
    sprintf(ctext,"z = %5.1f #mum", Time * skindepth / PUnits::um);
  else
    sprintf(ctext,"t = %5.1f #omega_{p}^{-1}",Time);
  
  TLatex *textTime = new TLatex(xMax - (xMax-xMin)/20.,yMax - (yMax-yMin)/10.,ctext);
  textTime->SetTextAlign(32);
  textTime->SetTextFont(fonttype);
  textTime->SetTextSize(fontsize-10);

  if(opt.Contains("units") && n0) 
    //    sprintf(ctext,"n_{0} = %5.1f x 10^{17} cm^{-3}", 1E-17*n0/(1./PUnits::cm3));
    sprintf(ctext,"n_{0} = %5.2f x 10^{18} cm^{-3}", 1E-18*n0/(1./PUnits::cm3));
   
  TLatex *textDen = new TLatex(xMin + (xMax-xMin)/20.,yMax - (yMax-yMin)/10.,ctext);
  textDen->SetTextAlign(12);
  textDen->SetTextFont(fonttype);
  //  textDen->SetTextColor(lineColor);
  textDen->SetTextColor(kGray+3);
  textDen->SetTextSize(fontsize-14);

  // Lines 
  TLine *lineYzero = new TLine(xMin,0.0,xMax,0.0);
  lineYzero->SetLineColor(kGray+2);
  lineYzero->SetLineStyle(2);

  zStartPlasma -= shiftz; 
  zStartNeutral -= shiftz; 
  zEndNeutral -= shiftz; 
  
  if(opt.Contains("units") && n0) {
    zStartPlasma *= skindepth / PUnits::um;
    zStartNeutral *= skindepth / PUnits::um;
    zEndNeutral *= skindepth / PUnits::um;
  }

  //  cout << "Start plasma = " << zStartPlasma << endl;
  TLine *lineStartPlasma = new TLine(zStartPlasma,yMin,zStartPlasma,yMax);
  lineStartPlasma->SetLineColor(kGray+2);
  lineStartPlasma->SetLineStyle(2);
  lineStartPlasma->SetLineWidth(3);

  //  cout << "Start neutral = " << zStartNeutral << endl;
  TLine *lineStartNeutral = new TLine(zStartNeutral,yMin,zStartNeutral,yMax);
  lineStartNeutral->SetLineColor(kGray+3);
  lineStartNeutral->SetLineStyle(3);
  lineStartNeutral->SetLineWidth(2);

  //  cout << "End neutral = " << zEndNeutral << endl;
  TLine *lineEndNeutral = new TLine(zEndNeutral,yMin,zEndNeutral,yMax);
  lineEndNeutral->SetLineColor(kGray+3);
  lineEndNeutral->SetLineStyle(3);
  lineEndNeutral->SetLineWidth(2);



  // Usable accelerating field (calculated from the real slope):
  Double_t EzMaxUse;
  {
    Float_t zshift = 0.1; 
    if(opt.Contains("units")) 
      zshift *= PUnits::um * kp;
    
    Int_t bin1 = hE1D[0]->FindBin(EzCross[0]-zshift);
    Int_t bin2 = hE1D[0]->FindBin(EzCross[0]+zshift);
    Float_t EzSlope = (hE1D[0]->GetBinContent(bin2)-hE1D[0]->GetBinContent(bin1))/
      (hE1D[0]->GetBinCenter(bin2)-hE1D[0]->GetBinCenter(bin1));

    EzMaxUse  = -EzSlope * (EzCross[1]-EzCross[0]);
  }

  // Get position of the injected beam:
  Float_t zInjBeam = -999;
  Float_t zrmsInjBeam = -999;
  Float_t EInjBeam = -999;
  Float_t ESlopeInjBeam = -999;


  if(Nspecies>2) {
    if(hCur1D[2] && noIndex!=2) {
      hCur1D[2]->ResetStats();
      zInjBeam = hCur1D[2]->GetMean();
      zrmsInjBeam = hCur1D[2]->GetRMS();
      
      // Find field value corresponding to this position
      Int_t ibin0 = hE1D[0]->FindBin(zInjBeam);
      Int_t ibin1 = hE1D[0]->FindBin(zInjBeam-zrmsInjBeam);
      Int_t ibin2 = hE1D[0]->FindBin(zInjBeam+zrmsInjBeam);
      EInjBeam = hE1D[0]->GetBinContent(ibin0);
      ESlopeInjBeam = (hE1D[0]->GetBinContent(ibin1)-hE1D[0]->GetBinContent(ibin2))/
	(hE1D[0]->GetBinCenter(ibin1)-hE1D[0]->GetBinCenter(ibin2));

      cout << endl << Form(" Witness bunch properties") << endl;
      if(opt.Contains("units")) {
	cout << Form("  Position at z bunch = %6.2f um ", zInjBeam) << endl ;
	cout << Form("  RMS of the bunch    = %6.2f um ", zrmsInjBeam) << endl ;
	cout << Form("  E_z at z bunch      = %6.2f GV/m ", EInjBeam) << endl ;
	cout << Form("  Slope at z bunch    = %6.2f (GV/m)um ", ESlopeInjBeam) << endl ;
      } else {
	cout << Form("  Position at z bunch = %6.2f ", zInjBeam) << endl ;
	cout << Form("  RMS of the bunch    = %6.2f ", zrmsInjBeam) << endl ;
	cout << Form("  E_z at z bunch      = %6.2f ", EInjBeam) << endl ;
	cout << Form("  Slope at z bunch    = %6.2f ", ESlopeInjBeam) << endl ;
      }
      cout << Form("  Transformer ratio   = %6.2f ", EInjBeam/EzExtr[0]) << endl ;
    } 
    
  }


  // Calculate blowout parameters
  Double_t maxcurr;
  Double_t rmslength;
  Double_t radius;
  Double_t EzMaxLo;
  Double_t radiusLu;
  Double_t EzMaxLuLinear;
  Double_t EzMaxBeamLo;
  Double_t DeltaPsiLo;
  Double_t RadiusPsiLo;

  // Get max current in norm units:
  if(Nspecies>1) {
    if(hCur1D[1] && noIndex!=1) {  
 
      hCur1D[1]->ResetStats();
      maxcurr = hCur1D[1]->GetMaximum();
      rmslength = hCur1D[1]->GetRMS();
      if(opt.Contains("units")) {
	maxcurr *= PUnits::kA / PConst::I0;
	rmslength *= PUnits::um * kp;
      }

      radius = PFunc::RadiusBO(maxcurr,rmslength);
      EzMaxLo = 0.5 * radius;
      radiusLu = 2.0 * TMath::Sqrt(maxcurr);
      EzMaxLuLinear = PFunc::EzMaxLu(maxcurr);
      EzMaxBeamLo = PFunc::EzMaxBeamLo(maxcurr);
      DeltaPsiLo = radius*radius/4.0;
      RadiusPsiLo = TMath::Sqrt(2*(DeltaPsiLo-1));

      cout << endl << Form(" Blowout parameters") << endl;
      if(opt.Contains("units")) {
	maxcurr /= PUnits::kA / PConst::I0;
	rmslength /= PUnits::um * kp;
	radius   *= skindepth / PUnits::um;
	radiusLu *= skindepth / PUnits::um;
	EzMaxLo  *= E0 / (PUnits::GV/PUnits::m);
	EzMaxLuLinear  *= E0 / (PUnits::GV/PUnits::m);
	EzMaxBeamLo  *= E0 / (PUnits::GV/PUnits::m);
	DeltaPsiLo   *= E0 * skindepth /(PUnits::MV);
	RadiusPsiLo  *= skindepth / PUnits::um;
	
	cout << Form("  Beam max. curr.  = %7.2f kA   RMS = %.2f um",maxcurr,rmslength) << endl;
	cout << Form("  Ez max beam      = %7.2f GV/m",EzExtr[0]) << endl;
	cout << Form("  Radius (Lotov) R = %7.2f um",radius) << endl;
	cout << Form("  Radius (Lu)    R = %7.2f um",radiusLu) << endl;
	cout << Form("  Ez max beam (Lo) = %7.2f GV/m",EzMaxBeamLo) << endl;
	cout << Form("  Ez max acc. (Lo) = %7.2f GV/m",-EzMaxLo) << endl;
	cout << Form("  Ez max (Lu)      = %7.2f GV/m",-EzMaxLuLinear) << endl;
	cout << Form("  E  max (ion)     = %7.2f GV/m",EmaxIon) << endl;
	cout << Form("  DPsi max. (Lo)   = %7.2f MV",DeltaPsiLo) << endl;
	cout << Form("  Trapping R (Lo)  = %7.2f um",RadiusPsiLo) << endl;
      } else {
	cout << Form("  Beam max. curr.  = %7.2f I0   RMS = %.2f kp^-1",maxcurr,rmslength) << endl;
	cout << Form("  Ez max beam      = %7.2f E0",EzExtr[0]) << endl;
	cout << Form("  Radius (Lotov) R = %7.2f kp^-1",radius) << endl;
	cout << Form("  Radius (Lu)    R = %7.2f kp^-1",radiusLu) << endl;
	cout << Form("  Ez max beam (Lo) = %7.2f E0",EzMaxBeamLo) << endl;
	cout << Form("  Ez max acc. (Lo) = %7.2f E0",-EzMaxLo) << endl;
	cout << Form("  Ez max acc. (Lu) = %7.2f E0",-EzMaxLuLinear) << endl;
	cout << Form("  E  max (ion)     = %7.2f E0",EmaxIon) << endl;
	cout << Form("  DPsi max.   (Lo) = %7.2f \\Psi0",DeltaPsiLo) << endl;
	cout << Form("  Trapping R  (Lo) = %7.2f kp^-1",RadiusPsiLo) << endl;
      }
    }  
  }
  
  cout << endl;

  // Access to color Palettes 
  TExec *exPlasma = new TExec("exPlasma","plasmaPalette->SetAlpha(0.8); plasmaPalette->cd();");
  //  TExec *exPlasma = new TExec("exPlasma","plasmaPalette->cd();");
  //  TExec *exElec   = new TExec("exElec","redelectronPalette->SetAlpha(0.8); redelectronPalette->cd();");
  TExec *exElec   = new TExec("exElec","redelectronPalette->cd();");
  //  TExec *exHot    = new TExec("exHot","hotPalette->SetAlpha(0.8); hotPalette->cd();");
  TExec *exHot    = new TExec("exHot","hotPalette->cd();");
  TExec *exField  = new TExec("exField","rbowwhitePalette->cd();");
  TExec *exFieldT = new TExec("exFieldT","red0Palette->cd();");
  //TExec *exIonP   = new TExec("exIonP","grayPalette->cd();");
  TExec *exIonP   = new TExec("exIonP","redelectron0Palette->cd();");
  TExec *exPot    = new TExec("exPot","rbowinvPalette->cd();");
     
  TString drawopt = "colz same";

  // Actual Plotting!
  // ------------------------------------------------------------

  C->cd(0);

  Float_t x1,x2,y1,y2;
  Float_t gap = 0.005;
  TPaletteAxis *palette = NULL;
  UInt_t ip = NPad-1;

  if(mask == 0) {
        
    pad[ip]->Draw();
    pad[ip]->cd(); 
    if(opt.Contains("logz")) {
      pad[ip]->SetLogz(1);
    } else {
      pad[ip]->SetLogz(0);
    }
    
    hFrame[ip]->Draw("col");

    if(opt.Contains("cont")) {

      // // ADK contours  
      // for(Int_t i=0;i<graphsI2D.GetEntriesFast();i++) {
      // 	TGraph *gr = (TGraph*) graphsI2D_main.At(i);
      // 	if(!gr) continue;
      //   gr->Draw("C");
      // }

      // PSI MAIN contour  
      for(Int_t i=0;i<graphsV2D_main.GetEntriesFast();i++) {
	TGraph *gr = (TGraph*) graphsV2D_main.At(i);
	if(!gr) continue;
	gr->Draw("C");
      }
      
    }
    
    
    ip--;  
    C->cd(0);
  }

  if(mask & 0x1) {
    
    pad[ip]->Draw();
    pad[ip]->cd(); 
    if(opt.Contains("logz")) {
      pad[ip]->SetLogz(1);
    } else {
      pad[ip]->SetLogz(0);
    }
    
    hFrame[ip]->Draw();   

    // cout << Form(" N species = %i",Nspecies) << endl;

    if(Nspecies==1) {
      // Just plasma
      if(hDen2D[0] && noIndex!=0) {
	hDen2D[0]->GetZaxis()->SetNdivisions(503);
	hDen2D[0]->GetZaxis()->SetTitleFont(fonttype);
	hDen2D[0]->GetZaxis()->SetTickLength(tylength);
	exPlasma->Draw();
	//exElec->Draw();
	
	hDen2D[0]->Draw(drawopt);
      }
    } else if(Nspecies==2) {
      // Plasma
      if(hDen2D[0] && noIndex!=0) {
	hDen2D[0]->GetZaxis()->SetNdivisions(503);
	hDen2D[0]->GetZaxis()->SetTitleFont(fonttype);
	exPlasma->Draw();
	hDen2D[0]->Draw(drawopt);
      }

      // Beam driver.
      if(hDen2D[1] && noIndex!=1) {
	hDen2D[1]->GetZaxis()->SetNdivisions(503);
	hDen2D[1]->GetZaxis()->SetTitleFont(fonttype);
	exElec->Draw();
	hDen2D[1]->Draw(drawopt);
      }
    } else if (Nspecies==3) {

      // Injected electrons ?
      if(hDen2D[2] && noIndex!=2) {
	exHot->Draw();
	//exElec->Draw();
	hDen2D[2]->GetZaxis()->SetNdivisions(503);
	hDen2D[2]->GetZaxis()->SetTitleFont(fonttype);
	hDen2D[2]->Draw(drawopt);
      }

      // Plasma
      if(hDen2D[0] && noIndex!=0) {
	hDen2D[0]->GetZaxis()->SetNdivisions(503);
	hDen2D[0]->GetZaxis()->SetTitleFont(fonttype);
	//plasmaPalette->SetAlpha(1.0);  
	exPlasma->Draw();
	hDen2D[0]->Draw("colz same");
      }
      
      // Beam driver.
      if(hDen2D[1] && noIndex!=1) {
	hDen2D[1]->GetZaxis()->SetNdivisions(503);
	hDen2D[1]->GetZaxis()->SetTitleFont(fonttype);
	exElec->Draw();
	//exPlasma->Draw();
	hDen2D[1]->Draw(drawopt);
      }
      
    }
    
 
    if(opt.Contains("contall")) {
      // ADK contours  
      for(Int_t i=0;i<graphsI2D.GetEntriesFast();i++) {
	TGraph *gr = (TGraph*) graphsI2D.At(i);
	if(!gr) continue;
	gr->Draw("C");
      }  

      // PSI contours  
      for(Int_t i=0;i<graphsV2D.GetEntriesFast();i++) {
	TGraph *gr = (TGraph*) graphsV2D.At(i);
	if(!gr) continue;
	gr->Draw("C");
      }
    } else if(opt.Contains("cont")) {
      // ADK MAIN contour  
      for(Int_t i=0;i<graphsI2D_main.GetEntriesFast();i++) {
	TGraph *gr = (TGraph*) graphsI2D_main.At(i);
	if(!gr) continue;
	gr->Draw("C");
      }  

      // PSI MAIN contour  
      for(Int_t i=0;i<graphsV2D_main.GetEntriesFast();i++) {
	TGraph *gr = (TGraph*) graphsV2D_main.At(i);
	if(!gr) continue;
	gr->Draw("C");
      }
    } 


    if(opt.Contains("bopar")) {
      TEllipse *ellipLo = new TEllipse(EzCross[0],0.0,radius);	
      ellipLo->SetFillStyle(0);
      ellipLo->SetLineStyle(3);
      ellipLo->SetLineColor(kGray+3);
      ellipLo->Draw();
	
      TEllipse *ellipLu = new TEllipse(EzCross[0],0.0,radiusLu);
      ellipLu->SetFillStyle(0);
      ellipLu->SetLineStyle(4);
      ellipLu->SetLineColor(kGray+3);
      ellipLu->Draw();
    }

    if(opt.Contains("sline")) {
      if(zStartPlasma>xMin && zStartPlasma<xMax)
	lineStartPlasma->Draw();
      if(zStartNeutral>xMin && zStartNeutral<xMax)
	lineStartNeutral->Draw();
      if(zEndNeutral>xMin && zEndNeutral<xMax)
	lineEndNeutral->Draw();
    }
    
    // Palettes re-arrangement
    pad[ip]->Update();
    y1 = pad[ip]->GetBottomMargin();
    y2 = 1 - pad[ip]->GetTopMargin();
    x1 = pad[ip]->GetLeftMargin();
    x2 = 1 - pad[ip]->GetRightMargin();
    Float_t pvsize = (y2-y1)/float(Nspecies);
    if(noIndex>=0 && Nspecies>2) pvsize = (y2-y1)/float(Nspecies-1);
    
    for(Int_t i=0;i<Nspecies;i++) {
      if(i==noIndex) continue;
      if(!hDen2D[i]) continue;
      palette = (TPaletteAxis*)hDen2D[i]->GetListOfFunctions()->FindObject("palette");
      if(palette) {
	palette->SetY1NDC(y1 + i*pvsize + gap);
	palette->SetY2NDC(y1 + (i+1)*pvsize - gap);
	palette->SetX1NDC(x2 + 0.005);
	palette->SetX2NDC(x2 + 0.03);
	palette->SetTitleOffset(tzoffset);
	//	palette->SetTitleOffset(999.9);
	palette->SetTitleSize(tzsize);
	palette->SetLabelFont(fonttype);
	palette->SetLabelSize(lzsize);
	if(opt.Contains("logz")) { 
	  palette->SetLabelOffset(-lyoffset * (250/ypadsize));
	} else
	  palette->SetLabelOffset(lyoffset);
	palette->SetBorderSize(2);
	palette->SetLineColor(1);
      } 
    }

    // 1D plots
    Float_t yaxismin  =  pad[ip]->GetUymin();
    Float_t yaxismax  =  pad[ip]->GetUymin() + 0.3*(pad[ip]->GetUymax() - pad[ip]->GetUymin()) - 0.00;
    TGaxis **axis = new TGaxis*[Nspecies];
    if(opt.Contains("den1d")) {

      Float_t denmin = Min[1];
      Float_t denmax = Max[1];
      if(opt.Contains("logz")) {
	denmin = TMath::Log10(denmin);
	denmax = TMath::Log10(denmax);
      }
      
      for(Int_t i=0;i<Nspecies;i++) {
	if(i==noIndex) continue;
	if(!hDen1D[i] || i!=1) continue;
    
	Float_t slope = (yaxismax - yaxismin)/(denmax - denmin);
    
	for(Int_t j=0;j<hDen1D[i]->GetNbinsX();j++) {
	  Float_t content = hDen1D[i]->GetBinContent(j+1);
	  if(opt.Contains("logz")) content = TMath::Log10(content); 
	  
	  if(content<denmin) 
	    hDen1D[i]->SetBinContent(j+1,yaxismin);
	  else 
	    hDen1D[i]->SetBinContent(j+1,(content - denmin) * slope + yaxismin);
	}

	hDen1D[i]->SetLineWidth(2);
	hDen1D[i]->SetLineColor(PGlobals::elecLine);
	hDen1D[i]->Draw("same C");
      }
    } else if(opt.Contains("cur1d")) {
      
      for(Int_t i=0;i<Nspecies;i++) {
	if(i==noIndex) continue;
	if(!hCur1D[i] || i==0 ) continue;
	
	Float_t curmin = 0.001;
	Float_t curmax = hCur1D[i]->GetMaximum();
	// Round for better plotting
	//curmax = TMath::CeilNint(10*curmax)/10.0;
	//curmax = TMath::CeilNint(curmax);

	if(curmax ==0) continue;
	if(curmax ==9) curmax = 10;

	Float_t slope = (yaxismax - yaxismin)/(curmax - curmin);
	Float_t zPos = hCur1D[i]->GetMean() + 1.5 * hCur1D[i]->GetRMS();
	if(i==2) zPos = hCur1D[i]->GetMean() + 5.0 * hCur1D[i]->GetRMS();
      
	for(Int_t j=0;j<hCur1D[i]->GetNbinsX();j++) {
	  Float_t content = hCur1D[i]->GetBinContent(j+1);
	  
	  hCur1D[i]->SetBinContent(j+1,(content - curmin) * slope + yaxismin);	
	}
	
	hCur1D[i]->SetLineWidth(2);
	if(i==1)
	  hCur1D[i]->SetLineColor(PGlobals::elecLine);
	else
	  hCur1D[i]->SetLineColor(kOrange+8);
	  
	hCur1D[i]->Draw("same C");

	if(!opt.Contains("noaxis")) {
	
	  // Current axis
	  axis[i] = new TGaxis(zPos,yMin,
			       zPos,yaxismax,
			       curmin,curmax,505,"+LS");
	  // axis[i]->SetNdivisions(0);
	  
	  axis[i]->SetLineWidth(1);
	  axis[i]->SetLineColor(kGray+3);//PGlobals::elecLine);
	  axis[i]->SetLabelColor(kGray+3);//PGlobals::elecLine);
	  axis[i]->SetLabelFont(fonttype);
	  axis[i]->SetLabelSize(lzsize-10);
	  axis[i]->SetLabelOffset(lyoffset);
	  axis[i]->SetTitleColor(kGray+3);//PGlobals::elecLine);
	  axis[i]->SetTitleFont(fonttype);
	  axis[i]->SetTitleSize(tzsize-14);
	  axis[i]->SetTitleOffset(tyoffset);
	  axis[i]->SetTickSize(0.03);
	  axis[i]->SetTitle(hCur1D[i]->GetYaxis()->GetTitle());

	  axis[i]->CenterTitle();
	  axis[i]->SetNdivisions(503);
    
	  axis[i]->Draw();
	}
      }
    } else if(opt.Contains("ez1d")) {
      // Fit the E1D in the pad:
      Float_t y1 = yMin + yRange/80.;
      Float_t y2 = y1 + yRange/2.8 - 2*yRange/80.;
      Float_t slope = (y2-y1)/(Emax[0]-Emin[0]);
      
      TLine *lineEzero = new TLine(xMin,(0.0-Emin[0])*slope + y1,xMax,(0.0-Emin[0])*slope + y1);
      lineEzero->SetLineColor(lineColor);
      lineEzero->SetLineStyle(2);
      lineEzero->SetLineWidth(1);
      lineEzero->Draw();

      TH1F *hE1Dclone = (TH1F*) hE1D[0]->Clone("hE1Dclone");

      for(Int_t j=0;j<hE1Dclone->GetNbinsX();j++) {
	hE1Dclone->SetBinContent(j+1,(hE1Dclone->GetBinContent(j+1)-Emin[0])*slope + y1);
      }
      hE1Dclone->SetLineStyle(1);
      // hE1Dclone->SetLineWidth(3);
      hE1Dclone->SetLineWidth(2);
        
      hE1Dclone->SetLineColor(lineColor);

      hE1Dclone->Draw("sameL");
    } else if (opt.Contains("et1d")) {
      // Fit the E1D in the E2D pad:
      Float_t HeIon  = 92.75;  // GV/m
      if(!opt.Contains("units") || !n0 ) {
	HeIon  /= ( E0 / (PUnits::GV/PUnits::m));
      }

      // Fit the E1D in the pad:
      Float_t ETmin =  0.00;//hETotal1D->GetMinimum();
      Float_t ETmax =  hETotal1D->GetMaximum();
      Float_t y1 = yMin + yRange/80.;
      Float_t y2 = y1 + yRange/2.8 - 2*yRange/80.;
      Float_t slope = (y2-y1)/(ETmax-ETmin);
  
      TLine *lineEHe = new TLine(xMin,(HeIon-ETmin)*slope + y1,xMax,(HeIon-ETmin)*slope + y1);
      lineEHe->SetLineColor(lineColor);
      lineEHe->SetLineStyle(2);
      lineEHe->SetLineWidth(1);
      lineEHe->Draw();

      TLine *lineEzero = new TLine(xMin,(0.0-ETmin)*slope + y1,xMax,(0.0-ETmin)*slope + y1);
      lineEzero->SetLineColor(kGray+1);
      lineEzero->SetLineStyle(2);
      lineEzero->SetLineWidth(1);
      lineEzero->Draw();
       
      for(Int_t j=0;j<hETotal1D->GetNbinsX();j++) {
	hETotal1D->SetBinContent(j+1,(hETotal1D->GetBinContent(j+1)-ETmin)*slope + y1);
      }
      hETotal1D->SetLineStyle(1);
      hETotal1D->SetLineWidth(2);
    
      hETotal1D->SetLineColor(lineColor);
    
      hETotal1D->Draw("sameC");
    } 

    if (opt.Contains("curInj")) {

      if(hCur1D[2]) {

	// Fit the 
	Float_t curInjmin =  0.00;
	Float_t curInjmax =  hCur1D[2]->GetMaximum();
	Float_t y1 = yMin + yRange/80.;
	Float_t y2 = y1 + yRange/2.8 - 2*yRange/80.;
	Float_t slope = (y2-y1)/(curInjmax-curInjmin);
  
	TLine *lineCurzero = new TLine(xMin,(0.0-curInjmin)*slope + y1,xMax,(0.0-curInjmin)*slope + y1);
	lineCurzero->SetLineColor(kGray+1);
	lineCurzero->SetLineStyle(2);
	lineCurzero->SetLineWidth(1);
	lineCurzero->Draw();
       
	for(Int_t j=0;j<hCur1D[2]->GetNbinsX();j++) {
	  hCur1D[2]->SetBinContent(j+1,(hCur1D[2]->GetBinContent(j+1)-curInjmin)*slope + y1);
	}
	hCur1D[2]->SetLineStyle(1);
	hCur1D[2]->SetLineWidth(2);
    
	hCur1D[2]->SetLineColor(kGray+2);
    
	hCur1D[2]->Draw("sameC");
      } 

    }

  
    if(!opt.Contains("notim"))
      textTime->Draw();
    
    if(!opt.Contains("noden") && opt.Contains("units") && n0) 
      textDen->Draw();

    pad[ip]->RedrawAxis("g"); 

    ip--;  
    C->cd(0);
  }    


  if(mask & 0x2) {

    pad[ip]->Draw();
    pad[ip]->cd(); 
    // if(opt.Contains("logz")) {
    //   pad[ip]->SetLogz(1);
    // } else {
    //   pad[ip]->SetLogz(0);
    // }

    hFrame[ip]->Draw("col");

    //  hE2D[0]->GetZaxis()->SetNdivisions(503);
    hE2D[0]->GetZaxis()->SetTitleFont(fonttype);
    Float_t xFactor = pad[0]->GetAbsWNDC()/pad[ip]->GetAbsWNDC();
    Float_t yFactor = pad[0]->GetAbsHNDC()/pad[ip]->GetAbsHNDC();
    hE2D[0]->GetZaxis()->SetTickLength(xFactor*tylength/yFactor);
  
    exField->Draw();
    hE2D[0]->Draw(drawopt);

    if(opt.Contains("1dline")) {
      lineYzero->Draw();
    }


    if(opt.Contains("sline")) {
      if(zStartPlasma>xMin && zStartPlasma<xMax)
	lineStartPlasma->Draw();
      if(zStartNeutral>xMin && zStartNeutral<xMax)
	lineStartNeutral->Draw();
      if(zEndNeutral>xMin && zEndNeutral<xMax)
	lineEndNeutral->Draw();
    }
   
    // Fit the 1D line in the 2D pad:
    TLine *lineE = NULL;
    TLine *lineEmaxLuLinear = NULL;
    TLine *lineEmaxUSE = NULL;
    TLine *lineEmaxLo = NULL;
    TLine *lineEmaxLoVert = NULL;
    TLine *lineEmaxBeam = NULL;
    {
      Float_t rightmin = Emin[0];
      Float_t rightmax = Emax[0];
      Float_t slope = (yMax-yMin)/(rightmax-rightmin);
    
      for(Int_t j=0;j<hE1D[0]->GetNbinsX();j++) {
	hE1D[0]->SetBinContent(j+1,slope*(hE1D[0]->GetBinContent(j+1)-rightmin)+yMin);
      }
    
      Double_t sl = EzMaxLo/radius;
      Double_t zmatch = (EzMaxBeamLo/sl)+EzCross[0];

      Double_t EzAdapted  = slope*(EzMaxLo-rightmin)+yMin;
      Double_t EzAdapted2 = slope*(EzMaxBeamLo-rightmin)+yMin;
      lineE = new TLine(EzCross[0]-radius,-EzAdapted,zmatch,EzAdapted2);

      Double_t linelength = 1;
      if(opt.Contains("units")) linelength *= skindepth / PUnits::um;
 
      EzAdapted = slope*(EzMaxLuLinear-rightmin)+yMin;
      lineEmaxLuLinear = new TLine(EzCross[0]-radius,-EzAdapted,EzCross[0]-radius-linelength,-EzAdapted);

      EzAdapted = slope*(EzMaxUse-rightmin)+yMin;
      lineEmaxUSE = new TLine(EzCross[0]-radius,EzAdapted,EzCross[0],EzAdapted);

      EzAdapted = slope*(EzMaxLo-rightmin)+yMin;
      lineEmaxLo = new TLine(EzCross[0]-radius,-EzAdapted,EzCross[0]-radius-linelength,-EzAdapted);
      lineEmaxLoVert = new TLine(EzCross[0]-radius,-EzAdapted,EzCross[0]-radius,0.0);

      EzAdapted = slope*(EzMaxBeamLo-rightmin)+yMin;
      lineEmaxBeam = new TLine(zmatch,EzAdapted,-zmatch,EzAdapted);
      
      if(opt.Contains("mark")) {
	Float_t e0 = slope*(0-rightmin)+yMin;

	Float_t zF = hE1D[0]->GetBinCenter(binFcross);
	Float_t EF = hE1D[0]->GetBinContent(binFcross);

	TLine *lineEF = new TLine(zF,0,zF,EF);
	lineEF->SetLineColor(kRed);
	lineEF->SetLineStyle(2);
	lineEF->Draw();

	TMarker *markEF = new TMarker(zF,EF,20);
	markEF->SetMarkerColor(kRed);
	markEF->SetMarkerSize(1.4);
	markEF->Draw();

	Float_t zB = zEbMax;
	Float_t EB = slope*(EzExtr[0]-rightmin)+yMin;
      
	TLine *lineEB = new TLine(zB,0,zB,EB);
	lineEB->SetLineColor(kRed);
	lineEB->SetLineStyle(2);
	lineEB->Draw();

	TMarker *markEB = new TMarker(zB,EB,20);
	markEB->SetMarkerColor(kRed);
	markEB->SetMarkerSize(1.4);
	markEB->Draw();

	Float_t z0 = EzCross[0];
      
	TMarker *markE0 = new TMarker(z0,e0,20);
	markE0->SetMarkerColor(kRed);
	markE0->SetMarkerSize(1.4);
	markE0->Draw();

	Float_t zI = hE1D[0]->GetBinCenter(binPotValueIni);
	Float_t EI = hE1D[0]->GetBinContent(binPotValueIni);
	TLine *lineEI = new TLine(zI,e0,zI,EI);
	lineEI->SetLineColor(kRed);
	lineEI->SetLineStyle(2);
	lineEI->Draw();

	TMarker *markEI = new TMarker(zI,EI,20);
	markEI->SetMarkerColor(kRed);
	markEI->SetMarkerSize(1.4);
	markEI->Draw();
      }


    }
    
    if(opt.Contains("zero")) {
      TLine *lineZero = new TLine(xMin,0,xMax,0);
      lineZero->SetLineColor(lineColor);
      lineZero->SetLineStyle(2);
      lineZero->Draw();
    }

    if(opt.Contains("off")) {
      TLine *lineOff = new TLine(xMin,0,xMax,0);
      lineOff->SetLineColor(lineColor);
      lineOff->SetLineStyle(2);
      lineOff->Draw();
    }


    if(opt.Contains("bopar")) {

      lineE->SetLineColor(kGray+3);
      lineE->SetLineStyle(3);
      lineE->Draw();    
      
      lineEmaxLuLinear->SetLineColor(kGray+3);
      lineEmaxLuLinear->SetLineStyle(4);
      //lineEmaxLuLinear->Draw();    
      
      // lineEmaxUSE->SetLineColor(kGray+3);
      // lineEmaxUSE->SetLineStyle(4);
      // lineEmaxUSE->Draw();    
      
      lineEmaxLo->SetLineColor(kGray+3);
      lineEmaxLo->SetLineStyle(3);
      lineEmaxLo->Draw();    
      
      lineEmaxLoVert->SetLineColor(kGray+3);
      lineEmaxLoVert->SetLineStyle(3);
      lineEmaxLoVert->Draw();    
      
      lineEmaxBeam->SetLineColor(kGray+3);
      lineEmaxBeam->SetLineStyle(3);
      lineEmaxBeam->Draw();    


      TEllipse *ellipLo = new TEllipse(EzCross[0],0.0,radius);	
      ellipLo->SetFillStyle(0);
      ellipLo->SetLineStyle(3);
      ellipLo->SetLineColor(kGray+3);
      ellipLo->Draw();
	
      TEllipse *ellipLu = new TEllipse(EzCross[0],0.0,radiusLu);
      ellipLu->SetFillStyle(0);
      ellipLu->SetLineStyle(4);
      ellipLu->SetLineColor(kGray+3);
      ellipLu->Draw();


    }

    hE1D[0]->SetLineStyle(1);
    hE1D[0]->SetLineWidth(2);
    hE1D[0]->SetLineColor(lineColor);

    hE1D[0]->Draw("sameL");    
  
    pad[ip]->Update();
  
    y1 = pad[ip]->GetBottomMargin();
    y2 = 1 - pad[ip]->GetTopMargin();
    x1 = pad[ip]->GetLeftMargin();
    x2 = 1 - pad[ip]->GetRightMargin();
  
    palette = (TPaletteAxis*) hE2D[0]->GetListOfFunctions()->FindObject("palette");  
    if(palette) {
      palette->SetY2NDC(y2 - gap);
      palette->SetY1NDC(y1 + gap);
      palette->SetX1NDC(x2 + 0.005);
      palette->SetX2NDC(x2 + 0.03);
      palette->SetTitleOffset(tzoffset);
      palette->SetTitleSize(tzsize);
      palette->SetLabelFont(fonttype);
      palette->SetLabelSize(lzsize);
      palette->SetLabelOffset(lyoffset);
      palette->SetBorderSize(2);
      palette->SetLineColor(1);
    }
   
    pad[ip]->RedrawAxis(); 

    if(opt.Contains("blines")) {
      TLine *lineEbunch = NULL;
      TLine *lineESlope = NULL;
      TLine *lineZbunch = NULL;
      if(zInjBeam>-999) { 
	Float_t slope = (yMax-yMin)/(Emax[0]-Emin[0]);

	lineEbunch = new TLine(xMin,EInjBeam*slope,xMax,EInjBeam*slope);
	lineEbunch->SetLineColor(kGray+2);
	lineEbunch->SetLineStyle(2);
	lineEbunch->SetLineWidth(1);
	lineEbunch->Draw();
      
	lineZbunch = new TLine(zInjBeam,yMin,zInjBeam,yMax);
	lineZbunch->SetLineColor(kGray+2);
	lineZbunch->SetLineStyle(2);
	lineZbunch->SetLineWidth(1);
	lineZbunch->Draw();
      
	Float_t z1 = zInjBeam - 10;
	Float_t z2 = zInjBeam + 10;
	Float_t E1 = EInjBeam + ESlopeInjBeam * (z1-zInjBeam);
	Float_t E2 = EInjBeam + ESlopeInjBeam * (z2-zInjBeam);
	lineESlope = new TLine(z1,E1*slope,z2,E2*slope);
	lineESlope->SetLineColor(lineColor);
	lineESlope->SetLineStyle(1);
	lineESlope->SetLineWidth(1);
	//lineESlope->Draw();    
      
      }
    }

  
  
    if(opt.Contains("1dline")) {
      lineYzero->Draw();
    }

    if(opt.Contains("sline")) {
      if(zStartPlasma>xMin && zStartPlasma<xMax)
	lineStartPlasma->Draw();
      if(zStartNeutral>xMin && zStartNeutral<xMax)
	lineStartNeutral->Draw();
      if(zEndNeutral>xMin && zEndNeutral<xMax)
	lineEndNeutral->Draw();
    }
 
    ip--;  
    C->cd(0);
  }



  if(mask & 0x4) {

    pad[ip]->Draw();
    pad[ip]->cd(); // <---------------------------------------------------------- 2nd panel
    // if(opt.Contains("logz")) {
    //   pad[ip]->SetLogz(1);
    // } else {
    //   pad[ip]->SetLogz(0);
    // }

    hFrame[ip]->Draw("col");

    //  hE2D[0]->GetZaxis()->SetNdivisions(503);
    hE2D[1]->GetZaxis()->SetTitleFont(fonttype);
    Float_t xFactor = pad[0]->GetAbsWNDC()/pad[ip]->GetAbsWNDC();
    Float_t yFactor = pad[0]->GetAbsHNDC()/pad[ip]->GetAbsHNDC();
    hE2D[1]->GetZaxis()->SetTickLength(xFactor*tylength/yFactor);
  
    exField->Draw();
    hE2D[1]->Draw("colz0 same");

    if(opt.Contains("bopar")) {
      
      TEllipse *ellipLo = new TEllipse(EzCross[0],0.0,radius);	
      ellipLo->SetFillStyle(0);
      ellipLo->SetLineStyle(3);
      ellipLo->SetLineColor(kGray+3);
      ellipLo->Draw();
	
      TEllipse *ellipLu = new TEllipse(EzCross[0],0.0,radiusLu);
      ellipLu->SetFillStyle(0);
      ellipLu->SetLineStyle(4);
      ellipLu->SetLineColor(kGray+3);
      ellipLu->Draw();

    }


    if(opt.Contains("1dline")) {
      lineYzero->Draw();
    }

    if(opt.Contains("sline")) {
      if(zStartPlasma>xMin && zStartPlasma<xMax)
	lineStartPlasma->Draw();
      if(zStartNeutral>xMin && zStartNeutral<xMax)
	lineStartNeutral->Draw();
      if(zEndNeutral>xMin && zEndNeutral<xMax)
	lineEndNeutral->Draw();
    }
   
    // Fit the 1D line in the 2D pad:
    {
      Float_t rightmin = Emin[1];
      Float_t rightmax = Emax[1];
      Float_t slope = (yMax-yMin)/(rightmax-rightmin);
    
      for(Int_t j=0;j<hE1D[1]->GetNbinsX();j++) {
	hE1D[1]->SetBinContent(j+1,slope*(hE1D[1]->GetBinContent(j+1)-rightmin)+yMin);
      }
    }


    if(opt.Contains("zero")) {
      TLine *lineZero = new TLine(xMin,0,xMax,0);
      lineZero->SetLineColor(lineColor);
      lineZero->SetLineStyle(2);
      lineZero->Draw();
    }

    if(opt.Contains("off")) {
      TLine *lineOff = new TLine(xMin,xoff,xMax,xoff);
      lineOff->SetLineColor(lineColor);
      lineOff->SetLineStyle(2);
      lineOff->Draw();
    }


    hE1D[1]->SetLineStyle(1);
    hE1D[1]->SetLineWidth(2);
    hE1D[1]->SetLineColor(lineColor);

    hE1D[1]->Draw("sameL");
  
    pad[ip]->Update();
  
    y1 = pad[ip]->GetBottomMargin();
    y2 = 1 - pad[ip]->GetTopMargin();
    x1 = pad[ip]->GetLeftMargin();
    x2 = 1 - pad[ip]->GetRightMargin();
  
    palette = (TPaletteAxis*) hE2D[1]->GetListOfFunctions()->FindObject("palette");  
    if(palette) {
      palette->SetY2NDC(y2 - gap);
      palette->SetY1NDC(y1 + gap);
      palette->SetX1NDC(x2 + 0.005);
      palette->SetX2NDC(x2 + 0.03);
      palette->SetTitleOffset(tzoffset);
      palette->SetTitleSize(tzsize);
      palette->SetLabelFont(fonttype);
      palette->SetLabelSize(lzsize);
      palette->SetLabelOffset(lyoffset);
      palette->SetBorderSize(2);
      palette->SetLineColor(1);
    }
   
    pad[ip]->RedrawAxis(); 
  
  
    if(opt.Contains("1dline")) {
      lineYzero->Draw();
    }

    if(opt.Contains("sline")) {
      if(zStartPlasma>xMin && zStartPlasma<xMax)
	lineStartPlasma->Draw();
      if(zStartNeutral>xMin && zStartNeutral<xMax)
	lineStartNeutral->Draw();
      if(zEndNeutral>xMin && zEndNeutral<xMax)
	lineEndNeutral->Draw();
    }
 
    ip--;  
    C->cd(0);
  }
 

  if(mask & 0x8) {

    pad[ip]->Draw();
    pad[ip]->cd(); // <---------------------------------------------------------- 2nd panel
    // if(opt.Contains("logz")) {
    //   pad[ip]->SetLogz(1);
    // } else {
    //   pad[ip]->SetLogz(0);
    // }

    hFrame[ip]->Draw("col");

    //  hE2D[2]->GetZaxis()->SetNdivisions(503);
    hE2D[2]->GetZaxis()->SetTitleFont(fonttype);
    Float_t xFactor = pad[0]->GetAbsWNDC()/pad[ip]->GetAbsWNDC();
    Float_t yFactor = pad[0]->GetAbsHNDC()/pad[ip]->GetAbsHNDC();
    hE2D[2]->GetZaxis()->SetTickLength(xFactor*tylength/yFactor);
  
    exField->Draw();
    hE2D[2]->Draw(drawopt);

    if(opt.Contains("1dline")) {
      lineYzero->Draw();
    }

    if(opt.Contains("sline")) {
      if(zStartPlasma>xMin && zStartPlasma<xMax)
	lineStartPlasma->Draw();
      if(zStartNeutral>xMin && zStartNeutral<xMax)
	lineStartNeutral->Draw();
      if(zEndNeutral>xMin && zEndNeutral<xMax)
	lineEndNeutral->Draw();
    }
   
    // Fit the 1D line in the 2D pad:
    {
      Float_t rightmin = Emin[2];
      Float_t rightmax = Emax[2];
      Float_t slope = (yMax-yMin)/(rightmax-rightmin);
    
      for(Int_t j=0;j<hE1D[2]->GetNbinsX();j++) {
	hE1D[2]->SetBinContent(j+1,slope*(hE1D[2]->GetBinContent(j+1)-rightmin)+yMin);
      }
    
    }


    if(opt.Contains("zero")) {
      TLine *lineZero = new TLine(xMin,0,xMax,0);
      lineZero->SetLineColor(lineColor);
      lineZero->SetLineStyle(2);
      lineZero->Draw();
    }

    if(opt.Contains("off")) {
      TLine *lineOff = new TLine(xMin,xoff,xMax,xoff);
      lineOff->SetLineColor(lineColor);
      lineOff->SetLineStyle(2);
      lineOff->Draw();
    }

    hE1D[2]->SetLineStyle(1);
    hE1D[2]->SetLineWidth(2);
    hE1D[2]->SetLineColor(lineColor);

    hE1D[2]->Draw("sameL");
  
    pad[ip]->Update();
  
    y1 = pad[ip]->GetBottomMargin();
    y2 = 1 - pad[ip]->GetTopMargin();
    x1 = pad[ip]->GetLeftMargin();
    x2 = 1 - pad[ip]->GetRightMargin();
  
    palette = (TPaletteAxis*) hE2D[2]->GetListOfFunctions()->FindObject("palette");  
    if(palette) {
      palette->SetY2NDC(y2 - gap);
      palette->SetY1NDC(y1 + gap);
      palette->SetX1NDC(x2 + 0.005);
      palette->SetX2NDC(x2 + 0.03);
      palette->SetTitleOffset(tzoffset);
      palette->SetTitleSize(tzsize);
      palette->SetLabelFont(fonttype);
      palette->SetLabelSize(lzsize);
      palette->SetLabelOffset(lyoffset);
      palette->SetBorderSize(2);
      palette->SetLineColor(1);
    }
   
    pad[ip]->RedrawAxis(); 
  
  
    if(opt.Contains("1dline")) {
      lineYzero->Draw();
    }

    if(opt.Contains("sline")) {
      if(zStartPlasma>xMin && zStartPlasma<xMax)
	lineStartPlasma->Draw();
      if(zStartNeutral>xMin && zStartNeutral<xMax)
	lineStartNeutral->Draw();
      if(zEndNeutral>xMin && zEndNeutral<xMax)
	lineEndNeutral->Draw();
    }
 
    ip--;
    C->cd(0);
  }

  if(mask & 0x10) {

    pad[ip]->Draw();
    pad[ip]->cd(); // <---------------------------------------------------------- 2nd panel
    // if(opt.Contains("logz")) {
    //   pad[ip]->SetLogz(1);
    // } else {
    //   pad[ip]->SetLogz(0);
    // }

    hFrame[ip]->Draw("col");

    //  hFocus2D->GetZaxis()->SetNdivisions(503);
    hFocus2D->GetZaxis()->SetTitleFont(fonttype);
    Float_t xFactor = pad[0]->GetAbsWNDC()/pad[ip]->GetAbsWNDC();
    Float_t yFactor = pad[0]->GetAbsHNDC()/pad[ip]->GetAbsHNDC();
    hFocus2D->GetZaxis()->SetTickLength(xFactor*tylength/yFactor);
  
    exField->Draw();
    hFocus2D->Draw(drawopt);


    if(opt.Contains("bopar")) {
      
      TEllipse *ellipLo = new TEllipse(EzCross[0],0.0,radius);	
      ellipLo->SetFillStyle(0);
      ellipLo->SetLineStyle(3);
      ellipLo->SetLineColor(kGray+3);
      ellipLo->Draw();
	
      TEllipse *ellipLu = new TEllipse(EzCross[0],0.0,radiusLu);
      ellipLu->SetFillStyle(0);
      ellipLu->SetLineStyle(4);
      ellipLu->SetLineColor(kGray+3);
      ellipLu->Draw();

    }


    if(opt.Contains("1dline")) {
      lineYzero->Draw();
    }

    if(opt.Contains("sline")) {
      if(zStartPlasma>xMin && zStartPlasma<xMax)
	lineStartPlasma->Draw();
      if(zStartNeutral>xMin && zStartNeutral<xMax)
	lineStartNeutral->Draw();
      if(zEndNeutral>xMin && zEndNeutral<xMax)
	lineEndNeutral->Draw();
    }
  
    if(opt.Contains("zero")) {
      TLine *lineZero = new TLine(xMin,0,xMax,0);
      //lineZero->SetLineColor(lineColor);
      lineZero->SetLineColor(kGray);
      lineZero->SetLineStyle(2);
      lineZero->Draw();
    }


    if(opt.Contains("off")) {
      TLine *lineOff = new TLine(xMin,xoff,xMax,xoff);
      lineOff->SetLineColor(lineColor);
      lineOff->SetLineStyle(2);
      lineOff->Draw();
    }

 
    // Fit the 1D line in the 2D pad:
    {
      Float_t rightmin = Fmin;
      Float_t rightmax = Fmax;
      Float_t slope = (yMax-yMin)/(rightmax-rightmin);
    
      for(Int_t j=0;j<hFocus1D->GetNbinsX();j++) {
	hFocus1D->SetBinContent(j+1,slope*(hFocus1D->GetBinContent(j+1)-rightmin)+yMin);
      }
      
      if(opt.Contains("mark")) {
	
	Float_t zF = hFocus1D->GetBinCenter(binFcross);
	Float_t FF = slope*(0.0-rightmin)+yMin;
	
	TMarker *markFF = new TMarker(zF,FF,20);
	markFF->SetMarkerColor(kRed);
	markFF->SetMarkerSize(1.4);
	markFF->Draw();
	
      }
      
    }
  
    hFocus1D->SetLineStyle(1);
    hFocus1D->SetLineWidth(2);
    hFocus1D->SetLineColor(lineColor);

    hFocus1D->Draw("sameL");
  
    pad[ip]->Update();
  
    y1 = pad[ip]->GetBottomMargin();
    y2 = 1 - pad[ip]->GetTopMargin();
    x1 = pad[ip]->GetLeftMargin();
    x2 = 1 - pad[ip]->GetRightMargin();
  
    palette = (TPaletteAxis*) hFocus2D->GetListOfFunctions()->FindObject("palette");  
    if(palette) {
      palette->SetY2NDC(y2 - gap);
      palette->SetY1NDC(y1 + gap);
      palette->SetX1NDC(x2 + 0.005);
      palette->SetX2NDC(x2 + 0.03);
      palette->SetTitleOffset(tzoffset);
      palette->SetTitleSize(tzsize);
      palette->SetLabelFont(fonttype);
      palette->SetLabelSize(lzsize);
      palette->SetLabelOffset(lyoffset);
      palette->SetBorderSize(2);
      palette->SetLineColor(1);
    }
   
    pad[ip]->RedrawAxis(); 
  
  
    if(opt.Contains("1dline")) {
      lineYzero->Draw();
    }

    if(opt.Contains("sline")) {
      if(zStartPlasma>xMin && zStartPlasma<xMax)
	lineStartPlasma->Draw();
      if(zStartNeutral>xMin && zStartNeutral<xMax)
	lineStartNeutral->Draw();
      if(zEndNeutral>xMin && zEndNeutral<xMax)
	lineEndNeutral->Draw();
    }
 
    ip--;
    C->cd(0);
  }

  


  if(mask & 0x20) {
    
    pad[ip]->Draw();
    pad[ip]->cd(); // <--------------------------------------------------- 3rd panel
    // if(opt.Contains("logz")) {
    //   pad[ip]->SetLogz(1);
    // } else {
    //   pad[ip]->SetLogz(0);
    // }

    hFrame[ip]->Draw("col");

    //  hETotal2D->GetZaxis()->SetNdivisions(503);
    hETotal2D->GetZaxis()->SetTitleFont(fonttype);
    Float_t xFactor = pad[0]->GetAbsWNDC()/pad[ip]->GetAbsWNDC();
    Float_t yFactor = pad[0]->GetAbsHNDC()/pad[ip]->GetAbsHNDC();
    hETotal2D->GetZaxis()->SetTickLength(xFactor*tylength/yFactor);

    exFieldT->Draw();
  
    hETotal2D->Draw(drawopt);

    // ADK contours  
    // for(Int_t i=0;i<graphsI2D.GetEntriesFast();i++) {
    //   TGraph *gr = (TGraph*) graphsI2D.At(i);
    //   if(!gr) continue;
    //   gr->Draw("C");
    // }

    
    // Fit the E1D in the E2D pad:
    Float_t rightmin = ETmin;
    Float_t rightmax = ETmax;
    Float_t slope = (yMax-yMin)/(rightmax-rightmin);
    for(Int_t j=0;j<hETotal1D->GetNbinsX();j++) {
      hETotal1D->SetBinContent(j+1,slope*(hETotal1D->GetBinContent(j+1)-rightmin)+yMin);
    }


    // if(opt.Contains("zero")) {
    //   TLine *lineZero = new TLine(xMin,0,xMax,0);
    //   lineZero->SetLineColor(lineColor);
    //   lineZero->SetLineStyle(2);
    //   lineZero->Draw();
    // }

    if(opt.Contains("off")) {
      TLine *lineOff = new TLine(xMin,xoff,xMax,xoff);
      lineOff->SetLineColor(lineColor);
      lineOff->SetLineStyle(2);
      lineOff->Draw();
    }

    
    if(opt.Contains("ionth")) {
      // const Int_t Nat = 6;
      // // char atNames[Nat][8] = {"H", "He","He^{+}", "Ne", "Ne+","Custom"};
      // Float_t IonTh[Nat] = {33.79,92.77,234.97,74.42,143.69,500.59};
      // TLine *lineTh[Nat];
      // Int_t atColor[Nat] = {kAzure+1,kOrange+8,kOrange+8,kSpring-1,kSpring-1,kCyan};
      // Int_t atStyle[Nat] = {1,1,2,1,2,2};
      const Int_t Nat = 5;
      // char atNames[Nat][8] = {"H", "He","He^{+}", "Ne", "Ne+"};
      Float_t IonTh[Nat] = {33.79,92.77,234.97,74.42,143.69};
      TLine *lineTh[Nat];
      Int_t atColor[Nat] = {kAzure+1,kOrange+8,kOrange+8,kSpring-1,kSpring-1};
      Int_t atStyle[Nat] = {1,1,2,1,2};
      
      for(Int_t i=0;i<Nat;i++) {
	if(!opt.Contains("units") || !n0 ) 
	  IonTh[i] /= ( E0 / (PUnits::GV/PUnits::m));
	
	lineTh[i] = new TLine(xMin,slope*(IonTh[i]-rightmin)+yMin,
			      xMax,slope*(IonTh[i]-rightmin)+yMin);

	lineTh[i]->SetLineColor(atColor[i]);
	lineTh[i]->SetLineStyle(atStyle[i]);
	lineTh[i]->Draw();


      }
    }

    if(opt.Contains("mark")) {

      Float_t zI = hETotal1D->GetBinCenter(binPotValueIni);
      Float_t EI = hETotal1D->GetBinContent(binPotValueIni);
      Float_t e0 = slope*(0-rightmin)+yMin;
 
      TLine *lineEI = new TLine(zI,e0,zI,EI);
      lineEI->SetLineColor(kRed);
      lineEI->SetLineStyle(2);
      lineEI->Draw();

      TMarker *markEI = new TMarker(zI,EI,20);
      markEI->SetMarkerColor(kRed);
      markEI->SetMarkerSize(1.4);
      markEI->Draw();

      Float_t zB = zEbMax;
      Float_t EB = hETotal1D->GetBinContent(binEbMax);
      
      TLine *lineEB = new TLine(zB,e0,zB,EB);
      lineEB->SetLineColor(kRed);
      lineEB->SetLineStyle(2);
      lineEB->Draw();

      TMarker *markEB = new TMarker(zB,EB,20);
      markEB->SetMarkerColor(kRed);
      markEB->SetMarkerSize(1.4);
      markEB->Draw();

    }

    if(opt.Contains("bopar")) {
      
      TEllipse *ellipLo = new TEllipse(EzCross[0],0.0,radius);	
      ellipLo->SetFillStyle(0);
      ellipLo->SetLineStyle(3);
      ellipLo->SetLineColor(kGray+3);
      ellipLo->Draw();
	
      TEllipse *ellipLu = new TEllipse(EzCross[0],0.0,radiusLu);
      ellipLu->SetFillStyle(0);
      ellipLu->SetLineStyle(4);
      ellipLu->SetLineColor(kGray+3);
      ellipLu->Draw();

    }

  
    hETotal1D->SetLineStyle(1);
    hETotal1D->SetLineWidth(2);
    hETotal1D->SetLineColor(lineColor);

    hETotal1D->Draw("sameL");
    
    if(opt.Contains("1dline")) {
      lineYzero->Draw();
    }

    if(opt.Contains("sline")) {
      if(zStartPlasma>xMin && zStartPlasma<xMax)
	lineStartPlasma->Draw();
      if(zStartNeutral>xMin && zStartNeutral<xMax)
	lineStartNeutral->Draw();
      if(zEndNeutral>xMin && zEndNeutral<xMax)
	lineEndNeutral->Draw();
    }
  
    
    if(opt.Contains("cont")) {
      // PSI MAIN contour  
      for(Int_t i=0;i<graphsV2D_main.GetEntriesFast();i++) {
	TGraph *gr = (TGraph*) graphsV2D_main.At(i);
	if(!gr) continue;
	// gr->Draw("C");
      }
    }
    
  
    pad[ip]->Update();
  
    y1 = pad[ip]->GetBottomMargin();
    y2 = 1 - pad[ip]->GetTopMargin();
    x1 = pad[ip]->GetLeftMargin();
    x2 = 1 - pad[ip]->GetRightMargin();
  
    palette = (TPaletteAxis*) hETotal2D->GetListOfFunctions()->FindObject("palette");  
    if(palette) {
      palette->SetY2NDC(y2 - gap);
      palette->SetY1NDC(y1 + gap);
      palette->SetX1NDC(x2 + 0.005);
      palette->SetX2NDC(x2 + 0.03);
      palette->SetTitleOffset(tzoffset);
      palette->SetTitleSize(tzsize);
      palette->SetLabelFont(fonttype);
      palette->SetLabelSize(lzsize);
      palette->SetLabelOffset(lyoffset);
      palette->SetBorderSize(2);
      palette->SetLineColor(1);
    }
   
    pad[ip]->RedrawAxis(); 

    ip--;
    C->cd(0);

  }

  if(mask & 0x40) {

    pad[ip]->Draw();
    pad[ip]->cd(); // <--------------------------------------------------------- 5th panel
 
    hFrame[ip]->Draw("col");

    hV2D->GetZaxis()->SetTitleFont(fonttype);
    Float_t xFactor = pad[0]->GetAbsWNDC()/pad[ip]->GetAbsWNDC();
    Float_t yFactor = pad[0]->GetAbsHNDC()/pad[ip]->GetAbsHNDC();
    hV2D->GetZaxis()->SetTickLength(xFactor*tylength/yFactor);
 
    exPot->Draw();
  
    hV2D->GetZaxis()->SetRangeUser(Vmin,Vmax);
    hV2D->Draw("col z same");
  
    for(Int_t i=0;i<graphsV2D.GetEntriesFast();i++) {
      TGraph *gr = (TGraph*) graphsV2D.At(i);
      if(!gr) continue;
      gr->Draw("C");
    }
  
    // PSI MAIN contour  
    // for(Int_t i=0;i<graphsV2D_main.GetEntriesFast();i++) {
    //   TGraph *gr = (TGraph*) graphsV2D_main.At(i);
    //   if(!gr) continue;
    //   gr->Draw("C");
    // }

    // Fit the V1D in the V2D pad:
    TLine *lineT = NULL;
    {
      Float_t rightmin = Vmin;
      Float_t rightmax = Vmax;
      Float_t slope = (yMax-yMin)/(rightmax-rightmin);
      Float_t psi0 = slope*(0.0-rightmin)+yMin;

      for(Int_t j=0;j<hV1D->GetNbinsX();j++) {
	hV1D->SetBinContent(j+1,slope*(hV1D->GetBinContent(j+1)-rightmin)+yMin);
      }

      if(opt.Contains("zero")) {
	lineT = new TLine(xMin,psi0,xMax,psi0);

	
	lineT->SetLineColor(lineColor);
	lineT->SetLineStyle(2);
	lineT->Draw();
      }
    
    

      if(opt.Contains("mark")) {
	Float_t zI = hV1D->GetBinCenter(binPotValueIni);
	Float_t psiI = hV1D->GetBinContent(binPotValueIni);
	Float_t zF = hV1D->GetBinCenter(binFcross);
	Float_t psiF = hV1D->GetBinContent(binFcross);
	
	TLine *linePsiI = new TLine(zI,psiF,zI,psiI);
	linePsiI->SetLineColor(kRed);
	linePsiI->SetLineStyle(2);
	linePsiI->Draw();
      
	TMarker *markPsiI = new TMarker(zI,psiI,20);
	markPsiI->SetMarkerColor(kRed);
	markPsiI->SetMarkerSize(1.4);
	markPsiI->Draw();

		
	TLine *linePsiF = new TLine(zF,psiF,zF,psiF);
	linePsiF->SetLineColor(kRed);
	linePsiF->SetLineStyle(2);
	linePsiF->Draw();
      
	TMarker *markPsiF = new TMarker(zF,psiF,20);
	markPsiF->SetMarkerColor(kRed);
	markPsiF->SetMarkerSize(1.4);
	markPsiF->Draw();
	
      }
    }

    hV1D->SetLineStyle(1);
    hV1D->SetLineWidth(2);
    //  hV1D->SetLineColor(PGlobals::elecLine);
    hV1D->SetLineColor(lineColor);

    hV1D->Draw("sameL");
    
    // ADK MAIN contour  
    // for(Int_t i=0;i<graphsI2D_main.GetEntriesFast();i++) {
    //   TGraph *gr = (TGraph*) graphsI2D_main.At(i);
    //   if(!gr) continue;
    //   gr->Draw("C");
    // }

    if(opt.Contains("bopar")) {
      TEllipse *ellipTrap = new TEllipse(EzCross[0],0.0,RadiusPsiLo);           
      ellipTrap->SetFillStyle(0);
      ellipTrap->SetLineStyle(3);
      ellipTrap->SetLineColor(kWhite);
      //ellipTrap->Draw();

      TEllipse *ellipLo = new TEllipse(EzCross[0],0.0,radius);	
      ellipLo->SetFillStyle(0);
      ellipLo->SetLineStyle(3);
      ellipLo->SetLineColor(kGray+3);
      ellipLo->Draw();

      TEllipse *ellipLu = new TEllipse(EzCross[0],0.0,radiusLu);
      ellipLu->SetFillStyle(0);
      ellipLu->SetLineStyle(4);
      ellipLu->SetLineColor(kGray+3);
      ellipLu->Draw();


    }

    if(opt.Contains("1dline")) {
      lineYzero->Draw();
    }

    if(opt.Contains("sline")) {
      if(zStartPlasma>xMin && zStartPlasma<xMax)
	lineStartPlasma->Draw();
      if(zStartNeutral>xMin && zStartNeutral<xMax)
	lineStartNeutral->Draw();
      if(zEndNeutral>xMin && zEndNeutral<xMax)
	lineEndNeutral->Draw();
    }
  
    
    pad[ip]->Update();

    y1 = pad[ip]->GetBottomMargin();
    y2 = 1 - pad[ip]->GetTopMargin();
    x1 = pad[ip]->GetLeftMargin();
    x2 = 1 - pad[ip]->GetRightMargin();
 
    palette = (TPaletteAxis*)hV2D->GetListOfFunctions()->FindObject("palette");
    if(palette) {
      palette->SetY2NDC(y2 - gap);
      palette->SetY1NDC(y1 + gap);
      palette->SetX1NDC(x2 + 0.005);
      palette->SetX2NDC(x2 + 0.03);
      palette->SetTitleOffset(tzoffset);
      // palette->SetTitleOffset(999.9);
      palette->SetTitleSize(tzsize);
      palette->SetLabelFont(fonttype);
      palette->SetLabelSize(lzsize);
      palette->SetLabelOffset(lyoffset);
      palette->SetBorderSize(2);
      palette->SetLineColor(1);
    }

    pad[ip]->RedrawAxis(); 
  
    ip--;
    C->cd();
  }


  if(mask & 0x80) {
    pad[ip]->Draw();
    pad[ip]->cd(); // <---------------------------------------------------------- 4th panel
    // if(opt.Contains("logz")) {
    //   pad[ip]->SetLogz(1);
    // } else {
    //   pad[ip]->SetLogz(0);
    // }

    hFrame[ip]->Draw("col");

    hIonProb2D[ii]->GetZaxis()->SetNdivisions(503);
    hIonProb2D[ii]->GetZaxis()->SetTitleFont(fonttype);
    Float_t xFactor = pad[0]->GetAbsWNDC()/pad[ip]->GetAbsWNDC();
    Float_t yFactor = pad[0]->GetAbsHNDC()/pad[ip]->GetAbsHNDC();
    hIonProb2D[ii]->GetZaxis()->SetTickLength(xFactor*tylength/yFactor);
  
    exIonP->Draw();
    hIonProb2D[ii]->Draw(drawopt);

    // ADK contours  
    for(Int_t i=0;i<graphsI2D.GetEntriesFast();i++) {
      TGraph *gr = (TGraph*) graphsI2D.At(i);
      if(!gr) continue;
      gr->Draw("C");
    }

    {
      Float_t rightmin = IPmin;
      Float_t rightmax = IPmax;
      Float_t slope = (yMax-yMin)/(rightmax-rightmin);
    
      for(Int_t j=0;j<hIonProb1D[ii]->GetNbinsX();j++) {
	hIonProb1D[ii]->SetBinContent(j+1,slope*(hIonProb1D[ii]->GetBinContent(j+1)-rightmin)+yMin);
      }
    }
  
    hIonProb1D[ii]->SetLineStyle(1);
    hIonProb1D[ii]->SetLineWidth(2);
    hIonProb1D[ii]->SetLineColor(lineColor);

    hIonProb1D[ii]->Draw("sameL");

    if(opt.Contains("bopar")) {
      
      TEllipse *ellipLo = new TEllipse(EzCross[0],0.0,radius);	
      ellipLo->SetFillStyle(0);
      ellipLo->SetLineStyle(3);
      ellipLo->SetLineColor(kGray+3);
      ellipLo->Draw();
	
      TEllipse *ellipLu = new TEllipse(EzCross[0],0.0,radiusLu);
      ellipLu->SetFillStyle(0);
      ellipLu->SetLineStyle(4);
      ellipLu->SetLineColor(kGray+3);
      ellipLu->Draw();

    }

    if(opt.Contains("cont")) {
      // PSI MAIN contour  
      // for(Int_t i=0;i<graphsV2D_main.GetEntriesFast();i++) {
      // 	TGraph *gr = (TGraph*) graphsV2D_main.At(i);
      // 	if(!gr) continue;
      // 	gr->Draw("C");
      // }

      

    }

    // PSI contours  
    for(Int_t i=0;i<graphsV2D.GetEntriesFast();i++) {
      TGraph *gr = (TGraph*) graphsV2D.At(i);
      if(!gr) continue;
      gr->Draw("C");
    }


    if(opt.Contains("1dline")) {
      lineYzero->Draw();
    }

    if(opt.Contains("sline")) {
      if(zStartPlasma>xMin && zStartPlasma<xMax)
	lineStartPlasma->Draw();
      if(zStartNeutral>xMin && zStartNeutral<xMax)
	lineStartNeutral->Draw();
      if(zEndNeutral>xMin && zEndNeutral<xMax)
	lineEndNeutral->Draw();
    }
  
    pad[ip]->Update();
  
    y1 = pad[ip]->GetBottomMargin();
    y2 = 1 - pad[ip]->GetTopMargin();
    x1 = pad[ip]->GetLeftMargin();
    x2 = 1 - pad[ip]->GetRightMargin();
  
    palette = (TPaletteAxis*) hIonProb2D[ii]->GetListOfFunctions()->FindObject("palette");  
    if(palette) {
      palette->SetY2NDC(y2 - gap);
      palette->SetY1NDC(y1 + gap);
      palette->SetX1NDC(x2 + 0.005);
      palette->SetX2NDC(x2 + 0.03);
      palette->SetTitleOffset(tzoffset);
      palette->SetTitleSize(tzsize);
      palette->SetLabelFont(fonttype);
      palette->SetLabelSize(lzsize);
      palette->SetLabelOffset(lyoffset);
      palette->SetBorderSize(2);
      palette->SetLineColor(1);
    }
   
    pad[ip]->RedrawAxis(); 

    ip--;
    C->cd(0);
  }


  // Writing the figure
  TString fOutName = Form("./%s/Plots/Snapshots/Snapshot%s-%s_%i",pData->GetPath().c_str(),imask.c_str(),pData->GetName(),timestep);
  PGlobals::imgconv(C,fOutName,opt);
  // ---------------------------------------------------------

  //  PGlobals::DestroyCanvases();
}
 
 
