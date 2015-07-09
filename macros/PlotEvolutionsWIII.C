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
#include <TMarker.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TMath.h>
#include <TPaveText.h>
#include <TExec.h>
#include <TGaxis.h>
#include <TPaletteAxis.h>

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

void PlotEvolutionsWIII(const TString &sim, UInt_t mask = 3, const TString &options="png") { 
  
#ifdef __CINT__  
  gSystem->Load("libptools.so");
#endif

  string imask = DecToBin(mask);
  cout << Form("\n Plotting Evolultion with mask: %s",imask.c_str()) << endl; 

  PGlobals::Initialize();
  
  // Palettes!
  gROOT->Macro("PPalettes.C");

  TString opt = options;
  
  // More makeup            
  Float_t margins[4] = {0.15,0.15,0.20,0.10};
  gStyle->SetPadLeftMargin(margins[0]);  // Margin left axis  
  gStyle->SetPadRightMargin(margins[2]);
  gStyle->SetPadTopMargin(margins[3]);  // Margin left axis  
  gStyle->SetPadBottomMargin(margins[1]);

  gStyle->SetPadTickX(0);
  gStyle->SetPadTickY(0);


  if(opt.Contains("grid")) {
    gStyle->SetPadGridX(1);
    gStyle->SetPadGridY(1);
  }

  // Load first simulation data (for instance)
  PData *pData = PData::Get(sim.Data());
  Double_t E0 = pData->GetPlasmaE0();
  
  Float_t maxEcross = -999.;
  Float_t minEcross = 999.;
  Float_t maxEextr  = -999.;
  Float_t minEextr  = 999.;
  Float_t maxEdephas  = -999.;
  Float_t minEdephas  = 999.;
 
  const Int_t Nfields = 2; // E_z, E_x 
  TH2F   *hEvsTime[Nfields];
  Int_t NCross[Nfields];
  TGraph **gEcross[Nfields]; 
  TGraph **gEextr[Nfields]; 
  TGraph **gEdephas[Nfields]; 
  TGraph *gTRatio;

  TH2F   *hFvsTime;
  TGraph **gFcross = NULL;
  TGraph **gFextr = NULL;

  TH2F   *hETvsTime;
  TGraph **gETcross = NULL;
  TGraph **gETextr = NULL;

  TH2F   *hVvsTime;
  TGraph **gVcross = NULL;
  TGraph **gVextr = NULL;

  const Int_t NAtoms = 3;
  char atNames[NAtoms][4] = {"H","He","He2"};
  TH2F   *hIonProbvsTime[NAtoms]; // For H, He and He+. 
  TGraph *gIonProb10[NAtoms];
  TGraph *gIonProb100[NAtoms];
  Float_t IonTh[NAtoms] = {33.8,92.75,234.96} ;  // GV/m
  // if(!opt.Contains("units")) 
  //   for(Int_t i=0;i<NAtoms;i++) IonTh[i]  /= ( E0 / (PUnits::GV/PUnits::m));

  char hName[24];
  char gName[24];
  
  TString filename;
  filename = Form("./%s/Plots/Evolutions/Evolutions-%s.root",sim.Data(),sim.Data());
  
  TFile  *ifile = (TFile*) gROOT->GetListOfFiles()->FindObject(filename.Data());
  if (!ifile) ifile = new TFile(filename,"READ");

  TH2F *hDen1DvsTime = NULL;
  sprintf(hName,"hDenvsTime_1"); 
  hDen1DvsTime = (TH2F*) ifile->Get(hName);

  TH2F *hRmsvsTime = NULL;
  sprintf(hName,"hRmsvsTime_1"); 
  hRmsvsTime = (TH2F*) ifile->Get(hName);
    
  for(Int_t i=0;i<Nfields;i++) {
    
    sprintf(hName,"hEvsTime_%i",i); 
    hEvsTime[i] = (TH2F*) ifile->Get(hName);
    if(!hEvsTime[i]) continue;
    
    cout << Form("ANALYZING FIELD %i ...",i) << endl;

    Int_t NTBins = hEvsTime[i]->GetNbinsX();
    for(Int_t it=NTBins;it>0;it--) {
      
      // 1D field at certain timestep "it".
      TH1F *hE1D = (TH1F*) hEvsTime[i]->ProjectionY("_py",it,it);
      
      Int_t MAXCROSS = 2;
      Float_t *Cross = new Float_t[MAXCROSS];
      Float_t *Extr = new Float_t[MAXCROSS];
      memset(Cross,0,sizeof(Float_t)*MAXCROSS);
      memset(Extr,0,sizeof(Float_t)*MAXCROSS);

      Int_t auxNcross = PGlobals::HCrossings(hE1D,Cross,Extr,MAXCROSS,0.,0.);
      // cout << Form("  -> Number of crossings for histogram \"%s\" : %i ",hE1D->GetName(),auxNcross) << endl;
      // for(Int_t ic=0;ic<auxNcross;ic++) {
      // 	cout << Form(" %2i:  cross = %6.4f  extreme = %6.4f", ic, Cross[ic], Extr[ic]) << endl; 
      // }
      
      if(it==NTBins) {
	NCross[i] = auxNcross;
	
	gEcross[i] = new TGraph*[NCross[i]];
	gEextr[i] = new TGraph*[NCross[i]];

	for(Int_t ic = 0;ic<NCross[i];ic++) {
	  gEcross[i][ic] = new TGraph(NTBins);
	  sprintf(gName,"gEcross_%i_%i",i,ic); 
	  gEcross[i][ic]->SetName(gName);
	  
	  gEextr[i][ic] = new TGraph(NTBins);
	  sprintf(gName,"gEextr_%i_%i",i,ic); 
	  gEextr[i][ic]->SetName(gName);
	}

      }

      Float_t time = hEvsTime[i]->GetXaxis()->GetBinCenter(it);
      // cout << Form("Time step %i (%.2f): %i crossings",it,time,NCross[i]) << endl;
     
      for(Int_t ic=0;ic<NCross[i];ic++) {
	// cout << Form("  - Adding %i crossing: cross = %6.4f extreme = %6.4f",ic,Cross[ic],Extr[ic]) << endl;
	
	gEcross[i][ic]->SetPoint(it-1,time,Cross[ic]);
	gEextr[i][ic]->SetPoint(it-1,time,Extr[ic]);
      }
      
    }

    
    // Calculate the max and min of every crossing.
    // Also calculates dephasing:
    gEdephas[i] = new TGraph*[NCross[i]];
    for(Int_t ic = 0;ic<NCross[i];ic++) {

      Int_t Npoints = gEcross[i][ic]->GetN();
      Double_t *yEcross = gEcross[i][ic]->GetY();
      Double_t *yEextr  = gEextr[i][ic]->GetY();
      Double_t *xEextr  = gEextr[i][ic]->GetX();
      
      Double_t *yEdephas = new Double_t[Npoints];
      for(Int_t j=0;j<Npoints;j++) {
	yEdephas[j] = yEcross[j] - yEcross[0];
      }
      gEdephas[i][ic] = new TGraph(Npoints,xEextr,yEdephas);
      sprintf(gName,"gEdephas_%i_%i",i,ic); 
      gEdephas[i][ic]->SetName(gName);
      
      for(Int_t j=0;j<Npoints;j++) {
	if(yEcross[j]>maxEcross)
	  maxEcross = yEcross[j];
	if(yEcross[j]<minEcross)
	  minEcross = yEcross[j];
	
	if(yEextr[j]>maxEextr)
	  maxEextr = yEextr[j];
	if(yEextr[j]<minEextr)
	  minEextr = yEextr[j];
      
	// Only takes into account the minimums of the accelerating field:
	if(ic%2 || i!=0) continue;
	if(yEdephas[j]>maxEdephas)
	  maxEdephas = yEdephas[j];
	if(yEdephas[j]<minEdephas)
	  minEdephas = yEdephas[j];
	
      }
    
    }
  }

  // Transformer ratio vs time:
  // Take the ratio of the minimum over the maximum of the first oscillation of the Ez.
  {
    Int_t Npoints = gEextr[0][1]->GetN();
    Double_t *TR = new Double_t[Npoints];
    Double_t *yExtrPrev = gEextr[0][0]->GetY();
    Double_t *yExtr = gEextr[0][1]->GetY();
    Double_t *xExtr = gEextr[0][1]->GetX();
    for(Int_t j=0;j<Npoints;j++) {
      TR[j] = TMath::Abs(yExtr[j]/yExtrPrev[j]);	
    }
    
    gTRatio = new TGraph(Npoints,xExtr,TR);
    sprintf(gName,"gTRatio"); 
    gTRatio->SetName(gName);
  }
  

  sprintf(hName,"hVvsTime"); 
  hVvsTime = (TH2F*) ifile->Get(hName);
  Int_t NVCross = 0;
  cout << Form("ANALYZING POTENTIAL") << endl;

  Int_t NTBins = hVvsTime->GetNbinsX();
  for(Int_t it=NTBins;it>0;it--) {
      
    // 1D field at certain timestep "it".
    TH1F *hV1D = (TH1F*) hVvsTime->ProjectionY("_py",it,it);
      
    Int_t MAXCROSS = 2;
    Float_t *Cross = new Float_t[MAXCROSS];
    Float_t *Extr = new Float_t[MAXCROSS];
    memset(Cross,0,sizeof(Float_t)*MAXCROSS);
    memset(Extr,0,sizeof(Float_t)*MAXCROSS);
    
    Int_t auxNcross = PGlobals::HCrossings(hV1D,Cross,Extr,MAXCROSS,0.,20.);
    // cout << Form("  -> Number of crossings for histogram \"%s\" : %i ",hV1D->GetName(),auxNcross) << endl;
    // for(Int_t ic=0;ic<auxNcross;ic++) {
    // 	cout << Form(" %2i:  cross = %6.4f  extreme = %6.4f", ic, Cross[ic], Extr[ic]) << endl; 
    // }
    
    if(it==NTBins) {
      NVCross = auxNcross;
	
      gVcross = new TGraph*[NVCross];
      gVextr = new TGraph*[NVCross];

      for(Int_t ic = 0;ic<NVCross;ic++) {
	gVcross[ic] = new TGraph(NTBins);
	sprintf(gName,"gVcross_%i",ic); 
	gVcross[ic]->SetName(gName);
	  
	gVextr[ic] = new TGraph(NTBins);
	sprintf(gName,"gVextr_%i",ic); 
	gVextr[ic]->SetName(gName);
      }

    }

    Float_t time = hVvsTime->GetXaxis()->GetBinCenter(it);
    // cout << Form("Time step %i (%.2f): %i crossings",it,time,NVCross) << endl;
     
    for(Int_t ic=0;ic<NVCross;ic++) {
      // cout << Form("  - Adding %i crossing: cross = %6.4f extreme = %6.4f",ic,Cross[ic],Extr[ic]) << endl;
	
      gVcross[ic]->SetPoint(it-1,time,Cross[ic]);
      gVextr[ic]->SetPoint(it-1,time,Extr[ic]);
    }
      
  }


  sprintf(hName,"hETotalvsTime"); 
  hETvsTime = (TH2F*) ifile->Get(hName);
  Int_t NETCross = 0;
  cout << Form("ANALYZING TOTAL ELECTRIC FIELD") << endl;

  NTBins = hETvsTime->GetNbinsX();
  for(Int_t it=NTBins;it>0;it--) {
      
    // 1D field at certain timestep "it".
    TH1F *hET1D = (TH1F*) hETvsTime->ProjectionY("_py",it,it);
      
    Int_t MAXCROSS = 2;
    Float_t *Cross = new Float_t[MAXCROSS];
    Float_t *Extr = new Float_t[MAXCROSS];
    memset(Cross,0,sizeof(Float_t)*MAXCROSS);
    memset(Extr,0,sizeof(Float_t)*MAXCROSS);

    Float_t baseline = 0.5;
    if(opt.Contains("units"))
      baseline *= E0 / (PUnits::GV/PUnits::m);
    
    Int_t auxNcross = PGlobals::HCrossings(hET1D,Cross,Extr,MAXCROSS,baseline,-999,-999,"cum");
    // cout << Form("  -> Number of crossings for histogram \"%s\" : %i ",hET1D->GetName(),auxNcross) << endl;
    // for(Int_t ic=0;ic<auxNcross;ic++) {
    // 	cout << Form(" %2i:  cross = %6.4f  extreme = %6.4f", ic, Cross[ic], Extr[ic]) << endl; 
    // }
    
    if(it==NTBins) {
      NETCross = auxNcross;
	
      gETcross = new TGraph*[NETCross];
      gETextr = new TGraph*[NETCross];

      for(Int_t ic = 0;ic<NETCross;ic++) {
	gETcross[ic] = new TGraph(NTBins);
	sprintf(gName,"gETcross_%i",ic); 
	gETcross[ic]->SetName(gName);
	  
	gETextr[ic] = new TGraph(NTBins);
	sprintf(gName,"gETextr_%i",ic); 
	gETextr[ic]->SetName(gName);
      }

    }

    Float_t time = hETvsTime->GetXaxis()->GetBinCenter(it);
    // if(it==1)
    //   cout << Form("Time step %i (%.2f): %i crossings",it,time,NETCross) << endl;
     
    for(Int_t ic=0;ic<NETCross;ic++) {
      //if(it==1)
      //	cout << Form("  - Adding %i crossing: cross = %6.4f extreme = %6.4f",ic,Cross[ic],Extr[ic]) << endl;
	
      gETcross[ic]->SetPoint(it-1,time,Cross[ic]);
      gETextr[ic]->SetPoint(it-1,time,Extr[ic]);
      

    }
    
  }

  sprintf(hName,"hFocusvsTime"); 
  hFvsTime = (TH2F*) ifile->Get(hName);
  Int_t NFCross = 0;
  cout << Form("ANALYZING FOCUSING") << endl;

  NTBins = hFvsTime->GetNbinsX();
  for(Int_t it=NTBins;it>0;it--) {
      
    // 1D field at certain timestep "it".
    TH1F *hF1D = (TH1F*) hFvsTime->ProjectionY("_py",it,it);
      
    Int_t MAXCROSS = 2;
    Float_t *Cross = new Float_t[MAXCROSS];
    Float_t *Extr = new Float_t[MAXCROSS];
    memset(Cross,0,sizeof(Float_t)*MAXCROSS);
    memset(Extr,0,sizeof(Float_t)*MAXCROSS);

    Int_t auxNcross = PGlobals::HCrossings(hF1D,Cross,Extr,MAXCROSS,0.,0.);
    // cout << Form("  -> Number of crossings for histogram \"%s\" : %i ",hF1D->GetName(),auxNcross) << endl;
    // for(Int_t ic=0;ic<auxNcross;ic++) {
    // 	cout << Form(" %2i:  cross = %6.4f  extreme = %6.4f", ic, Cross[ic], Extr[ic]) << endl; 
    // }
    
    if(it==NTBins) {
      NFCross = auxNcross;
	
      gFcross = new TGraph*[NFCross];
      gFextr = new TGraph*[NFCross];

      for(Int_t ic = 0;ic<NFCross;ic++) {
	gFcross[ic] = new TGraph(NTBins);
	sprintf(gName,"gFcross_%i",ic); 
	gFcross[ic]->SetName(gName);
	  
	gFextr[ic] = new TGraph(NTBins);
	sprintf(gName,"gFextr_%i",ic); 
	gFextr[ic]->SetName(gName);
      }

    }

    Float_t time = hFvsTime->GetXaxis()->GetBinCenter(it);
    // cout << Form("Time step %i (%.2f): %i crossings",it,time,NFCross) << endl;
     
    for(Int_t ic=0;ic<NFCross;ic++) {
      // cout << Form("  - Adding %i crossing: cross = %6.4f extreme = %6.4f",ic,Cross[ic],Extr[ic]) << endl;
	
      gFcross[ic]->SetPoint(it-1,time,Cross[ic]);
      gFextr[ic]->SetPoint(it-1,time,Extr[ic]);
      

    }
    
  }
  

  for(Int_t i=0;i<NAtoms;i++) {
    sprintf(hName,"hIonProbvsTime_%s",atNames[i]); 
    hIonProbvsTime[i] = (TH2F*) ifile->Get(hName);
    if(!hIonProbvsTime[i]) continue;
    
    cout << Form("ANALYZING Ionization probability %i ...",i) << endl;

    Int_t NTBins = hIonProbvsTime[i]->GetNbinsX();
    for(Int_t it=NTBins;it>0;it--) {
      
      // 1D field at certain timestep "it".
      TH1F *hIonProb1D = (TH1F*) hIonProbvsTime[i]->ProjectionY("_py",it,it);
      
      Int_t MAXCROSS = 2;
      Float_t *Cross = new Float_t[MAXCROSS];
      Float_t *Extr = new Float_t[MAXCROSS];
      memset(Cross,0,sizeof(Float_t)*MAXCROSS);
      memset(Extr,0,sizeof(Float_t)*MAXCROSS);
      
      Int_t auxNcross = PGlobals::HCrossings(hIonProb1D,Cross,Extr,MAXCROSS,10.0);

      if(it==NTBins) {
	gIonProb10[i] = new TGraph(NTBins);
	sprintf(gName,"gIonProb10_%i",i); 
	gIonProb10[i]->SetName(gName);
      }

      Float_t time = hIonProbvsTime[i]->GetXaxis()->GetBinCenter(it);
      gIonProb10[i]->SetPoint(it-1,time,Cross[0]);
     
      memset(Cross,0,sizeof(Float_t)*MAXCROSS);
      memset(Extr,0,sizeof(Float_t)*MAXCROSS);
      auxNcross = PGlobals::HCrossings(hIonProb1D,Cross,Extr,MAXCROSS,99.0);
      
      if(it==NTBins) {
	gIonProb100[i] = new TGraph(NTBins);
	sprintf(gName,"gIonProb100_%i",i); 
	gIonProb100[i]->SetName(gName);
      }
      
      gIonProb100[i]->SetPoint(it-1,time,Cross[0]);
    }
  }

  // Set the color of the different evolutions according to a palette
  // UInt_t np = 50;
  // PPalette * colorPalette = (PPalette*) gROOT->FindObject("colorPalette");
  // if(!colorPalette) {
  //   const UInt_t Number = 3;
  //   Double_t Red[Number] = { 1.00, 0.00, 0.00};
  //   Double_t Green[Number]  = { 0.00, 1.00, 0.00};
  //   Double_t Blue[Number]   = { 1.00, 0.00, 1.00};
  //   Double_t Length[Number] = { 0.00, 0.50, 1.00 };
  //   colorPalette = new PPalette("colorPalette");
  //   colorPalette->CreateGradientColorTable(Number,Length,Red,Green,Blue,np);
  // }
  
  // for(Int_t i=0;i<Nfields;i++) { 
  //   for(Int_t ic=0;ic<Ncross;ic++) {
  //     Float_t step = (np/Nosc);
  //     Int_t icolor = TMath::Nint( ((ic+1)/2) * step - 1 );
  //     gEextr[i][ic]->SetLineColor(colorPalette->GetColor(icolor));
  //     gEextr[i][ic]->SetLineWidth(2);
  //     gEdephas[i][ic]->SetLineColor(colorPalette->GetColor(icolor));
  //     gEdephas[i][ic]->SetLineWidth(2);
  //   }
  // }
  // --------------------------------------------------------------------------

  // Manual coloring:
  const Int_t NCOLORS = 5;
  //  Int_t colors[NCOLORS] = {kMagenta+2,kRed,kBlue,kYellow+2,kCyan+2};
  Int_t colors[NCOLORS] = {kGray+3,kGray+2,kGray+1,kGray};
  for(Int_t i=0;i<Nfields;i++) { 
    for(Int_t ic=0;ic<NCross[i];ic++) {

      if( !gEcross[i][ic] || !gEextr[i][ic] ) continue;
            
      Int_t index = ic/2;
      if(index>=NCOLORS) index = NCOLORS-1;
      gEcross[i][ic]->SetLineColor(colors[index]);
      gEextr[i][ic]->SetLineColor(colors[index]);
      gEextr[i][ic]->SetLineWidth(1);
      gEdephas[i][ic]->SetLineColor(colors[index]);
      gEdephas[i][ic]->SetLineWidth(1);

      // cout << "EEEOOO" << endl;
      // if(ic%2) { 
      // 	gEcross[i][ic]->SetLineStyle(2);
      // 	gEextr[i][ic]->SetLineStyle(2);
      // 	gEdephas[i][ic]->SetLineStyle(2);
      // } else {
      // 	gEcross[i][ic]->SetLineStyle(1);
      // 	gEextr[i][ic]->SetLineStyle(1);
      // 	gEdephas[i][ic]->SetLineStyle(1);
      // }
      
      
    }
  }


  for(Int_t ic = 0;ic<NFCross;ic++) {
    // Graph's attributes
    Int_t index = ic/2;
    if(index>=NCOLORS) index = NCOLORS-1;
    gFcross[ic]->SetLineColor(colors[index]);
    if(ic%2-1) { 
      gFcross[ic]->SetLineStyle(2);
      gFextr[ic]->SetLineStyle(2);
      
    } else {
      gFcross[ic]->SetLineStyle(1);
      gFextr[ic]->SetLineStyle(1);
    }
  }

  for(Int_t ic = 0;ic<NVCross;ic++) {
    // Graph's attributes
    Int_t index = ic/2;
    if(index>=NCOLORS) index = NCOLORS-1;
    gVcross[ic]->SetLineColor(colors[index]);
    gVcross[ic]->SetLineStyle(3);
    gVextr[ic]->SetLineStyle(1);
  }
  
  for(Int_t ic = 0;ic<NETCross;ic++) {
    // Graph's attributes
    Int_t index = ic/2;
    if(index>=NCOLORS) index = NCOLORS-1;
    gETcross[ic]->SetLineColor(colors[index]);
    if(ic%2-1) { 
      gETcross[ic]->SetLineStyle(2);
      gETextr[ic]->SetLineStyle(2);
      
    } else {
      gETcross[ic]->SetLineStyle(1);
      gETextr[ic]->SetLineStyle(1);
    }
  }
  
  for(Int_t i=0;i<NAtoms;i++) {
    gIonProb10[i]->SetLineStyle(2);
    gIonProb10[i]->SetLineColor(kGray+2);

    gIonProb100[i]->SetLineStyle(1);
    gIonProb100[i]->SetLineColor(kGray+2);
  }
  

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
  Float_t mMargin = 0.015 * (950/ysize);
  Float_t pfactor = 1.0;
  if(opt.Contains("nomar"))
    bMargin = tMargin = lMargin = rMargin = mMargin = 0.0;
  if(NPad==1)
    PGlobals::CanvasPartition(C,NPad,lMargin,rMargin,bMargin,tMargin,mMargin);
  else
    PGlobals::CanvasAsymPartition(C,NPad,lMargin,rMargin,bMargin,tMargin,pfactor,mMargin);
 
  // Define the frames for plotting
  Int_t fonttype = 43;
  Int_t fontsize = 32;
  Int_t tfontsize = 38;
  Int_t txsize = tfontsize+6;
  Int_t lxsize = fontsize;
  Int_t tysize = tfontsize;
  Int_t lysize = fontsize-2;
  Int_t tzsize = tfontsize-4;
  Int_t lzsize = fontsize-2;
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
    pad[i]->SetTickx(1);
    pad[i]->SetTicky(1);
    if(opt.Contains("trans"))
      pad[i]->SetFillStyle(4000);
    pad[i]->SetFrameFillStyle(4000);

    sprintf(name,"hFrame_%i",i);
    hFrame[i] = (TH2F*) gROOT->FindObject(name);
    if(hFrame[i]) delete hFrame[i];
    hFrame[i] = (TH2F*) hEvsTime[0]->Clone(name);
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

  // Access to color Palettes 
  TExec *exPlasma = new TExec("exPlasma","plasmaPalette->cd();");
  TExec *exElec   = new TExec("exElec","electron0Palette->cd();");
  TExec *exHot    = new TExec("exHot","hotPalette->cd();");
  TExec *exField  = new TExec("exField","rbowwhitePalette->cd();");
  TExec *exFieldT = new TExec("exFieldT","red0Palette->cd();");
  TExec *exIonP   = new TExec("exIonP","redelectron0Palette->cd();");
  TExec *exPot    = new TExec("exPot","rbowinvPalette->cd();");

  

  // Actual Plotting!
  // ------------------------------------------------------------

  C->cd(0);

  Float_t x1,x2,y1,y2;
  Float_t gap = 0.01;
  TPaletteAxis *palette = NULL;
  UInt_t ip = NPad-1;

  if(mask & 0x1) {

    pad[ip]->Draw();
    pad[ip]->cd(); 
    if(opt.Contains("logz")) {
      pad[ip]->SetLogz(1);
    } else {
      pad[ip]->SetLogz(0);
    }

    hFrame[ip]->Draw("col");

    //  hDen1DvsTime->GetZaxis()->SetNdivisions(503);
    hDen1DvsTime->GetZaxis()->SetTitleFont(fonttype);
    Float_t xFactor = pad[0]->GetAbsWNDC()/pad[ip]->GetAbsWNDC();
    Float_t yFactor = pad[0]->GetAbsHNDC()/pad[ip]->GetAbsHNDC();
    hDen1DvsTime->GetZaxis()->SetTickLength(xFactor*tylength/yFactor);
    
    exElec->Draw();
    hDen1DvsTime->Draw("colz same");

    pad[ip]->Update();
  
    y1 = pad[ip]->GetBottomMargin();
    y2 = 1 - pad[ip]->GetTopMargin();
    x1 = pad[ip]->GetLeftMargin();
    x2 = 1 - pad[ip]->GetRightMargin();
  
    palette = (TPaletteAxis*) hDen1DvsTime->GetListOfFunctions()->FindObject("palette");  
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

  if(mask & 0x2) {

    pad[ip]->Draw();
    pad[ip]->cd(); 
    if(opt.Contains("logz")) {
      pad[ip]->SetLogz(1);
    } else {
      pad[ip]->SetLogz(0);
    }

    hFrame[ip]->Draw("col");

    //  hRmsvsTime->GetZaxis()->SetNdivisions(503);
    hRmsvsTime->GetZaxis()->SetTitleFont(fonttype);
    Float_t xFactor = pad[0]->GetAbsWNDC()/pad[ip]->GetAbsWNDC();
    Float_t yFactor = pad[0]->GetAbsHNDC()/pad[ip]->GetAbsHNDC();
    hRmsvsTime->GetZaxis()->SetTickLength(xFactor*tylength/yFactor);
    
    exElec->Draw();
    hRmsvsTime->Draw("colz same");

    pad[ip]->Update();
  
    y1 = pad[ip]->GetBottomMargin();
    y2 = 1 - pad[ip]->GetTopMargin();
    x1 = pad[ip]->GetLeftMargin();
    x2 = 1 - pad[ip]->GetRightMargin();
  
    palette = (TPaletteAxis*) hRmsvsTime->GetListOfFunctions()->FindObject("palette");  
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

  if(mask & 0x4) {

    pad[ip]->Draw();
    pad[ip]->cd(); 
    if(opt.Contains("logz")) {
      pad[ip]->SetLogz(1);
    } else {
      pad[ip]->SetLogz(0);
    }

    hFrame[ip]->Draw("col");

    //  hEvsTime[0]->GetZaxis()->SetNdivisions(503);
    hEvsTime[0]->GetZaxis()->SetTitleFont(fonttype);
    Float_t xFactor = pad[0]->GetAbsWNDC()/pad[ip]->GetAbsWNDC();
    Float_t yFactor = pad[0]->GetAbsHNDC()/pad[ip]->GetAbsHNDC();
    hEvsTime[0]->GetZaxis()->SetTickLength(xFactor*tylength/yFactor);
    
    exField->Draw();
    hEvsTime[0]->Draw("colz same");

    pad[ip]->Update();
  
    y1 = pad[ip]->GetBottomMargin();
    y2 = 1 - pad[ip]->GetTopMargin();
    x1 = pad[ip]->GetLeftMargin();
    x2 = 1 - pad[ip]->GetRightMargin();
  
    palette = (TPaletteAxis*) hEvsTime[0]->GetListOfFunctions()->FindObject("palette");  
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

    if(!opt.Contains("nocross")) {
    
      for(Int_t ic=0;ic<NVCross;ic++) {
	if( gVcross[ic] )
	  gVcross[ic]->Draw("L");
      }

      if(gEcross[1][0]) {
	gEcross[1][0]->SetLineStyle(4);
	gEcross[1][0]->Draw("L");
      }
    }

    pad[ip]->RedrawAxis(); 

    ip--;
    C->cd(0);
    
  }

  if(mask & 0x8) {

    pad[ip]->Draw();
    pad[ip]->cd(); 
    if(opt.Contains("logz")) {
      pad[ip]->SetLogz(1);
    } else {
      pad[ip]->SetLogz(0);
    }

    hFrame[ip]->Draw("col");

    //  hEvsTime[1]->GetZaxis()->SetNdivisions(503);
    hEvsTime[1]->GetZaxis()->SetTitleFont(fonttype);
    Float_t xFactor = pad[0]->GetAbsWNDC()/pad[ip]->GetAbsWNDC();
    Float_t yFactor = pad[0]->GetAbsHNDC()/pad[ip]->GetAbsHNDC();
    hEvsTime[1]->GetZaxis()->SetTickLength(xFactor*tylength/yFactor);
    
    exField->Draw();
    hEvsTime[1]->Draw("colz same");

    pad[ip]->Update();
  
    y1 = pad[ip]->GetBottomMargin();
    y2 = 1 - pad[ip]->GetTopMargin();
    x1 = pad[ip]->GetLeftMargin();
    x2 = 1 - pad[ip]->GetRightMargin();
  
    palette = (TPaletteAxis*) hEvsTime[1]->GetListOfFunctions()->FindObject("palette");  
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

  if(mask & 0x10) {

    pad[ip]->Draw();
    pad[ip]->cd(); 
    if(opt.Contains("logz")) {
      pad[ip]->SetLogz(1);
    } else {
      pad[ip]->SetLogz(0);
    }

    hFrame[ip]->Draw("col");

    //  hETvsTime->GetZaxis()->SetNdivisions(503);
    hETvsTime->GetZaxis()->SetTitleFont(fonttype);
    Float_t xFactor = pad[0]->GetAbsWNDC()/pad[ip]->GetAbsWNDC();
    Float_t yFactor = pad[0]->GetAbsHNDC()/pad[ip]->GetAbsHNDC();
    hETvsTime->GetZaxis()->SetTickLength(xFactor*tylength/yFactor);
    
    exFieldT->Draw();
    hETvsTime->Draw("colz same");

    pad[ip]->Update();
  
    y1 = pad[ip]->GetBottomMargin();
    y2 = 1 - pad[ip]->GetTopMargin();
    x1 = pad[ip]->GetLeftMargin();
    x2 = 1 - pad[ip]->GetRightMargin();
  
    palette = (TPaletteAxis*) hETvsTime->GetListOfFunctions()->FindObject("palette");  
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

    if(!opt.Contains("nocross")) {
    
      for(Int_t ic=0;ic<NVCross;ic++) {
	if( gVcross[ic] )
	  gVcross[ic]->Draw("L");
      }

      if(gEcross[1][0]) {
	gEcross[1][0]->SetLineStyle(4);
	gEcross[1][0]->Draw("L");
      }
    }
     
    pad[ip]->RedrawAxis(); 

    ip--;
    C->cd(0);
    
  }

  if(mask & 0x20) {

    pad[ip]->Draw();
    pad[ip]->cd(); 
    if(opt.Contains("logz")) {
      pad[ip]->SetLogz(1);
    } else {
      pad[ip]->SetLogz(0);
    }

    hFrame[ip]->Draw("col");

    //  hVvsTime->GetZaxis()->SetNdivisions(503);
    hVvsTime->GetZaxis()->SetTitleFont(fonttype);
    Float_t xFactor = pad[0]->GetAbsWNDC()/pad[ip]->GetAbsWNDC();
    Float_t yFactor = pad[0]->GetAbsHNDC()/pad[ip]->GetAbsHNDC();
    hVvsTime->GetZaxis()->SetTickLength(xFactor*tylength/yFactor);
    
    exPot->Draw();
    hVvsTime->Draw("colz same");

    pad[ip]->Update();
  
    y1 = pad[ip]->GetBottomMargin();
    y2 = 1 - pad[ip]->GetTopMargin();
    x1 = pad[ip]->GetLeftMargin();
    x2 = 1 - pad[ip]->GetRightMargin();
  
    palette = (TPaletteAxis*) hVvsTime->GetListOfFunctions()->FindObject("palette");  
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
    pad[ip]->cd(); 
    if(opt.Contains("logz")) {
      pad[ip]->SetLogz(1);
    } else {
      pad[ip]->SetLogz(0);
    }

    hFrame[ip]->Draw("col");

    hIonProbvsTime[1]->GetZaxis()->SetNdivisions(503);
    hIonProbvsTime[1]->GetZaxis()->SetTitleFont(fonttype);
    Float_t xFactor = pad[0]->GetAbsWNDC()/pad[ip]->GetAbsWNDC();
    Float_t yFactor = pad[0]->GetAbsHNDC()/pad[ip]->GetAbsHNDC();
    hIonProbvsTime[1]->GetZaxis()->SetTickLength(xFactor*tylength/yFactor);
    hIonProbvsTime[1]->GetZaxis()->SetRangeUser(0.,0.12);
    
    // exFieldT->Draw();
    // exPlasma->Draw();
    exIonP->Draw();
    hIonProbvsTime[1]->Draw("colz same");

    pad[ip]->Update();
  
    y1 = pad[ip]->GetBottomMargin();
    y2 = 1 - pad[ip]->GetTopMargin();
    x1 = pad[ip]->GetLeftMargin();
    x2 = 1 - pad[ip]->GetRightMargin();
  
    palette = (TPaletteAxis*) hIonProbvsTime[1]->GetListOfFunctions()->FindObject("palette");  
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

    if(!opt.Contains("nocross")) {
      
      for(Int_t ic=0;ic<NVCross;ic++) {
	if( gVcross[ic] )
	  gVcross[ic]->Draw("L");
      }
      
      if(gEcross[1][0]) {
	gEcross[1][0]->SetLineStyle(4);
	gEcross[1][0]->Draw("L");
      }
    }
     
    pad[ip]->RedrawAxis(); 

    ip--;
    C->cd(0);
    
  }
  
  // Print to a file
  // Output file
  TString fOutName = Form("./%s/Plots/Evolutions/Evolutions%s-WII-%s",sim.Data(),imask.c_str(),sim.Data());
  PGlobals::imgconv(C,fOutName,opt);
  // ---------------------------------------------------------

  
  ifile->Close();
  cout << endl;
}
  
