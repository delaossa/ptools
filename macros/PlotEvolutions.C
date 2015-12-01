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

void PlotEvolutions(const TString &sim, const TString &options="png") { 
  
#ifdef __CINT__  
  gSystem->Load("libptools.so");
#endif

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
  Float_t maxPsiextr = -999.;
  Float_t minPsiextr = 999.;

  
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
  TGraph **gVextr_alt = NULL;
  TGraph **gVextr_avg = NULL;

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
    
    Int_t auxNcross = PGlobals::HCrossings(hV1D,Cross,Extr,MAXCROSS,0.,0.);
    // cout << Form("  -> Number of crossings for histogram \"%s\" : %i ",hV1D->GetName(),auxNcross) << endl;
    // for(Int_t ic=0;ic<auxNcross;ic++) {
    // 	cout << Form(" %2i:  cross = %6.4f  extreme = %6.4f", ic, Cross[ic], Extr[ic]) << endl; 
    // }
    
    if(it==NTBins) {
      NVCross = auxNcross;
	
      gVcross = new TGraph*[NVCross];
      gVextr = new TGraph*[NVCross];
      gVextr_alt = new TGraph*[NVCross];
      gVextr_avg = new TGraph*[NVCross];
      
      for(Int_t ic = 0;ic<NVCross;ic++) {
	gVcross[ic] = new TGraph(NTBins);
	sprintf(gName,"gVcross_%i",ic); 
	gVcross[ic]->SetName(gName);
	  
	gVextr[ic] = new TGraph(NTBins);
	sprintf(gName,"gVextr_%i",ic); 
	gVextr[ic]->SetName(gName);

	gVextr_alt[ic] = new TGraph(NTBins);
	sprintf(gName,"gVextr_alt_%i",ic); 
	gVextr_alt[ic]->SetName(gName);

	gVextr_avg[ic] = new TGraph(NTBins);
	sprintf(gName,"gVextr_avg_%i",ic); 
	gVextr_avg[ic]->SetName(gName);
      }

    }

    Float_t time = hVvsTime->GetXaxis()->GetBinCenter(it);
    // cout << Form("Time step %i (%.2f): %i crossings",it,time,NVCross) << endl;
     
    for(Int_t ic=0;ic<NVCross;ic++) {
      // cout << Form("  - Adding %i crossing: cross = %6.4f extreme = %6.4f",ic,Cross[ic],Extr[ic]) << endl;
	
      gVcross[ic]->SetPoint(it-1,time,Cross[ic]);
      gVextr[ic]->SetPoint(it-1,time,Extr[ic]);

      if(Extr[ic]>maxPsiextr)
	maxPsiextr = Extr[ic];
      if(Extr[ic]<minPsiextr)
	minPsiextr = Extr[ic];

      if(gEcross[0][ic]) {
	Double_t x,y;
	gEcross[0][ic]->GetPoint(it-1,x,y);      
	Double_t value = hV1D->GetBinContent(hV1D->FindBin(y));
	gVextr_alt[ic]->SetPoint(it-1,time,value);
	gVextr_avg[ic]->SetPoint(it-1,time,(value+Extr[ic])/2);
      }
            
    }
      
  }


  sprintf(hName,"hETotalvsTime"); 
  hETvsTime = (TH2F*) ifile->Get(hName);
  Int_t NETCross = 0;
  if(hETvsTime) {
    
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
    hIonProbvsTime[i] = NULL;
    gIonProb10[i] = NULL;
    gIonProb100[i] = NULL;
    
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
    gVextr[ic]->SetLineWidth(2);
    gVextr_alt[ic]->SetLineStyle(1);
    gVextr_alt[ic]->SetLineColor(2);
    gVextr_alt[ic]->SetLineWidth(2);
    gVextr_avg[ic]->SetLineStyle(1);
    gVextr_avg[ic]->SetLineColor(1);
    gVextr_avg[ic]->SetLineWidth(1);
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
    if(gIonProb10[i]) {
      gIonProb10[i]->SetLineStyle(2);
      gIonProb10[i]->SetLineColor(kGray+2);
    }
    
    if(gIonProb100[i]) {
      gIonProb100[i]->SetLineStyle(1);
      gIonProb100[i]->SetLineColor(kGray+2);
    }
  }
  

  // palettes for drawing
  PPalette * redPalette = (PPalette*) gROOT->FindObject("red0");
  PPalette * rbowwhitePalette = (PPalette*) gROOT->FindObject("rbowwhite");
  PPalette * electronPalette = (PPalette*) gROOT->FindObject("oli");
  PPalette * electroninvPalette = (PPalette*) gROOT->FindObject("oli");
  PPalette * barsaPalette = (PPalette*) gROOT->FindObject("barsa");
  
  // Canvas setup
  Int_t sizex = 1024;
  Int_t sizey = 640;
  TCanvas *C = new TCanvas("C","Evolution of Electric fields",sizex,sizey);
  
  C->cd(0);
  gPad->SetFrameLineWidth(2);  
  gPad->SetTickx(1);
  gPad->SetTicky(1);
  // gPad->SetFillStyle(4000);
  // gPad->SetFrameFillStyle(4000);

  PGlobals::SetH1LabelSize(hEvsTime[0]);
 
  // Change the range of z axis for the fields to be symmetric.
  Float_t Emax = hEvsTime[0]->GetMaximum();
  Float_t Emin = hEvsTime[0]->GetMinimum();
  if(Emax > TMath::Abs(Emin))
    Emin = -Emax;
  else
    Emax = -Emin;
  hEvsTime[0]->GetZaxis()->SetRangeUser(Emin,Emax); 
 
  hEvsTime[0]->GetYaxis()->SetLabelFont(43);
  hEvsTime[0]->GetYaxis()->SetLabelSize(32);
  //  hEvsTime[0]->GetYaxis()->SetLabelOffset(0.8);
  hEvsTime[0]->GetYaxis()->SetTitleFont(43);
  hEvsTime[0]->GetYaxis()->SetTitleSize(36);
  hEvsTime[0]->GetYaxis()->SetTitleOffset(0.9);

  hEvsTime[0]->GetXaxis()->SetLabelFont(43);
  hEvsTime[0]->GetXaxis()->SetLabelSize(32);
  //  hEvsTime[0]->GetXaxis()->SetLabelOffset(0.8);
  hEvsTime[0]->GetXaxis()->SetTitleFont(43);
  hEvsTime[0]->GetXaxis()->SetTitleSize(36);
  hEvsTime[0]->GetXaxis()->SetTitleOffset(0.95);


  hEvsTime[0]->GetZaxis()->SetLabelFont(43);
  hEvsTime[0]->GetZaxis()->SetLabelSize(32);
  // hEvsTime[0]->GetZaxis()->SetLabelOffset(0.8);
  hEvsTime[0]->GetZaxis()->SetTitleFont(43);
  hEvsTime[0]->GetZaxis()->SetTitleSize(36);
  hEvsTime[0]->GetZaxis()->SetTitleOffset(0.9);

  rbowwhitePalette->cd();
  
  hEvsTime[0]->Draw("colz");

  gPad->Update();
  TPaletteAxis *palette = (TPaletteAxis*)hEvsTime[0]->GetListOfFunctions()->FindObject("palette");
  palette->SetX1NDC(1-margins[2]+0.02);
  palette->SetX2NDC(1-margins[2]+0.05);
  palette->SetY1NDC(margins[1]);
  palette->SetY2NDC(1-margins[3]);
  C->Modified();
  C->Update();

  // Plot the crossings
  if(!opt.Contains("nocross")) {
    for(Int_t ic=0;ic<NCross[0];ic++) {
      if( gEcross[0][ic] )
	gEcross[0][ic]->Draw("L");
    }
  }
  
  
  gPad->RedrawAxis("G");

  // Print to a file
  // Output file
  TString fOutName = Form("./%s/Plots/Evolutions/Evolutions-Ez-%s",sim.Data(),sim.Data());
  PGlobals::imgconv(C,fOutName,opt);
  // ---------------------------------------------------------

  // Canvas setup
  TCanvas *C0 = new TCanvas("C0","Evolution of Electric fields",sizex,sizey);
  
  C0->cd(0);
  gPad->SetFrameLineWidth(2);  
  gPad->SetTickx(1);
  gPad->SetTicky(1);

  PGlobals::SetH1LabelSize(hEvsTime[1]);
 
  Emax = hEvsTime[1]->GetMaximum();
  Emin = hEvsTime[1]->GetMinimum();
  hEvsTime[1]->GetZaxis()->SetRangeUser(Emin,Emax); 
 
  // Dynamic Ex palette
  const Int_t exPNRGBs = 6;
  const Int_t exPNCont = 64;
  Float_t zeroPos = -Emin/(Emax-Emin);
  
  Double_t exPStops[exPNRGBs] = { 0.00, zeroPos-3.0/exPNCont,zeroPos-1.0/exPNCont, zeroPos, zeroPos+1.0/exPNCont, 1.00 };
  Double_t exPRed[exPNRGBs]   = { 0.106, 0.698, 0.90, 0.90, 0.965, 0.518 };
  Double_t exPGreen[exPNRGBs] = { 0.078, 0.818, 0.90, 0.90, 0.925, 0.078 };
  Double_t exPBlue[exPNRGBs]  = { 0.518, 0.880, 0.90, 0.90, 0.353, 0.106 };
   
  PPalette * exPalette = (PPalette*) gROOT->FindObject("rbowwhite");
  exPalette->CreateGradientColorTable(exPNRGBs, exPStops, 
				      exPRed, exPGreen, exPBlue, exPNCont);

  
  hEvsTime[1]->Draw("colz");

  gPad->Update();
  palette = (TPaletteAxis*)hEvsTime[1]->GetListOfFunctions()->FindObject("palette");
  palette->SetX1NDC(1-margins[2]+0.02);
  palette->SetX2NDC(1-margins[2]+0.05);
  palette->SetY1NDC(margins[1]);
  palette->SetY2NDC(1-margins[3]);
  C0->Modified();
  C0->Update();

  // Plot the crossings
  if(!opt.Contains("nocross")) {
    for(Int_t ic=0;ic<NCross[1];ic++) {
      if( gEcross[1][ic] )
	gEcross[1][ic]->Draw("L");
    }
  }
    
  gPad->RedrawAxis("G");

  // Print to a file
  // Output file
  fOutName = Form("./%s/Plots/Evolutions/Evolutions-Ex-%s",sim.Data(),sim.Data());
  PGlobals::imgconv(C0,fOutName,opt);
  // ---------------------------------------------------------

  if(hEvsTime[0]) {
    
    // Canvas setup
    TCanvas *C1 = new TCanvas("C1","Evolution of Electric fields",sizex,sizey);
    gPad->SetTickx(0);
    gPad->SetTicky(0);


    TLegend *Leg = new TLegend(0.20,0.20,0.35,0.40);
    PGlobals::SetPaveStyle(Leg);
    Leg->SetTextAlign(22);
    Leg->SetTextColor(kGray+3);
    Leg->SetLineColor(kGray+1);
    Leg->SetBorderSize(0);
    Leg->SetFillColor(0);
    Leg->SetFillStyle(1001);
    for(Int_t ic=0;ic<NCross[0];ic++) {
      // Only minimums:
      // if(ic%2-1) continue;
      Leg->AddEntry(gEextr[0][ic],Form("Minimum %2i",ic/2+1),"L");
    }
    Leg->SetTextColor(kGray+3);
  
  
    C1->cd(0);
  
    Float_t margin = (maxEextr - minEextr)/10;
    Float_t Ez0   = minEextr - margin;
    Float_t Ez1   = maxEextr + margin;
    TH1F *hFrame2 = new TH1F("hFrame2","",100,hEvsTime[0]->GetXaxis()->GetXmin(),hEvsTime[0]->GetXaxis()->GetXmax());
    // TH2F *hFrame2 = (TH2F*) hEvsTime[0]->Clone("hFrame2");
  
    hFrame2->GetYaxis()->SetRangeUser(Ez0,Ez1);
    hFrame2->GetXaxis()->SetTitle(hEvsTime[0]->GetXaxis()->GetTitle());
    hFrame2->GetYaxis()->SetTitle("E_{z}/E_{0}");
    PGlobals::SetH1LabelSize(hFrame2);
    
    hFrame2->Draw("axis");

    gPad->Update();

    if(!opt.Contains("nocross")) {

      for(Int_t ic=0;ic<NCross[0];ic++) {
      
	// Only minimums:
	// if(ic%2-1) continue;
	if( !gEextr[0][ic] ) continue;

	// Remove crazy points
	Double_t *yEextr  = gEextr[0][ic]->GetY();

	for(Int_t k=0; k < gEextr[0][ic]->GetN(); k++) {
	  //cout << "PEPEEEEEEEEEEE" << endl;
	  if(fabs(yEextr[k])>1000) {
	    gEextr[0][ic]->RemovePoint(k);
	    k--;
	  }
	}
      
	gEextr[0][ic]->Draw("L");
      
      }
    }


    // Transformer ratio
    Double_t rightmin = 0.0;
    Double_t rightmax = 12.0;
    Double_t slope = (gPad->GetUymax() - gPad->GetUymin())/(rightmax-rightmin); 
  
    Double_t *xTR = gTRatio->GetX();
    Double_t *yTR = gTRatio->GetY();
    for(Int_t j=0;j<gTRatio->GetN();j++) {
      gTRatio->SetPoint(j,xTR[j],yTR[j]*slope + gPad->GetUymin());
    }
  
    gTRatio->SetLineWidth(2);
    gTRatio->SetLineColor(kRed);
    gTRatio->Draw("L");
  
    //draw an axis on the right side
    TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),gPad->GetUxmax(),
			      gPad->GetUymax(),rightmin,rightmax,505,"+L");
  
    axis->SetLineWidth(1);
    axis->SetLineColor(gTRatio->GetLineColor());
    axis->SetLabelColor(gTRatio->GetLineColor());
    // axis->SetLabelSize(0.05);
    // axis->SetTitleSize(0.06);
    // axis->SetTitleOffset(0.7);
    axis->SetTitle("|E^{-} / E^{+}|");
    axis->CenterTitle();
    axis->SetTitleColor(gTRatio->GetLineColor());
  
    axis->Draw();


    C1->cd(0);
  
    // Print to a file
    // Output file
    fOutName = Form("./%s/Plots/Evolutions/Evolutions-EzMax-%s",sim.Data(),sim.Data());
    PGlobals::imgconv(C1,fOutName,opt);
  }

  if(hVvsTime) {
    // Canvas setup
    TCanvas *CV = new TCanvas("CV","Evolution of the trapping potential",sizex,sizey);
  
    CV->cd(0);
    gPad->SetFrameLineWidth(2);  
    gPad->SetTickx(1);
    gPad->SetTicky(1);

    PGlobals::SetH1LabelSize(hVvsTime);

    // Change the range of z axis for the fields to be symmetric.
    Float_t Vmax = hVvsTime->GetMaximum();
    Float_t Vmin = hVvsTime->GetMinimum();
    if(Vmax > TMath::Abs(Vmin))
      Vmin = -Vmax;
    else
      Vmax = -Vmin;
    hVvsTime->GetZaxis()->SetRangeUser(Emin,Emax); 
 
    hVvsTime->GetYaxis()->SetLabelFont(43);
    hVvsTime->GetYaxis()->SetLabelSize(32);
    //  hVvsTime->GetYaxis()->SetLabelOffset(0.8);
    hVvsTime->GetYaxis()->SetTitleFont(43);
    hVvsTime->GetYaxis()->SetTitleSize(36);
    hVvsTime->GetYaxis()->SetTitleOffset(0.9);

    hVvsTime->GetXaxis()->SetLabelFont(43);
    hVvsTime->GetXaxis()->SetLabelSize(32);
    //  hVvsTime->GetXaxis()->SetLabelOffset(0.8);
    hVvsTime->GetXaxis()->SetTitleFont(43);
    hVvsTime->GetXaxis()->SetTitleSize(36);
    hVvsTime->GetXaxis()->SetTitleOffset(0.95);


    hVvsTime->GetZaxis()->SetLabelFont(43);
    hVvsTime->GetZaxis()->SetLabelSize(32);
    // hVvsTime->GetZaxis()->SetLabelOffset(0.8);
    hVvsTime->GetZaxis()->SetTitleFont(43);
    hVvsTime->GetZaxis()->SetTitleSize(36);
    hVvsTime->GetZaxis()->SetTitleOffset(0.9);

    rbowwhitePalette->cd();
  
    // Float_t Vmax = hVvsTime->GetMaximum();
    // Float_t Vmin = hVvsTime->GetMinimum();
    // hVvsTime->GetZaxis()->SetRangeUser(Vmin,Vmax); 
    // if(Vmax<0.1) Vmax = 0.1;
  
    // // Dynamic potential palette
    // const Int_t potPNRGBs = 6;
    // const Int_t potPNCont = 64;
    // zeroPos = -Vmin/(Vmax-Vmin);
  
    // Double_t potPStops[potPNRGBs] = { 0.00, zeroPos-3.0/potPNCont,zeroPos-1.0/potPNCont, zeroPos, zeroPos+3.0/potPNCont, 1.00 };
    // Double_t potPRed[potPNRGBs]   = { 0.518, 0.965, 0.90,0.90, 0.498, 0.106 };
    // Double_t potPGreen[potPNRGBs] = { 0.078, 0.925, 0.90,0.90, 0.718, 0.078 };
    // Double_t potPBlue[potPNRGBs]  = { 0.106, 0.353, 0.90,0.90, 0.780, 0.518 };
   
    // PPalette * potentialPalette = (PPalette*) gROOT->FindObject("rbowinv");
    // potentialPalette->CreateGradientColorTable(potPNRGBs, potPStops, 
    // 					     potPRed, potPGreen, potPBlue, potPNCont);
  
    hVvsTime->Draw("colz");

    gPad->Update();
    palette = (TPaletteAxis*)hVvsTime->GetListOfFunctions()->FindObject("palette");
    palette->SetX1NDC(1-margins[2]+0.02);
    palette->SetX2NDC(1-margins[2]+0.05);
    palette->SetY1NDC(margins[1]);
    palette->SetY2NDC(1-margins[3]);
    CV->Modified();
    CV->Update();

    // Plot the crossings
    if(!opt.Contains("nocross")) {
    
      for(Int_t ic=0;ic<NVCross;ic++) {
	if( gVcross[ic] )
	  gVcross[ic]->Draw("L");
      }
    
      for(Int_t ic=0;ic<NCross[0];ic++) {
	if( gEcross[0][ic] ) {
	  gEcross[0][ic]->SetLineColor(kWhite);
	  gEcross[0][ic]->Draw("L");
	
	}
      
      }
      
    
    
    }
  
  
    gPad->RedrawAxis("G");

    // Print to a file
    // Output file
    fOutName = Form("./%s/Plots/Evolutions/Evolutions-Psi-%s",sim.Data(),sim.Data());
    PGlobals::imgconv(CV,fOutName,opt);
    // ---------------------------------------------------------

    
    // Canvas setup
    TCanvas *CV1 = new TCanvas("CV1","Evolution of potential extremes",sizex,sizey);
    gPad->SetTickx(0);
    gPad->SetTicky(0);


    TLegend *Leg = new TLegend(0.20,0.20,0.35,0.40);
    PGlobals::SetPaveStyle(Leg);
    Leg->SetTextAlign(22);
    Leg->SetTextColor(kGray+3);
    Leg->SetLineColor(kGray+1);
    Leg->SetBorderSize(0);
    Leg->SetFillColor(0);
    Leg->SetFillStyle(1001);
    for(Int_t ic=0;ic<NVCross;ic++) {
      // Only minimums:
      // if(ic%2-1) continue;
      Leg->AddEntry(gVextr[ic],Form("Minimum %2i",ic/2+1),"L");
    }
    Leg->SetTextColor(kGray+3);
  
  
    CV1->cd(0);
  
    Float_t margin = (maxPsiextr - minPsiextr)/10;
    Float_t Psi0   = minPsiextr - margin;
    Float_t Psi1   = maxPsiextr + margin;
    TH1F *hFrame2 = new TH1F("hFrame2","",100,hVvsTime->GetXaxis()->GetXmin(),hVvsTime->GetXaxis()->GetXmax());
    // TH2F *hFrame2 = (TH2F*) hVvsTime->Clone("hFrame2");
  
    hFrame2->GetYaxis()->SetRangeUser(Psi0,Psi1);
    hFrame2->GetXaxis()->SetTitle(hVvsTime->GetXaxis()->GetTitle());
    hFrame2->GetYaxis()->SetTitle("#psi");
    PGlobals::SetH1LabelSize(hFrame2);
    
    hFrame2->Draw("axis");

    gPad->Update();

    if(!opt.Contains("nocross")) {

      for(Int_t ic=0;ic<NVCross;ic++) {
      
	// Only minimums:
	// if(ic%2-1) continue;
	if( !gVextr[ic] ) continue;

	// Remove crazy points
	Double_t *yVextr  = gVextr[ic]->GetY();

	for(Int_t k=0; k < gVextr[ic]->GetN(); k++) {
	  if(fabs(yVextr[k])>1000) {
	    gVextr[ic]->RemovePoint(k);
	    k--;
	  }
	}
      
	gVextr[ic]->Draw("L");
	gVextr_alt[ic]->Draw("L");
	gVextr_avg[ic]->Draw("L");
      
      }
    }


    CV1->cd(0);
  
    // Print to a file
    // Output file
    fOutName = Form("./%s/Plots/Evolutions/Evolutions-PsiMax-%s",sim.Data(),sim.Data());
    PGlobals::imgconv(CV1,fOutName,opt);
  


    // Special plot
    
    // Canvas setup
    TCanvas *CVDen = new TCanvas("CVDen","Evolution of potential extremes",sizex,sizey);
    gPad->SetTickx(0);
    gPad->SetTicky(0);
  
    CVDen->cd(0);

    Float_t ntop = 15.63;
    Float_t sigmal = 37.64;

    TH1F *hFrame3 = new TH1F("hFrame3","",100,0.,16.);
  
    hFrame3->GetYaxis()->SetRangeUser(Psi0,Psi1);
    hFrame3->GetXaxis()->SetTitle("n/n_{0}");
    hFrame3->GetYaxis()->SetTitle("#psi");
    PGlobals::SetH1LabelSize(hFrame3);
    
    hFrame3->Draw("axis");

    gPad->Update();

    Int_t Np = gVextr_avg[1]->GetN();
    Double_t *yVextr  = gVextr_avg[1]->GetY();
    Double_t *xVextr  = gVextr_avg[1]->GetX();
    Double_t *den = new Double_t[Np];

    for(Int_t ip=0;ip<Np;ip++) {
      
      den[ip] = 1 + (ntop-1) * TMath::Exp(-(xVextr[ip]*xVextr[ip])/(2*sigmal*sigmal));
    }
    
    TGraph *gVvsDen = new TGraph(Np,den,yVextr);
    gVvsDen->SetLineStyle(1);
    gVvsDen->SetLineColor(kBlue);
    gVvsDen->SetLineWidth(2);

    gVvsDen->Draw("L");
    
    CVDen->cd(0);
  
    // Print to a file
    // Output file
    fOutName = Form("./%s/Plots/Evolutions/Evolutions-PsivsDen-%s",sim.Data(),sim.Data());
    PGlobals::imgconv(CVDen,fOutName,opt);
    
  }
  
  
  if(hFvsTime) {
    // Canvas setup
    TCanvas *CF = new TCanvas("CF","Evolution of the focusing field",sizex,sizey);
 
    CF->cd(0);
    gPad->SetFrameLineWidth(2);  
    gPad->SetTickx(1);
    gPad->SetTicky(1);
 
    PGlobals::SetH1LabelSize(hFvsTime);
 
    Float_t Fmax = hFvsTime->GetMaximum();
    Float_t Fmin = hFvsTime->GetMinimum();
    hFvsTime->GetZaxis()->SetRangeUser(Fmin,Fmax); 
    //if(Fmax<0.1) Fmax = 0.1;
  
    // Dynamic focusing palette
    const Int_t focPNRGBs = 6;
    const Int_t focPNCont = 64;
    zeroPos = -Fmin/(Fmax-Fmin);
  
    Double_t focPStops[focPNRGBs] = { 0.00, zeroPos-3.0/focPNCont,zeroPos-1.0/focPNCont, zeroPos, zeroPos+3.0/focPNCont, 1.00 };
    Double_t focPRed[focPNRGBs]   = { 0.106, 0.698, 0.90, 0.90, 0.965, 0.518 };
    Double_t focPGreen[focPNRGBs] = { 0.078, 0.818, 0.90, 0.90, 0.925, 0.078 };
    Double_t focPBlue[focPNRGBs]  = { 0.518, 0.880, 0.90, 0.90, 0.353, 0.106 };
   
    PPalette * focusPalette = (PPalette*) gROOT->FindObject("rbowwhite");
    focusPalette->CreateGradientColorTable(focPNRGBs, focPStops, 
					   focPRed, focPGreen, focPBlue, focPNCont);
  
    hFvsTime->Draw("colz");

    gPad->Update();
    palette = (TPaletteAxis*)hFvsTime->GetListOfFunctions()->FindObject("palette");
    palette->SetX1NDC(1-margins[2]+0.02);
    palette->SetX2NDC(1-margins[2]+0.05);
    palette->SetY1NDC(margins[1]);
    palette->SetY2NDC(1-margins[3]);
    CF->Modified();
    CF->Update();

    // Plot the crossings
    if(!opt.Contains("nocross")) {

      for(Int_t ic=0;ic<NFCross;ic++) {
	if( gFcross[ic] )
	  gFcross[ic]->Draw("L");
      }
    }
  
    gPad->RedrawAxis("G");

    // Print to a file
    // Output file
    fOutName = Form("./%s/Plots/Evolutions/Evolutions-Focus-%s",sim.Data(),sim.Data());
    PGlobals::imgconv(CF,fOutName,opt);
    // ---------------------------------------------------------
  }

  if(hETvsTime) {
    
    // Canvas setup
    TCanvas *CT = new TCanvas("CT","Evolution of the total electric field",sizex,sizey);
    
    CT->cd(0);
    gPad->SetFrameLineWidth(2);  
    gPad->SetTickx(1);
    gPad->SetTicky(1);
    // gPad->SetFillStyle(4000);
    // gPad->SetFrameFillStyle(4000);
    
    PGlobals::SetH1LabelSize(hETvsTime);
    
    // Change the range of z axis for the fields to be symmetric.
    Float_t ETmax = hETvsTime->GetMaximum();
    Float_t ETmin = hETvsTime->GetMinimum();
    hETvsTime->GetZaxis()->SetRangeUser(ETmin,ETmax); 
    
    hETvsTime->GetYaxis()->SetLabelFont(43);
    hETvsTime->GetYaxis()->SetLabelSize(32);
    //  hETvsTime->GetYaxis()->SetLabelOffset(0.8);
    hETvsTime->GetYaxis()->SetTitleFont(43);
    hETvsTime->GetYaxis()->SetTitleSize(36);
    hETvsTime->GetYaxis()->SetTitleOffset(0.9);
    
    hETvsTime->GetXaxis()->SetLabelFont(43);
    hETvsTime->GetXaxis()->SetLabelSize(32);
    //  hETvsTime->GetXaxis()->SetLabelOffset(0.8);
    hETvsTime->GetXaxis()->SetTitleFont(43);
    hETvsTime->GetXaxis()->SetTitleSize(36);
    hETvsTime->GetXaxis()->SetTitleOffset(0.95);
    
    
    hETvsTime->GetZaxis()->SetLabelFont(43);
    hETvsTime->GetZaxis()->SetLabelSize(32);
    // hETvsTime->GetZaxis()->SetLabelOffset(0.8);
    hETvsTime->GetZaxis()->SetTitleFont(43);
    hETvsTime->GetZaxis()->SetTitleSize(36);
    hETvsTime->GetZaxis()->SetTitleOffset(0.9);

    
    redPalette->cd();
    hETvsTime->Draw("colz");

    gPad->Update();
    palette = (TPaletteAxis*)hETvsTime->GetListOfFunctions()->FindObject("palette");
    palette->SetX1NDC(1-margins[2]+0.02);
    palette->SetX2NDC(1-margins[2]+0.05);
    palette->SetY1NDC(margins[1]);
    palette->SetY2NDC(1-margins[3]);
    CT->Modified();
    CT->Update();
    
    gPad->RedrawAxis("G");
    
    if(!opt.Contains("nocross")) {
      
      // if(NETCross>0)
      //   gETcross[0]->Draw("L");
      
      // if(NETCross>2)
      //   gETcross[2]->Draw("L");
    
      for(Int_t ic=0;ic<NVCross;ic++) {
	if( gVcross[ic] )
	  gVcross[ic]->Draw("L");
      }
    }
    
    // Print to a file
    // Output file
    fOutName = Form("./%s/Plots/Evolutions/Evolutions-ETotal-%s",sim.Data(),sim.Data());
    PGlobals::imgconv(CT,fOutName,opt);
  }
  // ---------------------------------------------------------

    // Canvas setup
  TCanvas *C2 = new TCanvas("C2","Evolution of the dephasing",sizex,sizey);
  C2->cd(0);
 
  gPad->SetTickx(0);
  gPad->SetTicky(0);

  Float_t margin = (maxEdephas - minEdephas)/10;
  Float_t D0   = minEdephas-margin;
  Float_t D1   = maxEdephas+margin;
  TH1F *hFrame4 = new TH1F("hFrame4","",100,hEvsTime[0]->GetXaxis()->GetXmin(),hEvsTime[0]->GetXaxis()->GetXmax());
  
  hFrame4->GetYaxis()->SetRangeUser(D0,D1);
  hFrame4->GetXaxis()->SetTitle(hEvsTime[0]->GetXaxis()->GetTitle());
  hFrame4->GetYaxis()->SetTitle("#Delta#zeta");

  PGlobals::SetH1LabelSize(hFrame4);
  hFrame4->Draw("axis");
  
  gPad->Update();
     
  for(Int_t i=0;i<Nfields;i++) {
    if(i>0) continue;
    
    for(Int_t ic=0;ic<NCross[i];ic++) {
      
      // Only minimums:
      if(ic%2) continue;
      if( !gEdephas[i][ic] ) continue;
      
      gEdephas[i][ic]->Draw("L");
    }
    
    // Leg->Draw();
  }
  
  gPad->RedrawAxis("G");
  
  C2->cd(0);
  
  // Print to a file
  // Output file
  fOutName = Form("./%s/Plots/Evolutions/Evolutions-EzDephas-%s",sim.Data(),sim.Data());
  PGlobals::imgconv(C2,fOutName,opt);
  // ---------------------------------------------------------


  // TH2F *hGvsTime = NULL;
  // sprintf(hName,"hGvsTime"); 
  // hGvsTime = (TH2F*) ifile->Get(hName);
  
  // if(hGvsTime) {
  //   // Canvas setup
  //   TCanvas *C3 = new TCanvas("C3","Evolution of the driver energy",sizex,sizey);
  //   C3->cd(0);
  //   gPad->SetFrameLineWidth(2);  
  //   gPad->SetTickx(1);
  //   gPad->SetTicky(1);

  //   Float_t minGamma =  43.07 - 0.5;
  //   Float_t maxGamma =  43.07 + 0.5;
  
  //   hGvsTime->GetZaxis()->SetRangeUser(minGamma,maxGamma); 
  //   PGlobals::SetH1LabelSize(hGvsTime);
  
  //   // barsaPalette->cd();
  //   rbowwhitePalette->cd();
    
  //   hGvsTime->Draw("colz");

  //   gPad->Update();
  
  //   palette = (TPaletteAxis*)hGvsTime->GetListOfFunctions()->FindObject("palette");
  //   palette->SetX1NDC(1-margins[2]+0.02);
  //   palette->SetX2NDC(1-margins[2]+0.05);
  //   palette->SetY1NDC(margins[1]);
  //   palette->SetY2NDC(1-margins[3]);
  //   C3->Modified();  
  //   gPad->RedrawAxis("G");
  
  //   C3->cd(0);
  
  //   // Print to a file
  //   // Output file
  //   fOutName = Form("./%s/Plots/Evolutions/Evolutions-GammaDriver-%s",sim.Data(),sim.Data());
  //   PGlobals::imgconv(C3,fOutName,opt);
  //   // ---------------------------------------------------------
  // }

  TH2F *hRmsvsTime = NULL;
  sprintf(hName,"hRmsvsTime_1"); 
  hRmsvsTime = (TH2F*) ifile->Get(hName);
  if(hRmsvsTime) {
    
    // Canvas setup
    TCanvas *C4 = new TCanvas("C4","Evolution of the transverse RMS",sizex,sizey);
    C4->cd(0);
    gPad->SetFrameLineWidth(2);  
    gPad->SetTickx(1);
    gPad->SetTicky(1);
  
    Float_t Rmsmax = hRmsvsTime->GetMaximum();
  
    electroninvPalette->cd();

    hRmsvsTime->GetZaxis()->SetRangeUser(0,Rmsmax);   
    PGlobals::SetH1LabelSize(hRmsvsTime);
    hRmsvsTime->GetZaxis()->SetTitle("#sigma_{x} [#mum]");
    hRmsvsTime->Draw("colz");

    gPad->Update();
  
    palette = (TPaletteAxis*)hRmsvsTime->GetListOfFunctions()->FindObject("palette");
    palette->SetX1NDC(1-margins[2]+0.02);
    palette->SetX2NDC(1-margins[2]+0.05);
    palette->SetY1NDC(margins[1]);
    palette->SetY2NDC(1-margins[3]);
    C4->Modified();  
    gPad->RedrawAxis("G");
  
    C4->cd(0);
  
    // Print to a file
    // Output file
    fOutName = Form("./%s/Plots/Evolutions/Evolutions-DriverRMS-%s",sim.Data(),sim.Data());
    PGlobals::imgconv(C4,fOutName,opt);
    // ---------------------------------------------------------
  }

  TH2F *hDen1DvsTime = NULL;
  sprintf(hName,"hDenvsTime_1"); 
  hDen1DvsTime = (TH2F*) ifile->Get(hName);
  if(hDen1DvsTime) {
    // Canvas setup
    TCanvas *C5 = new TCanvas("C5","Evolution of the on-axis density",sizex,sizey);
    C5->cd(0);
    gPad->SetFrameLineWidth(2);  
    gPad->SetTickx(1);
    gPad->SetTicky(1);

    Float_t Denmax = hDen1DvsTime->GetMaximum();
  
    electronPalette->cd();

    hDen1DvsTime->GetZaxis()->SetRangeUser(0,Denmax);   
    PGlobals::SetH1LabelSize(hDen1DvsTime); 
    hDen1DvsTime->Draw("colz");

    gPad->Update();
  
    palette = (TPaletteAxis*)hDen1DvsTime->GetListOfFunctions()->FindObject("palette");
    palette->SetX1NDC(1-margins[2]+0.02);
    palette->SetX2NDC(1-margins[2]+0.05);
    palette->SetY1NDC(margins[1]);
    palette->SetY2NDC(1-margins[3]);
    C5->Modified();  
    gPad->RedrawAxis("G");
  
    C5->cd(0);
  
    // Print to a file
    // Output file
    fOutName = Form("./%s/Plots/Evolutions/Evolutions-DriverDen-%s",sim.Data(),sim.Data());
    PGlobals::imgconv(C5,fOutName,opt);
    // ---------------------------------------------------------
  }

  // Canvas setup
  
  for(Int_t i=0;i<NAtoms;i++) {

    if(hIonProbvsTime[i]) {
      
      TCanvas *CP = new TCanvas("CP","Evolution of the ionization probability",sizex,sizey);
    
      CP->cd(0);
      CP->Clear();

      gPad->SetFrameLineWidth(2);  
      gPad->SetTickx(1);
      gPad->SetTicky(1);

      PGlobals::SetH1LabelSize(hIonProbvsTime[i]);
 
      redPalette->cd();
      hIonProbvsTime[i]->Draw("colz");
    
      // Plot ionzization limits
      if( gIonProb10[i] )
	gIonProb10[i]->Draw("L");

      if( gIonProb100[i] )
	gIonProb100[i]->Draw("L");
    
      for(Int_t ic=0;ic<NVCross;ic++) {
	if( gVcross[ic] )
	  gVcross[ic]->Draw("L");
      }
    
    
      if(gEcross[1][0]) {
	gEcross[1][0]->SetLineStyle(3);
	gEcross[1][0]->Draw("L");
      }
    
    
      gPad->Update();
      palette = (TPaletteAxis*)hIonProbvsTime[i]->GetListOfFunctions()->FindObject("palette");
      palette->SetX1NDC(1-margins[2]+0.02);
      palette->SetX2NDC(1-margins[2]+0.05);
      palette->SetY1NDC(margins[1]);
      palette->SetY2NDC(1-margins[3]);
      CP->Modified();
      CP->Update();
    
      gPad->RedrawAxis("G");
    
      // Print to a file
      // Output file
      fOutName = Form("./%s/Plots/Evolutions/Evolutions-IonProb-%s-%s",sim.Data(),atNames[i],sim.Data());
      PGlobals::imgconv(CP,fOutName,opt);
    }
    // ---------------------------------------------------------
    
  } // end Atoms loop
    
  
  ifile->Close();
  cout << endl;
}
  
