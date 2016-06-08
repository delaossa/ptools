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
#include <TPaveText.h>
#include <TPaletteAxis.h>
#include <TExec.h>
#include <TClonesArray.h>


#include "PData.hh"
#include "PGlobals.hh"
#include "PPalette.hh"

void PlotChargeFancy2D( const TString &sim, Int_t time, Float_t zoom=2, Int_t NonBin=2, const TString &options="") {
  
#ifdef __CINT__  
  gSystem->Load("libptools.so");
#endif

  PGlobals::Initialize();

  // Palettes!
  gROOT->Macro("PPalettes.C");

  // Init Units table
  PUnits::UnitsTable::Get();
  
  // Load PData
  PData *pData = PData::Get(sim.Data());
  pData->LoadFileNames(time);
  if(!pData->IsInit()) return;

  TString opt = options;

  gStyle->SetPadGridY(0);
  if(opt.Contains("gridx")) {
    gStyle->SetPadGridX(1);
  }
  if(opt.Contains("gridy")) {
    gStyle->SetPadGridY(1);
  }
  gStyle->SetNumberContours(255);
  
  // Some plasma constants
  Double_t n0 = pData->GetPlasmaDensity();
  Double_t kp = pData->GetPlasmaK();
  Double_t skindepth = 1.;
  if(kp!=0.0) skindepth = 1/kp;
  Double_t E0 = pData->GetPlasmaE0();
  
  // Some beam properties:
  Double_t Ebeam = pData->GetBeamEnergy();
  Double_t gamma = pData->GetBeamGamma();
  Double_t vbeam = pData->GetBeamVelocity();

  cout << Form(" - Bunch gamma      = %8.4f", gamma ) << endl;
  cout << Form(" - Bunch velocity   = %8.4f c", vbeam ) << endl;

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
  
  if(opt.Contains("center")) {
    Time -= zStartPlasma;
    if(opt.Contains("comov"))      // Centers on the head of the beam.
      Time += zStartBeam;
  }
  Float_t shiftz = pData->Shift(opt);
  //  cout << "Shift = " << shiftz << endl;
  

  // Calculate the "axis range" in number of bins. If NonBin==0 a RMS width is taken.
  Double_t rms0 = pData->GetBeamRmsY() * kp;
  if(pData->IsCyl())  rms0  = pData->GetBeamRmsR() * kp;

  // Calculate the "axis range" in number of bins. 
  if(NonBin==0) {  // If NonBin==0 a RMS width is taken.
      
    if(rms0>0.0)
      NonBin =  TMath::Nint(rms0 / pData->GetDX(1));
    else
      NonBin = 1;
    
  } else if(NonBin<0) { // If negative, the number is taken in units of 10/kp
    NonBin = -TMath::Nint(float(NonBin*0.1) / pData->GetDX(1));
  }
  
  Int_t FirstxBin = 0;
  Int_t LastxBin = 0;
  if(NonBin==0) { 
    if(rms0>0.0)
      NonBin =  TMath::Nint(rms0 / pData->GetDX(1));
    else
      NonBin = 1;
  }
  
  // Slice width limits.
  if(!pData->IsCyl()) {
    FirstxBin = pData->GetNX(1)/2 + 1 - NonBin;
    LastxBin =  pData->GetNX(1)/2 + NonBin;
  } else {
    FirstxBin = 1; 
    LastxBin  = NonBin;
  }

  // Zoom window:  
  Double_t xRange = (pData->GetXMax(1) - pData->GetXMin(1))/zoom;
  Double_t xMid   = (pData->GetXMax(1) + pData->GetXMin(1))/2.;
  Double_t xMin = xMid - xRange/2.0;
  Double_t xMax = xMid + xRange/2.0;
  if(pData->IsCyl()) {
    xMin = pData->GetXMin(1);
    xMax = xRange;
  }
  pData->SetX2Min(xMin);
  pData->SetX2Max(xMax);


  // ----------------------------------------------------------------------------------
  
  
  // Get charge density histos
  Int_t Nspecies = pData->NSpecies();
  TH2F **hDen2D = new TH2F*[Nspecies];
  // Get charge density on-axis
  TH1F **hDen1D = new TH1F*[Nspecies];
  for(Int_t i=0;i<Nspecies;i++) {
 
    
    hDen2D[i] = NULL;
    
    if(!pData->GetChargeFileName(i)) 
      continue;

    cout << Form(" Getting charge density of specie: ") << i << endl;

    
    char hName[24];
    sprintf(hName,"hDen2D_%i",i);
    hDen2D[i] = (TH2F*) gROOT->FindObject(hName);
    if(hDen2D[i]) delete hDen2D[i];
    
    if(!pData->Is3D())
      hDen2D[i] = pData->GetCharge(i,opt);
    else
      hDen2D[i] = pData->GetCharge2DSliceZX(i,-1,NonBin,opt+"avg");
    
    hDen2D[i]->SetName(hName);
    hDen2D[i]->GetXaxis()->CenterTitle();
    hDen2D[i]->GetYaxis()->CenterTitle();
    hDen2D[i]->GetZaxis()->CenterTitle();
    
    if(opt.Contains("comov"))
      hDen2D[i]->GetXaxis()->SetTitle("k_{p} #zeta");
    else
      hDen2D[i]->GetXaxis()->SetTitle("k_{p} z");
    
    if(pData->IsCyl()) 
      hDen2D[i]->GetYaxis()->SetTitle("k_{p} r");
    else
      hDen2D[i]->GetYaxis()->SetTitle("k_{p} y");
    
    hDen2D[i]->GetZaxis()->SetTitle("n [n_{0}]");
  }
  

  // Get electric fields 2D
  const Int_t Nfields = 3;
  TH1F **hE1D = new TH1F*[Nfields];
  TH2F **hE2D = new TH2F*[Nfields];
  for(Int_t i=0;i<Nfields;i++) {
    hE1D[i] = NULL;
    hE2D[i] = NULL;
    
    if(!opt.Contains("etotal") && i>0) 
      continue;

    if(!pData->GetEfieldFileName(i))
      continue;

    cout << Form(" Getting electric field number ") << i+1 << endl;
    
    char hName[24];
    sprintf(hName,"hE1D_%i",i);
    hE1D[i] = (TH1F*) gROOT->FindObject(hName);
    if(hE1D[i]) delete hE1D[i];
      
    // 1D histograms
    char nam[3]; sprintf(nam,"e%i",i+1);
    if(pData->Is3D()) {
      
      if(i==0) 
	hE1D[i] = pData->GetH1SliceZ3D(pData->GetEfieldFileName(i)->c_str(),nam,-1,NonBin,-1,NonBin,opt+"avg");
      else  
	hE1D[i] = pData->GetH1SliceZ3D(pData->GetEfieldFileName(i)->c_str(),nam,-NonBin,NonBin,-NonBin,NonBin,opt+"avg");
      
    } else if(pData->IsCyl()) { // Cylindrical: The first bin with r>0 is actually the number 1 (not the 0).
      
      hE1D[i] = pData->GetH1SliceZ(pData->GetEfieldFileName(i)->c_str(),nam,1,NonBin,opt+"avg");
      
    } else { // 2D cartesian
      
      if(i==0) 
	hE1D[i] = pData->GetH1SliceZ(pData->GetEfieldFileName(i)->c_str(),nam,-1,NonBin,opt+"avg");
      else 
	hE1D[i] = pData->GetH1SliceZ(pData->GetEfieldFileName(i)->c_str(),nam,-NonBin,NonBin,opt+"avg");    
    }
    
    hE1D[i]->SetName(hName);
    if(opt.Contains("comov"))
      hE1D[i]->GetXaxis()->SetTitle("#zeta [c/#omega_{p}]");
    else
      hE1D[i]->GetXaxis()->SetTitle("z [c/#omega_{p}]");
    
    if(i==0)
      hE1D[i]->GetYaxis()->SetTitle("E_{z} [E_{0}]");
    else if(i==1)
      hE1D[i]->GetYaxis()->SetTitle("E_{y} [E_{0}]");
    else if(i==2)
      hE1D[i]->GetYaxis()->SetTitle("E_{x} [E_{0}]");

    if(opt.Contains("pot") && i==0) {
      sprintf(hName,"hE2D_%i",i);
      hE2D[i] = (TH2F*) gROOT->FindObject(hName);
      if(hE2D[i]) delete hE2D[i];
      
      if(!pData->Is3D())
	hE2D[i] = pData->GetEField(i,opt);
      else
	hE2D[i] = pData->GetEField2DSliceZY(i,-1,NonBin,opt+"avg");
    } 
    
  }

  // Computing on-axis total field 
  TH1F *hETotal1D = (TH1F*) hE1D[0]->Clone("hETotal1D");
  hETotal1D->Reset();
  if(opt.Contains("etotal")) {
    Float_t integral = 0.0;
    Float_t dx = pData->GetDX(0);
    
    Int_t NbinsX = hE1D[0]->GetNbinsX();
    for(Int_t j=NbinsX;j>=0;j--) {
      Float_t E1 = hE1D[0]->GetBinContent(j);
      Float_t E2 = hE1D[1]->GetBinContent(j);
      Float_t E3 = hE1D[2]->GetBinContent(j);
      Float_t E  = TMath::Sqrt(E1*E1+E2*E2+E3*E3);
      
      hETotal1D->SetBinContent(j,E);
    }
  }
  
  // Computing trapping potential
  TH1F *hV1D = NULL;
  TH2F *hV2D = NULL;
  Float_t Vmin = 0.0;
  Float_t Vmax = 0.0;
  // define the contours
  const Int_t Ncontours = 25;
  Double_t contours[Ncontours];
  for(Int_t i=0; i<Ncontours; i++) {
    contours[i] = i*(trapPotential/5.0) - trapPotential; 
  }
  TClonesArray *graphsV2D = NULL;
  if(opt.Contains("pot")) {
    Int_t   NbinsX = hE2D[0]->GetNbinsX();
    Int_t   NbinsY = hE2D[0]->GetNbinsY();

    Float_t dx = pData->GetDX(0);
    
    char hName[24];  
    sprintf(hName,"hV2D");
    hV2D = (TH2F*) hE2D[0]->Clone(hName);
    hV2D->Reset();
    
    sprintf(hName,"hV1D");
    hV1D = (TH1F*) hE1D[0]->Clone(hName);
    hV1D->Reset();

    for(Int_t j=NbinsY;j>0;j--) {
      Float_t integral = 0.0;
      for(Int_t k=NbinsX;k>0;k--) {
	integral += hE2D[0]->GetBinContent(k,j) * dx;
	hV2D->SetBinContent(k,j,integral);
      }
    }

    Float_t integral = 0.0;
    for(Int_t k=NbinsX;k>0;k--) {
      integral += hE1D[0]->GetBinContent(k) * dx;
      hV1D->SetBinContent(k,integral);
    }
   
    Vmin = hV1D->GetMinimum();    
    { // Shift potential
      Int_t   NbinsX = hV2D->GetNbinsX(); 
      Int_t   NbinsY = hV2D->GetNbinsY(); 
      for(Int_t j=0;j<NbinsX;j++) {
	for(Int_t k=0;k<NbinsY;k++) {
	  hV2D->SetBinContent(j,k, hV2D->GetBinContent(j,k) - Vmin -trapPotential);
	}
	hV1D->SetBinContent(j, hV1D->GetBinContent(j) - Vmin -trapPotential);
      }
    }

    Vmin = hV1D->GetMinimum();   
    Vmax = hV1D->GetMaximum();    
    if(Vmax<0.1) Vmax = 0.1;
   

    // Extract contours
    TCanvas* c = new TCanvas("c","Contour List",0,0,600,600);
    c->cd();
  
    // Potential
    TH2F *hV2Dc = (TH2F*) hV2D->Clone("hV2Dc");
    hV2Dc->SetContour(Ncontours, contours);
    hV2Dc->Draw("cont list");
  
    c->Update();
    TObjArray *contsV2D = (TObjArray*) gROOT->GetListOfSpecials()->FindObject("contours");
    graphsV2D = new TClonesArray("TGraph",Ncontours);
    {
      Int_t ncontours = contsV2D->GetSize();
      TList* clist = NULL;
      Int_t nGraphs = 0;
      TGraph *gr = NULL;
      for(Int_t i = 0; i < ncontours; i++){
	if(i==0) continue;
	
	clist = (TList*) contsV2D->At(i);
	
	for(Int_t j = 0 ; j < clist->GetSize(); j++) {
	  gr = (TGraph*) clist->At(j);
	  if(!gr) continue;
	  
	  gr->SetLineWidth(2);
	  gr->SetLineColor(kGray+2);
	  gr->SetLineStyle(2);
	  
	  if(i-1==5 && j==1) {
	    new((*graphsV2D)[nGraphs]) TGraph(*gr);
	    //graphsV2D->AddLast(gr);
	    nGraphs++;
	  }
	  	  
	}
      }
    }
  }    
  

  
  // Tunning the Histograms
  // ---------------------
    
  Float_t zMin = hDen2D[0]->GetXaxis()->GetXmin();
  Float_t zMax = hDen2D[0]->GetXaxis()->GetXmax();
  // Float_t zRange = zMax - zMin;

  // ------------- z Zoom --------------------------------- Plasma palette -----------
  // Set the range of the plasma charge density histogram for maximum constrast 
  // using a dynamic palette wich adjust the nominal value to a certain color.
  

  Float_t density = 1;
  Float_t Base  = density;
  
  Float_t *Max = new Float_t[Nspecies];
  Float_t *Min = new Float_t[Nspecies];
  
  for(Int_t i=0;i<Nspecies;i++) {
    if(!hDen2D[i]) continue;
   
    Max[i] = hDen2D[i]->GetMaximum();
    Min[i] = 1.01E-1 * Base;
    if(i==0) {
      if(Max[i]<1) {
	Max[i] = 1.1*Base;
      } else if(Max[i]<2) {
	Min[i] = 2*Base - Max[i];
      } else {
	Max[i] = 0.4*hDen2D[i]->GetMaximum(); // enhance plasma contrast.
      }
    }
        
    if(i==1) Min[i] = 5.01E-1 * Base;
    if(i==2) Min[i] = 1.01E-1 * Base;

    if(pData->GetDenMax(i)>0)
      Max[i] = pData->GetDenMax(i);
    else if(pData->GetDenMax(i)==0) {
      delete hDen2D[i];
      hDen2D[i] = NULL;
    }
    
    if(pData->GetDenMin(i)>=0)
      Min[i] = pData->GetDenMin(i);
    
    hDen2D[i]->GetZaxis()->SetRangeUser(Min[i],Max[i]);
  }
  
  // Dynamic plasma palette
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

  // Palette for plasma
  const Int_t NRGBs = 3;
  const Int_t NCont = 255;
  Double_t Stops[NRGBs] = { 0.00, basePos, 1.00 };
  Int_t cindex[NRGBs];
  if(opt.Contains("dark")) {
    cindex[0] = TColor::GetColor("#141515"); // grey11
    cindex[1] = TColor::GetColor("#303F4B");   // grey31
    // cindex[1] = TColor::GetColor("#1E2D78"); // lighter blue but dark
    // cindex[1] = TColor::GetColor("#121F40"); // dark blue
    cindex[2] = TColor::GetColor("#FFFFFF"); // white
  } else {
    cindex[0] = TColor::GetColor("#FFFFFF"); // white
    cindex[1] = TColor::GetColor((Float_t) 0.90,(Float_t) 0.90,(Float_t) 0.90);  // light gray
    cindex[2] = TColor::GetColor((Float_t) 0.00,(Float_t) 0.00,(Float_t) 0.00);  // black
  }
  
  PPalette * plasmaPalette = (PPalette*) gROOT->FindObject("plasma");
  if(!plasmaPalette) {
    plasmaPalette = new PPalette("plasma");
    plasmaPalette->CreateGradientColorTable(NRGBs, Stops, cindex, NCont);
  } else {
    plasmaPalette->ChangeGradientColorTable(NRGBs, Stops, cindex);
  }
  
  if(opt.Contains("alpha")) {
    if(opt.Contains("dark")) 
      plasmaPalette->SetAlpha(0.7);
    else
      plasmaPalette->SetAlpha(0.9);
  }
  
  // Palette for beam
  PPalette * beamPalette = (PPalette*) gROOT->FindObject("beam");
  if(!beamPalette) {
    beamPalette = new PPalette("beam");
  }

  if(opt.Contains("dark")) {
    const Int_t bNRGBs = 4;
    const Int_t bNCont = 255;
    Double_t bStops[bNRGBs] = { 0.00, 0.20, 0.40, 1.00};
    Int_t bcindex[bNRGBs];
    
    bcindex[0] = TColor::GetColor("#380A3C"); // deep purple
    //bcindex[0] = cindex[1];
    bcindex[1] = TColor::GetColor((Float_t) 0.39, (Float_t) 0.05, (Float_t) 0.33); // dark magenta 
    bcindex[2] = TColor::GetColor((Float_t) 0.70, (Float_t) 0.20, (Float_t) 0.30); // pinky orange
    bcindex[3] = TColor::GetColor((Float_t) 1.00, (Float_t) 1.00, (Float_t) 0.20); // yellow
    beamPalette->ChangeGradientColorTable(bNRGBs, bStops, bcindex);
  } else {
    const Int_t bNRGBs = 5;
    const Int_t bNCont = 255;
    Double_t bStops[bNRGBs] = { 0.00, 0.30, 0.45, 0.55, 1.00};
    Int_t bcindex[bNRGBs];
    bcindex[0] = cindex[1];
    bcindex[1] = TColor::GetColor("#386EA5"); // steel blue
    bcindex[2] = TColor::GetColor((Float_t) 0.39, (Float_t) 0.05, (Float_t) 0.33); // dark magenta
    bcindex[3] = TColor::GetColor((Float_t) 0.70, (Float_t) 0.20, (Float_t) 0.30); // pinky orange
    bcindex[4] = TColor::GetColor((Float_t) 1.00, (Float_t) 1.00, (Float_t) 0.20); // yellow
    beamPalette->ChangeGradientColorTable(bNRGBs, bStops, bcindex);
    //    beamPalette->SetPalette("elec");
  }
  
  PPalette * beam2Palette = (PPalette*) gROOT->FindObject("beam2");
  if(!beam2Palette) {
    beam2Palette = new PPalette("beam2");
    beam2Palette->SetPalette("hot");
  }
  // if(opt.Contains("dark")) {
  //   bcindex[0] = bcindex[0];
  //   bcindex[1] = bcindex[2];    
  //   bStops[1] = bStops[2] = 0.1;
  //   beam2Palette->ChangeGradientColorTable(bNRGBs, bStops, bcindex);
  // } 
      
  // Change the range of z axis for the fields to be symmetric.
  Float_t *Emax = new Float_t[Nfields];
  Float_t *Emin = new Float_t[Nfields];
  for(Int_t i=0;i<Nfields;i++) {
    if(!hE1D[i]) continue;
    Emax[i] = hE1D[i]->GetMaximum();
    Emin[i] = hE1D[i]->GetMinimum();
    if(Emax[i] > TMath::Abs(Emin[i]))
      Emin[i] = -Emax[i];
    else
      Emax[i] = -Emin[i];
    hE1D[i]->GetZaxis()->SetRangeUser(Emin[i],Emax[i]); 
  }
  
  
  // "Axis range" in Osiris units:
  Double_t ylow  = hDen2D[0]->GetYaxis()->GetBinLowEdge(FirstxBin);
  Double_t yup = hDen2D[0]->GetYaxis()->GetBinUpEdge(LastxBin);

  zStartPlasma -= shiftz; 
  zStartNeutral -= shiftz; 
  zEndNeutral -= shiftz; 
  
  if(opt.Contains("units")) {
    zStartPlasma *= skindepth / PUnits::um;
    zStartNeutral *= skindepth / PUnits::um;
    zEndNeutral *= skindepth / PUnits::um;
  }

  //  cout << "Start plasma = " << zStartPlasma << endl;
  TLine *lineStartPlasma = new TLine(zStartPlasma,xMin,zStartPlasma,xMax);
  lineStartPlasma->SetLineColor(kGray+2);
  lineStartPlasma->SetLineStyle(2);
  lineStartPlasma->SetLineWidth(3);

  //  cout << "Start neutral = " << zStartNeutral << endl;
  TLine *lineStartNeutral = new TLine(zStartNeutral,xMin,zStartNeutral,xMax);
  lineStartNeutral->SetLineColor(kGray+3);
  lineStartNeutral->SetLineStyle(3);
  lineStartNeutral->SetLineWidth(2);

  //  cout << "End neutral = " << zEndNeutral << endl;
  TLine *lineEndNeutral = new TLine(zEndNeutral,xMin,zEndNeutral,xMax);
  lineEndNeutral->SetLineColor(kGray+3);
  lineEndNeutral->SetLineStyle(3);
  lineEndNeutral->SetLineWidth(2);
  

  // Plotting
  // -----------------------------------------------

  // Canvas setup
  TCanvas *C = new TCanvas("C","",1024,576);
 
  // Palettes setup
  TExec *exPlasma = new TExec("exPlasma","plasmaPalette->cd();");
  TExec *exDriver   = new TExec("exDriver","beamPalette->cd();");
  TExec *exWitness    = new TExec("exWitness","beam2Palette->cd();");

  if(opt.Contains("alpha0"))
    plasmaPalette->SetAlpha(0.8);
  else if(opt.Contains("alpha1"))
    beamPalette->SetAlpha(0.8);
  else if(opt.Contains("alpha2"))
    beam2Palette->SetAlpha(0.8);
  
  // Actual Plotting!
  // ------------------------------------------------------------

  // Output file
  TString fOutName = Form("./%s/Plots/ChargeFancy2D/ChargeFancy2D",pData->GetPath().c_str());
  fOutName += Form("-%s_%i",pData->GetName(),time);

  // Setup Pad layout:
  const Int_t NPad = 1;
  TPad *pad[NPad];
  TH2F *hFrame[NPad];

  // Text objects
  TPaveText **textLabel = new TPaveText*[NPad];
  Float_t lMargin = 0.0;
  Float_t rMargin = 0.0;
  Float_t bMargin = 0.0;
  Float_t tMargin = 0.0;
  PGlobals::CanvasPartition(C,NPad,lMargin,rMargin,bMargin,tMargin);
  
  // Define the frames for plotting
  Int_t fonttype = 43;
  Int_t fontsize = 24;
  Int_t tfontsize = 28;
  Float_t txoffset = 2.0;
  Float_t lxoffset = 0.02;
  Float_t tyoffset = 1.3;
  Float_t lyoffset = 0.01;
  Float_t tylength = 0.0;
  Float_t txlength = 0.0;
  for(Int_t i=0;i<NPad;i++) {
    char name[16];
    sprintf(name,"pad_%i",i);
    pad[i] = (TPad*) gROOT->FindObject(name);
    pad[i]->SetFrameLineWidth(2);  
    // pad[i]->SetTickx(1);
    // pad[i]->SetTicky(1);

    sprintf(name,"hFrame_%i",i);
    hFrame[i] = (TH2F*) gROOT->FindObject(name);
    if(hFrame[i]) delete hFrame[i];
    hFrame[i] = (TH2F*) hDen2D[0]->Clone(name);
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

    // Format for x axis
    hFrame[i]->GetXaxis()->SetTitleFont(fonttype);
    hFrame[i]->GetXaxis()->SetTitleSize(tfontsize+2);
    hFrame[i]->GetXaxis()->SetTitleOffset(txoffset);
    hFrame[i]->GetXaxis()->SetLabelFont(fonttype);
    hFrame[i]->GetXaxis()->SetLabelSize(fontsize+2);
    hFrame[i]->GetXaxis()->SetLabelOffset(lxoffset);
    
    hFrame[i]->GetXaxis()->SetTickLength(yFactor*txlength/xFactor);      
  }


  C->cd(0);
  pad[0]->Draw();
  pad[0]->cd(); // <---------------------------------------------- Top Plot ---------

  if(opt.Contains("dark")) {
    Int_t backcolor = TColor::GetColor(10,10,10);
    pad[0]->SetFillColor(backcolor);
  }
  
  if(opt.Contains("logz")) {
    pad[0]->SetLogz(1);
  } else {
    pad[0]->SetLogz(0);
  }

  hFrame[0]->Draw("AH");

  // Injected electrons if any
  if(Nspecies>=3) {
    if(hDen2D[2]) {
      exWitness->Draw();
      hDen2D[2]->Draw("col same");
    }
  }
    
  // Plasma
  hDen2D[0]->GetZaxis()->SetTitleFont(fonttype);
  exPlasma->Draw();
  hDen2D[0]->Draw("col same");
  
  // Beam driver.
  if(hDen2D[1]) {
    //    hDen2D[1]->GetZaxis()->SetNdivisions(505);
    exDriver->Draw();
    hDen2D[1]->Draw("col same");
  }

  if(!opt.Contains("etotal")) {
    // Fit the E1D in the pad:
    Float_t x1 = xMin + xRange/80.;
    Float_t x2 = x1 + xRange/2.8 - 2*xRange/80.;
    Float_t slope = (x2-x1)/(Emax[0]-Emin[0]);
  
    Int_t lineColor = TColor::GetColor(196,30,78);
    //Int_t lineColor = TColor::GetColor(236,78,35);

    TLine *lineEzero = new TLine(zMin,(0.0-Emin[0])*slope + x1,zMax,(0.0-Emin[0])*slope + x1);
    lineEzero->SetLineColor(lineColor);
    lineEzero->SetLineStyle(2);
    lineEzero->SetLineWidth(1);
    lineEzero->Draw();

    for(Int_t j=0;j<hE1D[0]->GetNbinsX();j++) {
      hE1D[0]->SetBinContent(j+1,(hE1D[0]->GetBinContent(j+1)-Emin[0])*slope + x1);
    }
    hE1D[0]->SetLineStyle(1);
    hE1D[0]->SetLineWidth(3);
        
    hE1D[0]->SetLineColor(lineColor);

    hE1D[0]->Draw("sameC");
  } else {
    // Fit the E1D in the E2D pad:
    Float_t HeIon  = 92.75;  // GV/m
    if(!opt.Contains("units") || !n0 ) {
      HeIon  /= ( E0 / (PUnits::GV/PUnits::m));
    }

    // Fit the E1D in the pad:
    Float_t ETmin =  0.00;//hETotal1D->GetMinimum();
    Float_t ETmax =  hETotal1D->GetMaximum();
    Float_t x1 = xMin + xRange/80.;
    Float_t x2 = x1 + xRange/2.8 - 2*xRange/80.;
    Float_t slope = (x2-x1)/(ETmax-ETmin);
  
    Int_t lineColor = TColor::GetColor(196,30,78);
    //Int_t lineColor = TColor::GetColor(236,78,35);

    TLine *lineEHe = new TLine(zMin,(HeIon-ETmin)*slope + x1,zMax,(HeIon-ETmin)*slope + x1);
    lineEHe->SetLineColor(lineColor);
    lineEHe->SetLineStyle(2);
    lineEHe->SetLineWidth(1);
    lineEHe->Draw();

    TLine *lineEzero = new TLine(zMin,(0.0-ETmin)*slope + x1,zMax,(0.0-ETmin)*slope + x1);
    lineEzero->SetLineColor(kGray+1);
    lineEzero->SetLineStyle(2);
    lineEzero->SetLineWidth(1);
    lineEzero->Draw();
       
    for(Int_t j=0;j<hETotal1D->GetNbinsX();j++) {
      hETotal1D->SetBinContent(j+1,(hETotal1D->GetBinContent(j+1)-ETmin)*slope + x1);
    }
    hETotal1D->SetLineStyle(1);
    hETotal1D->SetLineWidth(2);
    
    hETotal1D->SetLineColor(lineColor);
    
    hETotal1D->Draw("sameC");
  }

  if(opt.Contains("pot")) {
    //  PSI MAIN contour  
    for(Int_t i=0;i<graphsV2D->GetEntriesFast();i++) {
      TGraph *gr = (TGraph*) graphsV2D->At(i);
      if(!gr) continue;
      gr->Draw("C");
      
    }
  }

  if(opt.Contains("sline")) {
    if(zStartPlasma>zMin && zStartPlasma<zMax)
      lineStartPlasma->Draw();
    if(zStartNeutral>zMin && zStartNeutral<zMax)
      lineStartNeutral->Draw();
    if(zEndNeutral>zMin && zEndNeutral<zMax)
      lineEndNeutral->Draw();
  }
  
  pad[0]->RedrawAxis(); 
 
  C->cd(0);
  C->cd();

  // Print to a file
  PGlobals::imgconv(C,fOutName,opt);
  // ---------------------------------------------------------

  PGlobals::DestroyCanvases();
}
  
