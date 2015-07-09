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

void PlotChargeFancy2D( const TString &sim, Int_t time, Float_t zoom=2, Int_t Nbins=2, const TString &options="") {
  
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
  

  // Calculate the "axis range" in number of bins. If Nbins==0 a RMS width is taken.
  Double_t rms0 = pData->GetBeamRmsY() * kp;
  if(pData->IsCyl())  rms0  = pData->GetBeamRmsR() * kp;
  
  Int_t FirstyBin = 0;
  Int_t LastyBin = 0;
  if(Nbins==0) { 
    if(rms0>0.0)
      Nbins =  TMath::Nint(rms0 / pData->GetDX(1));
    else
      Nbins = 1;
  }
  
  // Slice width limits.
  if(!pData->IsCyl()) {
    FirstyBin = pData->GetNX(1)/2 + 1 - Nbins;
    LastyBin =  pData->GetNX(1)/2 + Nbins;
  } else {
    FirstyBin = 1; 
    LastyBin  = Nbins;
  }


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
      hDen2D[i] = pData->GetCharge2DSliceZY(i,-1,Nbins,opt+"avg");
    
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
	hE1D[i] = pData->GetH1SliceZ3D(pData->GetEfieldFileName(i)->c_str(),nam,-1,Nbins,-1,Nbins,opt+"avg");
      else  
	hE1D[i] = pData->GetH1SliceZ3D(pData->GetEfieldFileName(i)->c_str(),nam,-Nbins,Nbins,-Nbins,Nbins,opt+"avg");
      
    } else if(pData->IsCyl()) { // Cylindrical: The first bin with r>0 is actually the number 1 (not the 0).
      
      hE1D[i] = pData->GetH1SliceZ(pData->GetEfieldFileName(i)->c_str(),nam,1,Nbins,opt+"avg");
      
    } else { // 2D cartesian
      
      if(i==0) 
	hE1D[i] = pData->GetH1SliceZ(pData->GetEfieldFileName(i)->c_str(),nam,-1,Nbins,opt+"avg");
      else 
	hE1D[i] = pData->GetH1SliceZ(pData->GetEfieldFileName(i)->c_str(),nam,-Nbins,Nbins,opt+"avg");    
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
	hE2D[i] = pData->GetEField2DSliceZY(i,-1,Nbins,opt+"avg");
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
  

  // --------------------------------------------------- Vertical Zoom ------------
  
  Float_t yRange   = (hDen2D[0]->GetYaxis()->GetXmax() - hDen2D[0]->GetYaxis()->GetXmin())/zoom;
  Float_t midPoint = (hDen2D[0]->GetYaxis()->GetXmax() + hDen2D[0]->GetYaxis()->GetXmin())/2.;
  Float_t yMin = midPoint-yRange/2.0;
  Float_t yMax = midPoint+yRange/2.0;
  if(pData->IsCyl()) {
    yMin = pData->GetXMin(1);
    yMax = yRange;
  }

  //  cout << Form(" %f  %f   %f  %f",hDen2D[0]->GetYaxis()->GetXmin(),hDen2D[0]->GetYaxis()->GetXmax(),yMin,yMax) << endl;

  for(Int_t i=0;i<Nspecies;i++) {
    if(!hDen2D[i]) continue;
    hDen2D[i]->GetYaxis()->SetRangeUser(yMin,yMax);
  }
  
  Float_t xMin = hDen2D[0]->GetXaxis()->GetXmin();
  Float_t xMax = hDen2D[0]->GetXaxis()->GetXmax();
  Float_t xRange = xMax - xMin;

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
        
    if(i==1) Min[i] = 1.01E-1 * Base;
    if(i==2) Min[i] = 1.01E-3 * Base;
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

  // FANCY palette for plasma
  // Double_t plasmaDStops[plasmaDNRGBs] = { 0.00, basePos, 1.00 };
  // Double_t plasmaDRed[plasmaDNRGBs]   = { 0.04, 0.09, 1.00 };
  // Double_t plasmaDGreen[plasmaDNRGBs] = { 0.04, 0.17, 1.00 };
  // Double_t plasmaDBlue[plasmaDNRGBs]  = { 0.04, 0.32, 1.00 };
   
  // PPalette * plasmaPalette = (PPalette*) gROOT->FindObject("plasma");
  // plasmaPalette->CreateGradientColorTable(plasmaDNRGBs, plasmaDStops, 
  // 					  plasmaDRed, plasmaDGreen, plasmaDBlue, plasmaDNCont);

  Double_t plasmaDStops[plasmaDNRGBs] = { 0.00, basePos, 1.00 };
  Double_t plasmaDRed[plasmaDNRGBs]   = { 0.99, 0.90, 0.00 };
  Double_t plasmaDGreen[plasmaDNRGBs] = { 0.99, 0.90, 0.00 };
  Double_t plasmaDBlue[plasmaDNRGBs]  = { 0.99, 0.90, 0.00 };
   
  PPalette * plasmaPalette = (PPalette*) gROOT->FindObject("plasma");
  plasmaPalette->CreateGradientColorTable(plasmaDNRGBs, plasmaDStops, 
					  plasmaDRed, plasmaDGreen, plasmaDBlue, plasmaDNCont);


  // const Int_t redelectronNRGBs = 4;
  // const Int_t redelectronNCont = 64;
  // Double_t redelectronStops[redelectronNRGBs] = { 0.30, 0.40, 0.60, 1.00};
  // Double_t redelectronRed[redelectronNRGBs] =   { 0.09, 0.39, 0.70, 1.00};
  // Double_t redelectronGreen[redelectronNRGBs] = { 0.17, 0.05, 0.20, 1.00};
  // Double_t redelectronBlue[redelectronNRGBs] =  { 0.32, 0.33, 0.30, 0.20};

  // PPalette * redelectronPalette = (PPalette*) gROOT->FindObject("redelectron");
  // redelectronPalette->CreateGradientColorTable(redelectronNRGBs, redelectronStops, 
  // 					       redelectronRed, redelectronGreen, redelectronBlue, redelectronNCont);
  

  // const Int_t hotNRGBs = 3;
  // const Int_t hotNCont = 64;
  // Double_t hotStops[hotNRGBs] =  { 0.25, 0.5, 1.00 };
  // Double_t hotRed[hotNRGBs] =   { 1.00, 1.000, 1.000 };
  // Double_t hotGreen[hotNRGBs] = { 0.15, 0.984, 1.000 };
  // Double_t hotBlue[hotNRGBs] =  { 0.00, 0.000, 1.000 };
  // PPalette * hotPalette = (PPalette*) gROOT->FindObject("hot");
  // hotPalette->CreateGradientColorTable(hotNRGBs, hotStops,hotRed, hotGreen, hotBlue, hotNCont);

  
    
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
  Double_t ylow  = hDen2D[0]->GetYaxis()->GetBinLowEdge(FirstyBin);
  Double_t yup = hDen2D[0]->GetYaxis()->GetBinUpEdge(LastyBin);

  zStartPlasma -= shiftz; 
  zStartNeutral -= shiftz; 
  zEndNeutral -= shiftz; 
  
  if(opt.Contains("units")) {
    zStartPlasma *= skindepth / PUnits::um;
    zStartNeutral *= skindepth / PUnits::um;
    zEndNeutral *= skindepth / PUnits::um;
  }

  // Plotting
  // -----------------------------------------------

  // Canvas setup
  TCanvas *C = new TCanvas("C","",1024,576);
 
  // Palettes setup
  TExec *exPlasma = new TExec("exPlasma","plasmaPalette->cd();");
  TExec *exDriver   = new TExec("exDriver","electronPalette->cd();");
  TExec *exWitness    = new TExec("exWitness","hotPalette->cd();");
     
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
  
  // Int_t backcolor = TColor::GetColor(10,10,10);
  // pad[0]->SetFillColor(backcolor);

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
    Float_t y1 = yMin + yRange/80.;
    Float_t y2 = y1 + yRange/2.8 - 2*yRange/80.;
    Float_t slope = (y2-y1)/(Emax[0]-Emin[0]);
  
    Int_t lineColor = TColor::GetColor(196,30,78);
    //Int_t lineColor = TColor::GetColor(236,78,35);

    TLine *lineEzero = new TLine(xMin,(0.0-Emin[0])*slope + y1,xMax,(0.0-Emin[0])*slope + y1);
    lineEzero->SetLineColor(lineColor);
    lineEzero->SetLineStyle(2);
    lineEzero->SetLineWidth(1);
    lineEzero->Draw();

    for(Int_t j=0;j<hE1D[0]->GetNbinsX();j++) {
      hE1D[0]->SetBinContent(j+1,(hE1D[0]->GetBinContent(j+1)-Emin[0])*slope + y1);
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
    Float_t y1 = yMin + yRange/80.;
    Float_t y2 = y1 + yRange/2.8 - 2*yRange/80.;
    Float_t slope = (y2-y1)/(ETmax-ETmin);
  
    Int_t lineColor = TColor::GetColor(196,30,78);
    //Int_t lineColor = TColor::GetColor(236,78,35);

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

  if(opt.Contains("pot")) {
    //  PSI MAIN contour  
    for(Int_t i=0;i<graphsV2D->GetEntriesFast();i++) {
      TGraph *gr = (TGraph*) graphsV2D->At(i);
      if(!gr) continue;
      gr->Draw("C");
      
    }
  }
  
  pad[0]->RedrawAxis(); 
 
  C->cd(0);
  C->cd();

  // Print to a file
  PGlobals::imgconv(C,fOutName,opt);
  // ---------------------------------------------------------

  PGlobals::DestroyCanvases();
}
  
