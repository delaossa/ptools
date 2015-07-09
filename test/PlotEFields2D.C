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
#include "PlasmaGlob.hh"
#include "PPalette.hh"

void PlotEFields2D( const TString &sim, Int_t time, Int_t zoom=2, Int_t Nbins=2, const TString &options="") {
  
#ifdef __CINT__  
  gSystem->Load("libplasma.so");
#endif

  PlasmaGlob::Initialize();

  // Palettes!
  gROOT->Macro("PlasmaPalettes.C");

  // Init Units table
  PUnits::UnitsTable::Get();
  
  // Load PData
  PData *pData = PData::Get(sim.Data());
  pData->LoadFileNames(time);
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
  Double_t n0 = pData->GetPlasmaDensity();
  Double_t omegap = pData->GetPlasmaFrequency();
  Double_t timedepth = 1.;
  if(omegap!=0.0) timedepth = 1/omegap;
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
  cout << Form(" - Trap. potential  = %8.4f mc2/e",trapPotential) << endl;
  cout << endl;

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
  
  
 
  
  // Get electric fields 2D
  const Int_t Nfields = 3;
  TH2F **hE2D = new TH2F*[Nfields];
  TH1F **hE1D = new TH1F*[Nfields];
  TH2F *hV2D = NULL;
  TH1F *hV1D = NULL;
  for(Int_t i=0;i<Nfields;i++) {
    hE2D[i] = NULL;
    hE1D[i] = NULL;

    if(!pData->GetEfieldFileName(i))
      continue;

    cout << Form(" Getting electric field number ") << i+1 << endl;
    
    char hName[24];
    sprintf(hName,"hE2D_%i",i);
    hE2D[i] = (TH2F*) gROOT->FindObject(hName);
    if(hE2D[i]) delete hE2D[i];
    
    if(!pData->Is3D())
      hE2D[i] = pData->GetEField(i,opt);
    else
      hE2D[i] = pData->GetEField2DSliceZY(i,-1,Nbins,opt+"avg");
    
    hE2D[i]->SetName(hName);   
    hE2D[i]->GetXaxis()->CenterTitle();
    hE2D[i]->GetYaxis()->CenterTitle();
    hE2D[i]->GetZaxis()->CenterTitle();
    if(opt.Contains("comov"))
      hE2D[i]->GetXaxis()->SetTitle("k_{p} #zeta");
    else
      hE2D[i]->GetXaxis()->SetTitle("k_{p} z");
    
    if(pData->IsCyl()) 
      hE2D[i]->GetYaxis()->SetTitle("k_{p} r");
    else
      hE2D[i]->GetYaxis()->SetTitle("k_{p} y");
    
    if(i==0)
      hE2D[i]->GetZaxis()->SetTitle("E_{z}/E_{0}");
    else if(i==1)
      hE2D[i]->GetZaxis()->SetTitle("E_{y}/E_{0}");
    else if(i==2)
      hE2D[i]->GetZaxis()->SetTitle("E_{x}/E_{0}");
    
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
    
    // Alternative
    // if(hE1D[i]) delete hE1D[i];
    // hE1D[i] = (TH1F*) hE2D[i]->ProjectionX(hName,FirstyBin,LastyBin);
    // hE1D[i]->Scale(1.0/(LastyBin-FirstyBin+1));
    
    if(i==0) {
      Int_t   NbinsX = hE2D[i]->GetNbinsX();
      Int_t   NbinsY = hE2D[i]->GetNbinsY();

      Float_t dx = pData->GetDX(0);
     
      sprintf(hName,"hV2D");
      hV2D = (TH2F*) hE2D[i]->Clone(hName);
      hV2D->Reset();

      sprintf(hName,"hV1D");
      hV1D = (TH1F*) hE1D[i]->Clone(hName);
      hV1D->Reset();

      for(Int_t j=NbinsY;j>0;j--) {
	Double_t integral = 0.0;
	for(Int_t k=NbinsX;k>0;k--) {
	  integral += hE2D[i]->GetBinContent(k,j) * dx;
	  hV2D->SetBinContent(k,j,integral);
	}
      }

      Double_t integral = 0.0;
      for(Int_t k=NbinsX;k>0;k--) {
	integral += hE1D[i]->GetBinContent(k) * dx;
	hV1D->SetBinContent(k,integral);
      }

    }
  
  }

  // Now, combine the electric field components into the total |E|
  TH2F *hETotal2D = (TH2F*) hE2D[0]->Clone("hETotal2D");
  hETotal2D->Reset();
  {
    Int_t NbinsX = hE2D[0]->GetNbinsX();
    Int_t NbinsY = hE2D[0]->GetNbinsY();
    for(Int_t j=0;j<NbinsX;j++) {
      for(Int_t k=0;k<NbinsY;k++) {
	Double_t E1 = hE2D[0]->GetBinContent(j,k);
	Double_t E2 = hE2D[1]->GetBinContent(j,k);
	Double_t E3 = hE2D[2]->GetBinContent(j,k);
	Double_t E  = TMath::Sqrt(E1*E1+E2*E2+E3*E3);
    
	hETotal2D->SetBinContent(j,k,E);
      }
    }
  }
  hETotal2D->GetZaxis()->SetTitle("E [E_{0}]");
  
  TH1F *hETotal1D = (TH1F*) hE1D[0]->Clone("hETotal1D");
  hETotal1D->Reset();
  {
    Int_t NbinsX = hE1D[0]->GetNbinsX();
    for(Int_t j=1;j<=NbinsX;j++) {
      Double_t E1 = hE1D[0]->GetBinContent(j);
      Double_t E2 = hE1D[1]->GetBinContent(j);
      Double_t E3 = hE1D[2]->GetBinContent(j);
      Double_t E  = TMath::Sqrt(E1*E1+E2*E2+E3*E3);
    
      hETotal1D->SetBinContent(j,E);
    }
  } 
  hETotal1D->GetYaxis()->SetTitle("E [E_{0}]");

  // Tunning the Histograms
  // ---------------------
  
  // Chaning to user units: 
  // --------------------------
  
  if(opt.Contains("units") && n0) {
    
    for(Int_t i=0;i<Nfields;i++) {
      Int_t NbinsX = hE2D[i]->GetNbinsX();
      Float_t xMin = skindepth * hE2D[i]->GetXaxis()->GetXmin() / PUnits::um;
      Float_t xMax = skindepth * hE2D[i]->GetXaxis()->GetXmax() / PUnits::um;
      Int_t NbinsY = hE2D[i]->GetNbinsY();
      Float_t ymin = skindepth * hE2D[i]->GetYaxis()->GetXmin() / PUnits::um;
      Float_t ymax = skindepth * hE2D[i]->GetYaxis()->GetXmax() / PUnits::um;
      hE2D[i]->SetBins(NbinsX,xMin,xMax,NbinsY,ymin,ymax);
      hE1D[i]->SetBins(NbinsX,xMin,xMax);
            
      for(Int_t j=0;j<hE2D[i]->GetNbinsX();j++) {
	for(Int_t k=0;k<hE2D[i]->GetNbinsY();k++) {
	  hE2D[i]->SetBinContent(j,k, hE2D[i]->GetBinContent(j,k) * ( E0 / (PUnits::GV/PUnits::m) ) );
	}
	hE1D[i]->SetBinContent(j, hE1D[i]->GetBinContent(j) * ( E0 / (PUnits::GV/PUnits::m) ) );
      }
      
      if(pData->IsCyl())
	hE2D[i]->GetYaxis()->SetTitle("r [#mum]");      
      else
	hE2D[i]->GetYaxis()->SetTitle("y [#mum]");      
      
      if(opt.Contains("comov"))
	hE2D[i]->GetXaxis()->SetTitle("#zeta [#mum]");
      else
	hE2D[i]->GetXaxis()->SetTitle("z [#mum]");
      
      if(i==0)
	hE2D[i]->GetZaxis()->SetTitle("E_{z} [GV/m]");
      else if(i==1)
	hE2D[i]->GetZaxis()->SetTitle("E_{y} [GV/m]");
      else if(i==2)
	hE2D[i]->GetZaxis()->SetTitle("E_{x} [GV/m]");
      
      
      if(opt.Contains("comov"))
	hE1D[i]->GetXaxis()->SetTitle("#zeta [mm]");
      else
	hE1D[i]->GetXaxis()->SetTitle("z [mm]");
      
      if(i==0)
	hE1D[i]->GetYaxis()->SetTitle("E_{z} [GV/m]");
      else if(i==1)
	hE1D[i]->GetYaxis()->SetTitle("E_{y} [GV/m]");
      else if(i==2)
	hE1D[i]->GetYaxis()->SetTitle("E_{x} [GV/m]");
      
      
      if(i==0) {    
	hV2D->SetBins(NbinsX,xMin,xMax,NbinsY,ymin,ymax);
	hETotal2D->SetBins(NbinsX,xMin,xMax,NbinsY,ymin,ymax);

	hV1D->SetBins(NbinsX,xMin,xMax);
	hETotal1D->SetBins(NbinsX,xMin,xMax);

	for(Int_t j=0;j<NbinsX;j++) {
	  for(Int_t k=0;k<NbinsY;k++) {
	    hV2D->SetBinContent(j,k, hV2D->GetBinContent(j,k) * E0 * skindepth / (PUnits::MV));
	    hETotal2D->SetBinContent(j,k, hETotal2D->GetBinContent(j,k) * ( E0 / (PUnits::GV/PUnits::m) ) );
	  }
	  hV1D->SetBinContent(j, hV1D->GetBinContent(j) * ( E0 * skindepth / (PUnits::MV) ) );
	  hETotal1D->SetBinContent(j, hETotal1D->GetBinContent(j) * ( E0 / (PUnits::GV/PUnits::m) ) );
	}
	
	if(pData->IsCyl()) {
	  hV2D->GetYaxis()->SetTitle("r [#mum]");      
	  hETotal2D->GetYaxis()->SetTitle("r [#mum]");      
	} else {
	  hV2D->GetYaxis()->SetTitle("y [#mum]");      
	  hETotal2D->GetYaxis()->SetTitle("y [#mum]");      
	
	}
   
	if(opt.Contains("comov")) {
	  hV2D->GetXaxis()->SetTitle("#zeta [#mum]");
	  hV1D->GetXaxis()->SetTitle("#zeta [#mum]");
	  hETotal2D->GetXaxis()->SetTitle("#zeta [#mum]");
	  hETotal1D->GetXaxis()->SetTitle("#zeta [#mum]");
	} else {
	  hV2D->GetXaxis()->SetTitle("z [#mum]");
	  hV2D->GetXaxis()->SetTitle("z [#mum]");
	  hETotal2D->GetXaxis()->SetTitle("z [#mum]");
	  hETotal1D->GetXaxis()->SetTitle("z [#mum]");
	}
	 
	hV2D->GetZaxis()->SetTitle("#Psi [MV]");
	hV1D->GetYaxis()->SetTitle("#Psi [MV]");
	hETotal2D->GetZaxis()->SetTitle("E [GV/m]");
	hETotal1D->GetYaxis()->SetTitle("E [GV/m]");
      }
    }
  }


  // --------------------------------------------------- Vertical Zoom ------------
  
  Float_t yRange    = (hE2D[0]->GetYaxis()->GetXmax() - hE2D[0]->GetYaxis()->GetXmin())/zoom;
  Float_t midPoint = (hE2D[0]->GetYaxis()->GetXmax() + hE2D[0]->GetYaxis()->GetXmin())/2.;
  Float_t yMin = midPoint-yRange/2;
  Float_t yMax = midPoint+yRange/2;
  if(pData->IsCyl()) {
    yMin = pData->GetXMin(1);
    yMax = yRange;
  }

  for(Int_t i=0;i<Nfields;i++) {
    if(!hE2D[i]) continue;
    hE2D[i]->GetYaxis()->SetRangeUser(yMin,yMax);
  }

  hETotal2D->GetYaxis()->SetRangeUser(yMin,yMax);

  Float_t xMin = hE2D[0]->GetXaxis()->GetXmin();
  Float_t xMax = hE2D[0]->GetXaxis()->GetXmax();
  Float_t xRange = xMax - xMin;
  
  // Change the range of z axis for the fields to be symmetric.
  Float_t *Emax = new Float_t[Nfields];
  Float_t *Emin = new Float_t[Nfields];
  for(Int_t i=0;i<Nfields;i++) {
    Emax[i] = hE2D[i]->GetMaximum();
    Emin[i] = hE2D[i]->GetMinimum();
    if(Emax[i] > TMath::Abs(Emin[i]))
      Emin[i] = -Emax[i];
    else
      Emax[i] = -Emin[i];
    hE2D[i]->GetZaxis()->SetRangeUser(Emin[i],Emax[i]); 
    //hE2D[i]->GetZaxis()->SetRangeUser(Emin[0],Emax[0]); 
  }
  
  // Potential 
  if(opt.Contains("units")) {
    trapPotential *=  ( E0 * skindepth / (PUnits::MV) ); 
  }  
  
  Float_t Vmin = hV1D->GetMinimum();    
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
  Float_t Vmax = hV1D->GetMaximum();    
   
  // // Dynamic potential palette
  // const Int_t potPNRGBs = 5;
  // const Int_t potPNCont = 64;
  // Float_t zeroPos = -Vmin/(Vmax-Vmin);

  // Double_t potPStops[potPNRGBs] = { 0.00, zeroPos-3.0/potPNCont,zeroPos, zeroPos+3.0/potPNCont, 1.00 };
  // Double_t potPRed[potPNRGBs]   = { 0.518, 0.965, 0.90, 0.498, 0.106 };
  // Double_t potPGreen[potPNRGBs] = { 0.078, 0.925, 0.90, 0.718, 0.078 };
  // Double_t potPBlue[potPNRGBs]  = { 0.106, 0.353, 0.90, 0.780, 0.518 };
   
  // PPalette * potentialPalette = (PPalette*) gROOT->FindObject("rbow2inv");
  // potentialPalette->CreateGradientColorTable(potPNRGBs, potPStops, 
  // 					     potPRed, potPGreen, potPBlue, potPNCont);
  
  // // Extract contours
  // TCanvas* c = new TCanvas("c","Contour List",0,0,600,600);
  // c->cd();
  
  // // Potential
  // TH2F *hV2Dc = (TH2F*) hV2D->Clone("hV2Dc");
  // const Int_t Ncontours = 25;
  // Double_t contours[Ncontours];
  // for(Int_t i=0; i<Ncontours; i++) {
  //   contours[i] = i*(trapPotential/5.0) - trapPotential; 
  // }
  // hV2Dc->SetContour(Ncontours, contours);
  // hV2Dc->Draw("cont list");
  
  // c->Update();
  // TObjArray *contsV2D = (TObjArray*) gROOT->GetListOfSpecials()->FindObject("contours");
  // TClonesArray graphsV2D("TGraph",Ncontours);
  // {
  //   Int_t ncontours = contsV2D->GetSize();
  //   TList* clist = NULL;
  //   Int_t nGraphs = 0;
  //   TGraph *gr = NULL;
  //   for(Int_t i = 0; i < ncontours; i++){
  //     if(i==0) continue;
      
  //     clist = (TList*) contsV2D->At(i);
    
  //     for(Int_t j = 0 ; j < clist->GetSize(); j++) {
  // 	gr = (TGraph*) clist->At(j);
  // 	if(!gr) continue;
      
  // 	gr->SetLineWidth(1);
  // 	gr->SetLineColor(kGray+1);
      
  // 	if( !((i)%5) ) {
  // 	  gr->SetLineWidth(2);
  // 	  gr->SetLineColor(kGray+2);
  // 	} 
  //     	new(graphsV2D[nGraphs]) TGraph(*gr) ;
  // 	nGraphs++;
  //     }
  //   }
  // }

  // "Axis range" in Osiris units:
  Double_t ylow  = hE2D[1]->GetYaxis()->GetBinLowEdge(FirstyBin);
  Double_t yup = hE2D[1]->GetYaxis()->GetBinUpEdge(LastyBin);
  Double_t xmin = hE2D[1]->GetXaxis()->GetXmin();
  Double_t xmax = hE2D[1]->GetXaxis()->GetXmax();

  TLine *lineYzero = new TLine(xmin,0.0,xmax,0.0);
  lineYzero->SetLineColor(kGray+2);
  lineYzero->SetLineStyle(2);

  TLine *lineYup = new TLine(xmin,yup,xmax,yup);
  lineYup->SetLineColor(kGray+1);
  lineYup->SetLineStyle(2);
 
  TLine *lineYdown = new TLine(xmin,ylow,xmax,ylow);
  lineYdown->SetLineColor(kGray+1);
  lineYdown->SetLineStyle(2);

  zStartPlasma -= shiftz; 
  zStartNeutral -= shiftz; 
  zEndNeutral -= shiftz; 
  
  if(opt.Contains("units")) {
    zStartPlasma *= skindepth / PUnits::um;
    zStartNeutral *= skindepth / PUnits::um;
    zEndNeutral *= skindepth / PUnits::um;
  }

  //  cout << "Start plasma = " << zStartPlasma << endl;
  TLine *lineStartPlasma = new TLine(zStartPlasma,yMin,zStartPlasma,yMax);
  lineStartPlasma->SetLineColor(kGray+2);
  lineStartPlasma->SetLineStyle(2);
  lineStartPlasma->SetLineWidth(3);

  //  cout << "Start plasma = " << zStartNeutral << endl;
  TLine *lineStartNeutral = new TLine(zStartNeutral,yMin,zStartNeutral,yMax);
  lineStartNeutral->SetLineColor(kGray+1);
  lineStartNeutral->SetLineStyle(2);
  lineStartNeutral->SetLineWidth(3);

  //  cout << "End plasma = " << zEndNeutral << endl;
  TLine *lineEndNeutral = new TLine(zEndNeutral,yMin,zEndNeutral,yMax);
  lineEndNeutral->SetLineColor(kGray+1);
  lineEndNeutral->SetLineStyle(2);
  lineEndNeutral->SetLineWidth(3);


  // Plotting
  // -----------------------------------------------

  // Canvas setup
  TCanvas *C = new TCanvas("C","2D Charge density and Electric field",858,1000);
 
  // Palettes setup
  TExec *exField  = new TExec("exField","rbow2Palette->cd();");
  TExec *exFieldT = new TExec("exFieldT","redPalette->cd();");
  TExec *exIonP   = new TExec("exIonP","redPalette->cd();");
  TExec *exPot    = new TExec("exPot","rbow2invPalette->cd();");
     
  // Actual Plotting!
  // ------------------------------------------------------------

  // Output file
  TString fOutName = Form("./%s/Plots/EFields2D/EFields2D",pData->GetPath().c_str());
  fOutName += Form("-%s_%i",pData->GetName(),time);

  // Setup Pad layout:
  Float_t lMargin = 0.15;
  Float_t rMargin = 0.18;
  Float_t bMargin = 0.10;
  Float_t tMargin = 0.04;
  Float_t factor = 1.0;    
  PlasmaGlob::CanvasAsymPartition(C,3,lMargin,rMargin,bMargin,tMargin,factor);
  
  TPad *pad[3];
  TString sLabels[] = {"(a)","(b)","(c)"};
  // Text objects
  TPaveText **textLabel = new TPaveText*[3];

  C->cd(0);
  char pname[16];
  sprintf(pname,"pad_%i",2);
  pad[0] = (TPad*) gROOT->FindObject(pname);
  pad[0]->Draw();
  pad[0]->cd(); // <---------------------------------------------- Top Plot ---------
  // if(opt.Contains("logz")) {
  //   pad[0]->SetLogz(1);
  // } else {
  //   pad[0]->SetLogz(0);
  // }
  pad[0]->SetFrameLineWidth(3);  
  pad[0]->SetTickx(1);

  TH2F *hFrame = (TH2F*) gROOT->FindObject("hFrame1");
  if(hFrame) delete hFrame;
  hFrame = (TH2F*) hE2D[0]->Clone("hFrame1");
  hFrame->Reset();
  
  hFrame->SetLabelFont(42,"xyz");
  hFrame->SetTitleFont(42,"xyz");
  
  hFrame->GetYaxis()->SetNdivisions(505);
  hFrame->GetYaxis()->SetLabelSize(0.085);
  hFrame->GetYaxis()->SetTitleSize(0.09);
  hFrame->GetYaxis()->SetTitleOffset(0.7);
  hFrame->GetYaxis()->SetTickLength(0.02);
  
  hFrame->GetXaxis()->SetLabelOffset(999.);
  hFrame->GetXaxis()->SetTitleOffset(999.);
  hFrame->GetXaxis()->SetTickLength(0.04);
  
  // Frame asymmetry:
  hFrame->Draw("col");

    //  hE2D[0]->GetZaxis()->SetNdivisions(505);
  hE2D[0]->GetZaxis()->SetTitleFont(42);
  hE2D[0]->GetZaxis()->SetTickLength(0.02);
  
  exField->Draw();
  hE2D[0]->Draw("colz same");

  lineYzero->Draw();
  if(opt.Contains("1dline")) {
    lineYzero->Draw();
    lineYdown->Draw();
    lineYup->Draw();
  }

  if(opt.Contains("sline")) {
    if(zStartPlasma>xmin && zStartPlasma<xmax)
      lineStartPlasma->Draw();
    if(zStartNeutral>xmin && zStartNeutral<xmax)
      lineStartNeutral->Draw();
    if(zEndNeutral>xmin && zEndNeutral<xmax)
      lineEndNeutral->Draw();
  }
  
  // Fit the V1D in the E2D pad:
  // Float_t rightmin = Vmin;
  // Float_t rightmax = Vmax;
  Float_t rightmin = Vmin;
  Float_t rightmax = Vmax;
  Float_t slope = (yMax/2)/rightmax;
  
  for(Int_t j=0;j<hV1D->GetNbinsX();j++) {
    hV1D->SetBinContent(j+1,hV1D->GetBinContent(j+1)*slope);
  }
  hV1D->SetLineStyle(1);
  hV1D->SetLineWidth(2);
  hV1D->SetLineColor(PlasmaGlob::elecLine);

  hV1D->Draw("sameC");
  
  // Line for trapping potential
  TLine *lineTrap = new TLine(hV1D->GetXaxis()->GetXmin(),
			      -(trapPotential)*slope,
			      hV1D->GetXaxis()->GetXmax(),
			      -(trapPotential)*slope);
  lineTrap->SetLineColor(PlasmaGlob::elecLine);
  lineTrap->SetLineStyle(2);
  lineTrap->SetLineWidth(1);
  lineTrap->Draw();
  
  // Fit the E1D in the E2D pad:
  rightmin = Emin[0];
  rightmax = Emax[0];
  slope = yMax/rightmax;
  
  for(Int_t j=0;j<hE1D[0]->GetNbinsX();j++) {
    hE1D[0]->SetBinContent(j+1,hE1D[0]->GetBinContent(j+1)*slope);
  }
  hE1D[0]->SetLineStyle(1);
  hE1D[0]->SetLineWidth(2);
  hE1D[0]->SetLineColor(kOrange+10);

  hE1D[0]->Draw("sameC");
  
  //lineYdown->Draw();
  //lineYup->Draw();

  pad[0]->Update();
  
  Float_t y1 = pad[0]->GetBottomMargin();
  Float_t y2 = 1 - pad[0]->GetTopMargin();
  Float_t x1 = pad[0]->GetLeftMargin();
  Float_t x2 = 1 - pad[0]->GetRightMargin();
  
  TPaletteAxis *palette = (TPaletteAxis*) hE2D[0]->GetListOfFunctions()->FindObject("palette");  
  if(palette) {
    palette->SetY2NDC(y2 - 0.00);
    palette->SetY1NDC(y1 + 0.00);
    palette->SetX1NDC(x2 + 0.005);
    palette->SetX2NDC(x2 + 0.03);
    // palette->SetTitleFont(42);
    palette->SetTitleSize(0.075);
    palette->SetTitleOffset(0.80);
    palette->SetLabelSize(0.075);
    palette->SetLabelFont(42);
    palette->SetLabelOffset(0.01);
    palette->SetBorderSize(2);
    palette->SetLineColor(1);
  }
   
  pad[0]->RedrawAxis(); 


  
  TPaveText *textTime = new TPaveText(xMax - 0.3*xRange, yMax-0.15*yRange, xMax-0.1, yMax-0.05*yRange);
  //x2-0.17,y2-0.12,x2-0.02,y2-0.02,"NDC");
  PlasmaGlob::SetPaveTextStyle(textTime,32); 
  char ctext[128];
  if(opt.Contains("units") && n0) 
    sprintf(ctext,"z = %5.1f #mum", Time * skindepth / PUnits::um);
  else
    sprintf(ctext,"t = %5.1f #omega_{p}^{-1}",Time);
  textTime->SetTextFont(42);
  textTime->AddText(ctext);

  textTime->Draw();
  // textDen->Draw();
  // if(opt.Contains("units"))
  //   textWav->Draw();

  textLabel[0] = new TPaveText(xMin + 0.02*xRange, yMax-0.2*yRange, xMin+0.30*xRange, yMax-0.05*yRange);
  PlasmaGlob::SetPaveTextStyle(textLabel[0],12); 
  textLabel[0]->SetTextFont(42);
  textLabel[0]->AddText(sLabels[0]);
  textLabel[0]->Draw();
  
  
  pad[0]->RedrawAxis(); 
 
  C->cd(0);
  sprintf(pname,"pad_%i",1);
  pad[1] = (TPad*) gROOT->FindObject(pname);
  pad[1]->Draw();
  pad[1]->cd(); // <---------------------------------------------------------- Middle Plot
  // if(opt.Contains("logz")) {
  //   pad[1]->SetLogz(1);
  // } else {
  //   pad[1]->SetLogz(0);
  // }
  pad[1]->SetFrameLineWidth(3);  
  pad[1]->SetTickx(1);
   
  hFrame = (TH2F*) gROOT->FindObject("hFrame2");
  if(hFrame) delete hFrame;
  hFrame = (TH2F*) hE2D[0]->Clone("hFrame2");
  hFrame->Reset();
  
  Float_t yFactor = pad[0]->GetAbsHNDC()/pad[1]->GetAbsHNDC();

  hFrame->GetYaxis()->SetLabelSize(yFactor*0.085);
  hFrame->GetYaxis()->SetTitleSize(yFactor*0.09);
  hFrame->GetYaxis()->SetTitleOffset(0.7/yFactor);
  hFrame->GetYaxis()->SetTickLength(0.02/yFactor);

  hFrame->GetXaxis()->SetLabelOffset(999.);
  hFrame->GetXaxis()->SetTitleOffset(999.);
  hFrame->GetXaxis()->SetTickLength(0.04*yFactor);
  
  hFrame->SetLabelFont(42,"xyz");
  hFrame->SetTitleFont(42,"xyz");

  hFrame->Draw("col");

    //  hE2D[1]->GetZaxis()->SetNdivisions(505);
  hE2D[1]->GetZaxis()->SetTitleFont(42);
  hE2D[1]->GetZaxis()->SetTickLength(0.02/yFactor);

  exField->Draw();
  hE2D[1]->Draw("colz same");

  if(opt.Contains("1dline")) {
    lineYzero->Draw();
    lineYdown->Draw();
    lineYup->Draw();
  }

  if(opt.Contains("sline")) {
    if(zStartPlasma>xmin && zStartPlasma<xmax)
      lineStartPlasma->Draw();
    if(zStartNeutral>xmin && zStartNeutral<xmax)
      lineStartNeutral->Draw();
    if(zEndNeutral>xmin && zEndNeutral<xmax)
      lineEndNeutral->Draw();
  }
  
  pad[1]->Update();
  
  y1 = pad[1]->GetBottomMargin();
  y2 = 1 - pad[1]->GetTopMargin();
  x1 = pad[1]->GetLeftMargin();
  x2 = 1 - pad[1]->GetRightMargin();
  
  palette = (TPaletteAxis*) hE2D[1]->GetListOfFunctions()->FindObject("palette");  
  if(palette) {
    palette->SetY2NDC(y2 - 0.00);
    palette->SetY1NDC(y1 + 0.00);
    palette->SetX1NDC(x2 + 0.005);
    palette->SetX2NDC(x2 + 0.03);
    // palette->SetTitleFont(42);
    palette->SetTitleSize(yFactor*0.075);
    palette->SetTitleOffset(0.80/yFactor);
    palette->SetLabelSize(yFactor*0.075);
    palette->SetLabelFont(42);
    palette->SetLabelOffset(0.01/yFactor);
    palette->SetBorderSize(2);
    palette->SetLineColor(1);
  }
   
  pad[1]->RedrawAxis(); 

  
  if(opt.Contains("1dline")) {
    lineYzero->Draw();
  }

  if(opt.Contains("sline")) {
    if(zStartPlasma>xmin && zStartPlasma<xmax)
      lineStartPlasma->Draw();
    if(zStartNeutral>xmin && zStartNeutral<xmax)
      lineStartNeutral->Draw();
    if(zEndNeutral>xmin && zEndNeutral<xmax)
      lineEndNeutral->Draw();
  }
  
 
  pad[1]->Update();
  
  y1 = pad[1]->GetBottomMargin();
  y2 = 1 - pad[1]->GetTopMargin();
  x1 = pad[1]->GetLeftMargin();
  x2 = 1 - pad[1]->GetRightMargin();
  
  palette = (TPaletteAxis*) hE2D[1]->GetListOfFunctions()->FindObject("palette");  
  if(palette) {
    palette->SetY2NDC(y2 - 0.00);
    palette->SetY1NDC(y1 + 0.00);
    palette->SetX1NDC(x2 + 0.005);
    palette->SetX2NDC(x2 + 0.03);
    // palette->SetTitleFont(42);
    palette->SetTitleSize(yFactor*0.075);
    palette->SetTitleOffset(0.80/yFactor);
    palette->SetLabelSize(yFactor*0.075);
    palette->SetLabelFont(42);
    palette->SetLabelOffset(0.01/yFactor);
    palette->SetBorderSize(2);
    palette->SetLineColor(1);
  }
   
  pad[1]->RedrawAxis(); 

  textLabel[1] = new TPaveText(xMin + 0.02*xRange, yMax-0.2*yRange, xMin+0.30*xRange, yMax-0.05*yRange);
  PlasmaGlob::SetPaveTextStyle(textLabel[1],12); 
  textLabel[1]->SetTextFont(42);
  textLabel[1]->AddText(sLabels[1]);
  textLabel[1]->Draw();


  C->cd(0);
  sprintf(pname,"pad_%i",0);
  pad[2] = (TPad*) gROOT->FindObject(pname);
  pad[2]->Draw();
  pad[2]->cd(); // <--------------------------------------------------------- Bottom Plot
  pad[2]->SetFrameLineWidth(3);  
  pad[2]->SetTickx(1);

  hFrame = (TH2F*) gROOT->FindObject("hFrame3");
  if(hFrame) delete hFrame;
  hFrame = (TH2F*) hE2D[0]->Clone("hFrame3");
  hFrame->Reset();

  yFactor = pad[0]->GetAbsHNDC()/pad[2]->GetAbsHNDC();
  
  
  hFrame->GetYaxis()->SetLabelSize(yFactor*0.085);
  hFrame->GetYaxis()->SetTitleSize(yFactor*0.09);
  hFrame->GetYaxis()->SetTitleOffset(0.7/yFactor);
  hFrame->GetYaxis()->SetTickLength(0.02/yFactor);

  hFrame->GetXaxis()->SetTitleSize(0.10);
  hFrame->GetXaxis()->SetLabelSize(0.08);
  hFrame->GetXaxis()->SetLabelOffset(0.02);
  hFrame->GetXaxis()->SetTitleOffset(1.0);
  hFrame->GetXaxis()->SetTickLength(0.04*yFactor);
  
  hFrame->SetLabelFont(42,"xyz");
  hFrame->SetTitleFont(42,"xyz");
  
  hFrame->Draw("col");

  //  hE2D[2]->GetZaxis()->SetNdivisions(505);
  hE2D[2]->GetZaxis()->SetTitleFont(42);
  hE2D[2]->GetZaxis()->SetTickLength(0.02/yFactor);

  exField->Draw();
  hE2D[2]->Draw("colz same");


  if(opt.Contains("1dline")) {
    lineYzero->Draw();
    lineYdown->Draw();
    lineYup->Draw();
  }

  if(opt.Contains("sline")) {
    if(zStartPlasma>xmin && zStartPlasma<xmax)
      lineStartPlasma->Draw();
    if(zStartNeutral>xmin && zStartNeutral<xmax)
      lineStartNeutral->Draw();
    if(zEndNeutral>xmin && zEndNeutral<xmax)
      lineEndNeutral->Draw();
  }
  
    
  pad[2]->Update();

  y1 = pad[2]->GetBottomMargin();
  y2 = 1 - pad[2]->GetTopMargin();
  x1 = pad[2]->GetLeftMargin();
  x2 = 1 - pad[2]->GetRightMargin();
 
  palette = (TPaletteAxis*)hE2D[2]->GetListOfFunctions()->FindObject("palette");
  if(palette) {
    palette->SetY2NDC(y2 - 0.00);
    palette->SetY1NDC(y1 + 0.00);
    palette->SetX1NDC(x2 + 0.005);
    palette->SetX2NDC(x2 + 0.03);
    // palette->SetTitleFont(42);
    palette->SetTitleSize(yFactor*0.075);
    palette->SetTitleOffset(0.80/yFactor);
    palette->SetLabelSize(yFactor*0.075);
    palette->SetLabelFont(42);
    palette->SetLabelOffset(0.01/yFactor);
    palette->SetBorderSize(2);
    palette->SetLineColor(1);
  }



  pad[2]->RedrawAxis(); 

  textLabel[2] = new TPaveText(xMin + 0.02*xRange, yMax-0.2*yRange, xMin+0.30*xRange, yMax-0.05*yRange);
  PlasmaGlob::SetPaveTextStyle(textLabel[2],12); 
  textLabel[2]->SetTextFont(42);
  textLabel[2]->AddText(sLabels[2]);
  textLabel[2]->Draw();

  C->cd();

  // Print to a file
  PlasmaGlob::imgconv(C,fOutName,opt);
  // ---------------------------------------------------------

  PlasmaGlob::DestroyCanvases();
}
  
