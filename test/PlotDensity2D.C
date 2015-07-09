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
#include <TExec.h>
#include <TPaletteAxis.h>

#include "PData.hh"
#include "PlasmaGlob.hh"
#include "PPalette.hh"

void PlotDensity2D( const TString &sim, Int_t time, Int_t zoom=2, Int_t Nbins=1, const TString &options="") {
  
  PlasmaGlob::Initialize();

  // Palettes!
  gROOT->Macro("PlasmaPalettes.C");
  
  // Load PData
  PData *pData = PData::Get(sim.Data());
  pData->LoadFileNames(time);
  if(!pData->IsInit()) return;

  TString opt = options;
 
  // More makeup
  gStyle->SetPadTopMargin(0.07); 
  gStyle->SetPadLeftMargin(0.08);    // Margin left axis 
  gStyle->SetPadRightMargin(0.12);   // Margin for palettes in 2D histos
  gStyle->SetTitleSize(0.06, "x");
  gStyle->SetTitleSize(0.06, "y");
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetTitleOffset(0.6,"y");
  gStyle->SetTitleOffset(0.7,"z");
  gStyle->SetLabelSize(0.06, "x");
  gStyle->SetLabelSize(0.06, "y");
  gStyle->SetLabelOffset(0.012,"x");
  // gStyle->SetLabelOffset(0.6,"y");
  if(opt.Contains("grid")) {
    gStyle->SetPadGridX(1);
    gStyle->SetPadGridY(1);
  }

  Bool_t CYL = kFALSE;
  if(sim.Contains("cyl")) { CYL = kTRUE; opt += "cyl"; } 
    
  Bool_t ThreeD = kFALSE;
  if(sim.Contains("3D")) ThreeD = kTRUE; 
  
  // Some plasma constants
  Double_t n0 = pData->GetPlasmaDensity();
  Double_t kp = pData->GetPlasmaK();
  Double_t skindepth = 1/kp;
  Double_t E0 = pData->GetPlasmaE0();

  // Some beam properties:
  Float_t Ebeam = pData->GetBeamEnergy() * PUnits::MeV;
  Float_t gamma = Ebeam / PConst::ElectronMassE;
  Float_t vbeam = TMath::Sqrt(1 - 1/(gamma*gamma));
  
  // Time in OU
  Float_t Time = pData->GetRealTime();

  // z start of the plasma in normalized units.
  Float_t zStartPlasma = pData->GetPlasmaStart()*kp;

  // z start of the beam in normalized units.
  Float_t zStartBeam = pData->GetBeamStart()*kp;
  Time -= zStartPlasma - zStartBeam;
  

  // Get charge density histos
  Int_t Nspecies = pData->NSpecies();
  TH2F **hDen2D = new TH2F*[Nspecies];
  for(Int_t i=0;i<Nspecies;i++) {
    hDen2D[i] = NULL;

    if(!pData->GetChargeFileName(i)) 
      continue;
    
    if(!ThreeD)
      hDen2D[i] = pData->GetCharge(i,opt);
    else
      hDen2D[i] = pData->GetCharge2DSliceZY(i,-1,Nbins);
    
    char hName[24];
    sprintf(hName,"hDen_%i",i);
    hDen2D[i]->SetName(hName);
    hDen2D[i]->GetXaxis()->CenterTitle();
    hDen2D[i]->GetYaxis()->CenterTitle();
    hDen2D[i]->GetZaxis()->CenterTitle();
    hDen2D[i]->GetYaxis()->SetTitle("y [c/#omega_{p}]");
    if(opt.Contains("comov"))
      hDen2D[i]->GetXaxis()->SetTitle("#zeta [c/#omega_{p}]");
    else
      hDen2D[i]->GetXaxis()->SetTitle("z [c/#omega_{p}]");
    if(i==0)
      hDen2D[i]->GetZaxis()->SetTitle("#LTn_{e}#GT [n_{0}]");
    else
      hDen2D[i]->GetZaxis()->SetTitle("#LTn_{b}#GT [n_{0}]");
    
  }
  
  // Get electric fields
  const Int_t Nfields = 2;
  TH2F **hE2D = new TH2F*[Nfields];
  for(Int_t i=0;i<Nfields;i++) {
    hE2D[i] = NULL;

    if(!pData->GetEfieldFileName(i))
      continue;
    
    if(!ThreeD)
      hE2D[i] = pData->GetEField(i,opt);
    else
      hE2D[i] = pData->GetEField2DSliceZY(i,-1,Nbins);
    
    char hName[24];
    sprintf(hName,"hE2D_%i",i);
    hE2D[i]->SetName(hName);   
    hE2D[i]->GetXaxis()->CenterTitle();
    hE2D[i]->GetYaxis()->CenterTitle();
    hE2D[i]->GetZaxis()->CenterTitle();
    if(opt.Contains("comov"))
      hE2D[i]->GetXaxis()->SetTitle("#zeta [c/#omega_{p}]");
    else
      hE2D[i]->GetXaxis()->SetTitle("z [c/#omega_{p}]");
    hE2D[i]->GetYaxis()->SetTitle("y [c/#omega_{p}]");
    if(i==0)
      hE2D[i]->GetZaxis()->SetTitle("E_{z} [E_{0}]");
    else if(i==1)
      hE2D[i]->GetZaxis()->SetTitle("E_{y} [E_{0}]");
    else if(i==2)
      hE2D[i]->GetZaxis()->SetTitle("E_{x} [E_{0}]");
    
  }
  
  // Tunning the Histograms
  // ---------------------


  // Changing to user units: 
  // --------------------------

  if(opt.Contains("units") && n0) {
    
    for(Int_t i=0;i<Nspecies;i++) {
      Int_t NbinsX = hDen2D[i]->GetNbinsX();
      Float_t xMin = skindepth * hDen2D[i]->GetXaxis()->GetXmin() / PUnits::mm;
      Float_t xMax = skindepth * hDen2D[i]->GetXaxis()->GetXmax() / PUnits::mm;
      Int_t NbinsY = hDen2D[i]->GetNbinsY();
      Float_t yMin = skindepth * hDen2D[i]->GetYaxis()->GetXmin() / PUnits::mm;
      Float_t yMax = skindepth * hDen2D[i]->GetYaxis()->GetXmax() / PUnits::mm;
      hDen2D[i]->SetBins(NbinsX,xMin,xMax,NbinsY,yMin,yMax);
      for(Int_t j=0;j<hDen2D[i]->GetNbinsX();j++) {
	for(Int_t k=0;k<hDen2D[i]->GetNbinsY();k++) {
	  hDen2D[i]->SetBinContent(j,k, hDen2D[i]->GetBinContent(j,k) * n0 / (1e15/PUnits::cm3) );
	}
      }
      
      hDen2D[i]->GetYaxis()->SetTitle("y [mm]");      
      if(opt.Contains("comov"))
	hDen2D[i]->GetXaxis()->SetTitle("#zeta [mm]");
      else
	hDen2D[i]->GetXaxis()->SetTitle("z [mm]");
      
      if(i==0)
	hDen2D[i]->GetZaxis()->SetTitle("#LTn_{e}#GT [10^{15}/cm^{3}]");
      else
	hDen2D[i]->GetZaxis()->SetTitle("#LTn_{b}#GT [10^{15}/cm^{3}]"); 
    }

    for(Int_t i=0;i<Nfields;i++) {
      Int_t NbinsX = hE2D[i]->GetNbinsX();
      Float_t xMin = skindepth * hE2D[i]->GetXaxis()->GetXmin() / PUnits::mm;
      Float_t xMax = skindepth * hE2D[i]->GetXaxis()->GetXmax() / PUnits::mm;
      Int_t NbinsY = hE2D[i]->GetNbinsY();
      Float_t yMin = skindepth * hE2D[i]->GetYaxis()->GetXmin() / PUnits::mm;
      Float_t yMax = skindepth * hE2D[i]->GetYaxis()->GetXmax() / PUnits::mm;
      hE2D[i]->SetBins(NbinsX,xMin,xMax,NbinsY,yMin,yMax);
      for(Int_t j=0;j<hE2D[i]->GetNbinsX();j++) {
	for(Int_t k=0;k<hE2D[i]->GetNbinsY();k++) {
	  hE2D[i]->SetBinContent(j,k, hE2D[i]->GetBinContent(j,k) * ( E0 / (PUnits::GV/PUnits::m) ) );
	}
      }
      
      hE2D[i]->GetYaxis()->SetTitle("y [mm]");      
      if(opt.Contains("comov"))
	hE2D[i]->GetXaxis()->SetTitle("#zeta [mm]");
      else
	hE2D[i]->GetXaxis()->SetTitle("z [mm]");
      
      if(i==0)
	hE2D[i]->GetZaxis()->SetTitle("E_{z} [GV/m]");
      else if(i==1)
	hE2D[i]->GetZaxis()->SetTitle("E_{y} [GV/m]");
      else if(i==2)
	hE2D[i]->GetZaxis()->SetTitle("E_{x} [GV/m]");
    }
    
  }
  
  
  // Zoom !!
  Float_t range    = (hDen2D[0]->GetYaxis()->GetXmax() - hDen2D[0]->GetYaxis()->GetXmin())/zoom;
  Float_t midPoint = (hDen2D[0]->GetYaxis()->GetXmax() + hDen2D[0]->GetYaxis()->GetXmin())/2.;
  
  if(!CYL) {
    hDen2D[0]->GetYaxis()->SetRangeUser(midPoint-range/2,midPoint+range/2);  
    hDen2D[1]->GetYaxis()->SetRangeUser(midPoint-range/2,midPoint+range/2);  
    hE2D[0]->GetYaxis()->SetRangeUser(midPoint-range/2,midPoint+range/2);  
    hE2D[1]->GetYaxis()->SetRangeUser(midPoint-range/2,midPoint+range/2);  
  } else {
    hDen2D[0]->GetYaxis()->SetRangeUser(hDen2D[0]->GetYaxis()->GetXmin(),range);  
    hDen2D[1]->GetYaxis()->SetRangeUser(hDen2D[1]->GetYaxis()->GetXmin(),range);
    hE2D[0]->GetYaxis()->SetRangeUser(hE2D[0]->GetYaxis()->GetXmin(),range);
    hE2D[1]->GetYaxis()->SetRangeUser(hE2D[1]->GetYaxis()->GetXmin(),range);  
  }

  // Set the range of the plasma charge density histogram for maximum constrast 
  // using a dynamic palette wich adjust the nominal value to a certain color.
  Float_t density = 1;
  if(opt.Contains("units") && pData->GetPlasmaDensity())
    density = 1e-15 * 1e-6 * pData->GetPlasmaDensity();
  
  Float_t Max  = 1.1 * hDen2D[0]->GetMaximum();
  Float_t Base = density;
  Float_t Min  = 2.* Base - Max;
  if(Max >= 2. * Base) {
    Min = 0;
  } else if(Max<1.0 * Base) {
    Max = 1.1 * Base;
    Min = 0.;
  }
    
  if(hDen2D[0])
    hDen2D[0]->GetZaxis()->SetRangeUser(Min,Max); 
  if(hDen2D[1])
    hDen2D[1]->GetZaxis()->SetRangeUser(0.001,hDen2D[1]->GetMaximum());
    

  // Change the range of z axis for the fields to be symmetric.
  Float_t Emax = hE2D[0]->GetMaximum();
  Float_t Emin = hE2D[0]->GetMinimum();
  if(Emax > TMath::Abs(Emin))
    Emin = -Emax;
  else
    Emax = -Emin;
  hE2D[0]->GetZaxis()->SetRangeUser(Emin,Emax); 
  
  Emax = hE2D[1]->GetMaximum();
  Emin = hE2D[1]->GetMinimum();
  if(Emax > TMath::Abs(Emin))
    Emin = -Emax;
  else
    Emax = -Emin;
  hE2D[1]->GetZaxis()->SetRangeUser(Emin,Emax); 

  // Dynamic plasma palette
  const Int_t plasmaDNRGBs = 3;
  const Int_t plasmaDNCont = 64;
  Double_t plasmaDStops[plasmaDNRGBs] = { 0.00, 0.50, 1.00 };
  Double_t plasmaDRed[plasmaDNRGBs]   = { 0.99, 0.90, 0.00 };
  Double_t plasmaDGreen[plasmaDNRGBs] = { 0.99, 0.90, 0.00 };
  Double_t plasmaDBlue[plasmaDNRGBs]  = { 0.99, 0.90, 0.00 };
  Double_t basePos = (1.0/(Max-Min))*(Base - Min);
  plasmaDStops[1] = basePos;
 
  PPalette * plasmaPalette = (PPalette*) gROOT->FindObject("plasma");
  plasmaPalette->CreateGradientColorTable(plasmaDNRGBs, plasmaDStops, 
					  plasmaDRed, plasmaDGreen, plasmaDBlue, plasmaDNCont);
  // Plotting
  // -----------------------------------------------

  // Canvas setup
  TCanvas *C = new TCanvas("C","2D Charge density and Electric field",750,1000);
  C->Divide(1,3);

  TExec *exPlasma;
  TExec *exElec;
  if( (Max >= 2. * Base) || (Max<1.0 * Base) ) {
    exPlasma = new TExec("exPlasma","plasmaPalette->cd();"); 
    exElec   = new TExec("exElec","electronPalette->cd();");
  } else {
    exPlasma= new TExec("exPlasma","plasmaPalette->cd();");
    exElec   = new TExec("exElec","electronPalette->cd();");
  }
  TExec *exField  = new TExec("exField","rbowPalette->cd();");

  // Text objects
  TPaveText *text1 = new TPaveText(0.6,0.94,0.88,1.0,"NDC");
  PlasmaGlob::SetPaveTextStyle(text1);
  if(opt.Contains("units") && pData->GetPlasmaDensity()) 
    text1->AddText("Electron density [10^{15}/cm^{3}]");
  else
    text1->AddText("Electron density [n_{0}]");

  TPaveText *text2 = new TPaveText(0.6,0.94,0.88,1.0,"NDC");
  PlasmaGlob::SetPaveTextStyle(text2);  
  if(opt.Contains("units") && pData->GetPlasmaDensity()) 
    text2->AddText("Longitudinal E field [GV/m]");
  else
    text2->AddText("Longitudinal E field [E_{0}]");

  TPaveText *text3 = new TPaveText(0.6,0.94,0.88,1.0,"NDC");
  PlasmaGlob::SetPaveTextStyle(text3);  
  if(opt.Contains("units") && pData->GetPlasmaDensity()) 
    text3->AddText("Transverse E field [GV/m]");
  else
    text3->AddText("Transverse E field [E_{0}]");
 
  TPaveText *textTime = new TPaveText(0.7,0.82,0.85,0.88,"NDC");
  PlasmaGlob::SetPaveTextStyle(textTime); 
  char ctext[128];
  if(opt.Contains("units") && pData->GetPlasmaDensity()) 
    sprintf(ctext,"Z = %5.1f mm", 1e3 * skindepth * Time);
  else
    sprintf(ctext,"T = %5.1f #omega_{p}^{-1}",Time);
  textTime->AddText(ctext);
 
  TPaveText *textDen = new TPaveText(0.13,0.82,0.38,0.88,"NDC");
  PlasmaGlob::SetPaveTextStyle(textDen,12); 
  textDen->SetTextColor(kOrange+10);
   if(opt.Contains("units") && pData->GetPlasmaDensity())
    sprintf(ctext,"n_{0} = %5.2f x 10^{15} / cm^{3}", 1e-6 * 1e-15 * pData->GetPlasmaDensity());
  else if(pData->GetBeamDensity() && pData->GetPlasmaDensity())
    sprintf(ctext,"n_{b}/n_{0} = %5.2f", pData->GetBeamDensity()/pData->GetPlasmaDensity());
   textDen->AddText(ctext);
  
  // Actual Plotting!
  // ------------------------------------------------------------

  // Output file
  TString fOutName = Form("./%s/Plots/Density2D/Density2D",sim.Data());
  fOutName += Form("-%s_%i",sim.Data(),time);

  C->cd(1); // <--- Top Plot

  gPad->SetFrameLineWidth(3);  

  
  TH2F *hFrame = (TH2F*) gROOT->FindObject("hFrame1");
  if(hFrame) delete hFrame;
  hFrame = (TH2F*) hDen2D[0]->Clone("hFrame1");
  hFrame->Reset();
  hFrame->Draw("col");
  
  exPlasma->Draw();
  hDen2D[0]->Draw("colz same");
  
  if(hDen2D[1]) {
    exElec->Draw();
    hDen2D[1]->Draw("colz same");
  }
  
  if(Nspecies>=3) {
    if(hDen2D[2]) {
      exElec->Draw();
      hDen2D[2]->Draw("col same");
    }
  }

  gPad->Update();
  TPaletteAxis *palette = (TPaletteAxis*)hDen2D[1]->GetListOfFunctions()->FindObject("palette");
  
  Float_t y1 = gPad->GetBottomMargin();
  Float_t y2 = 1 - gPad->GetTopMargin();
  palette->SetY2NDC(y2);
  palette->SetY1NDC(0.5*(y1+y2));
  //palette->SetTitleOffset(0.8);
  palette->SetBorderSize(2);
  palette->SetLineColor(1);
  
  TPaletteAxis *palette2 = (TPaletteAxis*)hDen2D[0]->GetListOfFunctions()->FindObject("palette");
  palette2->SetY2NDC(0.5*(y1+y2));
  palette2->SetY1NDC(y1);

  text1->Draw();
  textTime->Draw();
  textDen->Draw();

  gPad->RedrawAxis(); 

  C->cd(2); // <--- Mid Plot
  gPad->SetFrameLineWidth(3);  
  
  hFrame = (TH2F*) gROOT->FindObject("hFrame2");
  if(hFrame) delete hFrame;
  hFrame = (TH2F*) hE2D[0]->Clone("hFrame2");
  hFrame->Reset();

  hFrame->Draw("col");
  exField->Draw();
  hE2D[0]->Draw("colz same");

  text2->Draw();
  textTime->Draw();

  gPad->RedrawAxis(); 

  C->cd(3); // <--- Bottom Plot
  gPad->SetFrameLineWidth(3); 
   
  hFrame = (TH2F*) gROOT->FindObject("hFrame3");
  if(hFrame) delete hFrame;
  hFrame = (TH2F*) hE2D[1]->Clone("hFrame3");
  hFrame->Reset();

  hFrame->Draw("col");
  exField->Draw();
  hE2D[1]->Draw("colz same");

  text3->Draw();
  textTime->Draw();

  gPad->RedrawAxis(); 

  C->cd();

  // Print to a file
  PlasmaGlob::imgconv(C,fOutName,opt);
  // ---------------------------------------------------------

}
