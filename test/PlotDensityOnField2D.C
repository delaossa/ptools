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

#include "PData.hh"
#include "PlasmaGlob.hh"
#include "PPalette.hh"

void PlotDensityOnField2D( const TString &sim, Int_t time, Int_t zoom=2, const TString &opt="") {
  
  PlasmaGlob::Initialize();

  // Palettes!
  gROOT->Macro("PlasmaPalettes.C");

  // More makeup
  gStyle->SetPadLeftMargin(0.08);    // Margin left axis 
  gStyle->SetPadRightMargin(0.12);   // Margin for palettes in 2D histos
  gStyle->SetTitleSize(0.06, "x");
  gStyle->SetTitleSize(0.06, "y");
  gStyle->SetTitleOffset(1.0,"x");
  gStyle->SetTitleOffset(0.6,"y");
  gStyle->SetPadGridY(1);
  gStyle->SetPadGridX(0);

  // Load PData
  PData *pData = PData::Get(sim.Data());
  pData->LoadFileNames(time);
  if(!pData->IsInit()) return;

  Bool_t ThreeD = kFALSE;
  if(sim.Contains("3D")) ThreeD = kTRUE; 
  
  // Get charge density histos
  Int_t Nspecies = pData->NSpecies();
  TH2F **hDen2D = new TH2F*[Nspecies];
  for(Int_t i=0;i<Nspecies;i++) {
    hDen2D[i] = NULL;

    if(!pData->GetChargeFileName(i)) 
      continue;

    if(pData->GetSpeciesName(i).find("beam")==string::npos)
      continue;

    if(!ThreeD)
      hDen2D[i] = pData->GetCharge(i);
    else
      hDen2D[i] = pData->GetCharge2DSliceZY(i);
    
    char hName[24];
    sprintf(hName,"hDen_%i",i);
    hDen2D[i]->SetName(hName);
    hDen2D[i]->GetXaxis()->CenterTitle();
    hDen2D[i]->GetYaxis()->CenterTitle();
    hDen2D[i]->GetXaxis()->SetTitle("z [c/#omega_{p}]");
    hDen2D[i]->GetYaxis()->SetTitle("y [c/#omega_{p}]");
    
    { // Center the y coordinate
      Float_t midY = 0.5 * (hDen2D[i]->GetYaxis()->GetXmin()+hDen2D[i]->GetYaxis()->GetXmax());
      Int_t NbinsX = hDen2D[i]->GetNbinsX();
      Float_t xMin = hDen2D[i]->GetXaxis()->GetXmin();
      Float_t xMax = hDen2D[i]->GetXaxis()->GetXmax();
      Int_t NbinsY = hDen2D[i]->GetNbinsY();
      Float_t yMin = hDen2D[i]->GetYaxis()->GetXmin() - midY;
      Float_t yMax = hDen2D[i]->GetYaxis()->GetXmax() - midY;
      hDen2D[i]->SetBins(NbinsX,xMin,xMax,NbinsY,yMin,yMax);
    }


    // Change to co-moving coordinate
    if(opt.Contains("comov")) {
      Int_t NbinsX = hDen2D[i]->GetNbinsX();
      Float_t xMin = hDen2D[i]->GetXaxis()->GetXmin()-hDen2D[i]->GetXaxis()->GetXmax();
      Float_t xMax = hDen2D[i]->GetXaxis()->GetXmax()-hDen2D[i]->GetXaxis()->GetXmax();
      Int_t NbinsY = hDen2D[i]->GetNbinsY();
      Float_t yMin = hDen2D[i]->GetYaxis()->GetXmin();
      Float_t yMax = hDen2D[i]->GetYaxis()->GetXmax();
      hDen2D[i]->SetBins(NbinsX,xMin,xMax,NbinsY,yMin,yMax);
      hDen2D[i]->GetXaxis()->SetTitle("#zeta [c/#omega_{p}]");
    }
    
    // Chaning to user units: 
    if(opt.Contains("units") && pData->GetPlasmaDensity()) {
      Int_t NbinsX = hDen2D[i]->GetNbinsX();
      Float_t xMin = 1e3 * pData->GetPlasmaSkinDepth() * hDen2D[i]->GetXaxis()->GetXmin();
      Float_t xMax = 1e3 * pData->GetPlasmaSkinDepth() * hDen2D[i]->GetXaxis()->GetXmax();
      Int_t NbinsY = hDen2D[i]->GetNbinsY();
      Float_t yMin = 1e3 * pData->GetPlasmaSkinDepth() * hDen2D[i]->GetYaxis()->GetXmin();
      Float_t yMax = 1e3 * pData->GetPlasmaSkinDepth() * hDen2D[i]->GetYaxis()->GetXmax();
      hDen2D[i]->SetBins(NbinsX,xMin,xMax,NbinsY,yMin,yMax);
      
      for(Int_t j=0;j<hDen2D[i]->GetNbinsX();j++) {
	for(Int_t k=0;k<hDen2D[i]->GetNbinsY();k++) {
	  hDen2D[i]->SetBinContent(j,k,1e-15 * 1e-6 * pData->GetPlasmaDensity() * hDen2D[i]->GetBinContent(j,k));
	}
      }

      hDen2D[i]->GetYaxis()->SetTitle("y [mm]");      
      if(opt.Contains("comov"))
	hDen2D[i]->GetXaxis()->SetTitle("#zeta [mm]");
      else
	hDen2D[i]->GetXaxis()->SetTitle("z [mm]");
      
    }
  }
  
  // Get electric fields
  const Int_t Nfields = 2;
  TH2F **hE2D = new TH2F*[Nfields];
  for(Int_t i=0;i<Nfields;i++) {
    hE2D[i] = NULL;

    if(!pData->GetEfieldFileName(i))
      continue;
    
    if(!ThreeD)
      hE2D[i] = pData->GetEField(i);
    else
      hE2D[i] = pData->GetEField2DSliceZY(i);
    
    char hName[24];
    sprintf(hName,"hE2D_%i",i);
    hE2D[i]->SetName(hName);   
    hE2D[i]->GetXaxis()->CenterTitle();
    hE2D[i]->GetYaxis()->CenterTitle();
    hE2D[i]->GetXaxis()->SetTitle("z [c/#omega_{p}]");
    hE2D[i]->GetYaxis()->SetTitle("y [c/#omega_{p}]");
    
    { // Center the y coordinate
      Float_t midY = 0.5 * (hE2D[i]->GetYaxis()->GetXmin()+hE2D[i]->GetYaxis()->GetXmax());
      Int_t NbinsX = hE2D[i]->GetNbinsX();
      Float_t xMin = hE2D[i]->GetXaxis()->GetXmin();
      Float_t xMax = hE2D[i]->GetXaxis()->GetXmax();
      Int_t NbinsY = hE2D[i]->GetNbinsY();
      Float_t yMin = hE2D[i]->GetYaxis()->GetXmin() - midY;
      Float_t yMax = hE2D[i]->GetYaxis()->GetXmax() - midY;
      hE2D[i]->SetBins(NbinsX,xMin,xMax,NbinsY,yMin,yMax);
    }

    // Change to co-moving coordinate
    if(opt.Contains("comov")) {
      Int_t NbinsX = hE2D[i]->GetNbinsX();
      Float_t xMin = hE2D[i]->GetXaxis()->GetXmin()-hE2D[i]->GetXaxis()->GetXmax();
      Float_t xMax = hE2D[i]->GetXaxis()->GetXmax()-hE2D[i]->GetXaxis()->GetXmax();
      Int_t NbinsY = hE2D[i]->GetNbinsY();
      Float_t yMin = hE2D[i]->GetYaxis()->GetXmin();
      Float_t yMax = hE2D[i]->GetYaxis()->GetXmax();
      hE2D[i]->SetBins(NbinsX,xMin,xMax,NbinsY,yMin,yMax);
      hE2D[i]->GetXaxis()->SetTitle("#zeta [c/#omega_{p}]");
    }
  
    // Chaning to user units: 
    if(opt.Contains("units") && pData->GetPlasmaDensity()) {
      Int_t NbinsX = hE2D[i]->GetNbinsX();
      Float_t xMin = 1e3 * pData->GetPlasmaSkinDepth() * hE2D[i]->GetXaxis()->GetXmin();
      Float_t xMax = 1e3 * pData->GetPlasmaSkinDepth() * hE2D[i]->GetXaxis()->GetXmax();
      Int_t NbinsY = hE2D[i]->GetNbinsY();
      Float_t yMin = 1e3 * pData->GetPlasmaSkinDepth() * hE2D[i]->GetYaxis()->GetXmin();
      Float_t yMax = 1e3 * pData->GetPlasmaSkinDepth() * hE2D[i]->GetYaxis()->GetXmax();
      hE2D[i]->SetBins(NbinsX,xMin,xMax,NbinsY,yMin,yMax);
      
      for(Int_t j=0;j<hE2D[i]->GetNbinsX();j++) {
	for(Int_t k=0;k<hE2D[i]->GetNbinsY();k++) {
	  hE2D[i]->SetBinContent(j,k,1e-9 * pData->GetPlasmaE0() * hE2D[i]->GetBinContent(j,k));
	}
      }

      hE2D[i]->GetYaxis()->SetTitle("y [mm]");	
      if(opt.Contains("comov"))
	hE2D[i]->GetXaxis()->SetTitle("#zeta [mm]");
      else
	hE2D[i]->GetXaxis()->SetTitle("z [mm]");
	
    }
    
  }
  
  // We need this histo!!
  if(!hDen2D[1])
    return;

  // Tunning the Histograms
  // ---------------------
  
  Float_t Time = pData->GetRealTime();
   
  // Zoom
  Float_t range    = (hDen2D[1]->GetYaxis()->GetXmax() - hDen2D[1]->GetYaxis()->GetXmin())/zoom;
  Float_t midPoint =  (hDen2D[1]->GetYaxis()->GetXmax() + hDen2D[1]->GetYaxis()->GetXmin())/2.;
  
  hDen2D[1]->GetYaxis()->SetRangeUser(midPoint-range/2,midPoint+range/2);  
  hDen2D[1]->GetZaxis()->SetRangeUser(0.002,1.05*hDen2D[1]->GetMaximum());
  
  hE2D[0]->GetYaxis()->SetRangeUser(midPoint-range/2,midPoint+range/2);  
  hE2D[1]->GetYaxis()->SetRangeUser(midPoint-range/2,midPoint+range/2);  
  
  // Plotting
  // -----------------------------------------------

  // Canvas setup
  TCanvas *C = new TCanvas("C","2D Charge density and Electric field",850,500);
 
  // Palettes setup
  TExec *exPlasma = new TExec("exPlasma","plasmaPalette->cd();"); 
  TExec *exElec   = new TExec("exElec","hot2Palette->cd();");
  TExec *exField  = new TExec("exPlasmaField","grayPalette->cd();");

  // Text objects
  TPaveText *text1 = new TPaveText(0.6,0.94,0.88,1.0,"NDC");
  PlasmaGlob::SetPaveTextStyle(text1);
  if(opt.Contains("units") && pData->GetPlasmaDensity()) 
    text1->AddText("dQ/dydz [10^{15}/cc]");
  else
    text1->AddText("dQ/dydz [a.u.]");

  TPaveText *text2 = new TPaveText(0.6,0.94,0.88,1.0,"NDC");
  PlasmaGlob::SetPaveTextStyle(text2);  
  if(opt.Contains("units") && pData->GetPlasmaDensity()) 
    text2->AddText("Longitudinal E field [GeV/m]");
  else
    text2->AddText("Longitudinal E field [E_{0}]");

  TPaveText *textTime = new TPaveText(0.7,0.82,0.85,0.88,"NDC");
  PlasmaGlob::SetPaveTextStyle(textTime); 
  char ctext[128];
  if(opt.Contains("units") && pData->GetPlasmaDensity()) 
    sprintf(ctext,"L = %5.1f mm", 1e3 * pData->GetPlasmaSkinDepth() * Time);
  else
    sprintf(ctext,"T = %5.1f 1/#omega_{p}",Time);
  textTime->AddText(ctext);
 

  
  // Actual Plotting!
  // ------------------------------------------------------------

  // Output file
  TString fOutName = Form("./%s/Plots/DensityOnField2D/DensityOnField2D",sim.Data());
  fOutName += Form("-%s_%i",sim.Data(),time);

  C->cd(0); // <--- Top Plot

  gPad->SetFrameLineWidth(2);  

  
  TH2F *hFrame = (TH2F*) gROOT->FindObject("hFrame1");
  if(hFrame) delete hFrame;
  hFrame = (TH2F*) hDen2D[1]->Clone("hFrame1");
  hFrame->Reset();
  hFrame->Draw("col");
  
  exField->Draw();
  hE2D[0]->Draw("col same");

  exElec->Draw();
  hDen2D[1]->Draw("colz same");
  
  text1->Draw();
  textTime->Draw();
 
  gPad->RedrawAxis(); 

  C->cd();

  // Print to a file
  PlasmaGlob::imgconv(C,fOutName,opt);
  // ---------------------------------------------------------

}
