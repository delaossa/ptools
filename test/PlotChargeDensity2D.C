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

#include "PData.hh"
#include "PlasmaGlob.hh"
#include "PPalette.hh"

void PlotChargeDensity2D( const TString &sim, Int_t time, Int_t zoom=2, Int_t Nbins=1, const TString &opt="") {
  
  PlasmaGlob::Initialize();

  // Palettes!
  gROOT->Macro("PlasmaPalettes.C");

  gStyle->SetPadTopMargin(0.07); 
  gStyle->SetPadLeftMargin(0.10); 
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetTitleOffset(0.95,"z");
  gStyle->SetLabelSize(0.042,"z");
  gStyle->SetTitleOffset(1.3,"x");
  gStyle->SetLabelOffset(0.01,"x");
  gStyle->SetTitleOffset(0.95,"y");
  if(opt.Contains("grid")) {
    gStyle->SetPadGridX(1);
    gStyle->SetPadGridY(1);
  }

  // Load PData
  PData *pData = PData::Get(sim.Data());
  pData->LoadFileNames(time);
  if(!pData->IsInit()) return;

  Bool_t CYL = kFALSE;
  if(sim.Contains("cyl")) CYL = kTRUE; 
    
  Bool_t ThreeD = kFALSE;
  if(sim.Contains("3D")) ThreeD = kTRUE; 

  // Some plasma constants
  Double_t n0 = pData->GetPlasmaDensity();
  Double_t kp = pData->GetPlasmaK();
  Double_t skindepth = 1/kp;
  
  // Some beam properties:
  Float_t Ebeam = pData->GetBeamEnergy() * PUnits::MeV;
  Float_t gamma = Ebeam / PConst::ElectronMassE;
  Float_t vbeam = TMath::Sqrt(1 - 1/(gamma*gamma));
  // cout << Form(" - Bunch gamma      = %8.4f", gamma ) << endl;
  // cout << Form(" - Bunch velocity   = %8.4f c", vbeam ) << endl;

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
    
    char hName[24];
    sprintf(hName,"hDen_%i",i);
    hDen2D[i] = (TH2F*) gROOT->FindObject(hName);
    if(hDen2D[i]) delete hDen2D[i];

    if(!ThreeD)
      hDen2D[i] = pData->GetCharge(i);
    else
      hDen2D[i] = pData->GetCharge2DSliceZY(i,-1,Nbins);
    
  
    hDen2D[i]->SetName(hName);
    hDen2D[i]->GetXaxis()->CenterTitle();
    hDen2D[i]->GetYaxis()->CenterTitle();
    hDen2D[i]->GetZaxis()->CenterTitle();
    hDen2D[i]->GetXaxis()->SetTitle("z [c/#omega_{p}]");
    hDen2D[i]->GetYaxis()->SetTitle("y [c/#omega_{p}]");
    if(i==0)
      hDen2D[i]->GetZaxis()->SetTitle("#LTn_{e}#GT [n_{0}]");
    else
      hDen2D[i]->GetZaxis()->SetTitle("#LTn_{b}#GT [n_{0}]");

    if(opt.Contains("center")) { // Center the y and z coordinate
      Int_t NbinsX = hDen2D[i]->GetNbinsX();
      Float_t xMin = hDen2D[i]->GetXaxis()->GetXmin() - zStartPlasma;
      Float_t xMax = hDen2D[i]->GetXaxis()->GetXmax() - zStartPlasma;
      Int_t NbinsY = hDen2D[i]->GetNbinsY();
      Float_t yMin = hDen2D[i]->GetYaxis()->GetXmin();
      Float_t yMax = hDen2D[i]->GetYaxis()->GetXmax();
      
      if(!CYL) {
	Float_t midY = 0.5 * (hDen2D[i]->GetYaxis()->GetXmin()+hDen2D[i]->GetYaxis()->GetXmax());
	yMin -= midY;
	yMax -= midY;
      }

      hDen2D[i]->SetBins(NbinsX,xMin,xMax,NbinsY,yMin,yMax);
    }


    // Change to co-moving coordinate
    if(opt.Contains("comov")) {
      Int_t NbinsX = hDen2D[i]->GetNbinsX();
      Float_t xMin = hDen2D[i]->GetXaxis()->GetXmin() - vbeam*Time;
      Float_t xMax = hDen2D[i]->GetXaxis()->GetXmax() - vbeam*Time;
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
     
      if(i==0)
	hDen2D[i]->GetZaxis()->SetTitle("#LTn_{e}#GT [10^{15}/cc]");
      else
	hDen2D[i]->GetZaxis()->SetTitle("#LTn_{b}#GT [10^{15}/cc]"); 
    }
  }
  
  // Tunning the Histograms
  // ---------------------
  
  // Zoom
  Float_t range    = (hDen2D[0]->GetYaxis()->GetXmax() - hDen2D[0]->GetYaxis()->GetXmin())/zoom;
  Float_t midPoint = (hDen2D[0]->GetYaxis()->GetXmax() + hDen2D[0]->GetYaxis()->GetXmin())/2.;
  
  if(!CYL) {
    hDen2D[0]->GetYaxis()->SetRangeUser(midPoint-range/2,midPoint+range/2);  
  } else {
    hDen2D[0]->GetYaxis()->SetRangeUser(hDen2D[0]->GetYaxis()->GetXmin(),range);  
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
    
  // Dynamic plasma palette
  const Int_t plasmaDNRGBs = 3;
  const Int_t plasmaDNCont = 64;
  Double_t plasmaDStops[plasmaDNRGBs] = { 0.00,  0.5, 1.00 };
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
  TCanvas *C;
  if(opt.Contains("hres") && !opt.Contains("pdf")) // high resolution for plain grahics output.
    C = new TCanvas("C","2D Charge density and Electric field",1000,625);
  else
    C = new TCanvas("C","2D Charge density and Electric field",800,500);

  // Palettes setup
  TExec *exPlasma = new TExec("exPlasma","plasmaPalette->cd();");
  TExec *exElec   = new TExec("exElec","electronPalette->cd();");
  TExec *exField  = new TExec("exField","rbowPalette->cd();");
  
  // Text objects
  TPaveText *textTime = new TPaveText(0.55,0.85,0.82,0.9,"NDC");
  PlasmaGlob::SetPaveTextStyle(textTime,32); 
  char ctext[128];
  if(opt.Contains("units") && pData->GetPlasmaDensity()) 
    sprintf(ctext,"Z = %5.1f mm", 1e3 * pData->GetPlasmaSkinDepth() * Time);
  else
    sprintf(ctext,"T = %5.1f #omega_{p}^{-1}",Time);
  textTime->AddText(ctext);
 
  TPaveText *textDen = new TPaveText(0.12,0.85,0.45,0.9,"NDC");
  PlasmaGlob::SetPaveTextStyle(textDen,12); 
  textDen->SetTextColor(kOrange+10);
  if(opt.Contains("units") && pData->GetPlasmaDensity())
    sprintf(ctext,"n_{0} = %5.2f x 10^{15} / cc", 1e-6 * 1e-15 * pData->GetPlasmaDensity());
  else if(pData->GetBeamDensity() && pData->GetPlasmaDensity())
    sprintf(ctext,"n_{b}/n_{0} = %5.2f", pData->GetBeamDensity()/pData->GetPlasmaDensity());
  textDen->AddText(ctext);

  TPaveText *textWav = new TPaveText(0.12,0.78,0.45,0.83,"NDC");
  PlasmaGlob::SetPaveTextStyle(textWav,12); 
  textWav->SetTextColor(kGray+2);
  sprintf(ctext,"#lambda_{p} = %5.3f mm", 1e3 * pData->GetPlasmaWaveLength());
  textWav->AddText(ctext);

  
  // Actual Plotting!
  // ------------------------------------------------------------

  // Output file
  TString fOutName = Form("./%s/Plots/ChargeDensity2D/ChargeDensity2D",sim.Data());
  fOutName += Form("-%s_%i",sim.Data(),time);

  C->cd();
 
  gPad->SetFrameLineWidth(3);  

  TH2F *hFrame = (TH2F*) gROOT->FindObject("hFrame1");
  if(hFrame) delete hFrame;
  hFrame = (TH2F*) hDen2D[0]->Clone("hFrame1");
  hFrame->Reset();
 
  // hFrame->GetYaxis()->SetLabelSize(0.075);
  // hFrame->GetYaxis()->SetTitleSize(0.075);
  // hFrame->GetYaxis()->SetTitleOffset(0.5);
  // hFrame->GetXaxis()->SetLabelOffset(0.10);
  hFrame->Draw("col");
   
  // hDen2D[0]->GetZaxis()->SetLabelSize(0.06);  
  // hDen2D[0]->GetZaxis()->SetTitleSize(0.06);
  // hDen2D[0]->GetZaxis()->SetTitleOffset(0.06);
 

  exPlasma->Draw();
  hDen2D[0]->Draw("colz same");
  
  if(hDen2D[1]) {
    // hDen2D[1]->GetZaxis()->SetLabelSize(0.06);  
    // hDen2D[1]->GetZaxis()->SetTitleSize(0.06);
 

    exElec->Draw();
    hDen2D[1]->Draw("colz same");

    // DeBUG
    // Int_t nx = hDen2D[1]->GetNbinsX();
    // Int_t ny = hDen2D[1]->GetNbinsY();
    // Double_t lowx = hDen2D[1]->GetXaxis()->GetXmin();
    // Double_t higx = hDen2D[1]->GetXaxis()->GetXmax();
    // Double_t lowy = hDen2D[1]->GetYaxis()->GetXmin();
    // Double_t higy = hDen2D[1]->GetYaxis()->GetXmax();
    
    // Int_t bini[4] = {700,700,700,700};
    // Int_t binj[4] = {ny/2,ny/2+1,1,ny};
    // Double_t binx[4]; 
    // Double_t biny[4];
    // Double_t val[4];
  
    // cout << Form(" NX = %3i (%6f,%6f)  NY = %3i (%6f,%6f)",nx,lowx,higx,ny,lowy,higy) << endl;  
    // for(Int_t k=0; k<4;k++) {
    //   binx[k] = hDen2D[1]->GetXaxis()->GetBinCenter(bini[k]);
    //   biny[k] = hDen2D[1]->GetYaxis()->GetBinCenter(binj[k]);
    //   val[k]  = hDen2D[1]->GetBinContent(bini[k],binj[k]);
    //   cout << Form("       (%3i,%3i) = (%6f,%6f) -> value = %6f",bini[k],binj[k],binx[k],biny[k],val[k]) << endl;   
      
    // }
  }
  
  gPad->Update();
  TPaletteAxis *palette = (TPaletteAxis*)hDen2D[1]->GetListOfFunctions()->FindObject("palette");
  
  Float_t y1 = gPad->GetBottomMargin();
  Float_t y2 = 1 - gPad->GetTopMargin();
  palette->SetY2NDC(y2);
  palette->SetY1NDC(0.5*(y1+y2));
  //  palette->SetTitleOffset(0.75);
  palette->SetBorderSize(2);
  palette->SetLineColor(1);
  
  TPaletteAxis *palette2 = (TPaletteAxis*)hDen2D[0]->GetListOfFunctions()->FindObject("palette");
  palette2->SetY2NDC(0.5*(y1+y2));
  palette2->SetY1NDC(y1);
  //  palette2->SetTitleOffset(0.75);
  palette2->SetBorderSize(2);
  palette2->SetLineColor(1);

  textTime->Draw();
  textDen->Draw();
  if(opt.Contains("units"))
    textWav->Draw();
  
  gPad->RedrawAxis(); 

  C->cd();

  // Print to a file
  PlasmaGlob::imgconv(C,fOutName,opt);
  // ---------------------------------------------------------

}
