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

void PlotChargeRange2D( const TString &sim, Int_t time, Int_t index = 0, Int_t zoom=2, const TString &options="") {
  //, Double_t x1Min = -1, Double_t x1Max = -1, Double_t x2Min = -1, Double_t x2Max = -1) {
  
  PlasmaGlob::Initialize();

  TString opt = options;
 
  // Palettes!
  gROOT->Macro("PlasmaPalettes.C");

  if(opt.Contains("grid")) {
    gStyle->SetPadGridX(1);
    gStyle->SetPadGridY(1);
  }

  gStyle->SetPadRightMargin(0.20);   // Margin for palettes in 2D histos
  gStyle->SetLabelFont(42,"xyz");

  // Load PData
  PData *pData = PData::Get(sim.Data());
  pData->LoadFileNames(time);
  if(!pData->IsInit()) return;

  Bool_t CYL = kFALSE;
  if(sim.Contains("cyl")) CYL = kTRUE; 
    
  Bool_t ThreeD = kFALSE;
  if(sim.Contains("3D")) ThreeD = kTRUE; 

  // Set range for the charge calculation
  Double_t x1Min = -1;
  Double_t x1Max = -1;
  Double_t x2Min = -1;
  Double_t x2Max = -1;

  // Some plasma constants
  Double_t n0 = pData->GetPlasmaDensity();
  Double_t kp = pData->GetPlasmaK();
  Double_t skindepth = 1.;
  if(kp!=0.0) skindepth = 1/kp;
  Double_t E0 = pData->GetPlasmaE0();

  // Some beam properties:
  Float_t Ebeam = pData->GetBeamEnergy() * PUnits::MeV;
  Float_t gamma = Ebeam / PConst::ElectronMassE;
  Float_t vbeam = TMath::Sqrt(1 - 1/(gamma*gamma));
  Double_t rms0  = pData->GetBeamRmsY() * kp;
  // cout << Form(" - Bunch gamma      = %8.4f", gamma ) << endl;
  // cout << Form(" - Bunch velocity   = %8.4f c", vbeam ) << endl;
  
  // Time in OU
  Float_t Time = pData->GetRealTime();
  // z start of the plasma in normalized units.
  Float_t zStartPlasma = pData->GetPlasmaStart()*kp;
  // z start of the beam in normalized units.
  Float_t zStartBeam = pData->GetBeamStart()*kp;
  
  if(opt.Contains("center")) {
    Time -= zStartPlasma;
    if(opt.Contains("comov"))      // Centers on the head of the beam.
      Time += zStartBeam;
  } 

  // Get charge density histos
  Int_t Nspecies = pData->NSpecies();
  TH2F **hDen2D = new TH2F*[Nspecies];
  Double_t *Q = new Double_t[Nspecies];
  for(Int_t i=0;i<Nspecies;i++) {
    hDen2D[i] = NULL;
    Q[i] = 0.0;

    if(i!=index) continue;

    if(!pData->GetChargeFileName(i)) 
      continue;
    
    char hName[24];
    sprintf(hName,"hDen_%i",i);
    hDen2D[i] = (TH2F*) gROOT->FindObject(hName);
    if(hDen2D[i]) delete hDen2D[i];

    if(!ThreeD)
      hDen2D[i] = pData->GetCharge(i,opt);
    else
      hDen2D[i] = pData->GetCharge2DSliceZY(i,-1,1,opt+"avg");
    
  
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

    

    // INTEGRATED Beam's Charge:

    Int_t x1Nbin = hDen2D[i]->GetNbinsX();
    Int_t x2Nbin = hDen2D[i]->GetNbinsY();
    Double_t ymean = 0.0;
    Double_t ymean2 = 0.0;
    for(Int_t j=1;j<=x1Nbin;j++) {
      Double_t x = hDen2D[i]->GetXaxis()->GetBinCenter(j);
      if(x<x1Min && x1Min!=-1) continue;
      if(x>x1Max && x1Max!=-1) continue;
	
      for(Int_t k=1;k<=x2Nbin;k++) {
	Double_t y = hDen2D[i]->GetYaxis()->GetBinCenter(k);
	if(y<x2Min && x2Min!=-1) continue;
	if(y>x2Max && x2Max!=-1) continue;
	Double_t value  = hDen2D[i]->GetBinContent(j,k);
	ymean += value*y;
	ymean2 += value*y*y;
	if(CYL) {
	  Q[i] += y * value;
	  // cout << Form(" (%i,%i) -> radius = %7.4f , value = %7.4f",j,k,radius,value) << endl;
	} else {
	  Q[i] += value;
	}
      }    
    }
    ymean  /= Q[i];
    ymean2 /= Q[i];
    Double_t yrms  = TMath::Sqrt(ymean2 - ymean*ymean);  
    
    Double_t xbinsize = hDen2D[i]->GetXaxis()->GetBinWidth(1);
    Double_t ybinsize = hDen2D[i]->GetYaxis()->GetBinWidth(1); 
    Q[i] *= xbinsize * ybinsize;
   
 
    if(!CYL && !ThreeD) {
      // Q[i] *= TMath::Sqrt(2*TMath::Pi()) * yrms; 
      Q[i] *= TMath::Sqrt(2*TMath::Pi()) * rms0; 
    } else if(CYL) {
      Q[i] *= 2*TMath::Pi();
    }
    
    if(opt.Contains("units")) {
      Double_t dV = skindepth * skindepth * skindepth;
      Q[i] *= n0 * dV;
      Q[i] *= (PConst::ElectronCharge/PUnits::picocoulomb); 
      cout << Form(" Integrated charge of specie %3i = %8i pC",i, TMath::Nint(Q[i])) << endl;
    } else {
      cout << Form(" Integrated charge of specie %3i = %8.4f n0 * kp^-3",i,Q[i]) << endl;
    }
    // ------------------------------------------------------------------------------------
  }


  // Tunning the Histograms
  // ---------------------

  // Chaning to user units: 
  // --------------------------
  
  if(opt.Contains("units") && n0) {
    
    for(Int_t i=0;i<Nspecies;i++) {
      
      if(i!=index) continue;

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

      if(CYL)
	hDen2D[i]->GetYaxis()->SetTitle("r [mm]");      
      else
	hDen2D[i]->GetYaxis()->SetTitle("y [mm]");      

      if(opt.Contains("comov"))
	hDen2D[i]->GetXaxis()->SetTitle("#zeta [mm]");
      else
	hDen2D[i]->GetXaxis()->SetTitle("z [mm]");
      
      if(i==0)
	hDen2D[i]->GetZaxis()->SetTitle("n [10^{15}/cm^{3}]");
      else
	hDen2D[i]->GetZaxis()->SetTitle("n_{b} [10^{15}/cm^{3}]"); 
    }

  }

  // Zoom
  Float_t range    = (hDen2D[index]->GetYaxis()->GetXmax() - hDen2D[index]->GetYaxis()->GetXmin())/zoom;
  Float_t midPoint = (hDen2D[index]->GetYaxis()->GetXmax() + hDen2D[index]->GetYaxis()->GetXmin())/2.;
  
  if(!CYL) {
    hDen2D[index]->GetYaxis()->SetRangeUser(midPoint-range/2,midPoint+range/2);  
  } else {
    hDen2D[index]->GetYaxis()->SetRangeUser(hDen2D[index]->GetYaxis()->GetXmin(),range);  
  }

  // Set palette:
  PPalette * pPalette = (PPalette*) gROOT->FindObject("electron0");
  pPalette->cd();
    
  Float_t Max  = hDen2D[index]->GetMaximum();
  Float_t Min  = hDen2D[index]->GetMinimum();


  if(pData->GetChargeFileName(index)->find("plasma") != string::npos) { 
    
    // ------- z Zoom --------------------------------- Plasma palette -----------
    // Set the range of the plasma charge density histogram for maximum constrast 
    // using a dynamic palette wich adjust the nominal value to a certain color.
    
    Float_t density = 1;
    if(opt.Contains("units") && n0)
      density = n0 / (1e15/PUnits::cm3);
    Float_t Base = density;
    if(Max>Base && Max<2*Base) {
      Min  = 2.* Base - Max;
    } else if(Max <= Base)  {
      Max = 1.01 * Base;
      Min = 0;
    } else if(Max >= 2.0*Base) {
      Min = 0;
    }
    
  
    // Dynamic plasma palette
    const Int_t plasmaDNRGBs = 3;
    const Int_t plasmaDNCont = 64;
    Double_t basePos = 0.5;
    if(Max!=Min)
      basePos = (1.0/(Max-Min))*(Base - Min);
 
    Double_t plasmaDStops[plasmaDNRGBs] = { 0.00, basePos, 1.00 };
    Double_t plasmaDRed[plasmaDNRGBs]   = { 0.99, 0.90, 0.00 };
    Double_t plasmaDGreen[plasmaDNRGBs] = { 0.99, 0.90, 0.00 };
    Double_t plasmaDBlue[plasmaDNRGBs]  = { 0.99, 0.90, 0.00 };
    
    pPalette->CreateGradientColorTable(plasmaDNRGBs, plasmaDStops, 
				       plasmaDRed, plasmaDGreen, plasmaDBlue, plasmaDNCont);
  }
  
  if(hDen2D[index])
    hDen2D[index]->GetZaxis()->SetRangeUser(Min,Max); 

  
  // Plotting
  // -----------------------------------------------

  // Canvas setup
  TCanvas *C;
  if(opt.Contains("hres") && !opt.Contains("pdf")) // high resolution for plain grahics output.
    C = new TCanvas("C","2D Charge",1000,625);
  else
    C = new TCanvas("C","2D Charge",800,500);

    
  // Text objects
  TPaveText *textTime = new TPaveText(0.50,0.85,0.77,0.9,"NDC");
  PlasmaGlob::SetPaveTextStyle(textTime,32); 
  char ctext[128];
  if(opt.Contains("units") && pData->GetPlasmaDensity()) 
    sprintf(ctext,"Z = %5.1f mm", 1e3 * pData->GetPlasmaSkinDepth() * Time);
  else
    sprintf(ctext,"T = %5.1f #omega_{p}^{-1}",Time);
  textTime->AddText(ctext);
 
  TPaveText *textDen = new TPaveText(0.15,0.85,0.48,0.9,"NDC");
  PlasmaGlob::SetPaveTextStyle(textDen,12); 
  textDen->SetTextColor(kOrange+10);
  if(opt.Contains("units") && pData->GetPlasmaDensity())
    sprintf(ctext,"n_{0} = %5.2f x 10^{15} / cc", 1e-6 * 1e-15 * pData->GetPlasmaDensity());
  else if(pData->GetBeamDensity() && pData->GetPlasmaDensity())
    sprintf(ctext,"n_{b}/n_{0} = %5.2f", pData->GetBeamDensity()/pData->GetPlasmaDensity());
  textDen->AddText(ctext);

  TPaveText *textWav = new TPaveText(0.15,0.78,0.48,0.83,"NDC");
  PlasmaGlob::SetPaveTextStyle(textWav,12); 
  textWav->SetTextColor(kGray+2);
  sprintf(ctext,"#lambda_{p} = %5.3f mm", 1e3 * pData->GetPlasmaWaveLength());
  textWav->AddText(ctext);

  TPaveText *textCharge = new TPaveText(0.50,0.20,0.77,0.25,"NDC");
  PlasmaGlob::SetPaveTextStyle(textCharge,32); 
  textCharge->SetTextColor(kGray+1);
  if(opt.Contains("units")) 
    sprintf(ctext,"Q = %5i pC",TMath::Nint(Q[index]));
  else
    sprintf(ctext,"Q = %8.4f (c/#omega_{p})^{3} #times n_{0}",Q[index]);
  textCharge->AddText(ctext);
  

  
  // Actual Plotting!
  // ------------------------------------------------------------

  // Output file
  TString fOutName = Form("./%s/Plots/ChargeRange2D/ChargeRange2D",sim.Data());
  fOutName += Form("-%s_%i",sim.Data(),time);

  C->cd();
 
  gPad->SetFrameLineWidth(3);  

  TH2F *hFrame = (TH2F*) gROOT->FindObject("hFrame1");
  if(hFrame) delete hFrame;
  hFrame = (TH2F*) hDen2D[index]->Clone("hFrame1");
  hFrame->Reset();
 
  hFrame->Draw("col");
   
  hDen2D[index]->Draw("colz same");
  
  if(Nspecies>1) {
    if(hDen2D[1]) {
        
        hDen2D[1]->Draw("colz same");
    }
  }
  
  gPad->Update();
  
  textTime->Draw();
  textDen->Draw();
  if(opt.Contains("units"))
    textWav->Draw();
  
  textCharge->Draw();
  gPad->RedrawAxis(); 

  C->cd();

  // Print to a file
  PlasmaGlob::imgconv(C,fOutName,opt);
  // ---------------------------------------------------------

}
