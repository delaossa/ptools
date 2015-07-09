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


#ifndef __CINT__
#include "PFunctions.hh"
#endif

void PlotInjProb2D( const TString &sim, Int_t time, Int_t zoom=2, Int_t Nbins=2, const TString &options="") {
  
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
  gStyle->SetTitleFont(42,"xyz");
 
  // Load PData
  PData *pData = PData::Get(sim.Data());
  pData->LoadFileNames(time);
  if(!pData->IsInit()) return;

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
  const Int_t Nfields = 3;
  TH2F **hE2D = new TH2F*[Nfields];
  for(Int_t i=0;i<Nfields;i++) {
    hE2D[i] = NULL;

    //if(i!=index) continue;

    if(!pData->GetEfieldFileName(i))
      continue;
    
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
      hE2D[i]->GetXaxis()->SetTitle("#zeta [c/#omega_{p}]");
    else
      hE2D[i]->GetXaxis()->SetTitle("z [c/#omega_{p}]");
    
    if(pData->IsCyl()) 
      hE2D[i]->GetYaxis()->SetTitle("r [c/#omega_{p}]");
    else
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
  
  // Chaning to user units: 
  // --------------------------
  
  if(opt.Contains("units") && n0) {
    
    for(Int_t i=0;i<Nfields;i++) {

      //if(i!=index) continue;

      Int_t NbinsX = hE2D[i]->GetNbinsX();
      Float_t xMin = skindepth * hE2D[i]->GetXaxis()->GetXmin() / PUnits::um;
      Float_t xMax = skindepth * hE2D[i]->GetXaxis()->GetXmax() / PUnits::um;
      Int_t NbinsY = hE2D[i]->GetNbinsY();
      Float_t yMin = skindepth * hE2D[i]->GetYaxis()->GetXmin() / PUnits::um;
      Float_t yMax = skindepth * hE2D[i]->GetYaxis()->GetXmax() / PUnits::um;
      hE2D[i]->SetBins(NbinsX,xMin,xMax,NbinsY,yMin,yMax);
      for(Int_t j=0;j<NbinsX;j++) {
	for(Int_t k=0;k<NbinsY;k++) {
	  hE2D[i]->SetBinContent(j,k, hE2D[i]->GetBinContent(j,k) * ( E0 / (PUnits::GV/PUnits::m) ) );
	}
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
    }
    
  }


  // Now, combine the electric field components into the total |E|
  // and calculate ionization probability for He:
  // Outter Helium electron
  Double_t Eion0 = 24.59 * PUnits::eV;
  Double_t Z     = 1;
  Double_t l     = 0;
  Double_t m     = 0;
  
  TH2F *hETotal2D = (TH2F*) hE2D[0]->Clone("hETotal2D");
  hETotal2D->Reset();
  TH2F *hIonProb2D = (TH2F*) hE2D[0]->Clone("hIonProb2D");
  hIonProb2D->Reset();
  Int_t NbinsX = hE2D[0]->GetNbinsX();
  Int_t NbinsY = hE2D[0]->GetNbinsY();
  for(Int_t j=0;j<NbinsX;j++) {
    for(Int_t k=0;k<NbinsY;k++) {
      Double_t E1 = hE2D[0]->GetBinContent(j,k);
      Double_t E2 = hE2D[1]->GetBinContent(j,k);
      Double_t E3 = hE2D[2]->GetBinContent(j,k);
      Double_t E  = TMath::Sqrt(E1*E1+E2*E2+E3*E3);
    
      hETotal2D->SetBinContent(j,k,E);
    
      if(opt.Contains("units") && n0) E *= PUnits::GV/PUnits::m;
      else E *= E0;
      
      Double_t IonProb = (PFunc::ADK(E,Eion0,Z,l,m)/PUnits::atomictime)*PUnits::femtosecond;
      //if(IonProb>10) cout << "Ion prob = " << IonProb << endl;
      hIonProb2D->SetBinContent(j,k,IonProb);
    }
  }
  
  if(opt.Contains("units") && n0)
    hETotal2D->GetZaxis()->SetTitle("E [GV/m]");
  else
    hETotal2D->GetZaxis()->SetTitle("E [E_{0}]");

  hIonProb2D->GetZaxis()->SetTitle("W_{ADK} [fs^{-1}]");
  
  
  // Zoom
  Float_t range    = (hETotal2D->GetYaxis()->GetXmax() - hETotal2D->GetYaxis()->GetXmin())/zoom;
  Float_t midPoint = (hETotal2D->GetYaxis()->GetXmax() + hETotal2D->GetYaxis()->GetXmin())/2.;
  
  Float_t yMin,yMax;
  if(!pData->IsCyl()) {
    yMin = midPoint-range/2;
    yMax = midPoint+range/2;
  } else {
    yMin = 0.;
    yMax = range;
  }
  hETotal2D->GetYaxis()->SetRangeUser(yMin,yMax);
  
  
  
  // // Change the range of z axis for the fields to be symmetric.
  // Float_t Emax = hETotal2D->GetMaximum();
  // Float_t Emin = hETotal2D->GetMinimum();
  // if(Emax > TMath::Abs(Emin))
  //   Emin = -Emax;
  // else
  //   Emax = -Emin;
  // hETotal2D->GetZaxis()->SetRangeUser(Emin,Emax); 

  // Plotting
  // -----------------------------------------------

  // Canvas setup
  TCanvas *C;
  C = new TCanvas("C","2D Electric field",600,800);
  
  // Palettes setup
  TExec *exField  = new TExec("exField","electronPalette->cd();");
  TExec *exProb   = new TExec("exProb","electronPalette->cd();");
  
  // Text objects
  TPaveText *textTime = new TPaveText(0.52,0.87,0.79,0.92,"NDC");
  PlasmaGlob::SetPaveTextStyle(textTime,32); 
  char ctext[128];
  if(opt.Contains("units") && pData->GetPlasmaDensity()) 
    sprintf(ctext,"z = %5.1f #mum", skindepth * Time / PUnits::um);
  else
    sprintf(ctext,"z = %5.1f #omega_{p}^{-1}",Time);
  textTime->AddText(ctext);
 
  TPaveText *textDen = new TPaveText(0.13,0.85,0.38,0.92,"NDC");
  PlasmaGlob::SetPaveTextStyle(textDen,12); 
  textDen->SetTextColor(kOrange+10);
  if(opt.Contains("units") && pData->GetPlasmaDensity())
    sprintf(ctext,"n_{0} = %5.2f x 10^{17} / cm^{3}", pData->GetPlasmaDensity() / (1e17/PUnits::cm3));
  else if(pData->GetBeamDensity() && pData->GetPlasmaDensity())
    sprintf(ctext,"n_{b}/n_{0} = %5.2f", pData->GetBeamDensity()/pData->GetPlasmaDensity());
  textDen->AddText(ctext);
  
  // Actual Plotting!
  // ------------------------------------------------------------

  // Output file
  TString fOutName = Form("./%s/Plots/InjProb2D/InjProb2D",sim.Data());
  fOutName += Form("-%s_%i",sim.Data(),time);

  C->Divide(1,2);

  C->cd(1);
  
  gPad->SetFrameLineWidth(3);  
  if(opt.Contains("logz")) {
    gPad->SetLogz(1);
  } else {
    gPad->SetLogz(0);
  }

  TH2F *hFrame = (TH2F*) gROOT->FindObject("hFrame1");
  if(hFrame) delete hFrame;
  hFrame = (TH2F*) hETotal2D->Clone("hFrame1");
  hFrame->Reset();
 
  hFrame->Draw("col");

  exField->Draw();
  hETotal2D->Draw("colz same");

  gPad->Update();

  textTime->Draw();
  //textDen->Draw();

  gPad->RedrawAxis(); 

  C->cd(2);
  
  gPad->SetFrameLineWidth(3);  
  if(opt.Contains("logz")) {
    gPad->SetLogz(1);
  } else {
    gPad->SetLogz(0);
  }

  hFrame = (TH2F*) gROOT->FindObject("hFrame2");
  if(hFrame) delete hFrame;
  hFrame = (TH2F*) hETotal2D->Clone("hFrame2");
  hFrame->Reset();
  
  hFrame->Draw("col");

  exProb->Draw();
  hIonProb2D->Draw("colz same");
  
  gPad->Update();
  gPad->RedrawAxis(); 

  C->cd();

  // Print to a file
  PlasmaGlob::imgconv(C,fOutName,opt);
  // ---------------------------------------------------------

}
