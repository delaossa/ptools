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

void PlotCFF2D( const TString &sim, Int_t time, Int_t zoom=2, Int_t Nbins=2, const TString &options="") {
  
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
   
  // -------------------------------------------------------------------------------
  
  // Get charge density histos
  Int_t Nspecies = pData->NSpecies();
  TH2F **hDen2D = new TH2F*[Nspecies];
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
  }
  
  // Get electric fields
  const Int_t Nfields = 2;
  TH2F **hE2D = new TH2F*[Nfields];
  for(Int_t i=0;i<Nfields;i++) {
    hE2D[i] = NULL;
  
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
  }

  
  // --------------------------------------------------- Vertical Zoom ------------
  
  Float_t range    = (pData->GetXMax(1)-pData->GetXMin(1))/zoom;
  Float_t midPoint = (pData->GetXMax(1)+pData->GetXMin(1))/2.;
  Double_t ymin = midPoint-range/2;
  Double_t ymax = midPoint+range/2;
  if(pData->IsCyl()) {
    ymin = pData->GetXMin(1);
    ymax = range;
  }

  for(Int_t i=0;i<Nspecies;i++) {
    if(!hDen2D[i]) continue;
    hDen2D[i]->GetYaxis()->SetRangeUser(ymin,ymax);
  }
  
  for(Int_t i=0;i<Nfields;i++) {
    if(!hE2D[i]) continue;
    hE2D[i]->GetYaxis()->SetRangeUser(ymin,ymax);
  }
  
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
    Min[i] = 1E-1 * Base;
    if(i==2) Min[i] = 1E-3 * Base;
    hDen2D[i]->GetZaxis()->SetRangeUser(Min[i],Max[i]);
  }
  
  // Dynamic plasma palette
  const Int_t plasmaDNRGBs = 3;
  const Int_t plasmaDNCont = 64;
  Double_t basePos = 0.5;
  if(Max[0]!=Min[0]) {
    if(opt.Contains("logz")) {
      Float_t a = 1.0/(TMath::Log10(Max[0])-TMath::Log10(Min[0]));
      Float_t b = TMath::Log10(Min[0]);
      basePos = a*(TMath::Log10(Base) - b);
      
    } else {
      basePos = (1.0/(Max[0]-Min[0]))*(Base - Min[0]);
    }
  }

  Double_t plasmaDStops[plasmaDNRGBs] = { 0.00, basePos, 1.00 };
  Double_t plasmaDRed[plasmaDNRGBs]   = { 0.99, 0.90, 0.00 };
  Double_t plasmaDGreen[plasmaDNRGBs] = { 0.99, 0.90, 0.00 };
  Double_t plasmaDBlue[plasmaDNRGBs]  = { 0.99, 0.90, 0.00 };
   
  PPalette * plasmaPalette = (PPalette*) gROOT->FindObject("plasma");
  plasmaPalette->CreateGradientColorTable(plasmaDNRGBs, plasmaDStops, 
					  plasmaDRed, plasmaDGreen, plasmaDBlue, plasmaDNCont);
  
  
  // Change the range of z axis for the fields to be symmetric.
  for(Int_t i=0;i<Nfields;i++) {
    if(!hE2D[i]) continue;
    Float_t Emax = hE2D[i]->GetMaximum();
    Float_t Emin = hE2D[i]->GetMinimum();
    if(Emax > TMath::Abs(Emin))
      Emin = -Emax;
    else
      Emax = -Emin;
    hE2D[i]->GetZaxis()->SetRangeUser(Emin,Emax); 
  }
  
  
  // Plotting
  // -----------------------------------------------

  // Canvas setup
  TCanvas *C = new TCanvas("C","2D Charge density, Accelerating and focusing fields",750,1000);
  
  Int_t NPADS = 3;
  TPad **pad = new TPad*[NPADS];
  PlasmaGlob::CanvasPartition(C,3);
  
  // Palettes setup
  TExec *exPlasma = new TExec("exPlasma","plasmaPalette->cd();");
  TExec *exHot    = new TExec("exHot","hotPalette->cd();");
  TExec *exElec   = new TExec("exElec","electronPalette->cd();");
  TExec *exField  = new TExec("exField","rbow2Palette->cd();");
  
  // Draw!

  // <---------------------------------------------- Top Plot ---------
  char name[16];
  sprintf(name,"pad_%i",NPADS-1);
  pad[0] = (TPad*) gROOT->FindObject(name);
  
  pad[0]->Draw();
  pad[0]->cd();
  
  if(opt.Contains("logz")) {
    pad[0]->SetLogz(1);
  } else {
    pad[0]->SetLogz(0);
  }
  pad[0]->SetFrameLineWidth(3);  
  
  TH2F *hFrame = (TH2F*) gROOT->FindObject("hFrame1");
  if(hFrame) delete hFrame;
  hFrame = (TH2F*) hDen2D[0]->Clone("hFrame1");
  hFrame->Reset();
  
  Float_t vfactor = 1.0;

  hFrame->GetXaxis()->SetLabelOffset(999);

  hFrame->GetYaxis()->SetLabelSize(0.065*vfactor);
  hFrame->GetYaxis()->SetLabelOffset(0.02/vfactor);
  hFrame->GetYaxis()->SetTickLength(0.02/vfactor);
  hFrame->GetYaxis()->SetTitleSize(0.075*vfactor);
  hFrame->GetYaxis()->SetTitleOffset(0.65/vfactor);

  hFrame->GetZaxis()->SetLabelSize(0.06*vfactor);  
  hFrame->GetZaxis()->SetTickLength(0.02/vfactor); 
  hFrame->GetZaxis()->SetTitleSize(0.06*vfactor);                        
  hFrame->GetZaxis()->SetTitleOffset(0.45/vfactor);
  
  hFrame->Draw("col");            

  if(Nspecies>=3) {
    if(hDen2D[2]) {
      exHot->Draw();
      hDen2D[2]->Draw("col same");
    }
  }
  exPlasma->Draw();
  hDen2D[0]->Draw("col same");
  
  if(hDen2D[1]) {
    exElec->Draw();
    hDen2D[1]->Draw("col same");
  }

  pad[0]->RedrawAxis();
  

  // <---------------------------------------------- Mid Plot ---------
  C->cd(0);

  sprintf(name,"pad_%i",NPADS-2);
  pad[1] = (TPad*) gROOT->FindObject(name);
  
  pad[1]->Draw();
  pad[1]->cd();
  
  pad[1]->SetFrameLineWidth(3);  
  
  hFrame = (TH2F*) gROOT->FindObject("hFrame2");
  if(hFrame) delete hFrame;
  hFrame = (TH2F*) hE2D[0]->Clone("hFrame2");
  hFrame->Reset();
  
  vfactor =  pad[0]->GetAbsHNDC()/pad[1]->GetAbsHNDC(); 

  hFrame->GetXaxis()->SetLabelOffset(999);

  hFrame->GetYaxis()->SetLabelSize(0.065*vfactor);
  hFrame->GetYaxis()->SetLabelOffset(0.02/vfactor);
  hFrame->GetYaxis()->SetTickLength(0.02/vfactor);
  hFrame->GetYaxis()->SetTitleSize(0.075*vfactor);
  hFrame->GetYaxis()->SetTitleOffset(0.65/vfactor);

  hFrame->GetZaxis()->SetLabelSize(0.06*vfactor);  
  hFrame->GetZaxis()->SetTickLength(0.02/vfactor); 
  hFrame->GetZaxis()->SetTitleSize(0.06*vfactor);                        
  hFrame->GetZaxis()->SetTitleOffset(0.45/vfactor);
  
  hFrame->Draw("col");            

  
  exField->Draw();
  hE2D[0]->Draw("col same");
 
  pad[1]->RedrawAxis();


  // <---------------------------------------------- Bottom Plot ---------
  C->cd(0);
  
  sprintf(name,"pad_%i",NPADS-3);
  pad[2] = (TPad*) gROOT->FindObject(name);
  
  pad[2]->Draw();
  pad[2]->cd();
  
  pad[2]->SetFrameLineWidth(3);  
  
  hFrame = (TH2F*) gROOT->FindObject("hFrame3");
  if(hFrame) delete hFrame;
  hFrame = (TH2F*) hE2D[1]->Clone("hFrame3");
  hFrame->Reset();
  
  vfactor =  pad[0]->GetAbsHNDC()/pad[2]->GetAbsHNDC(); 

  //  hFrame->GetXaxis()->SetLabelOffset(999);

  hFrame->GetYaxis()->SetLabelSize(0.065*vfactor);
  hFrame->GetYaxis()->SetLabelOffset(0.02/vfactor);
  hFrame->GetYaxis()->SetTickLength(0.02/vfactor);
  hFrame->GetYaxis()->SetTitleSize(0.075*vfactor);
  hFrame->GetYaxis()->SetTitleOffset(0.65/vfactor);

  hFrame->GetZaxis()->SetLabelSize(0.06*vfactor);  
  hFrame->GetZaxis()->SetTickLength(0.02/vfactor); 
  hFrame->GetZaxis()->SetTitleSize(0.06*vfactor);                        
  hFrame->GetZaxis()->SetTitleOffset(0.45/vfactor);

  hFrame->Draw("col");            
  
  exField->Draw();
  hE2D[1]->Draw("col same");
 
  pad[2]->RedrawAxis();

  C->cd();

  // Output file
  TString fOutName = Form("./%s/Plots/CFF2D/CFF2D",pData->GetPath().c_str());
  fOutName += Form("-%s_%i",pData->GetName(),time);
  // Print to a file
  PlasmaGlob::imgconv(C,fOutName,opt);
  // ---------------------------------------------------------


  
}
