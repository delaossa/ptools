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
#include "PGlobals.hh"
#include "PPalette.hh"

void PlotField2D( const TString &sim, Int_t time, Int_t index = 0, Float_t zoom=2, Int_t Nbins=2, const TString &options="") {
    
#ifdef __CINT__  
  gSystem->Load("libptools.so");
#endif

  // Load PData
  PData *pData = PData::Get(sim.Data());
  pData->LoadFileNames(time);
  if(!pData->IsInit()) return;

  // Refresh and Style
  PGlobals::Initialize();

  // Coloured palettes
  gROOT->Macro("PPalettes.C");

  TString opt = options;
 
  // More makeup
  if(opt.Contains("grid")) {
    gStyle->SetPadGridX(1);
    gStyle->SetPadGridY(1);
  }
 
  gStyle->SetTitleFont(42);
  gStyle->SetStatFont(42);
  gStyle->SetTextFont(42);
  gStyle->SetTitleFont(42,"xyz");
  gStyle->SetLabelFont(42,"xyz");

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
  Float_t nb = pData->GetBeamDensity();
  
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
    Nbins =  TMath::Nint(rms0 / pData->GetDX(1));
  }
  
  // Slice width limits.
  if(!pData->IsCyl()) {
    FirstyBin = pData->GetNX(1)/2 + 1 - Nbins;
    LastyBin =  pData->GetNX(1)/2 + Nbins;
  } else {
    FirstyBin = 1; 
    LastyBin  = Nbins;
  }
  // -------------------------------------------------------------------------------

  // Get field histos
  const Int_t Nfields = 3;
  TH2F **hE2D = new TH2F*[Nfields];
  for(Int_t i=0;i<Nfields;i++) {
    hE2D[i] = NULL;

    if(i!=index) continue;

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
      hE2D[i]->GetXaxis()->SetTitle("k_{p} #zeta");
    else
      hE2D[i]->GetXaxis()->SetTitle("k_{p} z");
    
    if(pData->IsCyl()) 
      hE2D[i]->GetYaxis()->SetTitle("k_{p} r");
    else
      hE2D[i]->GetYaxis()->SetTitle("k_{p} x");
    
    if(i==0)
      hE2D[i]->GetZaxis()->SetTitle("E_{z} / E_{0}");
    else if(i==1)
      hE2D[i]->GetZaxis()->SetTitle("E_{y} / E_{0}");
    else if(i==2)
      hE2D[i]->GetZaxis()->SetTitle("E_{x} / E_{0}");
  }


  // Tunning the Histograms
  // ---------------------

  // Chaning to user units: 
  // --------------------------
  
  if(opt.Contains("units") && n0) {
    
    for(Int_t i=0;i<Nfields;i++) {
      
      if(i!=index) continue;

      Int_t NbinsX = hE2D[i]->GetNbinsX();
      Float_t xMin = skindepth * hE2D[i]->GetXaxis()->GetXmin() / PUnits::um;
      Float_t xMax = skindepth * hE2D[i]->GetXaxis()->GetXmax() / PUnits::um;
      Int_t NbinsY = hE2D[i]->GetNbinsY();
      Float_t yMin = skindepth * hE2D[i]->GetYaxis()->GetXmin() / PUnits::um;
      Float_t yMax = skindepth * hE2D[i]->GetYaxis()->GetXmax() / PUnits::um;
      hE2D[i]->SetBins(NbinsX,xMin,xMax,NbinsY,yMin,yMax);
      for(Int_t j=0;j<hE2D[i]->GetNbinsX();j++) {
	for(Int_t k=0;k<hE2D[i]->GetNbinsY();k++) {
	  hE2D[i]->SetBinContent(j,k, hE2D[i]->GetBinContent(j,k) * ( E0 / (PUnits::GV/PUnits::m) ) );
	}
      }

      if(pData->IsCyl())
	hE2D[i]->GetYaxis()->SetTitle("r [#mum]");      
      else
	hE2D[i]->GetYaxis()->SetTitle("x [#mum]");      

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

  // --------------------------------------------------- Vertical Zoom ------------
  
  Float_t yRange    = (hE2D[index]->GetYaxis()->GetXmax() - hE2D[index]->GetYaxis()->GetXmin())/zoom;
  Float_t midPoint = (hE2D[index]->GetYaxis()->GetXmax() + hE2D[index]->GetYaxis()->GetXmin())/2.;
  Float_t yMin = midPoint-yRange/2;
  Float_t yMax = midPoint+yRange/2;
  if(pData->IsCyl()) {
    yMin = hE2D[index]->GetYaxis()->GetXmin();
    yMax = yRange;
  }
  
  hE2D[index]->GetYaxis()->SetRangeUser(yMin,yMax);

  Float_t xMin = hE2D[index]->GetXaxis()->GetXmin();
  Float_t xMax = hE2D[index]->GetXaxis()->GetXmax();
  Float_t xRange = xMax - xMin;
 

  // Change the range of z axis for the fields to be symmetric.
  Float_t Emax = hE2D[index]->GetMaximum();
  Float_t Emin = hE2D[index]->GetMinimum();
  if(Emax > TMath::Abs(Emin))
    Emin = -Emax;
  else
    Emax = -Emin;
  hE2D[index]->GetZaxis()->SetRangeUser(Emin,Emax); 

  
  PPalette * plasmaPalette = (PPalette*) gROOT->FindObject("rbowwhite");
  plasmaPalette->cd();

  // Plotting
  // -----------------------------------------------

  // Canvas setup
  TCanvas *C = new TCanvas("C","2D electric",800,500);
  
  // Actual Plotting!
  // ------------------------------------------------------------

  // Output file
  TString fname;
  if(index<3) fname = Form("E%i",index+1);
  TString fOutName = Form("./%s/Plots/Field2D/Field2D-%s",sim.Data(),fname.Data());
  fOutName += Form("-%s_%i",sim.Data(),time);
  
  C->cd();

  // Setup Pad layout:
  Float_t lMargin = 0.15;
  Float_t rMargin = 0.18;
  Float_t bMargin = 0.20;
  Float_t tMargin = 0.06;
  gPad->SetLeftMargin(lMargin);
  gPad->SetRightMargin(rMargin);
  gPad->SetBottomMargin(bMargin);
  gPad->SetTopMargin(tMargin);
  
  if(opt.Contains("logz")) {
    gPad->SetLogz(1);
  } else {
    gPad->SetLogz(0);
  }
  gPad->SetFrameLineWidth(2);  

  // Define the frames for plotting
  Int_t fonttype = 43;
  Int_t fontsize = 28;
  Int_t tfontsize = 30;
  Float_t txoffset = 1.3;
  Float_t lxoffset = 0.02;
  Float_t tyoffset = 1.0;
  Float_t lyoffset = 0.01;
  Float_t tylength = 0.02;
  Float_t txlength = 0.04;

  TH2F *hFrame = (TH2F*) gROOT->FindObject("hFrame");
  if(hFrame) delete hFrame;
  hFrame = (TH2F*) hE2D[index]->Clone("hFrame");
  hFrame->Reset();

  // Format for y axis
  hFrame->GetYaxis()->SetTitleFont(fonttype);
  hFrame->GetYaxis()->SetTitleSize(tfontsize);
  hFrame->GetYaxis()->SetTitleOffset(tyoffset);
  hFrame->GetYaxis()->SetLabelFont(fonttype);
  hFrame->GetYaxis()->SetLabelSize(fontsize);
  hFrame->GetYaxis()->SetLabelOffset(lyoffset);

  hFrame->GetYaxis()->SetTickLength(tylength);

  // Format for x axis
  hFrame->GetXaxis()->SetTitleFont(fonttype);
  hFrame->GetXaxis()->SetTitleSize(tfontsize+2);
  hFrame->GetXaxis()->SetTitleOffset(txoffset);
  hFrame->GetXaxis()->SetLabelFont(fonttype);
  hFrame->GetXaxis()->SetLabelSize(fontsize+2);
  hFrame->GetXaxis()->SetLabelOffset(lxoffset);
    
  hFrame->GetXaxis()->SetTickLength(txlength);      

  hFrame->Draw("col");            
 
  //  hE2D[index]->GetZaxis()->SetNdivisions(505);  
  hE2D[index]->GetZaxis()->SetTitleFont(fonttype);
  hE2D[index]->Draw("colz same");
    
  // Re-touchs
  gPad->Update();

  Float_t y1 = gPad->GetBottomMargin();
  Float_t y2 = 1 - gPad->GetTopMargin();
  Float_t x1 = gPad->GetLeftMargin();
  Float_t x2 = 1 - gPad->GetRightMargin();
  Float_t gap = 0.005;  

  TPaletteAxis *palette = (TPaletteAxis*)hE2D[index]->GetListOfFunctions()->FindObject("palette");
  if(palette) {
    palette->SetY2NDC(y2 - gap);
    palette->SetY1NDC(y1 + gap);
    palette->SetX1NDC(x2 + 0.005);
    palette->SetX2NDC(x2 + 0.03);
    palette->SetTitleOffset(tyoffset);
    palette->SetTitleSize(tfontsize);
    palette->SetLabelFont(fonttype);
    palette->SetLabelSize(fontsize);
    if(opt.Contains("logz")) 
      palette->SetLabelOffset(0);
    else
      palette->SetLabelOffset(lyoffset);
    palette->SetBorderSize(2);
    palette->SetLineColor(1);
  }
 
  
  // Text objects
  TPaveText *textTime = new TPaveText(xMax - 0.30*xRange, yMax-0.10*yRange, xMax-0.02*xRange, yMax-0.02*yRange);
  PGlobals::SetPaveTextStyle(textTime,32); 
  char ctext[128];
  if(opt.Contains("units") && n0) 
    sprintf(ctext,"z = %5.1f #mum", Time * skindepth / PUnits::um);
  else
    sprintf(ctext,"#omega_{p} t = %5.1f",Time);
  textTime->AddText(ctext);
  
  TPaveText *textDen = new TPaveText(xMin + 0.02*xRange, yMax-0.10*yRange, xMin + 0.40*xRange, yMax-0.02*yRange);
  PGlobals::SetPaveTextStyle(textDen,12); 
  textDen->SetTextColor(kOrange+10);
  if(opt.Contains("units") && n0) {
    sprintf(ctext,"n_{0} = %5.2f x 10^{17} / cm^{3}", 1e-17 * n0 * PUnits::cm3);
    textDen->AddText(ctext);
    textDen->Draw();
  }
  
  textTime->Draw();
  
  gPad->RedrawAxis(); 

  C->cd();

  // Print to a file
  PGlobals::imgconv(C,fOutName,opt);
  // ---------------------------------------------------------

}
