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

void PlotBFields2D( const TString &sim, Int_t time, Int_t zoom=2, Int_t Nbins=2, const TString &options="") {
  
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
  
  
 
  
  // Get magnetic fields 2D
  const Int_t Nfields = 3;
  TH2F **hB2D = new TH2F*[Nfields];
  for(Int_t i=0;i<Nfields;i++) {
    hB2D[i] = NULL;

    if(!pData->GetBfieldFileName(i))
      continue;

    cout << Form(" Getting magnetic field number ") << i+1 << endl;
    
    char hName[24];
    sprintf(hName,"hB2D_%i",i);
    hB2D[i] = (TH2F*) gROOT->FindObject(hName);
    if(hB2D[i]) delete hB2D[i];
    
    if(!pData->Is3D())
      hB2D[i] = pData->GetBField(i,opt);
    else
      hB2D[i] = pData->GetBField2DSliceZY(i,-1,Nbins,opt+"avg");
    
    hB2D[i]->SetName(hName);   
    hB2D[i]->GetXaxis()->CenterTitle();
    hB2D[i]->GetYaxis()->CenterTitle();
    hB2D[i]->GetZaxis()->CenterTitle();
    if(opt.Contains("comov"))
      hB2D[i]->GetXaxis()->SetTitle("k_{p} #zeta");
    else
      hB2D[i]->GetXaxis()->SetTitle("k_{p} z");
    
    if(pData->IsCyl()) 
      hB2D[i]->GetYaxis()->SetTitle("k_{p} r");
    else
      hB2D[i]->GetYaxis()->SetTitle("k_{p} y");
    
    if(i==0)
      hB2D[i]->GetZaxis()->SetTitle("B_{z}/E_{0}");
    else if(i==1)
      hB2D[i]->GetZaxis()->SetTitle("B_{y}/E_{0}");
    else if(i==2)
      hB2D[i]->GetZaxis()->SetTitle("B_{x}/E_{0}");
      
  }

  // Tunning the Histograms
  // ---------------------
  
  // Chaning to user units: 
  // --------------------------
  
  if(opt.Contains("units") && n0) {
    
    for(Int_t i=0;i<Nfields;i++) {
      Int_t NbinsX = hB2D[i]->GetNbinsX();
      Float_t xMin = skindepth * hB2D[i]->GetXaxis()->GetXmin() / PUnits::um;
      Float_t xMax = skindepth * hB2D[i]->GetXaxis()->GetXmax() / PUnits::um;
      Int_t NbinsY = hB2D[i]->GetNbinsY();
      Float_t ymin = skindepth * hB2D[i]->GetYaxis()->GetXmin() / PUnits::um;
      Float_t ymax = skindepth * hB2D[i]->GetYaxis()->GetXmax() / PUnits::um;
      hB2D[i]->SetBins(NbinsX,xMin,xMax,NbinsY,ymin,ymax);
            
      for(Int_t j=0;j<hB2D[i]->GetNbinsX();j++) {
	for(Int_t k=0;k<hB2D[i]->GetNbinsY();k++) {
	  hB2D[i]->SetBinContent(j,k, hB2D[i]->GetBinContent(j,k) * ( E0 / (PUnits::GV/PUnits::m) ) );
	}

      }
      
      if(pData->IsCyl())
	hB2D[i]->GetYaxis()->SetTitle("r [#mum]");      
      else
	hB2D[i]->GetYaxis()->SetTitle("y [#mum]");      
      
      if(opt.Contains("comov"))
	hB2D[i]->GetXaxis()->SetTitle("#zeta [#mum]");
      else
	hB2D[i]->GetXaxis()->SetTitle("z [#mum]");
      
      if(i==0)
	hB2D[i]->GetZaxis()->SetTitle("B_{z} [GV/m]");
      else if(i==1)
	hB2D[i]->GetZaxis()->SetTitle("B_{y} [GV/m]");
      else if(i==2)
	hB2D[i]->GetZaxis()->SetTitle("B_{x} [GV/m]");
      
      
    }
  }


  // --------------------------------------------------- Vertical Zoom ------------
  
  Float_t yRange    = (hB2D[0]->GetYaxis()->GetXmax() - hB2D[0]->GetYaxis()->GetXmin())/zoom;
  Float_t midPoint = (hB2D[0]->GetYaxis()->GetXmax() + hB2D[0]->GetYaxis()->GetXmin())/2.;
  Float_t yMin = midPoint-yRange/2;
  Float_t yMax = midPoint+yRange/2;
  if(pData->IsCyl()) {
    yMin = pData->GetXMin(1);
    yMax = yRange;
  }

  for(Int_t i=0;i<Nfields;i++) {
    if(!hB2D[i]) continue;
    hB2D[i]->GetYaxis()->SetRangeUser(yMin,yMax);
  }

  Float_t xMin = hB2D[0]->GetXaxis()->GetXmin();
  Float_t xMax = hB2D[0]->GetXaxis()->GetXmax();
  Float_t xRange = xMax - xMin;
  
  // Change the range of z axis for the fields to be symmetric.
  Float_t *Bmax = new Float_t[Nfields];
  Float_t *Bmin = new Float_t[Nfields];
  for(Int_t i=0;i<Nfields;i++) {
    Bmax[i] = hB2D[i]->GetMaximum();
    Bmin[i] = hB2D[i]->GetMinimum();
    if(Bmax[i] > TMath::Abs(Bmin[i]))
      Bmin[i] = -Bmax[i];
    else
      Bmax[i] = -Bmin[i];
    hB2D[i]->GetZaxis()->SetRangeUser(Bmin[i],Bmax[i]); 
    //hB2D[i]->GetZaxis()->SetRangeUser(Bmin[0],Bmax[0]); 
  }
  
  // "Axis range" in Osiris units:
  Double_t ylow  = hB2D[0]->GetYaxis()->GetBinLowEdge(FirstyBin);
  Double_t yup = hB2D[0]->GetYaxis()->GetBinUpEdge(LastyBin);
  Double_t xmin = hB2D[0]->GetXaxis()->GetXmin();
  Double_t xmax = hB2D[0]->GetXaxis()->GetXmax();

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
  TString fOutName = Form("./%s/Plots/BFields2D/BFields2D",pData->GetPath().c_str());
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
  hFrame = (TH2F*) hB2D[0]->Clone("hFrame1");
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

    //  hB2D[0]->GetZaxis()->SetNdivisions(505);
  hB2D[0]->GetZaxis()->SetTitleFont(42);
  hB2D[0]->GetZaxis()->SetTickLength(0.02);
  
  exField->Draw();
  hB2D[0]->Draw("colz same");

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
  
  pad[0]->Update();
  
  Float_t y1 = pad[0]->GetBottomMargin();
  Float_t y2 = 1 - pad[0]->GetTopMargin();
  Float_t x1 = pad[0]->GetLeftMargin();
  Float_t x2 = 1 - pad[0]->GetRightMargin();
  
  TPaletteAxis *palette = (TPaletteAxis*) hB2D[0]->GetListOfFunctions()->FindObject("palette");  
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
  hFrame = (TH2F*) hB2D[0]->Clone("hFrame2");
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

    //  hB2D[1]->GetZaxis()->SetNdivisions(505);
  hB2D[1]->GetZaxis()->SetTitleFont(42);
  hB2D[1]->GetZaxis()->SetTickLength(0.02/yFactor);

  exField->Draw();
  hB2D[1]->Draw("colz same");

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
  
  palette = (TPaletteAxis*) hB2D[1]->GetListOfFunctions()->FindObject("palette");  
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
  
  palette = (TPaletteAxis*) hB2D[1]->GetListOfFunctions()->FindObject("palette");  
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
  hFrame = (TH2F*) hB2D[0]->Clone("hFrame3");
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

  //  hB2D[2]->GetZaxis()->SetNdivisions(505);
  hB2D[2]->GetZaxis()->SetTitleFont(42);
  hB2D[2]->GetZaxis()->SetTickLength(0.02/yFactor);

  exField->Draw();
  hB2D[2]->Draw("colz same");


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
 
  palette = (TPaletteAxis*)hB2D[2]->GetListOfFunctions()->FindObject("palette");
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
  
