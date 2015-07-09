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
#include <TMarker.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TMath.h>
#include <TPaveText.h>
#include <TExec.h>
#include <TGaxis.h>

#include "PData.hh"
#include "PlasmaGlob.hh"

void PlotField1D( const TString &sim, Int_t time, Int_t Nbins=1, const TString &options="") { 
  
#ifdef __CINT__  
  gSystem->Load("libplasma.so");
#endif

  PlasmaGlob::Initialize();
    
  TString opt = options;

  gStyle->SetPadLeftMargin(0.10);  // Margin left axis  
  gStyle->SetPadRightMargin(0.12); // Margin right axis 
  if(opt.Contains("grid")) {
    gStyle->SetPadGridX(1);
    gStyle->SetPadGridY(1);
  }

  // Load PData
  PData *pData = PData::Get(sim.Data());
  pData->LoadFileNames(time);
  if(!pData->IsInit()) return;

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
  // cout << Form(" - Bunch gamma      = %8.4f", gamma ) << endl;
  // cout << Form(" - Bunch velocity   = %8.4f c", vbeam ) << endl;

  Float_t Time = pData->GetRealTime();
  // z start of the plasma in normalized units.
  Float_t zStartPlasma = pData->GetPlasmaStart()*kp;
  // z start of the beam in normalized units.
  Float_t zStartBeam = pData->GetBeamStart()*kp;
  Time -= zStartPlasma - zStartBeam;

  // 1D histograms
  TString opth1 = opt;
  opth1 += "avg";
  // Get electric fields
  const Int_t Nfields = 1;
  TH1F **hE1D = new TH1F*[Nfields];
  for(Int_t i=0;i<Nfields;i++) {
    hE1D[i] = NULL;
    if(!pData->GetEfieldFileName(i))
      continue;
    
    char nam[3]; sprintf(nam,"e%i",i+1);
    if(ThreeD) {
      if(i==0) 
	hE1D[i] = pData->GetH1SliceZ3D(pData->GetEfieldFileName(i)->c_str(),nam,-1,Nbins,-1,Nbins,opth1.Data());
      else  
	hE1D[i] = pData->GetH1SliceZ3D(pData->GetEfieldFileName(i)->c_str(),nam,-Nbins,Nbins,-Nbins,Nbins,opth1.Data());
    } else if(CYL) { // Cylindrical: The firt bin with r>0 is actually the number 1 (not the 0).
      if(i==0) 
	hE1D[i] = pData->GetH1SliceZ(pData->GetEfieldFileName(i)->c_str(),nam,1,Nbins,opth1.Data());
      else
	hE1D[i] = pData->GetH1SliceZ(pData->GetEfieldFileName(i)->c_str(),nam,1,Nbins,opth1.Data());
    } else { // 2D cartesian
      if(i==0) 
	hE1D[i] = pData->GetH1SliceZ(pData->GetEfieldFileName(i)->c_str(),nam,-1,Nbins,opth1.Data());
      else 
	hE1D[i] = pData->GetH1SliceZ(pData->GetEfieldFileName(i)->c_str(),nam,-Nbins,Nbins,opth1.Data());
    }
    
    char hName[24];
    sprintf(hName,"hE_%i",i);
    hE1D[i]->SetName(hName);
    if(opt.Contains("comov"))
      hE1D[i]->GetXaxis()->SetTitle("#zeta [c/#omega_{p}]");
    else
      hE1D[i]->GetXaxis()->SetTitle("z [c/#omega_{p}]");
   
    if(i==0)
      hE1D[i]->GetYaxis()->SetTitle("E_{z} [E_{0}]");
    else if(i==1)
      hE1D[i]->GetYaxis()->SetTitle("E_{y} [E_{0}}]");
    else if(i==2)
      hE1D[i]->GetYaxis()->SetTitle("E_{x} [E_{0}}]");
    
    
  }  
  
  // Chaning to user units: 
  if(opt.Contains("units") && n0) {
    for(Int_t i=0;i<Nfields;i++) {
      Int_t NbinsX = hE1D[i]->GetNbinsX();
      Float_t xMin = (skindepth/PUnits::mm) * hE1D[i]->GetXaxis()->GetXmin();
      Float_t xMax = (skindepth/PUnits::mm) * hE1D[i]->GetXaxis()->GetXmax();
      
      hE1D[i]->SetBins(NbinsX,xMin,xMax);
      
      for(Int_t j=0;j<NbinsX;j++) {
	hE1D[i]->SetBinContent(j, hE1D[i]->GetBinContent(j) * ( E0 / (PUnits::GV/PUnits::m) ) );
      }
      
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
    }
  }
  
  // Calculate wave positions:
  // ----------------------------------------------------------------
  
  // Retrieve the previous time TGraph if any;
  // Open TGraph
  TString filename = Form("./%s/Plots/Field1D/Field1D-%s.root",sim.Data(),sim.Data());
  TFile * ifile = (TFile*) gROOT->GetListOfFiles()->FindObject(filename);
  // if doesn't exist the directory should be created
  if (!ifile) {
    TString f = filename;
    TString dir2 = f.Remove( f.Last( '/' ), f.Length() - f.Last( '/' ) );
    TString dir1 = f.Remove( f.Last( '/' ), f.Length() - f.Last( '/' ) );
    gSystem->mkdir( dir1 );
    gSystem->mkdir( dir2 );
    ifile = new TFile(filename,"UPDATE");
  }  
  
  
  
  Float_t *EMaxPos = new Float_t[Nfields];
  Float_t *EMaxValue = new Float_t[Nfields];
  
  // Get the Graph with the x1 positions of the maximum E_1
  TGraph **graph = new TGraph*[Nfields];
  char gName[24];
  for(Int_t i=0;i<Nfields;i++) {
    if(!hE1D[i]) continue;
    
    EMaxPos[i] = EMaxValue[i] = -999;
    
    // Initial time search window:
    Float_t xCenter = (hE1D[i]->GetXaxis()->GetXmin()+hE1D[i]->GetXaxis()->GetXmax())/2.; 
    Float_t xs1min = 0.5*(hE1D[i]->GetXaxis()->GetXmin()+xCenter);
    Float_t xs1max = xCenter;

    // For focusing fields i==1,2 we use a narrower window based on the previosly found 
    // minimum for the accelerating fields i==0.
    if(i>0) {
      xCenter = EMaxPos[0];
      if(opt.Contains("units") && n0) {
	xs1min = xCenter - (0.25*TMath::Pi() * 1e3 * skindepth);
	xs1max = xCenter + (0.25*TMath::Pi() * 1e3 * skindepth);
      } else {
	xs1min = xCenter - 0.25*TMath::Pi();
	xs1max = xCenter + 0.25*TMath::Pi();
      }
    }

    sprintf(gName,"gEMaxPos_%i",i);
    graph[i] =  (TGraph*) ifile->Get(gName);
    if(graph[i]) {
      Double_t *y = graph[i]->GetY();
      
      // Setup the searching windows to +/- pi/2 respect the last found minimum.
      if(opt.Contains("units") && n0) {
	xs1min = y[graph[i]->GetN()-1] - (0.5*TMath::Pi() * 1e3 * skindepth);
	xs1max = y[graph[i]->GetN()-1] + (0.5*TMath::Pi() * 1e3 * skindepth);
      } else {
	xs1min = y[graph[i]->GetN()-1] - 0.5*TMath::Pi();
	xs1max = y[graph[i]->GetN()-1] + 0.5*TMath::Pi();
      }
      
      delete graph[i];
    }
    
    TH1F *htemp = (TH1F*) hE1D[i]->Clone("htemp");
    htemp->GetXaxis()->SetRangeUser(xs1min,xs1max);
    
    //htemp->Smooth(1,"R");
    Int_t binMax = htemp->GetMinimumBin();
    EMaxPos[i] = (Float_t) htemp->GetBinCenter(binMax);
    EMaxValue[i] = (Float_t) htemp->GetBinContent(binMax);
    delete htemp;
    
  }
  
  
  // Tunning the Histograms
  // ---------------------
  
  
  // Plotting
  // -----------------------------------------------
  
  // Output file
  TString fOutName = Form("./%s/Plots/Field1D/Field1D",sim.Data());
  fOutName += Form("-%s_%i",sim.Data(),time);
  
  // Canvas setup
  TCanvas *C = new TCanvas("C","Electric wakefield on axis",850,1000);
  C->Divide(1,2);

  // Draw objects
  TPaveText *textTime = new TPaveText(0.70,0.87,0.85,0.92,"NDC");
  PlasmaGlob::SetPaveTextStyle(textTime); 
  char ctext[128];
  if(opt.Contains("units") && n0) 
    sprintf(ctext,"Z = %5.1f mm", 1e3 * skindepth * Time);
  else
    sprintf(ctext,"T = %5.1f 1/#omega_{p}",Time);
  textTime->AddText(ctext);
  
  // Colors
  Int_t fieldC  = PlasmaGlob::fieldLine;
  Int_t phaseC  = kGray+1;
 
  // Actual Plotting!
  // ------------------------------------------------------------

  // More makeup
  C->cd(1);
  gPad->SetGridy(0);
  gPad->SetGridx(0);
  gPad->SetFrameLineWidth(2);  

  hE1D[0]->SetLineWidth(1);
  hE1D[0]->GetYaxis()->CenterTitle();
  hE1D[0]->GetXaxis()->CenterTitle();
  hE1D[0]->SetLineStyle(1);
  hE1D[0]->SetLineWidth(3);
  hE1D[0]->SetLineColor(fieldC);
  hE1D[0]->SetMarkerStyle(20);

  if(Nfields>1) {
    hE1D[1]->GetYaxis()->CenterTitle();
    hE1D[1]->GetXaxis()->CenterTitle();
    hE1D[1]->SetLineStyle(1);
    hE1D[1]->SetLineWidth(1);  
    hE1D[1]->SetLineColor(fieldC);
    hE1D[1]->SetMarkerStyle(24);
  }

  Float_t factor = 1.5;
  Float_t minimum = factor * hE1D[0]->GetMinimum();
  Float_t maximum = factor * hE1D[0]->GetMaximum();
    
  if(Nfields>1) {
    if(hE1D[1]->GetMaximum() > hE1D[0]->GetMaximum()) {
      maximum = factor * hE1D[1]->GetMaximum();
    } 
    
    if(hE1D[1]->GetMinimum() < hE1D[0]->GetMinimum()) {
    minimum = factor * hE1D[1]->GetMinimum();
    } 
  }
  
  if( maximum >= TMath::Abs(minimum)) minimum = -maximum;
  else maximum = - minimum;
  
  hE1D[0]->GetYaxis()->SetRangeUser(minimum,maximum);
  hE1D[0]->Draw("C");
  
  if(Nfields>1)
    hE1D[1]->Draw("C same");
  
  C->Update();
  
  TLine *line0 = new TLine(hE1D[0]->GetXaxis()->GetXmin(),
			   (gPad->GetUymin()+gPad->GetUymax())/2.,
			   hE1D[0]->GetXaxis()->GetXmax(),
			   (gPad->GetUymin()+gPad->GetUymax())/2.);
  line0->SetLineColor(kGray+1);
  line0->SetLineStyle(2);
  line0->Draw();
  
  TMarker *markEMax0 = new TMarker(EMaxPos[0],EMaxValue[0], 24);
  markEMax0->SetMarkerColor(fieldC);
  markEMax0->SetMarkerSize(1.6);
  markEMax0->Draw();

  if(Nfields>1) {
    TMarker *markEMax1 = new TMarker(EMaxPos[1],EMaxValue[1], 24);
    markEMax1->SetMarkerColor(fieldC);
    markEMax1->SetMarkerSize(1.4);
    markEMax1->Draw();
  }
  
  textTime->Draw();

  // ----
  
  // Define the TGraphs
  Int_t nPoints = 0;  
  TGraph **gEMaxPos = new TGraph*[Nfields];
  TGraph **gEMaxValue = new TGraph*[Nfields];
  for(Int_t i=0;i<Nfields;i++) {
    if(!hE1D[i]) continue;

    sprintf(gName,"gEMaxPos_%i",i);
    gEMaxPos[i] = (TGraph*) ifile->Get(gName);
    if(gEMaxPos[i]==NULL) {
      gEMaxPos[i] = new TGraph();
      gEMaxPos[i]->SetName(gName);
    } else {
      nPoints = gEMaxPos[i]->GetN(); 
    }  
    gEMaxPos[i]->Set(nPoints+1);
    if(opt.Contains("units") && n0) 
      gEMaxPos[i]->SetPoint(nPoints, 1e3 * skindepth * Time,EMaxPos[i]);
    else
      gEMaxPos[i]->SetPoint(nPoints,Time,EMaxPos[i]);
    
    if(opt.Contains("units") && n0) {
      gEMaxPos[i]->GetYaxis()->SetTitle("#zeta_{min} [mm]");
      gEMaxPos[i]->GetXaxis()->SetTitle("Z [mm]");
    } else {
      gEMaxPos[i]->GetYaxis()->SetTitle("#zeta_{min} [c/#omega_{p}]");
      gEMaxPos[i]->GetXaxis()->SetTitle("T [c/#omega_{p}]");
    }

    gEMaxPos[i]->Write(gEMaxPos[i]->GetName(),TObject::kOverwrite);

    sprintf(gName,"gEMaxValue_%i",i);
    gEMaxValue[i] = (TGraph*) ifile->Get(gName);
    if(gEMaxValue[i]==NULL) {
      gEMaxValue[i] = new TGraph();  
      gEMaxValue[i]->SetName(gName);
    } else {
      nPoints = gEMaxValue[i]->GetN(); 
    }  
    gEMaxValue[i]->Set(nPoints+1);
    if(opt.Contains("units") && n0) 
      gEMaxValue[i]->SetPoint(nPoints, 1e3 * skindepth * Time,EMaxValue[i]);
    else
      gEMaxValue[i]->SetPoint(nPoints,Time,EMaxValue[i]);
    
    if(opt.Contains("units") && n0) {
      gEMaxValue[i]->GetYaxis()->SetTitle("E_{min} [GV/m]");
      gEMaxValue[i]->GetXaxis()->SetTitle("Z [mm]");
    } else {
      gEMaxValue[i]->GetYaxis()->SetTitle("E_{min} [E_{0}]");
      gEMaxValue[i]->GetXaxis()->SetTitle("T [c/#omega_{p}]");
    }
    
    gEMaxValue[i]->Write(gEMaxValue[i]->GetName(),TObject::kOverwrite);    
  }


  C->cd(2);
  gPad->SetGridy(1);
  gPad->SetGridx(0);
  gPad->SetFrameLineWidth(2);  

  Float_t minPhase = 99.;
  Float_t maxPhase = -99.;
  Float_t minField = 99.;
  Float_t maxField = -99.;

  Double_t *yEMaxPos[Nfields];
  Double_t *yEMaxValue[Nfields];

  for(Int_t i=0;i<Nfields;i++) {
    yEMaxPos[i]   = gEMaxPos[i]->GetY();
    yEMaxValue[i] = gEMaxValue[i]->GetY();
    
    for(Int_t j=0;j<gEMaxPos[0]->GetN();j++) {
      if(yEMaxPos[i][j]>maxPhase)
	maxPhase = yEMaxPos[i][j];
      if(yEMaxPos[i][j]<minPhase)
	minPhase = yEMaxPos[i][j];
 
      if(yEMaxValue[i][j]>maxField)
	maxField = yEMaxValue[i][j];
      if(yEMaxValue[i][j]<minField)
	minField = yEMaxValue[i][j];
    }
  }

  Float_t margin = (maxPhase - minPhase)/10;
  gEMaxPos[0]->GetYaxis()->SetRangeUser(minPhase-margin,maxPhase+margin);
  gEMaxPos[0]->GetYaxis()->CenterTitle();
  gEMaxPos[0]->GetXaxis()->CenterTitle();
  gEMaxPos[0]->SetLineColor(phaseC);
  gEMaxPos[0]->SetMarkerColor(phaseC);
  gEMaxPos[0]->SetLineWidth(3);
  gEMaxPos[0]->SetMarkerStyle(20);
  gEMaxPos[0]->SetMarkerSize(1.4);
  gEMaxPos[0]->Draw("APC");

  if(Nfields>1) {
    gEMaxPos[1]->SetLineStyle(1);
    gEMaxPos[1]->SetLineColor(phaseC);
    gEMaxPos[1]->SetMarkerColor(phaseC);
    gEMaxPos[1]->SetLineWidth(1);
    gEMaxPos[1]->SetMarkerStyle(24);
    gEMaxPos[1]->SetMarkerSize(1.4);
    gEMaxPos[1]->Draw("PC");
  }

  // Emax value
  // New axis first:
  C->Update();  // Needed for the axis!

  margin = (maxField - minField)/10;
  if (margin==0) margin = 1; 
  Float_t rightmin = minField-margin;
  Float_t rightmax = maxField+margin;
  Float_t slope = (gPad->GetUymax() - gPad->GetUymin())/(rightmax-rightmin);
  TGaxis *axisEmax = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),gPad->GetUxmax(),
				gPad->GetUymax(),rightmin,rightmax,505,"+L");
  axisEmax->SetLineWidth(1);
  axisEmax->SetLineColor(fieldC);
  axisEmax->SetLabelColor(fieldC);
  axisEmax->SetTitleColor(fieldC);
  if(opt.Contains("units") && n0) 
    axisEmax->SetTitle("E_{min} [GV/m]");
  else
    axisEmax->SetTitle("E_{min} [E_{0}]");
  axisEmax->CenterTitle();
  axisEmax->SetTitleSize(0.05);
  axisEmax->SetTitleOffset(1.2);
  axisEmax->SetLabelSize(0.05);
  axisEmax->SetLabelOffset(0.006);
  
  axisEmax->Draw();
  
  // Adjust the TGraph
  Double_t *x = gEMaxValue[0]->GetX();
  Double_t *y = gEMaxValue[0]->GetY();
  for(Int_t i=0;i<gEMaxValue[0]->GetN();i++) {
    gEMaxValue[0]->SetPoint(i,x[i],(y[i]-rightmin)*slope + gPad->GetUymin());
  }  
  gEMaxValue[0]->SetLineColor(fieldC);
  gEMaxValue[0]->SetMarkerColor(fieldC);
  gEMaxValue[0]->SetLineWidth(3);
  gEMaxValue[0]->SetMarkerStyle(20);
  gEMaxValue[0]->SetMarkerSize(1.4);
  gEMaxValue[0]->Draw("PC");

   if(Nfields>1) {
    x = gEMaxValue[1]->GetX();
    y = gEMaxValue[1]->GetY();
    for(Int_t i=0;i<gEMaxValue[1]->GetN();i++) {
      gEMaxValue[1]->SetPoint(i,x[i],(y[i]-rightmin)*slope + gPad->GetUymin());
    }  
    gEMaxValue[1]->SetLineColor(fieldC);
    gEMaxValue[1]->SetMarkerColor(fieldC);
    gEMaxValue[1]->SetLineWidth(1);
    gEMaxValue[1]->SetMarkerStyle(24);
    gEMaxValue[1]->SetMarkerSize(1.4);
    gEMaxValue[1]->Draw("PC");
  }

  // Emax value
  // New axis first:
  C->Update();  // Needed for the axis!

  C->cd();

  ifile->Close();
  
  // Print to a file
  PlasmaGlob::imgconv(C,fOutName,opt);
  // ---------------------------------------------------------
  
  }
  
