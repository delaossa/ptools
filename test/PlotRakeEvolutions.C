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
#include <TPaletteAxis.h>

#include "PData.hh"
#include "PlasmaGlob.hh"
#include "PPalette.hh"

void PlotRakeEvolutions(const TString &sim, const TString &options="png") { 
  
#ifdef __CINT__  
  gSystem->Load("libplasma.so");
#endif

  PlasmaGlob::Initialize();
  
  // Palettes!
  gROOT->Macro("PlasmaPalettes.C");

  TString opt = options;
  
  gStyle->SetPadTickX(0);
  gStyle->SetPadTickY(0);
  
  if(opt.Contains("grid")) {
    gStyle->SetPadGridX(1);
    gStyle->SetPadGridY(1);
  }

  const Int_t Nspaces = 4;
  TString phaname[Nspaces] = {"p1x1","p2x2","p3x3","x2x1"};
  TGraph *gXmean[Nspaces]; 
  TGraph *gYmean[Nspaces]; 
  TGraph *gXrms[Nspaces]; 
  TGraph *gYrms[Nspaces]; 
  TGraph *gEmit[Nspaces];
  TGraph *gCharge = NULL;
  // Special graph with an Energy spread band:
  TGraphErrors *gEneRms = NULL;
  TGraphErrors *gYposRms = NULL;
  TGraphErrors *gPyRms = NULL;
 
  Float_t maxEmit[Nspaces] = { -999., -999., -999., -999.};
  Float_t minEmit[Nspaces] = { 999., 999.,999., 999.};
  Float_t maxXmean[Nspaces] = { -999., -999., -999., -999.};
  Float_t minXmean[Nspaces] = { 999., 999., 999., 999.};
  Float_t maxXrms[Nspaces] = { -999., -999., -999., -999.};
  Float_t minXrms[Nspaces] = { 999., 999., 999., 999.};
  Float_t maxYmean[Nspaces] = { -999., -999., -999., -999.};
  Float_t minYmean[Nspaces] = { 999., 999., 999.};
  Float_t maxYrms[Nspaces] = { -999., -999., -999., -999.};
  Float_t minYrms[Nspaces] = { 999., 999., 999., 999.};
  Float_t maxCharge =-999.;
  Float_t minCharge = 999.;
  
  Float_t XmeanShift = -30.775;   // um in comoving var.

  

  for(Int_t i=0;i<Nspaces;i++) {
    TString filename;
    filename = Form("./%s/Plots/EmittanceEvolution/Evolutions-%s-%s.root",sim.Data(),sim.Data(),phaname[i].Data());
    
    TFile  *ifile = (TFile*) gROOT->GetListOfFiles()->FindObject(filename.Data());
    if (!ifile) ifile = new TFile(filename,"READ");
    
    
    gEmit[i] = (TGraph*) ifile->Get("gEmitvsTime");  
    gXmean[i] = (TGraph*) ifile->Get("gXmeanvsTime");
    gXrms[i] = (TGraph*) ifile->Get("gXrmsvsTime");
    gYmean[i] = (TGraph*) ifile->Get("gYmeanvsTime");
    gYrms[i] = (TGraph*) ifile->Get("gYrmsvsTime");

    // Energy spread
    if(i==0) {
      gCharge = (TGraph*) ifile->Get("gChargevsTime");  
      Int_t Npoints = gCharge->GetN();
      Double_t *yCharge = gCharge->GetY();
      for(Int_t j=0;j<Npoints;j++) {
	if(yCharge[j]>maxCharge)
	  maxCharge = yCharge[j];
	if(yCharge[j]<minCharge)
	  minCharge = yCharge[j];
      }
      Double_t *xValues = gYmean[i]->GetX();
      Double_t *yMean = gYmean[i]->GetY();
      Double_t *yRms = gYrms[i]->GetY();
      Npoints =  gYmean[i]->GetN();

      gEneRms = new TGraphErrors(Npoints,xValues,yMean,0,yRms);
      // cout << "NPOints = " << gEneRms->GetN() << endl;
      // for(Int_t j=0;j<gEneRms->GetN();j++) {
      // 	gEneRms->SetPointError(j,0.0,yRms[j]);
      // 	cout << "eoooo" << endl;
      // }
    } else if(i==1) {
      Int_t Npoints = gYmean[i]->GetN();
      Double_t *xValues = gXmean[i]->GetX();
      Double_t *xMean = gXmean[i]->GetY();
      Double_t *xRms = gXrms[i]->GetY();
      Double_t *yMean = gYmean[i]->GetY();
      Double_t *yRms = gYrms[i]->GetY();
  
      gPyRms   = new TGraphErrors(Npoints,xValues,yMean,0,yRms);
      gYposRms = new TGraphErrors(Npoints,xValues,xMean,0,xRms);

    } 
    // Calculate the max and min of every set of graphs:
    
    Int_t Npoints = gEmit[i]->GetN();
    Double_t *yEmit = gEmit[i]->GetY();
    for(Int_t j=0;j<Npoints;j++) {
      if(yEmit[j]>maxEmit[i])
	maxEmit[i] = yEmit[j];
      if(yEmit[j]<minEmit[i])
	minEmit[i] = yEmit[j];
    }

    Npoints = gXmean[i]->GetN();
    Double_t *yXmean = gXmean[i]->GetY();
    for(Int_t j=0;j<Npoints;j++) {
      if(i==0) {
	yXmean[j] -= XmeanShift;
      }
      
      if(yXmean[j]>maxXmean[i])
	maxXmean[i] = yXmean[j];
      if(yXmean[j]<minXmean[i])
	minXmean[i] = yXmean[j];

    }

    Npoints = gXrms[i]->GetN();
    Double_t *yXrms = gXrms[i]->GetY();
    for(Int_t j=0;j<Npoints;j++) {
      if(yXrms[j]>maxXrms[i])
	maxXrms[i] = yXrms[j];
      if(yXrms[j]<minXrms[i])
	minXrms[i] = yXrms[j];
    }

    Npoints = gYmean[i]->GetN();
    Double_t *yYmean = gYmean[i]->GetY();
    for(Int_t j=0;j<Npoints;j++) {
      if(yYmean[j]>maxYmean[i])
	maxYmean[i] = yYmean[j];
      if(yYmean[j]<minYmean[i])
	minYmean[i] = yYmean[j];
    }

    Npoints = gYrms[i]->GetN();
    Double_t *yYrms = gYrms[i]->GetY();
    for(Int_t j=0;j<Npoints;j++) {
      if(yYrms[j]>maxYrms[i])
	maxYrms[i] = yYrms[j];
      if(yYrms[j]<minYrms[i])
	minYrms[i] = yYrms[j];
    }
    
  }

  // --------------------

  // Canvas setup

  // Resolution:
  Int_t sizex = 800;
  Int_t sizey = 600;
  if(opt.Contains("hres")) {
    Int_t sizex = 1600;
    Int_t sizey = 1200;    
  }

  TCanvas *C1 = new TCanvas("C1","Evolution of Emittance",sizex,sizey);
  C1->cd();
  
  // C1->Divide(1,Nspaces);
  // PlasmaGlob::CanvasPartition(C1,pad,Nspaces);
  
  const Int_t Npads = 3;
  TString sLabels[Npads] = {"(c)","(b)","(a)"};
 
  // Text objects
  TPaveText **textLabel = new TPaveText*[Npads];

  TPad **pad = new TPad*[Npads];
  TH1F *hFrame[Npads];
  TGaxis *axis[Npads];

  // Setup Pad layout:
  Double_t lMargin = 0.15;
  Double_t rMargin = 0.15;
  Double_t bMargin = 0.12;
  Double_t tMargin = 0.04;
  Double_t vSpacing = 0.0; 
  Double_t hStep = (1.-lMargin-rMargin);
  Double_t vStep = (1.- bMargin - tMargin - (Npads-1) * vSpacing) / Npads;
  
  Float_t vposd = 0.0;
  Float_t vposu = 0.0;
  Float_t vmard = 0.0;
  Float_t vmaru = 0.0;
  Float_t vfactor = 0.0;
  Float_t hposl = 0.0;
  Float_t hposr = 1.0;
  Float_t hmarl = lMargin;
  Float_t hmarr = rMargin;
  Float_t hfactor = 1.0;
  
  for(Int_t i=0;i<Npads;i++) {
    
    if(i==0) {
      vposd = 0.0;
      vposu = bMargin + vStep;
      vfactor = vposu-vposd;  
      vmard = bMargin / vfactor;
      vmaru = 0.0;
    } else if(i == Npads-1) {
      vposd = vposu + vSpacing;
      vposu = vposd + vStep + tMargin;
      vfactor = vposu-vposd;   
      vmard = 0.0;
      vmaru = tMargin / (vposu-vposd);
    } else {
      vposd = vposu + vSpacing;
      vposu = vposd + vStep; 
      vfactor = vposu-vposd;
      vmard = 0.0;
      vmaru = 0.0;
    } 
    hfactor = hposl-hposr;
    
    C1->cd();

    char name[16];
    sprintf(name,"pad_%i",i);
    pad[i] = new TPad(name,"",hposl,vposd,hposr,vposu);
    pad[i]->SetLeftMargin(hmarl);
    pad[i]->SetRightMargin(hmarr);  
    pad[i]->SetBottomMargin(vmard);
    pad[i]->SetTopMargin(vmaru);

    pad[i]->SetFrameLineWidth(3);  
    if(opt.Contains("gridx")) {
      pad[i]->SetGridx(1);
    } else if(opt.Contains("gridy")){
      pad[i]->SetGridy(1);
    }
   
    pad[i]->Draw();
    
    pad[i]->cd();
  
    sprintf(name,"hFrame_%i",i);  
    hFrame[i] = new TH1F(name,"",100,-2.999,84.999);

    hFrame[i]->GetXaxis()->CenterTitle();
    hFrame[i]->GetYaxis()->CenterTitle();
    hFrame[i]->GetZaxis()->CenterTitle();
    hFrame[i]->SetLabelFont(42,"xyz");
    hFrame[i]->SetTitleFont(42,"xyz");

    hFrame[i]->SetNdivisions(505,"xyz");

    hFrame[i]->SetTickLength(0.06,"xz");
    hFrame[i]->SetTickLength(0.06*vfactor,"y");
   

    hFrame[i]->GetYaxis()->SetLabelSize(0.04/vfactor);
    hFrame[i]->GetYaxis()->SetLabelOffset(0.02);
    hFrame[i]->GetYaxis()->SetTitleSize(0.045/vfactor);
    hFrame[i]->GetYaxis()->SetTitleOffset(0.1);
       
    hFrame[i]->GetXaxis()->SetLabelSize(0.10);
    hFrame[i]->GetXaxis()->SetLabelOffset(0.02);
    hFrame[i]->GetXaxis()->SetTitleSize(0.14);
    hFrame[i]->GetXaxis()->SetTitleOffset(0.9);
    
    hFrame[i]->Draw("axis");    

    Int_t lineColor = kOrange+8;

    Float_t y1 = pad[i]->GetBottomMargin();
    Float_t y2 = 1 - pad[i]->GetTopMargin();
    Float_t x1 = pad[i]->GetLeftMargin();
    Float_t x2 = 1 - pad[i]->GetRightMargin();

    // textLabel[i] = new TPaveText(x2+0.02,y1+(0.02/vfactor),x1+0.07,y1+(0.06/vfactor),"NDC"); 
    textLabel[i] = new TPaveText(x2-0.08,y1+(0.02/vfactor),x2-0.02,y1+(0.08/vfactor),"NDC"); 
    PlasmaGlob::SetPaveTextStyle(textLabel[i],32); 
    textLabel[i]->SetTextFont(42);
    textLabel[i]->AddText(sLabels[i]);
    textLabel[i]->Draw();
    

    if(i==Npads-1) {
      hFrame[i]->GetYaxis()->SetTitleOffset(0.4);

      Float_t yMin = minYmean[0] - (maxYmean[0]-minYmean[0])*0.2;
      Float_t yMax = maxYmean[0] + (maxYmean[0]-minYmean[0])*0.2;
      hFrame[i]->GetYaxis()->SetRangeUser(yMin,yMax);
      hFrame[i]->GetXaxis()->SetTitle("z [mm]");
      hFrame[i]->GetYaxis()->SetTitle("#LTp_{z}#GT [GeV/c]");

      //    PlasmaGlob::SetH1LabelSize(hFrame[i]);
      
      gEneRms->SetFillColor(kGray);
      gEneRms->Draw("3");
      gYmean[0]->SetLineColor(lineColor);//PlasmaGlob::elecLine);
      gYmean[0]->Draw("C");  

      // Charge on right axis:
      //   Float_t yrMin = minCharge - (maxCharge-minCharge)*0.2;
      Float_t yrMin = 0.0001;
      Float_t yrMax = maxCharge + (maxCharge-minCharge)*0.2;
      Float_t slope = (yMax-yMin)/(yrMax-yrMin);
      
      Double_t *x = gCharge->GetX();
      Double_t *y = gCharge->GetY();
      for(Int_t j=0;j<gCharge->GetN();j++) {
      	gCharge->SetPoint(j,x[j],slope*(y[j]-yrMin)+yMin);
      }

      //hFrame[i]->GetYaxis()->SetRangeUser(yrMin,yrMax);
      gCharge->SetLineStyle(1);
      gCharge->SetLineWidth(2);
      gCharge->SetLineColor(kGray+2);      
      gCharge->Draw("C");

      pad[i]->Update();

      // Charge axis
      axis[i] = new TGaxis(pad[i]->GetUxmax(),pad[i]->GetUymin(),pad[i]->GetUxmax(),pad[i]->GetUymax(),yrMin,yrMax,505,"+LS");
      
      axis[i]->SetLineWidth(1);
      axis[i]->SetLineColor(kGray+2);//PlasmaGlob::elecLine);
      axis[i]->SetLabelColor(kGray+2);//PlasmaGlob::elecLine);
      axis[i]->SetLabelSize(0.04/vfactor);
      axis[i]->SetLabelOffset(0.015);
      axis[i]->SetLabelFont(42);
      axis[i]->SetTitleSize(0.04/vfactor);
      axis[i]->SetTitleOffset(0.4);
      axis[i]->SetTitleFont(42);
      axis[i]->SetTitle("charge [pC]");
      axis[i]->CenterTitle();
      axis[i]->SetTickSize(0.01);
      axis[i]->SetTitleColor(kGray+2);//PlasmaGlob::elecLine);
      
      axis[i]->Draw();


      pad[i]->Update();

    } else if(i==1) {
      hFrame[i]->GetYaxis()->SetTitleOffset(0.35);

      Float_t yMin,yMax;
      yMin = minYmean[1];
      yMax = maxYmean[1];
      if(fabs(yMin)>fabs(yMax)) // centering
	yMax = -yMin;
      else
	yMin = -yMax;
      Float_t yDist = yMax - yMin;
      yMin -= 0.4*yDist; 
      yMax += 0.4*yDist;
      hFrame[i]->GetYaxis()->SetRangeUser(yMin,yMax);
      hFrame[i]->GetXaxis()->SetTitle("z [mm]");
      hFrame[i]->GetYaxis()->SetTitle("#LTp_{y}#GT [MeV/c]");
      //    PlasmaGlob::SetH1LabelSize(hFrame[i]);

      Int_t grayFill = TColor::GetColor(200,200,200);
      Int_t grayFill2 = TColor::GetColor(170,170,170);
      gYposRms->SetFillColor(grayFill);//);
      gYposRms->Draw("3");
      gPyRms->SetFillColor(grayFill2);//kGray);
      //gPyRms->SetFillColor(kOrange-4);
      gPyRms->Draw("3");
      
      // Right axis:
      Float_t yrMin = minXmean[1] - (maxXmean[1]-minXmean[1])*0.2;
      Float_t yrMax = maxXmean[1] + (maxXmean[1]-minXmean[1])*0.2;
      if(fabs(yrMin)>fabs(yrMax)) // centering
	yrMax = -yrMin;
      else
	yrMin = -yrMax;
	
      Float_t slope = (yMax-yMin)/(yrMax-yrMin);
      
      Double_t *x = gXmean[1]->GetX();
      Double_t *y = gXmean[1]->GetY();
      for(Int_t j=0;j<gXmean[1]->GetN();j++) {
      	gXmean[1]->SetPoint(j,x[j],slope*(y[j]-yrMin)+yMin);
       }
   
      x = gYposRms->GetX();
      y = gYposRms->GetY();
      Double_t *ex = gYposRms->GetEX();
      Double_t *ey = gYposRms->GetEY();
      for(Int_t j=0;j<gYposRms->GetN();j++) {
  	gYposRms->SetPoint(j,x[j],slope*(y[j]-yrMin)+yMin);
      	gYposRms->SetPointError(j,ex[j],slope*(ey[j]-yrMin)+yMin);
      }   
       //hFrame[i]->GetYaxis()->SetRangeUser(yrMin,yrMax);
      gXmean[1]->SetLineStyle(1);
      gXmean[1]->SetLineWidth(2);
      gXmean[1]->SetLineColor(kGray+2);//PlasmaGlob::elecLine);      
      gXmean[1]->Draw("L");
   
      gYmean[1]->SetLineStyle(1);
      gYmean[1]->SetLineWidth(2);
      gYmean[1]->SetLineColor(lineColor);//PlasmaGlob::elecLine);
      gYmean[1]->Draw("L");
  
      pad[i]->Update();

      // Right axis
      axis[i] = new TGaxis(pad[i]->GetUxmax(),pad[i]->GetUymin(),pad[i]->GetUxmax(),pad[i]->GetUymax(),yrMin,yrMax,505,"+LS");
      
      axis[i]->SetLineWidth(1);
      axis[i]->SetLineColor(kGray+2);//PlasmaGlob::elecLine);
      axis[i]->SetLabelColor(kGray+2);//PlasmaGlob::elecLine);
      axis[i]->SetLabelSize(0.04/vfactor);
      axis[i]->SetLabelOffset(0.015);
      axis[i]->SetLabelFont(42);
      axis[i]->SetTitleSize(0.04/vfactor);
      axis[i]->SetTitleOffset(0.35);
      axis[i]->SetTitleFont(42);
      axis[i]->SetTitle("#LTy#GT [#mum]");
      axis[i]->CenterTitle();
      axis[i]->SetTickSize(0.01);
      axis[i]->SetTitleColor(kGray+2);//PlasmaGlob::elecLine);
      
      axis[i]->Draw();
      pad[i]->RedrawAxis(); 


    } else if(i==0) {
      hFrame[i]->GetYaxis()->SetTitleOffset(0.5);

      Float_t yMin,yMax;

      yMin = minEmit[1];
      yMax = maxEmit[1];
      
      Float_t yDist = yMax - yMin;
      yMin -= 0.2*yDist; 
      yMax += 0.2*yDist;
      hFrame[i]->GetYaxis()->SetRangeUser(yMin,yMax);
      hFrame[i]->GetXaxis()->SetTitle("z [mm]");
      hFrame[i]->GetYaxis()->SetTitle("#epsilon_{y} [#mum]");
      // PlasmaGlob::SetH1LabelSize(hFrame[i]);
      
      gEmit[1]->SetLineStyle(1);
      gEmit[1]->SetLineWidth(2);
      gEmit[1]->SetLineColor(lineColor);//PlasmaGlob::elecLine);
      gEmit[1]->Draw("L");
 
      pad[i]->Update();
      pad[i]->RedrawAxis(); 

    } 
  }
    
  C1->cd();
  
  // Print to a file
  // Output file
  TString fOutName = Form("./%s/Plots/RakeEvolution/RakeEvolution-%s",sim.Data(),sim.Data());
  PlasmaGlob::imgconv(C1,fOutName,opt);
  
}
