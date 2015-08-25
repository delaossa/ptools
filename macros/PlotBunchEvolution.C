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
#include "PGlobals.hh"
#include "PPalette.hh"

void PlotBunchEvolution(const TString &sim, Int_t index = 2, const TString &options="png") { 
  
#ifdef __CINT__  
  gSystem->Load("libptools.so");
#endif

  PGlobals::Initialize();
  
  // Palettes!
  gROOT->Macro("PPalettes.C");

  TString opt = options;
  
  gStyle->SetPadTickX(0);
  gStyle->SetPadTickY(0);
  
  if(opt.Contains("grid")) {
    gStyle->SetPadGridX(1);
    gStyle->SetPadGridY(1);
  }

  // Load PData
  PData *pData = PData::Get(sim.Data());

  pData->LoadFileNames(210);
  if(!pData->IsInit()) return;

  Float_t maxEmitx = -999.;
  Float_t minEmitx = 999.;
  Float_t maxBeta = -999.;
  Float_t minBeta = 999.;
  Float_t maxZmean = -999.;
  Float_t minZmean = 999.;
  Float_t maxZrms = -999.;
  Float_t minZrms = 999.;
  Float_t maxPzmean = -999.;
  Float_t minPzmean = 999.;
  Float_t maxPzrms = -999.;
  Float_t minPzrms = 999.;
  Float_t maxXmean = -999.;
  Float_t minXmean = 999.;
  Float_t maxXrms = -999.;
  Float_t minXrms = 999.;
  Float_t maxPxmean = -999.;
  Float_t minPxmean = 999.;
  Float_t maxPxrms = -999.;
  Float_t minPxrms = 999.;
  Float_t maxXpmean = -999.;
  Float_t minXpmean = 999.;
  Float_t maxXprms = -999.;
  Float_t minXprms = 999.;
  Float_t maxCharge =-999.;
  Float_t minCharge = 999.;
  Float_t maxGrms =-999.;
  Float_t minGrms = 999.;
  Float_t maxTime =-999.;
  Float_t minTime = 999.;
  
  TString filename = Form("./%s/Plots/Bunch/%s/Bunch-Evolution-%s.root",sim.Data(),pData->GetSpeciesName(index).c_str(),sim.Data());
  
  cout << filename << endl;
    
  TFile  *ifile = (TFile*) gROOT->GetListOfFiles()->FindObject(filename.Data());
  if (!ifile) ifile = new TFile(filename,"READ");

      
  TGraph *gZmean = (TGraph*) ifile->Get("gZmeanvsTime");
  TGraph *gZrms = (TGraph*) ifile->Get("gZrmsvsTime");
  TGraph *gPzmean = (TGraph*) ifile->Get("gPzmeanvsTime");
  TGraph *gPzrms = (TGraph*) ifile->Get("gPzrmsvsTime");
  TGraph *gXmean = (TGraph*) ifile->Get("gXmeanvsTime");
  TGraph *gXrms = (TGraph*) ifile->Get("gXrmsvsTime");
  TGraph *gPxmean = (TGraph*) ifile->Get("gPxmeanvsTime");
  TGraph *gPxrms = (TGraph*) ifile->Get("gPxrmsvsTime");
  TGraph *gCharge = (TGraph*) ifile->Get("gChargevsTime");  
  TGraph *gEmitx = (TGraph*) ifile->Get("gEmitxvsTime");  
  
  Int_t Npoints =  gPzmean->GetN();
  Double_t *tValues = gPzmean->GetX();
  Double_t *zMean = gZmean->GetY();
  Double_t *zRms = gZrms->GetY();
  Double_t *pzMean = gPzmean->GetY();
  Double_t *pzRms = gPzrms->GetY();
  Double_t *xmean = gXmean->GetY();
  Double_t *xrms = gXrms->GetY();
  Double_t *pxmean = gPxmean->GetY();
  Double_t *pxrms = gPxrms->GetY();
  Double_t *charge = gCharge->GetY();
  Double_t *emitx = gEmitx->GetY();

  // Calculate other quantities
  Double_t *beta = new Double_t[Npoints];
  Double_t *xpmean = new Double_t[Npoints];
  Double_t *xprms  = new Double_t[Npoints];
  Double_t *Grms  = new Double_t[Npoints];
  for(Int_t i=0;i<Npoints;i++) {
    Double_t gamma = pzMean[i]/(PConst::ElectronMassE/PUnits::GeV);
    beta[i] = gamma * (xrms[i]*xrms[i]/emitx[i]) * PUnits::um/PUnits::mm;
    xpmean[i] = pxmean[i]/pzMean[i]; // mrad
    xprms[i] = pxrms[i]/pzMean[i];   // mrad
    Grms[i] = 100.0 * pzRms[i] / pzMean[i];
  }
  TGraph *gBeta = new TGraph(Npoints,tValues,beta);
  TGraph *gGrms =  new TGraph(Npoints,tValues,Grms);

  TGraphErrors *gPzRms   = new TGraphErrors(Npoints,tValues,pzMean,0,pzRms);
  TGraphErrors *gPxRms   = new TGraphErrors(Npoints,tValues,pxmean,0,pxrms);
  TGraphErrors *gXRms = new TGraphErrors(Npoints,tValues,xmean,0,xrms);
  
  TGraph *gXpmean =  new TGraph(Npoints,tValues,xpmean);
  TGraphErrors *gXprms = new TGraphErrors(Npoints,tValues,xpmean,0,xprms);
  
  // Calculate the max and min of every set of graphs:
  for(Int_t j=0;j<Npoints;j++) {
    if(charge[j]>maxCharge)
      maxCharge = charge[j];
    if(charge[j]<minCharge)
      minCharge = charge[j];

    if(Grms[j]>maxGrms)
      maxGrms = Grms[j];
    if(Grms[j]<minGrms)
      minGrms = Grms[j];
  
    if(emitx[j]>maxEmitx)
      maxEmitx = emitx[j];
    if(emitx[j]<minEmitx)
      minEmitx = emitx[j];

    if(beta[j]>maxBeta)
      maxBeta = beta[j];
    if(beta[j]<minBeta)
      minBeta = beta[j];
        
    if(zMean[j]>maxZmean)
      maxZmean = zMean[j];
    if(zMean[j]<minZmean)
      minZmean = zMean[j];

    if(zRms[j]>maxZrms)
      maxZrms = zRms[j];
    if(zRms[j]<minZrms)
      minZrms = zRms[j];

    if(pzMean[j]>maxPzmean)
      maxPzmean = pzMean[j];
    if(pzMean[j]<minPzmean)
      minPzmean = pzMean[j];

    if(pzRms[j]>maxPzrms)
      maxPzrms = pzRms[j];
    if(pzRms[j]<minPzrms)
      minPzrms = pzRms[j];
    
    if(xmean[j]>maxXmean)
      maxXmean = xmean[j];
    if(xmean[j]<minXmean)
      minXmean = xmean[j];
  
    if(xrms[j]>maxXrms)
      maxXrms = xrms[j];
    if(xrms[j]<minXrms)
      minXrms = xrms[j];

    if(pxmean[j]>maxPxmean)
      maxPxmean = pxmean[j];
    if(pxmean[j]<minPxmean)
      minPxmean = pxmean[j];
  
    if(pxrms[j]>maxPxrms)
      maxPxrms = pxrms[j];
    if(pxrms[j]<minPxrms)
      minPxrms = pxrms[j];

    if(xpmean[j]>maxXpmean)
      maxXpmean = xpmean[j];
    if(xpmean[j]<minXpmean)
      minXpmean = xpmean[j];
  
    if(xprms[j]>maxXprms)
      maxXprms = xprms[j];
    if(xprms[j]<minXprms)
      minXprms = xprms[j];

    if(tValues[j]>maxTime)
      maxTime = tValues[j];
    if(tValues[j]<minTime)
      minTime = tValues[j];

  }

  // --------------------

  // Canvas setup

  // Resolution:
  Int_t sizex = 800;
  Int_t sizey = 600;

  TCanvas *C = new TCanvas("C","Bunch evolution",sizex,sizey);
  C->SetFillStyle(4000);
  C->cd();
  
  const Int_t Npads = 3;
  TString sLabels[Npads] = {"(c)","(b)","(a)"};
 
  // Text objects
  // TPaveText **textLabel = new TPaveText*[Npads];

  TPad **pad = new TPad*[Npads];
  TH1F *hFrame[Npads];
  TGaxis *axis[Npads];

  // Setup Pad layout:
  Double_t lMargin = 0.15;
  Double_t rMargin = 0.15;
  Double_t bMargin = 0.17;
  Double_t tMargin = 0.04;
  Float_t vSpacing = 0.02;
  PGlobals::CanvasPartition(C,Npads,lMargin,rMargin,bMargin,tMargin,vSpacing);
  
  // Define the frames for plotting
  Int_t fonttype = 43;
  Int_t fontsize = 24;
  Int_t tfontsize = 28;
  Float_t txoffset = 2.7;
  Float_t lxoffset = 0.02;
  Float_t tyoffset = 1.1;
  Float_t lyoffset = 0.01;
  Float_t tylength = 0.02;
  Float_t txlength = 0.04;
  for(Int_t i=0;i<Npads;i++) {
    char name[16];
    sprintf(name,"pad_%i",i);
    pad[i] = (TPad*) gROOT->FindObject(name);
    pad[i]->SetFrameLineWidth(2);  
    pad[i]->SetFillStyle(4000);
    pad[i]->SetFrameFillStyle(4000);
    // pad[i]->SetTickx(1);
    // pad[i]->SetTicky(1);
    
    sprintf(name,"hFrame_%i",i);  
    hFrame[i] = (TH1F*) gROOT->FindObject(name);
    if(hFrame[i]) delete hFrame[i];
    hFrame[i] = new TH1F(name,"",100,minTime,maxTime);
    hFrame[i]->Reset();
    
    Float_t xFactor = pad[0]->GetAbsWNDC()/pad[i]->GetAbsWNDC();
    Float_t yFactor = pad[0]->GetAbsHNDC()/pad[i]->GetAbsHNDC();
    
    // Format for y axis
    hFrame[i]->GetYaxis()->SetTitleFont(fonttype);
    hFrame[i]->GetYaxis()->SetTitleSize(tfontsize);
    hFrame[i]->GetYaxis()->SetTitleOffset(tyoffset);
    hFrame[i]->GetYaxis()->SetLabelFont(fonttype);
    hFrame[i]->GetYaxis()->SetLabelSize(fontsize);
    hFrame[i]->GetYaxis()->SetLabelOffset(lyoffset);
    hFrame[i]->GetYaxis()->SetTickLength(xFactor*tylength/yFactor);
    hFrame[i]->GetYaxis()->CenterTitle();
    
    // Format for x axis
    hFrame[i]->GetXaxis()->SetTitleFont(fonttype);
    hFrame[i]->GetXaxis()->SetTitleSize(tfontsize+2);
    hFrame[i]->GetXaxis()->SetTitleOffset(txoffset);
    hFrame[i]->GetXaxis()->SetLabelFont(fonttype);
    hFrame[i]->GetXaxis()->SetLabelSize(fontsize+2);
    hFrame[i]->GetXaxis()->SetLabelOffset(lxoffset);
    hFrame[i]->GetXaxis()->CenterTitle();
    hFrame[i]->GetXaxis()->SetTickLength(yFactor*txlength/xFactor);      

    Float_t y1 = pad[i]->GetBottomMargin();
    Float_t y2 = 1 - pad[i]->GetTopMargin();
    Float_t x1 = pad[i]->GetLeftMargin();
    Float_t x2 = 1 - pad[i]->GetRightMargin();
  }
  
  Int_t leftLineColor = kGray+2;//kOrange+8;
  Int_t rightLineColor = kOrange+8;//kMagenta-2;//kGray+2;
  Int_t leftFillColor = TColor::GetColor(200,200,200);
  Int_t rightFillColor = TColor::GetColor(255,206,149);//kOrange+8;//TColor::GetColor(170,170,170);
  Int_t leftFillStyle = 0;//1001;//0;//3151;//3144;
  Int_t rightFillStyle = 0;//1001;//0;//3115;//3144;


  C->cd(0);
  pad[2]->Draw();
  pad[2]->cd(); // <---------------------------------------------- Top Plot ---------
  if(opt.Contains("logz")) {
    pad[2]->SetLogz(1);
  } else {
    pad[2]->SetLogz(0);
  }

  Float_t xFactor = pad[0]->GetAbsWNDC()/pad[2]->GetAbsWNDC();
  Float_t yFactor = pad[0]->GetAbsHNDC()/pad[2]->GetAbsHNDC();


  Float_t yMin = minPzmean - (maxPzmean-minPzmean)*0.2;
  Float_t yMax = maxPzmean + (maxPzmean-minPzmean)*0.2;
  
  hFrame[2]->GetYaxis()->SetRangeUser(yMin,yMax);
  hFrame[2]->GetYaxis()->SetTitle("#LTp_{z}#GT [GeV/c]");
  hFrame[2]->Draw("axis");    

  gPzRms->SetFillColor(kGray);
  gPzRms->SetFillStyle(leftFillStyle);
  gPzRms->Draw("3");

  gPzmean->SetLineWidth(2);
  gPzmean->SetLineColor(leftLineColor);//PGlobals::elecLine);
  gPzmean->Draw("L");  

  // // Charge on right axis:
  // //   Float_t yrMin = minCharge - (maxCharge-minCharge)*0.2;
  // Float_t yrMin = 0.0001;
  // //  Float_t yrMax = maxCharge + (maxCharge-minCharge)*0.2;
  // Float_t yrMax = maxCharge*1.2;
  // Float_t slope = (yMax-yMin)/(yrMax-yrMin);
  
  // Double_t *x = gCharge->GetX();
  // Double_t *y = gCharge->GetY();
  // for(Int_t j=0;j<gCharge->GetN();j++) {
  //   gCharge->SetPoint(j,x[j],slope*(y[j]-yrMin)+yMin);
  // }
  
  // //hFrame[i]->GetYaxis()->SetRangeUser(yrMin,yrMax);
  // gCharge->SetLineStyle(1);
  // gCharge->SetLineWidth(2);
  // gCharge->SetLineColor(rightLineColor);      
  // gCharge->Draw("L");

  // pad[2]->Update();

  // // Charge axis
  // axis[2] = new TGaxis(pad[2]->GetUxmax(),pad[2]->GetUymin(),pad[2]->GetUxmax(),pad[2]->GetUymax(),yrMin,yrMax,505,"+LS");

  // Grms on right axis:
  //   Float_t yrMin = minGrms - (maxGrms-minGrms)*0.2;
  Float_t yrMin = 0.0001;
  //  Float_t yrMax = maxGrms + (maxGrms-minGrms)*0.2;
  Float_t yrMax = maxGrms*1.2;
  Float_t slope = (yMax-yMin)/(yrMax-yrMin);
  
  Double_t *x = gGrms->GetX();
  Double_t *y = gGrms->GetY();
  for(Int_t j=0;j<gGrms->GetN();j++) {
    gGrms->SetPoint(j,x[j],slope*(y[j]-yrMin)+yMin);
  }
  
  //hFrame[i]->GetYaxis()->SetRangeUser(yrMin,yrMax);
  gGrms->SetLineStyle(1);
  gGrms->SetLineWidth(2);
  gGrms->SetLineColor(rightLineColor);      
  gGrms->Draw("L");
  
  pad[2]->Update();

  // Grms axis
  axis[2] = new TGaxis(pad[2]->GetUxmax(),pad[2]->GetUymin(),pad[2]->GetUxmax(),pad[2]->GetUymax(),yrMin,yrMax,505,"+LS");
  
  axis[2]->SetLineWidth(1);
  axis[2]->SetLineColor(rightLineColor);//PGlobals::elecLine);
  axis[2]->SetTitleColor(rightLineColor);//PGlobals::elecLine);
  axis[2]->SetTitleFont(fonttype);
  axis[2]->SetTitleSize(tfontsize);
  axis[2]->SetTitleOffset(tyoffset);
  axis[2]->SetTitle("#Delta#gamma/#gamma [%]");
  axis[2]->CenterTitle();
  axis[2]->SetLabelFont(fonttype);
  axis[2]->SetLabelColor(rightLineColor);//PGlobals::elecLine);
  axis[2]->SetLabelSize(fontsize);
  axis[2]->SetLabelOffset(lyoffset);
  axis[2]->SetTickSize(xFactor*tylength/yFactor);
  
  axis[2]->Draw();
  
  pad[2]->RedrawAxis();
    


  C->cd(0);
  pad[1]->Draw();
  pad[1]->cd(); // <---------------------------------------------- Mid Plot ---------
  if(opt.Contains("logz")) {
    pad[1]->SetLogz(1);
  } else {
    pad[1]->SetLogz(0);
  }

  xFactor = pad[0]->GetAbsWNDC()/pad[1]->GetAbsWNDC();
  yFactor = pad[0]->GetAbsHNDC()/pad[1]->GetAbsHNDC();
  
  yMin = minXpmean - (maxXpmean-minXpmean)*0.2;
  yMax = maxXpmean + (maxXpmean-minXpmean)*0.2;
  if(-yMin<yMax)
    yMin = -yMax;
  else
    yMax = -yMin;   
  if(yMax<maxXprms) {
    yMax = 1.2*maxXprms;
    yMin = -yMax;
  }

  hFrame[1]->GetYaxis()->SetRangeUser(yMin,yMax);
  hFrame[1]->GetYaxis()->SetTitle("#LTx'#GT [mrad]");
  hFrame[1]->Draw("axis");    

  // Trans. position on right axis:
  yrMin = minXmean - (maxXmean-minXmean)*0.2;
  yrMax = maxXmean + (maxXmean-minXmean)*0.2;
  if(-yrMin<yrMax)
    yrMin = -yrMax;
  else
    yrMax = -yrMin;   
  if(yrMax<maxXrms) {
    yrMax = 1.2*maxXrms;
    yrMin = -yrMax;
  }
  slope = (yMax-yMin)/(yrMax-yrMin);
  
  x = gXRms->GetX();
  y = gXRms->GetY();
  Double_t *ex = gXRms->GetEX();
  Double_t *ey = gXRms->GetEY();
  for(Int_t j=0;j<gXRms->GetN();j++) {
    gXmean->SetPoint(j,x[j],slope*(y[j]-yrMin)+yMin);
    gXRms->SetPoint(j,x[j],slope*(y[j]-yrMin)+yMin);
    gXRms->SetPointError(j,ex[j],slope*(ey[j]-yrMin)+yMin);
  }

  TLine *lineZero = new TLine(minTime,0.0,maxTime,0.0);
  lineZero->SetLineColor(kGray+2);
  lineZero->SetLineStyle(2);
  lineZero->Draw();

  gXprms->SetFillColor(leftFillColor);
  gXprms->SetFillStyle(leftFillStyle);
  gXprms->Draw("3");

  gXRms->SetFillColor(rightFillColor);      
  gXRms->SetFillStyle(rightFillStyle);
  gXRms->Draw("3");

  gXmean->SetLineColor(rightLineColor);
  gXmean->SetLineWidth(2);
  gXmean->Draw("L");

  gXpmean->SetLineWidth(2);
  gXpmean->SetLineColor(leftLineColor);
  gXpmean->Draw("L");  

  
  pad[1]->Update();

  // XRms axis
  axis[1] = new TGaxis(pad[1]->GetUxmax(),pad[1]->GetUymin(),pad[1]->GetUxmax(),pad[1]->GetUymax(),yrMin,yrMax,505,"+LS");
  
  axis[1]->SetLineWidth(1);
  axis[1]->SetLineColor(rightLineColor);//PGlobals::elecLine);
  axis[1]->SetTitleColor(rightLineColor);//PGlobals::elecLine);
  axis[1]->SetTitleFont(fonttype);
  axis[1]->SetTitleSize(tfontsize);
  axis[1]->SetTitleOffset(tyoffset);
  axis[1]->SetTitle("#LTx#GT [#mum]");
  axis[1]->CenterTitle();
  axis[1]->SetLabelFont(fonttype);
  axis[1]->SetLabelColor(rightLineColor);//PGlobals::elecLine);
  axis[1]->SetLabelSize(fontsize);
  axis[1]->SetLabelOffset(lyoffset);
  axis[1]->SetTickSize(xFactor*tylength/yFactor);
  
  axis[1]->Draw();
  
  pad[1]->RedrawAxis();

  C->cd(0);

  {
    pad[0]->Draw();
    pad[0]->cd(); // <---------------------------------------------- Bottom Plot ---------
    if(opt.Contains("logz")) {
      pad[0]->SetLogz(1);
    } else {
      pad[0]->SetLogz(0);
    }

    Float_t xFactor = pad[0]->GetAbsWNDC()/pad[0]->GetAbsWNDC();
    Float_t yFactor = pad[0]->GetAbsHNDC()/pad[0]->GetAbsHNDC();
  
    Float_t yMin = minEmitx - (maxEmitx-minEmitx)*0.2;
    Float_t yMax = maxEmitx + (maxEmitx-minEmitx)*0.2;

    if(yMin<0.0) 
      yMin = 0.00001;

    hFrame[0]->GetYaxis()->SetRangeUser(yMin,yMax);
    hFrame[0]->GetYaxis()->SetTitle("#varepsilon_{x} [#mum]");

    hFrame[0]->GetXaxis()->SetTitle("z [#mum]");
    hFrame[0]->Draw("axis");    

  
    gEmitx->SetLineWidth(2);//PGlobals::elecLine);
    gEmitx->SetLineColor(leftLineColor);//PGlobals::elecLine);
    gEmitx->Draw("L");  

    // Beta on right axis:
    Float_t yrMin = minBeta - (maxBeta-minBeta)*0.2;
    // Float_t yrMin = 0.0001;
    Float_t yrMax = maxBeta + (maxBeta-minBeta)*0.2;
    // Float_t yrMax = maxBeta*1.2;
    Float_t slope = (yMax-yMin)/(yrMax-yrMin);
  
    Double_t *x = gBeta->GetX();
    Double_t *y = gBeta->GetY();
    for(Int_t j=0;j<gBeta->GetN();j++) {
      gBeta->SetPoint(j,x[j],slope*(y[j]-yrMin)+yMin);
    }
  
    //hFrame[i]->GetYaxis()->SetRangeUser(yrMin,yrMax);
    gBeta->SetLineStyle(1);
    gBeta->SetLineWidth(2);
    gBeta->SetLineColor(rightLineColor);      
    gBeta->Draw("L");
  
    pad[0]->Update();

    // Beta axis
    axis[0] = new TGaxis(pad[0]->GetUxmax(),pad[0]->GetUymin(),pad[0]->GetUxmax(),pad[0]->GetUymax(),yrMin,yrMax,505,"+LS");
  
    axis[0]->SetLineWidth(1);
    axis[0]->SetLineColor(rightLineColor);//PGlobals::elecLine);
    axis[0]->SetTitleColor(rightLineColor);//PGlobals::elecLine);
    axis[0]->SetTitleFont(fonttype);
    axis[0]->SetTitleSize(tfontsize);
    axis[0]->SetTitleOffset(tyoffset);
    axis[0]->SetTitle("#beta_{x} [mm]");
    axis[0]->CenterTitle();
    axis[0]->SetLabelFont(fonttype);
    axis[0]->SetLabelColor(rightLineColor);//PGlobals::elecLine);
    axis[0]->SetLabelSize(fontsize);
    axis[0]->SetLabelOffset(lyoffset);
    axis[0]->SetTickSize(xFactor*tylength/yFactor);
  
    axis[0]->Draw();
  
    pad[0]->RedrawAxis();
  }

  C->cd();
  
  // Print to a file
  // Output file
  TString fOutName = Form("./%s/Plots/Bunch/%s/Bunch-Evolution-%s",sim.Data(),pData->GetSpeciesName(index).c_str(),sim.Data());
  PGlobals::imgconv(C,fOutName,opt);
  
}
