#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>

#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TCut.h>
#include <TChain.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
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

using namespace std;

void PlotBunchComp(Int_t time, const TString &options="") {

#ifdef __CINT__  
  gSystem->Load("libptools.so");
#endif

  PGlobals::Initialize();

  TString opt = options;

  const Int_t Nsim = 13;
  char sName[Nsim][36] = {"facet_v23kA.G.RI.c.n2.3D", 
  "facet_v23kA.G.RI.c.n4.3D", 
  "facet_v23kA.G.RI.c.n6.3D", 
  "facet_v23kA.G.RI.c.n8.3D", 
  "facet_v23kA.G.RI.c.n10.3D", 
  "facet_v23kA.G.RI.c.n20.3D", 
  "facet_v23kA.G.RI.c.n25.3D", 
  "facet_v23kA.G.RI.c.n30.3D", 
  "facet_v23kA.G.RI.c.n35.3D", 
  "facet_v23kA.G.RI.c.n40.3D", 
  "facet_v23kA.G.RI.c.n60.3D", 
  "facet_v23kA.G.RI.c.n80.3D", 
  "facet_v23kA.G.RI.c.n100.3D"};

  // Load first simulation data (for instance)
  PData *pData = PData::Get(sName[0]);
  pData->LoadFileNames(time);
  if(!pData->IsInit()) return;

  // Some plasma constants
  Double_t n0 = pData->GetPlasmaDensity();
  Double_t kp = pData->GetPlasmaK();
  Double_t skindepth = 1.;
  if(kp!=0.0) skindepth = 1/kp;
  Double_t E0 = pData->GetPlasmaE0();
  Double_t  lightspeed =  PConst::c_light;


  // Time in OU
  Float_t Time = pData->GetRealTime();
  // z start of the plasma in normalized units.
  Float_t zStartPlasma = pData->GetPlasmaStart()*kp;
  // z start of the beam in normalized units.
  Float_t zStartBeam = pData->GetBeamStart()*kp;
  
  // Compulsory options
  opt += "comovcenter";

  // Centering time and z position:
  Double_t shiftz = pData->Shift(opt);
  if(opt.Contains("center")) {
    Time -= zStartPlasma;
    if(opt.Contains("comov"))      // Centers on the head of the beam.
      Time += zStartBeam;
  } 

  Float_t x1Min = -7.75;
  Float_t x1Max = -7.0;
  
  // Float_t HeDen[Nsim] = {2.,4.,6.,8.,10.,20.,25.,30.,35.,40.,60.,80,100.};
  // for(Int_t i=0;i<Nsim;i++) 
  //   HeDen[i] *= 1E-3 * n0 / (1E15/PUnits::cm3);
  Float_t HeDen[Nsim] = {0.2,0.4,0.6,0.8,1.,2.,2.5,3.,3.5,4.,6.,8,10.};
 
  Float_t Imax[Nsim];
  Float_t Irms[Nsim];
  Float_t Charge[Nsim];

  Float_t Zrms[Nsim];
  Float_t Yrms[Nsim];

  Float_t Emit[Nsim];
  Float_t Pmean[Nsim];
  Float_t Prms[Nsim];
  Float_t PmeanFWHM[Nsim];
  Float_t PrmsFWHM[Nsim];
  
  TFile *sFile[Nsim];
  TString *sPath[Nsim];
  Int_t color[Nsim];

  // histos and graphs
  TH1 *hX1[Nsim];
  TH1 *hP1[Nsim];
  TH1 *hP1cut[Nsim];
  TH2 *hP1X1[Nsim];
  TH2 *hP2X2[Nsim];

  TGraph *gEmit[Nsim];
  TGraph *gErms[Nsim];
  TGraph *gYrms[Nsim];
  
  for(Int_t i=0;i<Nsim;i++) {
    
    TString filename = Form("./%s/Plots/Bunch/He-electrons/Bunch-%s_%i.root",sName[i],sName[i],time);
    sFile[i] = new TFile(filename,"READ");
    
    // Load histos and graphs
    hX1[i]   = (TH1*) sFile[i]->Get("hX1"); 
    hP1[i]   = (TH1*) sFile[i]->Get("hP1"); 
    hP1X1[i] = (TH2*) sFile[i]->Get("hP1X1"); 
    hP2X2[i] = (TH2*) sFile[i]->Get("hP2X2");
    hX1[i]->ResetStats();
    hP1[i]->ResetStats();
    hP1X1[i]->ResetStats();
    hP2X2[i]->ResetStats();
 

    gEmit[i] = (TGraph*) sFile[i]->Get("gEmitx"); 
    gYrms[i] = (TGraph*) sFile[i]->Get("gXrms"); 
    gErms[i] = (TGraph*) sFile[i]->Get("gErms"); 

    Imax[i]   = hX1[i]->GetMaximum(); // in kA .

    Float_t binsize =  hX1[i]->GetBinWidth(1) * PUnits::um / lightspeed;
    Charge[i] = hX1[i]->Integral() * binsize / PUnits::femtosecond; // pC
    Zrms[i]   = hX1[i]->GetRMS(); // um
    Irms[i]   = (Charge[i] / (TMath::Sqrt(2.0*PConst::pi) * Zrms[i] * PUnits::um / lightspeed) ) * PUnits::femtosecond;

    Pmean[i]  = hP1[i]->GetMean();
    Prms[i]   = 100*hP1[i]->GetRMS()/hP1[i]->GetMean();
    //Prms[i]   = hP1[i]->GetRMS()/PUnits::MeV;

    
    // Total relative energy spread within FWHM:
    char hName[16];
    sprintf(hName,"hP1cut_%i",i);
    hP1cut[i] = (TH1F*) hP1[i]->Clone(hName);
    hP1cut[i]->Reset();

    Float_t maxValue = hP1[i]->GetBinContent(hP1[i]->GetMaximumBin());
    Int_t   lBin = -1;
    for(Int_t k=1;k<=hP1[i]->GetNbinsX();k++) {
      Float_t binValue = hP1[i]->GetBinContent(k);
      if(binValue>maxValue/2) {
	lBin = k;
	break;
      }
    }

    Int_t rBin = -1;
    for(Int_t k=hP1[i]->GetNbinsX();k>0;k--) {
      Float_t binValue = hP1[i]->GetBinContent(k);
      if(binValue>maxValue/2) {
	rBin = k;
	break;
      }
    }
    
    for(Int_t k=lBin;k<=rBin;k++) {
      Float_t binValue = hP1[i]->GetBinContent(k);
      hP1cut[i]->SetBinContent(k,binValue);
    }
    //  hP1[i]cut->ResetStats();
    PmeanFWHM[i] = hP1cut[i]->GetMean();
    PrmsFWHM[i]  = 100*hP1cut[i]->GetRMS()/PmeanFWHM[i];
    

    
    // emittance
    Double_t xmean  = 0.0;
    Double_t ymean  = 0.0;
    Double_t x2mean = 0.0;
    Double_t y2mean = 0.0;
    Double_t xymean = 0.0;
    Double_t Ntotal = 0.0;
    
    for(Int_t k=1;k<=hP2X2[i]->GetNbinsX();k++) {
      Double_t x = hP2X2[i]->GetXaxis()->GetBinCenter(k);
      // if(x<xmin || x>xmax) continue;
      for(Int_t l=1;l<=hP2X2[i]->GetNbinsY();l++) {
	Double_t y = hP2X2[i]->GetYaxis()->GetBinCenter(l);
	// if(y<ymin || y>ymax) continue;
	Double_t value = TMath::Abs(hP2X2[i]->GetBinContent(k,l));
	xmean += x*value;
	ymean += y*value;
	x2mean += x*x*value;
	y2mean += y*y*value;
	xymean += x*y*value;
	
	Ntotal += value;
      }
    }
    
    xmean  /= Ntotal;
    ymean  /= Ntotal;
    x2mean /= Ntotal;
    y2mean /= Ntotal;
    xymean /= Ntotal;

    Double_t xrms2  = x2mean - xmean*xmean;
    Double_t yrms2  = y2mean - ymean*ymean;
    Double_t xrms   = TMath::Sqrt(xrms2);
    Double_t yrms   = TMath::Sqrt(yrms2);
    Double_t xyrms2 = xymean - xmean*ymean;

    Double_t emittance = TMath::Sqrt(xrms2*yrms2 - xyrms2*xyrms2);

    Yrms[i] = yrms;
    Emit[i] = emittance;

    // Cosmetics 
    color[i] = TColor::GetColor(100+20*i,200-20*i,200-20*i);
    
    hX1[i]->SetLineColor(color[i]);
    hP1[i]->SetLineColor(color[i]);

    gEmit[i]->SetLineColor(color[i]);
    gEmit[i]->SetMarkerColor(color[i]);
    gErms[i]->SetLineColor(color[i]);
    gErms[i]->SetMarkerColor(color[i]);
    gYrms[i]->SetLineColor(color[i]);
    gYrms[i]->SetMarkerColor(color[i]);
    
  } 

  // Graphs for He density dependence.
  TGraph *gImaxDen = new TGraph(Nsim,HeDen,Imax);
  TGraph *gIrmsDen = new TGraph(Nsim,HeDen,Irms);
  TGraph *gChargeDen = new TGraph(Nsim,HeDen,Charge);
  TGraph *gZrmsDen = new TGraph(Nsim,HeDen,Zrms);
  TGraph *gYrmsDen = new TGraph(Nsim,HeDen,Yrms);
  // TGraph *gPmeanDen = new TGraph(Nsim,HeDen,Pmean);
  // TGraph *gPrmsDen = new TGraph(Nsim,HeDen,Prms);
  TGraph *gPmeanDen = new TGraph(Nsim,HeDen,PmeanFWHM);
  TGraph *gPrmsDen = new TGraph(Nsim,HeDen,PrmsFWHM);
  TGraph *gEmitDen = new TGraph(Nsim,HeDen,Emit);
  
  gImaxDen->SetLineColor(kMagenta-2);
  gImaxDen->SetLineWidth(2);
  gImaxDen->SetMarkerColor(kMagenta-2);
  gImaxDen->SetMarkerStyle(20);
  gImaxDen->SetMarkerSize(1.2);

  gIrmsDen->SetLineColor(kMagenta-2);
  gIrmsDen->SetLineWidth(2);
  gIrmsDen->SetMarkerColor(kMagenta-2);
  gIrmsDen->SetMarkerStyle(20);
  gIrmsDen->SetMarkerSize(1.2);

  gChargeDen->SetLineColor(kMagenta-2);
  gChargeDen->SetLineWidth(2);
  gChargeDen->SetMarkerColor(kMagenta-2);
  gChargeDen->SetMarkerStyle(20);
  gChargeDen->SetMarkerSize(1.2);

  gZrmsDen->SetLineColor(kMagenta-2);
  gZrmsDen->SetLineWidth(2);
  gZrmsDen->SetMarkerColor(kMagenta-2);
  gZrmsDen->SetMarkerStyle(20);
  gZrmsDen->SetMarkerSize(1.2);

  gPrmsDen->SetLineColor(kMagenta-2);
  gPrmsDen->SetLineWidth(2);
  gPrmsDen->SetMarkerColor(kMagenta-2);
  gPrmsDen->SetMarkerStyle(20);
  gPrmsDen->SetMarkerSize(1.2);

  gPmeanDen->SetLineColor(kMagenta-2);
  gPmeanDen->SetLineWidth(2);
  gPmeanDen->SetMarkerColor(kMagenta-2);
  gPmeanDen->SetMarkerStyle(20);
  gPmeanDen->SetMarkerSize(1.2);

  gEmitDen->SetLineColor(kMagenta-2);
  gEmitDen->SetLineWidth(2);
  gEmitDen->SetMarkerColor(kMagenta-2);
  gEmitDen->SetMarkerStyle(20);
  gEmitDen->SetMarkerSize(1.2);


  // Canvas setup
  // Create the canvas and the pads before the Frame loop
  // Resolution:
  Int_t sizex = 600;
  Int_t sizey = 800;
  char cName[32];
  sprintf(cName,"C");     
  TCanvas *C = (TCanvas*) gROOT->FindObject(cName);
  if(C==NULL) C = new TCanvas("C","",sizex,sizey);
  C->cd();
  C->Clear();

  // Float_t nMin = 0.9;
  // Float_t nMax = 55.;
  Float_t nMin = 0.18;
  Float_t nMax = 11.;

  // Setup Pad layout: 
  const Int_t NPad = 6;
  TPad *pad[NPad];
  TH1F *hFrame[NPad];

  Float_t lMargin = 0.18;
  Float_t rMargin = 0.10;
  Float_t bMargin = 0.15;
  Float_t tMargin = 0.04;
  Float_t factor = 1.0;    
  PGlobals::CanvasAsymPartition(C,NPad,lMargin,rMargin,bMargin,tMargin,factor,0.01);

  // Define the frames for plotting
  Int_t fonttype = 43;
  Int_t fontsize = 18;
  Int_t tfontsize = 20;
  Float_t txoffset = 5.0;
  Float_t lxoffset = 0.02;
  Float_t tyoffset = 2.5;
  Float_t lyoffset = 0.01;
  Float_t tylength = 0.02;
  Float_t txlength = 0.04;
  for(Int_t i=0;i<NPad;i++) {
    char name[16];
    sprintf(name,"pad_%i",i);
    pad[i] = (TPad*) gROOT->FindObject(name);
    pad[i]->SetFrameLineWidth(2);  
    pad[i]->SetTickx(1);
    pad[i]->SetTicky(1);

    sprintf(name,"hFrame_%i",i);

    hFrame[i] = new TH1F(name,"",10,nMin,nMax);
    
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
    hFrame[i]->GetYaxis()->SetNdivisions(505);

    // Format for x axis
    hFrame[i]->GetXaxis()->SetTitleFont(fonttype);
    hFrame[i]->GetXaxis()->SetTitleSize(tfontsize+2);
    hFrame[i]->GetXaxis()->SetTitleOffset(txoffset);
    hFrame[i]->GetXaxis()->SetLabelFont(fonttype);
    hFrame[i]->GetXaxis()->SetLabelSize(fontsize+2);
    hFrame[i]->GetXaxis()->SetLabelOffset(lxoffset);
    hFrame[i]->GetXaxis()->CenterTitle();
    hFrame[i]->GetXaxis()->SetTickLength(yFactor*txlength/xFactor);      
  }

  Double_t Min,Max;
  Double_t mfactor = 0.3;

  TLine line0(nMin,0.,nMax,0.);
  line0.SetLineColor(kGray+2);
  line0.SetLineStyle(2);

  C->cd(0);
  pad[0]->Draw();
  pad[0]->cd();

  gPad->SetLogx(1);

  PGlobals::GetGraphMaxMin(gPmeanDen,Min,Max);
  hFrame[0]->GetYaxis()->SetRangeUser(Min-(Max-Min)*mfactor,Max+(Max-Min)*mfactor);  
  
  // hFrame[0]->GetXaxis()->SetTitle("n_{He} [10^{15} e/cm^{3}]");
  hFrame[0]->GetXaxis()->SetTitle("n_{He}/n_{0} [%]");
  hFrame[0]->GetYaxis()->SetTitle("#LTp_{z}#GT [GeV/c]");

  hFrame[0]->Draw();

  gPmeanDen->Draw("PL");

  C->cd(0);
  pad[1]->Draw();
  pad[1]->cd();

  gPad->SetLogx(1);

  PGlobals::GetGraphMaxMin(gPrmsDen,Min,Max);
  
  hFrame[1]->GetYaxis()->SetRangeUser(Min-(Max-Min)*mfactor,Max+(Max-Min)*mfactor);
  hFrame[1]->GetYaxis()->SetTitle("#Delta#gamma/#LT#gamma#GT [%]");
  hFrame[1]->Draw("AXIS");

  line0.Draw();
  gPrmsDen->Draw("PL");
  
  C->cd(0);
  pad[2]->Draw();
  pad[2]->cd();

  gPad->SetLogx(1);

  PGlobals::GetGraphMaxMin(gZrmsDen,Min,Max);
  hFrame[2]->GetYaxis()->SetRangeUser(Min-(Max-Min)*mfactor,Max+(Max-Min)*mfactor);
  hFrame[2]->GetYaxis()->SetTitle("#sigma_{z} [#mum]");
  hFrame[2]->Draw("AXIS");
  gZrmsDen->Draw("PL");

  C->cd(0);
  pad[3]->Draw();
  pad[3]->cd();

  gPad->SetLogx(1);

  PGlobals::GetGraphMaxMin(gEmitDen,Min,Max);
  hFrame[3]->GetYaxis()->SetRangeUser(Min-(Max-Min)*mfactor,Max+(Max-Min)*mfactor);
  hFrame[3]->GetYaxis()->SetTitle("#varepsilon_{n} [#mum]");
  hFrame[3]->Draw("AXIS");

  line0.Draw();
  gEmitDen->Draw("PL");

  C->cd(0);
  pad[4]->Draw();
  pad[4]->cd();

  gPad->SetLogx(1);

  PGlobals::GetGraphMaxMin(gIrmsDen,Min,Max);
  hFrame[4]->GetYaxis()->SetRangeUser(Min-(Max-Min)*mfactor,Max+(Max-Min)*mfactor);
  hFrame[4]->GetYaxis()->SetTitle("I_{rms} [kA]");
  hFrame[4]->Draw("AXIS");

  line0.Draw();
  gIrmsDen->Draw("PL");

  C->cd(0);
  pad[5]->Draw();
  pad[5]->cd();

  gPad->SetLogx(1);

  PGlobals::GetGraphMaxMin(gChargeDen,Min,Max);
  hFrame[5]->GetYaxis()->SetRangeUser(Min-(Max-Min)*mfactor,Max+(Max-Min)*mfactor);
  hFrame[5]->GetYaxis()->SetTitle("Q [pC]");
  hFrame[5]->Draw("AXIS");

  line0.Draw();
  gChargeDen->Draw("PL");


  // Print to file --------------------------------------
  
  C->cd();
  
  // Print to a file
  // Output file
  TString fOutName = Form("./BunchComp/BunchComp");
  fOutName += Form("_%i",time);
  
  PGlobals::imgconv(C,fOutName,opt);
  // ---------------------------------------------------------
  
  
}


