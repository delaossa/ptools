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
#include <TF1.h>
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

Double_t Density(Double_t *x, Double_t *par) {
  Double_t f;
  Double_t skd = par[0];
  Double_t length = (391.4 - 15.1) * skd;
  Double_t slope = 1.0/length;
  if(x[0] < -length) f = 0;
  else if(x[0] < 0.0) f = slope * (x[0]+length);
  else f = 1.0;
  return f;
}

Double_t DensityGauss(Double_t *x, Double_t *par) {
  Double_t f;
  Double_t sz = par[0];

  if(x[0] < -4 * sz) f = 0;
  else if(x[0] < 0.0) f = TMath::Exp(-x[0]*x[0]/(2*sz*sz));
  else f = 1.0;
  return f;
}

Double_t DensityTap(Double_t *x, Double_t *par) {
  Double_t f;
  Double_t Ltap = par[0];
  Double_t lambda = par[1];
  
  if(x[0] < -Ltap) f = 0;
  if(x[0] < 0.0) f = TMath::Power(1.0-(x[0]/lambda),-4.0);
  else f = 1.0;
  return f;
}


void PlotCentroidCompS2E(const TString &options="") {

#ifdef __CINT__  
  gSystem->Load("libptools.so");
#endif

  PGlobals::Initialize();

  TString opt = options;

  // const Int_t Nsim = 6;
  // char sName[Nsim][56] = {"flash-500pC-2kA-JZ-2016.06.3D", 
  // 			  "flash-500pC-2kA-JZ-2016.06.tap10.3D",
  // 			  "flash-500pC-2kA-JZ-2016.06.tap10.corr.3D",
  // 			  "flash-500pC-2kA-JZ-2016.06.tap10.R56max.3D",
  // 			  "flash-500pC-2kA-JZ-2016.06.tap10.R56zero.3D",
  // 			  "flash-500pC-2kA-JZ-2016.06.tap10.R56zero.ES.3D"};

  // char lName[Nsim][56] = {"No taper", 
  // 			  "k_{#beta,0}L = 10",
  // 			  "k_{#beta,0}L = 10 (corr.)",
  // 			  "k_{#beta,0}L = 10 (R_{56} max)",
  // 			  "k_{#beta,0}L = 10 (R_{56} zero)",
  // 			  "k_{#beta,0}L = 10 (R_{56} zero ES)"};

  const Int_t Nsim = 4;
  char sName[Nsim][56] = {//"flash-500pC-2kA-JZ-2016.06.tap10.3D",
   			  //"flash-500pC-2kA-JZ-2016.06.tap10.R56max.3D",
   			  "flash-500pC-2kA-JZ-2016.06.tap10.R56zero.3D",
  			  "flash-500pC-2kA-JZ-2016.06.tap10.R56zero.ES.3D",
 			  "flash-500pC-3kA-JZ-2016.06.tap10.R56zero.3D",
			  "flash-500pC-3kA-JZ-2016.06.tap10.R56zero.ES.3D"};

  // char lName[Nsim][56] = {//"2 kA R_{56} opt",
  //  			  //"2 kA R_{56} max",
  //  			  "2 kA R_{56} zero",
  // 			  //   			  "R_{56} zero ES (25 #mum Al)"};
  //  			  "2 kA R_{56} zero ES",
  //  			  "3 kA R_{56} zero",
  // 			  "3 kA R_{56} zero ES"};
  char lName[Nsim][56] = {"(a): 2 kA",
			  //   			  "R_{56} zero ES (25 #mum Al)"};
   			  "(b): 2 kA (emit. spoiler)",
   			  "(c): 3 kA",
			  "(d): 3 kA (emit. spoiler)"};
  
  
  // Load first simulation data (for instance)
  PData *pData = PData::Get(sName[0]);
  //if(!pData->IsInit()) return;
  
  Float_t skd = pData->GetPlasmaSkinDepth() / PUnits::mm;
  cout << Form("\n skd = %.6f mm",skd) << endl;
  
  TFile *sFile[Nsim];
  TGraph *gX[Nsim];
  TGraph *gXrms[Nsim];
  TGraph *gXup[Nsim];
  TGraph *gXdo[Nsim];
  TGraph *gY[Nsim];
  TGraph *gYrms[Nsim];
  TGraph *gYup[Nsim];
  TGraph *gYdo[Nsim];

  Int_t N[Nsim];
  Double_t *T[Nsim];
  Double_t *X[Nsim];
  Double_t *Xrms[Nsim];
  Double_t *Xup[Nsim]; 
  Double_t *Xdo[Nsim];
  Double_t *Y[Nsim];
  Double_t *Yrms[Nsim];
  Double_t *Yup[Nsim]; 
  Double_t *Ydo[Nsim];

  
  Int_t color[Nsim];

  PPalette *pal = new PPalette("pal");
  pal->SetPalette("oli");

  Double_t xMax = -99999;
  Double_t xMin = 99999;
  Double_t yMax = -99999;
  Double_t yMin = 99999;
  
  for(Int_t i=0;i<Nsim;i++) {

    Int_t index =  i * pal->GetNColors() / Nsim;
    color[i] = pal->GetColorIndex(index);

    // switch(i) {
    // case 0 :
    //   color[i] = kGray+2;
    //   break;
    // case 1 :
    //   color[i] = kRed-7;
    //   break;
    // case 2 :
    //   color[i] = kGray+3;
    //   break;
    // case 3 :
    //   color[i] = kRed;
    //   break;
    // }
    
    TString filename = Form("./%s/Plots/Bunch/beam-driver/Bunch-Evolution-%s.root",sName[i],sName[i]);
    sFile[i] = new TFile(filename,"READ");
    
    // Load histos and graphs
    gX[i]    = (TGraph*) sFile[i]->Get("gXmeanvsTime"); 
    gXrms[i] = (TGraph*) sFile[i]->Get("gXrmsvsTime");

    T[i] =  gX[i]->GetX();
    X[i] =  gX[i]->GetY();
    Xrms[i] = gXrms[i]->GetY();
    N[i] = gX[i]->GetN();

    Xup[i] = new Double_t[N[i]];
    Xdo[i] = new Double_t[N[i]];
  
    for(Int_t j=0;j<N[i];j++) {
      
      Xup[i][j] = X[i][j] + Xrms[i][j];
      Xdo[i][j] = X[i][j] - Xrms[i][j];

      if( Xup[i][j] > xMax) xMax =  Xup[i][j];
      if( Xdo[i][j] < xMin) xMin =  Xdo[i][j];
      
    }

    gXup[i] = new TGraph(N[i],T[i],Xup[i]);
    gXdo[i] = new TGraph(N[i],T[i],Xdo[i]);

    gXup[i]->SetLineWidth(1);
    gXup[i]->SetLineColor(color[i]);
    gXup[i]->SetLineStyle(3);
    
    gXdo[i]->SetLineWidth(1);
    gXdo[i]->SetLineColor(color[i]);
    gXdo[i]->SetLineStyle(3);

    gX[i]->SetLineWidth(3);
    gX[i]->SetLineColor(color[i]);

    gY[i]    = (TGraph*) sFile[i]->Get("gYmeanvsTime"); 
    gYrms[i] = (TGraph*) sFile[i]->Get("gYrmsvsTime");

    T[i] =  gY[i]->GetX();
    Y[i] =  gY[i]->GetY();
    Yrms[i] = gYrms[i]->GetY();
    N[i] = gY[i]->GetN();

    Yup[i] = new Double_t[N[i]];
    Ydo[i] = new Double_t[N[i]];
  
    for(Int_t j=0;j<N[i];j++) {
      
      Yup[i][j] = Y[i][j] + Yrms[i][j];
      Ydo[i][j] = Y[i][j] - Yrms[i][j];

      if( Yup[i][j] > yMax) yMax =  Yup[i][j];
      if( Ydo[i][j] < yMin) yMin =  Ydo[i][j];
      
    }

    gYup[i] = new TGraph(N[i],T[i],Yup[i]);
    gYdo[i] = new TGraph(N[i],T[i],Ydo[i]);

    gYup[i]->SetLineWidth(1);
    gYup[i]->SetLineColor(color[i]);
    gYup[i]->SetLineStyle(2);
    
    gYdo[i]->SetLineWidth(1);
    gYdo[i]->SetLineColor(color[i]);
    gYdo[i]->SetLineStyle(2);

    gY[i]->SetLineWidth(3);
    gY[i]->SetLineColor(color[i]);

  }
  
  // Canvas setup
  // Create the canvas and the pads before the Frame loop
  // Resolution:
  Int_t sizex = 800;
  Int_t sizey = 500;

  if(opt.Contains("xy")) {
    sizey = 800;
  }
  
  char cName[32];
  sprintf(cName,"C");     
  TCanvas *C = (TCanvas*) gROOT->FindObject(cName);
  if(C==NULL) C = new TCanvas("C","",sizex,sizey);
  C->SetFillStyle(4000);
  C->cd();
  C->Clear();

  Float_t zmin = 999.0;
  Float_t zmax = -999.0;
  for(Int_t i=0;i<Nsim;i++) {

    if(T[i][0]<zmin) zmin = T[i][0];
    if(T[i][N[i]-1]>zmax) zmax = T[i][N[i]-1];
  }

  zmax = 100;
  
  TLine *lineZero = new TLine(zmin,0.0,zmax,0.0);
  lineZero->SetLineColor(kGray+2);
  lineZero->SetLineStyle(2);
 
  // Setup Pad layout: 
  Int_t NPad = 2;

  if(opt.Contains("xy"))
    NPad = 3;
  
  TPad **pad = new TPad*[NPad];
  TH1F **hFrame = new TH1F*[NPad];
  
  Float_t lMargin = 0.12;
  Float_t rMargin = 0.04;
  Float_t bMargin = 0.15;
  Float_t tMargin = 0.04;
  Float_t factor = 0.2;    

  // Define the frames for plotting
  Int_t fonttype = 43;
  Int_t fontsize = 24;
  Int_t tfontsize = 26;
  Float_t txoffset = 1.4;
  Float_t lxoffset = 0.02;
  Float_t tyoffset = 1.0;
  Float_t lyoffset = 0.01;
  Float_t tylength = 0.015;
  Float_t txlength = 0.02;

  if(opt.Contains("xy")) {
    lMargin = 0.12;
    rMargin = 0.04;
    bMargin = 0.10;
    tMargin = 0.04;
    factor = 0.2;    

    fontsize = 24;
    tfontsize = 26;
    txoffset = 2.2;
    lxoffset = 0.02;
    tyoffset = 1.4;
    lyoffset = 0.01;
    tylength = 0.015;
    txlength = 0.02;
  }

  PGlobals::CanvasAsymPartition(C,NPad,lMargin,rMargin,bMargin,tMargin,factor,0.02);

  for(Int_t i=0;i<NPad;i++) {
    char name[16];
    sprintf(name,"pad_%i",i);
    pad[i] = (TPad*) gROOT->FindObject(name);
    pad[i]->SetFrameLineWidth(2);  
    pad[i]->SetTickx(1);
    pad[i]->SetTicky(1);

    sprintf(name,"hFrame_%i",i);

    hFrame[i] = new TH1F(name,"",10,zmin,zmax);
    
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

  Double_t mfactor = 0.1;
  if(xMax>-xMin) xMin = -xMax;
  if(-xMin>xMax) xMax = -xMin;
  if(yMax>-yMin) yMin = -yMax;
  if(-yMin>yMax) yMax = -yMin;

  Double_t xmin =  xMin-(xMax-xMin)*mfactor;
  Double_t xmax =  xMax+(xMax-xMin)*mfactor;
  Double_t xrange = xmax - xmin;
  Double_t ymin =  yMin-(yMax-yMin)*mfactor;
  Double_t ymax =  yMax+(yMax-yMin)*mfactor;
  Double_t yrange = ymax - ymin;
  Double_t zrange = zmax - zmin;

  TLegend *Leg = new TLegend(zmin + 0.02 * zrange, xmax - 0.3*xrange , zmin + 0.4 * zrange, xmax - 0.05 * xrange,"","tl");
  PGlobals::SetPaveStyle(Leg);
  Leg->SetTextAlign(12);
  Leg->SetTextColor(kGray+3);
  Leg->SetTextFont(43);
  Leg->SetTextSize(16);
  Leg->SetLineColor(1);
  Leg->SetBorderSize(0);
  Leg->SetFillColor(0);
  Leg->SetFillStyle(1001);
  Leg->SetFillStyle(0); // Hollow  

  Int_t ipad = 0;

  if(opt.Contains("xy")) {
    C->cd(0);
    pad[ipad]->Draw();
    pad[ipad]->cd();

    hFrame[ipad]->GetYaxis()->SetRangeUser(ymin,ymax);  
  
    hFrame[ipad]->GetXaxis()->SetTitle("z [mm]");
    hFrame[ipad]->GetYaxis()->SetTitle("y [#mum]");

    hFrame[ipad]->Draw("AXIS");

    lineZero->Draw();

    for(Int_t i=0;i<Nsim;i++) {
      gY[i]->Draw("L");
      gYup[i]->Draw("L");
      gYdo[i]->Draw("L");
    }
    ipad++;
  }
  
  C->cd(0);
  pad[ipad]->Draw();
  pad[ipad]->cd();

  hFrame[ipad]->GetYaxis()->SetRangeUser(xmin,xmax);
  
  hFrame[ipad]->GetXaxis()->SetTitle("z [mm]");
  hFrame[ipad]->GetYaxis()->SetTitle("x [#mum]");

  hFrame[ipad]->Draw("AXIS");

  lineZero->Draw();

  for(Int_t i=0;i<Nsim;i++) {
    gX[i]->Draw("L");
    gXup[i]->Draw("L");
    gXdo[i]->Draw("L");

    Leg->AddEntry(gX[i],lName[i],"L");
  }
  Leg->Draw();

  ipad++;
  
  C->cd(0);
  pad[ipad]->Draw();
  pad[ipad]->cd();

  TF1 **fDen = new TF1*[Nsim]; 
  
  Float_t kp = pData->GetPlasmaK();
  Float_t sz = 2.0 / kp ;
  fDen[0] = new TF1("fDen",DensityGauss,zmin,zmax,1);
  fDen[0]->SetParameter(0,sz/PUnits::mm);
  
  Float_t gamma = 1956.95;
  Float_t kbeta = kp / TMath::Sqrt(2*gamma);
  Float_t Ltap = 10.0 / kbeta;
  Float_t lambda = (Ltap/TMath::Sqrt(Ltap * kbeta));
  
  cout << Form("  L = %.4e     lambda = %.4e      kbeta = %.4e ",Ltap * kp,lambda * kp,kbeta / kp) << endl;
  cout << Form("  L = %.4e mm  lambda = %.4e mm   kbeta = %.4e mm^-1",Ltap/PUnits::mm,lambda/PUnits::mm,kbeta * PUnits::mm) << endl;
  cout << Form("  zmin = %.2f mm  zmax = %.2f mm",zmin,zmax) << endl;
  
  fDen[1] = new TF1("fDen2",DensityTap,zmin,zmax,2);

  fDen[1]->SetParameter(0,Ltap/PUnits::mm);
  fDen[1]->SetParameter(1,lambda/PUnits::mm);

  Float_t denmintap = TMath::Power(1.0-(-Ltap/lambda),-4.0);
  cout << Form("  denmin = %.4f",denmintap) << endl;
  hFrame[ipad]->GetYaxis()->SetRangeUser(denmintap,1.4);  
  
  // hFrame[ipad]->GetXaxis()->SetTitle("n_{He} [10^{15} e/cm^{3}]");
  hFrame[ipad]->GetXaxis()->SetTitle("");
  hFrame[ipad]->GetYaxis()->SetTitle("n/n_{0}");
  
  hFrame[ipad]->Draw("AXIS");

  lineZero->Draw();

  fDen[0]->SetLineColor(color[0]);
  fDen[0]->SetLineWidth(2);
  fDen[0]->SetNpx(1000);
  // fDen[0]->Draw("same C");

  //fDen[1]->SetLineColor(color[3]);
  fDen[1]->SetLineColor(kGray+2);
  fDen[1]->SetLineWidth(2);
  fDen[1]->SetNpx(1000);
  fDen[1]->Draw("same C");
 
  C->cd();
  
  // Print to file -------------------------------------------

  TString fOutName = Form("./CentroidCompS2E/CentroidCompS2E");
  PGlobals::imgconv(C,fOutName,opt);

  // ---------------------------------------------------------
  
  
}


