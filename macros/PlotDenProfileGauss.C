#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>

#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TMath.h>
#include <TPaveText.h>
#include <TBox.h>

#include "PData.hh"
#include "PGlobals.hh"
#include "PPalette.hh"


Double_t DensityGauss(Double_t *x, Double_t *par) {
  Double_t f = 0;

  Double_t n0 = par[0];
  Double_t z0 = par[1];
  Double_t sigma0 = par[2];
  Double_t n1 = par[3];
  Double_t z1 = par[4];
  Double_t sigma1 = par[5];
  Double_t n2 = par[6];
  Double_t z2 = par[7];
  Double_t sigma2 = par[8];
  
  if( x[0] < z0) f = n0 * TMath::Gaus(x[0],z0,sigma0);
  else if( x[0] < z1) f = n0;
  else if( x[0] >= z1 && x[0]<z2) f = (n0-n1) * TMath::Gaus(x[0],z1,sigma1) + n1;
  else if( x[0] >= z2) f = (n1-n2) *  TMath::Gaus(x[0],z2,sigma2) + n2;
  
  return f;
}

Double_t DensityGaussTap(Double_t *x, Double_t *par) {
  Double_t f = 0;

  Double_t n0 = par[0];
  Double_t z0 = par[1];
  Double_t sigma0 = par[2];
  Double_t n1 = par[3];
  Double_t z1 = par[4];
  Double_t sigma1 = par[5];
  Double_t nt = par[6];
  Double_t zt = par[7];
  Double_t sigmat = par[8];
  
  if( x[0] < z0) f = n0 * TMath::Gaus(x[0],z0,sigma0);
  else if( x[0] < z1) f = n0;
  else if( x[0] >= z1 && x[0]<zt) {
    f = (n0-n1) * TMath::Gaus(x[0],z1,sigma1) + n1
      + (nt-n1) *  TMath::Gaus(x[0],zt,sigmat);
  } else if( x[0] >= zt) f = nt;
  
  return f;
}


Double_t DensityGaussSum(Double_t *x, Double_t *par) {
  Double_t f = DensityGauss(x,par) +  DensityGauss(x,&par[9]);

  
  return f;
}


void PlotDenProfileGauss(const TString &sim, const Float_t zmin0 = 0, const Float_t zmax0 = 100, const TString &options="pdf") { 
  
#ifdef __CINT__  
  gSystem->Load("libptools.so");
#endif

  PGlobals::Initialize();
  

  TString opt = options;

  // More makeup            
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  //gStyle->SetLineWidth(3);  
  
  if(opt.Contains("grid")) {
    gStyle->SetPadGridX(1);
    gStyle->SetPadGridY(1);
  }

  Int_t myBlue = TColor::GetColor((Float_t) 0.16, (Float_t) 0.83, (Float_t) 0.5);
  Int_t myNaranja = TColor::GetColor((Float_t) 0.992157, (Float_t) 0.411765, (Float_t) 0.027451);

  // Load first simulation data (for instance)
  PData *pData = PData::Get(sim.Data());
  Double_t np  = pData->GetPlasmaDensity() / (1E17/PUnits::cm3);
  Double_t skd = pData->GetPlasmaSkinDepth() / PUnits::mm;

  Double_t shiftz = 60.0;
  Double_t  zmin = zmin0 - shiftz; 
  Double_t  zmax = zmax0 - shiftz; 

  if(opt.Contains("units")) {
    zmin *= skd;
    zmax *= skd;
  }

  // arrays for graphs
  Int_t NP = 1000;
  Double_t *zarray = new Double_t[NP];
  Double_t *darray = new Double_t[NP];
  
  cout << Form("  Plasma density = %.2f E17 cm^{-3}", np) << endl;
  
  Double_t nres = 2.26;

  Double_t n0 = 10.0;
  Double_t z0 = 60.0 - shiftz;
  Double_t sigma0 = 7.5;
  Double_t n1 = 1.0;
  Double_t z1 = z0;
  Double_t sigma1 = 2.5;
  Double_t n2 = 0;
  Double_t z2 = 10000 - shiftz;
  Double_t sigma2 = 0.0;
  Double_t nt = 1.5;
  Double_t zt = 100 - shiftz;
  Double_t sigmat = 5.0;

  if(opt.Contains("units")) {
    n0 *= np;
    z0 *= skd;
    sigma0 *= skd;
    n1 *= np;
    z1 *= skd;
    sigma1 *= skd;
    n2 *= np;
    z2 *= skd;
    sigma2 *= skd;
    nt *= np;
    zt *= skd;
    sigmat *= skd;

    nres *= np;
  }
  
  
  const Int_t Nsim = 7;
  char lName[Nsim][128] = {"#tilde{n}_{top} = 10,  #tilde{#sigma}_{l} = 1.25",
			   "(a) #tilde{n}_{top} = 10,  #tilde{#sigma}_{l} = 2.5",
			   "(b) #tilde{n}_{top} = 10,  #tilde{#sigma}_{l} = 5.0",
			   "(c) #tilde{n}_{top} = 10,  #tilde{#sigma}_{l} = 7.5",
			   "(d) #tilde{n}_{top} =  5,  #tilde{#sigma}_{l} = 2.5",
			   "(e) #tilde{n}_{top} = 2.5,  #tilde{#sigma}_{l} = 2.5",
			   "(f) #tilde{n}_{top} = 10,  #tilde{#sigma}_{l} = 2.5 [tap]"};
  
  TF1 **fDenProf = new TF1*[Nsim];
  TGraph **gDenProf = new TGraph*[Nsim];

  Int_t color[Nsim];
  PPalette *pal = new PPalette("pal");
  pal->SetPalette("oli");

  for(Int_t i=0;i<Nsim;i++) {

    if(i==6) 
      fDenProf[i] = new TF1(Form("fDenProf_%i",i),DensityGaussTap,zmin,zmax,9);
    else
      fDenProf[i] = new TF1(Form("fDenProf_%i",i),DensityGauss,zmin,zmax,9);
   
    fDenProf[i]->SetNpx(10000);
    
    Int_t index =  i * pal->GetNColors() / Nsim;
    color[i] = pal->GetColorIndex(index);
    
    if(i==0)
      color[i] = color[i];
    else if(i==1)
      color[i] = kGray+3;
    else if(i==2)
      color[i] = kRed-6;
    else if(i==3)
      color[i] = kRed-9;
    else if(i==4)
      color[i] = kBlue-6;
    else if(i==5)
      color[i] = kBlue-9;
    else if(i==6)
      color[i] = kOrange; // kGray+1;
    
    
    fDenProf[i]->SetLineColor(color[i]);
    fDenProf[i]->SetLineWidth(3);

    // if(i==4) 
    //   fDenProf[i]->SetLineStyle(2);
    
  }
    
  fDenProf[0]->SetParameters(n0,z0,sigma0,n1,z1,0.5*sigma1,n2,z2,sigma2);
  
  fDenProf[1]->SetParameters(n0,z0,sigma0,n1,z1,sigma1,n2,z2,sigma2);
  
  fDenProf[2]->SetParameters(n0,z0,sigma0,n1,z1,2*sigma1,n2,z2,sigma2);
  
  fDenProf[3]->SetParameters(n0,z0,sigma0,n1,z1,3*sigma1,n2,z2,sigma2);

  fDenProf[6]->SetParameters(n0,z0,sigma0,n1,z1,sigma1,nt,zt,sigmat);
  
  // Change ramp profile
  Float_t n05 = 2.2/1.8;
  Float_t kp5 = TMath::Sqrt(1.0/n05);
  fDenProf[4]->SetParameters(0.5*n0*n05,z0*kp5,sigma0*kp5,n1*n05,z1*kp5,sigma1*kp5,n2*n05,z2*kp5,sigma2*kp5);
  
  // Change ramp profile
  Float_t n025 = 2.2/1.56;
  Float_t kp25 = TMath::Sqrt(1.0/n025);
  fDenProf[5]->SetParameters(0.25*n0*n025,z0*kp25,sigma0*kp25,n1*n025,z1*kp25,sigma1*kp25,n2*n025,z2*kp25,sigma2*kp25);

  // GRAPHS
  for(Int_t i=0;i<Nsim;i++) {
    
    for(Int_t j=0;j<NP;j++) {
      zarray[j] = zmin + j * (zmax-zmin)/(NP-1);
      darray[j] = fDenProf[i]->Eval(zarray[j]);
    }
    gDenProf[i] = new TGraph(NP,zarray,darray);
    gDenProf[i]->SetLineColor(color[i]);
    gDenProf[i]->SetLineWidth(4);
  }

    // Canvas setup
  // Create the canvas and the pads before the Frame loop
  // Resolution:
  // Int_t sizex = 800;
  // Int_t sizey = 400;
  Int_t sizex = 1024;
  Int_t sizey =  380;

  Int_t font = 43;
  Int_t labelsize = 34;
  Int_t titlesize = 38;
  
  Float_t labeloffset = 0.01;
  Float_t titleoffsety = 0.5;
  Float_t titleoffsetx = 1.0;
  
  Float_t ticksizex = 0.03;
  Float_t ticksizey = 0.01;

  // Int_t NdivX = 505;
  // Int_t NdivY = 505;
  Int_t NdivX = 6;
  Int_t NdivY = 6;
  
  //  gStyle->SetJoinLinePS(2);

  gStyle->SetTitleFont(font,"xyz");
  gStyle->SetLabelFont(font,"xyz");
  
  gStyle->SetLabelSize(labelsize,"xyz");
  gStyle->SetTitleSize(titlesize,"xyz");
    
  gStyle->SetLabelOffset(labeloffset,"xyz");

  gStyle->SetTitleOffset(titleoffsetx,"x");
  gStyle->SetTitleOffset(titleoffsety,"yz");
  
  gStyle->SetTickLength(ticksizex,"x");
  gStyle->SetTickLength(ticksizey,"yz");

  Float_t lMargin = 0.15;
  Float_t rMargin = 0.10;
  Float_t bMargin = 0.25;
  Float_t tMargin = 0.05;

  gStyle->SetPadLeftMargin(lMargin);
  gStyle->SetPadRightMargin(rMargin);
  gStyle->SetPadBottomMargin(bMargin);
  gStyle->SetPadTopMargin(tMargin);

  Int_t frameWidth = 3;
  gStyle->SetLineWidth(frameWidth);

  if(opt.Contains("gridx")) {
    gStyle->SetPadGridX(1);
  }
  if(opt.Contains("gridy")) {
    gStyle->SetPadGridY(1);
  }
  
  char cName[32];
  sprintf(cName,"C");     
  TCanvas *C = (TCanvas*) gROOT->FindObject(cName);
  if(C==NULL) C = new TCanvas("C","",sizex,sizey);
  C->SetFillStyle(4000);
  C->cd();
  C->Clear();

  char name[16];
  sprintf(name,"hFrame");
  TH1F *hFrame = new TH1F(name,"",10,zmin,zmax);

  hFrame->GetXaxis()->CenterTitle();
  hFrame->GetXaxis()->SetNdivisions(NdivX);

  hFrame->GetYaxis()->CenterTitle();
  hFrame->GetYaxis()->SetNdivisions(NdivY);

  Float_t maxDen = fDenProf[0]->GetMaximum();
  hFrame->GetYaxis()->SetRangeUser(0.0, 1.2 * maxDen);  
  
  // hFrame->GetXaxis()->SetTitle("n_{He} [10^{15} e/cm^{3}]");
  hFrame->GetXaxis()->SetTitle("k_{p}^{0} z");
  hFrame->GetYaxis()->SetTitle("n/n_{0}");

  if(opt.Contains("units")) {
    hFrame->GetXaxis()->SetTitle("z [mm]");
    hFrame->GetYaxis()->SetTitle("n [10^{17} cm^{-3}]");
  }
  
  hFrame->Draw("AXIS");

  // Lines and guides
  Float_t xMin = hFrame->GetXaxis()->GetXmin();
  Float_t xMax = hFrame->GetXaxis()->GetXmax();
  Float_t yMin = hFrame->GetYaxis()->GetXmin();
  Float_t yMax = 1.2 * maxDen;

  TLine *linezero = new TLine(0.0,yMin,0.0,yMax);
  linezero->SetLineColor(kGray);
  linezero->SetLineWidth(2);
  linezero->SetLineStyle(3);
  linezero->Draw();

  TLine *linemax = new TLine(xMin,maxDen,xMax,maxDen);
  linemax->SetLineColor(kGray);
  linemax->SetLineWidth(2);
  linemax->SetLineStyle(3);
  linemax->Draw();

  TLine *lineplateau = new TLine(xMin,n1,xMax,n1);
  lineplateau->SetLineColor(kGray);
  lineplateau->SetLineWidth(2);
  lineplateau->SetLineStyle(3);
  lineplateau->Draw();

  TLine *linenres = new TLine(xMin,nres,xMax,nres);
  linenres->SetLineColor(kGray+2);
  linenres->SetLineWidth(2);
  linenres->SetLineStyle(2);
  linenres->Draw();

  // Functions


  for(Int_t i=Nsim-1;i>=0;i--) {
    //for(Int_t i=0;i<Nsim;i++) {
    if(i==0) continue;
    if(i==3) 
      gDenProf[6]->Draw("C same");
    gDenProf[i]->Draw("C same");
  }
  
  hFrame->Draw("AXIS same");

  gPad->Update();

  TBox *lFrame = new TBox(gPad->GetUxmin(), gPad->GetUymin(),
  			  gPad->GetUxmax(), gPad->GetUymax());
  lFrame->SetFillStyle(0);
  lFrame->SetLineColor(kBlack);
  lFrame->SetLineWidth(frameWidth);
  lFrame->Draw();

  TLegend *Leg = new TLegend(0.55,0.40,0.85,0.83);
  PGlobals::SetPaveStyle(Leg);
  Leg->SetTextAlign(12);
  Leg->SetTextColor(kGray+3);
  Leg->SetTextFont(53);
  Leg->SetTextSize(20);
  Leg->SetLineColor(kGray+1);
  Leg->SetBorderSize(0);
  Leg->SetLineWidth(1);
  Leg->SetFillColor(0);
  Leg->SetFillStyle(1001);
  Leg->SetFillStyle(0); // Hollow

  for(Int_t i=0;i<Nsim;i++) {
    if(i==0) continue;
    Leg->AddEntry(fDenProf[i],lName[i],"L");
  }
  
  Leg->Draw();
  
  
  // Print to file --------------------------------------
  
  C->cd();
  
  // Print to a file
  // Output file
  TString fOutName = Form("./DenProf/DenProf-%s",sim.Data());
  PGlobals::imgconv(C,fOutName,opt);
  // ---------------------------------------------------------

  
}
  
