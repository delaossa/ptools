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
#include "PlasmaGlob.hh"
#include "PPalette.hh"

void PlotRakeBunch( const TString &sim, Int_t time, Int_t index = 0, const TString &options="") {
  
#ifdef __CINT__  
  gSystem->Load("libplasma.so");
#endif

  PlasmaGlob::Initialize();

  TString opt = options;
 
  // Palettes!
  gROOT->Macro("PlasmaPalettes.C");

  if(opt.Contains("grid")) {
    gStyle->SetPadGridX(1);
    gStyle->SetPadGridY(1);
  }
  gStyle->SetLabelFont(42,"xyz");
  gStyle->SetTextFont(62);
 

  // Load PData
  PData *pData = PData::Get(sim.Data());
  pData->LoadFileNames(time);
  if(!pData->IsInit()) return;


  Int_t Nspecies = pData->NSpecies();
  if(index>Nspecies-1) {
    return;
  }
  if(!pData->GetRawFileName(index)) {
    return;    
  }


  Bool_t CYL = kFALSE;
  if(sim.Contains("cyl")) CYL = kTRUE; 
    
  Bool_t ThreeD = kFALSE;
  if(sim.Contains("3D")) ThreeD = kTRUE; 

  // Some plasma constants
  Double_t n0 = pData->GetPlasmaDensity();
  Double_t kp = pData->GetPlasmaK();
  Double_t skindepth = 1.;
  if(kp!=0.0) skindepth = 1/kp;
  Double_t E0 = pData->GetPlasmaE0();

  // Time in OU
  Float_t Time = pData->GetRealTime();
  // z start of the plasma in normalized units.
  Float_t zStartPlasma = pData->GetPlasmaStart()*kp;
  // z start of the beam in normalized units.
  Float_t zStartBeam = pData->GetBeamStart()*kp;
  
  opt += "comovcenter";

  // Centering time and z position:
  Double_t shiftz = pData->Shift(opt);
  TString sshiftz = Form("(x1-%f)",shiftz);

  if(opt.Contains("center")) {
    Time -= zStartPlasma;
    if(opt.Contains("comov"))      // Centers on the head of the beam.
      Time += zStartBeam;
  } 

  // Spatial resolution
  Float_t dx1 = pData->GetDX(0);
  Float_t dx2 = pData->GetDX(1);
  Float_t dx3 = pData->GetDX(2);
  
  // Spatial coordinates intervals:
  Float_t x1Min = -4.3;
  Float_t x1Max = -3.9;
  Float_t x2Min = -0.5;
  Float_t x2Max =  0.5;
  Float_t x3Min = -0.5;
  Float_t x3Max =  0.5;

  // Momentum coordinates intervals:
  Float_t p1Min =  6500.01;
  Float_t p1Max =  10099.99;
  Float_t p2Min = -2.0;
  Float_t p2Max =  2.0;
  Float_t p3Min = -2.0;
  Float_t p3Max =  2.0;

  // Bining, intervals, labels, etc.
  Int_t xNbin = 200;
  Int_t pNbin = 200;
  
  if(sim.Contains(".DDR.") || sim.Contains(".DR.")) {
    // xNbin = 150;
    // pNbin = 150;
   
    // t = 1195
    p1Min =  1700;
    p1Max =  4200;
    p2Min = -5;
    p2Max =  5;
    p3Min = -5;
    p3Max =  5;

    x1Min = -5.0;
    x1Max = -3.5;
    x2Min = -0.1;
    x2Max =  0.1;
    x3Min = -0.1;
    x3Max =  0.1;

  } else if(sim.Contains("flash")) {
    if(sim.Contains(".G.") ) {    
      x1Min = -6.3;
      x1Max = -5.0;
      
      //p1Min =  0.001;
      p1Min =  650.001;
      p1Max =  1499.99;
    } else if(sim.Contains("v10.RI") ) {
      x1Min = -4.3;
      x1Max = -3.9;
      
      p2Min = -15;
      p2Max =  15;
      p3Min = -15;
      p3Max =  15;
      // t=1510
      p1Min =  6000.;
      p1Max =  10500.;
    }
  } else if(sim.Contains("facet_v23kA.G.RI.e")) {
    x1Min = -7.8;
    x1Max = -7.2;

    x2Min = -0.5;
    x2Max =  0.5;
    x3Min = -0.5;
    x3Max =  0.5;
    p2Min =  -15;
    p2Max =  15;
    p3Min =  -15;
    p3Max =  15;
    // t = 38
    p1Min =  400.01;
    p1Max =  899.99;
  } else if(sim.Contains("facet_v23kA.G.RI")) {
    xNbin = 200;
    pNbin = 200;

    x1Min = -7.8;
    x1Max = -7.0;

    x2Min = -0.5;
    x2Max =  0.5;
    x3Min = -0.5;
    x3Max =  0.5;

    p2Min = -20;
    p2Max =  20;
    p3Min = -20;
    p3Max =  20;

    // t=40
    // p1Min =  140.01;
    // p1Max =  459.99;
    // t=60
    // p1Min =  350.01;
    // p1Max =  799.99;
    // t = 100
    // p1Min =  700.01;
    // p1Max =  1399.99;
    // t=105
    // p1Min =  800.001;
    // p1Max =  1399.99;
    // t = 120
    // p1Min =  900.01;
    // p1Max =  1699.99;
    // t = 140
    // p1Min =  1100.01;
    // p1Max =  2199.99;
    // t=180
    // p1Min =  1600.0;
    // p1Max =  2700.0;
    // t=220
    // p1Min =  2100.01;
    // p1Max =  3100.99; 
    // t=235
    // p1Min =  2000.01;
    // p1Max =  3499.99; 
    // t=280
    // p1Min =  2400.01;
    // p1Max =  4299.99;
    // t=310
    // p1Min =  2800.01;
    // p1Max =  4799.99; 
    // t=320
    // p1Min =  2800.01;
    // p1Max =  4799.99; 
    // t=390
    // p1Min =  3400.01;
    // p1Max =  5899.99; 
    // t=420
    p1Min =  4050.01;
    p1Max =  6199.99; 
    // t=520
    // p1Min =  5000.01;
    // p1Max =  8199.99; 
    // t=620
    // p1Min =  6000.01;
    // p1Max =  9199.99; 
    // t=720
    // p1Min =  7000.01;
    // p1Max =  10999.99; 
  } else if(sim.Contains("facet_v23kA.G")) {
    x1Min = -7.8;
    x1Max = -7.1;

    x2Min = -0.5;
    x2Max =  0.5;
    x3Min = -0.5;
    x3Max =  0.5;
    // t = 120
    // p1Min =  800.01;
    // p1Max =  1599.99;
    // t=150
    // p1Min =  1100.001;
    // p1Max =  1999.99;
    // t=183
    p1Min =  1500.;
    p1Max =  2700.;
  } 
  
  cout << endl;
  cout << Form("Ranges:  x1Min = %5.2f, x1Min = %5.2f, p1Min = %5.2f, p1Max = %5.2f",x1Min,x1Max,p1Min,p1Max) << endl; 
  cout << Form("         x2Min = %5.2f, x2Min = %5.2f, p2Min = %5.2f, p2Max = %5.2f",x2Min,x2Max,p2Min,p2Max) << endl; 

  char hName[24];
  char dCommand[128];
 
  TH1F *hX1 = NULL;
  TH1F *hP1 = NULL;

  cout << Form("\n1. Getting TTree... ") << endl; 

  char cutString[512];
  sprintf(cutString,"TMath::Abs(q)*(x1 > %.1f && x1 < %.1f && x2 > %.1f && x2 < %.1f && x3 > %.1f && x3 < %.1f && p1 > %.1f && p1 < %.1f && p2 > %.1f && p2 < %.1f && p3 > %.1f && p3 < %.1f)",x1Min,x1Max,x2Min,x2Max,x3Min,x3Max,p1Min,p1Max,p2Min,p2Max,p3Min,p3Max); 
  //sprintf(cutString,"TMath::Abs(q)"); 
  
  TCut Cut = cutString;
  //  cout << Form("   (applied cut: \n %s)",cutString) << endl;
  
  TTree *tree = pData->GetTreeRaw(pData->GetRawFileName(index)->c_str(),opt);

 
  // Auto ranges 
  if (opt.Contains("autoz")) {
    sprintf(hName,"hX1");
    hX1 = (TH1F*) gROOT->FindObject(hName);
    if(hX1) delete hX1;
    sprintf(dCommand,"x1>>%s",hName);
    tree->Draw(dCommand,Cut,"goff");
    hX1 = (TH1F*) gROOT->FindObject(hName);
    cout << Form("   - x1. ") << endl;    
    x1Min = hX1->GetXaxis()->GetXmin();
    x1Max = hX1->GetXaxis()->GetXmax();


    sprintf(hName,"hP1");
    hP1 = (TH1F*) gROOT->FindObject(hName);
    if(hP1) delete hP1;
    hP1 = new TH1F(hName,"",pNbin,p1Min,p1Max);
    
  } else if (opt.Contains("autop")) {
    sprintf(hName,"hP1");
    hP1 = (TH1F*) gROOT->FindObject(hName);
    if(hP1) delete hP1;
    sprintf(dCommand,"p1>>%s",hName);
    tree->Draw(dCommand,Cut,"goff");
    hP1 = (TH1F*) gROOT->FindObject(hName);
    cout << Form("   - p1. ") << endl;        
    p1Min = hP1->GetXaxis()->GetXmin();
    p1Max = hP1->GetXaxis()->GetXmax();

    sprintf(hName,"hX1");
    hX1 = (TH1F*) gROOT->FindObject(hName);
    if(hX1) delete hX1;
    hX1 = new TH1F(hName,"",xNbin,x1Min,x1Max);
  } else if(opt.Contains("auto")) {

    sprintf(hName,"hX1");
    sprintf(dCommand,"x1>>%s",hName);
    tree->Draw(dCommand,Cut,"goff");
    hX1 = (TH1F*) gROOT->FindObject(hName);
    cout << Form("   - x1. ") << endl;    
    x1Min = hX1->GetXaxis()->GetXmin();
    x1Max = hX1->GetXaxis()->GetXmax();

    sprintf(hName,"hP1");
    sprintf(dCommand,"p1>>%s",hName);
    tree->Draw(dCommand,Cut,"goff");
    hP1 = (TH1F*) gROOT->FindObject(hName);
    cout << Form("   - p1. ") << endl;        
    p1Min = hP1->GetXaxis()->GetXmin();
    p1Max = hP1->GetXaxis()->GetXmax();    
  } else {
    sprintf(hName,"hX1");
    hX1 = (TH1F*) gROOT->FindObject(hName);
    if(hX1) delete hX1;
    hX1 = new TH1F(hName,"",xNbin,x1Min,x1Max);

    sprintf(hName,"hP1");
    hP1 = (TH1F*) gROOT->FindObject(hName);
    if(hP1) delete hP1;
    hP1 = new TH1F(hName,"",pNbin,p1Min,p1Max);
  }

  sprintf(hName,"hP1X1");
  TH2F *hP1X1 = (TH2F*) gROOT->FindObject(hName);
  if(hP1X1) delete hP1X1;
  hP1X1 = new TH2F(hName,"",xNbin,x1Min,x1Max,pNbin,p1Min,p1Max);

  sprintf(hName,"hP2X2");
  TH2F *hP2X2 =  (TH2F*) gROOT->FindObject(hName);
  if(hP2X2) delete hP2X2;
  hP2X2 = new TH2F(hName,"",xNbin,x2Min,x2Max,pNbin,p2Min,p2Max);

  // Sliced quantities:
  // --------------------------------------------------------------------------

  // cout << Form("\n5. Slicing ") << endl;
  // Binning
  
  cout << Form("\n2. Setting the binning: ") << endl; 

  // Analize hX1 for the binning:
  
  Int_t SNbin = 10;
  Float_t x1BinMin = -4.16;
  Float_t x1BinMax = -4.02;
  if(sim.Contains("DDR.") || sim.Contains(".DR.")) {
    SNbin = 10;
    x1BinMin = -4.7;
    x1BinMax = -3.8;
  } else if(sim.Contains("flash")) {
    SNbin = 20;
    x1BinMin = -4.15;
    x1BinMax = -4.02;
  } else if(sim.Contains("facet_v23kA.G.RI.e")) {
    SNbin = 30;
    x1BinMin = -7.60;
    x1BinMax = -7.35;
  } else if(sim.Contains("facet_v23kA.G.RI")) {
    SNbin = 40;
    x1BinMin = -7.54;
    x1BinMax = -7.14;
  } else if(sim.Contains("facet_v23kA.G")) {
    SNbin = 40;
    x1BinMin = -7.55;
    x1BinMax = -7.15;
  } 
  
  // Set the binning
  Float_t *sBinLim = new Float_t[SNbin+1];
  sBinLim[0] = x1BinMin;
  sBinLim[SNbin] = x1BinMax;
  Float_t slbinSize = (sBinLim[SNbin] - sBinLim[0])/SNbin;
  for(Int_t i=1;i<SNbin;i++) {
    sBinLim[i] = sBinLim[i-1] + slbinSize;
  }
  
  TH1F **hP1sl = new TH1F*[SNbin];
  TH2F **hP2X2sl = new TH2F*[SNbin];
  for(Int_t k=0;k<SNbin;k++) {
    sprintf(hName,"hP2X2sl_%2i",k);
    hP2X2sl[k] = (TH2F*) gROOT->FindObject(hName);
    if(hP2X2sl[k]) delete hP2X2sl[k];
    hP2X2sl[k] = new TH2F(hName,"",xNbin,x2Min,x2Max,pNbin,p2Min,p2Max);

    hP2X2sl[k]->GetXaxis()->SetTitle("k_{p} y");
    hP2X2sl[k]->GetYaxis()->SetTitle("p_{y}/mc");


    sprintf(hName,"hP1sl_%2i",k);
    hP1sl[k] = (TH1F*) gROOT->FindObject(hName);
    if(hP1sl[k]) delete hP1sl[k];
    hP1sl[k] = new TH1F(hName,"",pNbin,p1Min,p1Max);

  }

  // Prepare the branches for reading from the tree:
  const Int_t Nvar = 8;
  Float_t var[Nvar];  
  char varname[Nvar][4] = {{"ene"},{"p1"},{"p2"},{"p3"},{"q"},{"x1"},{"x2"},{"x3"}};
  for(Int_t i=0;i<Nvar;i++) {
    tree->SetBranchAddress(varname[i],&var[i]);
  }

  // Main loop!
  cout << Form("\n3. Filling histograms: ") << endl; 

  hX1->Reset();
  hX1->SetBins(xNbin,x1Min,x1Max);
  hP1->Reset();
  hP1->SetBins(pNbin,p1Min,p1Max);

  Int_t nentries = (Int_t)tree->GetEntries();  
  for(Int_t i=0;i<nentries;i++) {
    tree->GetEntry(i);
    var[5] = var[5] - shiftz;
    
    if(var[5]<x1Min || var[5]>x1Max ) continue; 
    if(var[6]<x2Min || var[6]>x2Max ) continue; 
    if(var[7]<x3Min || var[7]>x3Max ) continue; 
    if(var[1]<p1Min || var[1]>p1Max ) continue; 
    if(var[2]<p2Min || var[2]>p2Max ) continue; 
    if(var[3]<p3Min || var[3]>p3Max ) continue; 
    
    hX1->Fill(var[5],TMath::Abs(var[4]));
    hP1->Fill(var[1],TMath::Abs(var[4]));
 
    hP1X1->Fill(var[5],var[1],TMath::Abs(var[4]));
        
    // Slices    
    if(var[5]<sBinLim[0] || var[5]>sBinLim[SNbin]) continue;
    Int_t iBin = -1;
    for(Int_t j=0; j<SNbin; j++) {
      
      if(var[5]<sBinLim[j+1]) {
	iBin = j;
	break;
      }
    }
    if(iBin<0) continue;

    // projected emittance in the bunch range. (skip the tails to "match" the sliced ones)
    hP2X2->Fill(var[6],var[2],TMath::Abs(var[4]));


    //    cout << Form(" entry:  %i , bin:  %i",i,iBin)<< endl;

    hP1sl[iBin]->Fill(var[1],TMath::Abs(var[4]));   
    hP2X2sl[iBin]->Fill(var[6],var[2],TMath::Abs(var[4]));
  
  
  }
  

  // Integrated long. emittance:
  cout << Form("\n4. Calculating integrated quantities.. ") << endl;

  // Longitudinal phasespace
  Double_t xmean = 0.0;
  Double_t ymean = 0.0;
  Double_t x2mean = 0.0;
  Double_t y2mean = 0.0;
  Double_t xymean = 0.0;
  Double_t Ntotal = 0.0;
  for(Int_t i=1;i<=xNbin;i++) {
    Double_t x = hP1X1->GetXaxis()->GetBinCenter(i);
    // if(x<xmin || x>xmax) continue;
    for(Int_t j=1;j<=pNbin;j++) {
      Double_t y = hP1X1->GetYaxis()->GetBinCenter(j);
      // if(y<ymin || y>ymax) continue;
      Double_t value = TMath::Abs(hP1X1->GetBinContent(i,j));
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

  cout << Form("  xMean = %7.3f   yMean = %7.3f",xmean,ymean) << endl;
  cout << Form("  xRms  = %7.3f   yRms  = %7.3f",xrms,yrms) << endl;
  cout << Form("  Emittance = %7.3f",emittance) << endl;

  Double_t emitz = emittance;
  Double_t zmean = xmean;
  Double_t zrms = xrms;
  Double_t pzmean = ymean;
  Double_t pzrms = yrms;

  // Transverse phasespace
  xmean = 0.0;
  ymean = 0.0;
  x2mean = 0.0;
  y2mean = 0.0;
  xymean = 0.0;
  Ntotal = 0.0;
  for(Int_t i=1;i<=xNbin;i++) {
    Double_t x = hP2X2->GetXaxis()->GetBinCenter(i);
    // if(x<xmin || x>xmax) continue;
    for(Int_t j=1;j<=pNbin;j++) {
      Double_t y = hP2X2->GetYaxis()->GetBinCenter(j);
      // if(y<ymin || y>ymax) continue;
      Double_t value = TMath::Abs(hP2X2->GetBinContent(i,j));
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

  xrms2  = x2mean - xmean*xmean;
  yrms2  = y2mean - ymean*ymean;
  xrms   = TMath::Sqrt(xrms2);
  yrms   = TMath::Sqrt(yrms2);
  xyrms2 = xymean - xmean*ymean;

  emittance = TMath::Sqrt(xrms2*yrms2 - xyrms2*xyrms2);

  cout << Form("  xMean = %7.3f   yMean = %7.3f",xmean,ymean) << endl;
  cout << Form("  xRms  = %7.3f   yRms  = %7.3f",xrms,yrms) << endl;
  cout << Form("  Emittance = %7.3f",emittance) << endl;

  Double_t emity = emittance;
  Double_t y_mean = xmean;
  Double_t y_rms = xrms;
  Double_t pymean = ymean;
  Double_t pyrms = yrms;

  
  // Charge  
  hX1->Scale(dx1*dx2*dx3);
  Double_t Charge = hX1->Integral();
  
  // Charge *= dx1*dx2*dx3;
 
  if(opt.Contains("units")) {
    Double_t dV = skindepth * skindepth * skindepth;
    Charge *= n0 * dV * (PConst::ElectronCharge/PUnits::picocoulomb);
    cout << Form(" Integrated charge (RAW) of specie %3i = %8f pC",index,Charge) << endl;
  } else {
    cout << Form(" Integrated charge (RAW) of specie %3i = %8.4f n0 * kp^-3",index,Charge) << endl;
  }
  


  cout << Form("\n6. Calculating sliced quantities.. ") << endl;

  TGraph *gemit = NULL;
  TGraph *gYrms = NULL;
  TGraph *gErms = NULL;
  TGraph *gErmsB = NULL;
 
  Double_t * sxmean = new Double_t[SNbin];
  Double_t * symean = new Double_t[SNbin];
  Double_t * sx2mean = new Double_t[SNbin];
  Double_t * sy2mean = new Double_t[SNbin];
  Double_t * sxymean = new Double_t[SNbin];
  Double_t * sNtotal = new Double_t[SNbin];
  Double_t * sxrms2 = new Double_t[SNbin];  
  Double_t * syrms2 = new Double_t[SNbin]; 
  Double_t * sxrms = new Double_t[SNbin];  
  Double_t * syrms = new Double_t[SNbin];  
  Double_t * sxyrms2 = new Double_t[SNbin];
 
  Double_t * xbin = new Double_t[SNbin];
  Double_t * semittance = new Double_t[SNbin];

  Double_t * sNEtotal = new Double_t[SNbin]; 
  Double_t * sEmean = new Double_t[SNbin];
  Double_t * sE2mean = new Double_t[SNbin];
  Double_t * sErms = new Double_t[SNbin];

  for(Int_t k=0;k<SNbin;k++) {
    sxmean[k] = symean[k] = sx2mean[k] = sy2mean[k] = sxymean[k] 
      = sNtotal[k] = sxrms2[k] = syrms2[k] = sxrms[k] = syrms[k]
      = sxyrms2[k] = xbin[k] = semittance[k] = 0.0;
    sNEtotal[k] = sEmean[k] = sE2mean[k] = sErms[k] = 0.0;
    
    xbin[k] = (sBinLim[k] + sBinLim[k+1])/2.;
    
    for(Int_t i=1;i<=xNbin;i++) {
      Double_t x = hP2X2sl[k]->GetXaxis()->GetBinCenter(i);
      // if(x<xmin || x>xmax) continue;
      for(Int_t j=1;j<=pNbin;j++) {
	Double_t y = hP2X2sl[k]->GetYaxis()->GetBinCenter(j);
	// if(y<ymin || y>ymax) continue;
	Double_t value = TMath::Abs(hP2X2sl[k]->GetBinContent(i,j));
	sxmean[k] += x*value;
	symean[k] += y*value;
	sx2mean[k] += x*x*value;
	sy2mean[k] += y*y*value;
	sxymean[k] += x*y*value;
	
	sNtotal[k] += value;
      }	
    }
    
    for(Int_t i=1;i<=pNbin;i++) {
      Double_t y = hP1sl[k]->GetXaxis()->GetBinCenter(i);
      Double_t value = TMath::Abs(hP1sl[k]->GetBinContent(i));
      sEmean[k] += y*value;
      sE2mean[k] += y*y*value;
      sNEtotal[k] += value;
    }
    
    sxmean[k]  /= sNtotal[k];
    symean[k]  /= sNtotal[k];
    sx2mean[k] /= sNtotal[k];
    sy2mean[k] /= sNtotal[k];
    sxymean[k] /= sNtotal[k];
      
    sxrms2[k]  = sx2mean[k] - sxmean[k]*sxmean[k];
    syrms2[k]  = sy2mean[k] - symean[k]*symean[k];
    sxrms[k]   = TMath::Sqrt(sxrms2[k]);
    syrms[k]   = TMath::Sqrt(syrms2[k]);
    sxyrms2[k] = sxymean[k] - sxmean[k]*symean[k];
      
    semittance[k] = TMath::Sqrt(sxrms2[k]*syrms2[k] - sxyrms2[k]*sxyrms2[k]);

    sEmean[k]  /= sNEtotal[k];
    sE2mean[k] /= sNEtotal[k];
    sErms[k]   =  TMath::Sqrt(sE2mean[k] - sEmean[k]*sEmean[k]);
    
    
    cout<< Form("\nk = %i : (x1 > %f && x1 < %f)",k,sBinLim[k],sBinLim[k+1]) << endl; 

    cout << Form("  xMean = %7.3f   yMean = %7.3f",sxmean[k],symean[k]) << endl;
    cout << Form("  xRms  = %7.3f   yRms  = %7.3f",sxrms[k],syrms[k]) << endl;
    cout << Form("  Emittance = %7.3f",semittance[k]) << endl;

    cout << Form("  Emean = %7.3f   Erms = %7.3f",sEmean[k],sErms[k]) << endl;
    

  }



  // Changing to user units: 
  // --------------------------
  
  if(opt.Contains("units") && n0) {
    

    x1Min *= skindepth / PUnits::um;
    x1Max *= skindepth / PUnits::um;
    p1Min *= PConst::ElectronMassE / PUnits::GeV;
    p1Max *= PConst::ElectronMassE / PUnits::GeV;

    hP1X1->SetBins(xNbin,x1Min,x1Max,pNbin,p1Min,p1Max);
    
    // Converting electron density
    Double_t dVb = skindepth * skindepth * skindepth;
    Double_t dX = (x1Max-x1Min)/xNbin; 
    Double_t dE = (p1Max-p1Min)/pNbin; 
    for(Int_t j=0;j<hP1X1->GetNbinsX();j++) {
      for(Int_t k=0;k<hP1X1->GetNbinsY();k++) {
	Double_t binValue =  fabs(hP1X1->GetBinContent(j,k) * dx1 * dx2 * dx3 * dVb * n0 *
				  (PConst::ElectronCharge/PUnits::picocoulomb));
     	//cout << Form(" value = %f",binValue) << endl;
	hP1X1->SetBinContent(j,k,binValue);
	
      }
    }
    
    if(opt.Contains("comov"))
      hP1X1->GetXaxis()->SetTitle("#zeta [#mum]");
    else
      hP1X1->GetXaxis()->SetTitle("z [#mum]");
    
    hP1X1->GetYaxis()->SetTitle("p_{z} [GeV/c]");
    
    hP1X1->GetZaxis()->SetTitle("|Charge| [pC]");
    hP1X1->GetZaxis()->CenterTitle();

    hP1->SetBins(pNbin,p1Min,p1Max);
    hP1->GetYaxis()->SetTitle("p_{z} [GeV/c]");

    hX1->SetBins(xNbin,x1Min,x1Max);
    Double_t binSize = (x1Max - x1Min)/xNbin;
    
    Double_t dV = skindepth * skindepth * skindepth;
    Double_t  lightspeed =  PConst::c_light / (PUnits::um/PUnits::femtosecond);
    hX1->Scale(TMath::Abs(n0 * dV * (PConst::ElectronCharge/PUnits::picocoulomb) * (lightspeed/binSize)));
    
    // hX1->Scale(TMath::Abs((PUnits::um/skindepth)*(PConst::ElectronCharge/PUnits::picocoulomb)*PConst::c_light));
    
    // hX1->GetYaxis()->SetTitle("I[kA]");
    hX1->GetYaxis()->SetTitle("");
    if(opt.Contains("comov"))
      hX1->GetXaxis()->SetTitle("#zeta [#mum]");
    else
      hX1->GetXaxis()->SetTitle("z [#mum]");
    
    hP1->SetBins(pNbin,p1Min,p1Max);
    hP1->GetYaxis()->SetTitle("p_{z} [GeV/c]");   
    
  
    zmean *= skindepth / PUnits::um;
    zrms  *= skindepth / PUnits::um;
  
    pzmean *= PConst::ElectronMassE / PUnits::GeV;
    pzrms  *= PConst::ElectronMassE / PUnits::GeV;
    
    emitz *= (skindepth / PUnits::um);

    // Transverse phase-space
    x2Min *= skindepth/PUnits::um;
    x2Max *= skindepth/PUnits::um;
    p2Min *= PConst::ElectronMassE / PUnits::MeV;
    p2Max *= PConst::ElectronMassE / PUnits::MeV;
    hP2X2->SetBins(xNbin,x2Min,x2Max,pNbin,p2Min,p2Max);
    for(Int_t j=0;j<hP2X2->GetNbinsX();j++) {
      for(Int_t k=0;k<hP2X2->GetNbinsY();k++) {
	Double_t binValue =  fabs(hP2X2->GetBinContent(j,k) * dx1 * dx2 * dx3 * dVb * n0 *
				  (PConst::ElectronCharge/PUnits::picocoulomb));
     	//cout << Form(" value = %f",binValue) << endl;
	hP2X2->SetBinContent(j,k,binValue);
	
      }
    }
    
    hP2X2->GetXaxis()->SetTitle("y [#mum]");
    hP2X2->GetYaxis()->SetTitle("p_{y} [MeV/c]");
    hP2X2->GetZaxis()->SetTitle("dQ/dydp_{y} [pC]");
    hP2X2->GetZaxis()->CenterTitle();

    emity *= (skindepth / PUnits::um);
 
    for(Int_t k=0;k<SNbin;k++) {
      
      hP2X2sl[k]->SetBins(xNbin,x2Min,x2Max,pNbin,p2Min,p2Max);
      // for(Int_t j=0;j<hP2X2sl[k]->GetNbinsX();j++) {
      // 	for(Int_t l=0;l<hP2X2sl[k]->GetNbinsY();l++) {
      // 	  Double_t binValue =  fabs(hP2X2sl[k]->GetBinContent(j,l) * dx1 * dx2 * dx3 * dVb * n0 *
      // 				    (PConst::ElectronCharge/PUnits::picocoulomb));
      // 	  //cout << Form(" value = %f",binValue) << endl;
      // 	  hP2X2sl[k]->SetBinContent(j,l,binValue);
	  
      // 	}
      // }
	
      hP2X2sl[k]->GetZaxis()->SetTitle("dQ/dydp_{y} [a.u.]");
      hP2X2sl[k]->GetXaxis()->SetTitle("y [#mum]");
      hP2X2sl[k]->GetYaxis()->SetTitle("p_{y} [MeV/c]");

      hP2X2sl[k]->GetZaxis()->CenterTitle();
      
      xbin[k] *= skindepth / PUnits::um;
      
      sxmean[k] *= skindepth / PUnits::um;
      sxrms[k]  *= skindepth / PUnits::um;
      symean[k] *= PConst::ElectronMassE / PUnits::MeV;
      syrms[k] *= PConst::ElectronMassE / PUnits::MeV;
      
      semittance[k] *= (skindepth / PUnits::um);
      
      sEmean[k] *= PConst::ElectronMassE / PUnits::GeV;
      sErms[k]  *= 100 * PConst::ElectronMassE / PUnits::GeV / pzmean; //sEmean[k];
      // sErms[k]  *= PConst::ElectronMassE / PUnits::GeV;
      
    }
  }
    
  // Centering in the x1 mean
  if(opt.Contains("zmean")) {
    hX1->SetBins(xNbin,x1Min-zmean,x1Max-zmean);
    hP1X1->SetBins(xNbin,x1Min-zmean,x1Max-zmean,pNbin,p1Min,p1Max);
    for(Int_t k=0;k<SNbin;k++) {
      xbin[k] -= zmean;
    }
    zmean = 0.0;
  }
  // ------

  // Create the graph with the sliced quantities:
  gemit = new TGraph(SNbin,xbin,semittance);
  gYrms = new TGraph(SNbin,xbin,sxrms);
  gErms = new TGraph(SNbin,xbin,sErms);
  
  
  // Profile energy for p1 vs x1:
  TString pname = hP1X1->GetName();
  pname += "_pfx";
  TProfile *hP1X1prof = (TProfile*) gROOT->FindObject(pname.Data());
  if(hP1X1prof) { delete hP1X1prof; hP1X1prof = NULL; }
  hP1X1prof = hP1X1->ProfileX("_pfx",1,-1,"s");

  // get the errors from the profile:
  Int_t NP1X1Bins = hP1X1prof->GetNbinsX();
  Double_t *x1bins = new Double_t[NP1X1Bins];
  Double_t *eRms   = new Double_t[NP1X1Bins];
  for(Int_t i=1;i<=hP1X1prof->GetNbinsX();i++) {
    x1bins[i] = hP1X1prof->GetBinCenter(i);
    eRms[i] = 100 * hP1X1prof->GetBinError(i) / hP1X1prof->GetBinContent(i);
  }
  gErmsB = new TGraph(NP1X1Bins,x1bins,eRms);
  
  // Vertical Energy histogram:
  // --------------------------------------------------------------------------------   
  TGraph *gP1left = NULL;
  if(hP1) {
    Double_t *yarray   = new Double_t[pNbin];
    Double_t *xarray   = new Double_t[pNbin];
    
    // This is for the right side:
    // Double_t xMax = x1Min + (x1Max-x1Min) * 0.9;
    // Double_t xMin = x1Max;
    // And this for left:
    Double_t xMin = hX1->GetXaxis()->GetXmin();
    Double_t xMax = hX1->GetXaxis()->GetXmin() + (hX1->GetXaxis()->GetXmax()
						  -hX1->GetXaxis()->GetXmin()) * 0.2;
    Double_t EneMax = hP1->GetMaximum();
    // cout << Form("  EneMax = %f ", EneMax) << endl;
 
    for(Int_t j=0; j<pNbin; j++) {
      yarray[j] = hP1->GetBinCenter(j+1);
      xarray[j] = ((xMax-xMin)/EneMax)*hP1->GetBinContent(j+1) + xMin;

      // cout << Form("  x = %f  y = %f ", xarray[j],yarray[j]) << endl;
    }

    gP1left = new TGraph(pNbin,xarray,yarray);
    gP1left->SetLineColor(PlasmaGlob::elecLine);
    gP1left->SetLineWidth(2);
    gP1left->SetFillStyle(1001);
    gP1left->SetFillColor(PlasmaGlob::elecFill);
       
  }

  
  // Ranges!!
  Double_t yMin =  999.9;
  Double_t yMax =  -999.9;
  for(Int_t k=0;k<SNbin;k++) {
    if(semittance[k]<yMin)
      yMin = semittance[k];
    
    if(semittance[k]>yMax)
      yMax = semittance[k];

    if(sErms[k]<yMin)
      yMin = sErms[k];
    
    if(sErms[k]>yMax)
      yMax = sErms[k];
  }

  for(Int_t k=1;k<=xNbin;k++) {
    Double_t value = hX1->GetBinContent(k);
    if(value<yMin)
      yMin = value;
    
    if(value>yMax)
      yMax = value;

  }

  yMax *= 1.1;

  // Plotting
  // -----------------------------------------------
    
  cout << "\n7. Plotting... " << endl;

  // Canvas setup
  // Create the canvas and the pads before the Frame loop
  // Resolution:
  Int_t sizex = 800;
  Int_t sizey = 600;
  if(opt.Contains("hres")) {
    Int_t sizex = 1600;
    Int_t sizey = 1200;    
  }
  
  TCanvas *C = new TCanvas("C1","Evolution of Injection",sizex,sizey);
  C->cd();

  // Set palette:
  PPalette * pPalette = (PPalette*) gROOT->FindObject("electron");
  pPalette->cd();

  // Float_t Max  = hP1X1->GetMaximum();
  // Float_t Min  = hP1X1->GetMinimum();
  
  // hP1X1->GetZaxis()->SetRangeUser(Min,Max); 

  // Text objects
  TPaveText *textTime =  new TPaveText(0.55,0.76,0.82,0.86,"NDC");
  PlasmaGlob::SetPaveTextStyle(textTime,32); 
  textTime->SetTextColor(kGray+2);
  char ctext[128];
  if(opt.Contains("units") && pData->GetPlasmaDensity()) 
    sprintf(ctext,"z = %5.1f mm", Time * skindepth / PUnits::mm);
  else
    sprintf(ctext,"t = %5.1f #omega_{p}^{-1}",Time);
  textTime->AddText(ctext);
 
  TPaveText *textDen = new TPaveText(0.15,0.85,0.48,0.9,"NDC");
  PlasmaGlob::SetPaveTextStyle(textDen,12); 
  textDen->SetTextColor(kOrange+10);
  if(opt.Contains("units") && pData->GetPlasmaDensity())
    sprintf(ctext,"n_{0} = %5.2f x 10^{17} / cc", n0 / (1e17/PUnits::cm3));
  else if(pData->GetBeamDensity() && pData->GetPlasmaDensity())
    sprintf(ctext,"n_{b}/n_{0} = %5.2f", pData->GetBeamDensity()/n0);
  textDen->AddText(ctext);

  TPaveText *textCharge = new TPaveText(0.15,0.25,0.48,0.3,"NDC");
  PlasmaGlob::SetPaveTextStyle(textCharge,12); 
  textCharge->SetTextColor(kGray+2);
  if(opt.Contains("units") && pData->GetPlasmaDensity())
    sprintf(ctext,"Q = %5.2f pC", Charge);
  else
    sprintf(ctext,"Q = %5.2f n0#timeskp^{-3}", Charge);    
  textCharge->AddText(ctext);

  TPaveText *textMom = new TPaveText(0.55,0.03,0.82,0.13,"NDC");
  PlasmaGlob::SetPaveTextStyle(textMom,32); 
  textMom->SetTextColor(kGray+3);
  textMom->SetTextFont(62);
  if(opt.Contains("units") && pData->GetPlasmaDensity())
    sprintf(ctext,"#LTp_{z}#GT = %5.2f GeV/c", pzmean);
  else
    sprintf(ctext,"#LTp_{z}#GT = %5.2f mc", pzmean);    
  textMom->AddText(ctext);


  TPaveText *textInfo = new TPaveText(0.55,0.45,0.82,0.75,"NDC");
  PlasmaGlob::SetPaveTextStyle(textInfo,32); 
  textInfo->SetTextColor(kGray+2);
  textInfo->SetTextFont(42);
  sprintf(ctext,"Q = %5.2f pC",Charge);
  textInfo->AddText(ctext);
  sprintf(ctext,"#LT#zeta#GT_{rms} = %5.2f #mum",zrms);
  textInfo->AddText(ctext);
  sprintf(ctext,"#LTp_{z}#GT_{rms} = %5.2f GeV/c",pzrms);
  textInfo->AddText(ctext);
  sprintf(ctext,"#epsilon_{y} = %5.2f #mum",emity);
  textInfo->AddText(ctext);
  
  // Setup Pad layout: 
  const Int_t NPad = 2;
  TPad *pad[NPad];
  TH1F *hFrame[NPad];
  TString sLabels[] = {"(b)","(a)"};
  TPaveText **textLabel = new TPaveText*[NPad];
  
  Float_t lMargin = 0.15;
  Float_t rMargin = 0.18;
  Float_t bMargin = 0.15;
  Float_t tMargin = 0.04;
  Float_t factor = 1.0;    
  PlasmaGlob::CanvasAsymPartition(C,NPad,lMargin,rMargin,bMargin,tMargin,factor);

    // Define the frames for plotting
  Int_t fonttype = 43;
  Int_t fontsize = 24;
  Int_t tfontsize = 28;
  Float_t txoffset = 2.0;
  Float_t lxoffset = 0.02;
  Float_t tyoffset = 1.3;
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
    hFrame[i] = (TH1F*) gROOT->FindObject(name);
    if(hFrame[i]) delete hFrame[i];
    hFrame[i] = (TH1F*) hX1->Clone(name);
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
  }

  C->cd(0);
  pad[1]->Draw();
  pad[1]->cd(); // <---------------------------------------------- Top Plot ---------
  if(opt.Contains("logz")) {
    pad[1]->SetLogz(1);
  } else {
    pad[1]->SetLogz(0);
  }
  
  hFrame[1]->GetYaxis()->SetRangeUser(hP1X1->GetYaxis()->GetXmin(),hP1X1->GetYaxis()->GetXmax());
  
  if(opt.Contains("units"))
    hFrame[1]->GetYaxis()->SetTitle("p_{z} [GeV/c]");
  
  hFrame[1]->Draw();
  
  gP1left->SetLineWidth(2);
  gP1left->Draw("F");
  gP1left->Draw("L");
  
  TLine lZmean(zmean,hP1X1->GetYaxis()->GetXmin(),zmean,hP1X1->GetYaxis()->GetXmax());
  lZmean.SetLineColor(kGray+2);
  lZmean.SetLineStyle(2);
  lZmean.Draw();

  TLine lPmean(hP1X1->GetXaxis()->GetXmin(),pzmean,hP1X1->GetXaxis()->GetXmax(),pzmean);
  lPmean.SetLineColor(kGray+2);
  lPmean.SetLineStyle(2);
  lPmean.Draw();

  // 2D histogram z range
  Float_t dmax = hP1X1->GetMaximum();
  Float_t dmin = 0.0;
  hP1X1->GetZaxis()->SetRangeUser(1.1*dmin,dmax);

  hP1X1->GetZaxis()->SetTitleFont(fonttype);

  hP1X1->Draw("colzsame");
  // hP1X1->SetContour(20);
  // hP1X1->Draw("contzsame");
  // hP1X1prof->SetMarkerStyle(1);
  // hP1X1prof->SetLineWidth(2);
  // hP1X1prof->Draw("zsame");

  //hP1->Draw("C");
  
  gPad->Update();

  TPaletteAxis *palette = (TPaletteAxis*)hP1X1->GetListOfFunctions()->FindObject("palette");
  if(palette) {
    Float_t y1 = gPad->GetBottomMargin();
    Float_t y2 = 1 - gPad->GetTopMargin();
    Float_t x1 = 1 - gPad->GetRightMargin();
    palette->SetY2NDC(y2 - 0.04);
    palette->SetY1NDC(y1 + 0.04);
    palette->SetX1NDC(x1 + 0.01);
    palette->SetX2NDC(x1 + 0.04);
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


  if(!opt.Contains("notext")) {  
    textTime->Draw();
    textInfo->Draw();
    // textCharge->Draw();
    textMom->Draw();
  }

  Float_t y1 = gPad->GetBottomMargin();
  Float_t y2 = 1 - gPad->GetTopMargin();
  Float_t x1 = gPad->GetLeftMargin();
  Float_t x2 = 1 - gPad->GetRightMargin();
  Float_t yrange = y2-y1; 
  Float_t xrange = x2-x1; 
  
  textLabel[1] = new TPaveText(x1 + 0.02*(x2-x1), y2-0.2*(y2-y1), x1+0.30*(x2-x1), y2-0.05*(y2-y1),"NDC");
  PlasmaGlob::SetPaveTextStyle(textLabel[1],12); 
  textLabel[1]->SetTextFont(42);
  textLabel[1]->AddText(sLabels[1]);
  textLabel[1]->Draw();

  gPad->RedrawAxis(); 


  // Bottom plot -----------------------------------------
  C->cd(0);
  pad[0]->Draw();
  pad[0]->cd(); // <---------------------------------------------------------- Bottom Plot
  // if(opt.Contains("logz")) {
  //   pad[0]->SetLogz(1);
  // } else {
  //   pad[0]->SetLogz(0);
  // }

  hFrame[0]->GetYaxis()->SetRangeUser(0.0,1.1*yMax);
  hFrame[0]->Draw();

  hX1->GetYaxis()->SetNdivisions(503);
  hX1->SetLineWidth(2);
  hX1->SetFillStyle(1001);
  hX1->SetFillColor(PlasmaGlob::elecFill);
  // hX1->SetLineColor(kBlue);
  
  //hX1->Smooth();
  hX1->Draw("FL same");
  //hX1->Draw("C");

  TLine lZmean2(zmean,0.0,zmean,1.1*yMax);
  lZmean2.SetLineColor(kGray+2);
  lZmean2.SetLineStyle(2);
  lZmean2.Draw();

  Int_t markerSize = 1.2; 
  Int_t lineWidth  = 2.0;   

  gYrms->SetMarkerStyle(20);
  gYrms->SetLineStyle(1);
  gYrms->SetMarkerColor(kGray+1);
  gYrms->SetMarkerSize(markerSize); 
  gYrms->SetLineColor(kGray+1);
  gYrms->SetLineWidth(lineWidth);
  // gYrms->Draw("PL");
  
  // hP2X2sl[0]->Draw("colz");
  gemit->SetMarkerStyle(20);
  // gemit->SetMarkerColor(kMagenta-2);
  gemit->SetMarkerColor(kGray+2);
  gemit->SetMarkerSize(markerSize);
  gemit->SetLineWidth(lineWidth);
  gemit->SetLineColor(kGray+2);
  gemit->Draw("PL");

  gErms->SetMarkerStyle(20);
  gErms->SetMarkerSize(markerSize);
  gErms->SetMarkerColor(kOrange+10);
  gErms->SetLineColor(kOrange+10);
  gErms->SetLineWidth(lineWidth);
  gErms->Draw("PL");


  TLegend *Leg;
  Leg=new TLegend(0.55,0.75,1 - gPad->GetRightMargin() - 0.02,0.95);
    
  PlasmaGlob::SetPaveStyle(Leg);
  Leg->SetTextAlign(12);
  Leg->SetTextColor(kGray+3);
  Leg->SetTextFont(42);
  Leg->SetLineColor(1);
  Leg->SetBorderSize(0);
  Leg->SetFillColor(0);
  Leg->SetFillStyle(1001);
  Leg->SetFillStyle(0); // Hollow
 
  Leg->AddEntry(hX1  ,"Current [kA]","L");
  // Leg->AddEntry(gErms,"Energy spread (GeV)","PL");
  Leg->AddEntry(gErms,"Energy spread [%]","PL");
  Leg->AddEntry(gemit,"Emittance [#mum]","PL");
  // Leg->AddEntry(gYrms,"Bunch width [#mum]","PL");
 
  Leg->Draw();

  y1 = gPad->GetBottomMargin();
  y2 = 1 - gPad->GetTopMargin();
  x1 = gPad->GetLeftMargin();
  x2 = 1 - gPad->GetRightMargin();
  yrange = y2-y1; 
  xrange = x2-x1; 

  textLabel[0] = new TPaveText(x1 + 0.02*(x2-x1), y2-0.2*(y2-y1), x1+0.30*(x2-x1), y2-0.05*(y2-y1),"NDC");
  PlasmaGlob::SetPaveTextStyle(textLabel[0],12); 
  textLabel[0]->SetTextFont(42);
  textLabel[0]->AddText(sLabels[0]);
  textLabel[0]->Draw();

  gPad->RedrawAxis(); 

  gPad->Update();
  

  // Print to file --------------------------------------
  
  C->cd();
  
  // Print to a file
  // Output file
  TString fOutName = Form("./%s/Plots/RakeBunch/RakeBunch",sim.Data());
  fOutName += Form("-%s_%i",sim.Data(),time);

  PlasmaGlob::imgconv(C,fOutName,opt);
  // ---------------------------------------------------------

  gStyle->SetOptStat(1);
  gStyle->cd();

  TCanvas *CA4 = new TCanvas("CA4","Sliced p2-x2",600,800);
  CA4->cd();

  Int_t ndiv = 4;
  CA4->Divide(1,ndiv);

  TString fOutName2 = Form("./%s/Plots/RakeBunch/RakeBunch-p2x2-%s_%i",sim.Data(),sim.Data(),time);

  CA4->Print(fOutName2 + ".ps[","Portrait");

  hP2X2->GetXaxis()->SetLabelSize(0.08);
  hP2X2->GetXaxis()->SetTitleSize(0.08);
  hP2X2->GetXaxis()->SetTitleOffset(1.0);
  hP2X2->GetXaxis()->CenterTitle();
  
  hP2X2->GetYaxis()->SetLabelSize(0.08);
  hP2X2->GetYaxis()->SetTitleSize(0.08);
  hP2X2->GetYaxis()->SetTitleOffset(0.8);
  hP2X2->GetYaxis()->CenterTitle();

  hP2X2->GetZaxis()->SetLabelSize(0.08);
  hP2X2->GetZaxis()->SetTitleSize(0.08);
  hP2X2->GetZaxis()->SetTitleOffset(0.8);
  hP2X2->GetZaxis()->CenterTitle();

  
  CA4->cd(1);
  hP2X2->Draw("colz");
    
  y1 = gPad->GetBottomMargin();
  y2 = 1 - gPad->GetTopMargin();
  x1 = gPad->GetLeftMargin();
  x2 = 1 - gPad->GetRightMargin();
  
  TPaveText *textStatInt = new TPaveText(x1+0.02,y2-0.40,x1+0.20,y2-0.05,"NDC");
  PlasmaGlob::SetPaveTextStyle(textStatInt,12); 
  textStatInt->SetTextColor(kGray+3);
  textStatInt->SetTextFont(42);
    
  char text[64];
  sprintf(text,"#LTy#GT_{rms} = %5.2f",xrms);
  textStatInt->AddText(text);
  sprintf(text,"#LTp_{y}#GT_{rms} = %5.2f",yrms);
  textStatInt->AddText(text);
  sprintf(text,"#epsilon_{y} = %5.2f",emity);
  textStatInt->AddText(text);
  textStatInt->Draw();

  TPaveText **textStat = new TPaveText*[SNbin];
  Int_t pnumber = 0;
  for(Int_t k=0;k<SNbin;k++) {
    pnumber++;
    Int_t ic = pnumber%ndiv;
    
    // new page
    if( ic==0 ) {
      CA4->cd(0);
      CA4->Clear();
      CA4->Divide(1,ndiv);
    }
    CA4->cd(ic+1);
    
    hP2X2sl[k]->GetXaxis()->SetLabelSize(0.08);
    hP2X2sl[k]->GetXaxis()->SetTitleSize(0.08);
    hP2X2sl[k]->GetXaxis()->SetTitleOffset(1.0);
    hP2X2sl[k]->GetXaxis()->CenterTitle();

    hP2X2sl[k]->GetYaxis()->SetLabelSize(0.08);
    hP2X2sl[k]->GetYaxis()->SetTitleSize(0.08);
    hP2X2sl[k]->GetYaxis()->SetTitleOffset(0.8);
    hP2X2sl[k]->GetYaxis()->CenterTitle();

    hP2X2sl[k]->GetZaxis()->SetLabelSize(0.08);
    hP2X2sl[k]->GetZaxis()->SetTitleSize(0.08);
    hP2X2sl[k]->GetZaxis()->SetTitleOffset(0.8);
    hP2X2sl[k]->GetZaxis()->CenterTitle();

    hP2X2sl[k]->Draw("colz");
    
    Float_t y1 = gPad->GetBottomMargin();
    Float_t y2 = 1 - gPad->GetTopMargin();
    Float_t x1 = gPad->GetLeftMargin();
    Float_t x2 = 1 - gPad->GetRightMargin();
    textStat[k] = new TPaveText(x1+0.02,y2-0.40,x1+0.20,y2-0.05,"NDC");
    PlasmaGlob::SetPaveTextStyle(textStat[k],12); 
    textStat[k]->SetTextColor(kGray+3);
    textStat[k]->SetTextFont(42);
    
    char text[64];
    sprintf(text,"#LTy#GT_{rms} = %5.2f",sxrms[k]);
    textStat[k]->AddText(text);
    sprintf(text,"#LTp_{y}#GT_{rms} = %5.2f",syrms[k]);
    textStat[k]->AddText(text);
    sprintf(text,"#epsilon_{y} = %5.2f",semittance[k]);
    textStat[k]->AddText(text);
    textStat[k]->Draw();
    
    if(ic+1==ndiv) {
      CA4->cd(0);
      CA4->Print(fOutName2 + ".ps");
    }
  }
  
  CA4->Print(fOutName2 + ".ps]");
  
  gSystem->Exec("ps2pdf " + fOutName2 + ".ps " + fOutName2 + ".pdf");
  gSystem->Exec("rm -rf " + fOutName2 + ".ps"); 
  
}
