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
  Float_t p2Min = -20.0;
  Float_t p2Max =  20.0;
  Float_t p3Min = -20.0;
  Float_t p3Max =  20.0;

  // Bining, intervals, labels, etc.
  Int_t xNbin = 200;
  Int_t yNbin = 200;

  if(sim.Contains("DR")) {
    xNbin = 200;
    yNbin = 200;
   
    // p1Min =  2200.01;
    // p1Max =  3399.99;
    // p1Min =  3500.01;
    // p1Max =  4999.99;
    p1Min =  1800.01;
    p1Max =  3999.99;

    // x1Min = -4.0;
    // x1Max = -3.2;
    x1Min = 8637.0;
    x1Max = 8639.0;
    x2Min = -0.5;
    x2Max =  0.5;
    x3Min = -0.5;
    x3Max =  0.5;
  } else if(sim.Contains("flash") && sim.Contains(".G.") ) {
    
    x1Min = -6.3;
    x1Max = -5.0;

    //p1Min =  0.001;
    p1Min =  650.001;
    p1Max =  1499.99;
    
  } else if(sim.Contains("facet_v23kA.G.A")) {
    x1Min = -7.8;
    x1Max = -7.1;

    x2Min = -1.0;
    x2Max =  1.0;
    x3Min = -1.0;
    x3Max =  1.0;
    // t=150
    p1Min =  1100.001;
    p1Max =  1999.99;
  } else if(sim.Contains("facet_v23kA.G")) {
    x1Min = -7.8;
    x1Max = -7.0;

    x2Min = -1.0;
    x2Max =  1.0;
    x3Min = -1.0;
    x3Max =  1.0;
    // t=60
    // p1Min =  350.01;
    // p1Max =  799.99;
    // t = 100
    p1Min =  700.01;
    p1Max =  1399.99;
    // t=105
    //p1Min =  800.001;
    //p1Max =  1399.99;
    // t=183
    // p1Min =  1500.001;
    // p1Max =  2499.99;
    // t=235
    // p1Min =  2000.01;
    // p1Max =  3499.99; 
    // t=310
    // p1Min =  2800.01;
    // p1Max =  4799.99; 
  }
  
  
  // Get phasespace histos
  Int_t Nspecies = pData->NSpecies();
  if(index>Nspecies-1) {
    return;
  }
  if(!pData->GetRawFileName(index)) {
    return;    
  }
  TH1F *hX1 = NULL;
  TH1F *hP1 = NULL;
  TH2F *hP1X1 = NULL;
  TH2F *hP2X2 = NULL;


  cout << Form("\n1. Getting data... ") << endl; 
  char cutString[512];
  sprintf(cutString,"TMath::Abs(q)*(%s > %.1f && %s < %.1f && x2 > %.1f && x2 < %.1f && x3 > %.1f && x3 < %.1f)",sshiftz.Data(),x1Min,sshiftz.Data(),x1Max,x2Min,x2Max,x3Min,x3Max); 
  TCut Cut = cutString;
  cout << Form("   (applied cut: \n %s)",cutString) << endl;
  TTree *tree = pData->GetTreeRaw(pData->GetRawFileName(index)->c_str(),opt);
  
  char hName[24];
  char dCommand[128];

  cout << Form("\n2. Dumping 1D histograms.. ") << endl;
  
  sprintf(hName,"hX1");
  hX1 = (TH1F*) gROOT->FindObject(hName);
  if(hX1) delete hX1;
  hX1 = new TH1F(hName,"",xNbin,x1Min,x1Max);
  sprintf(dCommand,"%s>>%s",sshiftz.Data(),hName);
  cout << Form("   - x1. ") << endl;    
  tree->Draw(dCommand,Cut,"goff");

  sprintf(hName,"hP1");
  hP1 = (TH1F*) gROOT->FindObject(hName);
  if(hP1) delete hP1;
  hP1 = new TH1F(hName,"",yNbin,p1Min,p1Max);
  sprintf(dCommand,"p1>>%s",hName);
  cout << Form("   - p1. ") << endl;    
  tree->Draw(dCommand,Cut,"goff");

  cout << Form("\n3. Dumping 2D histograms.. ") << endl;

  sprintf(hName,"hP1X1");
  hP1X1 = (TH2F*) gROOT->FindObject(hName);
  if(hP1X1) delete hP1X1;
  hP1X1 = new TH2F(hName,"",xNbin,x1Min,x1Max,yNbin,p1Min,p1Max);
  sprintf(dCommand,"p1:%s>>%s",sshiftz.Data(),hName);
  cout << Form("   - p1 vs. x1 ") << endl;    
  tree->Draw(dCommand,Cut,"goff");
  
  sprintf(hName,"hP2X2");
  hP2X2 = (TH2F*) gROOT->FindObject(hName);
  if(hP2X2) delete hP2X2;
  hP2X2 = new TH2F(hName,"",xNbin,x2Min,x2Max,yNbin,p2Min,p2Max);
  sprintf(dCommand,"p2:x2>>%s",hName);
  cout << Form("   - p2 vs. x2 ") << endl;    
  tree->Draw(dCommand,Cut,"goff");
  
  hX1->GetXaxis()->CenterTitle();
  hX1->GetYaxis()->CenterTitle();
  hX1->GetZaxis()->CenterTitle();
  hP1X1->GetXaxis()->CenterTitle();
  hP1X1->GetYaxis()->CenterTitle();
  hP1X1->GetZaxis()->CenterTitle();
  hP2X2->GetXaxis()->CenterTitle();
  hP2X2->GetYaxis()->CenterTitle();
  hP2X2->GetZaxis()->CenterTitle();
   
  // Integrated long. emittance:

  cout << Form("\n4. Calculating integrated quantities.. ") << endl;

  Double_t xmean = 0.0;
  Double_t ymean = 0.0;
  Double_t x2mean = 0.0;
  Double_t y2mean = 0.0;
  Double_t xymean = 0.0;
  Double_t Ntotal = 0.0;
  for(Int_t i=1;i<=xNbin;i++) {
    Double_t x = hP1X1->GetXaxis()->GetBinCenter(i);
    // if(x<xmin || x>xmax) continue;
    for(Int_t j=1;j<=yNbin;j++) {
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
  
  // cout << " xrms = " << xrms << endl;
  
  // Sliced emittance:
  // --------------------------------------------------------------------------

  cout << Form("\n5. Slicing ") << endl;


  // Bining for sliced quantities:
  // const Int_t SNbin = 7;
  // Float_t sBinLim[SNbin+1] = {-4.16,-4.14,-4.12,-4.10,-4.08,-4.06,-4.04,-4.02};
  // const Int_t SNbin = 7;
  // Float_t sBinLim[SNbin+1] = {-4.20,-4.17,-4.14,-4.12,-4.10,-4.08,-4.06,-4.00};
  // const Int_t SNbin = 8;
  // Float_t sBinLim[SNbin+1] = {-3.93,-3.87,-3.81,-3.75,-3.69,-3.63,-3.57,-3.51,-3.45};
  
  Int_t SNbin = 7;
  Float_t nsigma = 2;
  Float_t x1BinMin = -4.16;
  Float_t x1BinMax = -4.02;
  if(sim.Contains("DR")) {
    SNbin = 150; 
    nsigma = 1.4;
  } else if(sim.Contains("facet_v23kA.G.A")) {
    SNbin = 46;
    nsigma = 1;
    x1BinMin = -7.55;
    x1BinMax = -7.25;
  } else if(sim.Contains("facet_v23kA.G")) {
    SNbin = 10;
    nsigma = 1;
    x1BinMin = -7.55;
    x1BinMax = -7.15;
  }
  
  Float_t *sBinLim = new Float_t[SNbin+1];

  if(opt.Contains("rms")) {
    sBinLim[0] = xmean - nsigma*xrms;
    sBinLim[SNbin] = xmean + nsigma*xrms;
  } else {
    sBinLim[0] = x1BinMin;
    sBinLim[SNbin] = x1BinMax;
  }

  Float_t slbinSize = (sBinLim[SNbin] - sBinLim[0])/SNbin;
  
  for(Int_t i=1;i<SNbin;i++) {
    sBinLim[i] = sBinLim[i-1] + slbinSize;
  }
  
  TH1F **hP1sl = new TH1F*[SNbin];
  TH2F **hP2X2sl = new TH2F*[SNbin];

  cout << Form("\n   - Dumping in %i bins ",SNbin) << endl;

  
  for(Int_t k=0;k<SNbin;k++) {

    cout<< Form("k = %i : (x1 > %f && x1 < %f)",k,sBinLim[k],sBinLim[k+1]) << endl; 

    sprintf(hName,"hP2X2sl_%2i",k);
    hP2X2sl[k] = (TH2F*) gROOT->FindObject(hName);
    if(hP2X2sl[k]) delete hP2X2sl[k];
    hP2X2sl[k] = new TH2F(hName,"",xNbin,x2Min,x2Max,yNbin,p2Min,p2Max);

    char zCutString[128];
    sprintf(zCutString,"(%s > %f && %s < %f)",sshiftz.Data(),sBinLim[k],sshiftz.Data(),sBinLim[k+1]); 
    TCut zCut = zCutString;
    tree->Project(hName,"p2:x2",Cut + zCut);

    sprintf(hName,"hP1sl_%2i",k);
    hP1sl[k] = (TH1F*) gROOT->FindObject(hName);
    if(hP1sl[k]) delete hP1sl[k];
    hP1sl[k] = new TH1F(hName,"",yNbin,p1Min,p1Max);
    tree->Project(hName,"p1",Cut + zCut);

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
      for(Int_t j=1;j<=yNbin;j++) {
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
    
    for(Int_t i=1;i<=yNbin;i++) {
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
    
    
    cout << " Bunch properties: " << endl;
    cout << Form("  xMean = %7.3f   yMean = %7.3f",sxmean[k],symean[k]) << endl;
    cout << Form("  xRms  = %7.3f   yRms  = %7.3f",sxrms[k],syrms[k]) << endl;
    cout << Form("  Emittance = %7.3f",semittance[k]) << endl;

    cout << Form("  Emean = %7.3f   Erms = %7.3f",sEmean[k],sErms[k]) << endl;
    

  }

  // Charge
  Double_t dx1 = pData->GetDX(0);
  Double_t dx2 = pData->GetDX(1);
  Double_t dx3 = pData->GetDX(2);
  
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
  

  // Chaning to user units: 
  // --------------------------
  
  if(opt.Contains("units") && n0) {
    
    Int_t NbinsX = hP1X1->GetNbinsX();
    Double_t xMin = skindepth * hP1X1->GetXaxis()->GetXmin() / PUnits::um;
    Double_t xMax = skindepth * hP1X1->GetXaxis()->GetXmax() / PUnits::um;
    Int_t NbinsY = hP1X1->GetNbinsY();
    Double_t yMin = hP1X1->GetYaxis()->GetXmin() * pData->GetBeamMass() / PUnits::GeV;
    Double_t yMax = hP1X1->GetYaxis()->GetXmax() * pData->GetBeamMass() / PUnits::GeV;
    hP1X1->SetBins(NbinsX,xMin,xMax,NbinsY,yMin,yMax);
    // Converting electron density
    Double_t dVb = skindepth * skindepth * skindepth;
    Double_t dX = (xMax-xMin)/NbinsX; 
    Double_t dE = (yMax-yMin)/NbinsY; 
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
    
    hP1X1->GetZaxis()->SetTitle("dQ/d#zetadp_{z} [pC]");

    hP1->SetBins(NbinsY,yMin,yMax);
    hP1->GetYaxis()->SetTitle("p_{z} [GeV/c]");

    hX1->SetBins(NbinsX,xMin,xMax);
    Double_t binSize = (xMax - xMin)/NbinsX;

    Double_t dV = skindepth * skindepth * skindepth;
    Double_t  lightspeed =  PConst::c_light / (PUnits::um/PUnits::femtosecond);
    cout << Form("Speed of light = %f",lightspeed) << endl;
    hX1->Scale(TMath::Abs(n0 * dV * (PConst::ElectronCharge/PUnits::picocoulomb) * (lightspeed/binSize)));
    
    // hX1->Scale(TMath::Abs((PUnits::um/skindepth)*(PConst::ElectronCharge/PUnits::picocoulomb)*PConst::c_light));
    
    // hX1->GetYaxis()->SetTitle("I[kA]");
    hX1->GetYaxis()->SetTitle("");
    if(opt.Contains("comov"))
      hX1->GetXaxis()->SetTitle("#zeta [#mum]");
    else
      hX1->GetXaxis()->SetTitle("z [#mum]");
    
    
    xmean *= skindepth / PUnits::um;
    xrms  *= skindepth / PUnits::um;
    ymean *= pData->GetBeamMass() / PUnits::GeV;
    yrms  *= pData->GetBeamMass() / PUnits::GeV;
    
    emittance *= (skindepth / PUnits::um);
    
    for(Int_t k=0;k<SNbin;k++) {
      xbin[k] *= skindepth / PUnits::um;

      sxmean[k] *= skindepth / PUnits::um;
      sxrms[k]  *= skindepth / PUnits::um;
      symean[k] *= pData->GetBeamMass() / PUnits::MeV;
      syrms[k] *= pData->GetBeamMass() / PUnits::MeV;
      
      semittance[k] *= (skindepth / PUnits::um);

      sEmean[k] *= pData->GetBeamMass() / PUnits::GeV;
      sErms[k]  *= 100 * pData->GetBeamMass() / PUnits::GeV / ymean; //sEmean[k];
      // sErms[k]  *= pData->GetBeamMass() / PUnits::GeV;

    }

  }


  // Create the graph with the emittances:
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
    Double_t *yarray   = new Double_t[yNbin];
    Double_t *xarray   = new Double_t[yNbin];
    
    // This is for the right side:
    // Double_t xMax = x1Min + (x1Max-x1Min) * 0.9;
    // Double_t xMin = x1Max;
    // And this for left:
    Double_t xMin = hX1->GetXaxis()->GetXmin();
    Double_t xMax = hX1->GetXaxis()->GetXmin() + (hX1->GetXaxis()->GetXmax()
						  -hX1->GetXaxis()->GetXmin()) * 0.2;
    Double_t EneMax = hP1->GetMaximum();
    // cout << Form("  EneMax = %f ", EneMax) << endl;
 
    for(Int_t j=0; j<yNbin; j++) {
      yarray[j] = hP1->GetBinCenter(j+1);
      xarray[j] = ((xMax-xMin)/EneMax)*hP1->GetBinContent(j+1) + xMin;

      // cout << Form("  x = %f  y = %f ", xarray[j],yarray[j]) << endl;
    }

    gP1left = new TGraph(yNbin,xarray,yarray);
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


  // Plotting
  // -----------------------------------------------
    
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
  TPaveText *textTime = new TPaveText(0.55,0.8,0.82,0.9,"NDC");
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

  TPaveText *textWav = new TPaveText(0.15,0.2,0.48,0.25,"NDC");
  PlasmaGlob::SetPaveTextStyle(textWav,12); 
  textWav->SetTextColor(kGray+2);
  sprintf(ctext,"#lambda_{p} = %5.2f #mum", pData->GetPlasmaWaveLength() / PUnits::um);
  textWav->AddText(ctext);

  TPaveText *textCharge = new TPaveText(0.15,0.25,0.48,0.3,"NDC");
  PlasmaGlob::SetPaveTextStyle(textCharge,12); 
  textCharge->SetTextColor(kGray+2);
  if(opt.Contains("units") && pData->GetPlasmaDensity())
    sprintf(ctext,"Charge = %5.2f pC", Charge);
  else
    sprintf(ctext,"Charge = %5.2f n0#timeskp^{-3}", Charge);    
  textCharge->AddText(ctext);

  TPaveText *textMom = new TPaveText(0.55,0.03,0.82,0.13,"NDC");
  PlasmaGlob::SetPaveTextStyle(textMom,32); 
  textMom->SetTextColor(kGray+3);
  textMom->SetTextFont(62);
  if(opt.Contains("units") && pData->GetPlasmaDensity())
    sprintf(ctext,"#LTp_{z}#GT = %5.2f GeV/c", ymean);
  else
    sprintf(ctext,"Mom = %5.2f mc", ymean);    
  textMom->AddText(ctext);


  TPaveText *textInfo = new TPaveText(0.55,0.52,0.82,0.75,"NDC");
  PlasmaGlob::SetPaveTextStyle(textInfo,32); 
  textInfo->SetTextColor(kGray+2);
  textInfo->SetTextFont(42);
  sprintf(ctext,"Charge = %5.2f pC",Charge);
  textInfo->AddText(ctext);
  sprintf(ctext,"#LT#zeta#GT_{rms} = %5.2f #mum",xrms);
  textInfo->AddText(ctext);
  sprintf(ctext,"#LTp_{z}#GT_{rms} = %5.2f GeV/c",yrms);
  textInfo->AddText(ctext);
  // sprintf(ctext,"#epsilon_{N} = %5.2f #mum",emittance);
  // textInfo->AddText(ctext);
  
  // Setup Pad layout:
  const Int_t NFrames = 2;
  TPad *pad[NFrames];
  TH1F *hFrame[NFrames];
  TString sLabels[] = {"(a)","(b)"};
  TPaveText **textLabel = new TPaveText*[NFrames];
  
  Float_t lMargin = 0.15;
  Float_t rMargin = 0.18;
  Float_t bMargin = 0.15;
  Float_t tMargin = 0.04;
  Float_t factor = 1.0;    
  PlasmaGlob::CanvasAsymPartition(C,NFrames,lMargin,rMargin,bMargin,tMargin,factor);

  for(Int_t k=0;k<NFrames;k++) {

    char padname[16];
    sprintf(padname,"pad_%i",k);
    pad[k] = (TPad*) gROOT->FindObject(padname);

    Float_t vfactor = pad[0]->GetAbsHNDC()/pad[k]->GetAbsHNDC();

    char name[16];
    
    sprintf(name,"hFrame_%i",k);  
    hFrame[k] = (TH1F*) gROOT->FindObject(name);
    if(hFrame[k]) delete hFrame[k];
    hFrame[k] = (TH1F*) hX1->Clone(name);
    hFrame[k]->Reset();
    
    hFrame[k]->GetXaxis()->CenterTitle();
    hFrame[k]->GetYaxis()->CenterTitle();
    hFrame[k]->GetZaxis()->CenterTitle();
    hFrame[k]->SetLabelFont(42,"xyz");
    hFrame[k]->SetTitleFont(42,"xyz");
    
    hFrame[k]->SetNdivisions(505,"xyz");
    
    hFrame[k]->SetTickLength(0.04*vfactor,"x");
    hFrame[k]->SetTickLength(0.02/vfactor,"y");
    hFrame[k]->SetTickLength(0.02/vfactor,"z");
    
    hFrame[k]->GetYaxis()->SetLabelSize(0.08*vfactor);
    hFrame[k]->GetYaxis()->SetLabelOffset(0.02/vfactor);
    
    hFrame[k]->GetYaxis()->SetTitleSize(0.10*vfactor);
    hFrame[k]->GetYaxis()->SetTitleOffset(0.7/vfactor);

    if(k==0) {  
      hFrame[k]->GetXaxis()->SetLabelSize(0.08);
      hFrame[k]->GetXaxis()->SetLabelOffset(0.02);
      hFrame[k]->GetXaxis()->SetTitleSize(0.12);
      hFrame[k]->GetXaxis()->SetTitleOffset(0.7);
    } else {
      hFrame[k]->GetXaxis()->SetLabelSize(0.0);
      hFrame[k]->GetXaxis()->SetTitleSize(0.0);
    }
  }

  
  C->cd(0);
  char padname[16];
  sprintf(padname,"pad_%i",1);
  pad[1] = (TPad*) gROOT->FindObject(padname);
  pad[1]->Draw();
  pad[1]->cd(); // <---------------------------------------------- Top Plot ---------
  if(opt.Contains("logz")) {
    pad[1]->SetLogz(1);
  } else {
    pad[1]->SetLogz(0);
  }
  pad[1]->SetFrameLineWidth(3);  
  pad[1]->SetTickx(1);

  hFrame[1]->GetYaxis()->SetRangeUser(hP1X1->GetYaxis()->GetXmin(),hP1X1->GetYaxis()->GetXmax());
  
  if(opt.Contains("units"))
    hFrame[1]->GetYaxis()->SetTitle("p_{z} [GeV/c]");
  
  hFrame[1]->Draw();
  
  gP1left->SetLineWidth(2);
  gP1left->Draw("F");
  gP1left->Draw("L");

  TLine lZmean(xmean,hP1X1->GetYaxis()->GetXmin(),xmean,hP1X1->GetYaxis()->GetXmax());
  lZmean.SetLineColor(kGray+2);
  lZmean.SetLineStyle(2);
  lZmean.Draw();

  TLine lPmean(hP1X1->GetXaxis()->GetXmin(),ymean,hP1X1->GetXaxis()->GetXmax(),ymean);
  lPmean.SetLineColor(kGray+2);
  lPmean.SetLineStyle(2);
  lPmean.Draw();


  hP1X1->GetYaxis()->SetNdivisions(503);
  hP1X1->GetZaxis()->SetNdivisions(503);
  hP1X1->GetZaxis()->SetRangeUser(0.001*hP1X1->GetMaximum(),hP1X1->GetMaximum());
  hP1X1->GetZaxis()->SetLabelSize(0.05);
  hP1X1->GetZaxis()->SetTitleSize(0.04);
  hP1X1->GetZaxis()->SetTitleFont(42);

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
    
    palette->SetLabelFont(42);
    
    palette->SetLabelSize(0.08);
    //palette->SetLabelOffset(0.005/vfactor);
    palette->SetTitleSize(0.10);
    //    palette->SetTitleOffset(9999.0*vfactor);
    palette->SetTitleOffset(0.6);
    
    palette->SetBorderSize(2);
    palette->SetLineColor(1);
  }


  textTime->Draw();
  textInfo->Draw();
  // textCharge->Draw();
  textMom->Draw();
   
  gPad->RedrawAxis(); 

  // Bottom plot -----------------------------------------
  C->cd(0);
  sprintf(padname,"pad_%i",0);
  pad[0] = (TPad*) gROOT->FindObject(padname);
  pad[0]->Draw();
  pad[0]->cd(); // <---------------------------------------------------------- Middle Plot
  // if(opt.Contains("logz")) {
  //   pad[0]->SetLogz(1);
  // } else {
  //   pad[0]->SetLogz(0);
  // }
  pad[0]->SetFrameLineWidth(3);  
  pad[0]->SetTickx(1);

  hFrame[0]->GetYaxis()->SetRangeUser(0.0,1.1*yMax);
  hFrame[0]->Draw();

  hX1->SetLineWidth(2);
  hX1->SetFillStyle(1001);
  hX1->SetFillColor(PlasmaGlob::elecFill);
  // hX1->SetLineColor(kBlue);
  hX1->Draw("FL same");
  //hX1->Draw("C");

  TLine lZmean2(xmean,0.0,xmean,1.1*yMax);
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
  gYrms->Draw("PL");
  
  // hP2X2sl[0]->Draw("colz");
  gemit->SetMarkerStyle(20);
  //  gemit->SetMarkerColor(kMagenta-2);
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
  if(!sim.Contains("DR")) 
    Leg=new TLegend(0.55,0.60,1 - gPad->GetRightMargin() - 0.02,0.95);
  else
    Leg=new TLegend(gPad->GetLeftMargin() + 0.02, 1.0-gPad->GetTopMargin()-0.20,
		    gPad->GetLeftMargin() + 0.30, 1.0-gPad->GetTopMargin()-0.02);
  
  
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
  //  Leg->AddEntry(gErms,"Energy spread (GeV)","PL");
  Leg->AddEntry(gErms,"Energy spread [%]","PL");
  Leg->AddEntry(gemit,"Emittance [#mum]","PL");
  Leg->AddEntry(gYrms,"Bunch width [#mum]","PL");
 
  Leg->Draw();

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
  
}
