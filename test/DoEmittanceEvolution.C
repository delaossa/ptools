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
#include <TGraphErrors.h>
#include <TMarker.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TMath.h>
#include <TPaveText.h>
#include <TPaletteAxis.h>
#include <TExec.h>
#include <TGaxis.h>

#include "PData.hh"
#include "PlasmaGlob.hh"
#include "PPalette.hh"

void DoEmittanceEvolution( const TString &sim, Int_t time, Int_t index = 0, const TString &phaname = "p1x1", const TString &options="") { 
  
#ifdef __CINT__  
  gSystem->Load("libplasma.so");
#endif

  PlasmaGlob::Initialize();
    
  // Palettes!
  gROOT->Macro("PlasmaPalettes.C");

  TString opt = options;
  // cout << "options = " << opt << endl;

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
  Double_t skindepth = 1.0;
  if(kp!=0.0) skindepth = 1/kp;
  Double_t E0 = pData->GetPlasmaE0();

  // Time in OU
  Float_t Time = pData->GetRealTime();
  // z start of the plasma in normalized units.
  Float_t zStartPlasma = pData->GetPlasmaStart() * kp;
  // z start of the beam in normalized units.
  Float_t zStartBeam = pData->GetBeamStart() * kp;
  
 
  opt += "comovcenter";

  if(opt.Contains("center")) {
    Time -= zStartPlasma;
    if(opt.Contains("comov"))      // Centers on the head of the beam.
      Time += zStartBeam;
  } 
  
  // Spatial coordinates intervals:
  Double_t x1Min = -4.3;
  Double_t x1Max = -3.9;
  Double_t x2Min = -0.5;
  Double_t x2Max =  0.5;
  Double_t x3Min = -0.5;
  Double_t x3Max =  0.5;

  if(sim.Contains("DR")) {
    x1Min = -4.0;
    x1Max = -3.4;
  } else if(sim.Contains("flash") && sim.Contains(".G.") ) {
    x1Min = -6.0;
    x1Max = -5.0;
  } else if(sim.Contains("facet_v23kA.G")) {
    x1Min = -7.7;
    x1Max = -7.0;
  }

  char cutString[128];
  sprintf(cutString,"TMath::Abs(q)*(x1 > %f && x1 < %f && x2 > %f && x2 < %f && x3 > %f && x3 < %f)",x1Min,x1Max,x2Min,x2Max,x3Min,x3Max); 
  TCut Cut = cutString;

  // Bining, intervals, labels, etc.
  Int_t xNbin = 250;
  Int_t yNbin = 250;
  Double_t xMin,yMin,xMax,yMax;
  xMin = yMin = xMax = yMax = 0.0;

  // Labeling...
  char xname[24];
  char yname[24];
  char xlabel[24];
  char ylabel[24];
  char alabel[24];
  char tlabel[24];

  if(phaname.Contains("p1x1")) {
    xMin = x1Min;
    xMax = x1Max;
    yMin =  0.0;
    yMax =  10000.0; 
    sprintf(xname,"x1");
    sprintf(yname,"p1");
    sprintf(xlabel,"#zeta");
    sprintf(ylabel,"p_{z}");
    sprintf(alabel,"#epsilon_{z}");   
  } else if(phaname.Contains("p2x2")) {
    xMin =  x2Min;
    xMax =  x2Max;
    yMin = -20.0;
    yMax =  20.0;
    sprintf(xname,"x2");
    sprintf(yname,"p2");
    sprintf(xlabel,"y");
    sprintf(ylabel,"p_{y}");
    sprintf(alabel,"#epsilon_{y}");
  } else if(phaname.Contains("p3x3")) {
    xMin =  x3Min;
    xMax =  x3Max;
    yMin = -20.0;
    yMax =  20.0;
    sprintf(xname,"x3");
    sprintf(yname,"p3");
    sprintf(xlabel,"x");
    sprintf(ylabel,"p_{x}");
    sprintf(alabel,"#epsilon_{x}");
  } else if(phaname.Contains("x2x1")){
    xMin =  x1Min;
    xMax =  x1Max;
    yMin =  x2Min;
    yMax =  x3Max;
    sprintf(xname,"x1");
    sprintf(yname,"x2");
    sprintf(xlabel,"#zeta");
    sprintf(ylabel,"y");
  }
  sprintf(tlabel,"time");

  char xunits[24];
  char yunits[24];
  char aunits[24];
  char tunits[24];
  if(opt.Contains("units")) {
    sprintf(xunits,"#mum");
    sprintf(tunits,"mm");
    if(phaname.Contains("p1x1")) {
      sprintf(yunits,"MeV/c");
      //      sprintf(aunits,"(MeV/c) #mum");
      sprintf(aunits,"mm mrad");
    } else if(phaname.Contains("p2x2") || phaname.Contains("p3x3")) {
      sprintf(yunits,"MeV/c");
      //sprintf(aunits,"(MeV/c) #mum");
      sprintf(aunits,"mm mrad");
    } else if(phaname.Contains("x2x1")) {
      sprintf(yunits,"#mum");
      sprintf(aunits,"#mum^{2}");
    }
  } else {
    sprintf(xunits,"c/#omega_{p}");
    sprintf(tunits,"#omega_{p}^{-1}");
    if(phaname.Contains("p1x1") || phaname.Contains("p2x2") || phaname.Contains("p3x3")) {
      sprintf(yunits,"m_{e}c");
      sprintf(aunits,"m_{e}c^{2}/#omega_{p}");
    } else if(phaname.Contains("x2x1")) {
      sprintf(yunits,"c/#omega_{p}");
      sprintf(aunits,"(c/#omega_{p})^{2}");    
    }
  }
  
  char xaxisname[24];
  char yaxisname[24];
  sprintf(xaxisname,"%s [%s]",xlabel,xunits);
  sprintf(yaxisname,"%s [%s]",ylabel,yunits);

  // Get phasespace histos
  Int_t Nspecies = pData->NSpecies();
  TH1F *hPha1Dx = NULL;
  TH1F *hPha1Dy = NULL;
  TH2F *hPha2D = NULL;
  Double_t Charge = 0.;
  for(Int_t i=0;i<Nspecies;i++) {
    
    if(i!=index) continue;

    if(!pData->GetRawFileName(i))
      continue;    
    
    char hName[24];
    sprintf(hName,"hPha_%i_%s",i,yname);
    hPha1Dy = (TH1F*) gROOT->FindObject(hName);
    if(hPha1Dy) delete hPha1Dy;
    hPha1Dy = new TH1F(hName,"",yNbin,yMin,yMax);
    // cout << "filename = " << pData->GetRawFileName(i)->c_str() << endl;
    //    pData->GetH1Raw(pData->GetRawFileName(i)->c_str(),yname,hPha1Dy,opt);
    TTree *tree = pData->GetTreeRaw(pData->GetRawFileName(i)->c_str(),opt);
    tree->Project(hName,yname,Cut);
    
    // Explore the 1D histogram in the 2nd variable :
    // It gets the mean the rms and the boundaries for plotting
    Double_t yNtotal = 0.0;
    Double_t ymean = 0.0;
    Double_t ymean2 = 0.0;
    Double_t yLeft = -1;
    Double_t yRight = -1;
    Int_t ybinLeft = -1;
    Int_t ybinRight = -1;
    for(Int_t j=1;j<=yNbin;j++) 
      if(yLeft==-1 && hPha1Dy->GetBinContent(j)>0.001*hPha1Dy->GetMaximum()) {
	yLeft = hPha1Dy->GetBinCenter(j); 
	ybinLeft = j;
	break;
      } 
    
    for(Int_t j=yNbin;j>0;j--) 
      if(yRight==-1 && hPha1Dy->GetBinContent(j)>0.001*hPha1Dy->GetMaximum()) {
	yRight = hPha1Dy->GetBinCenter(j); 
	ybinRight = j;
	break;
      } 

    for(Int_t k=ybinLeft;k<=ybinRight;k++) {
      Double_t y = hPha1Dy->GetXaxis()->GetBinCenter(k);
      Double_t value = TMath::Abs(hPha1Dy->GetBinContent(k));
      ymean += y*value;
      ymean2 += y*y*value;
      yNtotal += value;
    }
    ymean  /= yNtotal;
    ymean2  /= yNtotal;
    Double_t yrms  = TMath::Sqrt(ymean2 - ymean*ymean);  
      
    Charge = hPha1Dy->Integral();

    // Redefine axis range according with the mean and the rms of the 1D distributions:
    TString YNAME = yname;
    if(YNAME.Contains("p1") || YNAME.Contains("p2") || YNAME.Contains("p3")) {
      yMin = yLeft  - TMath::Abs(yRight-yLeft)*0.5;
      yMax = yRight + TMath::Abs(yRight-yLeft)*0.5;
    }
    
    // Center the plot around 0
    if(phaname.Contains("p2x2") || phaname.Contains("p3x3")) {
      if(fabs(xMin)>fabs(xMax))
	xMax = -xMin;
      else
	xMin = -xMax;
      
      if(fabs(yMin)>fabs(yMax))
	yMax = -yMin;
      else
	yMin = -yMax;
    }

    sprintf(hName,"hPha_%i_%s",i,phaname.Data());
    hPha2D = (TH2F*) gROOT->FindObject(hName);
    if(hPha2D) delete hPha2D;
    hPha2D = new TH2F(hName,"",xNbin,xMin,xMax,yNbin,yMin,yMax);
    //pData->GetH2Raw(pData->GetRawFileName(i)->c_str(),xname,yname,hPha2D,opt);
    char dcommand[6];
    sprintf(dcommand,"%s:%s",yname,xname);
    tree->Project(hName,dcommand,Cut);
   
    hPha2D->GetXaxis()->CenterTitle();
    hPha2D->GetYaxis()->CenterTitle();
    hPha2D->GetZaxis()->CenterTitle();
    hPha2D->GetXaxis()->SetTitle(xaxisname);
    hPha2D->GetYaxis()->SetTitle(yaxisname);
    hPha2D->GetZaxis()->SetTitle("n [a.u.]");
  }

  // Emittance and others 
  Double_t xmean = 0.0;
  Double_t ymean = 0.0;
  Double_t x2mean = 0.0;
  Double_t y2mean = 0.0;
  Double_t xymean = 0.0;
  Double_t Ntotal = 0.0;
  for(Int_t i=1;i<=xNbin;i++) {
    Double_t x = hPha2D->GetXaxis()->GetBinCenter(i);
    // if(x<xmin || x>xmax) continue;
    for(Int_t j=1;j<=yNbin;j++) {
      Double_t y = hPha2D->GetYaxis()->GetBinCenter(j);
      // if(y<ymin || y>ymax) continue;
      Double_t value = TMath::Abs(hPha2D->GetBinContent(i,j));
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

  // Converts the charge from RAW data.
  Double_t dx1 = 0.05;
  Double_t dx2 = 0.06;
  Double_t dx3 = 0.06;
  Charge *= dx1*dx2*dx3;

  // Set palette:
  PPalette * pPalette = (PPalette*) gROOT->FindObject("electron");
  pPalette->cd();

  Float_t Max  = hPha2D->GetMaximum();
  Float_t Min  = hPha2D->GetMinimum();
  
  hPha2D->GetZaxis()->SetRangeUser(Min,Max); 

  
  // Chaning to user units: 
  // --------------------------
  
  if(opt.Contains("units") && n0) {

    Time *= skindepth / PUnits::mm;
    
    xMin *= skindepth / PUnits::um;
    xMax *= skindepth / PUnits::um;
    xmean *= skindepth / PUnits::um;
    xrms  *= skindepth / PUnits::um;
    if(phaname.Contains("p1x1")) {
      yMin *= pData->GetBeamMass() / PUnits::GeV;
      yMax *= pData->GetBeamMass() / PUnits::GeV;
      ymean *= pData->GetBeamMass() / PUnits::GeV;
      yrms  *= pData->GetBeamMass() / PUnits::GeV;
      emittance *= (skindepth / PUnits::um); // * (pData->GetBeamMass() / PUnits::MeV);
    } else if(phaname.Contains("p2x2")||phaname.Contains("p3x3")) {
      yMin *= pData->GetBeamMass() / PUnits::MeV;
      yMax *= pData->GetBeamMass() / PUnits::MeV;
      ymean *= pData->GetBeamMass() / PUnits::MeV;
      yrms  *= pData->GetBeamMass() / PUnits::MeV;
      emittance *= (skindepth / PUnits::um); // * (pData->GetBeamMass() / PUnits::MeV);
    } else if(phaname.Contains("x2x1")) {
      yMin *= skindepth / PUnits::um;
      yMax *= skindepth / PUnits::um;
      ymean *= skindepth / PUnits::um;
      yrms  *= skindepth / PUnits::um;
      emittance *= (skindepth / PUnits::um) * (skindepth / PUnits::um);
    }
    hPha2D->SetBins(xNbin,xMin,xMax,yNbin,yMin,yMax);
    // for(Int_t j=0;j<=xNbin;j++) {
    //   for(Int_t k=0;k<=yNbin;k++) {
    // 	hPha2D->SetBinContent(j,k, hPha2D->GetBinContent(j,k) * n0 / (1e15/PUnits::cm3) );
    //   }
    // }
        
    Double_t dV = skindepth * skindepth * skindepth;
    Charge *= n0 * dV * PConst::ElectronCharge / PUnits::picocoulomb; 
  }
  
  cout << " Bunch properties: " << endl;
  cout << Form("  xMean = %7.3f   yMean = %7.3f   xyMean = %7.3f",xmean,ymean,xymean) << endl;
  cout << Form("  xRms  = %7.3f   yRms  = %7.3f",xrms,yrms) << endl;
  cout << Form("  Emittance = %7.3f",emittance) << endl;
  cout << Form("  Charge = %7.3f",Charge) << endl;
  
  
  // OUTPUT ROOT FILE WITH THE PLOTS:
  TString filename = Form("./%s/Plots/EmittanceEvolution/Evolutions-%s-%s.root",sim.Data(),sim.Data(),phaname.Data());
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
  
  Int_t nPoints = 0;
  char gName[32];
  sprintf(gName,"gEmitvsTime");     
  TGraph *gEmitvsTime = NULL;
  gEmitvsTime = (TGraph*) ifile->Get(gName);
  if(gEmitvsTime==NULL) {
    gEmitvsTime = new TGraph();
    gEmitvsTime->SetName(gName);
    nPoints = 0;
    // Some cosmetics at creation time:
    gEmitvsTime->SetLineWidth(3);
    gEmitvsTime->SetLineColor(PlasmaGlob::fieldLine);
    gEmitvsTime->SetMarkerStyle(20);
    gEmitvsTime->SetMarkerSize(0.4);
    gEmitvsTime->SetMarkerColor(PlasmaGlob::fieldLine);	
  } else {
    nPoints = gEmitvsTime->GetN(); 
  }  
  
  gEmitvsTime->Set(nPoints+1);
  gEmitvsTime->SetPoint(nPoints,Time,emittance);
  gEmitvsTime->Write(gName,TObject::kOverwrite);
  
  
  sprintf(gName,"gYmeanvsTime");     
  TGraph *gYmeanvsTime = NULL;
  gYmeanvsTime = (TGraph*) ifile->Get(gName);
  if(gYmeanvsTime==NULL) {
    gYmeanvsTime = new TGraph();
    gYmeanvsTime->SetName(gName);
    nPoints = 0;
    // Some cosmetics at creation time:
    gYmeanvsTime->SetLineWidth(3);
    gYmeanvsTime->SetLineColor(PlasmaGlob::fieldLine);
    gYmeanvsTime->SetMarkerStyle(20);
    gYmeanvsTime->SetMarkerSize(0.4);
    gYmeanvsTime->SetMarkerColor(PlasmaGlob::fieldLine);	
  } else {
    nPoints = gYmeanvsTime->GetN(); 
  }  
  
  gYmeanvsTime->Set(nPoints+1);
  gYmeanvsTime->SetPoint(nPoints,Time,ymean);
  gYmeanvsTime->Write(gName,TObject::kOverwrite);
  
  sprintf(gName,"gYrmsvsTime");     
  TGraph *gYrmsvsTime = NULL;
  gYrmsvsTime = (TGraph*) ifile->Get(gName);
  if(gYrmsvsTime==NULL) {
    gYrmsvsTime = new TGraph();
    gYrmsvsTime->SetName(gName);
    nPoints = 0;
    // Some cosmetics at creation time:
    gYrmsvsTime->SetLineWidth(3);
    gYrmsvsTime->SetLineColor(PlasmaGlob::fieldLine);
    gYrmsvsTime->SetMarkerStyle(20);
    gYrmsvsTime->SetMarkerSize(0.4);
    gYrmsvsTime->SetMarkerColor(PlasmaGlob::fieldLine);	
  } else {
    nPoints = gYrmsvsTime->GetN(); 
  }  
  
  gYrmsvsTime->Set(nPoints+1);
  gYrmsvsTime->SetPoint(nPoints,Time,yrms);
  gYrmsvsTime->Write(gName,TObject::kOverwrite);
  

  
  sprintf(gName,"gXmeanvsTime");     
  TGraph *gXmeanvsTime = NULL;
  gXmeanvsTime = (TGraph*) ifile->Get(gName);
  if(gXmeanvsTime==NULL) {
    gXmeanvsTime = new TGraph();
    gXmeanvsTime->SetName(gName);
    nPoints = 0;
    // Some cosmetics at creation time:
    gXmeanvsTime->SetLineWidth(3);
    gXmeanvsTime->SetLineColor(PlasmaGlob::fieldLine);
    gXmeanvsTime->SetMarkerStyle(20);
    gXmeanvsTime->SetMarkerSize(0.4);
    gXmeanvsTime->SetMarkerColor(PlasmaGlob::fieldLine);	
  } else {
    nPoints = gXmeanvsTime->GetN(); 
  }  
  
  gXmeanvsTime->Set(nPoints+1);
  gXmeanvsTime->SetPoint(nPoints,Time,xmean);
  gXmeanvsTime->Write(gName,TObject::kOverwrite);
  
  sprintf(gName,"gXrmsvsTime");     
  TGraph *gXrmsvsTime = NULL;
  gXrmsvsTime = (TGraph*) ifile->Get(gName);
  if(gXrmsvsTime==NULL) {
    gXrmsvsTime = new TGraph();
    gXrmsvsTime->SetName(gName);
    nPoints = 0;
    // Some cosmetics at creation time:
    gXrmsvsTime->SetLineWidth(3);
    gXrmsvsTime->SetLineColor(PlasmaGlob::fieldLine);
    gXrmsvsTime->SetMarkerStyle(20);
    gXrmsvsTime->SetMarkerSize(0.4);
    gXrmsvsTime->SetMarkerColor(PlasmaGlob::fieldLine);	
  } else {
    nPoints = gXrmsvsTime->GetN(); 
  }  
  
  gXrmsvsTime->Set(nPoints+1);
  gXrmsvsTime->SetPoint(nPoints,Time,xrms);
  gXrmsvsTime->Write(gName,TObject::kOverwrite);
  

  sprintf(gName,"gChargevsTime");     
  TGraph *gChargevsTime = NULL;
  gChargevsTime = (TGraph*) ifile->Get(gName);
  if(gChargevsTime==NULL) {
    gChargevsTime = new TGraph();
    gChargevsTime->SetName(gName);
    nPoints = 0;
    // Some cosmetics at creation time:
    gChargevsTime->SetLineWidth(3);
    gChargevsTime->SetLineColor(PlasmaGlob::fieldLine);
    gChargevsTime->SetMarkerStyle(20);
    gChargevsTime->SetMarkerSize(0.4);
    gChargevsTime->SetMarkerColor(PlasmaGlob::fieldLine);	
  } else {
    nPoints = gChargevsTime->GetN(); 
  }  
  
  gChargevsTime->Set(nPoints+1);
  gChargevsTime->SetPoint(nPoints,Time,fabs(Charge));
  gChargevsTime->Write(gName,TObject::kOverwrite);
    
  // ------------------------------------------------------------------------------------
  

  
  ifile->Close();
  
}

