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
#include <TPaletteAxis.h>
#include <TExec.h>
#include <TGaxis.h>

#include "PData.hh"
#include "PGlobals.hh"
#include "PPalette.hh"

void DoEvolutions( const TString &sim, Int_t time, Int_t Nbins=1, const TString &options="") { 
  
#ifdef __CINT__  
  gSystem->Load("libptools.so");
#endif

  PGlobals::Initialize();
    
  // Palettes!
  gROOT->Macro("PPalettes.C");

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

  // Some initial beam properties:
  Float_t Ebeam = pData->GetBeamEnergy() * PUnits::MeV;
  Float_t gamma = pData->GetBeamGamma();
  Float_t vbeam = pData->GetBeamVelocity();
  
  Double_t rms0 = pData->GetBeamRmsY() * kp;
  if(CYL)  rms0 = pData->GetBeamRmsR() * kp;
  
  // Time in OU
  Float_t Time = pData->GetRealTime();
  // z start of the plasma in normalized units.
  Float_t zStartPlasma = pData->GetPlasmaStart() * kp;
  // z start of the beam in normalized units.
  Float_t zStartBeam = pData->GetBeamStart() * kp;
  
  if(opt.Contains("center")) {
    Time -= zStartPlasma;
    if(opt.Contains("comov"))      // Centers on the head of the beam.
      Time += zStartBeam;
  }

  // Beam charge 2D and 1D histogram (on-axis)
  // ------------------------------------------------------------------
  TH2F *hDen2D = NULL;
  if(pData->GetChargeFileName(1)) {
    char hName[24];
    sprintf(hName,"hDen2D");
    hDen2D = (TH2F*) gROOT->FindObject(hName);
    if(hDen2D) { delete hDen2D; hDen2D = NULL; }

    if(!ThreeD)
      hDen2D = pData->GetCharge(1,opt);
    else
      hDen2D = pData->GetCharge2DSliceZY(1,-1,1,opt+"avg");

    hDen2D->SetName(hName);
    hDen2D->GetXaxis()->CenterTitle();
    hDen2D->GetYaxis()->CenterTitle();
    hDen2D->GetZaxis()->CenterTitle();
    
    if(opt.Contains("comov"))
      hDen2D->GetXaxis()->SetTitle("k_{p}#zeta");
    else
      hDen2D->GetXaxis()->SetTitle("k_{p}z");
    
    if(CYL) 
      hDen2D->GetYaxis()->SetTitle("k_{p}r");
    else
      hDen2D->GetYaxis()->SetTitle("k_{p}y");

    hDen2D->GetZaxis()->SetTitle("n_{b}/n_{0}");

 
  }
  
  // Define ranges from the charge 2D histogram:
  // Binning for 2D histograms:
  // We get this values from the 2D density histogram.
  Int_t   x1Nbin    = hDen2D->GetNbinsX();
  Float_t x1Range   = (hDen2D->GetXaxis()->GetXmax() - hDen2D->GetXaxis()->GetXmin());
  Float_t x1Mid     = (hDen2D->GetXaxis()->GetXmax() + hDen2D->GetXaxis()->GetXmin())/2.;
  Float_t x1Min     = hDen2D->GetXaxis()->GetXmin();
  Float_t x1Max     = hDen2D->GetXaxis()->GetXmax();
  
  Int_t   x2Nbin    = hDen2D->GetNbinsY();      
  Float_t x2Range   = (hDen2D->GetYaxis()->GetXmax() - hDen2D->GetYaxis()->GetXmin());
  Float_t x2Mid     = (hDen2D->GetYaxis()->GetXmax() + hDen2D->GetYaxis()->GetXmin())/2.;
  Float_t x2Min     = x2Mid - x2Range/2;
  Float_t x2Max     = x2Mid + x2Range/2;
  
  if(Nbins==0) {
    Nbins = TMath::Nint(rms0 / hDen2D->GetYaxis()->GetBinWidth(1)) ;
    // cout << Form(" Rms0 = %6.2f  Dx = %6.2f  Nbins = %4i .", 
    // 	   rms0, hDen2D->GetYaxis()->GetBinWidth(1), Nbins) << endl;
  }
  
  // Slice width limits.
  Int_t FirstyBin = 0;
  Int_t LastyBin  = 0;
  if(!CYL) {
    FirstyBin = hDen2D->GetNbinsY()/2 + 1 - Nbins;
    LastyBin =  hDen2D->GetNbinsY()/2 + Nbins;
  } else {
    FirstyBin = 1; 
    LastyBin  = Nbins;
  }  


  // OUTPUT ROOT FILE WITH THE PLOTS:
  TString filename = Form("./%s/Plots/Evolutions/Evolutions-%s.root",sim.Data(),sim.Data());
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

  
  // Charge 1D histogram on axis
  TH1F *hDen1D = NULL;
  if(pData->GetChargeFileName(1)) {
    TString opth1 = opt;
    opth1 += "avg";
    
    char hName[24];
    sprintf(hName,"hDen1D");
    hDen1D = (TH1F*) gROOT->FindObject(hName);
    if(hDen1D) delete hDen1D;
    
    if(ThreeD) {
      hDen1D = pData->GetH1SliceZ3D(pData->GetChargeFileName(1)->c_str(),"charge",-1,Nbins,-1,Nbins,opth1.Data());
    } else if(CYL) { // Cylindrical: The first bin with r>0 is actually the number 1 (not the 0).
      hDen1D = pData->GetH1SliceZ(pData->GetChargeFileName(1)->c_str(),"charge",1,Nbins,opth1.Data());
    } else { // 2D cartesian
      hDen1D = pData->GetH1SliceZ(pData->GetChargeFileName(1)->c_str(),"charge",-1,Nbins,opth1.Data());
    }
    hDen1D->SetName(hName); 
    
    if(opt.Contains("comov"))
      hDen1D->GetXaxis()->SetTitle("k_{p}#zeta");
    else
      hDen1D->GetXaxis()->SetTitle("k_{p}z");
  
    hDen1D->GetYaxis()->SetTitle("n_{b}/n_{0}");
  }

  // On-axis beam density vs \zeta vs time! _________________________________
  TH2F *hDen1DvsTime = NULL; 
  if(hDen1D) {
    char hName[24];
    sprintf(hName,"hDen1DvsTime");
    TH2F *hDen1DvsTimeOld = (TH2F*) ifile->Get(hName);

    Int_t nBins   = 1;
    Float_t edge0 = Time-0.5;
    Float_t edge1 = Time+0.5;
    if(hDen1DvsTimeOld!=NULL) {
      nBins = hDen1DvsTimeOld->GetNbinsX()+1;
      Float_t binwidth =  (Time - hDen1DvsTimeOld->GetXaxis()->GetBinCenter(1))/(nBins-1);
      edge0 = hDen1DvsTimeOld->GetXaxis()->GetBinCenter(1) - binwidth/2.;
      edge1 = Time + binwidth/2.;
    }
    hDen1DvsTime = new TH2F("temp","",nBins,edge0,edge1,
		       	hDen1D->GetNbinsX(),
			hDen1D->GetBinLowEdge(1),
			hDen1D->GetBinLowEdge(hDen1D->GetNbinsX()+1));
    
    for(Int_t ix=1;ix<hDen1DvsTime->GetNbinsX();ix++) {
      for(Int_t iy=1;iy<hDen1DvsTime->GetNbinsY();iy++) {
	hDen1DvsTime->SetBinContent(ix,iy,hDen1DvsTimeOld->GetBinContent(ix,iy));
      }
    }  
    delete hDen1DvsTimeOld;
  
    // Fill last bin with the newest values.
    for(Int_t iy=1;iy<=hDen1D->GetNbinsX();iy++) {
      hDen1DvsTime->SetBinContent(nBins,iy,hDen1D->GetBinContent(iy));
    }   

    hDen1DvsTime->GetZaxis()->SetTitle("n_{b}/n_{0}");
    hDen1DvsTime->GetYaxis()->SetTitle("k_{p}#zeta");
    hDen1DvsTime->GetXaxis()->SetTitle("k_{p}z");
    hDen1DvsTime->GetZaxis()->CenterTitle();
    hDen1DvsTime->GetYaxis()->CenterTitle();
    hDen1DvsTime->GetXaxis()->CenterTitle();
    hDen1DvsTime->SetName(hName);

    // Change the range of z axis 
    Float_t Denmax = hDen1DvsTime->GetMaximum();
    hDen1DvsTime->GetZaxis()->SetRangeUser(0,Denmax); 
    hDen1DvsTime->Write(hName,TObject::kOverwrite);

  }

  // RMS (vs z) of the beam's charge distribution: 
  TProfile *hDen2Dprof = NULL;
  TH1F *hRms = NULL;
  Double_t axisPos = x2Mid;
  if(hDen2D) {
    TString pname = hDen2D->GetName();
    pname += "_pfx";
    
    hDen2Dprof =  (TProfile*) gROOT->FindObject(pname.Data());
    if(hDen2Dprof) { delete hDen2Dprof; hDen2Dprof = NULL; }
    hDen2Dprof = hDen2D->ProfileX("_pfx",1,-1,"s");
    
    hRms = (TH1F*) gROOT->FindObject("hRms");
    if(hRms) delete hRms;
    
    hRms = new TH1F("hRms","",x1Nbin,x1Min,x1Max);
    
    if(CYL) axisPos = 0.0;
    
    for(Int_t j=0;j<hRms->GetNbinsX();j++) {
      Double_t rms = 0;
      Double_t total = 0;
      for(Int_t k=1;k<=x2Nbin;k++) {
	Double_t value  = hDen2D->GetBinContent(j,k);
	Double_t radius = hDen2D->GetYaxis()->GetBinCenter(k) - axisPos;
	if(CYL) {
	  rms += radius*radius*radius*value;
	  total += radius*value;
	} else {
	  rms += radius*radius*value;
	  total += value;
	}
	// cout << Form(" (%i,%i) -> radius = %7.4f ,  density = %7.4f",j,k,radius,value) << endl;
      }
      
      rms /= total;
      rms = sqrt(rms);
      
      hRms->SetBinContent(j,rms); 
      
    }
    
    hRms->GetXaxis()->SetTitle("k_{p}z");
    if(opt.Contains("comov"))
      hRms->GetXaxis()->SetTitle("k_{p}#zeta");
    
    hRms->GetYaxis()->SetTitle("k_{p}#LTr#GT_{rms}");
  }
  
  // Transverse charge RMS vs \zeta vs time! _________________________________
  TH2F *hRmsvsTime = NULL; 
  if(hRms) {
    char hName[24];
    sprintf(hName,"hRmsvsTime");
    TH2F *hRmsvsTimeOld = (TH2F*) ifile->Get(hName);

    Int_t nBins   = 1;
    Float_t edge0 = Time-0.5;
    Float_t edge1 = Time+0.5;
    if(hRmsvsTimeOld!=NULL) {
      nBins = hRmsvsTimeOld->GetNbinsX()+1;
      Float_t binwidth =  (Time - hRmsvsTimeOld->GetXaxis()->GetBinCenter(1))/(nBins-1);
      edge0 = hRmsvsTimeOld->GetXaxis()->GetBinCenter(1) - binwidth/2.;
      edge1 = Time + binwidth/2.;
    }
    hRmsvsTime = new TH2F("temp","",nBins,edge0,edge1,
		       	hRms->GetNbinsX(),
			hRms->GetBinLowEdge(1),
			hRms->GetBinLowEdge(hRms->GetNbinsX()+1));
    
    for(Int_t ix=1;ix<hRmsvsTime->GetNbinsX();ix++) {
      for(Int_t iy=1;iy<hRmsvsTime->GetNbinsY();iy++) {
	hRmsvsTime->SetBinContent(ix,iy,hRmsvsTimeOld->GetBinContent(ix,iy));
      }
    }  
    delete hRmsvsTimeOld;
  
    // Fill last bin with the newest values.
    for(Int_t iy=1;iy<=hRms->GetNbinsX();iy++) {
      hRmsvsTime->SetBinContent(nBins,iy,hRms->GetBinContent(iy));
    }   

    hRmsvsTime->GetZaxis()->SetTitle("#LTr#GT_{rms}");
    hRmsvsTime->GetYaxis()->SetTitle("k_{p}#zeta");
    hRmsvsTime->GetXaxis()->SetTitle("k_{p}z");
    hRmsvsTime->GetZaxis()->CenterTitle();
    hRmsvsTime->GetYaxis()->CenterTitle();
    hRmsvsTime->GetXaxis()->CenterTitle();
    hRmsvsTime->SetName(hName);

    // Change the range of z axis
    Float_t Rmsmax = hRmsvsTime->GetMaximum();
    hRmsvsTime->GetZaxis()->SetRangeUser(0,Rmsmax); 
    hRmsvsTime->Write(hName,TObject::kOverwrite);

  }

  // INTEGRATED Beam's Charge:
  // Total charge vs time :
  TGraph *gQvsTime = NULL;
  if(hDen2D) {
    Double_t Q = 0;
    for(Int_t i=1;i<=x1Nbin;i++) {
      for(Int_t j=1;j<=x2Nbin;j++) {
	Double_t value  = hDen2D->GetBinContent(i,j);
	if(CYL) {
	  Double_t radius = hDen2D->GetYaxis()->GetBinCenter(j);
	  Q += radius * value;
	  // cout << Form(" (%i,%i) -> radius = %7.4f , value = %7.4f",i,j,radius,value) << endl;
	} else {
	  Q += value;
	}
      }    
    }
    Double_t xbinsize = hDen2D->GetXaxis()->GetBinWidth(1);
    Double_t ybinsize = hDen2D->GetYaxis()->GetBinWidth(1); 
    Q *= xbinsize * ybinsize;
    
    if(!CYL && !ThreeD) {
      Q *= TMath::Sqrt(2*TMath::Pi()) * rms0; 
    } else if(CYL) {
      Q *= 2*TMath::Pi();
    }
    
    if(opt.Contains("units")) {
      Double_t dV = skindepth * skindepth * skindepth;
      Q *= n0 * dV;
      Q *= (PConst::ElectronCharge/PUnits::picocoulomb); 
      cout << Form(" Integrated charge     = %8i pC", TMath::Nint(Q)) << endl;
    } else {
      cout << Form(" Integrated charge     = %8.4f n0 * kp^-3",Q) << endl;
    }
    
    Int_t nPoints = 0;
    char gName[32];
    sprintf(gName,"gQvsTime");     
    gQvsTime = (TGraph*) ifile->Get(gName);
    if(gQvsTime==NULL) {
      gQvsTime = new TGraph();
      gQvsTime->SetName(gName);
      nPoints = 0;
      // Some cosmetics at creation time:
      gQvsTime->SetLineWidth(3);
      gQvsTime->SetLineColor(PGlobals::fieldLine);
      gQvsTime->SetMarkerStyle(20);
      gQvsTime->SetMarkerSize(0.4);
      gQvsTime->SetMarkerColor(PGlobals::fieldLine);	
      gQvsTime->GetYaxis()->SetTitle("charge [n_{0}/k_{p}^{3}]");
      gQvsTime->GetXaxis()->SetTitle("k_{p}z");
    } else {
      nPoints = gQvsTime->GetN(); 
    }  
    
    gQvsTime->Set(nPoints+1);
    gQvsTime->SetPoint(nPoints,Time,Q);
    gQvsTime->Write(gName,TObject::kOverwrite);
  }
  
  // ------------------------------------------------------------------------------------
  

  // Longitudinal phasespace 
  Int_t  gNbin = 100;
  // Float_t gMin = 80;
  // Float_t gMax = 120;
  Float_t gMin = 43.07 - 1.2;
  Float_t gMax = 43.07 + 1.2;
  TH2F *hGvsZ = NULL;
  if(pData->GetRawFileName(1)) {
    char hName[24];
    sprintf(hName,"hGvsZ");
    hGvsZ = (TH2F*) gROOT->FindObject(hName);
    if(hGvsZ) { delete hGvsZ; hGvsZ = NULL; }
    hGvsZ = new TH2F(hName,"",x1Nbin,x1Min,x1Max,gNbin,gMin,gMax);
    pData->GetH2Raw(pData->GetRawFileName(1)->c_str(),"x1","gamma",hGvsZ,opt);
    
    hGvsZ->GetXaxis()->CenterTitle();
    hGvsZ->GetYaxis()->CenterTitle();
    hGvsZ->GetZaxis()->CenterTitle();
    hGvsZ->GetYaxis()->SetTitle("#gamma");
    if(opt.Contains("comov")) {
      hGvsZ->GetXaxis()->SetTitle("k_{p}#zeta");
      hGvsZ->GetZaxis()->SetTitle("dN/d#zetad#gamma [a.u.]");
    }  else {
      hGvsZ->GetXaxis()->SetTitle("k_{p}z");
      hGvsZ->GetZaxis()->SetTitle("dN/dzd#gamma [a.u.]");
    }    
  } else {
    cout << Form("--> No RAW data file is present for species 1") << endl;
  }

  TH2F *hGvsTime = NULL; 
  TProfile *hGvsZprof = NULL;
  TGraphErrors *gGvsZ = NULL;
  if(hGvsZ) {
    TString pname = hGvsZ->GetName();
    pname += "_pfx";
    hGvsZprof =  (TProfile*) gROOT->FindObject(pname.Data());
    if(hGvsZprof) delete hGvsZprof;

    hGvsZprof = hGvsZ->ProfileX("_pfx",1,-1,"s");

    gGvsZ = (TGraphErrors*) gROOT->FindObject("gGvsZ");
    if(gGvsZ) delete gGvsZ;

    Int_t Npoints = hGvsZprof->GetNbinsX();
    Double_t *x = new Double_t[Npoints];
    Double_t *y = new Double_t[Npoints];
    Double_t *ex = new Double_t[Npoints];
    Double_t *ey = new Double_t[Npoints];
    
    for(Int_t j=0;j<Npoints;j++) {
      x[j] = hGvsZprof->GetBinCenter(j);
      y[j] = hGvsZprof->GetBinContent(j);
      ex[j] = 0;
      ey[j] = hGvsZprof->GetBinError(j);   
    }
    
    gGvsZ = new TGraphErrors(Npoints,x,y,ex,ey);
    gGvsZ->SetName("gGvsZ");
        
    // PGlobals::SetH1Style((TH1*)gGvsZ,1);
    PGlobals::SetGraphStyle(gGvsZ,1);

   
    if(opt.Contains("comov")) 
      gGvsZ->GetXaxis()->SetTitle("k_{p}#zeta");
    else
      gGvsZ->GetXaxis()->SetTitle("k_{p}z");
    
    gGvsZ->GetYaxis()->SetTitle("#LT#gamma#GT [MeV]");

    char hName[24];
    sprintf(hName,"hGvsTime");
    TH2F *hGvsTimeOld = (TH2F*) ifile->Get(hName);

    Int_t nBins   = 1;
    Float_t edge0 = Time-0.5;
    Float_t edge1 = Time+0.5;
    if(hGvsTimeOld!=NULL) {
      nBins = hGvsTimeOld->GetNbinsX()+1;
      Float_t binwidth =  (Time - hGvsTimeOld->GetXaxis()->GetBinCenter(1))/(nBins-1);
      edge0 = hGvsTimeOld->GetXaxis()->GetBinCenter(1) - binwidth/2.;
      edge1 = Time + binwidth/2.;
    }
    hGvsTime = new TH2F("temp","",nBins,edge0,edge1,
		       	hGvsZprof->GetNbinsX(),
			hGvsZprof->GetBinLowEdge(1),
			hGvsZprof->GetBinLowEdge(hGvsZprof->GetNbinsX()+1));
    
    for(Int_t ix=1;ix<hGvsTime->GetNbinsX();ix++) {
      for(Int_t iy=1;iy<hGvsTime->GetNbinsY();iy++) {
	hGvsTime->SetBinContent(ix,iy,hGvsTimeOld->GetBinContent(ix,iy));
      }
    }  
    delete hGvsTimeOld;
  
    // Fill last bin with the newest values.
    for(Int_t iy=1;iy<=hGvsZprof->GetNbinsX();iy++) {
      hGvsTime->SetBinContent(nBins,iy,hGvsZprof->GetBinContent(iy));
    }   

    hGvsTime->GetZaxis()->SetTitle("#LT#gamma#GT");
    hGvsTime->GetYaxis()->SetTitle("k_{p}#zeta");
    hGvsTime->GetXaxis()->SetTitle("k_{p}z");
    hGvsTime->GetZaxis()->CenterTitle();
    hGvsTime->GetYaxis()->CenterTitle();
    hGvsTime->GetXaxis()->CenterTitle();
    hGvsTime->SetName(hName);

    // Change the range of z axis
    Float_t Gmax = hGvsTime->GetMaximum();
    Float_t Gmin = hGvsTime->GetMinimum();    
    hGvsTime->GetZaxis()->SetRangeUser(Gmin,Gmax); 
    hGvsTime->Write(hName,TObject::kOverwrite);
    
  }

  // ---------------------------------------------------------------------------------


  // EM fields on - axis :

  TString opth1 = opt;
  opth1 += "avg";
  // Get electric fields
  const Int_t Nfields = 2;
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
    } else if(CYL) { // Cylindrical: The first bin with r>0 is actually the number 1 (not the 0).
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
    sprintf(hName,"hE_%i_%i",i,time);
    hE1D[i]->SetName(hName);
    if(opt.Contains("comov"))
      hE1D[i]->GetXaxis()->SetTitle("k_{p}#zeta");
    else
      hE1D[i]->GetXaxis()->SetTitle("k_{p}z");
   
    if(i==0)
      hE1D[i]->GetYaxis()->SetTitle("E_{z}/E_{0}");
    else if(i==1)
      hE1D[i]->GetYaxis()->SetTitle("E_{y}/E_{0}");
    else if(i==2)
      hE1D[i]->GetYaxis()->SetTitle("E_{x}/E_{0}");
    
    hE1D[i]->GetYaxis()->CenterTitle();
    hE1D[i]->GetXaxis()->CenterTitle();
    
  }  
  
  // Calculate wave positions:
  // ----------------------------------------------------------------
    
  // Calculate the crossings and the extremes of the Electric fields
  Float_t Ecross[Nfields][100] = {{0.0}};
  Float_t Eextr[Nfields][100] = {{0.0}};
  Int_t Ncross[Nfields] = {0};

  for(Int_t i=0;i<Nfields;i++) {
    Ncross[i] = 0;
    
    if(!hE1D[i]) continue;
    
    // Only smooths the focusing if flag activated..
    if(i>0 && opt.Contains("smooth")) {
      // cout << " Smoothing fields on axis..." << endl;
      hE1D[i]->Smooth(10);
    } 

    Float_t maxZeta = zStartBeam;
    if(opt.Contains("center")) 
      maxZeta -= zStartBeam;
          
    for(Int_t ip=hE1D[i]->GetNbinsX();ip>1;ip--) {

      Float_t Z2 = hE1D[i]->GetBinCenter(ip-1);
      if(Z2 > maxZeta) continue;
      Float_t E1 = hE1D[i]->GetBinContent(ip);
      Float_t E2 = hE1D[i]->GetBinContent(ip-1);
      Float_t Z1 = hE1D[i]->GetBinCenter(ip);
      
      // cout << Form("Z1 = %6.4f  Z2 = %6.4f   E1 = %6.4f   E2 = %6.4f", Z1, Z2, E1, E2) << endl; 

      if(E1*E2 >= 0) { // No change of sign means we are in a side of the zero axis.
	if(fabs(E2)>fabs(Eextr[i][Ncross[i]])) {
	  Eextr[i][Ncross[i]] = E2;
	} 
      }
      
      if(E1*E2 < 0) { // change of sign means a crossing!
	
	// The next crossing has to be far enough from the previous one:
	Float_t zcross =  -E1 * ( (Z2-Z1)/(E2-E1) ) + Z1;
        if(Ncross[i]>0 && fabs(Ecross[i][Ncross[i]-1]-zcross)<TMath::PiOver2() ) continue;	
	// cout << " CROSS! " << endl;

	// add the point
	Ecross[i][Ncross[i]] = zcross;
	Ncross[i]++;
      }
    }
    
    cout << "  -> Number of crossings for field " << i << " : " << Ncross[i] << endl;
    for(Int_t ic=0;ic<Ncross[i];ic++) {
      //  cout << Form(" %2i:  zeta = %6.4f  E = %6.4f", ic, Ecross[i][ic], Eextr[i][ic]) << endl; 
    }  
    
    
    hE1D[i]->SetLineColor(kRed);
    hE1D[i]->Write(hE1D[i]->GetName(),TObject::kOverwrite);

  }
  
  // Get the Graphs and histos from file
  Int_t nPoints = 0;
  TGraph ***gEcross = new TGraph**[Nfields]; 
  TGraph ***gEextr  = new TGraph**[Nfields]; 
  TH2F **hEvsTime = new TH2F*[Nfields]; 
  for(Int_t i=0;i<Nfields;i++) {
    char hName[24];
    sprintf(hName,"hEvsTime_%i",i);
    TH2F *hEvsTimeOld = (TH2F*) ifile->Get(hName);
    Int_t nBins   = 1;
    Float_t edge0 = Time-0.5;
    Float_t edge1 = Time+0.5;
    if(hEvsTimeOld!=NULL) {
      nBins = hEvsTimeOld->GetNbinsX()+1;
      Float_t binwidth =  (Time - hEvsTimeOld->GetXaxis()->GetBinCenter(1))/(nBins-1);
      edge0 = hEvsTimeOld->GetXaxis()->GetBinCenter(1) - binwidth/2.;
      edge1 = Time + binwidth/2.;
    }
    hEvsTime[i] = new TH2F("temp","",nBins,edge0,edge1,
			hE1D[i]->GetNbinsX(),
			hE1D[i]->GetBinLowEdge(1),
			hE1D[i]->GetBinLowEdge(hE1D[i]->GetNbinsX()+1));
    
    for(Int_t ix=1;ix<hEvsTime[i]->GetNbinsX();ix++) {
      for(Int_t iy=1;iy<hEvsTime[i]->GetNbinsY();iy++) {
	hEvsTime[i]->SetBinContent(ix,iy,hEvsTimeOld->GetBinContent(ix,iy));
      }
    }  
    delete hEvsTimeOld;
  
    // Fill last bin with the newest values.
    for(Int_t iy=1;iy<=hE1D[i]->GetNbinsX();iy++) {
      hEvsTime[i]->SetBinContent(nBins,iy,hE1D[i]->GetBinContent(iy));
    }   

    if(i==0) 
      hEvsTime[i]->GetZaxis()->SetTitle("E_{z}/E_{0}");
    else if(i==1)
      hEvsTime[i]->GetZaxis()->SetTitle("E_{y}/E_{0}");
    else if(i==2)
      hEvsTime[i]->GetZaxis()->SetTitle("E_{x}/E_{0}");
  
    hEvsTime[i]->GetYaxis()->SetTitle("k_{p}#zeta");
    hEvsTime[i]->GetXaxis()->SetTitle("k_{p}z");
    hEvsTime[i]->GetZaxis()->CenterTitle();
    hEvsTime[i]->GetYaxis()->CenterTitle();
    hEvsTime[i]->GetXaxis()->CenterTitle();
    hEvsTime[i]->SetName(hName);

    // Change the range of z axis for the fields to be symmetric.
    Float_t Emax = hEvsTime[i]->GetMaximum();
    Float_t Emin = hEvsTime[i]->GetMinimum();
    if(Emax > TMath::Abs(Emin))
      Emin = -Emax;
    else
      Emax = -Emin;
    hEvsTime[i]->GetZaxis()->SetRangeUser(Emin,Emax); 
    
    hEvsTime[i]->Write(hName,TObject::kOverwrite);

    // ---

    gEcross[i] = new TGraph*[Ncross[i]];
    gEextr[i] = new TGraph*[Ncross[i]];
    char gName[24];
    Int_t ifail = 0;
    for(Int_t ic=0;ic<Ncross[i];ic++) {
      sprintf(gName,"gEcross_%i_%i",i,ic);     
      gEcross[i][ic] = (TGraph*) ifile->Get(gName);
      if(gEcross[i][ic]==NULL) {
	gEcross[i][ic] = new TGraph();
	gEcross[i][ic]->SetName(gName);
	nPoints = 0;
	// Some cosmetics at creation time:
	if(i==1) gEcross[i][ic]->SetLineStyle(2);
	else gEcross[i][ic]->SetLineStyle(1);
	gEcross[i][ic]->SetLineWidth(1);
	gEcross[i][ic]->SetLineColor(kGray+1);
	gEcross[i][ic]->SetMarkerStyle(20);
	gEcross[i][ic]->SetMarkerSize(0.4);
	gEcross[i][ic]->SetMarkerColor(kGray+1);	
	gEcross[i][ic]->GetYaxis()->SetTitle("k_{p}#zeta]");
	gEcross[i][ic]->GetXaxis()->SetTitle("k_{p}z");
      } else {
	nPoints = gEcross[i][ic]->GetN(); 
      }  
      
      // Check the new crossings respect the previous ones:
      // Double_t t,zeta;
      // if(nPoints>0) {
      // 	gEcross[i][ic]->GetPoint(nPoints-1,t,zeta);
      // 	if(fabs(zeta-Ecross[i][ic+ifail])>TMath::Pi()) {
      // 	  ic--;
      // 	  ifail++;
      // 	  continue;
      // 	}
      // }
      
      gEcross[i][ic]->Set(nPoints+1);
      gEcross[i][ic]->SetPoint(nPoints,Time,Ecross[i][ic+ifail]);
      gEcross[i][ic]->Write(gName,TObject::kOverwrite);
      
      // if(ic==Ncross[i]-1) continue;
      
      sprintf(gName,"gEextr_%i_%i",i,ic);     
      gEextr[i][ic] = (TGraph*) ifile->Get(gName);
      if(gEextr[i][ic]==NULL) {
	gEextr[i][ic] = new TGraph();
	gEextr[i][ic]->SetName(gName);
	nPoints = 0;
	// Some cosmetics at creation time:
	if(i==0) {
	  gEextr[i][ic]->SetLineWidth(3);
	  gEextr[i][ic]->SetLineColor(PGlobals::fieldLine);
	  gEextr[i][ic]->SetMarkerStyle(20);
	  gEextr[i][ic]->SetMarkerSize(0.4);
	  gEextr[i][ic]->SetMarkerColor(PGlobals::fieldLine);	
	  gEextr[i][ic]->GetYaxis()->SetTitle("E_{z}/E_{0}");
	  gEextr[i][ic]->GetXaxis()->SetTitle("k_{p}z");
	} else if(i==1) {
	  gEextr[i][ic]->SetLineWidth(1);
	  gEextr[i][ic]->SetLineColor(kGray+2);
	  gEextr[i][ic]->SetMarkerStyle(20);
	  gEextr[i][ic]->SetMarkerSize(0.4);
	  gEextr[i][ic]->SetMarkerColor(kGray+2);	
	  gEextr[i][ic]->GetYaxis()->SetTitle("E_{y}/E_{0}");
	  gEextr[i][ic]->GetXaxis()->SetTitle("k_{p}z");	  
	}
      } else {
	nPoints = gEextr[i][ic]->GetN(); 
      }  
      
      gEextr[i][ic]->Set(nPoints+1);
      gEextr[i][ic]->SetPoint(nPoints,Time,Eextr[i][ic]);
      gEextr[i][ic]->Write(gName,TObject::kOverwrite);
    }
  }

  
  ifile->Close();
  
}

