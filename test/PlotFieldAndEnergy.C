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

void PlotFieldAndEnergy( const TString &sim, Int_t time, Float_t Emin, Float_t Emax, Int_t NYbins=1, const TString &opt="") { 
  
  PlasmaGlob::Initialize();
    
  gStyle->SetPadLeftMargin(0.10);  // Margin left axis  
  gStyle->SetPadRightMargin(0.20); // Margin right axis 
  
  // Load PData
  PData *pData = PData::Get(sim.Data());
  pData->LoadFileNames(time);
  if(!pData->IsInit()) return;
  
  Bool_t ThreeD = kFALSE;
  if(sim.Contains("3D")) ThreeD = kTRUE; 
  
  // Histos for bunch specie
  TH2F *hDen2D  = NULL;   // Charge density
  TH1F *hEnergy = NULL;   // Energy spectrum
  TH2F *hEvsX1  = NULL;   // Energy vs x1
  Int_t Nspecies = pData->NSpecies();
  for(Int_t i=0;i<Nspecies;i++) {
    
    if(!pData->GetChargeFileName(i)) continue;
      
    // Get 2D charge density histogram of the bunch.
    // We use it to calculate the rms of the charge distribution
    // in x2 as a function of x1.
    if(pData->GetSpeciesName(i).find("beam")==string::npos) 
      continue;
     
    if(!ThreeD) 
      hDen2D = pData->GetCharge(i);
    else
      hDen2D = pData->GetCharge2DSliceZY(i);
    
    hDen2D->SetName("hDen2D");
    hDen2D->GetXaxis()->CenterTitle();
    hDen2D->GetYaxis()->CenterTitle();
    hDen2D->GetXaxis()->SetTitle("z [c/#omega_{p}]");
    hDen2D->GetYaxis()->SetTitle("y [c/#omega_{p}]");

    Float_t x1Min  = hDen2D->GetXaxis()->GetXmin();
    Float_t x1Max  = hDen2D->GetXaxis()->GetXmax();
    Int_t NbinX1 = hDen2D->GetNbinsX();
    
    // Get spectrum
    char hName[24];
    sprintf(hName,"hEnergy");
    hEnergy = (TH1F*) gROOT->FindObject(hName);
    if(hEnergy) { delete hEnergy; hEnergy = NULL; }
    hEnergy = new TH1F(hName,"",200,Emin,Emax);
    
    pData->GetH1Raw(pData->GetRawFileName(i)->c_str(),"ene",hEnergy);
    
    hEnergy->GetXaxis()->SetTitle("E [MeV]");
    hEnergy->GetYaxis()->SetTitle("dQ/dE [a.u.]");
    PlasmaGlob::SetH1Style(hEnergy,i);
    
    // Get Spectrum vs x1
    sprintf(hName,"hEvsX1");
    hEvsX1 = (TH2F*) gROOT->FindObject(hName);
    if(hEvsX1) { delete hEvsX1; hEvsX1 = NULL; }
    hEvsX1 = new TH2F(hName,"",NbinX1,x1Min,x1Max,200,Emin,Emax);
    
    pData->GetH2Raw(pData->GetRawFileName(i)->c_str(),"x1","ene",hEvsX1);
    hEvsX1->GetXaxis()->SetTitle("z [c/#omega_{p}]");
    hEvsX1->GetYaxis()->SetTitle("#LTE#GT [MeV]");
    PlasmaGlob::SetH1Style(hEvsX1,i);
    
    if(opt.Contains("comov")) {
      Int_t NbinsX = hDen2D->GetNbinsX();
      Float_t xMin = hDen2D->GetXaxis()->GetXmin()-hDen2D->GetXaxis()->GetXmax();
      Float_t xMax = hDen2D->GetXaxis()->GetXmax()-hDen2D->GetXaxis()->GetXmax();
      Int_t NbinsY = hDen2D->GetNbinsY();
      Float_t yMin = hDen2D->GetYaxis()->GetXmin();
      Float_t yMax = hDen2D->GetYaxis()->GetXmax();
      hDen2D->SetBins(NbinsX,xMin,xMax,NbinsY,yMin,yMax);
      
      NbinsX = hEvsX1->GetNbinsX();
      xMin = hEvsX1->GetXaxis()->GetXmin()-hEvsX1->GetXaxis()->GetXmax();
      xMax = hEvsX1->GetXaxis()->GetXmax()-hEvsX1->GetXaxis()->GetXmax();
      NbinsY = hEvsX1->GetNbinsY();
      yMin = hEvsX1->GetYaxis()->GetXmin();
      yMax = hEvsX1->GetYaxis()->GetXmax();
      hEvsX1->SetBins(NbinsX,xMin,xMax,NbinsY,yMin,yMax);    
    }

    // Chaning to user units: 
    if(opt.Contains("units") && pData->GetPlasmaDensity()) {
      Int_t NbinsX = hDen2D->GetNbinsX();
      Float_t xMin = 1e3 * pData->GetPlasmaSkinDepth() * hDen2D->GetXaxis()->GetXmin();
      Float_t xMax = 1e3 * pData->GetPlasmaSkinDepth() * hDen2D->GetXaxis()->GetXmax();
      Int_t NbinsY = hDen2D->GetNbinsY();
      Float_t yMin = 1e3 * pData->GetPlasmaSkinDepth() * hDen2D->GetYaxis()->GetXmin();
      Float_t yMax = 1e3 * pData->GetPlasmaSkinDepth() * hDen2D->GetYaxis()->GetXmax();
      hDen2D->SetBins(NbinsX,xMin,xMax,NbinsY,yMin,yMax);	
      
      hDen2D->GetYaxis()->SetTitle("y [mm]");	
      if(opt.Contains("comov"))
	hDen2D->GetXaxis()->SetTitle("#zeta [mm]");
      else
	hDen2D->GetXaxis()->SetTitle("z [mm]");
    
      NbinsX = hEvsX1->GetNbinsX();
      xMin = 1e3 * pData->GetPlasmaSkinDepth() * hEvsX1->GetXaxis()->GetXmin();
      xMax = 1e3 * pData->GetPlasmaSkinDepth() * hEvsX1->GetXaxis()->GetXmax();
      NbinsY = hEvsX1->GetNbinsY();
      yMin = hEvsX1->GetYaxis()->GetXmin();
      yMax = hEvsX1->GetYaxis()->GetXmax();
      hEvsX1->SetBins(NbinsX,xMin,xMax,NbinsY,yMin,yMax);  
      
      if(opt.Contains("comov"))
	hEvsX1->GetXaxis()->SetTitle("#zeta [mm]");
      else
	hEvsX1->GetXaxis()->SetTitle("z [mm]");
      
    }

  }
   
  TProfile *hDen2Dprof = NULL;
  TH1F *hRms = NULL;
  if(hDen2D) {
    TString pname = hDen2D->GetName();
    pname += "_pfx";
    hDen2Dprof =  (TProfile*) gROOT->FindObject(pname.Data());
    if(hDen2Dprof) delete hDen2Dprof;
    hDen2Dprof = hDen2D->ProfileX("_pfx",1,-1,"s");

    hRms = (TH1F*) gROOT->FindObject("hRms");
    if(hRms) delete hRms;

    hRms = new TH1F("hRms","",hDen2D->GetNbinsX(),hDen2D->GetXaxis()->GetXmin(),
		    hDen2D->GetXaxis()->GetXmax());

    for(Int_t j=0;j<hRms->GetNbinsX();j++) {
      if(opt.Contains("units") && pData->GetBeamRmsR())
	hRms->SetBinContent(j,hDen2Dprof->GetBinError(j)/(pData->GetBeamRmsR() * 1e-3) ); 
      else if(pData->GetBeamRmsZ())
	hRms->SetBinContent(j,hDen2Dprof->GetBinError(j) / pData->GetBeamRmsR() 
			    * pData->GetPlasmaSkinDepth() * 1e6 ); 
      
    }
    
    hRms->GetXaxis()->SetTitle("z [c/#omega_{p}]");
    if(opt.Contains("comov"))
      hRms->GetXaxis()->SetTitle("#zeta [c/#omega_{p}]");
    
    if(opt.Contains("units") && pData->GetBeamRmsR()) {
      hRms->GetXaxis()->SetTitle("z [mm]");
      if(opt.Contains("comov"))
	hRms->GetXaxis()->SetTitle("#zeta [mm]");

    }

    hRms->GetYaxis()->SetTitle("r_{b}/r_{0}");
    
  }
  
  
  TProfile *hEvsX1prof = NULL;
  TH1F *hEmean = NULL;
  if(hEvsX1) {
    TString pname = hEvsX1->GetName();
    pname += "_pfx";
    hEvsX1prof =  (TProfile*) gROOT->FindObject(pname.Data());
    if(hEvsX1prof) delete hEvsX1prof;
    hEvsX1prof = hEvsX1->ProfileX("_pfx",1,-1,"s");

    hEmean = (TH1F*) gROOT->FindObject("hEmean");
    if(hEmean) delete hEmean;

    hEmean = new TH1F("hEmean","",hEvsX1->GetNbinsX(),hEvsX1->GetXaxis()->GetXmin(),
		      hEvsX1->GetXaxis()->GetXmax());

    for(Int_t j=0;j<hEmean->GetNbinsX();j++) {
      hEmean->SetBinContent(j,hEvsX1prof->GetBinContent(j));
      // hEmean->SetBinError(j,hEvsX1prof->GetBinError(j));   
      hEmean->SetBinError(j,0);     
    }
    
    hEmean->GetXaxis()->SetTitle("z [c/#omega_{p}]");
    if(opt.Contains("comov"))
      hEmean->GetXaxis()->SetTitle("#zeta [c/#omega_{p}]");
    
    if(opt.Contains("units") && pData->GetBeamRmsR()) {
      hEmean->GetXaxis()->SetTitle("z [mm]");
      if(opt.Contains("comov"))
	hEmean->GetXaxis()->SetTitle("#zeta [mm]");
    }    

    hEmean->GetYaxis()->SetTitle("#LTE#GT [MeV]");
    
  }
  
  
  // Get electric fields
  const Int_t Nfields = 2;
  TH1F **hE1D = new TH1F*[Nfields];
  for(Int_t i=0;i<Nfields;i++) {
    hE1D[i] = NULL;
    if(!pData->GetEfieldFileName(i))
      continue;
    
    char nam[3]; sprintf(nam,"e%i",i+1);
    if(!ThreeD) {
      if(i==0) {
	hE1D[i] = pData->GetH1SliceZ(pData->GetEfieldFileName(i)->c_str(),nam,-1,NYbins);
      }
      else { 
	hE1D[i] = pData->GetH1SliceZ(pData->GetEfieldFileName(i)->c_str(),nam,-NYbins,NYbins);
      }
    } else {
      if(i==0) {
	hE1D[i] = pData->GetH1SliceZ3D(pData->GetEfieldFileName(i)->c_str(),nam,-1,NYbins,-1,NYbins);
	}
      else { 
	hE1D[i] = pData->GetH1SliceZ3D(pData->GetEfieldFileName(i)->c_str(),nam,-NYbins,NYbins,-NYbins,NYbins);
      }
      
    }
    
    char hName[24];
    sprintf(hName,"hE_%i",i);
    hE1D[i]->SetName(hName);
    hE1D[i]->GetXaxis()->SetTitle("z [c/#omega_{p}]");
     if(i==0)
      hE1D[i]->GetYaxis()->SetTitle("E_{z} [E_{0}]");
    else if(i==1)
      hE1D[i]->GetYaxis()->SetTitle("E_{y} [E_{0}}]");
    else if(i==2)
      hE1D[i]->GetYaxis()->SetTitle("E_{x} [E_{0}}]");
    
    // Change to co-moving coordinate
    if(opt.Contains("comov")) {
      Int_t NbinsX = hE1D[i]->GetNbinsX();
      Float_t xMin = hE1D[i]->GetXaxis()->GetXmin()-hE1D[i]->GetXaxis()->GetXmax();
      Float_t xMax = hE1D[i]->GetXaxis()->GetXmax()-hE1D[i]->GetXaxis()->GetXmax();
      hE1D[i]->SetBins(NbinsX,xMin,xMax);
      
      hE1D[i]->GetXaxis()->SetTitle("#zeta [c/#omega_{p}]");
    }
    
    // Chaning to user units: 
    if(opt.Contains("units") && pData->GetPlasmaDensity()) {
      Int_t NbinsX = hE1D[i]->GetNbinsX();
      Float_t xMin = 1e3 * pData->GetPlasmaSkinDepth() * hE1D[i]->GetXaxis()->GetXmin();
      Float_t xMax = 1e3 * pData->GetPlasmaSkinDepth() * hE1D[i]->GetXaxis()->GetXmax();
      hE1D[i]->SetBins(NbinsX,xMin,xMax);
      
      for(Int_t j=0;j<hE1D[i]->GetNbinsX();j++) {
	hE1D[i]->SetBinContent(j, 1e-9 * pData->GetPlasmaE0() * hE1D[i]->GetBinContent(j));
      }
      
      hE1D[i]->GetYaxis()->SetTitle("E [GV/m]");
      
      if(opt.Contains("comov"))
	hE1D[i]->GetXaxis()->SetTitle("#zeta [mm]");
      else
	hE1D[i]->GetXaxis()->SetTitle("z [mm]");
    }
    
  }

  
  // Calculate wave positions:
  // ----------------------------------------------------------------

  // Local Minimum search window
  Float_t xs1min = -99;
  Float_t xs1max = -99;
  
  // Retrieve the previous time TGraph if any;
  // Open TGraph
  TString filename = Form("./%s/Plots/FieldAndEnergy/FieldAndEnergy-%s.root",sim.Data(),sim.Data());
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

  // Get the Graph with the x1 positions of the maximum E_1
  TGraph *gEfieldMaxPos0 =  (TGraph*) ifile->Get("gEfieldMaxPos_0");
  if(gEfieldMaxPos0) {
    Double_t *y = gEfieldMaxPos0->GetY();
    // Setup the searching windows to +/- pi/2 respect the last found minimum.
    if(opt.Contains("units") && pData->GetPlasmaDensity()) {
      xs1min = y[gEfieldMaxPos0->GetN()-1] - (0.5*TMath::Pi() * 1e3 * pData->GetPlasmaSkinDepth());
      xs1max = y[gEfieldMaxPos0->GetN()-1] + (0.5*TMath::Pi() * 1e3 * pData->GetPlasmaSkinDepth());
    } else {
      xs1min = y[gEfieldMaxPos0->GetN()-1] - 0.5*TMath::Pi();
      xs1max = y[gEfieldMaxPos0->GetN()-1] + 0.5*TMath::Pi();
    }
  }

  Float_t *EfieldMaxPos = new Float_t[Nfields];
  Float_t *EfieldMaxValue = new Float_t[Nfields];
  for(Int_t i=0;i<Nfields;i++) {
    if(!hE1D[i]) continue;
    EfieldMaxPos[i] = EfieldMaxValue[i] = -999.;    
    Float_t xCenter = (hE1D[i]->GetXaxis()->GetXmin()+hE1D[i]->GetXaxis()->GetXmax())/2.;
    TH1F *htemp = (TH1F*) hE1D[i]->Clone("htemp");
    
    if(i==0 && (xs1min==-99) && (xs1max==-99)  ) { // Ez field: Find global minimum.
      htemp->GetXaxis()->SetRangeUser(hE1D[i]->GetXaxis()->GetXmin(),xCenter);
    } else {
      htemp->GetXaxis()->SetRangeUser(xs1min,xs1max);
    }
    
    //htemp->Smooth(1,"R");
    Int_t binMax = htemp->GetMinimumBin();
    EfieldMaxPos[i] = (Float_t) htemp->GetBinCenter(binMax);
    EfieldMaxValue[i] = (Float_t) htemp->GetBinContent(binMax);
    delete htemp;

    if(i==0) { // Ez field: Find global minimum.
      if(opt.Contains("units") && pData->GetPlasmaDensity()) {
	xs1min = EfieldMaxPos[i] - (0.5*TMath::Pi()  * 1e3 * pData->GetPlasmaSkinDepth());
	xs1max = EfieldMaxPos[i] + (0.5*TMath::Pi()  * 1e3 * pData->GetPlasmaSkinDepth());
      } else {
	xs1min = EfieldMaxPos[i] - 0.5*TMath::Pi();
	xs1max = EfieldMaxPos[i] + 0.5*TMath::Pi();
      }
    }
    
  }
  
  // Bunch RMS: Find local extreme
  // Trick for nicer plotting
  Float_t bin1 = -1;
  Float_t bin2 = -1;
  for(Int_t i=0;i<hRms->GetNbinsX();i++) 
    if(bin1==-1 && hRms->GetBinContent(i+1)>0) bin1 = hRms->GetBinCenter(i+1);  
  for(Int_t i=hRms->GetNbinsX()-1;i>=0;i--) 
    if(bin2==-1 && hRms->GetBinContent(i+1)>0) bin2 = hRms->GetBinCenter(i+1);
  hRms->GetXaxis()->SetRangeUser(bin1,bin2);
  Float_t hRmsMax = -999.;
  Float_t hRmsMaxValue = -999.;
  if(hRms) {
    Float_t xCenter = (hRms->GetXaxis()->GetXmin()+hRms->GetXaxis()->GetXmax())/2.;
    TH1F *htemp = (TH1F*) hRms->Clone("htemp");
    htemp->GetXaxis()->SetRangeUser(bin1,(bin1+bin2)/2.);
    //htemp->Smooth(1,"R");
    Int_t binMax = htemp->GetMaximumBin();
    hRmsMax = (Float_t) htemp->GetBinCenter(binMax);
    hRmsMaxValue = (Float_t) htemp->GetBinContent(binMax);
    delete htemp;
  }
  
  // Bunch Energy gain: Find local extreme
  // Trick for nicer plotting
  bin1 = -1;
  bin2 = -1;
  for(Int_t i=0;i<hEmean->GetNbinsX();i++) 
    if(bin1==-1 && hEmean->GetBinContent(i+1)>0) bin1 = hEmean->GetBinCenter(i+1);  
  for(Int_t i=hEmean->GetNbinsX()-1;i>=0;i--) 
    if(bin2==-1 && hEmean->GetBinContent(i+1)>0) bin2 = hEmean->GetBinCenter(i+1);
  hEmean->GetXaxis()->SetRangeUser(bin1,bin2);

  Float_t EmeanMaxPos = -999.;
  Float_t EmeanMaxValue = -999.;
  // Get the Graph with the x1 positions of the maximum Energy
  TGraph *gEmeanMaxPos0 =  (TGraph*) ifile->Get("gEmeanMaxPos");
  Float_t xCenter = 0;
  if(gEmeanMaxPos0) {
    Double_t *y = gEmeanMaxPos0->GetY();
    xCenter = y[gEmeanMaxPos0->GetN()-1];
  }
  if(hEmean) {
    if(xCenter==0) xCenter = (bin1+bin2)/2.;
    TH1F *htemp = (TH1F*) hEmean->Clone("htemp");
    Float_t x1,x2;
    if(opt.Contains("units") && pData->GetPlasmaDensity()) {
      x1 = xCenter - (0.5*TMath::Pi()  * 1e3 * pData->GetPlasmaSkinDepth());
      x2 = xCenter + (0.5*TMath::Pi()  * 1e3 * pData->GetPlasmaSkinDepth());
    } else {
      x1 = xCenter - pData->GetPlasmaSkinDepth();
      x2 = xCenter + pData->GetPlasmaSkinDepth();
    }
    htemp->GetXaxis()->SetRangeUser(x1,x2);
    //htemp->Smooth(1,"R");
    Int_t binMax = htemp->GetMaximumBin();
    EmeanMaxPos = (Float_t) htemp->GetBinCenter(binMax);
    EmeanMaxValue = (Float_t) htemp->GetBinContent(binMax);
    delete htemp;
  }
  
  // Calculate weighted energy gain:
  Float_t Egain = 0.;
  Float_t Eloss = 0.;
  if(hEnergy) {
    for(Int_t i=0;i<hEnergy->GetNbinsX();i++) {
      if(hEnergy->GetBinCenter(i+1) > pData->GetBeamEnergy()) {
	Egain += hEnergy->GetBinCenter(i+1)*hEnergy->GetBinContent(i+1);
      } else {
	Eloss += hEnergy->GetBinCenter(i+1)*hEnergy->GetBinContent(i+1);
      }
    } 
  }
  Float_t Asymm = (Egain-Eloss)/(Egain+Eloss);
  
  // Tunning the Histograms
  // ---------------------
  
  Float_t Time = pData->GetRealTime();
  
  // Plotting
  // -----------------------------------------------

  // Output file
  TString fOutName = Form("./%s/Plots/FieldAndEnergy/FieldAndEnergy",sim.Data());
  fOutName += Form("-%s_%i",sim.Data(),time);

  // Canvas setup
  TCanvas *C = new TCanvas("C","Charge density and Electric field",850,1000);
  C->Divide(1,2);

  // Draw objects
  TPaveText *textTime = new TPaveText(0.63,0.87,0.78,0.92,"NDC");
  PlasmaGlob::SetPaveTextStyle(textTime); 
  char ctext[128];
  if(opt.Contains("units") && pData->GetPlasmaDensity()) 
    sprintf(ctext,"L = %5.1f mm", 1e3 * pData->GetPlasmaSkinDepth() * Time);
  else
    sprintf(ctext,"T = %5.1f 1/#omega_{p}",Time);
  textTime->AddText(ctext);
  
  // Colors
  Int_t plasmaC = kGray+1;
  Int_t beamC   = kAzure-5;
  Int_t fieldC  = kGray+1; // PlasmaGlob::fieldLine;
  Int_t fieldCb = kGray+1; // PlasmaGlob::fieldLine;
  Int_t eneC = kOrange+10;
  Int_t fieldMaxC = kGray+2;

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
  hE1D[0]->SetLineColor(fieldC);
  hE1D[0]->SetMarkerStyle(20);

  hE1D[1]->GetYaxis()->CenterTitle();
  hE1D[1]->GetXaxis()->CenterTitle();
  hE1D[1]->SetLineStyle(3);
  hE1D[1]->SetLineColor(fieldCb);
  hE1D[1]->SetMarkerStyle(24);

  Float_t factor = 1.5;
  Float_t minimum = factor * hE1D[0]->GetMinimum();
  Float_t maximum = factor * hE1D[0]->GetMaximum();
    
  if(hE1D[1]->GetMaximum() > hE1D[0]->GetMaximum()) {
    maximum = factor * hE1D[1]->GetMaximum();
  } 

  if(hE1D[1]->GetMinimum() < hE1D[0]->GetMinimum()) {
    minimum = factor * hE1D[1]->GetMinimum();
  } 

  if( maximum >= TMath::Abs(minimum)) minimum = -maximum;
  else maximum = - minimum;
  
  hE1D[0]->GetYaxis()->SetRangeUser(minimum,maximum);
  hE1D[0]->Draw("C");
  hE1D[1]->Draw("C same");
  
  C->Update();
  
  TLine *line0 = new TLine(hE1D[0]->GetXaxis()->GetXmin(),
			   (gPad->GetUymin()+gPad->GetUymax())/2.,
			   hE1D[0]->GetXaxis()->GetXmax(),
			   (gPad->GetUymin()+gPad->GetUymax())/2.);
  line0->SetLineColor(kGray+1);
  line0->SetLineStyle(2);
  line0->Draw();
  
  TMarker *markEMax0 = new TMarker(EfieldMaxPos[0],EfieldMaxValue[0], 20);
  markEMax0->SetMarkerColor(fieldC);
  markEMax0->Draw();

  TMarker *markEMax1 = new TMarker(EfieldMaxPos[1],EfieldMaxValue[1], 24);
  markEMax1->SetMarkerColor(fieldC);
  markEMax1->Draw();

  Float_t rightmin = 0;
  Float_t rightmax = 2.5 * hRms->GetMaximum();
  Float_t slope = (gPad->GetUymax() - gPad->GetUymin())/(rightmax-rightmin);
  
  for(Int_t i=0;i<hRms->GetNbinsX();i++) {
    hRms->SetBinContent(i+1,(hRms->GetBinContent(i+1)-rightmin)*slope  + gPad->GetUymin());
  }
  
  hRms->SetLineStyle(1);
  hRms->SetLineColor(beamC);
  hRms->Draw("same C");

  TGaxis *axisRMS = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),gPad->GetUxmax(),
  			    gPad->GetUymax(),rightmin,rightmax,505,"+L");
  
  axisRMS->SetLineWidth(1);
  axisRMS->SetLineColor(beamC);
  axisRMS->SetLabelColor(beamC);

  axisRMS->SetTitle("r_{b}/r_{0}");
  axisRMS->CenterTitle();
  axisRMS->SetTitleColor(beamC);
  axisRMS->SetTitleOffset(0.8);
  axisRMS->Draw();

  TMarker *markRmsMax = new TMarker(hRmsMax,(hRmsMaxValue-rightmin)*slope  + gPad->GetUymin(), 20);
  markRmsMax->SetMarkerColor(beamC);
  // markRmsMax->Draw();

  // Line around 1
  Float_t value = (1-rightmin)*slope + gPad->GetUymin();
  TLine *line1 = new TLine(hE1D[0]->GetXaxis()->GetXmin(),value,
			   hE1D[0]->GetXaxis()->GetXmax(),value);
  line1->SetLineColor(beamC);
  line1->SetLineStyle(2);
  line1->Draw();
    
  // rightmax = 0.8*(hEmean->GetMaximum() - 25.0) + hEmean->GetMaximum();
  // rightmin = 25.0 - (rightmax-25.0);
  rightmin = Emin;
  rightmax = Emax;
  slope = (gPad->GetUymax() - gPad->GetUymin())/(rightmax-rightmin);
  
  for(Int_t i=0;i<hEmean->GetNbinsX();i++) {
    hEmean->SetBinContent(i+1,(hEmean->GetBinContent(i+1)-rightmin)*slope  + gPad->GetUymin());
  }
  
  hEmean->SetLineStyle(1);
  hEmean->SetLineColor(eneC);
  hEmean->Draw("same C");
  
  //draw an axis on the right side
  Float_t rightmargin = 0.08;
  Float_t ux = gPad->PixeltoX(gPad->UtoPixel(1-rightmargin));
  TGaxis *axisEmean = new TGaxis(ux,gPad->GetUymin(),
				 ux,gPad->GetUymax(),rightmin,rightmax,505,"+L");
  
  axisEmean->SetLineWidth(1);
  axisEmean->SetLineColor(eneC);
  axisEmean->SetLabelColor(eneC);
  axisEmean->SetTitleColor(eneC);
  axisEmean->SetTitle("#LTE#GT [MeV]");
  axisEmean->CenterTitle();
  axisEmean->SetTitleOffset(0.8);
  axisEmean->Draw();
  
  TMarker *markEmeanMax = new TMarker(EmeanMaxPos,(EmeanMaxValue-rightmin)*slope + gPad->GetUymin(), 20);
  markEmeanMax->SetMarkerColor(eneC);
  markEmeanMax->Draw();
  
  textTime->Draw();

  // Define the TGraphs
  Int_t nPoints = 0;  
  TGraph **gEfieldMaxPos = new TGraph*[Nfields];
  TGraph **gEfieldMaxValue = new TGraph*[Nfields];
  for(Int_t i=0;i<Nfields;i++) {
    if(!hE1D[i]) continue;

    char gName[24];
    sprintf(gName,"gEfieldMaxPos_%i",i);
    gEfieldMaxPos[i] = (TGraph*) ifile->Get(gName);
    if(gEfieldMaxPos[i]==NULL) {
      gEfieldMaxPos[i] = new TGraph();
      gEfieldMaxPos[i]->SetName(gName);
    } else {
      nPoints = gEfieldMaxPos[i]->GetN(); 
    }  
    gEfieldMaxPos[i]->Set(nPoints+1);
    if(opt.Contains("units") && pData->GetPlasmaDensity()) 
      gEfieldMaxPos[i]->SetPoint(nPoints, 1e3 * pData->GetPlasmaSkinDepth() * Time,EfieldMaxPos[i]);
    else
      gEfieldMaxPos[i]->SetPoint(nPoints,Time,EfieldMaxPos[i]);
    
    if(opt.Contains("units") && pData->GetPlasmaDensity()) {
      gEfieldMaxPos[i]->GetYaxis()->SetTitle("#zeta_{min} [mm]");
      gEfieldMaxPos[i]->GetXaxis()->SetTitle("L [mm]");
    } else {
      gEfieldMaxPos[i]->GetYaxis()->SetTitle("#zeta_{min} [c/#omega_{p}]");
      gEfieldMaxPos[i]->GetXaxis()->SetTitle("T [c/#omega_{p}]");
    }

    gEfieldMaxPos[i]->Write(gEfieldMaxPos[i]->GetName(),TObject::kOverwrite);

    sprintf(gName,"gEfieldMaxValue_%i",i);
    gEfieldMaxValue[i] = (TGraph*) ifile->Get(gName);
    if(gEfieldMaxValue[i]==NULL) {
      gEfieldMaxValue[i] = new TGraph();  
      gEfieldMaxValue[i]->SetName(gName);
    } else {
      nPoints = gEfieldMaxValue[i]->GetN(); 
    }  
    gEfieldMaxValue[i]->Set(nPoints+1);
    if(opt.Contains("units") && pData->GetPlasmaDensity()) 
      gEfieldMaxValue[i]->SetPoint(nPoints, 1e3 * pData->GetPlasmaSkinDepth() * Time,EfieldMaxValue[i]);
    else
      gEfieldMaxValue[i]->SetPoint(nPoints,Time,EfieldMaxValue[i]);
    
    if(opt.Contains("units") && pData->GetPlasmaDensity()) {
      gEfieldMaxValue[i]->GetYaxis()->SetTitle("E_{min} [GV/m]");
      gEfieldMaxValue[i]->GetXaxis()->SetTitle("L [mm]");
    } else {
      gEfieldMaxValue[i]->GetYaxis()->SetTitle("E_{min} [E_{0}]");
      gEfieldMaxValue[i]->GetXaxis()->SetTitle("T [c/#omega_{p}]");
    }
    
    gEfieldMaxValue[i]->Write(gEfieldMaxValue[i]->GetName(),TObject::kOverwrite);    
  }

  TGraph *gRmsMax = NULL;
  if(hRms) {
    char gName[24];
    sprintf(gName,"gRmsMax");
    gRmsMax = (TGraph*) ifile->Get(gName);
    
    if(gRmsMax==NULL) {
      gRmsMax = new TGraph();
      gRmsMax->SetName(gName);
    } else {
      nPoints = gRmsMax->GetN(); 
    }  
    gRmsMax->Set(nPoints+1);
     if(opt.Contains("units") && pData->GetPlasmaDensity()) 
      gRmsMax->SetPoint(nPoints, 1e3 * pData->GetPlasmaSkinDepth() * Time,hRmsMax);
    else
      gRmsMax->SetPoint(nPoints,Time,hRmsMax);

    if(opt.Contains("units") && pData->GetPlasmaDensity()) {
      gRmsMax->GetYaxis()->SetTitle("#zeta_{min} [mm]");
      gRmsMax->GetXaxis()->SetTitle("L [mm]");
    } else {
      gRmsMax->GetYaxis()->SetTitle("#zeta_{min} [c/#omega_{p}]");
      gRmsMax->GetXaxis()->SetTitle("T [c/#omega_{p}]");
    }


    gRmsMax->Write(gRmsMax->GetName(),TObject::kOverwrite);
    
  }

  TGraph *gEmeanMaxPos = NULL;
  TGraph *gEmeanMaxValue = NULL;
  if(hEmean) {
    char gName[24];
    sprintf(gName,"gEmeanMaxPos");
    gEmeanMaxPos = (TGraph*) ifile->Get(gName);
    if(gEmeanMaxPos==NULL) {
      gEmeanMaxPos = new TGraph();  
      gEmeanMaxPos->SetName(gName);
    } else {
      nPoints = gEmeanMaxPos->GetN(); 
    }  
    gEmeanMaxPos->Set(nPoints+1);
    if(opt.Contains("units") && pData->GetPlasmaDensity()) 
      gEmeanMaxPos->SetPoint(nPoints, 1e3 * pData->GetPlasmaSkinDepth() * Time,EmeanMaxPos);
    else
      gEmeanMaxPos->SetPoint(nPoints,Time,EmeanMaxPos);

    if(opt.Contains("units") && pData->GetPlasmaDensity()) {
      gEmeanMaxPos->GetYaxis()->SetTitle("#zeta_{min} [mm]");
      gEmeanMaxPos->GetXaxis()->SetTitle("L [mm]");
    } else {
      gEmeanMaxPos->GetYaxis()->SetTitle("#zeta_{min} [c/#omega_{p}]");
      gEmeanMaxPos->GetXaxis()->SetTitle("T [c/#omega_{p}]");
    }
    gEmeanMaxPos->Write(gEmeanMaxPos->GetName(),TObject::kOverwrite);

    sprintf(gName,"gEmeanMaxValue");
    gEmeanMaxValue = (TGraph*) ifile->Get(gName);
    if(gEmeanMaxValue==NULL) {
      gEmeanMaxValue = new TGraph();  
      gEmeanMaxValue->SetName(gName);
    } else {
      nPoints = gEmeanMaxValue->GetN(); 
    }  
    gEmeanMaxValue->Set(nPoints+1);
    if(opt.Contains("units") && pData->GetPlasmaDensity()) 
      gEmeanMaxValue->SetPoint(nPoints, 1e3 * pData->GetPlasmaSkinDepth() * Time,EmeanMaxValue);
    else
      gEmeanMaxValue->SetPoint(nPoints,Time,EmeanMaxValue);

    if(opt.Contains("units") && pData->GetPlasmaDensity()) {
      gEmeanMaxValue->GetYaxis()->SetTitle("#LTE#GT_{max} [MeV]");
      gEmeanMaxValue->GetXaxis()->SetTitle("L [mm]");
    } else {
      gEmeanMaxValue->GetYaxis()->SetTitle("#LTE#GT_{max} [MeV]");
      gEmeanMaxValue->GetXaxis()->SetTitle("T [c/#omega_{p}]");
    }
 
    gEmeanMaxValue->Write(gEmeanMaxValue->GetName(),TObject::kOverwrite);
    
  }

  TGraph *gEgain = NULL;
  if(hEnergy) {
    char gName[24];
    sprintf(gName,"gEgain");
    gEgain = (TGraph*) ifile->Get(gName);
    if(gEgain==NULL) {
      gEgain = new TGraph();  
      gEgain->SetName(gName);
    } else {
      nPoints = gEgain->GetN(); 
    }  
    gEgain->Set(nPoints+1);

    if(opt.Contains("units") && pData->GetPlasmaDensity()) 
      gEgain->SetPoint(nPoints, 1e3 * pData->GetPlasmaSkinDepth() * Time,Asymm);
    else
      gEgain->SetPoint(nPoints,Time,Asymm);

    if(opt.Contains("units") && pData->GetPlasmaDensity()) {
      gEgain->GetYaxis()->SetTitle("E^{gain} [MeV]");
      gEgain->GetXaxis()->SetTitle("L [mm]");
    } else {
      gEgain->GetYaxis()->SetTitle("E^{gain} [MeV]");
      gEgain->GetXaxis()->SetTitle("T [c/#omega_{p}]");
    }
   
    gEgain->Write(gEgain->GetName(),TObject::kOverwrite);
  }

  C->cd(2);
  gPad->SetGridy(1);
  gPad->SetGridx(0);
  gPad->SetFrameLineWidth(2);  

  xCenter = (hE1D[0]->GetXaxis()->GetXmin()+hE1D[0]->GetXaxis()->GetXmax())/2.;
  Float_t minPhase = 99.;
  Float_t maxPhase = -99.;
  Float_t minField = 99.;
  Float_t maxField = -99.;
  Float_t minEmean = 99.;
  Float_t maxEmean = -99.;

  Double_t *yEfieldMaxPos_0   = gEfieldMaxPos[0]->GetY();
  Double_t *yEfieldMaxPos_1   = gEfieldMaxPos[1]->GetY();
  Double_t *yEfieldMaxValue_0 = gEfieldMaxValue[0]->GetY();
  Double_t *yEmeanMaxValue    = gEmeanMaxValue->GetY();
  
  for(Int_t j=0;j<gEfieldMaxPos0->GetN();j++) {
    if(yEfieldMaxPos_0[j]>maxPhase)
      maxPhase = yEfieldMaxPos_0[j];
    if(yEfieldMaxPos_0[j]<minPhase)
      minPhase = yEfieldMaxPos_0[j];

    if(yEfieldMaxPos_1[j]>maxPhase)
      maxPhase = yEfieldMaxPos_1[j];
    if(yEfieldMaxPos_1[j]<minPhase)
      minPhase = yEfieldMaxPos_1[j];
 
    if(yEfieldMaxValue_0[j]>maxField)
      maxField = yEfieldMaxValue_0[j];
    if(yEfieldMaxValue_0[j]<minField)
      minField = yEfieldMaxValue_0[j];

    if(yEmeanMaxValue[j]>maxEmean)
      maxEmean = yEmeanMaxValue[j];
    if(yEmeanMaxValue[j]<minEmean)
      minEmean = yEmeanMaxValue[j];
  }
  
  Float_t margin = (maxPhase - minPhase)/10;
  gEfieldMaxPos[0]->GetYaxis()->SetRangeUser(minPhase-margin,maxPhase+margin);
  gEfieldMaxPos[0]->GetYaxis()->CenterTitle();
  gEfieldMaxPos[0]->GetXaxis()->CenterTitle();
  gEfieldMaxPos[0]->SetLineColor(fieldCb);
  gEfieldMaxPos[0]->SetMarkerColor(fieldCb);
  gEfieldMaxPos[0]->SetLineWidth(1);
  gEfieldMaxPos[0]->SetMarkerStyle(20);
  gEfieldMaxPos[0]->Draw("APC");

  gEfieldMaxPos[1]->SetLineStyle(3);
  gEfieldMaxPos[1]->SetLineColor(fieldCb);
  gEfieldMaxPos[1]->SetMarkerColor(fieldCb);
  gEfieldMaxPos[1]->SetLineWidth(1);
  gEfieldMaxPos[1]->SetMarkerStyle(24);
  gEfieldMaxPos[1]->Draw("PC");

  gRmsMax->SetLineColor(beamC);
  gRmsMax->SetMarkerColor(beamC);
  gRmsMax->SetLineWidth(1);
  gRmsMax->SetMarkerStyle(25);
  // gRmsMax->Draw("PC");

  // Emax value
  // New axis first:
  C->Update();  // Needed for the axis!

  margin = (maxField - minField)/10;
  if (margin==0) margin = 1; 
  rightmin = minField-margin;
  rightmax = maxField+margin;
  slope = (gPad->GetUymax() - gPad->GetUymin())/(rightmax-rightmin);
  TGaxis *axisEmax = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),gPad->GetUxmax(),
				gPad->GetUymax(),rightmin,rightmax,505,"+L");
  axisEmax->SetLineWidth(1);
  axisEmax->SetLineColor(kAzure-5);
  axisEmax->SetLabelColor(kAzure-5);
  axisEmax->SetTitleColor(kAzure-5);
  if(opt.Contains("units") && pData->GetPlasmaDensity()) 
    axisEmax->SetTitle("E_{min}^{z} [GV/m]");
  else
    axisEmax->SetTitle("E_{min}^{z} [E_{0}]");
  axisEmax->CenterTitle();
  axisEmax->SetTitleOffset(1.2);
  
  axisEmax->Draw();
  
  // Adjust the TGraph
  Double_t *x = gEfieldMaxValue[0]->GetX();
  Double_t *y = gEfieldMaxValue[0]->GetY();
  for(Int_t i=0;i<gEfieldMaxValue[0]->GetN();i++) {
    gEfieldMaxValue[0]->SetPoint(i,x[i],(y[i]-rightmin)*slope + gPad->GetUymin());
  }  
  gEfieldMaxValue[0]->SetLineColor(kAzure-5);
  gEfieldMaxValue[0]->SetMarkerColor(kAzure-5);
  gEfieldMaxValue[0]->SetLineWidth(1);
  gEfieldMaxValue[0]->SetMarkerStyle(21);
  gEfieldMaxValue[0]->Draw("PC");

  // Emax value
  // New axis first:
  C->Update();  // Needed for the axis!

  margin = (maxEmean - minEmean)/10;
  if (margin==0) margin = 1; 
  rightmin = minEmean-margin;
  rightmax = maxEmean+margin;
  slope = (gPad->GetUymax() - gPad->GetUymin())/(rightmax-rightmin);  
  rightmargin = 0.08;
  ux = gPad->PixeltoX(gPad->UtoPixel(1-rightmargin));
  TGaxis *axisEmeanB = new TGaxis(ux,gPad->GetUymin(),
				  ux,gPad->GetUymax(),rightmin,rightmax,505,"+L");
  axisEmeanB->SetLineWidth(1);
  axisEmeanB->SetLineColor(eneC);
  axisEmeanB->SetLabelColor(eneC);
  axisEmeanB->SetTitleColor(eneC);
  axisEmeanB->SetTitle("#LTE#GT_{max} [MeV]");
  axisEmeanB->CenterTitle();
  axisEmeanB->SetTitleOffset(0.8);
  axisEmeanB->Draw();
  
  // Adjust the TGraph for the right axis
  x = gEmeanMaxValue->GetX();
  y = gEmeanMaxValue->GetY();
  for(Int_t i=0;i<gEmeanMaxValue->GetN();i++) {
    gEmeanMaxValue->SetPoint(i,x[i],(y[i]-rightmin)*slope + gPad->GetUymin());
  }  
  gEmeanMaxValue->SetLineColor(eneC);
  gEmeanMaxValue->SetMarkerColor(eneC);
  gEmeanMaxValue->SetLineWidth(1);
  gEmeanMaxValue->SetMarkerStyle(21);
  gEmeanMaxValue->Draw("PC");

  rightmin = -1.;
  rightmax = 1.;
  slope = (gPad->GetUymax() - gPad->GetUymin())/(rightmax-rightmin);
  x = gEgain->GetX();
  y = gEgain->GetY();
  for(Int_t i=0;i<gEgain->GetN();i++) {
    gEgain->SetPoint(i,x[i],(y[i]-rightmin)*slope + gPad->GetUymin());
  }  
  gEgain->SetLineColor(eneC);
  gEgain->SetMarkerColor(eneC);
  gEgain->SetLineWidth(1);
  gEgain->SetMarkerStyle(24);
  gEgain->SetMarkerSize(0.6);
  gEgain->SetLineStyle(3);
  gEgain->Draw("PC");
 
  C->cd();

  ifile->Close();
  
  // Print to a file
  PlasmaGlob::imgconv(C,fOutName,opt);
  // ---------------------------------------------------------
  
}
