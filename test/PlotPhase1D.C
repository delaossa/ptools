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

void PlotPhase1D( const TString &sim, Int_t time, Int_t Nbins=1, const TString &opt="") { 
  
  PlasmaGlob::Initialize();

  gStyle->SetPadLeftMargin(0.10);  // Margin left axis 
  gStyle->SetPadRightMargin(0.20); // Margin right axis 

  // Load PData
  PData *pData = PData::Get(sim.Data());
  pData->LoadFileNames(time);
  if(!pData->IsInit()) return;

  Bool_t ThreeD = kFALSE;
  if(sim.Contains("3D")) ThreeD = kTRUE; 
  
  // Get charge density histos
  Int_t Nspecies = pData->NSpecies();
  TH1F **hDen1D = new TH1F*[Nspecies];
  TH2F *hDen2D = NULL;
  for(Int_t i=0;i<Nspecies;i++) {

    hDen1D[i] = NULL;
    if(!pData->GetChargeFileName(i)) 
      continue;
    
    // Get 2D charge density histogram of the bunch.
    // We use it to calculate the rms of the charge distribution
    // in x2 as a function of x1.
    if(pData->GetSpeciesName(i).find("beam")!=string::npos) {
      
      if(ThreeD)
	hDen2D = pData->GetCharge2DSliceZY(i,-1,Nbins);
      else
	hDen2D = pData->GetCharge(i,opt);
      
      hDen2D->SetName("hDen2D");
      hDen2D->GetXaxis()->SetTitle("z [c/#omega_{p}]");
      hDen2D->GetYaxis()->SetTitle("y [c/#omega_{p}]");
      
      if(opt.Contains("comov")) {
	Int_t NbinsX = hDen2D->GetNbinsX();
	Float_t xMin = hDen2D->GetXaxis()->GetXmin()-hDen2D->GetXaxis()->GetXmax();
	Float_t xMax = hDen2D->GetXaxis()->GetXmax()-hDen2D->GetXaxis()->GetXmax();
	Int_t NbinsY = hDen2D->GetNbinsY();
	Float_t yMin = hDen2D->GetYaxis()->GetXmin();
	Float_t yMax = hDen2D->GetYaxis()->GetXmax();
	hDen2D->SetBins(NbinsX,xMin,xMax,NbinsY,yMin,yMax);

	hDen2D->GetXaxis()->SetTitle("#zeta [c/#omega_{p}]");
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
	
      }
      
    }  
    
    // 1D charge density histos around the axis.
    if(!ThreeD) 
      hDen1D[i] = pData->GetH1SliceZ(pData->GetChargeFileName(i)->c_str(),"charge",-1,Nbins);
    else
      hDen1D[i] = pData->GetH1SliceZ3D(pData->GetChargeFileName(i)->c_str(),"charge",-1,Nbins,-1,Nbins);
    
    char hName[24];
    sprintf(hName,"hDen_%i",i);
    hDen1D[i]->SetName(hName);
    hDen1D[i]->GetXaxis()->CenterTitle();
    hDen1D[i]->GetYaxis()->CenterTitle();
    hDen1D[i]->GetXaxis()->SetTitle("z [c/#omega_{p}]");
    hDen1D[i]->GetYaxis()->SetTitle("dQ/dz [ou]");
    
    // Change to co-moving coordinate
    if(opt.Contains("comov")) {
      Int_t NbinsX = hDen1D[i]->GetNbinsX();
      Float_t xMin = hDen1D[i]->GetXaxis()->GetXmin()-hDen1D[i]->GetXaxis()->GetXmax();
      Float_t xMax = hDen1D[i]->GetXaxis()->GetXmax()-hDen1D[i]->GetXaxis()->GetXmax();
      hDen1D[i]->SetBins(NbinsX,xMin,xMax);
      hDen1D[i]->GetXaxis()->SetTitle("#zeta [c/#omega_{p}]");
      
    }
    
    // Chaning to user units: 
    if(opt.Contains("units") && pData->GetPlasmaDensity()) {
      Int_t NbinsX = hDen1D[i]->GetNbinsX();
      Float_t xMin = 1e3 * pData->GetPlasmaSkinDepth() * hDen1D[i]->GetXaxis()->GetXmin();
      Float_t xMax = 1e3 * pData->GetPlasmaSkinDepth() * hDen1D[i]->GetXaxis()->GetXmax();
      hDen1D[i]->SetBins(NbinsX,xMin,xMax);
      
      for(Int_t j=0;j<hDen1D[i]->GetNbinsX();j++) {
	hDen1D[i]->SetBinContent(j,1e-15 * 1e-6 * pData->GetPlasmaDensity() * hDen1D[i]->GetBinContent(j));
      }
      
      hDen1D[i]->GetYaxis()->SetTitle("dQ/dz [10^{15}/cc]");
      
      if(opt.Contains("comov"))
	hDen1D[i]->GetXaxis()->SetTitle("#zeta [mm]");
      else
	hDen1D[i]->GetXaxis()->SetTitle("z [mm]");
      
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
	hE1D[i] = pData->GetH1SliceZ(pData->GetEfieldFileName(i)->c_str(),nam,-1,Nbins);
      }
      else { 
	hE1D[i] = pData->GetH1SliceZ(pData->GetEfieldFileName(i)->c_str(),nam,-Nbins,Nbins);
      }
    } else {
      if(i==0) {
	hE1D[i] = pData->GetH1SliceZ3D(pData->GetEfieldFileName(i)->c_str(),nam,-1,Nbins,-1,Nbins);
      }
      else { 
	hE1D[i] = pData->GetH1SliceZ3D(pData->GetEfieldFileName(i)->c_str(),nam,-Nbins,Nbins,-Nbins,Nbins);
      }
      
    }
    
    char hName[24];
    sprintf(hName,"hE_%i",i);
    hE1D[i]->SetName(hName);
    hE1D[i]->GetXaxis()->SetTitle("z [c/#omega_{p}]");
    hE1D[i]->GetYaxis()->SetTitle("E [E_{0}]");

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

      hE1D[i]->GetYaxis()->SetTitle("E [GeV/m]");
      
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
  TString filename = Form("./%s/Plots/Phase1D/Phase1D-%s.root",sim.Data(),sim.Data());
  TFile * ifile = (TFile*) gROOT->GetListOfFiles()->FindObject(filename);
  if (!ifile) {
    TString f = filename;
    TString dir2 = f.Remove( f.Last( '/' ), f.Length() - f.Last( '/' ) );
    TString dir1 = f.Remove( f.Last( '/' ), f.Length() - f.Last( '/' ) );
    gSystem->mkdir( dir1 );
    gSystem->mkdir( dir2 );
    ifile = new TFile(filename,"UPDATE");
  }  
  TGraph *gEfieldMaxPos0 =  (TGraph*) ifile->Get("gEfieldMaxPos_0");
  if(gEfieldMaxPos0) {
    Double_t *y = gEfieldMaxPos0->GetY();
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
    
    htemp->Smooth(1,"R");
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
  
  
  // Charge density
  Float_t *hDenMax = new Float_t[Nspecies];
  Float_t *hDenMaxValue = new Float_t[Nspecies];  
  for(Int_t i=0;i<Nspecies;i++) {
    if(!hDen1D[i]) continue;
    hDenMax[i] = hDenMaxValue[i] = -999.;  
    Float_t xCenter = (hDen1D[i]->GetXaxis()->GetXmin()+hDen1D[i]->GetXaxis()->GetXmax())/2.;
    TH1F *htemp = (TH1F*) hDen1D[i]->Clone("htemp");
    htemp->GetXaxis()->SetRangeUser(xs1min,xs1max);
    htemp->Smooth(1,"R");
    Int_t binMax = htemp->GetMinimumBin();
    hDenMax[i] = (Float_t) htemp->GetBinCenter(binMax);
    hDenMaxValue[i] = (Float_t) htemp->GetBinContent(binMax);
    delete htemp;
  }
  
  Float_t hRmsMax = -999.;
  Float_t hRmsMaxValue = -999.;
  if(hRms) {
    Float_t xCenter = (hRms->GetXaxis()->GetXmin()+hRms->GetXaxis()->GetXmax())/2.;
    TH1F *htemp = (TH1F*) hRms->Clone("htemp");
    
    htemp->GetXaxis()->SetRangeUser(xs1min,xs1max);
    htemp->Smooth(1,"R");
    Int_t binMax = htemp->GetMaximumBin();
    hRmsMax = (Float_t) htemp->GetBinCenter(binMax);
    hRmsMaxValue = (Float_t) htemp->GetBinContent(binMax);
    delete htemp;
  }
  
  
  // Tunning the Histograms
  // ---------------------
  
  Float_t Time = pData->GetRealTime();
  
  // Set the range of the histogram for maximum constrast
  Float_t density = 1;
  if(opt.Contains("units") && pData->GetPlasmaDensity())
    density = 1e-15 * 1e-6 * pData->GetPlasmaDensity();
  
  Float_t Max  = 1.1 * hDen1D[0]->GetMaximum();
  Float_t Base = density;
  Float_t Min  = 2.* Base - Max;
  if(Max >= 2. * Base) {
    Max = 2. * Base;
    Min = 2. * Base  - Max;
  } else if(Max<1.0 * Base) {
    Max = 1.1 * Base;
    Min = 0.;
  }
  
  hDen1D[0]->GetYaxis()->SetRangeUser(Min,Max);  
  

  // Plotting
  // -----------------------------------------------

  // Output file
  TString fOutName = Form("./%s/Plots/Phase1D/Phase1D",sim.Data());
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
  Int_t fieldC  = kOrange+10;
  Int_t fieldCb = kGray+1;
  Int_t fieldMaxC = kOrange+10;

  // Actual Plotting!
  // ------------------------------------------------------------

  // More makeup
  C->cd(1);
  gPad->SetGridy(0);
  gPad->SetGridx(0);
  gPad->SetFrameLineWidth(2);  

  hDen1D[0]->SetLineColor(plasmaC);
  hDen1D[0]->SetLineWidth(1);
  hDen1D[0]->Draw("C");

  C->Update();

  TLine *line0 = new TLine(hDen1D[0]->GetXaxis()->GetXmin(),
  			  (gPad->GetUymin()+gPad->GetUymax())/2.,
  			  hDen1D[0]->GetXaxis()->GetXmax(),
  			  (gPad->GetUymin()+gPad->GetUymax())/2.);
  line0->SetLineColor(kGray+1);
  line0->SetLineStyle(2);
  line0->Draw();

  TLine *lineMax0 = new TLine(hDenMax[0],gPad->GetUymin(),hDenMax[0],hDenMaxValue[0]);
  lineMax0->SetLineStyle(1);
  lineMax0->SetLineColor(plasmaC);
  //  lineMax0->Draw();

  TMarker *markDenMax0 = new TMarker(hDenMax[0],hDenMaxValue[0],20);
  markDenMax0->SetMarkerColor(plasmaC);
  markDenMax0->Draw();

  // Bunch RMS
  // Trick for nicer plotting
  Float_t bin1 = -1;
  Float_t bin2 = -1;
  for(Int_t i=0;i<hRms->GetNbinsX();i++) 
    if(bin1==-1 && hRms->GetBinContent(i+1)>0) bin1 = hRms->GetBinCenter(i+1);  
  for(Int_t i=hRms->GetNbinsX()-1;i>=0;i--) 
    if(bin2==-1 && hRms->GetBinContent(i+1)>0) bin2 = hRms->GetBinCenter(i+1);
  hRms->GetXaxis()->SetRangeUser(bin1,bin2);

  Float_t rightmin = 0;
  Float_t rightmax = 2.5 * hRms->GetMaximum();
  Float_t slope = (gPad->GetUymax() - gPad->GetUymin())/(rightmax-rightmin);
  
  for(Int_t i=0;i<hRms->GetNbinsX();i++) {
    hRms->SetBinContent(i+1,(hRms->GetBinContent(i+1)-rightmin)*slope + gPad->GetUymin());
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

  TMarker *markRmsMax = new TMarker(hRmsMax,(hRmsMaxValue-rightmin)*slope + gPad->GetUymin(), 20);
  markRmsMax->SetMarkerColor(beamC);
  //markRmsMax->Draw();

  // Line around 1
  Float_t value = (1-rightmin)*slope + gPad->GetUymin();
  TLine *line1 = new TLine(hDen1D[0]->GetXaxis()->GetXmin(),value,
			   hDen1D[0]->GetXaxis()->GetXmax(),value);
  line1->SetLineColor(beamC);
  line1->SetLineStyle(2);
  line1->Draw();

  // Longitudinal Electric field
  Float_t factor = 1.5;
  if(hE1D[1]) {
    rightmin = factor * hE1D[0]->GetMinimum();
    rightmax = factor * hE1D[0]->GetMaximum();
    
    if(hE1D[1]->GetMaximum() > hE1D[0]->GetMaximum())
      rightmax = factor * hE1D[1]->GetMaximum();
  }
  
  if(rightmax > TMath::Abs(rightmin)) rightmin = -rightmax;
  else rightmax = - rightmin;
  slope = (gPad->GetUymax() - gPad->GetUymin())/(rightmax-rightmin);
  
  for(Int_t i=0;i<hE1D[0]->GetNbinsX();i++) {
    hE1D[0]->SetBinContent(i+1,(hE1D[0]->GetBinContent(i+1)-rightmin)*slope + gPad->GetUymin());
  }
  
  hE1D[0]->SetLineStyle(1);
  hE1D[0]->SetLineColor(fieldC);
  hE1D[0]->Draw("same C");

  TLine *lineEMax0 = new TLine(EfieldMaxPos[0],gPad->GetUymin(),EfieldMaxPos[0],(EfieldMaxValue[0]-rightmin)*slope + gPad->GetUymin());
  lineEMax0->SetLineStyle(1);
  lineEMax0->SetLineColor(fieldC);
  //  lineEMax0->Draw();

  TMarker *markEMax0 = new TMarker(EfieldMaxPos[0],(EfieldMaxValue[0]-rightmin)*slope + gPad->GetUymin(), 20);
  markEMax0->SetMarkerColor(fieldC);
  markEMax0->Draw();
  
  // Transverse field
  for(Int_t i=0;i<hE1D[1]->GetNbinsX();i++) {
    hE1D[1]->SetBinContent(i+1,(hE1D[1]->GetBinContent(i+1)-rightmin)*slope + gPad->GetUymin());
  }
  
  hE1D[1]->SetLineStyle(3);
  hE1D[1]->SetLineColor(fieldCb);
  hE1D[1]->Draw("same C");

  TLine *lineEMax1 = new TLine(EfieldMaxPos[1],gPad->GetUymin(),EfieldMaxPos[1],(EfieldMaxValue[1]-rightmin)*slope+Min);
  lineEMax1->SetLineStyle(3);
  lineEMax1->SetLineColor(fieldC);
  //  lineEMax1->Draw();

  TMarker *markEMax1 = new TMarker(EfieldMaxPos[1],(EfieldMaxValue[1]-rightmin)*slope + gPad->GetUymin(), 24);
  markEMax1->SetMarkerColor(fieldCb);
  markEMax1->Draw();

  //draw an axis on the right side
  Float_t rightmargin = 0.08;
  Float_t ux = gPad->PixeltoX(gPad->UtoPixel(1-rightmargin));
  TGaxis *axisE = new TGaxis(ux,gPad->GetUymin(),ux,
			     gPad->GetUymax(),rightmin,rightmax,505,"+L");
  
  axisE->SetLineWidth(1);
  axisE->SetLineColor(fieldC);
  axisE->SetLabelColor(fieldC);
  axisE->SetTitleColor(fieldC);
  if(opt.Contains("units") && pData->GetPlasmaDensity()) 
    axisE->SetTitle("E [GeV/m]");
  else
    axisE->SetTitle("E [E_{0}]");
  axisE->CenterTitle();
  axisE->SetTitleOffset(0.8);

  axisE->Draw();
  
  textTime->Draw();

  // Fitting...
  // TF1 *fitFunction = new TF1("fitFunction","[0]+[1]*(x-[2])*sin([3]*(x-[4]))",hDen1D[0]->GetXaxis()->GetXmin(),hDen1D[0]->GetXaxis()->GetXmax());
  // fitFunction->SetParameters(1.,hDen1D[0]->GetMaximum()-1.,-5.,1.,-5.);
  // hDen1D[0]->Fit("fitFunction","MIR0","",-35.,-10);
  // fitFunction->SetLineStyle(3);
  // fitFunction->SetLineWidth(1);
  // fitFunction->Draw("same");

  Int_t nPoints = 0;  
  TGraph **gDenMax = new TGraph*[Nspecies];

  for(Int_t i=0;i<Nspecies;i++) {
    if(!hDen1D[i]) continue;
    char gName[24];
    sprintf(gName,"gDenMax_%i",i);
    gDenMax[i] = (TGraph*) ifile->Get(gName);
    if(gDenMax[i]==NULL) {
      gDenMax[i] = new TGraph();
      gDenMax[i]->SetName(gName);
    } else {
      nPoints = gDenMax[i]->GetN(); 
    }  
    gDenMax[i]->Set(nPoints+1);
    if(opt.Contains("units") && pData->GetPlasmaDensity()) 
      gDenMax[i]->SetPoint(nPoints, 1e3 * pData->GetPlasmaSkinDepth() * Time,hDenMax[i]);
    else
      gDenMax[i]->SetPoint(nPoints,Time,hDenMax[i]);
    
    if(opt.Contains("units") && pData->GetPlasmaDensity()) {
      gDenMax[i]->GetYaxis()->SetTitle("#zeta_{min} [mm]");
      gDenMax[i]->GetXaxis()->SetTitle("L [mm]");
    } else {
      gDenMax[i]->GetYaxis()->SetTitle("#zeta_{min} [c/#omega_{p}]");
      gDenMax[i]->GetXaxis()->SetTitle("T [c/#omega_{p}]");
    }

    gDenMax[i]->Write(gDenMax[i]->GetName(),TObject::kOverwrite);
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
      gEfieldMaxValue[i]->GetYaxis()->SetTitle("E_{min} [GeV/m]");
      gEfieldMaxValue[i]->GetXaxis()->SetTitle("L [mm]");
    } else {
      gEfieldMaxValue[i]->GetYaxis()->SetTitle("E_{min} [E_{0}]");
      gEfieldMaxValue[i]->GetXaxis()->SetTitle("T [c/#omega_{p}]");
    }

    gEfieldMaxValue[i]->Write(gEfieldMaxValue[i]->GetName(),TObject::kOverwrite);
    
  }


  C->cd(2);
  gPad->SetGridy(0);
  gPad->SetGridx(0);
  gPad->SetFrameLineWidth(2);  

  Float_t xCenter = (hDen1D[0]->GetXaxis()->GetXmin()+hDen1D[0]->GetXaxis()->GetXmax())/2.; 
  Float_t minPhase = 99.;
  Float_t maxPhase = -99.;
  Float_t minField = 99.;
  Float_t maxField = -99.;

  Double_t *yEfieldMaxPos_0   = gEfieldMaxPos[0]->GetY();
  Double_t *yEfieldMaxPos_1   = gEfieldMaxPos[1]->GetY();
  Double_t *yEfieldMaxValue_0 = gEfieldMaxValue[0]->GetY();

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

  gDenMax[0]->SetLineColor(plasmaC);
  gDenMax[0]->SetMarkerColor(plasmaC);
  gDenMax[0]->SetLineWidth(1);
  gDenMax[0]->SetMarkerStyle(20);
  // gDenMax[0]->Draw("APC");
  
  
  gDenMax[1]->SetLineColor(beamC);
  gDenMax[1]->SetMarkerColor(beamC);
  // gDenMax[1]->Draw("PC");

  gRmsMax->SetLineColor(beamC);
  gRmsMax->SetMarkerColor(beamC);
  gRmsMax->SetLineWidth(1);
  gRmsMax->SetMarkerStyle(24);
  //gRmsMax->Draw("PC");

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
  axisEmax->SetLineColor(fieldMaxC);
  axisEmax->SetLabelColor(fieldMaxC);
  axisEmax->SetTitleColor(fieldMaxC);
  if(opt.Contains("units") && pData->GetPlasmaDensity()) 
    axisEmax->SetTitle("E_{min} [GeV/m]");
  else
    axisEmax->SetTitle("E_{min} [E_{0}]");
  axisEmax->CenterTitle();
  axisEmax->SetTitleOffset(1.4);

  axisEmax->Draw();

  // Adjust the TGraph
  Double_t *x = gEfieldMaxValue[0]->GetX();
  Double_t *y = gEfieldMaxValue[0]->GetY();
  for(Int_t i=0;i<gEfieldMaxValue[0]->GetN();i++) {
    gEfieldMaxValue[0]->SetPoint(i,x[i],(y[i]-rightmin)*slope + gPad->GetUymin());
  }  
  gEfieldMaxValue[0]->SetLineColor(fieldMaxC);
  gEfieldMaxValue[0]->SetMarkerColor(fieldMaxC);
  gEfieldMaxValue[0]->SetLineWidth(1);
  gEfieldMaxValue[0]->SetMarkerStyle(21);
  gEfieldMaxValue[0]->Draw("PC");
 

  ifile->Close();

  // Print to a file
  PlasmaGlob::imgconv(C,fOutName,opt);
  
  C->cd();
  // ---------------------------------------------------------
 
}
