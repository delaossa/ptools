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
#include <TPad.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>

#include <TPaveText.h>
#include <TExec.h>
#include <TGaxis.h>

#include "PData.hh"
#include "PlasmaGlob.hh"

void PlotDensity1D( const TString &sim, Int_t time, Int_t Nbins=1, const TString &options="") {
  
#ifdef __CINT__  
  gSystem->Load("libplasma.so");
#endif

  PlasmaGlob::Initialize();
  
  // Init Units table
  PUnits::UnitsTable::Get();

  // Load PData
  PData *pData = PData::Get(sim.Data());
  pData->LoadFileNames(time);
  if(!pData->IsInit()) return; 

  TString opt = options;
 
  gStyle->SetPadRightMargin(0.20); // Margin right axis 
  if(opt.Contains("grid")) {
    gStyle->SetPadGridX(1);
    gStyle->SetPadGridY(1);
  }

  Bool_t CYL = kFALSE;
  if(sim.Contains("cyl")) { CYL = kTRUE; opt += "cyl"; } 
  
  Bool_t ThreeD = kFALSE;
  if(sim.Contains("3D")) ThreeD = kTRUE; 

  Bool_t INT = kTRUE; // Integrate instead of averaging.
  
  // Some plasma constants
  Float_t n0 = pData->GetPlasmaDensity(); 
  Float_t kp = pData->GetPlasmaK();       
  Float_t skindepth = (1/kp);               
  Float_t E0 = pData->GetPlasmaE0();

  // Some beam properties:
  Double_t Ebeam = pData->GetBeamEnergy();
  Double_t gamma = pData->GetBeamGamma();
  Double_t vbeam = pData->GetBeamVelocity();
  Double_t kbeta = PFunc::BeamBetatronWavenumber(gamma,n0);
  Double_t rms0  = pData->GetBeamRmsY() * kp;
  if(CYL)  rms0  = pData->GetBeamRmsR() * kp;
  
  cout << Form(" - Bunch gamma      = %8.4f", gamma ) << endl;
  cout << Form(" - Bunch velocity   = %8.4f c", vbeam ) << endl;
  cout << Form(" - Bunch betatron k = %8.4f mm-1", kbeta * PUnits::mm) << endl;
  cout << Form(" - Bunch RMS_0      = %8.4f um", rms0 * skindepth / PUnits::um) << endl;
  cout << endl;

  // Time in OU
  Float_t Time = pData->GetRealTime();
  // z start of the plasma in normalized units.
  Float_t zStartPlasma = pData->GetPlasmaStart()*kp;
  // z start of the beam in normalized units.
  Float_t zStartBeam = pData->GetBeamStart()*kp;
  
  if(opt.Contains("center")) {
    Time -= zStartPlasma;
    if(opt.Contains("comov"))      // Centers on the head of the beam.
      Time += zStartBeam;
  }
  
  // Get charge density histos
  Int_t Nspecies = pData->NSpecies();
  TH1F **hDen1D = new TH1F*[Nspecies];
  TH2F *hDen2D = NULL;
  for(Int_t i=0;i<Nspecies;i++) {
    hDen1D[i] = NULL;
    
    if(!pData->GetChargeFileName(i)) 
      continue;
    
    if(i==0) {
      if(ThreeD)
	hDen2D = pData->GetCharge2DSliceZY(i,-1,Nbins);
      else
	hDen2D = pData->GetCharge(i,opt);
      
      char hName[24];
      sprintf(hName,"hDen_%i",i);
      hDen2D->SetName(hName);
      hDen2D->GetXaxis()->CenterTitle();
      hDen2D->GetYaxis()->CenterTitle();
      hDen2D->GetZaxis()->CenterTitle();
      hDen2D->GetXaxis()->SetTitle("z [c/#omega_{p}]");
      hDen2D->GetYaxis()->SetTitle("y [c/#omega_{p}]");
      if(i==0)
	hDen2D->GetZaxis()->SetTitle("#LTn_{e}#GT [n_{0}]");
      else
	hDen2D->GetZaxis()->SetTitle("#LTn_{b}#GT [n_{0}]");
    }
    
    if(Nbins==0) {
      Nbins = TMath::Nint(rms0 / hDen2D->GetYaxis()->GetBinWidth(1)) ;
      // cout << Form(" Rms0 = %6.2f  Dx = %6.2f  Nbins = %4i .", 
      // 	   rms0, hDen2D[i]->GetYaxis()->GetBinWidth(1), Nbins) << endl;
    }
    
    // 1D histograms
    TString opth1 = opt;
    opth1 += "avg";
    if(CYL) opth1 += "cyl";
    
    if(ThreeD) {
      hDen1D[i] = pData->GetH1SliceZ3D(pData->GetChargeFileName(i)->c_str(),"charge",-1,Nbins,-1,Nbins);
    } else if(CYL) { // Cylindrical: The first bin with r>0 is actually the number 1 (not the 0).
      hDen1D[i] = pData->GetH1SliceZ(pData->GetChargeFileName(i)->c_str(),"charge",1,Nbins,opth1.Data());
    } else {         // 2D cartesian
      hDen1D[i] = pData->GetH1SliceZ(pData->GetChargeFileName(i)->c_str(),"charge",-1,Nbins,opth1.Data());
    }
    
    char hName[24];
    sprintf(hName,"hDen_%i",i);
    hDen1D[i]->SetName(hName);
    hDen1D[i]->GetXaxis()->SetTitle("z [c/#omega_{p}]");
    hDen1D[i]->GetYaxis()->SetTitle("#LTn_{e}#GT [n_{0}]");
    
  }
  
  // Get electric fields
  const Int_t Nfields = 2;
  TH1F **hE1D = new TH1F*[Nfields];
  for(Int_t i=0;i<Nfields;i++) {
    hE1D[i] = NULL;
    
    if(!pData->GetEfieldFileName(i)) 
      continue;

    // 1D histograms
    TString opth1 = opt;
    opth1 += "avg";
    if(CYL) opth1 += "cyl";
    
    char nam[3]; sprintf(nam,"e%i",i+1);
    if(ThreeD) {
      if(i==0) 
	hE1D[i] = pData->GetH1SliceZ3D(pData->GetEfieldFileName(i)->c_str(),nam,-1,Nbins,-1,Nbins);
      else  
	hE1D[i] = pData->GetH1SliceZ3D(pData->GetEfieldFileName(i)->c_str(),nam,-Nbins,Nbins,-Nbins,Nbins);
    } else if(CYL) { // Cylindrical: The firt bin with r>0 is actually the number 1 (not the 0).
      if(i==0) 
	hE1D[i] = pData->GetH1SliceZ(pData->GetEfieldFileName(i)->c_str(),nam,1,Nbins,opth1.Data());
      else
	hE1D[i] = pData->GetH1SliceZ(pData->GetEfieldFileName(i)->c_str(),nam,1,Nbins,opth1.Data());
    } else {         // 2D cartesian
      if(i==0) 
	hE1D[i] = pData->GetH1SliceZ(pData->GetEfieldFileName(i)->c_str(),nam,-1,Nbins,opth1.Data());
      else 
	hE1D[i] = pData->GetH1SliceZ(pData->GetEfieldFileName(i)->c_str(),nam,-Nbins,Nbins,opth1.Data());
    }
    
    
    char hName[24];
    sprintf(hName,"hE_%i",i);
    hE1D[i]->SetName(hName);
    hE1D[i]->GetXaxis()->SetTitle("z [c/#omega_{p}]");
    hE1D[i]->GetYaxis()->SetTitle("E [E_{e}]");
  }
  
  // Y range:
  Int_t NbinsY  = (Int_t) hDen2D->GetNbinsY();
  Int_t midyBin = NbinsY/2;
  if(!CYL && (Nbins >  midyBin) ) Nbins = midyBin;
  if(CYL  && (Nbins >= NbinsY ) ) Nbins = NbinsY-1;
  
  Int_t FirstyBin = midyBin + 1 - Nbins;
  Int_t LastyBin  = midyBin + Nbins;
  if(CYL) {
    FirstyBin = 1; 
    LastyBin  = Nbins;
  }
  
  if(LastyBin>=NbinsY) LastyBin = NbinsY - 1;
  Float_t ymin = hDen2D->GetYaxis()->GetBinLowEdge(FirstyBin);
  Float_t ymax = hDen2D->GetYaxis()->GetBinLowEdge(LastyBin+1);

  cout << Form(" Nbins = %i. Firstbin = %i Lastbin = %i  ->  ymin = %7.2f , ymax = %7.2f",Nbins,FirstyBin,LastyBin,ymin,ymax) << endl;
  // ----

  // Tunning the Histograms
  // ---------------------

  // Chaning to user units: 

  // cout << Form(" n0 = %10e ", n0 * PUnits::cm3) << endl;
  
  if(opt.Contains("units") && n0) {
    
    Int_t NbinsX = hDen2D->GetNbinsX();
    Float_t xMin = skindepth * hDen2D->GetXaxis()->GetXmin() / PUnits::mm;
    Float_t xMax = skindepth * hDen2D->GetXaxis()->GetXmax() / PUnits::mm;
    Int_t NbinsY = hDen2D->GetNbinsY();
    Float_t yMin = skindepth * hDen2D->GetYaxis()->GetXmin() / PUnits::mm;
    Float_t yMax = skindepth * hDen2D->GetYaxis()->GetXmax() / PUnits::mm;
    hDen2D->SetBins(NbinsX,xMin,xMax,NbinsY,yMin,yMax);
    for(Int_t j=0;j<hDen2D->GetNbinsX();j++) {
      for(Int_t k=0;k<hDen2D->GetNbinsY();k++) {
	hDen2D->SetBinContent(j,k, hDen2D->GetBinContent(j,k) * n0 / (1e15/PUnits::cm3) );
      }
    }
    
    hDen2D->GetYaxis()->SetTitle("y [mm]");      
    if(opt.Contains("comov"))
      hDen2D->GetXaxis()->SetTitle("#zeta [mm]");
    else
      hDen2D->GetXaxis()->SetTitle("z [mm]");
    
    hDen2D->GetZaxis()->SetTitle("#LTn_{b}#GT [10^{15}/cm^{3}]"); 
  
    
    for(Int_t i=0;i<Nspecies;i++) {
    
      Int_t NbinsX = hDen1D[i]->GetNbinsX();
      Float_t xMin = skindepth * hDen1D[i]->GetXaxis()->GetXmin() / PUnits::mm;
      Float_t xMax = skindepth * hDen1D[i]->GetXaxis()->GetXmax() / PUnits::mm;
      hDen1D[i]->SetBins(NbinsX,xMin,xMax);
      for(Int_t j=0;j<hDen1D[i]->GetNbinsX();j++) {
	Float_t bincontent = (hDen1D[i]->GetBinContent(j) * n0 / (1e15/PUnits::cm3));
	hDen1D[i]->SetBinContent(j,bincontent);
      }
      
      hDen1D[i]->GetYaxis()->SetTitle("#LTn_{e}#GT [10^{15}/cm^{3}]");
      
      if(opt.Contains("comov"))
	hDen1D[i]->GetXaxis()->SetTitle("#zeta [mm]");
      else
	hDen1D[i]->GetXaxis()->SetTitle("z [mm]");
    
    }


    for(Int_t i=0;i<Nfields;i++) {
      Int_t NbinsX = hE1D[i]->GetNbinsX();
      Float_t xMin = skindepth * hE1D[i]->GetXaxis()->GetXmin() / PUnits::mm;
      Float_t xMax = skindepth * hE1D[i]->GetXaxis()->GetXmax() / PUnits::mm;
      hE1D[i]->SetBins(NbinsX,xMin,xMax);
      
      for(Int_t j=0;j<hE1D[i]->GetNbinsX();j++) {
	hE1D[i]->SetBinContent(j, hE1D[i]->GetBinContent(j) * ( E0 / (PUnits::GV/PUnits::m) ) );
      }
      
      hE1D[i]->GetYaxis()->SetTitle("E [GV/m]");
      
      if(opt.Contains("comov"))
	hE1D[i]->GetXaxis()->SetTitle("#zeta [mm]");
      else
	hE1D[i]->GetXaxis()->SetTitle("z [mm]");
    }
  }
  
  // Set the range of the histogram for maximum constrast
  Float_t density = 1;
  if(opt.Contains("units") && n0)
    density = 1e-15 * 1e-6 * n0;
  
  Float_t Max  = 1.1*hDen1D[0]->GetMaximum();
  Float_t Min  = 0;
  // Float_t Base = density;
  // Float_t Min  = 2.* Base - Max;
  // if(Max >= 2. * Base) {
  //   Max = 2. * Base;
  //   Min = 2. * Base  - Max;
  // } else if(Max<1.0 * Base) {
  //   Max = 1.1 * Base;
  //   Min = 0.;
  // }
 
  hDen1D[0]->GetYaxis()->SetRangeUser(0.,Max);  
  
  // Plotting
  // -----------------------------------------------

  // Canvas setup
  TCanvas *C = new TCanvas("C","Charge density and Electric field on axis",850,500);

  TPaveText *textTime = new TPaveText(0.63,0.87,0.78,0.92,"NDC");
  PlasmaGlob::SetPaveTextStyle(textTime); 
  textTime->SetTextColor(kGray+2);
  char ctext[128];

  if(opt.Contains("units") && n0) 
    sprintf(ctext,"Z = %5.1f mm", 1e3 * skindepth * Time);
  else
    sprintf(ctext,"T = %5.1f 1/#omega_{p}",Time);
  textTime->AddText(ctext);

  
  TPaveText *textRange = new TPaveText(0.13,0.87,0.38,0.92,"NDC");
  PlasmaGlob::SetPaveTextStyle(textRange,12); 
  textRange->SetTextColor(kGray+2);
  if(opt.Contains("units") && n0)
    sprintf(ctext,"%5.3f < y < %5.3f mm",ymin,ymax);
  else
    sprintf(ctext,"%5.3f < y < %5.3f c/#omega_{p}",ymin,ymax);
  textRange->AddText(ctext);

  // Actual Plotting!
  // ------------------------------------------------------------

  // Output file
  TString fOutName = Form("./%s/Plots/Density1D/Density1D",sim.Data());
  fOutName += Form("-%s_%i",sim.Data(),time);

  // Colors
  Int_t plasmaC = kGray+1;
  Int_t beamC   = kAzure-5;
  Int_t fieldC  = kOrange+10;
  Int_t fieldCb = kGray+1;

  C->cd(0);
  gPad->SetFrameLineWidth(2);  

  hDen1D[0]->SetLineColor(plasmaC);
  hDen1D[0]->SetLineWidth(1);
  hDen1D[0]->Draw("C");
  hDen1D[0]->GetYaxis()->CenterTitle();
  hDen1D[0]->GetXaxis()->CenterTitle();
  
  C->Update();

  TLine *line0 = new TLine(hDen1D[0]->GetXaxis()->GetXmin(),
			  (gPad->GetUymin()+gPad->GetUymax())/2.,
			  hDen1D[0]->GetXaxis()->GetXmax(),
			  (gPad->GetUymin()+gPad->GetUymax())/2.);
  line0->SetLineColor(kGray+1);
  line0->SetLineStyle(2);
  line0->Draw();

  Float_t rightmax = 2.5 * hDen1D[1]->GetMaximum();
  Float_t slope = (gPad->GetUymax() - gPad->GetUymin())/rightmax;

  for(Int_t i=0;i<hDen1D[1]->GetNbinsX();i++) {
    hDen1D[1]->SetBinContent(i+1,hDen1D[1]->GetBinContent(i+1)*slope + Min);
  }

  hDen1D[1]->SetLineWidth(2);
  hDen1D[1]->SetLineColor(beamC);
  hDen1D[1]->Draw("same C");
  // hTest->Draw("same");
 
  //draw an axis on the right side
  TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),gPad->GetUxmax(),
			    gPad->GetUymax(),0,rightmax,505,"+L");
  
  axis->SetLineWidth(1);
  axis->SetLineColor(beamC);
  axis->SetLabelColor(beamC);
  if(opt.Contains("units") && n0) 
    axis->SetTitle("#LTn_{b}#GT [10^{15}/cm^{3}]");
  else
    axis->SetTitle("#LTn_{b}#GT [n_{0}]");
  axis->CenterTitle();
  axis->SetTitleColor(beamC);
  axis->SetTitleOffset(1.2);
  
  axis->Draw();

  // Longitudinal Electric field
  Float_t factor = 1.5;
  Float_t rightmin = factor * hE1D[0]->GetMinimum();
  rightmax = factor * hE1D[0]->GetMaximum();
  if(hE1D[1]) {
    if(hE1D[1]->GetMaximum() > hE1D[0]->GetMaximum())
      rightmax = factor * hE1D[1]->GetMaximum();
  }
  
  if(rightmax > TMath::Abs(rightmin)) rightmin = -rightmax;
  else rightmax = - rightmin;
  slope = (gPad->GetUymax() - gPad->GetUymin())/(rightmax-rightmin);
  
  for(Int_t i=0;i<hE1D[0]->GetNbinsX();i++) {
    hE1D[0]->SetBinContent(i+1,(hE1D[0]->GetBinContent(i+1)-rightmin)*slope + Min);
  }
  
  hE1D[0]->SetLineStyle(1);
  hE1D[0]->SetLineWidth(2);
  hE1D[0]->SetLineColor(fieldC);
  hE1D[0]->Draw("same C");

  // Transverse field
  for(Int_t i=0;i<hE1D[1]->GetNbinsX();i++) {
    hE1D[1]->SetBinContent(i+1,(hE1D[1]->GetBinContent(i+1)-rightmin)*slope + Min);
  }
  
  hE1D[1]->SetLineStyle(2);
  hE1D[1]->SetLineWidth(1);
  hE1D[1]->SetLineColor(fieldC);
  hE1D[1]->Draw("same C");

  //draw an axis on the right side
  Float_t rightmargin = 0.08;
  Float_t ux = gPad->PixeltoX(gPad->UtoPixel(1-rightmargin));
  TGaxis *axisE = new TGaxis(ux,gPad->GetUymin(),ux,
			     gPad->GetUymax(),rightmin,rightmax,505,"+L");
  
  axisE->SetLineWidth(1);
  axisE->SetLineColor(fieldC);
  axisE->SetLabelColor(fieldC);
  axisE->SetTitleColor(fieldC);
  if(opt.Contains("units") && n0) 
    axisE->SetTitle("E [GV/m]");
  else
    axisE->SetTitle("E [E_{0}]");
  axisE->CenterTitle();
  axisE->SetTitleOffset(0.8);

  axisE->Draw();
  
  textTime->Draw();

  textRange->Draw();

  C->cd();

  // Print to a file
  PlasmaGlob::imgconv(C,fOutName,opt);
  // ---------------------------------------------------------
}
