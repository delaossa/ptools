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
#include <TProfile.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>

#include <TPaveText.h>
#include <TExec.h>
#include <TGaxis.h>

#include "PData.hh"
#include "PlasmaGlob.hh"

void PlotCharge1D( const TString &sim, Int_t time, Int_t Nbins=1, const TString &options="") {
  
  PlasmaGlob::Initialize();
  
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

  // Some plasma constants
  Double_t n0 = pData->GetPlasmaDensity(); 
  Double_t kp = pData->GetPlasmaK();       
  Double_t skindepth = (1/kp);               
  Double_t E0 = pData->GetPlasmaE0();

  // Some beam properties:
  Double_t Ebeam = pData->GetBeamEnergy();
  Double_t gamma = Ebeam / PConst::ElectronMassE;
  //  Double_t vbeam = TMath::Sqrt(1 - 1/(gamma*gamma));
  
  // Time in OU
  Double_t Time = pData->GetRealTime();
  // z start of the plasma in normalized units.
  Double_t zStartPlasma = pData->GetPlasmaStart()*kp;
  // z start of the beam in normalized units.
  Double_t zStartBeam = pData->GetBeamStart()*kp;
  Time -= zStartPlasma - zStartBeam;
  
  // Get charge density histos
  Int_t Nspecies = pData->NSpecies();
  TH1F  *hDen1D = NULL;
  TH2F  *hDen2D = NULL;
  for(Int_t i=0;i<Nspecies;i++) {
    hDen1D = NULL;
  
    if(!pData->GetChargeFileName(i)) 
      continue;
    
    // Only the beam.
    if(pData->GetSpeciesName(i).find("beam")!=string::npos) {
      
      if(ThreeD)
	hDen2D = pData->GetCharge2DSliceZY(i,-1,Nbins);
      else
	hDen2D = pData->GetCharge(i,opt);
      
      char hName[24];
      sprintf(hName,"hDen2D_%i",i);
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
      
      // 1D histograms
      TString opth1 = opt;
      
      // This option SUMS the values of the bins within the range given by Nbins.
      opth1 += "int"; 
      
      // This option AVERAGES the values of the bins within the range given by Nbins.
      // opth1 += "avg"; 
      
      // Flag for Cylindrical coordinates (needs special treatment).
      if(CYL) opth1 += "cyl";
      
      if(ThreeD) {
	hDen1D = pData->GetH1SliceZ3D(pData->GetChargeFileName(i)->c_str(),"charge",-1,Nbins,-1,Nbins);
      } else if(CYL) { // Cylindrical: The first bin with r>0 is actually the number 1 (not the 0).
	hDen1D = pData->GetH1SliceZ(pData->GetChargeFileName(i)->c_str(),"charge",1,Nbins,opth1.Data());
      } else { // 2D cartesian
	hDen1D = pData->GetH1SliceZ(pData->GetChargeFileName(i)->c_str(),"charge",-1,Nbins,opth1.Data());
      }
      
      sprintf(hName,"hDen1D_%i",i);
      hDen1D->SetName(hName);
      hDen1D->GetXaxis()->SetTitle("z [c/#omega_{p}]");

      if(CYL)
	hDen1D->GetYaxis()->SetTitle("N_{e} [n_{0} / (k_{p}^{2} rad)]");
      else
	hDen1D->GetYaxis()->SetTitle("N_{e} [n_{0} / k_{p}]");
    }
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
    
    // This option SUMS the values of the bins within the range given by Nbins.
    // opth1 += "int"; 
    
    // This option AVERAGES the values of the bins within the range given by Nbins.
    opth1 += "avg"; 
    
    // Flag for Cylindrical coordinates (needs special treatment).
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
    } else { // 2D cartesian
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
  

  // Manual binning range (as a function of Nbins).
  Int_t NbinsY  = (Int_t) hDen2D->GetNbinsY();
  Int_t midyBin = NbinsY/2;
  if(!CYL && (Nbins >  midyBin) ) Nbins = midyBin;
  if(CYL  && (Nbins >= NbinsY ) ) Nbins = NbinsY;
  
  Int_t FirstyBin = midyBin + 1 - Nbins;
  Int_t LastyBin  = midyBin + Nbins;
  if(CYL) {
    FirstyBin = 1; 
    LastyBin  = Nbins;
  }
  
  if(LastyBin>=NbinsY) LastyBin = NbinsY - 1;
   
  // Actual range in Osiris units:
  Double_t ymin = hDen2D->GetYaxis()->GetBinLowEdge(FirstyBin);
  Double_t ymax = hDen2D->GetYaxis()->GetBinUpEdge(LastyBin);
  cout << Form(" Nbins = %i (%3i,%3i) : (%6.3f,%6.3f)  ",LastyBin-FirstyBin+1,FirstyBin,LastyBin,ymin,ymax) << endl;

  
  // Integrating beam charge:
  // ---------------------------------
  Double_t Q = 0;
  for(Int_t i=1;i<=hDen1D->GetNbinsX();i++) {
    Q += hDen1D->GetBinContent(i);
  }
  Double_t xbinsize = hDen1D->GetBinWidth(0);
  Q *= xbinsize;
  // -----

  if(!CYL && !ThreeD) {
    Double_t sigmar = pData->GetBeamRmsY() * kp;
    Q *= TMath::Sqrt(2*TMath::Pi()) * sigmar;
  } else if(CYL) 
    Q *= 2*TMath::Pi();
  
  Q *= n0 / (kp * kp * kp);
  Q *= (PConst::ElectronCharge/PUnits::picocoulomb);
  
  cout << Form(" Integrated charge method 1 = %8.4f pC", Q) << endl;
  
  Q = 0;
  Double_t radius = 1;
  for(Int_t i=0;i<hDen2D->GetNbinsX();i++) {
    for(Int_t j=FirstyBin;j<=LastyBin;j++) {
      Double_t value  = hDen2D->GetBinContent(i,j);
      if(CYL) radius = hDen2D->GetYaxis()->GetBinCenter(j);
      Q += radius * value;
      // cout << Form(" (%i,%i) -> radius = %7.4f , value = %7.4f",i,j,radius,value) << endl;
    }    
  }
  
  xbinsize = hDen2D->GetXaxis()->GetBinWidth(0);
  Double_t ybinsize = hDen2D->GetYaxis()->GetBinWidth(0); 
  Q *= xbinsize * ybinsize;
  // -----

  if(!CYL && !ThreeD) {
    Double_t sigmar = pData->GetBeamRmsY() * kp;
    Q *= TMath::Sqrt(2*TMath::Pi()) * sigmar; 
  } else if(CYL) 
    Q *= 2*TMath::Pi();

  Q *= n0 / (kp * kp * kp);
  Q *= (PConst::ElectronCharge/PUnits::picocoulomb);
  
  cout << Form(" Integrated charge method 2 = %8.4f pC", Q) << endl;

  // ------------------------------------------------------------------------------------
  
  // Charge RMS calculation
  
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
      Double_t rms = 0;
      Double_t total = 0;
      for(Int_t k=1;k<=hDen2D->GetNbinsY();k++) {
	Double_t value  = hDen2D->GetBinContent(j,k);
	Double_t radius = hDen2D->GetYaxis()->GetBinCenter(k);
	if(CYL) {
	  rms += radius*radius*radius*value;
	  total += radius*value;
	} else {
	  rms += radius*radius*value;
	  total += value;
	}
	// cout << Form(" (%i,%i) -> radius = %7.4f  charge = %7.4f",j,k,radius,value) << endl;
      }

      if(total!=0.0)
	rms /= total;
      else
	rms = 0.;
      rms = sqrt(fabs(rms));

      // cout << Form(" j = %i  rms = %7.4f",j,rms) << endl;
      
      hRms->SetBinContent(j,rms); 
      
    }
    
    hRms->GetXaxis()->SetTitle("z [c/#omega_{p}]");
    if(opt.Contains("comov"))
      hRms->GetXaxis()->SetTitle("#zeta [c/#omega_{p}]");
    
    hRms->GetYaxis()->SetTitle("#LT r_{b} #GT_{rms} [c/#omega_{p}]");
    
  }
  
  // Tunning the Histograms
  // ---------------------

  // Chaning to user units: 

  cout << Form(" n0 = %10e ", n0 * PUnits::cm3) << endl;
  
  if(opt.Contains("units") && n0) {
    
    Int_t NbinsX = hDen2D->GetNbinsX();
    Double_t xMin = skindepth * hDen2D->GetXaxis()->GetXmin() / PUnits::mm;
    Double_t xMax = skindepth * hDen2D->GetXaxis()->GetXmax() / PUnits::mm;
    Int_t NbinsY = hDen2D->GetNbinsY();
    Double_t yMin = skindepth * hDen2D->GetYaxis()->GetXmin() / PUnits::mm;
    Double_t yMax = skindepth * hDen2D->GetYaxis()->GetXmax() / PUnits::mm;
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
  
    
    {
      Int_t NbinsX = hDen1D->GetNbinsX();
      Double_t xMin = skindepth * hDen1D->GetXaxis()->GetXmin() / PUnits::mm;
      Double_t xMax = skindepth * hDen1D->GetXaxis()->GetXmax() / PUnits::mm;
      hDen1D->SetBins(NbinsX,xMin,xMax);
      for(Int_t j=0;j<hDen1D->GetNbinsX();j++) {
	Double_t bincontent = (hDen1D->GetBinContent(j) * n0 / (1e15/PUnits::cm3)) * (skindepth/PUnits::cm);
	// if(CYL) bincontent *= skindepth/PUnits::cm; // yeeeeehaaaa
	hDen1D->SetBinContent(j,bincontent);
      }
      
      hDen1D->GetYaxis()->SetTitle("N_{e} [10^{15}/cm^{3}]");
      
      if(opt.Contains("comov"))
	hDen1D->GetXaxis()->SetTitle("#zeta [mm]");
      else
	hDen1D->GetXaxis()->SetTitle("z [mm]");
    }

    
    {
      Int_t NbinsX = hRms->GetNbinsX();
      Double_t xMin = skindepth * hRms->GetXaxis()->GetXmin() / PUnits::mm;
      Double_t xMax = skindepth * hRms->GetXaxis()->GetXmax() / PUnits::mm;
      hRms->SetBins(NbinsX,xMin,xMax);
      for(Int_t j=0;j<hRms->GetNbinsX();j++) {
	Double_t bincontent = (hRms->GetBinContent(j) * skindepth / PUnits::um);
	// if(CYL) bincontent *= skindepth/PUnits::cm; // yeeeeehaaaa
	hRms->SetBinContent(j,bincontent);
      }
      
      hRms->GetYaxis()->SetTitle("#LT r_{b} #GT_{rms} [#mum]");
      
      if(opt.Contains("comov"))
	hRms->GetXaxis()->SetTitle("#zeta [mm]");
      else
	hRms->GetXaxis()->SetTitle("z [mm]");
    }
    

    for(Int_t i=0;i<Nfields;i++) {
      Int_t NbinsX = hE1D[i]->GetNbinsX();
      Double_t xMin = skindepth * hE1D[i]->GetXaxis()->GetXmin() / PUnits::mm;
      Double_t xMax = skindepth * hE1D[i]->GetXaxis()->GetXmax() / PUnits::mm;
      hE1D[i]->SetBins(NbinsX,xMin,xMax);
      
      for(Int_t j=0;j<hE1D[i]->GetNbinsX();j++) {
	hE1D[i]->SetBinContent(j, hE1D[i]->GetBinContent(j) * ( E0 / (PUnits::GV/PUnits::m) ) * (skindepth/PUnits::m));
      }
      
      hE1D[i]->GetYaxis()->SetTitle("E [GV]");
      
      if(opt.Contains("comov"))
	hE1D[i]->GetXaxis()->SetTitle("#zeta [mm]");
      else
	hE1D[i]->GetXaxis()->SetTitle("z [mm]");
    }
  }
  
  // Plotting
  // -----------------------------------------------

  // Trick for nicer plotting
  // It uses hDen1D to setup the valid range
  Float_t bin1 = -1;
  Float_t bin2 = -1;
  for(Int_t i=0;i<hRms->GetNbinsX();i++) 
    if(bin1==-1 && hDen1D->GetBinContent(i+1)>0.0005) bin1 = hRms->GetBinCenter(i+1);  
  for(Int_t i=hRms->GetNbinsX()-1;i>=0;i--) 
    if(bin2==-1 && hDen1D->GetBinContent(i+1)>0.0005) bin2 = hRms->GetBinCenter(i+1);
  hRms->GetXaxis()->SetRangeUser(bin1,bin2);
  hRms->SetLineWidth(1);

  // Canvas setup
  TCanvas *C = new TCanvas("C","Charge density and Electric field",850,500);

  TPaveText *textTime = new TPaveText(0.63,0.87,0.78,0.92,"NDC");
  PlasmaGlob::SetPaveTextStyle(textTime); 
  textTime->SetTextColor(kGray+2);
  char ctext[128];

  if(opt.Contains("units") && n0) 
    sprintf(ctext,"Z = %5.1f mm", 1e3 * skindepth * Time);
  else
    sprintf(ctext,"T = %5.1f 1/#omega_{p}",Time);
  textTime->AddText(ctext);

  
  TPaveText *textRange = new TPaveText(0.13,0.85,0.38,0.92,"NDC");
  PlasmaGlob::SetPaveTextStyle(textRange,12); 
  textRange->SetTextColor(kGray+2);
  sprintf(ctext,"%5.3f < y < %5.3f",ymin,ymax);
  textRange->AddText(ctext);

  // Actual Plotting!
  // ------------------------------------------------------------

  // Output file
  TString fOutName = Form("./%s/Plots/Charge1D/Charge1D",sim.Data());
  fOutName += Form("-%s_%i",sim.Data(),time);

  // Colors
  Int_t plasmaC = kGray+1;
  Int_t beamC   = kAzure-5;
  Int_t fieldC  = kOrange+10;
 
  C->cd(0);
  gPad->SetFrameLineWidth(2);  

  Double_t Min = hDen1D->GetMinimum();
  Double_t Max = hDen1D->GetMaximum();
  hDen1D->GetYaxis()->SetRangeUser(Min,1.2*Max);

  hDen1D->SetLineColor(plasmaC);
  hDen1D->SetLineWidth(1);
  hDen1D->Draw("C");
  hDen1D->GetYaxis()->CenterTitle();
  hDen1D->GetXaxis()->CenterTitle();
  
  C->Update();

  TLine *line0 = new TLine(hDen1D->GetXaxis()->GetXmin(),
			   (gPad->GetUymin()+gPad->GetUymax())/2.,
			   hDen1D->GetXaxis()->GetXmax(),
			   (gPad->GetUymin()+gPad->GetUymax())/2.);
  line0->SetLineColor(kGray+1);
  line0->SetLineStyle(2);
  line0->Draw();


  Double_t rightmax = 2.5 * hRms->GetMaximum();
  Double_t slope = (gPad->GetUymax() - gPad->GetUymin())/rightmax;

  for(Int_t i=0;i<hRms->GetNbinsX();i++) {
    hRms->SetBinContent(i+1,hRms->GetBinContent(i+1)*slope + Min);
  }

  hRms->SetLineWidth(2);
  hRms->SetLineColor(beamC);
  hRms->Draw("same L");
  // hTest->Draw("same");
 
  //draw an axis on the right side
  TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),gPad->GetUxmax(),
  			    gPad->GetUymax(),0,rightmax,505,"+L");
  
  axis->SetLineWidth(1);
  axis->SetLineColor(beamC);
  axis->SetLabelColor(beamC);
  if(opt.Contains("units") && n0) 
    axis->SetTitle("#LT r_{b} #GT_{rms} [#mum]");
  else
    axis->SetTitle("#LT r_{b} #GT_{rms} [c/#omega_{p}]");
  axis->CenterTitle();
  axis->SetTitleColor(beamC);
  axis->SetTitleOffset(1.2);
  
  axis->Draw();

  // Longitudinal Electric field
  Double_t factor = 1.5;
  Double_t rightmin = factor * hE1D[0]->GetMinimum();
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
  Double_t rightmargin = 0.08;
  Double_t ux = gPad->PixeltoX(gPad->UtoPixel(1-rightmargin));
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
