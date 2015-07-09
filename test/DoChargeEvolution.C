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
#include "PlasmaGlob.hh"
#include "PPalette.hh"

void DoChargeEvolution( const TString &sim, Int_t time, Int_t index = 0, Int_t Nbins=1, const TString &options="") { 
  
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

  if(opt.Contains("center")) {
    Time -= zStartPlasma;
    if(opt.Contains("comov"))      // Centers on the head of the beam.
      Time += zStartBeam;
  } 


  // Integration range in z:
  Float_t zMin = -4.5;    // in norm units
  //Float_t zMax = -3.5;
  Float_t zMax = 10;
  opt += "comovcenter";   // Histo must be center for the validity of these limits.

  // Get charge density histos
  Int_t Nspecies = pData->NSpecies();
  TH2F **hDen2D = new TH2F*[Nspecies];
  TH2F **hDen2DRaw = new TH2F*[Nspecies];
  for(Int_t i=0;i<Nspecies;i++) {
    hDen2D[i] = NULL;

    if(i!=index) continue;

    if(!pData->GetChargeFileName(i) || !pData->GetRawFileName(i)) 
      continue;
    
    char hName[24];
    sprintf(hName,"hDen2D_%i",i);
    hDen2D[i] = (TH2F*) gROOT->FindObject(hName);
    if(hDen2D[i]) delete hDen2D[i];

    if(!ThreeD)
      hDen2D[i] = pData->GetCharge(i,opt);
    else // In 3D case, integrates X slice in between +-Nbins from the z axis:
      hDen2D[i] = pData->GetCharge2DSliceZY(i,-1,Nbins,opt+"int");
    
    hDen2D[i]->SetName(hName);
    hDen2D[i]->GetXaxis()->CenterTitle();
    hDen2D[i]->GetYaxis()->CenterTitle();
    hDen2D[i]->GetZaxis()->CenterTitle();
    
    if(opt.Contains("comov"))
      hDen2D[i]->GetXaxis()->SetTitle("#zeta [c/#omega_{p}]");
    else
      hDen2D[i]->GetXaxis()->SetTitle("z [c/#omega_{p}]");
    
    if(CYL) 
      hDen2D[i]->GetYaxis()->SetTitle("r [c/#omega_{p}]");
    else
      hDen2D[i]->GetYaxis()->SetTitle("y [c/#omega_{p}]");

    if(i==0)
      hDen2D[i]->GetZaxis()->SetTitle("n_{e} [n_{0}]");
    else
      hDen2D[i]->GetZaxis()->SetTitle("n_{b} [n_{0}]");


    if(opt.Contains("raw")) {
      sprintf(hName,"hDen2DRaw_%i",i);
      hDen2DRaw[i] = (TH2F*) gROOT->FindObject(hName);
      if(hDen2DRaw[i]) delete hDen2DRaw[i];
      
      // Define it like the charge density histogram
      Int_t Nx =  hDen2D[i]->GetXaxis()->GetNbins();
      Float_t xMin = hDen2D[i]->GetXaxis()->GetXmin();
      Float_t xMax = hDen2D[i]->GetXaxis()->GetXmax();
      Int_t Ny =  hDen2D[i]->GetYaxis()->GetNbins();
      Float_t yMin = hDen2D[i]->GetYaxis()->GetXmin();
      Float_t yMax = hDen2D[i]->GetYaxis()->GetXmax();
      hDen2DRaw[i] = new TH2F(hName,"",Nx,xMin,xMax,Ny,yMin,yMax);
      
      // TTree *treeRaw = pData->GetTreeRaw(pData->GetRawFileName(i)->c_str(),opt);
      // treeRaw->Project(hName,"x2:x1","ene>1");
      pData->GetH2Raw(pData->GetRawFileName(i)->c_str(),"x1","x2",hDen2DRaw[i],opt);
      
    }
  }

  // Find the range in z:
  Int_t FirstzBin = hDen2D[index]->GetXaxis()->FindBin(zMin); 
  Int_t LastzBin  = hDen2D[index]->GetXaxis()->FindBin(zMax); 

  // Calculate the "axis range" in number of bins. If Nbins==0 a RMS width is taken.
  Double_t rms0 = pData->GetBeamRmsY() * kp;
  if(CYL)  rms0  = pData->GetBeamRmsR() * kp;
  if(opt.Contains("units"))
    rms0 *= skindepth / PUnits::mm;
  
  Int_t FirstyBin = 0;
  Int_t LastyBin = 0;
  if(Nbins==0) { 
    Nbins =  TMath::Nint(rms0 / hDen2D[index]->GetYaxis()->GetBinWidth(1));
  }
  
  // Slice width limits.
  if(!CYL) {
    FirstyBin = hDen2D[index]->GetNbinsY()/2 + 1 - Nbins;
    LastyBin =  hDen2D[index]->GetNbinsY()/2 + Nbins;
  } else {
    FirstyBin = 1; 
    LastyBin  = Nbins;
  }  
  Float_t yMin = hDen2D[index]->GetYaxis()->GetBinLowEdge(FirstyBin); 
  Float_t yMax = hDen2D[index]->GetYaxis()->GetBinUpEdge(LastyBin); 

  
  Int_t yNbin = hDen2D[index]->GetNbinsY();
  Int_t zNbin = hDen2D[index]->GetNbinsX();

  // Calculate the charge:
  Double_t charge = 0.0;
  // INTEGRATED Beam's Charge:
  for(Int_t j=FirstzBin;j<=LastzBin;j++) {
    for(Int_t k=FirstyBin;k<=LastyBin;k++) {
      Double_t value  = hDen2D[index]->GetBinContent(j,k);
      if(CYL) {
	Double_t radius = hDen2D[index]->GetYaxis()->GetBinCenter(k);
	charge += radius * value;
	// cout << Form(" (%i,%i) -> radius = %7.4f , value = %7.4f",j,k,radius,value) << endl;
      } else {
	charge += value;
	// cout << Form(" (%i,%i) -> value = %7.4f",j,k,value) << endl;	
      }
    }    
  }
  Double_t xbinsize = hDen2D[index]->GetXaxis()->GetBinWidth(1);
  Double_t ybinsize = hDen2D[index]->GetYaxis()->GetBinWidth(1); 
  charge *= xbinsize * ybinsize;
  
  if(!CYL && !ThreeD) {
    charge *= TMath::Sqrt(2*TMath::Pi()) * rms0; 
  } else if(CYL) {
    charge *= 2*TMath::Pi();
  }
  
  if(opt.Contains("units")) {
  
    Time *= skindepth / PUnits::mm;

    Double_t dV = skindepth * skindepth * skindepth;
    charge *= n0 * dV;
    charge *= (PConst::ElectronCharge/PUnits::picocoulomb); 
    cout << Form(" Integrated charge of specie %3i = %8f pC",index,charge) << endl;
  } else {
    cout << Form(" Integrated charge of specie %3i = %8.4f n0 * kp^-3",index,charge) << endl;
  }

  // Calculate the charge:
  Double_t chargeRaw = 0.0;
  if(opt.Contains("raw")) {
    
    // INTEGRATED Beam's Charge:
    for(Int_t j=FirstzBin;j<=LastzBin;j++) {
      for(Int_t k=FirstyBin;k<=LastyBin;k++) {
	Double_t value  = hDen2DRaw[index]->GetBinContent(j,k);
	if(CYL) {
	  Double_t radius = hDen2DRaw[index]->GetYaxis()->GetBinCenter(k);
	  chargeRaw += radius * value;
	  // cout << Form(" (%i,%i) -> radius = %7.4f , value = %7.4f",j,k,radius,value) << endl;
	} else {
	  chargeRaw += value;
	  // cout << Form(" (%i,%i) -> value = %7.4f",j,k,value) << endl;	
	}
      }    
    }

    chargeRaw *=  xbinsize * ybinsize;
    if(!CYL && !ThreeD) {
      chargeRaw *= TMath::Sqrt(2*TMath::Pi()) * rms0; 
    } else if(CYL) {
      chargeRaw *= 2*TMath::Pi();
    } else {
      chargeRaw *= ybinsize;
    }

    //cout << Form("%f  %f",xbinsize,ybinsize);

    if(opt.Contains("units")) {
      Double_t dV = skindepth * skindepth * skindepth;
      chargeRaw *= n0 * dV;
      chargeRaw *= (PConst::ElectronCharge/PUnits::picocoulomb); 
      cout << Form(" Integrated charge (RAW) of specie %3i = %8f pC",index,chargeRaw) << endl;
    } else {
      cout << Form(" Integrated charge (RAW) of specie %3i = %8.4f n0 * kp^-3",index,chargeRaw) << endl;

      cout << Form(" MAGIC FACTOR = %.5f ",chargeRaw/charge) << endl; 
    }
  }
  
  // ------------------------------------------------------------------------------------
  

  // OUTPUT ROOT FILE WITH THE PLOTS:
  TString filename = Form("./%s/Plots/ChargeEvolution/ChargeEvolution-%s.root",sim.Data(),sim.Data());
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
  gChargevsTime->SetPoint(nPoints,Time,fabs(charge));
  gChargevsTime->Write(gName,TObject::kOverwrite);
  
  ifile->Close();
  
  // ------------------------------------------------------------------------------------
  
  // Output file
  TString fOutName = Form("./%s/Plots/ChargeEvolution/Charge",pData->GetPath().c_str());
  fOutName += Form("-%s_%i",pData->GetName(),time);

  // Canvas setup
  TCanvas *C = new TCanvas("C","Charge density and Electric field",850,500);
  C->cd(0);
  
  TExec *exElec   = new TExec("exElec","electronPalette->cd();");
  
  exElec->Draw();
  hDen2D[index]->Draw("colz");


  TLine *line[4];
  line[0] = new TLine(zMin,hDen2D[index]->GetYaxis()->GetXmin(),zMin,hDen2D[index]->GetYaxis()->GetXmax());
  line[1] = new TLine(zMax,hDen2D[index]->GetYaxis()->GetXmin(),zMax,hDen2D[index]->GetYaxis()->GetXmax());
  line[2] = new TLine(hDen2D[index]->GetXaxis()->GetXmin(),yMin,hDen2D[index]->GetXaxis()->GetXmax(),yMin);
  line[3] = new TLine(hDen2D[index]->GetXaxis()->GetXmin(),yMax,hDen2D[index]->GetXaxis()->GetXmax(),yMax);
  
  for(Int_t i=0;i<4;i++) {
    line[i]->SetLineColor(kGray+1);
    line[i]->SetLineStyle(2);
    line[i]->Draw();
  }
  
  C->cd();
  
  // Print to a file
  PlasmaGlob::imgconv(C,fOutName,opt);
  
  
}

