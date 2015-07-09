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
#include <TF1.h>
#include <TGraph.h>
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
#include <TFitResultPtr.h>
#include <TPaletteAxis.h>

#include "PData.hh"
#include "PlasmaGlob.hh"

void PlotEnergyModulation( const TString &sim, Int_t time, Double_t Emin, Double_t Emax, const TString &opt="") { 

  if(!opt.Contains("comov")) {
    cout << "PlotEnergyModulation:: This module needs \"comov\" option to work properly" << endl;
    return;
  }

  
  PlasmaGlob::Initialize();
  
  // Palettes!
  gROOT->Macro("PlasmaPalettes.C");

  // Init Units table
  PUnits::UnitsTable::Get();
  
  gStyle->SetPadTopMargin(0.05); 
  gStyle->SetPadLeftMargin(0.12);  // Margin left axis  
  gStyle->SetPadRightMargin(0.17); // Margin right axis 
  gStyle->SetGridStyle(2);
  gStyle->SetGridColor(kGray+1);
  gStyle->SetPadGridY(1);
  gStyle->SetPadGridX(0);

  // Load PData
  PData *pData = PData::Get(sim.Data());
  pData->LoadFileNames(time);
  if(!pData->IsInit()) return;
  
  Bool_t ThreeD = kFALSE;
  if(sim.Contains("3D")) ThreeD = kTRUE; 
  
  // Some plasma constants
  Double_t n0 = pData->GetPlasmaDensity();
  Double_t kp = pData->GetPlasmaK();
  Double_t skindepth = 1/kp;
  
  // Some beam properties:
  Double_t Ebeam = pData->GetBeamEnergy();
  Double_t gamma = Ebeam / PConst::ElectronMassE;
  Double_t vbeam = TMath::Sqrt(1 - 1/(gamma*gamma));
  // cout << Form(" - Bunch gamma      = %8.4f", gamma ) << endl;
  // cout << Form(" - Bunch velocity   = %8.4f c", vbeam ) << endl;

  // Time in OU
  Double_t Time = pData->GetRealTime();

  // z start of the plasma in normalized units.
  Double_t zStartPlasma = pData->GetPlasmaStart()*kp;

  // z start of the beam in normalized units.
  Double_t zStartBeam = pData->GetBeamStart()*kp;
  Time -= zStartPlasma - zStartBeam;
  
  // Energy and time resolutions
  // Double_t Eres = 0.01 * PUnits::Megaelectronvolt; 
  // Double_t Tres = 100 * PUnits::femtosecond; 
  Double_t Eres = 10.0 * PUnits::Kiloelectronvolt; 
  Double_t Tres = 0.10 * PUnits::picosecond; 
  Double_t Zres = Tres * PConst::c_light;
  // cout << PUnits::BestUnit(Zres, "Length") << endl;
  if(n0) {
    Zres = Zres/(skindepth*PUnits::meter);
    // cout << Zres << " norm units" << endl;
  }

  // Z axis limits and number of bins
  Int_t NbinZ; NbinZ = 0;
  Double_t ZMin,ZMax;
  ZMin = ZMax = 0.;
  
  // Histos for bunch specie
  TH2F *hDen2D     = NULL;   // Charge density
  TH1F *hEnergy    = NULL;   // Energy spectrum
  TH2F *hEvsZ     = NULL;   // Energy vs x1
  TH2F *hEvsZRes  = NULL;   // Energy vs x1
  Int_t Nspecies   = pData->NSpecies();
  for(Int_t i=0;i<Nspecies;i++) {
    
    if(!pData->GetChargeFileName(i)) continue;
      
    // Get 2D charge density histogram of the bunch.
    // We use it to calculate the rms of the charge distribution
    // in x2 as a function of x1.
    if(pData->GetSpeciesName(i).find("beam")==string::npos) 
      continue;
     
    char hName[24];
    sprintf(hName,"hDen2D");
    hDen2D = (TH2F*) gROOT->FindObject(hName);
    if(hDen2D) delete hDen2D;
    
    if(!ThreeD) 
      hDen2D = pData->GetCharge(i,opt);
    else
      hDen2D = pData->GetCharge2DSliceZY(i,-1,1,opt+"avg");
    
    hDen2D->SetName(hName);
    hDen2D->GetXaxis()->CenterTitle();
    hDen2D->GetYaxis()->CenterTitle();
    hDen2D->GetXaxis()->SetTitle("z [c/#omega_{p}]");
    hDen2D->GetYaxis()->SetTitle("y [c/#omega_{p}]");
    hDen2D->GetZaxis()->SetTitle("#LTn_{b}#GT [n_{0}]");

    // Get boundaries of the Z axis
    NbinZ = hDen2D->GetXaxis()->GetNbins();
    ZMin  = hDen2D->GetXaxis()->GetXmin();
    ZMax  = hDen2D->GetXaxis()->GetXmax();
    
    // Open z range for left spectrum
    ZMin  = ZMin - 0.2*(ZMax-ZMin);

    // Get spectrum
    char hName2[24];
    sprintf(hName2,"hEnergy");
    hEnergy = (TH1F*) gROOT->FindObject(hName2);
    if(hEnergy) { delete hEnergy; hEnergy = NULL; }

    Int_t NEbin = int((Emax - Emin) * PUnits::Megaelectronvolt/Eres);
    hEnergy = new TH1F(hName2,"",NEbin,Emin,Emax);
    pData->GetH1Raw(pData->GetRawFileName(i)->c_str(),"energy",hEnergy,opt);
    
    hEnergy->GetXaxis()->SetTitle("E [MeV]");
    hEnergy->GetYaxis()->SetTitle("dQ/dE [a.u.]");
    PlasmaGlob::SetH1Style(hEnergy,i);
    
    // Get Spectrum vs x1
    sprintf(hName2,"hEvsZ");
    hEvsZ = (TH2F*) gROOT->FindObject(hName2);
    if(hEvsZ) { delete hEvsZ; hEvsZ = NULL; }

    Int_t NZbin = int((ZMax - ZMin)/Zres);
    // hEvsZ = new TH2F(hName2,"",NbinZ,ZMin,ZMax,NEbin,Emin,Emax);
    hEvsZ = new TH2F(hName2,"",NZbin,ZMin,ZMax,NEbin,Emin,Emax);
    
    pData->GetH2Raw(pData->GetRawFileName(i)->c_str(),"x1","energy",hEvsZ,opt);

    hEvsZ->GetYaxis()->SetTitle("Energy [MeV]");
    if(opt.Contains("comov")) {
      hEvsZ->GetXaxis()->SetTitle("#zeta [c/#omega_{p}]");
      hEvsZ->GetZaxis()->SetTitle("dn_{b}/d#zetadE [a.u.]"); 
    } else {
      hEvsZ->GetXaxis()->SetTitle("z [c/#omega_{p}]");
      hEvsZ->GetZaxis()->SetTitle("dn_{b}/dzdE [a.u.]"); 
    }
   
    hEvsZ->GetZaxis()->CenterTitle();

    PlasmaGlob::SetH1Style(hEvsZ,i);
  }
  
  
  TProfile *hEvsZprof = NULL;
  TGraphErrors *gEvsZ = NULL;
  if(hEvsZ) {
    TString pname = hEvsZ->GetName();
    pname += "_pfx";
    hEvsZprof =  (TProfile*) gROOT->FindObject(pname.Data());
    if(hEvsZprof) delete hEvsZprof;

    hEvsZprof = hEvsZ->ProfileX("_pfx",1,-1,"s");

    gEvsZ = (TGraphErrors*) gROOT->FindObject("gEvsZ");
    if(gEvsZ) delete gEvsZ;

    Int_t Npoints = hEvsZprof->GetNbinsX();
    Double_t *x = new Double_t[Npoints];
    Double_t *y = new Double_t[Npoints];
    Double_t *ex = new Double_t[Npoints];
    Double_t *ey = new Double_t[Npoints];
    
    for(Int_t j=0;j<Npoints;j++) {
      x[j] = hEvsZprof->GetBinCenter(j);
      y[j] = hEvsZprof->GetBinContent(j);
      ex[j] = 0;
      ey[j] = hEvsZprof->GetBinError(j);   
    }
    
    gEvsZ = new TGraphErrors(Npoints,x,y,ex,ey);
    gEvsZ->SetName("gEvsZ");
        
    // PlasmaGlob::SetH1Style((TH1*)gEvsZ,1);
    PlasmaGlob::SetGraphStyle(gEvsZ,1);

   
    if(opt.Contains("comov")) 
      gEvsZ->GetXaxis()->SetTitle("#zeta [c/#omega_{p}]");
    else
      gEvsZ->GetXaxis()->SetTitle("z [c/#omega_{p}]");
    
    gEvsZ->GetYaxis()->SetTitle("#LTE#GT [MeV]");
    
  }


  // Tunning the Histograms
  // ---------------------


  // Chaning to user units: 
  // --------------------------
  
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

    
    NbinsX = hEvsZ->GetNbinsX();
    xMin = skindepth * hEvsZ->GetXaxis()->GetXmin() / PUnits::mm;
    xMax = skindepth * hEvsZ->GetXaxis()->GetXmax() / PUnits::mm;
    NbinsY = hEvsZ->GetNbinsY();
    yMin = hEvsZ->GetYaxis()->GetXmin();
    yMax = hEvsZ->GetYaxis()->GetXmax();
    hEvsZ->SetBins(NbinsX,xMin,xMax,NbinsY,yMin,yMax);
    for(Int_t j=0;j<hEvsZ->GetNbinsX();j++) {
      for(Int_t k=0;k<hEvsZ->GetNbinsY();k++) {
	hEvsZ->SetBinContent(j,k, hEvsZ->GetBinContent(j,k) * n0 / (1e15/PUnits::cm3) );
      }
    }
    
    if(opt.Contains("comov")) {
      hEvsZ->GetXaxis()->SetTitle("#zeta [mm]");
      hEvsZ->GetZaxis()->SetTitle("dn_{b}/d#zetadE [a.u.]"); 
    } else {
      hEvsZ->GetXaxis()->SetTitle("z [mm]");
      hEvsZ->GetZaxis()->SetTitle("dn_{b}/dzdE [a.u.]"); 
    }

    // New boundaries after change of units
    ZMin  *= skindepth / PUnits::mm;
    ZMax  *= skindepth / PUnits::mm;

    
    Int_t Npoints = gEvsZ->GetN();
    Double_t *x = gEvsZ->GetX();
    Double_t *y = gEvsZ->GetY();
    Double_t *ex = gEvsZ->GetEX();
    Double_t *ey = gEvsZ->GetEY();

    for(Int_t j=0;j<Npoints;j++) {
      // cout << Form(" i = %i  ->  x = %5f  y = %5f",j,x[j],y[j]) << endl;
      gEvsZ->SetPoint(j,x[j] * skindepth / PUnits::mm,y[j]);
      gEvsZ->SetPointError(j,ex[j] * skindepth / PUnits::mm,ey[j]);   
    }
    
    if(opt.Contains("comov"))
      gEvsZ->GetXaxis()->SetTitle("#zeta [mm]");
    else
      gEvsZ->GetXaxis()->SetTitle("z [mm]");
  }
  

  // Vertical graphs: Displayed at a side of 2D histograms
  char gName[24];
  sprintf(gName,"gEnergy");
  TGraph *gEnergy = (TGraph*) gROOT->FindObject(gName);
  if(gEnergy) { delete gEnergy; gEnergy = NULL; } 
  
  if(hEnergy) {
    Int_t Nbin   = hEnergy->GetNbinsX();
    Double_t *y   = new Double_t[Nbin];
    Double_t *x   = new Double_t[Nbin];
    
    Double_t xMin = ZMin;
    Double_t xMax = ZMin + (ZMax-ZMin) * 0.2;
    Double_t Emax = hEnergy->GetMaximum();
    for(Int_t i=0; i<Nbin; i++) {
      y[i] = hEnergy->GetBinCenter(i+1);
      x[i] = ((xMax-xMin)/Emax)*hEnergy->GetBinContent(i+1) + xMin;
    }

    gEnergy = new TGraph(Nbin,x,y);
    gEnergy->SetName(gName);
    PlasmaGlob::SetGraphStyle(gEnergy,1);
    // gEnergy->SetLineColor(PlasmaGlob::elecLine);
    // gEnergy->SetLineWidth(2);
    // gEnergy->SetFillStyle(1001);
    // gEnergy->SetFillColor(PlasmaGlob::elecFill);

    delete x;
    delete y;
  }
  
  
  // Fitting...
  // Range of fit in comoving coordinate (in mm).
  // It takes the range of the totally flat-top part.
  Double_t x1 = -5.5;
  Double_t x2 = -0.0;

  if( !opt.Contains("units") && n0) {
    x1 *= kp * PUnits::mm;
    x2 *= kp * PUnits::mm;
  }
  
  char fitName[64];
  sprintf(fitName,"fitFunction_%i",time);
  TF1 *fitFunction = (TF1*) gROOT->FindObject(fitName);
  if(!fitFunction) 
    fitFunction = new TF1(fitName,"[0]+[1]*sin([2]*(x-[3]))",x1,x2);
  
  Int_t Npoints = gEvsZ->GetN();
  Double_t *y = gEvsZ->GetY();
  Double_t EmeanMax = -99;
  for(Int_t j=0;j<Npoints;j++) {
    if(y[j]>EmeanMax) EmeanMax = y[j];
  }
  Double_t E0 = 25.;
  Double_t EmeanValue = hEnergy->GetMean();
  cout << Form(" Maximum energy = %5f MeV",EmeanMax) << endl;
  Double_t skdepth = 1.;
  
  if(opt.Contains("units") && n0) 
    skdepth = skindepth / PUnits::mm;    // in mm.
  
  Double_t b = x2 - skdepth * TMath::ASin(-1); 
  
  fitFunction->SetParameters(E0,EmeanMax-EmeanValue,1./skdepth,b);
  fitFunction->SetParLimits(1,0.,2.);
  
  Int_t res = gEvsZ->Fit(fitName,"R0EX0");
     
  // Retrieve parameters of the fit
  Double_t chi2 = fitFunction->GetChisquare();
  Double_t NDF = fitFunction->GetNDF();
  Double_t *par = fitFunction->GetParameters();
  Double_t *parErr = fitFunction->GetParErrors();
  
  // cout << " x1 = " << x1 << "  x2 = " << x2 << "  k = " << k << endl;
  if(res!=0) {
    cout << "Fit results: " << endl;
    cout << Form("chi2 = %5.2f : A = %5.2f  B = %5.2f  C = %5.2f  D  = %5.2f",chi2,par[0],par[1],par[2],par[3]) << endl;
  }

  // Calculate weighted energy gain:
  Double_t Egain = 0.;
  Double_t Eloss = 0.;
  if(hEnergy) {
    for(Int_t i=0;i<hEnergy->GetNbinsX();i++) {
      if(hEnergy->GetBinCenter(i+1) > E0) {
	Egain += hEnergy->GetBinCenter(i+1)*hEnergy->GetBinContent(i+1);
      } else {
	Eloss += hEnergy->GetBinCenter(i+1)*hEnergy->GetBinContent(i+1);
      }
    } 
  }
  Double_t Asym = (Egain-Eloss)/(Egain+Eloss);

  // Plots over time steps (looping)
  // Retrieve the previous time TGraph if any;
  // Open TGraph
  TString filename = Form("./%s/Plots/EnergyModulation/EnergyModulation-%s.root",sim.Data(),sim.Data());
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
  TGraphErrors *gEmean  = NULL;
  TGraphErrors *gEamp   = NULL;
  TGraphErrors *gLambda = NULL;
  TGraphErrors *gPhase = NULL;
  TGraphErrors *gChi2   = NULL;
  TGraphErrors *gAsym  = NULL;


  sprintf(gName,"gChi2");
  gChi2 = (TGraphErrors*) ifile->Get(gName);
  if(gChi2==NULL) {
    gChi2 = new TGraphErrors();
    gChi2->SetName(gName);  
  } else {
    nPoints = gChi2->GetN(); 
  }  

  Bool_t OK = kTRUE;
  if(nPoints>1) {
    Double_t chi2old;
    Double_t L;
    gChi2->GetPoint(nPoints-1,L,chi2old);
    if(fabs(chi2/NDF-chi2old)>20.) {
      OK = kFALSE;
      cout << " chi2old = " << chi2old << "  chi2 = " << chi2/NDF << endl;
    }
  }
  
  if(OK) {
    gChi2->Set(nPoints+1);
    if(opt.Contains("units") && n0)
      gChi2->SetPoint(nPoints,Time,chi2/NDF);  
    else
      gChi2->SetPoint(nPoints,Time,chi2/NDF);  
    
    gChi2->GetYaxis()->SetTitle("#chi^{2}/ndf");
    
    if(opt.Contains("units") && n0) {
      gChi2->GetXaxis()->SetTitle("Z [mm]");
    } else {
      gChi2->GetXaxis()->SetTitle("T [c/#omega_{p}]");
    }
    gChi2->Write(gChi2->GetName(),TObject::kOverwrite);
    
    sprintf(gName,"gEmean");
    gEmean = (TGraphErrors*) ifile->Get(gName);
    if(gEmean==NULL) {
      gEmean = new TGraphErrors();
      gEmean->SetName(gName);
    } else {
      nPoints = gEmean->GetN(); 
    }  
    gEmean->Set(nPoints+1);
    if(opt.Contains("units") && n0) {
      gEmean->SetPoint(nPoints,Time,par[0]);
    }
    else
      gEmean->SetPoint(nPoints,Time,par[0]);  
    
    gEmean->SetPointError(nPoints,0.,parErr[0]);
    
    gEmean->GetYaxis()->SetTitle("#LTE#GT [MeV]");
    if(opt.Contains("units") && n0) {
      gEmean->GetXaxis()->SetTitle("Z [mm]");
    } else {
      gEmean->GetXaxis()->SetTitle("T [c/#omega_{p}]");
    } 
    gEmean->Write(gEmean->GetName(),TObject::kOverwrite);
    
    sprintf(gName,"gEamp");
    gEamp = (TGraphErrors*) ifile->Get(gName);
    if(gEamp==NULL) {
      gEamp = new TGraphErrors();
      gEamp->SetName(gName);
    } else {
      nPoints = gEamp->GetN(); 
    }  
    gEamp->Set(nPoints+1);
    if(opt.Contains("units") && n0)
      gEamp->SetPoint(nPoints,Time,par[1]);  
    else
      gEamp->SetPoint(nPoints,Time,par[1]);  
    
    gEamp->SetPointError(nPoints,0,parErr[1]);
    
    gEamp->GetYaxis()->SetTitle("Amplitude [MeV]");
    
    if(opt.Contains("units") && n0) {
      gEamp->GetXaxis()->SetTitle("Z [mm]");
    } else {
      gEamp->GetXaxis()->SetTitle("T [c/#omega_{p}]");
    }
    
    gEamp->Write(gEamp->GetName(),TObject::kOverwrite);
    
    
    sprintf(gName,"gAsym");
    gAsym = (TGraphErrors*) ifile->Get(gName);
    if(gAsym==NULL) {
      gAsym = new TGraphErrors();
      gAsym->SetName(gName);   
    } else {
      nPoints = gAsym->GetN(); 
    }  
    gAsym->Set(nPoints+1);
    if(opt.Contains("units") && n0)
      gAsym->SetPoint(nPoints,Time,Asym);  
    else
      gAsym->SetPoint(nPoints,Time,Asym);  
    
    gAsym->GetYaxis()->SetTitle("Asymmetry");
    
    if(opt.Contains("units") && n0) {
      gAsym->GetXaxis()->SetTitle("Z [mm]");
    } else {
      gAsym->GetXaxis()->SetTitle("T [c/#omega_{p}]");
    }
    gAsym->Write(gAsym->GetName(),TObject::kOverwrite);
    
    
    sprintf(gName,"gLambda");
    gLambda = (TGraphErrors*) ifile->Get(gName);
    if(gLambda==NULL) {
      gLambda = new TGraphErrors();
      gLambda->SetName(gName);
    } else {
      nPoints = gLambda->GetN(); 
    }  
    gLambda->Set(nPoints+1);
    if(opt.Contains("units") && n0)
      gLambda->SetPoint(nPoints,Time,TMath::TwoPi()/par[2]);  
    else
      gLambda->SetPoint(nPoints,Time,TMath::TwoPi()/par[2]);  
    
    gLambda->SetPointError(nPoints,0,TMath::TwoPi()*parErr[2]/(par[2]*par[2]));
    
    if(opt.Contains("units") && n0) {
      gLambda->GetXaxis()->SetTitle("Z [mm]");
      gLambda->GetYaxis()->SetTitle("#lambda [mm]");
    } else {
      gLambda->GetXaxis()->SetTitle("T [c/#omega_{p}]");
      gLambda->GetYaxis()->SetTitle("#lambda [c/#omega_{p}]");
    }
    gLambda->Write(gLambda->GetName(),TObject::kOverwrite);
    

    sprintf(gName,"gPhase");
    gPhase = (TGraphErrors*) ifile->Get(gName);
    if(gPhase==NULL) {
      gPhase = new TGraphErrors();
      gPhase->SetName(gName);
    } else {
      nPoints = gPhase->GetN(); 
    }  
    gPhase->Set(nPoints+1);
    if(opt.Contains("units") && n0)
      gPhase->SetPoint(nPoints,Time,par[3]);  
    else
      gPhase->SetPoint(nPoints,Time,par[3]);  
    
    gPhase->SetPointError(nPoints,0,parErr[3]);
    
    if(opt.Contains("units") && n0) {
      gPhase->GetXaxis()->SetTitle("Z [mm]");
      gPhase->GetYaxis()->SetTitle("#psi [mm]");
    } else {
      gPhase->GetXaxis()->SetTitle("T [c/#omega_{p}]");
      gPhase->GetYaxis()->SetTitle("#psi [c/#omega_{p}]");
    }
    gPhase->Write(gPhase->GetName(),TObject::kOverwrite);
    
    fitFunction->Write(fitFunction->GetName(),TObject::kOverwrite);
  }
  ifile->Close();
  

  // Plotting
  // -----------------------------------------------

  // Output file
  TString fOutName = Form("./%s/Plots/EnergyModulation/EnergyModulation",sim.Data());
  fOutName += Form("-%s_%i",sim.Data(),time);

  // Canvas setup
  TCanvas *C;
  if(opt.Contains("hres") && !opt.Contains("pdf")) // high resolution for plain grahics output.
    C = new TCanvas("C","Energy Vs Z",1000,625);
  else
    C = new TCanvas("C","Energy Vs Z",800,500);
  
  // Palettes setup
  TExec *exElec = new TExec("exElec","electronPalette->cd();");

  // Draw objects
  TPaveText *textTime = new TPaveText(0.65,0.87,0.80,0.92,"NDC");
  PlasmaGlob::SetPaveTextStyle(textTime); 
  char ctext[128];
  if(opt.Contains("units") && n0) 
    sprintf(ctext,"z = %5.1f mm", 1e3 *skindepth * Time);
  else
    sprintf(ctext,"t = %5.1f 1/#omega_{p}",Time);
  textTime->SetTextFont(42);
  textTime->AddText(ctext);

  TPaveText *textChi = new TPaveText(0.15,0.87,0.30,0.92,"NDC");
  PlasmaGlob::SetPaveTextStyle(textChi,12); 
  textChi->SetTextColor(kMagenta+2);
  //textChi->SetTextColor(kGray+2);
  sprintf(ctext,"#chi^{2}/ndf = %5.2f",chi2/NDF);
  textChi->SetTextFont(42);
  textChi->AddText(ctext);

  
  // Colors

  // Actual Plotting!
  // ------------------------------------------------------------

  // More makeup
  C->cd(0);
  gPad->SetFrameLineWidth(2);  

  TH2F *hFrame = (TH2F*) hEvsZ->Clone("hFrame");
  hFrame->Reset();
  hFrame->Draw("col");
  
  hFrame->SetLabelFont(42,"xyz");
  hFrame->SetTitleFont(42,"xyz");

  hFrame->GetYaxis()->SetNdivisions(505);

  hFrame->GetXaxis()->SetNdivisions(510);
  hFrame->GetXaxis()->SetLabelSize(0.06);
  hFrame->GetXaxis()->SetTitleSize(0.065);
  hFrame->GetXaxis()->SetTitleOffset(1.0);

  hEvsZ->GetZaxis()->SetTitleFont(42);
  hEvsZ->GetZaxis()->SetLabelFont(42);
 
  exElec->Draw();
  if(hEvsZ) 
    hEvsZ->Draw("col same z");

  gPad->Update();
  
  Float_t y1 = gPad->GetBottomMargin();
  Float_t y2 = 1 - gPad->GetTopMargin();
  x1 = gPad->GetLeftMargin();
  x2 = 1 - gPad->GetRightMargin();
  
  TPaletteAxis *palette = (TPaletteAxis*) hEvsZ->GetListOfFunctions()->FindObject("palette");  
  if(palette) {
    palette->SetY2NDC(y2 - 0.00);
    palette->SetY1NDC(y1 + 0.00);
    palette->SetX1NDC(x2 + 0.005);
    palette->SetX2NDC(x2 + 0.03);
    // palette->SetLabelSize(yFactor*0.085);
    // palette->SetLabelFont(42);
    // palette->SetLabelOffset(0.01/yFactor);
    // // palette->SetTitleFont(42);
    // palette->SetTitleSize(yFactor*0.09);
    // palette->SetTitleOffset(0.7/yFactor);
    palette->SetBorderSize(2);
    palette->SetLineColor(1);
  }
   
  gPad->RedrawAxis(); 

  
  gEnergy->Draw("F");
  gEnergy->Draw("L");

  gEvsZ->SetMarkerSize(0.6);
  //gEvsZ->Draw("PZ same");
  
  textTime->Draw();
  textChi->Draw();

  gPad->RedrawAxis(); 
  gPad->RedrawAxis("G"); 

  fitFunction->SetLineStyle(2);
  fitFunction->SetLineWidth(1);
  fitFunction->SetLineColor(kMagenta+2);
  // fitFunction->SetLineColor(kGray+3);
  fitFunction->Draw("same");
 
  C->cd();

  // Print to a file
  PlasmaGlob::imgconv(C,fOutName,opt);
  // ---------------------------------------------------------
  
  delete fitFunction;
}
