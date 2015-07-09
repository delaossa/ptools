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

void DoField1DEvolution( const TString &sim, Int_t time, Int_t Nbins=1, const TString &options="") { 
  
#ifdef __CINT__  
  gSystem->Load("libplasma.so");
#endif

  PlasmaGlob::Initialize();

  // Palettes!
  gROOT->Macro("PlasmaPalettes.C");

  // Init Units table
  PUnits::UnitsTable::Get();
  
  // Load PData
  PData *pData = PData::Get(sim.Data());
  pData->LoadFileNames(time);
  if(!pData->IsInit()) return;

  TString opt = options;

  // Always comoving frame for this analysis
  opt += "comov";
 
  // More makeup
  gStyle->SetPadGridY(0);
  if(opt.Contains("gridx")) {
    gStyle->SetPadGridX(1);
  }
  if(opt.Contains("gridy")) {
    gStyle->SetPadGridY(1);
  }
 
  
  // Some plasma constants
  Float_t n0 = pData->GetPlasmaDensity();
  Float_t omegap = pData->GetPlasmaFrequency();
  Float_t timedepth = pData->GetPlasmaTimeDepth();
  Float_t kp = pData->GetPlasmaK();
  Float_t skindepth = pData->GetPlasmaSkinDepth();
  Float_t E0 = pData->GetPlasmaE0();
  
  // Some beam properties:
  Float_t Ebeam = pData->GetBeamEnergy();
  Float_t gamma = pData->GetBeamGamma();
  Float_t vbeam = pData->GetBeamVelocity();

  cout << Form(" - Bunch gamma      = %8.4f", gamma ) << endl;
  cout << Form(" - Bunch velocity   = %8.4f c", vbeam ) << endl;
  
  // Other parameters
  Float_t trapPotential = 1.0 - (1.0/gamma);
  cout << Form(" - Trap. potential  = %8.4f mc2/e",trapPotential) << endl;
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
  Float_t shiftz = pData->Shift(opt);

  
  // Calculate the "axis range" in number of bins. If Nbins==0 a RMS width is taken.
  Float_t rms0 = pData->GetBeamRmsY() * kp;
  if(pData->IsCyl())  rms0  = pData->GetBeamRmsR() * kp;
  
  Int_t FirstyBin = 0;
  Int_t LastyBin = 0;
  if(Nbins==0) { 
    if(rms0>0.0)
      Nbins =  TMath::Nint(rms0 / pData->GetDX(1));
    else
      Nbins = 1;
  }
  
  // Slice width limits.
  if(!pData->IsCyl()) {
    FirstyBin = pData->GetNX(1)/2 + 1 - Nbins;
    LastyBin =  pData->GetNX(1)/2 + Nbins;
  } else {
    FirstyBin = 1; 
    LastyBin  = Nbins;
  }


  // On-axis range

  Float_t zMin = pData->GetXMin(0);
  Float_t zMax = pData->GetXMax(0);
  //zMin = -9.0  + pData->Shift("centercomov");
  //zMax = -3.77 + pData->Shift("centercomov");
  Float_t zRange = zMax - zMin;  
  pData->SetX1Min(zMin);
  pData->SetX1Max(zMax);

  // ----------------------------------------------------------------------------------


  // OUTPUT ROOT FILE WITH THE PLOTS:
  TString filename = Form("./%s/Plots/Field1DEvolutions/Field1DEvolutions-%s.root",sim.Data(),sim.Data());
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

  
  // EM fields on - axis :
  // ------------------------------------------------------------------------

  const Int_t Nfields = 2;
  TH1F **hE1D = new TH1F*[Nfields];
  for(Int_t i=0;i<Nfields;i++) {
    hE1D[i] = NULL;

    if(!pData->GetEfieldFileName(i))
      continue;

    cout << Form(" Getting electric field number ") << i+1 << endl;
    
    char hName[24];
    sprintf(hName,"hE1D_%i",i);
    hE1D[i] = (TH1F*) gROOT->FindObject(hName);
    if(hE1D[i]) delete hE1D[i];
      
    // 1D histograms
    char nam[3]; sprintf(nam,"e%i",i+1);
    if(pData->Is3D()) {
      
      if(i==0) 
	hE1D[i] = pData->GetH1SliceZ3D(pData->GetEfieldFileName(i)->c_str(),nam,-1,Nbins,-1,Nbins,opt+"avg");
      else  
	hE1D[i] = pData->GetH1SliceZ3D(pData->GetEfieldFileName(i)->c_str(),nam,-Nbins,Nbins,-Nbins,Nbins,opt+"avg");
      
    } else if(pData->IsCyl()) { // Cylindrical: The first bin with r>0 is actually the number 1 (not the 0).
      
      hE1D[i] = pData->GetH1SliceZ(pData->GetEfieldFileName(i)->c_str(),nam,1,Nbins,opt+"avg");
      
    } else { // 2D cartesian
      
      if(i==0) 
	hE1D[i] = pData->GetH1SliceZ(pData->GetEfieldFileName(i)->c_str(),nam,-1,Nbins,opt+"avg");
      else 
	hE1D[i] = pData->GetH1SliceZ(pData->GetEfieldFileName(i)->c_str(),nam,-Nbins,Nbins,opt+"avg");    
    }
    
    hE1D[i]->SetName(hName);
    if(opt.Contains("comov"))
      hE1D[i]->GetXaxis()->SetTitle("#zeta [c/#omega_{p}]");
    else
      hE1D[i]->GetXaxis()->SetTitle("z [c/#omega_{p}]");
    
    if(i==0)
      hE1D[i]->GetYaxis()->SetTitle("E_{z} [E_{0}]");
    else if(i==1)
      hE1D[i]->GetYaxis()->SetTitle("E_{y} [E_{0}]");
    else if(i==2)
      hE1D[i]->GetYaxis()->SetTitle("E_{x} [E_{0}]");
    

  }  
  
  // Calculate wave positions:
  // ----------------------------------------------------------------
    
  // Calculate the crossings and the extremes of the Electric fields
  Float_t zEcross[Nfields][100] = {{0.0}};
  Float_t Eextr[Nfields][100] = {{0.0}};
  Float_t zEextr[Nfields][100] = {{0.0}};
  Int_t Ncross[Nfields] = {0};

  for(Int_t i=0;i<Nfields;i++) {
    Ncross[i] = 0;
    
    if(!hE1D[i]) continue;
    
    // Only smooths the focusing if flag activated..
    if(i>0 && opt.Contains("smooth")) {
      // cout << " Smoothing fields on axis..." << endl;
      hE1D[i]->Smooth(10);
    } 

    Float_t maxZeta = zStartBeam +5;
    if(opt.Contains("center")) 
      maxZeta -= zStartBeam;

    cout << "  Examining Field " << i << " : " << Ncross[i] << endl;
    //cout << "  Max Z =  " << maxZeta << endl;
          
    for(Int_t ip=hE1D[i]->GetNbinsX();ip>1;ip--) {

      Float_t Z2 = hE1D[i]->GetBinCenter(ip-1);
      if(Z2 > maxZeta) continue;
      Float_t E1 = hE1D[i]->GetBinContent(ip);
      Float_t E2 = hE1D[i]->GetBinContent(ip-1);
      Float_t Z1 = hE1D[i]->GetBinCenter(ip);
      
      //cout << Form("Z1 = %6.4f  Z2 = %6.4f   E1 = %6.4f   E2 = %6.4f", Z1, Z2, E1, E2) << endl; 

      if(E1*E2 >= 0) { // No change of sign means we are in a side of the zero axis.
	if(fabs(E2)>fabs(Eextr[i][Ncross[i]])) {
	  Eextr[i][Ncross[i]]  = E2;
	  zEextr[i][Ncross[i]] = Z2;
	  //cout << " --> Extreme found! " << endl;
	} 
      }
      
      if(E1*E2 < 0) { // change of sign means a crossing!
	
	// The next crossing has to be further away than a skindepth from the previous one:
	Float_t zcross =  -E1 * ( (Z2-Z1)/(E2-E1) ) + Z1;
        if(Ncross[i]>0 && fabs(zEcross[i][Ncross[i]-1]-zcross)<1 ) continue;	
	// cout << " CROSS! " << endl;

	// add the point
	zEcross[i][Ncross[i]] = zcross;
	Ncross[i]++;
      }
    }
    
    cout << "  -> Number of crossings for field " << i << " : " << Ncross[i] << endl;
    for(Int_t ic=0;ic<Ncross[i];ic++) {
      //  cout << Form(" %2i:  zeta = %6.4f  E = %6.4f", ic, zEcross[i][ic], Eextr[i][ic]) << endl; 
    }  
    
    
    hE1D[i]->SetLineColor(kRed);
    hE1D[i]->Write(hE1D[i]->GetName(),TObject::kOverwrite);

  }
  
  // Get the Graphs and histos from file
  Int_t nPoints = 0;
  TGraph ***gzEcross = new TGraph**[Nfields]; 
  TGraph ***gEextr  = new TGraph**[Nfields]; 
  TGraph ***gzEextr  = new TGraph**[Nfields]; 
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

    gzEcross[i] = new TGraph*[Ncross[i]];
    gEextr[i] = new TGraph*[Ncross[i]];
    gzEextr[i] = new TGraph*[Ncross[i]];
    char gName[24];
    for(Int_t ic=0;ic<Ncross[i];ic++) {
      sprintf(gName,"gzEcross_%i_%i",i,ic);     
      gzEcross[i][ic] = (TGraph*) ifile->Get(gName);
      if(gzEcross[i][ic]==NULL) {
	gzEcross[i][ic] = new TGraph();
	gzEcross[i][ic]->SetName(gName);
	nPoints = 0;
	// Some cosmetics at creation time:
	if(i==1) gzEcross[i][ic]->SetLineStyle(2);
	else gzEcross[i][ic]->SetLineStyle(1);
	gzEcross[i][ic]->SetLineWidth(1);
	gzEcross[i][ic]->SetLineColor(kGray+1);
	gzEcross[i][ic]->SetMarkerStyle(20);
	gzEcross[i][ic]->SetMarkerSize(0.4);
	gzEcross[i][ic]->SetMarkerColor(kGray+1);	
	gzEcross[i][ic]->GetYaxis()->SetTitle("k_{p}#zeta]");
	gzEcross[i][ic]->GetXaxis()->SetTitle("k_{p}z");
      } else {
	nPoints = gzEcross[i][ic]->GetN(); 
      }  
      
      gzEcross[i][ic]->Set(nPoints+1);
      gzEcross[i][ic]->SetPoint(nPoints,Time,zEcross[i][ic]);
      gzEcross[i][ic]->Write(gName,TObject::kOverwrite);
      
      // if(ic==Ncross[i]-1) continue;
      
      sprintf(gName,"gEextr_%i_%i",i,ic);     
      gEextr[i][ic] = (TGraph*) ifile->Get(gName);
      if(gEextr[i][ic]==NULL) {
	gEextr[i][ic] = new TGraph();
	gEextr[i][ic]->SetName(gName);
	nPoints = 0;
	// Some cosmetics at creation time:
	if(i==0) {
	  gEextr[i][ic]->SetLineWidth(2);
	  gEextr[i][ic]->SetLineColor(PlasmaGlob::fieldLine);
	  gEextr[i][ic]->SetMarkerStyle(16);
	  gEextr[i][ic]->SetMarkerSize(0.4);
	  gEextr[i][ic]->SetMarkerColor(PlasmaGlob::fieldLine);	
	  gEextr[i][ic]->GetYaxis()->SetTitle("E_{z}/E_{0}");
	  gEextr[i][ic]->GetXaxis()->SetTitle("k_{p}z");
	} else if(i==1) {
	  gEextr[i][ic]->SetLineWidth(2);
	  gEextr[i][ic]->SetLineColor(PlasmaGlob::fieldLine);
	  //gEextr[i][ic]->SetLineStyle(2);
	  gEextr[i][ic]->SetMarkerStyle(16);
	  gEextr[i][ic]->SetMarkerSize(0.4);
	  gEextr[i][ic]->SetMarkerColor(PlasmaGlob::fieldLine);	
	  gEextr[i][ic]->GetYaxis()->SetTitle("E_{y}/E_{0}");
	  gEextr[i][ic]->GetXaxis()->SetTitle("k_{p}z");	  
	}
      } else {
	nPoints = gEextr[i][ic]->GetN(); 
      }  
      
      gEextr[i][ic]->Set(nPoints+1);
      gEextr[i][ic]->SetPoint(nPoints,Time,Eextr[i][ic]);
      gEextr[i][ic]->Write(gName,TObject::kOverwrite);

      sprintf(gName,"gzEextr_%i_%i",i,ic);     
      gzEextr[i][ic] = (TGraph*) ifile->Get(gName);
      if(gzEextr[i][ic]==NULL) {
	gzEextr[i][ic] = new TGraph();
	gzEextr[i][ic]->SetName(gName);
	nPoints = 0;
	// Some cosmetics at creation time:
	if(i==1) gzEextr[i][ic]->SetLineStyle(2);
	else gzEextr[i][ic]->SetLineStyle(1);
	gzEextr[i][ic]->SetLineWidth(1);
	gzEextr[i][ic]->SetLineColor(kGray+1);
	gzEextr[i][ic]->SetMarkerStyle(20);
	gzEextr[i][ic]->SetMarkerSize(0.4);
	gzEextr[i][ic]->SetMarkerColor(kGray+1);	
	gzEextr[i][ic]->GetYaxis()->SetTitle("k_{p}#zeta]");
	gzEextr[i][ic]->GetXaxis()->SetTitle("k_{p}z");
      } else {
	nPoints = gzEextr[i][ic]->GetN(); 
      }  
      
      gzEextr[i][ic]->Set(nPoints+1);
      gzEextr[i][ic]->SetPoint(nPoints,Time,zEextr[i][ic]);
      gzEextr[i][ic]->Write(gName,TObject::kOverwrite);

    }
  }

  
  ifile->Close();
  
}

