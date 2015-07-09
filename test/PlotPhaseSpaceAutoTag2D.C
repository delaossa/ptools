#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <algorithm>    // std::sort

#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TCut.h>
#include <TChain.h>
#include <TH1F.h>
#include <TH2F.h>
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

void PlotPhaseSpaceAutoTag2D( const TString &sim, Int_t time, Int_t index = 0,TString dcom = "p1:x1", const TString &options="") {
  
#ifdef __CINT__  
  gSystem->Load("libplasma.so");
#endif

  PlasmaGlob::Initialize();

  TString opt = options;
 
  // Palettes!
  gROOT->Macro("PlasmaPalettes.C");
  
  if(opt.Contains("gridx")) {
    gStyle->SetPadGridX(1);
  } else if(opt.Contains("gridy")) {
    gStyle->SetPadGridY(1);
  }
  
  gStyle->SetTitleAlign(22);
  gStyle->SetPadRightMargin(0.17);   // Margin for palettes in 2D histos

  gStyle->SetLabelFont(42,"xyz");
  gStyle->SetTitleFont(42,"xyz");

  // Load PData
  PData *pData = PData::Get(sim.Data());
  pData->LoadFileNames(time);
  if(!pData->IsInit()) return;

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
  
  //opt += "comovcenter";

  // Centering time and z position:
  Double_t shiftz = 0.0;
  if(opt.Contains("center")) {
    Time -= zStartPlasma;
    if(opt.Contains("comov")) {     // Centers on the head of the beam.
      Time += zStartBeam;
      shiftz += zStartBeam;
    } else {
      shiftz += zStartPlasma;
    }
  } 
  if(opt.Contains("comov")) {
    Double_t v = pData->GetBeamVelocity();    
    if(v==0) v = 1.0; // If equals to 0 (default), then set to c.
    shiftz += v * pData->GetRealTime();
  }   

  TString sshiftz = Form("(x1-%f)",shiftz);
  
  //cout << Form("string = %s",sshiftz.Data()) << endl; 
  TString sx1 = "x1";
  dcom.ReplaceAll(sx1,sshiftz);
  //cout << Form("%s replaced by %s result = %s",sx1.Data(),sshiftz.Data(),dcom.Data()) << endl;
  // ----

  // Weighting by the macroparticle charge:
  char cutString[128];
  sprintf(cutString,"TMath::Abs(q)"); 
  
  TCut Cut = cutString;
  
  // Bining, intervals, labels, etc.
  Int_t xNbin = 200;
  Int_t yNbin = 200;
  
  // Pointer to data TTree
  TTree *tree = pData->GetTreeRaw(pData->GetRawFileName(index)->c_str(),opt);
  
  // Get phasespace histos
  TH2F *hPha2D = NULL;
  char hName[24]; 
  sprintf(hName,"hPha2D");
  char dCom[128];
  sprintf(dCom,"%s>>%s(%i,%i)",dcom.Data(),hName,xNbin,yNbin);

  tree->Draw(dCom,Cut,"goff");
  hPha2D = (TH2F*) gROOT->FindObject(hName);    
  hPha2D->GetXaxis()->CenterTitle();
  hPha2D->GetYaxis()->CenterTitle();
  hPha2D->GetZaxis()->CenterTitle();

  // hPha2D->GetXaxis()->SetTitle("#zeta [c/#omega_{p}]");
  // hPha2D->GetYaxis()->SetTitle("y [c/#omega_{p}]");
  // hPha2D->GetZaxis()->SetTitle("dn/d#zetady");
  
  Float_t xMean = hPha2D->GetMean(1);
  Float_t xRms  = hPha2D->GetRMS(1);
  Float_t yMean = hPha2D->GetMean(2);
  Float_t yRms  = hPha2D->GetRMS(2);

  cout << Form(" xMean = %f  xRms = %f   yMean = %f  yRms = %f",xMean,xRms,yMean,yRms) << endl;

  // Read array of tags:
  TString filename = Form("./%s/Plots/RakeDump/RakeDump-%s_%i.raw",sim.Data(),sim.Data(),100);
  ifstream fDataIn(filename.Data(),ios::in);
  Int_t Npart;
  fDataIn >> Npart;
  vector<Int_t*> tags;
  Int_t it = 0;
  while (!fDataIn.eof()) {
    tags.push_back(new Int_t[2]);
    fDataIn >> tags[it][0] >> tags[it][1];
    //  cout <<  Form("%10i  %10i",tags[it][0],tags[it][1]) << endl;
    it++;
  }
  
  // Redo the histogram with only tagged particles:
  TH2F *hPha2DTag = (TH2F*) hPha2D->Clone("hPha2DTag");
  hPha2DTag->Reset();
  
  // Loop for event selection:
  const Int_t Nvar = 8;
  Float_t var[Nvar];  
  char varname[Nvar][4] = {{"ene"},{"p1"},{"p2"},{"p3"},{"q"},{"x1"},{"x2"},{"x3"}};
  for(Int_t i=0;i<Nvar;i++) {
    tree->SetBranchAddress(varname[i],&var[i]);
  }
  Int_t tag[2];
  tree->SetBranchAddress("tag",&tag);
  
  Int_t nentries = (Int_t)tree->GetEntries();  
  for(Int_t i=0;i<nentries;i++) {
    tree->GetEntry(i);
    var[5] = var[5] - shiftz;
    Bool_t found = kFALSE;
    for(Int_t it = 0; it<(int)tags.size(); it++) {
      if((tag[0]==tags[it][0]) && (tag[1]==tags[it][1])) {
	hPha2DTag->Fill(var[5],var[6]);
	//tags.erase(tags.begin()+it);
	continue;
      }
    }
  }
  
  
  
  // Plotting
  // -----------------------------------------------
  
  // Canvas setup
  TCanvas *C;
  if(opt.Contains("hres") && !opt.Contains("pdf")) // high resolution for plain graphics output.
    C = new TCanvas("C","Phasespace",1200,750);
  else
    C = new TCanvas("C","Phasespace",960,600);

    
  // Text objects
  TPaveText *textTime = new TPaveText(0.55,0.85,0.82,0.9,"NDC");
  PlasmaGlob::SetPaveTextStyle(textTime,32); 
  char ctext[128];
  if(opt.Contains("units") && pData->GetPlasmaDensity()) 
    sprintf(ctext,"z = %5.1f mm", Time * skindepth / PUnits::mm);
  else
    sprintf(ctext,"t = %5.1f #omega_{p}^{-1}",Time);
  textTime->AddText(ctext);
 

  // Actual Plotting!
  // ------------------------------------------------------------

  // Output file
  TString fOutName = Form("./%s/Plots/PhaseSpaceAutoTag2D/PhaseSpaceAutoTag2D",sim.Data());
  fOutName += Form("-%s_%i",sim.Data(),time);

  C->cd();
 
  gPad->SetFrameLineWidth(2);  

  if(opt.Contains("logz")) {
    gPad->SetLogz(1);
  } else {
    gPad->SetLogz(0);
  }

  // Set palette:
  PPalette * pPalette = (PPalette*) gROOT->FindObject("electron");
  pPalette->cd();

  hPha2D->Draw("colz");

  hPha2DTag->SetMarkerStyle(6);
  hPha2DTag->SetMarkerColor(kGray+2);
  hPha2DTag->Draw("same");
  
  textTime->Draw();
 
  gPad->RedrawAxis(); 

  C->cd();

  // Print to a file
  PlasmaGlob::imgconv(C,fOutName,opt);
  // ---------------------------------------------------------

}
