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
#include <TEllipse.h>
#include <TExec.h>
#include <TRandom3.h>

#include "PData.hh"
#include "PlasmaGlob.hh"
#include "PPalette.hh"


bool myComp (int *x,int *y) {  
  if(x[0]<y[0]) return true;
  else if(x[0]>y[0]) return false;
  else return (x[1]<y[1]); 
}

void PlotRakeDump( const TString &sim, Int_t time, Int_t index = 0, const TString &options=""){
  
#ifdef __CINT__  
  gSystem->Load("libplasma.so");
#endif

  PlasmaGlob::Initialize();

  TString opt = options;
 
  // Palettes!
  gROOT->Macro("PlasmaPalettes.C");

  if(opt.Contains("grid")) {
    gStyle->SetPadGridX(1);
    gStyle->SetPadGridY(1);
  }
  gStyle->SetTitleAlign(22);
  gStyle->SetPadRightMargin(0.17);   // Margin for palettes in 2D histos
  //gStyle->SetTitleOffset(0.9,"z");
  gStyle->SetLabelFont(42,"xyz");
  gStyle->SetTitleFont(42,"xyz");

  // Load PData
  PData *pData = PData::Get(sim.Data());
  pData->LoadFileNames(time);
  if(!pData->IsInit()) return;

  Bool_t CYL = kFALSE;
  if(sim.Contains("cyl")) CYL = kTRUE; 
    
  Bool_t ThreeD = kFALSE;
  if(sim.Contains("3D")) ThreeD = kTRUE; 

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
  
  opt += "comovcenter";

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
  }

  char cutString[128];
  sprintf(cutString,"TMath::Abs(q)*(x1 > %f && x1 < %f && x2 > %f && x2 < %f && x3 > %f && x3 < %f)",x1Min+shiftz,x1Max+shiftz,x2Min,x2Max,x3Min,x3Max); 
  TCut Cut = cutString;

  // Bining, intervals, labels, etc.
  Int_t zNbin = 100;
  Int_t yNbin = 100;
  
  // Pointer to data TTree
  TTree *tree = pData->GetTreeRaw(pData->GetRawFileName(index)->c_str(),opt);
  // tree->Print();

  // Get phasespace histos
  TH2F *hYvsZ = NULL;
  TH2F *hYvsZcut = NULL;
  
  char hName[24];
  sprintf(hName,"hYvsZ");
  hYvsZ = new TH2F(hName,"",zNbin,x1Min,x1Max,yNbin,x2Min,x2Max);
  //  tree->Project(hName,"TMath::Sqrt(x2*x2+x3*x3):x1",Cut);
  TString dcom = Form("x2:(x1-%f)",shiftz);  
  tree->Project(hName,dcom,Cut);
  hYvsZ = (TH2F*) gROOT->FindObject(hName);    
  
  hYvsZ->GetXaxis()->CenterTitle();
  hYvsZ->GetYaxis()->CenterTitle();
  hYvsZ->GetZaxis()->CenterTitle();

  hYvsZ->GetXaxis()->SetTitle("#zeta [c/#omega_{p}]");
  hYvsZ->GetYaxis()->SetTitle("y [c/#omega_{p}]");
  hYvsZ->GetZaxis()->SetTitle("dn/d#zetady");
  
  Float_t x1Mean = hYvsZ->GetMean(1);
  Float_t x1Rms  = hYvsZ->GetRMS(1);
  Float_t x2Mean = hYvsZ->GetMean(2);
  Float_t x2Rms  = hYvsZ->GetRMS(2);

  //  cout << Form(" zmean = %f  zrms = %f",x1Mean,x1Rms) << endl;

  // Define ellipse for cutting
  Float_t zc = x1Mean;
  Float_t yc = (x2Min+x2Max)/2.0;
  Float_t ac = 2*x1Rms;
  Float_t bc = 2*x2Rms;

  // Loop for event selection:
  const Int_t Nvar = 8;
  Float_t var[Nvar];  
  char varname[Nvar][4] = {{"ene"},{"p1"},{"p2"},{"p3"},{"q"},{"x1"},{"x2"},{"x3"}};
  for(Int_t i=0;i<Nvar;i++) {
    tree->SetBranchAddress(varname[i],&var[i]);
  }
  Int_t tag[2];
  tree->SetBranchAddress("tag",&tag);

  // Main loop for selection
  Int_t nentries = (Int_t)tree->GetEntries();  
  vector<Int_t*> tags;
  Int_t it = 0;
  for(Int_t i=0;i<nentries;i++) {
    tree->GetEntry(i);
    var[5] = var[5] - shiftz - zc;
    var[6] -= yc;
    if((var[5]*var[5])/(ac*ac) + (var[6]*var[6])/(bc*bc) <= 1.0) {
      tags.push_back(new Int_t[2]);
      tags[it][0] = tag[0];
      tags[it][1] = tag[1];
      cout  << Form("%10i  %10i",tags[it][0],tags[it][1]) << endl ;
      it++;
    }
  }
  std::sort(tags.begin(),tags.end(),myComp);
  
  TString filename = Form("./%s/Plots/RakeDump/RakeDump-%s_%i.raw",sim.Data(),sim.Data(),time);
  ofstream fData(filename.Data(),ios::out);
  fData <<  Form("%10i",(int)tags.size()) << endl;  // writing the number of elements at first line overwrites first selected element.
  for(Int_t i=0;i<(int)tags.size();i++) {
    fData << Form("%10i  %10i",tags[i][0],tags[i][1]) << endl;
  }
  fData.close();
  
  // sprintf(cutString,"(x1 > %f && x1 < %f && TMath::Sqrt(x2*x2+x3*x3) > %f && TMath::Sqrt(x2*x2+x3*x3) < %f) ",zMean-zRms,zMean+zRms,x2Mean-x2Rms,x2Mean+x2Rms);
  sprintf(cutString,"(x1-%f)*(x1-%f)/(%f*%f) + (x2-%f)*(x2-%f)/(%f*%f) < 1.0",zc,zc,ac,ac,yc,yc,bc,bc);
  TCut Sel = cutString;
  sprintf(hName,"hYvsZcut");
  hYvsZcut = new TH2F(hName,"",zNbin,x1Min,x1Max,yNbin,x2Min,x2Max);
  // tree->Project(hName,"TMath::Sqrt(x2*x2+x3*x3):x1",Cut+Sel);
  tree->Project(hName,"x2:x1",Cut+Sel);
  hYvsZcut = (TH2F*) gROOT->FindObject(hName);
  
  hYvsZcut->GetXaxis()->CenterTitle();
  hYvsZcut->GetYaxis()->CenterTitle();
  hYvsZcut->GetZaxis()->CenterTitle();

  hYvsZcut->GetXaxis()->SetTitle("#zeta [c/#omega_{p}]");
  hYvsZcut->GetYaxis()->SetTitle("y [c/#omega_{p}]");
  hYvsZcut->GetZaxis()->SetTitle("dn/d#zetady");

  
  // Plotting
  // -----------------------------------------------
    
  // Canvas setup
  TCanvas *C;
  if(opt.Contains("hres") && !opt.Contains("pdf")) // high resolution for plain graphics output.
    C = new TCanvas("C","Phasespaces",1200,750);
  else
    C = new TCanvas("C","Phasespaces",960,600);

    
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
  TString fOutName = Form("./%s/Plots/RakeDump/RakeDump",sim.Data());
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

  hYvsZ->Draw("colz");

  TEllipse *ellip = new TEllipse(zc,yc,ac,bc);
  ellip->SetFillStyle(0);
  ellip->SetLineStyle(2);
  ellip->SetLineWidth(2);
  ellip->SetLineColor(kGray+2);
  ellip->Draw();

  hYvsZcut->SetMarkerStyle(6);
  hYvsZcut->SetMarkerColor(kGray+2);
  // hYvsZcut->Draw("same");

  gPad->Update();
  
  textTime->Draw();
 
  gPad->RedrawAxis(); 

  C->cd();

  // Print to a file
  PlasmaGlob::imgconv(C,fOutName,opt);
  // ---------------------------------------------------------

}
