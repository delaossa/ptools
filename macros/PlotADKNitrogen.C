#include <TF1.h>
#include <TF2.h>
#include <TPad.h>
#include <TLegend.h>
#include <TMarker.h>

#ifndef __CINT__
#include "PFunctions.hh"
#endif

#include "PData.hh"
#include "PGlobals.hh"

Double_t ADKfunction(Double_t *x, Double_t *par)
{
  Double_t E     = x[0]*(PUnits::GV/PUnits::m);  
  Double_t Eion0 = par[0]*PUnits::eV;           
  Int_t Z = (int)par[1];
  
  Double_t f = PFunc::ADK_ENG(E,Eion0,Z);
  
  // to fs-1:
  f *= PUnits::femtosecond;
  
  return f;
}


void PlotADKNitrogen( const TString &opt="")
{
#ifdef __CINT__
  gSystem->Load("libptools.so");
#endif

  PGlobals::Initialize();

  // Palettes!
  gROOT->Macro("PPalettes.C");
  
  gStyle->SetPadGridY(0);
  gStyle->SetPadGridX(0);

  gStyle->SetTextFont(42);
  gStyle->SetTitleFont(42,"xyz");
  gStyle->SetLabelFont(42,"xyz");

  gStyle->SetLabelOffset(0.014, "y");

  gStyle->SetTextFont(42);

  gStyle->SetPadRightMargin(0.10);

  PUnits::UnitsTable::Get();

  const Int_t Nat = 7;
  char atNames[Nat][8] = {   "N^{+}", "N^{2+}", "N^{3+}", "N^{4+}", "N^{5+}", "N^{6+}", "N^{7+}"};
  // Array of atomic species: http://en.wikipedia.org/wiki/Ionization_energies_of_the_elements_(data_page)
  // Double_t atEion[Nat] = { 14.53414,  29.6013, 47.44924, 77.4735, 97.8902, 552.0718, 667.046};
  // Ionization levels from NIST: http://physics.nist.gov/cgi-bin/ASD/ie.pl
  Double_t atEion[Nat] = { 14.53413,  29.60125, 47.4453, 77.4735, 97.89013, 552.06731, 667.04609};
  Double_t atZ[Nat]    = {      1.0,      2.0,      3.0,     4.0,     5.0,      6.0,     7.0};
  TF1 **fADKvsE = new TF1*[Nat]; 
  // Int_t atColor[Nat] = {kGray+2,kGray+2,kGray+2,kGray+2,kGray+2,kGray+2,kGray+2};
  Int_t atColor[Nat] = {kGray+2,kGray+2,kGray+2,kGray+2,kGray+2,kGray+2,kGray+2};
  Int_t atStyle[Nat] = {1,2,3,4,5,6,7};
  TLine *lineTh[Nat];
  TMarker *markTh[Nat];
  
  // Electric field range
  Double_t Emin[Nat] = {20,60,100,200,220,5000,6000}; //PUnits::GV/PUnits::m;
  Double_t Emax[Nat] = {100,150,250,450,550,12000,14000}; //PUnits::GV/PUnits::m;
  
  const Int_t NPAR = 2;
  Double_t par[Nat][NPAR];
  for(Int_t i=0; i<Nat; i++) {
    par[i][0] = atEion[i];
    par[i][1] = atZ[i];
    
    char name[24];
    sprintf(name,"fADKvsE_%s",atNames[i]);
    fADKvsE[i] = new TF1(name,ADKfunction,Emin[i],1*Emax[i],NPAR);
    fADKvsE[i]->SetParameters(par[i]);  
  }

  // Evaluate the function to find Ion thresholds
  Int_t Npoints = 10000;
  Float_t IonTh[Nat];
  Bool_t found[Nat] = {0};
  for(Int_t i=0; i<Npoints; i++) {
    for(Int_t j=0;j<Nat;j++) {
      Float_t E = (i+1)*(Emax[j]-Emin[j])/Npoints + Emin[j];
      
      Float_t value = fADKvsE[j]->Eval(E);
      if(value>0.1 && !found[j]) {
	IonTh[j] = E;
	found[j] = 1;
      }
    }
  }
  
  for(Int_t j=0;j<Nat;j++) {
    cout << Form("%s ion threshold: %.2f",atNames[j],IonTh[j]) << endl;
    
    // Compute OSIRIS parameters for Ionization
    Double_t Eion0 = atEion[j];
    Int_t Z = atZ[j];
    Double_t n =  3.69*Z/TMath::Sqrt(Eion0);
    Double_t A =  1.52E15 * TMath::Power(4,n) * Eion0 / ( n * TMath::Gamma(2*n) )
      * TMath::Power(20.5 * TMath::Power(Eion0,3./2.),2*n-1);
    Double_t B = 6.83 * TMath::Power(Eion0,3./2.);
    Double_t C = 2*n-1;
    
    cout << Form("OSIRIS coefficients") << endl;
    cout << Form("  %e   %e   %e",A,B,C) << endl;
    
  }
  
  TCanvas *Canv = new TCanvas("Canv","Tunnel-ionization probability",1024,640);

  Canv->cd();

  gPad->SetLogx(1);
  gPad->SetTickx(1);

  TH1F *hFrame = new TH1F("hFrame","",10,5.0,20000);
  hFrame->GetYaxis()->SetRangeUser(0.,10);
  hFrame->GetYaxis()->SetTitle("W_{ADK} [fs^{-1}]");
  hFrame->GetXaxis()->SetTitle("E [GV/m]"); 
  hFrame->GetXaxis()->CenterTitle();
  hFrame->GetYaxis()->CenterTitle();

  hFrame->Draw();

  for(Int_t i=0; i<Nat; i++) {
    fADKvsE[i]->GetYaxis()->SetRangeUser(0.,10);
    fADKvsE[i]->SetLineWidth(2);
    fADKvsE[i]->SetLineColor(atColor[i]);
    fADKvsE[i]->SetLineStyle(atStyle[i]);
    fADKvsE[i]->Draw("C same");    
  }

  for(Int_t i=0; i<Nat; i++) {
    lineTh[i] = new TLine(IonTh[i],0.0,IonTh[i],0.3);
    lineTh[i]->SetLineColor(atColor[i]);
    lineTh[i]->SetLineWidth(1);
    lineTh[i]->Draw();

    markTh[i] = new TMarker(IonTh[i],0.3,20);
    markTh[i]->SetMarkerColor(atColor[i]);
    markTh[i]->SetMarkerSize(1);
    markTh[i]->Draw();

  }
  
  gPad->Update();
  Float_t y1 = gPad->GetBottomMargin();
  Float_t y2 = 1 - gPad->GetTopMargin();
  Float_t x1 = gPad->GetLeftMargin();
  //  Float_t x2 = 1 - gPad->GetRightMargin();
  
  TLegend *Leg = new TLegend(x1+0.05,y1+0.30,x1+0.24,y2-0.07);
  PGlobals::SetPaveStyle(Leg);
  for(Int_t i=0; i<Nat; i++) {
    Leg->AddEntry(fADKvsE[i],Form("%s (%.1f eV)",atNames[i],atEion[i]),"L");
  }
  Leg->SetTextFont(82);
  Leg->SetTextColor(kGray+2);
  Leg->Draw();

  

  // Print to a file
  PGlobals::imgconv(Canv,"./ADK-probabilities-N",opt);
  // ---------------------------------------------------------

}
