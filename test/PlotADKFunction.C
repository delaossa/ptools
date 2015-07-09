#include <TF1.h>
#include <TF2.h>
#include <TPad.h>
#include <TLegend.h>

#ifndef __CINT__
#include "PFunctions.hh"
#endif

#include "PData.hh"
#include "PlasmaGlob.hh"

Double_t ADKfunction(Double_t *x, Double_t *par)
{
  Double_t E     = x[0]*(PUnits::GV/PUnits::m);  
  Double_t Eion0 = par[0]*PUnits::eV;           
  Int_t Z = (int)par[1];
  Int_t l = (int)par[2];
  Int_t m = (int)par[3];
    
  Double_t f = PFunc::ADK(E,Eion0,Z,l,m);

  // From atomic units to fs-1:
  f /= PUnits::atomictime/PUnits::femtosecond;
  
  return f;
}

Double_t ADK_ENGfunction(Double_t *x, Double_t *par)
{
  Double_t E     = x[0]*(PUnits::GV/PUnits::m);  
  Double_t Eion0 = par[0]*PUnits::eV;           
  Int_t Z = (int)par[1];
  
  Double_t f = PFunc::ADK_ENG(E,Eion0,Z);
  
  // to fs-1:
  f *= PUnits::femtosecond;
  
  return f;
}

Double_t PPTfunction(Double_t *x, Double_t *par)
{
  Double_t E     = x[0]*(PUnits::GV/PUnits::m);  
  Double_t Eion0 = par[0]*PUnits::eV;           
  Int_t Z = (int)par[1];
  Int_t l = (int)par[2];
  Int_t m = (int)par[3];
    
  Double_t f = PFunc::PPT(E,Eion0,Z,l,m);

  // From atomic units to fs-1:
  f /= PUnits::atomictime/PUnits::femtosecond;
    
  return f;
}

void PlotADKFunction( const TString &opt="")
{
#ifdef __CINT__
  gSystem->Load("libplasma.so");
#endif

  PlasmaGlob::Initialize();

  // Palettes!
  gROOT->Macro("PlasmaPalettes.C");
  
  gStyle->SetPadGridY(1);
  gStyle->SetPadGridX(1);

  PUnits::UnitsTable::Get();

  // Outter Helium electron
  Double_t Eion0 = 24.59;// * PUnits::eV;
  //Double_t Eion0 = 13.6;// * PUnits::eV;
  Double_t Z     = 1;
  Double_t l     = 0;
  Double_t m     = 0;
  
  // Electric field range
  Double_t Emin = 1; //PUnits::GV/PUnits::m;
  Double_t Emax = 1000; //PUnits::GV/PUnits::m;

  const Int_t NPAR = 4;
  Double_t par[NPAR] = {Eion0,Z,l,m};
  
  TF1 *fADKvsE = new TF1("fADKvsE",ADKfunction,Emin,Emax,NPAR);
  fADKvsE->SetParameters(par);  

  TF1 *fPPTvsE = new TF1("fPPTvsE",PPTfunction,Emin,Emax,NPAR);
  fPPTvsE->SetParameters(par);  

  TF1 *fADKENGvsE = new TF1("fADKENGvsE",ADK_ENGfunction,Emin,Emax,NPAR);
  fADKENGvsE->SetParameters(par);  

  
  TCanvas *C = new TCanvas("C","Tunnel-ionization probability",800,600);

  C->cd();

  gPad->SetLogx(1);

  fADKvsE->GetYaxis()->SetRangeUser(0.,10.);

  fADKvsE->GetYaxis()->SetTitle("W_{ADK} [fs^{-1}]");
  fADKvsE->GetXaxis()->SetTitle("E [GV/m]"); 
  fADKvsE->GetXaxis()->CenterTitle();
  fADKvsE->GetYaxis()->CenterTitle();
  
  fADKvsE->SetLineWidth(3);
  fADKvsE->Draw("L");

  fPPTvsE->SetLineColor(2); 
  fPPTvsE->SetLineWidth(1); 
  fPPTvsE->SetLineStyle(1);
  fPPTvsE->Draw("L same");

  fADKENGvsE->SetLineColor(kCyan); 
  fADKENGvsE->SetLineWidth(2); 
  fADKENGvsE->SetLineStyle(2);
  fADKENGvsE->Draw("L same");
    
  // Print to a file
  PlasmaGlob::imgconv(C,"./ADK-probability",opt);
  // ---------------------------------------------------------

}
