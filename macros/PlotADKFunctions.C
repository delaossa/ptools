#include <TF1.h>
#include <TF2.h>
#include <TPad.h>
#include <TLegend.h>
#include <TMarker.h>

#include "PFunctions.hh"

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

// Double_t ADKfunctionHe(Double_t *x, Double_t *par)
// {
//   Double_t E     = x[0];  
 
//   Double_t f = 159.07 * TMath::Power(13.16/E,0.488)*TMath::Exp(-4.38/E);
  
//   return f;
// }


void PlotADKFunctions( const TString &opt="")
{
#ifdef __CINT__
  gSystem->Load("libptools.so");
#endif

  PGlobals::Initialize();
  
  gStyle->SetPadGridY(0);
  gStyle->SetPadGridX(0);

  gStyle->SetPadRightMargin(0.10);

  PUnits::UnitsTable::Get();

  // Array of atomic species: http://en.wikipedia.org/wiki/Ionization_energies_of_the_elements_(data_page)
  // const Int_t Nat = 10;
  // char atNames[Nat][8] = { "H",  "He","He^{+}",  "Ar", "Ar^{+}","Ar^{++}", "Li","Li^{+}",   "N", "N^{+}"};
  // Double_t atEion[Nat] = {13.6, 24.59,    54.4, 15.76,    27.63,    40.74, 5.39,   75.64,  14.53,  29.60};
  // Double_t atZ[Nat]    = { 1.0,   1.0,     2.0,   1.0,      2.0,      3.0,  1.0,     2.0,    1.0,    2.0};
  // TF1 **fADKvsE = new TF1*[Nat]; 
  // Int_t atColor[Nat] = {kAzure+1,kOrange+8,kOrange+8,kMagenta-3,kMagenta-3,kMagenta-3,kGray+2,kGray+2,kGreen,kGreen};
  // Int_t atStyle[Nat] = {1,1,2,1,2,3,1,2,1,2};

  // // Electric field range
  // Double_t Emin[Nat] = {1,1,1,1,1,1,1,1,1,1}; //PUnits::GV/PUnits::m;
  // Double_t Emax[Nat] = {75.3,182.,602.,85.,85.,85.,12.,985.,30.,50.}; //PUnits::GV/PUnits::m;
  // const Int_t Nat = 12;
  // char atNames[Nat][8] = { "H",  "He","He^{+}",  "Ar", "Ar^{+}","Ar^{++}",   "N", "N^{+}", "Ne", "Ne+","Ne++","custom"};
  // Double_t atEion[Nat] = {13.6, 24.59,    54.4, 15.76,    27.63,    40.74, 14.53,  29.60, 21.56, 40.96, 63.45, 85};
  // Double_t atZ[Nat]    = { 1.0,   1.0,     2.0,   1.0,      2.0,      3.0,   1.0,    2.0,   1.0,   2.0, 3.0, 2};
  // TF1 **fADKvsE = new TF1*[Nat]; 
  // Int_t atColor[Nat] = {kAzure+1,kOrange+8,kOrange+8,kMagenta-3,kMagenta-3,kMagenta-3,kGray+1,kGray+1,kSpring-1,kSpring-1,kSpring-1,kRed};
  // Int_t atStyle[Nat] = {1,1,2,1,2,3,1,2,1,2,3,2};

  // Array of atomic species: http://en.wikipedia.org/wiki/Ionization_energies_of_the_elements_(data_page)
  // Double_t atEion[Nat] = { 14.53414,  29.6013, 47.44924, 77.4735, 97.8902, 552.0718, 667.046};
  // Ionization levels from NIST: http://physics.nist.gov/cgi-bin/ASD/ie.pl
  // const Int_t Nat = 11;
  // char atNames[Nat][8] = { "H",  "He","He^{+}",  "Ar", "Ar^{+}","Ar^{++}",   "N", "N^{+}", "Ne", "Ne^{+}","Ne^{++}"};
  // Double_t atEion[Nat] = {13.6, 24.59,    54.4, 15.76,    27.63,    40.74, 14.53,  29.60, 21.56, 40.96, 63.45};
  // Double_t atZ[Nat]    = { 1.0,   1.0,     2.0,   1.0,      2.0,      3.0,   1.0,    2.0,   1.0,   2.0, 3.0};
  // Int_t atColor[Nat] = {kAzure+1,kOrange+8,kOrange+8,kMagenta-3,kMagenta-3,kMagenta-3,kGray+1,kGray+1,kSpring-1,kSpring-1,kSpring-1};
  // Int_t atStyle[Nat] = {1,1,2,1,2,3,1,2,1,2,3};
  // // Electric field range
  // Double_t Emin[Nat] = {20,20,20,20,20,20,20,20,20,20,20}; //PUnits::GV/PUnits::m;
  // Double_t Emax[Nat] = {74,218.,450.,100.,150.,200.,90.,150.,180,390,390}; //PUnits::GV/PUnits::m;
  const Int_t Nat = 14;
  char atNames[Nat][16] = { "H", "H_{2}^{+}",  "He","He^{+}",      "N",   "N^{+}", "N^{2+}", "N^{3+}", "N^{4+}",  "Ne", "Ne^{+}","Ne^{2+}","Ne^{3+}","Ne^{4+}"};
  Double_t atEion[Nat] = {13.6, 30.0,  24.59,    54.4, 14.53413,  29.60125,  47.4453,  77.4735,  97.8901, 21.56,    40.96,    63.45,    97.19, 126.247};
  Double_t atZ[Nat]    = { 1.0,  2.0,    1.0,     2.0,      1.0,       2.0,      3.0,      4.0,      5.0,   1.0,      2.0,      3.0,      4.0,     5.0};
  Int_t atColor[Nat] = {kAzure+1,kAzure+2,kOrange+8,kOrange+8,kGray+2,kGray+2,kGray+2,kGray+2,kGray+2,kSpring-1,kSpring-1,kSpring-1,kSpring-1,kSpring-1};
  Int_t atStyle[Nat] = {1,2,1,2,1,2,3,4,5,1,2,3,4,5};
  // Electric field range
  Double_t Emin[Nat] = {20,50,50,150,20,60,100,200,220,50,100,150,300,400}; //PUnits::GV/PUnits::m;
  Double_t Emax[Nat] = {74,218.,218.,450.,100,150,250,450,550,180,390,390,700,1000}; //PUnits::GV/PUnits::m;

  TF1 **fADKvsE = new TF1*[Nat]; 
  TLine *lineTh[Nat];
  TMarker *markTh[Nat];

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
    Float_t E = (i+1)*(1000.0-20.0)/Npoints + 20.;
    
    for(Int_t j=0;j<Nat;j++) {
      Float_t value = fADKvsE[j]->Eval(E);
      if(value>0.1 && !found[j]) {
	IonTh[j] = E;
	found[j] = 1;
      }
    }
  }
  
  for(Int_t j=0;j<Nat;j++) {
    cout << Form("%s ion threshold: %.2f GV/m",atNames[j],IonTh[j]) << endl;
    
    // Compute OSIRIS parameters for Ionization
    Double_t Eion0 = atEion[j];
    Int_t Z = atZ[j];
    Double_t n =  3.69*Z/TMath::Sqrt(Eion0);
    Double_t A =  1.52E15 * TMath::Power(4,n) * Eion0 / ( n * TMath::Gamma(2*n) )
      * TMath::Power(20.5 * TMath::Power(Eion0,3./2.),2*n-1);
    Double_t B = 6.83 * TMath::Power(Eion0,3./2.);
    Double_t C = 2*n-1;
    
    cout << Form("OSIRIS parameters:  A = %e  B = %e  C = %e",A,B,C) << endl;
    
  }
  
  TCanvas *Canv = new TCanvas("Canv","Tunnel-ionization probability",1024,640);

  Canv->cd();

  gPad->SetLogx(1);

  TH1F *hFrame = new TH1F("hFrame","",10,5.0,1100);
  hFrame->GetYaxis()->SetRangeUser(0.,10);
  hFrame->GetYaxis()->SetTitle("W_{ADK} [fs^{-1}]");
  hFrame->GetXaxis()->SetTitle("E [GV/m]"); 
  hFrame->GetXaxis()->CenterTitle();
  hFrame->GetYaxis()->CenterTitle();

  hFrame->Draw();

  for(Int_t i=0; i<Nat; i++) {
    fADKvsE[i]->GetYaxis()->SetRangeUser(0.,10);
    fADKvsE[i]->SetLineWidth(3);
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
  
  TLegend *Leg = new TLegend(x1+0.05,y2-0.07,x1+0.35,y1+0.15);
  PGlobals::SetPaveStyle(Leg);
  for(Int_t i=0; i<Nat; i++) {
    Leg->AddEntry(fADKvsE[i],Form("%s ( %.1f )",atNames[i],IonTh[i]),"L");
  }
  // Leg->SetTextFont(42);
  Leg->SetTextColor(kGray+3);
  Leg->Draw();

  

  // Print to a file
  PGlobals::imgconv(Canv,"./ADK-probabilities",opt);
  // ---------------------------------------------------------

}
