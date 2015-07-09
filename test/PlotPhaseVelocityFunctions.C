#include <TF1.h>
#include <TF2.h>
#include <TPad.h>
#include <TLegend.h>

#ifndef __CINT__
#include "PFunctions.hh"
#endif

#include "PData.hh"
#include "PlasmaGlob.hh"

Double_t PhaseVelocityVsZ(Double_t *x, Double_t *par)
{
  Float_t z = x[0];
  Float_t zg = par[0];
  Float_t r0 = par[1];
  Float_t gamma = par[2];
  Float_t nb = par[3];
  Float_t n0 = par[4];
  
  Double_t f = PFunc::PhaseVelocity(z,zg,r0,gamma,nb,n0);
  return f - 1.0;
}

Double_t PhaseVelocityVsZg(Double_t *x, Double_t *par)
{
  Float_t zg = x[0];
  Float_t z = par[0];
  Float_t r0 = par[1];
  Float_t gamma = par[2];
  Float_t nb = par[3];
  Float_t n0 = par[4];
  
  Double_t f = PFunc::PhaseVelocity(z,zg,r0,gamma,nb,n0);
  return f - 1.0;
}

Double_t PhaseVelocityVsZVsZg(Double_t *x, Double_t *par)
{
  Float_t z = x[0];
  Float_t zg = x[1];
  Float_t r0 = par[0];
  Float_t gamma = par[1];
  Float_t nb = par[2];
  Float_t n0 = par[3];
  
  Double_t f = PFunc::PhaseVelocity(z,zg,r0,gamma,nb,n0);
  return f - 1.0;
}


Double_t PhaseVelocityVsZ2(Double_t *x, Double_t *par)
{
  Float_t z = x[0];
  Float_t zg = par[0];
  Float_t r0 = par[1];
  Float_t gamma = par[2];
  Float_t nb = par[3];
  Float_t n0 = par[4];
  
  Double_t f = PFunc::PhaseVelocity2(z,zg,gamma,nb,n0);
  return f - 1.0;
}

Double_t PhaseVelocityVsZg2(Double_t *x, Double_t *par)
{
  Float_t zg = x[0];
  Float_t z = par[0];
  Float_t r0 = par[1];
  Float_t gamma = par[2];
  Float_t nb = par[3];
  Float_t n0 = par[4];
  
  Double_t f = PFunc::PhaseVelocity2(z,zg,gamma,nb,n0);
  return f - 1.0;
}


Double_t PhaseVelocityVsZVsZg2(Double_t *x, Double_t *par)
{
  Float_t z = x[0];
  Float_t zg = x[1];
  Float_t r0 = par[0];
  Float_t gamma = par[1];
  Float_t nb = par[2];
  Float_t n0 = par[3];
  
  Double_t f = PFunc::PhaseVelocity2(z,zg,gamma,nb,n0);
  return f - 1.0;
}


void PlotPhaseVelocityFunctions( const TString &opt="")
{
#ifdef __CINT__
  gSystem->Load("libplasma.so");
#endif

  PlasmaGlob::Initialize();

  // Palettes!
  gROOT->Macro("PlasmaPalettes.C");
  
  gStyle->SetPadTopMargin(0.06); 
  gStyle->SetPadGridY(0);
  gStyle->SetPadGridX(0);
  gStyle->SetFrameLineWidth(2);
  gStyle->SetLabelSize(0.04, "xyz");
  gStyle->SetTitleOffset(1.2,"y");
  gStyle->SetTitleOffset(1.2,"z");
  gStyle->SetNdivisions(505,"xyz");

  PUnits::UnitsTable::Get();
  
  // // SPS parameters from Pukhov, Kumar et al. PRL107,145003(2011)
  // Double_t n0 = 7.76e14 / PUnits::cm3;
  // Double_t nb = 1.5e12 / PUnits::cm3;
  // Double_t lambda = PFunc::PlasmaWavelength(n0);
  // Double_t kp = PFunc::PlasmaWavenumber(n0);
  // Double_t skindepth =  PFunc::PlasmaSkindepth(n0);
  // Double_t r0 = 0.19 * PUnits::mm;
  // Double_t z  = 2.5 * PUnits::m;
  // Double_t zg = -28.6 * PUnits::mm;
  // Double_t E  = 450 * PUnits::GeV;
  // Double_t gamma = E / PConst::ProtonMassE ;

  // SPS parameters from Schroeder et al. PRL107,145002(2011)
  // Double_t n0 = 1.0;
  // Double_t nb = 0.002;
  // Double_t lambda = TMath::TwoPi();
  // Double_t kp = 1.0;
  // Double_t r0 = 1.0;
  // Double_t z  = 13105;
  // Double_t E  = 450 * PUnits::GeV;
  // Double_t gamma = E / PConst::ProtonMassE ;

  // PITZ parameters
  Double_t n0 = 1.0e15 / PUnits::cm3;
  Double_t nb = 1.05e13 / PUnits::cm3;
  Double_t lambda = PFunc::PlasmaWavelength(n0);
  Double_t kp = PFunc::PlasmaWavenumber(n0);
  Double_t skindepth =  PFunc::PlasmaSkindepth(n0);
  Double_t r0 = 34.259 * PUnits::um;
  Double_t z  = 150. * PUnits::mm;
  Double_t zg = -5. * PUnits::mm;
  Double_t E  = 25 * PUnits::MeV;
  Double_t gamma = (E / PConst::ElectronMassE) + 1.0;
  Double_t vb = TMath::Sqrt(1. - (1./(gamma*gamma)));

  cout << " n0 = " << n0 * PUnits::cm3 << "  e/cc" << endl;
  cout << " Wavelength = " << PUnits::BestUnit(lambda, "Length") << endl;
  cout << " Skindepth = " << PUnits::BestUnit(skindepth, "Length") << endl;
  cout << " Wavenumber = " << kp * PUnits::mm << "  mm^-1" << endl;

  // Normalized variables
  if(!opt.Contains("units")){
    r0 *= kp;
    nb /= n0;
    n0 = 1.0;
    z *= kp;
    zg *= kp;
  }
  

  Double_t N = PFunc::Nefoldings(z,zg,r0,gamma,nb,n0);
  Double_t r1 = (TMath::Power(3.,1./4.)/TMath::Power(8.*TMath::Pi()*N,1./2.)) * TMath::Exp(N) ;
  Double_t Gamma = PFunc::PhaseGamma(z,zg,r0,gamma,nb,n0);
  Double_t vph = PFunc::PhaseVelocity(z,zg,r0,gamma,nb,n0);
  Double_t vph2 = PFunc::PhaseVelocity2(z,zg,gamma,nb,n0);

  cout << " gamma beam = " << gamma << endl;
  cout << " zeta = " << z << endl;
  cout << " zeta comov. = " << zg << endl;
  cout << " r0 = " << r0 << endl;
  cout << " n_b0/n_0 = " << nb/n0 << endl;
  cout << " nu constant = " << PFunc::Nu(r0,n0) << endl;
  cout << " Number e-foldings = " << N << endl; 
  cout << " r1 = " << r1 << endl;
  cout << " Gamma wake = " << Gamma << endl;
  cout << " Wake Phase velocity  = " << vph << endl;
  cout << " Wake Phase velocity2 = " << vph2 << endl;
  cout << " Beam phase velocity  = " << vb << endl;

  const Int_t NPAR = 5;
  Double_t par[NPAR] = {zg,r0,gamma,nb,n0};
  TF1 *fPhaseVsZ2 = new TF1("fPhaseVsZ2",PhaseVelocityVsZ2,0,z,NPAR);
  fPhaseVsZ2->SetParameters(par);  
  TF1 *fPhaseVsZ = new TF1("fPhaseVsZ",PhaseVelocityVsZ,0,z,NPAR);
  fPhaseVsZ->SetParameters(par);

  Double_t par2[NPAR] = {z,r0,gamma,nb,n0};
  TF1 *fPhaseVsZg2 = new TF1("fPhaseVsZg2",PhaseVelocityVsZg2,zg,0.,NPAR);
  fPhaseVsZg2->SetParameters(par2);
  TF1 *fPhaseVsZg = new TF1("fPhaseVsZg",PhaseVelocityVsZg,zg,0.,NPAR);
  fPhaseVsZg->SetParameters(par2);
  
  const Int_t NPAR2D = 4;
  Double_t par3[NPAR2D] = {r0,gamma,nb,n0};
  TF2 *fPhaseVsZVsZg2 = new TF2("fPhaseVsZVsZg2",PhaseVelocityVsZVsZg2,0.,z,zg,0.,NPAR2D);
  fPhaseVsZVsZg2->SetParameters(par3);
  TF2 *fPhaseVsZVsZg = new TF2("fPhaseVsZVsZg",PhaseVelocityVsZVsZg,0.,z,zg,0.,NPAR2D);
  fPhaseVsZVsZg->SetParameters(par3);
  
  char ctext[64];
  TPaveText *textZetag = new TPaveText(0.13,0.85,0.38,0.92,"NDC");
  PlasmaGlob::SetPaveTextStyle(textZetag,12); 
  textZetag->SetTextColor(kGray+3);
  if(opt.Contains("units"))
    sprintf(ctext,"#zeta_{0} = %6.2f mm", zg / PUnits::mm);
  else
    sprintf(ctext,"#zeta_{0} = %6.2f c/#omega_{p}", zg);
  textZetag->AddText(ctext);
  
  TPaveText *textZeta = new TPaveText(0.13,0.85,0.38,0.92,"NDC");
  PlasmaGlob::SetPaveTextStyle(textZeta,12); 
  textZeta->SetTextColor(kGray+3);
  if(opt.Contains("units"))
    sprintf(ctext,"z_{0} = %6.2f mm", z / PUnits::mm);
  else
    sprintf(ctext,"z_{0} = %6.2f c/#omega_{p}", z);
  textZeta->AddText(ctext);

  // Graph for de-phasing:
  const Int_t NP = 100;
  TGraph *gPhase = new TGraph(NP);
  TGraph *gPhase2 = new TGraph(NP);
  Float_t phase  = zg;
  Float_t phase2 = zg;
  Float_t Dz = z / NP;
  for(Int_t i=0;i<NP;i++) {
    Float_t zp = (i+1)*Dz;
    Float_t v = PFunc::PhaseVelocity(zp,phase,r0,gamma,nb,n0);
    phase += (v - vb) * Dz;
    // cout << " z = " << zp << "  phase = " << phase << endl;
    if(opt.Contains("units"))
      gPhase->SetPoint(i,zp,(phase-zg)*kp);
    else
      gPhase->SetPoint(i,zp,(phase-zg));

    v = PFunc::PhaseVelocity2(zp,phase2,gamma,nb,n0);
    phase2 += (v - 1) * Dz;
    // cout << " z = " << zp << "  phase = " << phase << endl;
    if(opt.Contains("units"))
      gPhase2->SetPoint(i,zp,(phase2-zg)*kp);
    else
      gPhase2->SetPoint(i,zp,(phase2-zg));
  }
  
  TCanvas *C = new TCanvas("C","Wake phase velocity",1000,750);
  C->Divide(2,2);

  TLegend *Leg = new TLegend(0.6,0.20,0.85,0.35);
  PlasmaGlob::SetPaveStyle(Leg);
  Leg->SetTextAlign(22);
  Leg->SetTextColor(kGray+3);
  Leg->SetLineColor(1);
  Leg->SetBorderSize(1);
  Leg->SetFillColor(0);
  Leg->SetFillStyle(1001);
  //Leg-> SetNColumns(2);
  Leg->AddEntry(fPhaseVsZ,"PRL 107,145002","L");
  Leg->AddEntry(fPhaseVsZ2,"PRL 107,145003","L");
  Leg->SetTextColor(kGray+3);
  

  C->cd(1);
  fPhaseVsZ->GetYaxis()->SetTitle("(v_{p} - c)/c");
  if(opt.Contains("units")) 
    fPhaseVsZ->GetXaxis()->SetTitle("z [m]"); 
  else
    fPhaseVsZ->GetXaxis()->SetTitle("z [c/#omega_{p}]"); 
  fPhaseVsZ->GetYaxis()->SetRangeUser(-3e-4,0.);
  fPhaseVsZ->GetXaxis()->CenterTitle();
  fPhaseVsZ->GetYaxis()->CenterTitle();
  fPhaseVsZ->SetLineWidth(2);
  fPhaseVsZ->Draw("C");
  fPhaseVsZ2->SetLineWidth(2);
  fPhaseVsZ2->SetLineStyle(2);
  fPhaseVsZ2->Draw("C same");
  textZetag->Draw();
  Leg->Draw();

  C->cd(2);
  fPhaseVsZg->GetYaxis()->SetTitle("(v_{p} - c)/c");
  if(opt.Contains("units"))
    fPhaseVsZg->GetXaxis()->SetTitle("#zeta [m]"); 
  else
    fPhaseVsZg->GetXaxis()->SetTitle("#zeta [c/#omega_{p}]"); 
  fPhaseVsZg->GetYaxis()->SetRangeUser(-3e-4,0.);
  fPhaseVsZg->GetYaxis()->SetNdivisions(505);
  fPhaseVsZg->GetXaxis()->CenterTitle();
  fPhaseVsZg->SetLineWidth(2);
  fPhaseVsZg->GetYaxis()->CenterTitle();
  fPhaseVsZg->Draw("C");
  fPhaseVsZg2->SetLineWidth(2);
  fPhaseVsZg2->SetLineStyle(2);
  fPhaseVsZg2->Draw("C same");
  textZeta->Draw();
  
  C->cd(3);
  gPhase->GetYaxis()->SetTitle("#Delta#zeta [c/#omega_{p}]");
  
  if(opt.Contains("units")) {
    gPhase->GetXaxis()->SetTitle("z [m]"); 
    gPhase->GetXaxis()->SetNdivisions(510);
  } else
    gPhase->GetXaxis()->SetTitle("z [c/#omega_{p}]"); 
 
  
  // gPhase->GetXaxis()->SetNdivisions(510);
  gPhase->GetXaxis()->CenterTitle();
  gPhase->GetYaxis()->CenterTitle();
  gPhase->GetXaxis()->SetRangeUser(0.,z);
  gPhase->GetYaxis()->SetNdivisions(505);  
  gPhase->SetLineWidth(2);
  gPhase->Draw("AC");
  gPhase2->SetLineStyle(2);
  gPhase2->SetLineWidth(2);
  gPhase2->Draw("C");
  textZetag->Draw();

  C->cd(4);
  gPad->SetLeftMargin(0.18);

  if(opt.Contains("units")) {
    fPhaseVsZVsZg->GetYaxis()->SetTitle("#zeta [m]"); 
  } else
    fPhaseVsZVsZg->GetYaxis()->SetTitle("#zeta [c/#omega_{p}]"); 
  
  if(opt.Contains("units")) {
    fPhaseVsZVsZg->GetXaxis()->SetTitle("z [m]"); 
  } else
    fPhaseVsZVsZg->GetXaxis()->SetTitle("z [c/#omega_{p}]"); 
  
  fPhaseVsZVsZg->GetZaxis()->SetTitle("(v_{p} - c)/c"); 

  fPhaseVsZVsZg->GetYaxis()->SetNdivisions(505);
  // fPhaseVsZVsZg->GetXaxis()->SetNdivisions(510);
  fPhaseVsZVsZg->GetXaxis()->CenterTitle();
  fPhaseVsZVsZg->GetXaxis()->SetTitleOffset(1.6);
  fPhaseVsZVsZg->GetYaxis()->CenterTitle();
  fPhaseVsZVsZg->GetYaxis()->SetTitleOffset(2.0); 
  fPhaseVsZVsZg->GetZaxis()->CenterTitle();
  fPhaseVsZVsZg->GetZaxis()->SetTitleOffset(1.4); 
  
  fPhaseVsZVsZg->Draw("surf2");

  C->cd();

  // Print to a file
  PlasmaGlob::imgconv(C,"./WakePhaseVelocity",opt);
  // ---------------------------------------------------------

}
