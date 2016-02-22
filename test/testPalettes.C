#include "PPalette.hh"


void testPalettes(){
  // Here the Palettes are defined
  gROOT->Macro("PPalettes.C");

  gStyle->SetOptStat(0);
  gStyle->SetNumberContours(64);

  TCanvas *c = new TCanvas("c","c",0,700,1000,625);

  TH2F *h1 = new TH2F("h1","h1",100,-4,4,100,-4,4);
  TH2F *h2 = new TH2F("h2","h2",100,-4,4,100,-4,4);
  TH2F *h3 = new TH2F("h3","h3",100,-4,4,100,-4,4);
  TH2F *h4 = new TH2F("h4","h4",100,-4,4,100,-4,4);

  Double_t a,b;
  for (Int_t i=0;i<50000;i++) {
    gRandom->Rannor(a,b);
    h1->Fill(a-1.5,b-1.5);
    h2->Fill(a+1.5,b+1.5);
    h3->Fill(a-1.5,b+1.5);
    h4->Fill(a+1.5,b-1.5);
  }

  //  plasmaPalette->SetAlpha(0.7);

  TExec *ex1 = new TExec("ex1","plasmaPalette->SetAlpha(0.7); plasmaPalette->SetPalette(\"gray\");");
  TExec *ex2 = new TExec("ex2","plasmaPalette->SetPalette(\"rbow0\");");
  TExec *ex3 = new TExec("ex3","plasmaPalette->SetPalette(\"electron0\");");
  TExec *ex4 = new TExec("ex4","plasmaPalette->SetPalette(\"oli\");");

  h1->Draw("axis");
  ex1->Draw(); h1->Draw("col same");

  ex2->Draw(); h2->Draw("col same");

  ex3->Draw(); h3->Draw("col same");

  ex4->Draw(); h4->Draw("col same");

  c->Print("./testPalettes.pdf");
}
