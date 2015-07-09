#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>

#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>

using namespace std;

void GetHiRawData(const TString &sim, Float_t time, Int_t index = 0, const TString &options="") {
  
  TString filename = Form("%s/DATA/InjectedBeamPhaseSpace_time_%.3f",sim.Data(),time);
  cout << filename.Data() << endl;
  
  ifstream fDataIn(filename.Data(),ios::in);
  
  TTree *tree = new TTree("HiTree","oeoeoeo");
  const Int_t Nvar = 7;
  Double_t var[Nvar];  
  char varname[Nvar][4] = {{"x1"},{"x2"},{"x3"},{"p1"},{"p2"},{"p3"},{"q"}};
  for(Int_t i=0;i<Nvar;i++) {
    char vartype[8];
    sprintf(vartype,"%s/D",varname[i]);
    tree->Branch(varname[i],&var[i],vartype);
  }
  
  Int_t Npart = 0;
  while(!fDataIn.eof()) {
    for(Int_t i=0;i<Nvar;i++) {
      fDataIn >> var[i];
      cout << Form ("%e  ",var[i]);
    }
    cout << endl;
    tree->Fill();
    Npart++;
  }
  
  cout << Form("  %i  particles read! " , Npart) << endl;
  
  TCanvas *C = new TCanvas("C","",1024,640);
  C->cd();
  
  tree->Draw("x2:x1");
  
  
}



