#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>

#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>

using namespace std;

void GetHiRawData2(const TString &sim, Float_t time, Int_t index = 0, const TString &options="") {
  
  TString filename = Form("%s/DATA/InjectedBeamPhaseSpace_time_%.3f",sim.Data(),time);
  cout << filename.Data() << endl;
  
  FILE * pFile;
  long   lSize;
  double * buffer;
  size_t result;

  pFile = fopen (filename.Data(),"rb" );
  if (pFile==NULL) {fputs ("File error",stderr); exit (1);}
  
  // obtain file size:
  fseek (pFile , 0 , SEEK_END);
  lSize = ftell (pFile);
  rewind (pFile);

  // allocate memory to contain the whole file:
  buffer = (double*) malloc (lSize);
  if (buffer == NULL) {fputs ("Memory error",stderr); exit (2);}

  // copy the file into the buffer:
  result = fread (buffer,1,lSize,pFile);
  if (result != lSize) {fputs ("Reading error",stderr); exit (3);}
  
  /* the whole file is now loaded in the memory buffer. */

  // terminate
  fclose (pFile);


  TTree *tree = new TTree("HiTree","oeoeoeo");
  const Int_t Nvar = 9;
  Double_t var[Nvar];  
  char varname[Nvar][8] = {{"x1"},{"x2"},{"x3"},{"p1"},{"p2"},{"p3"},{"q"},{"ipart"},{"iproc"}};
  for(Int_t i=0;i<Nvar;i++) {
    char vartype[8];
    sprintf(vartype,"%s/D",varname[i]);
    tree->Branch(varname[i],&var[i],vartype);
  }

  Int_t Npart=lSize/(Nvar*sizeof(double));
  
  for(int i_part=0; i_part<Npart; i_part++) {
    for(Int_t i_var=0; i_var<Nvar; i_var++) {
      var[i_var]=buffer[i_var+i_part*Nvar];
      //cout << Form ("%e  ",var[i_var]);
    }
    // cout << endl;
    tree->Fill();
  }
  
  free (buffer);
    
  cout << Form("  %i  particles read! " , Npart) << endl;
  
  TCanvas *C = new TCanvas("C","",1024,640);
  C->cd();
  
  tree->Draw("p2:x2","q");
  
  
}



