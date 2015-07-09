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
#include <TCanvas.h>
#include <TStyle.h>
#include <TExec.h>

#include "PData.hh"

using namespace std;

int main(int argc,char *argv[]) {
  if(argc<=4) {
    printf("\n Usage: %s <source simulation> <target dir> <initial time> <final time> <time step>\n\n",argv[0]);
    return 0;
  }

  // General options
  TString  sim  = "";
  TString  simo = "";
  Int_t    time = 0;
  Int_t  iStart = 0;
  Int_t    iEnd = 0;
  Int_t   iStep = 1;
  TString  opt  = "";

  // Interfacing command line:
  for(int l=1;l<argc;l++){
    TString arg = argv[l];

    switch(l) {
    case 1 :
      sim = arg;
      break;
    case 2:
      simo = arg;
      break;
    case 3:
      iStart = arg.Atoi();
      break;
    case 4:
      iEnd = arg.Atoi();
      break;
    case 5:
      iStep = arg.Atoi();
      break;
    }
  }
  
  // The Data manager
  PData *pData = PData::Get(sim.Data());
  
  // Time looper
  for(Int_t i=iStart; i<iEnd+1; i+=iStep) {
    
    time = i;
    pData->LoadFileNames(time);
    
    if(pData->IsInit()) {
      cout << Form("\nCopying files for simulation %s at time step %i:\n",sim.Data(),time) << endl;
      
      pData->CopyData(simo.Data());
      
    } else {
      continue;
    }
    
  }
  
  if(pData)
    delete pData;

  // ---------------------------------------------------------
}
