#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>

#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TExec.h>

//#include "PData.hh"
//#include "PDataHiP.hh"

using namespace std;

int main(int argc,char *argv[]) {
  
  TString IncPath = gSystem->Getenv("PTOOLS");
  IncPath +=  "/inc";
  gROOT->ProcessLine(Form(".include %s",IncPath.Data()));
  
  if(argc<=1) {
    printf("\n Usage: %s <MACRO> <-t(time)>\n",argv[0]);
    printf("      <-i(initial time)> <-f(final time)> <-s(time step)>\n");
    printf("      <-m(plot mask)>\n");
    printf("      <--png> <--pdf> <--eps> <--units> <--comov> <--hres>\n\n");
    return 0;
  }
  
  TString   sim = "";
  Int_t    time = 0;
  Int_t  iStart = -1;
  Int_t    iEnd = -1;
  Int_t   iStep = 1;
  TString macro = "";
  
  // Interfacing command line:
  for(int l=1;l<argc;l++){
    TString arg = argv[l];
    
    if(arg.Contains(".C")) {
      macro = arg;
    } else if(arg.Contains("-t")) {
      char ss[2];
      sscanf(arg,"%2s%i",ss,&time);
    } else if(arg.Contains("-i")) {
      char ss[2];
      sscanf(arg,"%2s%i",ss,&iStart);
    } else if(arg.Contains("-f")) {
      char ss[2];
      sscanf(arg,"%2s%i",ss,&iEnd);
    } else if(arg.Contains("-s")) {
      char ss[2];
      sscanf(arg,"%2s%i",ss,&iStep);
    }
  }
  
  // Process macro string
  // Trick to get the macro and the simulation names
  char mname[50];
  char aux[2];
  char sname[50];
  sscanf(macro.Data(),"%[^\(]%2s%[^\"]",mname,aux,sname);
  
  // The Data manager
  // PData *pData = PData::Get(sname);
  // if(pData->isHiPACE()) {
  //   delete pData; pData = NULL;
  //   pData = PDataHiP::Get(sname);
  // }
  
  if(iStart<0) iStart = time;
  if(iEnd<=iStart) iEnd = iStart;

  // Time looper
  for(Int_t i=iStart; i<iEnd+1; i+=iStep) {

    time = i;

    cout << Form("\nLooping %s at time step %i:\n",sim.Data(),time) << endl;
    
    // Run macro
    TString command = Form(macro.Data(),time);
    command = Form(".x %s",command.Data());
    
    cout << command << endl;
    gROOT->ProcessLine(command);
    
  }


  // ---------------------------------------------------------
}
