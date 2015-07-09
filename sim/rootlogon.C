{
  TString IncPath = gSystem->Getenv("PTOOLS");
  IncPath +=  "/inc";
  gROOT->ProcessLine(Form(".include %s",IncPath.Data()));

#include "PGlobals.hh"

#ifdef __CINT__  
  gSystem->Load("libptools.so");
#endif
  PGlobals::Initialize();

  // Palettes!
  gROOT->Macro("PPalettes.C");

}

