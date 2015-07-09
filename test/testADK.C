
#ifndef __CINT__
#include "PFunctions.hh"
#endif

#include "PData.hh"
#include "PlasmaGlob.hh"

void testADK(const TString &opt=""){

#ifdef __CINT__
  gSystem->Load("libplasma.so");
#endif
  PlasmaGlob::Initialize();  

  // Helium
  Double_t Z = 1;       // First electron
  Double_t Ei = 24.59*PUnits::eV;  // eV
  Double_t E = 92.75*(PUnits::GV/PUnits::m);     // GV/m

  cout << Form(" First  He level at E = %.3f  -->  W = %8.6f fs^-1 ",E/(PUnits::GV/PUnits::m),PFunc::ADK_ENG(E,Ei,Z)*PUnits::femtosecond) << endl;


  // Helium (2nd level)
  Z = 2;       // First electron
  Ei = 54.4*PUnits::eV;  // eV
  E = 234.96*(PUnits::GV/PUnits::m);     // GV/m

  cout << Form(" Second He level at E = %.3f  -->  W = %8.6f fs^-1 ",E/(PUnits::GV/PUnits::m),PFunc::ADK_ENG(E,Ei,Z)*PUnits::femtosecond) << endl;

  // Hydrogen
  Z = 1;
  Ei = 13.6*PUnits::eV;  // eV
  E = 33.8*(PUnits::GV/PUnits::m);     // GV/m

  cout << Form(" Second He level at E = %.3f  -->  W = %8.6f fs^-1 ",E/(PUnits::GV/PUnits::m),PFunc::ADK_ENG(E,Ei,Z)*PUnits::femtosecond) << endl;


}
