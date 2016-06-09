R__LOAD_LIBRARY(libptools.so)
R__ADD_INCLUDE_PATH($PTOOLS/inc)

#include "PGlobals.hh"

void rootlogon()
{
  cout << "\n PTOOLS login script\n" << endl;
  PGlobals::Initialize();
}
