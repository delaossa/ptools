#ifndef PDATAHIP
#define PDATAHIP

#include "PData.hh"

// PData class definition
// -----------------------------------------------------------------

class PDataHiP : public PData 
{
public:

  PDataHiP(const char * name);
  virtual ~PDataHiP();

  // void Clear(Option_t *option="");
  static PDataHiP *  Get(const char * name = 0);


  void LoadFileNames(Int_t t);

  ClassDef(PDataHiP,1) // Singleton class for global configuration settings
};

#endif
