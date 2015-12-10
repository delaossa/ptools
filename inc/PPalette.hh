#ifndef PPalette_hh
#define PPalette_hh

#include <iostream>

#include <TNamed.h>

class PPalette : public TNamed
{
public:
  
  PPalette(const char * name);
  virtual ~PPalette();
  
  void                         CreateGradientColorTable(UInt_t Number, 
							Double_t* Stops, 
							Double_t* Red, 
							Double_t* Green,
							Double_t* Blue,
							UInt_t NColors,
							Float_t alpha=1);

  Int_t                        GetColor(UInt_t i=0) {
    if(i>=fNColors) {
      std::cout << "Color index out of range. Returning last color." << std::endl;
      return fColors[fNColors-1];
    }
    else if(i<0) {
      std::cout << "Color index out of range. Returning first color." << std::endl;
      return fColors[0];
    }
    return fColors[i];
  };

  UInt_t GetNColors() {
    return fNColors+1;
  }
  
  virtual void                 cd();
  static void                  SetPalette(const char * name);
  void                         SetAlpha(Float_t alpha=1);
  
protected:

  UInt_t                      fNColors;
  Int_t                      *fColors;

  ClassDef(PPalette,0);
};

#endif


