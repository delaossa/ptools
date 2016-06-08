#ifndef PPalette_hh
#define PPalette_hh

#include <iostream>

#include <TNamed.h>
#include <TColor.h>

class PPalette : public TNamed
{
public:
  
  PPalette(const char * name);
  virtual ~PPalette();
  
  Int_t   CreateGradientColorTable(UInt_t Number, 
				  Double_t* Stops, 
				  Double_t* Red, 
				  Double_t* Green,
				  Double_t* Blue,
				  UInt_t NColors,
				  Float_t alpha=1.0);

  Int_t   CreateGradientColorTable(UInt_t Number,
				   Double_t* Stops,
				   Int_t* cindex,
				   UInt_t NColors=64,
				   Float_t alpha=1.0);

  
  Int_t   ChangeGradientColorTable(UInt_t Number, 
				  Double_t* Stops, 
				  Double_t* Red, 
				  Double_t* Green,
				  Double_t* Blue,
				  Float_t alpha=1.0);

  Int_t   ChangeGradientColorTable(UInt_t Number,
				   Double_t* Stops,
				   Int_t* cindex,
				   Float_t alpha=1.0);

  
  Int_t  GetColor(UInt_t i=0) {
    if(i>=fNColors) {
      std::cout << "Color index out of range. Returning last color." << std::endl;
      return fColors[fNColors-1];
    }
    // else if(i<0) {
    //   std::cout << "Color index out of range. Returning first color." << std::endl;
    //   return fColors[0];
    // }
    return fColors[i];
  };

  UInt_t GetNColors() {
    return fNColors;
  }
  
  virtual void    cd();
  void    SetPalette(Int_t ncolors, Int_t *colors = 0,Float_t alpha=1.) {
    TColor::SetPalette(ncolors,colors,alpha);
    fNColors = TColor::GetNumberOfColors();
    Int_t ind = TColor::GetColorPalette(0);
    if (fColors) delete [] fColors;
    fColors = new Int_t[fNColors];
    for (UInt_t i = 0; i < fNColors; i++) fColors[i] = ind+i;
  }

  void            Invert();
  
  Int_t           SetPalette(const char * name);
  void            SetAlpha(Float_t alpha=1);

  Int_t           GetColorIndex(UInt_t i=0) {
    if(i<fNColors) {
      return fColors[i];
    } else
      return -1;
  }
  
protected:

  UInt_t                      fNColors;
  Int_t                      *fColors;
  Float_t                     fAlpha;         //Alpha (transparency)


  ClassDef(PPalette,0);
};

#endif


