#include <math.h>

#include <iostream>

#include <TROOT.h>
#include <TClass.h>
#include <TStyle.h>
#include <TColor.h>

#include "PPalette.hh"

ClassImp(PPalette);

PPalette::PPalette(const char * name)
  :TNamed(name, name)
{
  fNColors = 0;
  fColors = 0;
  fAlpha = 1.0;

  TObject * obj = gROOT->GetListOfSpecials()->First();
  while (obj) {
    if (strcmp(obj->GetName(), name)==0) {
      gROOT->GetListOfSpecials()->Remove(obj);
      delete obj;
      break;
    }
    obj = gROOT->GetListOfSpecials()->After(obj);
  }

  this->SetPalette("oli");

  gROOT->GetListOfSpecials()->Add(this);
}

PPalette::~PPalette()
{
  if (fColors) delete [] fColors;
}
  
Int_t PPalette::CreateGradientColorTable(UInt_t Number, 
					 Double_t* Stops, 
					 Double_t* Red, 
					 Double_t* Green,
					 Double_t* Blue,
					 UInt_t NColors,
					 Float_t alpha)
{

  fAlpha = alpha;
  
  if (fColors) delete [] fColors;
  fNColors = NColors;
  fColors = new Int_t[fNColors];
  
  Int_t ind = TColor::CreateGradientColorTable(Number,Stops,Red,Green,Blue,NColors,alpha);

  for (UInt_t i = 0; i < NColors; i++) fColors[i] = ind+i;

  return ind;
  
}

Int_t PPalette::CreateGradientColorTable(UInt_t Number, Double_t* Stops, Int_t *cindex, UInt_t NColors, Float_t alpha)
{

  Double_t *red = new Double_t[Number];
  Double_t *green = new Double_t[Number];
  Double_t *blue = new Double_t[Number];
  
  for (UInt_t i = 0; i < Number; i++) {
    Float_t r,g,b;
    gROOT->GetColor(cindex[i])->GetRGB(r,g,b);
    red[i] = r;
    green[i] = g;
    blue[i] = b;
  }
  
  Int_t ind = PPalette::CreateGradientColorTable(Number,Stops,red,green,blue,NColors,alpha);
  
  return ind;
}

Int_t PPalette::ChangeGradientColorTable(UInt_t Number, 
					 Double_t* Stops, 
					 Double_t* Red, 
					 Double_t* Green,
					 Double_t* Blue,
					 Float_t alpha)
{

  fAlpha = alpha;

  TColor::InitializeColors();
  
  if(Number < 2 || fNColors < 1){
    return -1;
  }
  
  UInt_t g, c;

  // Check if all RGB values are between 0.0 and 1.0 and
  // Stops goes from 0.0 to 1.0 in increasing order.
  for (c = 0; c < Number; c++) {
    if (Red[c] < 0 || Red[c] > 1.0 ||
	Green[c] < 0 || Green[c] > 1.0 ||
	Blue[c] < 0 || Blue[c] > 1.0 ||
	Stops[c] < 0 || Stops[c] > 1.0) {
      //Error("CreateGradientColorTable",
      //      "All RGB colors and stops have to be between 0.0 and 1.0");
      return -1;
    }
    if (c >= 1) {
      if (Stops[c-1] > Stops[c]) {
	//Error("CreateGradientColorTable",
	//      "Stops have to be in increasing order");
	return -1;
      }
    }
  }
  
  // Now create the new colors and replace the old ones
  UInt_t nPalette = 0;
  Int_t *palette = new Int_t[fNColors+1];
  UInt_t nColorsGradient;
  TColor *color;
  Int_t  index = fColors[0];

  for (g = 1; g < Number; g++) {
    // create the colors...
    nColorsGradient = (Int_t) (floor(fNColors*Stops[g]) - floor(fNColors*Stops[g-1]));
    for (c = 0; c < nColorsGradient; c++) {
      TColor *col = gROOT->GetColor(index);
      if(col) {
	col->SetRGB(Red[g-1] + c * (Red[g] - Red[g-1])/ nColorsGradient,
		    Green[g-1] + c * (Green[g] - Green[g-1])/ nColorsGradient,
		    Blue[g-1] + c * (Blue[g] - Blue[g-1])/ nColorsGradient);
	
	col->SetAlpha(alpha);

      } else {
	new TColor(index,
		   Red[g-1] + c * (Red[g] - Red[g-1])/ nColorsGradient,
		   Green[g-1] + c * (Green[g] - Green[g-1])/ nColorsGradient,
		   Blue[g-1] + c * (Blue[g] - Blue[g-1])/ nColorsGradient,
		   "  ",alpha);
      }
      
      palette[nPalette] = index;
      nPalette++;
      index++;
    }
  }
  
  TColor::SetPalette(nPalette, palette);  
  
  Int_t ind = index - fNColors; 
  for (UInt_t i = 0; i < fNColors; i++) fColors[i] = ind+i;

  return ind;

}

Int_t PPalette::ChangeGradientColorTable(UInt_t Number, Double_t* Stops, Int_t *cindex, Float_t alpha)
{

  Double_t *red = new Double_t[Number];
  Double_t *green = new Double_t[Number];
  Double_t *blue = new Double_t[Number];
  
  for (UInt_t i = 0; i < Number; i++) {
    Float_t r,g,b;
    gROOT->GetColor(cindex[i])->GetRGB(r,g,b);
    red[i] = r;
    green[i] = g;
    blue[i] = b;
  }
  
  Int_t ind = PPalette::ChangeGradientColorTable(Number,Stops,red,green,blue,alpha);
  
  return ind;
}


void PPalette::cd()
{
  gStyle->SetPalette(fNColors, fColors);
  //  gStyle->SetNumberContours(fNColors);
}
 
Int_t PPalette::SetPalette(const char * name)
{
  // TObject * obj = gROOT->GetListOfSpecials()->First();
  // while (obj) {
  //   if (strcmp(obj->GetName(), name)==0 && obj->IsA()->InheritsFrom(PPalette::Class())) {
  //     ((PPalette*)obj)->cd();
  //     break;
  //   }
  //   obj = gROOT->GetListOfSpecials()->After(obj);
  // }
  
  // if(obj) {
  //   if (fColors) delete [] fColors;
  //   fNColors = ((PPalette*)obj)->GetNColors();
  //   fColors = new Int_t[fNColors];
  //   for (UInt_t i = 0; i < fNColors; i++) fColors[i] = ((PPalette*)obj)->GetColorIndex(i);
    
  //   return 1;
  // } else {

    if(strcmp(name,"gray")==0) {
      const Int_t NRGBs = 2;
      const Int_t NCont = 255;
      Double_t Stops[NRGBs] = { 0.00, 1.00 };
      Double_t Red[NRGBs] =   { 0.99, 0.1 };
      Double_t Green[NRGBs] = { 0.99, 0.1 };
      Double_t Blue[NRGBs] =  { 0.99, 0.1 };
      if(!fColors)
	this->CreateGradientColorTable(NRGBs, Stops, Red, Green, Blue, NCont, fAlpha);
      else
	this->ChangeGradientColorTable(NRGBs, Stops, Red, Green, Blue, fAlpha);

      return 1;
      
    } else if(strcmp(name,"rbow0")==0) {
      const Int_t NRGBs = 6;
      const Int_t NCont = 255;
      Double_t Stops[NRGBs] = { 0.00, 0.40, 0.49, 0.51, 0.60, 1.00 };
      Double_t Red[NRGBs] =   { 0.106, 0.698, 1.0, 1.0, 0.965, 0.518 };
      Double_t Green[NRGBs] = { 0.078, 0.818, 1.0, 1.0, 0.925, 0.078 };
      Double_t Blue[NRGBs] =  { 0.518, 0.880, 1.0, 1.0, 0.353, 0.106 };
      if(!fColors)
	this->CreateGradientColorTable(NRGBs, Stops, Red, Green, Blue, NCont, fAlpha);
      else
	this->ChangeGradientColorTable(NRGBs, Stops, Red, Green, Blue, fAlpha);
      return 1;
      
    } else if(strcmp(name,"rbow")==0) {
      const Int_t NRGBs = 6;
      const Int_t NCont = 255;
      Double_t Stops[NRGBs] = { 0.00, 0.40, 0.49, 0.51, 0.60, 1.00 };
      Double_t Red[NRGBs] =   { 0.106, 0.698, 0.90, 0.90, 0.965, 0.518 };
      Double_t Green[NRGBs] = { 0.078, 0.818, 0.90, 0.90, 0.925, 0.078 };
      Double_t Blue[NRGBs] =  { 0.518, 0.880, 0.90, 0.90, 0.353, 0.106 };
      if(!fColors)
	this->CreateGradientColorTable(NRGBs, Stops, Red, Green, Blue, NCont, fAlpha);
      else
	this->ChangeGradientColorTable(NRGBs, Stops, Red, Green, Blue, fAlpha);
      return 1;
      
    } else if(strcmp(name,"electron0")==0) {
      const Int_t NRGBs = 6;
      const Int_t NCont = 255;
      Double_t Stops[NRGBs] = { 0.00, 0.10, 0.40, 0.55, 0.70, 1.00};
      Double_t Red[NRGBs] =   { 1.00, 0.52, 0.22, 0.39, 0.70, 1.00};
      Double_t Green[NRGBs] = { 1.00, 0.74, 0.34, 0.05, 0.20, 1.00};
      Double_t Blue[NRGBs] =  { 1.00, 0.80, 0.58, 0.33, 0.30, 0.20};
      if(!fColors)
	this->CreateGradientColorTable(NRGBs, Stops, Red, Green, Blue, NCont, fAlpha);
      else
	this->ChangeGradientColorTable(NRGBs, Stops, Red, Green, Blue, fAlpha);
      return 1;
      
    } else if(strcmp(name,"electron1")==0) {
      const Int_t NRGBs = 6;
      const Int_t NCont = 255;
      Double_t Stops[NRGBs] = { 0.00, 0.10, 0.40, 0.55, 0.70, 1.00};
      Double_t Red[NRGBs] =   { 0.00, 0.52, 0.22, 0.39, 0.70, 1.00};
      Double_t Green[NRGBs] = { 0.00, 0.74, 0.34, 0.05, 0.20, 1.00};
      Double_t Blue[NRGBs] =  { 0.00, 0.80, 0.58, 0.33, 0.30, 0.20};
      if(!fColors)
	this->CreateGradientColorTable(NRGBs, Stops, Red, Green, Blue, NCont, fAlpha);
      else
	this->ChangeGradientColorTable(NRGBs, Stops, Red, Green, Blue, fAlpha);
      return 1;
      
    } else if(strcmp(name,"electron")==0) {
      const Int_t NRGBs = 6;
      const Int_t NCont = 255;
      Double_t Stops[NRGBs] = { 0.00, 0.10, 0.40, 0.55, 0.70, 1.00};
      Double_t Red[NRGBs] =   { 0.90, 0.52, 0.22, 0.39, 0.70, 1.00};
      Double_t Green[NRGBs] = { 0.90, 0.74, 0.34, 0.05, 0.20, 1.00};
      Double_t Blue[NRGBs] =  { 0.90, 0.80, 0.58, 0.33, 0.30, 0.20};
      if(!fColors)
	this->CreateGradientColorTable(NRGBs, Stops, Red, Green, Blue, NCont, fAlpha);
      else
	this->ChangeGradientColorTable(NRGBs, Stops, Red, Green, Blue, fAlpha);
      return 1;
      
    } else if(strcmp(name,"elec0")==0) {
      const Int_t NRGBs = 5;
      const Int_t NCont = 255;
      Double_t Stops[NRGBs] = { 0.00, 0.40, 0.50, 0.60, 1.00};
      Double_t Red[NRGBs] =   { 1.00, 0.22, 0.39, 0.70, 1.00};
      Double_t Green[NRGBs] = { 1.00, 0.34, 0.05, 0.20, 1.00};
      Double_t Blue[NRGBs] =  { 1.00, 0.58, 0.33, 0.30, 0.20};
      if(!fColors)
	this->CreateGradientColorTable(NRGBs, Stops, Red, Green, Blue, NCont, fAlpha);
      else
	this->ChangeGradientColorTable(NRGBs, Stops, Red, Green, Blue, fAlpha);
      return 1;
      
    } else if(strcmp(name,"elec")==0) {
      const Int_t NRGBs = 5;
      const Int_t NCont = 255;
      Double_t Stops[NRGBs] = { 0.00, 0.40, 0.50, 0.60, 1.00};
      Double_t Red[NRGBs] =   { 0.90, 0.22, 0.39, 0.70, 1.00};
      Double_t Green[NRGBs] = { 0.90, 0.34, 0.05, 0.20, 1.00};
      Double_t Blue[NRGBs] =  { 0.90, 0.58, 0.33, 0.30, 0.20};
      if(!fColors)
	this->CreateGradientColorTable(NRGBs, Stops, Red, Green, Blue, NCont, fAlpha);
      else
	this->ChangeGradientColorTable(NRGBs, Stops, Red, Green, Blue, fAlpha);
      return 1;
      
    } else if(strcmp(name,"hot0")==0) {
      const Int_t NRGBs = 4;
      const Int_t NCont = 255;
      Double_t Stops[NRGBs] = { 0.00, 0.25, 0.75, 1.00 };
      Double_t Red[NRGBs] =   { 1.00, 1.000, 1.000, 1.000 };
      Double_t Green[NRGBs] = { 1.00, 0.149, 0.984, 1.000 };
      Double_t Blue[NRGBs] =  { 1.00, 0.000, 0.000, 1.000 };
      if(!fColors)
	this->CreateGradientColorTable(NRGBs, Stops, Red, Green, Blue, NCont, fAlpha);
      else
	this->ChangeGradientColorTable(NRGBs, Stops, Red, Green, Blue, fAlpha);
      return 1;
      
    } else if(strcmp(name,"hot1")==0) {
      const Int_t NRGBs = 4;
      const Int_t NCont = 255;
      Double_t Stops[NRGBs] = { 0.00, 0.25, 0.75, 1.00 };
      Double_t Red[NRGBs] =   { 0.90, 1.000, 1.000, 1.000 };
      Double_t Green[NRGBs] = { 0.90, 0.149, 0.984, 1.000 };
      Double_t Blue[NRGBs] =  { 0.90, 0.000, 0.000, 1.000 };
      if(!fColors)
	this->CreateGradientColorTable(NRGBs, Stops, Red, Green, Blue, NCont, fAlpha);
      else
	this->ChangeGradientColorTable(NRGBs, Stops, Red, Green, Blue, fAlpha);
      return 1;
      
    } else if(strcmp(name,"hot")==0) {
      const Int_t NRGBs = 3;
      const Int_t NCont = 255;
      Double_t Stops[NRGBs] = { 0.0, 0.75, 1.00 };
      Double_t Red[NRGBs] =   { 1.000, 1.000, 1.000 };
      Double_t Green[NRGBs] = { 0.149, 0.984, 1.000 };
      Double_t Blue[NRGBs] =  { 0.000, 0.000, 1.000 };
      if(!fColors)
	this->CreateGradientColorTable(NRGBs, Stops, Red, Green, Blue, NCont, fAlpha);
      else
	this->ChangeGradientColorTable(NRGBs, Stops, Red, Green, Blue, fAlpha);
      return 1;
      
    } else if(strcmp(name,"barsa")==0) {
      const Int_t NRGBs = 3;
      const Int_t NCont = 256;
      Double_t Stops[NRGBs] = {   0.00,  0.5, 1.00 }; 
      Double_t Red[NRGBs] =   { 0.0392, 1.00, 0.40 };
      Double_t Green[NRGBs] = { 0.1412, 1.00, 0.00 };
      Double_t Blue[NRGBs] =  { 0.4157, 1.00, 0.00 };
      if(!fColors)
	this->CreateGradientColorTable(NRGBs, Stops, Red, Green, Blue, NCont, fAlpha);
      else
	this->ChangeGradientColorTable(NRGBs, Stops, Red, Green, Blue, fAlpha);
      return 1;
      
    } else if(strcmp(name,"red0")==0) {
      const Int_t NRGBs = 3;
      const Int_t NCont = 255;
      Double_t Stops[NRGBs] = { 0.0, 0.10, 1.00 };
      Double_t Red[NRGBs] =   { 1.00, 0.965, 0.518 };
      Double_t Green[NRGBs] = { 1.00, 0.925, 0.078 };
      Double_t Blue[NRGBs] =  { 1.00, 0.353, 0.106 };
      if(!fColors)
	this->CreateGradientColorTable(NRGBs, Stops, Red, Green, Blue, NCont, fAlpha);
      else
	this->ChangeGradientColorTable(NRGBs, Stops, Red, Green, Blue, fAlpha);
      return 1;
      
    } else if(strcmp(name,"red")==0) {
      const Int_t NRGBs = 3;
      const Int_t NCont = 255;
      Double_t Stops[NRGBs] = { 0.0, 0.10, 1.00 };
      Double_t Red[NRGBs] =   { 0.90, 0.965, 0.518 };
      Double_t Green[NRGBs] = { 0.90, 0.925, 0.078 };
      Double_t Blue[NRGBs] =  { 0.90, 0.353, 0.106 };
      if(!fColors)
	this->CreateGradientColorTable(NRGBs, Stops, Red, Green, Blue, NCont, fAlpha);
      else
	this->ChangeGradientColorTable(NRGBs, Stops, Red, Green, Blue, fAlpha);
      return 1;
      
    } else if(strcmp(name,"blue")==0) {
      const Int_t NRGBs = 3;
      const Int_t NCont = 255;
      Double_t Stops[NRGBs] = { 0.0, 0.10, 1.00 };
      Double_t Red[NRGBs] =   { 0.90, 0.498, 0.106};
      Double_t Green[NRGBs] = { 0.90, 0.718, 0.078};
      Double_t Blue[NRGBs] =  { 0.90, 0.780, 0.518};
      if(!fColors)
	this->CreateGradientColorTable(NRGBs, Stops, Red, Green, Blue, NCont, fAlpha);
      else
	this->ChangeGradientColorTable(NRGBs, Stops, Red, Green, Blue, fAlpha);
      return 1;
      
    } else if(strcmp(name,"oli")==0) {
      const Int_t NRGBs = 9;
      const Int_t NCont = 255;
      Double_t Stops[NRGBs] = { 0.0000, 0.1250, 0.2500, 0.3750, 0.5000, 0.6250, 0.7500, 0.8750, 1.0000};
      Double_t Red[NRGBs]   = { 0.2081, 0.0591, 0.0779, 0.0231, 0.1801, 0.5300, 0.8185, 0.9955, 0.9763};
      Double_t Green[NRGBs] = { 0.1663, 0.3598, 0.5040, 0.6418, 0.7177, 0.7491, 0.7327, 0.7861, 0.9831};
      Double_t Blue[NRGBs]  = { 0.5292, 0.8683, 0.8384, 0.7913, 0.6424, 0.4661, 0.3498, 0.1967, 0.0538};
      if(!fColors)
	this->CreateGradientColorTable(NRGBs, Stops, Red, Green, Blue, NCont, fAlpha);
      else
	this->ChangeGradientColorTable(NRGBs, Stops, Red, Green, Blue, fAlpha);
      return 1;
      
    }  
    else {
      std::cout << Form(" PPalette :: palette name not recognised.") << std::endl;

      return 0;
    }
    
    //  }
  
}

  
void PPalette::Invert()
{
  if(!fColors) return; 

  Float_t *red = new Float_t[fNColors];
  Float_t *green = new Float_t[fNColors];
  Float_t *blue = new Float_t[fNColors];

  for(UInt_t i=0;i<fNColors;i++)
    gROOT->GetColor(fColors[i])->GetRGB(red[i],green[i],blue[i]);

  for(UInt_t i=0;i<fNColors;i++) {
    gROOT->GetColor(fColors[i])->SetRGB(red[fNColors-1-i],green[fNColors-1-i],blue[fNColors-1-i]);
  }
  
  
}

void PPalette::SetAlpha(Float_t alpha)
{

  fAlpha = alpha;

  if(fColors) {
    for(UInt_t i=0;i<fNColors;i++)
      gROOT->GetColor(fColors[i])->SetAlpha(fAlpha);
  }
  
}
