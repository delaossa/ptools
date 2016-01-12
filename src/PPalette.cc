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

  TObject * obj = gROOT->GetListOfSpecials()->First();
  while (obj) {
    if (strcmp(obj->GetName(), name)==0) {
      gROOT->GetListOfSpecials()->Remove(obj);
      delete obj;
      break;
    }
    obj = gROOT->GetListOfSpecials()->After(obj);
  }

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
 
  if (fColors) delete [] fColors;
  fNColors = NColors;
  fColors = new Int_t[fNColors];
  
  Int_t ind = TColor::CreateGradientColorTable(Number,Stops,Red,Green,Blue,NColors,alpha);

  for (UInt_t i = 0; i < NColors; i++) fColors[i] = ind+i;

  return ind;
  
}

Int_t PPalette::ChangeGradientColorTable(UInt_t Number, 
					Double_t* Stops, 
					Double_t* Red, 
					Double_t* Green,
					Double_t* Blue,
					Float_t alpha)
{

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

void PPalette::cd()
{
  gStyle->SetPalette(fNColors, fColors);
  // gStyle->SetNumberContours(fNColors);
}
 
Int_t PPalette::SetPalette(const char * name)
{
  TObject * obj = gROOT->GetListOfSpecials()->First();
  while (obj) {
    if (strcmp(obj->GetName(), name)==0 && obj->IsA()->InheritsFrom(PPalette::Class())) {
      ((PPalette*)obj)->cd();
      break;
    }
    obj = gROOT->GetListOfSpecials()->After(obj);
  }
  
  if(obj) {
    if (fColors) delete [] fColors;
    fNColors = ((PPalette*)obj)->GetNColors();
    fColors = new Int_t[fNColors];
    for (UInt_t i = 0; i < fNColors; i++) fColors[i] = ((PPalette*)obj)->GetColorIndex(i);
    
    return 1;
  } else {

    if(strcmp(name,"gray")==0) {
      const Int_t NRGBs = 2;
      const Int_t NCont = 64;
      Double_t Stops[NRGBs] = { 0.00, 1.00 };
      Double_t Red[NRGBs] =   { 0.99, 0.1 };
      Double_t Green[NRGBs] = { 0.99, 0.1 };
      Double_t Blue[NRGBs] =  { 0.99, 0.1 };
     
      this->CreateGradientColorTable(NRGBs, Stops, Red, Green, Blue, NCont);
      return 1;
      
    } else    
      return 0;
    
    
  }
  
}

  


void PPalette::SetAlpha(Float_t alpha)
{
  for(UInt_t i=0;i<fNColors;i++)
    gROOT->GetColor(fColors[i])->SetAlpha(alpha);
  
}
