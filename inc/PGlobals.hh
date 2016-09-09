// global PLASMA style settings
#ifndef PGLOBALS_HH
#define PGLOBALS_HH

#include <iostream>

#include "TMath.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TSystem.h"
#include "TImage.h"
#include "TKey.h"
#include "TH1.h"
#include "TGraph.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TObjArray.h"
#include "TPave.h"
#include "TPaveText.h"

namespace PGlobals {

  // static const Int_t c_Canvas         = TColor::GetColor( "#f0f0f0" );
  // static const Int_t c_FrameFill      = TColor::GetColor( "#fffffd" );
  // static const Int_t c_TitleBox       = TColor::GetColor( "#5D6B7D" );
  // static const Int_t c_TitleBorder    = TColor::GetColor( "#7D8B9D" );
  // static const Int_t c_TitleText      = TColor::GetColor( "#FFFFFF" );
  // static const Int_t c_SignalFill     = TColor::GetColor( "#7d99d1" );
  // static const Int_t c_BackgroundLine = TColor::GetColor( "#ff0000" );
  // static const Int_t c_NovelBlue      = TColor::GetColor( "#2244a5" );

  static const Int_t elecLine = TColor::GetColor(69,108,155);
  static const Int_t elecFill = TColor::GetColor(245,245,245);
  static const Int_t fieldLine = kOrange+10;
  static const Int_t energyLine = kOrange+7;
  static const Int_t energyFill = energyLine;
  // static const Int_t frameColor = kGray+3; 
  // static const Int_t fontColor = kGray+3; 
  static const Int_t frameColor = kBlack; 
  static const Int_t fontColor = kBlack; 
  static const Int_t fontType = 43; 
  static const Int_t titleSize = 24; 
  static const Int_t labelSize = 22; 
  static const Int_t frameWidth = 2; 


  void SetPlasmaStyle();
 
  void SetH1Style(TH1 *histo, Int_t index=0);
  void SetH1LabelSize(TH1 *histo);
  void SetGraphStyle(TGraph *graph, Int_t index=0);
  void GetGraphMaxMin(TGraph *graph, Double_t &Min, Double_t &Max);
  void SetMultiClassStyle( TObjArray* hists ); 
  void SetFrameStyle( TH1* frame, Float_t scale = 1.0 );
  void SetPaveStyle(TPave* pave);  
  void SetPaveTextStyle( TPaveText* text, Int_t align=12);
  void CanvasPartition(TCanvas *C, Int_t N=3, 
		       Float_t lMargin = 0.15,
		       Float_t rMargin = 0.18,
		       Float_t bMargin = 0.10,
		       Float_t tMargin = 0.04,
		       Float_t mMargin = 0.0);
  
  void CanvasAsymPartition(TCanvas *C, Int_t N=3, 
			   Float_t lMargin = 0.15,
			   Float_t rMargin = 0.18,
			   Float_t bMargin = 0.10,
			   Float_t tMargin = 0.04,
			   Float_t factor = 1.0,
			   Float_t mMargin = 0.0);
  
  void CanvasDoublePartition(TCanvas *C,const Int_t Nx = 2,const Int_t Ny = 2,
			     Float_t lMargin = 0.15, Float_t rMargin = 0.05,
			     Float_t bMargin = 0.15, Float_t tMargin = 0.05,
			     Float_t vSpacing = 0.0);

  void DestroyCanvases();

  void Initialize( Bool_t usePlasmaStyle = kTRUE );
  
  TFile* OpenFile( const TString& fin );

  void imgconv( TCanvas* c, const TString & fname, const TString & opt="png");

  Int_t HCrossings(TH1F *h, Float_t *Cross, Float_t *Extr, Int_t MAXCROSS = 10, 
		   Float_t baseline = 0.0, Float_t xmax = -999, Float_t xmin = -999,
		   TString opt="");


}

#endif
