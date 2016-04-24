// global PLASMA style settings
#ifndef PGLOBALS
#define PGLOBALS

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

using std::cout;
using std::endl;

namespace PGlobals {

  // static Int_t c_Canvas         = TColor::GetColor( "#f0f0f0" );
  // static Int_t c_FrameFill      = TColor::GetColor( "#fffffd" );
  // static Int_t c_TitleBox       = TColor::GetColor( "#5D6B7D" );
  // static Int_t c_TitleBorder    = TColor::GetColor( "#7D8B9D" );
  // static Int_t c_TitleText      = TColor::GetColor( "#FFFFFF" );
  // static Int_t c_SignalFill     = TColor::GetColor( "#7d99d1" );
  // static Int_t c_BackgroundLine = TColor::GetColor( "#ff0000" );
  // static Int_t c_NovelBlue      = TColor::GetColor( "#2244a5" );

  static Int_t elecLine = TColor::GetColor(69,108,155);
  static Int_t elecFill = TColor::GetColor(245,245,245);
  static Int_t fieldLine = kOrange+10;
  static Int_t energyLine = kOrange+7;
  static Int_t energyFill = energyLine;
  // static Int_t frameColor = kGray+3; 
  // static Int_t fontColor = kGray+3; 
  static Int_t frameColor = kBlack; 
  static Int_t fontColor = kBlack; 
  static Int_t fontType = 43; 
  static Int_t titleSize = 28; 
  static Int_t labelSize = 22; 
  static Int_t frameWidth = 2; 


  void SetPlasmaStyle() {
      
    TStyle *PlasmaStyle = gROOT->GetStyle("Plasma");
    if(PlasmaStyle!=0) {
      //      cout << " style already defined." << endl;
      gROOT->SetStyle("Plasma");
      return;
    }

    PlasmaStyle = new TStyle(*gROOT->GetStyle("Plain")); // our style is based on Plain
    PlasmaStyle->SetName("Plasma");
    PlasmaStyle->SetTitle("Plasma style based on \"Plain\" with modifications defined in PlasmaStyle.C");
    gROOT->GetListOfStyles()->Add(PlasmaStyle);
    gROOT->SetStyle("Plasma");
			
    // Line styles
    PlasmaStyle->SetLineColor(kMagenta+2);
    //PlasmaStyle->SetLineStyleString( 5, "[52 12]" );
    //PlasmaStyle->SetLineStyleString( 6, "[22 12]" );
    //PlasmaStyle->SetLineStyleString( 7, "[22 10 7 10]" );

    // use plain black on white colors
    PlasmaStyle->SetCanvasBorderMode(0);
    PlasmaStyle->SetPadBorderMode(0);
    PlasmaStyle->SetPadColor(0);
    PlasmaStyle->SetFrameBorderMode(0);
    PlasmaStyle->SetFrameFillColor(0);
    PlasmaStyle->SetFrameLineColor(kBlack);
    PlasmaStyle->SetFrameLineWidth(1);

    // Global line width (for TAxis)
    PlasmaStyle->SetLineWidth(1);

    PlasmaStyle->SetGridStyle(2);
    PlasmaStyle->SetGridColor(kGray+1);
    PlasmaStyle->SetPadGridX(0);
    PlasmaStyle->SetPadGridY(0);
    
    PlasmaStyle->SetLegendBorderSize(0);
      
    PlasmaStyle->SetCanvasDefH(500);
    PlasmaStyle->SetCanvasDefW(800);
    // set the paper & margin sizes
    PlasmaStyle->SetPaperSize(TStyle::kA4);
    PlasmaStyle->SetPadTopMargin(0.05);
    PlasmaStyle->SetPadRightMargin(0.20);
    PlasmaStyle->SetPadBottomMargin(0.17);
    PlasmaStyle->SetPadLeftMargin(0.15);

    // Line and markers
    PlasmaStyle->SetMarkerStyle(21);
    PlasmaStyle->SetMarkerSize(0.8);
    PlasmaStyle->SetHistLineWidth(1);
    PlasmaStyle->SetHistLineColor(elecLine);
    PlasmaStyle->SetFuncColor(kMagenta+2);

    // do not display any of the standard histogram decorations
    PlasmaStyle->SetOptTitle(0);
    PlasmaStyle->SetTitleH(0.052);
    PlasmaStyle->SetOptStat(0);
    PlasmaStyle->SetOptFit(0);

    // put/quit tick marks on top and right of plots
    PlasmaStyle->SetPadTickX(0);
    PlasmaStyle->SetPadTickY(0);

    // Fonts
    PlasmaStyle->SetTextColor(fontColor);
    PlasmaStyle->SetTitleFont(fontType);
    PlasmaStyle->SetStatFont(fontType);
    PlasmaStyle->SetTextFont(fontType);
    PlasmaStyle->SetTitleFont(fontType,"xyz");
    PlasmaStyle->SetLabelFont(fontType,"xyz");

    // axis's title and label sizes
    PlasmaStyle->SetLabelOffset(0.008, "x");
    PlasmaStyle->SetLabelOffset(0.008, "y");
    PlasmaStyle->SetLabelOffset(0.008, "z");
    PlasmaStyle->SetLabelSize(0.05, "xy");
    PlasmaStyle->SetLabelSize(0.045, "z");
    PlasmaStyle->SetLabelColor(fontColor, "xyz");
    PlasmaStyle->SetTitleFontSize(titleSize);
    PlasmaStyle->SetTitleOffset(1.4,"x");
    PlasmaStyle->SetTitleOffset(1.2,"y");
    PlasmaStyle->SetTitleOffset(1.2,"z");
    PlasmaStyle->SetTitleSize(0.05, "xyz");
    PlasmaStyle->SetTitleColor(fontColor, "xyz");
    PlasmaStyle->SetNdivisions(505,"x");
    PlasmaStyle->SetNdivisions(505,"y");
    PlasmaStyle->SetNdivisions(505,"z");
  }

  // More styling functions
  // ---------------------------------------------
 
  void SetH1Style(TH1 *histo, Int_t index=0) {
    // Lines and colors for 1D Histos
    UInt_t Colors1D[]  = {kGray+2,elecLine,kOrange+8,kGray};
    UInt_t Colors1Db[] = {kGray+1,elecFill,kOrange+8,kGray};
    UInt_t Lines1D[]   = {2,2,2,2};
    UInt_t Fill1D[]    = {3002,3001,3001,3001};

    histo->SetLineColor(Colors1D[index%4]);
    histo->SetLineWidth(Lines1D[index%4]);
    histo->SetFillColor(Colors1Db[index%4]);
    histo->SetFillStyle(Fill1D[index%4]);

    histo->SetMarkerColor(Colors1D[index%4]);
    histo->SetMarkerStyle( 20 );
    histo->SetMarkerSize(0.8);

    histo->GetXaxis()->CenterTitle();
    histo->GetYaxis()->CenterTitle();
    
  }

  void SetH1LabelSize(TH1 *histo) {
    // Lines and colors for 1D Histos
       
    histo->SetLabelSize(0.04,"xyz");
    histo->SetLabelOffset(0.01,"xyz");
 
    histo->SetTitleSize(0.05,"xyz");
    histo->SetTitleOffset(1.2,"xy");
    histo->SetTitleOffset(1.3,"z");
    
    histo->SetTickLength(0.018,"xyz");
    
    histo->GetXaxis()->CenterTitle();
    histo->GetYaxis()->CenterTitle();
    histo->GetZaxis()->CenterTitle();

    histo->SetTitleFont(fontType,"xyz");
    histo->SetLabelFont(fontType,"xyz");
    
    histo->SetNdivisions(505,"xyz");
  }

  void SetGraphStyle(TGraph *graph, Int_t index=0) {
    // Lines and colors for 1D Histos
    UInt_t Colors1D[]  = {kGray+2,elecLine,kOrange+8,kGray};
    UInt_t Colors1Db[] = {kGray+1,elecFill,kOrange+8,kGray};
    UInt_t Lines1D[]   = {2,2,2,2};
    UInt_t Fill1D[]    = {3002,1001,3001,3001};

    graph->SetLineColor(Colors1D[index%4]);
    graph->SetLineWidth(Lines1D[index%4]);
    graph->SetFillColor(Colors1Db[index%4]);
    graph->SetFillStyle(Fill1D[index%4]);

    graph->SetMarkerColor(Colors1D[index%4]);
    graph->SetMarkerStyle( 20 );
    graph->SetMarkerSize(0.8);

    graph->GetXaxis()->CenterTitle();
    graph->GetYaxis()->CenterTitle();
    
  }

  void GetGraphMaxMin(TGraph *graph, Double_t &Min, Double_t &Max) {
    
    Min = 99999.;
    Max = -99999.;

    Int_t Npoints = graph->GetN();
    //  Double_t *x = graph->GetX();
    Double_t *y = graph->GetY();
    
    for(Int_t i=0;i<Npoints;i++) {
      if(y[i]>Max) Max = y[i];
      if(y[i]<Min) Min = y[i];
    }
    
  }


  void SetMultiClassStyle( TObjArray* hists ) 
  {
    Int_t FillColors[10] = {38,2,3,6,7,8,9,11};
    Int_t LineColors[10] = {4,2,3,6,7,8,9,11};
    Int_t FillStyles[5] = {1001,3554,3003,3545,0};

    for(Int_t i=0; i<hists->GetEntriesFast(); ++i){
      ((TH1*)(*hists)[i])->SetFillColor(FillColors[i%10]);
      ((TH1*)(*hists)[i])->SetFillStyle(FillStyles[i%5]);
      ((TH1*)(*hists)[i])->SetLineColor(LineColors[i%10]);
      ((TH1*)(*hists)[i])->SetLineWidth(2);
    }
  }

  // set frame styles
  void SetFrameStyle( TH1* frame, Float_t scale = 1.0 )
  {
    frame->SetLabelOffset( 0.012, "X" );// label offset on x axis
    frame->SetLabelOffset( 0.012, "Y" );// label offset on x axis
    frame->GetXaxis()->SetTitleOffset( 1.25 );
    frame->GetYaxis()->SetTitleOffset( 1.22 );
    frame->GetXaxis()->SetTitleSize( 0.045*scale );
    frame->GetYaxis()->SetTitleSize( 0.045*scale );
    Float_t labelSize = 0.04*scale;
    frame->GetXaxis()->SetLabelSize( labelSize );
    frame->GetYaxis()->SetLabelSize( labelSize );

    // global style settings
    gPad->SetTicks();
    gPad->SetLeftMargin  ( 0.108*scale );
    gPad->SetRightMargin ( 0.050*scale );
    gPad->SetBottomMargin( 0.120*scale  );
  }

  // set frame styles
  void SetPaveStyle(TPave* pave)
  {
    pave->SetBorderSize(0);
    pave->SetLineColor(0);
    pave->SetFillColor(0);
    pave->SetFillStyle(0);
    pave->SetShadowColor(0);
  }
  
  void SetPaveTextStyle( TPaveText* text, Int_t align=12)
  {
    SetPaveStyle(text);
    text->SetTextAlign(align);
  }

  void CanvasPartition(TCanvas *C, Int_t N=3, 
		       Float_t lMargin = 0.15,
		       Float_t rMargin = 0.18,
		       Float_t bMargin = 0.10,
		       Float_t tMargin = 0.04,
		       Float_t mMargin = 0.0)
  {
    if(!C) return;
    
    // Setup Pad layout:
    Float_t vStep = (1.- bMargin - tMargin - (N-1) * mMargin) / N;
    
    Float_t vposd = 0.0;
    Float_t vposu = 0.0;
    Float_t vmard = 0.0;
    Float_t vmaru = 0.0;
    Float_t vfactor = 0.0;
    Float_t hposl = 0.0;
    Float_t hposr = 1.0;
    Float_t hmarl = lMargin;
    Float_t hmarr = rMargin;
    
    for(Int_t i=0;i<N;i++) {
      C->cd(0);
      
      if(i==0) {
	vposd = 0.0;
	vposu = bMargin + vStep + 0.5 * mMargin;
	vfactor = vposu-vposd;  
	vmard = bMargin / vfactor;
	vmaru = 0.5 * mMargin / vfactor;
      } else if(i == N-1) {
	vposd = vposu;
	vposu = vposd + vStep + 0.5 * mMargin + tMargin;
	vfactor = vposu-vposd;   
	vmard = 0.5 * mMargin / vfactor;
	vmaru = tMargin / vfactor;
      } else {
	vposd = vposu;
	vposu = vposd + vStep + mMargin; 
	vfactor = vposu-vposd;
	vmard = 0.5 * mMargin / vfactor;
	vmaru = 0.5 * mMargin / vfactor;
      }    
            
      // cout << Form("pad %i : vpos = %.2f",i,vposu) << endl;

      char name[16];
      sprintf(name,"pad_%i",i);
      TPad *pad = (TPad*) gROOT->FindObject(name);
      if(pad) delete pad;
      pad = new TPad(name,"",hposl,vposd,hposr,vposu);
      pad->SetLeftMargin(hmarl);
      pad->SetRightMargin(hmarr);  
      pad->SetBottomMargin(vmard);
      pad->SetTopMargin(vmaru);
      pad->Draw();
    }
    
    // cout << "Canvas partitioned" <<endl;

  }

  void CanvasAsymPartition(TCanvas *C, Int_t N=3, 
			   Float_t lMargin = 0.15,
			   Float_t rMargin = 0.18,
			   Float_t bMargin = 0.10,
			   Float_t tMargin = 0.04,
			   Float_t factor = 1.0,
			   Float_t mMargin = 0.0)
  {
    if(!C) return;
    
    // Setup Pad layout:
    Float_t vStep = (1.- bMargin - tMargin - (N-1) * mMargin) / N;
    Float_t vStepB = N*vStep/(N-1+factor);
    Float_t vStepA = factor*vStepB;
    
    Float_t vposd = 0.0;
    Float_t vposu = 0.0;
    Float_t vmard = 0.0;
    Float_t vmaru = 0.0;
    Float_t vfactor = 0.0;
    Float_t hposl = 0.0;
    Float_t hposr = 1.0;
    Float_t hmarl = lMargin;
    Float_t hmarr = rMargin;
    
    for(Int_t i=0;i<N;i++) {
      C->cd(0);
      
      if(i==0) {
	vposd = 0.0;
	vposu = bMargin + vStepB + 0.5 * mMargin;
	vfactor = vposu-vposd;  
	vmard = bMargin / vfactor;
	vmaru = 0.5 * mMargin / vfactor;
      } else if(i == N-1) {
	vposd = vposu;
	vposu = vposd + vStepA + 0.5 * mMargin + tMargin;
	if(vposu>1) vposu = 1.0;
	vfactor = vposu-vposd;   
	vmard = 0.5 * mMargin / vfactor;
	vmaru = tMargin / vfactor;
      } else {
	vposd = vposu;
	vposu = vposd + vStepB + mMargin; 
	vfactor = vposu-vposd;
	vmard = 0.5 * mMargin / vfactor;
	vmaru = 0.5 * mMargin / vfactor;
      }    
            
      char name[16];
      sprintf(name,"pad_%i",i);
      TPad *pad = (TPad*) gROOT->FindObject(name);
      if(pad) delete pad;
      pad = new TPad(name,"",hposl,vposd,hposr,vposu);
      pad->SetLeftMargin(hmarl);
      pad->SetRightMargin(hmarr);  
      pad->SetBottomMargin(vmard);
      pad->SetTopMargin(vmaru);
      pad->Draw();
    }
    
  }
  
  void CanvasDoublePartition(TCanvas *C,const Int_t Nx = 2,const Int_t Ny = 2,
			     Float_t lMargin = 0.15, Float_t rMargin = 0.05,
			     Float_t bMargin = 0.15, Float_t tMargin = 0.05,
			     Float_t vSpacing = 0.0)
  {
    if(!C) return;
    
    // Setup Pad layout:
    Float_t vStep  = (1.- bMargin - tMargin - (Ny-1) * vSpacing) / Ny;
    Float_t hStep  = (1.- Nx * (lMargin + rMargin) ) / Nx;
    
    Float_t vposd = 0.0;
    Float_t vposu = 0.0;
    Float_t vmard = 0.0;
    Float_t vmaru = 0.0;
    Float_t vfactor = 0.0;
    Float_t hposl = 0.0;
    Float_t hposr = 1.0;
    Float_t hmarl = lMargin;
    Float_t hmarr = rMargin;
    Float_t hfactor = 0.0;

    for(Int_t i=0;i<Nx;i++) {

      if(i==0) {
	hposl = 0.0;
	hposr = lMargin + hStep + rMargin;
	hfactor = hposr-hposl;  
	hmarl = lMargin / hfactor;
	hmarr = rMargin / hfactor;
      } else {
	hposl = hposr;
	hposr = hposl + lMargin + hStep + rMargin;
	hfactor = hposr-hposl;   
	hmarl = lMargin / hfactor;;
	hmarr = rMargin / hfactor;
      }        
    
      for(Int_t j=0;j<Ny;j++) {
      
	if(j==0) {
	  vposd = 0.0;
	  vposu = bMargin + vStep;
	  vfactor = vposu-vposd;  
	  vmard = bMargin / vfactor;
	  vmaru = 0.0;
	} else if(j == Ny-1) {
	  vposd = vposu + vSpacing;
	  vposu = vposd + vStep + tMargin;
	  vfactor = vposu-vposd;   
	  vmard = 0.0;
	  vmaru = tMargin / (vposu-vposd);
	} else {
	  vposd = vposu + vSpacing;
	  vposu = vposd + vStep; 
	  vfactor = vposu-vposd;
	  vmard = 0.0;
	  vmaru = 0.0;
	}    
         
	C->cd(0);
      
	//      cout << Form(" %i  %i : (%.2f,%.2f,%.2f,%.2f)",i,j,hposl,vposd,hposr,vposu) << endl;

	char name[16];
	sprintf(name,"pad_%i_%i",i,j);
	TPad *pad = (TPad*) gROOT->FindObject(name);
	if(pad) delete pad;
	pad = new TPad(name,"",hposl,vposd,hposr,vposu);
	pad->SetLeftMargin(hmarl);
	pad->SetRightMargin(hmarr);  
	pad->SetBottomMargin(vmard);
	pad->SetTopMargin(vmaru);

	pad->SetFrameBorderMode(0);
	pad->SetBorderMode(0);
	pad->SetBorderSize(0);
     
	pad->Draw();
      }
    }
  }

  void DestroyCanvases()
  {

    TList* loc = (TList*)gROOT->GetListOfCanvases();
    TListIter itc(loc);
    TObject *o(0);
    while ((o = itc())) delete o;
  }

  // set style and remove existing canvas'
  void Initialize( Bool_t usePlasmaStyle = kTRUE )
  {
    // destroy canvas'
    DestroyCanvases();

    // set style
    if (!usePlasmaStyle) {
      gROOT->SetStyle("Plain");
      gStyle->SetOptStat(0);
      return;
    }

    SetPlasmaStyle();
  }

  // checks if file with name "fin" is already open, and if not opens one
  TFile* OpenFile( const TString& fin )
  {
    TFile* file = gDirectory->GetFile();
    if (file==0 || fin != file->GetName()) {
      if (file != 0) {
	gROOT->cd();
	file->Close();
      }
      cout << "--- Opening root file " << fin << " in read mode" << endl;
      file = TFile::Open( fin, "READ" );
    }
    else {
      file = gDirectory->GetFile();
    }

    file->cd();
    return file;
  }

  // used to create output file for canvas
  void imgconv( TCanvas* c, const TString & fname, const TString & opt="png")
  {
    // return;
    if (NULL == c) {
      cout << "*** Error in PlasmaGlob::imgconv: canvas is NULL" << endl;
    }
    else {
      // create directory if not existing
      TString dir1,dir2;
      TString f = fname;
      if(f.Contains('/'))
	dir2 = f.Remove( f.Last( '/' ), f.Length() - f.Last( '/' ) );
      if(f.Contains('/')) {
	dir1 = f.Remove( f.Last( '/' ), f.Length() - f.Last( '/' ) );	
      }
      if(!dir1.IsNull())
	gSystem->mkdir( dir1 );
      if(!dir2.IsNull())
	gSystem->mkdir( dir2 );
      
      c->cd();
      
      if(opt.Contains("png"))
	c->Print(fname + ".png");

      if(opt.Contains("jpg"))
	c->Print(fname + ".jpg");

      if(opt.Contains("tiff"))
	c->Print(fname + ".tiff","tiff");
      
      if(opt.Contains("eps"))
	c->Print(fname + ".eps");
      
      if(opt.Contains("pdf"))  
	c->Print(fname + ".pdf");
   
      if(opt.Contains("svg"))  
	c->Print(fname + ".svg");
    }
  }

  Int_t HCrossings(TH1F *h, Float_t *Cross, Float_t *Extr, Int_t MAXCROSS = 10, 
		    Float_t baseline = 0.0, Float_t xmax = -999, Float_t xmin = -999, TString opt="") {
  
    // Calculate the crossings and the extremes of an histogram.
    for(Int_t ic=0;ic<MAXCROSS;ic++) {
      Cross[ic] = 0.0;
      Extr[ic] = 0.0;
    }
    Int_t Ncross = 0;

    if(!h) return Ncross;
  
    Float_t XMAX = h->GetXaxis()->GetXmax();
    Float_t XMIN = h->GetXaxis()->GetXmin();
  
    // Finds crossings from right to left, between xmax and xmin 
    if( (xmin == -999) || (xmin<XMIN) ) xmin = XMIN;
    if( (xmax == -999) || (xmax>XMAX) ) xmax = XMAX;
  
    for(Int_t ip=h->GetNbinsX();ip>1;ip--) {
    
      Float_t X1 = h->GetBinCenter(ip);
      if(X1 > xmax) continue;
      if(X1 < xmin) continue;
      Float_t X2 = h->GetBinCenter(ip-1);
      Float_t E1 = h->GetBinContent(ip)-baseline;
      Float_t E2 = h->GetBinContent(ip-1)-baseline;
    
      if(E1*E2 >= 0) { // No change of sign means we are in a side of the zero axis.
	if(fabs(E2 + baseline)>fabs(Extr[Ncross])) {
	  Extr[Ncross] = E1 + baseline;
	} 
      } else {  // change of sign means a crossing!
      
	Float_t zcross =  -E1 * ( (X2-X1)/(E2-E1) ) + X1;
            
	// add the point
	Cross[Ncross] = zcross;
	Ncross++;
 
	if(opt.Contains("cum") && Ncross==2) baseline = Extr[1];
      
	if(Ncross>=MAXCROSS) break;
      }
    }
    return Ncross;
  } 

  
}

#endif
