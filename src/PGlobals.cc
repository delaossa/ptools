#include "PGlobals.hh"
#include "TPad.h"

using std::cout;
using std::endl;

void PGlobals::SetPlasmaStyle() {
      
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
  PlasmaStyle->SetLineColor(kBlack);
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

  PlasmaStyle->SetLineWidth(1);

  // PlasmaStyle->SetGridStyle(2);
  PlasmaStyle->SetGridColor(kGray+1);
  PlasmaStyle->SetPadGridX(0);
  PlasmaStyle->SetPadGridY(0);
    
  PlasmaStyle->SetLegendBorderSize(0);
      
  PlasmaStyle->SetCanvasDefW(800);
  PlasmaStyle->SetCanvasDefH(600);
  // set the paper & margin sizes
  PlasmaStyle->SetPaperSize(TStyle::kA4);
  PlasmaStyle->SetPadTopMargin(0.05);
  // PlasmaStyle->SetPadRightMargin(0.20);
  PlasmaStyle->SetPadRightMargin(0.10);
  PlasmaStyle->SetPadBottomMargin(0.17);
  PlasmaStyle->SetPadLeftMargin(0.15);

  // Line and markers
  PlasmaStyle->SetMarkerStyle(20);
  PlasmaStyle->SetMarkerSize(0.4);
  PlasmaStyle->SetMarkerColor(kAzure+2);
  PlasmaStyle->SetHistLineWidth(1);
  PlasmaStyle->SetHistLineColor(elecLine);
  PlasmaStyle->SetFuncColor(kMagenta+2);

  // do not display any of the standard histogram decorations
  PlasmaStyle->SetOptTitle(0);
  PlasmaStyle->SetOptStat(0);
  PlasmaStyle->SetOptFit(0);

  // add/remove tick marks on top and right of plots
  PlasmaStyle->SetPadTickX(0);
  PlasmaStyle->SetPadTickY(0);

  // Fonts
  PlasmaStyle->SetStatFont(fontType);
  PlasmaStyle->SetTextFont(fontType);
    
  PlasmaStyle->SetTitleFont(fontType,"xyz");
  PlasmaStyle->SetLabelFont(fontType,"xyz");

  PlasmaStyle->SetTextColor(fontColor);
  PlasmaStyle->SetTitleColor(fontColor, "xyz");
  PlasmaStyle->SetLabelColor(fontColor, "xyz");

  // axis's title and label sizes
  PlasmaStyle->SetLabelSize(labelSize, "xyz");
  PlasmaStyle->SetTitleSize(titleSize, "xyz");
    
  PlasmaStyle->SetLabelOffset(0.01, "xyz");
  PlasmaStyle->SetTitleOffset(1.1,"xy");
  PlasmaStyle->SetTitleOffset(1.2,"z");
  PlasmaStyle->SetTickLength(0.018,"x");
  PlasmaStyle->SetTickLength(0.01,"yz");
  PlasmaStyle->SetNdivisions(505,"xyz");


  // palettes
  gROOT->Macro("PPalettes.C");
}
  
// More styling functions
// ---------------------------------------------
 
void PGlobals::SetH1Style(TH1 *histo, Int_t index) {
  // Lines and colors for 1D Histos
  Int_t Colors1D[4]  = {kGray+2,elecLine,kOrange+8,kGray};
  Int_t Colors1Db[4] = {kGray+1,elecFill,kOrange+8,kGray};
  Int_t Lines1D[4]   = {2,2,2,2};
  Int_t Fill1D[4]    = {3002,3001,3001,3001};

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

void PGlobals::SetH1LabelSize(TH1 *histo) {
  // Lines and colors for 1D Histos
       
  histo->SetLabelSize(labelSize,"xyz");
  histo->SetLabelOffset(0.01,"xyz");
 
  histo->SetTitleSize(titleSize,"xyz");
  histo->SetTitleOffset(1.4,"xyz");
    
  histo->SetTickLength(0.018,"x");
  histo->SetTickLength(0.01,"yz");
  
  histo->GetXaxis()->CenterTitle();
  histo->GetYaxis()->CenterTitle();
  histo->GetZaxis()->CenterTitle();

  histo->SetTitleFont(fontType,"xyz");
  histo->SetLabelFont(fontType,"xyz");
    
  histo->SetNdivisions(505,"xyz");
}

void PGlobals::SetGraphStyle(TGraph *graph, Int_t index) {
  // Lines and colors for 1D Histos
  Int_t Colors1D[]  = {kGray+2,elecLine,kOrange+8,kGray};
  Int_t Colors1Db[] = {kGray+1,elecFill,kOrange+8,kGray};
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

void PGlobals::GetGraphMaxMin(TGraph *graph, Double_t &Min, Double_t &Max) {
    
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


void PGlobals::SetMultiClassStyle( TObjArray* hists ) 
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
void PGlobals::SetFrameStyle( TH1* frame, Float_t scale)
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
void PGlobals::SetPaveStyle(TPave* pave)
{
  pave->SetBorderSize(0);
  pave->SetLineColor(0);
  pave->SetFillColor(0);
  pave->SetFillStyle(0);
  pave->SetShadowColor(0);
}
  
void PGlobals::SetPaveTextStyle( TPaveText* text, Int_t align)
{
  SetPaveStyle(text);
  text->SetTextAlign(align);
}

void PGlobals::CanvasPartition(TCanvas *C, Int_t N, 
			       Double_t lMargin,
			       Double_t rMargin,
			       Double_t bMargin,
			       Double_t tMargin,
			       Double_t mMargin)
{
  if(!C) return;
    
  // Setup Pad layout:
  Double_t vStep = (1.- bMargin - tMargin - (N-1) * mMargin) / N;
    
  Double_t vposd = 0.0;
  Double_t vposu = 0.0;
  Double_t vmard = 0.0;
  Double_t vmaru = 0.0;
  Double_t vfactor = 0.0;
  Double_t hposl = 0.0;
  Double_t hposr = 1.0;
  Double_t hmarl = lMargin;
  Double_t hmarr = rMargin;
    
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

    cout << Form(" Pad: %s created",name) << endl;
  }
    
  // cout << "Canvas partitioned" <<endl;

}

void PGlobals::CanvasAsymPartition(TCanvas *C, Int_t N, 
				   Double_t lMargin,
				   Double_t rMargin,
				   Double_t bMargin,
				   Double_t tMargin,
				   Double_t factor,
				   Double_t mMargin)
{
  if(!C) return;
    
  // Setup Pad layout:
  Double_t vStep = (1.- bMargin - tMargin - (N-1) * mMargin) / N;
  Double_t vStepB = N*vStep/(N-1+factor);
  Double_t vStepA = factor*vStepB;
    
  Double_t vposd = 0.0;
  Double_t vposu = 0.0;
  Double_t vmard = 0.0;
  Double_t vmaru = 0.0;
  Double_t vfactor = 0.0;
  Double_t hposl = 0.0;
  Double_t hposr = 1.0;
  Double_t hmarl = lMargin;
  Double_t hmarr = rMargin;
    
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
  
void PGlobals::CanvasDoublePartition(TCanvas *C,const Int_t Nx,const Int_t Ny,
				     Double_t lMargin, Double_t rMargin,
				     Double_t bMargin, Double_t tMargin,
				     Double_t vSpacing)
{
  if(!C) return;
    
  // Setup Pad layout:
  Double_t vStep  = (1.- bMargin - tMargin - (Ny-1) * vSpacing) / Ny;
  Double_t hStep  = (1.- Nx * (lMargin + rMargin) ) / Nx;
    
  Double_t vposd = 0.0;
  Double_t vposu = 0.0;
  Double_t vmard = 0.0;
  Double_t vmaru = 0.0;
  Double_t vfactor = 0.0;
  Double_t hposl = 0.0;
  Double_t hposr = 1.0;
  Double_t hmarl = lMargin;
  Double_t hmarr = rMargin;
  Double_t hfactor = 0.0;

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

void PGlobals::DestroyCanvases()
{

  TList* loc = (TList*)gROOT->GetListOfCanvases();
  TListIter itc(loc);
  TObject *o(0);
  while ((o = itc())) delete o;
}

// set style and remove existing canvas'
void PGlobals::Initialize( Bool_t usePlasmaStyle )
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
TFile* PGlobals::OpenFile( const TString& fin )
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

// create directory if not existing
void PGlobals::mkdir(const TString & fname) {
  TString dir1,dir2,dir3;
  TString f = fname;
  if(f.Contains('/'))
    dir3 = f.Remove( f.Last( '/' ), f.Length() - f.Last( '/' ) );
  if(f.Contains('/'))
    dir2 = f.Remove( f.Last( '/' ), f.Length() - f.Last( '/' ) );
  if(f.Contains('/')) 
    dir1 = f.Remove( f.Last( '/' ), f.Length() - f.Last( '/' ) );	
  
  if(!dir1.IsNull())
    gSystem->mkdir( dir1 );
  if(!dir2.IsNull())
    gSystem->mkdir( dir2 );
  if(!dir3.IsNull())
    gSystem->mkdir( dir3 );
}


// used to create output file for canvas
void PGlobals::imgconv( TCanvas* c, const TString & fname, const TString & opt)
{
  // return;
  if (NULL == c) {
    cout << "*** Error in PlasmaGlob::imgconv: canvas is NULL" << endl;
  }
  else {

    PGlobals::mkdir(fname);
    
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

Int_t PGlobals::HCrossings(TH1F *h, Float_t *Cross, Float_t *Extr, Int_t MAXCROSS, 
			   Float_t baseline, Float_t xmax, Float_t xmin, TString opt) {
  
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
	Extr[Ncross] = E2 + baseline;
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
