#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>

#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TCut.h>
#include <TChain.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TPaletteAxis.h>
#include <TExec.h>

#include "PData.hh"
#include "PGlobals.hh"
#include "PPalette.hh"
#include "H5Cpp.h"

using namespace std;
using namespace H5;


int main(int argc,char *argv[]) {
  if(argc<=2) {
    printf("\n Usage: %s <simulation name> <-t(time)>\n",argv[0]);
    printf("      <-i(initial time)> <-f(final time)> <-s(time step)>\n");
    printf("      <--center> <--comov> <--units> <--logz>\n");
    printf("      <--png> <--pdf> <--eps> \n");
    printf("      <--file> <--loop>\n");
    return 0;
  }

  // General options
  TString   sim = "";
  Int_t    time = 0;
  Int_t  iStart = -1;
  Int_t    iEnd = -1;
  Int_t   iStep = 1;
  TString   opt = "";
  
  // Options for Spectrum
  Float_t Pmin =  99999.;
  Float_t Pmax = -99999.;

  // Interfacing command line:
  for(int l=1;l<argc;l++){
    TString arg = argv[l];

    if(arg.Contains("--pdf")) {
      opt += "pdf";
    } else if(arg.Contains("--eps")) {
      opt += "eps";
    } else if(arg.Contains("--png")) {
      opt += "png";
    } else if(arg.Contains("--tiff")) {
      opt += "tiff";
    } else if(arg.Contains("--comov")){
      opt += "comov";
    } else if(arg.Contains("--units")){
      opt += "units";
    } else if(arg.Contains("--center")){
      opt += "center";
    } else if(arg.Contains("--grid")){
      opt += "grid"; 
    } else if(arg.Contains("--logz")){
      opt += "logz"; 
    } else if(arg.Contains("--split")){
      opt += "split"; 
    } else if(arg.Contains("--loop")){
      opt += "loop"; 
    } else if(arg.Contains("--file")){
      opt += "file"; 
    } else if(arg.Contains("--notext")){
      opt += "notext"; 
    } else if(arg.Contains("--noinfo")){
      opt += "noinfo"; 
    } else if(arg.Contains("-t")) {
      char ss[2];
      sscanf(arg,"%2s%i",ss,&time);
    } else if(arg.Contains("-i")) {
      char ss[2];
      sscanf(arg,"%2s%i",ss,&iStart);
    } else if(arg.Contains("-f")) {
      char ss[2];
      sscanf(arg,"%2s%i",ss,&iEnd);
    } else if(arg.Contains("-s")) {
      char ss[2];
      sscanf(arg,"%2s%i",ss,&iStep);
    } else if(arg.Contains("-pmin")) {
      char ss[5];
      sscanf(arg,"%5s%f",ss,&Pmin);
    } else if(arg.Contains("-pmax")) {
      char ss[5];
      sscanf(arg,"%5s%f",ss,&Pmax);
    } else {
      sim = arg;
    }
  }
  
  PGlobals::Initialize();

  if(opt.Contains("grid")) {
    gStyle->SetPadGridX(1);
    gStyle->SetPadGridY(1);
  }

  // Load PData
  PData *pData = PData::Get(sim.Data());

  if(iStart<0) iStart = time;
  if(iEnd<=iStart) iEnd = iStart;
  
  // Some plasma constants
  Float_t n0 = pData->GetPlasmaDensity();
  Float_t kp = pData->GetPlasmaK();
  Float_t skindepth = 1.;
  if(kp!=0.0) skindepth = 1/kp;
  Float_t E0 = pData->GetPlasmaE0();

  // z start of the plasma in normalized units.
  Float_t zStartPlasma = pData->GetPlasmaStart()*kp;
  // z start of the beam in normalized units.
  Float_t zStartBeam = pData->GetBeamStart()*kp;
  
  opt += "comovcenter";
  
  // Time looper
  for(Int_t i=iStart; i<iEnd+1; i+=iStep) {

    time = i;
    pData->LoadFileNames(time);    
    
    if(!pData->IsInit()) continue;
    
    // Time in OU
    Float_t Time = pData->GetRealTime();

    // Centering time and z position:
    Float_t shiftz = pData->Shift(opt);

    if(opt.Contains("center")) {
      Time -= zStartPlasma;
      if(opt.Contains("comov"))      // Centers on the head of the beam.
	Time += zStartBeam;
    }

    // ----------------------------------------------------------------------------------
  
    Int_t Nspecies = pData->NSpecies();
    Float_t **ene = new Float_t*[Nspecies];
    Float_t **x1 = new Float_t*[Nspecies];
    Float_t **q = new Float_t*[Nspecies];
    UInt_t *NP = new UInt_t[Nspecies];
    Float_t eneMin = 99999;
    Float_t eneMax = -99999;
    for(Int_t i=0;i<Nspecies;i++) {
      ene[i] = NULL;
      x1[i] = NULL;
      NP[i] = 0;
    
      if(!pData->GetRawFileName(i)) continue;
      NP[i] = pData->GetRawSingleArray(pData->GetRawFileName(i)->c_str(),&x1[i],"x1");
      NP[i] = pData->GetRawSingleArray(pData->GetRawFileName(i)->c_str(),&ene[i],"ene");
      NP[i] = pData->GetRawSingleArray(pData->GetRawFileName(i)->c_str(),&q[i],"q");
    
      // for(UInt_t j=0;j<1000;j++)
      // 	cout << Form(" %6f  %6f  %6f",x1[i][j],ene[i][j],q[i][j]) << endl;
    
      for(UInt_t j=0;j<NP[i];j++) {
      	if(ene[i][j]<eneMin) eneMin = ene[i][j];
      	if(ene[i][j]>eneMax) eneMax = ene[i][j];
      }
    }    

    // BOX limits
    Double_t X1MIN = pData->GetXMin(0) - shiftz;
    Double_t X1MAX = pData->GetXMax(0) - shiftz;
    Double_t dx1 = pData->GetDX(0);    
    Float_t x1Min = pData->GetX1Min() - shiftz;
    Float_t x1Max = pData->GetX1Max() - shiftz;
    Int_t Nx1Bin  = pData->GetX1N();
    x1Min -= 0.4*(x1Max - x1Min);

    // Adjust binning
    x1Min = floor((x1Min-X1MIN)/dx1) * dx1 + X1MIN;  
    x1Max = floor((x1Max-X1MIN)/dx1) * dx1 + X1MIN;
    Nx1Bin  = ceil ((x1Max - x1Min)/(dx1));
    
    Float_t eMin = eneMin - (eneMax-eneMin) * 0.4;
    Float_t eMax = eneMax + (eneMax-eneMin) * 0.1;
    Int_t NeBin  = Nx1Bin;
    
    char hName[24];
    TH2F **hEneVsZ  = new TH2F*[Nspecies];
    TH2F *hEneVsZjoint = NULL;
    TH1D **hCurr = new TH1D*[Nspecies];
    TH1D **hSpec = new TH1D*[Nspecies];
    TGraph **gSpec = new TGraph*[Nspecies];   
    Float_t maxCur = -9999;
    Float_t maxEne = -9999;
    for(Int_t i=0;i<Nspecies;i++) {
    
      hEneVsZ[i] = NULL;
      hSpec[i] = NULL;
      hCurr[i] = NULL;
      gSpec[i] = NULL;
      
      if(ene[i]==NULL) continue;
    
      sprintf(hName,"EneVsZ_%i",i); 
      hEneVsZ[i] = new TH2F(hName,"",Nx1Bin,x1Min,x1Max,NeBin,eMin,eMax);

      for(UInt_t j=0;j<NP[i];j++) {
	hEneVsZ[i]->Fill(x1[i][j]-shiftz,ene[i][j],-q[i][j]);
      }

      // SI units
      if(opt.Contains("units")) {
	Double_t spaUnit, eneUnit;
	string spaSUnit, eneSUnit;
	
	PUnits::BestUnit bspaSUnit(TMath::TwoPi()*skindepth,"Length");
	bspaSUnit.GetBestUnits(spaUnit,spaSUnit);
	
	eneUnit = PUnits::MeV;
	eneSUnit = "MeV";

	Int_t NbinsX = hEneVsZ[i]->GetNbinsX();
	Double_t xMin = skindepth * hEneVsZ[i]->GetXaxis()->GetXmin() / spaUnit;
	Double_t xMax = skindepth * hEneVsZ[i]->GetXaxis()->GetXmax() / spaUnit;
	Int_t NbinsY = hEneVsZ[i]->GetNbinsY();
	Double_t ymin = PConst::ElectronMassE * hEneVsZ[i]->GetYaxis()->GetXmin() / eneUnit;
	Double_t ymax = PConst::ElectronMassE * hEneVsZ[i]->GetYaxis()->GetXmax() / eneUnit;
	hEneVsZ[i]->SetBins(NbinsX,xMin,xMax,NbinsY,ymin,ymax);

	hEneVsZ[i]->GetXaxis()->SetTitle(Form("#zeta [%s]",spaSUnit.c_str()));   
	hEneVsZ[i]->GetYaxis()->SetTitle(Form("Energy [%s]",eneSUnit.c_str()));   
	hEneVsZ[i]->GetZaxis()->SetTitle(Form("Charge [arb. unit]"));
      } else {
	hEneVsZ[i]->GetXaxis()->SetTitle(Form("k_{p} #zeta"));   
	hEneVsZ[i]->GetYaxis()->SetTitle(Form("#gamma"));   
	hEneVsZ[i]->GetZaxis()->SetTitle(Form("Charge [arb. unit]"));
      }

      
      hCurr[i] = hEneVsZ[i]->ProjectionX(Form("hCurr_%1i",i),1,Nx1Bin);
      if(hCurr[i]->GetMaximum()>maxCur) maxCur = hCurr[i]->GetMaximum();
      hSpec[i] = hEneVsZ[i]->ProjectionY(Form("hSpec_%1i",i),1,NeBin);
      if(hSpec[i]->GetMaximum()>maxEne) maxEne = hSpec[i]->GetMaximum();

      if(!hEneVsZjoint) {
	hEneVsZjoint = (TH2F*) hEneVsZ[i]->Clone("hEneVsZjoint");
      } else {
	hEneVsZjoint->Add(hEneVsZ[i]);
      }
      
    }
    
    // Colors and palettes        
    PPalette * beamPalette = (PPalette*) gROOT->FindObject("beam");
    if(!beamPalette) {
      beamPalette = new PPalette("beam");
    }
    beamPalette->SetPalette("elec0");

    PPalette * beam2Palette = (PPalette*) gROOT->FindObject("beam2");
    if(!beam2Palette) {
      beam2Palette = new PPalette("beam2");
    }
    beam2Palette->SetPalette("elec0");
    
    PPalette * beam3Palette = (PPalette*) gROOT->FindObject("beam3");
    if(!beam3Palette) {
      beam3Palette = new PPalette("beam3");
    }
    beam3Palette->SetPalette("elec0");
    
    TExec **exPal = new TExec*[Nspecies];
    Int_t Npal = 0;
    Int_t *colors = new Int_t[Nspecies];
    Int_t *fillstyle = new Int_t[Nspecies];
    for(Int_t i=0;i<Nspecies;i++) {
      if(!hEneVsZ[i]) continue;
      
      if(i==0) {
	exPal[i] = new TExec("exPlasma","plasmaPalette->cd();");
	colors[i] = kGray;
	fillstyle[i] = 0;
      } else if(i==1) {
	exPal[i] = new TExec("exBeam","beamPalette->cd();");
	//colors[i] = TColor::GetColor(69,108,155);
	colors[i] = kGray+1;
	fillstyle[i] = 0;
      } else if(i==2) {
	exPal[i] = new TExec(Form("exBeam%1i",i),Form("beam%1iPalette->cd();",i));  
	colors[i] = TColor::GetColor(83,109,161); // Bluish
	fillstyle[i] = 1001;
      } else {
	exPal[i] = new TExec(Form("exBeam%1i",i),Form("beam%1iPalette->cd();",i));  
	colors[i] = TColor::GetColor(231,99,51);
      	fillstyle[i] = 1001;
      }
      
      Npal++;
    }

    TPaletteAxis *palette = NULL;
        
    // Canvas setup
    Int_t sizex = 1024;
    Int_t sizey = 640;
    TCanvas *C = new TCanvas("C","Beams energy",sizex,sizey);
    C->SetFillStyle(4000);
    gPad->SetTickx(1);
    gPad->SetTicky(1);
    if(opt.Contains("logz"))
      gPad->SetLogz(1);

    TH2F *hFrame = (TH2F*) hEneVsZ[1]->Clone("hFrame");
    hFrame->Reset();
    PGlobals::SetH1LabelSize(hFrame);
    //hFrame->GetYaxis()->SetTickLength();

    hFrame->GetYaxis()->SetLabelFont(43);
    hFrame->GetYaxis()->SetLabelSize(32);
    // hFrame->GetYaxis()->SetLabelOffset(0.8);
    hFrame->GetYaxis()->SetTitleFont(43);
    hFrame->GetYaxis()->SetTitleSize(36);
    hFrame->GetYaxis()->SetTitleOffset(0.9);
    hFrame->GetYaxis()->SetTickLength(0.01);

    hFrame->GetXaxis()->SetLabelFont(43);
    hFrame->GetXaxis()->SetLabelSize(32);
    // hFrame->GetXaxis()->SetLabelOffset(0.8);
    hFrame->GetXaxis()->SetTitleFont(43);
    hFrame->GetXaxis()->SetTitleSize(36);
    hFrame->GetXaxis()->SetTitleOffset(0.95);
   
    hFrame->Draw("axis");

    gPad->Update();
    Float_t y1 = gPad->GetBottomMargin();
    Float_t y2 = 1 - gPad->GetTopMargin();
    Float_t x2 = 1 - gPad->GetRightMargin();
    Float_t gap = 0.02;
    Float_t pvsize = (y2-y1)/float(Npal);
    
    Int_t j = 0;
    for(Int_t i=0;i<Nspecies;i++) {
      if(!hEneVsZ[i])
	continue;
      
      if(opt.Contains("split")) {
	exPal[i]->Draw();
	hEneVsZ[i]->Draw("colz same");
	
	gPad->Update();
	palette = (TPaletteAxis*)hEneVsZ[i]->GetListOfFunctions()->FindObject("palette");
	if(palette) {
	  Float_t y1b,y2b;
	  
	  if(j==0)
	    y1b = y1 + j*pvsize;
	  else
	    y1b = y1 + j*pvsize + gap/2.0;
	  
	  if(j==Npal-1)  
	    y2b = y1 + (j+1)*pvsize;
	  else
	    y2b = y1 + (j+1)*pvsize - gap/2.0;
	  
	  palette->SetY1NDC(y1b);
	  palette->SetY2NDC(y2b);
	  
	  palette->SetX1NDC(x2 + 0.01);
	  palette->SetX2NDC(x2 + 0.03);
	  palette->SetBorderSize(2);
	  palette->SetLineColor(1);
	  
	  TPave *pFrame = new TPave((x2 + 0.01),y1b,(x2 + 0.03),y2b,1,"NDCL");
	  pFrame->SetFillStyle(0);
	  pFrame->SetLineColor(kBlack);
	  pFrame->SetShadowColor(0);
	  pFrame->Draw();
	  
	  j++;
	}
      } 

      // Plot 1D current
      Float_t yaxismin  =  gPad->GetUymin();
      Float_t yaxismax  =  gPad->GetUymin() + 0.3*(gPad->GetUymax() - gPad->GetUymin());

      Float_t curmin = 0.0;
      Float_t curmax = maxCur;
      // Round for better axis
      curmax = 0.1*TMath::Ceil(10*curmax);
      
      Float_t slope = (yaxismax - yaxismin)/(curmax - curmin);
      TH1D *hClone = (TH1D*) hCurr[i]->Clone(Form("hClone_%1i",i));
      for(Int_t k=0;k<hClone->GetNbinsX();k++) {
	Float_t content = hClone->GetBinContent(k+1);
	hClone->SetBinContent(k+1,(content - curmin) * slope + yaxismin);	
      }
      hClone->SetLineWidth(2);
      hClone->SetLineColor(colors[i]);

      hClone->Draw("same L HIST");


      // Plot 1D spectrum
      Float_t xaxismin  =  gPad->GetUxmin();
      Float_t xaxismax  =  gPad->GetUxmin() + 0.2*(gPad->GetUxmax() - gPad->GetUxmin());
      
      Float_t enemin = 0.0;
      Float_t enemax = maxEne;
      // Round for better axis
      enemax = 0.1*TMath::Ceil(10*enemax);

      Double_t *yarray   = new Double_t[NeBin];
      Double_t *xarray   = new Double_t[NeBin];
      for(Int_t k=0; k<NeBin; k++) {
	yarray[k] = hSpec[i]->GetBinCenter(k+1);
	xarray[k] = ((xaxismax-xaxismin)/enemax)*hSpec[i]->GetBinContent(k+1) + xaxismin;
      }
      gSpec[i] = new TGraph(NeBin,xarray,yarray);
      gSpec[i]->SetLineColor(colors[i]);
      gSpec[i]->SetLineWidth(2);
      gSpec[i]->SetFillStyle(fillstyle[i]);
      gSpec[i]->SetFillColor(PGlobals::elecFill);
      gSpec[i]->Draw("F");
      gSpec[i]->Draw("L");
      
    }

    Float_t xaxismin  =  gPad->GetUxmin();
    Double_t *yarray   = new Double_t[NeBin];
    Double_t *xarray   = new Double_t[NeBin];
    for(Int_t k=0;k<NeBin;k++) {
      xarray[k] = 0.0;
      yarray[k] = 0.0;
    }
    
    for(Int_t i=0;i<Nspecies;i++) {
      if(!gSpec[i])
	continue;
      // cout << Form("\n Specie %i ----------- ",i) << endl;
      Double_t *x = gSpec[i]->GetX();
      Double_t *y = gSpec[i]->GetY();

      for(Int_t k=0;k<NeBin;k++) {
	yarray[k] = y[k];
	xarray[k] += x[k]-xaxismin;
	//cout << Form("(%.2f, %2f)",xarray[k], yarray[k]) << endl;
      }
    }
    for(Int_t k=0;k<NeBin;k++) xarray[k] += xaxismin;
    
    TGraph *gSpecjoint = new TGraph(NeBin,xarray,yarray);
    //TGraph *gSpecjoint = (TGraph*) gSpec[2]->Clone("kk");
    gSpecjoint->SetLineColor(kGray+2);
    gSpecjoint->SetLineWidth(1);
    gSpecjoint->SetLineStyle(2);
    // gSpecjoint->SetFillStyle(1001);
    // gSpecjoint->SetFillColor(PGlobals::elecFill);
    // gSpecjoint->Draw("F");
    gSpecjoint->Draw("L");

    
    for(Int_t i=0;i<Nspecies;i++) {
      if(!gSpec[i])
	continue;
      gSpec[i]->Draw("L");
    }
      
    if(!opt.Contains("split")) {
      exPal[1]->Draw();

      hEneVsZjoint->GetZaxis()->SetLabelFont(43);
      hEneVsZjoint->GetZaxis()->SetLabelSize(32);
      // hEneVsZjoint->GetZaxis()->SetLabelOffset(0.8);
      hEneVsZjoint->GetZaxis()->SetTitleFont(43);
      hEneVsZjoint->GetZaxis()->SetTitleSize(36);
      hEneVsZjoint->GetZaxis()->SetTitleOffset(0.9);
      hEneVsZjoint->GetZaxis()->SetTickLength(0.01);
      hEneVsZjoint->GetZaxis()->CenterTitle();
      hEneVsZjoint->Draw("colz same");

      gPad->Update();      
      palette = (TPaletteAxis*)hEneVsZjoint->GetListOfFunctions()->FindObject("palette");
      if(palette) {
	palette->SetY2NDC(y2 - 0.0);
	palette->SetY1NDC(y1 + 0.0);
	palette->SetX1NDC(x2 + 0.01);
	palette->SetX2NDC(x2 + 0.03);
	palette->SetBorderSize(2);
	palette->SetLineColor(1);
	
	TPave *pFrame = new TPave((x2 + 0.01),y1 + 0.0,(x2 + 0.03),y2 - 0.0,2,"NDC");
	pFrame->SetFillStyle(0);
	pFrame->SetLineColor(kBlack);
	pFrame->SetLineWidth(1);
	pFrame->SetShadowColor(0);
	pFrame->Draw();
      }
    }
    
    TLine *linezero = new TLine(gPad->GetUxmin(),0.0,gPad->GetUxmax(),0.0);
    linezero->SetLineColor(kGray+2);
    linezero->SetLineStyle(2);
    linezero->Draw();
    
    
    TBox *lFrame = new TBox(gPad->GetUxmin(), gPad->GetUymin(),
			    gPad->GetUxmax(), gPad->GetUymax());
    lFrame->SetFillStyle(0);
    lFrame->SetLineColor(kBlack);
    lFrame->SetLineWidth(2);
    lFrame->Draw();
    
    gPad->RedrawAxis("g");

    
    // TH2F *hEneVsZjoint = (TH2F*) hEneVsZ[1]->Clone("hEneVsZjoint");
    // hEneVsZjoint->Add(hEneVsZ[2]);
    // // beamPalette->cd();
    // hEneVsZjoint->Draw("colz");
  
    // Print to a file
    // Output file
    TString fOutName = Form("./%s/Plots/BeamsEnergy/BeamsEnergy-%s_%i",pData->GetPath().c_str(),pData->GetName(),time);
  
    PGlobals::imgconv(C,fOutName,opt);
    // ---------------------------------------------------------
  

  }
}
