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
#include <TGaxis.h>
#include <TLatex.h>

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
    printf("      <-cmax(value)> <-imax(value)> <-gmax(value)> <-dmax(value)> \n");
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

  Float_t cMax = 0.0;
  Float_t iMax = -99999.;
  Float_t gMax = -99999.;
  Float_t dMax = -99999.;

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
    } else if(arg.Contains("--spec")){
      opt += "spec"; 
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
    } else if(arg.Contains("-cmax")) {
      char ss[5];
      sscanf(arg,"%5s%f",ss,&cMax);
    } else if(arg.Contains("-curmax")) {
      char ss[7];
      sscanf(arg,"%7s%f",ss,&iMax);
    } else if(arg.Contains("-gmax")) {
      char ss[5];
      sscanf(arg,"%5s%f",ss,&gMax);
    } else if(arg.Contains("-dmax")) {
      char ss[5];
      sscanf(arg,"%5s%f",ss,&dMax);
    } else {
      sim = arg;
    }
  }
  
  PGlobals::Initialize();

  if(opt.Contains("grid")) {
    gStyle->SetPadGridX(1);
    gStyle->SetPadGridY(1);
  }

  gStyle->SetJoinLinePS(2);
  gStyle->SetPadRightMargin(0.20);
  //  gStyle->SetFrameLineWidth(3);

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

    cout << Form("\n Reading data arrays") << endl;
    
    Int_t Nspecies = pData->NSpecies();
    Float_t **ene = new Float_t*[Nspecies];
    Float_t **x1 = new Float_t*[Nspecies];
    Float_t **q = new Float_t*[Nspecies];
    Float_t **p1 = new Float_t*[Nspecies];
    Float_t **p2 = new Float_t*[Nspecies];
    Float_t **p3 = new Float_t*[Nspecies];
    UInt_t *NP = new UInt_t[Nspecies];
    Float_t eneMin = 99999;
    Float_t eneMax = -99999;
    Float_t divMin = 99999;
    Float_t divMax = -99999;
    
    for(Int_t i=0;i<Nspecies;i++) {
      ene[i] = NULL;
      x1[i] = NULL;
      NP[i] = 0;
    
      if(!pData->GetRawFileName(i)) continue;

      cout << Form(" - Species %i ",i) << endl;


      NP[i] = pData->GetRawSingleArray(pData->GetRawFileName(i)->c_str(),&x1[i],"x1");
      NP[i] = pData->GetRawSingleArray(pData->GetRawFileName(i)->c_str(),&ene[i],"ene");
      NP[i] = pData->GetRawSingleArray(pData->GetRawFileName(i)->c_str(),&q[i],"q");

      if(opt.Contains("spec")) {
	pData->GetRawSingleArray(pData->GetRawFileName(i)->c_str(),&p1[i],"p1");
	pData->GetRawSingleArray(pData->GetRawFileName(i)->c_str(),&p2[i],"p2");
	pData->GetRawSingleArray(pData->GetRawFileName(i)->c_str(),&p3[i],"p3");
      }
      // for(UInt_t j=0;j<1000;j++)
      // 	cout << Form(" %6f  %6f  %6f",x1[i][j],ene[i][j],q[i][j]) << endl;
    
      for(UInt_t j=0;j<NP[i];j++) {
      	if(ene[i][j]<eneMin) eneMin = ene[i][j];
      	if(ene[i][j]>eneMax) eneMax = ene[i][j];

	//	Float_t div = TMath::Sqrt(p2[i][j]*p2[i][j]+p3[i][j]*p3[i][j])/p1[i][j];
	Float_t div = p2[i][j]/p1[i][j];
	if(div<divMin) divMin = div;
      	if(div>divMax) divMax = div;

      }
    }    

    // Histogram limits
    Double_t X1MIN = pData->GetXMin(0) - shiftz;
    Double_t X1MAX = pData->GetXMax(0) - shiftz;
    Double_t dx1 = pData->GetDX(0);    
    Float_t x1Min = pData->GetX1Min() - shiftz;
    Float_t x1Max = pData->GetX1Max() - shiftz;
    Int_t Nx1Bin  = pData->GetX1N();
    x1Min -= 0.1*(x1Max - x1Min);

    // Adjust binning
    x1Min = floor((x1Min-X1MIN)/dx1) * dx1 + X1MIN;  
    x1Max = floor((x1Max-X1MIN)/dx1) * dx1 + X1MIN;
    Nx1Bin  = ceil ((x1Max - x1Min)/(dx1));
    
    Float_t eMin = eneMin - (eneMax-eneMin) * 0.1;
    // if(eMin<0) eMin = 0.0;
    eMin = 0.0;
     
    Float_t eMax = eneMax + (eneMax-eneMin) * 0.1;
    if(gMax >= 0.0) eMax = gMax;
    
    Int_t NeBin  = Nx1Bin;
    Double_t dene = (eMax-eMin)/NeBin;

    Int_t NdivBin = 100;

    Float_t x1Minu(x1Min), x1Maxu(x1Max), eMinu(eMin), eMaxu(eMax);
    Double_t denUnit, propUnit, spaUnit, eneUnit, charUnit, curUnit;
    string denSUnit, propSUnit, spaSUnit, eneSUnit, charSUnit, curSUnit;
    if(opt.Contains("units")) {
      PUnits::BestUnit bdenSUnit(n0,"PartDensity");
      bdenSUnit.GetBestUnits(denUnit,denSUnit);
      
      propUnit = PUnits::mm;
      propSUnit = "mm";
      
      PUnits::BestUnit bspaSUnit(TMath::TwoPi()*skindepth,"Length");
      bspaSUnit.GetBestUnits(spaUnit,spaSUnit);
      
      eneUnit = PUnits::MeV;
      eneSUnit = "MeV";
      
      x1Minu = x1Min * skindepth / spaUnit;
      x1Maxu = x1Max * skindepth / spaUnit;
      eMinu = eMin * PConst::ElectronMassE / eneUnit;
      eMaxu = eMax * PConst::ElectronMassE / eneUnit;
    }

    // Charge normalisation
    Double_t dV =  pData->GetDX(0)*pData->GetDX(1)*pData->GetDX(2) * skindepth * skindepth * skindepth;
    Double_t Q0 = fabs(n0 * dV);
    charUnit = PUnits::picocoulomb;
    charSUnit = "pC";

    // Current units
    curUnit = PUnits::kA;
    curSUnit = "kA";
    
    cout << Form("\n Filling histograms ") << endl;
    
    char hName[24];
    TH2F **hEneVsZ  = new TH2F*[Nspecies];
    TH2F *hEneVsZjoint = NULL;
    TH2F **hEneVsDiv  = new TH2F*[Nspecies];
    TH2F *hEneVsDivjoint = NULL;
    TH1D **hCurr = new TH1D*[Nspecies];
    TH1D **hClone = new TH1D*[Nspecies];
    TH1D **hSpec = new TH1D*[Nspecies];
    TGraph **gSpec = new TGraph*[Nspecies];   
    Float_t maxCur = -9999;
    Float_t maxEne = -9999;
    for(Int_t i=0;i<Nspecies;i++) {
    
      hEneVsZ[i] = NULL;
      hEneVsDiv[i] = NULL;
      hSpec[i] = NULL;
      hCurr[i] = NULL;
      hClone[i] = NULL;
      gSpec[i] = NULL;
      
      if(ene[i]==NULL) continue;

      cout << Form(" - Species %i",i) << endl;

      // 2D histogram: #gamma vs \zeta  
      sprintf(hName,"EneVsZ_%i",i);
      hEneVsZ[i] = (TH2F*) gROOT->FindObject(hName);
      if(hEneVsZ[i]) delete hEneVsZ[i];
      hEneVsZ[i] = new TH2F(hName,"",Nx1Bin,x1Min,x1Max,NeBin,eMin,eMax);
      
      // 2D histogram: #gamma vs \zeta
      divMax = 19.99E-3;
      divMin = -divMax;
      sprintf(hName,"EneVsDiv_%i",i);
      hEneVsDiv[i] = (TH2F*) gROOT->FindObject(hName);
      if(hEneVsDiv[i]) delete hEneVsDiv[i];
      hEneVsDiv[i] = new TH2F(hName,"",NeBin,eMin,eMax,NdivBin,divMin,divMax);
      
      for(UInt_t j=0;j<NP[i];j++) {
	hEneVsZ[i]->Fill(x1[i][j]-shiftz,ene[i][j],-q[i][j]);
	// Float_t div = TMath::Sqrt(p2[i][j]*p2[i][j]+p3[i][j]*p3[i][j])/p1[i][j];
	Float_t div = p2[i][j]/p1[i][j];
	hEneVsDiv[i]->Fill(ene[i][j]*0.511,div,-q[i][j]);
      }

      // SI units
      if(opt.Contains("units")) {
	hEneVsZ[i]->SetBins(Nx1Bin,x1Minu,x1Maxu,NeBin,eMinu,eMaxu);
	hEneVsZ[i]->GetXaxis()->SetTitle(Form("#zeta [%s]",spaSUnit.c_str()));   
	hEneVsZ[i]->GetYaxis()->SetTitle(Form("Energy [%s]",eneSUnit.c_str()));   
	hEneVsZ[i]->GetZaxis()->SetTitle(Form("Charge [%s]",charSUnit.c_str()));
	//Double_t dx = dx1 * skindepth / spaUnit;
	//Double_t de = dene * PConst::ElectronMassE / eneUnit;
	//hEneVsZ[i]->Scale(((Q0/charUnit)/(dx*dene)) * );
	hEneVsZ[i]->Scale(Q0*(-PConst::ElectronCharge)/charUnit);
	//hEneVsZ[i]->ResetStats();
      } else {
	hEneVsZ[i]->GetXaxis()->SetTitle(Form("k_{p} #zeta"));   
	hEneVsZ[i]->GetYaxis()->SetTitle(Form("#gamma"));   
	hEneVsZ[i]->GetZaxis()->SetTitle(Form("Charge [e^{-}]"));
      }
      
      if(cMax>0.0) {
	hEneVsZ[i]->GetZaxis()->SetRangeUser(0.0,cMax);
      }
      
      hCurr[i] = hEneVsZ[i]->ProjectionX(Form("hCurr_%1i",i),1,Nx1Bin);
      if(opt.Contains("units")) {
	Double_t dt = dx1 * skindepth / PConst::c_light;
	hCurr[i]->Scale((PUnits::picocoulomb/dt)/curUnit);
	//hCurr[i]->ResetStats();
      }
      if(hCurr[i]->GetMaximum()>maxCur && i<3) maxCur = hCurr[i]->GetMaximum();

      hSpec[i] = hEneVsZ[i]->ProjectionY(Form("hSpec_%1i",i),1,NeBin);
      if(hSpec[i]->GetMaximum()>maxEne) maxEne = hSpec[i]->GetMaximum();
      
      if(!hEneVsZjoint) {
	hEneVsZjoint = (TH2F*) hEneVsZ[i]->Clone("hEneVsZjoint");
      } else {
	hEneVsZjoint->Add(hEneVsZ[i]);
      }

      if(!hEneVsDivjoint) {
	hEneVsDivjoint = (TH2F*) hEneVsDiv[i]->Clone("hEneVsDivjoint");
      } else {
	hEneVsDivjoint->Add(hEneVsDiv[i]);
      }

    }

    if(cMax>0.0) {
      hEneVsZjoint->GetZaxis()->SetRangeUser(0.0,cMax);
    }

    // Round for better axis
    maxCur = 0.1*TMath::Ceil(10*maxCur);

    // SI units
    if(opt.Contains("units")) {
      x1Min = x1Minu;
      x1Max = x1Maxu;
      eMin = eMinu;
      eMax = eMaxu;

      Double_t dt = dx1 * spaUnit / PConst::c_light;
      
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
    Int_t *fillcolor = new Int_t[Nspecies];
    for(Int_t i=0;i<Nspecies;i++) {
      if(!hEneVsZ[i]) continue;
      
      if(i==0) {
	exPal[i] = new TExec("exPlasma","plasmaPalette->cd();");
	colors[i] = kGray;
	fillstyle[i] = 0;
	fillcolor[i] = TColor::GetColor(250,250,250);
      } else if(i==1) {
	exPal[i] = new TExec("exBeam","beamPalette->cd();");
	//colors[i] = TColor::GetColor(69,108,155);
	colors[i] = kGray+1;
	fillstyle[i] = 1001;
	fillcolor[i] = TColor::GetColor(230,230,230);
      } else if(i==2) {
	exPal[i] = new TExec(Form("exBeam%1i",i),Form("beam%1iPalette->cd();",i));  
	colors[i] = TColor::GetColor(83,109,161); // Bluish
	fillstyle[i] = 1001;
	fillcolor[i] = PGlobals::elecFill;
      } else {
	exPal[i] = new TExec(Form("exBeam%1i",i),Form("beam%1iPalette->cd();",i));  
	colors[i] = TColor::GetColor(231,99,51);
      	fillstyle[i] = 1001;
	fillcolor[i] = PGlobals::elecFill;
      }
      
      Npal++;
    }

    TPaletteAxis *palette = NULL;
        
    // Canvas setup
    Int_t sizex = 1024;
    Int_t sizey = 640;
    TCanvas *C = new TCanvas("C","Beams energy",sizex,sizey);
    // C->SetFillStyle(4000);
    gPad->SetTickx(1);
    gPad->SetTicky(1);
    if(opt.Contains("logz"))
      gPad->SetLogz(1);

    // Plot frame
    // Float_t xmin = x1Min - 0.4*(x1Max - x1Min);
    Float_t xmin = x1Min;
    Float_t xmax = x1Max;
    Float_t ymin = eMin - 0.3*(eMax - eMin);
    Float_t ymax = eMax;

    Float_t curmin = 0.0;
    Float_t curmax = maxCur;

    TH2F *hFrame = (TH2F*) gROOT->FindObject("hFrame");
    if(hFrame) delete hFrame;
    hFrame = new TH2F("hFrame","",10,xmin,xmax,10,ymin,ymax);
    if(opt.Contains("units")) {
      hFrame->GetXaxis()->SetTitle(Form("#zeta [%s]",spaSUnit.c_str()));   
      hFrame->GetYaxis()->SetTitle(Form("Energy [%s]",eneSUnit.c_str()));   
    } else {
      hFrame->GetXaxis()->SetTitle(Form("k_{p} #zeta"));   
      hFrame->GetYaxis()->SetTitle(Form("#gamma"));   
    }
    hFrame->GetXaxis()->CenterTitle();
    hFrame->GetYaxis()->CenterTitle();
    hFrame->Draw("axis");

    gPad->Update();
    Float_t y1 = gPad->GetBottomMargin();
    Float_t y2 = 1 - gPad->GetTopMargin();
    Float_t x2 = 1 - gPad->GetRightMargin();

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
	  Float_t gap = 0.02;
	  Float_t pvsize = (y2-y1)/float(Npal);
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
      // Float_t yaxismax  =  gPad->GetUymin() + 0.3*(gPad->GetUymax() - gPad->GetUymin());
      Float_t yaxismax  =  eMin;
      
      if(iMax >= 0.0) curmax = iMax;
      
      Float_t slope = (yaxismax - yaxismin)/(curmax - curmin);
      hClone[i] = (TH1D*) hCurr[i]->Clone(Form("hClone_%1i",i));
      for(Int_t k=0;k<hClone[i]->GetNbinsX();k++) {
	Float_t content = hClone[i]->GetBinContent(k+1);
	hClone[i]->SetBinContent(k+1,(content - curmin) * slope + yaxismin);	
      }
      hClone[i]->SetLineWidth(2);
      hClone[i]->SetLineColor(colors[i]);
      hClone[i]->SetFillStyle(fillstyle[i]);
      hClone[i]->SetFillColor(fillcolor[i]);

      hClone[i]->Draw("same LF2 HIST");

      // Plot 1D spectrum
      Float_t xaxismin  =  gPad->GetUxmin();
      // Float_t xaxismax  =  gPad->GetUxmin() + 0.2*(gPad->GetUxmax() - gPad->GetUxmin());
      //  Float_t xaxismax  = x1Min;
      Float_t xaxismax  = xmin + (xmax-xmin) * 0.2;
      Float_t enemin = 0.0;
      Float_t enemax = maxEne;
      // Round for better axis
      enemax = 0.1*TMath::Ceil(10*enemax);

      Double_t *yarray   = new Double_t[NeBin];
      Double_t *xarray   = new Double_t[NeBin];
      for(Int_t k=0; k<NeBin; k++) {
	yarray[k] = hSpec[i]->GetBinCenter(k+1);
	xarray[k] = ((xaxismax-xaxismin)/enemax)*hSpec[i]->GetBinContent(k+1) + xaxismin;
	if(k==0) xarray[k] = xaxismin;
      }
      gSpec[i] = new TGraph(NeBin,xarray,yarray);
      gSpec[i]->SetLineColor(colors[i]);
      gSpec[i]->SetLineWidth(2);
      gSpec[i]->SetFillStyle(fillstyle[i]);
      gSpec[i]->SetFillColor(fillcolor[i]);
      gSpec[i]->Draw("F");
      gSpec[i]->Draw("L");
      
    }

    // Current axis
    Float_t x1pos = xmin + (xmax-xmin) * 0.2;
    TGaxis *axis = new TGaxis(x1pos,ymin,
			      x1pos,eMin,
			      0.0,curmax,503,"S");
    axis->SetLineWidth(1);
    axis->SetLineColor(kGray+3);//PGlobals::elecLine);
    axis->SetLabelColor(kGray+3);//PGlobals::elecLine);
    axis->SetLabelFont(PGlobals::fontType);
    axis->SetLabelSize(PGlobals::labelSize-10);
    axis->SetLabelOffset(0.005);
    axis->SetTitleColor(kGray+3);//PGlobals::elecLine);
    axis->SetTitleFont(PGlobals::fontType);
    axis->SetTitleSize(PGlobals::titleSize-10);
    axis->SetTitleOffset(0.8);
    axis->SetTickSize(0.03);
    axis->ChangeLabel(1,-1,0.);
    if(opt.Contains("units"))
      axis->SetTitle(Form("I [%s]",curSUnit.c_str()));
    axis->CenterTitle();
    //axis->ChangeLabel(1,-1,-1,-1,-1,-1,"");
    // axis->SetMaxDigits(2);
    axis->Draw();
    
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

      // if (hClone[i]) {
      // 	hClone[i]->SetFillStyle(0);
      // 	hClone[i]->SetFillColor(0);
      // 	hClone[i]->Draw("same L HIST");
      // }
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

    char ctext[128];
    if(opt.Contains("units") && n0) 
      sprintf(ctext,"z = %5.2f %s", Time * skindepth / propUnit, propSUnit.c_str());
    else
      sprintf(ctext,"k_{p}z = %5.1f",Time);
    TLatex *textTime = new TLatex(xmax - (xmax-xmin)/20.,ymax - (ymax-ymin)/15.,ctext);
    textTime->SetTextAlign(32);
    textTime->SetTextFont(PGlobals::fontType);
    textTime->SetTextSize(PGlobals::titleSize-10);
    textTime->Draw();
    
    if(opt.Contains("units") && n0) 
      sprintf(ctext,"n_{0} = %5.2f x %s", n0/denUnit,denSUnit.c_str());
    TLatex *textDen = new TLatex(xmin + (xmax-xmin)/20.,ymax - (ymax-ymin)/15.,ctext);
    textDen->SetTextAlign(12);
    textDen->SetTextFont(PGlobals::fontType);
    //  textDen->SetTextColor(lineColor);
    textDen->SetTextColor(kGray+3);
    textDen->SetTextSize(PGlobals::titleSize-14);
    textDen->Draw();
    
    TLine *linezero = new TLine(gPad->GetUxmin(),0.0,gPad->GetUxmax(),0.0);
    linezero->SetLineColor(kGray+3);
    linezero->SetLineStyle(3);
    //linezero->SetLineWidth(2);
    linezero->Draw();
    
    TLine *lineY = new TLine(x1Min,gPad->GetUymin(),x1Min,gPad->GetUymax());
    lineY->SetLineColor(kGray+3);
    lineY->SetLineStyle(3);
    //lineY->SetLineWidth(2);
    lineY->Draw();

    
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


    // Canvas setup
    sizex = 1024;
    sizey = 320;
    TCanvas *C2 = new TCanvas("C","Spectrum",sizex,sizey);
    // C->SetFillStyle(4000);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.20);
    gPad->SetTickx(1);
    gPad->SetTicky(1);
    if(opt.Contains("logz"))
      gPad->SetLogz(1);
      
    //  gPad->SetLogx(1);
    // gPad->SetFrameLineWidth(3);
    // hEneVsDivjoint->GetXaxis()->SetMoreLogLabels(1);
    // hEneVsDivjoint->GetXaxis()->SetNdivisions(510);
    hEneVsDivjoint->GetXaxis()->SetAxisColor(0);
    hEneVsDivjoint->GetYaxis()->SetAxisColor(0);
    hEneVsDivjoint->GetXaxis()->SetTickLength(0.03);
    
    PPalette * specPalette = (PPalette*) gROOT->FindObject("spec");
    if(!specPalette) {
      specPalette = new PPalette("spec");
    }
    //    specPalette->SetPalette(60);
    specPalette->SetPalette("spectrum");
    hEneVsDivjoint->SetMinimum(-0.00001);

    // cout << Form("dmax = %f",hEneVsDivjoint->GetMaximum()) << endl;
    if(dMax >= 0.0)
      hEneVsDivjoint->SetMaximum(dMax);

    if(gMax >= 0.0)
      // hEneVsDivjoint->GetXaxis()->SetRangeUser(39.999,gMax);
      hEneVsDivjoint->GetXaxis()->SetRangeUser(0.0,gMax);
    
    hEneVsDivjoint->Draw("col2");
    gPad->SetFrameLineWidth(5);
    // gPad->SetFrameLineColor(0);
    gPad->Update();
    gPad->RedrawAxis();
      
    
    // Print to a file
    // Output file
    fOutName = Form("./%s/Plots/BeamsEnergy/Spectrum-%s_%i",pData->GetPath().c_str(),pData->GetName(),time);
  
    PGlobals::imgconv(C2,fOutName,opt);
    // ---------------------------------------------------------
    
    // delete objects:
    if(hEneVsZjoint)
      delete hEneVsZjoint;
    if(hEneVsDivjoint)
      delete hEneVsDivjoint;
    for(Int_t i=0;i<Nspecies;i++) {
      if(hEneVsZ[i]) delete hEneVsZ[i];
      if(hEneVsDiv[i]) delete hEneVsDiv[i];
      if(hCurr[i]) delete hCurr[i];
      if(hClone[i]) delete hClone[i];
      if(hSpec[i]) delete hSpec[i];
      if(gSpec[i]) delete gSpec[i];
       
    }
    
    
    
    PGlobals::DestroyCanvases();
  }
}
