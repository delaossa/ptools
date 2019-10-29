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
#include "PDataHiP.hh"
#include "PGlobals.hh"
#include "PPalette.hh"
#include "H5Cpp.h"

using namespace std;
using namespace H5;

void redrawBorder()
{
   // this little macro redraws the axis tick marks and the pad border lines.
   gPad->Update();
   gPad->RedrawAxis();
   TLine l;
   l.DrawLine(gPad->GetUxmin(), gPad->GetUymax(), gPad->GetUxmax(), gPad->GetUymax());
   l.DrawLine(gPad->GetUxmax(), gPad->GetUymin(), gPad->GetUxmax(), gPad->GetUymax());
}

int main(int argc,char *argv[]) {
  if(argc<=2) {
    printf("\n Usage: %s <simulation name> <-t(time)>\n",argv[0]);
    printf("      <-i(initial time)> <-f(final time)> <-s(time step)>\n");
    printf("      <-cmax(value)> <-imax(value)> <-gmax(value)> <-dmin(value)> <-dmax(value)> <-smax(value)> \n");
    printf("      <--center> <--comov> <--units> <--logz> <--outl>\n");
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
  Float_t gMin = -99999.;
  Float_t gMax = -99999.;
  Float_t dMin = -99999.;
  Float_t dMax = -99999.;
  Float_t divMin = -99999.;
  Float_t divMax = -99999.;
  Float_t sMax = -99999.;

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
    } else if(arg.Contains("--nopal")){
      opt += "nopal"; 
    } else if(arg.Contains("--wpal")){
      opt += "wpal"; 
    } else if(arg.Contains("--hzdr")){
      opt += "hzdr"; 
    } else if(arg.Contains("--ffwd")){
      opt += "ffwd"; 
    } else if(arg.Contains("--outl")){
      opt += "outl"; 
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
    } else if(arg.Contains("-gmin")) {
      char ss[5];
      sscanf(arg,"%5s%f",ss,&gMin);
    } else if(arg.Contains("-gmax")) {
      char ss[5];
      sscanf(arg,"%5s%f",ss,&gMax);
    } else if(arg.Contains("-dmin")) {
      char ss[5];
      sscanf(arg,"%5s%f",ss,&dMin);
    } else if(arg.Contains("-dmax")) {
      char ss[5];
      sscanf(arg,"%5s%f",ss,&dMax);
    } else if(arg.Contains("-divmin")) {
      char ss[7];
      sscanf(arg,"%7s%f",ss,&divMin);
    } else if(arg.Contains("-divmax")) {
      char ss[7];
      sscanf(arg,"%7s%f",ss,&divMax);
    } else if(arg.Contains("-smax")) {
      char ss[5];
      sscanf(arg,"%5s%f",ss,&sMax);
    } else if(arg.Contains("-i")) {
      char ss[2];
      sscanf(arg,"%2s%i",ss,&iStart);
    } else if(arg.Contains("-f")) {
      char ss[2];
      sscanf(arg,"%2s%i",ss,&iEnd);
    } else if(arg.Contains("-s")) {
      char ss[2];
      sscanf(arg,"%2s%i",ss,&iStep);
    } else {
      sim = arg;
    }
  }
  
  PGlobals::Initialize();

  if(opt.Contains("grid")) {
    gStyle->SetPadGridX(1);
    gStyle->SetPadGridY(1);
  }

  gStyle->SetJoinLinePS(1);
  gStyle->SetPadRightMargin(0.20);
  //  gStyle->SetFrameLineWidth(3);
  gStyle->SetLineWidth(2);

  // Load PData
  PData *pData = PData::Get(sim.Data());
  if(pData->isHiPACE()) {
    delete pData; pData = NULL;
    pData = PDataHiP::Get(sim.Data());
    if(!opt.Contains("comov"))
      opt += "comov";
  }

  if(iStart<0) iStart = time;
  if(iEnd<=iStart) iEnd = iStart;
  
  // Some plasma constants
  Float_t n0 = pData->GetPlasmaDensity();
  Float_t kp = pData->GetPlasmaK();
  Float_t skindepth = pData->GetPlasmaSkinDepth();
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
    
    Int_t Nspecies = pData->NRawSpecies();
    Float_t **ene = new Float_t*[Nspecies];
    Float_t **x1 = new Float_t*[Nspecies];
    Float_t **q = new Float_t*[Nspecies];
    Float_t **p1 = new Float_t*[Nspecies];
    Float_t **p2 = new Float_t*[Nspecies];
    Float_t **p3 = new Float_t*[Nspecies];
    UInt_t *NP = new UInt_t[Nspecies];
    Float_t eneMin = 99999;
    Float_t eneMax = -99999;
    Float_t divmin = 99999;
    Float_t divmax = -99999;

    cout << Form("\n Number of RAW species: %i", Nspecies) << endl;

    for(Int_t i=0;i<Nspecies;i++) {
      ene[i] = NULL;
      x1[i] = NULL;
      NP[i] = 0;
    
      if(!pData->GetRawFileName(i)) continue;

      cout << Form(" - Species %i ",i) << endl;


      NP[i] = pData->GetRawSingleArray(pData->GetRawFileName(i)->c_str(),&x1[i],"x1");
      // NP[i] = pData->GetRawSingleArray(pData->GetRawFileName(i)->c_str(),&ene[i],"ene");
      NP[i] = pData->GetRawSingleArray(pData->GetRawFileName(i)->c_str(),&q[i],"q");

      if(opt.Contains("spec")) {
	pData->GetRawSingleArray(pData->GetRawFileName(i)->c_str(),&p1[i],"p1");
	pData->GetRawSingleArray(pData->GetRawFileName(i)->c_str(),&p2[i],"p2");
	pData->GetRawSingleArray(pData->GetRawFileName(i)->c_str(),&p3[i],"p3");
      }

      ene[i] = new Float_t[NP[i]];
      
      // for(UInt_t j=0;j<1000;j++)
      // 	cout << Form(" %6f  %6f  %6f",x1[i][j],ene[i][j],q[i][j]) << endl;
    
      for(UInt_t j=0;j<NP[i];j++) {
	ene[i][j] = TMath::Sqrt(p1[i][j]*p1[i][j] + p2[i][j]*p2[i][j] + p3[i][j]*p3[i][j]) - 1.0;
	
      	if(ene[i][j]<eneMin) eneMin = ene[i][j];
      	if(ene[i][j]>eneMax) eneMax = ene[i][j];

	//	Float_t div = TMath::Sqrt(p2[i][j]*p2[i][j]+p3[i][j]*p3[i][j])/p1[i][j];
	Float_t div = 1e3 * p2[i][j]/p1[i][j];
	if(div<divmin) divmin = div;
      	if(div>divmax) divmax = div;

      }
    }    

    
    // Histogram limits
    Double_t X1MIN = pData->GetXMin(0) - shiftz;
    Double_t X1MAX = pData->GetXMax(0) - shiftz;
    Double_t dx1 = pData->GetDX(0);    
    Float_t x1Min = pData->GetX1Min() - shiftz;
    Float_t x1Max = pData->GetX1Max() - shiftz;
    Int_t Nx1Bin  = pData->GetX1N();
    x1Min -= 0.2*(x1Max - x1Min);
    
    // Adjust binning
    x1Min = floor((x1Min-X1MIN)/dx1) * dx1 + X1MIN;  
    x1Max = floor((x1Max-X1MIN)/dx1) * dx1 + X1MIN;
    Nx1Bin  = ceil ((x1Max - x1Min)/(dx1));
    
    Float_t eMin = eneMin - (eneMax-eneMin) * 0.1;
    if(gMin >= 0.0) eMin = gMin;
     
    Float_t eMax = eneMax + (eneMax-eneMin) * 0.1;
    if(gMax >= 0.0) eMax = gMax;
    
    Int_t NeBin  = Nx1Bin;
    Double_t de = (eMax-eMin)/NeBin;
    Double_t dx = (x1Max-x1Min)/Nx1Bin;
    //Int_t NdivBin = NeBin;
    Int_t NdivBin = 100;
    
    Float_t dxu(dx), deu(de);
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
      dxu = dx * skindepth / spaUnit;
      
      eMinu = eMin * PConst::ElectronMassE / eneUnit;
      eMaxu = eMax * PConst::ElectronMassE / eneUnit;
      deu = de *  PConst::ElectronMassE / eneUnit;
      
    }

    // Charge normalisation
    Double_t dV =  pData->GetDX(0) * pData->GetDX(1) * pData->GetDX(2) * skindepth * skindepth * skindepth;
    Double_t Q0 = fabs(n0 * dV * PConst::ElectronCharge);
    charUnit = PUnits::picocoulomb;
    charSUnit = "pC";
    Double_t Q0u = Q0 / charUnit;

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

    Float_t diveMin = divmin - (divmax-divmin) * 0.1;
    if(divMin != -99999.) diveMin = divMin;
     
    Float_t diveMax = divmax + (divmax-divmin) * 0.1;
    if(divMax != -99999.) diveMax = divMax;

    Double_t dvu = (diveMax-diveMin)/NdivBin;

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
      sprintf(hName,"EneVsDiv_%i",i);
      hEneVsDiv[i] = (TH2F*) gROOT->FindObject(hName);
      if(hEneVsDiv[i]) delete hEneVsDiv[i];
      hEneVsDiv[i] = new TH2F(hName,"",NeBin,eMin,eMax,NdivBin,diveMin,diveMax);
      
      for(UInt_t j=0;j<NP[i];j++) {
	Float_t zeta = x1[i][j]-shiftz;
	Float_t ener = ene[i][j];
	Float_t div  = 1000*(p2[i][j]/p1[i][j]);

	// if (ener<eMin || ener>eMax)
	//   continue;
	if (div<diveMin || div>diveMax)
	  continue;
	
	hEneVsZ[i]->Fill(zeta,ener,TMath::Abs(q[i][j]));
	hEneVsDiv[i]->Fill(ener,div,TMath::Abs(q[i][j]));
      }

      // SI units
      if(opt.Contains("units")) {
	hEneVsZ[i]->SetBins(Nx1Bin,x1Minu,x1Maxu,NeBin,eMinu,eMaxu);
	hEneVsZ[i]->GetXaxis()->SetTitle(Form("#zeta [%s]",spaSUnit.c_str()));   
	hEneVsZ[i]->GetYaxis()->SetTitle(Form("Energy [%s]",eneSUnit.c_str()));   
	hEneVsZ[i]->GetZaxis()->SetTitle(Form("Charge [%s]",charSUnit.c_str()));
	hEneVsZ[i]->Scale(Q0u/(deu*dxu));
	hEneVsZ[i]->ResetStats();

	hEneVsDiv[i]->SetBins(NeBin,eMinu,eMaxu,NdivBin,diveMin,diveMax);
	hEneVsDiv[i]->Scale(Q0u/(deu*dvu));
	hEneVsDiv[i]->ResetStats();
	
      } else {
	hEneVsZ[i]->GetXaxis()->SetTitle(Form("k_{p} #zeta"));   
	hEneVsZ[i]->GetYaxis()->SetTitle(Form("#gamma"));   
	hEneVsZ[i]->GetZaxis()->SetTitle(Form("Charge [e^{-}]"));
      }
      
      if(cMax>0.0) {
	hEneVsZ[i]->GetZaxis()->SetRangeUser(0.0,cMax);
      }
      
      //hCurr[i] = hEneVsZ[i]->ProjectionX(Form("hCurr_%1i",i),1,NeBin);
      hCurr[i] = hEneVsZ[i]->ProjectionX(Form("hCurr_%1i",i));
      if(opt.Contains("units")) {
	Double_t c = PConst::c_light / (spaUnit / PUnits::fs);
	hCurr[i]->Scale((deu*c));
	hCurr[i]->ResetStats();
      }
      if(hCurr[i]->GetMaximum()>maxCur && i<3) maxCur = hCurr[i]->GetMaximum();

      // hSpec[i] = hEneVsZ[i]->ProjectionY(Form("hSpec_%1i",i),1,Nx1Bin);
      hSpec[i] = hEneVsZ[i]->ProjectionY(Form("hSpec_%1i",i));
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

    if(Nspecies>1) {
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
    } else {
      exPal[0] = new TExec("exBeam","beamPalette->cd();");
      colors[0] = kGray+1;
      fillstyle[0] = 1001;
      fillcolor[0] = TColor::GetColor(230,230,230); 
    }
    
    TPaletteAxis *palette = NULL;
        
    // Canvas setup
    Int_t sizex = 1024;
    Int_t sizey = 640;
    TCanvas *C = new TCanvas("C","Beams energy",sizex,sizey);
    C->SetFillStyle(4000);
    gPad->SetTickx(1);
    gPad->SetTicky(1);
    if(opt.Contains("logz")) {
      gPad->SetLogz(1);
      hEneVsZjoint->GetZaxis()->SetLabelOffset(-0.002);
      hEneVsDivjoint->GetZaxis()->SetLabelOffset(-0.01);
    }
      
    gPad->SetFrameFillColor(beamPalette->GetColor(0));

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
      } else {
	if (Nspecies==1) 
	  exPal[0]->Draw();
	else if (Nspecies>1)
	  exPal[1]->Draw();
	

	hEneVsZjoint->GetZaxis()->SetLabelFont(43);
	hEneVsZjoint->GetZaxis()->SetLabelSize(32);
	// hEneVsZjoint->GetZaxis()->SetLabelOffset(0.8);
	hEneVsZjoint->GetZaxis()->SetTitleFont(43);
	hEneVsZjoint->GetZaxis()->SetTitleSize(36);
	hEneVsZjoint->GetZaxis()->SetTitleOffset(1.2);
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
    
    // Current axis
    Float_t x1pos = xmax - (xmax-xmin) * 0.1;
    TGaxis *axis = new TGaxis(x1pos,ymin,
			      x1pos,eMin,
			      0.0,curmax,503,"S+L");
    axis->SetLineWidth(1);
    axis->SetLineColor(kGray+3);//PGlobals::elecLine);
    axis->SetLabelColor(kGray+3);//PGlobals::elecLine);
    axis->SetLabelFont(PGlobals::fontType);
    axis->SetLabelSize(PGlobals::labelSize-10);
    axis->SetLabelOffset(0.005);
    axis->SetTitleColor(kGray+3);//PGlobals::elecLine);
    axis->SetTitleFont(PGlobals::fontType);
    axis->SetTitleSize(PGlobals::titleSize-10);
    axis->SetTitleOffset(0.6);
    axis->SetTickSize(0.03);
    axis->ChangeLabel(1,-1,0.);
    if(opt.Contains("units"))
      axis->SetTitle(Form("I [%s]",curSUnit.c_str()));
    axis->CenterTitle();
    //axis->ChangeLabel(1,-1,-1,-1,-1,-1,"");
    // axis->SetMaxDigits(2);
    axis->Draw();

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
    TString fOutName = Form("./%s/Plots/Spectrum/Phasespace-%s_%i",pData->GetPath().c_str(),pData->GetName(),time);
  
    PGlobals::imgconv(C,fOutName,opt);
    // ---------------------------------------------------------


    // Canvas setup
    sizex = 1024;
    sizey = 320;
    TCanvas *C2 = new TCanvas("C2","Spectrum",sizex,sizey);
    C->SetFillStyle(4000);
    gPad->SetFillStyle(4000);

    if(opt.Contains("nopal"))
      gPad->SetRightMargin(0.02);
    else
      gPad->SetRightMargin(0.14);
      
    gPad->SetLeftMargin(0.10);    
    gPad->SetBottomMargin(0.25);
    gPad->SetTickx(1);
    gPad->SetTicky(1);
    if(opt.Contains("logz"))
      gPad->SetLogz(1);

    Int_t linewidth = 5;
    gPad->SetFrameLineWidth(linewidth);
    Int_t color = kWhite;
    Int_t color2 = kWhite;
    if (opt.Contains("hzdr"))  {
      color  = kBlack;
      color2 = 46;
      linewidth = 3;
      gPad->SetFrameLineWidth(linewidth);
      // gPad->SetFrameLineColor(kWhite);
    } else if (opt.Contains("ffwd"))  {
      color  = kBlack;
      color2 = kOrange-3;
      linewidth = 3;
      gPad->SetFrameLineWidth(linewidth);
      // gPad->SetFrameLineColor(kWhite);
    }

 
    // gPad->SetLogx(1);
    // gPad->SetFrameLineWidth(3);
    // hEneVsDivjoint->GetXaxis()->SetMoreLogLabels(1);
    hEneVsDivjoint->GetXaxis()->SetNdivisions(510);
    hEneVsDivjoint->GetXaxis()->SetAxisColor(color);
    hEneVsDivjoint->GetYaxis()->SetAxisColor(color);
    hEneVsDivjoint->GetXaxis()->SetTickLength(0.03);
    
    hEneVsDivjoint->GetXaxis()->CenterTitle();
    hEneVsDivjoint->GetXaxis()->SetTitleOffset(1.0);
    hEneVsDivjoint->GetXaxis()->SetTitleSize(32);
    hEneVsDivjoint->GetXaxis()->SetTitle("Energy [MeV]");
    hEneVsDivjoint->GetYaxis()->CenterTitle();
    hEneVsDivjoint->GetYaxis()->SetTitleOffset(0.4);
    hEneVsDivjoint->GetYaxis()->SetTitleSize(32);
    hEneVsDivjoint->GetYaxis()->SetTitle("#theta [mrad]");
    hEneVsDivjoint->GetZaxis()->CenterTitle();
    hEneVsDivjoint->GetZaxis()->SetTitleOffset(0.5);
    hEneVsDivjoint->GetZaxis()->SetTitle("pC/MeV/mrad");
    hEneVsDivjoint->GetZaxis()->SetTitleSize(28);

    UInt_t colorz = kBlack;
    if (opt.Contains("wpal")) 
      colorz = kWhite;
    
    hEneVsDivjoint->GetZaxis()->SetAxisColor(colorz);
    hEneVsDivjoint->GetZaxis()->SetLabelColor(colorz);
    hEneVsDivjoint->GetZaxis()->SetTitleColor(colorz);

    PPalette * specPalette = (PPalette*) gROOT->FindObject("spec");
    if(!specPalette) {
      specPalette = new PPalette("spec");
    }
    // specPalette->SetPalette(60);
    specPalette->SetPalette("spectrum");
    // specPalette->Invert();
    if (opt.Contains("hzdr")) 
      specPalette->SetPalette("spectrum1");
    else if (opt.Contains("ffwd")) 
      specPalette->SetPalette("screen");
    
    
    gPad->SetFrameFillColor(specPalette->GetColor(0));
    //    hEneVsDivjoint->SetMinimum(-1.0000);
    
    // cout << Form("dmax = %f",hEneVsDivjoint->GetMaximum()) << endl;
    if(dMax >= 0.0) {
      hEneVsDivjoint->SetMaximum(dMax);
      // hEneVsDivjoint->GetZaxis()->SetRangeUser(30,dMax);
      // hEneVsDivjoint->GetZaxis()->SetRangeUser(10,dMax);
    }
    cout << Form(" Spectrum maximum = %f",hEneVsDivjoint->GetMaximum()) << endl;
  
    
    if(dMin >= 0.0) {
      hEneVsDivjoint->SetMinimum(dMin);
    }
        
    hEneVsDivjoint->GetXaxis()->SetRangeUser(0.0,eMax);

    if(opt.Contains("nopal"))
      hEneVsDivjoint->Draw("col 0");
    else
      hEneVsDivjoint->Draw("colz 0");
      
    // Save histograms
    TString ofilename = Form("./%s/Plots/Spectrum/Spectrum-%s_%i.root",pData->GetPath().c_str(),pData->GetName(),time);
    TFile *ofile = new TFile(ofilename,"RECREATE");
    hEneVsDivjoint->Write("hSpectrum2D",TObject::kOverwrite);
    ofile->Close();
    // ----------------

    gPad->Update();

    y1 = gPad->GetBottomMargin();
    y2 = 1 - gPad->GetTopMargin();
    x2 = 1 - gPad->GetRightMargin();
    TPave *pFrame = new TPave((x2 + 0.01),y1,(x2 + 0.03),y2,1,"NDCL");
    pFrame->SetFillStyle(0);
    pFrame->SetLineColor(colorz);
    pFrame->SetShadowColor(0);

    palette = (TPaletteAxis*) hEneVsDivjoint->GetListOfFunctions()->FindObject("palette");
    if(palette) {
      palette->SetY1NDC(y1);
      palette->SetY2NDC(y2);
      palette->SetX1NDC(x2 + 0.01);
      palette->SetX2NDC(x2 + 0.03);
      palette->SetBorderSize(2);
      palette->SetLineColor(color);
      
      pFrame->Draw();
    }
      
    TH1D *hProj = NULL;
    if (opt.Contains("outl")) {
      hProj = hEneVsDivjoint->ProjectionX("hEneProj",1,hEneVsDivjoint->GetNbinsY());
      hProj->Scale(dvu);
      hProj->ResetStats();
      hProj->SetLineWidth(linewidth);
      hProj->SetLineColor(color2);
      Float_t smax = hProj->GetMaximum();
      if(sMax >= 0.0) 
	smax = sMax;
      cout << Form(" Maximum = %f",smax) << endl;

      Float_t Peak = hProj->GetMaximum();
      Float_t Emean = hProj->GetMean();
      Float_t Erms  = hProj->GetRMS();
      Float_t Qgauss = Peak * TMath::Sqrt(TMath::TwoPi()) * Erms ;
      Float_t Qtotal = hProj->Integral("width");
      
      cout << Form(" Peak   = %f pC/MeV",Peak) << endl;
      cout << Form(" Mean   = %f MeV",Emean) << endl;
      cout << Form(" Rms    = %f MeV",Erms) << endl;
      cout << Form(" Charge = %f pC",Qtotal) << endl;
      cout << Form(" C (g)  = %f pC",Qgauss) << endl;
      
      Float_t axismax = diveMin+(diveMax-diveMin)/3.0;
      Float_t slope = (axismax - diveMin)/smax;
      for(Int_t j=1;j<=hProj->GetNbinsX();j++) {
	Float_t value = hProj->GetBinContent(j);
	hProj->SetBinContent(j,slope*value+diveMin);
      }
      
      // Current axis
      Float_t x1pos = eMax - (eMax-eMin) * 0.10;
      TGaxis *axis = new TGaxis(x1pos,diveMin,x1pos,axismax,0.0,smax,203,"S+L");
      axis->SetLineWidth(2);
      axis->SetLineColor(color2);//PGlobals::elecLine);
      axis->SetLabelColor(color2);//PGlobals::elecLine);
      axis->SetLabelFont(PGlobals::fontType);
      axis->SetLabelSize(PGlobals::labelSize-12);
      axis->SetLabelOffset(0.005);
      axis->SetTitleColor(color2);//PGlobals::elecLine);
      axis->SetTitleFont(PGlobals::fontType);
      axis->SetTitleSize(PGlobals::titleSize-18);
      axis->SetTitleOffset(0.4);
      axis->SetTickSize(0.03);
      axis->ChangeLabel(1,-1,0.);
      if(opt.Contains("units"))
	axis->SetTitle(Form("pC/MeV"));
      //axis->CenterTitle();
      //axis->ChangeLabel(1,-1,-1,-1,-1,-1,"");
      //axis->SetMaxDigits(2);
      axis->Draw();
      
      hProj->GetXaxis()->SetRangeUser(eMinu,eMaxu);
      hProj->Draw("hist L same");
    }

    gPad->Update();
    lFrame = new TBox(gPad->GetUxmin(), gPad->GetUymin(),
		      gPad->GetUxmax(), gPad->GetUymax());
    lFrame->SetFillStyle(0);
    lFrame->SetLineColor(color);
    if (opt.Contains("hzdr") || opt.Contains("ffwd"))  
      lFrame->SetLineWidth(linewidth);
    else
      lFrame->SetLineWidth(linewidth-3);
      
    lFrame->Draw();

    xmin = eMinu;
    xmax = eMaxu;
    ymin = diveMin;
    ymax = diveMax;
    if(opt.Contains("units") && n0) 
      sprintf(ctext,"z = %5.1f %s", Time * skindepth / propUnit, propSUnit.c_str());
    else
      sprintf(ctext,"k_{p}z = %5.1f",Time);
    TLatex *texTime2 = new TLatex(xmax - (xmax-xmin)/20.,ymax - (ymax-ymin)/15.,ctext);
    texTime2->SetTextFont(PGlobals::fontType);
    texTime2->SetTextSize(PGlobals::titleSize-10);
    texTime2->SetTextAlign(33);
    texTime2->SetTextColor(color);
    texTime2->Draw();
    
    gPad->RedrawAxis();
    gPad->RedrawAxis("g");
    
    // Print to a file
    // Output file
    fOutName = Form("./%s/Plots/Spectrum/Spectrum-%s_%i",pData->GetPath().c_str(),pData->GetName(),time);
  
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
