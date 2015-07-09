#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>

#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TLatex.h>
#include <TPaletteAxis.h>
#include <TExec.h>

#include "PData.hh"
#include "PlasmaGlob.hh"
#include "PPalette.hh"

void PlotRakeInjection( const TString &sim, Int_t zoom=2, Int_t Nbins=2, const TString &options="") {
  
#ifdef __CINT__  
  gSystem->Load("libplasma.so");
#endif

  PlasmaGlob::Initialize();

  // Palettes!
  gROOT->Macro("PlasmaPalettes.C");

  // Init Units table
  PUnits::UnitsTable::Get();
  
  TString opt = options;
  
  // More makeup
  gStyle->SetPadGridY(0);
  if(opt.Contains("grid")) {
    gStyle->SetPadGridX(1);
    gStyle->SetPadGridY(1);
  }
  
  Bool_t CYL = kFALSE;
  if(sim.Contains("cyl")) { CYL = kTRUE; opt += "cyl"; } 
  
  Bool_t ThreeD = kFALSE;
  if(sim.Contains("3D")) ThreeD = kTRUE; 
  
  // Load PData
  PData *pData = PData::Get(sim.Data());

  // Create the canvas and the pads before the Frame loop
  // Resolution:
  Int_t sizex = 600;
  Int_t sizey = 640;
  if(opt.Contains("hres")) {
    Int_t sizex = 1200;
    Int_t sizey = 1280;    
  }
  
  TCanvas *C1 = new TCanvas("C1","Evolution of Injection",sizex,sizey);
  C1->cd();
  
  // Palettes
  TExec *exPlasma = new TExec("exPlasma","plasmaPalette->cd();");
  TExec *exElec   = new TExec("exElec","electronPalette->cd();");
  TExec *exHot    = new TExec("exHot","hotPalette->cd();");
  
  const Int_t NFrames = 4;
  // Int_t timestep[NFrames] = {18,22,26,34};
  // TString sLabels[NFrames] = {"(a)","(b)","(c)","(d)"};
  //  Int_t timestep[NFrames] = {34,26,22,18};
  Int_t timestep[NFrames] = {500,60,36,31};
  TString sLabels[NFrames] = {"(d)","(c)","(b)","(a)"};
  Float_t Time[NFrames] = {0.0};
  TH2F ***hDen2D = new TH2F**[NFrames];
  Float_t **zDenMax = new Float_t*[NFrames];
  Float_t *zDenMaxTotal = NULL;

  // Text objects
  TPaveText **textTime  = new TPaveText*[NFrames];
  TPaveText **textLabel = new TPaveText*[NFrames];

  TPad **pad = new TPad*[NFrames];
  TH2F *hFrame[NFrames];
 
  //  C1->Divide(1,NFrames);

  // Setup Pad layout:
  Double_t lMargin = 0.15;
  Double_t rMargin = 0.18;
  Double_t bMargin = 0.10;
  Double_t tMargin = 0.04;
  Double_t vSpacing = 0.00; 
  Double_t hStep = (1.-lMargin-rMargin);
  Double_t vStep = (1.- bMargin - tMargin - (NFrames-1) * vSpacing) / NFrames;
  
  Float_t vposd = 0.0;
  Float_t vposu = 0.0;
  Float_t vmard = 0.0;
  Float_t vmaru = 0.0;
  Float_t vfactor = 0.0;
  Float_t hposl = 0.0;
  Float_t hposr = 1.0;
  Float_t hmarl = lMargin;
  Float_t hmarr = rMargin;
  Float_t hfactor = 1.0;

  Double_t n0,kp,skindepth,E0;
  Int_t Nspecies;
  
  for(Int_t k=0;k<NFrames;k++) {
    
    pData->LoadFileNames(timestep[k]);
    if(!pData->IsInit()) continue;
    
    // Some plasma constants
    n0 = pData->GetPlasmaDensity();
    kp = pData->GetPlasmaK();
    skindepth = 1.;
    if(kp!=0.0) skindepth = 1/kp;
    E0 = pData->GetPlasmaE0();
    
    // Time in OU
    Time[k] = pData->GetRealTime();
    // z start of the plasma in normalized units.
    Float_t zStartPlasma = pData->GetPlasmaStart()*kp;
    // z start of the beam in normalized units.
    Float_t zStartBeam = pData->GetBeamStart()*kp;
    
    if(opt.Contains("center")) {
      Time[k] -= zStartPlasma;
      if(opt.Contains("comov"))      // Centers on the head of the beam.
	Time[k] += zStartBeam;
    }
    
    // Get charge density histos
    Nspecies = pData->NSpecies();
    hDen2D[k] = new TH2F*[Nspecies];
    zDenMax[k] = new Float_t[Nspecies];
    
    if(k==0) zDenMaxTotal = new Float_t[Nspecies];
   
    for(Int_t i=0;i<Nspecies;i++) {
      hDen2D[k][i] = NULL;
      if(k==0) zDenMaxTotal[i] = -999.;
      if(!pData->GetChargeFileName(i)) 
	continue;
      
      char hName[24];
      sprintf(hName,"hDen_%i_%i",k,i);
      hDen2D[k][i] = (TH2F*) gROOT->FindObject(hName);
      if(hDen2D[k][i]) delete hDen2D[k][i];
      
      if(!ThreeD)
	hDen2D[k][i] = pData->GetCharge(i,opt);
      else
	hDen2D[k][i] = pData->GetCharge2DSliceZY(i,-1,Nbins,opt+"avg");
      
      hDen2D[k][i]->SetName(hName);
      hDen2D[k][i]->GetXaxis()->CenterTitle();
      hDen2D[k][i]->GetYaxis()->CenterTitle();
      hDen2D[k][i]->GetZaxis()->CenterTitle();

      zDenMax[k][i] = hDen2D[k][i]->GetMaximum();
      if(zDenMax[k][i]>zDenMaxTotal[i])
	zDenMaxTotal[i] = zDenMax[k][i]; 
            
      // Changing to user units: 
      // --------------------------
      if(opt.Contains("units") && n0) {
	Int_t NbinsX = hDen2D[k][i]->GetNbinsX();
	Float_t xMin = skindepth * hDen2D[k][i]->GetXaxis()->GetXmin() / PUnits::um;
	Float_t xMax = skindepth * hDen2D[k][i]->GetXaxis()->GetXmax() / PUnits::um;
	Int_t NbinsY = hDen2D[k][i]->GetNbinsY();
	Float_t yMin = skindepth * hDen2D[k][i]->GetYaxis()->GetXmin() / PUnits::um;
	Float_t yMax = skindepth * hDen2D[k][i]->GetYaxis()->GetXmax() / PUnits::um;
	hDen2D[k][i]->SetBins(NbinsX,xMin,xMax,NbinsY,yMin,yMax);
	
	hDen2D[k][i]->GetXaxis()->SetTitle("#zeta [#mum]");
      } else {
	hDen2D[k][i]->GetXaxis()->SetTitle("k_{p} #zeta");
      }

      // Zoom
      Float_t range    = (hDen2D[k][i]->GetYaxis()->GetXmax() - hDen2D[k][i]->GetYaxis()->GetXmin())/zoom;
      Float_t midPoint = (hDen2D[k][i]->GetYaxis()->GetXmax() + hDen2D[k][i]->GetYaxis()->GetXmin())/2.;
      
      Float_t yMin,yMax;
      if(!CYL) {
	yMin = midPoint-range/2;
	yMax = midPoint+range/2;
      } else {
	yMin = 0.;
	yMax = range;
      }
      hDen2D[k][i]->GetYaxis()->SetRangeUser(yMin,yMax);
      

    }
  }

  Float_t density = 1;
  Float_t Base = density;

  Float_t *Max = new Float_t[Nspecies];
  Float_t *Min = new Float_t[Nspecies];
  
  for(Int_t k=0;k<NFrames;k++) {
    
    // Set Z ranges:
    for(Int_t i=0;i<Nspecies;i++) {
  
      Max[i] = zDenMaxTotal[i];
      Min[i] =  1.01E-1 * Base;
      if(i==2) {
	Min[i] = 1.01E-3 * Base;
	Max[i] = zDenMax[k][i];
      }
      
      hDen2D[k][i]->GetZaxis()->SetRangeUser(Min[i],Max[i]);
      
      
      if(i==0 && k==0) {
	// Dynamic plasma palette
	const Int_t plasmaDNRGBs = 3;
	const Int_t plasmaDNCont = 128;
	Double_t basePos = 0.5;
	if(Max[i]!=Min[i]) {
	  if(opt.Contains("logz")) {
	    Float_t a = 1.0/(TMath::Log10(Max[i])-TMath::Log10(Min[i]));
	    Float_t b = TMath::Log10(Min[i]);
	    basePos = a*(TMath::Log10(Base) - b);
	    //  cout << Form("Min = %f Base = %f Max = %f",a*(TMath::Log10(Min) - b),basePos,a*(TMath::Log10(Max) - b)) << endl;     
	  } else {
	    basePos = (1.0/(Max[i]-Min[i]))*(Base - Min[i]);
	  }
	}
	
	Double_t plasmaDStops[plasmaDNRGBs] = { 0.00, basePos, 1.00 };
	Double_t plasmaDRed[plasmaDNRGBs]   = { 0.99, 0.90, 0.00 };
	Double_t plasmaDGreen[plasmaDNRGBs] = { 0.99, 0.90, 0.00 };
	Double_t plasmaDBlue[plasmaDNRGBs]  = { 0.99, 0.90, 0.00 };
	
	PPalette * plasmaPalette = (PPalette*) gROOT->FindObject("plasma");
	plasmaPalette->CreateGradientColorTable(plasmaDNRGBs, plasmaDStops, 
						plasmaDRed, plasmaDGreen, plasmaDBlue, plasmaDNCont);
      }
    }
    
    // PLOTTING!
    if(k==0) {
      vposd = 0.0;
      vposu = bMargin + vStep;
      vfactor = vposu-vposd;  
      vmard = bMargin / vfactor;
      vmaru = 0.0;
    } else if(k == NFrames-1) {
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
    hfactor = hposl-hposr;
  
    C1->cd();
    
    char name[16];
    sprintf(name,"pad_%i",k);
    pad[k] = new TPad(name,"",hposl,vposd,hposr,vposu);
    // // cout << Form("%f %f %f %f",hposl,vposd,hposr,vposu) << endl;
    // // cout << Form("%f %f %f %f",hmarl,vmard,hmarr,vmaru) << endl;
    pad[k]->SetLeftMargin(hmarl);
    pad[k]->SetRightMargin(hmarr);  
    pad[k]->SetBottomMargin(vmard);
    pad[k]->SetTopMargin(vmaru);
  
    pad[k]->SetFrameLineWidth(3);  
    if(opt.Contains("logz")) {
      pad[k]->SetLogz(1);
    } else {
      pad[k]->SetLogz(0);
    }
    pad[k]->Draw();
    pad[k]->cd();
    
    sprintf(name,"hFrame_%i",k);  
    hFrame[k] = (TH2F*) gROOT->FindObject(name);
    if(hFrame[k]) delete hFrame[k];
    hFrame[k] = (TH2F*) hDen2D[k][0]->Clone(name);
    hFrame[k]->Reset();
  
    hFrame[k]->GetXaxis()->CenterTitle();
    hFrame[k]->GetYaxis()->CenterTitle();
    hFrame[k]->GetZaxis()->CenterTitle();

    hFrame[k]->SetLabelFont(42,"xyz");
    hFrame[k]->SetTitleFont(42,"xyz");
  
    hFrame[k]->SetNdivisions(505,"xyz");
  
    hFrame[k]->SetTickLength(0.04,"xyz");
    hFrame[k]->SetTickLength(0.04*vfactor,"y");
  
    hFrame[k]->GetYaxis()->SetLabelSize(0.03/vfactor);
    hFrame[k]->GetYaxis()->SetLabelOffset(0.02);
  
    hFrame[k]->GetYaxis()->SetTitleSize(0.03/vfactor);
    hFrame[k]->GetYaxis()->SetTitleOffset(999.0*vfactor);
  
    if(k==0) {  
      hFrame[k]->GetXaxis()->SetLabelSize(0.10);
      hFrame[k]->GetXaxis()->SetLabelOffset(0.02);
      hFrame[k]->GetXaxis()->SetTitleSize(0.14);
      hFrame[k]->GetXaxis()->SetTitleOffset(1.0);
    } else {
      hFrame[k]->GetXaxis()->SetLabelSize(0.0);
      hFrame[k]->GetXaxis()->SetTitleSize(0.0);
    }
  
    hFrame[k]->Draw("axis");    
  


    exPlasma->Draw();

    // Sum of histograms!
    hDen2D[k][0]->Add(hDen2D[k][1]);
        
    //    hDen2D[k][0]->GetZaxis()->SetRangeUser(Min[1],Max[1]);
 
    hDen2D[k][0]->Draw("colz same");
    
    pad[k]->Update();
    TPaletteAxis *palette = (TPaletteAxis*)hDen2D[k][0]->GetListOfFunctions()->FindObject("palette");
    
    Float_t y1 = pad[k]->GetBottomMargin();
    Float_t y2 = 1 - pad[k]->GetTopMargin();
    Float_t x1 = pad[k]->GetLeftMargin();
    Float_t x2 = 1 - pad[k]->GetRightMargin();
    palette->SetY2NDC(y2 - 1*(y2-y1)/2.0 - 0.00);
    palette->SetY1NDC(y1 + 0*(y2-y1)/2.0 + 0.00);
    palette->SetX1NDC(x2 + 0.005);
    palette->SetX2NDC(x2 + 0.03);

    palette->SetLabelFont(42);

    palette->SetLabelSize(0.03/vfactor);
    palette->SetLabelOffset(-0.004);
    palette->SetTitleSize(0.03/vfactor);
    palette->SetTitleOffset(9999.0*vfactor);
  
    palette->SetBorderSize(2);
    palette->SetLineColor(1);
  
    exHot->Draw();
    hDen2D[k][2]->Draw("colz same");

    pad[k]->Update();
    palette = (TPaletteAxis*)hDen2D[k][2]->GetListOfFunctions()->FindObject("palette");
    
    palette->SetY2NDC(y2 - 0*(y2-y1)/2.0 - 0.00);
    palette->SetY1NDC(y1 + 1*(y2-y1)/2.0 + 0.00);
    palette->SetX1NDC(x2 + 0.005);
    palette->SetX2NDC(x2 + 0.03);

    palette->SetLabelFont(42);

    palette->SetLabelSize(0.03/vfactor);
    palette->SetLabelOffset(-0.004);
    palette->SetTitleSize(0.03/vfactor);
    palette->SetTitleOffset(9999.0*vfactor);

    palette->SetBorderSize(2);
    palette->SetLineColor(1);

    textTime[k] = new TPaveText(x2-0.17,y2-(0.05/vfactor),x2-0.02,y2-0.04,"NDC"); 
    PlasmaGlob::SetPaveTextStyle(textTime[k],32); 
    char ctext[128];
    if(opt.Contains("units")) 
      sprintf(ctext,"z = %5.0f #mum", skindepth * Time[k] / PUnits::um);
    else
      sprintf(ctext,"z = %5.1f #omega_{p}^{-1}",Time[k]);
    textTime[k]->AddText(ctext);
    textTime[k]->Draw();

    textLabel[k] = new TPaveText(x1+0.02,y1+(0.02/vfactor),x1+0.07,y1+(0.06/vfactor),"NDC"); 
    PlasmaGlob::SetPaveTextStyle(textLabel[k],32); 
    textLabel[k]->SetTextFont(42);
    textLabel[k]->AddText(sLabels[k]);
    textLabel[k]->Draw();
    

    pad[k]->RedrawAxis(); 

  }
    
  C1->cd();

  // PlasmaGlob::SetPaveTextStyle(textYaxis,11); 
  char ctext[128];
  if(opt.Contains("units")) 
    sprintf(ctext,"y [#mum]");
  else
    sprintf(ctext,"k_{p} y");
  TLatex *textYaxis = new TLatex(0.05,0.6,ctext);
  textYaxis->SetTextAlign(22);
  textYaxis->SetTextFont(42);
  textYaxis->SetTextAngle(90);
  textYaxis->Draw();
  
  sprintf(ctext,"n [n_{0}]");
  TLatex *textZaxis = new TLatex(0.95,0.6,ctext);
  textZaxis->SetTextAlign(22);
  textZaxis->SetTextFont(42);
  textZaxis->SetTextAngle(90);
  textZaxis->Draw();

  
  // Output file
  TString fOutName = Form("./%s/Plots/RakeInjection/RakeInjection",pData->GetPath().c_str());
  fOutName += Form("-%s",pData->GetName());
  
  // Print to a file
  PlasmaGlob::imgconv(C1,fOutName,opt);
  // ---------------------------------------------------------

}
