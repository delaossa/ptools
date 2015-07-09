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
#include <TPaletteAxis.h>
#include <TExec.h>

#include "PData.hh"
#include "PlasmaGlob.hh"
#include "PPalette.hh"

void PlotPhaseSpace2D( const TString &sim, Int_t time, UInt_t index = 0, const TString &phaname = "p1x1", const TString &options="") {
  
  PlasmaGlob::Initialize();

  TString opt = options;
 
  // Palettes!
  gROOT->Macro("PlasmaPalettes.C");

  if(opt.Contains("grid")) {
    gStyle->SetPadGridX(1);
    gStyle->SetPadGridY(1);
  }

  gStyle->SetPadRightMargin(0.20);   // Margin for palettes in 2D histos
  gStyle->SetLabelFont(42,"xyz");

  // Load PData
  PData *pData = PData::Get(sim.Data());
  pData->LoadFileNames(time);
  if(!pData->IsInit()) return;

  Bool_t CYL = kFALSE;
  if(sim.Contains("cyl")) CYL = kTRUE; 
    
  Bool_t ThreeD = kFALSE;
  if(sim.Contains("3D")) ThreeD = kTRUE; 

  // Some plasma constants
  Double_t n0 = pData->GetPlasmaDensity();
  Double_t kp = pData->GetPlasmaK();
  Double_t skindepth = 1.;
  if(kp!=0.0) skindepth = 1/kp;
  Double_t E0 = pData->GetPlasmaE0();

  // Some beam properties:
  Float_t Ebeam = pData->GetBeamEnergy() * PUnits::MeV;
  Float_t gamma = Ebeam / PConst::ElectronMassE;
  Float_t vbeam = TMath::Sqrt(1 - 1/(gamma*gamma));
  // cout << Form(" - Bunch gamma      = %8.4f", gamma ) << endl;
  // cout << Form(" - Bunch velocity   = %8.4f c", vbeam ) << endl;
  
  // Time in OU
  Float_t Time = pData->GetRealTime();
  // z start of the plasma in normalized units.
  Float_t zStartPlasma = pData->GetPlasmaStart()*kp;
  // z start of the beam in normalized units.
  Float_t zStartBeam = pData->GetBeamStart()*kp;
  
  if(opt.Contains("center")) {
    Time -= zStartPlasma;
    if(opt.Contains("comov"))      // Centers on the head of the beam.
      Time += zStartBeam;
  } 

  // Labeling...
  char xname[24];
  char yname[24];
  if(phaname.Contains("p1x1")) {
    sprintf(xname,"z");
    sprintf(yname,"p_{z}");
  } else if(phaname.Contains("p2x2")) {
    sprintf(xname,"y");
    sprintf(yname,"p_{y}");
  }

  char xunits[24];
  char yunits[24];
  if(opt.Contains("units")) {
    sprintf(xunits,"#mum");
    sprintf(yunits,"MeV/c");
  } else {
    sprintf(xunits,"c/#omega_{p}");
    sprintf(yunits,"m_{e}/c");
  }
  
  char xaxisname[24];
  char yaxisname[24];
  sprintf(xaxisname,"%s [%s]",xname,xunits);
  sprintf(yaxisname,"%s [%s]",yname,yunits);

  // Get charge density histos
  UInt_t Nspecies = pData->NSpecies();
  TH2F *hPha2D = NULL;
  for(UInt_t i=0;i<Nspecies;i++) {
    
    if(i!=index) continue;
    
    for(UInt_t j=0;j<pData->NPhaseSpaces();j++) {
      if(!pData->GetPhasespaceFileName(i,j)) 
	continue;
      
      if(!strcmp(phaname.Data(),pData->GetPhasespaceName(j).c_str())) {
	char hName[24];
	sprintf(hName,"hPha_%i_%s",i,phaname.Data());
	hPha2D = (TH2F*) gROOT->FindObject(hName);
	if(hPha2D) delete hPha2D;
	hPha2D =  pData->GetPhasespace(i,j,opt);
	hPha2D->SetName(hName);
  	continue;
      }
    }

    hPha2D->GetXaxis()->CenterTitle();
    hPha2D->GetYaxis()->CenterTitle();
    hPha2D->GetZaxis()->CenterTitle();
    hPha2D->GetXaxis()->SetTitle(xaxisname);
    hPha2D->GetYaxis()->SetTitle(yaxisname);
    if(i==0)
      hPha2D->GetZaxis()->SetTitle("#LTn_{e}#GT [n_{0}]");
    else
      hPha2D->GetZaxis()->SetTitle("#LTn_{b}#GT [n_{0}]");
  }

  // Emittance and others
  // we assume that longitudinal variable is comoving for the range.
  Double_t xMin = 15.0;
  Double_t xMax = 35.0;
  Double_t yMin =  1500.0;
  Double_t yMax =  2400.0;
  
  UInt_t NBINX = hPha2D->GetNbinsX();
  UInt_t NBINY = hPha2D->GetNbinsY();
  Double_t dx = hPha2D->GetXaxis()->GetBinWidth(1);
  Double_t dy = hPha2D->GetYaxis()->GetBinWidth(1);
  Double_t xmean = 0.0;
  Double_t ymean = 0.0;
  Double_t xmean2 = 0.0;
  Double_t ymean2 = 0.0;
  Double_t xymean = 0.0;
  Double_t Ntotal = 0.0;
  for(UInt_t i=1;i<=NBINX;i++) {
    Double_t x = hPha2D->GetXaxis()->GetBinCenter(i);
    // if(x<xMin || x>xMax) continue;
    for(UInt_t j=1;j<=NBINY;j++) {
      Double_t y = hPha2D->GetYaxis()->GetBinCenter(j);
      //  if(y<yMin || y>yMax) continue;
      Double_t value = TMath::Abs(hPha2D->GetBinContent(i,j));
      xmean += x*value;
      ymean += y*value;
      xmean2 += x*x*value;
      ymean2 += y*y*value;
      xymean += x*y*value;

      Ntotal += value;
    }
  }
  // hPha2D->GetXaxis()->SetRangeUser(xMin,xMax);
  // hPha2D->GetYaxis()->SetRangeUser(yMin,yMax);
  
  xmean  /= Ntotal;
  ymean  /= Ntotal;
  xmean2 /= Ntotal;
  ymean2 /= Ntotal;
  xymean /= Ntotal;

  Double_t xrms2  = xmean2 - xmean*xmean;
  Double_t yrms2  = ymean2 - ymean*ymean;
  Double_t xyrms2 = xymean - xmean*ymean;
  Double_t xrms  = TMath::Sqrt(xrms2);
  Double_t yrms  = TMath::Sqrt(yrms2);
  Double_t xyrms = TMath::Sqrt(xyrms2); 

  Double_t emittance = TMath::Sqrt( xrms2*yrms2 - xyrms2);

  cout << " Bunch properties: " << endl;
  cout << Form("  xMean = %7.3f   yMean = %7.3f",xmean,ymean) << endl;
  cout << Form("  xRms  = %7.3f   yRms  = %7.3f",xrms,yrms) << endl;
  cout << Form("  Emittance = %7.3f",emittance) << endl;

  // Set palette:
  PPalette * pPalette = (PPalette*) gROOT->FindObject("electron0");
  
  Float_t Max  = hPha2D->GetMaximum();
  Float_t Min  = hPha2D->GetMinimum();
  
  hPha2D->GetZaxis()->SetRangeUser(Min,Max); 
  
  
  // Plotting
  // -----------------------------------------------

  // Canvas setup
  TCanvas *C;
  if(opt.Contains("hres") && !opt.Contains("pdf")) // high resolution for plain grahics output.
    C = new TCanvas("C","2D Phasespace",1000,625);
  else
    C = new TCanvas("C","2D Phasespace",800,500);

    
  // Text objects
  TPaveText *textTime = new TPaveText(0.50,0.85,0.77,0.9,"NDC");
  PlasmaGlob::SetPaveTextStyle(textTime,32); 
  char ctext[128];
  if(opt.Contains("units") && pData->GetPlasmaDensity()) 
    sprintf(ctext,"Z = %5.1f mm", 1e3 * pData->GetPlasmaSkinDepth() * Time);
  else
    sprintf(ctext,"T = %5.1f #omega_{p}^{-1}",Time);
  textTime->AddText(ctext);
 
  TPaveText *textDen = new TPaveText(0.15,0.85,0.48,0.9,"NDC");
  PlasmaGlob::SetPaveTextStyle(textDen,12); 
  textDen->SetTextColor(kOrange+10);
  if(opt.Contains("units") && pData->GetPlasmaDensity())
    sprintf(ctext,"n_{0} = %5.2f x 10^{15} / cc", 1e-6 * 1e-15 * pData->GetPlasmaDensity());
  else if(pData->GetBeamDensity() && pData->GetPlasmaDensity())
    sprintf(ctext,"n_{b}/n_{0} = %5.2f", pData->GetBeamDensity()/pData->GetPlasmaDensity());
  textDen->AddText(ctext);

  TPaveText *textWav = new TPaveText(0.15,0.85,0.48,0.9,"NDC");
  PlasmaGlob::SetPaveTextStyle(textWav,12); 
  textWav->SetTextColor(kGray+2);
  sprintf(ctext,"#lambda_{p} = %5.3f mm", 1e3 * pData->GetPlasmaWaveLength());
  textWav->AddText(ctext);

  TPaveText *textEmit = new TPaveText(0.50,0.65,0.77,0.80,"NDC");
  PlasmaGlob::SetPaveTextStyle(textEmit,32); 
  textEmit->SetTextColor(kGray+1);
  sprintf(ctext,"#LT%s#GT_{rms} = %5.3f %s",xname,xrms,xunits);
  textEmit->AddText(ctext);
  sprintf(ctext,"#LT%s#GT_{rms} = %5.3f %s",yname,yrms,yunits);
  textEmit->AddText(ctext);
  sprintf(ctext,"#epsilon_{N} = %5.3f #mum",emittance);
  textEmit->AddText(ctext);
  
  
  // Actual Plotting!
  // ------------------------------------------------------------

  // Output file
  TString fOutName = Form("./%s/Plots/Phase2D-%s/Phase2D-%s",sim.Data(),phaname.Data(),phaname.Data());
  fOutName += Form("-%s_%i",sim.Data(),time);

  C->cd();
 
  gPad->SetFrameLineWidth(3);  

  TH2F *hFrame = (TH2F*) gROOT->FindObject("hFrame1");
  if(hFrame) delete hFrame;
  hFrame = (TH2F*) hPha2D->Clone("hFrame1");
  hFrame->Reset();
 
  hFrame->Draw("col");
   
  hPha2D->Draw("colz same");
    
  gPad->Update();
  
  textTime->Draw();
  textDen->Draw();
  if(opt.Contains("units"))
    textWav->Draw();
  textEmit->Draw();
 
  gPad->RedrawAxis(); 

  C->cd();

  // Print to a file
  PlasmaGlob::imgconv(C,fOutName,opt);
  // ---------------------------------------------------------

}
