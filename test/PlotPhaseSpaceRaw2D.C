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

void PlotPhaseSpaceRaw2D( const TString &sim, Int_t time, UInt_t index = 0, const TString &phaname = "p1x1", const TString &options="") {
  
#ifdef __CINT__  
  gSystem->Load("libplasma.so");
#endif

  PlasmaGlob::Initialize();

  TString opt = options;
 
  // Palettes!
  gROOT->Macro("PlasmaPalettes.C");

  if(opt.Contains("grid")) {
    gStyle->SetPadGridX(1);
    gStyle->SetPadGridY(1);
  }

  gStyle->SetPadRightMargin(0.15);   // Margin for palettes in 2D histos
  gStyle->SetTitleOffset(0.9,"z");
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
  
  TString opthist = opt;

  // Initial bining and intervals:
  UInt_t xNbin = 250;
  UInt_t yNbin = 250;
  Double_t xMin,yMin,xMax,yMax;
  
  if(phaname.Contains("p1x1")) {
    opthist += "comovcenter";
    if(sim.Contains("TH") || sim.Contains("regae")){
      xMin = -5.5;
      xMax = -3.0;
    } else if(sim.Contains("DP") || sim.Contains("RI")){
      // xMin = -6.0;
      // xMax = -3.0;
      xMin = -4.5;
      xMax = -3.5;
    }
    yMin =  0.0;
    yMax =  6000.0; 
    if(sim.Contains("flash_v0"))
      yMax =  2000.0;   // OJO!!
  } else if(phaname.Contains("p2x2") || phaname.Contains("p3x3")){
    xMin = -1.0;
    xMax =  1.0;
    yMin = -20.0;
    yMax =  20.0;
  } else if(phaname.Contains("x2x1")){
    opthist += "comovcenter";
    if(sim.Contains("TH")){
      xMin = -5.5;
      xMax = -3.0;
    } else if(sim.Contains("DP") || sim.Contains("RI")){
      // xMin = -6.0;
      // xMax = -3.0;
      xMin = -4.5;
      xMax = -3.5;
    }
    yMin = -20.0;
    yMax =  20.0;
  }


  if(opthist.Contains("center")) {
    Time -= zStartPlasma;
    if(opthist.Contains("comov"))      // Centers on the head of the beam.
      Time += zStartBeam;
  } 

  // Labeling...
  char xname[24];
  char yname[24];
  char xlabel[24];
  char ylabel[24];
  char alabel[24];
  char tlabel[24];
  if(phaname.Contains("p1x1")) {
    sprintf(xname,"x1");
    sprintf(yname,"p1");
    sprintf(xlabel,"#zeta");
    sprintf(ylabel,"p_{z}");
  } else if(phaname.Contains("p2x2")) {
    sprintf(xname,"x2");
    sprintf(yname,"p2");
    sprintf(xlabel,"y");
    sprintf(ylabel,"p_{y}");
  } else if(phaname.Contains("p3x3")) {
    sprintf(xname,"x3");
    sprintf(yname,"p3");
    sprintf(xlabel,"x");
    sprintf(ylabel,"p_{x}");
    sprintf(alabel,"#epsilon_{x}");
  } else if(phaname.Contains("x2x1")) {
    sprintf(xname,"x1");
    sprintf(yname,"x2");
    sprintf(xlabel,"#zeta");
    sprintf(ylabel,"y");
  }
  sprintf(tlabel,"time");

  char xunits[24];
  char yunits[24];
  char aunits[24];
  char tunits[24];
  if(opt.Contains("units")) {
    sprintf(xunits,"#mum");
    sprintf(tunits,"mm");
    if(phaname.Contains("p1x1")) {
      sprintf(yunits,"GeV/c");
      sprintf(aunits,"(GeV/c) #mum");
    } else if(phaname.Contains("p2x2")||phaname.Contains("p3x3")) {
      sprintf(yunits,"MeV/c");
      sprintf(aunits,"(MeV/c) #mum");
    } else if(phaname.Contains("x2x1")) {
      sprintf(yunits,"#mum");
      sprintf(aunits,"#mum^{2}");
    }
  } else {
    sprintf(xunits,"c/#omega_{p}");
    if(phaname.Contains("p1x1") || phaname.Contains("p2x2") || phaname.Contains("p3x3")) {
      sprintf(yunits,"m_{e}c");
      sprintf(aunits,"m_{e}c^{2}/#omega_{p}");
    } else if(phaname.Contains("x2x1")) {
      sprintf(yunits,"c/#omega_{p}");
      sprintf(aunits,"(c/#omega_{p})^{2}");    
    }
  }
  
  char xaxisname[24];
  char yaxisname[24];
  sprintf(xaxisname,"%s [%s]",xlabel,xunits);
  sprintf(yaxisname,"%s [%s]",ylabel,yunits);

  // Get phasespace histos
  UInt_t Nspecies = pData->NSpecies();
  TH1F *hPha1Dx = NULL;
  TH1F *hPha1Dy = NULL;
  TH2F *hPha2D = NULL;
  for(UInt_t i=0;i<Nspecies;i++) {
    
    if(i!=index) continue;

    if(!pData->GetRawFileName(i))
      continue;    
    
    char hName[24];
    sprintf(hName,"hPha_%i_%s",i,xname);
    hPha1Dx = (TH1F*) gROOT->FindObject(hName);
    if(hPha1Dx) delete hPha1Dx;
    hPha1Dx = new TH1F(hName,"",xNbin,xMin,xMax);
    // cout << "filename = " << pData->GetRawFileName(i)->c_str() << endl;
    pData->GetH1Raw(pData->GetRawFileName(i)->c_str(),xname,hPha1Dx,opthist);

    // Explore the 1D histogram in the 1st variable:
    // It gets the mean, the rms and the boundaries for plotting
    Double_t xNtotal = 0.0;
    Double_t xmean = 0.0;
    Double_t xmean2 = 0.0;
    Double_t xLeft = -1;
    Double_t xRight = -1;
    Int_t xbinLeft = -1;
    Int_t xbinRight = -1;
    for(UInt_t j=1;j<=xNbin;j++) 
      if(xLeft==-1 && hPha1Dx->GetBinContent(j)>0.1*hPha1Dx->GetMaximum()) {
	xLeft = hPha1Dx->GetBinCenter(j); 
	xbinLeft = j;
	break;
      } 
    
    for(UInt_t j=xNbin;j>0;j--) 
      if(xRight==-1 && hPha1Dx->GetBinContent(j)>0.1*hPha1Dx->GetMaximum()) {
	xRight = hPha1Dx->GetBinCenter(j); 
	xbinRight = j;
	break;
      } 
    
    for(UInt_t k=xbinLeft;k<=xbinRight;k++) {
      Double_t x = hPha1Dx->GetXaxis()->GetBinCenter(k);
      Double_t value = TMath::Abs(hPha1Dx->GetBinContent(k));
      xmean += x*value;
      xmean2 += x*x*value;
      xNtotal += value;
    }
    xmean  /= xNtotal;
    xmean2  /= xNtotal;
    Double_t xrms  = TMath::Sqrt(xmean2 - xmean*xmean);  

    sprintf(hName,"hPha_%i_%s",i,yname);
    hPha1Dy = (TH1F*) gROOT->FindObject(hName);
    if(hPha1Dy) delete hPha1Dy;
    hPha1Dy = new TH1F(hName,"",yNbin,yMin,yMax);
    pData->GetH1Raw(pData->GetRawFileName(i)->c_str(),yname,hPha1Dy,opthist);

    // Explore the 1D histogram in the 2nd variable:
    // It gets the mean the rms and the boundaries for plotting
    Double_t yNtotal = 0.0;
    Double_t ymean = 0.0;
    Double_t ymean2 = 0.0;
    Double_t yLeft = -1;
    Double_t yRight = -1;
    Int_t ybinLeft = -1;
    Int_t ybinRight = -1;
    for(UInt_t j=1;j<=yNbin;j++) 
      if(yLeft==-1 && hPha1Dy->GetBinContent(j)>0.01*hPha1Dy->GetMaximum()) {
	yLeft = hPha1Dy->GetBinCenter(j); 
	ybinLeft = j;
	break;
      } 
    
    for(UInt_t j=yNbin;j>0;j--) 
      if(yRight==-1 && hPha1Dy->GetBinContent(j)>0.01*hPha1Dy->GetMaximum()) {
	yRight = hPha1Dy->GetBinCenter(j); 
	ybinRight = j;
	break;
      } 

    for(UInt_t k=ybinLeft;k<=ybinRight;k++) {
      Double_t y = hPha1Dy->GetXaxis()->GetBinCenter(k);
      Double_t value = TMath::Abs(hPha1Dy->GetBinContent(k));
      ymean += y*value;
      ymean2 += y*y*value;
      yNtotal += value;
    }
    ymean  /= yNtotal;
    ymean2  /= yNtotal;
    Double_t yrms  = TMath::Sqrt(ymean2 - ymean*ymean);  

    
    // Redefine axis range according with the mean and the rms of the 1D distributions:
    TString XNAME = xname;
    // if(!XNAME.Contains("x1")) {
    //  xMin = xLeft  - TMath::Abs(xRight-xLeft)*0.5;
    //  xMax = xRight + TMath::Abs(xRight-xLeft)*0.5;
    // } 

    yMin = yLeft  - TMath::Abs(yRight-yLeft)*0.5;
    yMax = yRight + TMath::Abs(yRight-yLeft)*0.5;
    
    if(phaname.Contains("p2x2") || phaname.Contains("p3x3")){
      if(fabs(xMin)>fabs(xMax))
	xMax = -xMin;
      else
	xMin = -xMax;
      
      if(fabs(yMin)>fabs(yMax))
	yMax = -yMin;
      else
	yMin = -yMax;
      
    }

    // Since the energy is increasing we dynamically set the range of the histogram
    // based on the energy average value. The binning is kept constant.
    // if(phaname.Contains("p1x1")) {
    //   yMin = ymean - 25;
    //   yMax = ymean + 25;
    //   if(yMin<0.0) {
    // 	yMin = 0.0;
    // 	yMax = 50.0;
    //   }
    // }

    sprintf(hName,"hPha_%i_%s",i,phaname.Data());
    hPha2D = (TH2F*) gROOT->FindObject(hName);
    if(hPha2D) delete hPha2D;
    hPha2D = new TH2F(hName,"",xNbin,xMin,xMax,yNbin,yMin,yMax);
    pData->GetH2Raw(pData->GetRawFileName(i)->c_str(),xname,yname,hPha2D,opthist);
   
    hPha2D->GetXaxis()->CenterTitle();
    hPha2D->GetYaxis()->CenterTitle();
    hPha2D->GetZaxis()->CenterTitle();
    hPha2D->GetXaxis()->SetTitle(xaxisname);
    hPha2D->GetYaxis()->SetTitle(yaxisname);
    hPha2D->GetZaxis()->SetTitle("n [a.u.]");
  }

  // Emittance and others 
  Double_t xmean = 0.0;
  Double_t ymean = 0.0;
  Double_t x2mean = 0.0;
  Double_t y2mean = 0.0;
  Double_t xymean = 0.0;
  Double_t Ntotal = 0.0;
  for(UInt_t i=1;i<=xNbin;i++) {
    Double_t x = hPha2D->GetXaxis()->GetBinCenter(i);
    // if(x<xmin || x>xmax) continue;
    for(UInt_t j=1;j<=yNbin;j++) {
      Double_t y = hPha2D->GetYaxis()->GetBinCenter(j);
      // if(y<ymin || y>ymax) continue;
      Double_t value = TMath::Abs(hPha2D->GetBinContent(i,j));
      xmean += x*value;
      ymean += y*value;
      x2mean += x*x*value;
      y2mean += y*y*value;
      xymean += x*y*value;

      Ntotal += value;
    }
  }
  
  xmean  /= Ntotal;
  ymean  /= Ntotal;
  x2mean /= Ntotal;
  y2mean /= Ntotal;
  xymean /= Ntotal;

  Double_t xrms2  = x2mean - xmean*xmean;
  Double_t yrms2  = y2mean - ymean*ymean;
  Double_t xrms   = TMath::Sqrt(xrms2);
  Double_t yrms   = TMath::Sqrt(yrms2);
  Double_t xyrms2 = xymean - xmean*ymean;

  Double_t emittance = TMath::Sqrt(xrms2*yrms2 - xyrms2*xyrms2);

  // Set palette:
  PPalette * pPalette = (PPalette*) gROOT->FindObject("electron");
  pPalette->cd();

  Float_t Max  = hPha2D->GetMaximum();
  Float_t Min  = hPha2D->GetMinimum();
  
  hPha2D->GetZaxis()->SetRangeUser(Min,Max); 

  
  // Chaning to user units: 
  // --------------------------
  
  if(opt.Contains("units") && n0) {
    
    xMin *= skindepth / PUnits::um;
    xMax *= skindepth / PUnits::um;
    if(phaname.Contains("p1x1")) {
      yMin *= pData->GetBeamMass() / PUnits::GeV;
      yMax *= pData->GetBeamMass() / PUnits::GeV;
    } else if(phaname.Contains("p2x2")||phaname.Contains("p3x3")) {
      yMin *= pData->GetBeamMass() / PUnits::MeV;
      yMax *= pData->GetBeamMass() / PUnits::MeV;
    } else if(phaname.Contains("x2x1")) {
      yMin *= skindepth / PUnits::um;
      yMax *= skindepth / PUnits::um;
    }
    hPha2D->SetBins(xNbin,xMin,xMax,yNbin,yMin,yMax);
    // for(UInt_t j=0;j<=xNbin;j++) {
    //   for(UInt_t k=0;k<=yNbin;k++) {
    // 	hPha2D->SetBinContent(j,k, hPha2D->GetBinContent(j,k) * n0 / (1e15/PUnits::cm3) );
    //   }
    // }
    
    xmean *= skindepth / PUnits::um;
    xrms  *= skindepth / PUnits::um;
    if(phaname.Contains("p1x1")) {
      ymean *= pData->GetBeamMass() / PUnits::GeV;
      yrms  *= pData->GetBeamMass() / PUnits::GeV;
      emittance *= (skindepth / PUnits::um) * (pData->GetBeamMass() / PUnits::GeV);
    } else if(phaname.Contains("p2x2")||phaname.Contains("p3x3")) {
      ymean *= pData->GetBeamMass() / PUnits::MeV;
      yrms  *= pData->GetBeamMass() / PUnits::MeV;
      emittance *= (skindepth / PUnits::um) * (pData->GetBeamMass() / PUnits::MeV);
    } else if(phaname.Contains("x2x1")) {
      ymean *= skindepth / PUnits::um;
      yrms  *= skindepth / PUnits::um;
      emittance *= (skindepth / PUnits::um) * (skindepth / PUnits::um);
    }
    
  }

  cout << " Bunch properties: " << endl;
  cout << Form("  xMean = %7.3f   yMean = %7.3f   xyMean = %7.3f",xmean,ymean,xymean) << endl;
  cout << Form("  xRms  = %7.3f   yRms  = %7.3f",xrms,yrms) << endl;
  cout << Form("  Emittance = %7.3f",emittance) << endl;
  
  // Plotting
  // -----------------------------------------------
    
  // Canvas setup
  TCanvas *C;
  if(opt.Contains("hres") && !opt.Contains("pdf")) // high resolution for plain graphics output.
    C = new TCanvas("C","2D Phasespace",1000,625);
  else
    C = new TCanvas("C","2D Phasespace",800,500);

    
  // Text objects
  TPaveText *textTime = new TPaveText(0.55,0.85,0.82,0.9,"NDC");
  PlasmaGlob::SetPaveTextStyle(textTime,32); 
  char ctext[128];
  if(opt.Contains("units") && pData->GetPlasmaDensity()) 
    sprintf(ctext,"Z = %5.1f mm", Time * skindepth / PUnits::mm);
  else
    sprintf(ctext,"T = %5.1f #omega_{p}^{-1}",Time);
  textTime->AddText(ctext);
 
  TPaveText *textDen = new TPaveText(0.15,0.85,0.48,0.9,"NDC");
  PlasmaGlob::SetPaveTextStyle(textDen,12); 
  textDen->SetTextColor(kOrange+10);
  if(opt.Contains("units") && pData->GetPlasmaDensity())
    sprintf(ctext,"n_{0} = %5.2f x 10^{17} / cc", n0 / (1e17/PUnits::cm3));
  else if(pData->GetBeamDensity() && pData->GetPlasmaDensity())
    sprintf(ctext,"n_{b}/n_{0} = %5.2f", pData->GetBeamDensity()/n0);
  textDen->AddText(ctext);

  TPaveText *textWav = new TPaveText(0.15,0.2,0.48,0.25,"NDC");
  PlasmaGlob::SetPaveTextStyle(textWav,12); 
  textWav->SetTextColor(kGray+2);
  sprintf(ctext,"#lambda_{p} = %5.3f #mum", pData->GetPlasmaWaveLength() / PUnits::um);
  textWav->AddText(ctext);

  TPaveText *textEmit = new TPaveText(0.55,0.2,0.82,0.35,"NDC");
  PlasmaGlob::SetPaveTextStyle(textEmit,32); 
  textEmit->SetTextColor(kGray+1);
  sprintf(ctext,"#LT%s#GT_{rms} = %5.3f %s",xlabel,xrms,xunits);
  textEmit->AddText(ctext);
  sprintf(ctext,"#LT%s#GT_{rms} = %5.3f %s",ylabel,yrms,yunits);
  textEmit->AddText(ctext);
  if(!phaname.Contains("x2x1")) {
    sprintf(ctext,"#epsilon_{N} = %5.3f %s",emittance,aunits);
    textEmit->AddText(ctext);
  }
  
  
  // Actual Plotting!
  // ------------------------------------------------------------

  // Output file
  TString fOutName = Form("./%s/Plots/PhaseRaw2D-%s/PhaseRaw2D-%s",sim.Data(),phaname.Data(),phaname.Data());
  fOutName += Form("-%s_%i",sim.Data(),time);

  C->cd();
 
  gPad->SetFrameLineWidth(3);  

  // TH2F *hFrame = (TH2F*) gROOT->FindObject("hFrame1");
  // if(hFrame) delete hFrame;
  // hFrame = (TH2F*) hPha2D->Clone("hFrame1");
  // hFrame->Reset();
 
  // hFrame->Draw("col");
   
  // hPha2D->Draw("colz same");
  hPha2D->Draw("colz");
    
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
