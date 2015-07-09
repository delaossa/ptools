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

void PlotFieldInjection( const TString &sim, Int_t zoom=2, Int_t Nbins=2, const TString &options="") {
  
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
  TExec *exrbow2    = new TExec("exrbow2","rbow2Palette->cd();");
  
  const Int_t NFrames = 4;
  // Int_t timestep[NFrames] = {18,22,26,34};
  // TString sLabels[NFrames] = {"(a)","(b)","(c)","(d)"};
  Int_t timestep[NFrames] = {34,26,22,18};
  //Int_t timestep[NFrames] = {21,20,19,18};
  TString sLabels[NFrames] = {"(d)","(c)","(b)","(a)"};
  Float_t Time[NFrames] = {0.0};
  TH2F ***hE2D = new TH2F**[NFrames];
  TH2F **hETotal2D = new TH2F*[NFrames];
  TH2F **hDenTag2D = new TH2F*[NFrames];
  Float_t *zEMax = new Float_t[NFrames];
  Float_t zEMaxTotal = -99;
  Float_t *zDenMaxTag = new Float_t[NFrames];
  Float_t zDenMaxTagTotal = -99;
 
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

  // Some plasma constants
  pData->LoadFileNames(timestep[0]);
  Double_t n0,kp,skindepth,E0;  
  n0 = pData->GetPlasmaDensity();
  kp = pData->GetPlasmaK();
  skindepth = 1.;
  if(kp!=0.0) skindepth = 1/kp;
  E0 = pData->GetPlasmaE0();
  
  // z start of the plasma in normalized units.
  Float_t zStartPlasma = pData->GetPlasmaStart()*kp;
  // z start of the beam in normalized units.
  Float_t zStartBeam = pData->GetBeamStart()*kp;

  // Read array of selected tags:
  cout << endl;
  cout << Form(" Reading tags ...") << endl;
  TString filename = Form("./%s/Plots/RakeDump/RakeDump-%s_%i.raw",sim.Data(),sim.Data(),100); // time step when selection.
  ifstream fDataIn(filename.Data(),ios::in);
  Int_t Npart;
  fDataIn >> Npart;
  vector<Int_t*> tags;
  Int_t it = 0;
  for(Int_t i = 0; i<Npart; i++) {
    tags.push_back(new Int_t[2]);
    fDataIn >> tags[it][0] >> tags[it][1];
    // cout <<  Form("%10i  %10i",tags[it][0],tags[it][1]) << endl;
    it++;
  } //while (!fDataIn.eof()) ;
  fDataIn.close();
  

  const Int_t Nfields = 3;
  for(Int_t k=0;k<NFrames;k++) {
    
    pData->LoadFileNames(timestep[k]);
    if(!pData->IsInit()) continue;
        
    // Time in OU
    Time[k] = pData->GetRealTime();
    Double_t shiftz = 0.0;    
    if(opt.Contains("center")) {
      Time[k] -= zStartPlasma;
      if(opt.Contains("comov")) {      // Centers on the head of the beam.
	Time[k] += zStartBeam;
	shiftz += zStartBeam;
      } else {
	shiftz += zStartPlasma;
      }
    }
    if(opt.Contains("comov")) {
      Double_t v = pData->GetBeamVelocity();    
      if(v==0) v = 1.0; // If equals to 0 (default), then set to c.
      shiftz += v * pData->GetRealTime();
    }
    
    // Get charge density histos
   
    hE2D[k] = new TH2F*[Nfields];
    for(Int_t i=0;i<Nfields;i++) {
      hE2D[k][i] = NULL;
      
      if(!pData->GetChargeFileName(i)) 
	continue;
      
      char hName[24];
      sprintf(hName,"hE_%i_%i",k,i);
      hE2D[k][i] = (TH2F*) gROOT->FindObject(hName);
      if(hE2D[k][i]) delete hE2D[k][i];
      
      if(!ThreeD)
	hE2D[k][i] = pData->GetEField(i,opt);
      else
	hE2D[k][i] = pData->GetEField2DSliceZY(i,-1,Nbins,opt+"avg");
      
      hE2D[k][i]->SetName(hName);
      hE2D[k][i]->GetXaxis()->CenterTitle();
      hE2D[k][i]->GetYaxis()->CenterTitle();
      hE2D[k][i]->GetZaxis()->CenterTitle();

      // Get just the tagged electrons:
      if(i==2) {
	
	cout << " Getting tagged electrons for time step = " << timestep[k] << endl;
	TTree *tree = pData->GetTreeRaw(pData->GetRawFileName(i)->c_str(),opt);
	
	// Redo the histogram with only tagged particles:
	sprintf(hName,"hDenTag2D_%i",k);
	hDenTag2D[k] = (TH2F*) hE2D[k][i]->Clone("hDenTag2D");
	hDenTag2D[k]->Reset();

	// We need to dump tagged electrons in selected slice.
	Float_t dx1 = hE2D[k][i]->GetXaxis()->GetBinWidth(1);
	Float_t dx2 = hE2D[k][i]->GetYaxis()->GetBinWidth(1);
	Float_t Dx3 = hE2D[k][i]->GetYaxis()->GetBinWidth(1) * Nbins;
	// Volume element
	Float_t dV = dx1 * dx2 * dx2; 
	  
	// Loop for event selection:
	const Int_t Nvar = 8;
	Float_t var[Nvar];  
	char varname[Nvar][4] = {{"ene"},{"p1"},{"p2"},{"p3"},{"q"},{"x1"},{"x2"},{"x3"}};
	for(Int_t i=0;i<Nvar;i++) {
	  tree->SetBranchAddress(varname[i],&var[i]);
	}
	Int_t tag[2];
	tree->SetBranchAddress("tag",&tag);
  

	Int_t nentries = (Int_t)tree->GetEntries();  
	for(Int_t i=0;i<nentries;i++) {
	  tree->GetEntry(i);

	  if(TMath::Abs(var[7])>Dx3) continue; // Select central slice.

	  var[5] = var[5] - shiftz;
	  
	  for(Int_t it = 0; it<(int)tags.size(); it++) {
	    if((tag[0]==tags[it][0]) && (tag[1]==tags[it][1])) {
	      hDenTag2D[k]->Fill(var[5],var[6],TMath::Abs(var[4]));
	      //tags.erase(tags.begin()+it);
	      continue;
	    }
	  }
	}
	
	zDenMaxTag[k] = hDenTag2D[k]->GetMaximum();
	if(zDenMaxTag[k]>zDenMaxTagTotal)
	  zDenMaxTagTotal = zDenMaxTag[k];
      }

      
      // Changing to user units: 
      // --------------------------
      if(opt.Contains("units") && n0) {
	Int_t NbinsX = hE2D[k][i]->GetNbinsX();
	Float_t xMin = skindepth * hE2D[k][i]->GetXaxis()->GetXmin() / PUnits::um;
	Float_t xMax = skindepth * hE2D[k][i]->GetXaxis()->GetXmax() / PUnits::um;
	Int_t NbinsY = hE2D[k][i]->GetNbinsY();
	Float_t yMin = skindepth * hE2D[k][i]->GetYaxis()->GetXmin() / PUnits::um;
	Float_t yMax = skindepth * hE2D[k][i]->GetYaxis()->GetXmax() / PUnits::um;
	hE2D[k][i]->SetBins(NbinsX,xMin,xMax,NbinsY,yMin,yMax);
	for(Int_t j=0;j<NbinsX;j++) {
	  for(Int_t l=0;l<NbinsY;l++) {
	    hE2D[k][i]->SetBinContent(j,l, hE2D[k][i]->GetBinContent(j,l) * ( E0 / (PUnits::GV/PUnits::m) ) );
	  }
	}
	hE2D[k][i]->GetXaxis()->SetTitle("#zeta [#mum]");
      
	if(i==2) {
	  hDenTag2D[k]->SetBins(NbinsX,xMin,xMax,NbinsY,yMin,yMax);
	  hDenTag2D[k]->GetXaxis()->SetTitle("#zeta [#mum]");
	}
	
      } else {
	hE2D[k][i]->GetXaxis()->SetTitle("k_{p} #zeta");
	if(i==2) {
	  hDenTag2D[k]->GetXaxis()->SetTitle("k_{p} #zeta");
	}
      }

      // Zoom
      Float_t range    = (hE2D[k][i]->GetYaxis()->GetXmax() - hE2D[k][i]->GetYaxis()->GetXmin())/zoom;
      Float_t midPoint = (hE2D[k][i]->GetYaxis()->GetXmax() + hE2D[k][i]->GetYaxis()->GetXmin())/2.;
      
      Float_t yMin,yMax;
      if(!CYL) {
	yMin = midPoint-range/2;
	yMax = midPoint+range/2;
      } else {
	yMin = 0.;
	yMax = range;
      }
      hE2D[k][i]->GetYaxis()->SetRangeUser(yMin,yMax);
      

    }
    
    
    // Now, combine the electric field components into the total |E| :
    char hName[24];
    sprintf(hName,"hETotal2D_%i",k);
    hETotal2D[k] = (TH2F*) hE2D[k][0]->Clone("hETotal2D");
    hETotal2D[k]->Reset();
    Int_t NbinsX = hE2D[k][0]->GetNbinsX();
    Int_t NbinsY = hE2D[k][0]->GetNbinsY();
    for(Int_t j=0;j<NbinsX;j++) {
      for(Int_t l=0;l<NbinsY;l++) {
	Float_t E1 = hE2D[k][0]->GetBinContent(j,l);
	Float_t E2 = hE2D[k][1]->GetBinContent(j,l);
	Float_t E3 = hE2D[k][2]->GetBinContent(j,l);
	hETotal2D[k]->SetBinContent(j,l,TMath::Sqrt(E1*E1+E2*E2+E3*E3));
      }
    }
    
    zEMax[k] = hETotal2D[k]->GetMaximum();
    if(zEMax[k]>zEMaxTotal)
      zEMaxTotal = zEMax[k];
  }
  
  

  
  Float_t density = 1;
  // if(opt.Contains("units") && n0)
  //   density = n0 / (1e17/PUnits::cm3);
  Float_t Base = density;
  
  for(Int_t k=0;k<NFrames;k++) {
    
    // Set Z ranges:
    hDenTag2D[k]->GetZaxis()->SetRangeUser(0.01001 * density,zDenMaxTagTotal);
    hETotal2D[k]->GetZaxis()->SetRangeUser(0.01001,zEMaxTotal);
    
    
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
    hFrame[k] = (TH2F*) hE2D[k][0]->Clone(name);
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

    exrbow2->Draw();
    hETotal2D[k]->Draw("colz same");
    
    pad[k]->Update();
    TPaletteAxis *palette = (TPaletteAxis*)hETotal2D[k]->GetListOfFunctions()->FindObject("palette");
    
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
    hDenTag2D[k]->Draw("colz same");
    
    pad[k]->Update();
    // palette = (TPaletteAxis*)hDen2D[k][2]->GetListOfFunctions()->FindObject("palette");
    palette = (TPaletteAxis*)hDenTag2D[k]->GetListOfFunctions()->FindObject("palette");
    
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
  TString fOutName = Form("./%s/Plots/FieldInjection/FieldInjection",pData->GetPath().c_str());
  fOutName += Form("-%s",pData->GetName());
  
  // Print to a file
  PlasmaGlob::imgconv(C1,fOutName,opt);
  // ---------------------------------------------------------

}
