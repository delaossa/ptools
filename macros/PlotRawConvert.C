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
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TPaletteAxis.h>
#include <TExec.h>
#include <TRandom3.h>

#include "PData.hh"
#include "PGlobals.hh"
#include "PPalette.hh"

void PlotRawConvert( const TString &sim, Int_t time, Int_t index = 0, Float_t factor = 1.0, const TString &options="") {
  
#ifdef __CINT__  
  gSystem->Load("libptools.so");
#endif

  PGlobals::Initialize();

  TString opt = options;
 
  // Palettes!
  // gROOT->Macro("PPalettes.C");

  if(opt.Contains("grid")) {
    gStyle->SetPadGridX(1);
    gStyle->SetPadGridY(1);
  }
  gStyle->SetTitleAlign(22);
  gStyle->SetPadRightMargin(0.17);   // Margin for palettes in 2D histos
  //gStyle->SetTitleOffset(0.9,"z");
  gStyle->SetLabelFont(43,"xyz");
  gStyle->SetTitleFont(43,"xyz");

  gStyle->SetLabelSize(24, "xyz");
  gStyle->SetTitleSize(36, "xyz");

  gStyle->SetLabelOffset(0.01, "xyz");
  gStyle->SetTitleOffset(1.6,"xyz");

  // Load PData
  PData *pData = PData::Get(sim.Data());
  pData->LoadFileNames(time);
  if(!pData->IsInit()) return;

  // Some plasma constants
  Double_t n0 = pData->GetPlasmaDensity();
  Double_t kp = pData->GetPlasmaK();
  Double_t skindepth = 1.;
  if(kp!=0.0) skindepth = 1/kp;
  Double_t E0 = pData->GetPlasmaE0();

  // Time in OU
  Float_t Time = pData->GetRealTime();
  // z start of the plasma in normalized units.
  Float_t zStartPlasma = pData->GetPlasmaStart()*kp;
  // z start of the beam in normalized units.
  Float_t zStartBeam = pData->GetBeamStart()*kp;
  
  // opt += "comovcenter";

  // Centering time and z position:
  Double_t shiftz = 0.0;
  if(opt.Contains("center")) {
    Time -= zStartPlasma;
    if(opt.Contains("comov")) {     // Centers on the head of the beam.
      Time += zStartBeam;
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
  

  // Pointer to data TTree
  TTree *tree = pData->GetRawTree(pData->GetRawFileName(index)->c_str(),opt);
  
  // Get phasespace histos
  TH2F *hPzvsZ = NULL;
 
  // Bining, intervals, labels, etc.
  Int_t xNbin = 200;
  Int_t yNbin = 200;
  
  TCut Cut = "TMath::Abs(q)";
  
  char hName[24];
  sprintf(hName,"hPzvsZ");
  //hPzvsZ = new TH2F(hName,"",zNbin,x1Min,x1Max,yNbin,0,x2Max);
  char dCom[128];
  sprintf(dCom,"p1:(x1-%f)>>%s(%i,%i)",shiftz,hName,xNbin,yNbin);
  //  tree->Project(hName,dCom,Cut);
  tree->Draw(dCom,Cut,"goff");
  hPzvsZ = (TH2F*) gROOT->FindObject(hName);    
  
  hPzvsZ->GetXaxis()->CenterTitle();
  hPzvsZ->GetYaxis()->CenterTitle();
  hPzvsZ->GetZaxis()->CenterTitle();

  hPzvsZ->GetXaxis()->SetTitle("#zeta [k_{p}^{-1}]");
  hPzvsZ->GetYaxis()->SetTitle("p_{z} [mc]");
  hPzvsZ->GetZaxis()->SetTitle("dn/d#zetadp_{z}");
  
  Float_t integral = hPzvsZ->Integral();
  Float_t zMean = hPzvsZ->GetMean(1);
  Float_t zRms  = hPzvsZ->GetRMS(1);
  Float_t yMean = hPzvsZ->GetMean(2);
  Float_t yRms  = hPzvsZ->GetRMS(2);

  cout << Form(" Original distribution: ") << endl; 
  cout << Form(" N = %10i  Integral = %7.2f  Mean = %7.2f  Rms = %f",(Int_t)tree->GetEntries(),integral,zMean,zRms) << endl;
  
  // Info for the plots
  char ctext[128];
  TPaveText *textInfo = new TPaveText(0.55,0.55,0.82,0.8,"NDC");
  PGlobals::SetPaveTextStyle(textInfo,33); 
  textInfo->SetTextColor(kGray+2);
  textInfo->SetTextFont(43);
  textInfo->SetTextSize(18);
  sprintf(ctext,"Integral = %5.2f",integral);
  textInfo->AddText(ctext);
  sprintf(ctext,"#LT#zeta#GT = %5.2f k_{p}^{-1}",zMean);
  textInfo->AddText(ctext);
  sprintf(ctext,"#LT#zeta#GT_{rms} = %5.2f k_{p}^{-1}",zRms);
  textInfo->AddText(ctext);
  sprintf(ctext,"#LTp_{z}#GT = %5.2f mc",yMean);
  textInfo->AddText(ctext);
  sprintf(ctext,"#LTp_{z}#GT_{rms} = %5.2f mc",yRms);
  textInfo->AddText(ctext);


  // Get the maximum weight of the sample
  sprintf(hName,"hWeight");
  tree->Project(hName,"TMath::Abs(q)");
  TH1F *hWeight = (TH1F*) gROOT->FindObject(hName);      
  Int_t bin = -1;
  for(Int_t j=hWeight->GetNbinsX()-1;j>=0;j--) 
    if(bin==-1 && hWeight->GetBinContent(j+1)>1e-8) {
      bin = hWeight->GetBin(j+1); break;
    }
  Float_t maxW = hWeight->GetBinLowEdge(bin+1);
  
  // Convert the TTree
  TString filename = Form("./%s/Plots/RawConvert/RawConvert-%s",sim.Data(),sim.Data());
  filename += Form("_%i.root",time);

  TFile * ifile = (TFile*) gROOT->GetListOfFiles()->FindObject(filename);
  // if doesn't exist the directory should be created
  if (!ifile) {
    TString f = filename;
    TString dir2 = f.Remove( f.Last( '/' ), f.Length() - f.Last( '/' ) );
    TString dir1 = f.Remove( f.Last( '/' ), f.Length() - f.Last( '/' ) );
    gSystem->mkdir( dir1 );
    gSystem->mkdir( dir2 );
    ifile = new TFile(filename,"RECREATE");
  }  
  
  UInt_t Nvar = 7;
  if(!pData->Is3D()) Nvar = 6;

  Float_t *var = new Float_t[Nvar];  
  char varname[7][4] = {{"p1"},{"p2"},{"p3"},{"q"},{"x1"},{"x2"},{"x3"}};
  
  TTree *treeNew = new TTree("RawTreeNew","");
  for(UInt_t i=0;i<Nvar;i++) {
    tree->SetBranchAddress(varname[i],&var[i]);
    

    char varinfo[6];
    sprintf(varinfo,"%s/F",varname[i]);
    treeNew->Branch(varname[i],&var[i],varinfo);
  } 
  
  
  // Main loop for randomizing and converting
  filename = Form("./%s/Plots/RawConvert/RawConvert-%s",sim.Data(),sim.Data());
  filename += Form("_%i.raw",time);
  ofstream fData;
  fData.open(filename.Data(),ios::out);

  if(!opt.Contains("astra")) {
    for(UInt_t i=0;i<Nvar;i++) {
      fData << Form("%10s ",varname[i]); 
    }
    fData << endl << Form(" -------------------------------------------------------------------------------------- ") << endl;
  }

  // Converting charge. Base volume element is needed:
  Double_t dx1, dx2, dx3;
    
  dx1 = pData->GetDX(0);
  dx2 = pData->GetDX(1);    
  if(pData->Is3D())
    dx3 = pData->GetDX(2);
  else
    dx3 = 1.0;

  Double_t dV = (dx1*skindepth)*(dx2*skindepth)*(dx3*skindepth);
  Double_t Q0 = n0 * dV;
  Double_t QTotal = 0.00;
  Double_t QTotalNew = 0.00;
  
  Int_t nentries = (Int_t)tree->GetEntries();  
  TRandom3 *ran = new TRandom3(0);
  for(Int_t i=0;i<nentries;i++) {
    tree->GetEntry(i);

    if( ((i % 100000) == 0) && (i!=0)  )
      cout << i << " macroparticles processed " << endl;

    QTotal += var[3] *  Q0 * PConst::ElectronCharge/PUnits::picocoulomb;

    // Random number between 0. and "q".
    Float_t rw = ran->Uniform(maxW/factor);
    
    //    cout << Form(" wi = %f    rw = %f   maxW = %f ",var[3],rw,maxW) << endl;

    if(rw>fabs(var[3])) continue;
    
    //    cout << "eeooo" << endl;

    var[3] = -maxW/factor;

    QTotalNew += var[3] *  Q0 * PConst::ElectronCharge/PUnits::picocoulomb;
    
    var[5] -= shiftz;
    treeNew->Fill();
    
    // Damping to ASCII file:
    if(!opt.Contains("astra")) {
      for(UInt_t j=0;j<Nvar;j++) {
	fData << Form("%10.4f ",var[j]);
      }
      fData << endl;
    } else {
      fData << Form("%12.6f   %12.6f   %12.6f   %12.6f   %12.6f   %12.6f   %12.6f   %12.6f   %12.6f ",
		    var[6] * skindepth / PUnits::um,
		    var[5] * skindepth / PUnits::um,
		    (var[4]-zMean) * skindepth / PUnits::um,
		    var[2] * pData->GetBeamMass() / PUnits::MeV,
		    var[1] * pData->GetBeamMass() / PUnits::MeV,
		    var[0] * pData->GetBeamMass() / PUnits::MeV,
		    -999.0,
		    var[3] * Q0 * PConst::ElectronCharge/PUnits::femptocoulomb,
		    0.0) << endl;
    }
  }
  
  fData.close();

  cout << Form(" Total charge = %10.4f pC ", QTotal) << endl;

  sprintf(hName,"hPzvsZnew");

  TH2F *hPzvsZnew = (TH2F*) hPzvsZ->Clone(hName);
  hPzvsZnew->Reset();

  sprintf(dCom,"p1:(x1-%f)>>%s",shiftz,hName);
  treeNew->Draw(dCom,Cut,"goff");
  hPzvsZnew = (TH2F*) ifile->FindObject(hName);    
  
  hPzvsZnew->GetXaxis()->CenterTitle();
  hPzvsZnew->GetYaxis()->CenterTitle();
  hPzvsZnew->GetZaxis()->CenterTitle();

  hPzvsZnew->GetXaxis()->SetTitle("#zeta [k_{p}^{-1}]");
  hPzvsZnew->GetYaxis()->SetTitle("p_{z} [mc]");
  hPzvsZnew->GetZaxis()->SetTitle("dn/d#zetadp_{z}");
  
  integral = hPzvsZnew->Integral();
  zMean = hPzvsZnew->GetMean(1);
  zRms  = hPzvsZnew->GetRMS(1);
  yMean = hPzvsZnew->GetMean(2);
  yRms  = hPzvsZnew->GetRMS(2);

  cout << Form(" New distribution: ") << endl; 
  cout << Form(" N = %10i  Integral = %7.2f  Mean = %7.2f  Rms = %f",(Int_t)treeNew->GetEntries(),integral,zMean,zRms) << endl;
  cout << Form(" Total charge = %10.4f pC ", QTotalNew) << endl;
  
  
  TPaveText *textInfoNew = new TPaveText(0.55,0.55,0.82,0.8,"NDC");
  PGlobals::SetPaveTextStyle(textInfoNew,33); 
  textInfoNew->SetTextColor(kGray+2);
  textInfoNew->SetTextFont(43);
  textInfoNew->SetTextSize(18);
  sprintf(ctext,"Integral = %5.2f",integral);
  textInfoNew->AddText(ctext);
  sprintf(ctext,"#LT#zeta#GT = %5.2f k_{p}^{-1}",zMean);
  textInfoNew->AddText(ctext);
  sprintf(ctext,"#LT#zeta#GT_{rms} = %5.2f k_{p}^{-1}",zRms);
  textInfoNew->AddText(ctext);
  sprintf(ctext,"#LTp_{z}#GT = %5.2f mc",yMean);
  textInfoNew->AddText(ctext);
  sprintf(ctext,"#LTp_{z}#GT_{rms} = %5.2f mc",yRms);
  textInfoNew->AddText(ctext);

  // Plotting
  // -----------------------------------------------
    
  // Canvas setup
  TCanvas *C = new TCanvas("C","Phasespaces",800,1200);
      
  // Text objects
  TPaveText *textTime = new TPaveText(0.55,0.85,0.82,0.9,"NDC");
  PGlobals::SetPaveTextStyle(textTime,32);
  textTime->SetTextSize(18);
  if(opt.Contains("units") && pData->GetPlasmaDensity()) 
    sprintf(ctext,"z = %5.1f mm", Time * skindepth / PUnits::mm);
  else
    sprintf(ctext,"t = %5.1f #omega_{p}^{-1}",Time);
  textTime->AddText(ctext);


  TPaveText *textLabel = new TPaveText(0.20,0.20,0.40,0.30,"NDC");
  PGlobals::SetPaveTextStyle(textLabel,32); 
  textLabel->SetTextSize(18);
  sprintf(ctext,"Original distribution.");
  textLabel->AddText(ctext);
  

  TPaveText *textLabelNew = new TPaveText(0.20,0.20,0.40,0.30,"NDC");
  PGlobals::SetPaveTextStyle(textLabelNew,32); 
  textLabelNew->SetTextSize(18);
  textLabelNew->SetTextColor(kRed);
  sprintf(ctext,"New distribution.");
  textLabelNew->AddText(ctext);
 

  // Actual Plotting!
  // ------------------------------------------------------------

  // Output file
  TString fOutName = Form("./%s/Plots/RawConvert/RawConvert",sim.Data());
  fOutName += Form("-%s_%i",sim.Data(),time);

  C->cd();

  C->Divide(1,2);

  C->cd(1);
 
  gPad->SetFrameLineWidth(2);  

  if(opt.Contains("logz")) {
    gPad->SetLogz(1);
  } else {
    gPad->SetLogz(0);
  }

  if(opt.Contains("logy")) {
    gPad->SetLogy(1);
  } else {
    gPad->SetLogy(0);
  }


  // Set palette:
  PPalette * pPalette = (PPalette*) gROOT->FindObject("beam");
  if(!pPalette) {
    pPalette = new PPalette("beam");
    pPalette->SetPalette("elec");
  }
  pPalette->cd();
  
  hPzvsZ->Draw("colz");
  
  gPad->Update();
  
  textTime->Draw();
  textInfo->Draw();
  textLabel->Draw();
 
  gPad->RedrawAxis(); 

  C->cd(2);

  hPzvsZnew->Draw("colz");

  gPad->Update();
  
  textTime->Draw();
  textInfoNew->Draw();
  textLabelNew->Draw();
  
  gPad->RedrawAxis(); 


  C->cd();

  // Print to a file
  PGlobals::imgconv(C,fOutName,opt);
  // ---------------------------------------------------------


  treeNew->Write("RawTreeNew",TObject::kOverwrite);

  
  
  ifile->Close();

}
