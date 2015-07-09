#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>

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
#include "PlasmaGlob.hh"
#include "PPalette.hh"

void PlotRawConvertHiPace( const TString &sim, Double_t time, Int_t index = 0, Float_t factor = 1.0, const TString &options="") {
  
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
  gStyle->SetTitleAlign(22);
  gStyle->SetPadRightMargin(0.17);   // Margin for palettes in 2D histos
  //gStyle->SetTitleOffset(0.9,"z");
  gStyle->SetLabelFont(42,"xyz");
  gStyle->SetTitleFont(42,"xyz");


  
  // Some plasma constants
  Double_t n0 = 5E16 * (1./PUnits::cm3);
  if(sim.Contains("DRI6")) {
    n0 = 1E17 * (1./PUnits::cm3);
  }
  Double_t kp = PFunc::PlasmaWavenumber(n0);
  Double_t skindepth = 1.;
  if(kp!=0.0) skindepth = 1/kp;
  Double_t E0 = PFunc::PlasmaWBField(n0);
  Double_t BeamMass = 0.511 * PUnits::MeV;
  
  cout << endl;
  cout << Form(" - Plasma density   = %8.4e cm3",n0 / (1./PUnits::cm3)) << endl;
  cout << Form(" - Skin depth       = %8.4f um", skindepth / PUnits::um) << endl;
  cout << Form(" - E_0              = %8.4f GV/m", E0 * PUnits::meter / PUnits::GV ) << endl;
  cout << Form(" - Electron mass    = %8.4f MeV", BeamMass/PUnits::MeV ) << endl;
  
  // Time in OU
  Float_t Time = time;
 
  // Spatial resolution
  Float_t dx1 = 0.005;
  Float_t dx2 = 0.04;
  Float_t dx3 = 0.04;
  
  // Spatial coordinates intervals:
  Float_t x1Min = 12.2;
  Float_t x1Max = 14.2;
  Float_t x2Min = -0.1;
  Float_t x2Max =  0.1;
  Float_t x3Min = -0.1;
  Float_t x3Max =  0.1;

  // Momentum coordinates intervals:
  Float_t p1Min =  7000.01;
  Float_t p1Max =  20999.99;
  Float_t p2Min = -5.0;
  Float_t p2Max =  5.0;
  Float_t p3Min = -5.0;
  Float_t p3Max =  5.0;

  if(sim.Contains("DRI6")) {
    // t = 20667.711
    p1Min =  15000.01;
    p1Max =  37999.99;
    x1Min = 11.8;
    x1Max = 14.0;
  }

    // Bining, intervals, labels, etc.
  Int_t xNbin = 120;
  Int_t pNbin = 120;
  Float_t dp1 = (p1Max-p1Min)/pNbin;
  
  cout << Form(" xNbin = %i , dx1 = %.4f  pNbin = %i , dp1 = %4.f",xNbin,dx1,pNbin,dp1) << endl;


  // READ FROM HiPACE
  TString filename = Form("%s/DATA/InjectedBeamPhaseSpace_time_%.3f",sim.Data(),time);
  cout << filename.Data() << endl;
  
  FILE * pFile;
  size_t lSize;
  double * buffer;
  size_t result;
  
  pFile = fopen (filename.Data(),"rb" );
  if (pFile==NULL) {fputs ("File error",stderr); exit (1);}
  
  // obtain file size:
  fseek (pFile , 0 , SEEK_END);
  lSize = ftell (pFile);
  rewind (pFile);

  // allocate memory to contain the whole file:
  buffer = (double*) malloc (lSize);
  if (buffer == NULL) {fputs ("Memory error",stderr); exit (2);}

  // copy the file into the buffer:
  result = fread (buffer,1,lSize,pFile);
  if (result != lSize) {fputs ("Reading error",stderr); exit (3);}
  
  /* the whole file is now loaded in the memory buffer. */

  // terminate
  fclose (pFile);

  TTree *tree = new TTree("HiTree","RAW tree");
  const Int_t Nvar = 9;
  Double_t var[Nvar];  
  char varname[Nvar][8] = {{"x1"},{"x2"},{"x3"},{"p1"},{"p2"},{"p3"},{"q"},{"ipart"},{"iproc"}};
  for(Int_t i=0;i<Nvar;i++) {
    char vartype[8];
    sprintf(vartype,"%s/D",varname[i]);
    tree->Branch(varname[i],&var[i],vartype);
  }

  Int_t Npart=lSize/(Nvar*sizeof(double));
  
  for(int i_part=0; i_part<Npart; i_part++) {
    for(Int_t i_var=0; i_var<Nvar; i_var++) {
      var[i_var]=buffer[i_var+i_part*Nvar];
      //cout << Form ("%e  ",var[i_var]);
    }
    // cout << endl;
    tree->Fill();
  }
  
  free (buffer);
    
  cout << Form("  %i  particles read! " , Npart) << endl;
  

  // -------------------------------------------------------------------------------------------

  
  // Get phasespace histos
 
  char cutString[512];
  sprintf(cutString,"x1 > %.1f && x1 < %.1f && x2 > %.1f && x2 < %.1f && x3 > %.1f && x3 < %.1f",x1Min,x1Max,x2Min,x2Max,x3Min,x3Max); 
  sprintf(cutString,"TMath::Abs(q)"); 

  TCut Cut = cutString;
  //  cout << Form("   (applied cut: \n %s)",cutString) << endl;
 
  TH2F *hPzvsZ = NULL;   
  char hName[24];
  sprintf(hName,"hPzvsZ");
  //hPzvsZ = new TH2F(hName,"",zNbin,x1Min,x1Max,pNbin,0,x2Max);
  char dCom[128];
  sprintf(dCom,"p1:x1>>%s(%i,%i)",hName,xNbin,pNbin);
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
  PlasmaGlob::SetPaveTextStyle(textInfo,32); 
  textInfo->SetTextColor(kGray+2);
  textInfo->SetTextFont(42);
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
  filename = Form("./%s/Plots/RawConvert/RawConvert-%s",sim.Data(),sim.Data());
  filename += Form("_%.3f.root",time);

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
  
  const Int_t NvarOut = 7;
  Double_t varOut[NvarOut];  
  char varOutname[NvarOut][4] = {{"p1"},{"p2"},{"p3"},{"q"},{"x1"},{"x2"},{"x3"}};
  
  TTree *treeNew = new TTree("RawTreeNew","");
  for(Int_t i=0;i<NvarOut;i++) {
    tree->SetBranchAddress(varOutname[i],&varOut[i]);


    char varOutinfo[6];
    sprintf(varOutinfo,"%s/D",varOutname[i]);
    treeNew->Branch(varOutname[i],&varOut[i],varOutinfo);
  } 


  // Main loop for randomizing and converting
  filename = Form("./%s/Plots/RawConvert/RawConvert-%s",sim.Data(),sim.Data());
  filename += Form("_%.3f.raw",time);
  ofstream fData;
  fData.open(filename.Data(),ios::out);

  if(!opt.Contains("astra")) {
    for(Int_t i=0;i<NvarOut;i++) {
      fData << Form("%10s ",varOutname[i]); 
    }
    fData << endl << Form(" -------------------------------------------------------------------------------------- ") << endl;
  }


  Double_t dV = (dx1*skindepth)*(dx2*skindepth)*(dx3*skindepth);
  Double_t Q0 = n0 * dV;
  Double_t QTotal = 0.00;
  Double_t QTotalNew = 0.00;
  
  Int_t nentries = (Int_t)tree->GetEntries();  
  TRandom3 *ran = new TRandom3(0);
  for(Int_t i=0;i<nentries;i++) {
    tree->GetEntry(i);

    if( (i % 10000) == 0 )
      cout << i << " macroparticles processed " << endl;

    QTotal += varOut[3] *  Q0 * PConst::ElectronCharge/PUnits::picocoulomb;

    // Random number between 0. and "q".
    Float_t rw = ran->Uniform(maxW/factor);
    
    //    cout << Form(" wi = %f    rw = %f   maxW = %f ",varOut[4],rw,maxW) << endl;

    if(rw>fabs(varOut[3])) continue;
    
    //    cout << "eeooo" << endl;

    varOut[3] = -maxW/factor;

    QTotalNew += varOut[3] *  Q0 * PConst::ElectronCharge/PUnits::picocoulomb;
    
    treeNew->Fill();
    
    // Damping to ASCII file:
    if(!opt.Contains("astra")) {
      for(Int_t j=0;j<NvarOut;j++) {
	fData << Form("%10.4f ",varOut[j]);
      }
      fData << endl;
    } else {
      fData << Form("%12.6f   %12.6f   %12.6f   %12.6f   %12.6f   %12.6f   %12.6f   %12.6f   %12.6f ",
		    varOut[6] * skindepth / PUnits::um,
		    varOut[5] * skindepth / PUnits::um,
		    varOut[4] * skindepth / PUnits::um,
		    varOut[2] * BeamMass / PUnits::MeV,
		    varOut[1] * BeamMass / PUnits::MeV,
		    varOut[0] * BeamMass / PUnits::MeV,
		    -999.0,
		    varOut[3] * Q0 * PConst::ElectronCharge/PUnits::femptocoulomb,
		    0.0) << endl;
    }
  }
  
  fData.close();

  cout << Form(" Total charge = %10.4f pC ", QTotal) << endl;

  sprintf(hName,"hPzvsZnew");

  TH2F *hPzvsZnew = (TH2F*) hPzvsZ->Clone(hName);
  hPzvsZnew->Reset();

  sprintf(dCom,"p1:x1>>%s",hName);
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
  PlasmaGlob::SetPaveTextStyle(textInfoNew,32); 
  textInfoNew->SetTextColor(kGray+2);
  textInfoNew->SetTextFont(42);
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
  TCanvas *C;
  if(opt.Contains("hres") && !opt.Contains("pdf")) // high resolution for plain graphics output.
    C = new TCanvas("C","Phasespaces",1000,750);
  else
    C = new TCanvas("C","Phasespaces",800,1000);

    
  // Text objects
  TPaveText *textTime = new TPaveText(0.55,0.85,0.82,0.9,"NDC");
  PlasmaGlob::SetPaveTextStyle(textTime,32); 
  if(opt.Contains("units") && n0) 
    sprintf(ctext,"z = %5.1f mm", Time * skindepth / PUnits::mm);
  else
    sprintf(ctext,"t = %5.1f #omega_{p}^{-1}",Time);
  textTime->AddText(ctext);


  TPaveText *textLabel = new TPaveText(0.20,0.20,0.40,0.30,"NDC");
  PlasmaGlob::SetPaveTextStyle(textLabel,32); 
  sprintf(ctext,"Original distribution.");
  textLabel->AddText(ctext);


  TPaveText *textLabelNew = new TPaveText(0.20,0.20,0.40,0.30,"NDC");
  PlasmaGlob::SetPaveTextStyle(textLabelNew,32); 
  textLabelNew->SetTextColor(kRed);
  sprintf(ctext,"New distribution.");
  textLabelNew->AddText(ctext);
 

  // Actual Plotting!
  // ------------------------------------------------------------

  // Output file
  TString fOutName = Form("./%s/Plots/RawConvert/RawConvert",sim.Data());
  fOutName += Form("-%s_%.3f",sim.Data(),time);

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
  PPalette * pPalette = (PPalette*) gROOT->FindObject("electron");
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
  PlasmaGlob::imgconv(C,fOutName,opt);
  // ---------------------------------------------------------


  treeNew->Write("RawTreeNew",TObject::kOverwrite);  
  ifile->Close();

}
