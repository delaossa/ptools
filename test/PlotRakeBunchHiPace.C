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
#include <TProfile.h>
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

void PlotRakeBunchHiPace( const TString &sim, Double_t time, Int_t index = 0, const TString &options="") {
  
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
  gStyle->SetLabelFont(42,"xyz");
  gStyle->SetTextFont(62);
 
  // Some plasma constants
  Double_t n0 = 5E16 * (1./PUnits::cm3);
  if(sim.Contains("DRI6")) {
    n0 = 1E17 * (1./PUnits::cm3);
  } else if(sim.Contains("DRI7")) {
    n0 = 5E16 * (1./PUnits::cm3);
  } else if(sim.Contains("facet_DDR.HiPACE.3D")) {
    n0 = 5E16 * (1./PUnits::cm3);
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

  if(sim.Contains("facet_DDR.HiPACE.3D")) {
    dx1 = 15.0/1024;
    dx2 = 10.0/256;
    dx3 = 10.0/256;
  }
  
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
  } else if(sim.Contains("DRI7")) {
    // t = 71061.387
    p1Min =  30000.01;
    p1Max =  99999.99;
    x1Min = 12.2;
    x1Max = 14.0;
  }

  // Bining, intervals, labels, etc.
  Int_t xNbin = 200;
  Int_t pNbin = 200;
  Float_t dp1 = (p1Max-p1Min)/pNbin;
  
  cout << Form(" xNbin = %i , dx1 = %.4f  pNbin = %i , dp1 = %4.f",xNbin,dx1,pNbin,dp1) << endl;
  char hName[24];
  char dCommand[128];
 
  TH1F *hX1 = NULL;
  TH1F *hP1 = NULL;

  cout << Form("\n1. Getting TTree... ") << endl; 

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
  

  // ----------------------------------------------------------------------------------------------


  char cutString[512];
  sprintf(cutString,"x1 > %.1f && x1 < %.1f && x2 > %.1f && x2 < %.1f && x3 > %.1f && x3 < %.1f",x1Min,x1Max,x2Min,x2Max,x3Min,x3Max); 
  sprintf(cutString,"TMath::Abs(q)"); 

  TCut Cut = cutString;
  //  cout << Form("   (applied cut: \n %s)",cutString) << endl;
  

  // Auto ranges 
  if (opt.Contains("autoz")) {
    sprintf(hName,"hX1");
    hX1 = (TH1F*) gROOT->FindObject(hName);
    if(hX1) delete hX1;
    sprintf(dCommand,"x1>>%s",hName);
    tree->Draw(dCommand,Cut,"goff");
    hX1 = (TH1F*) gROOT->FindObject(hName);
    cout << Form("   - x1. ") << endl;    
    x1Min = hX1->GetXaxis()->GetXmin();
    x1Max = hX1->GetXaxis()->GetXmax();


    sprintf(hName,"hP1");
    hP1 = (TH1F*) gROOT->FindObject(hName);
    if(hP1) delete hP1;
    hP1 = new TH1F(hName,"",pNbin,p1Min,p1Max);
    
  } else if (opt.Contains("autop")) {
    sprintf(hName,"hP1");
    hP1 = (TH1F*) gROOT->FindObject(hName);
    if(hP1) delete hP1;
    sprintf(dCommand,"p1>>%s",hName);
    tree->Draw(dCommand,Cut,"goff");
    hP1 = (TH1F*) gROOT->FindObject(hName);
    cout << Form("   - p1. ") << endl;        
    p1Min = hP1->GetXaxis()->GetXmin();
    p1Max = hP1->GetXaxis()->GetXmax();

    sprintf(hName,"hX1");
    hX1 = (TH1F*) gROOT->FindObject(hName);
    if(hX1) delete hX1;
    hX1 = new TH1F(hName,"",xNbin,x1Min,x1Max);
  } else if(opt.Contains("auto")) {

    sprintf(hName,"hX1");
    sprintf(dCommand,"x1>>%s",hName);
    tree->Draw(dCommand,Cut,"goff");
    hX1 = (TH1F*) gROOT->FindObject(hName);
    cout << Form("   - x1. ") << endl;    
    x1Min = hX1->GetXaxis()->GetXmin();
    x1Max = hX1->GetXaxis()->GetXmax();

    sprintf(hName,"hP1");
    sprintf(dCommand,"p1>>%s",hName);
    tree->Draw(dCommand,Cut,"goff");
    hP1 = (TH1F*) gROOT->FindObject(hName);
    cout << Form("   - p1. ") << endl;        
    p1Min = hP1->GetXaxis()->GetXmin();
    p1Max = hP1->GetXaxis()->GetXmax();    
  } else {
    sprintf(hName,"hX1");
    hX1 = (TH1F*) gROOT->FindObject(hName);
    if(hX1) delete hX1;
    hX1 = new TH1F(hName,"",xNbin,x1Min,x1Max);

    sprintf(hName,"hP1");
    hP1 = (TH1F*) gROOT->FindObject(hName);
    if(hP1) delete hP1;
    hP1 = new TH1F(hName,"",pNbin,p1Min,p1Max);
  }

  sprintf(hName,"hP1X1");
  TH2F *hP1X1 = (TH2F*) gROOT->FindObject(hName);
  if(hP1X1) delete hP1X1;
  hP1X1 = new TH2F(hName,"",xNbin,x1Min,x1Max,pNbin,p1Min,p1Max);

  sprintf(hName,"hP2X2");
  TH2F *hP2X2 =  (TH2F*) gROOT->FindObject(hName);
  if(hP2X2) delete hP2X2;
  hP2X2 = new TH2F(hName,"",xNbin,x2Min,x2Max,pNbin,p2Min,p2Max);

  // Sliced quantities:
  // --------------------------------------------------------------------------

  // cout << Form("\n5. Slicing ") << endl;
  // Binning
  
  cout << Form("\n2. Setting the binning: ") << endl; 

  // Analize hX1 for the binning:
  
  Int_t SNbin = 40;
  Float_t x1BinMin = 12.55;
  Float_t x1BinMax = 13.35;

  if(sim.Contains("DRI6")) {
    x1BinMin = 12.0;
    x1BinMax = 13.1;
  } else if(sim.Contains("DRI7")) {
    x1BinMin = 12.55;
    x1BinMax = 13.35;
  }
  
  // Set the binning
  Float_t *sBinLim = new Float_t[SNbin+1];
  sBinLim[0] = x1BinMin;
  sBinLim[SNbin] = x1BinMax;
  Float_t slbinSize = (sBinLim[SNbin] - sBinLim[0])/SNbin;
  for(Int_t i=1;i<SNbin;i++) {
    sBinLim[i] = sBinLim[i-1] + slbinSize;
  }
  
  TH1F **hP1sl = new TH1F*[SNbin];
  TH2F **hP2X2sl = new TH2F*[SNbin];
  for(Int_t k=0;k<SNbin;k++) {    
    sprintf(hName,"hP2X2sl_%2i",k);
    hP2X2sl[k] = (TH2F*) gROOT->FindObject(hName);
    if(hP2X2sl[k]) delete hP2X2sl[k];
    hP2X2sl[k] = new TH2F(hName,"",xNbin,x2Min,x2Max,pNbin,p2Min,p2Max);

    sprintf(hName,"hP1sl_%2i",k);
    hP1sl[k] = (TH1F*) gROOT->FindObject(hName);
    if(hP1sl[k]) delete hP1sl[k];
    hP1sl[k] = new TH1F(hName,"",pNbin,p1Min,p1Max);

  }


  // Main loop!
  cout << Form("\n3. Filling histograms: ") << endl; 

  hX1->Reset();
  hX1->SetBins(xNbin,x1Min,x1Max);
  hP1->Reset();
  hP1->SetBins(pNbin,p1Min,p1Max);

  Int_t nentries = (Int_t)tree->GetEntries();  
  for(Int_t i=0;i<nentries;i++) {
    tree->GetEntry(i);

    if(var[0]<x1Min || var[0]>x1Max ) continue; 
    if(var[1]<x2Min || var[1]>x2Max ) continue; 
    if(var[2]<x3Min || var[2]>x3Max ) continue; 
       
    hX1->Fill(var[0],TMath::Abs(var[6]));
    hP1->Fill(var[3],TMath::Abs(var[6]));
 
    hP1X1->Fill(var[0],var[3],TMath::Abs(var[6]));
    hP2X2->Fill(var[1],var[4],TMath::Abs(var[6]));
    
    // Slices    
    if(var[0]<sBinLim[0] || var[0]>sBinLim[SNbin]) continue;
    Int_t iBin = -1;
    for(Int_t j=0; j<SNbin; j++) {
      
      if(var[0]<sBinLim[j+1]) {
	iBin = j;
	break;
      }
    }
    if(iBin<0) continue;

    //    cout << Form(" entry:  %i , bin:  %i",i,iBin)<< endl;

    hP1sl[iBin]->Fill(var[3],TMath::Abs(var[6]));   
    hP2X2sl[iBin]->Fill(var[1],var[4],TMath::Abs(var[6]));
  
  
  }
  

  // Integrated long. emittance:

  cout << Form("\n4. Calculating integrated quantities.. ") << endl;

  // Longitudinal phasespace
  Double_t xmean = 0.0;
  Double_t ymean = 0.0;
  Double_t x2mean = 0.0;
  Double_t y2mean = 0.0;
  Double_t xymean = 0.0;
  Double_t Ntotal = 0.0;
  for(Int_t i=1;i<=xNbin;i++) {
    Double_t x = hP1X1->GetXaxis()->GetBinCenter(i);
    // if(x<xmin || x>xmax) continue;
    for(Int_t j=1;j<=pNbin;j++) {
      Double_t y = hP1X1->GetYaxis()->GetBinCenter(j);
      // if(y<ymin || y>ymax) continue;
      Double_t value = TMath::Abs(hP1X1->GetBinContent(i,j));
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

  cout << Form("  xMean = %7.3f   yMean = %7.3f",xmean,ymean) << endl;
  cout << Form("  xRms  = %7.3f   yRms  = %7.3f",xrms,yrms) << endl;
  cout << Form("  Emittance = %7.3f",emittance) << endl;

  Double_t emitz = emittance;
  Double_t zmean = xmean;
  Double_t zrms = xrms;
  Double_t pzmean = ymean;
  Double_t pzrms = yrms;

  // Transverse phasespace
  xmean = 0.0;
  ymean = 0.0;
  x2mean = 0.0;
  y2mean = 0.0;
  xymean = 0.0;
  Ntotal = 0.0;
  for(Int_t i=1;i<=xNbin;i++) {
    Double_t x = hP2X2->GetXaxis()->GetBinCenter(i);
    // if(x<xmin || x>xmax) continue;
    for(Int_t j=1;j<=pNbin;j++) {
      Double_t y = hP2X2->GetYaxis()->GetBinCenter(j);
      // if(y<ymin || y>ymax) continue;
      Double_t value = TMath::Abs(hP2X2->GetBinContent(i,j));
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

  xrms2  = x2mean - xmean*xmean;
  yrms2  = y2mean - ymean*ymean;
  xrms   = TMath::Sqrt(xrms2);
  yrms   = TMath::Sqrt(yrms2);
  xyrms2 = xymean - xmean*ymean;

  emittance = TMath::Sqrt(xrms2*yrms2 - xyrms2*xyrms2);

  cout << Form("  xMean = %7.3f   yMean = %7.3f",xmean,ymean) << endl;
  cout << Form("  xRms  = %7.3f   yRms  = %7.3f",xrms,yrms) << endl;
  cout << Form("  Emittance = %7.3f",emittance) << endl;

  Double_t emity = emittance;
  Double_t y_mean = xmean;
  Double_t y_rms = xrms;
  Double_t pymean = ymean;
  Double_t pyrms = yrms;

  
  // Charge  
  hX1->Scale(dx1*dx2*dx3);
  Double_t Charge = hX1->Integral();
  
  // Charge *= dx1*dx2*dx3;
 
  if(opt.Contains("units")) {
    Double_t dV = skindepth * skindepth * skindepth;
    Charge *= n0 * dV * (PConst::ElectronCharge/PUnits::picocoulomb);
    cout << Form(" Integrated charge (RAW) of specie %3i = %8f pC",index,Charge) << endl;
  } else {
    cout << Form(" Integrated charge (RAW) of specie %3i = %8.4f n0 * kp^-3",index,Charge) << endl;
  }
  

  cout << Form("\n6. Calculating sliced quantities.. ") << endl;

  TGraph *gemit = NULL;
  TGraph *gYrms = NULL;
  TGraph *gErms = NULL;
  TGraph *gErmsB = NULL;
 
  Double_t * sxmean = new Double_t[SNbin];
  Double_t * symean = new Double_t[SNbin];
  Double_t * sx2mean = new Double_t[SNbin];
  Double_t * sy2mean = new Double_t[SNbin];
  Double_t * sxymean = new Double_t[SNbin];
  Double_t * sNtotal = new Double_t[SNbin];
  Double_t * sxrms2 = new Double_t[SNbin];  
  Double_t * syrms2 = new Double_t[SNbin]; 
  Double_t * sxrms = new Double_t[SNbin];  
  Double_t * syrms = new Double_t[SNbin];  
  Double_t * sxyrms2 = new Double_t[SNbin];
 
  Double_t * xbin = new Double_t[SNbin];
  Double_t * semittance = new Double_t[SNbin];

  Double_t * sNEtotal = new Double_t[SNbin]; 
  Double_t * sEmean = new Double_t[SNbin];
  Double_t * sE2mean = new Double_t[SNbin];
  Double_t * sErms = new Double_t[SNbin];

  for(Int_t k=0;k<SNbin;k++) {
    sxmean[k] = symean[k] = sx2mean[k] = sy2mean[k] = sxymean[k] 
      = sNtotal[k] = sxrms2[k] = syrms2[k] = sxrms[k] = syrms[k]
      = sxyrms2[k] = xbin[k] = semittance[k] = 0.0;
    sNEtotal[k] = sEmean[k] = sE2mean[k] = sErms[k] = 0.0;
    
    xbin[k] = (sBinLim[k] + sBinLim[k+1])/2.;
    
    for(Int_t i=1;i<=xNbin;i++) {
      Double_t x = hP2X2sl[k]->GetXaxis()->GetBinCenter(i);
      // if(x<xmin || x>xmax) continue;
      for(Int_t j=1;j<=pNbin;j++) {
	Double_t y = hP2X2sl[k]->GetYaxis()->GetBinCenter(j);
	// if(y<ymin || y>ymax) continue;
	Double_t value = TMath::Abs(hP2X2sl[k]->GetBinContent(i,j));
	sxmean[k] += x*value;
	symean[k] += y*value;
	sx2mean[k] += x*x*value;
	sy2mean[k] += y*y*value;
	sxymean[k] += x*y*value;
	
	sNtotal[k] += value;
      }	
    }
    
    for(Int_t i=1;i<=pNbin;i++) {
      Double_t y = hP1sl[k]->GetXaxis()->GetBinCenter(i);
      Double_t value = TMath::Abs(hP1sl[k]->GetBinContent(i));
      sEmean[k] += y*value;
      sE2mean[k] += y*y*value;
      sNEtotal[k] += value;
    }
    
    sxmean[k]  /= sNtotal[k];
    symean[k]  /= sNtotal[k];
    sx2mean[k] /= sNtotal[k];
    sy2mean[k] /= sNtotal[k];
    sxymean[k] /= sNtotal[k];
      
    sxrms2[k]  = sx2mean[k] - sxmean[k]*sxmean[k];
    syrms2[k]  = sy2mean[k] - symean[k]*symean[k];
    sxrms[k]   = TMath::Sqrt(sxrms2[k]);
    syrms[k]   = TMath::Sqrt(syrms2[k]);
    sxyrms2[k] = sxymean[k] - sxmean[k]*symean[k];
      
    semittance[k] = TMath::Sqrt(sxrms2[k]*syrms2[k] - sxyrms2[k]*sxyrms2[k]);

    sEmean[k]  /= sNEtotal[k];
    sE2mean[k] /= sNEtotal[k];
    sErms[k]   =  TMath::Sqrt(sE2mean[k] - sEmean[k]*sEmean[k]);
    
    

    cout<< Form("\nk = %i : (x1 > %f && x1 < %f)",k,sBinLim[k],sBinLim[k+1]) << endl; 

    cout << Form("  xMean = %7.3f   yMean = %7.3f",sxmean[k],symean[k]) << endl;
    cout << Form("  xRms  = %7.3f   yRms  = %7.3f",sxrms[k],syrms[k]) << endl;
    cout << Form("  Emittance = %7.3f",semittance[k]) << endl;

    cout << Form("  Emean = %7.3f   Erms = %7.3f",sEmean[k],sErms[k]) << endl;
    

  }


  

  // Chaning to user units: 
  // --------------------------
  
  if(opt.Contains("units") && n0) {
    
    Int_t NbinsX = hP1X1->GetNbinsX();
    Double_t xMin = skindepth * hP1X1->GetXaxis()->GetXmin() / PUnits::um;
    Double_t xMax = skindepth * hP1X1->GetXaxis()->GetXmax() / PUnits::um;
    Int_t NbinsY = hP1X1->GetNbinsY();
    Double_t yMin = hP1X1->GetYaxis()->GetXmin() * BeamMass / PUnits::GeV;
    Double_t yMax = hP1X1->GetYaxis()->GetXmax() * BeamMass / PUnits::GeV;
    hP1X1->SetBins(NbinsX,xMin,xMax,NbinsY,yMin,yMax);
    // Converting electron density
    Double_t dVb = skindepth * skindepth * skindepth;
    Double_t dX = (xMax-xMin)/NbinsX; 
    Double_t dE = (yMax-yMin)/NbinsY; 
    for(Int_t j=0;j<hP1X1->GetNbinsX();j++) {
      for(Int_t k=0;k<hP1X1->GetNbinsY();k++) {
	Double_t binValue =  fabs(hP1X1->GetBinContent(j,k) * dx1 * dx2 * dx3 * dVb * n0 *
				  (PConst::ElectronCharge/PUnits::picocoulomb));
     	//cout << Form(" value = %f",binValue) << endl;
	hP1X1->SetBinContent(j,k,binValue);
	
      }
    }
    
    if(opt.Contains("comov"))
      hP1X1->GetXaxis()->SetTitle("#zeta [#mum]");
    else
      hP1X1->GetXaxis()->SetTitle("z [#mum]");
    
    hP1X1->GetYaxis()->SetTitle("p_{z} [GeV/c]");
    
    hP1X1->GetZaxis()->SetTitle("dQ/d#zetadp_{z} [pC]");
    hP1X1->GetZaxis()->CenterTitle();

    hP1->SetBins(NbinsY,yMin,yMax);
    hP1->GetYaxis()->SetTitle("p_{z} [GeV/c]");

    hX1->SetBins(NbinsX,xMin,xMax);
    Double_t binSize = (xMax - xMin)/NbinsX;

    Double_t dV = skindepth * skindepth * skindepth;
    Double_t  lightspeed =  PConst::c_light / (PUnits::um/PUnits::femtosecond);
    // cout << Form("Speed of light = %f",lightspeed) << endl;
    hX1->Scale(TMath::Abs(n0 * dV * (PConst::ElectronCharge/PUnits::picocoulomb) * (lightspeed/binSize)));
    
    // hX1->Scale(TMath::Abs((PUnits::um/skindepth)*(PConst::ElectronCharge/PUnits::picocoulomb)*PConst::c_light));
    
    // hX1->GetYaxis()->SetTitle("I[kA]");
    hX1->GetYaxis()->SetTitle("");
    if(opt.Contains("comov"))
      hX1->GetXaxis()->SetTitle("#zeta [#mum]");
    else
      hX1->GetXaxis()->SetTitle("z [#mum]");
    
    x1Min *= skindepth / PUnits::um;
    x1Max *= skindepth / PUnits::um;
    zmean *= skindepth / PUnits::um;
    zrms  *= skindepth / PUnits::um;

    p1Min *= BeamMass / PUnits::GeV;
    p1Max *= BeamMass / PUnits::GeV;
    pzmean *= BeamMass / PUnits::GeV;
    pzrms  *= BeamMass / PUnits::GeV;
    
    emitz *= (skindepth / PUnits::um);
    emity *= (skindepth / PUnits::um);
    
    for(Int_t k=0;k<SNbin;k++) {
      xbin[k] *= skindepth / PUnits::um;

      sxmean[k] *= skindepth / PUnits::um;
      sxrms[k]  *= skindepth / PUnits::um;
      symean[k] *= BeamMass / PUnits::MeV;
      syrms[k] *= BeamMass / PUnits::MeV;
      
      semittance[k] *= (skindepth / PUnits::um);

      sEmean[k] *= BeamMass / PUnits::GeV;
      sErms[k]  *= 100 * BeamMass / PUnits::GeV / pzmean; //sEmean[k];
      // sErms[k]  *= BeamMass / PUnits::GeV;

    }

  }

  // Centering in the x1 mean
  if(opt.Contains("zmean")) {
    hX1->SetBins(xNbin,x1Min-zmean,x1Max-zmean);
    hP1X1->SetBins(xNbin,x1Min-zmean,x1Max-zmean,pNbin,p1Min,p1Max);
    for(Int_t k=0;k<SNbin;k++) {
      xbin[k] -= zmean;
    }
    zmean = 0.0;
  }
  // ------

  // Create the graph with the sliced quantities:
  gemit = new TGraph(SNbin,xbin,semittance);
  gYrms = new TGraph(SNbin,xbin,sxrms);
  gErms = new TGraph(SNbin,xbin,sErms);
  
  
  // Profile energy for p1 vs x1:
  TString pname = hP1X1->GetName();
  pname += "_pfx";
  TProfile *hP1X1prof = (TProfile*) gROOT->FindObject(pname.Data());
  if(hP1X1prof) { delete hP1X1prof; hP1X1prof = NULL; }
  hP1X1prof = hP1X1->ProfileX("_pfx",1,-1,"s");

  // get the errors from the profile:
  Int_t NP1X1Bins = hP1X1prof->GetNbinsX();
  Double_t *x1bins = new Double_t[NP1X1Bins];
  Double_t *eRms   = new Double_t[NP1X1Bins];
  for(Int_t i=1;i<=hP1X1prof->GetNbinsX();i++) {
    x1bins[i] = hP1X1prof->GetBinCenter(i);
    eRms[i] = 100 * hP1X1prof->GetBinError(i) / hP1X1prof->GetBinContent(i);
  }
  gErmsB = new TGraph(NP1X1Bins,x1bins,eRms);
  
  // Vertical Energy histogram:
  // --------------------------------------------------------------------------------   
  TGraph *gP1left = NULL;
  if(hP1) {
    Double_t *yarray   = new Double_t[pNbin];
    Double_t *xarray   = new Double_t[pNbin];
    
    // This is for the right side:
    // Double_t xMax = x1Min + (x1Max-x1Min) * 0.9;
    // Double_t xMin = x1Max;
    // And this for left:
    Double_t xMin = hX1->GetXaxis()->GetXmin();
    Double_t xMax = hX1->GetXaxis()->GetXmin() + (hX1->GetXaxis()->GetXmax()
						  -hX1->GetXaxis()->GetXmin()) * 0.2;
    Double_t EneMax = hP1->GetMaximum();
    // cout << Form("  EneMax = %f ", EneMax) << endl;
 
    for(Int_t j=0; j<pNbin; j++) {
      yarray[j] = hP1->GetBinCenter(j+1);
      xarray[j] = ((xMax-xMin)/EneMax)*hP1->GetBinContent(j+1) + xMin;

      // cout << Form("  x = %f  y = %f ", xarray[j],yarray[j]) << endl;
    }

    gP1left = new TGraph(pNbin,xarray,yarray);
    gP1left->SetLineColor(PlasmaGlob::elecLine);
    gP1left->SetLineWidth(2);
    gP1left->SetFillStyle(1001);
    gP1left->SetFillColor(PlasmaGlob::elecFill);
       
  }

  
  // Ranges!!
  Double_t yMin =  999.9;
  Double_t yMax =  -999.9;
  for(Int_t k=0;k<SNbin;k++) {
    if(semittance[k]<yMin)
      yMin = semittance[k];
    
    if(semittance[k]>yMax)
      yMax = semittance[k];

    if(sErms[k]<yMin)
      yMin = sErms[k];
    
    if(sErms[k]>yMax)
      yMax = sErms[k];
  }

  for(Int_t k=1;k<=xNbin;k++) {
    Double_t value = hX1->GetBinContent(k);
    if(value<yMin)
      yMin = value;
    
    if(value>yMax)
      yMax = value;

  }


  // Plotting
  // -----------------------------------------------
    
  cout << "\n7. Plotting... " << endl;

  // Canvas setup
  // Create the canvas and the pads before the Frame loop
  // Resolution:
  Int_t sizex = 800;
  Int_t sizey = 600;
  if(opt.Contains("hres")) {
    Int_t sizex = 1600;
    Int_t sizey = 1200;    
  }
  
  TCanvas *C = new TCanvas("C1","Evolution of Injection",sizex,sizey);
  C->cd();

  // Set palette:
  PPalette * pPalette = (PPalette*) gROOT->FindObject("electron");
  pPalette->cd();

  // Float_t Max  = hP1X1->GetMaximum();
  // Float_t Min  = hP1X1->GetMinimum();
  
  // hP1X1->GetZaxis()->SetRangeUser(Min,Max); 

  // Text objects
  TPaveText *textTime = new TPaveText(0.55,0.76,0.82,0.86,"NDC");
  PlasmaGlob::SetPaveTextStyle(textTime,32); 
  textTime->SetTextColor(kGray+2);
  char ctext[128];
  if(opt.Contains("units") && n0) 
    sprintf(ctext,"z = %5.1f mm", Time * skindepth / PUnits::mm);
  else
    sprintf(ctext,"t = %5.1f #omega_{p}^{-1}",Time);
  textTime->AddText(ctext);
 
  TPaveText *textCharge = new TPaveText(0.15,0.25,0.48,0.3,"NDC");
  PlasmaGlob::SetPaveTextStyle(textCharge,12); 
  textCharge->SetTextColor(kGray+2);
  if(opt.Contains("units") && n0)
    sprintf(ctext,"Q = %5.1f pC", Charge);
  else
    sprintf(ctext,"Q = %5.1f n0#timeskp^{-3}", Charge);    
  textCharge->AddText(ctext);

  TPaveText *textMom = new TPaveText(0.55,0.03,0.82,0.13,"NDC");
  PlasmaGlob::SetPaveTextStyle(textMom,32); 
  textMom->SetTextColor(kGray+3);
  textMom->SetTextFont(62);
  if(opt.Contains("units") && n0)
    sprintf(ctext,"#LTp_{z}#GT = %5.1f GeV/c", pzmean);
  else
    sprintf(ctext,"#LTp_{z}#GT = %5.1f mc", pzmean);    
  textMom->AddText(ctext);


  TPaveText *textInfo = new TPaveText(0.55,0.40,0.82,0.75,"NDC");
  PlasmaGlob::SetPaveTextStyle(textInfo,32); 
  textInfo->SetTextColor(kGray+2);
  textInfo->SetTextFont(42);
  sprintf(ctext,"Q = %5.1f pC",Charge);
  textInfo->AddText(ctext);
  //  sprintf(ctext,"#LT#zeta#GT_{rms} = %5.1f #mum",zrms);
  sprintf(ctext,"#sigma_{#zeta} = %5.1f #mum",zrms);
  textInfo->AddText(ctext);
  //  sprintf(ctext,"#LTp_{z}#GT_{rms} = %5.1f GeV/c",pzrms);
  sprintf(ctext,"#sigma_{p_{z}} = %5.1f GeV/c",pzrms);
  textInfo->AddText(ctext);
  sprintf(ctext,"#epsilon_{y} = %5.2f nm",emity*1000);
  textInfo->AddText(ctext);
  
  // Setup Pad layout:
  const Int_t NFrames = 2;
  TPad *pad[NFrames];
  TH1F *hFrame[NFrames];
  TString sLabels[] = {"(a)","(b)"};
  TPaveText **textLabel = new TPaveText*[NFrames];
  
  Float_t lMargin = 0.15;
  Float_t rMargin = 0.18;
  Float_t bMargin = 0.15;
  Float_t tMargin = 0.04;
  Float_t factor = 1.0;    
  PlasmaGlob::CanvasAsymPartition(C,NFrames,lMargin,rMargin,bMargin,tMargin,factor);

  for(Int_t k=0;k<NFrames;k++) {

    char padname[16];
    sprintf(padname,"pad_%i",k);
    pad[k] = (TPad*) gROOT->FindObject(padname);

    Float_t vfactor = pad[0]->GetAbsHNDC()/pad[k]->GetAbsHNDC();

    char name[16];
    
    sprintf(name,"hFrame_%i",k);  
    hFrame[k] = (TH1F*) gROOT->FindObject(name);
    if(hFrame[k]) delete hFrame[k];
    hFrame[k] = (TH1F*) hX1->Clone(name);
    hFrame[k]->Reset();
    
    hFrame[k]->GetXaxis()->CenterTitle();
    hFrame[k]->GetYaxis()->CenterTitle();
    hFrame[k]->GetZaxis()->CenterTitle();
    hFrame[k]->SetLabelFont(42,"xyz");
    hFrame[k]->SetTitleFont(42,"xyz");
    
    hFrame[k]->SetNdivisions(505,"xyz");
    
    hFrame[k]->SetTickLength(0.04*vfactor,"x");
    hFrame[k]->SetTickLength(0.02/vfactor,"y");
    hFrame[k]->SetTickLength(0.02/vfactor,"z");
    
    hFrame[k]->GetYaxis()->SetLabelSize(0.08*vfactor);
    hFrame[k]->GetYaxis()->SetLabelOffset(0.02/vfactor);
    
    hFrame[k]->GetYaxis()->SetTitleSize(0.10*vfactor);
    hFrame[k]->GetYaxis()->SetTitleOffset(0.7/vfactor);

    if(k==0) {  
      hFrame[k]->GetXaxis()->SetLabelSize(0.08);
      hFrame[k]->GetXaxis()->SetLabelOffset(0.02);
      hFrame[k]->GetXaxis()->SetTitleSize(0.10);
      hFrame[k]->GetXaxis()->SetTitleOffset(1.1);
    } else {
      hFrame[k]->GetXaxis()->SetLabelSize(0.0);
      hFrame[k]->GetXaxis()->SetTitleSize(0.0);
    }
  }

  
  cout <<  "  - Top plot" << endl;
  C->cd(0);
  char padname[16];
  sprintf(padname,"pad_%i",1);
  pad[1] = (TPad*) gROOT->FindObject(padname);
  pad[1]->Draw();
  pad[1]->cd(); // <---------------------------------------------- Top Plot ---------
  if(opt.Contains("logz")) {
    pad[1]->SetLogz(1);
  } else {
    pad[1]->SetLogz(0);
  }
  pad[1]->SetFrameLineWidth(1);  
  pad[1]->SetTickx(1);

  hFrame[1]->GetYaxis()->SetRangeUser(hP1X1->GetYaxis()->GetXmin(),hP1X1->GetYaxis()->GetXmax());
  
  if(opt.Contains("units"))
    hFrame[1]->GetYaxis()->SetTitle("p_{z} [GeV/c]");
  
  hFrame[1]->Draw();
  
  gP1left->SetLineWidth(2);
  gP1left->Draw("F");
  gP1left->Draw("L");

  TLine lZmean(zmean,hP1X1->GetYaxis()->GetXmin(),zmean,hP1X1->GetYaxis()->GetXmax());
  lZmean.SetLineColor(kGray+2);
  lZmean.SetLineStyle(2);
  lZmean.Draw();

  TLine lPmean(hP1X1->GetXaxis()->GetXmin(),pzmean,hP1X1->GetXaxis()->GetXmax(),pzmean);
  lPmean.SetLineColor(kGray+2);
  lPmean.SetLineStyle(2);
  lPmean.Draw();

  // 2D histogram z range
  Float_t dmax = 1E-1*TMath::Nint(hP1X1->GetMaximum()*1E1);
  Float_t dmin = 0.0;
  hP1X1->GetZaxis()->SetRangeUser(dmin,dmax);

  hP1X1->GetYaxis()->SetNdivisions(503);
  hP1X1->GetZaxis()->SetNdivisions(503);
  hP1X1->GetZaxis()->SetLabelSize(0.05);
  hP1X1->GetZaxis()->SetTitleSize(0.04);
  hP1X1->GetZaxis()->SetTitleFont(42);

  hP1X1->Draw("colzsame");
  // hP1X1->SetContour(20);
  // hP1X1->Draw("contzsame");
  // hP1X1prof->SetMarkerStyle(1);
  // hP1X1prof->SetLineWidth(2);
  // hP1X1prof->Draw("zsame");

  //hP1->Draw("C");
  
  gPad->Update();

  TPaletteAxis *palette = (TPaletteAxis*)hP1X1->GetListOfFunctions()->FindObject("palette");
  if(palette) {
    Float_t y1 = gPad->GetBottomMargin();
    Float_t y2 = 1 - gPad->GetTopMargin();
    Float_t x1 = 1 - gPad->GetRightMargin();
    palette->SetY2NDC(y2 - 0.04);
    palette->SetY1NDC(y1 + 0.04);
    palette->SetX1NDC(x1 + 0.01);
    palette->SetX2NDC(x1 + 0.04);
    
    palette->SetLabelFont(42);
    
    palette->SetLabelSize(0.08);
    //palette->SetLabelOffset(0.005/vfactor);
    palette->SetTitleSize(0.10);
    //    palette->SetTitleOffset(9999.0*vfactor);
    palette->SetTitleOffset(0.6);
    
    palette->SetBorderSize(2);
    palette->SetLineColor(1);
  }
  // cout << " DRAW !!" << endl;


  textTime->Draw();
  textInfo->Draw();
  // textCharge->Draw();
  textMom->Draw();
   
  gPad->RedrawAxis(); 

  // Bottom plot -----------------------------------------
  cout << "  - Bottom plot" << endl;
  C->cd(0);
  sprintf(padname,"pad_%i",0);
  pad[0] = (TPad*) gROOT->FindObject(padname);
  pad[0]->Draw();
  pad[0]->cd(); // <---------------------------------------------------------- Bottom Plot
  // if(opt.Contains("logz")) {
  //   pad[0]->SetLogz(1);
  // } else {
  //   pad[0]->SetLogz(0);
  // }
  pad[0]->SetFrameLineWidth(1);  
  pad[0]->SetTickx(1);

  hFrame[0]->GetYaxis()->SetRangeUser(0.0,1.1*yMax);
  hFrame[0]->Draw();

  hX1->GetYaxis()->SetNdivisions(503);
  hX1->SetLineWidth(2);
  hX1->SetFillStyle(1001);
  hX1->SetFillColor(PlasmaGlob::elecFill);
  // hX1->SetLineColor(kBlue);
  
  //hX1->Smooth();
  hX1->Draw("FL same");
  //hX1->Draw("C");

  TLine lZmean2(zmean,0.0,zmean,1.1*yMax);
  lZmean2.SetLineColor(kGray+2);
  lZmean2.SetLineStyle(2);
  lZmean2.Draw();

  Int_t markerSize = 1.2; 
  Int_t lineWidth  = 2.0;   

  gYrms->SetMarkerStyle(20);
  gYrms->SetLineStyle(1);
  gYrms->SetMarkerColor(kGray+1);
  gYrms->SetMarkerSize(markerSize); 
  gYrms->SetLineColor(kGray+1);
  gYrms->SetLineWidth(lineWidth);
  // gYrms->Draw("PL");
  
  // hP2X2sl[0]->Draw("colz");
  gemit->SetMarkerStyle(20);
  //  gemit->SetMarkerColor(kMagenta-2);
  gemit->SetMarkerColor(kGray+2);
  gemit->SetMarkerSize(markerSize);
  gemit->SetLineWidth(lineWidth);
  gemit->SetLineColor(kGray+2);
  gemit->Draw("PL");

  gErms->SetMarkerStyle(20);
  gErms->SetMarkerSize(markerSize);
  gErms->SetMarkerColor(kOrange+10);
  gErms->SetLineColor(kOrange+10);
  gErms->SetLineWidth(lineWidth);
  gErms->Draw("PL");


  TLegend *Leg;
  Leg=new TLegend(0.55,0.75,1 - gPad->GetRightMargin() - 0.02,0.95);
    
  
  PlasmaGlob::SetPaveStyle(Leg);
  Leg->SetTextAlign(12);
  Leg->SetTextColor(kGray+3);
  Leg->SetTextFont(42);
  Leg->SetLineColor(1);
  Leg->SetBorderSize(0);
  Leg->SetFillColor(0);
  Leg->SetFillStyle(1001);
  Leg->SetFillStyle(0); // Hollow
 
  Leg->AddEntry(hX1  ,"Current [kA]","L");
  //  Leg->AddEntry(gErms,"Energy spread (GeV)","PL");
  Leg->AddEntry(gErms,"Energy spread [%]","PL");
  Leg->AddEntry(gemit,"Emittance [#mum]","PL");
  //  Leg->AddEntry(gYrms,"Bunch width [#mum]","PL");
 
  Leg->Draw();

  gPad->RedrawAxis(); 

  gPad->Update();
  
  // Print to file --------------------------------------
  
  C->cd();
  
  // Print to a file
  // Output file
  TString fOutName = Form("./%s/Plots/RakeBunch/RakeBunch",sim.Data());
  fOutName += Form("-%s_%.3f",sim.Data(),time);

  PlasmaGlob::imgconv(C,fOutName,opt);
  // ---------------------------------------------------------
  
}
