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
#include <TLatex.h>
#include <TEllipse.h>
#include <TMarker.h>
#include <TPaveText.h>
#include <TPaletteAxis.h>
#include <TExec.h>
#include <TClonesArray.h>


#include "PData.hh"
#include "PDataHiP.hh"
#include "PGlobals.hh"
#include "PConst.hh"
#include "PPalette.hh"


UInt_t BitCounter(UInt_t v) {
  UInt_t c;

  for (c = 0; v; c++) {
    v &= v - 1; // clear the least significant bit set
  }

  return c;
  
}

string DecToBin(int number)
{
  if ( number == 0 ) return "0";
  if ( number == 1 ) return "1";

  if ( number % 2 == 0 )
    return DecToBin(number / 2) + "0";
  else
    return DecToBin(number / 2) + "1";
}

int BinToDec(string number)
{
  int result = 0, pow = 1;
  for ( int i = number.length() - 1; i >= 0; --i, pow <<= 1 )
    result += (number[i] - '0') * pow;

  return result;
}


void PlotSnapshot( const TString &sim, Int_t timestep, UInt_t mask = 3, const TString &options="") {
  
#ifdef __CINT__  
  gSystem->Load("libptools.so");
#endif

  PData *pData = PData::Get(sim.Data());
  if(pData->isHiPACE()) {
    delete pData; pData = NULL;
    pData = PDataHiP::Get(sim.Data());
  }
  
  pData->LoadFileNames(timestep);
  if(!pData->IsInit()) return;
  
  string imask = DecToBin(mask);
  cout << Form("\n Plotting Snapshot with mask: %s",imask.c_str()) << endl; 
  
  PGlobals::Initialize();
  
  TString opt = options;
  
  // More makeup
  gStyle->SetPadGridY(0);
  if(opt.Contains("gridx")) {
    gStyle->SetPadGridX(1);
  }
  if(opt.Contains("gridy")) {
    gStyle->SetPadGridY(1);
  }
  gStyle->SetNumberContours(255);
  gStyle->SetJoinLinePS(2);
  
  // Some plasma constants
  Float_t n0 = pData->GetPlasmaDensity();
  Float_t omegap = pData->GetPlasmaFrequency();
  Float_t timedepth = pData->GetPlasmaTimeDepth();
  Float_t kp = pData->GetPlasmaK();
  Float_t skindepth = pData->GetPlasmaSkinDepth();
  Float_t E0 = pData->GetPlasmaE0();

  // Some beam properties:
  // Float_t Ebeam = pData->GetBeamEnergy();
  // Float_t gamma = pData->GetBeamGamma();
  // Float_t vbeam = pData->GetBeamVelocity();

  // Other parameters
  //  Float_t trapPotential = 1.0 - (1.0/gamma);
  Float_t trapPotential = 1.0;  
  
  // ----------------------------------------------------------------------------------
  
  // Open snapshot file and get the histograms
  TString filename;
  filename = Form("./%s/Plots/Snapshots/Snapshot-%s_%i.root",sim.Data(),sim.Data(),timestep);
  if(opt.Contains("zyslc"))
    filename = Form("./%s/Plots/Snapshots/SnapshotZY-%s_%i.root",sim.Data(),sim.Data(),timestep);
  
  TFile  *ifile = (TFile*) gROOT->GetListOfFiles()->FindObject(filename.Data());
  if (!ifile) ifile = new TFile(filename,"READ");

  ifile->cd();

  // Get options
  TTree *infotree = (TTree*) ifile->Get("infotree");
  TString *opttree = 0;
  infotree->SetBranchAddress("options",&opttree);

  Double_t xon = 0.1;
  infotree->SetBranchAddress("xon",&xon);
  Double_t xoff = 0.1;
  infotree->SetBranchAddress("xoff",&xoff);
  Double_t xint = 0.1;
  infotree->SetBranchAddress("xint",&xint);
  Double_t denUnit = 1;
  infotree->SetBranchAddress("denUnit",&denUnit);  
  string   *denSUnit = 0;
  infotree->SetBranchAddress("denSUnit",&denSUnit);  
  Double_t spaUnit = 1;
  infotree->SetBranchAddress("spaUnit",&spaUnit);  
  string   *spaSUnit = 0;
  infotree->SetBranchAddress("spaSUnit",&spaSUnit);  
  Double_t propUnit = 1;
  infotree->SetBranchAddress("propUnit",&propUnit);  
  string   *propSUnit = 0;
  infotree->SetBranchAddress("propSUnit",&propSUnit);  
  Double_t timUnit = 1;
  infotree->SetBranchAddress("timUnit",&timUnit);  
  string   *timSUnit = 0;
  infotree->SetBranchAddress("timSUnit",&timSUnit);  
  Double_t eUnit = 1;
  infotree->SetBranchAddress("eUnit",&eUnit);  
  string   *eSUnit = 0;
  infotree->SetBranchAddress("eSUnit",&eSUnit);  
  Double_t curUnit = 1;
  infotree->SetBranchAddress("curUnit",&curUnit);  
  string   *curSUnit = 0;
  infotree->SetBranchAddress("curSUnit",&curSUnit);  
    
  infotree->GetEntry(0);

  if(propSUnit==0) {
    propUnit = spaUnit;
    propSUnit = spaSUnit;
  }

  cout << Form(" Options from doSnapshot = %s \n", opttree->Data());
  cout << Form(" Options in PlotSnapshot = %s \n", opt.Data());
  opt += *opttree;
  
  // Off-axis
  if(opt.Contains("units")) {
    xoff *= skindepth/spaUnit;
    xon *= skindepth/spaUnit;
    xint *= skindepth/spaUnit;
  }  

  // Time in OU
  Double_t Time = pData->GetRealTime();
  // cout << Form(" Real time = %.2f  ",Time);
  Time += pData->ShiftT(opt);
  // cout << Form(" Shifted time = %.2f  ",Time) << endl;

  Float_t shiftz = pData->Shift(opt);

  // z start of the plasma in normalized units.
  Float_t zStartPlasma = pData->GetPlasmaStart()*kp - shiftz;
  // z start of the neutral in normalized units.
  Float_t zStartNeutral = pData->GetNeutralStart()*kp - shiftz;
  // z end of the neutral in normalized units.
  Float_t zEndNeutral = pData->GetNeutralEnd()*kp - shiftz;

  
  cout << Form("\n Getting the histograms ... ") ; 

  // Skip one of the species
  Int_t noIndex = -1;
    
  char hName[36];

  Int_t Nspecies = pData->NSpecies();
  // Get charge density histos
  TH2F **hDen2D = new TH2F*[Nspecies];
  // Get charge density on-axis
  TH1F **hDen1D = new TH1F*[Nspecies];
  // And electric current (integrated)
  TH1F **hCur1D = new TH1F*[Nspecies];
  // Get Jz histos
  TH2F **hJz2D = new TH2F*[Nspecies];
  TH1F **hJz1D = new TH1F*[Nspecies];
  for(Int_t i=0;i<Nspecies;i++) {
    hCur1D[i] = NULL;
    if(i==noIndex) continue;

    sprintf(hName,"hDen2D_%i",i); 
    hDen2D[i] = (TH2F*) ifile->Get(hName);
    if(hDen2D[i]) {
      if(hDen2D[i]->GetMaximum() < 1E-4) {
	delete hDen2D[i];
	hDen2D[i] = NULL;
      }
    }
    
    sprintf(hName,"hDen1D_%i",i); 
    hDen1D[i] = (TH1F*) ifile->Get(hName);
    if(hDen1D[i]) {
      if(hDen1D[i]->GetMaximum() < 1E-4) {
	delete hDen1D[i];
	hDen1D[i] = NULL;
      }
    }

    sprintf(hName,"hCur1D_%i",i); 
    hCur1D[i] = (TH1F*) ifile->Get(hName);
    if(hCur1D[i]) {
      if(hCur1D[i]->GetMaximum() < 1E-6) {
	delete hCur1D[i];
	hCur1D[i] = NULL;
      }
    }

    sprintf(hName,"hJz2D_%i",i); 
    hJz2D[i] = (TH2F*) ifile->Get(hName);
    if(hJz2D[i]) {
      if(hJz2D[i]->GetMaximum() < 1E-4) {
	delete hJz2D[i];
	hJz2D[i] = NULL;
      }
    }

    sprintf(hName,"hJz1D_%i",i); 
    hJz1D[i] = (TH1F*) ifile->Get(hName);
    if(hJz1D[i]) {
      if(hJz1D[i]->GetMaximum() < 1E-4) {
	delete hJz1D[i];
	hJz1D[i] = NULL;
      }
    }

  }
  
  
  // Get electric fields 2D
  const Int_t Nfields = 3;
  TH2F **hE2D = new TH2F*[Nfields];
  TH2F **hB2D = new TH2F*[Nfields];
  TH1F **hE1D = new TH1F*[Nfields];
  TH1F **hB1D = new TH1F*[Nfields];
  for(Int_t i=0;i<Nfields;i++) {
    hE2D[i] = hB2D[i] = NULL;
    hE1D[i] = hB1D[i] = NULL;

    sprintf(hName,"hE2D_%i",i); 
    hE2D[i] = (TH2F*) ifile->Get(hName);

    sprintf(hName,"hB2D_%i",i); 
    hB2D[i] = (TH2F*) ifile->Get(hName);

    sprintf(hName,"hE1D_%i",i); 
    hE1D[i] = (TH1F*) ifile->Get(hName);

    sprintf(hName,"hB1D_%i",i); 
    hB1D[i] = (TH1F*) ifile->Get(hName);
    
  }

  TH2F *hFocus2D = NULL;
  TH1F *hFocus1D = NULL;
  if(opt.Contains("zyslc")) {
    sprintf(hName,"hW2D_1"); 
    hFocus2D = (TH2F*) ifile->Get(hName);
    sprintf(hName,"hW1D_1"); 
    hFocus1D = (TH1F*) ifile->Get(hName);
  } else {
    sprintf(hName,"hW2D_0"); 
    hFocus2D = (TH2F*) ifile->Get(hName);
    sprintf(hName,"hW1D_0"); 
    hFocus1D = (TH1F*) ifile->Get(hName);
  }
  sprintf(hName,"hETotal2D");   
  TH2F *hETotal2D = (TH2F*) ifile->Get(hName);
  sprintf(hName,"hETotal1D");   
  TH1F *hETotal1D = (TH1F*) ifile->Get(hName);

  sprintf(hName,"hV2D"); 
  TH2F *hV2D = (TH2F*) ifile->Get(hName);
  sprintf(hName,"hV1D"); 
  TH1F *hV1D = (TH1F*) ifile->Get(hName);

  sprintf(hName,"hA2D"); 
  TH2F *hA2D = (TH2F*) ifile->Get(hName);
  sprintf(hName,"hA1D"); 
  TH1D *hA1D = (TH1D*) ifile->Get(hName);
  
  cout << Form(" done. ") << endl;
  
  // Rescaling box limits:
  if(opt.Contains("rscal") && !opt.Contains("units") && hDen2D[0]) {

    Float_t xMin, xMax, yMin, yMax;
    Int_t NbinsX, NbinsY;
    NbinsX = hDen2D[0]->GetNbinsX();
    xMin = hDen2D[0]->GetXaxis()->GetXmin();
    xMax = hDen2D[0]->GetXaxis()->GetXmax();
    NbinsY = hDen2D[0]->GetNbinsY();
    yMin = hDen2D[0]->GetYaxis()->GetXmin();
    yMax = hDen2D[0]->GetYaxis()->GetXmax();

    // Local value of the density at \zeta = 0 and y = a quarter of the yRange.
    Int_t binx = hDen2D[0]->GetXaxis()->FindBin(0.);
    Int_t biny = hDen2D[0]->GetYaxis()->FindBin(yMax - ((yMax-yMin)/8.0));
    
    Float_t localden = hDen2D[0]->GetBinContent(binx,biny);
    Float_t kfactor = TMath::Sqrt(localden);

    cout << Form(" - Local density = %.2f",localden) << endl;

    xMin *= kfactor;
    xMax *= kfactor;
    yMin *= kfactor;
    yMax *= kfactor;

    xint *= kfactor;
    
    for(Int_t i=0;i<Nspecies;i++) {
      if(i==noIndex) continue;
      
      if(hDen2D[i] && hDen1D[i]) {
	hDen2D[i]->SetBins(NbinsX,xMin,xMax,NbinsY,yMin,yMax);
	hDen1D[i]->SetBins(NbinsX,xMin,xMax);

	for(Int_t j=0;j<=hDen2D[i]->GetNbinsX();j++) {
	  for(Int_t k=0;k<=hDen2D[i]->GetNbinsY();k++) {
	    hDen2D[i]->SetBinContent(j,k, hDen2D[i]->GetBinContent(j,k) / localden );
	  }
	  hDen1D[i]->SetBinContent(j, hDen1D[i]->GetBinContent(j) / localden);
	}

	if(i==0)
	  hDen2D[i]->GetZaxis()->SetTitle("n/n_{0}");
	else if(i==1)
	  hDen2D[i]->GetZaxis()->SetTitle("n_{b}/n_{0}");

      }

      if(hCur1D[i]) {
	hCur1D[i]->SetBins(NbinsX,xMin,xMax);

	hCur1D[i]->Scale(PConst::I0/PUnits::kA);
	hCur1D[i]->GetYaxis()->SetTitle("I_{b}[kA]");
      }
      
      if(hJz2D[i] && hJz1D[i]) {
	hJz2D[i]->SetBins(NbinsX,xMin,xMax,NbinsY,yMin,yMax);
	hJz1D[i]->SetBins(NbinsX,xMin,xMax);

	for(Int_t j=0;j<=hJz2D[i]->GetNbinsX();j++) {
	  for(Int_t k=0;k<=hJz2D[i]->GetNbinsY();k++) {
	    hJz2D[i]->SetBinContent(j,k, hJz2D[i]->GetBinContent(j,k) / localden );
	  }
	  hJz1D[i]->SetBinContent(j, hJz1D[i]->GetBinContent(j) / localden);
	}
      }
      
    }

    for(Int_t i=0;i<Nfields;i++) {
      
      if(hE2D[i] && hE1D[i]) {
	hE2D[i]->SetBins(NbinsX,xMin,xMax,NbinsY,yMin,yMax);
	hE1D[i]->SetBins(NbinsX,xMin,xMax);
	
      	for(Int_t j=0;j<=hE2D[i]->GetNbinsX();j++) {
	  for(Int_t k=0;k<=hE2D[i]->GetNbinsY();k++) {
	    hE2D[i]->SetBinContent(j,k, hE2D[i]->GetBinContent(j,k) / kfactor );
	  }
	  hE1D[i]->SetBinContent(j, hE1D[i]->GetBinContent(j) / kfactor);
	}

	hE2D[i]->GetZaxis()->SetTitle("E_{z}/E_{0}");
      }
	
      if(hB2D[i] && hB1D[i]) {
	hB2D[i]->SetBins(NbinsX,xMin,xMax,NbinsY,yMin,yMax);
	
      	for(Int_t j=0;j<=hB2D[i]->GetNbinsX();j++) {
	  for(Int_t k=0;k<=hB2D[i]->GetNbinsY();k++) {
	    hB2D[i]->SetBinContent(j,k, hB2D[i]->GetBinContent(j,k) / kfactor );
	  }
	  hB1D[i]->SetBinContent(j, hB1D[i]->GetBinContent(j) / kfactor);
	}
      }

      
    }

    if(hFocus2D && hFocus1D) {
      hFocus2D->SetBins(NbinsX,xMin,xMax,NbinsY,yMin,yMax);
      hFocus1D->SetBins(NbinsX,xMin,xMax);

      for(Int_t j=0;j<=hFocus2D->GetNbinsX();j++) {
	for(Int_t k=0;k<=hFocus2D->GetNbinsY();k++) {
	  hFocus2D->SetBinContent(j,k, hFocus2D->GetBinContent(j,k) / kfactor );
	}
	hFocus1D->SetBinContent(j, hFocus1D->GetBinContent(j) / kfactor);
      }
      
    }
    
    if(hV2D)
      hV2D->SetBins(NbinsX,xMin,xMax,NbinsY,yMin,yMax);

    if(hV1D)
      hV1D->SetBins(NbinsX,xMin,xMax);
  }
  
  // Ionization probability rates (ADK)
  // Calculates from the total E the ionization prob. rate for a given species.
  const UInt_t NAtoms = 7;
  char atNames[NAtoms][8] = {"H","He","He2","Ne","Ne2","Ne5","HIT"};
  char atAxNames[NAtoms][8] = {"H","He","He^{+}","Ne","Ne^{+}","Ne^{4+}","HIT"};
  Double_t Eion0[NAtoms] = {13.6 * PUnits::eV, 24.59 * PUnits::eV, 54.4 * PUnits::eV, 21.56 * PUnits::eV, 40.96 * PUnits::eV, 126.247 * PUnits::eV, 85.00 *  PUnits::eV};
  Float_t Z[NAtoms]     = { 1.0,  1.0,  2.0,  1.0, 2.0, 5.0, 2.0};
  // Float_t z10[NAtoms]   = {999.0, 999.0, 999.0};
  // Float_t z100[NAtoms]  = {999.0, 999.0, 999.0};

  TH2F *hIonRate2D[NAtoms] = {NULL};
  TH2F *hIonProb2D[NAtoms] = {NULL};
  TH1F *hIonRate1D[NAtoms] = {NULL};
  TH1F *hIonProb1D[NAtoms] = {NULL};

  // Ion probability: Select species index.
  UInt_t ii = 0; // Hydrogen
  // UInt_t ii = 1; // Helium
  // UInt_t ii = 2; // Helium2
  // UInt_t ii = 5; // Ne4
  // UInt_t ii = 6; // Custom

  if(opt.Contains("Ne")) {
    ii = 3;
  } else if(opt.Contains("He2")) {
    ii = 2;
  } else if(opt.Contains("He")) {
    ii = 1;
  }
  
  if( ((mask & 0x80) || (mask == 0) ) && hETotal2D) { // only if ionization bit is selected
    cout << Form("\n Calculating ionization probability rates (ADK) ... ") ; 
    
    for(UInt_t iat=0;iat<NAtoms;iat++) {
      
      if(iat!=ii) continue;
    
      sprintf(hName,"hIonRate2D_%s",atNames[iat]);   
      hIonRate2D[iat] = (TH2F*) hETotal2D->Clone(hName);
      hIonRate2D[iat]->Reset();
      sprintf(hName,"hIonProb2D_%s",atNames[iat]);
      hIonProb2D[iat] = (TH2F*) hETotal2D->Clone(hName);
      hIonProb2D[iat]->Reset();
    
      Int_t NbinsX = hETotal2D->GetNbinsX();
      Int_t NbinsY = hETotal2D->GetNbinsY();
      Float_t dx = hETotal2D->GetXaxis()->GetBinWidth(1);

      if(opt.Contains("units")) dx *= spaUnit / PConst::c_light;
      else dx *= timedepth;
      
      // if(!opt.Contains("units")) dx *=  timedepth / PUnits::femtosecond; // to fs
      // else dx /= PConst::c_light / (PUnits::um/PUnits::femtosecond); // to fs
    
      for(Int_t j=1;j<=NbinsY;j++) {
      
	Float_t integral[NAtoms] = {0.0};
      
	for(Int_t k=NbinsX;k>0;k--) {

	  Float_t E = hETotal2D->GetBinContent(k,j);
	  if(opt.Contains("units")) E *= eUnit;
	  else E *= E0;
	  if(E/(PUnits::GV/PUnits::m)<10) continue;
		
	  Float_t IonRate = PFunc::ADK_ENG(E,Eion0[iat],Z[iat]);	
	  hIonRate2D[iat]->SetBinContent(k,j,IonRate);
	  integral[iat] += IonRate * dx;	
	  if(integral[iat]>=1) integral[iat] = 1;

	  if(opt.Contains("units"))
	    IonRate /= 1.0/PUnits::femtosecond;
	  else
	    IonRate /= omegap;
	
	  hIonProb2D[iat]->SetBinContent(k,j,(1-integral[iat])*IonRate);
	
	}
      }

      {
	sprintf(hName,"hIonRate1D_%s",atNames[iat]);
	hIonRate1D[iat] = (TH1F*) hETotal1D->Clone(hName);
	hIonRate1D[iat]->Reset();
	sprintf(hName,"hIonProb1D_%s",atNames[iat]);
	hIonProb1D[iat] = (TH1F*) hETotal1D->Clone(hName);
	hIonProb1D[iat]->Reset();

	Int_t NbinsX = hETotal1D->GetNbinsX();
	Float_t dx = hETotal1D->GetBinWidth(1);
      
	if(opt.Contains("units")) dx *= spaUnit / PConst::c_light;
	else dx *= timedepth;

	Float_t integral[NAtoms] = {0.0};
	for(Int_t j=NbinsX;j>0;j--) {
	
	  Float_t E = hETotal1D->GetBinContent(j);
	  if(opt.Contains("units")) E *= eUnit;
	  else E *= E0;
	  if(E/(PUnits::GV/PUnits::m)<10) continue;

	  Float_t IonRate = PFunc::ADK_ENG(E,Eion0[iat],Z[iat]);
	  hIonRate1D[iat]->SetBinContent(j,IonRate);
	  integral[iat] += IonRate * dx;	
	  if(integral[iat]>=1) integral[iat] = 1;
	
	  if(opt.Contains("units"))
	    IonRate /= 1.0/PUnits::femtosecond;
	  else
	    IonRate /= omegap;

	  hIonProb1D[iat]->SetBinContent(j,(1-integral[iat])*IonRate);	
	}
      
      }

      cout << Form(" done. ") << endl; 
    }  
    
    // axis labels
    if(opt.Contains("units")) {
      char axName[36];
      sprintf(axName,"W_{%s} [fs^{-1}]",atAxNames[ii]);
      hIonRate2D[ii]->GetZaxis()->SetTitle(axName);
      sprintf(axName,"#Gamma_{%s} [fs^{-1}]",atAxNames[ii]);
      hIonProb2D[ii]->GetZaxis()->SetTitle(axName);

      sprintf(axName,"W_{%s} [fs^{-1}]",atAxNames[ii]);
      hIonRate1D[ii]->GetZaxis()->SetTitle(axName);
      sprintf(axName,"#Gamma_{%s} [fs^{-1}]",atAxNames[ii]);
      hIonProb1D[ii]->GetZaxis()->SetTitle(axName);
      
    } else {
      char axName[36];
      sprintf(axName,"W_{%s} #omega_{p}",atAxNames[ii]);
      hIonRate2D[ii]->GetZaxis()->SetTitle(axName);
      sprintf(axName,"#Gamma_{%s} #omega_{p}",atAxNames[ii]);
      hIonProb2D[ii]->GetZaxis()->SetTitle(axName);

      sprintf(axName,"W_{%s} #omega_{p}",atAxNames[ii]);
      hIonRate1D[ii]->GetZaxis()->SetTitle(axName);
      sprintf(axName,"#Gamma_{%s} #omega_{p}",atAxNames[ii]);
      hIonProb1D[ii]->GetZaxis()->SetTitle(axName);
    }
    
    
    
  }
  // End of ionization histograms.
  
  // ----------- Plotting Range (ZOOM) ------------

  TH2F *h2D = NULL;
  for(Int_t i=0;i<Nspecies;i++) {
    if(hDen2D[i]) {
      h2D = hDen2D[i];
      break;
    }
  }
  if(!h2D) {
    for(Int_t i=0;i<Nfields;i++) {
      if(hE2D[i]) {
	h2D = hE2D[i];
	break;
      }
    }
  }
  if(!h2D) {
    cout <<  Form("\n Error: base histogram could not be retrieved! -> exit!\n") << endl;
    return;
  }
    
  
  Float_t yRange   = (h2D->GetYaxis()->GetXmax() - h2D->GetYaxis()->GetXmin());
  Float_t midPoint = (h2D->GetYaxis()->GetXmax() + h2D->GetYaxis()->GetXmin())/2.;
  Float_t yMin = midPoint-yRange/2;
  Float_t yMax = midPoint+yRange/2;
  if(pData->IsCyl()) {
    yMin = pData->GetXMin(1);
    yMax = yRange;
  }
  Float_t xMin = h2D->GetXaxis()->GetXmin();
  Float_t xMax = h2D->GetXaxis()->GetXmax();
  Float_t xRange = xMax - xMin;
  
  
  for(Int_t i=0;i<Nspecies;i++) {
    if(i==noIndex) continue;
    if(hDen2D[i]) {
      hDen2D[i]->GetYaxis()->SetRangeUser(yMin,yMax);
      hDen2D[i]->GetXaxis()->SetRangeUser(xMin,xMax);
    }
    if(hJz2D[i]) {
      hJz2D[i]->GetYaxis()->SetRangeUser(yMin,yMax);
      hJz2D[i]->GetXaxis()->SetRangeUser(xMin,xMax);
    }
  }
  
  for(Int_t i=0;i<Nfields;i++) {
    if(!hE2D[i]) continue;
    hE2D[i]->GetYaxis()->SetRangeUser(yMin,yMax);
    hE2D[i]->GetXaxis()->SetRangeUser(xMin,xMax);
  }

  for(Int_t i=0;i<Nfields;i++) {
    if(!hB2D[i]) continue;
    hB2D[i]->GetYaxis()->SetRangeUser(yMin,yMax);
    hB2D[i]->GetXaxis()->SetRangeUser(xMin,xMax);
  }

  if(hETotal2D) {
    hETotal2D->GetYaxis()->SetRangeUser(yMin,yMax);
    hETotal2D->GetXaxis()->SetRangeUser(xMin,xMax);
  }

  
  // ----- z Zoom ---------- Plasma palette -----------
  // Set the range of the plasma charge density histogram for maximum constrast 
  // using a dynamic palette that adjusts the base value to a certain color.
  
  Float_t Base  = 1;
  Float_t baseden = Base;
  Float_t localden = baseden;
  
  if(hDen2D[0]) {
    
    // Local value of the density at \zeta = 0 and y = a quarter of the yRange.
    Int_t binx = hDen2D[0]->GetXaxis()->FindBin(0.);
    Int_t biny = hDen2D[0]->GetYaxis()->FindBin(yMax - ((yMax-yMin)/8.0));
    localden = hDen2D[0]->GetBinContent(binx,biny);
    
    cout << Form(" Local density = %.4f n0",localden) << endl;
  }
  
  if(opt.Contains("lden"))
    baseden = localden;
  
  Float_t *Max = new Float_t[Nspecies];
  Float_t *Min = new Float_t[Nspecies];
  
  for(Int_t i=0;i<Nspecies;i++) {
    if(i==noIndex) continue;
    if(!hDen2D[i]) continue;
    
    Max[i] = hDen2D[i]->GetMaximum();
    Min[i] = 1.01E-1 * baseden;
    if(Max[i]==0) {
      Min[i]=0;
      continue;
    }
    else if(Min[i]>Max[i]) Min[i] = 1.01E-1 * Max[i];
    
    if(i==0) {
      if(Max[i]<baseden) { 
	Max[i] = 1.1*baseden;
      } else if(Max[i]<2*baseden) {
	Min[i] = 2*baseden - Max[i];
      } else {
	Max[i] *= 1.0; // enhance plasma contrast.
        // Max[i] = hDen2D[i]->GetMaximum(); 
	if(Max[i]<baseden) Max[i] = 1.01 * baseden;
	if(Min[i]>baseden) Min[i] = 0.9*baseden;
      }
      // cout << Form("Base = %f, Max = %f, Min = %f",Base,Max[i], Min[i]) << endl;
    } 
    
    if(i==1) {
      Min[i] = 1.01E-1 * Base;
      if(Max[i]<Min[i])
	Min[i] = 1.01E-3 * Base;

    }

    if(i==2) {
      Min[i] = 1.01E-3 * Base;
      //   //Min[i] = 5.01E-3 * Base;
      //   //Min[i] = 50.01E-2 * Base;
    }

    if(i==3) {
      Min[i] = 1.01E-3 * Base;
    }

    if(pData->GetDenMax(i)>0)
      Max[i] = pData->GetDenMax(i);
    else if(pData->GetDenMax(i)==0) {
      // here there was a problem when using <= instead of ==
      // the reason is that denMax parameters are negative (-999.) by default
      if(hDen2D[i]) {
	delete hDen2D[i];
	hDen2D[i] = NULL;
      }

      if(hDen1D[i]) {
	delete hDen1D[i];
	hDen1D[i] = NULL;
      }

      if(hCur1D[i]) {
	delete hCur1D[i];
	hCur1D[i] = NULL;
      }

    }

    if(hDen2D[i]->GetMaximum() < pData->GetDenMin(i)) {
      if(hDen2D[i]) {
	delete hDen2D[i];
	hDen2D[i] = NULL;
      }

      if(hDen1D[i]) {
	delete hDen1D[i];
	hDen1D[i] = NULL;
      }

      if(hCur1D[i]) {
	delete hCur1D[i];
	hCur1D[i] = NULL;
      }
    }

    if(pData->GetDenMin(i)>0) {
      Min[i] = pData->GetDenMin(i);
      if(pData->GetDenMax(i)==0) {
	if(opt.Contains("logz") && (Max[i]<10*Min[i])) Max[i] = 1.01E1 * Min[i];
      }
    }
    
    if(pData->GetDenLoc(i)>0) {
      localden = pData->GetDenLoc(i);
      baseden = localden;
      cout << Form(" Species %2i  denMin = %5f  denMax = %5f  denLoc = %5f",i,Min[i],Max[i],localden) << endl;
    } else {
      cout << Form(" Species %2i  denMin = %5f  denMax = %5f",i,Min[i],Max[i]) << endl;
    }
    
    if(hDen2D[i])
      hDen2D[i]->GetZaxis()->SetRangeUser(Min[i],Max[i]);  

  }

  // Dynamic plasma palette
  PPalette * plasmaPalette = (PPalette*) gROOT->FindObject("plasma");
  if(!plasmaPalette) {
    plasmaPalette = new PPalette("plasma");
    plasmaPalette->SetPalette("gray");
  }
  
  Float_t basePos = 0.5;
  Float_t localPos = basePos;
  Float_t onePos = basePos;
  if( hDen2D[0] ) {
    const Int_t NRGBs = 3;
    const Int_t NCont = 255;
    if(Max[0]!=Min[0]) {
      if(opt.Contains("logz")) {
	Float_t a = 1.0/(TMath::Log10(Max[0])-TMath::Log10(Min[0]));
	Float_t b = TMath::Log10(Min[0]);
	basePos = a*(TMath::Log10(baseden) - b);
	localPos = a*(TMath::Log10(localden) - b);
	onePos = a*(TMath::Log10(1.0) - b);
      } else {
	basePos = (1.0/(Max[0]-Min[0]))*(baseden - Min[0]);
	localPos = (1.0/(Max[0]-Min[0]))*(localden - Min[0]);
	onePos = (1.0/(Max[0]-Min[0]))*(1.0 - Min[0]);
      }
      // cout << Form(" Local density = %.2f    Local position  = %f ",localden,localPos) << endl;
    }
    
    Double_t Stops[NRGBs] = { 0.00, basePos, 1.00 };
    Double_t Red[NRGBs]   = { 0.99, 0.90, 0.00 };
    Double_t Green[NRGBs] = { 0.99, 0.90, 0.00 };
    Double_t Blue[NRGBs]  = { 0.99, 0.90, 0.00 };
    
    plasmaPalette->CreateGradientColorTable(NRGBs, Stops, Red, Green, Blue, NCont, 1.0);
  }
  
	       
  // Redefines the ground color of the electron palette to match background plasma.
  PPalette * beamPalette = (PPalette*) gROOT->FindObject("beam");
  if(!beamPalette) {
    beamPalette = new PPalette("beam");
    beamPalette->SetPalette("elec");
  }
  
  if(opt.Contains("mbeam") && hDen2D[0]) {
    //  Int_t localcolorindex =  TMath::Nint(localPos * plasmaPalette->GetNColors());
    Int_t localcolorindex =  TMath::Nint(onePos * plasmaPalette->GetNColors());
    Int_t rootcolorindex = plasmaPalette->GetColor(localcolorindex);
    TColor *localcolor = gROOT->GetColor(rootcolorindex);
    Float_t r,g,b;
    localcolor->GetRGB(r,g,b);

    cout << Form(" N colors = %i   Local position = %f  Color index  = %i  RGB = (%.2f,%.2f,%.2f)", plasmaPalette->GetNColors(), localPos, rootcolorindex, r, g, b) << endl;

    const Int_t elecNRGBs = 5;
    const Int_t elecNCont = 255;
    Double_t elecStops[elecNRGBs] = { 0.00, 0.40, 0.50, 0.60, 1.00};
    Double_t elecRed[elecNRGBs] =   { r, 0.22, 0.39, 0.70, 1.00};
    Double_t elecGreen[elecNRGBs] = { g, 0.34, 0.05, 0.20, 1.00};
    Double_t elecBlue[elecNRGBs] =  { b, 0.58, 0.33, 0.30, 0.20};
    
    beamPalette->ChangeGradientColorTable(elecNRGBs, elecStops, elecRed, elecGreen, elecBlue);
  }
  
  PPalette * beam2Palette = (PPalette*) gROOT->FindObject("beam2");
  if(!beam2Palette) {
    beam2Palette = new PPalette("beam2");
  }

  PPalette * beam3Palette = (PPalette*) gROOT->FindObject("beam3");
  if(!beam3Palette) {
    beam3Palette = new PPalette("beam3");
  }

  // Change the range of z axis for the fields to be symmetric.
  Float_t *Emax = new Float_t[Nfields];
  Float_t *Emin = new Float_t[Nfields];
  for(Int_t i=0;i<Nfields;i++) {
    if(!hE2D[i]) continue;
    
    Emax[i] = hE2D[i]->GetBinContent(hE2D[i]->GetMaximumBin());
    Emin[i] = hE2D[i]->GetBinContent(hE2D[i]->GetMinimumBin());
    if(Emax[i] > TMath::Abs(Emin[i]))
      Emin[i] = -Emax[i];
    else
      Emax[i] = -Emin[i];

    if(i==0) {
      if(pData->GetE1Min() != 999.) {
        Emin[0] = pData->GetE1Min();
	if(opt.Contains("units"))
	  Emin[0] *= E0 / eUnit; 
      }

      if(pData->GetE1Max() != -999.) {
        Emax[0] = pData->GetE1Max();
	if(opt.Contains("units"))
	  Emax[0] *= E0 / eUnit; 
      }
      
    } else {
      Emin[i] = Emin[0];
      Emax[i] = Emax[0];
    }
    
    hE2D[i]->GetZaxis()->SetRangeUser(Emin[i],Emax[i]); 
  }

  Float_t ETmin = 0.0001;  
  Float_t ETmax = 1;
  if(hETotal2D) {
    ETmax = hETotal2D->GetBinContent(hETotal2D->GetMaximumBin());
    hETotal2D->GetZaxis()->SetRangeUser(ETmin,ETmax);
  }
  
  Float_t Fmin, Fmax;
  if(hFocus2D) {
    Fmax = hFocus2D->GetBinContent(hFocus2D->GetMaximumBin());
    Fmin = hFocus2D->GetBinContent(hFocus2D->GetMinimumBin());
  
    if(Fmax > TMath::Abs(Fmin))
      Fmin = -Fmax;
    else
      Fmax = -Fmin;

    Fmin = Emin[0];
    Fmax = Emax[0];
    hFocus2D->GetZaxis()->SetRangeUser(Fmin,Fmax);
  }
  
  // calculate \partial_x W_x
  TH2F *hdW2D = NULL;
  TH1F *hdW1D = NULL;

  Float_t kmax = -999.;
  Float_t kmin = 999.;
  if( (mask & 0x200) && hFocus2D) {

    hdW2D = (TH2F*) hFocus2D->Clone("hdW2D");
    hdW2D->Reset();

    TH2F *hField2D = hFocus2D;
    //TH2F *hField2D = hE2D[1];
    
    Int_t NBinsZ = hField2D->GetXaxis()->GetNbins();
    Int_t NBinsX = hField2D->GetYaxis()->GetNbins();
    Double_t dx = hField2D->GetYaxis()->GetBinWidth(1);
    for(Int_t i=1; i<=NBinsZ; i++) {
      for(Int_t j=1; j<=NBinsX; j++) {
	
	Float_t der = 0.;
	if(j>2 && j<NBinsX-2) {
	  der =  (4.0 / 3.0) * (hField2D->GetBinContent(i,j+1) - hField2D->GetBinContent(i,j-1)) / (2.0 * dx)
	    - (1.0 / 3.0) * (hField2D->GetBinContent(i,j+2) - hField2D->GetBinContent(i,j-2)) / (4.0 * dx) ;	  
	}

	// cout << Form(" Bin (%i,%i) = %f",i,j,der) << endl;
	if(opt.Contains("units"))
	  der *= (eUnit / spaUnit) / (E0 / skindepth );  
	hdW2D->SetBinContent(i,j,der);
	
	if(der>kmax) kmax = der;
	if(der<kmin) kmin = der;
      }
    }

    //    hdW2D->ResetStats();
    Int_t firstBin = NBinsX/2 - 1;
    Int_t lastBin = NBinsX/2 + 1;
    hdW1D = (TH1F*) hdW2D->ProjectionX(hName,firstBin,lastBin);
    hdW1D->Scale(1.0/(lastBin-firstBin+1));
  }

  // CROSSINGS
  
  // Find the first point on-axis where Ez changes from positive to negative:
  Int_t MAXCROSS = 2;
  Float_t *EzCross = new Float_t[MAXCROSS];
  Float_t *EzExtr = new Float_t[MAXCROSS];
  memset(EzCross,0,sizeof(Float_t)*MAXCROSS);
  memset(EzExtr,0,sizeof(Float_t)*MAXCROSS);

  Int_t auxNcross;
  if(hE1D[0])
    auxNcross = PGlobals::HCrossings(hE1D[0],EzCross,EzExtr,MAXCROSS,0.,0.);
  
  // for(Int_t i=0;i<auxNcross;i++) {
  //   if(opt.Contains("units"))
  //     cout << Form(" %i Ez crossing found at \\zeta = %.2f -> Ez maximum = %.2f GV/m",i,EzCross[i],EzExtr[i]) << endl;
  //   else
  //     cout << Form(" %i Ez crossing found at \\zeta = %.2f -> Ez maximum = %.2f E0",i,EzCross[i],EzExtr[i]) << endl;
      
  // }

  // Find maximum in beams region
  TH1F *hEzclone = (TH1F*) hE1D[0]->Clone("hEzclone");
  hEzclone->GetXaxis()->SetRangeUser(EzCross[0],hEzclone->GetXaxis()->GetXmax());
  Int_t binEbMax = hEzclone->GetMaximumBin();
  Float_t zEbMax = hEzclone->GetBinCenter(binEbMax);
  
  // Potential 
  if(opt.Contains("units") && n0) {
    trapPotential *=  ( E0 * skindepth / (PUnits::MV) ); 
  }  

  // Find the first point on-axis where focusing goes defocusing.
  MAXCROSS = 1;
  Float_t *FocusCross = new Float_t[MAXCROSS];
  Float_t *FocusExtr = new Float_t[MAXCROSS];
  memset(FocusCross,0,sizeof(Float_t)*MAXCROSS);
  memset(FocusExtr,0,sizeof(Float_t)*MAXCROSS);

  Int_t binFcross = -1;
  if(hFocus1D) {
    auxNcross = PGlobals::HCrossings(hFocus1D,FocusCross,FocusExtr,MAXCROSS,0.0,EzCross[0]);
    binFcross = hFocus1D->FindBin(FocusCross[0]);
  }
  
  Float_t potValueFin = hV1D->GetBinContent(binFcross);
  Float_t EmaxIon = hE1D[0]->GetBinContent(binFcross);
  
  Int_t binVmin = hV1D->GetMinimumBin();  
  Int_t binVmax = hV1D->GetMaximumBin();  
  Float_t Vmin = hV1D->GetMinimum();   
  Float_t Vmax = hV1D->GetMaximum();    
 

  Int_t binVfin = binVmin;
  Float_t Vfin = Vmin;
  if(opt.Contains("vfoc")) { // minimum in focusing
    Vfin = potValueFin;
    binVfin = binFcross;
  }
  
  // Find first bin over the trapping threshold (if any)
  Int_t binPotValueIni = -99;  
  Float_t potValueIni = -99;
  for(Int_t j=binVfin;j<hV1D->GetNbinsX();j++) {
    if(hV1D->GetBinContent(j) >= Vfin + trapPotential) {
      binPotValueIni = j;
      potValueIni = hV1D->GetBinContent(j);
      break;
    }
  }

  PPalette * potPalette = (PPalette*) gROOT->FindObject("pot");
   
  // Dynamic potential palette (blue values indicate trapping volume in respect to the minimum.
  if(binPotValueIni>0 && opt.Contains("trap")) {
    
    { // Shift potential value in respect to the minimum
      Int_t NbinsX = hV2D->GetNbinsX(); 
      Int_t NbinsY = hV2D->GetNbinsY();
      for(Int_t j=0;j<=NbinsX;j++) {
	for(Int_t k=0;k<=NbinsY;k++) {
	  hV2D->SetBinContent(j,k, hV2D->GetBinContent(j,k) - Vfin);
	}
	hV1D->SetBinContent(j, hV1D->GetBinContent(j) - Vfin);
      }
      Vmin -= Vfin;
      Vmax -= Vfin;
    }
    
    Float_t Vzero = hV1D->GetBinContent(binPotValueIni);
    // cout << Form(" VZERO = %f",Vzero) << endl;
    
    const Int_t potPNRGBs = 5;
    const Int_t potPNCont = 255;
    Float_t zeroPos = (Vzero-Vmin)/(Vmax-Vmin);

    Double_t potPStops[potPNRGBs] = { 0.00,zeroPos-0.05, zeroPos, zeroPos+0.05, 1.00 };
    Double_t potPRed[potPNRGBs]   = { 0.518, 0.965, 0.90, 0.498, 0.106 };
    Double_t potPGreen[potPNRGBs] = { 0.078, 0.925, 0.90, 0.718, 0.078 };
    Double_t potPBlue[potPNRGBs]  = { 0.106, 0.353, 0.90, 0.780, 0.518 };
   
    potPalette->CreateGradientColorTable(potPNRGBs, potPStops, 
					 potPRed, potPGreen, potPBlue, potPNCont);

  } else {
    
    // const Int_t potPNRGBs = 4;
    // const Int_t potPNCont = 64;

    // Float_t Vzero = 0.;
    // Float_t zeroPos = (Vzero-Vmin)/(Vmax-Vmin);
    
    // Double_t potPStops[potPNRGBs] = { 0.00, zeroPos, zeroPos + (1.-zeroPos)/2. ,1.00 };
    // Double_t potPRed[potPNRGBs]   = {  1.0, 0.9, 0.965, 0.518};
    // Double_t potPGreen[potPNRGBs] = {  1.0, 0.9, 0.925, 0.078};
    // Double_t potPBlue[potPNRGBs]  = {  1.0, 0.9, 0.353, 0.106};
   
    // potPalette->CreateGradientColorTable(potPNRGBs, potPStops, 
    // 					 potPRed, potPGreen, potPBlue, potPNCont);

    if(Vmax > TMath::Abs(Vmin))
      Vmin = -Vmax;
    else
      Vmax = -Vmin;
    
  }

  // Extract contours from 2D histos
  TCanvas* c = new TCanvas("contours","Contour List",0,0,600,600);
  c->cd();
  
  // Potential 
  TH2F *hV2Dc = (TH2F*) hV2D->Clone("hV2Dc");
  const Int_t Ncontours = 25;
  Double_t contours[Ncontours];

  if(opt.Contains("trap")) {
    for(Int_t i=0; i<Ncontours; i++) {
      contours[i] = i*(trapPotential/5.0);// - trapPotential; 
    }
  } else {
    for(Int_t i=0; i<Ncontours; i++) {
      contours[i] = Vmin + i*(trapPotential/10.0);// - trapPotential; 
    }
  }
    
  hV2Dc->SetContour(Ncontours, contours);
  hV2Dc->Draw("cont list 0");
  
  c->Update();
  TObjArray *contsV2D = (TObjArray*) gROOT->GetListOfSpecials()->FindObject("contours");
  TClonesArray graphsV2D("TGraph",Ncontours);
  TClonesArray graphsV2D_main("TGraph",Ncontours);
  {
    Int_t ncontours = contsV2D->GetSize();
    TList* clist = NULL;
    Int_t nGraphs = 0;
    Int_t nGraphs_main = 0;
    TGraph *gr = NULL;
    for(Int_t i = 0; i < ncontours; i++){  
      // negative contours are returned first from low to high absolute values. 
      // e.g. -1, -2, -3, ..., 0, 1, 2, 3,...
            
      clist = (TList*) contsV2D->At(i);
    
      for(Int_t j = 0 ; j < clist->GetSize(); j++) {
	//	cout << Form(" psi contour (%i,%i)",i,j) << endl;

	gr = (TGraph*) clist->At(j);
	if(!gr) continue;
      
	gr->SetLineWidth(1);
	gr->SetLineColor(kGray+1);

	if(opt.Contains("trap")) {
	  if( i==0 || i==5 || i==10 ) {
	    gr->SetLineWidth(2);
	    gr->SetLineColor(kGray+2);
	  }
	}
	
	new(graphsV2D[nGraphs]) TGraph(*gr) ;
	nGraphs++;
	
	//	if(i==0 && j==0) {
	if(i==0) {
	  TGraph *grm = new(graphsV2D_main[nGraphs_main]) TGraph(*gr);
	  nGraphs_main++;
	  grm->SetLineWidth(1);
	  grm->SetLineStyle(3);
	  grm->SetLineColor(kGray+3);	  
	}

	if(i==5 && j==0) {
	  TGraph *grm = new(graphsV2D_main[nGraphs_main]) TGraph(*gr);
	  nGraphs_main++;
	  grm->SetLineWidth(2);
	  grm->SetLineStyle(3);
	  grm->SetLineColor(kGray+3);	  
	}

      }
    }
  }
  
  // Get ionization probablity rates contours
  TObjArray *graphsI2D = NULL;
  TObjArray *graphsI2D_main = NULL;
  Float_t IPmin, IPmax;
  if(hIonProb2D[ii]) {
    IPmin = 0.00;
    IPmax = hIonProb2D[ii]->GetMaximum();

    const Int_t NcontI = 1;
    //  Double_t contI[NcontI] = {0.1*IPmax};
    Double_t contI[NcontI] = {0.2*IPmax};
    // Double_t contI[NcontI] = {0.1*IPmax};
    //const Int_t NcontI = 1;
    //Float_t contI[NcontI] = {0.1};
  
    IPmax = hIonProb1D[ii]->GetMaximum();
    hIonProb2D[ii]->GetZaxis()->SetRangeUser(IPmin,IPmax);
    
    TH2F *hIonProb2Dc = (TH2F*) hIonProb2D[ii]->Clone("hIonProb2Dc");
  
    hIonProb2Dc->SetContour(NcontI, contI);
    hIonProb2Dc->Draw("cont list 0");
  
    c->Update();

    TObjArray *contsI2D = (TObjArray*) gROOT->GetListOfSpecials()->FindObject("contours");
    graphsI2D = new TObjArray();
    graphsI2D_main = new TObjArray();
    {
      Int_t ncontours = contsI2D->GetSize();
      TList* clist = NULL;
      TGraph *gr = NULL;
      for(Int_t i = 0; i < ncontours; i++){
	clist = (TList*) contsI2D->At(i);
      
	for(Int_t j = 0 ; j < clist->GetSize(); j++) {
	  gr = (TGraph*) clist->At(j);
	  if(!gr) continue;
	
	  if( i==0 ) {
	    TGraph *grm = new TGraph(*gr);
	    grm->SetLineWidth(1);
	    grm->SetLineStyle(2);
	    grm->SetLineColor(kGray+2);

	    graphsI2D_main->Add(grm);
	    
	    TGraph *gr2 = new TGraph(*gr);
	    if(j==0) {
	      gr2->SetLineWidth(1);
	      gr2->SetLineStyle(2);
	      //	    gr2->SetLineColor(PGlobals::elecLine);
	      gr2->SetLineColor(kGray+2);
	    } else {
	      gr2->SetLineWidth(1);
	      gr2->SetLineStyle(1);
	      gr2->SetLineColor(kGray+2);
	    }
	    graphsI2D->Add(gr2);
	    
	  } 
	
	}
      }
    }
  }
  
  delete c;


  // Plotting
  // ------------------------------------------------------------
 
  // Canvas setup
  Int_t NPad = BitCounter(mask);
  if(NPad==0) NPad = 1;
  Float_t ypadsize = 250;
  Float_t ymarginsize = 200;
  if(NPad==1) ypadsize = 300;
  Float_t ysize = ypadsize * NPad + ymarginsize; 
  Float_t boom = 1.2;
  if(opt.Contains("boom"))
    ysize *= boom;
  TCanvas *C = new TCanvas("C","Snapshot",1050,ysize);
  C->SetFillStyle(4000);

  UInt_t lineColor = kOrange+10;
  //  UInt_t lineColor2 =  TColor::GetColor(196,30,78);
  //  UInt_t lineColor2 = kGray+2;
  UInt_t lineColor2 = lineColor;
  
  // Setup Pad layout:
  TPad **pad = new TPad*[NPad];
  TH2F **hFrame = new TH2F*[NPad];
  Float_t bMargin = 0.12 * (950/ysize);
  Float_t tMargin = 0.04 * (950/ysize);
  Float_t lMargin = 0.14;
  Float_t rMargin = 0.18;
  Float_t mMargin = 0.02 * (950/ysize);
  Float_t pfactor = 1.6;
  if(opt.Contains("nomar"))
    bMargin = tMargin = lMargin = rMargin = mMargin = 0.0;
  if(NPad==1)
    PGlobals::CanvasPartition(C,NPad,lMargin,rMargin,bMargin,tMargin,mMargin);
  else
    PGlobals::CanvasAsymPartition(C,NPad,lMargin,rMargin,bMargin,tMargin,pfactor,mMargin);


  // Define the frames for plotting
  Int_t fonttype = 43;
  Int_t fontsize = 38;
  Int_t tfontsize = 42;
  Int_t txsize = tfontsize+6;
  Int_t lxsize = fontsize+2;
  Int_t tysize = tfontsize;
  Int_t lysize = fontsize;
  Int_t tzsize = tfontsize-4;
  Int_t lzsize = fontsize-6;
  Float_t txoffset = (250/ypadsize) * 2.4 / (950/ysize);
  Float_t lxoffset = 0.015;
  Float_t tyoffset = 1.2 / (950/ysize);
  Float_t lyoffset = 0.01;
  Float_t tzoffset = 1.45 / (950/ysize);
  Float_t lzoffset = 0.01;
  Float_t tylength = 0.015;
  Float_t txlength = 0.04;
  for(Int_t i=NPad-1;i>=0;i--) {
    char name[16];
    sprintf(name,"pad_%i",i);
    pad[i] = (TPad*) gROOT->FindObject(name);
    //pad[i]->SetFrameLineWidth(2);
    if(opt.Contains("tick")) {
      pad[i]->SetTickx(1);
      pad[i]->SetTicky(1);
    }
    if(opt.Contains("trans"))
      pad[i]->SetFillStyle(4000);
    pad[i]->SetFrameFillStyle(4000);
    
    sprintf(name,"hFrame_%i",i);
    hFrame[i] = (TH2F*) gROOT->FindObject(name);
    if(hFrame[i]) {
      delete hFrame[i];
      hFrame[i] = NULL;
    }

    hFrame[i] = (TH2F*) h2D->Clone(name);
    hFrame[i]->Reset();
    

    Float_t xFactor = pad[NPad-1]->GetAbsWNDC()/pad[i]->GetAbsWNDC();
    Float_t yFactor = pad[NPad-1]->GetAbsHNDC()/pad[i]->GetAbsHNDC();

    // Format for y axis
    hFrame[i]->GetYaxis()->SetLabelFont(fonttype);
    hFrame[i]->GetYaxis()->SetLabelSize(lysize);
    hFrame[i]->GetYaxis()->SetLabelOffset(lyoffset);
    hFrame[i]->GetYaxis()->SetTitleFont(fonttype);
    hFrame[i]->GetYaxis()->SetTitleSize(tysize);
    hFrame[i]->GetYaxis()->SetTitleOffset(tyoffset);

    if(i==NPad-1) hFrame[i]->GetYaxis()->SetTickLength(xFactor*tylength/(yFactor*pfactor));
    else  hFrame[i]->GetYaxis()->SetTickLength(xFactor*tylength/yFactor);

    // Format for x axis
    hFrame[i]->GetXaxis()->SetLabelFont(fonttype);
    hFrame[i]->GetXaxis()->SetLabelSize(lxsize);
    hFrame[i]->GetXaxis()->SetLabelOffset(lxoffset);
    hFrame[i]->GetXaxis()->SetTitleFont(fonttype);
    hFrame[i]->GetXaxis()->SetTitleSize(txsize);
    hFrame[i]->GetXaxis()->SetTitleOffset(txoffset);
    
    hFrame[i]->GetXaxis()->SetTickLength(yFactor*txlength/xFactor);      

    if(i>0) { // skip x axis labels except for the lowest one
      hFrame[i]->GetXaxis()->SetLabelSize(0.0);
      hFrame[i]->GetXaxis()->SetTitleSize(0.0);
    }

    if(opt.Contains("nomar")) {
      hFrame[i]->GetYaxis()->SetTickLength(0.0);
      hFrame[i]->GetXaxis()->SetTickLength(0.0);      
    }

    // Labels for the frames

  }

  // Text objects
  TPaveText **textLabel = new TPaveText*[NPad];
  TString pNames[] = {"Den","Ez","Ex","Ey","F","ET","I","V"};
  TString sLabels[] = {"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)"};
    
  char ctext[128];
  if(opt.Contains("units") && n0) 
    sprintf(ctext,"z = %5.2f %s", Time * skindepth / propUnit,propSUnit->c_str());
  else
    sprintf(ctext,"k_{p}z = %5.1f",Time);
  
  TLatex *textTime = new TLatex(xMax - (xMax-xMin)/20.,yMax - (yMax-yMin)/10.,ctext);
  textTime->SetTextAlign(32);
  textTime->SetTextFont(fonttype);
  textTime->SetTextSize(fontsize-10);

  if(opt.Contains("units") && n0) 
    sprintf(ctext,"n_{0} = %5.2f x %s", n0/denUnit,denSUnit->c_str());
   
  TLatex *textDen = new TLatex(xMin + (xMax-xMin)/20.,yMax - (yMax-yMin)/10.,ctext);
  textDen->SetTextAlign(12);
  textDen->SetTextFont(fonttype);
  //  textDen->SetTextColor(lineColor);
  textDen->SetTextColor(kGray+3);
  textDen->SetTextSize(fontsize-14);

  // Lines 
  TLine *lineYzero = new TLine(xMin,0.0,xMax,0.0);
  lineYzero->SetLineColor(kGray+2);
  lineYzero->SetLineStyle(2);

  if(opt.Contains("units") && n0) {
    zStartPlasma *= skindepth /spaUnit;
    zStartNeutral *= skindepth /spaUnit;
    zEndNeutral *= skindepth /spaUnit;
  }

  //  cout << "Start plasma = " << zStartPlasma << endl;
  TLine *lineStartPlasma = new TLine(zStartPlasma,yMin,zStartPlasma,yMax);
  lineStartPlasma->SetLineColor(kGray+2);
  lineStartPlasma->SetLineStyle(2);
  lineStartPlasma->SetLineWidth(3);

  //  cout << "Start neutral = " << zStartNeutral << endl;
  TLine *lineStartNeutral = new TLine(zStartNeutral,yMin,zStartNeutral,yMax);
  lineStartNeutral->SetLineColor(kGray+3);
  lineStartNeutral->SetLineStyle(3);
  lineStartNeutral->SetLineWidth(2);

  //  cout << "End neutral = " << zEndNeutral << endl;
  TLine *lineEndNeutral = new TLine(zEndNeutral,yMin,zEndNeutral,yMax);
  lineEndNeutral->SetLineColor(kGray+3);
  lineEndNeutral->SetLineStyle(3);
  lineEndNeutral->SetLineWidth(2);

  // Usable accelerating field (calculated from the real slope):
  Double_t EzMaxUse;
  {
    Float_t zshift = 0.1; 
    if(opt.Contains("units")) 
      zshift *=spaUnit * kp;
    
    Int_t bin1 = hE1D[0]->FindBin(EzCross[0]-zshift);
    Int_t bin2 = hE1D[0]->FindBin(EzCross[0]+zshift);
    Float_t EzSlope = (hE1D[0]->GetBinContent(bin2)-hE1D[0]->GetBinContent(bin1))/
      (hE1D[0]->GetBinCenter(bin2)-hE1D[0]->GetBinCenter(bin1));

    EzMaxUse  = -EzSlope * (EzCross[1]-EzCross[0]);
  }

  // Data about Ez
  //  cout << Form("Ez max (driver) = %.2f",EzExtr[0]) << endl;
  
  // Get position of the injected beam:
  Float_t zInjBeam = -999;
  Float_t zrmsInjBeam = -999;
  Float_t EInjBeam = -999;
  Float_t ESlopeInjBeam = -999;

  if(Nspecies>2) {
    if(hCur1D[2] && noIndex!=2) {
      hCur1D[2]->ResetStats();
      zInjBeam = hCur1D[2]->GetMean();
      zrmsInjBeam = hCur1D[2]->GetRMS();
      
      // Find field value corresponding to this position
      Int_t ibin0 = hE1D[0]->FindBin(zInjBeam);
      Int_t ibin1 = hE1D[0]->FindBin(zInjBeam-zrmsInjBeam);
      Int_t ibin2 = hE1D[0]->FindBin(zInjBeam+zrmsInjBeam);
      EInjBeam = hE1D[0]->GetBinContent(ibin0);
      ESlopeInjBeam = (hE1D[0]->GetBinContent(ibin1)-hE1D[0]->GetBinContent(ibin2))/
	(hE1D[0]->GetBinCenter(ibin1)-hE1D[0]->GetBinCenter(ibin2));

      cout << endl << Form(" Witness bunch properties") << endl;
      if(opt.Contains("units")) {
	cout << Form("  Position at z bunch = %6.2f um ", zInjBeam) << endl ;
	cout << Form("  RMS of the bunch    = %6.2f um ", zrmsInjBeam) << endl ;
	cout << Form("  E_z at z bunch      = %6.2f GV/m ", EInjBeam) << endl ;
	cout << Form("  Slope at z bunch    = %6.2f (GV/m)um ", ESlopeInjBeam) << endl ;
      } else {
	cout << Form("  Position at z bunch = %6.2f ", zInjBeam) << endl ;
	cout << Form("  RMS of the bunch    = %6.2f ", zrmsInjBeam) << endl ;
	cout << Form("  E_z at z bunch      = %6.2f ", EInjBeam) << endl ;
	cout << Form("  Slope at z bunch    = %6.2f ", ESlopeInjBeam) << endl ;
      }
      cout << Form("  Transformer ratio   = %6.2f ", EInjBeam/EzExtr[0]) << endl ;
    } 
    
  }

  // Calculate blowout parameters
  Double_t maxcurr;
  Double_t rmslength;
  Double_t radius;
  Double_t EzMaxLo;
  Double_t radiusLu;
  Double_t EzMaxLuLinear;
  Double_t EzMaxBeamLo;
  Double_t DeltaPsiLo;
  Double_t RadiusPsiLo;

  // Get max current in norm units:
  if(Nspecies>1 && opt.Contains("bopar")) {
    if(hCur1D[1] && noIndex!=1) {  
 
      hCur1D[1]->ResetStats();
      maxcurr = hCur1D[1]->GetMaximum();
      rmslength = hCur1D[1]->GetRMS();
      if(opt.Contains("units")) {
	maxcurr *= curUnit / PConst::I0;
	rmslength *= spaUnit * kp;
      }

      radius = PFunc::RadiusBO(maxcurr,rmslength);
      EzMaxLo = 0.5 * radius;
      radiusLu = 2.0 * TMath::Sqrt(maxcurr);
      EzMaxLuLinear = PFunc::EzMaxLu(maxcurr);
      EzMaxBeamLo = PFunc::EzMaxBeamLo(maxcurr);
      DeltaPsiLo = radius*radius/4.0;
      RadiusPsiLo = TMath::Sqrt(2*(DeltaPsiLo-1));

      cout << endl << Form(" Blowout parameters") << endl;
      if(opt.Contains("units")) {
	maxcurr *= PConst::I0 / curUnit;
	rmslength *= 1.0/(spaUnit * kp);
	radius   *= skindepth / spaUnit;
	radiusLu *= skindepth / spaUnit;
	EzMaxLo  *= E0 / eUnit;
	EzMaxLuLinear  *= E0 / eUnit;
	EzMaxBeamLo  *= E0 / eUnit;
	DeltaPsiLo   *= E0 * skindepth / (PUnits::MV);
	RadiusPsiLo  *= skindepth / spaUnit;
	
	cout << Form("  Beam max. curr.  = %7.2f %s   RMS = %.2f um",maxcurr,curSUnit->c_str(),rmslength) << endl;
	cout << Form("  Ez max beam      = %7.2f %s",EzExtr[0],eSUnit->c_str()) << endl;
	cout << Form("  Radius (Lotov) R = %7.2f %s",radius,spaSUnit->c_str()) << endl;
	cout << Form("  Radius (Lu)    R = %7.2f %s",radiusLu,spaSUnit->c_str()) << endl;
	cout << Form("  Ez max beam (Lo) = %7.2f %s",EzMaxBeamLo,eSUnit->c_str()) << endl;
	cout << Form("  Ez max acc. (Lo) = %7.2f %s",-EzMaxLo,eSUnit->c_str()) << endl;
	cout << Form("  Ez max (Lu)      = %7.2f %s",-EzMaxLuLinear,eSUnit->c_str()) << endl;
	cout << Form("  E  max (ion)     = %7.2f %s",EmaxIon,eSUnit->c_str()) << endl;
	cout << Form("  DPsi max. (Lo)   = %7.2f MV",DeltaPsiLo) << endl;
	cout << Form("  Trapping R (Lo)  = %7.2f %s",RadiusPsiLo,spaSUnit->c_str()) << endl;
      } else {
	cout << Form("  Beam max. curr.  = %7.2f I0   RMS = %.2f kp^-1",maxcurr,rmslength) << endl;
	cout << Form("  Ez max beam      = %7.2f E0",EzExtr[0]) << endl;
	cout << Form("  Radius (Lotov) R = %7.2f kp^-1",radius) << endl;
	cout << Form("  Radius (Lu)    R = %7.2f kp^-1",radiusLu) << endl;
	cout << Form("  Ez max beam (Lo) = %7.2f E0",EzMaxBeamLo) << endl;
	cout << Form("  Ez max acc. (Lo) = %7.2f E0",-EzMaxLo) << endl;
	cout << Form("  Ez max acc. (Lu) = %7.2f E0",-EzMaxLuLinear) << endl;
	cout << Form("  E  max (ion)     = %7.2f E0",EmaxIon) << endl;
	cout << Form("  DPsi max.   (Lo) = %7.2f \\Psi0",DeltaPsiLo) << endl;
	cout << Form("  Trapping R  (Lo) = %7.2f kp^-1",RadiusPsiLo) << endl;
      }
    }  
  }
  
  cout << endl;

  // Access to color Palettes
  PPalette * laserPalette = (PPalette*) gROOT->FindObject("laser");
  if(!laserPalette) {
    laserPalette = new PPalette("laser");
    // laserPalette->SetPalette(kBird);
    laserPalette->SetPalette(kColorPrintableOnGrey);
    laserPalette->Invert();
  }
  
  TExec *exPlasma = new TExec("exPlasma","plasmaPalette->cd();");
  TExec *exBeam   = new TExec("exBeam","beamPalette->cd();");
  TExec *exBeam2  = new TExec("exBeam2","beam2Palette->cd();");
  TExec *exBeam3  = new TExec("exBeam3","beam3Palette->cd();");
  TExec *exLaser  = new TExec("exLaser","laserPalette->cd();");
  TExec *exField  = new TExec("exField","fieldPalette->cd();");
  TExec *exFieldT = new TExec("exFieldT","fieldTPalette->cd();");
  TExec *exPot    = new TExec("exPot","potPalette->cd();");
  TExec *exIonP   = new TExec("exIonP","ionPalette->cd();");

  if(opt.Contains("alpha0"))
    plasmaPalette->SetAlpha(0.8);
  if(opt.Contains("alpha1"))
    beamPalette->SetAlpha(0.8);
  if(opt.Contains("alpha2"))
    beam2Palette->SetAlpha(0.8);
  if(opt.Contains("alpha3"))
    beam3Palette->SetAlpha(0.8);
  if(opt.Contains("alphal"))
    laserPalette->SetAlpha(0.8);
  
  // Current lines colors
  //Int_t cBeam = TColor::GetColor(69,108,155);
  Int_t cBeam = kGray+2;
  Int_t cBeam1 = kOrange+8;
  Int_t cBeam2 = kOrange+8;

  if(Nspecies==4) {
    //    beam2Palette->SetPalette("red");
    // cBeam1 = TColor::GetColor("#FFCC66"); 
    // beam2Palette->SetPalette("blue");
    // cBeam1 = TColor::GetColor("#84B1C4");
    // beam2Palette->SetPalette("greengray");
    // cBeam1 = TColor::GetColor("#5C7B2F"); // Green
    beam2Palette->SetPalette("elec0");
    cBeam1 = TColor::GetColor(83,109,161); // Bluish
    //cBeam1 = TColor::GetColor(192,114,49); // dark orange
    //cBeam2 = TColor::GetColor(137,158,107); // dark green
    beam3Palette->SetPalette("hot");
  }


  
  TString drawopt = "colz same";

  // Actual Plotting!
  // ------------------------------------------------------------

  C->cd(0);

  Float_t x1,x2,y1,y2;
  Float_t gap = 0.02;
  TPaletteAxis *palette = NULL;
  UInt_t ip = NPad-1;

  if(mask == 0) {
        
    pad[ip]->Draw();
    pad[ip]->cd();  // <---------------------------------------------------------- Just frame 
    if(opt.Contains("logz")) {
      pad[ip]->SetLogz(1);
    } else {
      pad[ip]->SetLogz(0);
    }
    
    hFrame[ip]->Draw("col");

    if(opt.Contains("contall")) {
      // ADK contours
      if(graphsI2D) {
	for(Int_t i=0;i<graphsI2D->GetEntriesFast();i++) {
	  TGraph *gr = (TGraph*) graphsI2D->At(i);
	  if(!gr) continue;
	  gr->Draw("C");
	}
      }
      
      // PSI contours  
      for(Int_t i=0;i<graphsV2D.GetEntriesFast();i++) {
	TGraph *gr = (TGraph*) graphsV2D.At(i);
	if(!gr) continue;
	gr->Draw("C");
      }
    } else if(opt.Contains("cont")) {
	
      // ADK contours  
      for(Int_t i=0;i<graphsI2D_main->GetEntriesFast();i++) {
       	TGraph *gr = (TGraph*) graphsI2D_main->At(i);
       	if(!gr) continue;
	gr->Draw("C");
      }
      
      // PSI MAIN contour  
      for(Int_t i=0;i<graphsV2D_main.GetEntriesFast();i++) {
	TGraph *gr = (TGraph*) graphsV2D_main.At(i);
	if(!gr) continue;
	gr->Draw("C");
      }
      
    } 
    
    
    ip--;  
    C->cd(0);
  }

  if(mask & 0x1) {
    
    pad[ip]->Draw();
    pad[ip]->cd();  // <---------------------------------------------------------- 1st panel
    if(opt.Contains("logz")) {
      pad[ip]->SetLogz(1);
    } else {
      pad[ip]->SetLogz(0);
    }
    
    hFrame[ip]->Draw();   

    // Set the z-axis properties

    // Count the species to plot
    Int_t Ns = 0;
    for(Int_t i=0;i<Nspecies;i++) {
      if(i==noIndex) continue;
      if(!hDen2D[i]) continue;
      Ns++;
    }
    
    for(Int_t i=0;i<Nspecies;i++) {
      if(i==noIndex) continue;
      if(!hDen2D[i]) continue;
      
      //  hDen2D[i]->GetZaxis()->SetNdivisions(503);
      hDen2D[i]->GetZaxis()->SetTitleFont(fonttype);
      hDen2D[i]->GetZaxis()->SetTitleOffset(tzoffset);
      hDen2D[i]->GetZaxis()->SetTitleSize(tzsize);
      hDen2D[i]->GetZaxis()->SetLabelFont(fonttype);
      hDen2D[i]->GetZaxis()->SetLabelSize(lzsize);
      hDen2D[i]->GetZaxis()->SetTickLength(Ns*tylength/(pfactor));     
      if(opt.Contains("logz")) { 
	hDen2D[i]->GetZaxis()->SetLabelOffset(-lyoffset * (250/ypadsize));
      } else
	hDen2D[i]->GetZaxis()->SetLabelOffset(lyoffset);
      
    }

    if(opt.Contains("laser") && hA2D) {
      hA2D->GetZaxis()->SetTitleFont(fonttype);
      hA2D->GetZaxis()->SetTitleOffset(tzoffset);
      hA2D->GetZaxis()->SetTitleSize(tzsize);
      hA2D->GetZaxis()->SetLabelFont(fonttype);
      hA2D->GetZaxis()->SetLabelSize(lzsize);
      hA2D->GetZaxis()->SetTickLength(Ns*tylength/(pfactor));     
      if(opt.Contains("logz")) { 
	hA2D->GetZaxis()->SetLabelOffset(-lyoffset * (250/ypadsize));
      } else
	hA2D->GetZaxis()->SetLabelOffset(lyoffset);
      
    }
        
    // cout << Form(" N species = %i",Nspecies) << endl;

    if(Nspecies==1) {
      // Just plasma
      if(hDen2D[0] && noIndex!=0) {
	exPlasma->Draw();
	hDen2D[0]->Draw(drawopt);
      }
    } else if(Nspecies==2) {
      // Plasma
      if(hDen2D[0] && noIndex!=0) {
	exPlasma->Draw();
	hDen2D[0]->Draw(drawopt);
      }

      // Beam driver.
      if(hDen2D[1] && noIndex!=1) {
	exBeam->Draw();
	hDen2D[1]->Draw(drawopt);
      }
    } else if (Nspecies==3) {
  
      // Injected electrons ?
      if(hDen2D[2] && noIndex!=2) {
	exBeam2->Draw();
	//exBeam->Draw();
	hDen2D[2]->Draw(drawopt);
      }

      // Plasma
      if(hDen2D[0] && noIndex!=0) {
	exPlasma->Draw();
	hDen2D[0]->Draw(drawopt);
      }

      // Beam driver.
      if(hDen2D[1] && noIndex!=1) {
	exBeam->Draw();
	//exPlasma->Draw();
	hDen2D[1]->Draw(drawopt);
      }

    } else if (Nspecies==4) {
      
      // Witness 2
      if(hDen2D[3] && noIndex!=3) {
	exBeam3->Draw();
	hDen2D[3]->Draw(drawopt);
      }

      // Witness
      if(hDen2D[2] && noIndex!=2) {
	exBeam2->Draw();
	//exBeam->Draw();
	hDen2D[2]->Draw(drawopt);
      }

      // Plasma
      if(hDen2D[0] && noIndex!=0) {
	exPlasma->Draw();
	hDen2D[0]->Draw(drawopt);
      }

      // Beam driver.
      if(hDen2D[1] && noIndex!=1) {
	exBeam->Draw();
	//exPlasma->Draw();
	hDen2D[1]->Draw(drawopt);
      }
    }

    if(hA2D && opt.Contains("laser")) {
      exLaser->Draw();

      Float_t aMax = hA2D->GetMaximum();
      Float_t aMin = 0.1*aMax;
      if(pData->GetaMax()>-999) { 
	aMax = pData->GetaMax();
	aMin = 0.1*aMax;
      }
      
      if(pData->GetaMin()<999) 
	aMin = pData->GetaMin();
      
      cout << Form(" Laser  aMin = %5f  aMax = %5f",aMin,aMax) << endl;
      
      hA2D->GetZaxis()->SetRangeUser(aMin,aMax);
      hA2D->Draw(drawopt);
    }
    
    if(opt.Contains("contall")) {
      // ADK contours
      if(graphsI2D) {
	for(Int_t i=0;i<graphsI2D->GetEntriesFast();i++) {
	  TGraph *gr = (TGraph*) graphsI2D->At(i);
	  if(!gr) continue;
	  gr->Draw("C");
	}
      }

      // PSI contours  
      for(Int_t i=0;i<graphsV2D.GetEntriesFast();i++) {
	TGraph *gr = (TGraph*) graphsV2D.At(i);
	if(!gr) continue;
	gr->Draw("C");
      }
    } else if(opt.Contains("conti")) {
      // ADK MAIN contour
      if(graphsI2D_main) {
	for(Int_t i=0;i<graphsI2D_main->GetEntriesFast();i++) {
	  TGraph *gr = (TGraph*) graphsI2D_main->At(i);
	  if(!gr) continue;
	  gr->Draw("C");
	}  
      }
      
      // PSI MAIN contour  
      for(Int_t i=0;i<graphsV2D_main.GetEntriesFast();i++) {
	TGraph *gr = (TGraph*) graphsV2D_main.At(i);
	if(!gr) continue;
	gr->Draw("C");
      }
    }  else if(opt.Contains("cont")) {
      // PSI MAIN contour  
      for(Int_t i=0;i<graphsV2D_main.GetEntriesFast();i++) {
	TGraph *gr = (TGraph*) graphsV2D_main.At(i);
	if(!gr) continue;
	gr->Draw("C");
      }
    } 



    if(opt.Contains("bopar")) {
      TEllipse *ellipLo = new TEllipse(EzCross[0],0.0,radius);	
      ellipLo->SetFillStyle(0);
      ellipLo->SetLineStyle(3);
      ellipLo->SetLineColor(kGray+3);
      ellipLo->Draw();
	
      TEllipse *ellipLu = new TEllipse(EzCross[0],0.0,radiusLu);
      ellipLu->SetFillStyle(0);
      ellipLu->SetLineStyle(4);
      ellipLu->SetLineColor(kGray+3);
      ellipLu->Draw();
    }

    if(opt.Contains("sline")) {
      if(zStartPlasma>xMin && zStartPlasma<xMax)
	lineStartPlasma->Draw();
      if(zStartNeutral>xMin && zStartNeutral<xMax)
	lineStartNeutral->Draw();
      if(zEndNeutral>xMin && zEndNeutral<xMax)
	lineEndNeutral->Draw();
    }
    
    // Palettes re-arrangement
    // Count the species to be plot
    Ns = 0;
    for(Int_t i=0;i<Nspecies;i++) {
      if(i==noIndex) continue;
      if(!hDen2D[i]) continue;
      Ns++;
    }
    if(hA2D && opt.Contains("laser")) {
      Ns++;
    }
   
    pad[ip]->Update();
    y1 = pad[ip]->GetBottomMargin();
    y2 = 1 - pad[ip]->GetTopMargin();
    x1 = pad[ip]->GetLeftMargin();
    x2 = 1 - pad[ip]->GetRightMargin();
    Float_t pvsize = (y2-y1)/float(Ns);

    Int_t j = 0;
    for(Int_t i=0;i<Nspecies;i++) {
      if(i==noIndex) continue;
      if(!hDen2D[i]) continue;

      palette = (TPaletteAxis*)hDen2D[i]->GetListOfFunctions()->FindObject("palette");
      if(palette) {
	Float_t y1b,y2b;
	
	if(j==0)
	  y1b = y1 + j*pvsize;
	else
	  y1b = y1 + j*pvsize + gap/2.0;

	if(j==Ns-1)  
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

    if(j<Ns && opt.Contains("laser")) {
      palette = (TPaletteAxis*)hA2D->GetListOfFunctions()->FindObject("palette");
      
      if(palette) {
	Float_t y1b,y2b;
	
	if(j==0)
	  y1b = y1 + j*pvsize;
	else
	  y1b = y1 + j*pvsize + gap/2.0;

	if(j==Ns-1)  
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

    // 1D plots
    Float_t yaxismin  =  pad[ip]->GetUymin();
    Float_t yaxismax  =  pad[ip]->GetUymin() + 0.3*(pad[ip]->GetUymax() - pad[ip]->GetUymin()) - 0.00;
    //    TGaxis **axis = new TGaxis*[Nspecies];
    TGaxis *axis = NULL;

    if(opt.Contains("den1d")) {
      for(Int_t i=0;i<Nspecies;i++) {

	if(i==noIndex) continue;
	if(!hDen1D[i] || !hDen2D[i]) continue;
	//if(!hDen1D[i] || i!=1) continue;

	Float_t denmin = Min[i];
	Float_t denmax = Max[i];
	if(opt.Contains("logz")) {
	  denmin = TMath::Log10(denmin);
	  denmax = TMath::Log10(denmax);
	}

    
	Float_t slope = (yaxismax - yaxismin)/(denmax - denmin);
    
	for(Int_t j=0;j<hDen1D[i]->GetNbinsX();j++) {
	  Float_t content = hDen1D[i]->GetBinContent(j+1);
	  if(opt.Contains("logz")) content = TMath::Log10(content); 
	  
	  if(content<denmin) 
	    hDen1D[i]->SetBinContent(j+1,yaxismin);
	  else 
	    hDen1D[i]->SetBinContent(j+1,(content - denmin) * slope + yaxismin);
	}

	hDen1D[i]->SetLineWidth(2);
	if(i==0)
	  hDen1D[i]->SetLineColor(kGray+2);
	else 
	  hDen1D[i]->SetLineColor(PGlobals::elecLine);

	hDen1D[i]->Draw("same L");
      }
    } else if(opt.Contains("cur1d")) {

      TLine *lineCurUp = new TLine(xMin,xint,xMax,xint);
      lineCurUp->SetLineColor(kGray+1);
      lineCurUp->SetLineStyle(2);

      TLine *lineCurDo = new TLine(xMin,-xint,xMax,-xint);
      lineCurDo->SetLineColor(kGray+1);
      lineCurDo->SetLineStyle(2);

      if(opt.Contains("cur1dl")) {
	lineCurUp->Draw();
	lineCurDo->Draw();
      }

      Float_t maxCur = -999.0;
      Int_t imaxCur = -1;

      for(Int_t i=0;i<Nspecies;i++) {
	if(i==noIndex) continue;
	if(!hCur1D[i] || !hDen2D[i] || i==0 ) continue;

	hCur1D[i]->ResetStats();

	if(hCur1D[i]->GetMaximum()>maxCur) {
	  maxCur = hCur1D[i]->GetMaximum();
	  imaxCur = i;
	}
	
      }

      Float_t curmin = 0.0;
      Float_t curmax = maxCur;
      // Round for better axis
      curmax = 0.1*TMath::Ceil(10*curmax);
      
      Float_t slope = (yaxismax - yaxismin)/(curmax - curmin);
      Float_t zPos = xMax - (xMax-xMin) * 0.14;
	
      for(Int_t i=0;i<Nspecies;i++) {
	if(i==noIndex) continue;
	if(!hCur1D[i] || !hDen2D[i] || i==0 ) continue;
	hCur1D[i]->ResetStats();

	for(Int_t j=0;j<hCur1D[i]->GetNbinsX();j++) {
	  Float_t content = hCur1D[i]->GetBinContent(j+1);
	  
	  hCur1D[i]->SetBinContent(j+1,(content - curmin) * slope + yaxismin);	
	}
	
	hCur1D[i]->SetLineWidth(2);
	if(i==1)
	  hCur1D[i]->SetLineColor(cBeam);
	else if(i==2)
	  hCur1D[i]->SetLineColor(cBeam1);
	else if(i==3)
	  hCur1D[i]->SetLineColor(cBeam2);
	  
	hCur1D[i]->Draw("hist LF2 same");
	
      }	    

      if(!opt.Contains("noaxis") && maxCur>0.0) {
	
	// Current axis
	axis = new TGaxis(zPos,yaxismin,
			  zPos,yaxismax,
			  curmin,curmax,2,"+LS");
	// axis->SetNdivisions(0);
	axis->SetNdivisions(502);
	
	axis->SetLineWidth(2);
	axis->SetLineColor(kGray+3);//PGlobals::elecLine);
	axis->SetLabelColor(kGray+3);//PGlobals::elecLine);
	axis->SetLabelFont(fonttype);
	axis->SetLabelSize(lzsize-10);
	axis->SetLabelOffset(lyoffset);
	axis->SetTitleColor(kGray+3);//PGlobals::elecLine);
	axis->SetTitleFont(fonttype);
	axis->SetTitleSize(tzsize-14);
	axis->SetTitleOffset(tyoffset);
	axis->SetTickSize(0.03);
	axis->SetTitle(hCur1D[imaxCur]->GetYaxis()->GetTitle());

	axis->ChangeLabel(1,-1,0.);

	axis->CenterTitle();
	// axis->SetMaxDigits(2);
	axis->Draw();
      }

      
    }

         
    if(opt.Contains("ez1d")) {
      // Fit the E1D in the pad:
      Float_t y1 = yMin + yRange/80.;
      Float_t y2 = y1 + yRange/2.8 - 2*yRange/80.;
      Float_t slope = (y2-y1)/(Emax[0]-Emin[0]);
      
      TLine *lineEzero = new TLine(xMin,(0.0-Emin[0])*slope + y1,xMax,(0.0-Emin[0])*slope + y1);
      lineEzero->SetLineColor(lineColor2);
      lineEzero->SetLineStyle(2);
      lineEzero->SetLineWidth(1);
      lineEzero->Draw();

      TH1F *hE1Dclone = (TH1F*) hE1D[0]->Clone("hE1Dclone");

      for(Int_t j=0;j<hE1Dclone->GetNbinsX();j++) {
	hE1Dclone->SetBinContent(j+1,(hE1Dclone->GetBinContent(j+1)-Emin[0])*slope + y1);
      }
      hE1Dclone->SetLineStyle(1);
      hE1Dclone->SetLineWidth(2);
      //hE1Dclone->SetLineWidth(2);

      hE1Dclone->SetLineColor(lineColor2);

      hE1Dclone->Draw("sameL");
    } else if (opt.Contains("et1d")) {
      // Fit the E1D in the E2D pad:
      Float_t HeIon  = 92.75;  // GV/m
      if(!opt.Contains("units") || !n0 ) {
	HeIon  /= ( E0 / eUnit );
      }

      // Fit the E1D in the pad:
      Float_t ETmin =  0.00;//hETotal1D->GetMinimum();
      Float_t ETmax =  hETotal1D->GetMaximum();
      Float_t y1 = yMin + yRange/80.;
      Float_t y2 = y1 + yRange/2.8 - 2*yRange/80.;
      Float_t slope = (y2-y1)/(ETmax-ETmin);
  
      TLine *lineEHe = new TLine(xMin,(HeIon-ETmin)*slope + y1,xMax,(HeIon-ETmin)*slope + y1);
      lineEHe->SetLineColor(lineColor);
      lineEHe->SetLineStyle(2);
      lineEHe->SetLineWidth(1);
      lineEHe->Draw();

      TLine *lineEzero = new TLine(xMin,(0.0-ETmin)*slope + y1,xMax,(0.0-ETmin)*slope + y1);
      lineEzero->SetLineColor(kGray+1);
      lineEzero->SetLineStyle(2);
      lineEzero->SetLineWidth(1);
      lineEzero->Draw();
       
      for(Int_t j=0;j<hETotal1D->GetNbinsX();j++) {
	hETotal1D->SetBinContent(j+1,(hETotal1D->GetBinContent(j+1)-ETmin)*slope + y1);
      }
      hETotal1D->SetLineStyle(1);
      hETotal1D->SetLineWidth(2);
    
      hETotal1D->SetLineColor(lineColor);
    
      hETotal1D->Draw("sameC");
    } 

    if (opt.Contains("curInj")) {

      if(hCur1D[2]) {

	// Fit the 
	Float_t curInjmin =  0.00;
	Float_t curInjmax =  hCur1D[2]->GetMaximum();
	
	// Float_t y1 = yMin + yRange/80.;
	// Float_t y2 = y1 + yRange/2.8 - 2*yRange/80.;
	// Float_t slope = (y2-y1)/(curInjmax-curInjmin);
	Float_t slope = (yaxismax - yaxismin)/(curInjmax - curInjmin);
   
	// TLine *lineCurzero = new TLine(xMin,(0.0-curInjmin)*slope + y1,xMax,(0.0-curInjmin)*slope + y1);
	// lineCurzero->SetLineColor(kGray+1);
	// lineCurzero->SetLineStyle(2);
	// lineCurzero->SetLineWidth(1);
	// lineCurzero->Draw();
       
	for(Int_t j=0;j<hCur1D[2]->GetNbinsX();j++) {
	  Float_t content = hCur1D[2]->GetBinContent(j+1);
	  
	  hCur1D[2]->SetBinContent(j+1,(content - curInjmin) * slope + yaxismin);
	}
	hCur1D[2]->SetLineStyle(1);
	hCur1D[2]->SetLineWidth(2);
    
	hCur1D[2]->SetLineColor(kGray+2);
    
	hCur1D[2]->Draw("sameC");
      } 

    }

  
    if(!opt.Contains("notim"))
      textTime->Draw();
    
    if(!opt.Contains("noden") && opt.Contains("units") && n0) 
      textDen->Draw();


    TBox *lFrame = new TBox(gPad->GetUxmin(), gPad->GetUymin(),
			    gPad->GetUxmax(), gPad->GetUymax());
    lFrame->SetFillStyle(0);
    lFrame->SetLineColor(kBlack);
    lFrame->SetLineWidth(2);
    lFrame->Draw();
    
    pad[ip]->RedrawAxis("g");
    

    ip--;  
    C->cd(0);
  }    


  if(mask & 0x2) {

    pad[ip]->Draw();
    pad[ip]->cd();  // <---------------------------------------------------------- 2nd panel
    // if(opt.Contains("logz")) {
    //   pad[ip]->SetLogz(1);
    // } else {
    //   pad[ip]->SetLogz(0);
    // }

    hFrame[ip]->Draw("col");

    Float_t xFactor = pad[NPad-1]->GetAbsWNDC()/pad[ip]->GetAbsWNDC();
    Float_t yFactor = pad[NPad-1]->GetAbsHNDC()/pad[ip]->GetAbsHNDC();
    // hE2D[0]->GetZaxis()->SetNdivisions(503);
    hE2D[0]->GetZaxis()->SetTitleFont(fonttype);
    hE2D[0]->GetZaxis()->SetTitleOffset(tzoffset);
    hE2D[0]->GetZaxis()->SetTitleSize(tzsize);
    hE2D[0]->GetZaxis()->SetLabelFont(fonttype);
    hE2D[0]->GetZaxis()->SetLabelSize(lzsize);
    hE2D[0]->GetZaxis()->SetLabelOffset(lyoffset);
    hE2D[0]->GetZaxis()->SetTickLength(xFactor*tylength/yFactor);
    
    exField->Draw();
    hE2D[0]->Draw(drawopt + "0");

    if(opt.Contains("1dline")) {
      lineYzero->Draw();
    }


    if(opt.Contains("sline")) {
      if(zStartPlasma>xMin && zStartPlasma<xMax)
	lineStartPlasma->Draw();
      if(zStartNeutral>xMin && zStartNeutral<xMax)
	lineStartNeutral->Draw();
      if(zEndNeutral>xMin && zEndNeutral<xMax)
	lineEndNeutral->Draw();
    }
   
    // Fit the 1D line in the 2D pad:
    TH1F *hE1Dcl2 = (TH1F*) hE1D[0]->Clone("hE1Dcl2");
    
    TLine *lineE = NULL;
    TLine *lineEmaxLuLinear = NULL;
    TLine *lineEmaxUSE = NULL;
    TLine *lineEmaxLo = NULL;
    TLine *lineEmaxLoVert = NULL;
    TLine *lineEmaxBeam = NULL;
    {
      Float_t rightmin = Emin[0];
      Float_t rightmax = Emax[0];
      Float_t slope = (yMax-yMin)/(rightmax-rightmin);
    
      for(Int_t j=0;j<hE1Dcl2->GetNbinsX();j++) {
	hE1Dcl2->SetBinContent(j+1,slope*(hE1Dcl2->GetBinContent(j+1)-rightmin)+yMin);
      }
    
      Double_t sl = EzMaxLo/radius;
      Double_t zmatch = (EzMaxBeamLo/sl)+EzCross[0];

      Double_t EzAdapted  = slope*(EzMaxLo-rightmin)+yMin;
      Double_t EzAdapted2 = slope*(EzMaxBeamLo-rightmin)+yMin;
      lineE = new TLine(EzCross[0]-radius,-EzAdapted,zmatch,EzAdapted2);

      Double_t linelength = 1;
      if(opt.Contains("units")) linelength *= skindepth / spaUnit;
 
      EzAdapted = slope*(EzMaxLuLinear-rightmin)+yMin;
      lineEmaxLuLinear = new TLine(EzCross[0]-radius,-EzAdapted,EzCross[0]-radius-linelength,-EzAdapted);

      EzAdapted = slope*(EzMaxUse-rightmin)+yMin;
      lineEmaxUSE = new TLine(EzCross[0]-radius,EzAdapted,EzCross[0],EzAdapted);

      EzAdapted = slope*(EzMaxLo-rightmin)+yMin;
      lineEmaxLo = new TLine(EzCross[0]-radius,-EzAdapted,EzCross[0]-radius-linelength,-EzAdapted);
      lineEmaxLoVert = new TLine(EzCross[0]-radius,-EzAdapted,EzCross[0]-radius,0.0);

      EzAdapted = slope*(EzMaxBeamLo-rightmin)+yMin;
      lineEmaxBeam = new TLine(zmatch,EzAdapted,-zmatch,EzAdapted);
      
      if(opt.Contains("mark")) {
	Float_t e0 = slope*(0-rightmin)+yMin;

	Float_t zF = hE1Dcl2->GetBinCenter(binFcross);
	Float_t EF = hE1Dcl2->GetBinContent(binFcross);

	TLine *lineEF = new TLine(zF,0,zF,EF);
	lineEF->SetLineColor(kRed);
	lineEF->SetLineStyle(2);
	lineEF->Draw();

	TMarker *markEF = new TMarker(zF,EF,20);
	markEF->SetMarkerColor(kRed);
	markEF->SetMarkerSize(1.4);
	markEF->Draw();

	Float_t zB = zEbMax;
	Float_t EB = slope*(EzExtr[0]-rightmin)+yMin;
      
	TLine *lineEB = new TLine(zB,0,zB,EB);
	lineEB->SetLineColor(kRed);
	lineEB->SetLineStyle(2);
	lineEB->Draw();

	TMarker *markEB = new TMarker(zB,EB,20);
	markEB->SetMarkerColor(kRed);
	markEB->SetMarkerSize(1.4);
	markEB->Draw();

	Float_t z0 = EzCross[0];
      
	TMarker *markE0 = new TMarker(z0,e0,20);
	markE0->SetMarkerColor(kRed);
	markE0->SetMarkerSize(1.4);
	markE0->Draw();

	Float_t zI = hE1Dcl2->GetBinCenter(binPotValueIni);
	Float_t EI = hE1Dcl2->GetBinContent(binPotValueIni);
	TLine *lineEI = new TLine(zI,e0,zI,EI);
	lineEI->SetLineColor(kRed);
	lineEI->SetLineStyle(2);
	lineEI->Draw();

	TMarker *markEI = new TMarker(zI,EI,20);
	markEI->SetMarkerColor(kRed);
	markEI->SetMarkerSize(1.4);
	markEI->Draw();
      }


    }
    
    if(opt.Contains("zero")) {
      TLine *lineZero = new TLine(xMin,0,xMax,0);
      lineZero->SetLineColor(lineColor);
      lineZero->SetLineStyle(2);
      lineZero->Draw();
    }

    // if(opt.Contains("off")) {
    //   TLine *lineOff = new TLine(xMin,xoff,xMax,xoff);
    //   lineOff->SetLineColor(lineColor);
    //   lineOff->SetLineStyle(2);
    //   lineOff->Draw();
    // }


    if(opt.Contains("bopar")) {

      lineE->SetLineColor(kGray+3);
      lineE->SetLineStyle(3);
      lineE->Draw();    
      
      lineEmaxLuLinear->SetLineColor(kGray+3);
      lineEmaxLuLinear->SetLineStyle(4);
      //lineEmaxLuLinear->Draw();    
      
      // lineEmaxUSE->SetLineColor(kGray+3);
      // lineEmaxUSE->SetLineStyle(4);
      // lineEmaxUSE->Draw();    
      
      lineEmaxLo->SetLineColor(kGray+3);
      lineEmaxLo->SetLineStyle(3);
      lineEmaxLo->Draw();    
      
      lineEmaxLoVert->SetLineColor(kGray+3);
      lineEmaxLoVert->SetLineStyle(3);
      lineEmaxLoVert->Draw();    
      
      lineEmaxBeam->SetLineColor(kGray+3);
      lineEmaxBeam->SetLineStyle(3);
      lineEmaxBeam->Draw();    


      TEllipse *ellipLo = new TEllipse(EzCross[0],0.0,radius);	
      ellipLo->SetFillStyle(0);
      ellipLo->SetLineStyle(3);
      ellipLo->SetLineColor(kGray+3);
      ellipLo->Draw();
	
      TEllipse *ellipLu = new TEllipse(EzCross[0],0.0,radiusLu);
      ellipLu->SetFillStyle(0);
      ellipLu->SetLineStyle(4);
      ellipLu->SetLineColor(kGray+3);
      ellipLu->Draw();


    }

    hE1Dcl2->SetLineStyle(1);
    hE1Dcl2->SetLineWidth(2);
    hE1Dcl2->SetLineColor(lineColor);

    if(!opt.Contains("no1d"))
      hE1Dcl2->Draw("sameL");    
    
    pad[ip]->Update();
  
    y1 = pad[ip]->GetBottomMargin();
    y2 = 1 - pad[ip]->GetTopMargin();
    x1 = pad[ip]->GetLeftMargin();
    x2 = 1 - pad[ip]->GetRightMargin();
  
    palette = (TPaletteAxis*) hE2D[0]->GetListOfFunctions()->FindObject("palette");  
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
      pFrame->SetShadowColor(0);
      pFrame->Draw();

    }

    TBox *lFrame = new TBox(gPad->GetUxmin(), gPad->GetUymin(),
			    gPad->GetUxmax(), gPad->GetUymax());
    lFrame->SetFillStyle(0);
    lFrame->SetLineColor(kBlack);
    lFrame->SetLineWidth(2);
    lFrame->Draw();
    
    pad[ip]->RedrawAxis("g"); 

    if(opt.Contains("blines")) {
      TLine *lineEbunch = NULL;
      TLine *lineESlope = NULL;
      TLine *lineZbunch = NULL;
      if(zInjBeam>-999) { 
	Float_t slope = (yMax-yMin)/(Emax[0]-Emin[0]);

	lineEbunch = new TLine(xMin,EInjBeam*slope,xMax,EInjBeam*slope);
	lineEbunch->SetLineColor(kGray+2);
	lineEbunch->SetLineStyle(2);
	lineEbunch->SetLineWidth(1);
	lineEbunch->Draw();
      
	lineZbunch = new TLine(zInjBeam,yMin,zInjBeam,yMax);
	lineZbunch->SetLineColor(kGray+2);
	lineZbunch->SetLineStyle(2);
	lineZbunch->SetLineWidth(1);
	lineZbunch->Draw();
      
	Float_t z1 = zInjBeam - 10;
	Float_t z2 = zInjBeam + 10;
	Float_t E1 = EInjBeam + ESlopeInjBeam * (z1-zInjBeam);
	Float_t E2 = EInjBeam + ESlopeInjBeam * (z2-zInjBeam);
	lineESlope = new TLine(z1,E1*slope,z2,E2*slope);
	lineESlope->SetLineColor(lineColor);
	lineESlope->SetLineStyle(1);
	lineESlope->SetLineWidth(1);
	//lineESlope->Draw();    
      
      }
    }

  
  
    if(opt.Contains("1dline")) {
      lineYzero->Draw();
    }

    if(opt.Contains("sline")) {
      if(zStartPlasma>xMin && zStartPlasma<xMax)
	lineStartPlasma->Draw();
      if(zStartNeutral>xMin && zStartNeutral<xMax)
	lineStartNeutral->Draw();
      if(zEndNeutral>xMin && zEndNeutral<xMax)
	lineEndNeutral->Draw();
    }
 
    ip--;  
    C->cd(0);
  }



  if(mask & 0x4) {

    pad[ip]->Draw();
    pad[ip]->cd(); // <---------------------------------------------------------- 3rd panel
    // if(opt.Contains("logz")) {
    //   pad[ip]->SetLogz(1);
    // } else {
    //   pad[ip]->SetLogz(0);
    // }

    hFrame[ip]->Draw("col");

    Float_t xFactor = pad[NPad-1]->GetAbsWNDC()/pad[ip]->GetAbsWNDC();
    Float_t yFactor = pad[NPad-1]->GetAbsHNDC()/pad[ip]->GetAbsHNDC();
    // hE2D[1]->GetZaxis()->SetNdivisions(503);
    hE2D[1]->GetZaxis()->SetTitleFont(fonttype);
    hE2D[1]->GetZaxis()->SetTitleOffset(tzoffset);
    hE2D[1]->GetZaxis()->SetTitleSize(tzsize);
    hE2D[1]->GetZaxis()->SetLabelFont(fonttype);
    hE2D[1]->GetZaxis()->SetLabelSize(lzsize);
    hE2D[1]->GetZaxis()->SetLabelOffset(lyoffset);
    hE2D[1]->GetZaxis()->SetTickLength(xFactor*tylength/yFactor);
  
    exField->Draw();
    hE2D[1]->Draw("colz0 same");

    if(opt.Contains("bopar")) {
      
      TEllipse *ellipLo = new TEllipse(EzCross[0],0.0,radius);	
      ellipLo->SetFillStyle(0);
      ellipLo->SetLineStyle(3);
      ellipLo->SetLineColor(kGray+3);
      ellipLo->Draw();
	
      TEllipse *ellipLu = new TEllipse(EzCross[0],0.0,radiusLu);
      ellipLu->SetFillStyle(0);
      ellipLu->SetLineStyle(4);
      ellipLu->SetLineColor(kGray+3);
      ellipLu->Draw();

    }


    if(opt.Contains("1dline")) {
      lineYzero->Draw();
    }

    if(opt.Contains("sline")) {
      if(zStartPlasma>xMin && zStartPlasma<xMax)
	lineStartPlasma->Draw();
      if(zStartNeutral>xMin && zStartNeutral<xMax)
	lineStartNeutral->Draw();
      if(zEndNeutral>xMin && zEndNeutral<xMax)
	lineEndNeutral->Draw();
    }
   
    // Fit the 1D line in the 2D pad:
    {
      Float_t rightmin = Emin[1];
      Float_t rightmax = Emax[1];
      Float_t slope = (yMax-yMin)/(rightmax-rightmin);
    
      for(Int_t j=0;j<hE1D[1]->GetNbinsX();j++) {
	hE1D[1]->SetBinContent(j+1,slope*(hE1D[1]->GetBinContent(j+1)-rightmin)+yMin);
      }
    }


    if(opt.Contains("zero")) {
      TLine *lineZero = new TLine(xMin,0,xMax,0);
      //lineZero->SetLineColor(lineColor);
      lineZero->SetLineColor(kGray);
      lineZero->SetLineStyle(2);
      lineZero->Draw();
    }

    if(opt.Contains("off")) {
      TLine *lineOff = new TLine(xMin,xoff,xMax,xoff);
      lineOff->SetLineColor(lineColor);
      lineOff->SetLineStyle(2);
      lineOff->Draw();
    }


    hE1D[1]->SetLineStyle(1);
    hE1D[1]->SetLineWidth(2);
    hE1D[1]->SetLineColor(lineColor);

    if(!opt.Contains("no1d"))
      hE1D[1]->Draw("sameL");
  
    pad[ip]->Update();
  
    y1 = pad[ip]->GetBottomMargin();
    y2 = 1 - pad[ip]->GetTopMargin();
    x1 = pad[ip]->GetLeftMargin();
    x2 = 1 - pad[ip]->GetRightMargin();
  
    palette = (TPaletteAxis*) hE2D[1]->GetListOfFunctions()->FindObject("palette");  
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
      pFrame->SetShadowColor(0);
      pFrame->Draw();

    }

    TBox *lFrame = new TBox(gPad->GetUxmin(), gPad->GetUymin(),
			    gPad->GetUxmax(), gPad->GetUymax());
    lFrame->SetFillStyle(0);
    lFrame->SetLineColor(PGlobals::frameColor);
    lFrame->SetLineWidth(PGlobals::frameWidth);
    lFrame->Draw();
    
    pad[ip]->RedrawAxis(); 
  
  
    if(opt.Contains("1dline")) {
      lineYzero->Draw();
    }

    if(opt.Contains("sline")) {
      if(zStartPlasma>xMin && zStartPlasma<xMax)
	lineStartPlasma->Draw();
      if(zStartNeutral>xMin && zStartNeutral<xMax)
	lineStartNeutral->Draw();
      if(zEndNeutral>xMin && zEndNeutral<xMax)
	lineEndNeutral->Draw();
    }
 
    ip--;  
    C->cd(0);
  }
 

  if(mask & 0x8) {

    pad[ip]->Draw();
    pad[ip]->cd(); // <---------------------------------------------------------- 4th panel
    // if(opt.Contains("logz")) {
    //   pad[ip]->SetLogz(1);
    // } else {
    //   pad[ip]->SetLogz(0);
    // }

    hFrame[ip]->Draw("col");

    Float_t xFactor = pad[NPad-1]->GetAbsWNDC()/pad[ip]->GetAbsWNDC();
    Float_t yFactor = pad[NPad-1]->GetAbsHNDC()/pad[ip]->GetAbsHNDC();
    // hE2D[2]->GetZaxis()->SetNdivisions(503);
    hE2D[2]->GetZaxis()->SetTitleFont(fonttype);
    hE2D[2]->GetZaxis()->SetTitleOffset(tzoffset);
    hE2D[2]->GetZaxis()->SetTitleSize(tzsize);
    hE2D[2]->GetZaxis()->SetLabelFont(fonttype);
    hE2D[2]->GetZaxis()->SetLabelSize(lzsize);
    hE2D[2]->GetZaxis()->SetLabelOffset(lyoffset);
    hE2D[2]->GetZaxis()->SetTickLength(xFactor*tylength/yFactor);  
  
    exField->Draw();
    hE2D[2]->Draw("colz0 same");

    if(opt.Contains("1dline")) {
      lineYzero->Draw();
    }

    if(opt.Contains("sline")) {
      if(zStartPlasma>xMin && zStartPlasma<xMax)
	lineStartPlasma->Draw();
      if(zStartNeutral>xMin && zStartNeutral<xMax)
	lineStartNeutral->Draw();
      if(zEndNeutral>xMin && zEndNeutral<xMax)
	lineEndNeutral->Draw();
    }
   
    // Fit the 1D line in the 2D pad:
    {
      Float_t rightmin = Emin[2];
      Float_t rightmax = Emax[2];
      Float_t slope = (yMax-yMin)/(rightmax-rightmin);
    
      for(Int_t j=0;j<hE1D[2]->GetNbinsX();j++) {
	hE1D[2]->SetBinContent(j+1,slope*(hE1D[2]->GetBinContent(j+1)-rightmin)+yMin);
      }
    
    }


    if(opt.Contains("zero")) {
      TLine *lineZero = new TLine(xMin,0,xMax,0);
      //lineZero->SetLineColor(lineColor);
      lineZero->SetLineColor(kGray);
      lineZero->SetLineStyle(2);
      lineZero->Draw();
    }

    if(opt.Contains("off")) {
      TLine *lineOff = new TLine(xMin,xoff,xMax,xoff);
      lineOff->SetLineColor(lineColor);
      lineOff->SetLineStyle(2);
      lineOff->Draw();
    }

    hE1D[2]->SetLineStyle(1);
    hE1D[2]->SetLineWidth(2);
    hE1D[2]->SetLineColor(lineColor);

    if(!opt.Contains("no1d"))
      hE1D[2]->Draw("sameL");
  
    pad[ip]->Update();
  
    y1 = pad[ip]->GetBottomMargin();
    y2 = 1 - pad[ip]->GetTopMargin();
    x1 = pad[ip]->GetLeftMargin();
    x2 = 1 - pad[ip]->GetRightMargin();
  
    palette = (TPaletteAxis*) hE2D[2]->GetListOfFunctions()->FindObject("palette");  
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
      pFrame->SetShadowColor(0);
      pFrame->Draw();

    }

    TBox *lFrame = new TBox(gPad->GetUxmin(), gPad->GetUymin(),
			    gPad->GetUxmax(), gPad->GetUymax());
    lFrame->SetFillStyle(0);
    lFrame->SetLineColor(PGlobals::frameColor);
    lFrame->SetLineWidth(PGlobals::frameWidth);
    lFrame->Draw();

    
    pad[ip]->RedrawAxis(); 
  
  
    if(opt.Contains("1dline")) {
      lineYzero->Draw();
    }

    if(opt.Contains("sline")) {
      if(zStartPlasma>xMin && zStartPlasma<xMax)
	lineStartPlasma->Draw();
      if(zStartNeutral>xMin && zStartNeutral<xMax)
	lineStartNeutral->Draw();
      if(zEndNeutral>xMin && zEndNeutral<xMax)
	lineEndNeutral->Draw();
    }
 
    ip--;
    C->cd(0);
  }

  if( (mask & 0x10) && (hFocus2D) ) {

    pad[ip]->Draw();
    pad[ip]->cd(); // <---------------------------------------------------------- 5th panel
    // if(opt.Contains("logz")) {
    //   pad[ip]->SetLogz(1);
    // } else {
    //   pad[ip]->SetLogz(0);
    // }

    hFrame[ip]->Draw("col");

    Float_t xFactor = pad[NPad-1]->GetAbsWNDC()/pad[ip]->GetAbsWNDC();
    Float_t yFactor = pad[NPad-1]->GetAbsHNDC()/pad[ip]->GetAbsHNDC();
    //   hFocus2D->GetZaxis()->SetNdivisions(503);
    hFocus2D->GetZaxis()->SetTitleFont(fonttype);
    hFocus2D->GetZaxis()->SetTitleOffset(tzoffset);
    hFocus2D->GetZaxis()->SetTitleSize(tzsize);
    hFocus2D->GetZaxis()->SetLabelFont(fonttype);
    hFocus2D->GetZaxis()->SetLabelSize(lzsize);
    hFocus2D->GetZaxis()->SetLabelOffset(lyoffset);
    hFocus2D->GetZaxis()->SetTickLength(xFactor*tylength/yFactor);

    //    hFocus2D->GetZaxis()->SetTitle("E_{z} [GV/m]");
    
    exField->Draw();
    hFocus2D->Draw(drawopt + "0");


    if(opt.Contains("bopar")) {
      
      TEllipse *ellipLo = new TEllipse(EzCross[0],0.0,radius);	
      ellipLo->SetFillStyle(0);
      ellipLo->SetLineStyle(3);
      ellipLo->SetLineColor(kGray+3);
      ellipLo->Draw();
	
      TEllipse *ellipLu = new TEllipse(EzCross[0],0.0,radiusLu);
      ellipLu->SetFillStyle(0);
      ellipLu->SetLineStyle(4);
      ellipLu->SetLineColor(kGray+3);
      ellipLu->Draw();

    }


    if(opt.Contains("1dline")) {
      lineYzero->Draw();
    }

    if(opt.Contains("sline")) {
      if(zStartPlasma>xMin && zStartPlasma<xMax)
	lineStartPlasma->Draw();
      if(zStartNeutral>xMin && zStartNeutral<xMax)
	lineStartNeutral->Draw();
      if(zEndNeutral>xMin && zEndNeutral<xMax)
	lineEndNeutral->Draw();
    }
  
    if(opt.Contains("zero")) {
      TLine *lineZero = new TLine(xMin,0,xMax,0);
      //lineZero->SetLineColor(lineColor);
      lineZero->SetLineColor(kGray);
      lineZero->SetLineStyle(2);
      lineZero->Draw();
    }


    if(opt.Contains("off") && !opt.Contains("no1d")) {
      TLine *lineOff = new TLine(xMin,xoff,xMax,xoff);
      lineOff->SetLineColor(lineColor);
      lineOff->SetLineStyle(2);
      lineOff->Draw();
    }

 
    // Fit the 1D line in the 2D pad:
    {
      Float_t rightmin = Fmin;
      Float_t rightmax = Fmax;
      Float_t slope = (yMax-yMin)/(rightmax-rightmin);
    
      for(Int_t j=0;j<hFocus1D->GetNbinsX();j++) {
	hFocus1D->SetBinContent(j+1,slope*(hFocus1D->GetBinContent(j+1)-rightmin)+yMin);
      }
      
      if(opt.Contains("mark")) {
	
	Float_t zF = hFocus1D->GetBinCenter(binFcross);
	Float_t FF = slope*(0.0-rightmin)+yMin;
	
	TMarker *markFF = new TMarker(zF,FF,20);
	markFF->SetMarkerColor(kRed);
	markFF->SetMarkerSize(1.4);
	markFF->Draw();
	
      }
      
    }
  
    hFocus1D->SetLineStyle(1);
    hFocus1D->SetLineWidth(2);
    hFocus1D->SetLineColor(lineColor);

    if(!opt.Contains("no1d"))
      hFocus1D->Draw("sameL");
  
    pad[ip]->Update();
  
    y1 = pad[ip]->GetBottomMargin();
    y2 = 1 - pad[ip]->GetTopMargin();
    x1 = pad[ip]->GetLeftMargin();
    x2 = 1 - pad[ip]->GetRightMargin();
  
    palette = (TPaletteAxis*) hFocus2D->GetListOfFunctions()->FindObject("palette");  
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
      pFrame->SetShadowColor(0);
      pFrame->Draw();

    }

    TBox *lFrame = new TBox(gPad->GetUxmin(), gPad->GetUymin(),
			    gPad->GetUxmax(), gPad->GetUymax());
    lFrame->SetFillStyle(0);
    lFrame->SetLineColor(PGlobals::frameColor);
    lFrame->SetLineWidth(PGlobals::frameWidth);
    lFrame->Draw();

    pad[ip]->RedrawAxis(); 
  
  
    if(opt.Contains("1dline")) {
      lineYzero->Draw();
    }

    if(opt.Contains("sline")) {
      if(zStartPlasma>xMin && zStartPlasma<xMax)
	lineStartPlasma->Draw();
      if(zStartNeutral>xMin && zStartNeutral<xMax)
	lineStartNeutral->Draw();
      if(zEndNeutral>xMin && zEndNeutral<xMax)
	lineEndNeutral->Draw();
    }
 
    ip--;
    C->cd(0);
  }

  


  if(mask & 0x20) {
    
    pad[ip]->Draw();
    pad[ip]->cd(); // <--------------------------------------------------- 6th panel
    // if(opt.Contains("logz")) {
    //   pad[ip]->SetLogz(1);
    // } else {
    //   pad[ip]->SetLogz(0);
    // }

    hFrame[ip]->Draw("col");

    Float_t xFactor = pad[NPad-1]->GetAbsWNDC()/pad[ip]->GetAbsWNDC();
    Float_t yFactor = pad[NPad-1]->GetAbsHNDC()/pad[ip]->GetAbsHNDC();
    // hETotal2D->GetZaxis()->SetNdivisions(503);
    hETotal2D->GetZaxis()->SetTitleFont(fonttype);
    hETotal2D->GetZaxis()->SetTitleOffset(tzoffset);
    hETotal2D->GetZaxis()->SetTitleSize(tzsize);
    hETotal2D->GetZaxis()->SetLabelFont(fonttype);
    hETotal2D->GetZaxis()->SetLabelSize(lzsize);
    hETotal2D->GetZaxis()->SetLabelOffset(lyoffset);
    hETotal2D->GetZaxis()->SetTickLength(xFactor*tylength/yFactor);


    exFieldT->Draw();
  
    hETotal2D->Draw(drawopt);

    // ADK contours  
    // for(Int_t i=0;i<graphsI2D->GetEntriesFast();i++) {
    //   TGraph *gr = (TGraph*) graphsI2D->At(i);
    //   if(!gr) continue;
    //   gr->Draw("C");
    // }

    
    // Fit the E1D in the E2D pad:
    Float_t rightmin = ETmin;
    Float_t rightmax = ETmax;
    Float_t slope = (yMax-yMin)/(rightmax-rightmin);
    for(Int_t j=0;j<hETotal1D->GetNbinsX();j++) {
      hETotal1D->SetBinContent(j+1,slope*(hETotal1D->GetBinContent(j+1)-rightmin)+yMin);
    }


    // if(opt.Contains("zero")) {
    //   TLine *lineZero = new TLine(xMin,0,xMax,0);
    //   lineZero->SetLineColor(lineColor);
    //   lineZero->SetLineStyle(2);
    //   lineZero->Draw();
    // }

    if(opt.Contains("off")) {
      TLine *lineOff = new TLine(xMin,xoff,xMax,xoff);
      lineOff->SetLineColor(lineColor);
      lineOff->SetLineStyle(2);
      lineOff->Draw();
    }

    
    if(opt.Contains("ionth")) {
      // const Int_t Nat = 6;
      // // char atNames[Nat][8] = {"H", "He","He^{+}", "Ne", "Ne+","Custom"};
      // Float_t IonTh[Nat] = {33.79,92.77,234.97,74.42,143.69,500.59};
      // TLine *lineTh[Nat];
      // Int_t atColor[Nat] = {kAzure+1,kOrange+8,kOrange+8,kSpring-1,kSpring-1,kCyan};
      // Int_t atStyle[Nat] = {1,1,2,1,2,2};
      const Int_t Nat = 5;
      // char atNames[Nat][8] = {"H", "He","He^{+}", "Ne", "Ne+"};
      Float_t IonTh[Nat] = {33.79,92.77,234.97,74.42,143.69};
      TLine *lineTh[Nat];
      Int_t atColor[Nat] = {kAzure+1,kOrange+8,kOrange+8,kSpring-1,kSpring-1};
      Int_t atStyle[Nat] = {1,1,2,1,2};
      
      for(Int_t i=0;i<Nat;i++) {
	if(!opt.Contains("units") || !n0 ) 
	  IonTh[i] /= ( E0 / eUnit );
	
	lineTh[i] = new TLine(xMin,slope*(IonTh[i]-rightmin)+yMin,
			      xMax,slope*(IonTh[i]-rightmin)+yMin);

	lineTh[i]->SetLineColor(atColor[i]);
	lineTh[i]->SetLineStyle(atStyle[i]);
	lineTh[i]->Draw();


      }
    }

    if(opt.Contains("mark")) {

      Float_t zI = hETotal1D->GetBinCenter(binPotValueIni);
      Float_t EI = hETotal1D->GetBinContent(binPotValueIni);
      Float_t e0 = slope*(0-rightmin)+yMin;
 
      TLine *lineEI = new TLine(zI,e0,zI,EI);
      lineEI->SetLineColor(kRed);
      lineEI->SetLineStyle(2);
      lineEI->Draw();

      TMarker *markEI = new TMarker(zI,EI,20);
      markEI->SetMarkerColor(kRed);
      markEI->SetMarkerSize(1.4);
      markEI->Draw();

      Float_t zB = zEbMax;
      Float_t EB = hETotal1D->GetBinContent(binEbMax);
      
      TLine *lineEB = new TLine(zB,e0,zB,EB);
      lineEB->SetLineColor(kRed);
      lineEB->SetLineStyle(2);
      lineEB->Draw();

      TMarker *markEB = new TMarker(zB,EB,20);
      markEB->SetMarkerColor(kRed);
      markEB->SetMarkerSize(1.4);
      markEB->Draw();

    }

    if(opt.Contains("bopar")) {
      
      TEllipse *ellipLo = new TEllipse(EzCross[0],0.0,radius);	
      ellipLo->SetFillStyle(0);
      ellipLo->SetLineStyle(3);
      ellipLo->SetLineColor(kGray+3);
      ellipLo->Draw();
	
      TEllipse *ellipLu = new TEllipse(EzCross[0],0.0,radiusLu);
      ellipLu->SetFillStyle(0);
      ellipLu->SetLineStyle(4);
      ellipLu->SetLineColor(kGray+3);
      ellipLu->Draw();

    }

  
    hETotal1D->SetLineStyle(1);
    hETotal1D->SetLineWidth(2);
    hETotal1D->SetLineColor(lineColor);

    hETotal1D->Draw("sameL");
    
    if(opt.Contains("1dline")) {
      lineYzero->Draw();
    }

    if(opt.Contains("sline")) {
      if(zStartPlasma>xMin && zStartPlasma<xMax)
	lineStartPlasma->Draw();
      if(zStartNeutral>xMin && zStartNeutral<xMax)
	lineStartNeutral->Draw();
      if(zEndNeutral>xMin && zEndNeutral<xMax)
	lineEndNeutral->Draw();
    }
  
    
    if(opt.Contains("cont")) {
      // PSI MAIN contour  
      for(Int_t i=0;i<graphsV2D_main.GetEntriesFast();i++) {
	TGraph *gr = (TGraph*) graphsV2D_main.At(i);
	if(!gr) continue;
	// gr->Draw("C");
      }
    }
    
  
    pad[ip]->Update();
  
    y1 = pad[ip]->GetBottomMargin();
    y2 = 1 - pad[ip]->GetTopMargin();
    x1 = pad[ip]->GetLeftMargin();
    x2 = 1 - pad[ip]->GetRightMargin();
  
    palette = (TPaletteAxis*) hETotal2D->GetListOfFunctions()->FindObject("palette");  
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
      pFrame->SetShadowColor(0);
      pFrame->Draw();

    }

    TBox *lFrame = new TBox(gPad->GetUxmin(), gPad->GetUymin(),
			    gPad->GetUxmax(), gPad->GetUymax());
    lFrame->SetFillStyle(0);
    lFrame->SetLineColor(PGlobals::frameColor);
    lFrame->SetLineWidth(PGlobals::frameWidth);
    lFrame->Draw();

    pad[ip]->RedrawAxis(); 

    ip--;
    C->cd(0);

  }

  if(mask & 0x40) {

    pad[ip]->Draw();
    pad[ip]->cd(); // <--------------------------------------------------------- 7th panel
 
    hFrame[ip]->Draw("col");

    Float_t xFactor = pad[NPad-1]->GetAbsWNDC()/pad[ip]->GetAbsWNDC();
    Float_t yFactor = pad[NPad-1]->GetAbsHNDC()/pad[ip]->GetAbsHNDC();
    // hV2D->GetZaxis()->SetNdivisions(503);
    hV2D->GetZaxis()->SetTitleFont(fonttype);
    hV2D->GetZaxis()->SetTitleOffset(tzoffset);
    hV2D->GetZaxis()->SetTitleSize(tzsize);
    hV2D->GetZaxis()->SetLabelFont(fonttype);
    hV2D->GetZaxis()->SetLabelSize(lzsize);
    hV2D->GetZaxis()->SetLabelOffset(lyoffset);
    hV2D->GetZaxis()->SetTickLength(xFactor*tylength/yFactor);

 
    exPot->Draw();
  
    hV2D->GetZaxis()->SetRangeUser(Vmin,Vmax);
    hV2D->Draw("colz0 same");
    
    if(opt.Contains("cont")) {
      for(Int_t i=0;i<graphsV2D.GetEntriesFast();i++) {
	TGraph *gr = (TGraph*) graphsV2D.At(i);
	if(!gr) continue;
	gr->Draw("C");
      }
      
      // PSI MAIN contour  
      // for(Int_t i=0;i<graphsV2D_main.GetEntriesFast();i++) {
      //   TGraph *gr = (TGraph*) graphsV2D_main.At(i);
      //   if(!gr) continue;
      //   gr->Draw("C");
      // }
    }
    
    // Fit the V1D in the V2D pad:
    TH1F *hV1Dcl = (TH1F*) hV1D->Clone("hV1Dcl");

    TLine *lineT = NULL;
    {
      Float_t rightmin = Vmin;
      Float_t rightmax = Vmax;
      Float_t slope = (yMax-yMin)/(rightmax-rightmin);
      Float_t psi0 = slope*(0.0-rightmin)+yMin;

      for(Int_t j=0;j<hV1Dcl->GetNbinsX();j++) {
	hV1Dcl->SetBinContent(j+1,slope*(hV1Dcl->GetBinContent(j+1)-rightmin)+yMin);
      }

      if(opt.Contains("zero")) {
	lineT = new TLine(xMin,psi0,xMax,psi0);

	
	lineT->SetLineColor(lineColor);
	lineT->SetLineStyle(2);
	lineT->Draw();
      }
    
    

      if(opt.Contains("mark")) {
	Float_t zI = hV1Dcl->GetBinCenter(binPotValueIni);
	Float_t psiI = hV1Dcl->GetBinContent(binPotValueIni);
	Float_t zF = hV1Dcl->GetBinCenter(binFcross);
	Float_t psiF = hV1Dcl->GetBinContent(binFcross);
	
	TLine *linePsiI = new TLine(zI,psiF,zI,psiI);
	linePsiI->SetLineColor(kRed);
	linePsiI->SetLineStyle(2);
	linePsiI->Draw();
      
	TMarker *markPsiI = new TMarker(zI,psiI,20);
	markPsiI->SetMarkerColor(kRed);
	markPsiI->SetMarkerSize(1.4);
	markPsiI->Draw();

		
	TLine *linePsiF = new TLine(zF,psiF,zF,psiF);
	linePsiF->SetLineColor(kRed);
	linePsiF->SetLineStyle(2);
	linePsiF->Draw();
      
	TMarker *markPsiF = new TMarker(zF,psiF,20);
	markPsiF->SetMarkerColor(kRed);
	markPsiF->SetMarkerSize(1.4);
	markPsiF->Draw();
	
      }
    }

    hV1Dcl->SetLineStyle(1);
    hV1Dcl->SetLineWidth(2);
    //  hV1Dcl->SetLineColor(PGlobals::elecLine);
    hV1Dcl->SetLineColor(lineColor);

    if(!opt.Contains("no1d"))
      hV1Dcl->Draw("sameL");
    
    // ADK MAIN contour  
    // for(Int_t i=0;i<graphsI2D_main->GetEntriesFast();i++) {
    //   TGraph *gr = (TGraph*) graphsI2D_main->At(i);
    //   if(!gr) continue;
    //   gr->Draw("C");
    // }

    if(opt.Contains("bopar")) {
      TEllipse *ellipTrap = new TEllipse(EzCross[0],0.0,RadiusPsiLo);           
      ellipTrap->SetFillStyle(0);
      ellipTrap->SetLineStyle(3);
      ellipTrap->SetLineColor(kWhite);
      //ellipTrap->Draw();

      TEllipse *ellipLo = new TEllipse(EzCross[0],0.0,radius);	
      ellipLo->SetFillStyle(0);
      ellipLo->SetLineStyle(3);
      ellipLo->SetLineColor(kGray+3);
      ellipLo->Draw();

      TEllipse *ellipLu = new TEllipse(EzCross[0],0.0,radiusLu);
      ellipLu->SetFillStyle(0);
      ellipLu->SetLineStyle(4);
      ellipLu->SetLineColor(kGray+3);
      ellipLu->Draw();


    }

    if(opt.Contains("1dline")) {
      lineYzero->Draw();
    }

    if(opt.Contains("zero")) {
      TLine *lineZero = new TLine(xMin,0,xMax,0);
      lineZero->SetLineColor(lineColor);
      lineZero->SetLineStyle(2);
      lineZero->Draw();
    }
    

    if(opt.Contains("sline")) {
      if(zStartPlasma>xMin && zStartPlasma<xMax)
	lineStartPlasma->Draw();
      if(zStartNeutral>xMin && zStartNeutral<xMax)
	lineStartNeutral->Draw();
      if(zEndNeutral>xMin && zEndNeutral<xMax)
	lineEndNeutral->Draw();
    }
  
    
    pad[ip]->Update();

    y1 = pad[ip]->GetBottomMargin();
    y2 = 1 - pad[ip]->GetTopMargin();
    x1 = pad[ip]->GetLeftMargin();
    x2 = 1 - pad[ip]->GetRightMargin();
 
    palette = (TPaletteAxis*)hV2D->GetListOfFunctions()->FindObject("palette");
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
      pFrame->SetShadowColor(0);
      pFrame->Draw();

    }

    TBox *lFrame = new TBox(gPad->GetUxmin(), gPad->GetUymin(),
			    gPad->GetUxmax(), gPad->GetUymax());
    lFrame->SetFillStyle(0);
    lFrame->SetLineColor(PGlobals::frameColor);
    lFrame->SetLineWidth(PGlobals::frameWidth);
    lFrame->Draw();

    pad[ip]->RedrawAxis(); 
  
    ip--;
    C->cd();
  }


  if(mask & 0x80) {
    pad[ip]->Draw();
    pad[ip]->cd(); // <---------------------------------------------------------- 8th panel
    // if(opt.Contains("logz")) {
    //   pad[ip]->SetLogz(1);
    // } else {
    //   pad[ip]->SetLogz(0);
    // }

    hFrame[ip]->Draw("col");

    Float_t xFactor = pad[NPad-1]->GetAbsWNDC()/pad[ip]->GetAbsWNDC();
    Float_t yFactor = pad[NPad-1]->GetAbsHNDC()/pad[ip]->GetAbsHNDC();
    // hIonProb2D[ii]->GetZaxis()->SetNdivisions(503);
    hIonProb2D[ii]->GetZaxis()->SetTitleFont(fonttype);
    hIonProb2D[ii]->GetZaxis()->SetTitleOffset(tzoffset);
    hIonProb2D[ii]->GetZaxis()->SetTitleSize(tzsize);
    hIonProb2D[ii]->GetZaxis()->SetLabelFont(fonttype);
    hIonProb2D[ii]->GetZaxis()->SetLabelSize(lzsize);
    hIonProb2D[ii]->GetZaxis()->SetLabelOffset(lyoffset);
    hIonProb2D[ii]->GetZaxis()->SetTickLength(xFactor*tylength/yFactor);

  
    exIonP->Draw();
    hIonProb2D[ii]->Draw(drawopt);

    // ADK contours
    if(graphsI2D) {
      for(Int_t i=0;i<graphsI2D->GetEntriesFast();i++) {
	TGraph *gr = (TGraph*) graphsI2D->At(i);
	if(!gr) continue;
	gr->Draw("C");
      }
    }
    
    {
      Float_t rightmin = IPmin;
      Float_t rightmax = IPmax;
      Float_t slope = (yMax-yMin)/(rightmax-rightmin);
      
      for(Int_t j=0;j<hIonProb1D[ii]->GetNbinsX();j++) {
	hIonProb1D[ii]->SetBinContent(j+1,slope*(hIonProb1D[ii]->GetBinContent(j+1)-rightmin)+yMin);
      }
    }

    if(opt.Contains("off")) {
      TLine *lineOff = new TLine(xMin,xoff,xMax,xoff);
      lineOff->SetLineColor(lineColor);
      lineOff->SetLineStyle(2);
      lineOff->Draw();
    }

    
    hIonProb1D[ii]->SetLineStyle(1);
    hIonProb1D[ii]->SetLineWidth(2);
    hIonProb1D[ii]->SetLineColor(lineColor);

    hIonProb1D[ii]->Draw("sameL");
    
    if(opt.Contains("bopar")) {
      
      TEllipse *ellipLo = new TEllipse(EzCross[0],0.0,radius);	
      ellipLo->SetFillStyle(0);
      ellipLo->SetLineStyle(3);
      ellipLo->SetLineColor(kGray+3);
      ellipLo->Draw();
	
      TEllipse *ellipLu = new TEllipse(EzCross[0],0.0,radiusLu);
      ellipLu->SetFillStyle(0);
      ellipLu->SetLineStyle(4);
      ellipLu->SetLineColor(kGray+3);
      ellipLu->Draw();

    }

    if(opt.Contains("cont")) {
      // PSI MAIN contour  
      // for(Int_t i=0;i<graphsV2D_main.GetEntriesFast();i++) {
      // 	TGraph *gr = (TGraph*) graphsV2D_main.At(i);
      // 	if(!gr) continue;
      // 	gr->Draw("C");
      // }
      
      
      // PSI contours
      for(Int_t i=0;i<graphsV2D.GetEntriesFast();i++) {
	TGraph *gr = (TGraph*) graphsV2D.At(i);
	if(!gr) continue;
	gr->Draw("C");
      }
      
    }



    if(opt.Contains("1dline")) {
      lineYzero->Draw();
    }

    if(opt.Contains("sline")) {
      if(zStartPlasma>xMin && zStartPlasma<xMax)
	lineStartPlasma->Draw();
      if(zStartNeutral>xMin && zStartNeutral<xMax)
	lineStartNeutral->Draw();
      if(zEndNeutral>xMin && zEndNeutral<xMax)
	lineEndNeutral->Draw();
    }
  
    pad[ip]->Update();
  
    y1 = pad[ip]->GetBottomMargin();
    y2 = 1 - pad[ip]->GetTopMargin();
    x1 = pad[ip]->GetLeftMargin();
    x2 = 1 - pad[ip]->GetRightMargin();
  
    palette = (TPaletteAxis*) hIonProb2D[ii]->GetListOfFunctions()->FindObject("palette");  
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
      pFrame->SetShadowColor(0);
      pFrame->Draw();

    }

    TBox *lFrame = new TBox(gPad->GetUxmin(), gPad->GetUymin(),
			    gPad->GetUxmax(), gPad->GetUymax());
    lFrame->SetFillStyle(0);
    lFrame->SetLineColor(PGlobals::frameColor);
    lFrame->SetLineWidth(PGlobals::frameWidth);
    lFrame->Draw();

    pad[ip]->RedrawAxis(); 

    ip--;
    C->cd(0);
  }

  if(mask & 0x100) {
    pad[ip]->Draw();
    pad[ip]->cd(); // <---------------------------------------------------------- 9th panel
    // if(opt.Contains("logz")) {
    //   pad[ip]->SetLogz(1);
    // } else {
    //   pad[ip]->SetLogz(0);
    // }

    TH2F *hClone2D = (TH2F*) hJz2D[0]->Clone("hClone2D");
    TH1F *hClone1D = NULL;
    if(hJz1D[0])
      hClone1D = (TH1F*) hJz1D[0]->Clone("hClone1D");
    
    Float_t jzmax = hJz2D[0]->GetMaximum();
    Float_t jzmin = hJz2D[0]->GetMinimum();
    
    if(opt.Contains("betaz")) {
      hClone2D->Divide(hDen2D[0]);
      hClone2D->Scale(-1);

      if(hClone1D && hDen1D[0]) {
	hClone1D->Divide(hDen1D[0]);
	hClone1D->Scale(-1);
      }
      
      Float_t betamax = hClone2D->GetMaximum();
      Float_t betamin = hClone2D->GetMinimum();
      if(betamax > TMath::Abs(betamin))
	betamin = -betamax;
      else
	betamax = -betamin;
      hClone2D->GetZaxis()->SetRangeUser(betamin,betamax); 
      hClone2D->GetZaxis()->SetTitle("#beta_{z}");

      hClone1D->GetYaxis()->SetRangeUser(betamin,betamax); 
           
      jzmax = betamax;
      jzmin = betamin;

      
    }
    
    hFrame[ip]->Draw("col");
    
    Float_t xFactor = pad[NPad-1]->GetAbsWNDC()/pad[ip]->GetAbsWNDC();
    Float_t yFactor = pad[NPad-1]->GetAbsHNDC()/pad[ip]->GetAbsHNDC();
    //hClone2D->GetZaxis()->SetNdivisions(503);
    hClone2D->GetZaxis()->SetTitleFont(fonttype);
    hClone2D->GetZaxis()->SetTitleOffset(tzoffset);
    hClone2D->GetZaxis()->SetTitleSize(tzsize);
    hClone2D->GetZaxis()->SetLabelFont(fonttype);
    hClone2D->GetZaxis()->SetLabelSize(lzsize);
    hClone2D->GetZaxis()->SetLabelOffset(lyoffset);
    hClone2D->GetZaxis()->SetTickLength(xFactor*tylength/yFactor);

  
    exField->Draw();
    hClone2D->Draw(drawopt);

    {
       Float_t rightmin = jzmin;
       Float_t rightmax = jzmax;
       Float_t slope = (yMax-yMin)/(rightmax-rightmin);
      
       for(Int_t j=0;j<hClone1D->GetNbinsX();j++) {
	 hClone1D->SetBinContent(j+1,slope*(hClone1D->GetBinContent(j+1)-rightmin)+yMin);
       }
    }
    
    // if(opt.Contains("off")) {
    //   TLine *lineOff = new TLine(xMin,xoff,xMax,xoff);
    //   lineOff->SetLineColor(lineColor);
    //   lineOff->SetLineStyle(2);
    //   lineOff->Draw();
    // }
    
    if(opt.Contains("zero")) {
      TLine *lineZero = new TLine(xMin,0,xMax,0);
      lineZero->SetLineColor(lineColor);
      lineZero->SetLineStyle(2);
      lineZero->Draw();
    }
    
    if(opt.Contains("sline")) {
      if(zStartPlasma>xMin && zStartPlasma<xMax)
	lineStartPlasma->Draw();
      if(zStartNeutral>xMin && zStartNeutral<xMax)
	lineStartNeutral->Draw();
      if(zEndNeutral>xMin && zEndNeutral<xMax)
	lineEndNeutral->Draw();
    }

    hClone1D->SetLineStyle(1);
    hClone1D->SetLineWidth(2);
    hClone1D->SetLineColor(lineColor);
    hClone1D->Draw("sameL");
  
    pad[ip]->Update();
  
    y1 = pad[ip]->GetBottomMargin();
    y2 = 1 - pad[ip]->GetTopMargin();
    x1 = pad[ip]->GetLeftMargin();
    x2 = 1 - pad[ip]->GetRightMargin();
  
    palette = (TPaletteAxis*)hClone2D->GetListOfFunctions()->FindObject("palette");  
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
      pFrame->SetShadowColor(0);
      pFrame->Draw();

    }

    TBox *lFrame = new TBox(gPad->GetUxmin(), gPad->GetUymin(),
			    gPad->GetUxmax(), gPad->GetUymax());
    lFrame->SetFillStyle(0);
    lFrame->SetLineColor(PGlobals::frameColor);
    lFrame->SetLineWidth(PGlobals::frameWidth);
    lFrame->Draw();

    pad[ip]->RedrawAxis(); 

    ip--;
    C->cd(0);
  }


  if(mask & 0x200) {
    pad[ip]->Draw();
    pad[ip]->cd(); // <---------------------------------------------------------- 9th panel
    // if(opt.Contains("logz")) {
    //   pad[ip]->SetLogz(1);
    // } else {
    //   pad[ip]->SetLogz(0);
    // }

    cout << Form(" kmax = %f  kmin = %f",kmax,kmin) << endl;
    
    // if(kmax > TMath::Abs(kmin))
    //   kmin = -kmax;
    // else
    //   kmax = -kmin;

    kmax = 0.8;
    kmin = - kmax;
        
    hdW2D->GetZaxis()->SetRangeUser(kmin,kmax); 
    hdW2D->GetZaxis()->SetTitle("K_{x} [m#omega_{p}^{2}]");
        
    hFrame[ip]->Draw("axis");
    
    Float_t xFactor = pad[NPad-1]->GetAbsWNDC()/pad[ip]->GetAbsWNDC();
    Float_t yFactor = pad[NPad-1]->GetAbsHNDC()/pad[ip]->GetAbsHNDC();
    //hdW2D->GetZaxis()->SetNdivisions(503);
    hdW2D->GetZaxis()->SetTitleFont(fonttype);
    hdW2D->GetZaxis()->SetTitleOffset(tzoffset);
    hdW2D->GetZaxis()->SetTitleSize(tzsize);
    hdW2D->GetZaxis()->SetLabelFont(fonttype);
    hdW2D->GetZaxis()->SetLabelSize(lzsize);
    hdW2D->GetZaxis()->SetLabelOffset(lyoffset);
    hdW2D->GetZaxis()->SetTickLength(xFactor*tylength/yFactor);

  
    exField->Draw();
    hdW2D->Draw(drawopt + "0");

    if(hdW1D) {
      Float_t rightmin = kmin;
      Float_t rightmax = kmax;
      Float_t slope = (yMax-yMin)/(rightmax-rightmin);
      
      for(Int_t j=0;j<hdW1D->GetNbinsX();j++) {
	hdW1D->SetBinContent(j+1,slope*(hdW1D->GetBinContent(j+1)-rightmin)+yMin);
      }
    }
    
    // if(opt.Contains("off")) {
    //   TLine *lineOff = new TLine(xMin,xoff,xMax,xoff);
    //   lineOff->SetLineColor(lineColor);
    //   lineOff->SetLineStyle(2);
    //   lineOff->Draw();
    // }
    
    if(opt.Contains("zero")) {
      TLine *lineZero = new TLine(xMin,0,xMax,0);
      lineZero->SetLineColor(lineColor);
      lineZero->SetLineStyle(2);
      lineZero->Draw();
    }
    
    if(opt.Contains("sline")) {
      if(zStartPlasma>xMin && zStartPlasma<xMax)
	lineStartPlasma->Draw();
      if(zStartNeutral>xMin && zStartNeutral<xMax)
	lineStartNeutral->Draw();
      if(zEndNeutral>xMin && zEndNeutral<xMax)
	lineEndNeutral->Draw();
    }

    if(hdW1D) {
      hdW1D->SetLineStyle(1);
      hdW1D->SetLineWidth(2);
      hdW1D->SetLineColor(lineColor);
      hdW1D->Draw("same L");
    }
    
    pad[ip]->Update();
  
    y1 = pad[ip]->GetBottomMargin();
    y2 = 1 - pad[ip]->GetTopMargin();
    x1 = pad[ip]->GetLeftMargin();
    x2 = 1 - pad[ip]->GetRightMargin();
  
    palette = (TPaletteAxis*) hdW2D->GetListOfFunctions()->FindObject("palette");  
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
      pFrame->SetShadowColor(0);
      pFrame->Draw();

    }

    TBox *lFrame = new TBox(gPad->GetUxmin(), gPad->GetUymin(),
			    gPad->GetUxmax(), gPad->GetUymax());
    lFrame->SetFillStyle(0);
    lFrame->SetLineColor(PGlobals::frameColor);
    lFrame->SetLineWidth(PGlobals::frameWidth);
    lFrame->Draw();

    pad[ip]->RedrawAxis(); 

    ip--;
    C->cd(0);
  }

  // Writing the figure
  TString fOutName = Form("./%s/Plots/Snapshots/Snapshot%s-%s_%i",pData->GetPath().c_str(),imask.c_str(),pData->GetName(),timestep);
  if(opt.Contains("zyslc"))
    fOutName = Form("./%s/Plots/Snapshots/SnapshotZY%s-%s_%i",pData->GetPath().c_str(),imask.c_str(),pData->GetName(),timestep);
  PGlobals::imgconv(C,fOutName,opt);
  // ---------------------------------------------------------

  // PGlobals::DestroyCanvases();
}
 
 
