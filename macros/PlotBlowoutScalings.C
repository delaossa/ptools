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

#include "PData.hh"
#include "PGlobals.hh"
#include "PPalette.hh"

using namespace std;

void PlotBlowoutScalings(Int_t time, const TString &options="") {

#ifdef __CINT__  
  gSystem->Load("libptools.so");
#endif

  PGlobals::Initialize();

  TString opt = options;

  const Int_t Nsim = 11;
  char sName[Nsim][36] = {"rake-v1.0kA.G.SR2.MC01.3D", 
			  "rake-v2.5kA.G.SR2.MC01.3D",
			  "rake-v04kA.G.SR2.MC01.3D",
			  "rake-v06kA.G.SR2.MC01.3D",
			  "rake-v08kA.G.SR2.MC01.3D",
			  "rake-v10kA.G.SR2.MC01.3D",
			  "rake-v15kA.G.SR2.MC01.3D",
			  "rake-v20kA.G.MC01.3D",
  			  "rake-v30kA.G.SR2.MC01.3D",
  			  "rake-v45kA.G.SR2.MC01.3D",
  			  "rake-v60kA.G.SR2.MC01.3D"};

  // Load first simulation data (for instance)
  PData *pData = PData::Get(sName[0]);
  pData->LoadFileNames(time);
  if(!pData->IsInit()) return;

  // Some plasma constants
  Double_t n0 = pData->GetPlasmaDensity();
  Double_t kp = pData->GetPlasmaK();
  Double_t skindepth = 1.;
  if(kp!=0.0) skindepth = 1/kp;
  Double_t E0 = pData->GetPlasmaE0();
  Double_t  lightspeed =  PConst::c_light;

  // Some beam properties:
  Float_t gamma = pData->GetBeamGamma();

  // Trapping potential
  Float_t trapPotential = 1.0 - (1.0/gamma);

  // Time in OU
  Float_t Time = pData->GetRealTime();
  // z start of the plasma in normalized units.
  Float_t zStartPlasma = pData->GetPlasmaStart()*kp;
  // z start of the beam in normalized units.
  Float_t zStartBeam = pData->GetBeamStart()*kp;
  
  // Compulsory options
  opt += "comovcenter";

  // Centering time and z position:
  Double_t shiftz = pData->Shift(opt);
  if(opt.Contains("center")) {
    Time -= zStartPlasma;
    if(opt.Contains("comov"))      // Centers on the head of the beam.
      Time += zStartBeam;
  } 

  TFile *sFile[Nsim];
  TString *sPath[Nsim];
  Int_t color[Nsim];

  // histos and graphs
  TH2F *hEz2D[Nsim];
  TH2F *hDen2D[Nsim];
  TH1F *hCur1D[Nsim];
  TH1F *hEz1D[Nsim];
  TH1F *hV1D[Nsim];
  TH1F *hFocus1D[Nsim];
  TH1F *hDenx1D[Nsim];

  // Arrays with the quantities
  Float_t Current[Nsim];
  Float_t Zeta0[Nsim];
  Float_t RadiusX[Nsim];
  Float_t RadiusZ[Nsim];
  Float_t RadiusPsi[Nsim];
  Float_t EzMax[Nsim]; 
  Float_t EzMaxUse[Nsim];
  Float_t EzLimit[Nsim];
  Float_t EzMaxFocus[Nsim];
  Float_t EzSlope[Nsim];
  Float_t EzSlopeMaxFocus[Nsim];
  Float_t EzMaxBeam[Nsim];  
  Float_t DeltaPsi[Nsim];
  Float_t DeltaPsiFocus[Nsim];
  Float_t TransRatio[Nsim];
  Float_t TransRatioFocus[Nsim];


  // Theory
  Float_t EzMaxLU[Nsim];  
  Float_t EzMaxLO[Nsim];  
  Float_t EzMaxUseBO[Nsim];  
  Float_t EzMaxBeamLO[Nsim];  
  Float_t EzSlopeLO[Nsim];
  Float_t RadiusLO[Nsim];
  Float_t RadiusLU[Nsim];
  Float_t RadiusPsiBO[Nsim];
  Float_t PsiMaxLO[Nsim];
  Float_t DeltaPsiLO[Nsim];
  Float_t DeltaPsiLU[Nsim];

  // Maximums and minimums for plotting
  Float_t maxCur = -999.;
  Float_t minCur = 999.;
  Float_t maxRadius = -999.;
  Float_t minRadius = 999.;
  Float_t maxEz = -999.;
  Float_t minEz = 999.;
  Float_t maxEzSlope = -999.;
  Float_t minEzSlope = 999.;
  Float_t maxPsi = -999.;
  Float_t minPsi = 999.;
  Float_t maxTR = -999.;
  Float_t minTR = 999.;


  for(Int_t i=0;i<Nsim;i++) {
    
    TString filename = Form("./%s/Plots/Snapshots/Snapshot-%s_%i.root",sName[i],sName[i],time);
    sFile[i] = new TFile(filename,"READ");
    
    // Load histos and graphs
    hEz2D[i]  = (TH2F*) sFile[i]->Get("hE2D_0"); 
    hDen2D[i] = (TH2F*) sFile[i]->Get("hDen2D_0"); 
    hEz1D[i]  = (TH1F*) sFile[i]->Get("hE1D_0"); 
    hCur1D[i] = (TH1F*) sFile[i]->Get("hCur1D_1"); 
    hV1D[i]   = (TH1F*) sFile[i]->Get("hV1D"); 
    hFocus1D[i]   = (TH1F*) sFile[i]->Get("hFocus1D"); 

    hDen2D[i]->ResetStats();
    hEz1D[i]->ResetStats();
    hCur1D[i]->ResetStats();
    hV1D[i]->ResetStats();
    
    cout << Form("\n Analyzing simulation %s : \n",sName[i]) << endl;
   
    // Driver RMS
    Double_t rmslength = hCur1D[i]->GetRMS();
  
    // Find the first point on-axis where Ez changes from positive to negative:
    Int_t MAXCROSS = 2;
    Float_t *EzCross = new Float_t[MAXCROSS];
    Float_t *EzExtr = new Float_t[MAXCROSS];
    memset(EzCross,0,sizeof(Float_t)*MAXCROSS);
    memset(EzExtr,0,sizeof(Float_t)*MAXCROSS);
  
    Int_t auxNcross = PGlobals::HCrossings(hEz1D[i],EzCross,EzExtr,MAXCROSS,0.0,2*rmslength);
  
    // for(Int_t j=0;j<auxNcross;j++) {
    //   if(opt.Contains("units"))
    // 	cout << Form(" %i Ez crossing found at zeta = %.2f um -> Ez maximum = %.2f GV/m",j,EzCross[j],EzExtr[j]) << endl;
    //   else
    // 	cout << Form(" %i Ez crossing found at zeta = %.2f kp^-1 -> Ez maximum = %.2f E0",j,EzCross[j],EzExtr[j]) << endl;      
    // }

    Zeta0[i] = EzCross[0];

    Int_t bin0 = hEz1D[i]->FindBin(Zeta0[i]);
    hDenx1D[i] = (TH1F*) hDen2D[i]->ProjectionY("_py",bin0);

    RadiusX[i] = fabs(hDenx1D[i]->GetBinCenter(hDenx1D[i]->GetMaximumBin()));
    RadiusZ[i] = fabs(hV1D[i]->GetBinCenter(hV1D[i]->GetMinimumBin()) - Zeta0[i]);

    DeltaPsi[i] =  hV1D[i]->GetBinContent(hV1D[i]->GetMaximumBin()) -  hV1D[i]->GetBinContent(hV1D[i]->GetMinimumBin());
    

    Float_t zshift = 0.5; 
    if(opt.Contains("units")) 
      zshift *= PUnits::um * kp;
    
    Int_t bin1 = hEz1D[i]->FindBin(Zeta0[i]-zshift);
    Int_t bin2 = hEz1D[i]->FindBin(Zeta0[i]+zshift);
    EzSlope[i] = (hEz1D[i]->GetBinContent(bin1)-hEz1D[i]->GetBinContent(bin2))/
	(hEz1D[i]->GetBinCenter(bin1)-hEz1D[i]->GetBinCenter(bin2));

    EzMaxUse[i]  = EzSlope[i] * RadiusZ[i];
    EzMaxBeam[i] = EzExtr[0];
    EzMax[i] = fabs(EzExtr[1]);

    Current[i] = hCur1D[i]->GetMaximum();
    if(opt.Contains("units")) {
      Current[i] *= PUnits::kA / PConst::I0;
      rmslength *= PUnits::um * kp;
    }
    RadiusLO[i] = PFunc::RadiusBO(Current[i],rmslength);
    RadiusLU[i] = 2.0 * TMath::Sqrt(Current[i]);
    EzMaxLU[i] = PFunc::EzMaxLu(Current[i]);
    EzMaxBeamLO[i] = PFunc::EzMaxBeamLo(Current[i]);
    EzMaxLO[i] = RadiusLO[i]/2.0;
    Float_t RadiusLO2 = RadiusLO[i]*RadiusLO[i];
    Float_t RadiusLU2 = RadiusLU[i]*RadiusLU[i];
    // EzSlopeLO[i] = RadiusLO2 * (RadiusLO2 + 8.0) / (2.0 * (RadiusLO2 + 4.0) * (RadiusLO2 + 4.0));
    EzSlopeLO[i] = RadiusLU2 * (RadiusLU2 + 8.0) / (2.0 * (RadiusLU2 + 4.0) * (RadiusLU2 + 4.0));
    PsiMaxLO[i] = RadiusLO[i]*RadiusLO[i]/4.0;
    DeltaPsiLO[i] = RadiusLO[i]*RadiusLO[i]/4.0;
    DeltaPsiLU[i] = RadiusLU[i]*RadiusLU[i]/4.0;
    RadiusPsiBO[i] = TMath::Sqrt(2*(DeltaPsiLO[i]-1)); 
    EzMaxUseBO[i]  = RadiusPsiBO[i]/2.0;
    if(opt.Contains("units")) {
      Current[i] /= PUnits::kA / PConst::I0;
      rmslength /= PUnits::um * kp;
      RadiusLO[i] *= skindepth / PUnits::um;
      RadiusLU[i] *= skindepth / PUnits::um;
      RadiusPsiBO[i] *= skindepth / PUnits::um;
      EzMaxLU[i]  *= E0 / (PUnits::GV/PUnits::m);
      EzMaxBeamLO[i] *= E0 / (PUnits::GV/PUnits::m);
      EzMaxUseBO[i] *= E0 / (PUnits::GV/PUnits::m);
      EzSlope[i] *= (E0 / (PUnits::GV/PUnits::m)) * skindepth / PUnits::um;
      cout << Form(" Beam max. current = %.2f kA   RMS = %.2f um",Current[i],rmslength) << endl;
      cout << Form(" Radius of the bubble (Lotov): R = %.2f um",RadiusLO[i]) << endl;
      cout << Form(" Radius of the bubble (Lu)   : R = %.2f um",RadiusLU[i]) << endl;
      cout << Form(" Ez max beam (Lo) = %.2f GV/m",EzMaxBeamLO[i]) << endl;
      cout << Form(" Ez max (Lu) = %.2f GV/m",-EzMaxLU[i]) << endl;
      cout << Form(" Ez slope at crossing = %.2f (GV/m)/um ",EzSlope[i]) << endl;
    } else {
      cout << Form(" Beam max. current = %.2f I0   RMS = %.2f kp^-1",Current[i],rmslength) << endl;
      cout << Form(" Radius of the bubble (Lotov): R = %.2f kp^-1",RadiusLO[i]) << endl;
      cout << Form(" Radius of the bubble (Lu)   : R = %.2f kp^-1",RadiusLU[i]) << endl;
      cout << Form(" Ez max beam (Lo) = %.2f E0",EzMaxBeamLO[i]) << endl;
      cout << Form(" Ez max (Lu) = %.2f E0",-EzMaxLU[i]) << endl;
      cout << Form(" Ez slope at crossing = %.2f ",EzSlope[i]) << endl;
    }

    if(opt.Contains("kA")) {
      Current[i] /= PUnits::kA / PConst::I0;
    }

    // Find the first Focusing crossing after the maximum 
    MAXCROSS = 1;
    Float_t *FocusCross = new Float_t[MAXCROSS];
    Float_t *FocusExtr = new Float_t[MAXCROSS];
    memset(FocusCross,0,sizeof(Float_t)*MAXCROSS);
    memset(FocusExtr,0,sizeof(Float_t)*MAXCROSS);
    auxNcross = PGlobals::HCrossings(hFocus1D[i],FocusCross,FocusExtr,MAXCROSS,0.0,Zeta0[i]);

    if(opt.Contains("units") && n0) {
      trapPotential *=  ( E0 * skindepth / (PUnits::MV) ); 
    }  

    Int_t binFcross = hFocus1D[i]->FindBin(FocusCross[0]);
    Float_t potValueFin = hV1D[i]->GetBinContent(binFcross);
    EzMaxFocus[i] = fabs(hEz1D[i]->GetBinContent(binFcross));


    TransRatio[i] = EzMax[i]/EzMaxBeam[i];
    TransRatioFocus[i] = EzMaxFocus[i]/EzMaxBeam[i];
   
    cout << Form(" Max. R       = %.2f ",TransRatio[i]) << endl;
    cout << Form(" Max. R focus = %.2f ",TransRatioFocus[i]) << endl;


    Int_t binPotValueIni = -99;
    Float_t potValueIni = -99;
    for(Int_t j=binFcross;j<hV1D[i]->GetNbinsX();j++) {
      if(hV1D[i]->GetBinContent(j)>=potValueFin+trapPotential) {
	binPotValueIni = j;
	potValueIni = hV1D[i]->GetBinContent(j);
	//break;
      }
    }

    // Ez slope at max focusing Ez
    zshift = 0.1;
    bin1 = hEz1D[i]->FindBin(FocusCross[0]-zshift);
    bin2 = hEz1D[i]->FindBin(FocusCross[0]+zshift);
    EzSlopeMaxFocus[i] = (hEz1D[i]->GetBinContent(bin1)-hEz1D[i]->GetBinContent(bin2))/
	(hEz1D[i]->GetBinCenter(bin1)-hEz1D[i]->GetBinCenter(bin2));
    cout << Form(" Ez slope at max focus = %.2f ",EzSlopeMaxFocus[i]) << endl;


    // Difference between psi maximum and psi final focus
    DeltaPsiFocus[i] = hV1D[i]->GetBinContent(hV1D[i]->GetMaximumBin()) - potValueFin;

    // Find the first Psi crossing after the maximum 
    MAXCROSS = 1;
    Float_t *PsiCross = new Float_t[MAXCROSS];
    Float_t *PsiExtr = new Float_t[MAXCROSS];
    memset(PsiCross,0,sizeof(Float_t)*MAXCROSS);
    memset(PsiExtr,0,sizeof(Float_t)*MAXCROSS);
    auxNcross = PGlobals::HCrossings(hV1D[i],PsiCross,PsiExtr,MAXCROSS,potValueIni,Zeta0[i]);

    Int_t binVcross = -1;
    EzLimit[i] = 0.0;
    RadiusPsi[i] = 0.0;
    if(auxNcross>0) {
      binVcross = hEz1D[i]->FindBin(PsiCross[0]);
      EzLimit[i] = fabs(hEz1D[i]->GetBinContent(binVcross));
      RadiusPsi[i] = fabs(PsiCross[0]- Zeta0[i]);
    }
        
    cout << Form(" Trapping zone radius = %.2f ",RadiusPsi[i]) << endl;

    if(RadiusX[i]>maxRadius) maxRadius = RadiusX[i];
    if(RadiusZ[i]>maxRadius) maxRadius = RadiusZ[i];
    if(RadiusLO[i]>maxRadius) maxRadius = RadiusLO[i];
    if(RadiusLU[i]>maxRadius) maxRadius = RadiusLU[i];

    if(RadiusX[i]<minRadius) minRadius = RadiusX[i];
    if(RadiusZ[i]<minRadius) minRadius = RadiusZ[i];
    if(RadiusLO[i]<minRadius) minRadius = RadiusLO[i];
    if(RadiusLU[i]<minRadius) minRadius = RadiusLU[i];

    if(EzMax[i]>maxEz) maxEz = EzMax[i];
    if(EzMaxUse[i]>maxEz) maxEz = EzMaxUse[i];
    if(EzLimit[i]>maxEz) maxEz = EzLimit[i];
    if(EzMaxFocus[i]>maxEz) maxEz = EzMaxFocus[i];
    if(EzMaxBeam[i]>maxEz) maxEz = EzMaxBeam[i];
    if(EzMaxLO[i]>maxEz) maxEz = EzMaxLO[i];
    if(EzMaxLU[i]>maxEz) maxEz = EzMaxLU[i];
    if(EzMaxBeamLO[i]>maxEz) maxEz = EzMaxBeamLO[i];

    if(EzMax[i]<minEz) minEz = EzMax[i];
    if(EzMaxUse[i]<minEz) minEz = EzMaxUse[i];
    if(EzLimit[i]<minEz) minEz = EzLimit[i];
    if(EzMaxFocus[i]<minEz) minEz = EzMaxFocus[i];
    if(EzMaxBeam[i]<minEz) minEz = EzMaxBeam[i];
    if(EzMaxLO[i]<minEz) minEz = EzMaxLO[i];
    if(EzMaxLU[i]<minEz) minEz = EzMaxLU[i];
    if(EzMaxBeamLO[i]<minEz) minEz = EzMaxBeamLO[i];
    
    if(EzSlope[i]>maxEzSlope) maxEzSlope = EzSlope[i]; 
    if(EzSlope[i]<minEzSlope) minEzSlope = EzSlope[i]; 
    // if(EzSlopeMaxFocus[i]>maxEzSlope) maxEzSlope = EzSlopeMaxFocus[i]; 
    // if(EzSlopeMaxFocus[i]<minEzSlope) minEzSlope = EzSlopeMaxFocus[i]; 

    if(DeltaPsi[i]>maxPsi) maxPsi = DeltaPsi[i];
    if(DeltaPsiFocus[i]>maxPsi) maxPsi = DeltaPsiFocus[i];
    if(DeltaPsiLO[i]>maxPsi) maxPsi = DeltaPsiLO[i];
    if(DeltaPsiLU[i]>maxPsi) maxPsi = DeltaPsiLU[i];
    if(DeltaPsi[i]<minPsi) minPsi = DeltaPsi[i];
    if(DeltaPsiFocus[i]<minPsi) minPsi = DeltaPsiFocus[i];
    if(DeltaPsiLO[i]<minPsi) minPsi = DeltaPsiLO[i];
    if(DeltaPsiLU[i]<minPsi) minPsi = DeltaPsiLU[i];
    
    if(TransRatio[i]<minTR) minTR = TransRatio[i];
    if(TransRatio[i]>maxTR) maxTR = TransRatio[i];
    if(TransRatioFocus[i]<minTR) minTR = TransRatioFocus[i];
    if(TransRatioFocus[i]>maxTR) maxTR = TransRatioFocus[i];

    if(Current[i]<minCur) minCur = Current[i];
    if(Current[i]>maxCur) maxCur = Current[i];

  }

  Float_t xMin = minCur - (maxCur-minCur)*0.1;
  Float_t xMax = maxCur + (maxCur-minCur)*0.1;


  // Graphs for the max. current dependence
  TGraph *gRadiusX  = new TGraph(Nsim,Current,RadiusX);
  TGraph *gRadiusZ  = new TGraph(Nsim,Current,RadiusZ);
  TGraph *gRadiusPsi  = new TGraph(Nsim,Current,RadiusPsi);
  TGraph *gRadiusLO  = new TGraph(Nsim,Current,RadiusLO);
  TGraph *gRadiusLU  = new TGraph(Nsim,Current,RadiusLU);
  TGraph *gRadiusPsiBO  = new TGraph(Nsim,Current,RadiusPsiBO);

  TGraph *gEzMax    = new TGraph(Nsim,Current,EzMax);
  TGraph *gEzMaxUse = new TGraph(Nsim,Current,EzMaxUse);
  TGraph *gEzMaxBeam = new TGraph(Nsim,Current,EzMaxBeam);
  TGraph *gEzLimit = new TGraph(Nsim,Current,EzLimit);
  TGraph *gEzMaxFocus = new TGraph(Nsim,Current,EzMaxFocus);
  TGraph *gEzMaxLU   = new TGraph(Nsim,Current,EzMaxLU);
  TGraph *gEzMaxLO  = new TGraph(Nsim,Current,EzMaxLO);
  TGraph *gEzMaxUseBO  = new TGraph(Nsim,Current,EzMaxUseBO);
  TGraph *gEzMaxBeamLO = new TGraph(Nsim,Current,EzMaxBeamLO);
 
  TGraph *gDeltaPsi = new TGraph(Nsim,Current,DeltaPsi);
  TGraph *gDeltaPsiFocus = new TGraph(Nsim,Current,DeltaPsiFocus);
  TGraph *gDeltaPsiLO = new TGraph(Nsim,Current,DeltaPsiLO);
  TGraph *gDeltaPsiLU = new TGraph(Nsim,Current,DeltaPsiLU);

  TGraph *gEzSlope = new TGraph(Nsim,Current,EzSlope);
  TGraph *gEzSlopeMaxFocus = new TGraph(Nsim,Current,EzSlopeMaxFocus);
  TGraph *gEzSlopeLO = new TGraph(Nsim,Current,EzSlopeLO);

  TGraph *gTransRatio = new TGraph(Nsim,Current,TransRatio);
  TGraph *gTransRatioFocus = new TGraph(Nsim,Current,TransRatioFocus);

  // Int_t color1 = kMagenta-8;
  // Int_t color2 = kMagenta-2;
  // Int_t color3 = kMagenta-4;
  // Int_t color4 = kAzure+1;
  // Int_t color5 = kGray+2;
  // Int_t markersty = 20;
  // Float_t markersiz = 1.0;
  // Int_t linewit = 2;

  Int_t color1 = kGray+3;
  Int_t color2 = kGray+2;
  Int_t color3 = kAzure+1;
  Int_t color4 = kRed;
  Int_t color5 = kGray+1;
  Int_t markersty = 20;
  Float_t markersiz = 1.0;
  Int_t linewit = 1;

  gRadiusX->SetLineColor(color1);
  gRadiusX->SetLineWidth(linewit);
  gRadiusX->SetMarkerColor(color1);
  gRadiusX->SetMarkerStyle(markersty);
  gRadiusX->SetMarkerSize(markersiz);

  gRadiusLU->SetLineColor(color1);
  gRadiusLU->SetLineWidth(1);
  gRadiusLU->SetLineStyle(4);
  gRadiusLU->SetMarkerColor(color1);
  gRadiusLU->SetMarkerStyle(markersty);
  gRadiusLU->SetMarkerSize(0.0);

  gRadiusZ->SetLineColor(color2);
  gRadiusZ->SetLineWidth(linewit);
  gRadiusZ->SetMarkerColor(color2);
  gRadiusZ->SetMarkerStyle(markersty);
  gRadiusZ->SetMarkerSize(markersiz);

  gRadiusLO->SetLineColor(color2);
  gRadiusLO->SetLineWidth(1);
  gRadiusLO->SetLineStyle(3);
  gRadiusLO->SetMarkerColor(color2);
  gRadiusLO->SetMarkerStyle(markersty);
  gRadiusLO->SetMarkerSize(0.0);

  gRadiusPsi->SetLineColor(color3);
  gRadiusPsi->SetLineWidth(linewit);
  gRadiusPsi->SetMarkerColor(color3);
  gRadiusPsi->SetMarkerStyle(markersty);
  gRadiusPsi->SetMarkerSize(markersiz);

  gEzMax->SetLineColor(color2);
  gEzMax->SetLineWidth(linewit);
  gEzMax->SetMarkerColor(color2);
  gEzMax->SetMarkerStyle(markersty);
  gEzMax->SetMarkerSize(markersiz);

  gEzMaxUse->SetLineColor(color2);
  gEzMaxUse->SetLineWidth(linewit);
  gEzMaxUse->SetMarkerColor(color2);
  gEzMaxUse->SetMarkerStyle(markersty);
  gEzMaxUse->SetMarkerSize(markersiz);

  gEzLimit->SetLineColor(color3);
  gEzLimit->SetLineWidth(linewit);
  gEzLimit->SetMarkerColor(color3);
  gEzLimit->SetMarkerStyle(markersty);
  gEzLimit->SetMarkerSize(markersiz);

  gEzMaxFocus->SetLineColor(color1);
  gEzMaxFocus->SetLineWidth(linewit);
  gEzMaxFocus->SetMarkerColor(color1);
  gEzMaxFocus->SetMarkerStyle(markersty);
  gEzMaxFocus->SetMarkerSize(markersiz);

  gEzSlope->SetLineColor(color2);
  gEzSlope->SetLineWidth(linewit);
  gEzSlope->SetMarkerColor(color2);
  gEzSlope->SetMarkerStyle(markersty);
  gEzSlope->SetMarkerSize(markersiz);

  gEzSlopeMaxFocus->SetLineColor(color1);
  gEzSlopeMaxFocus->SetLineWidth(linewit);
  gEzSlopeMaxFocus->SetMarkerColor(color1);
  gEzSlopeMaxFocus->SetMarkerStyle(markersty);
  gEzSlopeMaxFocus->SetMarkerSize(markersiz);

  gEzSlopeLO->SetLineColor(color2);
  gEzSlopeLO->SetLineWidth(1);
  gEzSlopeLO->SetLineStyle(3);
  gEzSlopeLO->SetMarkerColor(color2);
  gEzSlopeLO->SetMarkerStyle(markersty);
  gEzSlopeLO->SetMarkerSize(0);


  gDeltaPsi->SetLineColor(color2);
  gDeltaPsi->SetLineWidth(linewit);
  gDeltaPsi->SetMarkerColor(color2);
  gDeltaPsi->SetMarkerStyle(markersty);
  gDeltaPsi->SetMarkerSize(markersiz);

  gDeltaPsiFocus->SetLineColor(color1);
  gDeltaPsiFocus->SetLineWidth(linewit);
  gDeltaPsiFocus->SetMarkerColor(color1);
  gDeltaPsiFocus->SetMarkerStyle(markersty);
  gDeltaPsiFocus->SetMarkerSize(markersiz);

  gEzMaxBeam->SetLineColor(color4);
  gEzMaxBeam->SetLineWidth(linewit);
  gEzMaxBeam->SetMarkerColor(color4);
  gEzMaxBeam->SetMarkerStyle(markersty);
  gEzMaxBeam->SetMarkerSize(markersiz);

  gRadiusPsiBO->SetLineColor(color3);
  gRadiusPsiBO->SetLineWidth(1);
  gRadiusPsiBO->SetLineStyle(3);
  gRadiusPsiBO->SetMarkerColor(color3);
  gRadiusPsiBO->SetMarkerStyle(markersty);
  gRadiusPsiBO->SetMarkerSize(0.0);

  gEzMaxLU->SetLineColor(color5);
  gEzMaxLU->SetLineWidth(1);
  gEzMaxLU->SetLineStyle(4);
  gEzMaxLU->SetMarkerColor(color5);
  gEzMaxLU->SetMarkerStyle(markersty);
  gEzMaxLU->SetMarkerSize(0.0);

  gEzMaxLO->SetLineColor(color2);
  gEzMaxLO->SetLineWidth(1);
  gEzMaxLO->SetLineStyle(3);
  gEzMaxLO->SetMarkerColor(color2);
  gEzMaxLO->SetMarkerStyle(markersty);
  gEzMaxLO->SetMarkerSize(0.0);

  gEzMaxUseBO->SetLineColor(color3);
  gEzMaxUseBO->SetLineWidth(1);
  gEzMaxUseBO->SetLineStyle(3);
  gEzMaxUseBO->SetMarkerColor(color3);
  gEzMaxUseBO->SetMarkerStyle(markersty);
  gEzMaxUseBO->SetMarkerSize(0.0);

  gEzMaxBeamLO->SetLineColor(color4);
  gEzMaxBeamLO->SetLineWidth(1);
  gEzMaxBeamLO->SetLineStyle(3);
  gEzMaxBeamLO->SetMarkerColor(color4);
  gEzMaxBeamLO->SetMarkerStyle(markersty);
  gEzMaxBeamLO->SetMarkerSize(0.0);

  gDeltaPsiLO->SetLineColor(color2);
  gDeltaPsiLO->SetLineWidth(1);
  gDeltaPsiLO->SetLineStyle(3);
  gDeltaPsiLO->SetMarkerColor(color2);
  gDeltaPsiLO->SetMarkerStyle(markersty);
  gDeltaPsiLO->SetMarkerSize(0.0);

  gDeltaPsiLU->SetLineColor(color2);
  gDeltaPsiLU->SetLineWidth(1);
  gDeltaPsiLU->SetLineStyle(4);
  gDeltaPsiLU->SetMarkerColor(color2);
  gDeltaPsiLU->SetMarkerStyle(markersty);
  gDeltaPsiLU->SetMarkerSize(0.0);


  gTransRatio->SetLineColor(color2);
  gTransRatio->SetLineWidth(linewit);
  gTransRatio->SetMarkerColor(color2);
  gTransRatio->SetMarkerStyle(markersty);
  gTransRatio->SetMarkerSize(markersiz);

  gTransRatioFocus->SetLineColor(color1);
  gTransRatioFocus->SetLineWidth(linewit);
  gTransRatioFocus->SetMarkerColor(color1);
  gTransRatioFocus->SetMarkerStyle(markersty);
  gTransRatioFocus->SetMarkerSize(markersiz);


  // Canvas setup
  // Create the canvas and the pads before the Frame loop
  // Resolution:
  Int_t sizex = 800;
  Int_t sizey = 1000;
  char cName[32];

  sprintf(cName,"C");     
  TCanvas *C = (TCanvas*) gROOT->FindObject(cName);
  if(C==NULL) C = new TCanvas("C","",sizex,sizey);
  C->cd();

  C->Clear();

  // Setup Pad layout: 
  const Int_t Npads = 5;
 
  // Text objects
  // TPaveText **textLabel = new TPaveText*[Npads];

  TPad **pad = new TPad*[Npads];
  TH1F *hFrame[Npads];
  TGaxis *axis[Npads];
  TLegend *Leg[Npads];

  // Setup Pad layout:
  Double_t lMargin = 0.14;
  Double_t rMargin = 0.05;
  Double_t bMargin = 0.10;
  Double_t tMargin = 0.04;
  Float_t mMargin = 0.01;
  PGlobals::CanvasPartition(C,Npads,lMargin,rMargin,bMargin,tMargin,mMargin);
  //PGlobals::CanvasAsymPartition(C,Npads,lMargin,rMargin,bMargin,tMargin,1.,mMargin);
  
  // Define the frames for plotting
  Int_t fonttype = 43;
  Int_t fontsize = 24;
  Int_t tfontsize = 28;
  Float_t txoffset = 3.8;
  Float_t lxoffset = 0.02;
  Float_t tyoffset = 2.0;
  Float_t lyoffset = 0.01;
  Float_t tylength = 0.02;
  Float_t txlength = 0.04;
  for(Int_t i=0;i<Npads;i++) {
    char name[16];
    sprintf(name,"pad_%i",i);
    pad[i] = (TPad*) gROOT->FindObject(name);
    pad[i]->SetFrameLineWidth(2);  
    // pad[i]->SetTickx(1);
    // pad[i]->SetTicky(1);
    
    sprintf(name,"hFrame_%i",i);  
    hFrame[i] = (TH1F*) gROOT->FindObject(name);
    if(hFrame[i]) delete hFrame[i];

    hFrame[i] = new TH1F(name,"",100,xMin,xMax);
    hFrame[i]->Reset();
    
    Float_t xFactor = pad[0]->GetAbsWNDC()/pad[i]->GetAbsWNDC();
    Float_t yFactor = pad[0]->GetAbsHNDC()/pad[i]->GetAbsHNDC();
    
    // Format for y axis
    hFrame[i]->GetYaxis()->SetTitleFont(fonttype);
    hFrame[i]->GetYaxis()->SetTitleSize(tfontsize);
    hFrame[i]->GetYaxis()->SetTitleOffset(tyoffset);
    hFrame[i]->GetYaxis()->SetLabelFont(fonttype);
    hFrame[i]->GetYaxis()->SetLabelSize(fontsize);
    hFrame[i]->GetYaxis()->SetLabelOffset(lyoffset);
    hFrame[i]->GetYaxis()->SetTickLength(xFactor*tylength/yFactor);
    hFrame[i]->GetYaxis()->CenterTitle();
    
    // Format for x axis
    hFrame[i]->GetXaxis()->SetTitleFont(fonttype);
    hFrame[i]->GetXaxis()->SetTitleSize(tfontsize+2);
    hFrame[i]->GetXaxis()->SetTitleOffset(txoffset);
    hFrame[i]->GetXaxis()->SetLabelFont(fonttype);
    hFrame[i]->GetXaxis()->SetLabelSize(fontsize+2);
    hFrame[i]->GetXaxis()->SetLabelOffset(lxoffset);
    hFrame[i]->GetXaxis()->CenterTitle();
    hFrame[i]->GetXaxis()->SetTickLength(yFactor*txlength/xFactor);      
  }


  C->cd(0);

  pad[0]->Draw();
  pad[0]->cd();

  //  gPad->SetLogx(1);
  Float_t yMin = minTR - (maxTR-minTR)*0.2;
  Float_t yMax = maxTR + (maxTR-minTR)*0.2;
  if(yMin<0.0) yMin = 0.001;
  hFrame[0]->GetYaxis()->SetRangeUser(yMin,yMax);  
  
  hFrame[0]->GetXaxis()->SetTitle("#Lambda_{b}");
  if(opt.Contains("kA"))
    hFrame[0]->GetXaxis()->SetTitle("I_{b}^{ 0} [kA]");

  hFrame[0]->GetYaxis()->SetTitle("R");
  hFrame[0]->Draw("AXIS");

  // if(opt.Contains("theo")) {
  //   gDeltaPsiLO->Draw("PL");
  //   gDeltaPsiLU->Draw("PL");
  // }

  gTransRatio->Draw("PL");
  gTransRatioFocus->Draw("PL");

  Int_t fontsizeleg = 20;
  if(opt.Contains("leg")) { 
    Leg[0]=new TLegend(1 - gPad->GetRightMargin() - 0.26, gPad->GetBottomMargin() + 0.40,
		       1 - gPad->GetRightMargin() - 0.02, gPad->GetBottomMargin() + 0.55); 
    PGlobals::SetPaveStyle(Leg[0]);
    Leg[0]->SetTextAlign(13);
    Leg[0]->SetTextColor(kGray+3);
    Leg[0]->SetTextFont(43);
    Leg[0]->SetTextSizePixels(fontsizeleg);
    Leg[0]->SetLineColor(1);
    Leg[0]->SetBorderSize(0);
    Leg[0]->SetFillColor(0);
    Leg[0]->SetFillStyle(1001);
    Leg[0]->SetFillStyle(0); // Hollow

    Leg[0]->AddEntry(gTransRatio,"Max R","LP");
    Leg[0]->AddEntry(gTransRatioFocus,"Max R in focus","LP");

    Leg[0]->Draw();
  }    

  
  C->cd(0);
  pad[1]->Draw();
  pad[1]->cd();

  //  gPad->SetLogx(1);
  yMin = minPsi - (maxPsi-minPsi)*0.2;
  yMax = maxPsi + (maxPsi-minPsi)*0.2;
  //if(yMin<0.0) yMin = 0.001;
  hFrame[1]->GetYaxis()->SetRangeUser(yMin,yMax);  
  
  hFrame[1]->GetXaxis()->SetTitle("#Lambda_{b}");
  if(opt.Contains("kA"))
    hFrame[1]->GetXaxis()->SetTitle("I_{b} [kA]");

  hFrame[1]->GetYaxis()->SetTitle("#Delta#psi");
  hFrame[1]->Draw("AXIS");

  TLine *line1 = new TLine(xMin,1.0,xMax,1.0);
  line1->SetLineColor(kGray+2);
  line1->SetLineStyle(2);
  line1->Draw();

  if(opt.Contains("theo")) {
    gDeltaPsiLO->Draw("PL");
    gDeltaPsiLU->Draw("PL");
  }
  gDeltaPsi->Draw("PL");
  gDeltaPsiFocus->Draw("PL");

  if(opt.Contains("leg")) { 
    Leg[1]=new TLegend(gPad->GetLeftMargin() + 0.03, 1 - gPad->GetTopMargin() - 0.32,
		       gPad->GetLeftMargin() + 0.24, 1 - gPad->GetTopMargin() - 0.08); 
    PGlobals::SetPaveStyle(Leg[1]);
    Leg[1]->SetTextAlign(13);
    Leg[1]->SetTextColor(kGray+3);
    Leg[1]->SetTextFont(43);
    Leg[1]->SetTextSizePixels(fontsizeleg);
    Leg[1]->SetLineColor(1);
    Leg[1]->SetBorderSize(0);
    Leg[1]->SetFillColor(0);
    Leg[1]->SetFillStyle(1001);
    Leg[1]->SetFillStyle(0); // Hollow

    Leg[1]->AddEntry(gDeltaPsi,"Max #Delta#psi","LP");
    Leg[1]->AddEntry(gDeltaPsiFocus,"Max #Delta#psi in focus","LP");
    // Leg[1]->AddEntry(gDeltaPsiLO,"Max #Psi Lotov","L");
    // Leg[1]->AddEntry(gDeltaPsiLU,"Max #Psi Lu","L");

    Leg[1]->Draw();
  }    

  C->cd(0);
  pad[2]->Draw();
  pad[2]->cd();

  // gPad->SetLogx(1);
  yMin = minEz - (maxEz-minEz)*0.2;
  yMax = maxEz + (maxEz-minEz)*0.2;
  if(yMin<=0.0) yMin = -0.999;
  hFrame[2]->GetYaxis()->SetRangeUser(yMin,yMax);  
  hFrame[2]->GetYaxis()->SetTitle("|E_{z}/E_{0}|");
  hFrame[2]->Draw("AXIS");

  TLine *line0 = new TLine(xMin,0.0,xMax,0.0);
  line0->SetLineColor(kGray+2);
  line0->SetLineStyle(2);
  line0->Draw();

  if(opt.Contains("theo")) {
    gEzMaxLO->Draw("PL");
    //gEzMaxUseBO->Draw("PL");
    gEzMaxBeamLO->Draw("PL");
    //   gEzMaxLU->Draw("PL");
  }

  gEzMax->Draw("PL");
  // gEzMaxUse->Draw("PL");
  gEzMaxBeam->Draw("PL");
  gEzLimit->Draw("PL");
  gEzMaxFocus->Draw("PL");

  if(opt.Contains("leg")) { 
    Leg[2]=new TLegend(gPad->GetLeftMargin() + 0.03, 1 - gPad->GetTopMargin() - 0.4, gPad->GetLeftMargin() + 0.24, 1 - gPad->GetTopMargin() - 0.04); 
    PGlobals::SetPaveStyle(Leg[2]);
    Leg[2]->SetTextAlign(13);
    Leg[2]->SetTextColor(kGray+3);
    Leg[2]->SetTextFont(43);
    Leg[2]->SetTextSizePixels(fontsizeleg-2);
    Leg[2]->SetLineColor(1);
    Leg[2]->SetBorderSize(0);
    Leg[2]->SetFillColor(0);
    Leg[2]->SetFillStyle(1001);
    Leg[2]->SetFillStyle(0); // Hollow

    Leg[2]->AddEntry(gEzMax,"Min E_{z}","LP");
    Leg[2]->AddEntry(gEzMaxFocus,"Min E_{z} in focus","LP");
    Leg[2]->AddEntry(gEzLimit,"Min E_{z} trapping","LP");
    Leg[2]->AddEntry(gEzMaxBeam,"Max E_{z} beam","LP");

    Leg[2]->Draw();
  }

  C->cd(0);
  pad[3]->Draw();
  pad[3]->cd();

  // gPad->SetLogx(1);
  yMin = minEzSlope - (maxEzSlope-minEzSlope)*0.25;
  yMax = maxEzSlope + (maxEzSlope-minEzSlope)*0.25;
  if(yMin<=0.0) yMin = -0.999;
  //  hFrame[3]->GetYaxis()->SetRangeUser(yMin,yMax);
  hFrame[3]->GetYaxis()->SetRangeUser(-0.0499,0.5499);
  hFrame[3]->GetYaxis()->SetNdivisions(510);
  hFrame[3]->GetYaxis()->SetTitle("#partial_{#zeta}E_{z}^{0}/k_{p}E_{0}");
  hFrame[3]->Draw("AXIS");

  TLine *line05 = new TLine(xMin,0.5,xMax,0.5);
  line05->SetLineColor(kGray+2);
  line05->SetLineStyle(2);
  line05->Draw();

  line0->Draw();
  
  if(opt.Contains("theo")) {
    //   gEzSlopeLO->Draw("PL");
  }

  gEzSlope->Draw("PL");
  // gEzSlopeMaxFocus->Draw("PL");

  if(opt.Contains("leg")) { 
    Leg[3]=new TLegend(1 - gPad->GetRightMargin() - 0.26, gPad->GetBottomMargin() + 0.14,1- gPad->GetRightMargin() - 0.02, gPad->GetBottomMargin() + 0.24); 
    PGlobals::SetPaveStyle(Leg[3]);
    Leg[3]->SetTextAlign(12);
    Leg[3]->SetTextColor(kGray+3);
    Leg[3]->SetTextFont(43);
    Leg[3]->SetTextSizePixels(fontsizeleg);
    Leg[3]->SetLineColor(1);
    Leg[3]->SetBorderSize(0);
    Leg[3]->SetFillColor(0);
    Leg[3]->SetFillStyle(1001);
    Leg[3]->SetFillStyle(0); // Hollow

    Leg[3]->AddEntry(gEzSlope,"E_{z} slope at #zeta_{m}","LP");
 
    Leg[3]->Draw();
  }

  C->cd(0);

  pad[4]->Draw();
  pad[4]->cd();

  //  gPad->SetLogx(1);

  yMin = minRadius - (maxRadius-minRadius)*0.2;
  yMax = maxRadius + (maxRadius-minRadius)*0.2;
  if(yMin<=0.0) yMin = -0.999;
  hFrame[4]->GetYaxis()->SetRangeUser(yMin,yMax);  
  hFrame[4]->GetYaxis()->SetTitle("k_{p} x");
  hFrame[4]->Draw("AXIS");
  
  line0->Draw();

  if(opt.Contains("theo")) {
    gRadiusLO->Draw("PL");
    gRadiusLU->Draw("PL");
    //    gRadiusPsiBO->Draw("PL");
  }

  gRadiusPsi->Draw("PL");
  gRadiusX->Draw("PL");
  gRadiusZ->Draw("PL");

  if(opt.Contains("leg")) { 
    Leg[4]=new TLegend(1 - gPad->GetRightMargin() - 0.15, gPad->GetBottomMargin() + 0.15,
		       1 - gPad->GetRightMargin() + 0.09, gPad->GetBottomMargin() + 0.45); 
    PGlobals::SetPaveStyle(Leg[4]);
    Leg[4]->SetTextAlign(12);
    Leg[4]->SetTextColor(kGray+3);
    Leg[4]->SetTextFont(43);
    Leg[4]->SetTextSizePixels(fontsizeleg);
    Leg[4]->SetLineColor(1);
    Leg[4]->SetBorderSize(0);
    Leg[4]->SetFillColor(0);
    Leg[4]->SetFillStyle(1001);
    Leg[4]->SetFillStyle(0); // Hollow

    // Leg[4]->AddEntry(gRadiusZ,"Longitudinal radius","LP");
    // Leg[4]->AddEntry(gRadiusX,"Transverse radius","LP");
    // Leg[4]->AddEntry(gRadiusPsi,"Trapping radius","LP");
    Leg[4]->AddEntry(gRadiusZ,"l_{m}","LP");
    Leg[4]->AddEntry(gRadiusX,"r_{m}","LP");
    Leg[4]->AddEntry(gRadiusPsi,"l_{t}","LP");
 
    Leg[4]->Draw();
  }
  C->cd();
  
  // Output file
  TString fOutName = Form("./Blowout/Blowout-Scalings");
  fOutName += Form("_%i",time);
  
  PGlobals::imgconv(C,fOutName,opt);
  // ---------------------------------------------------------
  
  
}


