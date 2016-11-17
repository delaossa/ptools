#include <TH2F.h>
#include <TClonesArray.h>
#include <TPaletteAxis.h>
#include <TGraphSmooth.h>

#include "PGlobals.hh"
#include "PData.hh"
#include "PPalette.hh"


Double_t DenGauss(Double_t z, Double_t sigmal, Double_t ntop) {
  
  Double_t nnorm =  1.0 + (ntop - 1.0) * TMath::Exp(-(z*z)/(2*sigmal*sigmal)) ;
  if(z<0) nnorm = ntop;
  return nnorm;
}

Double_t DerDenGauss(Double_t z, Double_t sigmal, Double_t ntop) {
  
  Double_t dernnorm =  -(z/(sigmal*sigmal)) * (ntop - 1) * TMath::Exp(-(z*z)/(2*sigmal*sigmal));
  if(z<0) dernnorm = 0.0;
  return dernnorm;
}

Double_t DenLin(Double_t z, Double_t sigmal, Double_t ntop) {

  Double_t zmin  = (2.0 - TMath::Sqrt(TMath::Exp(1))) * sigmal;
  Double_t zmax  = 2.0 * sigmal;
  
  Double_t nnorm = ntop;
  Double_t slop  = (1-ntop)/(TMath::Sqrt(TMath::Exp(1)) * sigmal);
  
  if(z>=zmin && z<=zmax)
    nnorm = 1 + slop * ( z - sigmal ) - slop * sigmal;
  else if(z>zmax) nnorm = 1.0;
  
  return nnorm;
  
}

Double_t DerDenLin(Double_t z, Double_t sigmal, Double_t ntop) {

  Double_t zmin  = (2.0 - TMath::Sqrt(TMath::Exp(1))) * sigmal;
  Double_t zmax  = 2.0 * sigmal;
  
  Double_t dernnorm = 0.0;
  Double_t slop  = (1-ntop)/(TMath::Sqrt(TMath::Exp(1)) * sigmal);
  
  if(z>=zmin && z<=zmax) {
    dernnorm = slop;
  }
  
  return dernnorm;
  
}

Double_t betaph(Double_t z, Double_t phase, Double_t sigmal, Double_t ntop) {

  Double_t beta = 1.0 / ( 1.0 + (phase/2) * TMath::Power(DenGauss(z,sigmal,ntop),-3/2.) * DerDenGauss(z,sigmal,ntop)  ) ;
  return beta;
    
}

Double_t betaphlin(Double_t z, Double_t phase, Double_t sigmal, Double_t ntop) {

  Double_t beta = 1.0 / ( 1.0 + (phase/2) * TMath::Power(DenLin(z,sigmal,ntop),-3/2.) * DerDenLin(z,sigmal,ntop)  ) ;
  return beta;
    
}


Double_t betapsi(Double_t psi) {
  Double_t beta = - ( 2*psi + psi*psi ) / (2*psi + psi*psi + 2) ;
  return beta;
}

void PlotDDR2D( const TString &options="" ){
  
  PGlobals::Initialize();

  gStyle->SetNumberContours(255);

  Int_t font = 43;
  Int_t labelsize = 34;
  Int_t titlesize = 38;
  Float_t labeloffset = 0.01;
  Float_t titleoffset = 1.0;
  Float_t ticksize = 0.01;

  
  gStyle->SetTitleFont(43,"xyz");
  gStyle->SetLabelFont(43,"xyz");
  
  gStyle->SetLabelSize(labelsize,"xyz");
  gStyle->SetTitleSize(titlesize,"xyz");
    
  gStyle->SetLabelOffset(labeloffset,"xyz");
  gStyle->SetTitleOffset(titleoffset,"yz");
  gStyle->SetTitleOffset(titleoffset,"x");
  
  gStyle->SetTickLength(ticksize,"yz");
  gStyle->SetTickLength(ticksize+0.01,"x");

  Int_t frameWidth = 3;
  gStyle->SetLineWidth(frameWidth);  
  gStyle->SetNdivisions(006,"xyz");

  Int_t ndivx = 6;
  Int_t ndivy = 6;
  
  // Colors
  Int_t myBlue = TColor::GetColor((Float_t) 0.16, (Float_t) 0.83, (Float_t) 0.5);
  Int_t myNaranja = TColor::GetColor((Float_t) 0.992157, (Float_t) 0.411765, (Float_t) 0.027451);

  //  Int_t colorDen = myNaranja;
  //  Int_t colorBeta = kBlack;
  Int_t colorDen = kBlack;
  Int_t colorDenLin = kGray+2;
  Int_t colorBeta = kOrange+10;
  //  Int_t colorBetaLin = TColor::GetColor("#FCBFC1");
  Int_t colorBetaLin = kOrange+10;
  Int_t lineWidth = 6;
  Int_t lineWidthLin = 2;
  Int_t lineStyleLin = 8;

  // 
  Int_t ylabelsize = 16;
  Int_t ytitlesize = 20;
  Float_t ylabeloffset = 0.016;
  Float_t ytitleoffset = 1.0;
  Float_t yticksize = 0.02;
  
  // DDR 1D plot
  Int_t sizex = 300;
  Int_t sizey = 400;
    
  // Canvas setup
  TCanvas *CDDR = new TCanvas("CDDR","DDR vs phase velocity",sizex,sizey);
  CDDR->cd(0);
 
  gPad->SetTickx(0);
  gPad->SetTicky(0);
  gPad->SetFillStyle(4000);
  gPad->SetFrameFillStyle(4000);

  gPad->SetRightMargin(0.18);
  gPad->SetLeftMargin(0.18);
  

  Int_t Np = 200;
  Double_t *zarray = new Double_t[Np];
  Double_t *denarray = new Double_t[Np];
  Double_t *denarraylin = new Double_t[Np];
  Double_t *betapharray = new Double_t[Np];
  Double_t *betapharray2 = new Double_t[Np];
  Double_t *betapharraylin = new Double_t[Np];
  Double_t ntop   = 10;
  Double_t sigmal = 2.0;
  Double_t phase = -TMath::TwoPi();

  Float_t dMin = -0.5;
  Float_t dMax = ntop+0.5;

  //  Float_t zmin = - 0.49 * sigmal;
  Float_t zmin = - 0.0;
  Float_t zmax =  5 * sigmal;
  
  for(Int_t i=0; i<Np; i++) {
    zarray[i] = (i + 0.5) * (zmax-zmin)/Np + zmin;
    denarray[i] = DenGauss(zarray[i],sigmal,ntop);
    denarraylin[i] = DenLin(zarray[i],sigmal,ntop);
    betapharray[i] = betaph(zarray[i],phase,sigmal,ntop);
    betapharray2[i] = betaph(zarray[i],2.0*phase,sigmal,ntop);
    betapharraylin[i] = betaphlin(zarray[i],phase,sigmal,ntop);
  }
    
  TGraph *denvsz = new TGraph(Np,zarray,denarray);
  denvsz->SetLineColor(colorDen);
  denvsz->SetLineWidth(lineWidth);

  TGraph *denlinvsz = new TGraph(Np,zarray,denarraylin);
  denlinvsz->SetLineColor(colorDenLin);
  denlinvsz->SetLineWidth(lineWidthLin);
  denlinvsz->SetLineStyle(lineStyleLin);
    
    
  TH1F *hFrameD = new TH1F("hFrameD","",100,zmin,zmax);
  hFrameD->GetYaxis()->SetRangeUser(dMin,dMax);
  hFrameD->GetXaxis()->SetTitle("k_{p}^{0} z");
  hFrameD->GetYaxis()->SetTitle("n/n_{0}");

  // PGlobals::SetH1LabelSize(hFrameD);

  hFrameD->GetXaxis()->SetNdivisions(ndivx);
  hFrameD->GetXaxis()->SetLabelSize(ylabelsize);
  hFrameD->GetXaxis()->SetTitleSize(ytitlesize);
  hFrameD->GetXaxis()->CenterTitle();
  hFrameD->GetYaxis()->SetNdivisions(ndivy);
  hFrameD->GetYaxis()->SetAxisColor(colorDen);
  hFrameD->GetYaxis()->SetLabelColor(colorDen);
  hFrameD->GetYaxis()->SetTitleColor(colorDen);
  hFrameD->GetYaxis()->SetLabelSize(ylabelsize);
  hFrameD->GetYaxis()->SetLabelOffset(ylabeloffset);
  hFrameD->GetYaxis()->SetTitleSize(ytitlesize);
  hFrameD->GetYaxis()->SetTitleOffset(ytitleoffset);
  hFrameD->GetYaxis()->SetTickLength(yticksize);
  hFrameD->GetYaxis()->CenterTitle();
    
  hFrameD->Draw("axis");

  gPad->Update();

  // Draw velocities
  Double_t rightmin = -0.05;
  Double_t rightmax = 1.05;
  Double_t slope = (gPad->GetUymax() - gPad->GetUymin())/(rightmax-rightmin); 
  for(Int_t i=0; i<Np; i++) {
    betapharray[i] = (betapharray[i]-rightmin) * slope +  gPad->GetUymin();
    betapharray2[i] = (betapharray2[i]-rightmin) * slope +  gPad->GetUymin();
    betapharraylin[i] = (betapharraylin[i]-rightmin) * slope +  gPad->GetUymin();
  }
    
  TGraph *betaphvsz = new TGraph(Np,zarray,betapharray); 
  //   betaphvsz->SetLineColor(myBlue);
  betaphvsz->SetLineColor(colorBeta);
  betaphvsz->SetLineWidth(lineWidth);
  betaphvsz->SetLineStyle(1);

  TGraph *betaphvsz2 = new TGraph(Np,zarray,betapharray2); 
  //   betaphvsz2->SetLineColor(myBlue);
  betaphvsz2->SetLineColor(colorBeta);
  betaphvsz2->SetLineWidth(lineWidth);
  betaphvsz2->SetLineStyle(2);
// betaphvsz2->Draw("C");

  TGraph *betaphvszlin = new TGraph(Np,zarray,betapharraylin); 
  betaphvszlin->SetLineColor(colorBetaLin);
  betaphvszlin->SetLineWidth(lineWidthLin);
  betaphvszlin->SetLineStyle(lineStyleLin);


  betaphvszlin->Draw("L");
  denlinvsz->Draw("L");
  denvsz->Draw("C");
    
  //draw an axis on the right side
  TGaxis *axisb = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),gPad->GetUxmax(),
			     gPad->GetUymax(),rightmin,rightmax,ndivy,"+LS");
  
  axisb->SetLineWidth(frameWidth);
  axisb->SetLineColor(betaphvsz->GetLineColor());
  axisb->SetLabelColor(betaphvsz->GetLineColor());
  axisb->SetTitle("#beta");
  axisb->CenterTitle();
  axisb->SetTitleColor(betaphvsz->GetLineColor());
  axisb->SetLabelFont(font);
  axisb->SetLabelColor(colorBeta);
  axisb->SetLabelSize(ylabelsize);
  axisb->SetLabelOffset(ylabeloffset);
  axisb->SetTitleFont(font);   
  axisb->SetTitleSize(ytitlesize);
  axisb->SetTitleOffset(ytitleoffset);
  axisb->SetTickSize(yticksize);
  axisb->Draw();

  TBox *lFrame = new TBox(gPad->GetUxmin(), gPad->GetUymin(),
			  gPad->GetUxmax(), gPad->GetUymax());
  lFrame->SetFillStyle(0);
  lFrame->SetLineColor(PGlobals::frameColor);
  lFrame->SetLineWidth(lineWidth);
  lFrame->Draw();

  TLine xaxis(gPad->GetUxmin(), 0.0,
	      gPad->GetUxmax(), 0.0);
  xaxis.SetLineColor(kGray+2);
  xaxis.SetLineWidth(2);
  xaxis.SetLineStyle(2);
  xaxis.Draw();
    
  TLine laxis(gPad->GetUxmin(), gPad->GetUymin(),
	      gPad->GetUxmin(), gPad->GetUymax());
  laxis.SetLineColor(colorDen);
  laxis.SetLineWidth(lineWidth);
  laxis.Draw();
    
  TLine raxis(gPad->GetUxmax(), gPad->GetUymin(),
	      gPad->GetUxmax(), gPad->GetUymax());
  raxis.SetLineColor(colorBeta);
  raxis.SetLineWidth(lineWidth);
  raxis.Draw();
    
  gPad->Update();
  gPad->RedrawAxis("G");

  //betaphvsz2->Draw("C");
  betaphvsz->Draw("C");

  // Get Data from PIC simulations with test ramp
  // --------------------------------------------
  
  Float_t ntop0 = 15.63;
  Float_t nres0 = 2.98;
  Float_t sigmal0 = 68.49;

  // const Int_t Nsim = 5;
  // char sName[Nsim][36] = {"flash_v1kA.200pC.G.DDR-long2s.3D",
  // 			  "flash_v2.0kA.400pC.G.DDR-long2s.3D",
  // 			  "flash_v2.5kA.500pC.G.DDR-long2s.3D",
  // 			  "flash_v3.0kA.600pC.G.DDR-long2s.3D",
  // 			  "flash_v5.0kA.1000pC.G.DDR-long2s.3D"};
  
  // char lName[Nsim][36] = {"1.0 kA",
  // 			  "2.0 kA",
  // 			  "2.5 kA", 
  // 			  "3.0 kA", 
  // 			  "5.0 kA"};

  const Int_t Nsim = 3;
  char sName[Nsim][36] = {"flash_v1kA.200pC.G.DDR-long2s.3D",
			  "flash_v2.5kA.500pC.G.DDR-long2s.3D",
			  "flash_v5.0kA.1000pC.G.DDR-long2s.3D"};
  
  char lName[Nsim][36] = {"1.0 kA",
			  "2.5 kA", 
			  "5.0 kA"};

  // simualtions colors
  Int_t color[Nsim] = {};
  Int_t color2[Nsim] = {};

  // color[0] =  TColor::GetColor("#000054");
  // color[1] =  TColor::GetColor("#000089");
  // color[2] =  TColor::GetColor("#0000C0");
  color2[0] = kBlue-9;
  color2[1] = kBlue-6;
  color2[2] = kBlue-2;
  // color[3] = kBlue-3;
  // color[4] = kBlue-4;

  color[0] = kRed-9;
  color[1] = kRed-6;
  color[2] = kRed-2;
  // color2[3] = kRed-3;
  // color2[4] = kRed-4;
  
  
  const Int_t NRGBs = 2;
  const Int_t NCont = 64;
  Double_t Stops[NRGBs] = { 0.000, 1.00 };
  Int_t cindex[NRGBs];
  Int_t cindex2[NRGBs];
  PPalette *palmin = new PPalette("palmin");
  PPalette *palmax = new PPalette("palmax");
  cindex[0] = kBlue-10;
  cindex[1] = kBlue+2;
  palmin->CreateGradientColorTable(NRGBs, Stops, cindex, NCont);
  cindex2[0] = kRed-10;
  cindex2[1] = kRed+2;
  palmax->CreateGradientColorTable(NRGBs, Stops, cindex2, NCont);

  // Double_t Red[NRGBs] =   { 0.698, 0.106 };
  // Double_t Green[NRGBs] = { 0.818, 0.078 };
  // Double_t Blue[NRGBs] =  { 0.880, 0.518 };
  // palmin->CreateGradientColorTable(NRGBs, Stops, Red, Green, Blue, NCont);
  
  // Double_t Red2[NRGBs] =   { 0.965, 0.518 };
  // Double_t Green2[NRGBs] = { 0.925, 0.078 };
  // Double_t Blue2[NRGBs] =  { 0.353, 0.106 };
  // palmax->CreateGradientColorTable(NRGBs, Stops, Red2, Green2, Blue2, NCont);
  
  // for(Int_t i=0;i<Nsim;i++) {
  //   Int_t index = ((palmin->GetNColors()-2)/(Nsim-1)) * i + 1;
  //   color[i] = palmin->GetColorIndex(index);
  //   index =  ((palmax->GetNColors()-2)/(Nsim-1)) * i + 1;
  //   color2[i] = palmax->GetColorIndex(index);
  // }
  
  TFile *sFile[Nsim];

  TH2F **hVvsTime = new TH2F*[Nsim];
  TH2F **hEzvsTime = new TH2F*[Nsim];
  TH2F **hDenvsTime = new TH2F*[Nsim];
  TH2F **hJzvsTime = new TH2F*[Nsim];
  TH1F **hV1D = new TH1F*[Nsim];
  TH1F **hEz1D = new TH1F*[Nsim];
  TH1F **hDen1D = new TH1F*[Nsim];
  TH1F **hJz1D = new TH1F*[Nsim];
  TH1F **hBeta1D = new TH1F*[Nsim];
  TGraph ***gVcross = new TGraph**[Nsim];
  TGraph ***gVextr = new TGraph**[Nsim];
  TGraph ***gEzcross = new TGraph**[Nsim];
  TGraph ***gEzextr = new TGraph**[Nsim];
  TGraph **gBeta = new TGraph*[Nsim];

  TGraph **gBetasimvsz = new TGraph*[Nsim];
  TGraph **gBetasimvszcorr = new TGraph*[Nsim];
  TGraph **gBetasimvszsmooth = new TGraph*[Nsim];

  Float_t *betamax = new Float_t[Nsim];

  // Kernel Smoother
  // create new kernel smoother and smooth data with bandwidth = 2.0
  TGraphSmooth **gsmooth = new TGraphSmooth*[Nsim];
     
  // Load first simulation data (for instance)
  PData *pData = PData::Get(sName[0]);
  //if(!pData->IsInit()) return;

  for(Int_t i=0;i<Nsim;i++) {

    gsmooth[i] = new TGraphSmooth("normal");
     
    betamax[i] = -9999.0;

    TString filename = Form("./%s/Plots/Evolutions/Evolutions-%s.root",sName[i],sName[i]);
    sFile[i] = new TFile(filename,"READ");

    // Get evolution of the wakefield potential (on axis)
    hVvsTime[i] = (TH2F*) sFile[i]->Get("hVvsTime");
    hEzvsTime[i] = (TH2F*) sFile[i]->Get("hEvsTime_0");
    hDenvsTime[i] = (TH2F*) sFile[i]->Get("hDenvsTime_0");
    hJzvsTime[i] = (TH2F*) sFile[i]->Get("hJzvsTime_0");
    
    Int_t NTBins = hVvsTime[i]->GetNbinsX();
    for(Int_t it=NTBins;it>0;it--) {
      
      // 1D field at timestep "it".
      hV1D[i] = (TH1F*) hVvsTime[i]->ProjectionY("_py",it,it);
      hEz1D[i] = (TH1F*) hEzvsTime[i]->ProjectionY("_py",it,it);
      hJz1D[i] = (TH1F*) hJzvsTime[i]->ProjectionY("_py",it,it);
      hDen1D[i] = (TH1F*) hDenvsTime[i]->ProjectionY("_py",it,it);
      
      // calculate the crossings
      Int_t MAXCROSS = 2;
      Float_t *VCross = new Float_t[MAXCROSS];
      Float_t *VExtr = new Float_t[MAXCROSS];
      
      memset(VCross,0,sizeof(Float_t)*MAXCROSS);
      memset(VExtr,0,sizeof(Float_t)*MAXCROSS);    
      Int_t NVCross = PGlobals::HCrossings(hV1D[i],VCross,VExtr,MAXCROSS,0.,0.);

      Float_t *EzCross = new Float_t[MAXCROSS];
      Float_t *EzExtr = new Float_t[MAXCROSS];
      
      memset(EzCross,0,sizeof(Float_t)*MAXCROSS);
      memset(EzExtr,0,sizeof(Float_t)*MAXCROSS);    
      Int_t NEzCross = PGlobals::HCrossings(hEz1D[i],EzCross,EzExtr,MAXCROSS,0.,0.);

      
      if(it==NTBins) {
	gVcross[i] = new TGraph*[NVCross];
	gVextr[i] = new TGraph*[NVCross]; 

	for(Int_t ic = 0;ic<NVCross;ic++) {
	  gVcross[i][ic] = new TGraph(NTBins);
	  gVcross[i][ic]->SetName(Form("gVcross_%i_%i",i,ic)); 
	  
	  gVextr[i][ic] = new TGraph(NTBins);
	  gVextr[i][ic]->SetName(Form("gVextr_%i_%i",i,ic)); 
	}

	gEzcross[i] = new TGraph*[NEzCross];
	gEzextr[i] = new TGraph*[NEzCross]; 
	
	for(Int_t ic = 0;ic<NEzCross;ic++) {
	  gEzcross[i][ic] = new TGraph(NTBins);
	  gEzcross[i][ic]->SetName(Form("gEzcross_%i_%i",i,ic)); 
	  
	  gEzextr[i][ic] = new TGraph(NTBins);
	  gEzextr[i][ic]->SetName(Form("gEzextr_%i_%i",i,ic)); 
	}

	gBeta[i] = new TGraph(NTBins);
	gBeta[i]->SetName(Form("gBeta_%i",i));
	
      }
      
      Float_t time = hVvsTime[i]->GetXaxis()->GetBinCenter(it);
      for(Int_t ic=0;ic<NVCross;ic++) {
	// cout << Form("  - Adding %i crossing: cross = %6.4f extreme = %6.4f",ic,VCross[ic],VExtr[ic]) << endl;
	
	gVcross[i][ic]->SetPoint(it-1,time,VCross[ic]);
	
	// Use the 2nd crossing of Ez as the position of the minimum \phi
	Double_t value = VExtr[ic];
	// if(EzCross[ic]) {
	//   value = (value + hV1D[i]->GetBinContent(hV1D[i]->FindBin(EzCross[ic])))/2.0;
	// }
	gVextr[i][ic]->SetPoint(it-1,time,value);
      }
      
      for(Int_t ic=0;ic<NEzCross;ic++) {
	gEzcross[i][ic]->SetPoint(it-1,time,EzCross[ic]);
	gEzextr[i][ic]->SetPoint(it-1,time,EzExtr[ic]);
	
      }

      Float_t jzmax = hJz1D[i]->GetMinimum();
      Int_t ijzmax = hJz1D[i]->GetMinimumBin();
      jzmax = (jzmax + hJz1D[i]->GetBinContent(ijzmax-1) + hJz1D[i]->GetBinContent(ijzmax+1) )/ 3.0;
      
      
      Float_t dmax =  (hDen1D[i]->GetBinContent(ijzmax-1)
		       + hDen1D[i]->GetBinContent(ijzmax)
		       + hDen1D[i]->GetBinContent(ijzmax+1) )/ 3.0;
      Float_t betamax = - jzmax/dmax; 
      gBeta[i]->SetPoint(it-1,time,betamax);
    }    
    
    // Draw beta from simulation 
    Double_t *betaarray = new Double_t[Np];
    Double_t *betaarraycorr = new Double_t[Np];
    Double_t *den0array = new Double_t[Np];
    Double_t nfactor = (ntop0 - 1 + TMath::Exp(2.)) / (ntop - 1 + TMath::Exp(2.));


    Int_t Npp = gVextr[i][1]->GetN();
    Double_t *psival = gVextr[i][1]->GetY();
    Double_t *zval = gVextr[i][1]->GetX();
    Double_t *nval = new Double_t[Npp];

    // Get density array from array of z (assuming a function n(z))
    for(Int_t j=0; j<Npp; j++) {
      Float_t density =  1 + (ntop0-1) * TMath::Exp(-(zval[j]*zval[j])/(2*sigmal0*sigmal0));
      nval[j] = density;
    }
    
    for(Int_t ip=0; ip<Np; ip++) {
      den0array[ip] = denarray[ip] * nfactor;
      Bool_t found = kFALSE;
      for(Int_t j=1; j<Npp; j++) {
	if( den0array[ip]<nval[j-1] && den0array[ip]>nval[j]) {
	  // interpolate
	  Double_t psi_int = (psival[j]-psival[j-1])/(nval[j]-nval[j-1]) * (den0array[ip] - nval[j-1]) + psival[j-1];

	  betaarray[ip] = betapsi(psi_int);
	  betaarraycorr[ip] = (betaarray[ip]-rightmin) * slope +  gPad->GetUymin();
	 
	  if(betaarray[ip]>betamax[i]) betamax[i] = betaarray[ip];
	  found = kTRUE;
      	}
	if(found) break;
      }
      if(found) continue;
      //  cout << Form("  %f ", den0array[ip]) << endl;
    }

    gBetasimvsz[i] = new TGraph(Np,zarray,betaarray); 
    gBetasimvsz[i]->SetLineColor(color[i]);
    gBetasimvsz[i]->SetLineWidth(lineWidth);
    gBetasimvsz[i]->SetLineStyle(1);
    // gBetasimvsz[i]->Draw("C");

    gBetasimvszcorr[i] = new TGraph(Np,zarray,betaarraycorr); 
    gBetasimvszcorr[i]->SetLineColor(color[i]);
    gBetasimvszcorr[i]->SetLineWidth(lineWidth);
    gBetasimvszcorr[i]->SetLineStyle(1);
    //    gBetasimvszcorr[i]->Draw("C");

    gBetasimvszsmooth[i] = gsmooth[i]->SmoothSuper(gBetasimvszcorr[i],"",0.0);
    gBetasimvszsmooth[i]->SetLineColor(color[i]);
    gBetasimvszsmooth[i]->SetLineWidth(lineWidth);
    gBetasimvszsmooth[i]->SetLineStyle(1);
    gBetasimvszsmooth[i]->Draw("C");
  }



  // draw last    
  //  betaphvsz->Draw("C");

  laxis.Draw();
  raxis.Draw();

  
  CDDR->cd(0);
  
  // Print to a file
  // Output file
  PGlobals::imgconv(CDDR,Form("./DDR/DDR-DenVsVel"),"pdf");
  // ---------------------------------------------------------


  PPalette *pal = new PPalette("pal");
  pal->SetPalette(kBird);

  
  // Canvas setup
  sizex = 1024;
  sizey = 640;
  TCanvas *C = new TCanvas("C","Minimum phase velocity",sizex,sizey);
  
  Float_t sigmin = 0.0;
  Float_t sigmax = 10.0;
  Float_t denmin = 1.0;
  Float_t denmax = 10.0;

  //  Float_t phase = -TMath::TwoPi();
  //  Int_t Np = 100;
  
  TH2F *hBetamin2D = new TH2F("hBetamin2D","",Np,sigmin,sigmax,Np,denmin,denmax);
  TH2F *hBetaminpos2D = new TH2F("hBetaminpos2D","",Np,sigmin,sigmax,Np,denmin,denmax);
  TGraph ***gBetavsz = new TGraph**[Np]; 
  
  // Double_t *zarray = new Double_t[Np];
  // Double_t *betapharray = new Double_t[Np];
  
  for(Int_t i=0; i<Np; i++) {

    gBetavsz[i] = new TGraph*[Np]; 
    Float_t sig  = (i + 0.5) * (sigmax-sigmin)/Np + sigmin;
    
    for(Int_t j=0; j<Np; j++) {
      
      Float_t ntop = (j + 0.5) * (denmax-denmin)/Np + denmin;
      
      Float_t betamin = 999.;
      Float_t betaminpos = 999.;
      
      for(Int_t k=0; k<Np; k++) {
	zarray[k] = (k + 0.5) * 5 * sig/Np;
	betapharray[k] = betaph(zarray[k],phase,sig,ntop);
	//	if(betapharray[k]<0.0) continue;
	if(betapharray[k]<betamin) {
	  betamin = betapharray[k];
	  betaminpos = zarray[k]/sig;
	}
      }
      hBetamin2D->SetBinContent(i+1,j+1,betamin); 
      hBetaminpos2D->SetBinContent(i+1,j+1,betaminpos); 
      gBetavsz[i][j] = new TGraph(Np,zarray,betapharray);
      //  gBetavsz[i][j]->Write();
    }
    
  }
  
  hBetamin2D->GetXaxis()->SetTitle("k_{p}^{0} #sigma_{l}");
  hBetamin2D->GetXaxis()->CenterTitle();
  hBetamin2D->GetYaxis()->SetTitle("n_{top}/n_{0}");
  hBetamin2D->GetYaxis()->CenterTitle();
  hBetamin2D->GetZaxis()->SetTitle("#beta_{#chi,min}");
  hBetamin2D->GetZaxis()->CenterTitle();

  TH2F *hClone = (TH2F*) hBetamin2D->Clone("clone");
  
  const Int_t Ncontours = 10;
  Double_t contours[Ncontours];
  for(Int_t i=0; i<Ncontours; i++) {
    contours[i] = (i+1) * 0.1;
  }

  hClone->SetContour(Ncontours,contours);
  hClone->Draw("cont list 0");

  C->Update();
  TObjArray *conts = (TObjArray*) gROOT->GetListOfSpecials()->FindObject("contours");
  TClonesArray graphs("TGraph",Ncontours);
  Int_t nGraphs = 0;
  TGraph *gr = NULL;
  TList* clist = NULL;
  Int_t ncont = conts->GetSize();
  for(Int_t i = 0; i < ncont; i++){  
    clist = (TList*) conts->At(i);
    for(Int_t j = 0 ; j < clist->GetSize(); j++) {
      gr = (TGraph*) clist->At(j);
      if(!gr) continue;

      gr->SetLineWidth(3);
      gr->SetLineStyle(7);
      gr->SetLineColor(kWhite);
      //gr->SetLineColor(kGray);
      // if(i==4) {
      // 	gr->SetLineColor(kRed);
      // 	gr->SetLineWidth(2);
      // } 
      
      new(graphs[nGraphs]) TGraph(*gr) ;
      nGraphs++;
    }
  }

  // Simulation contours
  Double_t contsim[Nsim];
  for(Int_t i=0; i<Nsim; i++) {
    contsim[i] = betamax[i];
    //   cout << Form("sim %i:  beta_max = %.4f",i,contsim[i]) << endl;
  }
  
  hClone->SetContour(Nsim,contsim);
  hClone->Draw("cont list 0");

  C->Update();
  TObjArray *contssim = (TObjArray*) gROOT->GetListOfSpecials()->FindObject("contours");
  TClonesArray graphssim("TGraph",Nsim);
  nGraphs = 0;
  gr = NULL;
  clist = NULL;
  ncont = contssim->GetSize();
  for(Int_t i = 0; i < ncont; i++){  
    clist = (TList*) contssim->At(i);
    for(Int_t j = 0 ; j < clist->GetSize(); j++) {
      gr = (TGraph*) clist->At(j);
      if(!gr) continue;
      
      gr->SetLineWidth(lineWidth);
      gr->SetLineStyle(1);
      gr->SetLineColor(color[i]);
      //gr->SetLineColor(myNaranja);
      
      new(graphssim[nGraphs]) TGraph(*gr) ;
      nGraphs++;
    }
  }
  

  
  hBetamin2D->Draw("colz");
  
  // Contours  
  for(Int_t i=0;i<graphs.GetEntriesFast();i++) {
    TGraph *gr = (TGraph*) graphs.At(i);
    if(!gr) continue;
    gr->Draw("C");
  }

  // Simulations contours  
  for(Int_t i=0;i<graphssim.GetEntriesFast();i++) {
    TGraph *gr = (TGraph*) graphssim.At(i);
    if(!gr) continue;
    gr->Draw("C");
  }
  
  gPad->Update();
  
  TPaletteAxis *palette = (TPaletteAxis*)hBetamin2D->GetListOfFunctions()->FindObject("palette");
  if(palette) {
    Double_t y1 = gPad->GetBottomMargin();
    Double_t y2 = 1 - gPad->GetTopMargin();
    Double_t x1 = 1 - gPad->GetRightMargin();
    
    Double_t x1b = x1 + 0.01;
    Double_t x2b = x1 + 0.035;
    Double_t y1b = y1 + 0.03;
    Double_t y2b = y2 - 0.03;
    palette->SetX1NDC(x1b);
    palette->SetY1NDC(y1b);
    palette->SetX2NDC(x2b);
    palette->SetY2NDC(y2b);
    palette->SetBorderSize(2);
    palette->SetLineColor(1);

    TPave *pFrame = new TPave(x1b,y1b,x2b,y2b,1,"NDCL");
    pFrame->SetFillStyle(0);
    pFrame->SetLineColor(kBlack);
    pFrame->SetLineWidth(frameWidth);
    pFrame->SetShadowColor(0);
    pFrame->Draw();
  }

  TBox *lFrame1 = new TBox(gPad->GetUxmin(), gPad->GetUymin(),
			   gPad->GetUxmax(), gPad->GetUymax());
  lFrame1->SetFillStyle(0);
  lFrame1->SetLineColor(kBlack);
  lFrame1->SetLineWidth(frameWidth);
  lFrame1->Draw();
  
  gPad->RedrawAxis(); 

  PGlobals::imgconv(C,Form("./DDR/DDR2D-BetaMin"),"pdf");

  // ------------------------
  
  TCanvas *C2 = new TCanvas("C2","Minimum phase velocity position",sizex,sizey);

  hBetaminpos2D->GetXaxis()->SetTitle("k_{p}^{0} #sigma_{l}");
  hBetaminpos2D->GetXaxis()->CenterTitle();
  hBetaminpos2D->GetYaxis()->SetTitle("n_{top}/n_{0}");
  hBetaminpos2D->GetYaxis()->CenterTitle();
  hBetaminpos2D->GetZaxis()->SetTitle("z_{min}/#sigma_{l}");
  hBetaminpos2D->GetZaxis()->CenterTitle();

  //  hBetaminpos2D->GetZaxis()->SetRangeUser(1.0,2.0);
  
  TH2F *hClone2 = (TH2F*) hBetaminpos2D->Clone("clone2");
  
  const Int_t Ncontours2 = 10;
  Double_t contours2[Ncontours2];
  for(Int_t i=0; i<Ncontours2; i++) {
    contours2[i] = (i+1) * 0.1 + 1.0;
  }

  hClone2->SetContour(Ncontours2,contours2);
  hClone2->Draw("cont list 0");

  C2->Update();
  TObjArray *conts2 = (TObjArray*) gROOT->GetListOfSpecials()->FindObject("contours");
  TClonesArray graphs2("TGraph",Ncontours2);
  Int_t nGraphs2 = 0;
  TGraph *gr2 = NULL;
  TList* clist2 = NULL;
  Int_t ncont2 = conts2->GetSize();
  for(Int_t i = 0; i < ncont2; i++){  
    clist2 = (TList*) conts2->At(i);
    for(Int_t j = 0 ; j < clist2->GetSize(); j++) {
      gr2 = (TGraph*) clist2->At(j);
      if(!gr2) continue;

      gr2->SetLineWidth(3);
      gr2->SetLineStyle(7);
      gr2->SetLineColor(kWhite);
      // if(i==4) {
      // 	gr2->SetLineColor(kRed);
      // 	gr2->SetLineWidth(2);
      // } 
      
      new(graphs2[nGraphs2]) TGraph(*gr2) ;
      nGraphs2++;
    }
  }

  hBetaminpos2D->Draw("colz");
  
  // Contours  
  for(Int_t i=0;i<graphs2.GetEntriesFast();i++) {
    TGraph *gr = (TGraph*) graphs2.At(i);
    if(!gr) continue;
    gr->Draw("C");
  }

  gPad->Update();

  TPaletteAxis *palette2 = (TPaletteAxis*)hBetaminpos2D->GetListOfFunctions()->FindObject("palette");
  if(palette2) {
    Double_t y1 = gPad->GetBottomMargin();
    Double_t y2 = 1 - gPad->GetTopMargin();
    Double_t x1 = 1 - gPad->GetRightMargin();
    
    Double_t x1b = x1 + 0.01;
    Double_t x2b = x1 + 0.035;
    Double_t y1b = y1 + 0.03;
    Double_t y2b = y2 - 0.03;
    palette2->SetX1NDC(x1b);
    palette2->SetY1NDC(y1b);
    palette2->SetX2NDC(x2b);
    palette2->SetY2NDC(y2b);
    palette2->SetBorderSize(2);
    palette2->SetLineColor(1);

    TPave *pFrame = new TPave(x1b,y1b,x2b,y2b,1,"NDCL");
    pFrame->SetFillStyle(0);
    pFrame->SetLineColor(kBlack);
    pFrame->SetShadowColor(0);
    pFrame->Draw();
  }

  TBox *lFrame2 = new TBox(gPad->GetUxmin(), gPad->GetUymin(),
			   gPad->GetUxmax(), gPad->GetUymax());
  lFrame2->SetFillStyle(0);
  lFrame2->SetLineColor(kBlack);
  lFrame2->SetLineWidth(frameWidth);
  lFrame2->Draw();
  
  gPad->RedrawAxis();

  PGlobals::imgconv(C2,Form("./DDR/DDR2D-BetaMinPos"),"pdf");


  // ------------------------
  sizex = 1024;
  sizey = 320;

  TCanvas *C3 = new TCanvas("C3","Minimum phase velocity position",sizex,sizey);
  C3->cd();

  gPad->SetTopMargin(0.05);
  gPad->SetRightMargin(0.20);
  gPad->SetBottomMargin(0.30);
  gPad->SetLeftMargin(0.15);
  
  // gPad->SetTickx(1);
  // gPad->SetTicky(1);
  gPad->SetFillStyle(4000);
  gPad->SetFrameFillStyle(4000);

  Int_t NBins = hBetaminpos2D->GetNbinsX();
  TH1F *hBetaminpos1D = (TH1F*) hBetaminpos2D->ProjectionY("_py",NBins/2.0,NBins/2.0);

  hBetaminpos1D->GetXaxis()->SetTitle("n_{top}/n_{0}");
  hBetaminpos1D->GetXaxis()->CenterTitle();
  hBetaminpos1D->GetXaxis()->SetTickLength(0.04);
  hBetaminpos1D->GetYaxis()->SetTitle("z_{min}/#sigma_{l}");
  hBetaminpos1D->GetYaxis()->CenterTitle();
  hBetaminpos1D->GetYaxis()->SetTitleOffset(0.5);

  hBetaminpos1D->Draw("axis");
  
  Int_t Npp = hBetaminpos1D->GetXaxis()->GetNbins();
  Float_t *xval = new Float_t[Npp];
  Float_t *yval = new Float_t[Npp];
  for(Int_t j=0;j<Npp;j++) {
    xval[j] = hBetaminpos1D->GetBinCenter(j+1);
    yval[j] = hBetaminpos1D->GetBinContent(j+1);
  }
  
  TGraph *gBetaminpos1D = new TGraph(Npp,xval,yval);

  TGraphSmooth *gsm = new TGraphSmooth("normal");
  TGraph *gBetaminpos1Dsmooth = gsm->SmoothSuper(gBetaminpos1D,"",0.0);
  
  gBetaminpos1Dsmooth->SetLineWidth(lineWidth);
  gBetaminpos1Dsmooth->SetLineColor(kGray+3);
  gBetaminpos1Dsmooth->Draw("L");
  
  gPad->Update();

  TBox *lFrame3 = new TBox(gPad->GetUxmin(), gPad->GetUymin(),
			   gPad->GetUxmax(), gPad->GetUymax());
  lFrame3->SetFillStyle(0);
  lFrame3->SetLineColor(kBlack);
  lFrame3->SetLineWidth(frameWidth);
  lFrame3->Draw();
  
  gPad->RedrawAxis();

  PGlobals::imgconv(C3,Form("./DDR/DDR1D-BetaMinPos"),"pdf");
  

  // Psi potential vs time (2D)
  // For one of the simulations  
  sizex = 1024;
  sizey = 640;
  gStyle->SetLineWidth(2);  
  Int_t isim = 1;

  // Canvas setup
  TCanvas *CV = new TCanvas("CV","Evolution of the wakefield potential",sizex,sizey);
  
  CV->cd(0);
  gPad->SetFrameLineWidth(2);  
  gPad->SetTickx(0);
  gPad->SetTicky(0);
  
  //  PGlobals::SetH1LabelSize(hVvsTime[isim]);
  hVvsTime[isim]->GetXaxis()->SetNdivisions(ndivx);
  hVvsTime[isim]->SetLabelSize(labelsize,"xyz");
  hVvsTime[isim]->SetTitleSize(titlesize,"xyz");
  hVvsTime[isim]->SetLabelOffset(labeloffset,"xyz");
  hVvsTime[isim]->SetTitleOffset(titleoffset,"yz");
  hVvsTime[isim]->SetTitleOffset(titleoffset+0.1,"x");
  
  hVvsTime[isim]->SetTickLength(ticksize,"yz");
  hVvsTime[isim]->SetTickLength(ticksize+0.01,"x");

  // Change the range of z axis for the fields to be symmetric.
  Float_t Vmax = hVvsTime[isim]->GetMaximum();
  Float_t Vmin = hVvsTime[isim]->GetMinimum();
  if(Vmax > TMath::Abs(Vmin))
    Vmin = -Vmax;
  else
    Vmax = -Vmin;

  
  hVvsTime[isim]->GetZaxis()->SetRangeUser(Vmin,Vmax);

  pal->SetPalette("rbow0");
  
  hVvsTime[isim]->Draw("colz");
  
  gPad->Update();
  TPaletteAxis *palette3 = (TPaletteAxis*)hVvsTime[isim]->GetListOfFunctions()->FindObject("palette");
  if(palette3) {
    Double_t y1 = gPad->GetBottomMargin();
    Double_t y2 = 1 - gPad->GetTopMargin();
    Double_t x1 = 1 - gPad->GetRightMargin();
    
    Double_t x1b = x1 + 0.01;
    Double_t x2b = x1 + 0.035;
    Double_t y1b = y1 + 0.03;
    Double_t y2b = y2 - 0.03;
    palette3->SetX1NDC(x1b);
    palette3->SetY1NDC(y1b);
    palette3->SetX2NDC(x2b);
    palette3->SetY2NDC(y2b);
    palette3->SetBorderSize(2);
    palette3->SetLineColor(1);

    TPave *pFrame = new TPave(x1b,y1b,x2b,y2b,1,"NDCL");
    pFrame->SetFillStyle(0);
    pFrame->SetLineColor(kBlack);
    pFrame->SetShadowColor(0);
    pFrame->Draw();
  }

  TBox *lFrame4 = new TBox(gPad->GetUxmin(), gPad->GetUymin(),
			   gPad->GetUxmax(), gPad->GetUymax());
  lFrame4->SetFillStyle(0);
  lFrame4->SetLineColor(kBlack);
  lFrame4->SetLineWidth(2);
  lFrame4->Draw();
  
  TGraph *gChiminvsz = gEzcross[isim][1];
  gChiminvsz->SetLineWidth(4);
  gChiminvsz->SetLineColor(myNaranja);
  gChiminvsz->Draw("C");

  TGraph *gChimaxvsz = gEzcross[isim][0];
  gChimaxvsz->SetLineWidth(4);
  gChimaxvsz->SetLineColor(myBlue);
  gChimaxvsz->Draw("C");

  gVcross[isim][0]->SetLineWidth(3);
  gVcross[isim][0]->SetLineColor(kGray+2);
  gVcross[isim][0]->SetLineStyle(7);
  gVcross[isim][0]->Draw("C");
  
  gVcross[isim][1]->SetLineWidth(3);
  gVcross[isim][1]->SetLineColor(kGray+2);
  gVcross[isim][1]->SetLineStyle(7);
  gVcross[isim][1]->Draw("C");
  
  gPad->RedrawAxis();

  // Print to a file
  // Output file
  PGlobals::imgconv(CV,"./DDR/DDR2D-PsiVsZ","pdf");
  // ---------------------------------------------------------

  
  TLine xaxis2(gPad->GetUxmin(), nres0,
	       gPad->GetUxmax(), nres0);
  xaxis2.SetLineColor(kBlack);
  xaxis2.SetLineWidth(3);
  xaxis2.SetLineStyle(3);

  TLine yaxis2(2*sigmal0,gPad->GetUymin(),
	       2*sigmal0,gPad->GetUymax());
  yaxis2.SetLineColor(kBlack);
  yaxis2.SetLineWidth(3);
  yaxis2.SetLineStyle(3);
    
  {
    // Canvas setup
    sizex = 1024;
    sizey = 640;
    TCanvas *CVR = new TCanvas("CVR","Evolution of the wakefield potential",sizex,sizey);

    const Int_t NPad = 2;
    TPad *pad[NPad];
    TH1F *hFrameP[NPad];
    
    Float_t lMargin = gStyle->GetPadLeftMargin();
    Float_t rMargin = gStyle->GetPadRightMargin();
    Float_t bMargin = gStyle->GetPadBottomMargin();
    Float_t tMargin = gStyle->GetPadTopMargin();
    Float_t pfactor = 0.2;    
    PGlobals::CanvasAsymPartition(C,NPad,lMargin,rMargin,bMargin,tMargin,pfactor,0.02);

    for(Int_t i=0;i<NPad;i++) {
      char name[16];
      sprintf(name,"pad_%i",i);
      pad[i] = (TPad*) gROOT->FindObject(name);
      pad[i]->SetFrameLineWidth(2);  
      pad[i]->SetTickx(0);
      pad[i]->SetTicky(0);

      sprintf(name,"hFrameP2_%i",i);

      hFrameP[i] = new TH1F(name,"",10,hVvsTime[isim]->GetXaxis()->GetXmin(),hVvsTime[isim]->GetXaxis()->GetXmax());
    
      Float_t xFactor = pad[0]->GetAbsWNDC()/pad[i]->GetAbsWNDC();
      Float_t yFactor = pad[0]->GetAbsHNDC()/pad[i]->GetAbsHNDC();

      hFrameP[i]->GetXaxis()->SetNdivisions(ndivx);
      hFrameP[i]->SetLabelSize(labelsize,"xyz");
      hFrameP[i]->SetTitleSize(titlesize,"xyz");
      hFrameP[i]->SetLabelOffset(labeloffset,"xyz");
      hFrameP[i]->SetTitleOffset(titleoffset,"yz");
      hFrameP[i]->SetTitleOffset(titleoffset+0.3,"x");
      
      if(i==NPad-1) hFrameP[i]->SetTickLength(ticksize*xFactor/(yFactor*pfactor),"yz");
      else  hFrameP[i]->SetTickLength(ticksize*xFactor/yFactor,"yz");
      
      hFrameP[i]->SetTickLength((ticksize+0.01)*yFactor/xFactor,"x");
      
      if(i!=0) {
	hFrameP[i]->GetXaxis()->SetLabelOffset(999);
	hFrameP[i]->GetXaxis()->SetTitleOffset(999);
      }
      
    }

    //----------------------------------  
    CVR->cd(0);
    pad[1]->Draw();
    pad[1]->cd();
    
    hFrameP[1]->GetXaxis()->SetTitle("k_{p}^{0} z");
    hFrameP[1]->GetXaxis()->CenterTitle();

    hFrameP[1]->GetYaxis()->SetTitle("n/n_{0}");
    hFrameP[1]->GetYaxis()->CenterTitle();
    hFrameP[1]->GetYaxis()->SetNdivisions(2);
 
    hFrameP[1]->GetYaxis()->SetRangeUser(0.01,20);

    hFrameP[1]->Draw("axis");
    gPad->Update();


    Int_t Npp = hVvsTime[isim]->GetXaxis()->GetNbins();
    Double_t *zval = new Double_t[Npp];
    Double_t *nval = new Double_t[Npp];
    for(Int_t j=0; j<Npp; j++) {
      zval[j] = hVvsTime[isim]->GetXaxis()->GetBinCenter(j+1);
      Float_t density =  1 + (ntop0-1) * TMath::Exp(-(zval[j]*zval[j])/(2*sigmal0*sigmal0));
      nval[j] = density;
    }

    xaxis2.DrawLine(gPad->GetUxmin(), nres0,
		    gPad->GetUxmax(), nres0);
    yaxis2.DrawLine(2*sigmal0,gPad->GetUymin(),
		    2*sigmal0,gPad->GetUymax());
    
    TGraph *gDenvsz = new TGraph(Npp,zval,nval);
    gDenvsz->SetLineWidth(4);
    gDenvsz->SetLineColor(kGray+3);
    gDenvsz->Draw("L");
    
    //----------------------------------  
    CVR->cd(0);
    pad[0]->Draw();
    pad[0]->cd();
    
    //   Float_t Zmax = hVvsTime[isim]->GetYaxis()->GetXmax();
    Float_t Zmax = 1.0;
    Float_t Zmin = hVvsTime[isim]->GetYaxis()->GetXmin();

    hFrameP[0]->GetYaxis()->SetRangeUser(Zmin,Zmax);
    hFrameP[0]->GetXaxis()->SetTitle("k_{p}^{0} z");
    hFrameP[0]->GetXaxis()->CenterTitle();
    hFrameP[0]->GetYaxis()->SetTitle("k_{p}^{0} #zeta");
    hFrameP[0]->GetYaxis()->CenterTitle();
    
    hFrameP[0]->Draw("axis");
    gPad->Update();

    // Change the range of z axis for the fields to be symmetric.
    Float_t Vmax = hVvsTime[isim]->GetMaximum();
    Float_t Vmin = hVvsTime[isim]->GetMinimum();
    if(Vmax > TMath::Abs(Vmin))
      Vmin = -Vmax;
    else
      Vmax = -Vmin;
    hVvsTime[isim]->GetZaxis()->SetRangeUser(Vmin,Vmax);
    
    pal->SetPalette("rbow0");
  
    hVvsTime[isim]->Draw("colz same");

    yaxis2.DrawLine(2*sigmal0,gPad->GetUymin(),
		    2*sigmal0,gPad->GetUymax());
    
    gPad->Update();
    TPaletteAxis *palette3 = (TPaletteAxis*)hVvsTime[isim]->GetListOfFunctions()->FindObject("palette");
    if(palette3) {
      Double_t y1 = gPad->GetBottomMargin();
      Double_t y2 = 1 - gPad->GetTopMargin();
      Double_t x1 = 1 - gPad->GetRightMargin();
    
      Double_t x1b = x1 + 0.01;
      Double_t x2b = x1 + 0.035;
      Double_t y1b = y1 + 0.03;
      Double_t y2b = y2 - 0.03;
      palette3->SetX1NDC(x1b);
      palette3->SetY1NDC(y1b);
      palette3->SetX2NDC(x2b);
      palette3->SetY2NDC(y2b);
      palette3->SetBorderSize(2);
      palette3->SetLineColor(1);

      TPave *pFrame = new TPave(x1b,y1b,x2b,y2b,1,"NDCL");
      pFrame->SetFillStyle(0);
      pFrame->SetLineColor(kBlack);
      pFrame->SetShadowColor(0);
      pFrame->Draw();
    }

    TBox *lFrame4 = new TBox(gPad->GetUxmin(), gPad->GetUymin(),
			     gPad->GetUxmax(), gPad->GetUymax());
    lFrame4->SetFillStyle(0);
    lFrame4->SetLineColor(kBlack);
    lFrame4->SetLineWidth(2);
    lFrame4->Draw();

    gChiminvsz->Draw("C");
    gChimaxvsz->Draw("C");

    gVcross[isim][0]->Draw("C");
    
    gVcross[isim][1]->Draw("C");
  
    gPad->RedrawAxis();

    // Print to a file
    // Output file
    PGlobals::imgconv(CVR,"./DDR/DDR2D-PsiVsZ-withramp","pdf");
    // ---------------------------------------------------------
  }
  

  // Psi extremes vs density (1D)
  // ------------------
  
  sizex = 1024;
  sizey = 380;

  gStyle->SetLineWidth(frameWidth);  

  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadRightMargin(0.10);
  gStyle->SetPadBottomMargin(0.25);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetTitleOffset(titleoffset-0.5,"yz");
  gStyle->SetTitleOffset(titleoffset-0.2,"x");

  // Canvas setup
  TCanvas *CV1D = new TCanvas("CV1D","Evolution of the wakefield potential 1D",sizex,sizey);
  
  CV1D->cd(0);
  
  // gPad->SetFillStyle(4000);
  // gPad->SetFrameFillStyle(4000);

  // gPad->SetLogx(1);
  
  TGraph *gDenvsZ = NULL;
  
  TGraph **gVminvsd = new TGraph*[Nsim];
  TGraph **gVmaxvsd = new TGraph*[Nsim];
  
  for(Int_t i=0;i<Nsim;i++) {
    if(hVvsTime[i]->GetMaximum()>Vmax) Vmax = hVvsTime[i]->GetMaximum();
    if(hVvsTime[i]->GetMinimum()<Vmin) Vmin = hVvsTime[i]->GetMinimum();
  }

  // if(Vmax > TMath::Abs(Vmin))
  //   Vmin = -Vmax;
  // else
  //   Vmax = -Vmin;

  Float_t Vrange = Vmax - Vmin;
  Vmin -= Vrange * 0.3;
  Vmax += Vrange * 0.3;
  Vrange = Vmax-Vmin;

  Vmin = -1.4;
  Vmax =  1.4;

  TH1F *hFrame = new TH1F("hFrame","",100,0.7,16.);
  hFrame->GetXaxis()->SetTitle("n/n_{0}");
  hFrame->GetXaxis()->CenterTitle();

  hFrame->GetYaxis()->SetTitle("#psi");
  hFrame->GetYaxis()->CenterTitle();
  hFrame->GetYaxis()->SetNdivisions(005);
  hFrame->GetYaxis()->SetRangeUser(Vmin,Vmax);

  hFrame->Draw("axis");

  gPad->Update();  
  
  for(Int_t i=0;i<Nsim;i++) {
    
    //  cout << Form("Sim %i: Npoints = %i ",i,gVextr[i][1]->GetN()) << endl;
    
    Int_t Npp = gVextr[i][1]->GetN();
    Double_t *psimin = gVextr[i][1]->GetY();
    Double_t *psimax = gVextr[i][0]->GetY();
    Double_t *zval = gVextr[i][1]->GetX();
    Double_t *nval = new Double_t[Npp];
    
    // Get density array from array of z (assuming a function n(z))
    for(Int_t j=0; j<Npp; j++) {
      Float_t density =  1 + (ntop0-1) * TMath::Exp(-(zval[j]*zval[j])/(2*sigmal0*sigmal0));
      nval[j] = density;

      //  cout << Form(" %i:  z = %.4f   n = %.4f   psimin = %.4f   psimax  = %.4f",j,zval[j],nval[j],psimin[j],psimax[j]) << endl;
      
    }

    gVminvsd[i] = new TGraph(Npp,nval,psimin);
    gVminvsd[i]->SetLineWidth(4);
    gVminvsd[i]->SetLineColor(color[i]);
    gVminvsd[i]->Draw("C");

    gVmaxvsd[i] = new TGraph(Npp,nval,psimax);
    gVmaxvsd[i]->SetLineWidth(4);
    gVmaxvsd[i]->SetLineColor(color2[i]);
    gVmaxvsd[i]->Draw("C");
  }

  xaxis2.DrawLine(gPad->GetUxmin(), 0.0,
		  gPad->GetUxmax(), 0.0);
  
  yaxis2.DrawLine(nres0,gPad->GetUymin(),
		  nres0,gPad->GetUymax());

  gPad->RedrawAxis();
  
  lFrame->SetFillStyle(0);
  lFrame->SetLineColor(PGlobals::frameColor);
  lFrame->SetLineWidth(frameWidth);
  lFrame->DrawBox(gPad->GetUxmin(), gPad->GetUymin(),
		  gPad->GetUxmax(), gPad->GetUymax());
  

  // Print to a file
  // Output file
  PGlobals::imgconv(CV1D,"./DDR/DDR1D-PsiVsDen","pdf");
  // ---------------------------------------------------------


  // Canvas setup
  TCanvas *CB1D = new TCanvas("CB1D","Evolution of the electron velocity",sizex,sizey);
  
  CB1D->cd(0);
  gPad->SetFrameLineWidth(2);  
  gPad->SetTickx(0);
  gPad->SetTicky(0);

  // gPad->SetLogx(1);
  
  TGraph **gBminvsd = new TGraph*[Nsim];
  TGraph **gBmaxvsd = new TGraph*[Nsim];
  TGraph **gBmaxrealvsd = new TGraph*[Nsim];
  TGraph **gBmaxrealvsdsmooth = new TGraph*[Nsim];

  Float_t Bmin = -1.4;
  Float_t Bmax = 1.4;
  
  hFrame->GetYaxis()->SetTitle("#beta");
  hFrame->GetYaxis()->SetRangeUser(Bmin,Bmax);
  hFrame->Draw("axis");

  gPad->Update();
  
  for(Int_t i=0;i<Nsim;i++) {
    
    //  cout << Form("Sim %i: Npoints = %i ",i,gVextr[i][1]->GetN()) << endl;
    
    Int_t Npp = gVextr[i][1]->GetN();
    Double_t *psimin = gVextr[i][1]->GetY();
    Double_t *psimax = gVextr[i][0]->GetY();
    Double_t *zval = gVextr[i][1]->GetX();
    Double_t *nval = new Double_t[Npp];
    Double_t *betamin = new Double_t[Npp];
    Double_t *betamax = new Double_t[Npp];
    Double_t *betamaxreal = gBeta[i]->GetY();
    
    // Get density array from array of z (assuming a function n(z))
    for(Int_t j=0; j<Npp; j++) {
      Float_t density =  1 + (ntop0-1) * TMath::Exp(-(zval[j]*zval[j])/(2*sigmal0*sigmal0));
      nval[j] = density;

      betamin[j] = betapsi(psimin[j]);
      betamax[j] = betapsi(psimax[j]);
    }

    gBminvsd[i] = new TGraph(Npp,nval,betamin);
    gBminvsd[i]->SetLineWidth(4);
    gBminvsd[i]->SetLineColor(color[i]);
    gBminvsd[i]->Draw("C");

    gBmaxvsd[i] = new TGraph(Npp,nval,betamax);
    gBmaxvsd[i]->SetLineWidth(4);
    gBmaxvsd[i]->SetLineColor(color2[i]);
    gBmaxvsd[i]->Draw("C");

    gBmaxrealvsd[i] = new TGraph(Npp,nval,betamaxreal);
    // gBmaxrealvsd[i]->SetLineWidth(1);
    // gBmaxrealvsd[i]->SetLineStyle(1);
    // gBmaxrealvsd[i]->SetLineColor(color[i]);
    // gBmaxrealvsd[i]->Draw("C");

    gBmaxrealvsdsmooth[i] = gsmooth[i]->SmoothSuper(gBmaxrealvsd[i],"",0.0);
    gBmaxrealvsdsmooth[i]->SetLineWidth(3);
    gBmaxrealvsdsmooth[i]->SetLineStyle(7);
    gBmaxrealvsdsmooth[i]->SetLineColor(color[i]);    
    gBmaxrealvsdsmooth[i]->Draw("C");
  }

  yaxis2.DrawLine(nres0,gPad->GetUymin(),
		  nres0,gPad->GetUymax());
  xaxis2.DrawLine(gPad->GetUxmin(), 0.0,
		  gPad->GetUxmax(), 0.0);

  gPad->RedrawAxis();

  lFrame->DrawBox(gPad->GetUxmin(), gPad->GetUymin(),
		  gPad->GetUxmax(), gPad->GetUymax());

  // Print to a file
  // Output file
  PGlobals::imgconv(CB1D,"./DDR/DDR1D-BetaVsDen","pdf");
  // ---------------------------------------------------------

  // Canvas setup
  TCanvas *CC1D = new TCanvas("CC1D","Evolution of the electron velocity",sizex,sizey);
  
  CC1D->cd(0);
  gPad->SetFrameLineWidth(2);  
  gPad->SetTickx(0);
  gPad->SetTicky(0);

  // gPad->SetLogx(1);
  
  TGraph **gCminvsd = new TGraph*[Nsim];
  TGraph **gCmaxvsd = new TGraph*[Nsim];

  //  Float_t Cmin = -TMath::TwoPi();
  Float_t Cmin = -8;
  Float_t Cmax = 1.0;
  
  hFrame->GetYaxis()->SetRangeUser(Cmin,Cmax);
  hFrame->GetYaxis()->SetTitle("#chi");

  hFrame->Draw("axis");

  gPad->Update();

  for(Int_t i=0;i<Nsim;i++) {
    
    //  cout << Form("Sim %i: Npoints = %i ",i,gVextr[i][1]->GetN()) << endl;
    
    Int_t Npp = gVextr[i][1]->GetN();
    Double_t *ezcross1 = gEzcross[i][1]->GetY();
    Double_t *ezcross0 = gEzcross[i][0]->GetY();
    Double_t *zval = gVextr[i][1]->GetX();
    Double_t *nval = new Double_t[Npp];
    Double_t *chimin = new Double_t[Npp];
    Double_t *chimax = new Double_t[Npp];
    
    // Get density array from array of z (assuming a function n(z))
    for(Int_t j=0; j<Npp; j++) {
      Float_t density =  1 + (ntop0-1) * TMath::Exp(-(zval[j]*zval[j])/(2*sigmal0*sigmal0));
      nval[j] = density;
      
      chimin[j] = ezcross1[j] * TMath::Sqrt(nval[j]);
      chimax[j] = ezcross0[j] * TMath::Sqrt(nval[j]);
    }
    
    gCminvsd[i] = new TGraph(Npp,nval,chimin);
    gCminvsd[i]->SetLineWidth(4);
    gCminvsd[i]->SetLineColor(color[i]);
    gCminvsd[i]->Draw("C");

    gCmaxvsd[i] = new TGraph(Npp,nval,chimax);
    gCmaxvsd[i]->SetLineWidth(4);
    gCmaxvsd[i]->SetLineColor(color2[i]);
    gCmaxvsd[i]->Draw("C");
  }
  
  gPad->RedrawAxis();
  
  yaxis2.DrawLine(nres0,gPad->GetUymin(),
		  nres0,gPad->GetUymax());

  xaxis2.DrawLine(gPad->GetUxmin(),0,
		  gPad->GetUxmax(),0);

  xaxis2.DrawLine(gPad->GetUxmin(),-TMath::Pi(),
		  gPad->GetUxmax(),-TMath::Pi());

  xaxis2.DrawLine(gPad->GetUxmin(),-TMath::TwoPi(),
		  gPad->GetUxmax(),-TMath::TwoPi());

  lFrame->DrawBox(gPad->GetUxmin(), gPad->GetUymin(),
		  gPad->GetUxmax(), gPad->GetUymax());

  // Print to a file
  // Output file
  PGlobals::imgconv(CC1D,"./DDR/DDR1D-ChiVsDen","pdf");
  // ---------------------------------------------------------

  // Canvas setup
  TCanvas *CE1D = new TCanvas("CE1D","Evolution of the electric field",sizex,sizey);
  
  CE1D->cd(0);
  gPad->SetFrameLineWidth(2);  
  gPad->SetTickx(0);
  gPad->SetTicky(0);

  // gPad->SetLogx(1);
  TGraph **gEminvsd = new TGraph*[Nsim];
  TGraph **gEmaxvsd = new TGraph*[Nsim];
  TGraph **gRmaxvsd = new TGraph*[Nsim];

  TGraph **gEminvsdsmooth = new TGraph*[Nsim];
    
  Float_t Emin = -5.5;
  Float_t Emax = 1.5;
  
  hFrame->GetYaxis()->SetRangeUser(Emin,Emax);
  hFrame->GetYaxis()->SetTitle("E_{z}/E_{0}");

  hFrame->Draw("axis");

  gPad->Update();

  for(Int_t i=0;i<Nsim;i++) {
    
    Int_t Npp = gVextr[i][1]->GetN();
    Double_t *ezmin = gEzextr[i][1]->GetY();
    Double_t *ezmax = gEzextr[i][0]->GetY();
    Double_t *zval = gEzextr[i][1]->GetX();
    Double_t *nval = new Double_t[Npp];
    Double_t *rmax = new Double_t[Npp];
    
    // Get density array from array of z (assuming a function n(z))
    for(Int_t j=0; j<Npp; j++) {
      Float_t density =  1 + (ntop0-1) * TMath::Exp(-(zval[j]*zval[j])/(2*sigmal0*sigmal0));
      nval[j] = density;

      // if(j==0 ) {
      // 	ezmin[j] = (ezmin[j] + ezmin[j+1])/2.0;
      // } else if(j==Npp-1) {
      // 	ezmin[j] = (ezmin[j] + ezmin[j-1])/2.0;
      // } else {
      // 	ezmin[j] = (ezmin[j] + ezmin[j-1] + ezmin[j+1])/3.0;
      // }

      rmax[j] = fabs(ezmin[j]/ezmax[j]);
      
    }
    
    gEminvsd[i] = new TGraph(Npp,nval,ezmin);
    // gEminvsd[i]->SetLineWidth(2);
    gEminvsd[i]->SetLineColor(color[i]);
    // gEminvsd[i]->Draw("L");

    gEminvsdsmooth[i] = gsmooth[i]->SmoothSuper(gEminvsd[i],"",0.0);
    gEminvsdsmooth[i]->SetLineWidth(4);
    // gEminvsdsmooth[i]->SetLineStyle(7);
    gEminvsdsmooth[i]->SetLineColor(color[i]);    
    gEminvsdsmooth[i]->Draw("L");
    
    gEmaxvsd[i] = new TGraph(Npp,nval,ezmax);
    gEmaxvsd[i]->SetLineWidth(4);
    gEmaxvsd[i]->SetLineColor(color2[i]);
    gEmaxvsd[i]->Draw("L");

    gRmaxvsd[i] = new TGraph(Npp,nval,rmax);
    gRmaxvsd[i]->SetLineWidth(4);
    gRmaxvsd[i]->SetLineColor(color[i]);
    //gRmaxvsd[i]->Draw("L");

  }

  yaxis2.DrawLine(nres0,gPad->GetUymin(),
	      nres0,gPad->GetUymax());
  xaxis2.DrawLine(gPad->GetUxmin(), 0.0,
		  gPad->GetUxmax(), 0.0);

  gPad->RedrawAxis();

  lFrame->DrawBox(gPad->GetUxmin(), gPad->GetUymin(),
		  gPad->GetUxmax(), gPad->GetUymax());
  
  // Print to a file
  // Output file
  PGlobals::imgconv(CE1D,"./DDR/DDR1D-EzVsDen","pdf");
  // ---------------------------------------------------------


  
  // Canvas setup
  sizex = 1024;
  sizey = 1638;
  TCanvas *CA1D = new TCanvas("CA1D","versus density",sizex,sizey);

  const Int_t NPad = 4;
  TPad *pad[NPad];
  TH1F *hFrameP[NPad];

  Float_t lMargin = 0.12;
  Float_t rMargin = 0.08;
  Float_t bMargin = 0.06;
  Float_t tMargin = 0.04;
  Float_t factor = 1.0;    
  PGlobals::CanvasAsymPartition(C,NPad,lMargin,rMargin,bMargin,tMargin,factor,0.02);

    // Define the frames for plotting
  Int_t fonttype = 43;
  Int_t fontsize = 34;
  Int_t tfontsize = 40;
  Float_t txoffset = 2.8;
  Float_t lxoffset = 0.02;
  Float_t tyoffset = 1.8;
  Float_t lyoffset = 0.01;
  Float_t tylength = 0.01;
  Float_t txlength = 0.02;
  for(Int_t i=0;i<NPad;i++) {
    char name[16];
    sprintf(name,"pad_%i",i);
    pad[i] = (TPad*) gROOT->FindObject(name);
    pad[i]->SetFrameLineWidth(2);  
    pad[i]->SetTickx(0);
    pad[i]->SetTicky(0);

    sprintf(name,"hFrameP_%i",i);

    hFrameP[i] = new TH1F(name,"",10,0.8,16.);
    
    Float_t xFactor = pad[0]->GetAbsWNDC()/pad[i]->GetAbsWNDC();
    Float_t yFactor = pad[0]->GetAbsHNDC()/pad[i]->GetAbsHNDC();

    // Format for y axis
    hFrameP[i]->GetYaxis()->SetTitleFont(fonttype);
    hFrameP[i]->GetYaxis()->SetTitleSize(tfontsize);
    hFrameP[i]->GetYaxis()->SetTitleOffset(tyoffset);
    hFrameP[i]->GetYaxis()->SetLabelFont(fonttype);
    hFrameP[i]->GetYaxis()->SetLabelSize(fontsize);
    hFrameP[i]->GetYaxis()->SetLabelOffset(lyoffset);
    hFrameP[i]->GetYaxis()->SetTickLength(xFactor*tylength/yFactor);
    hFrameP[i]->GetYaxis()->CenterTitle();
    hFrameP[i]->GetYaxis()->SetNdivisions(ndivy);

    // Format for x axis
    hFrameP[i]->GetXaxis()->SetTitleFont(fonttype);
    hFrameP[i]->GetXaxis()->SetTitleSize(tfontsize+2);
    hFrameP[i]->GetXaxis()->SetTitleOffset(txoffset);
    hFrameP[i]->GetXaxis()->SetLabelFont(fonttype);
    hFrameP[i]->GetXaxis()->SetLabelSize(fontsize+2);
    hFrameP[i]->GetXaxis()->SetLabelOffset(lxoffset);
    hFrameP[i]->GetXaxis()->CenterTitle();
    hFrameP[i]->GetXaxis()->SetTickLength(yFactor*txlength/xFactor);

    if(i!=0) {
      hFrameP[i]->GetXaxis()->SetLabelOffset(999);
    }
    
  }


  //----------------------------------  
  CA1D->cd(0);
  pad[3]->Draw();
  pad[3]->cd();

  Float_t Rmin = 1.0;
  Float_t Rmax = 9.0;
  
  hFrameP[3]->GetYaxis()->SetRangeUser(Rmin,Rmax);
  hFrameP[3]->GetXaxis()->SetTitle("n/n_{0}");
  hFrameP[3]->GetXaxis()->CenterTitle();
  hFrameP[3]->GetYaxis()->SetTitle("R_{max}");
  hFrameP[3]->GetYaxis()->CenterTitle();

  hFrameP[3]->Draw("axis");
  gPad->Update();
  yaxis2.DrawLine(nres0,gPad->GetUymin(),
	      nres0,gPad->GetUymax());
  xaxis2.DrawLine(gPad->GetUxmin(), 0.0,
		  gPad->GetUxmax(), 0.0);
  
  for(Int_t i=0;i<Nsim;i++) {
    gRmaxvsd[i]->Draw("L");
  }
  
  gPad->RedrawAxis();
  
  
  //----------------------------------  
  CA1D->cd(0);
  pad[2]->Draw();
  pad[2]->cd();

  Emin = -5.99;
  Emax = 1.99;
  
  hFrameP[2]->GetYaxis()->SetRangeUser(Emin,Emax);
  hFrameP[2]->GetXaxis()->SetTitle("n/n_{0}");
  hFrameP[2]->GetXaxis()->CenterTitle();
  hFrameP[2]->GetYaxis()->SetTitle("E_{z}/E_{0}");
  hFrameP[2]->GetYaxis()->CenterTitle();

  hFrameP[2]->Draw("axis");
  gPad->Update();
  yaxis2.DrawLine(nres0,gPad->GetUymin(),
	      nres0,gPad->GetUymax());
  xaxis2.DrawLine(gPad->GetUxmin(), 0.0,
		  gPad->GetUxmax(), 0.0);
  
  for(Int_t i=0;i<Nsim;i++) {
    gEminvsd[i]->Draw("L");
    gEmaxvsd[i]->Draw("L");
  }
  
  gPad->RedrawAxis();


  //----------------------------------
  CA1D->cd(0);
  
  pad[1]->Draw();
  pad[1]->cd();

  Vmin = -1.499;
  Vmax =  1.499;
  
  hFrameP[1]->GetYaxis()->SetRangeUser(Vmin,Vmax);
  hFrameP[1]->GetXaxis()->SetTitle("n/n_{0}");
  hFrameP[1]->GetXaxis()->CenterTitle();
  hFrameP[1]->GetYaxis()->SetTitle("#psi");
  hFrameP[1]->GetYaxis()->CenterTitle();

  hFrameP[1]->Draw("axis");
  gPad->Update();
  yaxis2.DrawLine(nres0,gPad->GetUymin(),
	      nres0,gPad->GetUymax());
  xaxis2.DrawLine(gPad->GetUxmin(), 0.0,
		  gPad->GetUxmax(), 0.0);
  
  for(Int_t i=0;i<Nsim;i++) {
    gVminvsd[i]->Draw("C");
    gVmaxvsd[i]->Draw("C");
  }
  
  gPad->RedrawAxis();

  //----------------------------------
  CA1D->cd(0);
  
  pad[0]->Draw();
  pad[0]->cd();

  Bmin = -1.4;
  Bmax =  1.4;
  
  hFrameP[0]->GetYaxis()->SetRangeUser(Bmin,Bmax);
  hFrameP[0]->GetXaxis()->SetTitle("n/n_{0}");
  hFrameP[0]->GetXaxis()->CenterTitle();
  hFrameP[0]->GetYaxis()->SetTitle("#beta");
  hFrameP[0]->GetYaxis()->CenterTitle();

  hFrameP[0]->Draw("axis");
  gPad->Update();
  yaxis2.DrawLine(nres0,gPad->GetUymin(),
	      nres0,gPad->GetUymax());
  xaxis2.DrawLine(gPad->GetUxmin(), 0.0,
		  gPad->GetUxmax(), 0.0);
  
  for(Int_t i=0;i<Nsim;i++) {
    gBminvsd[i]->Draw("C");
    gBmaxvsd[i]->Draw("C");
    gBmaxrealvsd[i]->Draw("C");
  }
  
  gPad->RedrawAxis();

  
  CA1D->cd(0);
  
  // Print to a file
  // Output file
  PGlobals::imgconv(CA1D,"./DDR/DDR1D-All","pdf");
  // ---------------------------------------------------------



  // Saving objects
  // TString filename = Form("./DDR/DDR.root");
  // TFile *ofile = new TFile(filename,"RECREATE");
  
  // for(Int_t i=0; i<Np; i++) {
  //   for(Int_t j=0; j<Np; j++) {
      
  //     gBetavsz[i][j]->Write(Form("gBetavsz_%i_%i",i,j));
      
  //   }
  // }
  // ofile->Close();
  
  
}
