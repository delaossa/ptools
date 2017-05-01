#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <stdlib.h>

#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TCut.h>
#include <TChain.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TProfile.h>
#include <TGraph.h>
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

void DrawErrors(TProfile *hProf,Float_t xmin = -999.0, Float_t xmax = -999.0,
		Int_t color1 = kWhite, Int_t color2 = kWhite) {
  
  Int_t Npoints = hProf->GetNbinsX();
  Int_t Np = 0;
  Double_t *x = new Double_t[Npoints];
  Double_t *y = new Double_t[Npoints];
  Double_t *yup = new Double_t[Npoints];
  Double_t *ydo = new Double_t[Npoints];
  Double_t *ex = new Double_t[Npoints];
  Double_t *ey = new Double_t[Npoints];
  
  if(xmin==-999.0)
    xmin = hProf->GetXaxis()->GetXmin();

  if(xmax==-999.0)
    xmax = hProf->GetXaxis()->GetXmax();

  for(Int_t j=0;j<Npoints;j++) {

    if(hProf->GetBinCenter(j)<xmin || hProf->GetBinCenter(j)>xmax) {
      continue;
    }
    
    x[Np] = hProf->GetBinCenter(j);
    y[Np] = hProf->GetBinContent(j);
    ex[Np] = 0;
    ey[Np] = hProf->GetBinError(j);
    yup[Np] = y[Np]+ey[Np];
    ydo[Np] = y[Np]-ey[Np];
    
    Np++;
  }
  
  // grapherr = new TGraphErrors(Np,x,y,ex,ey);
  // grapherr->SetName("grapherr");
  // grapherr->Draw("3");

  TGraph *grapherravg = new TGraph(Np,x,y); 
  TGraph *grapherrup = new TGraph(Np,x,yup); 
  TGraph *grapherrdo = new TGraph(Np,x,ydo); 
        
  // grapherravg->GetXaxis()->SetRangeUser(xmin,xmax);
  // grapherrup->GetXaxis()->SetRangeUser(xmin,xmax);
  // grapherrdo->GetXaxis()->SetRangeUser(xmin,xmax);

  grapherravg->SetMarkerStyle(20);
  grapherravg->SetMarkerSize(1.0);
  grapherravg->SetMarkerColor(color1);
  grapherravg->SetLineStyle(1);
  grapherravg->SetLineWidth(1);
  grapherravg->SetLineColor(color1);  
  grapherravg->Draw("PL");
  grapherrup->SetMarkerStyle(20);
  grapherrup->SetMarkerSize(0.2);
  grapherrup->SetMarkerColor(color2);
  grapherrup->SetLineStyle(2);
  grapherrup->SetLineWidth(1);
  grapherrup->SetLineColor(color2);
  grapherrup->Draw("L");
  grapherrdo->SetMarkerStyle(20);
  grapherrdo->SetMarkerSize(0.2);
  grapherrdo->SetMarkerColor(color2);
  grapherrdo->SetLineStyle(2);
  grapherrdo->SetLineWidth(1);
  grapherrdo->SetLineColor(color2);
  grapherrdo->Draw("L");

  return;
}

Double_t DensityGauss(Double_t *x, Double_t *par) {
  Double_t f = 0;

  Double_t n0 = par[0];
  Double_t z0 = par[1];
  Double_t sigma0 = par[2];
  Double_t n1 = par[3];
  Double_t z1 = par[4];
  Double_t sigma1 = par[5];
  Double_t n2 = par[6];
  Double_t z2 = par[7];
  Double_t sigma2 = par[8];
  
  if( x[0] < z0) f = n0 * TMath::Gaus(x[0],z0,sigma0);
  else if( x[0] < z1) f = n0;
  else if( x[0] >= z1 && x[0]<z2) f = (n0-n1) * TMath::Gaus(x[0],z1,sigma1) + n1;
  else if( x[0] >= z2) f = (n1-n2) *  TMath::Gaus(x[0],z2,sigma2) + n2;
  
  return f;
}

bool Comp(TTree *a, TTree *b)
{
  if(b == NULL) 
    return true;

  if(a == NULL)
    return false;

  Double_t xa,xb;
  Double_t ya,yb;
  Double_t ra,rb;

  a->SetBranchAddress("x2",&xa);
  a->SetBranchAddress("x3",&ya);
  a->GetEntry(0);

  b->SetBranchAddress("x2",&xb);
  b->SetBranchAddress("x3",&yb);
  b->GetEntry(0);

  ra = sqrt(xa*xa + ya*ya);
  rb = sqrt(xb*xb + yb*yb);
  
  return ra < rb;
}

void PlotTracks( const TString &sim, Int_t index = 0, const TString &options="") {
  
#ifdef __CINT__  
  gSystem->Load("libptools.so");
#endif

  PGlobals::Initialize();

  TString opt = options;
 
  if(opt.Contains("grid")) {
    gStyle->SetPadGridX(1);
    gStyle->SetPadGridY(1);
  }
  
  gStyle->SetJoinLinePS(2);
  Int_t frameWidth = 2;
  gStyle->SetLineWidth(frameWidth);
    
  // Load PData
  Int_t time = 820;
  PData *pData = PData::Get(sim.Data());
  pData->LoadFileNames(time);
  pData->PrintData();
  
  string *trackfile = pData->GetTrackFileName(index);
  
  // if(trackfile != NULL) {
  //   cout << Form("\n Track file name : %s \n",trackfile->c_str()) << endl;
  // } else {
  //   return 0;
  // }
  
  // Some plasma constants
  Double_t np = pData->GetPlasmaDensity();
  Double_t kp = pData->GetPlasmaK();
  Double_t skd = 1.;
  if(kp!=0.0) skd = 1/kp;
  Double_t E0 = pData->GetPlasmaE0();

  // Time in OU
  Float_t Time = pData->GetRealTime();
  // z start of the plasma in normalized units.
  Float_t zStartPlasma = pData->GetPlasmaStart()*kp;
  // z start of the beam in normalized units.
  Float_t zStartBeam = pData->GetBeamStart()*kp;
  
  opt += "comov";

  // Centering shifts
  Double_t shiftzeta = 0.0;
  if(opt.Contains("center")) {
    shiftzeta -= zStartBeam;
  }
  
  Double_t  zmin = 0; 
  Double_t  zmax = 1000; 
  if(opt.Contains("units")) {
    zmin *= skd;
    zmax *= skd;
  }
  
  Double_t zetamin = 2 + shiftzeta; 
  Double_t zetamax = 10.5 + shiftzeta; 

  Double_t zetamin2 = 2.5 + shiftzeta; 
  Double_t zetamax2 = 5.0 + shiftzeta; 

  Double_t szetamin = 2.8 + shiftzeta; 
  Double_t szetamax = 4.5 + shiftzeta; 
  
  if(sim.Contains("DDR.20")) {
    szetamin = 3.0 + shiftzeta; 
    szetamax = 4.3 + shiftzeta; 
  }
  
  Double_t rb0min = 0.001;
  Double_t rb0max = 1.499;

  // Define density profile
  Double_t n0 = 10.0;
  Double_t z0 = 40.0;
  Double_t sigma0 = 7.5;
  Double_t n1 = 1.0;
  Double_t z1 = z0;
  Double_t sigma1 = 2.5;
  //Double_t sigma1 = 5.0;

  if(sim.Contains("DDR.10"))
    sigma1 = 2.5;
  else if(sim.Contains("DDR.20"))
    sigma1 = 5.0;
    
  Double_t n2 = 0;
  Double_t z2 = 10000;
  Double_t sigma2 = 0.0;
  
  if(opt.Contains("units")) {
    n0 *= np;
    z0 *= skd;
    sigma0 *= skd;
    n1 *= np;
    z1 *= skd;
    sigma1 *= skd;
    n2 *= np;
    z2 *= skd;
    sigma2 *= skd;
  }

  TF1 *fDenProf = new TF1("fDenProf",DensityGauss,zmin,zmax,9);
  fDenProf->SetParameters(n0,z0,sigma0,n1,z1,sigma1,n2,z2,sigma2);
  fDenProf->SetNpx(10000);
    
  PPalette *pal = new PPalette("pal");
  pal->SetPalette(kBird);
  // pal->SetPalette("oli");
  // pal->SetPalette(kColorPrintableOnGrey);
  // pal->Invert();
 
  // Read tracks
  // --------------------------------
  Int_t NT = -1;
  TTree **tracktree = pData->GetTrackTree(trackfile->c_str(),NT);

  cout << Form(" Number of tracks = %i : ",NT) << endl;

  // sort tracks
  std::sort(&tracktree[0],&tracktree[NT-1],Comp);
  
  TGraph **grvszeta = new TGraph*[NT];
  TGraph **gvrvszeta = new TGraph*[NT];
  TGraph **gvzvszeta = new TGraph*[NT];
  
  //  TH2D *hRbvszeta = new TH2D("hRbvszeta","",50,zetamin,zetamax,50,rb0min,rb0max);
  Int_t NBinsY = 100;

  // BOX limits
  Double_t X1MIN = pData->GetXMin(0);
  Double_t X1MAX = pData->GetXMax(0);

  Double_t dx1 = pData->GetDX(0);
  Double_t dzf = 2.0;
  Double_t ddx1 = dzf * dx1;
  zetamin2 = floor((zetamin2-X1MIN)/ddx1) * ddx1 + X1MIN;  
  zetamax2 = floor((zetamax2-X1MIN)/ddx1) * ddx1 + X1MIN;  
  
  Int_t NBinsX = ceil ((zetamax2 - zetamin2)/(ddx1));
  
  TH2D *hRbvszeta = new TH2D("hRbvszeta","",NBinsX,zetamin2,zetamax2,NBinsY,0.301,1.099); 
  TH2D *hn0vszeta = new TH2D("hn0vszeta","",NBinsX,zetamin2,zetamax2,NBinsY,1.0,5.99); 
  
  Double_t rmin =  999;
  Double_t rmax = -999;
  Double_t zetamint =  999;
  Double_t zetamaxt = -999;
  Double_t z;
  Double_t t;
  Double_t x;
  Double_t y;
  Double_t pz,px,py;
  Double_t vz,vx,vy;
  Double_t gamma;
  Double_t r;
  Double_t q;
  Int_t NSel = 0;
  for(Int_t i=0;i<NT;i++) {
        
    tracktree[i]->SetBranchAddress("x1",&z);
    tracktree[i]->SetBranchAddress("x2",&x);
    tracktree[i]->SetBranchAddress("x3",&y);
    tracktree[i]->SetBranchAddress("p1",&pz);
    tracktree[i]->SetBranchAddress("p2",&px);
    tracktree[i]->SetBranchAddress("p3",&py);
    tracktree[i]->SetBranchAddress("t",&t);
    tracktree[i]->SetBranchAddress("q",&q);

    tracktree[i]->GetEntry(0);

    Float_t den = fDenProf->Eval(z);
    Float_t kpr = TMath::Sqrt(den) * TMath::Sqrt(x*x + y*y);
    Float_t r = TMath::Sqrt(x*x + y*y);

    Float_t den0 = den;
    Double_t rb0 = r;

    tracktree[i]->GetEntry(tracktree[i]->GetEntries()-1);
    Double_t zetaf = z-t + shiftzeta;
    
    // cout << Form(" n = %.2f rb0 = %.2f  zetaf = %.2f  q = %.2f ",den,rb0,zetaf,q) << endl;
    
    //hRbvszeta->Fill(zetaf,rb0,TMath::Abs(q));
    hRbvszeta->Fill(zetaf,kpr,TMath::Abs(q));
    hn0vszeta->Fill(zetaf,den0,TMath::Abs(q));
      
    grvszeta[i] = NULL;
    gvrvszeta[i] = NULL;
    gvzvszeta[i] = NULL;
    
    // Cut!!
    // if(z>45.0 && z<45.05)
    if(den>2.0 && den<2.03) {
    // if(den>2.24 && den<2.28)
    // if(den>3.0 && den<3.03)
      grvszeta[i] = new TGraph(tracktree[i]->GetEntries());
      gvrvszeta[i] = new TGraph(tracktree[i]->GetEntries());
      gvzvszeta[i] = new TGraph(tracktree[i]->GetEntries());

    } else
      continue;

    if(rb0>rmax) rmax = rb0;
    if(rb0<rmin) rmin = rb0;

    if(zetaf>zetamaxt) zetamaxt = zetaf;
    if(zetaf<zetamint) zetamint = zetaf;
       
    NSel++;
    
    for(Int_t j=0;j<tracktree[i]->GetEntries();j++) {
      
      tracktree[i]->GetEntry(j);

      den = fDenProf->Eval(z);
      //r = TMath::Sqrt(den) * TMath::Sqrt(x*x + y*y);
      r = TMath::Sqrt(x*x + y*y);

      Double_t p2 = pz*pz+py*py+px*px;
      gamma = TMath::Sqrt(p2+1);

      vz = pz/gamma;
      vx = px/gamma;
      vy = py/gamma;

      Double_t zeta = z-t + shiftzeta;
      grvszeta[i]->SetPoint(j,zeta,r);

      gvzvszeta[i]->SetPoint(j,zeta,vz);

      Double_t vr = TMath::Sqrt(vx*vx+vy*vy);

      // Double_t cphi = x/r;
      // Double_t sphi = y/r;
      // Double_t vr;
      // if(vx>vy)
      // 	vr = vx/cphi;
      // else
      // 	vr = vy/sphi;
      
      gvrvszeta[i]->SetPoint(j,zeta,vr);
      
      
    }
    
    grvszeta[i]->SetName(Form("track_rvsz_%i",i));
    gvrvszeta[i]->SetName(Form("track_vrvsz_%i",i));
    gvzvszeta[i]->SetName(Form("track_vzvsz_%i",i));
    
  }
  
  cout << Form(" Number of tracks selected = %i : ",NSel) << endl;;

  // Plot selected tracks
  // --------------------------

  // Set color code
  for(Int_t i=0;i<NT;i++) {
    if(!grvszeta[i]) continue;
    
    tracktree[i]->SetBranchAddress("x1",&z);
    tracktree[i]->SetBranchAddress("x2",&x);
    tracktree[i]->SetBranchAddress("x3",&y);
    tracktree[i]->SetBranchAddress("t",&t);
    
    tracktree[i]->GetEntry(0);
    
    Float_t den = fDenProf->Eval(z);
    Float_t kpr = TMath::Sqrt(den) * TMath::Sqrt(x*x + y*y);
    Float_t r = TMath::Sqrt(x*x + y*y);

    Int_t index =  (r - rmin) * ((pal->GetNColors()-1) / (rmax-rmin));
    Int_t color = pal->GetColorIndex(index);
    
    grvszeta[i]->SetLineColor(color);
    // grvszeta[i]->SetLineColorAlpha(color,0.8);
    grvszeta[i]->SetLineWidth(2);

    gvrvszeta[i]->SetLineColor(color);
    // gvrvszeta[i]->SetLineColorAlpha(color,0.8);
    gvrvszeta[i]->SetLineWidth(2);

    gvzvszeta[i]->SetLineColor(color);
    // gvzvszeta[i]->SetLineColorAlpha(color,0.8);
    gvzvszeta[i]->SetLineWidth(2);

    grvszeta[i]->Draw("C same");
  }

  Int_t sizex = 800;
  Int_t sizey = 900;

  char cName[32];
  sprintf(cName,"C");     
  TCanvas *C = (TCanvas*) gROOT->FindObject(cName);
  if(C==NULL) C = new TCanvas("C","",sizex,sizey);
  C->SetFillStyle(4000);
  C->cd();
  C->Clear();

  // Setup Pad layout: 
  Int_t NPad = 3;
  TPad **pad = new TPad*[NPad];
  TH1F **hFrame = new TH1F*[NPad];
  
  Float_t lMargin = 0.16;
  Float_t rMargin = 0.04;
  Float_t bMargin = 0.12;
  Float_t tMargin = 0.04;
  Float_t mMargin = 0.02;

  // Define the frames for plotting
  Int_t fonttype = 43;
  Int_t fontsize = 32;
  Int_t tfontsize = 38;
  Float_t txoffset = 2.8;
  Float_t lxoffset = 0.01;
  Float_t tyoffset = 1.6;
  Float_t lyoffset = 0.01;
  Float_t tylength = 0.01;
  Float_t txlength = 0.02;
  Float_t pfactor = 1.4;
  
  // Int_t NdivX = 505;
  // Int_t NdivY = 505;
  Int_t NdivX = 405;
  Int_t NdivY = 405;

  // PGlobals::CanvasAsymPartition(C,NPad,lMargin,rMargin,bMargin,tMargin);
  PGlobals::CanvasAsymPartition(C,NPad,lMargin,rMargin,bMargin,tMargin,pfactor,mMargin);
  
  for(Int_t i=0;i<NPad;i++) {
    char name[16];
    sprintf(name,"pad_%i",i);
    pad[i] = (TPad*) gROOT->FindObject(name);
    pad[i]->SetFrameLineWidth(2);  
    pad[i]->SetTickx(1);
    pad[i]->SetTicky(1);

    sprintf(name,"hFrame_%i",i);
    hFrame[i] = new TH1F(name,"",10,zetamin,zetamax);
    
    Float_t xFactor = pad[0]->GetAbsWNDC()/pad[i]->GetAbsWNDC();
    Float_t yFactor = pad[0]->GetAbsHNDC()/pad[i]->GetAbsHNDC();
    
    // Format for y axis
    hFrame[i]->GetYaxis()->SetTitleFont(fonttype);
    hFrame[i]->GetYaxis()->SetTitleSize(tfontsize);
    hFrame[i]->GetYaxis()->SetTitleOffset(tyoffset);
    hFrame[i]->GetYaxis()->SetLabelFont(fonttype);
    hFrame[i]->GetYaxis()->SetLabelSize(fontsize-2);
    hFrame[i]->GetYaxis()->SetLabelOffset(lyoffset);
    hFrame[i]->GetYaxis()->SetTickLength(xFactor*tylength/yFactor);
    hFrame[i]->GetYaxis()->CenterTitle();
    hFrame[i]->GetYaxis()->SetNdivisions(NdivY);

    // Format for x axis
    hFrame[i]->GetXaxis()->SetTitleFont(fonttype);
    hFrame[i]->GetXaxis()->SetTitleSize(tfontsize+2);
    hFrame[i]->GetXaxis()->SetTitleOffset(txoffset);
    hFrame[i]->GetXaxis()->SetLabelFont(fonttype);
    hFrame[i]->GetXaxis()->SetLabelSize(fontsize+2);
    hFrame[i]->GetXaxis()->SetLabelOffset(lxoffset);
    hFrame[i]->GetXaxis()->SetTickLength(yFactor*txlength/xFactor);      
    hFrame[i]->GetXaxis()->CenterTitle();
    hFrame[i]->GetXaxis()->SetNdivisions(NdivX);

    if(i>0) {
      hFrame[i]->GetXaxis()->SetLabelOffset(99);
      hFrame[i]->GetXaxis()->SetTitleOffset(99);
    }

  }

  UInt_t ipad = NPad-1;

  C->cd(0);
  pad[ipad]->Draw();
  pad[ipad]->cd();

  hFrame[ipad]->GetYaxis()->SetRangeUser(rb0min,rb0max);
  
  hFrame[ipad]->GetXaxis()->SetTitle("k_{p}^{0} #zeta");
  hFrame[ipad]->GetYaxis()->SetTitle("k_{p}^{0} r");

  
  // Plot frame
  hFrame[ipad]->Draw("axis");

  // plasma density
  TH2F *hDen2D = pData->GetCharge2DSliceZX(0,-1,2,opt+"avg");	  
  
  Float_t baseden = 1;
  Float_t localden = baseden;
  Float_t basePos = 0.5;
  Float_t localPos = basePos;

  Float_t Max = hDen2D->GetMaximum();
  Float_t Min = 1.01E-1 * baseden;

  if(pData->GetDenMax(0)>0)
    Max = pData->GetDenMax(0);

  if(pData->GetDenMin(0)>0) 
    Min = pData->GetDenMin(0);

  hDen2D->GetZaxis()->SetRangeUser(Min,Max);  

  // Dynamic plasma palette
  PPalette * plasmaPalette = (PPalette*) gROOT->FindObject("plasma");
  if(!plasmaPalette) {
    plasmaPalette = new PPalette("plasma");
    plasmaPalette->SetPalette("gray");
  }

  if( hDen2D ) {
    const Int_t NRGBs = 3;
    const Int_t NCont = 64;
    if(Max!=Min) {
      if(opt.Contains("logz")) {
	Float_t a = 1.0/(TMath::Log10(Max)-TMath::Log10(Min));
	Float_t b = TMath::Log10(Min);
	basePos = a*(TMath::Log10(baseden) - b);
	localPos = a*(TMath::Log10(localden) - b);
      } else {
	basePos = (1.0/(Max-Min))*(baseden - Min);
	localPos = (1.0/(Max-Min))*(localden - Min);
      }
      // cout << Form(" Local density = %.2f    Local position  = %f ",localden,localPos) << endl;
    }
    
    Double_t Stops[NRGBs] = { 0.00, basePos, 1.00 };
    Double_t Red[NRGBs]   = { 0.99, 0.90, 0.00 };
    Double_t Green[NRGBs] = { 0.99, 0.90, 0.00 };
    Double_t Blue[NRGBs]  = { 0.99, 0.90, 0.00 };
    
    plasmaPalette->CreateGradientColorTable(NRGBs, Stops, Red, Green, Blue, NCont, 1.0);
  }

  TExec *exPlasma = new TExec("exPlasma","plasmaPalette->cd();");

  if(opt.Contains("logz")) {
    gPad->SetLogz(1);
  } else {
    gPad->SetLogz(0);
  }

  exPlasma->Draw();
  hDen2D->Draw("col same");

  gPad->RedrawAxis("g");

  for(Int_t i=0;i<NT;i++) {
    if(!grvszeta[i]) continue;
    
    grvszeta[i]->Draw("C same");
  }

  gPad->Update();

  TLine zetalinemin(zetamint,gPad->GetUymin(),zetamint,gPad->GetUymax());
  zetalinemin.SetLineColor(kGray+3);
  zetalinemin.SetLineStyle(3);
  zetalinemin.SetLineWidth(2);
  zetalinemin.Draw();

  TLine zetalinemax(zetamaxt,gPad->GetUymin(),zetamaxt,gPad->GetUymax());
  zetalinemax.SetLineColor(kGray+3);
  zetalinemax.SetLineStyle(3);
  zetalinemax.SetLineWidth(2);
  zetalinemax.Draw();
  
  TBox *lFrame = new TBox(gPad->GetUxmin(), gPad->GetUymin(),
			  gPad->GetUxmax(), gPad->GetUymax());
  
  lFrame->SetFillStyle(0);
  lFrame->SetLineColor(kBlack);
  lFrame->SetLineWidth(frameWidth);
  lFrame->Draw();

  hFrame[ipad]->Draw("axis same");

  C->cd(0);

  ipad--;
  
  C->cd(0);
  pad[ipad]->Draw();
  pad[ipad]->cd();
  
  hFrame[ipad]->GetYaxis()->SetRangeUser(0.001,0.599);
  //hFrame[ipad]->GetYaxis()->SetRangeUser(-0.999,0.999);
  
  hFrame[ipad]->GetXaxis()->SetTitle("k_{p}^{0} #zeta");
  hFrame[ipad]->GetYaxis()->SetTitle("#beta_{r}");
  
  // Plot frame
  hFrame[ipad]->Draw("axis");

  gPad->RedrawAxis("g");
  
  for(Int_t i=0;i<NT;i++) {
    if(!gvrvszeta[i]) continue;
    
    gvrvszeta[i]->Draw("C same");
  }
  
  gPad->Update();
  
  zetalinemin.DrawLine(zetamint,gPad->GetUymin(),zetamint,gPad->GetUymax());
  zetalinemax.DrawLine(zetamaxt,gPad->GetUymin(),zetamaxt,gPad->GetUymax());
  
  lFrame->DrawBox(gPad->GetUxmin(), gPad->GetUymin(),
		  gPad->GetUxmax(), gPad->GetUymax());
  
  hFrame[ipad]->Draw("axis same");


  C->cd(0);

  ipad--;
  
  C->cd(0);
  pad[ipad]->Draw();
  pad[ipad]->cd();
  
  hFrame[ipad]->GetYaxis()->SetRangeUser(-0.599,1.0);
  
  hFrame[ipad]->GetXaxis()->SetTitle("k_{p}^{0} #zeta");
  hFrame[ipad]->GetYaxis()->SetTitle("#beta_{z}");
  
  // Plot frame
  hFrame[ipad]->Draw("axis");

  gPad->Update();

  TLine zeroline(gPad->GetUxmin(),0.0,gPad->GetUxmax(),0.0);
  zeroline.SetLineColor(kBlack);
  zeroline.SetLineStyle(2);
  zeroline.SetLineWidth(frameWidth);
  zeroline.Draw();

  gPad->RedrawAxis("g");

  for(Int_t i=0;i<NT;i++) {
    if(!gvzvszeta[i]) continue;
    
    gvzvszeta[i]->Draw("C same");
  }
  
  
  zetalinemin.DrawLine(zetamint,gPad->GetUymin(),zetamint,gPad->GetUymax());
  zetalinemax.DrawLine(zetamaxt,gPad->GetUymin(),zetamaxt,gPad->GetUymax());
  
  lFrame->DrawBox(gPad->GetUxmin(), gPad->GetUymin(),
		  gPad->GetUxmax(), gPad->GetUymax());
  
  hFrame[ipad]->Draw("axis same");

  C->cd(0);

  // C->cd();
  
  // Print to file -------------------------------------------

  TString fOutName = Form("./%s/Plots/Tracks/%s/Tracks-rvszeta-%s-%s",sim.Data(),pData->GetRawSpeciesName(index).c_str(),pData->GetRawSpeciesName(index).c_str(),sim.Data());
  PGlobals::imgconv(C,fOutName,opt);
  
  // ---------------------------------------------------------

  // New plot over the same canvas
  ipad = NPad-1;
  
  PPalette *pal2 = new PPalette("pal2");
  pal2->SetPalette(kColorPrintableOnGrey);
  //pal2->SetPalette(kGreyScale);
  //pal2->Invert();
  //pal2->SetPalette(kBird);
  //pal2->SetPalette("oli");
  //pal2->Invert();
  pal2->cd();

  TExec *exPal2 = new TExec("exPal2","pal2->cd();");
  // -----  

  C->cd(0);
  pad[ipad]->Draw();
  pad[ipad]->cd();
  gPad->SetLogz(0);

  hFrame[ipad]->SetBins(10,hRbvszeta->GetXaxis()->GetXmin(),hRbvszeta->GetXaxis()->GetXmax());
  hFrame[ipad]->GetYaxis()->SetRangeUser(hRbvszeta->GetYaxis()->GetXmin(),hRbvszeta->GetYaxis()->GetXmax());

  hFrame[ipad]->GetXaxis()->SetTitle("k_{p}^{0} #zeta_{f}");
  hFrame[ipad]->GetYaxis()->SetTitle("k_{p} r_{i}");
  
  hFrame[ipad]->Draw("axis");

  gPad->RedrawAxis("g");
  
  TProfile *hRbvszetaProf = hRbvszeta->ProfileX("hRbvszetaProf",1,-1,"s");

  exPal2->Draw();
  hRbvszeta->Draw("col same");

  Int_t pcolor = kGray+3;
  Int_t mstyle = 21;
  
  hRbvszetaProf->SetLineColor(pcolor);
  hRbvszetaProf->SetMarkerColor(pcolor);
  hRbvszetaProf->SetMarkerStyle(mstyle);
  hRbvszetaProf->SetMarkerSize(0.8);
  //hRbvszetaProf->Draw("same");

  Int_t linecolor = kGray+3;
  DrawErrors(hRbvszetaProf,szetamin,szetamax,linecolor,linecolor);

  hFrame[ipad]->Draw("axis same");

  // -----
  C->cd(0);
  ipad--;
  
  pad[ipad]->Draw();
  pad[ipad]->cd();

  hFrame[ipad]->SetBins(10,hn0vszeta->GetXaxis()->GetXmin(),hn0vszeta->GetXaxis()->GetXmax());
  hFrame[ipad]->GetYaxis()->SetRangeUser(hn0vszeta->GetYaxis()->GetXmin(),hn0vszeta->GetYaxis()->GetXmax());

  hFrame[ipad]->GetXaxis()->SetTitle("k_{p}^{0} #zeta_{f}");
  hFrame[ipad]->GetYaxis()->SetTitle("n_{i}/n_{0}");

  hFrame[ipad]->GetYaxis()->LabelsOption("h");
  hFrame[ipad]->GetYaxis()->SetDrawOption("L");
  hFrame[ipad]->GetYaxis()->ChangeLabel(1,-1,0);//-1,-1,3,-1,"6th label");
  hFrame[ipad]->Draw("axis");

  gPad->RedrawAxis("g");

  TProfile *hn0vszetaProf = hn0vszeta->ProfileX("hn0vszetaProf",1,-1,"s");

  exPal2->Draw();
  hn0vszeta->Draw("col same");

  // pcolor = kGray+2;
  hn0vszetaProf->SetLineColor(pcolor);
  hn0vszetaProf->SetMarkerColor(pcolor);
  hn0vszetaProf->SetMarkerStyle(mstyle);
  //hn0vszetaProf->Draw("same");

  DrawErrors(hn0vszetaProf,szetamin,szetamax,linecolor,linecolor);

  hFrame[ipad]->Draw("axis same");

  
  // -----
  C->cd(0);
  ipad--;
  
  pad[ipad]->Draw();
  pad[ipad]->cd();

  hFrame[ipad]->SetBins(10,hn0vszeta->GetXaxis()->GetXmin(),hn0vszeta->GetXaxis()->GetXmax());
  hFrame[ipad]->GetYaxis()->SetRangeUser(0.001,0.0399);
  
  hFrame[ipad]->GetXaxis()->SetTitle("k_{p}^{0} #zeta_{f}");
  hFrame[ipad]->GetYaxis()->SetTitle("k_{p}^{0} #varepsilon_{n}");

  hFrame[ipad]->Draw("axis");

  //
  TString filename = Form("./%s/Plots/Bunch/%s/Bunch-%s-%s_%i.root",sim.Data(),pData->GetRawSpeciesName(index).c_str(),pData->GetRawSpeciesName(index).c_str(),sim.Data(),time);
  TFile *sFile = new TFile(filename,"READ");

  cout << Form("\n Reading bunch file : %s \n",filename.Data()) << endl;
  TGraph *gEmitx = NULL;
  TGraph *gEmity = NULL;
  if(sFile) {
    gEmitx = (TGraph*) sFile->Get("gEmitx");
    gEmity = (TGraph*) sFile->Get("gEmity");
  } else {
    cout << Form("  [no file]") << endl;
  }

  gPad->RedrawAxis("g");

  gEmity->Draw("LP");
  gEmitx->Draw("LP");

  hFrame[ipad]->Draw("axis same");


  fOutName = Form("./%s/Plots/Tracks/%s/Tracks-hn0vszeta-%s-%s",sim.Data(),pData->GetRawSpeciesName(index).c_str(),pData->GetRawSpeciesName(index).c_str(),sim.Data());
  PGlobals::imgconv(C,fOutName,opt);
  
}
