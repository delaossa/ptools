#include "PGlobals.hh"

/**
 * \brief Save a canvas as .tex, then convert it to .pdf via pdflatex or lualatex
 * \param c pointer to the TCanvas
 * \param output the output name (.pdf)
 * \note the intermediate tex figure and skeleton are also saved (.tex) in /tmp. Adapt to your needs this function.
 */
void SaveAsTikz(TCanvas* const c, const TString output)
{
   if(!c) return;
   if(!output.EndsWith(".pdf",TString::kIgnoreCase)) return;
   
   const TString texName = "/tmp/"+TString(output).ReplaceAll(".pdf",".tex");
   const TString skelName= "/tmp/"+TString(output).ReplaceAll(".pdf","_skeleton.tex");
   
   c->SaveAs(texName);

   //Now we create the skeleton tex file
   TString s = 
   "\\documentclass[crop,tikz]{standalone}\n"
   "\\usepackage{tikz}\n"
   "\\usetikzlibrary{patterns,plotmarks}\n"
   "\\begin{document}\n"
   "\\input{";
   s+=texName+
   "}\n"
   "\\end{document}";

   cout << texName << " " << skelName << endl;
   ofstream g(skelName);//or "/tmp/" + tmpName
   g << s << endl;
   g.close();
   
   //gSystem->ChangeDirectory("/tmp");//To use sth like that if you want to change the subdirectory structure
   gSystem->Exec("lualatex -job-name="+output+" "+skelName);
}


void plotLotov( const TString &options="" ){
  PGlobals::Initialize();

  // Palettes!
  gROOT->Macro("PPalettes.C");
  
  TString opt = options;
  
  TFile *temp=new TFile("cav.root", "RECREATE");
  TNtuple *ntuple = new TNtuple("ntuple","NTUPLE","ib:sz:rm:xi0:ximin:flux0:fluxback:fluxezmax:fluxend:fluxmax:rbox:lbox:ezmax");
  ntuple->ReadFile("cav.dat");
  ntuple->Write();

  const Int_t XBINS = 50;
  Double_t xEdges[XBINS + 1] = {0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7,1.9,
				2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,
				4.1,4.3,4.5,4.7,4.9,5.1,5.3,5.5,5.7,5.9,
				6.1,6.3,6.5,6.7,6.9,7.1,7.3,7.5,7.7,7.9,
				8.1,8.3,8.5,8.7,8.9,9.1,9.3,9.5,9.7,9.9,10.1};
  const Int_t YBINS = 50; 
  Double_t yEdges[YBINS + 1] = {0.01,0.03,0.05,0.07,0.09,0.11,0.13,0.15,0.17,0.19,0.21,
				0.275,0.325,0.375,0.425,0.475,0.525,0.575,0.625,0.675,
				0.75,0.85,0.95,1.05,1.15,1.25,1.35,1.45,1.55,1.65,1.75,1.85,
				2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,
				4.25,4.75,5.25,5.75,6.25,6.75,7.25,8.5,9.5,10.5};
  // Double_t yEdges[YBINS + 1] = {0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95,
  // 				 1.05,1.15,1.25,1.35,1.45,1.55,1.65,1.75,1.85,1.95,
  // 				 2.05,2.15,2.25,2.35,2.45,2.55,2.65,2.75,2.85,2.95,
  // 				 3.05,3.15,3.25,3.35,3.45,3.55,3.65,3.75,3.85,3.95,
  // 				 4.05,4.15,4.25,4.35,4.45,4.55,4.65,4.75,4.85,4.95,5.05};
   				 

  const Int_t Nvar = 3;
  TString varName[Nvar] = {"r_{m}","E_{z,max}","#Psi_{flux,back}"};
  TH2* h2[Nvar]; 
  for(Int_t i=0;i<Nvar;i++) {
    char hName[8];
    sprintf(hName,"h2_%i",i);
    h2[i] = new TH2D(hName, "", XBINS, xEdges, YBINS, yEdges);			 
    //  h2[i]->GetXaxis()->SetTitle("I_{b,max}/8.5kA");
    h2[i]->GetXaxis()->SetTitle("#Lambda_{b}");
    h2[i]->GetXaxis()->CenterTitle();
    h2[i]->GetYaxis()->SetTitle("k_{p}#sigma_{z}");
    h2[i]->GetYaxis()->CenterTitle();
    h2[i]->GetZaxis()->SetTitle(varName[i].Data());
    h2[i]->GetZaxis()->CenterTitle();
  }
   
  Float_t Ib, sz, rm, ezmax, fluxback;
  ntuple->SetBranchAddress("ib",&Ib);
  ntuple->SetBranchAddress("sz",&sz);
  ntuple->SetBranchAddress("rm",&rm);
  ntuple->SetBranchAddress("ezmax",&ezmax);
  ntuple->SetBranchAddress("fluxback",&fluxback);

  Int_t nentries = (Int_t)ntuple->GetEntries();
  for (Int_t i=0; i<nentries; i++) {
    ntuple->GetEntry(i);
     
    Int_t i1=0;
    for(Int_t j=0;j<XBINS+1;j++)
      if(2*Ib<xEdges[j]) { 
	i1 = j;
	break;
      }

    Int_t i2=0;
    for(Int_t j=0;j<YBINS+1;j++)
      if(sz<yEdges[j]) { 
	i2 = j;
	break;
      } 	 

    h2[0]->SetBinContent(i1,i2,rm);
    h2[1]->SetBinContent(i1,i2,ezmax);
    if(fluxback<1000)
      h2[2]->SetBinContent(i1,i2,fluxback);
  }
   
  Float_t Ibmax[XBINS];
  Float_t szmax[XBINS];
  Float_t rmmax[XBINS];
  Float_t szmaxerr[XBINS];
  Float_t Ibmaxerr[XBINS];
  TH1F *hRmax = new TH1F("hRmax","",XBINS,xEdges[0],xEdges[XBINS]);
  for(Int_t i=0;i<XBINS;i++) {
     
    TH1F *h = (TH1F*) h2[0]->ProjectionY("szmax",i+1,i+1);
     
    Ibmax[i] = h2[0]->GetXaxis()->GetBinCenter(i+1);
    rmmax[i] = h->GetMaximum();
    szmax[i] = h->GetBinCenter(h->GetMaximumBin());
    szmaxerr[i] = h->GetBinWidth(h->GetMaximumBin())/2.0;
    Ibmaxerr[i] = h2[0]->GetXaxis()->GetBinWidth(i+1)/2.0;
    hRmax->SetBinContent(i+1,szmax[i]);
  }
   
   
  hRmax->Smooth(1000);
  Float_t szmaxsmooth[XBINS];
  for(Int_t i=0;i<XBINS;i++) {
    szmaxsmooth[i] = hRmax->GetBinContent(i+1);
  }
   
  TGraph *gZmax = new TGraph(XBINS,Ibmax,szmax);
  TGraph *gZmaxSmooth = new TGraph(XBINS,Ibmax,szmaxsmooth);
  TGraphErrors *gZmaxError = new TGraphErrors(XBINS,Ibmax,szmaxsmooth,Ibmaxerr,szmaxerr);

  TGraph *gRmax = new TGraph(XBINS,Ibmax,rmmax);
   
  // FIT
  Float_t xMin = h2[0]->GetXaxis()->GetXmin();
  Float_t xMax = h2[0]->GetXaxis()->GetXmax();

  // model 1
  TF1 *fRsqrt2 = new TF1("fRsqrt2","2.0*sqrt(x)",xMin,xMax);   
  
  char fitName[64];
  sprintf(fitName,"sqrtsqrt");
  TF1 *fitSqrSqr = (TF1*) gROOT->FindObject(fitName);
  if(!fitSqrSqr) 
    fitSqrSqr = new TF1(fitName,"[0]*sqrt(sqrt(x))",xMin,xMax);   
  fitSqrSqr->SetParameter(0,2.0);

  Int_t res = gZmaxError->Fit(fitName,"R0EX0");
   
  // Retrieve parameters of the fit
  Double_t chi2 = fitSqrSqr->GetChisquare();
  Double_t NDF = fitSqrSqr->GetNDF();
  Double_t *par = fitSqrSqr->GetParameters();
  Double_t *parErr = fitSqrSqr->GetParErrors();
   
  cout << Form(" Sqrt fit :  %.2f/%.2f = %.2f",chi2,NDF,chi2/NDF) << endl;

  sprintf(fitName,"pol3");
  TF1 *fitPol3 = (TF1*) gROOT->FindObject(fitName);
  if(!fitPol3) 
    fitPol3 = new TF1(fitName,"[0] + [1]*x + [2]*x*x + [3]*x*x*x",xMin,xMax);   
  fitPol3->SetParameters(0.0,0.0,0.0,0.0);
   
  res = gZmaxError->Fit(fitName,"R0EX0");
   
  // Retrieve parameters of the fit
  chi2 = fitPol3->GetChisquare();
  NDF = fitPol3->GetNDF();
  par = fitPol3->GetParameters();
  parErr = fitPol3->GetParErrors();
   
  cout << endl;
  cout << Form(" Pol3 fit :  chi2/NDF (%.2f/%.2f) = %.2f",chi2,NDF,chi2/NDF) << endl;
  cout << Form(" a0 = %.3f +/- %.3f",par[0],parErr[0]) << endl;
  cout << Form(" a1 = %.3f +/- %.3f",par[1],parErr[1]) << endl;
  cout << Form(" a2 = %.4f +/- %.4f",par[2],parErr[2]) << endl;
  cout << Form(" a3 = %.4f +/- %.4f",par[3],parErr[3]) << endl;
   

  sprintf(fitName,"pol2");
  TF1 *fitPol2 = (TF1*) gROOT->FindObject(fitName);
  if(!fitPol2) 
    fitPol2 = new TF1(fitName,"[0] + [1]*x + [2]*x*x",xMin,xMax);   
  fitPol2->SetParameters(0.0,0.0,0.0,0.0);
   
  res = gZmaxError->Fit(fitName,"R0EX0");
   
  // Retrieve parameters of the fit
  chi2 = fitPol2->GetChisquare();
  NDF = fitPol2->GetNDF();
  par = fitPol2->GetParameters();
  parErr = fitPol2->GetParErrors();
   
  cout << endl;
  cout << Form(" Pol2 fit :  chi2/NDF (%.2f/%.2f) = %.2f",chi2,NDF,chi2/NDF) << endl;
  cout << Form(" a0 = %.3f +/- %.3f",par[0],parErr[0]) << endl;
  cout << Form(" a1 = %.3f +/- %.3f",par[1],parErr[1]) << endl;
  cout << Form(" a2 = %.4f +/- %.4f",par[2],parErr[2]) << endl;


  sprintf(fitName,"sqrt");
  TF1 *fitSqr = (TF1*) gROOT->FindObject(fitName);
  if(!fitSqr) 
    fitSqr = new TF1(fitName,"[0]*sqrt(x)",xMin,xMax);   
  fitSqr->SetParameter(0,2.0);
  res = gRmax->Fit(fitName,"R0EX0");
   
  // Retrieve parameters of the fit
  chi2 = fitSqr->GetChisquare();
  NDF = fitSqr->GetNDF();
  par = fitSqr->GetParameters();
  parErr = fitSqr->GetParErrors();

  cout << endl;
  cout << Form(" Sqrt fit :  chi2/NDF (%.2f/%.2f) = %.2f",chi2,NDF,chi2/NDF) << endl;
  cout << Form(" a0 = %.3f +/- %.3f",par[0],parErr[0]) << endl;
      
  Int_t color = kWhite;
  Int_t markersty = 20;
  Float_t markersiz = 0.6;
  Int_t linewit = 1;

  gZmax->SetLineColor(color);
  gZmax->SetLineWidth(linewit);
  gZmax->SetMarkerColor(color);
  gZmax->SetMarkerStyle(markersty);
  gZmax->SetMarkerSize(markersiz);

  color = kRed;
  gZmaxSmooth->SetLineColor(color);
  gZmaxSmooth->SetLineWidth(linewit);
  gZmaxSmooth->SetMarkerColor(color);
  gZmaxSmooth->SetMarkerStyle(markersty);
  gZmaxSmooth->SetMarkerSize(markersiz);

  color = kRed;
  gZmaxError->SetLineColor(color);
  gZmaxError->SetLineWidth(linewit);
  gZmaxError->SetMarkerColor(color);
  gZmaxError->SetMarkerStyle(markersty);
  gZmaxError->SetMarkerSize(markersiz);

  gRmax->SetLineColor(color);
  gRmax->SetLineWidth(linewit);
  gRmax->SetMarkerColor(color);
  gRmax->SetMarkerStyle(markersty);
  gRmax->SetMarkerSize(markersiz);

  fitSqrSqr->SetLineStyle(2);
  fitSqrSqr->SetLineWidth(1);
  fitSqrSqr->SetLineColor(kWhite);
  
  fitSqr->SetLineStyle(1);
  fitSqr->SetLineWidth(1);
  fitSqr->SetLineColor(kGray+2);
  
  fitPol3->SetLineStyle(1);
  fitPol3->SetLineWidth(1);
  fitPol3->SetLineColor(kWhite);
   
  fitPol2->SetLineStyle(3);
  fitPol2->SetLineWidth(1);
  fitPol2->SetLineColor(kWhite);
   
  fRsqrt2->SetLineStyle(2);
  fRsqrt2->SetLineWidth(1);
  fRsqrt2->SetLineColor(kGray+2);

  // Canvas
  Int_t sizex = 800;
  Int_t sizey = 600;
  char cName[32];
  sprintf(cName,"C");     
  TCanvas *C = (TCanvas*) gROOT->FindObject(cName);
  if(C==NULL) C = new TCanvas("C","Blowout radius",sizex,sizey);
  C->SetFillStyle(4000);
  C->cd();
  gPad->SetTopMargin(0.10);

  C->Clear();

  // Set palette:
  // PPalette * pPalette = (PPalette*) gROOT->FindObject("electron0");
  PPalette * pPalette = (PPalette*) gROOT->FindObject("oli");
  pPalette->cd();
  // const Int_t nRGBs = 9;
  // Double_t stops[nRGBs] = { 0.0000, 0.1250, 0.2500, 0.3750, 0.5000, 0.6250, 0.7500, 0.8750, 1.0000};
  // Double_t red[nRGBs]   = { 0.2081, 0.0591, 0.0779, 0.0231, 0.1801, 0.5300, 0.8185, 0.9955, 0.9763};
  // Double_t green[nRGBs] = { 0.1663, 0.3598, 0.5040, 0.6418, 0.7177, 0.7491, 0.7327, 0.7861, 0.9831};
  // Double_t blue[nRGBs]  = { 0.5292, 0.8683, 0.8384, 0.7913, 0.6424, 0.4661, 0.3498, 0.1967, 0.0538};
  // TColor::CreateGradientColorTable(nRGBs, stops, red, green, blue, 255, 1);
  // gStyle->SetNumberContours(255);

  h2[0]->GetYaxis()->SetRangeUser(0.1,4.99);
  h2[0]->Draw("colz");
  h2[0]->GetZaxis()->SetTickLength(0.02);
  
  // gZmax->Draw("PL");
  // gZmaxSmooth->Draw("C");
  fitSqrSqr->Draw("same");
  fitPol3->Draw("same");
  fitPol2->Draw("same");
  gZmaxError->Draw("PEZ");
  
  gPad->RedrawAxis();
  gPad->Update();

  TPaletteAxis *palAxis = (TPaletteAxis*)h2[0]->GetListOfFunctions()->FindObject("palette");
  if(palAxis) {
    Double_t y1 = gPad->GetBottomMargin();
    Double_t y2 = 1 - gPad->GetTopMargin();
    Double_t x1 = 1 - gPad->GetRightMargin();
    palAxis->SetY2NDC(y2 - 0.04);
    palAxis->SetY1NDC(y1 + 0.04);
    palAxis->SetX1NDC(x1 + 0.01);
    palAxis->SetX2NDC(x1 + 0.04);
  }
   
   
  // Writing the figure
  TString fOutName = Form("./LotovRadius");
  if(opt.Contains("tex")) {
    fOutName += ".tex";
    C->SaveAs(fOutName);
  } else if (opt.Contains("tikz")) {
    fOutName += ".pdf";
    SaveAsTikz(C,fOutName);
  } else
    PGlobals::imgconv(C,fOutName,opt);
  
  // ---------------------------------------------------------
  
  C->Clear();
  
  TH2F *hFrame = (TH2F*) h2[0]->Clone("hFrame");
  hFrame->Reset();
  hFrame->GetYaxis()->SetTitle("k_{p} r_{m}");
  hFrame->GetYaxis()->SetRangeUser(0.001,8.0);
   
  hFrame->Draw("AXIS");
   
  gRmax->Draw("P");
  fRsqrt2->Draw("same");
  fitSqr->Draw("same");
   
  fOutName = Form("./LotovRadius2");
  if(opt.Contains("tex")) {
    fOutName += ".tex";
    C->SaveAs(fOutName);
  } else if (opt.Contains("tikz")) {
    fOutName += ".pdf";
    SaveAsTikz(C,fOutName);
  } else
    PGlobals::imgconv(C,fOutName,opt);

    // ---------------------------------------------------------
   
}
