#include "PGlobals.hh"


Double_t DenGauss(Double_t z, Double_t sigmal, Double_t ntop) {
  
  Double_t nnorm =  1.0 + (ntop - 1.0) * TMath::Exp(-(z*z)/(2*sigmal*sigmal)) ;
  return nnorm;
}

Double_t DerDenGauss(Double_t z, Double_t sigmal, Double_t ntop) {
  
  Double_t dernnorm =  -(z/(sigmal*sigmal)) * (ntop - 1) * TMath::Exp(-(z*z)/(2*sigmal*sigmal));
  return dernnorm;
}


Double_t betaph(Double_t z, Double_t phase, Double_t sigmal, Double_t ntop) {

  Double_t beta = 1.0 / ( 1.0 + (phase/2) * TMath::Power(DenGauss(z,sigmal,ntop),-3/2.) * DerDenGauss(z,sigmal,ntop)  ) ;
  return beta;
    
}

void PlotDDR2D( const TString &options="" ){
  PGlobals::Initialize();

  gStyle->SetNumberContours(255);

  gStyle->SetTitleFont(43,"xyz");
  gStyle->SetLabelFont(43,"xyz");
  
  gStyle->SetLabelSize(32, "xyz");
  gStyle->SetTitleSize(36, "xyz");
  
  gStyle->SetLabelOffset(0.01, "xyz");
  gStyle->SetTitleOffset(0.9,"yz");
  gStyle->SetTitleOffset(1.1,"x");

  gStyle->SetTickLength(0.01,"yz");
  gStyle->SetTickLength(0.02,"x");
  
  //  gStyle->SetNdivisions(505,"xyz");

  
  Float_t sigmin = 0.0;
  Float_t sigmax = 10.0;
  Float_t denmin = 1.0;
  Float_t denmax = 10.0;

  Float_t phase = -TMath::TwoPi();

  Int_t Np = 100;
  
  TH2F *hBetamin2D = new TH2F("hBetamin2D","",Np,sigmin,sigmax,Np,denmin,denmax);
  TH2F *hBetaminpos2D = new TH2F("hBetaminpos2D","",Np,sigmin,sigmax,Np,denmin,denmax);
  TGraph ***gBetavsz = new TGraph**[Np]; 
  
  Double_t *zarray = new Double_t[Np];
  Double_t *betapharray = new Double_t[Np];
  
  for(Int_t i=0; i<Np; i++) {

    gBetavsz[i] = new TGraph*[Np]; 
    
    for(Int_t j=0; j<Np; j++) {
      
      Float_t sig  = (i + 0.5) * (sigmax-sigmin)/Np + sigmin;
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
  

  // Canvas setup
  Int_t sizex = 1024;
  Int_t sizey = 640;
  TCanvas *C = new TCanvas("C","Minimum phase velocity",sizex,sizey);

  hBetamin2D->GetXaxis()->SetTitle("k_{p}^{0} #sigma_{ramp}");
  hBetamin2D->GetXaxis()->CenterTitle();
  hBetamin2D->GetYaxis()->SetTitle("n/n_{0}");
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

      gr->SetLineWidth(1);
      gr->SetLineStyle(3);
      gr->SetLineColor(kWhite);
      if(i==4) {
	gr->SetLineColor(kRed);
	gr->SetLineWidth(2);
      } 
      
      new(graphs[nGraphs]) TGraph(*gr) ;
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
    pFrame->SetShadowColor(0);
    pFrame->Draw();
  }

  TBox *lFrame = new TBox(gPad->GetUxmin(), gPad->GetUymin(),
			      gPad->GetUxmax(), gPad->GetUymax());
  lFrame->SetFillStyle(0);
  lFrame->SetLineColor(kBlack);
  lFrame->SetLineWidth(1);
  lFrame->Draw();
  
  gPad->RedrawAxis(); 

  PGlobals::imgconv(C,Form("./DDR2D-BetaMin"),"pdf");

  // ------------------------
  
  TCanvas *C2 = new TCanvas("C2","Minimum phase velocity position",sizex,sizey);

  hBetaminpos2D->GetXaxis()->SetTitle("k_{p}^{0} #sigma_{ramp}");
  hBetaminpos2D->GetXaxis()->CenterTitle();
  hBetaminpos2D->GetYaxis()->SetTitle("n/n_{0}");
  hBetaminpos2D->GetYaxis()->CenterTitle();
  hBetaminpos2D->GetZaxis()->SetTitle("z_{min}/#sigma_{ramp}");
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

      gr2->SetLineWidth(1);
      gr2->SetLineStyle(3);
      gr2->SetLineColor(kWhite);
      if(i==4) {
	gr2->SetLineColor(kRed);
	gr2->SetLineWidth(2);
      } 
      
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
  lFrame2->SetLineWidth(1);
  lFrame2->Draw();
  
  gPad->RedrawAxis();

  PGlobals::imgconv(C2,Form("./DDR2D-BetaMinPos"),"pdf");
  
}
