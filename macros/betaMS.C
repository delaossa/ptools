#include "PGlobals.hh"

Double_t theta(Double_t *x, Double_t *par) // mrad
{
  Double_t xx  = x[0];
  Double_t E0  = par[0];
  Double_t X0  = par[1];

  Double_t f = 0.0;
  if(xx>0) f = 1000. * (13.6/E0) * TMath::Sqrt(xx/X0) * (1.0 - 0.038 * TMath::Log(xx/X0)) ;
          
  return f;
}

Double_t Sigmax(Double_t *x, Double_t *par)
{
  Double_t xx  = x[0];
  Double_t E0  = par[0];
  Double_t X0  = par[1];
  Double_t sx0  = par[2];

  Double_t theta0 = theta(x,par) / 1000.; // rad.
  Double_t f = TMath::Sqrt( (sx0*sx0 + (xx*xx/3.0)*theta0*theta0) );
        
  return f;
}

Double_t Sigmaxp(Double_t *x, Double_t *par) // mrad
{
  Double_t xx  = x[0];
  Double_t E0  = par[0];
  Double_t X0  = par[1];
  Double_t sd0  = par[3];

  Double_t theta0 = theta(x,par);
  Double_t f = TMath::Sqrt( (sd0*sd0 + theta0*theta0) ); // mrad
        
  return f;
}

Double_t emittance(Double_t *x, Double_t *par)
{
  Double_t E0  = par[0];
  Double_t gamma = E0/0.511; 
  Double_t f = gamma * Sigmax(x,par) * (Sigmaxp(x,par) / 1000.); // um
        
  return f;
}

Double_t beta(Double_t *x, Double_t *par)
{
  Double_t f = Sigmax(x,par) / (Sigmaxp(x,par)); // mm
        
  return f;
}

Double_t finalemittance(Double_t *x, Double_t *par)
{
  Double_t iniemit = emittance(x,par);
  Double_t betamatch = par[4]; 
  Double_t inibetap = beta(x,par)/betamatch;

  Double_t f = iniemit * 0.5 * ( 1.0/inibetap + inibetap);
        
  return f;
}


void betaMS(Double_t z = 1000, const TString &opt="pdf") {

  gROOT->Reset();
  PGlobals::Initialize();
  
  // gStyle->SetPadTopMargin(0.10);
  // gStyle->SetPadRightMargin(0.05);
  // gStyle->SetPadBottomMargin(0.15);
  // gStyle->SetPadLeftMargin(0.15);

  Double_t E0  = 1000;  // MeV
  Double_t gamma = E0/0.511; 
  Double_t X0  = 35.0E4; // um (Berylium)
  Double_t sx0 = 10.0;   // um
  Double_t emitxn  = 1.00; // um
  Double_t emitx   = emitxn/gamma; // um
  Double_t sd0 = 1000. * (emitx/sx0); // mrad
  Double_t kp  = 0.0595; // um^-1  for n0 = 10^17 cm^-3
  Double_t kbeta  = kp/TMath::Sqrt(2 * gamma);
  Double_t betamatch = (1/kbeta) / 1000.; // mm

  TF1 *ftheta = new TF1("ftheta",theta,-0.01*z,z,4);
  ftheta->SetParameters(E0,X0,sx0,sd0);
  ftheta->SetParNames("Energy","X0","sx0","sd0");
  ftheta->SetTitle("");
  ftheta->GetXaxis()->SetTitle("#Deltaz [#mum]");
  ftheta->GetYaxis()->SetTitle("#Theta_{0} [mrad]");

  TF1 *fsigmax = new TF1("fsigmax",Sigmax,-0.01*z,z,4);
  fsigmax->SetParameters(E0,X0,sx0,sd0);
  fsigmax->SetParNames("Energy","X0","sx0","sd0");
  fsigmax->SetTitle("");
  fsigmax->GetXaxis()->SetTitle("#Deltaz [#mum]");
  fsigmax->GetYaxis()->SetTitle("#sigma_{x} [#mum]");

  TF1 *fsigmaxp = new TF1("fsigmaxp",Sigmaxp,-0.01*z,z,4);
  fsigmaxp->SetParameters(E0,X0,sx0,sd0);
  fsigmaxp->SetParNames("Energy","X0","sx0","sd0");
  fsigmaxp->SetTitle("");
  fsigmaxp->GetXaxis()->SetTitle("#Deltaz [#mum]");
  fsigmaxp->GetYaxis()->SetTitle("#sigma_{x'} [mrad]");

  TF1 *femit = new TF1("femit",emittance,-0.01*z,z,4);
  femit->SetParameters(E0,X0,sx0,sd0);
  femit->SetParNames("Energy","X0","sx0","sd0");
  femit->SetTitle("");
  femit->GetXaxis()->SetTitle("#Deltaz [#mum]");
  femit->GetYaxis()->SetTitle("#epsilon_{x,n} [#mum]");

  TF1 *fbeta = new TF1("fbeta",beta,-0.01*z,z,4);
  fbeta->SetParameters(E0,X0,sx0,sd0);
  fbeta->SetParNames("Energy","X0","sx0","sd0");
  fbeta->SetTitle("");
  fbeta->GetXaxis()->SetTitle("#Deltaz [#mum]");
  fbeta->GetYaxis()->SetTitle("#beta [mm]");
  
  TF1 *femitfin = new TF1("femitfin",finalemittance,-0.01*z,z,5);
  femitfin->SetParameters(E0,X0,sx0,sd0,betamatch);
  femitfin->SetParNames("Energy","X0","sx0","sd0","betam");
  femitfin->SetTitle("");
  femitfin->GetXaxis()->SetTitle("#Deltaz [#mum]");
  femitfin->GetYaxis()->SetTitle("#epsilon_{x,n,F} [#mum]");

  const UInt_t NPad = 6;
  TCanvas *C = new TCanvas("C","",640,1024);
  // Setup Pad layout:
  TPad **pad = new TPad*[NPad];
  Float_t bMargin = 0.10;
  Float_t tMargin = 0.04;
  Float_t lMargin = 0.20;
  Float_t rMargin = 0.08;
  Float_t vSpacing = 0.017;
  PGlobals::CanvasPartition(C,NPad,lMargin,rMargin,bMargin,tMargin,vSpacing);
 
  // Define the frames for plotting
  Int_t fonttype = 43;
  Int_t fontsize = 22;
  Int_t tfontsize = 26;
  Float_t txoffset = 5.0;
  Float_t lxoffset = 0.02;
  Float_t tyoffset = 3.0;
  Float_t lyoffset = 0.01;
  Float_t tylength = 0.02;
  Float_t txlength = 0.04;

  TF1 *func[NPad];
  func[5] = ftheta;
  func[4] = fsigmax;
  func[3] = fsigmaxp;
  func[2] = femit;
  func[1] = fbeta;
  func[0] = femitfin;

  for(Int_t i=0;i<NPad;i++) {
    char name[16];
    sprintf(name,"pad_%i",i);
    pad[i] = (TPad*) gROOT->FindObject(name);
    pad[i]->SetFrameLineWidth(2);  
    pad[i]->SetTickx(1);
    pad[i]->SetTicky(1);
    pad[i]->SetFillStyle(4000);
    pad[i]->SetFrameFillStyle(4000);

    Float_t xFactor = pad[0]->GetAbsWNDC()/pad[i]->GetAbsWNDC();
    Float_t yFactor = pad[0]->GetAbsHNDC()/pad[i]->GetAbsHNDC();

    // Format for y axis
    func[i]->GetYaxis()->SetTitleFont(fonttype);
    func[i]->GetYaxis()->SetTitleSize(tfontsize);
    func[i]->GetYaxis()->SetTitleOffset(tyoffset);
    func[i]->GetYaxis()->SetLabelFont(fonttype);
    func[i]->GetYaxis()->SetLabelSize(fontsize);
    func[i]->GetYaxis()->SetLabelOffset(lyoffset);

    func[i]->GetYaxis()->SetTickLength(xFactor*tylength/yFactor);
    func[i]->GetYaxis()->CenterTitle();

    // Format for x axis
    func[i]->GetXaxis()->SetTitleFont(fonttype);
    func[i]->GetXaxis()->SetTitleSize(tfontsize+2);
    func[i]->GetXaxis()->SetTitleOffset(txoffset);
    func[i]->GetXaxis()->SetLabelFont(fonttype);
    func[i]->GetXaxis()->SetLabelSize(fontsize);
    func[i]->GetXaxis()->SetLabelOffset(lxoffset);
    
    func[i]->GetXaxis()->SetTickLength(yFactor*txlength/xFactor);      
    func[i]->GetXaxis()->CenterTitle();

    if(i>0) {
      func[i]->GetXaxis()->SetTitleSize(0.0);
      func[i]->GetXaxis()->SetLabelSize(0.0);
    }


    pad[i]->cd();
    //    pad[i]->SetLogy(1);
    func[i]->Draw("L");

  }

  // Calculations:
  cout << Form(" Matched beta = %.4f mm", betamatch) << endl;
  
  Double_t iniemit0  = femit->Eval(0.0);
  Double_t inibeta0 = fbeta->Eval(0.0);
  Double_t inibetap0 = inibeta0/betamatch;
  Double_t finemit0  = iniemit0 * 0.5 * ( 1.0/inibetap0 + inibetap0);
  cout << Form(" Initial beta = %6.2f mm, emit = %6.2f um --> final emittance = %6.3f um %6.3f um ",
	       inibeta0,iniemit0,finemit0,femitfin->Eval(0)) << endl;

  Double_t iniemitZ  = femit->Eval(z);
  Double_t inibetaZ = fbeta->Eval(z);
  Double_t inibetapZ = inibetaZ/betamatch;
  Double_t finemitZ  = iniemitZ * 0.5 * ( 1.0/inibetapZ + inibetapZ);
  cout << Form(" Initial beta = %6.2f mm, emit = %6.2f um --> final emittance = %6.3f um %6.3f um",
	       inibetaZ,iniemitZ,finemitZ,femitfin->Eval(z) ) << endl;

  

  // Writing the figure
  TString fOutName = Form("./betaMS");
  PGlobals::imgconv(C,fOutName,opt);
  // ---------------------------------------------------------
  
  PGlobals::DestroyCanvases();

  

}

