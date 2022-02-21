#include "PUtils.hh"

void PUtils::FindLimits(TH1 *h,Double_t &xmin, Double_t &xmax, Double_t factor) {
  Double_t maxValue = h->GetBinContent(h->GetMaximumBin());
  for(Int_t i=1;i<=h->GetNbinsX();i++) {
    Double_t binValue = h->GetBinContent(i);
    if(binValue>maxValue*factor) {
      xmin = h->GetBinCenter(i);
      break;
    }
  }
  
  for(Int_t i=h->GetNbinsX();i>0;i--) {
    Double_t binValue = h->GetBinContent(i);
    if(binValue>maxValue*factor) {
      xmax = h->GetBinCenter(i);
      break;
    }
  } 
}

void PUtils::FillHisto(TH1 *histo, Double_t *data, Int_t NP, Double_t *weights){
  for(Int_t ip=0; ip<NP; ip++){
    histo->Fill(data[ip], weights[ip]);
  }
}

void PUtils::FillHisto2D(TH2 *histo, Double_t *datax, Double_t *datay, Int_t NP, Double_t *weights){
  for(Int_t ip=0; ip<NP; ip++){
    histo->Fill(datax[ip], datay[ip], weights[ip]);
  }
}

void PUtils::SetBinsHisto2D(TH2 *histo, Double_t *data){
  UInt_t NX = histo->GetNbinsX();
  UInt_t NY = histo->GetNbinsY();
  for(UInt_t i=0; i<NX; i++){
    for(UInt_t j=0; j<NY; j++){
      UInt_t index = (long)j*(long)NX + (long)i;
      histo->SetBinContent(i + 1, j + 1, data[index]);
    }
  }
}

