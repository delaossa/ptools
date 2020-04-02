#ifndef PUTILS_HH
#define PUTILS_HH

#include "TH1.h"
#include "TH2.h"

namespace PUtils
{
  void FindLimits(TH1 *h, Double_t &xmin, Double_t &xmax,
		  Double_t factor = 0.11);

  void FillHisto(TH1 *histo, Double_t *data, Int_t NP, Double_t *weights);

  void FillHisto2D(TH2 *histo, Double_t *datax, Double_t *datay, Int_t NP, Double_t *weights);

  void SetBinsHisto2D(TH2 *histo, Double_t *data);

}

#endif
