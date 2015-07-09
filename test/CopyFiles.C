#include "PData.hh"

void CopyFiles(const TString &sim, Int_t time,const TString &opath) {
  
  // Load PData
  PData *pData = PData::Get(sim.Data());
  pData->LoadFileNames(time);
  pData->Copy(opath);


}
