{
  TString IncPath = gSystem->Getenv("PPLASMA");
  IncPath +=  "/inc";
  gROOT->ProcessLine(Form(".include %s",IncPath.Data()));
}
