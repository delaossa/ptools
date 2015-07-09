{
  TString curIncludePath(gSystem->GetIncludePath());
  gSystem->SetIncludePath( " -I/data/netapp/fla/plasma/software/pPlasmaAlt/include " + curIncludePath );

  cout << gSystem->GetIncludePath() << endl;
}
