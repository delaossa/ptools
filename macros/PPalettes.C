{    
  // Global palettes
  
  PPalette *plasmaPalette = (PPalette*) gROOT->FindObject("plasma");
  if(!plasmaPalette) {
    plasmaPalette = new PPalette("plasma");
    plasmaPalette->SetPalette("gray");
  }

  PPalette *beamPalette = (PPalette*) gROOT->FindObject("beam");
  if(!beamPalette) {
    beamPalette = new PPalette("beam");
    beamPalette->SetPalette("elec");
  }

  PPalette *beam2Palette = (PPalette*) gROOT->FindObject("beam2");
  if(!beam2Palette) {
    beam2Palette = new PPalette("beam2");
    beam2Palette->SetPalette("hot");
    //beam2Palette->SetPalette("red");
  }

  PPalette *fieldPalette = (PPalette*) gROOT->FindObject("field");
  if(!fieldPalette) {
    fieldPalette = new PPalette("field");
    fieldPalette->SetPalette("rbow0");
  }

  PPalette *fieldTPalette = (PPalette*) gROOT->FindObject("fieldT");
  if(!fieldTPalette) {
    fieldTPalette = new PPalette("fieldT");
    fieldTPalette->SetPalette("red0");
  }

  PPalette *potPalette = (PPalette*) gROOT->FindObject("pot");
  if(!potPalette) {
    potPalette = new PPalette("pot");
    potPalette->SetPalette("rbow0");
  }

  PPalette *ionPalette = (PPalette*) gROOT->FindObject("ion");
  if(!ionPalette) {
    ionPalette = new PPalette("ion");
    ionPalette->SetPalette("elec0");
  }

  PPalette *defPalette = (PPalette*) gROOT->FindObject("def");
  if(!defPalette) {
    defPalette = new PPalette("def");
    defPalette->SetPalette("electron");
  }
  
  gStyle->SetNumberContours(255);
}
