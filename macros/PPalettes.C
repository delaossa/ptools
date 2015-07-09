{    

  // Palettes !
  PPalette * grayPalette =  (PPalette*) gROOT->FindObject("gray");
  if(!grayPalette) {

    //    cout << "Loading palettes for plasma...  " ;

    const Int_t grayNRGBs = 2;
    const Int_t grayNCont = 64;
    Double_t grayStops[grayNRGBs] = { 0.00, 1.00 };
    Double_t grayRed[grayNRGBs] =   { 0.99, 0.1 };
    Double_t grayGreen[grayNRGBs] = { 0.99, 0.1 };
    Double_t grayBlue[grayNRGBs] =  { 0.99, 0.1 };
    grayPalette = new PPalette("gray");
    grayPalette->CreateGradientColorTable(grayNRGBs, grayStops,grayRed, grayGreen, grayBlue, grayNCont);
  }  

  PPalette * yellowPalette = (PPalette*) gROOT->FindObject("yellow");
  if(!yellowPalette) {
    const Int_t yellowNRGBs = 4;
    const Int_t yellowNCont = 64;
    Double_t yellowStops[yellowNRGBs] = { 0.00, 0.25, 0.75, 1.00 };
    Double_t yellowRed[yellowNRGBs] =   { 1.00, 1.00, 0.60, 0.40 };
    Double_t yellowGreen[yellowNRGBs] = { 0.95, 0.80, 0.30, 0.20 };
    Double_t yellowBlue[yellowNRGBs] =  { 0.80, 0.00, 0.00, 0.00 };
    yellowPalette = new PPalette("yellow");
    yellowPalette->CreateGradientColorTable(yellowNRGBs, yellowStops, 
					    yellowRed, yellowGreen, yellowBlue, yellowNCont);
  }
  
  PPalette * plasmaPalette = (PPalette*) gROOT->FindObject("plasma");
  if(!plasmaPalette) {
    const Int_t plasmaNRGBs = 3;
    const Int_t plasmaNCont = 64;
    Double_t plasmaStops[plasmaNRGBs] = { 0.00,  0.5, 1.00 };
    Double_t plasmaRed[plasmaNRGBs] =   { 0.99, 0.90, 0.00 };
    Double_t plasmaGreen[plasmaNRGBs] = { 0.99, 0.90, 0.00 };
    Double_t plasmaBlue[plasmaNRGBs] =  { 0.99, 0.90, 0.00 };
    plasmaPalette = new PPalette("plasma");
    plasmaPalette->CreateGradientColorTable(plasmaNRGBs, plasmaStops, 
					    plasmaRed, plasmaGreen, plasmaBlue, plasmaNCont);
  }

  PPalette * plasmaHPalette =  (PPalette*) gROOT->FindObject("plasmaH");
  if(!plasmaHPalette) {
    const Int_t plasmaHNRGBs = 3;
    const Int_t plasmaHNCont = 64;
    Double_t plasmaHStops[plasmaHNRGBs] = { 0.00,  0.1, 1.00 };
    Double_t plasmaHRed[plasmaHNRGBs] =   { 0.99, 0.90, 0.00 };
    Double_t plasmaHGreen[plasmaHNRGBs] = { 0.99, 0.90, 0.00 };
    Double_t plasmaHBlue[plasmaHNRGBs] =  { 0.99, 0.90, 0.00 };
    plasmaHPalette = new PPalette("plasmaH");
    plasmaHPalette->CreateGradientColorTable(plasmaHNRGBs, plasmaHStops, 
					     plasmaHRed, plasmaHGreen, plasmaHBlue, plasmaHNCont);
  }

  PPalette * electronPalette =  (PPalette*) gROOT->FindObject("electron");
  if(!electronPalette) {
    const Int_t electronNRGBs = 6;
    const Int_t electronNCont = 64;
    Double_t electronStops[electronNRGBs] = { 0.00, 0.10, 0.40, 0.55, 0.70, 1.00};
    Double_t electronRed[electronNRGBs] =   { 0.90, 0.52, 0.22, 0.39, 0.70, 1.00};
    Double_t electronGreen[electronNRGBs] = { 0.90, 0.74, 0.34, 0.05, 0.20, 1.00};
    Double_t electronBlue[electronNRGBs] =  { 0.90, 0.80, 0.58, 0.33, 0.30, 0.20};
    electronPalette = new PPalette("electron");
    electronPalette->CreateGradientColorTable(electronNRGBs, electronStops, 
					      electronRed, electronGreen, electronBlue, electronNCont);
  }

  PPalette * electroninvPalette =  (PPalette*) gROOT->FindObject("electroninv");
  if(!electronPalette) {
    const Int_t electronNRGBs = 6;
    const Int_t electronNCont = 64;
    Double_t electronStops[electronNRGBs] = { 0.00, 0.30, 0.45, 0.60, 0.90, 1.00};
    Double_t electronRed[electronNRGBs] =   { 1.00, 0.70, 0.39, 0.22, 0.52, 0.90};
    Double_t electronGreen[electronNRGBs] = { 1.00, 0.20, 0.05, 0.34, 0.74, 0.90};
    Double_t electronBlue[electronNRGBs] =  { 0.20, 0.30, 0.33, 0.58, 0.80, 0.90};
    electroninvPalette = new PPalette("electroninv");
    electroninvPalette->CreateGradientColorTable(electronNRGBs, electronStops, 
					      electronRed, electronGreen, electronBlue, electronNCont);
  }

  PPalette * electron0Palette =  (PPalette*) gROOT->FindObject("electron0");
  if(!electron0Palette) {
    const Int_t electron0NRGBs = 6;
    const Int_t electron0NCont = 64;
    Double_t electron0Stops[electron0NRGBs] = { 0.00, 0.10, 0.40, 0.55, 0.70, 1.00};
    Double_t electron0Red[electron0NRGBs] =   { 1.00, 0.52, 0.22, 0.39, 0.70, 1.00};
    Double_t electron0Green[electron0NRGBs] = { 1.00, 0.74, 0.34, 0.05, 0.20, 1.00};
    Double_t electron0Blue[electron0NRGBs] =  { 1.00, 0.80, 0.58, 0.33, 0.30, 0.20};
    electron0Palette = new PPalette("electron0");
    electron0Palette->CreateGradientColorTable(electron0NRGBs, electron0Stops, 
					      electron0Red, electron0Green, electron0Blue, electron0NCont);
  }

  PPalette * electron0invPalette =  (PPalette*) gROOT->FindObject("electron0inv");
  if(!electron0invPalette) {
    const Int_t electron0invNRGBs = 6;
    const Int_t electron0invNCont = 64;
    Double_t electron0invStops[electron0invNRGBs] = { 0.00, 0.30, 0.45, 0.60, 0.90, 1.00};
    Double_t electron0invRed[electron0invNRGBs] =   { 1.00, 0.70, 0.39, 0.22, 0.52, 1.00};
    Double_t electron0invGreen[electron0invNRGBs] = { 1.00, 0.20, 0.05, 0.34, 0.74, 1.00};
    Double_t electron0invBlue[electron0invNRGBs] =  { 0.20, 0.30, 0.33, 0.58, 0.80, 1.00};
    electron0invPalette = new PPalette("electron0inv");
    electron0invPalette->CreateGradientColorTable(electron0invNRGBs, electron0invStops, 
					      electron0invRed, electron0invGreen, electron0invBlue, electron0invNCont);
  }


  PPalette * redelectronPalette =  (PPalette*) gROOT->FindObject("redelectron");
  if(!redelectronPalette) {
    const Int_t redelectronNRGBs = 6;
    const Int_t redelectronNCont = 64;
    Double_t redelectronStops[redelectronNRGBs] = { 0.00, 0.20, 0.40, 0.50, 0.60, 1.00};
    Double_t redelectronRed[redelectronNRGBs] =   { 0.90, 0.52, 0.22, 0.39, 0.70, 1.00};
    Double_t redelectronGreen[redelectronNRGBs] = { 0.90, 0.74, 0.34, 0.05, 0.20, 1.00};
    Double_t redelectronBlue[redelectronNRGBs] =  { 0.90, 0.80, 0.58, 0.33, 0.30, 0.20};
    redelectronPalette = new PPalette("redelectron");
    redelectronPalette->CreateGradientColorTable(redelectronNRGBs, redelectronStops, 
					      redelectronRed, redelectronGreen, redelectronBlue, redelectronNCont);
  }

  PPalette * redelectron0Palette =  (PPalette*) gROOT->FindObject("redelectron0");
  if(!redelectron0Palette) {
    const Int_t redelectronNRGBs = 6;
    const Int_t redelectronNCont = 64;
    Double_t redelectronStops[redelectronNRGBs] = { 0.00, 0.20, 0.40, 0.50, 0.60, 1.00};
    Double_t redelectronRed[redelectronNRGBs] =   { 1.00, 0.52, 0.22, 0.39, 0.70, 1.00};
    Double_t redelectronGreen[redelectronNRGBs] = { 1.00, 0.74, 0.34, 0.05, 0.20, 1.00};
    Double_t redelectronBlue[redelectronNRGBs] =  { 1.00, 0.80, 0.58, 0.33, 0.30, 0.20};
    redelectron0Palette = new PPalette("redelectron0");
    redelectron0Palette->CreateGradientColorTable(redelectronNRGBs, redelectronStops, 
					      redelectronRed, redelectronGreen, redelectronBlue, redelectronNCont);
  }

  


  PPalette * electron1Palette =  (PPalette*) gROOT->FindObject("electron1");
  if(!electron1Palette) {
    const Int_t electron1NRGBs = 5;
    const Int_t electron1NCont = 64;
    Double_t electron1Stops[electron1NRGBs] = { 0.00, 0.40, 0.55, 0.70, 1.00};
    Double_t electron1Red[electron1NRGBs] =   { 0.52, 0.22, 0.39, 0.70, 1.00};
    Double_t electron1Green[electron1NRGBs] = { 0.74, 0.34, 0.05, 0.20, 1.00};
    Double_t electron1Blue[electron1NRGBs] =  { 0.80, 0.58, 0.33, 0.30, 0.20};
    electron1Palette = new PPalette("electron1");
    electron1Palette->CreateGradientColorTable(electron1NRGBs, electron1Stops, 
					      electron1Red, electron1Green, electron1Blue, electron1NCont);
  }

  PPalette *hotPalette = (PPalette*) gROOT->FindObject("hot");
  if(!hotPalette) {
    const Int_t hotNRGBs = 3;
    const Int_t hotNCont = 64;
    Double_t hotStops[hotNRGBs] = { 0.25, 0.75, 1.00 };
    Double_t hotRed[hotNRGBs] =   { 1.000, 1.000, 1.000 };
    Double_t hotGreen[hotNRGBs] = { 0.149, 0.984, 1.000 };
    Double_t hotBlue[hotNRGBs] =  { 0.000, 0.000, 1.000 };
    hotPalette = new PPalette("hot");
    hotPalette->CreateGradientColorTable(hotNRGBs, hotStops,hotRed, hotGreen, hotBlue, hotNCont);
  }    

  PPalette *hot2Palette = (PPalette*) gROOT->FindObject("hot2");
  if(!hot2Palette) {
    const Int_t hot2NRGBs = 4;
    const Int_t hot2NCont = 64;
    Double_t hot2Stops[hot2NRGBs] = { 0.00, 0.25, 0.75, 1.00 };
    Double_t hot2Red[hot2NRGBs] =   { 0.90, 1.000, 1.000, 1.000 };
    Double_t hot2Green[hot2NRGBs] = { 0.90, 0.149, 0.984, 1.000 };
    Double_t hot2Blue[hot2NRGBs] =  { 0.90, 0.000, 0.000, 1.000 };
    hot2Palette = new PPalette("hot2");
    hot2Palette->CreateGradientColorTable(hot2NRGBs, hot2Stops,hot2Red, hot2Green, hot2Blue, hot2NCont);
  }

  PPalette *hot3Palette = (PPalette*) gROOT->FindObject("hot3");
  if(!hot3Palette) {
    const Int_t hot3NRGBs = 4;
    const Int_t hot3NCont = 64;
    Double_t hot3Stops[hot3NRGBs] = { 0.00, 0.25, 0.75, 1.00 }; 
    Double_t hot3Red[hot3NRGBs] =   { 0.91, 1.000, 1.000, 0.160 }; 
    Double_t hot3Green[hot3NRGBs] = { 0.91, 0.7, 0.149, 0.160 };
    Double_t hot3Blue[hot3NRGBs] =  { 0.91, 0.000, 0.000, 0.160 };
    hot3Palette = new PPalette("hot3");
    hot3Palette->CreateGradientColorTable(hot3NRGBs, hot3Stops,hot3Red, hot3Green, hot3Blue, hot3NCont);                                       
  }

  PPalette *rbowPalette = (PPalette*) gROOT->FindObject("rbow");
  if(!rbowPalette) {
    const Int_t rbowNRGBs = 6;
    const Int_t rbowNCont = 64;
    Double_t rbowStops[rbowNRGBs] = { 0.00, 0.40, 0.48, 0.5, 0.60, 1.00 };
    Double_t rbowRed[rbowNRGBs] =   { 0.106, 0.698, 0.90, 0.90, 0.965, 0.518 };
    Double_t rbowGreen[rbowNRGBs] = { 0.078, 0.818, 0.90, 0.90, 0.925, 0.078 };
    Double_t rbowBlue[rbowNRGBs] =  { 0.518, 0.880, 0.90, 0.90, 0.353, 0.106 };
    PPalette *rbowPalette = new PPalette("rbow");
    rbowPalette->CreateGradientColorTable(rbowNRGBs, rbowStops,rbowRed, rbowGreen, rbowBlue, rbowNCont);
  }

  PPalette *rbowinvPalette = (PPalette*) gROOT->FindObject("rbowinv");
  if(!rbowinvPalette) {
    const Int_t rbowinvNRGBs = 6;
    const Int_t rbowinvNCont = 64;
    Double_t rbowinvStops[rbowinvNRGBs] = { 0.00, 0.40, 0.48, 0.50, 0.60, 1.00 };
    Double_t rbowinvRed[rbowinvNRGBs] =   { 0.518, 0.965, 0.90, 0.90, 0.698, 0.106 };
    Double_t rbowinvGreen[rbowinvNRGBs] = { 0.078, 0.925, 0.90, 0.90, 0.818, 0.078 };
    Double_t rbowinvBlue[rbowinvNRGBs] =  { 0.106, 0.353, 0.90, 0.90, 0.880, 0.518 };
    PPalette *rbowinvPalette = new PPalette("rbowinv");
    rbowinvPalette->CreateGradientColorTable(rbowinvNRGBs, rbowinvStops,rbowinvRed, rbowinvGreen, rbowinvBlue, rbowinvNCont);
  }


  PPalette *rbowwhitePalette = (PPalette*) gROOT->FindObject("rbowwhite");
  if(!rbowwhitePalette) {
    const Int_t rbowwhiteNRGBs = 6;
    const Int_t rbowwhiteNCont = 64;
    Double_t rbowwhiteStops[rbowwhiteNRGBs] = { 0.00, 0.40, 0.48, 0.5, 0.60, 1.00 };
    Double_t rbowwhiteRed[rbowwhiteNRGBs] =   { 0.106, 0.698, 1.0, 1.0, 0.965, 0.518 };
    Double_t rbowwhiteGreen[rbowwhiteNRGBs] = { 0.078, 0.818, 1.0, 1.0, 0.925, 0.078 };
    Double_t rbowwhiteBlue[rbowwhiteNRGBs] =  { 0.518, 0.880, 1.0, 1.0, 0.353, 0.106 };
    PPalette *rbowwhitePalette = new PPalette("rbowwhite");
    rbowwhitePalette->CreateGradientColorTable(rbowwhiteNRGBs, rbowwhiteStops,rbowwhiteRed, rbowwhiteGreen, rbowwhiteBlue, rbowwhiteNCont);
  }
  
  PPalette * barsaPalette = (PPalette*) gROOT->FindObject("barsa");
  if(!barsaPalette) {
    const Int_t barsaNRGBs = 3;
    const Int_t barsaNCont = 256;
    Double_t barsaStops[barsaNRGBs] = { 0.00,  0.5, 1.00 }; // { 0.00,  0.5, 1.00 }; { 0.00,  0.05, 1.00 }; 
    Double_t barsaRed[barsaNRGBs] =   { 0.0392, 1.00, 0.40 };
    Double_t barsaGreen[barsaNRGBs] = { 0.1412, 1.00, 0.00 };
    Double_t barsaBlue[barsaNRGBs] =  { 0.4157, 1.00, 0.00 };
    barsaPalette = new PPalette("barsa");
    barsaPalette->CreateGradientColorTable(barsaNRGBs, barsaStops, 
					     barsaRed, barsaGreen, barsaBlue, barsaNCont);
  }    


  PPalette *redPalette = (PPalette*) gROOT->FindObject("red");
  if(!redPalette) {
    const Int_t redNRGBs = 3;
    const Int_t redNCont = 64;
    Double_t redStops[redNRGBs] = { 0.0, 0.10, 1.00 };
    Double_t redRed[redNRGBs] =   { 0.90, 0.965, 0.518 };
    Double_t redGreen[redNRGBs] = { 0.90, 0.925, 0.078 };
    Double_t redBlue[redNRGBs] =  { 0.90, 0.353, 0.106 };
    PPalette *redPalette = new PPalette("red");
    redPalette->CreateGradientColorTable(redNRGBs, redStops,redRed, redGreen, redBlue, redNCont);
  }

  PPalette *red0Palette = (PPalette*) gROOT->FindObject("red0");
  if(!red0Palette) {
    const Int_t red0NRGBs = 3;
    const Int_t red0NCont = 64;
    Double_t red0Stops[red0NRGBs] = { 0.0, 0.10, 1.00 };
    Double_t red0Red0[red0NRGBs] =   { 1.00, 0.965, 0.518 };
    Double_t red0Green[red0NRGBs] = { 1.00, 0.925, 0.078 };
    Double_t red0Blue[red0NRGBs] =  { 1.00, 0.353, 0.106 };
    PPalette *red0Palette = new PPalette("red0");
    red0Palette->CreateGradientColorTable(red0NRGBs, red0Stops,red0Red0, red0Green, red0Blue, red0NCont);
  }

  PPalette *bluePalette = (PPalette*) gROOT->FindObject("blue");
  if(!bluePalette) {
    const Int_t blueNRGBs = 3;
    const Int_t blueNCont = 64;
    Double_t blueStops[blueNRGBs] = { 0.0, 0.10, 1.00 };
    Double_t blueRed[blueNRGBs] =   { 0.90, 0.498, 0.106};
    Double_t blueGreen[blueNRGBs] = { 0.90, 0.718, 0.078};
    Double_t blueBlue[blueNRGBs] =  { 0.90, 0.780, 0.518};
    PPalette *bluePalette = new PPalette("blue");
    bluePalette->CreateGradientColorTable(blueNRGBs, blueStops, blueRed, blueGreen, blueBlue, blueNCont);
  }

  
  // Default palette
  electron0Palette->cd();
  // cout << "Done! (default palette : electron)." << endl;
}
