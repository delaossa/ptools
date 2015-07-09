#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>


#include <TROOT.h>
#include <TSystem.h>


using namespace std;

void PlotTest2(Int_t Nx, Int_t Ny) {

  
  Float_t **array = NULL;
  
  array = new Float_t*[Nx];
  for(Int_t i=0;i<Nx;i++) {
    array[i] = new Float_t[Ny];
  }

  
  for(Int_t i=0;i<Nx;i++) {
    for(Int_t j=0;j<Ny;j++) {
      array[i][j] = (i+1)*(j+1);
      cout << Form(" (%i,%i) element = %f",i,j,array[i][j]) << endl;
    }
  }
  
  
}
