#include <stdlib.h>
#include <stdio.h>     
#include <math.h>      
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>

#include <TFile.h>
#include <TTree.h>

#include "H5Cpp.h"

using namespace std;
using namespace H5;

int main(int argc,char *argv[]) {
  if(argc<2) {
    printf("\n Usage: %s <input file> \n",argv[0]);
    return 0;
  }

  // Input file
  string filename = "";

  // General options
  string opt = "";

  const float c_light = 299792458; // speed of light [m/s]
  const float e_mass = 0.510998910; // MeV
    
  // Interfacing command line:
  for(int l=1;l<argc;l++){
    string arg = argv[l];

    if(arg.find("--inv") != string::npos) {
      opt += "inv";
    } else if(arg.find("--astra") != string::npos) {
      opt += "astra";
    } else {
      filename = arg;
    }
  }

  // --------------------------------------------------
  // READ FROM TEXT ELEGANT OR ASTRA FILE
  ifstream file(filename.c_str());
    
  double Q = 0.0;
  unsigned int Np = 0;  
  const int Nvar = 6;
  double *var[Nvar]; 
  double varMean[Nvar]; 
  double varRms[Nvar]; 
  double varMin[Nvar]; 
  double varMax[Nvar]; 
  char varname[Nvar][8] = {{"x1"},{"x2"},{"x3"},{"p1"},{"p2"},{"p3"}};
  for(int i=0;i<Nvar;i++) {
    varMean[i] = 0.0;
    varRms[i] = 0.0;
    varMin[i] = 1E20;
    varMax[i] = -1E20;
  }

  printf("\n 1. Reading ELEGANT file (ascii mode) .. \n");
      
  string str; 
  int irow = 0;
  while (std::getline(file, str)) {
    std::stringstream stream(str);
      
    // Process str    
    if(irow==0) {
      stream >> Q;
      Q *= 1E12;   // pC
    } else if(irow==1) {
      stream >> Np;
      for(int i=0;i<Nvar;i++)
	var[i] = new double[Np];
    } else {
      for(int i=0;i<Nvar;i++) {
	stream >> var[i][irow-2];	
      }
      
      var[0][irow-2] *= -c_light; // m
      var[1][irow-2] *= 1;     // m
      var[2][irow-2] *= 1;     // m	
      var[4][irow-2] *= var[3][irow-2]; // mc
      var[5][irow-2] *= var[3][irow-2]; // mc
      
      for(int i=0;i<Nvar;i++) {
	varMean[i] += var[i][irow-2];
	varRms[i]  += var[i][irow-2] * var[i][irow-2];	
	
	if(var[i][irow-2]<varMin[i]) varMin[i] = var[i][irow-2];
	if(var[i][irow-2]>varMax[i]) varMax[i] = var[i][irow-2];
	
      }
    }
    irow++;
  }    
  
  for(int i=0;i<Nvar;i++) {
    varMean[i] /= Np;
    varRms[i]  /= Np;
    varRms[i] = sqrt( varRms[i] - varMean[i]*varMean[i] );
  }
    
  printf("\n  %i  particles read!    Charge = %.1f pC" , Np, Q);
  printf("\n  x1 = %e +/- %e \n  x2 = %e +/- %e \n  x3 = %e +/- %e \n  p1 = %e +/- %e \n  p2 = %e +/- %e \n  p3 = %e +/- %e",
	 varMean[0], varRms[0], varMean[1], varRms[1], varMean[2], varRms[2], varMean[3], varRms[3], varMean[4], varRms[4], varMean[5], varRms[5]);      
  
    
  printf("\n\n 2. Writing to TTree.. \n");
  
  char iext[5] = ".txt";
  char oext[6] = ".root";
  string ofilename = filename;
  ofilename.replace(ofilename.find(iext),5,oext);

  TFile *ofile = new TFile(ofilename.c_str(),"RECREATE");
  TTree *tree = new TTree("tree","");

  // Define the branches
  Double_t *darray = new Double_t[Nvar];
  for(UInt_t i=0;i<Nvar;i++) {
    tree->Branch(varname[i],&darray[i],Form("%s/D",varname[i]));
  }


  for(UInt_t j=0;j<Np;j++) {
    for(UInt_t i=0;i<Nvar;i++)
      darray[i] = var[i][j];
    
    tree->Fill();
  }
  
  ofile->Write();
  ofile->Close();
  
  return 0;
  
  // ----------------------------------------------------------------------------------------------
  
}
