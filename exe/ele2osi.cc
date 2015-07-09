#include <stdlib.h>
#include <stdio.h>     
#include <math.h>      
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>

#include <TRandom3.h>

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
  
  // Option for emittance spoiler
  // const Float_t X0 = 35.28E4; // um (Beryllium)
  const Float_t X0 = 8.897E4; // um (Aluminum)
  Float_t X = 1000.; // um default value
  TRandom3 *rndEngine = new TRandom3();
  
  // Interfacing command line:
  for(int l=1;l<argc;l++){
    string arg = argv[l];

    if(arg.find("--hdf") != string::npos) {
      opt += "hdf";
    } else if(arg.find("-spoil") != string::npos) {
      char ss[6];
      sscanf(arg.c_str(),"%6s%f",ss,&X);
      opt += "spoil";
    } else if(arg.find("--inv") != string::npos) {
      opt += "inv";
    } else {
      filename = arg;
    }
  }

  if(filename.find(".h5")==string::npos) {
    
    printf("\n 1. Reading ELEGANT file (ascii mode) .. \n");
    if(opt.find("spoil")!= string::npos) 
      printf("\n  -> Spoiling emittance with %.1f um of Aluminum foil ...\n",X);
    
    // --------------------------------------------------
    // READ FROM TEXT ELEGANT FILE
    ifstream file(filename.c_str());
    
    double Q;
    unsigned int Np;  
    const int Nvar = 6;
    double *var[Nvar]; 
    double varMean[Nvar]; 
    double varRms[Nvar]; 
    double varMin[Nvar]; 
    double varMax[Nvar]; 
    //  char varname[Nvar][8] = {{"x1"},{"x2"},{"x3"},{"p1"},{"p2"},{"p3"}};
    for(int i=0;i<Nvar;i++) {
      varMean[i] = 0.0;
      varRms[i] = 0.0;
      varMin[i] = 1E20;
      varMax[i] = -1E20;
    }
    
    string str; 
    int irow = 0;
    while (std::getline(file, str)) {
      std::stringstream stream(str);
      
      // Process str    
      if(irow==0) {
	stream >> Q;
	// Q *= 1E12;   // pC
      } else if(irow==1) {
	stream >> Np;
	for(int i=0;i<Nvar;i++)
	  var[i] = new double[Np];
      } else {
	for(int i=0;i<Nvar;i++) {
	  stream >> var[i][irow-2];	
	}
	
	// transform variables
	if(opt.find("inv")!= string::npos)
	  var[0][irow-2] *= c_light; // m
	else
	  var[0][irow-2] *= -c_light; // m
	  
	var[1][irow-2] *= 1;     // m
	var[2][irow-2] *= 1;     // m

	// Emittance spoiler:
	if(opt.find("spoil")!= string::npos) {
	  // \Theta_0 from the multiple scattering model:
	  // http://pdg.lbl.gov/2014/reviews/rpp2014-rev-passage-particles-matter.pdf
	  
	  Float_t E0 = var[3][irow-2] * e_mass;
	  Float_t Theta0 = (13.6/E0) * sqrt(X/X0) * (1.0 - 0.038 * log(X/X0)); // rad
	  Float_t xr1, xr2;
	  rndEngine->Rannor(xr1,xr2);
	  Float_t xplane  = ((xr1 * X * Theta0 / sqrt(12.)) + (xr2 * X * Theta0 / 2)) * 1E-6; // m
	  Float_t txplane = xr2 * Theta0; // rad
	  var[1][irow-2] += xplane;
	  var[4][irow-2] += txplane;
	  
	  Float_t yr1, yr2;
	  rndEngine->Rannor(yr1,yr2);
	  Float_t yplane  = (yr1 * X * Theta0 / sqrt(12.) + yr2 * X * Theta0 / 2) * 1E-6; // m
	  Float_t typlane = yr2 * Theta0; // rad
	  var[2][irow-2] += yplane;
	  var[5][irow-2] += typlane;
	}
	
	var[4][irow-2] *= var[3][irow-2]; // mc
	var[5][irow-2] *= var[3][irow-2]; // mc
	
	for(int i=0;i<Nvar;i++) {
	  varMean[i] += var[i][irow-2];
	  varRms[i]  += var[i][irow-2]*var[i][irow-2];	
	  
	  if(var[i][irow-2]<varMin[i]) varMin[i] = var[i][irow-2];
	  if(var[i][irow-2]>varMax[i]) varMax[i] = var[i][irow-2];
	  
	}
      }
      irow++;
    }
    
    file.close();
    
    for(int i=0;i<Nvar;i++) {
      varMean[i] /= Np;
      varRms[i]  /= Np;
      varRms[i] = sqrt( varRms[i] -  varMean[i]*varMean[i] );
    }
    
    printf("\n  %i  particles read!    Charge = %.1f pC" , Np, Q*1E12);
    printf("\n  x1 = %e +/- %e \n  x2 = %e +/- %e \n  x3 = %e +/- %e \n  p1 = %e +/- %e \n  p2 = %e +/- %e \n  p3 = %e +/- %e",
	   varMean[0], varRms[0], varMean[1], varRms[1], varMean[2], varRms[2], varMean[3], varRms[3], varMean[4], varRms[4], varMean[5], varRms[5]);      
    
    
    printf("\n\n 2. Writing OSIRIS file (ascii mode) .. \n");
    
    // Normalized charge per particle
    Q /= -Q;
    
    char iext[5] = ".txt";
    char oext[5] = ".osi";
    string ofilename = filename;
    ofilename.replace(ofilename.find(iext),5,oext);
    
    ofstream outfile(ofilename.c_str(),ios::out | ios::trunc);
    char ocstr[256];
    for(unsigned int i=0;i<Np;i++) {
      
      sprintf(ocstr,"%15.8e %15.8e %15.8e %10.5f %15.8e %15.8e %2i",var[0][i] - varMean[0],var[1][i] - varMean[1],var[2][i] - varMean[2],var[3][i],var[4][i],var[5][i],(int)Q);
      
      outfile << ocstr << endl;
      
    }
    
    outfile.close();

  } else {

    printf("\n 1. Reading ELEGANT file (hdf5 mode) .. \n");

    // --------------------------------------------------
    // READ FROM H5 ELEGANT FILE
    H5File h5file = H5File(filename,H5F_ACC_RDONLY);
    
    // Open "parameters" group 
    Group *paraGroup = new Group(h5file.openGroup("/page1/parameters"));
    
    // Total charge
    DataSet *chDataSet = new DataSet(paraGroup->openDataSet("Charge"));
    DataSpace chDataSpace = chDataSet->getSpace();
    int rank = chDataSpace.getSimpleExtentNdims();  
    hsize_t *dims = new hsize_t[rank];
    chDataSpace.getSimpleExtentDims(dims,NULL);
    DataSpace chmemSpace(rank,dims);
    
    // Preparing buffer
    const DataType chType = chDataSet->getDataType();
    double Q;
    
    // Reading
    chDataSet->read(&Q,chType,chmemSpace,chDataSpace);
    chDataSet->close();
    delete [] dims;
    
    printf("\n Total charge = %e",Q);
    printf("\n");

    // Open "columns" group 
    Group *dataGroup = new Group(h5file.openGroup("/page1/columns"));

    unsigned int Np;  
    const int Nvar = 6;
    DataSet *varDataSet[Nvar];
    double *var[Nvar];
    varDataSet[0] = new DataSet(dataGroup->openDataSet("t"));
    varDataSet[1] = new DataSet(dataGroup->openDataSet("x"));
    varDataSet[2] = new DataSet(dataGroup->openDataSet("y"));
    varDataSet[3] = new DataSet(dataGroup->openDataSet("p"));
    varDataSet[4] = new DataSet(dataGroup->openDataSet("xp"));
    varDataSet[5] = new DataSet(dataGroup->openDataSet("yp"));
    for(int i=0;i<Nvar;i++) {
      DataSpace varDataSpace = varDataSet[i]->getSpace();
      int rank = varDataSpace.getSimpleExtentNdims();  
      hsize_t *dims = new hsize_t[rank];
      varDataSpace.getSimpleExtentDims(dims,NULL);
      DataSpace memSpace(rank,dims);
      
      // Preparing buffer
      const DataType Type = varDataSet[i]->getDataType();
      Np = dims[0];

      var[i] = new double[Np];
      
      // Reading
      varDataSet[i]->read(var[i],Type,memSpace,varDataSpace);
      varDataSet[i]->close();

      delete [] dims; 
    }

    // Process data and dump into file:
    double varMean[Nvar]; 
    double varRms[Nvar]; 
    double varMin[Nvar]; 
    double varMax[Nvar]; 
    for(int i=0;i<Nvar;i++) {
      varMean[i] = 0.0;
      varRms[i] = 0.0;
      varMin[i] = 1E20;
      varMax[i] = -1E20;
    }
    
    for(unsigned int ip=0;ip<Np;ip++) {
      // transform variables
      var[0][ip] *= c_light; // m
      var[1][ip] *= 1;     // m
      var[2][ip] *= 1;     // m
      
      var[4][ip] *= var[3][ip]; // mc
      var[5][ip] *= var[3][ip]; // mc
      
      for(int i=0;i<Nvar;i++) {
	varMean[i] += var[i][ip];
	varRms[i]  += var[i][ip]*var[i][ip];	
	
	if(var[i][ip]<varMin[i]) varMin[i] = var[i][ip];
	if(var[i][ip]>varMax[i]) varMax[i] = var[i][ip];
	
      }
      
    }

    for(int i=0;i<Nvar;i++) {
      varMean[i] /= Np;
      varRms[i]  /= Np;
      varRms[i] = sqrt( varRms[i] -  varMean[i]*varMean[i] );
    }
    
    printf("\n  %i  particles read!    Charge = %.1f pC" , Np, Q*1E12);
    printf("\n  x1 = %e +/- %e \n  x2 = %e +/- %e \n  x3 = %e +/- %e \n  p1 = %e +/- %e \n  p2 = %e +/- %e \n  p3 = %e +/- %e",
	   varMean[0], varRms[0], varMean[1], varRms[1], varMean[2], varRms[2], varMean[3], varRms[3], varMean[4], varRms[4], varMean[5], varRms[5]);      

    printf("\n\n 2. Writing OSIRIS file (ascii mode) .. \n");
    
    // Normalized charge per particle
    Q /= -Q;
    
    char iext[5] = ".h5";
    char oext[5] = ".osi";
    string ofilename = filename;
    ofilename.replace(ofilename.find(iext),5,oext);
    
    ofstream outfile(ofilename.c_str(),ios::out | ios::trunc);
    char ocstr[256];
    for(unsigned int i=0;i<Np;i++) {
      
      sprintf(ocstr,"%15.8e %15.8e %15.8e %10.5f %15.8e %15.8e %2i",var[0][i] - varMean[0],var[1][i] - varMean[1],var[2][i] - varMean[2],var[3][i],var[4][i],var[5][i],(int)Q);
      
      outfile << ocstr << endl;
      
    }
    
    outfile.close();
        
  }
  
    return 0;
  
  // ----------------------------------------------------------------------------------------------
    
}
