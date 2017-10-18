#include <stdlib.h>
#include <stdio.h>     
#include <math.h>      
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>

#include <TRandom3.h>
#include <PUnits.hh>
#include <PFunctions.hh>

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

  const double c_light = 299792458; // speed of light [m/s]
  const double e_mass = 0.510998910; // MeV
  
  // Option for emittance spoiler
  // const double X0 = 35.28E4; // um (Beryllium)
  const float X0 = 8.897E4; // um (Aluminum)
  float X = 1000.; // um default value
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
    } else if(arg.find("--astra") != string::npos) {
      opt += "astra";
    } else {
      filename = arg;
    }
  }

    
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


  // --------------------------------------------------
  // READ FROM TEXT ELEGANT OR ASTRA FILE (ASCII mode)
  if(filename.find(".h5")==string::npos) {

    ifstream file(filename.c_str());

    if(opt.find("astra")==string::npos) {
      printf("\n 1. Reading ELEGANT file (ascii mode) .. \n");
      if(opt.find("spoil")!= string::npos) 
	printf("\n  -> Spoiling emittance with %.1f um of Aluminum foil ...\n",X);
      
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
	  
	    double E0 = var[3][irow-2] * e_mass;
	    double Theta0 = (13.6/E0) * sqrt(X/X0) * (1.0 - 0.038 * log(X/X0)); // rad
	    double xr1, xr2;
	    rndEngine->Rannor(xr1,xr2);
	    double xplane  = ((xr1 * X * Theta0 / sqrt(12.)) + (xr2 * X * Theta0 / 2)) * 1E-6; // m
	    double txplane = xr2 * Theta0; // rad
	    var[1][irow-2] += xplane;
	    var[4][irow-2] += txplane;
	  
	    double yr1, yr2;
	    rndEngine->Rannor(yr1,yr2);
	    double yplane  = (yr1 * X * Theta0 / sqrt(12.) + yr2 * X * Theta0 / 2) * 1E-6; // m
	    double typlane = yr2 * Theta0; // rad
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
    } else {
      printf("\n 1. Reading ASTRA file (ascii mode) .. \n");
          
      string str; 
      Int_t irow = 0;
      double t,q;
      Int_t pi,ps;

      // Count the lines
      while (std::getline(file, str)) ++Np;
      for(Int_t i=0;i<Nvar;i++)
	var[i] = new double[Np];

      // Rewind
      file.clear();
      file.seekg(0);
      while (std::getline(file, str)) {
	istringstream stream(str);

	//      cout << str << endl;
      
	for(Int_t i=0;i<Nvar;i++) {
	  stream >> var[i][irow];	
	}
      
	stream >> t;
	stream >> q;
	q *= -1;
	stream >> pi;
	stream >> ps;
      
	if(pi!=1) continue;
	if(ps!=5) continue;
	if(ps!=5) continue;
	if(var[2][irow]<-1) continue;
	
	// transform spatial coordinates
	double aux0 = var[0][irow];
	double aux1 = var[1][irow];
	var[0][irow] = var[2][irow];  // m;
	var[1][irow] = aux0;  // m;
	var[2][irow] = aux1;  // m;

	// transform momentum coordinates
	aux0 = var[3][irow];
	aux1 = var[4][irow];
	var[3][irow] = var[5][irow] / (e_mass*1E6); // mc
	var[4][irow] = aux0 / (e_mass*1E6); // mc
	var[5][irow] = aux1 / (e_mass*1E6); // mc
      
	if(irow>0) {
	  var[0][irow] += var[0][0];
	  var[3][irow] += var[3][0];
	}
      
	Q += q;

	for(Int_t i=0;i<Nvar;i++) {
	  varMean[i] += var[i][irow];
	  varRms[i]  += var[i][irow]*var[i][irow];	
	
	  if(var[i][irow]<varMin[i]) varMin[i] = var[i][irow];
	  if(var[i][irow]>varMax[i]) varMax[i] = var[i][irow];
	}
      
	irow++;
      }

      Q *= 1E3; // pC
      Np = irow;
    }
  } else {  // Reading HDF5 file (ELEGANT ONLY)
    
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

    DataSet *varDataSet[Nvar];
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

  }
    
  for(int i=0;i<Nvar;i++) {
    varMean[i] /= Np;
    varRms[i]  /= Np;
    varRms[i] = sqrt( varRms[i] -  varMean[i]*varMean[i] );
  }
    
  printf("\n  %i  particles read!    Charge = %.1f pC" , Np, Q);
  printf("\n  x1 = %e +/- %e \n  x2 = %e +/- %e \n  x3 = %e +/- %e \n  p1 = %e +/- %e \n  p2 = %e +/- %e \n  p3 = %e +/- %e",
	 varMean[0], varRms[0], varMean[1], varRms[1], varMean[2], varRms[2], varMean[3], varRms[3], varMean[4], varRms[4], varMean[5], varRms[5]);      
  
  
  if(opt.find("hdf")!=string::npos) {
    printf("\n\n 2. Writing OSIRIS RAW particle file (HDF5) .. \n");
    
    double n0  = 1E22;
    double kp  = PFunc::PlasmaWavenumber(n0);
    cout << Form(" Plasma density = %.4e cm^-3  --> skindepth = %f um",n0*PUnits::cm3,
    		 (1/kp)/PUnits::um) << endl;

    // Normalized charge 
    double Q0 = PUnits::echarge * n0 / (kp*kp*kp);

    // Charge per particle
    double q = (Q*PUnits::picocoulomb/Np)/Q0;
    
    double *charge = new double[Np];
    varMean[0] *= kp; 
    varMean[1] *= kp; 
    varMean[2] *= kp;
    for(unsigned int j=0;j<Np;j++) {
      var[0][j] *= kp;
      var[1][j] *= kp;
      var[2][j] *= kp;

      var[0][j] -= varMean[0];
      var[1][j] -= varMean[1];
      var[2][j] -= varMean[2];

      charge[j] = -q;
    }
    
    char iext[16] = ".txt";  
    char oext[16] = ".osi.h5";

    string ofilename = filename;
    ofilename.replace(ofilename.find(iext),string::npos,oext);

    H5File file( ofilename, H5F_ACC_TRUNC);

    hsize_t dimsf[1]; // dataset dimensions
    dimsf[0] = Np;
    DataSpace dataspace( 1, dimsf );
    
    FloatType datatype( PredType::NATIVE_DOUBLE );

    DataSet varDataSet[Nvar];
    for(int i=0;i<Nvar;i++) {
      varDataSet[i] = file.createDataSet(varname[i], datatype, dataspace);      
      varDataSet[i].write(var[i], PredType::NATIVE_DOUBLE);   
    }
    DataSet cDataSet = file.createDataSet("q", datatype, dataspace);
    cDataSet.write(charge,PredType::NATIVE_DOUBLE);
    
    file.close();
    
  } else {
    printf("\n\n 2. Writing OSIRIS RAW particle file (text mode) .. \n");
    
    // Normalized charge per particle
    Q /= -Q;
    
    char iext[16] = ".txt";  
    char oext[16] = ".osi";
 
    string ofilename = filename;
    ofilename.replace(ofilename.find(iext),string::npos,oext);
    
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
