#ifndef PDATA
#define PDATA

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <dirent.h>
#include <errno.h>
#include <vector>

#include <Rtypes.h>
#include <TMath.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TGraph.h>

#include <TNamed.h>

#include "PFunctions.hh"

#define DATA(i,j,k) (data[i * dim[2] * dim[1] + j * dim[2] + k])

using namespace std;

// Method for properly deleting a vector or list of pointers.
template <class C> void FreeClear( C & cntr ) {
    for ( typename C::iterator it = cntr.begin(); 
              it != cntr.end(); ++it ) {
        delete * it;
    }
    cntr.clear();
}

// External set of parameters
struct pparam {
  // Plasma parameters
  Float_t pDensity;  // Plasma density
  Float_t pStart;    // Start of plasma channel (first point of constant density).

  // Neutral specie parameters
  Float_t nDensity;  // Density of Neutral specie.
  Float_t nStart;    // Start of neutral channel (first point of constant density).
  Float_t nEnd;      // End of neutral channel (first point of constant density).
  
  // Bunch parameters
  Float_t bDensity;  // Peak density of the bunch = nb .
  Float_t bStart;    // Start of the bunch (first point at nb).
  Float_t bMass;     // Mass of the beam particles.
  Float_t bGamma;    // Mean energy of the beam in gamma units.
  Float_t bRmsZ;     // Characteristic length of the bunch on Z.
  Float_t bRmsX;     // Charasteristic transverse size of the bunch on X.
  Float_t bRmsY;     // Charasteristic transverse size of the bunch on Y.
  Float_t bRmsR;     // Charasteristic transverse size of the bunch on R.
  
  // Box parameters: These are 4D box ranges for plotting. 
  Float_t x1Min;
  Float_t x1Max;
  Float_t x2Min;
  Float_t x2Max;
  Float_t x3Min;
  Float_t x3Max;
  Float_t EMin;
  Float_t EMax;
};

// PData class definition
// -----------------------------------------------------------------

class PData : public TNamed 
{
public:
  PData(const char * name, const char * title);
  PData(const char * name);
  PData(const char * name, UInt_t time);
  
  static PData *  Get(const char * name = 0);

  virtual ~PData();

  virtual void    Clear(Option_t *option="");
  void    PrintData(Option_t *option="");
  void    SetPath(const char *path) { simPath = path; }
  void    SetTime(Int_t t) { time = t; LoadFileNames(time); }
  void    ReadParameters(const char *pfile="");
  void    LoadFileNames(Int_t t);
  void    CopyData(const char *opath, const char *cpcmd="cp -v");
  void    Delete(const char *dlcmd="rm -f");

  void    GetBoxDimensionsFromFile(const char *filename);  
  string  GetPath() { return simPath; };
  Int_t   GetTime() {  return time; } ;
  Double_t GetRealTimeFromFile(const char *filename);
  Double_t GetRealTime() { return rtime; };
  Bool_t  IsInit() { return Init; }
  Bool_t  Is3D() { return ThreeD; }
  Bool_t  IsCyl() { return Cyl; }
  
  UInt_t  NSpecies() { return species.size(); }
  string  GetSpeciesName(UInt_t i) { return species.at(i); }
  string *GetChargeFileName(UInt_t i) { return sCHG->at(i); }
  string *GetEfieldFileName(UInt_t i) { return sEF->at(i); }
  string *GetBfieldFileName(UInt_t i) { return sBF->at(i); }
  string *GetRawFileName(UInt_t i) { return sRAW->at(i); }

  UInt_t  NPhaseSpaces() { return pspaces.size(); }
  string  GetPhasespaceName(UInt_t i) { return pspaces.at(i); }
  string *GetPhasespaceFileName(UInt_t i, UInt_t j) { return sPHA->at(i).at(j); }

  inline Int_t ListDir(string dir, string pattern, vector<string> &files,string option="");
  inline void ResetParameters();
  inline Double_t Shift(TString option="");

  // Give access to exeternal parameters
  Float_t  GetPlasmaDensity() { return pParam.pDensity * (1/PUnits::cm3);  }  
  Float_t  GetPlasmaStart()   { return pParam.pStart * GetPlasmaSkinDepth(); }    
 
  // Other quantities
  Float_t  GetPlasmaFrequency() { return PFunc::PlasmaFrequency(GetPlasmaDensity()); }
  Float_t  GetPlasmaTimeDepth() { return PFunc::PlasmaTimeDepth(GetPlasmaDensity()); }
  Float_t  GetPlasmaPeriod() { return PFunc::PlasmaPeriod(GetPlasmaDensity()); }
  Float_t  GetPlasmaK() { return PFunc::PlasmaWavenumber(GetPlasmaDensity());}
  Float_t  GetPlasmaSkinDepth() { return PFunc::PlasmaSkindepth(GetPlasmaDensity());}
  Float_t  GetPlasmaWaveLength() { return PFunc::PlasmaWavelength(GetPlasmaDensity()); }
  Float_t  GetPlasmaE0() { return PFunc::PlasmaWBField(GetPlasmaDensity()); }

  // Parameters neutral specie
  Float_t  GetNeutralDensity() { return pParam.nDensity * (1/PUnits::cm3); }  
  Float_t  GetNeutralStart()   { return pParam.nStart * GetPlasmaSkinDepth(); }    
  Float_t  GetNeutralEnd()     { return pParam.nEnd * GetPlasmaSkinDepth(); }    

  // Parameters of the beam
  Float_t  GetBeamDensity()   { return pParam.bDensity * GetPlasmaDensity(); }
  Float_t  GetBeamStart()     { return pParam.bStart * GetPlasmaSkinDepth(); }      
  Float_t  GetBeamEnergy()    { return GetBeamMass() * pParam.bGamma; }
  Float_t  GetBeamMass()      { return pParam.bMass * PUnits::MeV; }
  Float_t  GetBeamGamma()     { return pParam.bGamma; }
  Float_t  GetBeamVelocity()  { return TMath::Sqrt(1. - 1./(GetBeamGamma()*GetBeamGamma()) ); } // units of c.
  Float_t  GetBeamRmsZ()    { return pParam.bRmsZ * GetPlasmaSkinDepth() ;  }
  Float_t  GetBeamRmsX()    { return pParam.bRmsX * GetPlasmaSkinDepth() ;  }
  Float_t  GetBeamRmsY()    { return pParam.bRmsY * GetPlasmaSkinDepth() ;  }
  Float_t  GetBeamRmsR()    { return pParam.bRmsR * GetPlasmaSkinDepth() ;  }

  Float_t  GetBeamBetatronK() { return PFunc::BeamBetatronWavenumber(GetPlasmaDensity(),GetBeamGamma()); }

  Float_t  GetEMin()  { return pParam.EMin; }
  Float_t  GetEMax()  { return pParam.EMax; }
 
  // Hyperslab access methods
  Double_t  GetX1Min() { return XMINR[0]; }
  Double_t  GetX1Max() { return XMAXR[0]; }
  Double_t  GetX2Min() { return XMINR[1]; }
  Double_t  GetX2Max() { return XMAXR[1]; }
  Double_t  GetX3Min() { return XMINR[2]; }
  Double_t  GetX3Max() { return XMAXR[2]; }

  Int_t  GetX1iMin() { return round((GetX1Min()-GetXMin(0))/GetDX(0)); }
  Int_t  GetX2iMin() { return round((GetX2Min()-GetXMin(1))/GetDX(1)); }
  Int_t  GetX3iMin() { return round((GetX3Min()-GetXMin(2))/GetDX(2)); }

  Int_t  GetX1iMax() { return round((GetX1Max()-GetXMin(0))/GetDX(0)); }
  Int_t  GetX2iMax() { return round((GetX2Max()-GetXMin(1))/GetDX(1)); }
  Int_t  GetX3iMax() { return round((GetX3Max()-GetXMin(2))/GetDX(2)); }

  Int_t  GetX1N() { return round((GetX1Max()-GetX1Min())/GetDX(0)); }
  Int_t  GetX2N() { return round((GetX2Max()-GetX2Min())/GetDX(1)); }
  Int_t  GetX3N() { return round((GetX3Max()-GetX3Min())/GetDX(2)); }

  void  SetX1Min(Double_t x) { XMINR[0] = x; }
  void  SetX1Max(Double_t x) { XMAXR[0] = x; }
  void  SetX2Min(Double_t x) { XMINR[1] = x; }
  void  SetX2Max(Double_t x) { XMAXR[1] = x; }
  void  SetX3Min(Double_t x) { XMINR[2] = x; }
  void  SetX3Max(Double_t x) { XMAXR[2] = x; }



  // Simulation parameters:
  Int_t    GetNDim()  { return NDIM; } 
  Int_t    GetNX(Int_t i) {return NX[i];}
  Double_t GetXMin(Int_t i) {return XMIN[i];}
  Double_t GetXMax(Int_t i) {return XMAX[i];}
  Double_t GetDX(Int_t i) {return (XMAX[i]-XMIN[i])/NX[i];}
  
  // General access to Histos in a certain H5 file.
  
  // 2D simulations (defined in z-x plane)
  TH1F* GetH1SliceZ(const char *filename, const char *dataname, 
		    Int_t Firstx2Bin = -1, Int_t Lastx2Bin = 1, const char *options="avg");
  
  TH1F* GetH1SliceX(const char *filename, const char *dataname, 
		    Int_t Firstx1Bin = -1, Int_t Lastx1Bin = 1, const char *options="avg");
  
  TH2F* GetH2(const char *filename, const char *dataname, const char *options="");
  
  // 3D simulations
  TH1F* GetH1SliceZ3D(const char *filename, const char *dataname, 
		      Int_t Firstx2Bin = -1, Int_t Lastx2Bin = 1, 
		      Int_t Firstx3Bin = -1, Int_t Lastx3Bin = 1, const char *options="avg");

  TH1F* GetH1CylSliceZ3D(const char *filename, const char *dataname, 
			 Int_t FirstrBin = -1, Int_t LastrBin = 1,
			 const char *options="avg");
  
  TH2F* GetH2SliceZX(const char *filename, const char *dataname, 
		     Int_t Firstx3Bin = -1, Int_t Lastx3Bin = 1, const char *options="avg");

  TH2F* GetH2SliceZY(const char *filename, const char *dataname, 
		     Int_t Firstx2Bin = -1, Int_t Lastx2Bin = 1, const char *options="avg");
  
  TH2F* GetH2SliceXY(const char *filename, const char *dataname, 
		     Int_t Firstx1Bin = -1, Int_t Lastx1Bin = 1, const char *options="avg");

  TH2F* GetH2ZR(const char *filename, const char *dataname, const char *options="");
  
  TH3F* GetH3(const char *filename, const char *dataname) ;

  Float_t* Get3Darray(const char *filename, const char *dataname, UInt_t dim[3]) ;
  
  // RAW data
  TTree* GetTreeRaw(const char *filename, const char *options="");

  UInt_t GetRawArray(const char *filename, Float_t **var = NULL, const char *options="");

  UInt_t GetRawSingleArray(const char *filename, Float_t **var,const char *varname = "", const char *options="");

  void  GetH1Raw(const char *filename, const char *dataname, 
		 TH1F *h1D, const char *options="");
  
  void  GetH1RawCut(const char *filename,const char *dataname, 
		    const char *cutname, Float_t cutMin, Float_t cutMax,
		    TH1F *h1D, const char *options="" );

  void  GetH1RawCut2(const char *filename,const char *dataname, 
		    const char *cutname1, Float_t cutMin1, Float_t cutMax1,
		    const char *cutname2, Float_t cutMin2, Float_t cutMax2,
		    TH1F *h1D, const char *options="" );

  void  GetH2Raw(const char *filename, const char *dataname1, const char *dataname2,
		 TH2F *h2D, const char *options="");
  
  void  GetH2RawCut(const char *filename, const char *dataname1, const char *dataname2,
		    const char *cutname, Float_t cutMin, Float_t cutMax,
		    TH2F *h2D, const char *options="");
  

  // Specific access to Histos

  // 1D Histos
  // TH1F*  GetCharge

  // 2D Histos
  TH2F*  GetCharge(UInt_t i, const char *options="") 
  { 
    return GetH2(GetChargeFileName(i)->c_str(),"charge",options); 
  } ;
  
  TH2F*  GetEField(UInt_t i, const char *options="") 
  { char nam[3]; sprintf(nam,"e%i",i+1); 
    return GetH2(GetEfieldFileName(i)->c_str(),nam,options); 
  } ;
  
  TH2F*  GetBField(UInt_t i, const char *options="") 
  { char nam[3]; sprintf(nam,"b%i",i+1); 
    return GetH2(GetBfieldFileName(i)->c_str(),nam,options); 
  } ;

  TH2F*  GetPhasespace(UInt_t i, UInt_t j,const char *options="") 
  { 
    return GetH2(GetPhasespaceFileName(i,j)->c_str(),GetPhasespaceName(j).c_str(),options); 
  } ;
  
  // 3D Histos
  TH3F* GetCharge3D(UInt_t i) 
  { 
    return GetH3(GetChargeFileName(i)->c_str(),"charge"); 
  } ;
  
  TH3F* GetEField3D(UInt_t i) 
  { char nam[3]; sprintf(nam,"e%i",i+1); 
    return GetH3(GetEfieldFileName(i)->c_str(),nam); 
  } ;
  
  TH3F* GetBField3D(UInt_t i) 
  { char nam[3]; sprintf(nam,"b%i",i+1); 
    return GetH3(GetBfieldFileName(i)->c_str(),nam); 
  } ;
  

  // ZX Slices
  TH2F* GetCharge2DSliceZX(UInt_t i, Int_t Firstx3Bin = -1, Int_t Lastx3Bin = 1, const char *options="" ) 
  { 
    return GetH2SliceZX(GetChargeFileName(i)->c_str(),"charge",Firstx3Bin,Lastx3Bin,options); 
  } ;
  
  TH2F* GetEField2DSliceZX(UInt_t i, Int_t Firstx3Bin = -1, Int_t Lastx3Bin = 1, const char *options="" ) 
  { char nam[3]; sprintf(nam,"e%i",i+1); 
    return GetH2SliceZX(GetEfieldFileName(i)->c_str(),nam,Firstx3Bin,Lastx3Bin,options); 
  } ;
  
  TH2F* GetBField2DSliceZX(UInt_t i, Int_t Firstx3Bin = -1, Int_t Lastx3Bin = 1, const char *options="" ) 
  { char nam[3]; sprintf(nam,"b%i",i+1); 
    return GetH2SliceZX(GetBfieldFileName(i)->c_str(),nam,Firstx3Bin,Lastx3Bin,options); 
  } ;

  
  // ZY Slices 
  TH2F* GetCharge2DSliceZY(UInt_t i, Int_t Firstx2Bin = -1, Int_t Lastx2Bin = 1, const char *options="" ) 
  { 
    return GetH2SliceZY(GetChargeFileName(i)->c_str(),"charge",Firstx2Bin,Lastx2Bin,options); 
  } ;
  
  TH2F* GetEField2DSliceZY(UInt_t i, Int_t Firstx2Bin = -1, Int_t Lastx2Bin = 1, const char *options="" ) 
  { char nam[3]; sprintf(nam,"e%i",i+1); 
    return GetH2SliceZY(GetEfieldFileName(i)->c_str(),nam,Firstx2Bin,Lastx2Bin,options); 
  } ;
  
  TH2F* GetBField2DSliceZY(UInt_t i, Int_t Firstx2Bin = -1, Int_t Lastx2Bin = 1, const char *options="" ) 
  { char nam[3]; sprintf(nam,"b%i",i+1); 
    return GetH2SliceZY(GetBfieldFileName(i)->c_str(),nam,Firstx2Bin,Lastx2Bin,options); 
  } ;

  // XY Slices
  TH2F* GetCharge2DSliceXY(UInt_t i, Int_t Firstx1Bin = -1, Int_t Lastx1Bin = 1, const char *options="" ) 
  { 
    return GetH2SliceXY(GetChargeFileName(i)->c_str(),"charge",Firstx1Bin,Lastx1Bin,options); 
  } ;
  
  TH2F* GetEField2DSliceXY(UInt_t i, Int_t Firstx1Bin = -1, Int_t Lastx1Bin = 1, const char *options="" ) 
  { char nam[3]; sprintf(nam,"e%i",i+1); 
    return GetH2SliceXY(GetEfieldFileName(i)->c_str(),nam,Firstx1Bin,Lastx1Bin,options); 
  } ;
  
  TH2F* GetBField2DSliceXY(UInt_t i, Int_t Firstx1Bin = -1, Int_t Lastx1Bin = 1, const char *options="" ) 
  { char nam[3]; sprintf(nam,"b%i",i+1); 
    return GetH2SliceXY(GetBfieldFileName(i)->c_str(),nam,Firstx1Bin,Lastx1Bin,options); 
  } ;

  // ZR Slices

  TH2F* GetChargeR(UInt_t i, const char *option="")  
  { 
    return GetH2ZR(GetChargeFileName(i)->c_str(),"charge",option);
  } ;
  
  TH2F* GetEFieldR(UInt_t i, const char *option="") 
  { char nam[3]; sprintf(nam,"e%i",i+1); 
    return GetH2ZR(GetEfieldFileName(i)->c_str(),nam,option); 
  } ;
  
protected:

  static PData             *fgData;

  Int_t                       time;   // Current file time step.
  Double_t                   rtime;   // Current time.
  Int_t                       NDIM;   // Dimension of the simulation.
  Int_t                        *NX;   // Number of cells in every dimension.
  Double_t                   *XMIN;   // Low edges of simulation box.
  Double_t                   *XMAX;   // Up edges of simulation box.
  Double_t                  *XMINR;   // Low edges of simulation box (zoom).
  Double_t                  *XMAXR;   // Up edges of simulation box (zoom).
  string                   simPath;   // Path to the simulation directory.
  Bool_t                      Init;   // Flag for initialization.
  Bool_t                       Cyl;   // Flag for cylindrical coordinates.
  Bool_t                    ThreeD;   // Flag for 3D sim.

  vector<string>           species;   // vector of species names.
  vector<string*>            *sCHG;   // vector of files with the Charge density.
  vector<string*>             *sEF;   // vector of files with the Electric field.
  vector<string*>             *sBF;   // vector of files with the Magnetic field. 
  vector<string*>            *sRAW;   // vector of files with the RAW data.
  vector<string>           pspaces;   // vector of phase spaces names.
  vector<vector<string*> >   *sPHA;   // vector of files with PHASESPACE info.

  pparam                    pParam;   // Struct with simulation parameters (see above).

  ClassDef(PData,1) // Singleton class for global configuration settings
};

extern PData *gData;


//______________________________________________________________________________________
Int_t PData::ListDir(string dir, string pattern, vector<string> &files,string option)
{
  //  unsigned char isFile =0x8;
  
  DIR *dp;
  struct dirent *dirp;
  if((dp  = opendir(dir.c_str())) == NULL) {
    //    cout << "Error(" << errno << ") opening " << dir << endl;
    return errno;
  }  

  dir += "/"; 

  while ((dirp = readdir(dp)) != NULL) {
    string filename = string(dirp->d_name);

    if(filename == "." || filename == "..") continue;

    //    if( (dirp->d_type != isFile) && (option.find("rec") != string::npos) ) {
    if( (filename.find(".h5") == string::npos) && (option.find("rec") != string::npos) ) {
      ListDir(dir + filename,pattern,files,option);
    }

    size_t found = filename.find(pattern);
    if(found != string::npos) {
      
      if(option.find("nam") != string::npos) 
	files.push_back(filename);
      else 
	files.push_back(dir + filename);
      
    }
  }

  closedir(dp);
  return 0;
}

//______________________________________________________________________________________
void PData::ResetParameters() {
  pParam.pDensity = pParam.nDensity = 1;

  pParam.pStart = pParam.nStart = pParam.nEnd
    = pParam.bDensity = pParam.bStart 
    = pParam.bRmsZ = pParam.bRmsX  = pParam.bRmsY    = pParam.bRmsR 
    = pParam.EMin  = pParam.EMax   =  0.;
  
  pParam.x1Min = pParam.x1Max  = pParam.x2Min = pParam.x2Max 
    = pParam.x3Min = pParam.x3Max = -999.;
  
  // Default plasma density
  pParam.pDensity = 1e17;  // particles per cm3. 
  
  // Default beam mass
  pParam.bMass = 0.511;    // In MeV.

  // Default beam gamma
  pParam.bGamma = 1.0;
}

//______________________________________________________________________________________
Double_t PData::Shift(TString option) {
  TString opt = option;

  Double_t shiftx1 = 0;
  if(opt.Contains("center")) {
    Double_t kp = GetPlasmaK();
    if(opt.Contains("comov"))        // Centers on the head of the beam.
      shiftx1 += GetBeamStart()*kp;   // Centers on the start of the plasma channel.
    else
      shiftx1 += GetPlasmaStart()*kp;
  }
  
  if(opt.Contains("comov")) {
    Double_t v = GetBeamVelocity();    
    if(v==0) v = 1.0; // If equals to 0 (default), then set to c.
    shiftx1 += v * GetRealTime();
  }   

  return shiftx1;
}

#endif