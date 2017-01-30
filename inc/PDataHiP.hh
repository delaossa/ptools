#ifndef PDATAHIP
#define PDATAHIP

#include "PData.hh"

// PData class definition
// -----------------------------------------------------------------

class PDataHiP : public PData 
{
public:

  PDataHiP(const char * name, const char * title);
  PDataHiP(const char * name);
  PDataHiP(const char * name, UInt_t time);

  virtual ~PDataHiP();

  static PDataHiP *  Get(const char * name = 0);
  void     Clear(Option_t *option="");
  void     LoadFileNames(Int_t t);
  void     PrintData(Option_t *option="");
  void     ReadOutputSummary(const char *pfile="");
  Double_t GetRealTimeFromFile(const char *filename);

  UInt_t   NRawSpecies();
  string   GetRawSpeciesName(UInt_t i);

  string  *GetWfieldFileName(UInt_t i);

  TH1F* GetH1SliceZ3D(const char *filename, const char *dataname, 
		      Int_t Firstx2Bin = -1, Int_t Lastx2Bin = 1, 
		      Int_t Firstx3Bin = -1, Int_t Lastx3Bin = 1, const char *options="avg");
  
  TH2F*    GetH2SliceZX(const char *filename, const char *dataname, 
		     Int_t Firstx3Bin = -1, Int_t Lastx3Bin = 1, const char *options="avg");

  TH2F*    GetH2SliceZY(const char *filename, const char *dataname, 
			Int_t Firstx2Bin = -1, Int_t Lastx2Bin = 1, const char *options="avg");
  
  Double_t Shift(TString option="");
  Double_t ShiftT(TString option="");
  
  TH2F* GetEField2DSliceZX(UInt_t i, Int_t Firstx3Bin = -1, Int_t Lastx3Bin = 1, const char *options="avg" ) 
  { char nam[3];
    if(i==0) sprintf(nam,"Ez");
    else return NULL;
    return GetH2SliceZX(GetEfieldFileName(i)->c_str(),nam,Firstx3Bin,Lastx3Bin,options); 
  } 
  
  TH2F* GetBField2DSliceZX(UInt_t i, Int_t Firstx3Bin = -1, Int_t Lastx3Bin = 1, const char *options="avg" ) 
  { char nam[3];
    if(i==1) sprintf(nam,"Bx");
    else if(i==2) sprintf(nam,"By");
    return GetH2SliceZX(GetBfieldFileName(i)->c_str(),nam,Firstx3Bin,Lastx3Bin,options); 
  } 
  
  TH2F*  GetWField2DSliceZX(UInt_t i, Int_t Firstx3Bin = -1, Int_t Lastx3Bin = 1, const char *options="avg" ) 
  { char nam[3];
    if(i==0) sprintf(nam,"ExmBy");
    else if(i==1) sprintf(nam,"EypBx");
    return GetH2SliceZX(GetWfieldFileName(i)->c_str(),nam,Firstx3Bin,Lastx3Bin,options); 
  } 

  TH2F*  GetWField(UInt_t i, const char *options="") 
  { char nam[3]; 
    if(i==0) sprintf(nam,"ExmBy");
    else if(i==1) sprintf(nam,"EypBx");
    return GetH2(GetWfieldFileName(i)->c_str(),nam,options); 
  }

  TH2F* GetEField2DSliceZY(UInt_t i, Int_t Firstx3Bin = -1, Int_t Lastx3Bin = 1, const char *options="avg" ) 
  { char nam[3];
    if(i==0) sprintf(nam,"Ez");
    else return NULL;
    return GetH2SliceZY(GetEfieldFileName(i)->c_str(),nam,Firstx3Bin,Lastx3Bin,options); 
  } 
  
  TH2F* GetBField2DSliceZY(UInt_t i, Int_t Firstx3Bin = -1, Int_t Lastx3Bin = 1, const char *options="avg" ) 
  { char nam[3];
    if(i==1) sprintf(nam,"Bx");
    else if(i==2) sprintf(nam,"By");
    return GetH2SliceZY(GetBfieldFileName(i)->c_str(),nam,Firstx3Bin,Lastx3Bin,options); 
  } 
  
  TH2F*  GetWField2DSliceZY(UInt_t i, Int_t Firstx3Bin = -1, Int_t Lastx3Bin = 1, const char *options="avg" ) 
  { char nam[3];
    if(i==0) sprintf(nam,"ExmBy");
    else if(i==1) sprintf(nam,"EypBx");
    return GetH2SliceZY(GetWfieldFileName(i)->c_str(),nam,Firstx3Bin,Lastx3Bin,options); 
  } 
  
  
protected:

  vector<string>      rawspecies;   // vector of RAW species names.
  vector<string*>           *sWF;   // vector of files with the transverse wakefields.

  ClassDef(PDataHiP,1) 
};



#endif
