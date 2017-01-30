#include "PDataHiP.hh"

#include <H5Cpp.h>

ClassImp(PDataHiP);

using namespace H5; 

std::string between(std::string const &in,
                    std::string const &before, std::string const &after) {
  size_t beg = in.find(before);
  beg += before.size();
  size_t end = in.find(after, beg);
  return in.substr(beg, end-beg);
}

//_______________________________________________________________________
PDataHiP::PDataHiP(const char * name, const char * title) : PData(name,title) {
  rawspecies.clear();

  sWF = NULL;
}

//_______________________________________________________________________
PDataHiP::PDataHiP(const char * name) : PData(name) {
  rawspecies.clear();

  sWF = NULL;  
}

//_______________________________________________________________________
PDataHiP::PDataHiP(const char * name, UInt_t time) : PData(name, time) {
  rawspecies.clear();

  sWF = NULL;  
}

//_______________________________________________________________________
PDataHiP::~PDataHiP()
{
  Clear();
  
  fgData = 0;
  gData  = 0;
}


//_______________________________________________________________________
PDataHiP * PDataHiP::Get(const char * name)
{
  if (!fgData) {
    if (name) {
      fgData = new PDataHiP(name);
    } else {
      fgData = new PDataHiP("");
    }
  }

  if (name && strcmp(fgData->GetName(), name)!=0) {
    cout << "PDataHiP: Changing analysis to " << name << endl;

    if(fgData) delete fgData;
    fgData = new PDataHiP(name);
  }
  
  return dynamic_cast<PDataHiP *>(fgData);
}

//_______________________________________________________________________
void PDataHiP::Clear(Option_t *option) {
  PData::Clear(option);

  if(sWF) {
    FreeClear(*sWF);
    delete sWF;
    sWF = NULL;
  }

  rawspecies.clear();
  
}
//_______________________________________________________________________
void PDataHiP::LoadFileNames(Int_t t) {
  
  if(time == t && Init) {
    // cout << "Time step already loaded: doing nothing." << endl;
    cout << " ... " << endl << endl;
    return;  
  }
  
  Clear();
  
  time = t;

  species.push_back("plasma");
  species.push_back("beam");

  // Time
  char stime[10];
  if(time<1e1) sprintf(stime,"00000%1i",time);
  else if(time<1e2) sprintf(stime,"0000%2i",time);
  else if(time<1e3) sprintf(stime,"000%3i",time);
  else if(time<1e4) sprintf(stime,"00%4i",time);
  else if(time<1e5) sprintf(stime,"0%5i",time);
  else if(time<1e6) sprintf(stime,"%6i",time);

  // Get the vector of files to be analized at certain time value.
  vector<string> files;
  string Dir = simPath + "/DATA";
  ListDir(Dir,string(stime),files,"rec");
  
  if(files.size()==0) {
    cout << "PDataHiP:: No files found in time step ( " << time << " ) : doing nothing!" << endl;
    Init = kFALSE;
    return;
  }

  
  // Initialize the pointers to NULLS
  sCHG = new vector<string*>(NSpecies(),NULL);
  sEF  = new vector<string*>(3,NULL);
  sBF  = new vector<string*>(3,NULL);
  sWF  = new vector<string*>(2,NULL);
  sPHA = new vector<vector<string*> >(NSpecies(),vector<string*>(NPhaseSpaces(),NULL));

  for(UInt_t i=0;i<3;i++) {
    sJ[i] = new vector<string*>(NSpecies(),NULL);    
  }
  
  // Fishing the pieces ...
  for(UInt_t i=0;i<files.size();i++) {
        
    // Get species files:
    for(UInt_t j=0;j<NSpecies();j++) {
      string filename = between(files[i],simPath,".h5");
      if(filename.find(species[j]) != string::npos) {
	
	if( (filename.find("density") != string::npos) ) {
	  sCHG->at(j) = new string(files[i]);
	  continue;
	} else if(filename.find("pha") != string::npos) {
	  // Loop over Phase spaces
	  for(UInt_t ip=0;ip<NPhaseSpaces();ip++) {
	    if( filename.find(pspaces[ip].c_str()) != string::npos ) {
	      sPHA->at(j).at(ip) = new string(files[i]);
	      continue;
	    }
	  }
	}
      }
    }

    // RAW data: First extract number of RAW species
    if((files[i].find("raw") != string::npos) && (files[i].find(".h5") != string::npos)) {
      string s_beg = "raw_";
      string s_end = Form("_%s.h5",stime);
      
      string rawname = between(files[i],s_beg,s_end);
      rawspecies.push_back(rawname);
      // cout << Form(" Raw species name = %s",rawname.c_str()) << endl;
    }
    
    
    // Get Electromagnetic fields files:
    for(UInt_t j=0;j<3;j++) {
      char eName[16];
      if(j==0) {
	sprintf(eName,"Ez");
	if(files[i].find(eName) != string::npos) 
	  sEF->at(j) = new string(files[i]);
      } else if(j==1) {
	sprintf(eName,"Bx");
	if( (files[i].find(eName) != string::npos) && (files[i].find("Ey") == string::npos)) 
	  sBF->at(j) = new string(files[i]);
      } else if(j==2) {
	sprintf(eName,"By");
	if( (files[i].find(eName) != string::npos) && (files[i].find("Ex") == string::npos)) 
	  sBF->at(j) = new string(files[i]);
      }
      
      if(j>=2) continue;
      
      if(j==0)
	sprintf(eName,"ExmBy");
      else if(j==1)
	sprintf(eName,"EypBx");
      
      if(files[i].find(eName) != string::npos) {
	sWF->at(j) = new string(files[i]);
      }
    }
  }

  sRAW = new vector<string*>(NRawSpecies(),NULL);
  sTrack = new vector<string*>(NRawSpecies(),NULL);
  
  // Put the driver first
  for(UInt_t i=0;i<NRawSpecies();i++) {
    if(rawspecies[i].find("driver") != string::npos)
      if(i>0) {
	string temp = rawspecies[0];
	rawspecies[0] = rawspecies[i];
	rawspecies[i] = temp;
      }
  }
  
  // Fishing just RAW data
  for(UInt_t i=0;i<files.size();i++) {    
    // Get raw species files:
    for(UInt_t j=0;j<NRawSpecies();j++) {
      if(files[i].find(rawspecies[j]) != string::npos) {
	if((files[i].find("raw") != string::npos) && (files[i].find(".h5") != string::npos)) {
	  sRAW->at(j) = new string(files[i]);
	}	
      }
    }
  }

  ReadOutputSummary();
  ThreeD = kTRUE;
  HiP = kTRUE;

  Bool_t tfound = kFALSE;
  if(NSpecies()) {
    for(UInt_t j=0;j<NSpecies();j++) {
      if(GetChargeFileName(j)) {
	rtime = GetRealTimeFromFile(GetChargeFileName(j)->c_str());
	tfound = kTRUE;
    	continue;
      }
    }
  }
  
  if(NRawSpecies() && !tfound) {
    for(UInt_t j=0;j<NRawSpecies();j++) {
      if(GetRawFileName(j)) {
	rtime = GetRealTimeFromFile(GetRawFileName(j)->c_str());
    	continue;
      }
    }
  }
  
  // Defines the sub-range for the analysis.
  // Here at initialization, it is set to the whole simulation range.
  Double_t shiftx1 = Shift("comov");  
  
  if(pParam.x1Min == -999.) {
    XMINR[0] = GetXMin(0);
  } else {
    XMINR[0] = pParam.x1Min + shiftx1;
    cout << Form(" x1Min = %f",XMINR[0]) << endl;
  }
  if(pParam.x1Max == -999.) {
    XMAXR[0] = GetXMax(0);
  } else {
    XMAXR[0] = pParam.x1Max + shiftx1;
    cout << Form(" x1Max = %f",XMAXR[0]) << endl;
  }
  if(pParam.x2Min == -999.) {
    XMINR[1] = GetXMin(1);
  } else {
    XMINR[1] = pParam.x2Min;
  }
  if(pParam.x2Max == -999.) {
    XMAXR[1] = GetXMax(1);
  } else {
    XMAXR[1] = pParam.x2Max;
  }
  if(pParam.x3Min == -999.) {
    XMINR[2] = GetXMin(2);
  } else {
    XMINR[2] = pParam.x3Min;
  }
  if(pParam.x3Max == -999.) {
    XMAXR[2] = GetXMax(2);
  } else {
    XMAXR[2] = pParam.x3Max;
  }
  
  Init = kTRUE;  
}

//_______________________________________________________________________
void PDataHiP::PrintData(Option_t *option) {
  PData::PrintData(option);

  cout << "Data for transverse Wakefields: " << endl;
  for(UInt_t iwf=0;iwf<2;iwf++) {
    if(sWF) 
      if(sWF->at(iwf))
	cout << " - " << sWF->at(iwf)->c_str() << endl;
    
  }

  cout << endl;
  
}


//_______________________________________________________________________
void PDataHiP::ReadOutputSummary(const char * pfile)
{
  ifstream ifile;
  string ifilename;
  string word;


  // Allocating box boundaries
  NDIM = 3;
  NX   = new Int_t[NDIM];
  XMIN = new Double_t[NDIM];
  XMAX = new Double_t[NDIM];
  XMINR = new Double_t[NDIM];
  XMAXR = new Double_t[NDIM];
  Float_t TMIN = -999, TMAX = -999;
  Int_t NT = -999;
  
  if(strcmp(pfile,"") == 0) {
    ifilename  = simPath + "/DATA/output_summary.txt";
  }
  
  ifile.open(ifilename.c_str(),ios::in);
  if (ifile.is_open()) {
    cout << "PDataHiP:: Reading output summary: " << ifilename.c_str() << endl;

    UInt_t iline = 0;
    string line;
    while ( getline(ifile,line) )
      {
	istringstream iss(line);
	
	if(iline==0) {
	  iss >> XMIN[0];
	  iss >> XMAX[0];
	  iss >> NX[0];
	  cout << Form(" - x1 range: [%8.2f, %8.2f]  N = %5i",XMIN[0],XMAX[0],NX[0]) << endl;
	} else if(iline==1) {
	  iss >> XMIN[1];
	  iss >> XMAX[1];
	  iss >> NX[1];
	  cout << Form(" - x1 range: [%8.2f, %8.2f]  N = %5i",XMIN[1],XMAX[1],NX[1]) << endl;
	} else if(iline==2) {
	  iss >> XMIN[2];
	  iss >> XMAX[2];
	  iss >> NX[2];
	  cout << Form(" - x3 range: [%8.2f, %8.2f]  N = %5i",XMIN[2],XMAX[2],NX[2]) << endl;
	} else if(iline==3) {
	  iss >> TMIN;
	  iss >> TMAX;
	  iss >> NT;
	  cout << Form(" - t  range: [%8.2f, %8.2f]  N = %5i --> dt = %8.2f",TMIN,TMAX,NT,fabs(TMAX-TMIN)/NT) << endl;
	} 
	
	iline++;
	
      }
    
    ifile.close();
    
  } else {
    cout << "PDataHiP:: No output summary file: " << ifilename.c_str() << endl;
  }
  
  
}

//_______________________________________________________________________
Double_t PDataHiP::GetRealTimeFromFile(const char *filename) {
  
  // Open input HDF5 file
  H5File h5 = H5File(filename,H5F_ACC_RDONLY);

  //  cout << Form("Reading time from %s", filename) << endl;
  
  // Open main group for data reading
  Group *root = new Group(h5.openGroup("/"));
  
  // Get time info from its attributes
  // ---------------------------------
  Attribute *att = new Attribute(root->openAttribute("TIME"));
  //Attribute *att = new Attribute(root->openAttribute("T"));
  Float_t rtime;
  att->read(PredType::NATIVE_FLOAT,&rtime); 
  att->close();
  delete att;
  delete root;

  return (Double_t) rtime;
}

//_______________________________________________________________________
TH1F* PDataHiP::GetH1SliceZ3D(const char *filename,const char *datanameold, 
			   Int_t Firstx2Bin, Int_t Lastx2Bin, 
			   Int_t Firstx3Bin, Int_t Lastx3Bin, const char *options) {
  
  // Save memory and time with an specific range selection.

  // Check for valid HDF5 file 
  if (!H5File::isHdf5(filename))    
    return NULL;
  
  // Options
  TString opt = options;
  
  // Open input HDF5 file
  H5File h5 = H5File(filename,H5F_ACC_RDONLY);
  // Open main group for data reading
  Group *root = new Group(h5.openGroup("/"));

  // Read data
  // ----------
  H5std_string strreadbuf ("");
  StrType strdatatype(PredType::C_S1, 32); 
  Attribute *att = new Attribute(root->openAttribute("NAME"));
  att->read(strdatatype,strreadbuf); 
  att->close();
  delete att;
  const char *dataname = strreadbuf.c_str();

  DataSet  *dataSet = new DataSet(root->openDataSet(dataname));
  DataSpace dataSpace = dataSet->getSpace();
  const DataType type = dataSet->getDataType();
    
  Int_t rank = dataSpace.getSimpleExtentNdims();
  hsize_t  *dataDims = new hsize_t[rank];
  dataSpace.getSimpleExtentDims(dataDims,NULL);

  // x3 slice
  if(Lastx3Bin>=(Int_t)dataDims[2]) Lastx3Bin = dataDims[2] - 1;
  
  if(Firstx3Bin<0)
    DoSlice(dataDims[2],Firstx3Bin,Lastx3Bin);
  
  if(Lastx3Bin<Firstx3Bin) 
    Lastx3Bin = Firstx3Bin;
  
  UInt_t x3Dim = Lastx3Bin - Firstx3Bin + 1; 
  
  //----------

  // x2 slice
  if(Lastx2Bin>=(Int_t)dataDims[1]) Lastx2Bin = dataDims[1] - 1;
  
  if(Firstx2Bin<0) {
    DoSlice(dataDims[1],Firstx2Bin,Lastx2Bin);
  }
  
  if(Lastx2Bin<Firstx2Bin) 
    Lastx2Bin = Firstx2Bin;

  UInt_t x2Dim = Lastx2Bin - Firstx2Bin + 1; 
 
  //----------
  UInt_t x1AvgFactor = GetNX(0)/dataDims[0];
  UInt_t x1Dim = GetX1N()/x1AvgFactor;
  Double_t x1Min = GetX1Min();
  Double_t x1Max = GetX1Max();
  
  hsize_t  *count  = new hsize_t[rank];
  count[0] = x1Dim;
  count[1] = x2Dim;
  count[2] = x3Dim;

  hsize_t  *offset = new hsize_t[rank];
  offset[0] = GetX1iMin()/x1AvgFactor;
  offset[1] = Firstx2Bin;
  offset[2] = Firstx3Bin;
  
  Float_t *data = new Float_t[x1Dim * x2Dim * x3Dim];

  DataSpace memSpace(rank,count);
  dataSpace.selectHyperslab(H5S_SELECT_SET,count,offset);
  dataSet->read(data,PredType::NATIVE_FLOAT,memSpace,dataSpace);
  dataSet->close();
  
  root->close();

  Float_t x1Array[x1Dim];
  // Set all elements to 0.
  memset(x1Array,0,sizeof(Float_t)*x1Dim);
  
  string sdata = dataname;
  // Sum the values 
  for(UInt_t i=0;i<x1Dim;i++) {
    for(UInt_t j=0;j<x2Dim;j++) {
      for(UInt_t k=0;k<x3Dim;k++) {
	UInt_t index = (long)i*(long)x2Dim*(long)x3Dim + (long)j*(long)x3Dim + (long)k;
		
	if(sdata.find("charge") != string::npos || sdata.find("p1x1") != string::npos || sdata.find("p2x2") != string::npos)
	  x1Array[i] += -data[index];
	else
	  x1Array[i] +=  data[index]; 
	
      }
    }
  }
  
  // Histogram centering
  Double_t shiftx1 = Shift(opt);  

  x1Min -= shiftx1;
  x1Max -= shiftx1;

  Float_t x2AvgFactor = GetNX(1)/dataDims[1];
  Float_t dx2   = GetDX(1) * x2AvgFactor;

  Float_t x3AvgFactor = GetNX(2)/dataDims[2];
  Float_t dx3   = GetDX(2) * x3AvgFactor;
  
  TH1F *h1D = new TH1F(); 
  h1D->SetBins(x1Dim,x1Min,x1Max);
  for(UInt_t k=0;k<x1Dim;k++) { 
    Double_t content = x1Array[k];
    if(opt.Contains("avg")) content /= x3Dim * x2Dim;
    else if(opt.Contains("int")) content *= dx2 * dx3;
    h1D->SetBinContent(k+1,content);
  }
  root->close();
  
  delete count;
  delete offset;
  delete dataSet;
  delete root;
  delete data;
  
  return h1D;
}

//_______________________________________________________________________
TH2F* PDataHiP::GetH2SliceZX(const char *filename,const char *datanameold, Int_t Firstx3Bin, Int_t Lastx3Bin, const char *options) {
  
  // Check for valid HDF5 file 
  if (!H5File::isHdf5(filename))    
    return NULL;
  
  // Save memory and time with an specific range selection.

  // Options
  TString opt = options;

  // Open input HDF5 file
  H5File h5 = H5File(filename,H5F_ACC_RDONLY);
  
  // Open main group for data reading
  Group *root = new Group(h5.openGroup("/"));

  // Read data
  // ----------
  H5std_string strreadbuf ("");
  StrType strdatatype(PredType::C_S1, 32); 
  Attribute *att = new Attribute(root->openAttribute("NAME"));
  att->read(strdatatype,strreadbuf); 
  att->close();
  delete att;
  const char *dataname = strreadbuf.c_str();
  
  //  cout << Form(" Dataname = %s ", dataname) << endl;

  DataSet  *dataSet = new DataSet(root->openDataSet(dataname));
  DataSpace dataSpace = dataSet->getSpace();
  const DataType type = dataSet->getDataType();
  
  Int_t rank = dataSpace.getSimpleExtentNdims();
  hsize_t  *dataDims = new hsize_t[rank];
  dataSpace.getSimpleExtentDims(dataDims,NULL);

  // x3 slice
  if(Lastx3Bin>=(Int_t)dataDims[2]) Lastx3Bin = dataDims[2] - 1;
  
  if(Firstx3Bin<0)
    DoSlice(dataDims[2],Firstx3Bin,Lastx3Bin);
  
  if(Lastx3Bin<Firstx3Bin) 
    Lastx3Bin = Firstx3Bin;

  UInt_t x3Dim = Lastx3Bin - Firstx3Bin + 1; 
 
  //----------
  // Fix average problem.
  UInt_t x1AvgFactor = GetNX(0)/dataDims[0];
  UInt_t x1Dim = GetX1N()/x1AvgFactor;
  Double_t x1Min = GetX1Min();
  Double_t x1Max = GetX1Max();

  UInt_t x2AvgFactor = GetNX(1)/dataDims[1];
  UInt_t x2Dim = GetX2N()/x2AvgFactor;
  Double_t x2Min = GetX2Min();
  Double_t x2Max = GetX2Max();
   
  hsize_t  *count  = new hsize_t[rank];
  count[0] = x1Dim;
  count[1] = x2Dim;
  count[2] = x3Dim;

  hsize_t  *offset = new hsize_t[rank];  
  offset[0] = GetX1iMin()/x1AvgFactor;
  offset[1] = GetX2iMin()/x2AvgFactor;
  offset[2] = Firstx3Bin; 

  //  cout << Form(" x1Dim = %i  offset = %i", count[0], offset[0]) << endl;
  //  cout << Form(" x2Dim = %i  offset = %i", count[1], offset[1]) << endl;
  //  cout << Form(" x3Dim = %i  offset = %i", count[2], offset[2]) << endl;

  Float_t *data = new Float_t[x1Dim * x2Dim * x3Dim];

  DataSpace memSpace(rank,count);
  dataSpace.selectHyperslab(H5S_SELECT_SET,count,offset);
  dataSet->read(data,PredType::NATIVE_FLOAT,memSpace,dataSpace);
  dataSet->close();

  root->close();
  
  Float_t x1x2Array[x1Dim][x2Dim];
  memset(x1x2Array,0,sizeof(Float_t)*x1Dim*x2Dim);
   
  // Sum the values 
  string sdata = dataname;
  for(UInt_t i=0;i<x1Dim;i++) {
    for(UInt_t j=0;j<x2Dim;j++) {
      for(UInt_t k=0;k<x3Dim;k++) {
	UInt_t index = (long)i*(long)x2Dim*(long)x3Dim + (long)j*(long)x3Dim + (long)k;
	if(sdata.find("charge") != string::npos || sdata.find("p1x1") != string::npos || sdata.find("p2x2") != string::npos)
	  x1x2Array[i][j] += -data[index];
	else
	  x1x2Array[i][j] +=  data[index]; 
      }
    }
  }
  
  // Histogram centering
  Double_t shiftx1 = Shift(opt);  
  Double_t shiftx2 = 0;
  if(opt.Contains("center")) 
    if(!opt.Contains("cyl"))
      shiftx2 += (x2Min + x2Max)/2.;
  
  x1Min -= shiftx1;
  x1Max -= shiftx1;
  
  x2Min -= shiftx2;
  x2Max -= shiftx2;

  UInt_t x3AvgFactor = GetNX(2)/dataDims[2];
  Double_t dx3 = GetDX(2) * x3AvgFactor;

  TH2F *h2D = new TH2F(); 
  h2D->SetBins(x1Dim,x1Min,x1Max,x2Dim,x2Min,x2Max); 
  for(UInt_t i=0;i<x1Dim;i++) 
    for(UInt_t j=0;j<x2Dim;j++) { 
      Double_t content = x1x2Array[i][j];
      if(opt.Contains("avg")) content /= x3Dim;
      else if(opt.Contains("int")) content *= dx3;
      h2D->SetBinContent(i+1,j+1,content);
    }
    
  delete count;
  delete offset;
  delete dataSet;
  delete root;
  delete data;
  
  return h2D;
}

//_______________________________________________________________________
TH2F* PDataHiP::GetH2SliceZY(const char *filename,const char *datanameold, Int_t Firstx2Bin, Int_t Lastx2Bin, const char *options) {
  
  // Check for valid HDF5 file 
  if (!H5File::isHdf5(filename))    
    return NULL;
  
  // Save memory and time with an specific range selection.

  // Options
  TString opt = options;

  // Open input HDF5 file
  H5File h5 = H5File(filename,H5F_ACC_RDONLY);
  
  // Open main group for data reading
  Group *root = new Group(h5.openGroup("/"));

  // Read data
  // ----------
  H5std_string strreadbuf ("");
  StrType strdatatype(PredType::C_S1, 32); 
  Attribute *att = new Attribute(root->openAttribute("NAME"));
  att->read(strdatatype,strreadbuf); 
  att->close();
  delete att;
  const char *dataname = strreadbuf.c_str();
  
  //  cout << Form(" Dataname = %s ", dataname) << endl;

  DataSet  *dataSet = new DataSet(root->openDataSet(dataname));
  DataSpace dataSpace = dataSet->getSpace();
  const DataType type = dataSet->getDataType();
  
  Int_t rank = dataSpace.getSimpleExtentNdims();
  hsize_t  *dataDims = new hsize_t[rank];
  dataSpace.getSimpleExtentDims(dataDims,NULL);

  // x2 slice
  if(Lastx2Bin>=(Int_t)dataDims[1]) Lastx2Bin = dataDims[1] - 1;
  
  if(Firstx2Bin<0)
    DoSlice(dataDims[1],Firstx2Bin,Lastx2Bin);
  
  if(Lastx2Bin<Firstx2Bin) 
    Lastx2Bin = Firstx2Bin;

  UInt_t x2Dim = Lastx2Bin - Firstx2Bin + 1; 
 
  //----------
  // Fix average problem.
  UInt_t x1AvgFactor = GetNX(0)/dataDims[0];
  UInt_t x1Dim = GetX1N()/x1AvgFactor;
  Double_t x1Min = GetX1Min();
  Double_t x1Max = GetX1Max();

  UInt_t x3AvgFactor = GetNX(2)/dataDims[2];
  UInt_t x3Dim = GetX3N()/x3AvgFactor;
  Double_t x3Min = GetX3Min();
  Double_t x3Max = GetX3Max();
   
  hsize_t  *count  = new hsize_t[rank];
  count[0] = x1Dim;
  count[1] = x2Dim;
  count[2] = x3Dim;

  hsize_t  *offset = new hsize_t[rank];  
  offset[0] = GetX1iMin()/x1AvgFactor;
  offset[1] = Firstx2Bin;
  offset[2] = GetX3iMin()/x3AvgFactor; 

  //  cout << Form(" x1Dim = %i  offset = %i", count[0], offset[0]) << endl;
  //  cout << Form(" x2Dim = %i  offset = %i", count[1], offset[1]) << endl;
  //  cout << Form(" x3Dim = %i  offset = %i", count[2], offset[2]) << endl;

  Float_t *data = new Float_t[x1Dim * x2Dim * x3Dim];

  DataSpace memSpace(rank,count);
  dataSpace.selectHyperslab(H5S_SELECT_SET,count,offset);
  dataSet->read(data,PredType::NATIVE_FLOAT,memSpace,dataSpace);
  dataSet->close();

  root->close();
  
  Float_t x1x3Array[x1Dim][x3Dim];
  memset(x1x3Array,0,sizeof(Float_t)*x1Dim*x3Dim);
   
  // Sum the values 
  string sdata = dataname;
  for(UInt_t i=0;i<x1Dim;i++) {
    for(UInt_t j=0;j<x2Dim;j++) {
      for(UInt_t k=0;k<x3Dim;k++) {
	UInt_t index = (long)i*(long)x2Dim*(long)x3Dim + (long)j*(long)x3Dim + (long)k;
	if(sdata.find("charge") != string::npos || sdata.find("p1x1") != string::npos || sdata.find("p2x2") != string::npos)
	  x1x3Array[i][k] += -data[index];
	else
	  x1x3Array[i][k] +=  data[index]; 
      }
    }
  }
  
  // Histogram centering
  Double_t shiftx1 = Shift(opt);  
  Double_t shiftx3 = 0;
  if(opt.Contains("center")) 
    if(!opt.Contains("cyl"))
      shiftx3 += (x3Min + x3Max)/2.;
  
  x1Min -= shiftx1;
  x1Max -= shiftx1;
  
  x3Min -= shiftx3;
  x3Max -= shiftx3;

  UInt_t x2AvgFactor = GetNX(1)/dataDims[1];
  Double_t dx2 = GetDX(1) * x2AvgFactor;

  TH2F *h2D = new TH2F(); 
  h2D->SetBins(x1Dim,x1Min,x1Max,x3Dim,x3Min,x3Max); 
  for(UInt_t i=0;i<x1Dim;i++) 
    for(UInt_t j=0;j<x3Dim;j++) { 
      Double_t content = x1x3Array[i][j];
      if(opt.Contains("avg")) content /= x2Dim;
      else if(opt.Contains("int")) content *= dx2;
      h2D->SetBinContent(i+1,j+1,content);
    }
    
  delete count;
  delete offset;
  delete dataSet;
  delete root;
  delete data;
  
  return h2D;
}

//_______________________________________________________________________
string*  PDataHiP::GetWfieldFileName(UInt_t i) { return sWF->at(i); }

//_______________________________________________________________________
UInt_t  PDataHiP::NRawSpecies() { return rawspecies.size(); }

//_______________________________________________________________________
string  PDataHiP::GetRawSpeciesName(UInt_t i) { return rawspecies.at(i); }

//______________________________________________________________________________________
Double_t PDataHiP::Shift(TString option) {
  TString opt = option;

  Double_t shiftx1 = 0;
  if(opt.Contains("center")) {
    Double_t kp = GetPlasmaK();
    shiftx1 += GetXMin(0) + GetBeamStart()*kp;   // Centers on the beam
  }
  
  return shiftx1;
}

//______________________________________________________________________________________
Double_t PDataHiP::ShiftT(TString option) {
  TString opt = option;
  
  Double_t shiftt = 0;
  Double_t kp = GetPlasmaK();
  if(opt.Contains("center")) {
    shiftt += GetXMin(0) + GetBeamStart()*kp - GetPlasmaStart()*kp;
  }
  
  return shiftt;
}
