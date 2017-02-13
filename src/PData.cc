#include "PData.hh"
#include "hdf2root.hh"

#include <H5Cpp.h>

ClassImp(PData)

using namespace H5; 

PData * PData::fgData = NULL;
PData * gData = NULL;

//_______________________________________________________________________
PData::PData(const char * name, const char * title) : TNamed(name,title) {
  time = 0;
  rtime = 0;
  simPath = name;

  species.clear(); 
  pspaces.clear();

  sCHG = sEF = sBF = sRAW = sTrack = NULL;
  sJ[0] = sJ[1] = sJ[2] = NULL;
  sPHA = NULL;
  NX = NULL; XMIN = XMAX = NULL;
  XMINR = XMAXR = NULL;
  NDIM = 0;

  ResetParameters();
  ReadParameters();

  Init = kFALSE;
  Cyl = kFALSE;
  ThreeD = kFALSE;
  Osi = kFALSE;
  HiP = kFALSE;
  Navg = 1;

  TString NAME = name;
  if(NAME.Contains("cyl")) Cyl = kTRUE;

  gData = this;
}

//_______________________________________________________________________
PData::PData(const char * name) : TNamed(name,name) {
  time = 0;
  rtime = 0;
  simPath = name;

  species.clear(); 
  pspaces.clear();

  sCHG = sEF = sBF = sRAW = sTrack = NULL;
  sJ[0] = sJ[1] = sJ[2] = NULL;
  sPHA = NULL;
  NX = NULL; XMIN = XMAX = NULL;
  XMINR = XMAXR = NULL;
  NDIM = 0;
  
  ResetParameters();
  ReadParameters();

  Init = kFALSE;
  Cyl = kFALSE;
  ThreeD = kFALSE;
  Osi = kFALSE;
  HiP = kFALSE;
  Navg = 1;

  TString NAME = name;
  if(NAME.Contains("cyl")) Cyl = kTRUE;
 
  gData = this;
}

//_______________________________________________________________________
PData::PData(const char* name, UInt_t t) : TNamed(name,name), time(t)  {
  time = 0;
  rtime = 0;
  simPath = name;

  species.clear(); 
  pspaces.clear();
  
  sCHG = sEF = sBF = sRAW = sTrack = NULL;
  sJ[0] = sJ[1] = sJ[2] = NULL;
  sPHA = NULL;
  NX = NULL; XMIN = XMAX = NULL;
  XMINR = XMAXR = NULL;
  NDIM = 0;

  ResetParameters();
  ReadParameters();

  gData = this;
  
  Init = kFALSE;
  Cyl = kFALSE;
  ThreeD = kFALSE;
  Osi = kFALSE;
  HiP = kFALSE;
  Navg = 1;

  TString NAME = name;
  if(NAME.Contains("cyl")) Cyl = kTRUE;
 
  LoadFileNames(time);
}

//_______________________________________________________________________
PData * PData::Get(const char * name)
{
  if (!fgData) {
    if (name) {
      fgData = new PData(name);
    } else {
      fgData = new PData("");
    }
  }

  if (name && strcmp(fgData->GetName(), name)!=0) {
    cout << "PData: Changing analysis to " << name << endl;

    if(fgData) delete fgData;
    fgData = new PData(name);
    // fgData->SetName(name);
    // fgData->SetTitle(name);
    // TString path = "./";
    // path += name;
    // fgData->SetPath(path.Data());
  }
  
  return fgData;
}

//_______________________________________________________________________
void PData::ReadParameters(const char * pfile)
{
  ifstream ifile;
  string ifilename;
  string word;

  if(strcmp(pfile,"") == 0) {
    ifilename  = simPath + "/";
    ifilename += GetName();
    ifilename += ".in";
  }
  
  ifile.open(ifilename.c_str(),ios::in);
  if (ifile.is_open()) {
    cout << "\n PData:: Reading parameters file: " << ifilename.c_str() << endl;

    string line;
    while ( getline(ifile,line) )
      {
	istringstream iss(line);
	
	iss >> word;
	if(word.find("!") != string::npos) continue;
	if(word.empty()) continue;
	
	if(word.find("pDensity") != string::npos) {  // Plasma density is read here in part/cm^3 .
	  iss >> pParam.pDensity;
	  cout << endl;
	  cout << Form(" - Plasma density   = %8.4e cm3",GetPlasmaDensity() / (1./PUnits::cm3)) << endl;
	  cout << Form(" - Plasma frequency = %8.4f fs^-1", GetPlasmaFrequency() * PUnits::femtosecond) << endl;
	  cout << Form(" - Time depth       = %8.4f fs", GetPlasmaTimeDepth() / PUnits::femtosecond) << endl;
	  cout << Form(" - Wavenumber       = %8.4f um^-1", GetPlasmaK() * PUnits::um) << endl;
	  cout << Form(" - Skin depth       = %8.4f um", GetPlasmaSkinDepth() / PUnits::um) << endl;
	  cout << Form(" - E_0              = %8.4f GV/m", GetPlasmaE0() * PUnits::meter / PUnits::GV ) << endl;
	}
	else if(word.find("pStart") != string::npos) {
	  iss >> pParam.pStart;
	  cout << Form(" - Plasma start     = %8.4f mm", GetPlasmaStart() / PUnits::mm ) << endl;
	}
	else if(word.find("nDensity") != string::npos) {  // Neutral density is read here in part/cm^3 .
	  iss >> pParam.nDensity;
	  cout << endl;
	  cout << Form(" - Neutral density  = %8.4e cm3",GetNeutralDensity() / (1./PUnits::cm3)) << endl;
	}
	else if(word.find("nStart") != string::npos) {
	  iss >> pParam.nStart;
	  cout << Form(" - Neutral start    = %8.4f mm", GetNeutralStart() / PUnits::mm ) << endl;
	}
	else if(word.find("nEnd") != string::npos) {
	  iss >> pParam.nEnd;
	  cout << Form(" - Neutral end      = %8.4f mm", GetNeutralEnd() / PUnits::mm ) << endl;
	}
	else if(word.find("bMass") != string::npos) {
	  iss >> pParam.bMass;
	  cout << endl;
	  cout << Form(" - Beam mass        = %8.4f MeV/c2", GetBeamMass() / PUnits::MeV) << endl;
	} 
	else if(word.find("bDensity") != string::npos) {
	  iss >> pParam.bDensity;
	  cout << Form(" - Bunch density    = %8.4e cm3", GetBeamDensity() / (1./PUnits::cm3)) << endl;
	}
	else if(word.find("bGamma") != string::npos) {
	  iss >> pParam.bGamma;
	  cout << Form(" - Bunch gamma      = %8.4f", GetBeamGamma() ) << endl;
	  cout << Form(" - Bunch energy     = %8.4f MeV", GetBeamEnergy() / PUnits::MeV) << endl;
	  cout << Form(" - Bunch velocity   = %8.4f c",  GetBeamVelocity()) << endl;
	}  
	else if(word.find("bStart") != string::npos) {
	  iss >> pParam.bStart;
	  cout << Form(" - Beam start       = %8.4f mm", GetBeamStart() / PUnits::mm ) << endl;
	}
	else if(word.find("bRmsZ") != string::npos)
	  iss >> pParam.bRmsZ;
	else if(word.find("bRmsX") != string::npos) {
	  iss >> pParam.bRmsX;
	  cout << Form(" - Beam X RMS       = %8.4f um", GetBeamRmsX() / PUnits::um) << endl;
	} else if(word.find("bRmsY") != string::npos) {
	  iss >> pParam.bRmsY;
	  cout << Form(" - Beam Y RMS       = %8.4f um", GetBeamRmsY() / PUnits::um) << endl;
	} else if(word.find("bRmsR") != string::npos)
	  iss >> pParam.bRmsR;
	else if(word.find("lOmega") != string::npos) {
	  iss >> pParam.lOmega;
	  cout << Form(" - Laser Omega      = %8.4f fs^-1", GetLaserOmega() * PUnits::fs) << endl;
	} else if(word.find("x1Min") != string::npos) {
	  iss >> pParam.x1Min; 
	  cout << Form(" - x1Min            = %8.4f", pParam.x1Min) << endl;
	} else if(word.find("x1Max") != string::npos) {
	  iss >> pParam.x1Max;
	  cout << Form(" - x1Max            = %8.4f", pParam.x1Max) << endl;
	} else if(word.find("x2Min") != string::npos)
	  iss >> pParam.x2Min;
	else if(word.find("x2Max") != string::npos)
	  iss >> pParam.x2Max;
	else if(word.find("x3Min") != string::npos)
	  iss >> pParam.x3Min;
	else if(word.find("x3Max") != string::npos)
	  iss >> pParam.x3Max;
	else if(word.find("EMin") != string::npos)
	  iss >> pParam.EMin;
	else if(word.find("EMax") != string::npos)
	  iss >> pParam.EMax;
	else if(word.find("denMin0") != string::npos)
	  iss >> pParam.denMin;
	else if(word.find("denMax0") != string::npos)
	  iss >> pParam.denMax;
	else if(word.find("denLoc0") != string::npos)
	  iss >> pParam.denLoc;
	else if(word.find("denMin1") != string::npos)
	  iss >> pParam.denMin1;
	else if(word.find("denMax1") != string::npos)
	  iss >> pParam.denMax1;
	else if(word.find("denMin2") != string::npos)
	  iss >> pParam.denMin2;
	else if(word.find("denMax2") != string::npos)
	  iss >> pParam.denMax2;
	else if(word.find("denMin3") != string::npos)
	  iss >> pParam.denMin3;
	else if(word.find("denMax3") != string::npos)
	  iss >> pParam.denMax3;

	word.clear();
      }

    cout << endl;
    ifile.close();
    
  } else {
    cout << "PData:: No input parameters file: " << ifilename.c_str() << endl;
  }
  
  
}


//_______________________________________________________________________
void PData::LoadFileNames(Int_t t) {
  
  if(time == t && Init) {
    // cout << "Time step already loaded: doing nothing." << endl;
    cout << " ... " << endl << endl;
    return;  
  }
  
  Clear();
  
  time = t;
  
  string Dir;
  
  // Get the list of different species
  vector<string> sptemp;
  Dir = simPath + "/MS/DENSITY";
  ListDir(Dir,string(""),sptemp,"nam");

  // Get the list of RAW species
  Dir = simPath + "/MS/RAW";
  ListDir(Dir,string(""),sptemp,"nam");

  // Eliminate redundancies
  for(UInt_t i=0;i<sptemp.size();i++) {
    for(UInt_t j=0;j<i;j++) {
      if(sptemp[j] == sptemp[i]) {
	sptemp.erase(sptemp.begin()+i);
	i--;
	break;
      }
    }
  }
  
  // It is very important the order of the species:
  // The analisys macros always assume that the plasma specie is the first one of the list.
  // The electron bunch, for instance, uses to be the second.
  // We swap here the order (if necessary) to put the "plasma" species at first.
  for(UInt_t i=0;i<sptemp.size();i++) {
    if(sptemp[i].find("plasma") != string::npos)
      if(i!=0 && sptemp.size()>0) {
	string temp = sptemp[0];
	sptemp[0] = sptemp[i];
	sptemp[i] = temp;
	i--;
	continue;
      }
    if(sptemp[i].find("driver") != string::npos)
      if(i!=1 && sptemp.size()>1) {
	string temp = sptemp[1];
	sptemp[1] = sptemp[i];
	sptemp[i] = temp;
	i--;
	continue;
      }
    if( (sptemp[i].find("high") != string::npos) ||
	(sptemp[i].find("1") != string::npos))
      if(i!=2 && sptemp.size()>2) {
	string temp = sptemp[2];
	sptemp[2] = sptemp[i];
	sptemp[i] = temp;
	i--;
	continue;
      }
  }

  for(UInt_t i=0;i<sptemp.size();i++) {
    species.push_back(sptemp[i]);
    //cout << Form(" Species = %s", species[i].c_str()) << endl;
  }
  
  // cout << Form("\n Number species = %i",species.size()) << endl;
  // for(UInt_t i=0;i<species.size();i++) {
  //   cout << Form(" Species = %s", species[i].c_str()) << endl;
  // }
  
  // Get the list of different phase spaces
  Dir = simPath + "/MS/PHA";
  ListDir(Dir,string(""),pspaces,"nam");
    
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
  Dir = simPath + "/MS";
  ListDir(Dir,string(stime),files,"recursive");
  
  // File containers for the different species
  // For every kind we need the following:
  // 1. Charge density.
  // 2. RAW data (discrete variables).
  // 3. Phase space files (not implemented yet).

  if(files.size()==0) {
    cout << "PData:: No files found in time step ( " << time << " ) : doing nothing!" << endl;
    Init = kFALSE;
    return;
  }
  
  // Initialize the pointers to NULLS
  sCHG = new vector<string*>(NSpecies(),NULL);
  sPHA = new vector<vector<string*> >(NSpecies(),vector<string*>(NPhaseSpaces(),NULL));
  sEF  = new vector<string*>(3,NULL);
  sBF  = new vector<string*>(3,NULL);
  sRAW = new vector<string*>(NRawSpecies(),NULL);
  sTrack = new vector<string*>(NRawSpecies(),NULL);

  for(UInt_t i=0;i<3;i++) {
    sJ[i] = new vector<string*>(NSpecies(),NULL);    
  }
  
  // Fishing the pieces ...
  for(UInt_t i=0;i<files.size();i++) {

    if(files[i].find("slice") != string::npos)
      continue;
    
    // Get species files:
    for(UInt_t j=0;j<NSpecies();j++) {

      if(files[i].find(species[j]) != string::npos) {

	if( (files[i].find("DENSITY") != string::npos) && ( files[i].find("charge") != string::npos)) {
	  if(files[i].find("savg") == string::npos)
	    sCHG->at(j) = new string(files[i]);
	  else if(!sCHG->at(j))
	    sCHG->at(j) = new string(files[i]);
	  
	} else if( (files[i].find("DENSITY") != string::npos) && ( files[i].find("j1") != string::npos)) {
	  if(files[i].find("savg") == string::npos)
	    sJ[0]->at(j) = new string(files[i]);
	  else if(!sCHG->at(j))
	    sJ[0]->at(j) = new string(files[i]);
	  
	} else if( (files[i].find("DENSITY") != string::npos) && ( files[i].find("j2") != string::npos)) {
	  if(files[i].find("savg") == string::npos)
	    sJ[1]->at(j) = new string(files[i]);
	  else if(!sCHG->at(j))
	    sJ[1]->at(j) = new string(files[i]);
	  
	} else if( (files[i].find("DENSITY") != string::npos) && ( files[i].find("j3") != string::npos)) {
	  if(files[i].find("savg") == string::npos)
	    sJ[2]->at(j) = new string(files[i]);
	  else if(!sCHG->at(j))
	    sJ[2]->at(j) = new string(files[i]);
	  
	} else if(files[i].find("PHA") != string::npos) {
	  // Loop over Phase spaces
	  for(UInt_t ip=0;ip<NPhaseSpaces();ip++) {
	    if( files[i].find(pspaces[ip].c_str()) != string::npos ) {
	      sPHA->at(j).at(ip) = new string(files[i]);
	      continue;
	    }
	  }
	}
      }
    }

    // Get RAW species files:
    for(UInt_t j=0;j<NRawSpecies();j++) {

      if(files[i].find(species[j]) != string::npos) {
	
	if((files[i].find("RAW") != string::npos) && (files[i].find(".h5") != string::npos)) {
	  sRAW->at(j) = new string(files[i]);
	}

	continue;
      }
    }
        
    // Get Electromagnetic fields files:
    for(UInt_t j=0;j<3;j++) {
      char eName[16];

      sprintf(eName,"e%1i-",j+1);   
      if(files[i].find(eName) != string::npos) {
	if(files[i].find("savg") == string::npos)
	  sEF->at(j) = new string(files[i]);
	else if(!sEF->at(j))
	  sEF->at(j) = new string(files[i]);
	
	continue;
      }

      sprintf(eName,"b%1i-",j+1);
      if(files[i].find(eName) != string::npos) {
	if(files[i].find("savg") == string::npos)
	  sBF->at(j) = new string(files[i]);
	else if(!sBF->at(j))
	  sBF->at(j) = new string(files[i]);
	
	continue;
      }
    }
  
  }
  
  // Get the vector of particle tracking files.
  vector<string> tfiles;
  Dir = simPath + "/MS/TRACKS";
  ListDir(Dir,".h5",tfiles,"recursive");
  
  // Fishing the pieces ...
  for(UInt_t i=0;i<tfiles.size();i++) {

    // Get RAW species files:
    for(UInt_t j=0;j<NRawSpecies();j++) {
      
      if(tfiles[i].find(species[j]) != string::npos) {
	if((tfiles[i].find("tracks") != string::npos) && (tfiles[i].find(".h5") != string::npos)) {
	  sTrack->at(j) = new string(tfiles[i]);
	}
	
	continue;
      }
    }
  }
  
  
  if(species.size()) {
    if(sCHG->at(0)) {
      rtime = GetRealTimeFromFile(GetChargeFileName(0)->c_str());
      GetBoxDimensionsFromFile(GetChargeFileName(0)->c_str());
    } else if(sRAW->at(0)) {
      rtime = GetRealTimeFromFile(GetRawFileName(0)->c_str());
      GetBoxDimensionsFromFile(GetRawFileName(0)->c_str());
    }
  } else {
    // cout << "PData:: No species folders in this simulation: " << GetPath() << "." << endl;
    // cout << "Checking the fields..." << endl;
    
    for(UInt_t i=0;i<3;i++) { 
      if(sEF->at(i)) {
	rtime = GetRealTimeFromFile(GetEfieldFileName(i)->c_str());
	GetBoxDimensionsFromFile(GetEfieldFileName(i)->c_str());
	break;
      } else
	continue;
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
void PData::PrintData(Option_t *option) {

  cout << endl <<  "PData:: Print available data in simulation:  " << GetName() << endl;
  
  if(NSpecies()==0) {
    cout << "No species loaded!!" << endl;
  } else {
    cout << "Species to analize (" << NSpecies() << "):  ";
    for(UInt_t i=0;i<NSpecies();i++) {
      cout << species[i];
      if(i<NSpecies()-1) cout << ", " ;
      else cout << endl;
    }
  }

  if(NPhaseSpaces()==0) {
    cout << "No phasespaces loaded!!" << endl;
  }  else {
    cout << "Phase spaces to analize (" << NPhaseSpaces() << "):  ";
    for(UInt_t i=0;i<NPhaseSpaces();i++) {
      cout << pspaces[i];
      if(i<NPhaseSpaces()-1) cout << ", " ;
      else cout << endl;
    }
  }
  cout << endl;

  // Print selected files list to stdout 
  for(UInt_t is=0;is<NSpecies();is++) {
    cout << "Data for '" << species[is] << "' : " << endl;

    if(sCHG->at(is)) 
      cout << " - " << sCHG->at(is)->c_str() << endl;
    
    for(UInt_t i=0;i<3;i++) {
      if(!sJ[i]) continue;
      if(sJ[i]->at(is)) 
	cout << " - " << sJ[i]->at(is)->c_str() << endl;
    }
    
    for(UInt_t ip=0;ip<NPhaseSpaces();ip++) {
      if(sPHA->at(is).at(ip))
	cout << " - " << sPHA->at(is).at(ip)->c_str() << endl;
    }

    cout << endl;
  }
  
  cout << "Data for RAW macroparticles:" << endl;
  for(UInt_t is=0;is<NRawSpecies();is++) {
    if(sRAW->at(is)) 
      cout << " - " << sRAW->at(is)->c_str() << endl;
  }
  cout << endl;

  cout << "Data for TRACKING macroparticles:" << endl;
  for(UInt_t is=0;is<NRawSpecies();is++) {
    if(sTrack->at(is)) 
      cout << " - " << sTrack->at(is)->c_str() << endl;
  }
  cout << endl;
  
  cout << "Data for Electromagnetic fields: " << endl;
  for(UInt_t ief=0;ief<3;ief++) {
    if(sEF) 
      if(sEF->at(ief))
	cout << " - " << sEF->at(ief)->c_str() << endl;
    if(sBF)
      if(sBF->at(ief))
	cout << " - " << sBF->at(ief)->c_str() << endl;
  }
  cout << endl;
  
}

//_______________________________________________________________________
void PData::CopyData(const char *opath, const char *cpcmd) {
 
  string command;
  string odir;
   
  // Print selected files list to stdout 
  for(UInt_t is=0;is<NSpecies();is++) {
    cout << "Copying data for '" << species[is] << "' : " << endl;
    if(sCHG->at(is)) {
      odir = *sCHG->at(is);
      Int_t pos = odir.find("/");
      odir.replace(0,pos,opath);
      pos = odir.find_last_of("/");
      odir.erase(pos,string::npos);

      command = Form("%s %s","mkdir -p",odir.c_str());
      system(command.c_str());
      
      command = Form("%s %s %s",cpcmd,sCHG->at(is)->c_str(),odir.c_str());
      system(command.c_str());
    }

    if(sRAW->at(is)) {
      odir = *sRAW->at(is);
      Int_t pos = odir.find("/");
      odir.replace(0,pos,opath);
      pos = odir.find_last_of("/");
      odir.erase(pos,string::npos);

      command = Form("%s %s","mkdir -p",odir.c_str());
      system(command.c_str());
      
      command = Form("%s %s %s",cpcmd,sRAW->at(is)->c_str(),odir.c_str());
      system(command.c_str());
    }

    
    for(UInt_t ip=0;ip<NPhaseSpaces();ip++) {
      if(sPHA->at(is).at(ip)) {
	odir = *sPHA->at(is).at(ip);
	Int_t pos = odir.find("/");
	odir.replace(0,pos,opath);
	pos = odir.find_last_of("/");
	odir.erase(pos,string::npos);
	
	command = Form("%s %s","mkdir -p",odir.c_str());
	system(command.c_str());
	
	command = Form("%s %s %s",cpcmd,sPHA->at(is).at(ip)->c_str(),odir.c_str());
	system(command.c_str());
	
      }
    }
    
    cout << endl;
  }
 
  cout << "Copying data for Electromagnetic fields: " << endl;
  for(UInt_t ief=0;ief<3;ief++) {
    if(sEF) 
      if(sEF->at(ief)) {
	odir = *sEF->at(ief);
	Int_t pos = odir.find("/");
	odir.replace(0,pos,opath);
	pos = odir.find_last_of("/");
	odir.erase(pos,string::npos);
	
	command = Form("%s %s","mkdir -p",odir.c_str());
	system(command.c_str());
	
	command = Form("%s %s %s",cpcmd,sEF->at(ief)->c_str(),odir.c_str());
	system(command.c_str());	
      }
    if(sBF)
      if(sBF->at(ief)) {
	odir = *sBF->at(ief);
	Int_t pos = odir.find("/");
	odir.replace(0,pos,opath);
	pos = odir.find_last_of("/");
	odir.erase(pos,string::npos);
	
	command = Form("%s %s","mkdir -p",odir.c_str());
	system(command.c_str());
	
	command = Form("%s %s %s",cpcmd,sBF->at(ief)->c_str(),odir.c_str());
	system(command.c_str());
      }
  }
  cout << endl;
  
}

//_______________________________________________________________________
void PData::Delete(const char *dlcmd) {
 
  string command;
  string odir;
   
  // Print selected files list to stdout 
  for(UInt_t is=0;is<NSpecies();is++) {
    cout << "Deleting data for '" << species[is] << "' : " << endl;
    if(sCHG->at(is)) {
      command = Form("%s %s",dlcmd,sCHG->at(is)->c_str());
      cout << command.c_str() << endl;
      system(command.c_str());
    }

    if(sRAW->at(is)) {
      command = Form("%s %s",dlcmd,sRAW->at(is)->c_str());
      cout << command.c_str() << endl;
      system(command.c_str());
    }

    
    for(UInt_t ip=0;ip<NPhaseSpaces();ip++) {
      if(sPHA->at(is).at(ip)) {
	command = Form("%s %s",dlcmd,sPHA->at(is).at(ip)->c_str());
      cout << command.c_str() << endl;
	system(command.c_str());
      }
    }
    
    cout << endl;
  }
 
  cout << "Deleting data for Electromagnetic fields: " << endl;
  for(UInt_t ief=0;ief<3;ief++) {
    if(sEF) 
      if(sEF->at(ief)) {
	command = Form("%s %s",dlcmd,sEF->at(ief)->c_str());
	cout << command.c_str() << endl;
	system(command.c_str());	
      }
    if(sBF)
      if(sBF->at(ief)) {
	command = Form("%s %s",dlcmd,sBF->at(ief)->c_str());
	cout << command.c_str() << endl;
	system(command.c_str());
      }
  }
  cout << endl;
  
}

//_______________________________________________________________________
void PData::Clear(Option_t *option)
{
  time = 0;

  if(sCHG) {
    FreeClear(*sCHG);
    delete sCHG;
    sCHG = NULL;
  }

  for(UInt_t i=0;i<3;i++) {
    if(sJ[i]) {
      FreeClear(*sJ[i]);
      delete sJ[i];
      sJ[i] = NULL;
    }
  }
  
  if(sEF) {
    FreeClear(*sEF);
    delete sEF;
    sEF = NULL;
  }
  
  if(sBF) {
    FreeClear(*sBF);
    delete sBF;
    sBF = NULL;
  }
  
  if(sRAW) {
    FreeClear(*sRAW);
    delete sRAW;
    sRAW = NULL;
  }

  if(sTrack) {
    FreeClear(*sTrack);
    delete sTrack;
    sTrack = NULL;
  }

  if(sPHA) {
    for(UInt_t i=0;i<NSpecies();i++) {
      FreeClear(sPHA->at(i));
    }
    delete sPHA;  
    sPHA = NULL;
  }
  
  species.clear(); 
  pspaces.clear();  
  
}

//_______________________________________________________________________
PData::~PData()
{
  Clear();

  fgData = 0;
  gData  = 0;
}


//_______________________________________________________________________
void PData::GetBoxDimensionsFromFile(const char *filename) {
  
  // Open input HDF5 file
  H5File h5 = H5File(filename,H5F_ACC_RDONLY);
  
  // Open main group for data reading
  Group *root = new Group(h5.openGroup("/"));
  
  // Get box info from files attributes
  // ---------------------------------

  Attribute *att = new Attribute(root->openAttribute("NX"));
  DataSpace attSpace = att->getSpace();
  Int_t rank = attSpace.getSimpleExtentNdims();
  hsize_t *attDim = new hsize_t[rank];
  attSpace.getSimpleExtentDims(attDim,NULL);
  NDIM = attDim[0];

  NX   = new Int_t[NDIM];
  XMIN = new Double_t[NDIM];
  XMAX = new Double_t[NDIM];
  XMINR = new Double_t[NDIM];
  XMAXR = new Double_t[NDIM];
  
  att->read(PredType::NATIVE_INT,NX); 
  att->close();  delete att;
  
  att = new Attribute(root->openAttribute("XMIN"));
  att->read(PredType::NATIVE_DOUBLE,XMIN); 
  att->close();  delete att;
  
  att = new Attribute(root->openAttribute("XMAX"));
  att->read(PredType::NATIVE_DOUBLE,XMAX); 
  att->close();  delete att;  

  // Correct time
  // rtime += XMIN[0]-rtime;
  // Correct x1
  // Double_t X1RANGE =  XMAX[0] - XMIN[0];
  // XMIN[0] += rtime-XMIN[0];
  // XMAX[0] = XMIN[0] + X1RANGE;
  // cout << Form(" Time = %.10e",rtime) << endl;
  
  
  // for(Int_t i=0;i<NDIM;i++) 
  //   cout << Form(" Dimension %i:  N = %i  XMIN = %.10e  XMIN = %.10e -> DX = %.10e",i,NX[i],XMIN[i],XMAX[i], GetDX(i)) << endl; 

  // cout << Form(" x1-t = %.5e",XMIN[0]-rtime) << endl;
  
  if(NDIM==3) ThreeD = kTRUE;
  
  delete root;
}

//_______________________________________________________________________
Double_t PData::GetRealTimeFromFile(const char *filename) {
  
  // Open input HDF5 file
  H5File h5 = H5File(filename,H5F_ACC_RDONLY);
  
  // Open main group for data reading
  Group *root = new Group(h5.openGroup("/"));
  
  // Get time info from its attributes
  // ---------------------------------
  Attribute *att = new Attribute(root->openAttribute("TIME"));
  Double_t rtime;
  att->read(PredType::NATIVE_DOUBLE,&rtime); 

  // Write time info to a string
  // char stime[48];
  // sprintf(stime," at t=%6.5f [%s]",time,"#omega_{p}^{-1}"); 
  // cout << Form("\n Time = %.10e",rtime) << endl;

  att->close();
  delete att;
  delete root;

  return rtime;
}

//_______________________________________________________________________
TH1F* PData::GetH1SliceZ(const char *filename,const char *dataname,Int_t Firstx2Bin, Int_t Lastx2Bin, const char *options) {

  // Options
  TString opt = options;

  // Open input HDF5 file
  H5File h5 = H5File(filename,H5F_ACC_RDONLY);

  // Open main group for data reading
  Group *root = new Group(h5.openGroup("/"));

  // Read data
  // ----------
  DataSet    *dataSet = new DataSet(root->openDataSet(dataname));
  DataSpace dataSpace = dataSet->getSpace();
  const DataType type = dataSet->getDataType();

  UInt_t       rank = dataSpace.getSimpleExtentNdims();
  hsize_t *dataDims = new hsize_t[rank];
  dataSpace.getSimpleExtentDims(dataDims,NULL);

  // x2 slice 
  if(Lastx2Bin>=(Int_t)dataDims[0]) Lastx2Bin = dataDims[0] - 1;
  if(Firstx2Bin<0) {
    DoSlice(dataDims[0],Firstx2Bin,Lastx2Bin);
  }
  if(Lastx2Bin<Firstx2Bin) 
    Lastx2Bin = Firstx2Bin;
  
  UInt_t x2DimSl = Lastx2Bin - Firstx2Bin + 1; 

  UInt_t x1AvgFactor = GetNX(0)/dataDims[1];
  UInt_t x1Dim = GetX1N()/x1AvgFactor;
  Double_t x1Min = GetX1Min();
  Double_t x1Max = GetX1Max();

  UInt_t x2AvgFactor = GetNX(1)/dataDims[0];
  UInt_t x2Dim = GetX2N()/x2AvgFactor;
  Double_t x2Min = GetX2Min();
  Double_t x2Max = GetX2Max();

  hsize_t  *count  = new hsize_t[rank];
  count[0] = x2DimSl;
  count[1] = x1Dim;

  hsize_t  *offset = new hsize_t[rank];
  offset[0] = Firstx2Bin;
  offset[1] = GetX1iMin()/x1AvgFactor;
  
  Float_t *data = new Float_t[x2DimSl*x1Dim];   
  memset(data,0,sizeof(Float_t)*x2DimSl*x1Dim);
 
  DataSpace memSpace(rank,count);
  dataSpace.selectHyperslab(H5S_SELECT_SET,count,offset);
  dataSet->read(data,PredType::NATIVE_FLOAT,memSpace,dataSpace);
  dataSet->close();
  
  root->close();

  Float_t *x1Array = new Float_t[x1Dim];
  memset(x1Array,0,sizeof(Float_t)*x1Dim);

  // Sum the values 
  string sdata = dataname;
  Double_t x2 = 1;
  Double_t x2binsize = (x2Max-x2Min)/x2Dim;
   
  for(UInt_t i=0;i<x1Dim;i++) {
    for(UInt_t j=0;j<x2DimSl;j++) {
     
      if(opt.Contains("cyl") && (opt.Contains("sum") || opt.Contains("int") ) ) {
      	// For cylindrical coordinates, this is the radius in Osiris units:
	x2 = x2binsize * (j-0.5) + x2Min;
      }
      
      if(sdata.find("charge") != string::npos || sdata.find("p1x1") != string::npos || sdata.find("p2x2") != string::npos) {
	x1Array[i] += -data[j*x1Dim + i] * x2;
      } else
	x1Array[i] +=  data[j*x1Dim + i] * x2;
    }
  }
  
  // Fill the array for the TH1F histo
  // Note: The Y values of the slice are averaged.
  
  // Histogram centering
  Double_t shiftx1 = Shift(opt);  
  x1Min -= shiftx1;
  x1Max -= shiftx1;
  
  TH1F* h1D = new TH1F();
  UInt_t x1DimAvg = x1Dim/Navg;
  h1D->SetBins(x1DimAvg,x1Min,x1Max);
  
  for(UInt_t i=0;i<x1Dim;i++) {
    UInt_t iavg = i/Navg;
    if(i%Navg) continue;
    Double_t content = x1Array[i];
    if(opt.Contains("avg")) content /= (Lastx2Bin-Firstx2Bin+1);
    else if(opt.Contains("int")) content *= x2binsize;
    h1D->SetBinContent(iavg,content);
  }
  
  delete x1Array;
  delete count;
  delete offset;
  delete dataDims;
  delete dataSet;
  delete root;

  return h1D;
}


//_______________________________________________________________________
TH1F* PData::GetH1SliceX(const char *filename,const char *dataname,Int_t Firstx1Bin, Int_t Lastx1Bin, const char *options) {

  // Options
  TString opt = options;

  // Open input HDF5 file
  H5File h5 = H5File(filename,H5F_ACC_RDONLY);
  
  // Open "axis" group 
  Group *axis = new Group(h5.openGroup("/AXIS"));

  // Get axes
  Double_t ax2lims[2];
  DataSet  *ax2 = new DataSet(axis->openDataSet("AXIS2"));
  hsize_t  *ax2Dims = new hsize_t[1];
  ax2->getSpace().getSimpleExtentDims(ax2Dims,NULL);
  DataSpace mem2(1,ax2Dims);
  ax2->read(ax2lims,PredType::NATIVE_DOUBLE,mem2,ax2->getSpace());
  ax2->close(); 
  delete ax2Dims;
  delete ax2;

  // The "int" option integrates the in the range of bins specified by
  // Firstx2Bin and Lastx2Bin. 
  // We need then the binning in axis2 to properly weight the values.
  Double_t ax1lims[2] = {0.,0.};
  if(opt.Contains("int")) {
    DataSet  *ax1 = new DataSet(axis->openDataSet("AXIS1"));
    hsize_t  *ax1Dims = new hsize_t[1];
    ax1->getSpace().getSimpleExtentDims(ax1Dims,NULL);
    DataSpace mem1(1,ax1Dims);
    ax1->read(ax1lims,PredType::NATIVE_DOUBLE,mem1,ax1->getSpace());
    ax1->close(); 
    delete ax1Dims;
    delete ax1;
  }
  axis->close();
  delete axis;

  // Open main group for data reading
  Group *root = new Group(h5.openGroup("/"));

  // Read data
  // ----------
  DataSet  *dataSet   = new DataSet(root->openDataSet(dataname));
  DataSpace dataSpace = dataSet->getSpace();
  const DataType type = dataSet->getDataType();

  UInt_t   rank     = dataSpace.getSimpleExtentNdims();
  hsize_t *dataDims = new hsize_t[rank];
  dataSpace.getSimpleExtentDims(dataDims,NULL);

  // x1 slice 
  if(Lastx1Bin>=(Int_t)dataDims[1]) Lastx1Bin = dataDims[1] - 1;
  if(Firstx1Bin<0) {
    DoSlice(dataDims[1],Firstx1Bin,Lastx1Bin);
  }
  if(Lastx1Bin<Firstx1Bin) 
    Lastx1Bin = Firstx1Bin;

  UInt_t x1Dim = Lastx1Bin - Firstx1Bin + 1; 

  UInt_t x2Dim = GetX2N();
  Double_t x2Min = GetX2Min();
  Double_t x2Max = GetX2Max();

  hsize_t  *count  = new hsize_t[rank];
  count[0] = x2Dim;
  count[1] = x1Dim;

  hsize_t  *offset = new hsize_t[rank];
  offset[0] = GetX2iMin();
  offset[1] = Firstx1Bin;

  Float_t *data = new Float_t[x2Dim*x1Dim];  
  memset(data,0,sizeof(Float_t)*x2Dim*x1Dim);
  
  DataSpace memSpace(rank,dataDims);
  dataSpace.selectHyperslab(H5S_SELECT_SET,count,offset);
  dataSet->read(data,type,memSpace,dataSpace);
  dataSet->close();

  root->close();

  Float_t *x2Array = new Float_t[x2Dim];
  memset(x2Array,0,sizeof(Float_t)*x2Dim);

  // Sum the values 
  string sdata = dataname;
  for(UInt_t i=0;i<x1Dim;i++) {
    for(UInt_t j=0;j<x2Dim;j++) {
      if(sdata.find("charge") != string::npos || sdata.find("p1x1") != string::npos || sdata.find("p2x2") != string::npos)
	x2Array[j] += -data[j*x1Dim + i];
      else
	x2Array[j] += data[j*x1Dim + i];
    }
  }

  // Fill the array for the TH1F histo

  // Histogram centering
  Double_t shiftx2 = 0;
  if(opt.Contains("center")) 
    if(!opt.Contains("cyl"))
      shiftx2 += (ax2lims[0] + ax2lims[1])/2.;
  
  x2Min -= shiftx2;
  x2Max -= shiftx2;
  
  TH1F* h1D = new TH1F();
  h1D->SetBins(x2Dim,x2Min,x2Max);
  
  Double_t x1binsize = (ax1lims[1]-ax1lims[0])/x1Dim;
  for(UInt_t j=0;j<x2Dim;j++) {
    Double_t content = x2Array[j];
    if(opt.Contains("avg")) content /= (Lastx1Bin-Firstx1Bin+1);
    else if(opt.Contains("int")) content *= x1binsize;
    h1D->SetBinContent(j+1,content);
  }
  
  root->close();
  
  delete x2Array;
  delete count;
  delete offset;
  delete dataDims;
  delete dataSet;
  delete root;

  return h1D;
}
//_______________________________________________________________________
TTree* PData::GetRawTree(const char *filename, const char *options) {
  
  // Options
  TString opt = options;

  // Open input HDF5 file
  // H5File h5 = H5File(filename,H5F_ACC_RDONLY);
  
  // Open main group for data reading
  // Group *root = new Group(h5.openGroup("/"));
  
  // TTree *tree = new TTree("RawTree","");
  // - not working -
  // GroupToTree(root,tree,kTRUE,kTRUE);  // kTRUE for sequential reading.

  UInt_t Nvar = 7;
  if(!Is3D()) Nvar = 6;
  char varname[7][4] = {{"p1"},{"p2"},{"p3"},{"q"},{"x1"},{"x2"},{"x3"}};

  
  Float_t **var;
  var = new Float_t*[Nvar];
  UInt_t Np = GetRawArray(filename,var);

  TTree *tree = new TTree("RawTree","");
  
  // Define the branches
  Float_t *darray = new Float_t[Nvar];
  for(UInt_t i=0;i<Nvar;i++) {
    tree->Branch(varname[i],&darray[i],Form("%s/F",varname[i]));
  }

  for(UInt_t j=0;j<Np;j++) {
    for(UInt_t i=0;i<Nvar;i++)
      darray[i] = var[i][j];
    
    tree->Fill();
  }
  
  // tree->Print();

  return tree;
}

//_______________________________________________________________________
UInt_t PData::GetRawArray(const char *filename, Float_t **var, const char *options) {
  
  // Options
  TString opt = options;

  // Open input HDF5 file
  H5File h5 = H5File(filename,H5F_ACC_RDONLY);
  
  // Open main group for data reading
  Group *root = new Group(h5.openGroup("/"));

  //UInt_t Ndata = root->getNumObjs();
  //cout << Form(" %i  datasets", Ndata) << endl;
  
  UInt_t Nvar = 7;
  if(!Is3D()) Nvar = 6;
  
  char varname[8][4] = {{"p1"},{"p2"},{"p3"},{"q"},{"x1"},{"x2"},{"x3"}};

  DataSet *dataSet[Nvar];
  UInt_t Np = 0;
  
  for(UInt_t i=0;i<Nvar;i++) {
    
    dataSet[i] = new DataSet(root->openDataSet(varname[i]));
    Int_t rank = dataSet[i]->getSpace().getSimpleExtentNdims();
    hsize_t  *dataDims = new hsize_t[rank];
    dataSet[i]->getSpace().getSimpleExtentDims(dataDims,NULL);
    if(i==0) Np = dataDims[0];  
    
    var[i] = new Float_t[Np];
    
    DataSpace memspace(rank,dataDims);
    dataSet[i]->read(var[i],PredType::NATIVE_FLOAT,memspace,dataSet[i]->getSpace());
    dataSet[i]->close();

  }
  root->close();

  return Np;
  
}

//_______________________________________________________________________
UInt_t PData::GetRawSingleArray(const char *filename, Float_t **var,const char *varname,const char *options) {
  
  // Options
  TString opt = options;

  // Open input HDF5 file
  H5File h5 = H5File(filename,H5F_ACC_RDONLY);
  
  // Open main group for data reading
  Group *root = new Group(h5.openGroup("/"));
 
  DataSet *dataSet = new DataSet(root->openDataSet(varname));
  Int_t rank = dataSet->getSpace().getSimpleExtentNdims();
  hsize_t  *dataDims = new hsize_t[rank];
  dataSet->getSpace().getSimpleExtentDims(dataDims,NULL);
  UInt_t Np = dataDims[0];  
  
  *var = new Float_t[Np];
  
  DataSpace memspace(rank,dataDims);
  dataSet->read(*var,PredType::NATIVE_FLOAT,memspace,dataSet->getSpace());
  dataSet->close();
  root->close();

  return Np;
  
}

//_______________________________________________________________________
void PData::GetH1Raw(const char *filename,const char *dataname, TH1F *h1D, const char *options) {
 
  // Options
  TString opt = options;

  // Open input HDF5 file
  H5File h5 = H5File(filename,H5F_ACC_RDONLY);
  
  // Open main group for data reading
  Group *root = new Group(h5.openGroup("/"));

  // Read data
  // ----------
  char dname[10];
  strcpy(dname,dataname);
  string sdata = dataname;
  if(sdata.find("energy") != string::npos) {
    strcpy(dname,"ene");
  }
  DataSet  *dataSet = new DataSet(root->openDataSet(dname));
  DataSpace dataSpace = dataSet->getSpace();
  const DataType type = dataSet->getDataType();

  Int_t     rank = dataSpace.getSimpleExtentNdims();
  hsize_t  *dataDims = new hsize_t[rank];
  dataSpace.getSimpleExtentDims(dataDims,NULL);

  UInt_t    nDim = dataDims[0];
  if(nDim == 0) return;

  hsize_t  *count  = new hsize_t[rank];
  hsize_t  *offset = new hsize_t[rank];
  for(Int_t j=0;j<rank;j++) {
    offset[j] = 0;
    count[j]  = dataDims[j];
  }

  Float_t *data = new Float_t[nDim];
  DataSpace memspace(rank,dataDims);
  dataSpace.selectHyperslab(H5S_SELECT_SET,count,offset);
  dataSet->read(data,type,memspace,dataSpace);
  dataSet->close();

  // The weights are also needed in any case.
  DataSet  *weightSet = new DataSet(root->openDataSet("q"));
  DataSpace weightSpace = weightSet->getSpace();
  DataSpace wmemspace(rank,dataDims);
  Float_t  *weight = new Float_t[nDim];
  weightSpace.selectHyperslab(H5S_SELECT_SET,count,offset);
  weightSet->read(weight,type,wmemspace,weightSpace);
  weightSet->close();

  // And the x2 coordinate in case of cylindrical coordinates
  Float_t *radius = NULL;
  if(opt.Contains("cyl")) {
    radius = new Float_t[nDim];
    DataSet  *radiusSet = new DataSet(root->openDataSet("x2"));
    DataSpace radiusSpace = radiusSet->getSpace();
    DataSpace rmemSpace(rank,dataDims);
    radiusSpace.selectHyperslab(H5S_SELECT_SET,count,offset);
    radiusSet->read(radius,type,rmemSpace,radiusSpace);
    radiusSet->close();
    delete radiusSet;
  }
  
  // Histogram centering
  Double_t shiftx1 = Shift(opt);   
  for(UInt_t i=0;i<nDim;i++) {
    if(weight[i]<0.) weight[i] *= -1;
    // if(opt.Contains("cyl")) weight[i] *= radius[i];
    if(sdata.find("energy") != string::npos) data[i] = (data[i]+1) * GetBeamMass() / PUnits::MeV;
    if(sdata.find("x1") != string::npos) data[i] -= shiftx1;
    h1D->Fill(data[i],weight[i]);
  }

  // Free memory:
  delete weight;
  delete weightSet;
  delete data;
  delete offset;
  delete count;
  delete dataDims;
  delete dataSet;
  if(radius) delete radius;
  delete root;
}


//_______________________________________________________________________
void PData::GetH1RawCut(const char *filename,const char *dataname, 
			const char *cutname, Float_t cutMin, Float_t cutMax,
			TH1F *h1D,const char *options) {
 
  // Options
  TString opt = options;

  // Open input HDF5 file
  H5File h5 = H5File(filename,H5F_ACC_RDONLY);
  
  // Open main group for data reading
  Group *root = new Group(h5.openGroup("/"));

  // Read data
  // ----------
  char dname[10];
  strcpy(dname,dataname);
  string sdata = dataname;
  if(sdata.find("energy") != string::npos) {
    strcpy(dname,"ene");
  }
  DataSet  *dataSet = new DataSet(root->openDataSet(dname));
  DataSpace dataSpace = dataSet->getSpace();
  const DataType type = dataSet->getDataType();

  Int_t     rank = dataSpace.getSimpleExtentNdims();
  hsize_t  *dataDims = new hsize_t[rank];
  dataSpace.getSimpleExtentDims(dataDims,NULL);

  UInt_t    nDim = dataDims[0];
  hsize_t  *count  = new hsize_t[rank];
  hsize_t  *offset = new hsize_t[rank];
  for(Int_t j=0;j<rank;j++) {
    offset[j] = 0;
    count[j]  = dataDims[j];
  }

  Float_t  *data = new Float_t[nDim];
  DataSpace memspace(rank,dataDims);
  dataSpace.selectHyperslab(H5S_SELECT_SET,count,offset);
  dataSet->read(data,type,memspace,dataSpace);
  dataSet->close();

  // The weights are also needed in any case.
  DataSet  *weightSet = new DataSet(root->openDataSet("q"));
  DataSpace weightSpace = weightSet->getSpace();
  DataSpace wmemspace(rank,dataDims);
  Float_t  *weight = new Float_t[nDim];
  weightSpace.selectHyperslab(H5S_SELECT_SET,count,offset);
  weightSet->read(weight,type,wmemspace,weightSet->getSpace());
  weightSet->close();

  // And the x2 coordinate in case of cylindrical coordinates
  Float_t *radius = NULL;
  if(opt.Contains("cyl")) {
    radius = new Float_t[nDim];
    DataSet  *radiusSet = new DataSet(root->openDataSet("x2"));
    DataSpace radiusSpace = radiusSet->getSpace();
    DataSpace rmemSpace(rank,dataDims);
    radiusSpace.selectHyperslab(H5S_SELECT_SET,count,offset);
    radiusSet->read(radius,type,rmemSpace,radiusSpace);
    radiusSet->close(); 
    delete radiusSet;
  }

  // And the variable where to cut.
  char cname[10];
  strcpy(cname,cutname);
  string scut = cutname;
  if(scut.find("energy") != string::npos) {
    strcpy(cname,"ene");
  }
  DataSet   *cutSet = new DataSet(root->openDataSet(cname));
  DataSpace  cutSpace = cutSet->getSpace();
  DataSpace  cmemspace(rank,dataDims);
  Float_t   *cut = new Float_t[nDim];
  cutSpace.selectHyperslab(H5S_SELECT_SET,count,offset);
  cutSet->read(cut,type,cmemspace,cutSpace);
  cutSet->close();

  // Histogram centering
  Double_t shiftx1 = Shift(opt);
  
  for(UInt_t i=0;i<nDim;i++) {
    if(scut.find("energy") != string::npos) cut[i] = (cut[i]+1) * GetBeamMass() / PUnits::MeV;
    if(scut.find("x1") != string::npos) cut[i] -= shiftx1;
    if(cut[i]<cutMin || cut[i]>cutMax) continue;
    
    if(weight[i]<0.) weight[i] *= -1;
    // if(opt.Contains("cyl")) weight[i] *= radius[i];
    if(sdata.find("energy") != string::npos) data[i] = (data[i]+1) * GetBeamMass() / PUnits::MeV;
    if(sdata.find("x1") != string::npos) data[i] -= shiftx1;
    
    h1D->Fill(data[i],weight[i]);
  }

  // Free memory:
  delete weight;
  delete weightSet;
  delete data;
  delete offset;
  delete count;
  delete dataDims;
  delete dataSet;
  delete cut;
  delete cutSet;
  if(radius) delete radius;
  delete root;

}

//_______________________________________________________________________
void PData::GetH1RawCut2(const char *filename,const char *dataname, 
			 const char *cutname1, Float_t cutMin1, Float_t cutMax1,
			 const char *cutname2, Float_t cutMin2, Float_t cutMax2,
			 TH1F *h1D,const char *options) {
 
  // Options
  TString opt = options;

  // Open input HDF5 file
  H5File h5 = H5File(filename,H5F_ACC_RDONLY);
  
  // Open main group for data reading
  Group *root = new Group(h5.openGroup("/"));

  // Read data
  // ----------
  char dname[10];
  strcpy(dname,dataname);
  string sdata = dataname;
  if(sdata.find("energy") != string::npos) {
    strcpy(dname,"ene");
  }
  DataSet  *dataSet = new DataSet(root->openDataSet(dname));
  DataSpace dataSpace = dataSet->getSpace();
  const DataType type = dataSet->getDataType();

  Int_t     rank = dataSpace.getSimpleExtentNdims();
  hsize_t  *dataDims = new hsize_t[rank];
  dataSpace.getSimpleExtentDims(dataDims,NULL);

  UInt_t    nDim = dataDims[0];
  hsize_t  *count  = new hsize_t[rank];
  hsize_t  *offset = new hsize_t[rank];
  for(Int_t j=0;j<rank;j++) {
    offset[j] = 0;
    count[j]  = dataDims[j];
  }

  Float_t  *data = new Float_t[nDim];
  DataSpace memspace(rank,dataDims);
  dataSpace.selectHyperslab(H5S_SELECT_SET,count,offset);
  dataSet->read(data,type,memspace,dataSpace);
  dataSet->close();

  // The weights are also needed in any case.
  DataSet  *weightSet = new DataSet(root->openDataSet("q"));
  DataSpace weightSpace = weightSet->getSpace();
  DataSpace wmemspace(rank,dataDims);
  Float_t  *weight = new Float_t[nDim];
  weightSpace.selectHyperslab(H5S_SELECT_SET,count,offset);
  weightSet->read(weight,type,wmemspace,weightSet->getSpace());
  weightSet->close();

  // And the x2 coordinate in case of cylindrical coordinates
  Float_t *radius = NULL;
  if(opt.Contains("cyl")) {
    radius = new Float_t[nDim];
    DataSet  *radiusSet = new DataSet(root->openDataSet("x2"));
    DataSpace radiusSpace = radiusSet->getSpace();
    DataSpace rmemSpace(rank,dataDims);
    radiusSpace.selectHyperslab(H5S_SELECT_SET,count,offset);
    radiusSet->read(radius,type,rmemSpace,radiusSpace);
    radiusSet->close(); 
    delete radiusSet;
  }

  // And the variable where to cut.
  char cname[10];
  strcpy(cname,cutname1);
  string scut1 = cutname1;
  if(scut1.find("energy") != string::npos) {
    strcpy(cname,"ene");
  }
  DataSet   *cutSet1 = new DataSet(root->openDataSet(cname));
  DataSpace  cutSpace1 = cutSet1->getSpace();
  DataSpace  cmemspace1(rank,dataDims);
  Float_t   *cut1 = new Float_t[nDim];
  cutSpace1.selectHyperslab(H5S_SELECT_SET,count,offset);
  cutSet1->read(cut1,type,cmemspace1,cutSpace1);
  cutSet1->close();

  strcpy(cname,cutname2);
  string scut2 = cutname2;
  if(scut2.find("energy") != string::npos) {
    strcpy(cname,"ene");
  }
  DataSet   *cutSet2 = new DataSet(root->openDataSet(cname));
  DataSpace  cutSpace2 = cutSet2->getSpace();
  DataSpace  cmemspace2(rank,dataDims);
  Float_t   *cut2 = new Float_t[nDim];
  cutSpace2.selectHyperslab(H5S_SELECT_SET,count,offset);
  cutSet2->read(cut2,type,cmemspace2,cutSpace2);
  cutSet2->close();

  // Histogram centering
  Double_t shiftx1 = Shift(opt);
  
  for(UInt_t i=0;i<nDim;i++) {
    if(scut1.find("energy") != string::npos) cut1[i] = (cut1[i]+1) * GetBeamMass() / PUnits::MeV;
    if(scut1.find("x1") != string::npos) cut1[i] -= shiftx1;
    if(cut1[i]<cutMin1 || cut1[i]>cutMax1) continue;

    if(scut2.find("energy") != string::npos) cut2[i] = (cut2[i]+1) * GetBeamMass() / PUnits::MeV;
    if(scut2.find("x1") != string::npos) cut2[i] -= shiftx1;
    if(cut2[i]<cutMin2 || cut2[i]>cutMax2) continue;
    
    if(weight[i]<0.) weight[i] *= -1;
    // if(opt.Contains("cyl")) weight[i] *= radius[i];
    if(sdata.find("energy") != string::npos) data[i] = (data[i]+1) * GetBeamMass() / PUnits::MeV;
    if(sdata.find("x1") != string::npos) data[i] -= shiftx1;
    
    h1D->Fill(data[i],weight[i]);
  }

  // Free memory:
  delete weight;
  delete weightSet;
  delete data;
  delete offset;
  delete count;
  delete dataDims;
  delete dataSet;
  delete cut1;
  delete cutSet1;
  delete cut2;
  delete cutSet2;
  if(radius) delete radius;
  delete root;

}

//_______________________________________________________________________
TTree** PData::GetTrackTree(const char *filename,Int_t &NTracks, const char *options) {
  
  // Options
  TString opt = options;

  // Open input HDF5 file
  H5File h5 = H5File(filename,H5F_ACC_RDONLY);
  
  // Open main group for data reading
  Group *root = new Group(h5.openGroup("/"));

  UInt_t Nvar = 8;
  if(!Is3D()) Nvar = 7;
  char varname[8][4] = {{"t"},{"p1"},{"p2"},{"p3"},{"q"},{"x1"},{"x2"},{"x3"}};

  Int_t NOBJ = root->getNumObjs();
  TTree **tree = new TTree*[NOBJ];
  Int_t NT = 0;
  for(Int_t i=0; (i<NOBJ) && ( (NTracks<0) || (NT<NTracks) ); i++) {
    char gname[128];
    root->getObjnameByIdx( i, gname, 128 );
    //cout << Form(" Group name = %s",gname) << endl;
    Group *track = new Group(h5.openGroup(Form("/%s",gname)));
    if(track == NULL) continue;

    Double_t **var = new Double_t*[Nvar];
    
    DataSet *dataSet = NULL;
    UInt_t Np = 0;
    
    for(UInt_t i=0;i<Nvar;i++) {
      
      dataSet = new DataSet(track->openDataSet(varname[i]));
      Int_t rank = dataSet->getSpace().getSimpleExtentNdims();
      hsize_t  *dataDims = new hsize_t[rank];
      dataSet->getSpace().getSimpleExtentDims(dataDims,NULL);
      if(i==0) Np = dataDims[0];  
    
      var[i] = new Double_t[Np];
      
      DataSpace memspace(rank,dataDims);
      dataSet->read(var[i],PredType::NATIVE_DOUBLE,memspace,dataSet->getSpace());
      dataSet->close();

      delete dataSet;
    }
    track->close(); 

    tree[NT] = new TTree(Form("TrackTree_%i",NT),"");
    // Define the branches
    Double_t *darray = new Double_t[Nvar];
    for(UInt_t i=0;i<Nvar;i++) {
      tree[NT]->Branch(varname[i],&darray[i],Form("%s/D",varname[i]));
    }

    // Fill the tree
    for(UInt_t j=0;j<Np;j++) {
      for(UInt_t i=0;i<Nvar;i++)
	darray[i] = var[i][j];
      
      tree[NT]->Fill();
    }
    
    NT++;
    
    delete track;
  }
  
  NTracks = NT;
  return tree;
}


//_______________________________________________________________________
UInt_t PData::GetTrackArray(const char *filename, Double_t ***var, Int_t *np, const char *options) {
  
  // Options
  TString opt = options;

  // Open input HDF5 file
  H5File h5 = H5File(filename,H5F_ACC_RDONLY);
  
  // Open main group for data reading
  Group *root = new Group(h5.openGroup("/"));
  Int_t NOBJ = root->getNumObjs();
  var = new Double_t**[NOBJ];
  np = new Int_t[NOBJ];
  for(Int_t i=0; (i<NOBJ); i++) {
    char gname[128];
    root->getObjnameByIdx( i, gname, 128 );
    Group *track = new Group(h5.openGroup(Form("/%s",gname)));
    if(track == NULL) continue;
   
    UInt_t Nvar = 8;
    if(!Is3D()) Nvar = 7;
    char varname[8][4] = {{"t"},{"p1"},{"p2"},{"p3"},{"q"},{"x1"},{"x2"},{"x3"}};

    var[i] = new Double_t*[Nvar];
    
    DataSet *dataSet[Nvar];
    UInt_t Ntime = 0;
    for(UInt_t j=0;j<Nvar;j++) {
    
      dataSet[j] = new DataSet(track->openDataSet(varname[j]));
      Int_t rank = dataSet[j]->getSpace().getSimpleExtentNdims();
      hsize_t  *dataDims = new hsize_t[rank];
      dataSet[j]->getSpace().getSimpleExtentDims(dataDims,NULL);
      if(j==0) {
	Ntime = dataDims[0];  
	np[i] = Ntime;
      }
      var[i][j] = new Double_t[Ntime];
      
      DataSpace memspace(rank,dataDims);
      dataSet[j]->read(var[i][j],PredType::NATIVE_DOUBLE,memspace,dataSet[j]->getSpace());

      dataSet[j]->close();
      
    }
    track->close();
  }
  root->close();

  return NOBJ;
  
}

//_______________________________________________________________________
TH2F* PData::GetH2(const char *filename,const char *dataname, const char *options) {
  
  // Options
  TString opt = options;
  
  // Open input HDF5 file
  H5File h5 = H5File(filename,H5F_ACC_RDONLY);
  
  // Open main group for data reading
  Group *root = new Group(h5.openGroup("/"));

  // Read data
  // ----------
  DataSet  *dataSet = new DataSet(root->openDataSet(dataname));
  DataSpace dataSpace = dataSet->getSpace();
  const DataType type = dataSet->getDataType();
  
  Int_t rank = dataSpace.getSimpleExtentNdims();
  hsize_t  *dataDims = new hsize_t[rank];
  dataSpace.getSimpleExtentDims(dataDims,NULL);
  
  UInt_t x1AvgFactor = GetNX(0)/dataDims[1];
  UInt_t x1Dim = GetX1N()/x1AvgFactor;
  Double_t x1Min = GetX1Min();
  Double_t x1Max = GetX1Max();

  // cout << Form(" NX1 = %i,   NX1range = %i,   NX1avg = %i,   NX1rangeAvg = %i",
  //  	       GetNX(0), GetX1N(), dataDims[1], x1Dim) << endl;
  
  UInt_t x2AvgFactor = GetNX(1)/dataDims[0];
  UInt_t x2Dim = GetX2N()/x2AvgFactor;
  Double_t x2Min = GetX2Min();
  Double_t x2Max = GetX2Max();

  // cout << Form(" NX2 = %i,   NX2range = %i,   NX2avg = %i,   NX2rangeAvg = %i",
  //  	       GetNX(1), GetX2N(), dataDims[0], x2Dim) << endl;

  hsize_t  *count  = new hsize_t[rank];
  count[0] = x2Dim;
  count[1] = x1Dim;
  
  hsize_t  *offset = new hsize_t[rank];  
  offset[0] = GetX2iMin()/x2AvgFactor;
  offset[1] = GetX1iMin()/x1AvgFactor;

  // cout << Form(" x1Dim = %i  offset = %i", count[1], offset[1]) << endl;
  // cout << Form(" x2Dim = %i  offset = %i", count[0], offset[0]) << endl;
  
  Float_t *data = new Float_t[x2Dim*x1Dim];   // <-- Solution: 1D arrays
  memset(data,0,sizeof(Float_t)*x2Dim*x1Dim);  
  DataSpace memSpace(rank,count);

  dataSpace.selectHyperslab(H5S_SELECT_SET,count,offset);
  dataSet->read(data,PredType::NATIVE_FLOAT,memSpace,dataSpace);
  dataSet->close();

  root->close();

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
 
  string sdata = dataname;
  TH2F *h2D = new TH2F();
  UInt_t x1DimAvg = x1Dim/Navg;
  h2D->SetBins(x1DimAvg,x1Min,x1Max,x2Dim,x2Min,x2Max);

  for(UInt_t i=0;i<x1Dim;i++) {
    UInt_t iavg = i/Navg;
    if(i%Navg) continue;
    for(UInt_t j=0;j<x2Dim;j++) {
      UInt_t index = (long)j*(long)x1Dim + (long)i;
      if(sdata.find("charge") != string::npos || sdata.find("p1x1") != string::npos || sdata.find("p2x2") != string::npos)
      	h2D->SetBinContent(iavg+1,j+1,-data[index]);
      else
	h2D->SetBinContent(iavg+1,j+1,data[index]);
      
    } 
  }

  // Free memory:
  delete dataDims;
  delete dataSet;
  delete root;
 
  return h2D;
}

//_______________________________________________________________________
void PData::GetH2Raw(const char *filename,const char *dataname1,const char *dataname2,TH2F *h2D,const char *options) {
 
  // Options
  TString opt = options;

  // Open input HDF5 file
  H5File h5 = H5File(filename,H5F_ACC_RDONLY);
  
  // Open main group for data reading
  Group *root = new Group(h5.openGroup("/"));

  // Read data
  // ----------
  char dname1[10];
  strcpy(dname1,dataname1);
  char dname2[10];
  strcpy(dname2,dataname2);  
  
  string sdata1 = dataname1;
  if( (sdata1.find("energy") != string::npos) ||
      (sdata1.find("gamma") != string::npos) ) {
    strcpy(dname1,"ene");
  }
 

  string sdata2 = dataname2;
  if( (sdata2.find("energy") != string::npos) ||
      (sdata2.find("gamma") != string::npos) ) {
    strcpy(dname2,"ene");
  }
  
  // 1st variable
  DataSet  *dataSet = new DataSet(root->openDataSet(dname1));
  DataSpace dataSpace = dataSet->getSpace();
  const DataType type = dataSet->getDataType();
  
  Int_t     rank = dataSpace.getSimpleExtentNdims();
  hsize_t  *dataDims = new hsize_t[rank];
  dataSpace.getSimpleExtentDims(dataDims,NULL);

  UInt_t    nDim = dataDims[0];
  hsize_t  *count  = new hsize_t[rank];
  hsize_t  *offset = new hsize_t[rank];
  for(Int_t j=0;j<rank;j++) {
    offset[j] = 0;
    count[j]  = dataDims[j];
  }

  Float_t *data1 = new Float_t[nDim];
  DataSpace memspace(rank,dataDims);
  dataSpace.selectHyperslab(H5S_SELECT_SET,count,offset);
  dataSet->read(data1,type,memspace,dataSpace);
  dataSet->close();

  // 2nd variable
  DataSet  *dataSet2 = new DataSet(root->openDataSet(dname2));
  DataSpace dataSpace2 = dataSet2->getSpace();
  Float_t  *data2 = new Float_t[nDim];
  DataSpace memSpace2(rank,dataDims);
  dataSpace2.selectHyperslab(H5S_SELECT_SET,count,offset);
  dataSet2->read(data2,type,memSpace2,dataSpace2);
  dataSet2->close();

  // The weights are also needed in any case.
  DataSet  *weightSet = new DataSet(root->openDataSet("q"));
  DataSpace weightSpace = weightSet->getSpace();
  DataSpace wmemspace(rank,dataDims);
  Float_t  *weight = new Float_t[nDim];
  weightSpace.selectHyperslab(H5S_SELECT_SET,count,offset);
  weightSet->read(weight,type,wmemspace,weightSpace);
  weightSet->close();

  // And the x2 coordinate in case of cylindrical coordenates
  Float_t *radius = NULL;
  if(opt.Contains("cyl")) {
    radius = new Float_t[nDim];
    DataSet  *radiusSet = new DataSet(root->openDataSet("x2"));
    DataSpace radiusSpace = radiusSet->getSpace();
    DataSpace rmemSpace(rank,dataDims);
    radiusSpace.selectHyperslab(H5S_SELECT_SET,count,offset);
    radiusSet->read(radius,type,rmemSpace,radiusSpace);
    radiusSet->close();
    delete radiusSet;
  }

  // Histogram centering
  Double_t shiftx1 = Shift(opt);  
  for(UInt_t i=0;i<nDim;i++) {
    if(weight[i]<0.) weight[i] *= -1;
    // if(opt.Contains("cyl")) weight[i] *= radius[i];
    if(sdata1.find("energy") != string::npos) data1[i] = (data1[i]+1) * GetBeamMass() / PUnits::MeV;
    if(sdata2.find("energy") != string::npos) data2[i] = (data2[i]+1) * GetBeamMass() / PUnits::MeV;
    if(sdata1.find("gamma") != string::npos) data1[i] = (data1[i]+1);
    if(sdata2.find("gamma") != string::npos) data2[i] = (data2[i]+1);
    if(sdata1.find("x1") != string::npos) data1[i] -= shiftx1;
    if(sdata2.find("x1") != string::npos) data2[i] -= shiftx1;

    //  cout << sdata1 << " = " << data1[i] << "  " << sdata2 << " = " << data2[i] << endl;
    
    h2D->Fill(data1[i],data2[i],weight[i]);
  }

  // Free memory:
  delete weight;
  delete weightSet;
  delete data2;
  delete dataSet2;
  delete data1;
  delete offset;
  delete count;
  delete dataDims;
  delete dataSet;
  if(radius) delete radius;
  delete root;
}

//_______________________________________________________________________
void PData::GetH2RawCut(const char *filename,const char *dataname1,const char *dataname2,
			const char *cutname, Float_t cutMin, Float_t cutMax,
			TH2F *h2D,const char *options) {
 
  // Options
  TString opt = options;

  // Open input HDF5 file
  H5File h5 = H5File(filename,H5F_ACC_RDONLY);
  
  // Open main group for data reading
  Group *root = new Group(h5.openGroup("/"));

  // Read data
  // ----------
  char dname1[10];
  strcpy(dname1,dataname1);
  char dname2[10];
  strcpy(dname2,dataname2);  
  
  string sdata1 = dataname1;
  if(sdata1.find("energy") != string::npos) {
    strcpy(dname1,"ene");
  }
  string sdata2 = dataname2;
  cout << dataname2 << endl;
  if(sdata2.find("energy") != string::npos) {
    strcpy(dname2,"ene");
  }

  // 1st variable
  DataSet  *dataSet = new DataSet(root->openDataSet(dname1));
  DataSpace dataSpace = dataSet->getSpace();
  const DataType type = dataSet->getDataType();

  Int_t     rank = dataSpace.getSimpleExtentNdims();
  hsize_t  *dataDims = new hsize_t[rank];
  dataSpace.getSimpleExtentDims(dataDims,NULL);

  UInt_t    nDim = dataDims[0];
  hsize_t  *count  = new hsize_t[rank];
  hsize_t  *offset = new hsize_t[rank];
  for(Int_t j=0;j<rank;j++) {
    offset[j] = 0;
    count[j]  = dataDims[j];
  }

  Float_t *data1 = new Float_t[nDim];
  DataSpace memspace(rank,dataDims);
  dataSpace.selectHyperslab(H5S_SELECT_SET,count,offset);
  dataSet->read(data1,type,memspace,dataSpace);
  dataSet->close();

  // 2nd variable
  DataSet  *dataSet2 = new DataSet(root->openDataSet(dname2));
  DataSpace dataSpace2 = dataSet2->getSpace();
  Float_t  *data2 = new Float_t[nDim];
  DataSpace memSpace2(rank,dataDims);
  dataSpace2.selectHyperslab(H5S_SELECT_SET,count,offset);
  dataSet2->read(data2,type,memSpace2,dataSpace2);
  dataSet2->close();

  // The weights are also needed in any case.
  DataSet  *weightSet = new DataSet(root->openDataSet("q"));
  DataSpace weightSpace = weightSet->getSpace();
  DataSpace wmemspace(rank,dataDims);
  Float_t  *weight = new Float_t[nDim];
  weightSpace.selectHyperslab(H5S_SELECT_SET,count,offset);
  weightSet->read(weight,type,wmemspace,weightSpace);
  weightSet->close();

  // And the x2 coordinate in case of cylindrical coordenates
  Float_t *radius = NULL;
  if(opt.Contains("cyl")) {
    radius = new Float_t[nDim];
    DataSet  *radiusSet = new DataSet(root->openDataSet("x2"));
    DataSpace radiusSpace = radiusSet->getSpace();
    DataSpace rmemSpace(rank,dataDims);
    radiusSpace.selectHyperslab(H5S_SELECT_SET,count,offset);
    radiusSet->read(radius,type,rmemSpace,radiusSpace);
    radiusSet->close();
    delete radiusSet;
  }

  // And the variable where to cut.
  char cname[10];
  strcpy(cname,cutname);
  string scut = cutname;
  if(scut.find("energy") != string::npos) {
    strcpy(cname,"ene");
  }
  DataSet   *cutSet = new DataSet(root->openDataSet(cname));
  DataSpace  cutSpace = cutSet->getSpace();
  DataSpace  cmemspace(rank,dataDims);
  Float_t   *cut = new Float_t[nDim];
  cutSpace.selectHyperslab(H5S_SELECT_SET,count,offset);
  cutSet->read(cut,type,cmemspace,cutSpace);
  cutSet->close();

  // Histogram centering
  Double_t shiftx1 = Shift(opt);
  for(UInt_t i=0;i<nDim;i++) {
    if(scut.find("energy") != string::npos) cut[i] = (cut[i]+1) * GetBeamMass() / PUnits::MeV;
    if(scut.find("x1") != string::npos) cut[i] -= shiftx1;
    if(cut[i]<cutMin || cut[i]>cutMax) continue;
  
    if(weight[i]<0.) weight[i] *= -1;
    // if(opt.Contains("cyl")) weight[i] *= radius[i];
    if(sdata1.find("energy") != string::npos) data1[i] = (data1[i]+1) * GetBeamMass() / PUnits::MeV;
    if(sdata2.find("energy") != string::npos) data2[i] = (data2[i]+1) * GetBeamMass() / PUnits::MeV;
    if(sdata1.find("x1") != string::npos) data1[i] -= shiftx1;
    if(sdata2.find("x1") != string::npos) data2[i] -= shiftx1;

    h2D->Fill(data1[i],data2[i],weight[i]);
  }

  // Free memory:
  delete weight;
  delete weightSet;
  delete data1;
  delete data2;
  delete offset;
  delete count;
  delete dataDims;
  delete dataSet;
  delete dataSet2;
  delete cut;
  delete cutSet;
  if(radius) delete radius;
  delete root;

}


//_______________________________________________________________________
TH1F* PData::GetH1SliceZ3D(const char *filename,const char *dataname, 
			   Int_t Firstx2Bin, Int_t Lastx2Bin, 
			   Int_t Firstx3Bin, Int_t Lastx3Bin, const char *options) {
  
  // Save memory and time with an specific range selection.

  // Options
  TString opt = options;
  
  // Open input HDF5 file
  H5File h5 = H5File(filename,H5F_ACC_RDONLY);
  
  // Open main group for data reading
  Group *root = new Group(h5.openGroup("/"));

  // Read data
  // ----------
  DataSet  *dataSet = new DataSet(root->openDataSet(dataname));
  DataSpace dataSpace = dataSet->getSpace();
  const DataType type = dataSet->getDataType();
    
  Int_t rank = dataSpace.getSimpleExtentNdims();
  hsize_t  *dataDims = new hsize_t[rank];
  dataSpace.getSimpleExtentDims(dataDims,NULL);

  // x3 slice
  if(Lastx3Bin>=(Int_t)dataDims[0]) Lastx3Bin = dataDims[0] - 1;
  
  if(Firstx3Bin<0)
    DoSlice(dataDims[0],Firstx3Bin,Lastx3Bin);
  
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
  UInt_t x1AvgFactor = GetNX(0)/dataDims[2];
  UInt_t x1Dim = GetX1N()/x1AvgFactor;
  Double_t x1Min = GetX1Min();
  Double_t x1Max = GetX1Max();
  
  hsize_t  *count  = new hsize_t[rank];
  count[0] = x3Dim;
  count[1] = x2Dim;
  count[2] = x1Dim;

  hsize_t  *offset = new hsize_t[rank];
  offset[0] = Firstx3Bin;
  offset[1] = Firstx2Bin;
  offset[2] = GetX1iMin()/x1AvgFactor;
  
  Float_t *data = new Float_t[x3Dim * x2Dim * x1Dim];

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
  for(UInt_t i=0;i<x3Dim;i++) {
    for(UInt_t j=0;j<x2Dim;j++) {
      for(UInt_t k=0;k<x1Dim;k++) {
	UInt_t index = (long)i*(long)x2Dim*(long)x1Dim + (long)j*(long)x1Dim + (long)k;
		
	if(sdata.find("charge") != string::npos || sdata.find("p1x1") != string::npos || sdata.find("p2x2") != string::npos)
	  x1Array[k] += -data[index];
	else
	  x1Array[k] +=  data[index]; 
	
      }
    }
  }
  
  // Histogram centering
  Double_t shiftx1 = Shift(opt);  

  x1Min -= shiftx1;
  x1Max -= shiftx1;

  Float_t x2AvgFactor = GetNX(1)/dataDims[1];
  Float_t dx2   = GetDX(1) * x2AvgFactor;

  Float_t x3AvgFactor = GetNX(2)/dataDims[0];
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
TH1F* PData::GetH1SliceX3D(const char *filename,const char *dataname, 
			   Int_t Firstx1Bin, Int_t Lastx1Bin, 
			   Int_t Firstx3Bin, Int_t Lastx3Bin, const char *options) {
  
  // Save memory and time with an specific range selection.

  // Options
  TString opt = options;
  
  // Open input HDF5 file
  H5File h5 = H5File(filename,H5F_ACC_RDONLY);
  
  // Open main group for data reading
  Group *root = new Group(h5.openGroup("/"));

  // Read data
  // ----------
  DataSet  *dataSet = new DataSet(root->openDataSet(dataname));
  DataSpace dataSpace = dataSet->getSpace();
  const DataType type = dataSet->getDataType();
    
  Int_t rank = dataSpace.getSimpleExtentNdims();
  hsize_t  *dataDims = new hsize_t[rank];
  dataSpace.getSimpleExtentDims(dataDims,NULL);

  // x3 slice
  if(Lastx3Bin>=(Int_t)dataDims[0]) Lastx3Bin = dataDims[0] - 1;
  
  if(Firstx3Bin<0)
    DoSlice(dataDims[0],Firstx3Bin,Lastx3Bin);
  
  if(Lastx3Bin<Firstx3Bin) 
    Lastx3Bin = Firstx3Bin;
  
  UInt_t x3Dim = Lastx3Bin - Firstx3Bin + 1; 
  
  //----------

  // x1 slice
  if(Lastx1Bin>=(Int_t)dataDims[2]) Lastx1Bin = dataDims[1] - 1;
  
  if(Firstx1Bin<0) {
    DoSlice(dataDims[1],Firstx1Bin,Lastx1Bin);
  }
  
  if(Lastx1Bin<Firstx1Bin) 
    Lastx1Bin = Firstx1Bin;

  UInt_t x1Dim = Lastx1Bin - Firstx1Bin + 1; 
 
  //----------
  UInt_t x2AvgFactor = GetNX(1)/dataDims[1];
  UInt_t x2Dim = GetX2N()/x2AvgFactor;
  Double_t x2Min = GetX2Min();
  Double_t x2Max = GetX2Max();
  
  hsize_t  *count  = new hsize_t[rank];
  count[0] = x3Dim;
  count[1] = x2Dim;
  count[2] = x1Dim;

  hsize_t  *offset = new hsize_t[rank];
  offset[0] = Firstx3Bin;
  offset[1] = GetX2iMin()/x2AvgFactor;
  offset[2] = Firstx1Bin;
  
  Float_t *data = new Float_t[x3Dim * x2Dim * x1Dim];

  DataSpace memSpace(rank,count);
  dataSpace.selectHyperslab(H5S_SELECT_SET,count,offset);
  dataSet->read(data,PredType::NATIVE_FLOAT,memSpace,dataSpace);
  dataSet->close();
  
  root->close();

  Float_t x2Array[x2Dim];
  // Set all elements to 0.
  memset(x2Array,0,sizeof(Float_t)*x2Dim);
  
  string sdata = dataname;
  // Sum the values 
  for(UInt_t i=0;i<x3Dim;i++) {
    for(UInt_t j=0;j<x2Dim;j++) {
      for(UInt_t k=0;k<x1Dim;k++) {
	UInt_t index = (long)i*(long)x2Dim*(long)x1Dim + (long)j*(long)x1Dim + (long)k;
		
	if(sdata.find("charge") != string::npos || sdata.find("p1x1") != string::npos || sdata.find("p2x2") != string::npos)
	  x2Array[j] += -data[index];
	else
	  x2Array[j] +=  data[index]; 
	
      }
    }
  }
  
  Float_t x1AvgFactor = GetNX(0)/dataDims[2];
  Float_t dx1   = GetDX(0) * x1AvgFactor;

  Float_t x3AvgFactor = GetNX(2)/dataDims[0];
  Float_t dx3   = GetDX(2) * x3AvgFactor;
  
  TH1F *h1D = new TH1F(); 
  h1D->SetBins(x2Dim,x2Min,x2Max);
  for(UInt_t k=0;k<x2Dim;k++) { 
    Double_t content = x2Array[k];
    if(opt.Contains("avg")) content /= x3Dim * x1Dim;
    else if(opt.Contains("int")) content *= dx1 * dx3;
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
TH1F* PData::GetH1CylSliceZ3D(const char *filename,const char *dataname, 
			      Int_t FirstrBin, Int_t LastrBin,
			      const char *options) {

  // Options
  TString opt = options;

  // Open input HDF5 file
  H5File h5 = H5File(filename,H5F_ACC_RDONLY);
  
  // Open "axis" group 
  Group *axis = new Group(h5.openGroup("/AXIS"));

  // Get axes
  Double_t ax1lims[2];
  DataSet  *ax1 = new DataSet(axis->openDataSet("AXIS1"));
  hsize_t  *ax1Dims = new hsize_t[1];
  ax1->getSpace().getSimpleExtentDims(ax1Dims,NULL);
  DataSpace mem1(1,ax1Dims);
  ax1->read(ax1lims,PredType::NATIVE_DOUBLE,mem1,ax1->getSpace());
  ax1->close(); 
  delete ax1Dims;
  delete ax1;

  Double_t ax2lims[2];
  DataSet  *ax2 = new DataSet(axis->openDataSet("AXIS2"));
  hsize_t  *ax2Dims = new hsize_t[1];
  ax2->getSpace().getSimpleExtentDims(ax2Dims,NULL);
  DataSpace mem2(1,ax2Dims);
  ax2->read(ax2lims,PredType::NATIVE_DOUBLE,mem2,ax2->getSpace());
  ax2->close();
  delete ax2Dims;
  delete ax2;

  Double_t ax3lims[2];
  DataSet  *ax3 = new DataSet(axis->openDataSet("AXIS3"));
  hsize_t  *ax3Dims = new hsize_t[1];
  ax3->getSpace().getSimpleExtentDims(ax3Dims,NULL);
  DataSpace mem3(1,ax3Dims);
  ax3->read(ax3lims,PredType::NATIVE_DOUBLE,mem3,ax3->getSpace());
  ax3->close();
  delete ax3Dims;
  delete ax3;

  axis->close();
  delete axis;

  // Open main group for data reading
  Group *root = new Group(h5.openGroup("/"));

  // Read data
  // ----------
  DataSet    *dataSet = new DataSet(root->openDataSet(dataname));
  DataSpace dataSpace = dataSet->getSpace();
  const DataType type = dataSet->getDataType();

  Int_t        rank = dataSpace.getSimpleExtentNdims();
  hsize_t *dataDims = new hsize_t[rank];
  dataSpace.getSimpleExtentDims(dataDims,NULL);

  UInt_t x1Dim = dataDims[2];
  UInt_t x2Dim = dataDims[1];
  UInt_t x3Dim = dataDims[0];

  // This is for sequential access to the 3D array.
  dataDims[0] = 1; 
  hsize_t  *count  = new hsize_t[rank];
  hsize_t  *offset = new hsize_t[rank];
  for(Int_t j=0;j<rank;j++) {
    offset[j] = 0;
    count[j]  = dataDims[j];
  }
  
  // Define radial bins
  UInt_t rDim = x2Dim; 
  // Define radial axis limits
  // Double_t axrlims[2] = {-(ax2lims[1]+ax2lims[0])/2.,(ax2lims[1]+ax2lims[0])/2.};

  if(LastrBin>=(Int_t)rDim) LastrBin = rDim - 1;

  if(FirstrBin<0)
    DoSlice(rDim,FirstrBin,LastrBin);
  
  if(LastrBin<FirstrBin) 
    LastrBin = FirstrBin;
  
  Float_t   data[rDim][x1Dim];
  DataSpace memSpace(rank,dataDims);

  Float_t x1Array[x1Dim];
  memset(x1Array,0,sizeof(Float_t)*x1Dim);
   
  Int_t rBins[x1Dim];
  memset(rBins,0,sizeof(Int_t)*x1Dim);
 
  Float_t xmin = ax3lims[0];
  Float_t xmax = ax3lims[1];
  Float_t xbinsize = (xmax-xmin)/x3Dim;
  
  Float_t x2min = ax2lims[0];
  Float_t x2max = ax2lims[1];
  Float_t x2binsize = (x2max-x2min)/x2Dim;

  // This is the max value of r in the range specified.
  Float_t rmax =  x2min + (LastrBin+0.5)*x2binsize - (x2max+x2min)/2.;
  
  string sdata = dataname;
  // Sum the values 
  offset[0] = 0;
  for(Int_t i=FirstrBin;i<=LastrBin;i++) {
    dataSpace.selectHyperslab(H5S_SELECT_SET,count,offset);
    dataSet->read(data,type,memSpace,dataSpace);
    
    Float_t x = xmin + (i+0.5)*xbinsize - (xmax+xmin)/2.;
    
    for(Int_t j=FirstrBin;j<LastrBin;j++) {
      
      // Coordinates of current bin center:
      Float_t y = x2min + (j+0.5)*x2binsize - (x2min+x2max)/2.;
      Float_t r = sqrt(x*x + y*y);

      if(r>rmax) continue; // Out of range.
      // if(r>=axrlims[1]) continue; // out of the limits
      if(y<0) r *= -1;
      
      // Calculates radial bin number
      // Int_t ir = 0;
      // Float_t rstep = (axrlims[1]-axrlims[0])/rDim;
      // for(UInt_t ii=0;ii<rDim;ii++) {
      // 	if(r < axrlims[0] + (ii+1)*rstep) {
      // 	  ir = ii; break;
      // 	}
      // }

      // if(ir<FirstrBin || ir>LastrBin) continue;
      
      for(UInt_t k=0;k<x1Dim;k++) {
	if(sdata.find("charge") != string::npos || sdata.find("p1x1") != string::npos || sdata.find("p2x2") != string::npos)
	  x1Array[k] += -data[j][k];
	else
	  x1Array[k] +=  data[j][k];
	
	rBins[k]++;
      }
    }
    
    offset[0]++;
  }
  dataSet->close();

  // Histogram centering
  Double_t shiftx1 = Shift(opt);
  Double_t x1min = ax1lims[0] - shiftx1;
  Double_t x1max = ax1lims[1] - shiftx1;
  
  TH1F *h1D = new TH1F(); 
  h1D->SetBins(x1Dim,x1min,x1max);
  for(UInt_t k=0;k<x1Dim;k++) {
    Double_t content = x1Array[k];
    if(opt.Contains("avg")) content /= rBins[k];
    else if(opt.Contains("int")) content *= x2binsize * xbinsize;
    h1D->SetBinContent(k+1,content);
  }
  
  root->close();
  
  delete count;
  delete offset;
  delete dataDims;
  delete dataSet;
  delete root;

  return h1D;

}


//_______________________________________________________________________
TH2F* PData::GetH2SliceZX(const char *filename,const char *dataname, Int_t Firstx3Bin, Int_t Lastx3Bin, const char *options) {
  
  // Save memory and time with an specific range selection.

  // Options
  TString opt = options;

  // Open input HDF5 file
  H5File h5 = H5File(filename,H5F_ACC_RDONLY);
  
  // Open main group for data reading
  Group *root = new Group(h5.openGroup("/"));

  // Read data
  // ----------
  DataSet  *dataSet = new DataSet(root->openDataSet(dataname));
  DataSpace dataSpace = dataSet->getSpace();
  const DataType type = dataSet->getDataType();
  
  Int_t rank = dataSpace.getSimpleExtentNdims();
  hsize_t  *dataDims = new hsize_t[rank];
  dataSpace.getSimpleExtentDims(dataDims,NULL);

  // x3 slice
  if(Lastx3Bin>=(Int_t)dataDims[0]) Lastx3Bin = dataDims[0] - 1;
  
  if(Firstx3Bin<0)
    DoSlice(dataDims[0],Firstx3Bin,Lastx3Bin);
  
  if(Lastx3Bin<Firstx3Bin) 
    Lastx3Bin = Firstx3Bin;

  UInt_t x3Dim = Lastx3Bin - Firstx3Bin + 1; 
 
  //----------
  // Fix average problem.
  UInt_t x1AvgFactor = GetNX(0)/dataDims[2];
  UInt_t x1Dim = GetX1N()/x1AvgFactor;
  Double_t x1Min = GetX1Min();
  Double_t x1Max = GetX1Max();

  // cout << Form(" NX1 = %i,   NX1range = %i,   NX1avg = %i,   NX1rangeAvg = %i",
  // 	       GetNX(0), GetX1N(), dataDims[2], x1Dim) << endl;
  
  UInt_t x2AvgFactor = GetNX(1)/dataDims[1];
  UInt_t x2Dim = GetX2N()/x2AvgFactor;
  Double_t x2Min = GetX2Min();
  Double_t x2Max = GetX2Max();
   
  hsize_t  *count  = new hsize_t[rank];
  count[0] = x3Dim;
  count[1] = x2Dim;
  count[2] = x1Dim;

  hsize_t  *offset = new hsize_t[rank];  
  offset[0] = Firstx3Bin;
  offset[1] = GetX2iMin()/x2AvgFactor;
  offset[2] = GetX1iMin()/x1AvgFactor;

  // cout << Form(" x1Dim = %i  offset = %i", count[2], offset[2]) << endl;
  // cout << Form(" x2Dim = %i  offset = %i", count[1], offset[1]) << endl;
  // cout << Form(" x3Dim = %i  offset = %i", count[0], offset[0]) << endl;

  Float_t *data = new Float_t[x3Dim * x2Dim * x1Dim];

  DataSpace memSpace(rank,count);
  dataSpace.selectHyperslab(H5S_SELECT_SET,count,offset);
  dataSet->read(data,PredType::NATIVE_FLOAT,memSpace,dataSpace);
  dataSet->close();

  root->close();
  
  Float_t x2x1Array[x2Dim][x1Dim];
  memset(x2x1Array,0,sizeof(Float_t)*x1Dim*x2Dim);
   
  // Sum the values 
  string sdata = dataname;
  for(UInt_t i=0;i<x3Dim;i++) {
    for(UInt_t j=0;j<x2Dim;j++) {
      for(UInt_t k=0;k<x1Dim;k++) {
	UInt_t index = (long)i*(long)x2Dim*(long)x1Dim + (long)j*(long)x1Dim + (long)k;
	
	if(sdata.find("charge") != string::npos || sdata.find("p1x1") != string::npos || sdata.find("p2x2") != string::npos)
	  x2x1Array[j][k] += -data[index];
	else
	  x2x1Array[j][k] +=  data[index]; 
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

  UInt_t x3AvgFactor = GetNX(2)/dataDims[0];
  Double_t dx3 = GetDX(2) * x3AvgFactor;

  TH2F *h2D = new TH2F(); 
  h2D->SetBins(x1Dim,x1Min,x1Max,x2Dim,x2Min,x2Max); 
  for(UInt_t j=0;j<x2Dim;j++) 
    for(UInt_t k=0;k<x1Dim;k++) { 
      Double_t content = x2x1Array[j][k];
      if(opt.Contains("avg")) content /= x3Dim;
      else if(opt.Contains("int")) content *= dx3;
      h2D->SetBinContent(k+1,j+1,content);
    }
    
  delete count;
  delete offset;
  delete dataSet;
  delete root;
  delete data;
  
  return h2D;

}


//_______________________________________________________________________
TH2F* PData::GetH2SliceZY(const char *filename,const char *dataname, Int_t Firstx2Bin, Int_t Lastx2Bin, const char *options) {
  
  // Save memory and time with an specific range selection.

  // Options
  TString opt = options;

  // Open input HDF5 file
  H5File h5 = H5File(filename,H5F_ACC_RDONLY);
  
  // Open main group for data reading
  Group *root = new Group(h5.openGroup("/"));

  // Read data
  // ----------
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
  UInt_t x1AvgFactor = GetNX(0)/dataDims[2];
  UInt_t x1Dim = GetX1N()/x1AvgFactor;
  Double_t x1Min = GetX1Min();
  Double_t x1Max = GetX1Max();

  UInt_t x3AvgFactor = GetNX(2)/dataDims[0];
  UInt_t x3Dim = GetX3N()/x3AvgFactor;
  Double_t x3Min = GetX3Min();
  Double_t x3Max = GetX3Max();
   
  hsize_t  *count  = new hsize_t[rank];
  count[0] = x3Dim;
  count[1] = x2Dim;
  count[2] = x1Dim;

  hsize_t  *offset = new hsize_t[rank];  
  offset[0] = GetX3iMin()/x3AvgFactor;
  offset[1] = Firstx2Bin;
  offset[2] = GetX1iMin()/x1AvgFactor;

  Float_t *data = new Float_t[x3Dim * x2Dim * x1Dim];

  DataSpace memSpace(rank,count);
  dataSpace.selectHyperslab(H5S_SELECT_SET,count,offset);
  dataSet->read(data,PredType::NATIVE_FLOAT,memSpace,dataSpace);
  dataSet->close();

  root->close();
  
  Float_t x3x1Array[x3Dim][x1Dim];
  memset(x3x1Array,0,sizeof(Float_t)*x1Dim*x3Dim);
   
  // Sum the values 
  string sdata = dataname;
  for(UInt_t i=0;i<x3Dim;i++) {
    for(UInt_t j=0;j<x2Dim;j++) {
      for(UInt_t k=0;k<x1Dim;k++) {
	UInt_t index = (long)i*(long)x2Dim*(long)x1Dim + (long)j*(long)x1Dim + (long)k;
	
	if(sdata.find("charge") != string::npos || sdata.find("p1x1") != string::npos || sdata.find("p2x2") != string::npos)
	  x3x1Array[i][k] += -data[index];
	else
	  x3x1Array[i][k] +=  data[index]; 
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
  Double_t dx2   = GetDX(1) * x2AvgFactor;
  
  TH2F *h2D = new TH2F(); 
  h2D->SetBins(x1Dim,x1Min,x1Max,x3Dim,x3Min,x3Max); 
  for(UInt_t j=0;j<x3Dim;j++) 
    for(UInt_t k=0;k<x1Dim;k++) { 
      Double_t content = x3x1Array[j][k];
      if(opt.Contains("avg")) content /= x2Dim;
      else if(opt.Contains("int")) content *= dx2;
      h2D->SetBinContent(k+1,j+1,content);
    }
    
  delete count;
  delete offset;
  delete dataSet;
  delete root;
  delete data;
  
  return h2D;

}


//_______________________________________________________________________
TH2F* PData::GetH2SliceXY(const char *filename,const char *dataname, Int_t Firstx1Bin, Int_t Lastx1Bin, const char *options) {
  
  // Save memory and time with an specific range selection.

  // Options
  TString opt = options;

  // Open input HDF5 file
  H5File h5 = H5File(filename,H5F_ACC_RDONLY);
  
  // Open main group for data reading
  Group *root = new Group(h5.openGroup("/"));

  // Read data
  // ----------
  DataSet  *dataSet = new DataSet(root->openDataSet(dataname));
  DataSpace dataSpace = dataSet->getSpace();
  const DataType type = dataSet->getDataType();
  
  Int_t rank = dataSpace.getSimpleExtentNdims();
  hsize_t  *dataDims = new hsize_t[rank];
  dataSpace.getSimpleExtentDims(dataDims,NULL);

  // x1 slice
  if(Lastx1Bin>=(Int_t)dataDims[2]) Lastx1Bin = dataDims[2] - 1;
  
  if(Firstx1Bin<0)
    DoSlice(dataDims[2],Firstx1Bin,Lastx1Bin);
  
  if(Lastx1Bin<Firstx1Bin) 
    Lastx1Bin = Firstx1Bin;

  UInt_t x1Dim = Lastx1Bin - Firstx1Bin + 1; 
 
  //----------
  // Fix average problem.
  UInt_t x3AvgFactor = GetNX(2)/dataDims[0];
  UInt_t x3Dim = GetX3N()/x3AvgFactor;
  Double_t x3Min = GetX3Min();
  Double_t x3Max = GetX3Max();

  UInt_t x2AvgFactor = GetNX(1)/dataDims[1];
  UInt_t x2Dim = GetX2N()/x2AvgFactor;
  Double_t x2Min = GetX2Min();
  Double_t x2Max = GetX2Max();
   
  hsize_t  *count  = new hsize_t[rank];
  count[0] = x3Dim;
  count[1] = x2Dim;
  count[2] = x1Dim;

  hsize_t  *offset = new hsize_t[rank];  
  offset[0] = GetX3iMin()/x3AvgFactor;
  offset[1] = GetX2iMin()/x2AvgFactor;
  offset[2] = Firstx1Bin;

  Float_t *data = new Float_t[x3Dim * x2Dim * x1Dim];

  DataSpace memSpace(rank,count);
  dataSpace.selectHyperslab(H5S_SELECT_SET,count,offset);
  dataSet->read(data,PredType::NATIVE_FLOAT,memSpace,dataSpace);
  dataSet->close();

  root->close();
  
  Float_t x3x2Array[x3Dim][x2Dim];
  memset(x3x2Array,0,sizeof(Float_t)*x3Dim*x2Dim);
   
  // Sum the values 
  string sdata = dataname;
  for(UInt_t i=0;i<x3Dim;i++) {
    for(UInt_t j=0;j<x2Dim;j++) {
      for(UInt_t k=0;k<x1Dim;k++) {
	UInt_t index = (long)i*(long)x2Dim*(long)x1Dim + (long)j*(long)x1Dim + (long)k;
	
	if(sdata.find("charge") != string::npos || sdata.find("p1x1") != string::npos || sdata.find("p2x2") != string::npos)
	  x3x2Array[i][j] += -data[index];
	else
	  x3x2Array[i][j] +=  data[index]; 
      }
    }
  }
  
  // Histogram centering
  Double_t shiftx2 = 0;
  Double_t shiftx3 = 0;
  if(opt.Contains("center")) 
    if(!opt.Contains("cyl")) {
      shiftx2 += (x2Min + x2Max)/2.;
      shiftx3 += (x3Min + x3Max)/2.;
    }
  
  x2Min -= shiftx2;
  x2Max -= shiftx2;

  x3Min -= shiftx3;
  x3Max -= shiftx3;

  UInt_t x1AvgFactor = GetNX(0)/dataDims[2];
  Double_t dx1 = GetDX(0) * x1AvgFactor;

  TH2F *h2D = new TH2F(); 
  h2D->SetBins(x2Dim,x2Min,x2Max,x3Dim,x3Min,x3Max); 
  for(UInt_t j=0;j<x3Dim;j++) 
    for(UInt_t k=0;k<x2Dim;k++) { 
      Double_t content = x3x2Array[j][k];
      if(opt.Contains("avg")) content /= x1Dim;
      else if(opt.Contains("int")) content *= dx1;
      h2D->SetBinContent(k+1,j+1,content);
    }
  
  delete count;
  delete offset;
  delete dataSet;
  delete root;
  delete data;
  
  return h2D;

}


//_______________________________________________________________________
TH2F* PData::GetH2ZR(const char *filename,const char *dataname, const char *options) {
  
  // Options
  TString opt = options;

  // Open input HDF5 file
  H5File h5 = H5File(filename,H5F_ACC_RDONLY);
  
  // Open "axis" group 
  Group *axis = new Group(h5.openGroup("/AXIS"));

  // Get axes
  Double_t ax1lims[2];
  DataSet  *ax1 = new DataSet(axis->openDataSet("AXIS1"));
  hsize_t  *ax1Dims = new hsize_t[1];
  ax1->getSpace().getSimpleExtentDims(ax1Dims,NULL);
  DataSpace mem1(1,ax1Dims);
  ax1->read(ax1lims,PredType::NATIVE_DOUBLE,mem1,ax1->getSpace());
  ax1->close(); 
  delete ax1Dims;
  delete ax1;

  Double_t ax2lims[2];
  DataSet  *ax2 = new DataSet(axis->openDataSet("AXIS2"));
  hsize_t  *ax2Dims = new hsize_t[1];
  ax2->getSpace().getSimpleExtentDims(ax2Dims,NULL);
  DataSpace mem2(1,ax2Dims);
  ax2->read(ax2lims,PredType::NATIVE_DOUBLE,mem2,ax2->getSpace());
  ax2->close();
  delete ax2Dims;
  delete ax2;

  Double_t ax3lims[2];
  DataSet  *ax3 = new DataSet(axis->openDataSet("AXIS3"));
  hsize_t  *ax3Dims = new hsize_t[1];
  ax3->getSpace().getSimpleExtentDims(ax3Dims,NULL);
  DataSpace mem3(1,ax3Dims);
  ax3->read(ax3lims,PredType::NATIVE_DOUBLE,mem3,ax3->getSpace());
  ax3->close();
  delete ax3Dims;
  delete ax3;

  axis->close();
  delete axis;

  // Open main group for data reading
  Group *root = new Group(h5.openGroup("/"));

  // Read data
  // ----------
  DataSet    *dataSet = new DataSet(root->openDataSet(dataname));
  DataSpace dataSpace = dataSet->getSpace();
  const DataType type = dataSet->getDataType();

  Int_t        rank = dataSpace.getSimpleExtentNdims();
  hsize_t *dataDims = new hsize_t[rank];
  dataSpace.getSimpleExtentDims(dataDims,NULL);

  UInt_t x1Dim = dataDims[2];
  UInt_t x2Dim = dataDims[1];
  UInt_t x3Dim = dataDims[0];

  // This is for sequential access to the 3D array.
  dataDims[0] = 1; 
  hsize_t  *count  = new hsize_t[rank];
  hsize_t  *offset = new hsize_t[rank];
  for(Int_t j=0;j<rank;j++) {
    offset[j] = 0;
    count[j]  = dataDims[j];
  }
  
  // Define radial bins
  UInt_t rDim = x2Dim; 
  // Define radial axis limits
  Double_t axrlims[2] = {-(ax2lims[1]+ax2lims[0])/2.,(ax2lims[1]+ax2lims[0])/2.};
  
  Float_t   data[rDim][x1Dim];
  DataSpace memSpace(rank,dataDims);

  Float_t rx1Array[rDim][x1Dim];
  memset(rx1Array,0,sizeof(Float_t)*rDim*x1Dim);
   
  Int_t rzBins[rDim][x1Dim];
  memset(rzBins,0,sizeof(Int_t)*rDim*x1Dim);

  string sdata = dataname;
  // Sum the values 
  offset[0] = 0;
  for(UInt_t i=0;i<x3Dim;i++) {
    dataSpace.selectHyperslab(H5S_SELECT_SET,count,offset);
    dataSet->read(data,type,memSpace,dataSpace);

    Double_t x = ax3lims[0] + (i+0.5)*(ax3lims[1]-ax3lims[0])/x3Dim - (ax3lims[1]+ax3lims[0])/2.;
    
    for(UInt_t j=0;j<x2Dim;j++) {
      
      // Coordinates of current bin center:
      Float_t y = ax2lims[0] + (j+0.5)*(ax2lims[1]-ax2lims[0])/x2Dim - (ax2lims[1]+ax2lims[0])/2.;
      Float_t r = sqrt(x*x + y*y);
      
      if(r>=axrlims[1]) continue; // out of the limits
      if(y<0) r *= -1;
      
      // Calculates radial bin number
      Int_t ir = 0;
      Float_t rstep = (axrlims[1]-axrlims[0])/rDim;
      for(UInt_t ii=0;ii<rDim;ii++) {
	if(r < axrlims[0] + (ii+1)*rstep) {
	  ir = ii; break;
	}
      }
      for(UInt_t k=0;k<x1Dim;k++) {
	if(sdata.find("charge") != string::npos || sdata.find("p1x1") != string::npos || sdata.find("p2x2") != string::npos)
	  rx1Array[ir][k] += -data[j][k];
	else
	  rx1Array[ir][k] +=  data[j][k];
	
	rzBins[ir][k]++;
      }
    }
    
    offset[0]++;
  }
  dataSet->close();
  
  
  TH2F *h2D = new TH2F(); 
  h2D->SetBins(x1Dim,ax1lims[0],ax1lims[1],rDim,axrlims[0],axrlims[1]);
  
  for(UInt_t ir=0;ir<rDim;ir++) 
    for(UInt_t k=0;k<x1Dim;k++) {
      if(!opt.Contains("sum"))
	h2D->SetBinContent(k+1,ir+1,rx1Array[ir][k]/rzBins[ir][k]);
      else
       	h2D->SetBinContent(k+1,ir+1,rx1Array[ir][k]);
    }
  
  root->close();
  
  delete count;
  delete offset;
  delete dataDims;
  delete dataSet;
  delete root;

  return h2D;
}


//_______________________________________________________________________
TH3F* PData::GetH3(const char *filename,const char *dataname, const char *options) {

  // Save memory and time with an specific range selection.

  // Options
  TString opt = options;
  
  // Open input HDF5 file
  H5File h5 = H5File(filename,H5F_ACC_RDONLY);
  
  // Open main group for data reading
  Group *root = new Group(h5.openGroup("/"));

  // Read data
  // ----------
  DataSet  *dataSet = new DataSet(root->openDataSet(dataname));
  DataSpace dataSpace = dataSet->getSpace();
  const DataType type = dataSet->getDataType();
  
  Int_t rank = dataSpace.getSimpleExtentNdims();
  hsize_t  *dataDims = new hsize_t[rank];
  dataSpace.getSimpleExtentDims(dataDims,NULL);
  
  //----------
  UInt_t x1AvgFactor = GetNX(0)/dataDims[2];
  UInt_t x1Dim = GetX1N()/x1AvgFactor;
  Double_t x1Min = GetX1Min();
  Double_t x1Max = GetX1Max();
  
  UInt_t x2AvgFactor = GetNX(1)/dataDims[1];
  UInt_t x2Dim = GetX2N()/x2AvgFactor;
  Double_t x2Min = GetX2Min();
  Double_t x2Max = GetX2Max();

  UInt_t x3AvgFactor = GetNX(2)/dataDims[0];
  UInt_t x3Dim = GetX3N()/x3AvgFactor;
  Double_t x3Min = GetX3Min();
  Double_t x3Max = GetX3Max();

  hsize_t  *count  = new hsize_t[rank];
  count[0] = x3Dim;
  count[1] = x2Dim;
  count[2] = x1Dim;
  
  hsize_t  *offset = new hsize_t[rank];  
  offset[0] = GetX3iMin()/x3AvgFactor;
  offset[1] = GetX2iMin()/x2AvgFactor;
  offset[2] = GetX1iMin()/x1AvgFactor;

  Float_t *data = new Float_t[x3Dim * x2Dim * x1Dim];

  DataSpace memSpace(rank,count);
  dataSpace.selectHyperslab(H5S_SELECT_SET,count,offset);
  dataSet->read(data,PredType::NATIVE_FLOAT,memSpace,dataSpace);
  dataSet->close();

  root->close();
  
  // Histogram centering
  Double_t shiftx1 = Shift(opt);  
  Double_t shiftx2(0), shiftx3(0);
  if(opt.Contains("center")) 
    if(!opt.Contains("cyl")) {
      shiftx2 += (x2Min + x2Max)/2.;
      shiftx3 += (x3Min + x3Max)/2.;
    }
  
  x1Min -= shiftx1;
  x1Max -= shiftx1;
  
  x2Min -= shiftx2;
  x2Max -= shiftx2;

  x3Min -= shiftx2;
  x3Max -= shiftx2;
  
  string sdata = dataname;
 
  TH3F *h3D = new TH3F();
  h3D->SetBins(x1Dim,x1Min,x1Max,x2Dim,x2Min,x2Max,x3Dim,x3Min,x3Max); 
  
  for(UInt_t i=0; i<x3Dim; i++) {
    for(UInt_t j=0;j<x2Dim;j++) {
      for(UInt_t k=0;k<x1Dim;k++) {
	UInt_t index = (long)i*(long)x2Dim*(long)x1Dim + (long)j*(long)x1Dim + (long)k;
	if(sdata.find("charge") != string::npos || sdata.find("p1x1") != string::npos || sdata.find("p2x2") != string::npos)
	  h3D->SetBinContent(k+1,j+1,i+1,-data[index]);
	else
	  h3D->SetBinContent(k+1,j+1,i+1,data[index]);
      }
    }
  }
  
  delete count;
  delete offset;
  delete dataDims;
  delete dataSet;
  delete root;

  return h3D;
}


//_______________________________________________________________________
Float_t* PData::Get3Darray(const char *filename,const char *dataname, UInt_t dim[3]) {

  // Open input HDF5 file
  H5File h5 = H5File(filename,H5F_ACC_RDONLY);
  
  // Open main group for data reading
  Group *root = new Group(h5.openGroup("/"));

  // Read data
  // ----------
  DataSet  *dataSet = new DataSet(root->openDataSet(dataname));
  DataSpace dataSpace = dataSet->getSpace();
  const DataType type = dataSet->getDataType();
  
  Int_t rank = dataSpace.getSimpleExtentNdims();
  hsize_t  *dataDims = new hsize_t[rank];
  dataSpace.getSimpleExtentDims(dataDims,NULL);
  
  // Defines the limits of the hyperslab.
  hsize_t  *count  = new hsize_t[rank];
  count[0] = GetX3N();
  count[1] = GetX2N();
  count[2] = GetX1N();
  dim[0] = GetX3N();
  dim[1] = GetX2N();
  dim[2] = GetX1N();

  hsize_t  *offset = new hsize_t[rank];
  offset[0] = GetX3iMin();
  offset[1] = GetX2iMin();
  offset[2] = GetX1iMin();

  // cout << Form(" Hyperslab: x1 (%i: %.2f, %.2f)  x2 (%i: %.2f, %.2f)  x3 (%i: %.2f, %.2f)",x1Dim,x1Min,x1Max,x2Dim,x2Min,x2Max,x3Dim,x3Min,x3Max) << endl;
     
  dataSpace.selectHyperslab(H5S_SELECT_SET,count,offset);  
  
  // It doesn't work I do not know why.
  // data = new Float_t**[x3Dim];
  // for(UInt_t i=0;i<x3Dim;i++) {
  //   data[i] = new Float_t*[x2Dim];
  //   for(UInt_t j=0;j<x2Dim;j++) {
  //     data[i][j] = new Float_t[x1Dim];
  //     for(UInt_t k=0;k<x1Dim;k++) {
  // 	data[i][j][k] = 0.0;
  //     }
  //   }
  // }

  Float_t *data = (Float_t *) malloc(count[0]*count[1]*count[2] * sizeof(Float_t));
  // 
  // Float_t data[x3Dim][x2Dim][x1Dim];
  
  DataSpace memSpace(rank,count);   
  dataSet->read(data,PredType::NATIVE_FLOAT,memSpace,dataSpace);
  dataSet->close();
  
  // if(data) {
  //    for(UInt_t i=0;i<dim[0];i++)
  //      for(UInt_t j=0;j<dim[1];j++)
  // 	 for(UInt_t k=0;k<dim[2];k++)
  // 	   cout << Form("%i, %i, %i : %.2f",i,j,k,DATA(i,j,k)) << endl; 
  // }
  
  root->close();

  delete count;
  delete offset;
  delete dataDims;
  delete dataSet;
  delete root;

  // return data
  return data;

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

//______________________________________________________________________________________
Double_t PData::ShiftT(TString option) {
  TString opt = option;
  
  Double_t shiftt = 0;
  Double_t kp = GetPlasmaK();
  if(opt.Contains("center")) {
    shiftt -= GetPlasmaStart()*kp;
    if(opt.Contains("comov"))      // Centers on the head of the beam.
      shiftt += GetBeamStart()*kp;
  }
  
  return shiftt;
}
