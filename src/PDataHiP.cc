#include "PDataHiP.hh"

ClassImp(PDataHiP);


//_______________________________________________________________________
PDataHiP::PDataHiP(const char * name, const char * title) : PData(name,title) {
  
}

//_______________________________________________________________________
PDataHiP::PDataHiP(const char * name) : PData(name) {
  
}

//_______________________________________________________________________
PDataHiP::PDataHiP(const char * name, UInt_t time) : PData(name, time) {
  
}

//_______________________________________________________________________
PDataHiP::~PDataHiP()
{

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
    // fgData->SetName(name);
    // fgData->SetTitle(name);
    // TString path = "./";
    // path += name;
    // fgData->SetPath(path.Data());
  }
  
  return dynamic_cast<PDataHiP *>(fgData);
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

  cout << Form(" Loading filenames timestep = %i ",time) << endl;
  
}
