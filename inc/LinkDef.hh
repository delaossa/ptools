#ifdef __CINT__

#pragma link off all global;
#pragma link off all class;
#pragma link off all function;

#pragma link C++ nestedclasses;
#pragma link C++ nestedtypedefs;

#pragma link C++ namespace PFunc;

#pragma link C++ namespace PConst;
#pragma link C++ namespace PUnits;
#pragma link C++ class PUnits::UnitsTable-!;
#pragma link C++ class PUnits::BestUnit-!;
#pragma link C++ function operator<<(std::ostream&, const PUnits::BestUnit&);

#pragma link C++ class PData+;
#pragma link C++ global gData;

#pragma link C++ class PDataHiP+;

#pragma link C++ class PPalette+;

#endif
