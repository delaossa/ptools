
#include <iomanip>
#include <sstream>

#include "PConst.hh"

#include "PUnits.hh"

/*!
 */
PUnits::UnitDefinition::UnitDefinition(const char * name, const char * shortname, Double_t unit)
{
  fName = name;
  fShortName = shortname;
  fUnit = unit;
}

/*!
 */
PUnits::UnitCategory::UnitCategory(const char * name)
{
  fName = name;
  fShortNameMaxLength = 0;
}

/*!
 */
PUnits::UnitCategory::~UnitCategory()
{
  
}

/*!
 */
void PUnits::UnitCategory::AddUnit(const char * name, const char * shortname,
				   Double_t unit)
{
  UnitDefinition * def = 0;

  for (UInt_t i=0;i<fUnits.size();i++) {
    if (fUnits[i]->GetName()==name) {
      def = fUnits[i];
      break;
    }
  }
  
  if (!def) {
    fUnits.push_back(new UnitDefinition(name, shortname, unit));
    if (strlen(shortname)>fShortNameMaxLength) fShortNameMaxLength = strlen(shortname);
  } else {
    def->SetShortName(shortname);
    def->SetUnit(unit);
  }
}

ClassImp(PUnits::UnitsTable);

PUnits::UnitsTable * PUnits::UnitsTable::fgUnitsTable = 0;

/*!
 */
PUnits::UnitsTable::UnitsTable()
{
  // Length
  AddUnit("Femtometer", "fm", "Length", femtometer);
  AddUnit("Nanometer",  "nm", "Length", nanometer);
  AddUnit("Decananometer", "10xnm", "Length", 10*nanometer);
  AddUnit("Hectonanometer", "100 nm", "Length", 100*nanometer);
  AddUnit("Micrometer", "#mum", "Length", micrometer);
  AddUnit("Decamicrometer", "10x#mum", "Length", decamicrometer);
  AddUnit("Millimeter", "mm", "Length", millimeter);
  AddUnit("Centimeter", "cm", "Length", centimeter);
  AddUnit("Meter",      "m",  "Length", meter);
  AddUnit("Kilometer",  "km", "Length", kilometer);

  // Time
  AddUnit("Femtosecond", "fs", "Time", femtosecond);
  AddUnit("Picosecond",  "ps", "Time", picosecond);
  AddUnit("Nanosecond",  "ns", "Time", nanosecond);
  AddUnit("Microsecond", "us", "Time", microsecond);
  AddUnit("Millisecond", "ms", "Time", millisecond);
  AddUnit("Second",      "s",  "Time", second);
  
  // Area
  AddUnit("Nanobarn",    "nbarn", "Area", nanobarn);
  AddUnit("Microbarn",   "ubarn", "Area", microbarn);
  AddUnit("Millibarn",   "mbarn", "Area", millibarn);
  AddUnit("Barn",        "barn",  "Area", barn);
  AddUnit("Millimeter2", "mm^2",  "Area", millimeter2);
  AddUnit("Centimeter2", "cm^2",  "Area", centimeter2);
  AddUnit("Meter2",      "m^2",   "Area", meter2);
  AddUnit("Kilometer2",  "km^2",  "Area", kilometer2);

  // Volume
  AddUnit("Millimeter3", "mm^3",  "Volume", millimeter3);
  AddUnit("Centimeter3", "cm^3",  "Volume", centimeter3);
  AddUnit("Liter",       "l",     "Volume", liter);
  AddUnit("Meter3",      "m^3",   "Volume", meter3);
  AddUnit("Kilometer3",  "km^3",  "Volume", kilometer3);

  // Mass
  AddUnit("Milligram",        "mg",      "Mass", milligram);
  AddUnit("Gram",             "g",       "Mass", gram);
  AddUnit("Kilogram",         "kg",      "Mass", kilogram);
  
  // Energy
  AddUnit("Electronvolt",     "eV",  "Energy", eV);
  AddUnit("Kiloelectronvolt", "keV", "Energy", keV);
  AddUnit("Megaelectronvolt", "MeV", "Energy", MeV);
  AddUnit("Gigaelectronvolt", "GeV", "Energy", GeV);
  AddUnit("Teraelectronvolt", "TeV", "Energy", TeV);
  AddUnit("Joule",            "J",   "Energy", joule);

  // Momentum
  AddUnit("Electronvolt",     "eV/c",  "Momentum", eV);
  AddUnit("Kiloelectronvolt", "keV/c", "Momentum", keV);
  AddUnit("Megaelectronvolt", "MeV/c", "Momentum", MeV);
  AddUnit("Gigaelectronvolt", "GeV/c", "Momentum", GeV);
  AddUnit("Teraelectronvolt", "TeV/c", "Momentum", TeV);

  // Density
  AddUnit("MilligramPerCm3",  "mg/cm^3", "Density", milligram/cm3);
  AddUnit("GramPerCm3",       "g/cm^3",  "Density", g/cm3);
  AddUnit("KilogramPerCm3",   "kg/cm^3", "Density", kg/cm3);

  // Charge
  AddUnit("Electrons",  "e", "Charge", echarge);
  AddUnit("Picocoulomb",  "pC", "Charge", picocoulomb);
  AddUnit("Nanocoulomb",  "nC", "Charge", nanocoulomb);
  AddUnit("Microcoulomb", "uC", "Charge", microcoulomb);
  AddUnit("Millicoulomb", "mC", "Charge", millicoulomb);
  AddUnit("Coulomb",      "C",  "Charge", coulomb);

  // Current
  AddUnit("Milliampere", "mA", "Current", milliampere);
  AddUnit("Centiampere", "cA", "Current", centiampere);
  AddUnit("Deciampere", "dA", "Current", deciampere);
  AddUnit("Ampere",  "A", "Current", ampere);
  AddUnit("Decaampere", "10xA", "Current", decaampere);
  AddUnit("Hectoampere","100xA", "Current", hectoampere);
  AddUnit("Kiloampere", "kA", "Current", kA);
  AddUnit("Decakiloampere", "10xkA", "Current", decakA);

  // Percentages
  AddUnit("Percent", "%", "Percentage", perCent);
  AddUnit("Perten", "10x%", "Percentage", perTen);
  AddUnit("Permille", "0.1x%", "Percentage", perThousand);
  AddUnit("Permillion", "0.01x%", "Percentage", perMillion);

}

/*!
 */
PUnits::UnitsTable::~UnitsTable()
{

}

/*!
 */
PUnits::UnitsTable * PUnits::UnitsTable::Get()
{
  if (!fgUnitsTable) {
    fgUnitsTable = new UnitsTable();
  }
  
  return fgUnitsTable;
}

/*!
 */
void PUnits::UnitsTable::Print(Option_t * option) const
{
  std::cout << "UnitsTable:\n";

  for (UInt_t c=0;c<fCategories.size();c++) {
    std::cout << "  " << fCategories[c]->GetName() << "\n";

    for (UInt_t u=0;u<fCategories[c]->GetUnits().size();u++) {
      std::cout << "    "
		<< std::setw(20) << std::left
		<< fCategories[c]->GetUnits()[u]->GetName()
		<< "  "
		<< std::setw(8) << std::right
		<< fCategories[c]->GetUnits()[u]->GetShortName()
		<< std::endl;
    }
  }
}

/*!
 */
void PUnits::UnitsTable::AddUnit(const char * name, const char * shortname,
				const char * category, Double_t unit)
{
  UnitCategory * cat = 0;

  for (UInt_t i=0;i<fCategories.size();i++) {
    if (fCategories[i]->GetName()==category) {
      cat = fCategories[i];
      break;
    }
  }

  if (!cat) {
    cat = new UnitCategory(category);
    fCategories.push_back(cat);
  }

  cat->AddUnit(name, shortname, unit);
}

ClassImp(PUnits::BestUnit);

/*!
 */
PUnits::BestUnit::BestUnit(const char * category, Int_t precision)
{
  fValue = 0;
  fCategory = 0;

  UnitsTable * table = UnitsTable::Get();

  for (UInt_t i=0;i<table->GetCategories().size();i++) {
    if (table->GetCategories()[i]->GetName()==category) {
      fCategory = table->GetCategories()[i];
      break;
    }
  }

  fPrecision = precision;
}

/*!
 */
PUnits::BestUnit::BestUnit(Double_t value, const char * category, Int_t precision)
{
  fValue = value;
  fCategory = 0;

  UnitsTable * table = UnitsTable::Get();

  for (UInt_t i=0;i<table->GetCategories().size();i++) {
    if (table->GetCategories()[i]->GetName()==category) {
      fCategory = table->GetCategories()[i];
      break;
    }
  }

  fPrecision = precision;
}

/*!
 */
PUnits::BestUnit::~BestUnit()
{

}

/*!
 */
const char * PUnits::BestUnit::AsCString() const
{
  std::ostringstream oss;
  oss << *this;
  return oss.str().c_str();
}

/*!
 */
const TString PUnits::BestUnit::AsString() const
{
  std::ostringstream oss;
  oss << *this;
  return TString(oss.str().c_str());
}

/*!
 */
PUnits::BestUnit::operator const char * () const
{
  std::ostringstream oss;
  oss << *this;
  return oss.str().c_str();
}

/*!
 */
std::ostream & operator<<(std::ostream & os, const PUnits::BestUnit & bu)
{
  Int_t ksup = -1;
  Int_t kinf = -1;
  Double_t umax = 0.;
  Double_t umin = DBL_MAX;
  Double_t rsup = DBL_MAX;
  Double_t rinf = 0.;
  Double_t value, unit, ratio;
  
  PUnits::UnitDefinitions defs = bu.GetCategory()->GetUnits();

  value = bu.GetValue();

  for (UInt_t k=0;k<defs.size();k++) {

    unit = defs[k]->GetUnit();

    if (!(value!=DBL_MAX)) {
      if(unit>umax) {
	umax = unit; 
	ksup = k;
      }
    } else if (value<=DBL_MIN) {
      if(unit<umin) {
	umin = unit; 
	kinf = k;
      }
    } else {
      ratio = value/unit;
      if ((ratio>=1.)&&(ratio<rsup)) {
	rsup = ratio; 
	ksup = k;
      }
      if ((ratio<1.)&&(ratio>rinf)) {
	rinf = ratio;
	kinf = k;
      }
    }
  }

  Int_t idx = ksup;
  if (idx==-1) idx = kinf;
  if (idx==-1) idx = 0;

  os << std::setprecision(bu.GetPrecision())
     << std::fixed
     << std::setw(bu.GetPrecision()+4)
     << value/defs[idx]->GetUnit() 
     << std::resetiosflags(std::ios::fixed) 
     << " " 
     << std::setw(bu.GetCategory()->GetShortNameMaxLength())
     << defs[idx]->GetShortName();

  return os;
}
