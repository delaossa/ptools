#ifndef _PUNITS_HH_
#define _PUNITS_HH_

#include <iostream>
//#include <fstream>
#include <sstream>
//#include <iomanip>
//#include <stdlib.h>

#include <vector>

#include <Rtypes.h>
#include <TString.h>

/*! \namespace PUnits
 *  \ingroup common
 *  \brief Helper classes and functions to deal with units
 *
 * The namespace PUnits offers classes and constants that allow to deal
 * with units. For the common user it should be enough to use the class
 * PUnits::BestUnit that will give the best unit for a given value and
 * a unit category. To use this class one has to initialize the units
 * table via a call to PUnits::UnitsTable::Get().
 *
 * The following is an example ROOT session demonstrating the use of
 * the namespace PUnits:
\verbatim
 root [0] PUnits::UnitsTable::Get()
 (class PUnits::UnitsTable*)0x8312ea0
 root [1] PUnits::UnitsTable::Get()->Print()
 UnitsTable:
   Length
     Femtometer                  fm
     Nanometer                   nm
     Micrometer                  um
     Millimeter                  mm
     Centimeter                  cm
     Meter                        m
     Kilometer                   km
   Area
     Nanobarn                 nbarn
     Microbarn                ubarn
     Millibarn                mbarn
     Barn                      barn
     Millimeter2               mm^2
 ...
 root [2] Double_t length = 0.01234 * PUnits::millimeter;
 root [3] length
 (Double_t)1.23400000000000007e-03
 root [4] cout << PUnits::BestUnit(length, "Length") << endl;
  12.340000 um
 root [5] cout << PUnits::BestUnit(1234.567 * PUnits::millimeter2, "Area") << endl;
  12.345670  cm^2
\endverbatim
 *
 */
namespace PUnits
{
  class UnitCategory;
  class UnitsTable;

  /*! \class UnitDefinition
   *  \ingroup common
   *  \brief Definition of a physical unit
   */
  class UnitDefinition
  {
    friend class UnitCategory;

  public:

    /*! \brief The destructor
     */
    virtual ~UnitDefinition() {};

    /*! \brief Get name of the unit
     */
    const TString &           GetName() { return fName; }

    /*! \brief Get short name of the unit
     */
    const TString &           GetShortName() { return fShortName; }

    /*! \brief Get the value of the unit
     */
    Double_t                  GetUnit() { return fUnit; }

  protected:

    /*! \brief The constructor
     */
    UnitDefinition(const char * name, const char * shortname, Double_t unit);

    /*! \brief Set short name of the unit
     */    
    void                      SetShortName(const char * n) { fShortName = n; }

    /*! \brief Set the value of the unit
     */
    void                      SetUnit(Double_t unit) { fUnit = unit; }

    TString                   fName;        ///< The name of the unit
    TString                   fShortName;   ///< The short name of the unit
    Double_t                  fUnit;        ///< The value of the unit
  };

  /*! \var std::vector<UnitDefinition*> UnitDefinitions;
   * \brief A vector of UnitDefinition objects
   */
  typedef std::vector<UnitDefinition*> UnitDefinitions;

  /*! \class UnitCategory
   *  \ingroup common
   *  \brief Definition of a category of a physical unit
   */
  class UnitCategory
  {
    friend class UnitsTable;

  public:

    /*! \brief The destructor
     */
    virtual ~UnitCategory();

    /*! \brief Get the name of the unit category
     */
    const TString &           GetName() { return fName; }

    /*! \brief Get vector of unit definitions (UnitDefinition)
     */
    const UnitDefinitions &   GetUnits() { return fUnits; }

    /*! \brief Get vector of unit definitions (UnitDefinition)
     */
    const UnitDefinitions &   GetUnits() const { return fUnits; }

    /*! \brief Get the maximum length of all short names in this category
     */
    Int_t                     GetShortNameMaxLength() const { return fShortNameMaxLength; }

  protected:

    /*! \brief The constructor
     */
    UnitCategory(const char * name);

    /*! \brief Add a new unit to this unit category
     */
    void                      AddUnit(const char * name, const char * shortname,
				      Double_t unit);
    
    TString                   fName;               ///< The name of the unit category
    UnitDefinitions           fUnits;              ///< The vector of units in this category
    UInt_t                    fShortNameMaxLength; ///< The maximum length of all short names in this category
  };
  
  /*! \var std::vector<UnitCategories*> UnitCategories;
   * \brief A vector of UnitCategory objects
   */
  typedef std::vector<UnitCategory*> UnitCategories;
  
  /*! \class UnitsTable
   *  \ingroup common
   *  \brief Table of units
   */
  class UnitsTable
  {
  public:

    /*! \brief Returns the pointer to the singleton unit table object
     */
    static UnitsTable *       Get();

    /*! \brief The destructor
     */
    virtual ~UnitsTable();

    /*! \brief Returns a vector with all unit categories (UnitCategory)
     */
    UnitCategories &          GetCategories() { return fCategories; }

    /*! \brief Add a new unit to the table
     */
    void                      AddUnit(const char * name, const char * shortname,
				      const char * category, Double_t unit);
    
    /*! \brief Print a list of known units
     */
    void                      Print(Option_t * option = "") const;

  protected:

    UnitCategories            fCategories;   ///< The vector of unit categories

    /*! \brief The constructor
     */
    UnitsTable();
    static UnitsTable        *fgUnitsTable;  ///< The pointer to the singleton object
  };
  
  class BestUnit
  {
  public:

    /*! \brief The constructor taking a category name
     */
    BestUnit(const char * category, Int_t precision = 6);
    
    /*! \brief The constructor taking a category name and value
     */
    BestUnit(Double_t value, const char * category, Int_t precision = 6);

    /*! \brief The destructor
     */
    ~BestUnit();

    /*! \brief Set the value
     */
    void                      SetValue(Double_t v) { fValue = v; }

    /*! \brief Assignment operator
     */
    void                      operator=(Double_t v) { fValue = v; }

    /*! \brief Set the floating point precision for output
     */
    void                      SetPrecision(Int_t p) { fPrecision = p; }

    /*! \brief Get the value
     */
    Double_t                  GetValue() const { return fValue; }

    /*! \brief Get a pointer to the unit category
     */
    const UnitCategory *      GetCategory() const { return fCategory; }

    /*! \brief Get the floating point precision for output
     */
    Int_t                     GetPrecision() const { return fPrecision; }

    /*! \brief Get the value with best unit as a c string
     */
    const char *              AsCString() const;

    /*! \brief Get the value with best unit as a TString
     */
    const TString             AsString() const;

    /*! \brief Get the value with best unit as a c string
     */
    operator const char * () const;

    /*! \brief Get best unit and its name
     */
    void GetBestUnits(Double_t &unit, std::string &sunit) ;


  protected:

    Double_t                  fValue;      ///< The actual value
    Int_t                     fPrecision;  ///< The precission
    UnitCategory             *fCategory;   ///< The unit category
  };

  /** @name Electric charge [Q]
   */
  //@{
  static const Double_t coulomb = 1;                        ///< Coulomb
  static const Double_t millicoulomb = 1.e-3*coulomb;       ///< Millicoulomb
  static const Double_t microcoulomb = 1.e-3*millicoulomb;  ///< Microcoulomb
  static const Double_t nanocoulomb = 1.e-3*microcoulomb;  ///< Nanocoulomb
  static const Double_t picocoulomb = 1.e-3*nanocoulomb;   ///< Picocoulomb
  static const Double_t femptocoulomb = 1.e-3*picocoulomb;  ///< Femptocoulomb
  static const Double_t echarge  = 1.602176565e-19*coulomb; ///< Positron charge

  //@}

  /** @name Length [L], Area [L*L] and Volume [L*L*L]
   */
  //@{

  static const Double_t m           = 1.; ///< Meter
  static const Double_t m2          = m*m; ///< Squaremeter
  static const Double_t m3          = m*m*m; ///< Cubicmeter
  static const Double_t meter       = 1.; ///< Meter
  static const Double_t meter2      = meter*meter; ///< Squaremeter
  static const Double_t meter3      = meter*meter*meter; ///< Cubicmeter
  static const Double_t liter       = 1.e-3*meter3; ///< Liter

  static const Double_t cm          = 1e-2*m; ///< Centimeter
  static const Double_t cm2         = cm*cm; ///< Squarecentimeter
  static const Double_t cm3         = cm*cm*cm; ///< Cubiccentimeter
  static const Double_t centimeter  = 1e-2*meter; ///< Centimeter
  static const Double_t centimeter2 = centimeter*centimeter; ///< Squarecentimeter
  static const Double_t centimeter3 = centimeter*centimeter*centimeter; ///< Cubiccentimeter

  static const Double_t mm          = 1e-1*cm; ///< Millimeter
  static const Double_t mm2         = mm*mm; ///< Squaremillimeter
  static const Double_t mm3         = mm*mm*mm; ///< Cubicillimeter
  static const Double_t millimeter  = 1e-1*centimeter; ///< Millimeter
  static const Double_t millimeter2 = millimeter*millimeter; ///< Squaremillimeter
  static const Double_t millimeter3 = millimeter*millimeter*millimeter; ///< Cubicmillimeter
  
  static const Double_t kilometer   = 1000.*meter; ///< Kilometer
  static const Double_t kilometer2  = kilometer*kilometer; ///< Squarekilometer
  static const Double_t kilometer3  = kilometer*kilometer*kilometer; ///< Cubickilometer
  static const Double_t      barn   = 1.e-28*meter2; ///< Barn
  static const Double_t millibarn   = 1.e-3*barn; ///< Millibarn
  static const Double_t microbarn   = 1.e-6*barn; ///< Microbarn
  static const Double_t  nanobarn   = 1.e-9*barn; ///< Nanobarn

  static const Double_t um          = 1.e-6*meter; ///< Micrometer
  static const Double_t micrometer  = 1.e-6*meter; ///< Micrometer
  static const Double_t nm          = 1.e-9*meter; ///< Nanometer
  static const Double_t nanometer   = 1.e-9*meter; ///< Nanometer
  static const Double_t femtometer  = 1.e-15*meter; ///< Femtometer
  static const Double_t fermi       = 1*femtometer; ///< Femtometer
  static const Double_t inch        = 2.54*centimeter; ///< Inch
  static const Double_t Angstrom    = 0.1*nanometer; ///< Angstrom
  //@}
  
  /** @name Time [T]
   */
  //@{
  static const Double_t      second = 1; ///< Second
  static const Double_t millisecond = 1.e-3*second; ///< Millisecond
  static const Double_t microsecond = 1.e-3*millisecond; ///< Microsecond
  static const Double_t  nanosecond = 1.e-3*microsecond; ///< Nanosecond
  static const Double_t  picosecond = 1.e-3*nanosecond; ///< Picosecond
  static const Double_t femtosecond = 1.e-3*picosecond; ///< Femtosecond
  static const Double_t  atomictime = 0.02419*femtosecond; ///< Atomictime
  //@}
  
  /** @name Frequency [T^-1]
   */
  //@{
  static const Double_t     hertz   = 1./second; ///< Hertz
  static const Double_t kilohertz   = 1.e+3*hertz; ///< Kilohertz
  static const Double_t Megahertz   = 1.e+6*hertz; ///< Megahertz
  
  static const Double_t  Hz         = 1*hertz; ///< Hertz
  static const Double_t kHz         = 1*kilohertz; ///< Kilohertz
  static const Double_t MHz         = 1*Megahertz; ///< Megahertz
  //@}
  
  /** @name Energy [E]
   */
  //@{
  static const Double_t joule   = 1.0;  ///< Joule
  static const Double_t  eV     = echarge*joule; ///< eV
  static const Double_t Gigaelectronvolt = 1.e9*eV; ///< GeV
  static const Double_t Megaelectronvolt = 1.e-3*Gigaelectronvolt; ///< MeV
  static const Double_t     electronvolt = 1.e-6*Megaelectronvolt; ///< eV
  static const Double_t Kiloelectronvolt = 1.e+3*electronvolt; ///< keV
  static const Double_t Teraelectronvolt = 1.e+3*Gigaelectronvolt; ///< TeV
  static const Double_t MeV     = Megaelectronvolt; ///< MeV
 
  static const Double_t keV     = Kiloelectronvolt; ///< keV
  static const Double_t GeV     = Gigaelectronvolt; ///< GeV
  static const Double_t TeV     = Teraelectronvolt; ///< TeV

  static const Double_t EA      = 27.211*eV;   ///< EA
  //@}

  /** @name current [Q][T^-1]
   */
  //@{
  static const double      ampere = coulomb/second; ///< Ampere
  static const double microampere = 1.e-6*ampere; ///< Microampere
  static const double milliampere = 1.e-3*ampere; ///< Milliampere
  static const double centiampere = 1.e-2*ampere; ///< Centiampere
  static const double  deciampere = 1.e-1*ampere; ///< Deciampere
  static const double  nanoampere = 1.e-9*ampere; ///< Nanoampere
  static const double  decaampere = 1.e+1*ampere; ///< Decaampere
  static const double hectoampere = 1.e+2*ampere; ///< Hectoampere
  static const double  kiloampere = 1.e+3*ampere; ///< Kiloampere
  static const double          kA = kiloampere;   ///< kA
  static const double      decakA = 1.e+1*kiloampere;   ///< DecakA
  //@}

  /** @name Electric potential [E][Q^-1]
   */
  //@{
  static const double     volt = joule/coulomb; ///< Volt
  static const double teravolt = 1e12*volt; ///< Tegavolt
  static const double gigavolt = 1e9*volt; ///< Megavolt
  static const double megavolt = 1e6*volt; ///< Megavolt
  static const double kilovolt = 1e3*volt; ///< Kilovolt
  static const Double_t MV     = megavolt; ///< MV
  static const Double_t kV     = kilovolt; ///< kV
  static const Double_t GV     = gigavolt; ///< GV
  static const Double_t TV     = teravolt; ///< TV
  //@}
  
  /** @name Electric resistance [E][T][Q^-2]
   */
  //@{
  static const double ohm = volt/ampere; ///< Ohm
  //@}

  /** @name Electric capacitance [Q^2][E^-1]
   */
  //@{
  static const double farad = coulomb/volt; ///< Farad
  static const double millifarad = 1.e-3*farad; ///< Millifarad
  static const double microfarad = 1.e-6*farad; ///< Microfarad
  static const double  nanofarad = 1.e-9*farad; ///< Nanofarad
  static const double  picofarad = 1.e-12*farad; ///< Picofarad
  //@}

  /** @name Magnetic Flux [T][E][Q^-1]
   */
  //@{
  static const double weber = volt*second; ///< Weber
  //@}

  /** @name Magnetic Field [T][E][Q^-1][L^-2]
   */
  //@{
  static const double tesla     = volt*second/meter2; ///< Tesla
  static const double gauss     = 1.e-4*tesla; ///< Gauss
  static const double kilogauss = 1.e-1*tesla; ///< Kilogauss
  //@}

  /** @name Inductance [T^2][E][Q^-2]
   */
  //@{
  static const double henry = weber/ampere; ///< Henry
  //@}
  
  /** @name Mass [E][T^2][L^-2]
   */
  //@{
  static const Double_t  kilogram = joule*second*second/(meter*meter); ///< Kilogram
  static const Double_t        kg = joule*second*second/(meter*meter); ///< Kilogram
  static const Double_t      gram = 1.e-3*kilogram; ///< Gram
  static const Double_t      g    = 1.e-3*kilogram; ///< Gram
  static const Double_t milligram = 1.e-3*gram; ///< Milligram
  //@}

  /** @name Power [E][T^-1]
   */
  //@{
  static const Double_t watt    = joule/second; ///< Watt
  //@}
  
  /** @name Force [E][L^-1]
   */
  //@{
  static const Double_t newton  = joule/meter; ///< Newton
  //@}

  /** @name Pressure [E][L^-3]
   */
  //@{
  static const Double_t pascal     = newton/meter2; ///< Pascal
  static const Double_t bar        = 100000*pascal; ///< Bar
  static const Double_t atmosphere = 101325*pascal; ///< Atmosphere
  //@}

  /** @name Temperature
   */
  //@{
  static const Double_t kelvin = 1.; ///< Kelvin
  //@}

  /** @name Amount of substance
   */
  //@{
  static const Double_t mole = 1.; ///< Mole
  //@}
  
  /** @name Activity [T^-1]
   */
  //@{
  static const Double_t becquerel = 1./second; ///< Becquerel
  static const Double_t curie = 3.7e+10 * becquerel; ///< Curie
  //@}
  
  /** @name Absorbed dose [L^2][T^-2]
   */
  //@{
  static const Double_t gray = joule/kilogram; ///< Gray
  //@}

  /** @name Miscellaneous constants
   */
  //@{
  static const Double_t perTen      = 0.1;    ///< Perten
  static const Double_t perCent     = 0.01;   ///< Percent
  static const Double_t perThousand = 0.001;  ///< Permile
  static const Double_t perMillion  = 0.0001; ///< PPB
  //@}
};

std::ostream & operator<<(std::ostream & os, const PUnits::BestUnit & bu);

#endif
