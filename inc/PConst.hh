#ifndef PCONST_HH
#define PCONST_HH

//#include <Rtypes.h>
#include <TMath.h>

#include "PUnits.hh"

/*! \namespace PConst
 *  \ingroup common
 *  \brief Some useful constants
 */
namespace PConst
{
  /*! \brief pi */
  static const Double_t pi = 3.14159265358979323846;

  /*! \brief Constant to convert angles from radians to degree */
  static const Double_t kRad2Deg = 180.0/pi;

  /*! \brief Constant to convert angles from degree to radians */
  static const Double_t kDeg2Rad = pi/180.0;

  /*! \brief The speed of light */
  static const Double_t c_light   = 299792458*PUnits::meter/PUnits::second;

  /*! \brief The speed of light squared */
  static const Double_t c_squared = c_light * c_light;

  /*! \brief The electron mass in SI */
  static const Double_t ElectronMass =  9.10938291e-31 * PUnits::kilogram;
  
  /*! \brief The electron mass in energy units */
  static const Double_t ElectronMassE = 0.510998910 * PUnits::MeV;

  /*! \brief The electron charge */
  static const Double_t ElectronCharge = -PUnits::echarge;

  /*! \brief The electron mass in SI */
  static const Double_t ProtonMass =  1.672621777e-27 * PUnits::kilogram;
  
  /*! \brief The electron mass in GeV */
  static const Double_t ProtonMassE = 938.272046 * PUnits::MeV;

  /*! \brief The electron charge */
  static const Double_t ProtonCharge = PUnits::echarge;

  /*! \brief The electron charge squared */
  static const Double_t ESquared = ElectronCharge * ElectronCharge;

  /*! \brief The electron radius */
  static const Double_t ElectronRadius = 2.817940285e-15 * PUnits::meter;

  /*! \brief Vacuum permitivity */
  static const Double_t eps0 = 8.854187817620e-12 * PUnits::farad / PUnits::meter;

  /*! \brief Coulomb constant */
  static const Double_t k0 = 1.0 / (4.0 * pi * eps0) ;

  /*! \brief Plank constant over 2pi */
  //static const Double_t hbar = 1.054571628e-34 * PUnits::joule * PUnits::second;
  static const Double_t hbar = 1.054571726e-34 * PUnits::joule * PUnits::second;

  /*! \brief Electric field strength of first Bohr orbit */
  static const Double_t EF0 = 514.221 * PUnits::GV / PUnits::meter;

  /*! \brief Atomic time interval  */
  static const Double_t AT = 2.41888e-17 * PUnits::second;
 
  /*! \brief Atomic energy scale  */
  static const Double_t XiA = k0 * k0 * ElectronMass * ESquared * ESquared / ( hbar * hbar) ;
  
  /*! \brief Ionization energy of Hydrogen atom */
  static const Double_t EionH = 13.598434005136 * PUnits::eV;
  
  /*! \brief Plasma current */
  static const Double_t I0 = 2.0 * pi * eps0 * ElectronMassE * c_light / fabs(ElectronCharge);
  
  
};


#endif
