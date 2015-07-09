#ifndef PFUNCTIONS_HH
#define PFUNCTIONS_HH

#include <TLorentzVector.h>

#include "PConst.hh"

/*! \namespace PFunc
 *  \ingroup common
 *  \brief Some useful functions.
 */
namespace PFunc
{

  /** @name Characteristic functions depending on the plasma density
   */
  //@{
  
  /*! \brief Plasma frequency
   */
  Double_t                    PlasmaFrequency(Double_t n0); 

  /*! \brief Plasma time depth: inverse of the frequency
   */
  Double_t                    PlasmaTimeDepth(Double_t n0) ; 

  /*! \brief Plasma time period
   */
  Double_t                    PlasmaPeriod(Double_t n0); 

  /*! \brief Plasma wave number
   */
  Double_t                    PlasmaWavenumber(Double_t n0); 
 
  /*! \brief Plasma skin depth
   */
  Double_t                    PlasmaSkindepth(Double_t n0); 
  
  /*! \brief Plasma wave length
   */
  Double_t                    PlasmaWavelength(Double_t n0); 

  /*! \brief Plasma wave breaking field (E0)
   */
  Double_t                    PlasmaWBField(Double_t n0); 

  //@}

  /** @name Characteristic functions of the beam
   */
  //@{
  
  /*! \brief Beam betatron wavenumber
   */
  Double_t                    BeamBetatronWavenumber(Double_t gamma, Double_t n0 = 1); 

  //@}


  /** @name Expressions for long beam self-modulation: See C. Schroeder et al. PRL107,145002(2011)
   */
  //@{
  
  /*! \brief nu constant
   */
  Double_t                   Nu(Double_t r0, Double_t n0 = 1.0);

  /*! \brief Number of e-foldings
   */
  Double_t                   Nefoldings(Double_t z, Double_t zg, Double_t r0, Double_t gamma, Double_t nb, Double_t n0 = 1.0);
  
  /*! \brief Lorenz factor of the phase velocity
   */
  Double_t                   PhaseGamma(Double_t z, Double_t zg, Double_t r0, Double_t gamma, Double_t nb, Double_t n0 = 1.0);
  
  /*! \brief Phase velocity from function above
   */
  Double_t                   PhaseVelocity(Double_t z, Double_t zg, Double_t r0, Double_t gamma, Double_t nb, Double_t n0 = 1.0) ;
  
  
  /*! \brief Phase velocity (Pukhov, Kumar et al. PRL107,145003(2011) ) 
   */
  Double_t                   PhaseVelocity2(Double_t z, Double_t zg, Double_t gamma, Double_t nb, Double_t n0 = 1.0);
  
  
  //@}


  /** @name Expressions for ionization probability rate:
   */
  //@{
  
  /*! \brief Ionization probability from ADK paper: 
    M.V. Ammosov, N.B. Delone, and V.P. Krainov, Sov. Phys. JETP 64, 1191 (1986).
    
  */
  Double_t                   ADK(Double_t F, Double_t Eion0, Int_t Z=1, Int_t l=0, Int_t m=0);
  
  /*! \brief Ionization probability from original paper (PPT):
    A. M. PERELOMOV, V. S. POPOV, and M. V. TERENT'EV J. Exptl. Theoret. Phys. (U.S.S.R.) 50, 1393-1409 (May, 1966).
  */
  Double_t                   PPT(Double_t F, Double_t Eion0, Int_t Z=1, Int_t l=0, Int_t m=0);
  
  /*! \brief Ionization probability from Phys. Plasmas paper (engineering formula) : 
    Bruhwiler D et al 2003 Phys. Plasmas 10 2022.
  */
  Double_t                   ADK_ENG(Double_t F, Double_t Eion0, Int_t Z=1);
  
  //@}


  /** @name Expressions for Blowout regime
   */
  //@{
  
  /*! \brief Blowout radius
    K.V. Lotov Phys. Rev. E 69, 046405 (2004)    
  */
  Double_t                   RadiusBO(Double_t curr, Double_t rmslength);
  
  /*! \brief Maximum accelerating gradient
    W. Lu et al. Phys. Plasmas 12, 063101 (2005)
    W. Lu et al. Phys. Plasmas 13, 056709 (2006)
    Narrow beam approximation for \lambda < 1
  */
  Double_t                   EzMaxLu(Double_t lambda);
  
  /*! \brief Maximum Ez beam
    K.V. Lotov Phys. Rev. E 69, 046405 (2004) 
  */
  Double_t                   EzMaxBeamLo(Double_t lambda);
    

  //@}


  /** @name More technical functions including special fitting functions.
   */
  //@{

  /*! \brief The square root of value with sign
   */
  Double_t                    SignSqrt(Double_t value);

  /*! \brief The scalar product of two TLorentzVectors
   */
  Double_t                    ScalarProduct(const TLorentzVector & a, const TLorentzVector & b);

  /*! \brief The cross product of two TLorentzVectors
   */
  void                        CrossProduct(TLorentzVector & c,
                                           const TLorentzVector & a, const TLorentzVector & b);

  /*! \brief Convolution of a Landau with a Gaussian distribution
   */
  Double_t                    LandauGaussConvolution(Double_t * x, Double_t * par);

  /*! \brief Get the real most probable value together with the full width at half max values
   *         for the Landau-Gauss-Convolution function
   */
  Bool_t                      GetLandauGaussParameter(Double_t * par,
						      Double_t &xmax, Double_t &HWHMlow, Double_t &HWHMup);

  /*! \brief An asymmetric Gaussian distribution
   */
  Double_t                    AsymGaussian(Double_t * x, Double_t * par);

  /*! \brief The Weibull distribution
   */
  Double_t                    Weibull(Double_t * x, Double_t * par);

  /*! \brief A triangular distribution
   */
  Double_t                    Triangular(Double_t * x, Double_t * par);
  
  //@}
}

#endif
