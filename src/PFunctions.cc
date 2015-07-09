#include <iostream>

#include <PConst.hh>

#include "PFunctions.hh"

using namespace PConst;

/*! Plasma frequency
 */
Double_t PFunc::PlasmaFrequency(Double_t n0) {
  return TMath::Sqrt(n0 * ESquared / (eps0 * ElectronMass)) ;
}				   
      
/*! Plasma time depth
 */
Double_t PFunc::PlasmaTimeDepth(Double_t n0) {
  if(n0==0) 
    return 0;
  else
    return 1.0/PlasmaFrequency(n0);
}
     
/*! Plasma time period
 */
Double_t PFunc::PlasmaPeriod(Double_t n0) {
  if(n0==0)
    return 0;
  else
    return TMath::TwoPi()/PlasmaFrequency(n0);
}

/*! Plasma wavenumber
 */
Double_t PFunc::PlasmaWavenumber(Double_t n0) {
  return  PlasmaFrequency(n0) / c_light; 
}

/*! Plasma skindepth
 */
Double_t PFunc::PlasmaSkindepth(Double_t n0) {
  if(n0==0) 
    return 0;
  else
    return   c_light / PlasmaFrequency(n0);
}

/*! Plasma wavelength
 */
Double_t PFunc::PlasmaWavelength(Double_t n0) {
  if(n0==0)
    return 0;
  else
    return  TMath::TwoPi() * c_light / PlasmaFrequency(n0); 
}

/*! Plasma wave breaking field (E0)
 */
Double_t PFunc::PlasmaWBField(Double_t n0) {
  return  fabs(ElectronMass * c_light * PlasmaFrequency(n0) / ElectronCharge);
}


/*! Beam betatron wavenumber
 */
Double_t PFunc::BeamBetatronWavenumber(Double_t gamma, Double_t n0) {
  return PFunc::PlasmaFrequency(n0) / (c_light * sqrt(2) * gamma);
}


// Expressions for long beam self-modulation: See C. Schroeder et al. PRL107,145002(2011)
// --------------------------------------------------------------------------------------

/*! nu constant
 */
Double_t PFunc::Nu(Double_t r0, Double_t n0) {
  Double_t kp = 1.0;
  if(n0!=1.0) kp = PlasmaWavenumber(n0);
  return  4.0 * TMath::BesselI(2,kp*r0) * TMath::BesselK(2,kp*r0);
}

/*! Number of e-foldings (we assume that the beam is an electron beam)
 */
Double_t PFunc::Nefoldings(Double_t z, Double_t zg, Double_t r0, Double_t gamma, Double_t nb, Double_t n0) {
  
  Double_t kp = 1.0;
  if(n0!=1.0) { 
    nb /= n0;
    kp = PlasmaWavenumber(n0);
  }
  Double_t nu = Nu(r0,n0);
  Double_t massratio = ElectronMassE/ElectronMassE;
  Double_t factor = massratio * nu * nb * kp * kp * kp * fabs(zg) * z * z / gamma;
  return (TMath::Power(3.,3./2.)/4.) * TMath::Power(factor,1./3.);
}

/*! Lorenz factor of the phase velocity
 */
Double_t PFunc::PhaseGamma(Double_t z, Double_t zg, Double_t r0, Double_t gamma, Double_t nb, Double_t n0) {
  
  Double_t kp = 1.0;
  if(n0!=1.0) { 
    nb /= n0;
    kp = PlasmaWavenumber(n0);
  }
  Double_t nu = Nu(r0,n0);
  Double_t massratio = ElectronMassE/ElectronMassE;
  return TMath::Power((z/fabs(zg)) * gamma/(massratio * nu * nb),1./6.);
}

/*! Phase velocity from function above
 */
Double_t PFunc::PhaseVelocity(Double_t z, Double_t zg, Double_t r0, Double_t gamma, Double_t nb, Double_t n0) { 
  Double_t x = PhaseGamma(z,zg,r0,gamma,nb,n0);
  return TMath::Sqrt(1. - (1./(x*x)) );
}

/*! Phase velocity  (Pukhov, Kumar et al. PRL107,145003(2011) ) 
 */
Double_t PFunc::PhaseVelocity2(Double_t z, Double_t zg, Double_t gamma, Double_t nb, Double_t n0) {
  
  if(n0!=1.0) { 
    nb /= n0;
  }
  Float_t vb = TMath::Sqrt(1. - (1./(gamma*gamma))); 
  Double_t massratio = ElectronMassE/ElectronMassE;
  return vb * ( 1. - 0.5 * TMath::Power( (fabs(zg)/z) * massratio * nb / (2. * gamma) , 1./3.));
}



// Expressions for ionization probability rate according to ADK theory:
// M.V. Ammosov, N.B. Delone, and V.P. Krainov, Sov. Phys. JETP 64, 1191 (1986).
// --------------------------------------------------------------------------------------

/*! Ionization probability from ADK paper: 
  M.V. Ammosov, N.B. Delone, and V.P. Krainov, Sov. Phys. JETP 64, 1191 (1986).
*/
Double_t PFunc::ADK(Double_t F, Double_t Eion0, Int_t Z, Int_t l, Int_t m) {
  
  // Change to atomic units:
  F      = F/PConst::EF0;
  Eion0  = Eion0/(2*PConst::EionH);
  // ----

  if(F<0.001) return 0; 

  
  Double_t n  = Z/TMath::Sqrt(2*Eion0);
  Double_t F0 = TMath::Power(2*Eion0,3./2.);

  Double_t C = TMath::Power(2*TMath::E()/n,n) / TMath::Sqrt(TMath::TwoPi()*n);
  
  m = TMath::Abs(m);
  Double_t W = Eion0 * (C*C) //* TMath::Sqrt( (3/(TMath::Pi()) * (F/F0))) 
    * ( (2*l+1) * TMath::Factorial(l+m) ) / (TMath::Power(2,m)*TMath::Factorial(m)*TMath::Factorial(l-m))
    * TMath::Power( 2*F0/F , 2*n-m-1 ) 
    * TMath::Exp(-(2./3.)*(F0/F));
  
  return W;  
}

/*!  Ionization probability from original paper (PPT):
  A. M. PERELOMOV, V. S. POPOV, and M. V. TERENT'EV J. Exptl. Theoret. Phys. (U.S.S.R.) 50, 1393-1409 (May, 1966).
*/
Double_t PFunc::PPT(Double_t F, Double_t Eion0, Int_t Z, Int_t l, Int_t m) {
  
  Double_t F0 = TMath::Power(Eion0/PConst::EionH,3./2.) * PConst::EF0;
  Double_t n  = Z/TMath::Sqrt(Eion0/PConst::EionH);

  Double_t C = TMath::Power(2*TMath::E()/n,n) / TMath::Sqrt(TMath::TwoPi()*n);
  
  m = TMath::Abs(m);
  Double_t W = Eion0 * (C*C) //* TMath::Sqrt( (3/(TMath::Pi()) * (F/F0))) 
    * ( (2*l+1) * TMath::Factorial(l+m) ) / (TMath::Power(2,m)*TMath::Factorial(m)*TMath::Factorial(l-m))
    * TMath::Power( 2*F0/F , 2*n-m-1 ) 
    * TMath::Exp(-(2./3.)*(F0/F));
 
  // Return value in atomic units
  return W/(2*PConst::EionH);  
}

/*! Ionization probability from Phys. Plasmas paper (engineering formula) : 
  Bruhwiler D et al 2003 Phys. Plasmas 10 2022.
*/
Double_t PFunc::ADK_ENG(Double_t E, Double_t Eion0, Int_t Z) {
  
  E      = E / (PUnits::GV/PUnits::m);
  Eion0  = Eion0/PUnits::eV;
  Double_t n  = 3.69*Z/TMath::Sqrt(Eion0);
  
  Double_t W = 1.52 * TMath::Power(4,n) * Eion0 / ( n * TMath::Gamma(2*n) )
    * TMath::Power(20.5 * TMath::Power(Eion0,3./2.) / E ,2*n-1)
    * TMath::Exp( -6.83 * TMath::Power(Eion0,3./2.) / E);

  // In terms of the OSIRIS parameters
  // Double_t A =  1.52 * TMath::Power(4,n) * Eion0 / ( n * TMath::Gamma(2*n) )
  //   * TMath::Power(20.5 * TMath::Power(Eion0,3./2.),2*n-1);
  // Double_t B = 6.83 * TMath::Power(Eion0,3./2.);
  // Double_t C = 2*n-1;

  // Double_t WO = A * TMath::Power(E,-C) * TMath::Exp(-B/E);

  // std::cout << Form(" W = %e  WO = %e",W,WO) << std::endl;
  
  return W / PUnits::femtosecond;
}

// Expressions for Blowout regime
// --------------------------------------------------------------------------------------

/*! Blowout radius
    K.V. Lotov Phys. Rev. E 69, 046405 (2004)  
*/
Double_t PFunc::RadiusBO(Double_t lambda, Double_t rmslength) {
  // Both quantities given in normalized units
  //Double_t radiusbo = TMath::Sqrt( 10.0265 * rmslength * TMath::Sqrt(lambda));
  //Double_t radiusbo = 3.8 * TMath::Sqrt( rmslength * TMath::Sqrt(lambda/2.0));
  Double_t radiusbo = 3.17 * TMath::Sqrt( rmslength * TMath::Sqrt(lambda));

  return radiusbo;
  
}


/*! Maximum accelerating gradient
    W. Lu et al. Phys. Plasmas 12, 063101 (2005)
    W. Lu et al. Phys. Plasmas 13, 056709 (2006)
    Narrow beam approximation for \lambda < 1
*/
Double_t PFunc::EzMaxLu(Double_t lambda) {
  // Normalized units
  Double_t EzMax = 1.3 * lambda * TMath::Log(1.0/TMath::Sqrt(lambda/10.0));

  return EzMax;
  
}

/*! Maximum Ez beam
    K.V. Lotov Phys. Rev. E 69, 046405 (2004) 
*/
Double_t PFunc::EzMaxBeamLo(Double_t lambda) {
  // Normalized units
  Double_t EzMax = TMath::Sqrt(lambda/2.0);
  
  return EzMax;
  
}


// Utilities:
// -----------------

/*! Propagating the sign of a root's argument (e.g. for the missing mass case)
 */
Double_t PFunc::SignSqrt(Double_t value)
{
  return value < 0.0 ? -TMath::Sqrt(-value) : TMath::Sqrt(value);
}

/*! The scalar product of the room-coordinates of two TLorentzVectors. c = a * b
 */
Double_t PFunc::ScalarProduct(const TLorentzVector & a, const TLorentzVector & b)
{
  return a.X()*b.X()+a.Y()*b.Y()+a.Z()*b.Z();
}

/*! The cross product of two TLorentzVectors. c = a x b
 */
void PFunc::CrossProduct(TLorentzVector & c,
                         const TLorentzVector & a, const TLorentzVector & b)
{
  c.SetT(0.);
  c.SetX(a.Y()*b.Z()-a.Z()*b.Y());
  c.SetY(a.Z()*b.X()-a.X()*b.Z());
  c.SetZ(a.X()*b.Y()-a.Y()*b.X());
}

/*! Convolution of a Landau with a Gaussian distribution
 *  \param x       The x value where the function is evaluated
 *  \param par     Array of function parameter
\verbatim
 0        Width of the Landau distribution
 1        Most probable value (MPV)
 2        Total area
 3        Width (sigma) of the convoluted Gaussian distribution
\endverbatim
 *
 * In the Landau distribution (represented by the CERNLIB approximation), 
 * the maximum is located at x=-0.22278298 with the location parameter=0.
 * This shift is corrected within this function, so that the actual
 * maximum is identical to the MP parameter. To retrieve the actual peak
 * position one can use the function PFunc::GetLandauGaussParameter()
 */
Double_t PFunc::LandauGaussConvolution(Double_t * x, Double_t * par)
{
  // Numeric constants
  Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  Double_t mpshift  = -0.22278298;       // Landau maximum location
  
  // Control constants
  Double_t np = 100.0;      // number of convolution steps
  Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas
  
  // Variables
  Double_t xx;
  Double_t mpc;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t xlow,xupp;
  Double_t step;
  Double_t i;
  
  // MP shift correction
  mpc = par[1] - mpshift * par[0]; 

  // Range of convolution integral
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];

  step = (xupp-xlow) / np;

  // Convolution integral of Landau and Gaussian by sum
  for(i=1.0; i<=np/2; i++) {
    xx = xlow + (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);

    xx = xupp - (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
  }

  return (par[2] * step * sum * invsq2pi / par[3]);
}

/*! Get the real most probable and full width at half max values for the
 *  Landau-Gauss-Convolution function
 *  \param par       Array of function parameter
\verbatim
 0       Width of the Landau distribution
 1       Most probable value (MPV)
 2       Total area
 3       Width (sigma) of the convoluted Gaussian distribution
\endverbatim
 *  \param xmax      Position of the peak (the most probable value)
 *  \param HWHMlow   Half width at half maximum to the left of the peak
 *  \param HWHMup    Half width at half maximum to the right of the peak
 */
Bool_t PFunc::GetLandauGaussParameter(Double_t * par,
				      Double_t & xmax, Double_t & HWHMlow, Double_t & HWHMup)
{
  Double_t p,x,fy,fxr,fxl;
  Double_t step;
  Double_t l,lold;
  Int_t i = 0;
  Int_t MAXCALLS = 10000;

  // Search for maximum
  
  p = par[1] - 0.5 * par[0];
  step = 0.05 * par[0];
  lold = -2.0;
  l    = -1.0;

  while ( (l != lold) && (i < MAXCALLS) ) {
    i++;
    
    lold = l;
    x = p + step;
    l = LandauGaussConvolution(&x,par);
    
    if (l < lold)
      step = -step/10;
    
    p += step;
  }
  
  if (i==MAXCALLS)
    return kFALSE;
  
  xmax = x;
  
  fy = l/2;
  
  // Search for right x location of fy
  
  p = xmax + 10.0 * par[0];
  step = 4.0 * par[0];
  lold = -2.0;
  l    = -1e300;
  i    = 0;
  
  while ((l!=lold) && (i<MAXCALLS)) {
    i++;
    
    lold = l;
    x = p + step;
    l = TMath::Abs(LandauGaussConvolution(&x,par) - fy);
    
    //std::cout << xmax << "\t" << x << "\t" << fy << "\t" << l << std::endl;

    if (l > lold)
      step = -step/10;
    
    p += step;
  }
  
  if (i == MAXCALLS)
    return kFALSE;
  
  fxr = x;
  
  // Search for left x location of fy

  p = xmax - 10.0 * par[0];
  step = -4.0 * par[0];
  lold = -2.0;
  l    = -1e300;
  i    = 0;
  
  while ( (l != lold) && (i < MAXCALLS) ) {
    i++;
    
    lold = l;
    x = p + step;
    l = TMath::Abs(LandauGaussConvolution(&x,par) - fy);
    
    if (l > lold)
      step = -step/10;
    
    p += step;
  }
  
  if (i==MAXCALLS)
    return kTRUE;

  fxl = x;

  HWHMlow = xmax - fxl;
  HWHMup = fxr - xmax;

  return kTRUE;
}

/*! An asymmetric Gaussian distribution
 *  \param x       The x value where the function is evaluated
 *  \param par     Array of function parameter
\verbatim
 0       Constant
 1       Peak position
 2       Sigma to the left of the peak
 3       Sigma to the right of the peak
\endverbatim
 */
Double_t PFunc::AsymGaussian(Double_t * x, Double_t * par)
{
  Double_t value;
  
  if (x[0]<par[1]) {
    value = par[0]*TMath::Gaus(x[0], par[1], par[2]);
  } else {
    value = par[0]*TMath::Gaus(x[0], par[1], par[3]);
  }
  
  return value;
}

/*! The Weibull distribution
 *
\verbatim
    f(x) = k/lambda (x/lambda)^(k-1) e^((x/lambda)^k)      [x>=0]
    f(x) = 0                                               [x<0]
\endverbatim
 *  \param x       The x value where the function is evaluated
 *  \param par     Parameter array
\verbatim
 0       k
 1       lambda
\endverbatim
 */
Double_t PFunc::Weibull(Double_t * x, Double_t * par)
{
  if (x[0]<0) return 0;

  Double_t k = par[0];
  Double_t lambda = par[1];
  Double_t xn = x[0]/lambda;

  return k/lambda * TMath::Power(xn, k) * TMath::Exp(-TMath::Power(xn, k));
}

/*! A triangular distribution
 *
 *  \param x       The x value where the function is evaluated
 *  \param par     Array of function parameter
\verbatim
 0       x position of the maximum
 1       width to the left of the maximum
 2       width to the right of the maximum
\endverbatim
 */
Double_t PFunc::Triangular(Double_t * x, Double_t * par)
{
  if (par[1]<=0) par[1] = 1.0;
  if (par[2]<=0) par[2] = 1.0;

  if (x[0]<par[0]-par[1]) return 0.0;
  if (x[0]>par[0]+par[2]) return 0.0;

  if (x[0]<=par[0]) return (x[0]-(par[0]-par[1]))/par[1];
  if (x[0]>par[0]) return ((par[0]+par[2])-x[0])/par[2];

  return 0.0;
}
