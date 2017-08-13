import scipy.special as sp
import numpy as np
import scipy.integrate as integrate
from scipy.misc import derivative as der

# Wakefields: Electromagnectic fields definition

# Blowout wakefields
class WakeBlowout:
    def __init__(self,Ez0=0,K=0.5,S=0.5,Sk=0.1,rb=5.0,Drho=1.0):
        self.Ez0 = Ez0
        self.K = K
        self.S = S
        self.Sk = Sk
        self.rb = rb
        self.Drho = Drho

    def Ez(self, z, x) :
        # return 0
        z0 = -np.pi/2
        if abs(x) < self.rb :
            return self.S * (z-z0) + self.Ez0
        else :
            return (self.S * (z-z0) + self.Ez0) * np.exp(-(abs(x)-self.rb)/self.Drho)

    def dEz(self, z, x) :
        return self.S
        
    def Ex(self, z, x) :
        if abs(x) < self.rb :
            return (1-self.S) * ( x/2.0 )
        else :
            return (1-self.S) * ( (np.sign(x) * self.rb)/2.0 ) * np.exp(-(np.abs(x)-self.rb)/self.Drho)

        
    def Kx(self, z, x) :
        z0 = -np.pi/2
        return self.K + self.Sk * (z-z0)
        
    def Wx(self, z, x) :
         if abs(x) < self.rb :
             return self.Kx(z,x) * x 
         else :
             return self.Kx(z,x) * np.sign(x) * self.rb * np.exp(-(np.abs(x)-self.rb)/self.Drho)
        
    def dKx(self, z, x) :
        return self.Sk


         

# Linear wakefields (analytic version)
class WakeLinearSimple:
    def __init__(self,nb0=0.1,sz=1.0,sx=1.0):

        # Driver paramerers
        self.nb0 = nb0
        self.sz = sz
        self.sx = sx

    def g(self,z) :
        return np.exp(-z**2/(2*self.sz**2))

    def f(self,x) :
        return np.exp(-x**2/(2*self.sx**2))
  
    def nb(self,z,x) :
        return  self.nb0 * g(z) * f(x)

    def n1(self,z,x) :
        return self.nb0 * np.sqrt(2*np.pi) * self.sz * np.exp(-self.sz**2/2) * self.f(np.abs(x)) * np.sin(z) 

    def dEz(self,z,x) :
        return -self.nb0 * np.sqrt(2*np.pi) * self.sz * np.exp(-self.sz**2/2) * self.f(np.abs(x)) * np.sin(z) 

    def Ez(self,z,x) :
        return self.nb0 * np.sqrt(2*np.pi) * self.sz * np.exp(-self.sz**2/2) * self.f(np.abs(x)) * np.cos(z) 

    def Wx(self,z,x) :
        return -self.nb0 * np.sqrt(2*np.pi) * self.sz * np.exp(-self.sz**2/2) * (x * self.f(np.abs(x)) / self.sx**2) * np.sin(z) 

    def Ex(self,z,x) :
        return 0
               #self.nb0 * np.sqrt(2*np.pi) * self.sz * np.exp(-self.sz**2/2) * np.sin(z) * (self.sx**2 * (1 - self.f(np.abs(x)))/np.abs(x) ) * np.sign(x)

    def Kx(self,z,x) :
        return -self.nb0 * np.sqrt(2*np.pi) * self.sz * np.exp(-self.sz**2/2) * ( (1 - (x**2/self.sx**2)) * self.f(np.abs(x)) / self.sx**2) * np.sin(z) 

    def dKx(self,z,x) :
        return -self.nb0 * np.sqrt(2*np.pi) * self.sz * np.exp(-self.sz**2/2) * ( (1 - (x**2/self.sx**2)) * self.f(np.abs(x)) / self.sx**2) * np.cos(z) 

    def Psi(self,z,x) :
        return -self.nb0 * np.sqrt(2*np.pi) * self.sz * np.exp(-self.sz**2/2) * self.f(np.abs(x)) * np.sin(z) 


class WakeLinearTimon:
    def __init__(self,a0=0.1,w0=1.0):

        # Driver paramerers
        self.a0 = a0
        self.w0 = w0
      
    def Ez(self,z,x) :
        return self.a0**2 * np.sqrt(np.pi/(2*np.e)) * np.cos(z) 

    def Kx(self,z,x) :
        K = np.sqrt(8*np.pi/np.e)*self.a0**2/self.w0**2
        return - K * np.sin(z)
    
    def Wx(self,z,x) :
        return self.Kx(z,x) * x 

    def Ex(self,z,x) :
        return 0



# Linear wakefields
class WakeLinear:
    def __init__(self,nb0=0.1,sz=1.0,sx=1.0):

        # Driver paramerers
        self.nb0 = nb0
        self.sz = sz
        self.sx = sx

       
    def g(self,z) :
        return np.exp(-z**2/(2*self.sz**2))

    def f(self,x) :
        return np.exp(-x**2/(2*self.sx**2))

    def nb(self,z,x) :
        return self.nb0 * self.g(z) * self.f(x)
    
    def G(self,z) :
        def integrand(z,z0) :
            return self.g(z) * np.sin(z-z0)
        
        return integrate.quad(integrand,z,np.inf,args=(z))[0]

    def F(self,x) :
        def integrand(x,x0) :
            if x<=x0 :
                return x * self.f(x) * sp.i0(x)  * sp.k0(x0)
            else :
                return x * self.f(x) * sp.i0(x0) * sp.k0(x)
        
        return integrate.quad(integrand,0,np.inf,args=(x))[0]

    def dG(self,z) :
        def integrand(z,z0) :
            return - self.g(z) * np.cos(z-z0)
        
        return integrate.quad(integrand,z,np.inf,args=(z))[0]

    def ddG(self,z) :
        def integrand(z,z0) :
            return - self.g(z) * np.cos(z-z0)
        
        return integrate.quad(integrand,z,np.inf,args=(z))[0]


    def dF(self,x) :
        def integrand(x,x0) :
            if x<=x0 :
                return - x * self.f(x) * sp.i0(x)  * sp.k1(x0)
            else :
                return   x * self.f(x) * sp.i1(x0) * sp.k0(x)
        
        return integrate.quad(integrand,0,np.inf,args=(x))[0]

    def ddF(self,x) :
        def integrand(x,x0) :
            if x<=x0 :
                return (x/2) * self.f(x) * sp.i0(x) * (sp.k0(x0) + sp.kv(2,x0)) 
            else :
                return (x/2) * self.f(x) * sp.k0(x) * (sp.i0(x0) + sp.iv(2,x0))
            
        
        return - self.f(x) + integrate.quad(integrand,0,np.inf,args=(x))[0]

    def n1(self,z,x) :
        return -self.nb0 * self.G(z) * self.f(x)

    def Psi(self,z,x) :
        return self.nb0 * self.G(z) * self.F(np.abs(x))
    
    def Ez(self,z,x) :
        return -self.nb0 * self.dG(z) * self.F(np.abs(x))

    def Ex(self,z,x) :
        return 0

    def Wx(self,z,x) :
        return -self.nb0 * np.sign(x) * self.G(z) * self.dF(np.abs(x))

    def Kx(self,z,x) :
        return -self.nb0 * self.G(z) * self.ddF(np.abs(x))

    
