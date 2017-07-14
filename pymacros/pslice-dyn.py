#!/usr/bin/env python

import scipy.integrate as integrate
import matplotlib.pyplot as plt
import numpy as np
import os
import types
import collections
from subprocess import call
from ROOT import PData, PGlobals, PPalette, TH2F, TCanvas

class Wake:
	def __init__(self, rb, Delta_rho):
                # Parameters
		self.rb = rb
		self.Delta_rho = Delta_rho
                
	# Force definition Fx(x)
	def Fx(self, *args):
		exp = np.exp
		rb = self.rb
		Delta_rho = self.Delta_rho
		
		if len(args) == 1 and isinstance(args[0], float):
			# case: argument is float
			x = args[0]  
			if x<=(-rb):
				f = exp((x+rb)/Delta_rho)/2
			elif (x>(-rb)) & (x<rb):
				f = -x/2
			elif x>=rb:
				f = -exp(-(x-rb)/Delta_rho)/2
		
		elif len(args) == 1 and isinstance(args[0], collections.Iterable):
				# case: argument is iterable
				x = args[0]
				f = np.piecewise(x, [x<=(-rb), (x>(-rb)) & (x<rb), x>=rb], \
				[lambda x: exp((x+rb)/Delta_rho)/2, lambda x: -x/2, lambda x: -exp(-(x-rb)/Delta_rho)/2])
	
		else:
				# case: wrong arguments
				print('Wrong arguments for Fx(x)')
				
		return f


def dY(Y, t):
	global wake
	return [Y[1], wake.Fx(Y[0])]


def main():

        # Sample definition
	N_t = 10
	N_part = 20000
 	a_t = np.linspace(0, 100, N_t)

        # Distribution function f(x,px)
 	sigma_x = 1.0
 	sigma_p = 0.1
 	mean = [0.5, 0]
 	cov = [[np.power(sigma_x, 2), 0.0], [0.0, np.power(sigma_p, 2)]]
 	X0, P0 = np.random.multivariate_normal(mean, cov, N_part).T

        # Matrix storing the entire solution:
        # (x,px) for each particle and each time. 
 	Yall = np.empty((N_t,2,N_part))
 	
        # Focusing force definition
	rb = 1.0
	Delta_rho = 1.0
 	global wake 
 	wake = Wake(rb, Delta_rho)

        # Plot Fx(x)
        figFx = plt.figure(6)
 	a_x = np.linspace(-rb-9*Delta_rho, rb+9*Delta_rho, 1000)
	plt.plot(a_x, wake.Fx(a_x), color='#1E90FF')
	#plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	plt.xlabel('x (k_p)')
	plt.ylabel('F_x')
	figFx.savefig('figs/Fx.pdf', format='pdf')
	plt.close(figFx)

        
        print('Computing ODEs...')
 	for i in range(0, N_part):
                # Initial values
                x0 = X0[i]
 		p0 = P0[i]
                # Integrate differential equation of motion
 		a_sol = integrate.odeint(dY, [x0, p0], a_t)
 		Yall[:,:,i] = a_sol
	print('done!')

 	
 	if not os.path.exists('./figs'):
		os.makedirs('./figs')
                
        # Plot initital distribution
	finit = plt.figure(1)
	plt.scatter(Yall[0,0,:], Yall[0,1,:], color='#1E90FF', marker='.', alpha=0.8)
 	plt.axis('equal')
 	finit.savefig('figs/finit.pdf', format='pdf')
	plt.close(finit)
	
        # Plot final distribution
	ffin = plt.figure(2)
	plt.scatter(Yall[-1,0,:], Yall[-1,1,:], color='#1E90FF', marker='.', alpha=0.8)
 	plt.axis('equal')
 	ffin.savefig('figs/ffin.pdf', format='pdf')
	plt.close(ffin)

        # Statistics vs time
	cov_all = np.empty((2,2,N_t))
	mean_all = np.empty((2,N_t))

	C = TCanvas('C','Transverse phasespace',800,600)
        C.SetFillStyle(4000);
        PGlobals.SetPlasmaStyle()
       # gStyle.SetNumberContours(255);
        palette = PPalette('electron')
        palette.SetPalette('electron0')

        for i in range(0, N_t):	
		mean_all[:,i] = [sum(Yall[i,0,:])/N_part, sum(Yall[i,1,:])/N_part]
		cov_all[:,:,i] = np.cov(Yall[i,0,:],Yall[i,1,:])
                sx = np.sqrt(cov_all[0,0,i])
                spx = np.sqrt(cov_all[1,1,i])
                hname = 'hpxvsx-%i' % i 
                hpxvsx = TH2F(hname,'',100,-4*sigma_x,4*sigma_x,
                              100,-2*sigma_x,2*sigma_x)
                #                                      mean_all[0,i]-3*sx,
                #                                      mean_all[0,i]+3*sx,
                #                                      100,
                #                                      mean_all[1,i]-3*spx,
                #                                      mean_all[1,i]+3*spx)

                for j in range(0,N_part):
                        hpxvsx.Fill(Yall[i,0,j],Yall[i,1,j])                  

                hpxvsx.Draw('colz')

                if i == 0 :
                        C.Print('figs/hpxvsx.pdf(')
                elif i == N_t-1 :
                        C.Print('figs/hpxvsx.pdf)')
                else :
                        C.Print('figs/hpxvsx.pdf')

                C.Print('figs/hpxvsx-%i.png' % i)
                                      

        

if __name__ == '__main__':
    main()
