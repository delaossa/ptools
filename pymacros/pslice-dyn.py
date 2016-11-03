#!/usr/bin/env python

import scipy.integrate as integrate
import matplotlib.pyplot as plt
import numpy as np
import os
import types
import collections
from subprocess import call

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
	N_t = 1000
	N_part = 5000
 	a_t = np.linspace(0, 100, N_t)

        # Distribution function f(x,px)
 	sigma_x = 1.0
 	sigma_p = 0.1
 	mean = [0.2, 0]
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

	for i in range(0, N_t):	
		mean_all[:,i] = [sum(Yall[i,0,:])/N_part, sum(Yall[i,1,:])/N_part]
		cov_all[:,:,i] = np.cov(Yall[i,0,:],Yall[i,1,:])
	

 	fmeanx = plt.figure(4)	
 	plt.plot(a_t, mean_all[0,:].T, color='#1E90FF')	
 	fmeanx.savefig('figs/meanx.pdf', format='pdf')	
 	plt.close(fmeanx)

 	fmeanpx = plt.figure(4)	
 	plt.plot(a_t, mean_all[1,:].T, color='#1E90FF')	
 	fmeanpx.savefig('figs/meanpx.pdf', format='pdf')	
 	plt.close(fmeanpx)

        # Randomly pick a subset of particles to track their trajectories 
 	lineplot_stride = int(N_part/10)
 	
 	figps = plt.figure(4)
	for i in range(0, N_part, lineplot_stride):
		plt.plot(Yall[:,0,i], Yall[:,1,i])
	figps.savefig('figs/p_vs_x.pdf', format='pdf')	
	plt.close(figps)

	figx = plt.figure(5)
	for i in range(0, N_part, lineplot_stride):
		plt.plot(a_t, Yall[:,0,i])
	figx.savefig('figs/x_vs_t.pdf', format='pdf')
	plt.close(figx)

        # Phase-space vs time 
	if not os.path.exists('./frames'):
		os.makedirs('./frames')

	frame_stride = 5
	
 	for i in range(0, N_t, frame_stride):	
 		fi = plt.figure(3)
 		print('Saving frame %0.3d of %0.3d' % (i/frame_stride+1, N_t/frame_stride))
		plt.scatter(Yall[i,0,:], Yall[i,1,:], color='#1E90FF', marker='.', alpha=0.1)
		plt.scatter(sum(Yall[i,0,:])/N_part, sum(Yall[i,1,:])/N_part, s=20, color='#000000', marker='+')
 		axes = plt.gca()
 		axes.set_xlim([-4*sigma_x,4*sigma_x])
		axes.set_ylim([-2*sigma_x,2*sigma_x])
 		framenostr = '%0.6d' % (i/frame_stride)
 		fi.savefig('frames/f_' + framenostr + '.png', format='png')
 		plt.close(fi)
 		
	
	#call(["ffmpeg", "-framerate 20", "-i frames/f_%06d.png", "-c:v libx264", "-pix_fmt yuv420p", "beam-ps.mp4"])
	

if __name__ == '__main__':
    main()
