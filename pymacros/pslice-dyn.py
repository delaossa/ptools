#!/usr/bin/env python

import scipy.integrate as integrate
import numpy as np
import os
import types
import collections
import ROOT
from ROOT import PData, PGlobals, PPalette, TF1, TH2F, TH1F, TGraph, TCanvas, TColor, TStyle, gStyle, gPad

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
		f = 0
                
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
	N_t = 20
	N_part = 1000
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
	rb = 2.0
	Delta_rho = 1.0
 	global wake 
 	wake = Wake(rb, Delta_rho)
              
        print('Computing ODEs...')
 	for i in range(0, N_part):
                # Initial values
                x0 = X0[i]
 		p0 = P0[i]
                # Integrate differential equation of motion
 		a_sol = integrate.odeint(dY, [x0, p0], a_t)
 		Yall[:,:,i] = a_sol
	print('done!')
	                
        # Statistics vs time
	cov_all = np.empty((2,2,N_t))
	mean_all = np.empty((2,N_t))

        for i in range(0, N_t):	
	        mean_all[:,i] = [sum(Yall[i,0,:])/N_part, sum(Yall[i,1,:])/N_part]
	        cov_all[:,:,i] = np.cov(Yall[i,0,:],Yall[i,1,:])

        # Plotting with ROOT
        
        # Style
        PGlobals.SetPlasmaStyle()
        gStyle.SetPadRightMargin(0.05);
        gStyle.SetNumberContours(255);

        # Canvas
        C = TCanvas('C','',800,600)
        C.SetFillStyle(4000);

        # Plot Fx(x)
        # def myWake(x):
        #        f = wake.Fx(x[0])
        #        return f
                
        # Fxfunc = TF1('Fxfunc',myWake,-4*sigma_x,4*sigma_x)
        
        Fxfunc = TF1('Fxfunc',lambda x : wake.Fx(x[0]),-4*sigma_x,4*sigma_x)
        Fxfunc.Draw()

        if not os.path.exists('./figs'):
                os.makedirs('./figs')

        C.Print('figs/Fx.pdf')
        
        fmeanx = TGraph(N_t,a_t,mean_all[0,:].T)
        fmeanx.SetLineWidth(2)
        fmeanx.SetLineColor(ROOT.kMagenta+2)
        fmeanx.SetMarkerSize(1)
        fmeanx.SetMarkerStyle(20)
        fmeanx.SetMarkerColor(ROOT.kMagenta+2)
        fmeanx.Draw("apc")
        C.Print('figs/meanx.pdf')

        fmeanpx = TGraph(N_t,a_t,mean_all[1,:].T)
        fmeanpx.SetLineWidth(2)
        fmeanpx.SetLineColor(ROOT.kMagenta+2)
        fmeanpx.SetMarkerSize(1)
        fmeanpx.SetMarkerStyle(20)
        fmeanpx.SetMarkerColor(ROOT.kMagenta+2)
        fmeanpx.Draw("apc")
        C.Print('figs/meanpx.pdf')

        # Randomly pick a subset of particles to track their trajectories 
 	lineplot_stride = int(N_part/10)

        # Palette
        palette = PPalette('bird')
        palette.SetPalette(ROOT.kBird)
        palette.cd()
        
        frame = TH1F('frame','',10,a_t[0],a_t[-1])
        frame.GetYaxis().SetRangeUser(-2*sigma_x,2*sigma_x)
        frame.Draw('axis')

        gxvst = []
	for i in range(0, N_part, lineplot_stride):
                gxvst.append(TGraph(N_t,a_t,Yall[:,0,i].flatten()))
                gxvst[-1].SetLineWidth(2)
                cindex = int((abs(Yall[0,0,i])/(3*sigma_x)) * (palette.GetNColors()-1))
                if cindex>palette.GetNColors()-1 :
                        cindex = palette.GetNColors()-1
                gxvst[-1].SetLineColor(palette.GetColorIndex(cindex))
                gxvst[-1].SetMarkerSize(0.2)
                gxvst[-1].SetMarkerStyle(20)
                gxvst[-1].SetMarkerColor(palette.GetColorIndex(cindex))
                gxvst[-1].Draw('pc')

        C.Print('figs/xvst.pdf')

        frame2 = TH1F('frame2','',10,-2*sigma_x,2*sigma_x)
        frame2.GetYaxis().SetRangeUser(-2*sigma_x/np.sqrt(2),2*sigma_x/np.sqrt(2))
        frame2.Draw('axis')
                
        gpvsx = []
	for i in range(0, N_part, lineplot_stride):
                gpvsx.append(TGraph(N_t,Yall[:,0,i].flatten(),Yall[:,1,i].flatten()))
                gpvsx[-1].SetLineWidth(2)
                cindex = int((abs(Yall[0,0,i])/(3*sigma_x)) * (palette.GetNColors()-1))
                gpvsx[-1].SetLineColor(palette.GetColorIndex(cindex))
                gpvsx[-1].SetMarkerSize(0.2)
                gpvsx[-1].SetMarkerStyle(20)
                gpvsx[-1].SetMarkerColor(palette.GetColorIndex(cindex))
                gpvsx[-1].Draw('pc')
                
        C.Print('figs/pvsx.pdf')
        
        if not os.path.exists('./figs/phasespace'):
		os.makedirs('./figs/phasespace')

        gPad.SetRightMargin(0.20)
        
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

#                if i == 0 :
#                        C.Print('figs/hpxvsx.pdf(')
#                elif i == N_t-1 :
#                        C.Print('figs/hpxvsx.pdf)')
#                else :
#                        C.Print('figs/hpxvsx.pdf')

                C.Print('figs/phasespace/hpxvsx-%i.png' % i)
                                      

        

if __name__ == '__main__':
    main()
