#!/usr/bin/env python

import scipy.integrate as integrate
import numpy as np
import os
import sys, argparse
import time as clock
import ROOT
from ROOT import PData, PGlobals, PPalette, TF1, TH2F, TH1F, TGraph, TCanvas, TColor, TStyle, gStyle, gPad

# Relativistic factor
def Gamma(pz,px):
    return np.sqrt(1 + (pz*pz + px*px) )


def nplasma(t,tout=2000,sout=100) :
    if t<tout :
        return 1.0
    elif t<tout+4*sout :
        return np.exp(-np.power(t - tout, 2.) / (2 * np.power(sout, 2.)))
    else :
        return 0.0 

# Electromagnectic fields definition
def Ez(zeta,s,t):
    return s * nplasma(t) * zeta

def Ex(x,s,t):
    return (1-s) * nplasma(t) * (x/2.0)

def Wx(x,t):
    return nplasma(t) * x/2.0

# Equations of motion:
def dY(Y, t):
    
    # Y is the vector storing the particle coordinates : Y = (z, pz, x, px)
    # dY returns the time diferential of that vector.

    return [ Y[1]/Gamma(Y[1],Y[3]),
             
             -( Ez(Y[0]-t,0.5,t) + ( Ex(Y[2],0.5,t) - Wx(Y[2],t) ) * Y[3] / Gamma(Y[1],Y[3]) ),
             
             Y[3]/Gamma(Y[1],Y[3]),
             
             -( Ex(Y[2],0.5,t) * (1 - Y[1]/Gamma(Y[1],Y[3]) ) + Wx(Y[2],t) * Y[1]/Gamma(Y[1],Y[3])) ]
            

def main():

    # Command argument line parser 
    parser = argparse.ArgumentParser(description='Betatron.py: A particle-in-wakefield trajectory solver.')

    parser.add_argument('--psnaps', action='store_true', default=0,help='plot phasespace snapshots')

    try:
        args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)

    
    # Bunch definition
    # -------------------------------
    # Number of particles
    NP = 500

    # Initial istribution function f(z,pz,x,px)
    sz = 1.0
    spz = 0.1
    sx = 1.0
    spx = 0.5
    meanY = [-1.0, 1000, 0.5, 0.0]   # mean(Y0) 
    covY  = [[sz*sz, 0.0    ,   0.0,     0.0],
             [  0.0, spz*spz,   0.0,     0.0],
             [  0.0, 0.0    , sx*sx,     0.0],
             [  0.0, 0.0    ,   0.0, spx*spx]]

    Z0, PZ0, X0, PX0 = np.random.multivariate_normal(meanY, covY, NP).T
    
    # Time grid
    # --------------------------------
    tmin = 0.0
    dt   = 0.5 * np.sqrt(2 * meanY[1])
    tmax = 200 * dt
    time = np.arange(tmin,tmax,dt)
    NT = len(time)
    # print(' tmax = %f ' % tmax ) 

    # Phase space (Particles and coordinates) 
    # ---------------------------------------
    
    # Matrix storing particle trajectories :
    # (z,pz,x,px) for each particle at each time. 
    Yall = np.empty((NT,4,NP))
    
    # Focusing force definition
    #
    # Wakefield parameters
    # rb = 1.0
    # Delta_rho = 1.0

    # Physics
    # --------------------------------

    t0 = clock.clock()
    print('Solving equations of motions (ODEs) for %i particles...' % NP)

    print('Time grid : t = (%.2f,%.2f) NT = %i  dt = %.2f' % (time[0],time[-1],NT,dt)) 

    for i in range(0, NP):
        # Initial values
        Y0 = [Z0[i], PZ0[i], X0[i], PX0[i]]
        
        # Integrate differential equation of motion
 	Y = integrate.odeint(dY, Y0, time)

        # Array storting the particle trajectories in phase space (z,pz,x,px)
        Yall[:,:,i] = Y

    t1 = clock.clock()
    print('Done in %.3f s. ' % (t1-t0))

    # ---------------------------------
	                
    print('Now analysing and plotting...')

    # Statistics vs time
    # ---------------------------------

    # Statistical moments (Twiss parameters)
    covz_all = np.empty((2,2,NT))
    meanz_all = np.empty((2,NT))
    covx_all = np.empty((2,2,NT))
    meanx_all = np.empty((2,NT))
    
    for i in range(0, NT):	
	meanz_all[:,i] = [sum(Yall[i,0,:])/NP, sum(Yall[i,1,:])/NP]
	covz_all[:,:,i] = np.cov(Yall[i,0,:],Yall[i,1,:])
	meanx_all[:,i] = [sum(Yall[i,2,:])/NP, sum(Yall[i,3,:])/NP]
	covx_all[:,:,i] = np.cov(Yall[i,2,:],Yall[i,3,:])
    
        
    
    # Plotting with ROOT
    # ----------------------------------------------------------------- 

    # Ranges
    xmin = -meanY[2]-6*sx
    xmax = meanY[2]+6*sx
    pxmin = xmin * np.sqrt(1.0/(2.0 * meanY[1])) * meanY[1]
    pxmax = xmax * np.sqrt(1.0/(2.0 * meanY[1])) * meanY[1]

    zetamin = meanY[0]-6*sz-3.0
    zetamax = meanY[0]+6*sz
    pzmin = 0.0
    pzmax = 2000.

    # print(' xmin  = %f  xmax  = %f' % (xmin,xmax) ) 
    # print(' pxmin = %f  pxmax = %f' % (pxmin,pxmax) ) 
    
    # Style
    PGlobals.SetPlasmaStyle()
    gStyle.SetPadRightMargin(0.05);
    gStyle.SetNumberContours(255);
    gStyle.SetPadGridX(1);
    gStyle.SetPadGridY(1);
    gStyle.SetJoinLinePS(2);

    # Canvas
    C = TCanvas('C',' -- O -- ',800,600)
    C.SetFillStyle(4000);
    
    # Plot Fields
    Fxfunc = TF1('Wx',lambda x : Wx(x[0]),xmin,xmax)
    Fxfunc.Draw()

    if not os.path.exists('./beta2figs'):
        os.makedirs('./beta2figs')

    # C.Print('beta2figs/Wx.pdf')
    
    fmeanx = TGraph(NT,time,meanx_all[0,:].flatten())
    fmeanx.SetLineWidth(2)
    fmeanx.SetLineColor(ROOT.kMagenta+2)
    fmeanx.SetMarkerSize(0.4)
    fmeanx.SetMarkerStyle(20)
    fmeanx.SetMarkerColor(ROOT.kMagenta+2)
    fmeanx.Draw("apl")
    C.Print('beta2figs/meanx.pdf')
    
    fmeanpx = TGraph(NT,time,meanx_all[1,:].flatten())
    fmeanpx.SetLineWidth(2)
    fmeanpx.SetLineColor(ROOT.kMagenta+2)
    fmeanpx.SetMarkerSize(0.4)
    fmeanpx.SetMarkerStyle(20)
    fmeanpx.SetMarkerColor(ROOT.kMagenta+2)
    fmeanpx.Draw("apl")
    C.Print('beta2figs/meanpx.pdf')

    frmsx = TGraph(NT,time,np.sqrt(covx_all[0,0,:]).flatten())
    frmsx.SetLineWidth(2)
    frmsx.SetLineColor(ROOT.kMagenta+2)
    frmsx.SetMarkerSize(0.4)
    frmsx.SetMarkerStyle(20)
    frmsx.SetMarkerColor(ROOT.kMagenta+2)
    frmsx.Draw("apl")
    C.Print('beta2figs/rmsx.pdf')
    
    frmspx = TGraph(NT,time,np.sqrt(covx_all[1,1,:]).flatten())
    frmspx.SetLineWidth(2)
    frmspx.SetLineColor(ROOT.kMagenta+2)
    frmspx.SetMarkerSize(0.4)
    frmspx.SetMarkerStyle(20)
    frmspx.SetMarkerColor(ROOT.kMagenta+2)
    frmspx.Draw("apl")
    C.Print('beta2figs/rmspx.pdf')

    fmeanz = TGraph(NT,time,meanz_all[0,:].flatten()-time)
    fmeanz.SetLineWidth(2)
    fmeanz.SetLineColor(ROOT.kMagenta+2)
    fmeanz.SetMarkerSize(0.4)
    fmeanz.SetMarkerStyle(20)
    fmeanz.SetMarkerColor(ROOT.kMagenta+2)
    fmeanz.Draw("apl")
    C.Print('beta2figs/meanz.pdf')
    
    fmeanpz = TGraph(NT,time,meanz_all[1,:].flatten())
    fmeanpz.SetLineWidth(2)
    fmeanpz.SetLineColor(ROOT.kMagenta+2)
    fmeanpz.SetMarkerSize(0.4)
    fmeanpz.SetMarkerStyle(20)
    fmeanpz.SetMarkerColor(ROOT.kMagenta+2)
    fmeanpz.Draw("apl")
    C.Print('beta2figs/meanpz.pdf')

    frmsz = TGraph(NT,time,np.sqrt(covz_all[0,0,:]).flatten())
    frmsz.SetLineWidth(2)
    frmsz.SetLineColor(ROOT.kMagenta+2)
    frmsz.SetMarkerSize(0.4)
    frmsz.SetMarkerStyle(20)
    frmsz.SetMarkerColor(ROOT.kMagenta+2)
    frmsz.Draw("apl")
    C.Print('beta2figs/rmsz.pdf')
    
    frmspz = TGraph(NT,time,np.sqrt(covz_all[1,1,:]).flatten())
    frmspz.SetLineWidth(2)
    frmspz.SetLineColor(ROOT.kMagenta+2)
    frmspz.SetMarkerSize(0.4)
    frmspz.SetMarkerStyle(20)
    frmspz.SetMarkerColor(ROOT.kMagenta+2)
    frmspz.Draw("apl")
    C.Print('beta2figs/rmspz.pdf')

    emit = np.sqrt(covx_all[0,0,:]*covx_all[1,1,:] - covx_all[0,1,:]*covx_all[0,1,:]).flatten()
    femit = TGraph(NT,time,emit)
    femit.SetLineWidth(2)
    femit.SetLineColor(ROOT.kMagenta+2)
    femit.SetMarkerSize(0.4)
    femit.SetMarkerStyle(20)
    femit.SetMarkerColor(ROOT.kMagenta+2)
    femit.GetYaxis().SetRangeUser(emit[0]*(0.8),emit[0]*(300.0))
    femit.Draw("apl")
    C.Print('beta2figs/emit.pdf')

    t2 = clock.clock()
    print('time elapsed %.3f s. ' % (t2-t1))
    
    # Draw trajectories of subset of particles

    Nsample = 100
    if Nsample > NP :
        Nsample = NP
        
    DeltaP = int(NP/Nsample)

    # Color palette
    palette = PPalette('bird')
    palette.SetPalette(ROOT.kBird)
    palette.cd()
    
    frame = TH1F('frame','',10,time[0],time[-1])
    frame.GetYaxis().SetRangeUser(xmin,xmax)
    frame.Draw()
    
    gxvst = []
    for i in range(0, NP, DeltaP):
        gxvst.append(TGraph(NT,time,Yall[:,2,i].flatten()))
        gxvst[-1].SetLineWidth(2)
        cindex = int((abs(Yall[0,2,i])/(meanY[2]+2*sx)) * (palette.GetNColors()-1))
        if cindex>palette.GetNColors()-1 :
            cindex = palette.GetNColors()-1
        gxvst[-1].SetLineColor(palette.GetColorIndex(cindex))
        gxvst[-1].SetMarkerSize(0.4)
        gxvst[-1].SetMarkerStyle(20)
        gxvst[-1].SetMarkerColor(palette.GetColorIndex(cindex))
        gxvst[-1].Draw('pl')

    C.Print('beta2figs/xvst.pdf')

    frame2 = TH1F('frame2','',10,xmin,xmax)
    frame2.GetYaxis().SetRangeUser(pxmin,pxmax)
    frame2.Draw()
    
    gpxvsx = []
    for i in range(0, NP, DeltaP):
        gpxvsx.append(TGraph(NT,Yall[:,2,i].flatten(),Yall[:,3,i].flatten()))
        gpxvsx[-1].SetLineWidth(2)
        cindex = int((abs(Yall[0,2,i])/(2*sx+meanY[2])) * (palette.GetNColors()-1))
        if cindex>palette.GetNColors()-1 :
            cindex = palette.GetNColors()-1

        gpxvsx[-1].SetLineColor(palette.GetColorIndex(cindex))
        gpxvsx[-1].SetMarkerSize(0.2)
        gpxvsx[-1].SetMarkerStyle(20)
        gpxvsx[-1].SetMarkerColor(palette.GetColorIndex(cindex))
        gpxvsx[-1].Draw('pl')
                
    C.Print('beta2figs/pxvsx.pdf')

    frame3 = TH1F('frame3','',10,zetamin,zetamax)
    frame3.GetYaxis().SetRangeUser(xmin,xmax)
    frame3.Draw()

    gxvszeta = []
    for i in range(0, NP, DeltaP):
        gxvszeta.append(TGraph(NT,Yall[:,0,i].flatten()-time[:],Yall[:,2,i].flatten()))
        gxvszeta[-1].SetLineWidth(2)
        cindex = int((abs(Yall[0,2,i])/(2*sx+meanY[2])) * (palette.GetNColors()-1))
        if cindex>palette.GetNColors()-1 :
            cindex = palette.GetNColors()-1

        gxvszeta[-1].SetLineColor(palette.GetColorIndex(cindex))
        gxvszeta[-1].SetMarkerSize(0.2)
        gxvszeta[-1].SetMarkerStyle(20)
        gxvszeta[-1].SetMarkerColor(palette.GetColorIndex(cindex))
        gxvszeta[-1].Draw('pl')
                
    C.Print('beta2figs/xvszeta.pdf')

    frame4 = TH1F('frame4','',10,time[0],time[-1])
    frame4.GetYaxis().SetRangeUser(0.99,1.0)
    frame4.Draw()

    gvzvsz = []
    for i in range(0, NP, DeltaP):
        vgamma = Gamma(Yall[:,1,i].flatten(),Yall[:,3,i].flatten())
        gvzvsz.append(TGraph(NT,Yall[:,0,i].flatten(),Yall[:,1,i].flatten()/vgamma))
        gvzvsz[-1].SetLineWidth(2)
        cindex = int((abs(Yall[0,2,i])/(2*sx+meanY[2])) * (palette.GetNColors()-1))
        if cindex>palette.GetNColors()-1 :
            cindex = palette.GetNColors()-1

        gvzvsz[-1].SetLineColor(palette.GetColorIndex(cindex))
        gvzvsz[-1].SetMarkerSize(0.2)
        gvzvsz[-1].SetMarkerStyle(20)
        gvzvsz[-1].SetMarkerColor(palette.GetColorIndex(cindex))
        gvzvsz[-1].Draw('pl')
                
    C.Print('beta2figs/vzvsz.pdf')

    t3 = clock.clock()
    print('time elapsed %.3f s. ' % (t3-t2))


    if args.psnaps :
        if not os.path.exists('./beta2figs/phasespace'):
	    os.makedirs('./beta2figs/phasespace')

        # gPad.SetRightMargin(0.20)
        palette.SetPalette('electron0')

        print('Plotting phasespace distributions...')

        Nsample = 50
        if Nsample > NT :
            Nsample = NT            
        DeltaT = int(NT/Nsample)

        for i in range(0, NT, DeltaT):	
	    # meanx_all[:,i] = [sum(Yall[i,0,:])/NP, sum(Yall[i,1,:])/NP]
	    # covx_all[:,:,i] = np.cov(Yall[i,0,:],Yall[i,1,:])
            # sx_t = np.sqrt(covx_all[0,0,i])
            # spx_t = np.sqrt(covx_all[1,1,i])
                
            hname = 'hpzvsz-%i' % i 
            hpzvsz = TH2F(hname,'',100,zetamin,zetamax,100,pzmin,pzmax)
            hname = 'hpxvsx-%i' % i 
            hpxvsx = TH2F(hname,'',100,xmin,xmax,100,pxmin,pxmax)
            hname = 'hxvszeta-%i' % i 
            hxvszeta = TH2F(hname,'',100,zetamin,zetamax,100,xmin,xmax)

            for j in range(0,NP):
                hpzvsz.Fill(Yall[i,0,j]-time[i],Yall[i,1,j])                  
                hpxvsx.Fill(Yall[i,2,j],Yall[i,3,j])                  
                hxvszeta.Fill(Yall[i,0,j]-time[i],Yall[i,2,j])                  

            # subs = ''
            # if i == 0 :
            #    subs = '('
            #elif i == NT-1 :
            #    subs = ')'
                        
            #hpxvsx.Draw('col')
            #C.Print('beta2figs/phasespace/hpxvsx.pdf'+subs)
            #hpzvsz.Draw('col')        
            #C.Print('beta2figs/phasespace/hpzvsz.pdf'+subs)
            #hxvszeta.Draw('col')        
            #C.Print('beta2figs/phasespace/hxvszeta.pdf'+subs)
                    
            hpxvsx.Draw('col')        
            C.Print('beta2figs/phasespace/hpxvsx-%i.png' % i)
            hpzvsz.Draw('col')        
            C.Print('beta2figs/phasespace/hpzvsz-%i.png' % i)
            hxvszeta.Draw('col')        
            C.Print('beta2figs/phasespace/hxvszeta-%i.png' % i)


    t4 = clock.clock()
    print('time elapsed %.3f s. ' % (t4-t3))

    print('Output at : %s' % 'beta2figs/' )
    print('Ciao.')


if __name__ == '__main__':
    main()
