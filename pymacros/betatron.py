#!/usr/bin/env python

import scipy.integrate as integrate
import scipy.constants as ct
import numpy as np
import os, sys, argparse
import time as clock
from ROOT import PData, PGlobals, PPalette, PUnits, TF1, TH2F, TH1F, TGraph, TCanvas, TColor, TStyle, gStyle, gPad, kMagenta, kBird


def parse_args():

    # Command argument line parser 
    parser = argparse.ArgumentParser(description='Betatron.py: A particle-in-wakefield trajectory solver.')

    parser.add_argument('name', nargs='?', default='betaout', help='Simulation name (output folder)')
    parser.add_argument('-NP', type=int, dest='NP', default=1000, help='Number of macroparticles')
    parser.add_argument('--ps', action='store_true', default=1, help='Plot phasespace snapshots')
    parser.add_argument('--png', action='store_true', default=0, help='Plot phasespace snapshots in png')
    parser.add_argument('--si', action='store_true', default=0, help='SI units')
    parser.add_argument('--man', action='store_true', default=0, help='Print help')
    
    args = parser.parse_args()

    if args.man :
        parser.print_help()
        sys.exit(0)
        
    return args


# Relativistic factor
def Gamma(pz,px):
    return np.sqrt( 1 + (pz*pz + px*px) )


# Plasma density profile
class Plasma:
    def __init__(self,zout,sout) :
        self.zout = zout
        self.sout = sout
        
    def n(self,z) :
        if z<self.zout :
            return 1.0
        elif z<self.zout+4*self.sout :
            return np.exp(-np.power(z - self.zout, 2.) / (2 * np.power(self.sout, 2.)))
        else :
            return 0.0 

        
# Wakefields: Electromagnectic fields definition
class Wake:
    def __init__(self,slope=0.5,rb=10.0,Drho=1.0):
        self.slope = slope
        self.rb = rb
        self.Drho = Drho

    def Ez(self, zeta, x) :
        # return 0
        if abs(x) < self.rb :
            return self.slope * zeta 
        else :
            return self.slope * zeta * np.exp(-(abs(x)-self.rb)/self.Drho)
      
    def Ex(self, x) :
        if abs(x) < self.rb :
            return (1-self.slope) * ( x/2.0 )
        else :
            return (1-self.slope) * ( (np.sign(x) * self.rb)/2.0 ) * np.exp(-(np.abs(x)-self.rb)/self.Drho)
        
    def Wx(self, x) :
         if abs(x) < self.rb :
             return x/2.0
         else :
             return (np.sign(x) * self.rb)/2.0 * np.exp(-(np.abs(x)-self.rb)/self.Drho)
        

        
# Equations of motion:
def dY(Y, t, wake, plasma):
    
    # Y is the vector of particle coordinates : Y = (z, pz, x, px)
    # dY returns the time diferential of that vector.

    kp = np.sqrt(plasma.n(Y[0]))
    
    return [ Y[1]/Gamma(Y[1],Y[3]),
             
             - kp * ( wake.Ez(kp*(Y[0]-t),Y[2]) + ( wake.Ex(kp*Y[2]) - wake.Wx(kp*Y[2]) ) * Y[3] / Gamma(Y[1],Y[3]) ),
             
             Y[3]/Gamma(Y[1],Y[3]),
             
             - kp * ( wake.Ex(kp*Y[2]) * (1 - Y[1]/Gamma(Y[1],Y[3]) ) + wake.Wx(kp*Y[2]) * Y[1]/Gamma(Y[1],Y[3])) ]



def main():

    args = parse_args()

        
    # Wakefields definition
    # -------------------------------
    rb = 1.0
    Drho = 1.0
    slope = 0.1
    wake = Wake(slope,rb,Drho)

    # -------------------------------

    
    # Bunch definition
    # -------------------------------

    # Number of particles
    NP = args.NP
    print('Bunch definition : N = %i particles' % NP)

    # Initial istribution function f(z,pz,x,px)

    # Initial statistical moments (means and rms')
    #       [    z,   pz,   x,  px ]
    meanY = [ -1.0,  100, 0.0, 0.0 ]   # mean(Y0) 
    rmsY  = [  0.1,  0.1, 1.0, 1.0 ]   # rms(Y0)         

    emit0 = 0.1
    rmsY[3]  = emit0 / rmsY[2]

    if args.si :
        n0 = 1.0E23 # m^-3
        kp = np.sqrt( n0 * ct.e * ct.e / (ct.epsilon_0 * ct.m_e) ) / ct.c
        mc2 = ct.m_e * ct.c * ct.c

        # print('Speed of light %f' % ct.c)
        # print('Skindepth = %f' % ((1/kp)/ct.micro))
        # print('Electron mass = %f MeV' % (mc2 / (ct.eV * ct.mega)))

        meanY_SI = np.array([-10 * ct.micro, 100 * ct.mega * ct.eV, 0.0, 0.0 ])
        rmsY_SI  = np.array([2.0 * ct.micro, 1.0 * ct.mega * ct.eV, 20.0 * ct.micro, 1.0 * ct.kilo * ct.eV])
        
        normarray = np.array([kp,1/mc2,kp,1/mc2])
        
        meanY = meanY_SI * normarray
        rmsY  = rmsY_SI * normarray      
        print(meanY)
        print(rmsY)
        
    # Blowout betatron frequency
    K_BO = 0.5
    omegaBO = np.sqrt(K_BO/meanY[1])
    print('Betatron wavelength = %.2f' % (2*np.pi/omegaBO) )

    # Second moments, covariance matrix:
    sz  = rmsY[0]
    spz = rmsY[1] 
    sx  = rmsY[2]
    spx = rmsY[3] 
        
    covY  = [[sz*sz, 0.0    ,   0.0,     0.0],
             [  0.0, spz*spz,   0.0,     0.0],
             [  0.0, 0.0    , sx*sx,     0.0],
             [  0.0, 0.0    ,   0.0, spx*spx]]

    Z0, PZ0, X0, PX0 = np.random.multivariate_normal(meanY, covY, NP).T

    # Matched emittance at <xpx> = 0
    emitM = meanY[1] * sx * sx * omegaBO

    # Time grid
    # --------------------------------
    tmin = 0.0
    dt   = 1.0 / omegaBO
    tmax = 150 * dt
    time = np.arange(tmin,tmax,dt)
    NT = len(time)
    # print(' tmax = %f ' % tmax ) 
    print('Time grid : t = (%.2f,%.2f) NT = %i  dt = %.2f' % (time[0],time[-1],NT,dt)) 

    
    # Plasma profile definition
    plasma = Plasma(100 * dt,10 * dt)
        
    # ODE solver
    # ----------

    tclock = clock.clock()
    print('Solving equations of motions (ODEs) ...')

    # Phase space (Particles and coordinates) 
    # ---------------------------------------
    
    # Matrix storing particle trajectories :
    # (z,pz,x,px) for each particle at each time. 
    Yall = np.empty((NT,4,NP))

    for i in range(0, NP):
        # Initial values
        Y0 = [Z0[i], PZ0[i], X0[i], PX0[i]]
        
        # Integrate differential equation of motion
 	Y = integrate.odeint(dY, Y0, time, args=(wake,plasma))

        # Array storting the particle trajectories in phase space (z,pz,x,px)
        Yall[:,:,i] = Y

    print('Done in %.3f s. ' % (clock.clock() - tclock))
    tclock = clock.clock()

    # ---------------------------------
	                
    print('Now analysing...')
    
    # Statistics vs time
    # ---------------------------------

    print('Extracting statistical moments...')

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
    

    print('%.3f seconds ' % (clock.clock() - tclock))
    tclock = clock.clock()
    
    
    # Plotting with ROOT
    # ----------------------------------------------------------------- 

    print('Plotting beam moments...')

    # Ranges
    xmin = -meanY[2]-6*sx
    xmax = meanY[2]+6*sx
    pxmin = xmin * np.sqrt(1.0/(2.0 * meanY[1])) * meanY[1]
    pxmax = xmax * np.sqrt(1.0/(2.0 * meanY[1])) * meanY[1]

    zetamin = meanY[0]-6*sz-2.0
    zetamax = meanY[0]+6*sz
    pzmin = 0.0
    pzmax = 1000.

    # print(' xmin  = %f  xmax  = %f' % (xmin,xmax) ) 
    # print(' pxmin = %f  pxmax = %f' % (pxmin,pxmax) ) 
    
    # Style
    PGlobals.SetPlasmaStyle()
    gStyle.SetPadRightMargin(0.05)
    gStyle.SetNumberContours(255)
    gStyle.SetPadGridX(1)
    gStyle.SetPadGridY(1)
    gStyle.SetJoinLinePS(2)
    gStyle.SetFrameLineWidth(2)

    # Canvas
    C = TCanvas('C',' -- O -- ',800,600)
    C.SetFillStyle(4000);

    opath = 'betaout'
    if args.name :
        opath = args.name
    
    if not os.path.exists(opath):
        os.makedirs(opath)

    # Plot Fields
    fWx = TF1('Wx',lambda x : wake.Wx(x[0]),xmin,xmax)
    fWx.SetNpx(1000)
    fWx.Draw()

    fEx = TF1('Ex',lambda x : wake.Ex(x[0]),xmin,xmax)
    fEx.SetNpx(1000)
    fEx.SetLineStyle(2)
    fEx.Draw("same")

    if not os.path.exists(opath + '/fields'):
        os.makedirs(opath + '/fields')

    C.Print(opath +'/fields/Wx&Ex.pdf')

    fEz = TF1('Ez',lambda x : wake.Ez(x[0],0.0),zetamin,zetamax)
    fEz.Draw()
    C.Print(opath +'/fields/Ez.pdf')

    fProf = TF1('Prof',lambda x : plasma.n(x[0]),time[0],time[-1])
    fProf.SetNpx(1000)
    fProf.Draw()
    C.Print(opath +'/fields/DenProf.pdf')

    lcolor = kMagenta + 2
    
    fmeanx = TGraph(NT,time,meanx_all[0,:].flatten())
    fmeanx.SetLineWidth(2)
    fmeanx.SetLineColor(lcolor)
    fmeanx.SetMarkerSize(0.4)
    fmeanx.SetMarkerStyle(20)
    fmeanx.SetMarkerColor(lcolor)
    fmeanx.Draw("apl")
    C.Print(opath +'/meanx.pdf')
    
    fmeanpx = TGraph(NT,time,meanx_all[1,:].flatten())
    fmeanpx.SetLineWidth(2)
    fmeanpx.SetLineColor(lcolor)
    fmeanpx.SetMarkerSize(0.4)
    fmeanpx.SetMarkerStyle(20)
    fmeanpx.SetMarkerColor(lcolor)
    fmeanpx.Draw("apl")
    C.Print(opath +'/meanpx.pdf')

    frmsx = TGraph(NT,time,np.sqrt(covx_all[0,0,:]).flatten())
    frmsx.SetLineWidth(2)
    frmsx.SetLineColor(lcolor)
    frmsx.SetMarkerSize(0.4)
    frmsx.SetMarkerStyle(20)
    frmsx.SetMarkerColor(lcolor)
    frmsx.Draw("apl")
    C.Print(opath +'/rmsx.pdf')
    
    frmspx = TGraph(NT,time,np.sqrt(covx_all[1,1,:]).flatten())
    frmspx.SetLineWidth(2)
    frmspx.SetLineColor(lcolor)
    frmspx.SetMarkerSize(0.4)
    frmspx.SetMarkerStyle(20)
    frmspx.SetMarkerColor(lcolor)
    frmspx.Draw("apl")
    C.Print(opath +'/rmspx.pdf')

    fmeanz = TGraph(NT,time,meanz_all[0,:].flatten()-time)
    fmeanz.SetLineWidth(2)
    fmeanz.SetLineColor(lcolor)
    fmeanz.SetMarkerSize(0.4)
    fmeanz.SetMarkerStyle(20)
    fmeanz.SetMarkerColor(lcolor)
    fmeanz.Draw("apl")
    C.Print(opath +'/meanz.pdf')
    
    fmeanpz = TGraph(NT,time,meanz_all[1,:].flatten())
    fmeanpz.SetLineWidth(2)
    fmeanpz.SetLineColor(lcolor)
    fmeanpz.SetMarkerSize(0.4)
    fmeanpz.SetMarkerStyle(20)
    fmeanpz.SetMarkerColor(lcolor)
    fmeanpz.Draw("apl")
    C.Print(opath +'/meanpz.pdf')

    frmsz = TGraph(NT,time,np.sqrt(covz_all[0,0,:]).flatten())
    frmsz.SetLineWidth(2)
    frmsz.SetLineColor(lcolor)
    frmsz.SetMarkerSize(0.4)
    frmsz.SetMarkerStyle(20)
    frmsz.SetMarkerColor(lcolor)
    frmsz.Draw("apl")
    C.Print(opath +'/rmsz.pdf')
    
    frmspz = TGraph(NT,time,np.sqrt(covz_all[1,1,:]).flatten())
    frmspz.SetLineWidth(2)
    frmspz.SetLineColor(lcolor)
    frmspz.SetMarkerSize(0.4)
    frmspz.SetMarkerStyle(20)
    frmspz.SetMarkerColor(lcolor)
    frmspz.Draw("apl")
    C.Print(opath +'/rmspz.pdf')

    emit = np.sqrt(covx_all[0,0,:]*covx_all[1,1,:] - covx_all[0,1,:]*covx_all[0,1,:]).flatten()
    femit = TGraph(NT,time,emit)
    femit.SetLineWidth(2)
    femit.SetLineColor(lcolor)
    femit.SetMarkerSize(0.4)
    femit.SetMarkerStyle(20)
    femit.SetMarkerColor(lcolor)
    femit.GetYaxis().SetRangeUser(emit[0]*(0.8),emit[0]*(300.0))
    femit.Draw("apl")
    C.Print(opath +'/emit.pdf')

    print('%.3f seconds ' % (clock.clock() - tclock))
    tclock = clock.clock()


    # Draw trajectories of subset of particles
    Nsample = 100
    if Nsample > NP :
        Nsample = NP
        
    DeltaP = int(NP/Nsample)

    if not os.path.exists(opath + '/tracks'):
	os.makedirs(opath + '/tracks')

    print('Plotting particle trajectories (Sample = %i part.)' % Nsample)
    
    # Color palette
    palette = PPalette('bird')
    palette.SetPalette(kBird)
    palette.cd()
    palette.SetAlpha(0.6)
    
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

    gPad.RedrawAxis()
    gPad.RedrawAxis('g')
    C.Print(opath +'/tracks/xvst.pdf')

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
                
    gPad.RedrawAxis()
    gPad.RedrawAxis('g')
    C.Print(opath +'/tracks/pxvsx.pdf')

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
                
    gPad.RedrawAxis()
    gPad.RedrawAxis('g')
    C.Print(opath +'/tracks/xvszeta.pdf')

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
                
    gPad.RedrawAxis()
    gPad.RedrawAxis('g')
    C.Print(opath +'/tracks/vzvsz.pdf')
        

    frame5 = TH1F('frame5','',10,time[0],time[-1])
    frame5.GetYaxis().SetRangeUser(1.0,pzmax)
    frame5.Draw()

    ggammavst = []
    for i in range(0, NP, DeltaP):
        vgamma = Gamma(Yall[:,1,i].flatten(),Yall[:,3,i].flatten())
        ggammavst.append(TGraph(NT,time,vgamma))
        ggammavst[-1].SetLineWidth(2)
        cindex = int((abs(Yall[0,2,i])/(2*sx+meanY[2])) * (palette.GetNColors()-1))
        if cindex>palette.GetNColors()-1 :
            cindex = palette.GetNColors()-1

        ggammavst[-1].SetLineColor(palette.GetColorIndex(cindex))
        ggammavst[-1].SetMarkerSize(0.2)
        ggammavst[-1].SetMarkerStyle(20)
        ggammavst[-1].SetMarkerColor(palette.GetColorIndex(cindex))
        ggammavst[-1].Draw('pl')
                
    gPad.RedrawAxis()
    gPad.RedrawAxis('g')
    C.Print(opath +'/tracks/gammavst.pdf')
        

    print('%.3f seconds ' % (clock.clock() - tclock))
    tclock = clock.clock()


    if args.ps :
        if not os.path.exists(opath + '/phasespace'):
	    os.makedirs(opath + '/phasespace')

        # gPad.SetRightMargin(0.20)
        palette.SetPalette('electron0')
        palette.SetAlpha(1.0)

        print('Plotting phasespace distributions...')

        Nsample = 50
        if Nsample > NT :
            Nsample = NT            
        DeltaT = int(NT/Nsample)

        tindex = range(0, NT, DeltaT)
        for i in tindex:	
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

            if args.png :
                hpxvsx.Draw('col')        
                C.Print(opath +'/phasespace/hpxvsx-%i.png' % i)
                hpzvsz.Draw('col')        
                C.Print(opath +'/phasespace/hpzvsz-%i.png' % i)
                hxvszeta.Draw('col')        
                C.Print(opath +'/phasespace/hxvszeta-%i.png' % i)
            else :
                subs = ''
                if i == 0 : subs = '('
                elif i == tindex[-1] : subs = ')'
                   
                hpxvsx.Draw('col')
                C.Print(opath +'/phasespace/hpxvsx.pdf'+subs)
                hpzvsz.Draw('col')        
                C.Print(opath +'/phasespace/hpzvsz.pdf'+subs)
                hxvszeta.Draw('col')        
                C.Print(opath +'/phasespace/hxvszeta.pdf'+subs)


        print('%.3f seconds ' % (clock.clock() - tclock))
        tclock = clock.clock()

    print('Output at : %s' % opath +'/' )
    print('Ciao.')


if __name__ == '__main__':
    main()
