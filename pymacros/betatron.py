#!/usr/bin/env python

import scipy.integrate as integrate
import scipy.constants as ct
import numpy as np
import os, sys, argparse
import time as clock
import multiprocessing
from functools import partial
import ROOT
from ROOT import PData, PGlobals, PPalette, PUnits, TF1, TH2F, TH1F, TGraph, TCanvas, TBox, TLine, TColor, TStyle, gStyle, gPad, kMagenta, kBird, kBlack
import wakefields as wk

def parse_args():

    # Command argument line parser 
    parser = argparse.ArgumentParser(description='Betatron.py: A particle-in-wakefield trajectory solver.')

    parser.add_argument('name', nargs='?', default='betaout', help='Simulation name (output folder)')
    parser.add_argument('-NP', type=int, dest='NP', default=1000, help='Number of macroparticles')
    parser.add_argument('--ps', action='store_true', default=1, help='Plot phasespace snapshots')
    parser.add_argument('--png', action='store_true', default=0, help='Plot phasespace snapshots in png (slow!)')
    parser.add_argument('--si', action='store_true', default=0, help='SI units')
    parser.add_argument('--mp', action='store_true', default=0, help='For parallel processing')
    parser.add_argument('--man', action='store_true', default=0, help='Print help')
    
    args = parser.parse_args()

    if args.man :
        parser.print_help()
        sys.exit(0)
        
    return args

def DrawFrame() :

    gPad.Update()
    gPad.RedrawAxis('g')
    gPad.RedrawAxis()

    l = TLine()
    l.SetLineColor(kBlack)
    l.SetLineWidth(2)
    l.DrawLine(gPad.GetUxmin(), gPad.GetUymax(), gPad.GetUxmax(), gPad.GetUymax())
    l.DrawLine(gPad.GetUxmax(), gPad.GetUymin(), gPad.GetUxmax(), gPad.GetUymax())
    l.DrawLine(gPad.GetUxmin(), gPad.GetUymin(), gPad.GetUxmin(), gPad.GetUymax())
    l.DrawLine(gPad.GetUxmin(), gPad.GetUymin(), gPad.GetUxmax(), gPad.GetUymin())
        
    
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

   
        
# Equations of motion:
def dY(Y, t, wake, plasma):
    
    # Y is the vector of particle coordinates : Y = (z, pz, x, px)
    # dY returns the time diferential of that vector.

    kp = np.sqrt(plasma.n(Y[0]))
    
    return [ Y[1]/Gamma(Y[1],Y[3]),
             
             - kp * ( wake.Ez(kp*(Y[0]-t),kp*Y[2]) + ( wake.Ex(kp*(Y[0]-t),kp*Y[2]) - wake.Wx(kp*(Y[0]-t),kp*Y[2]) ) * Y[3] / Gamma(Y[1],Y[3]) ),
             
             Y[3]/Gamma(Y[1],Y[3]),
             
             - kp * ( wake.Wx(kp*(Y[0]-t),kp*Y[2]) * Y[1]/Gamma(Y[1],Y[3]) + wake.Ex(kp*(Y[0]-t),kp*Y[2]) * (1 - Y[1]/Gamma(Y[1],Y[3]) ) )  ]


def integra(wake,plasma,time,dY,Y) :
    return integrate.odeint(dY, Y, time, args=(wake,plasma))

def statistics(Y) :
    NP = len(Y[0])

    meanz = [sum(Y[0])/NP, sum(Y[1])/NP]
    meanx = [sum(Y[2])/NP, sum(Y[3])/NP]
    
    covz = np.cov(Y[0],Y[1])
    covx = np.cov(Y[2],Y[3])

    return meanz,meanx,covz,covx


def main():

    args = parse_args()
    
    # Wakefields definition
    # -------------------------------
    Ez0 = -0.1
    K = 0.174
    S = 0.02
    Ks = 0.1
    rb = 10.0
    Drho = 1.0
    wake = wk.WakeBlowout(Ez0,K,S,Ks,rb,Drho)
    # -------------------------------

    a0 = 0.8
    nb0 = a0**2/2
    L  = np.sqrt(2)
    #L  = 1.49
    sbz = L/np.sqrt(2)
    w0 = 2.976
    sbx = w0/2

    print('Setting linear wakefields with : nb0 = %.2f, sigma_z = %.3f, sigma_r = %.3f' % (nb0,sbz,sbx) ) 
    # wake = wk.WakeLinearSimple(nb0,sbz,sbx)
    # -------------------------------------

    # Bunch definition
    # -------------------------------

    # Number of particles
    NP = args.NP
    print('Bunch definition : N = %i particles' % NP)

    # Initial istribution function f(z,pz,x,px)

    # Initial statistical moments (means and rms')
    #       [    z,   pz,   x,  px ]
    meanY = [ -1.0,  100, 0.0, 0.0 ]   # mean(Y0) 
    rmsY  = [  1.0,  0.1, 0.1, 0.1 ]   # rms(Y0)         

    # Blowout betatron frequency
    omegaB = np.sqrt(K/meanY[1])
    lambdaB = 2*np.pi/omegaB
    # print('Betatron wavelength = %.2f' % lambdaB )

    # Matched conditions at <xpx> = 0
    betaM = np.sqrt(meanY[1]/K) # matched beta
    emitM = meanY[1] * np.power(rmsY[2],2) / betaM  # matched emittance
    emit0 = emitM
    # emit0 = 0.1
    rmsY[3]  = emit0 / rmsY[2]

    # Bunch definition in SI units    
    if args.si :

        # Normalization constants
        n0 = 1.0E23 # m^-3
        kp = np.sqrt( n0 * ct.e * ct.e / (ct.epsilon_0 * ct.m_e) ) / ct.c
        mc2 = ct.m_e * ct.c * ct.c
        MeV = ct.mega * ct.eV
        um  = ct.micron

        print('Speed of light %f m/s' % ct.c)
        print('Skindepth = %f um' % ((1/kp)/um))
        print('Electron mass = %f MeV' % (mc2 / MeV))
        
        # Initial statistical moments (mean and rms)   
        #                   [        z,        pz,         x,           px  ]
        meanY_SI = np.array([  0.0 * um,     5 * MeV,  0.0 * um,  0.0 * MeV ])
        rmsY_SI  = np.array([  3.0 * um,0.0325 * MeV,  0.0 * um,  0.0 * MeV ])

        meangamma = meanY_SI[1]/mc2
        print('Gamma = %f' % meangamma)

        # Blowout betatron frequency
        K = wake.Kx(kp*meanY_SI[0],0)
        omegaB = np.sqrt(K/meangamma)
        lambdaB = 2*np.pi/omegaB
        print('Betatron wavelength = %.3f mm' % (lambdaB/(kp*1000*um)) )

        # Matched beta 
        betaM = (1/omegaB)/kp       
        print('Matched beta = %.3f mm' % (betaM/(1000*um)) )

        # Set norm. emittance 
        emit0 = 0.30 * um
        print('Emittance = %.3f um' % (emit0/um) )

        # Matched beam width
        sxM = np.sqrt(betaM * emit0 / meangamma)
        beta0 = 5.38 * betaM 
        #beta0 = betaM 
        sx0 = np.sqrt(beta0 * emit0 / meangamma)
        
        # Beam width
        rmsY_SI[2] = sx0
        print('Sigmax    = %.3f um' % (rmsY_SI[2]/um))

        # emitM = meangamma * np.power(rmsY_SI[2],2) / betaM  # matched emittance
        # print('Matched emittance = %f um' % (emitM/um) )        
        # emit0 = emitM                
        beta0 =  meangamma * np.power(rmsY_SI[2],2) / emit0 
        print('Beta      = %.3f mm' % (beta0/(1000*um)) )
        rmsY_SI[3]  = mc2 * emit0 / rmsY_SI[2]
        
        # Normalizes
        norma = np.array([kp,1/mc2,kp,1/mc2])
        # norma = np.array([1,1,1,1])
        
        meanY = meanY_SI * norma
        rmsY  = rmsY_SI * norma      
        

    # Second moments, covariance matrix initialization.
    sz  = rmsY[0]
    spz = rmsY[1] 
    sx  = rmsY[2]
    spx = rmsY[3] 
        
    covY  = [[sz*sz, 0.0    ,   0.0,     0.0],
             [  0.0, spz*spz,   0.0,     0.0],
             [  0.0, 0.0    , sx*sx,     0.0],
             [  0.0, 0.0    ,   0.0, spx*spx]]

    Y0 = np.random.multivariate_normal(meanY, covY, NP)

    # Plasma profile definition
    # ----------------------------------------- 
    plasma = Plasma(5 * lambdaB , 0.05 * lambdaB)
    # plasma = Plasma(2400 , 2 * lambdaB)

    # Time grid
    # --------------------------------
    tmin = 0.0
    dt   = 0.2 / omegaB
    tmax = 5.7 * lambdaB
    # tmax = 2400
    time = np.arange(tmin,tmax,dt)
    NT = len(time)
    # print(' tmax = %f ' % tmax ) 
    print('Time grid : t = (%.2f,%.2f) NT = %i  dt = %.2f' % (time[0],time[-1],NT,dt))   
        
    # ODE solver
    # ----------

    tclock = clock.clock()
    print('Solving equations of motion (ODEs) ...')

    # Phase space (Particles and coordinates) 
    # ---------------------------------------
    
    # Matrix storing particle trajectories :
    # (z,pz,x,px) for each particle at each time. 
    Yall = np.empty((NP,NT,4))

    # Serial version:
    if args.mp == 0 :
        for i in range(NP):        
            # Integrate differential equation of motion
 	    # Y = integrate.odeint(dY, Y0[i], time, args=(wake,plasma))
            Yall[i] = integra(wake,plasma,time,dY,Y0[i])

            if args.si :
                units = np.array([(1/kp)/um,1,(1/kp)/um,1])
                Yall[i,:] = Yall[i,:] * units    

    else :
    # Parallel version:
        pintegra = partial(integra,wake,plasma,time,dY)
        pool = multiprocessing.Pool()
        NCPU = multiprocessing.cpu_count()
        print('Number of CPUs = %i' % NCPU)
    
        Ylist = pool.map(pintegra,(Y0[i] for i in range(NP)) )
        pool.close()
        pool.join()

        # This is needed to avoid broadcasting problems 
        for i in range(NP):
            Yall[i] = Ylist[i]

            if args.si :
                units = np.array([(1/kp)/um,1,(1/kp)/um,1])
                Yall[i,:] = Yall[i,:] * units    
                
    # ----------------------------------------------------
        
    print('Done in %.3f s. ' % ((clock.clock() - tclock)))
    tclock = clock.clock()

    # ---------------------------------
	                
    print('Now analysing...')
    
    # Statistics vs time
    # ---------------------------------

    print('Extracting statistical moments...')

    # Statistical moments 
    meanz_all = np.empty((NT,2))
    covz_all  = np.empty((NT,2,2))
    meanx_all = np.empty((NT,2))
    covx_all  = np.empty((NT,2,2))
    
    # Twiss parameters
    emitx_all = np.empty(NT)
    betax_all = np.empty(NT)
    alphax_all = np.empty(NT)
                 
    # Serial version:
    if args.mp == 0 :
        for i in range(NT):

            # This bit can be quite lengthy for big NP --> parallelize on t?
	    # meanz_all[i,:] = [sum(Yall[:,i,0])/NP, sum(Yall[:,i,1])/NP]
	    # meanx_all[i,:] = [sum(Yall[:,i,2])/NP, sum(Yall[:,i,3])/NP]
        
            # covz_all[i,:,:] = np.cov(Yall[:,i,0],Yall[:,i,1])
            # covx_all[i,:,:] = np.cov(Yall[:,i,2],Yall[:,i,3])
            # --        
            meanz_all[i], meanx_all[i], covz_all[i], covx_all[i] = statistics(Yall[:,i].T)
                 
            emitx_all[i] = np.sqrt(covx_all[i,0,0]*covx_all[i,1,1] - covx_all[i,0,1]*covx_all[i,0,1])
            betax_all[i] = meanz_all[i,1] * covx_all[i,0,0] / emitx_all[i]
            # alphax_all[i] = - meanz_all[i,1] * covx_all[i,0,1] / emitx_all[i]
            alphax_all[i] = - covx_all[i,0,1] / emitx_all[i]

    else:
    # Parallel version:
        pool = multiprocessing.Pool()
        stats = pool.map(statistics,(Yall[:,i].T for i in range(NT)))
        pool.close()
        pool.join()

        # Assign values (quick)
        for i in range(NT):            
            meanz_all[i], meanx_all[i], covz_all[i], covx_all[i] = stats[i]

            emitx_all[i] = np.sqrt(covx_all[i,0,0]*covx_all[i,1,1] - covx_all[i,0,1]*covx_all[i,0,1])
            betax_all[i] = meanz_all[i,1] * covx_all[i,0,0] / emitx_all[i]
            # alphax_all[i] = - meanz_all[i,1] * covx_all[i,0,1] / emitx_all[i]
            alphax_all[i] = - covx_all[i,0,1] / emitx_all[i]

            
    print('%.3f seconds ' % (clock.clock() - tclock))
    tclock = clock.clock()

    # --
    
    # Plotting with ROOT
    # ----------------------------------------------------------------- 
    
    ROOT.gROOT.SetBatch(1)

    print('Plotting beam moments...')
        
    # Ranges
    xmin  = -meanY[2]-6*sx
    xmax  = meanY[2]+6*sx
    pxmin = xmin * np.sqrt(K/meanY[1]) * meanY[1]
    pxmax = xmax * np.sqrt(K/meanY[1]) * meanY[1]
    # print(' xmin  = %f  xmax  = %f' % (xmin,xmax) ) 
    # print(' pxmin = %f  pxmax = %f' % (pxmin,pxmax) ) 

    if args.si :
        units = np.array([(1/kp)/um,1,(1/kp)/um,1/mc2])
        meanY = meanY * units
        rmsY  = rmsY * units
        sx = sx * (1/kp) / um
        sz = sz * (1/kp) / um

        xmin = xmin * (1/kp) / um
        xmax = xmax * (1/kp) / um

        time = time * (1/kp) / um

    zetamin = meanY[0] - 10 * np.sqrt(np.max(covz_all[:,0,0]))
    zetamax = meanY[0] + 4 * sz

    maxDpz = np.sqrt(np.max(covz_all[:,1,1]))
    pzmin = 1.0 * np.min(meanz_all[:,1]) - 4 * maxDpz 
    pzmax = 1.0 * np.max(meanz_all[:,1]) + 4 * maxDpz

        
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
    fWx = TF1('Wx',lambda x : wake.Wx(meanY[0],x[0]),xmin,xmax)
    fWx.SetNpx(1000)
    fWx.Draw()

    fEx = TF1('Ex',lambda x : wake.Ex(meanY[0],x[0]),xmin,xmax)
    fEx.SetNpx(1000)
    fEx.SetLineStyle(2)
    # fEx.Draw("same")

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

    # Plot bunch parameters
    lcolor = kMagenta + 2
    def DrawGraph(graph) :
        graph.SetLineWidth(2)
        graph.SetLineColor(lcolor)
        graph.SetMarkerSize(0.4)
        graph.SetMarkerStyle(20)
        graph.SetMarkerColor(lcolor)
        graph.Draw("apl")
        
    
    fmeanx = TGraph(NT,time,meanx_all[:,0].flatten())
    DrawGraph(fmeanx)
    C.Print(opath +'/meanx.pdf')
    
    fmeanpx = TGraph(NT,time,meanx_all[:,1].flatten())
    DrawGraph(fmeanpx)
    C.Print(opath +'/meanpx.pdf')

    frmsx = TGraph(NT,time,np.sqrt(covx_all[:,0,0]).flatten())
    DrawGraph(frmsx)
    C.Print(opath +'/rmsx.pdf')
    
    frmspx = TGraph(NT,time,np.sqrt(covx_all[:,1,1]).flatten())
    DrawGraph(frmspx)
    C.Print(opath +'/rmspx.pdf')

    fmeanz = TGraph(NT,time,meanz_all[:,0].flatten()-time)
    DrawGraph(fmeanz)
    C.Print(opath +'/meanz.pdf')
    
    fmeanpz = TGraph(NT,time,meanz_all[:,1].flatten())
    DrawGraph(fmeanpz)
    C.Print(opath +'/meanpz.pdf')

    frmsz = TGraph(NT,time,np.sqrt(covz_all[:,0,0]).flatten())
    DrawGraph(frmsz)
    C.Print(opath +'/rmsz.pdf')
    
    frmspz = TGraph(NT,time,np.sqrt(covz_all[:,1,1]).flatten())
    DrawGraph(frmspz)
    C.Print(opath +'/rmspz.pdf')

    frmspzrel = TGraph(NT,time,(100*np.sqrt(covz_all[:,1,1])/meanz_all[:,1]).flatten())
    DrawGraph(frmspzrel)
    C.Print(opath +'/rmspzrel.pdf')

    
    femitx = TGraph(NT,time,emitx_all[:].flatten())
    DrawGraph(femitx)
    C.Print(opath +'/emitx.pdf')

    fbetax = TGraph(NT,time,betax_all[:].flatten())
    DrawGraph(fbetax)
    C.Print(opath +'/betax.pdf')

    falphax = TGraph(NT,time,alphax_all[:].flatten())
    DrawGraph(falphax)
    C.Print(opath +'/alphax.pdf')

    
    print('%.3f seconds ' % (clock.clock() - tclock))
    tclock = clock.clock()

    # Draw trajectories of a subset of particles
    # -----------------------------------------------------
    
    Nsample = 100
    if Nsample > NP :
        Nsample = NP
        
    DP = int(NP/Nsample)
    # List of selected particle indexes
    plist = range(0, NP, DP)

    print('Plotting particle trajectories (Sample = %i part.)' % Nsample)

    if not os.path.exists(opath + '/tracks'):
	os.makedirs(opath + '/tracks')
    
    # Color palette
    palette = PPalette('bird')
    palette.SetPalette(kBird)
    palette.SetAlpha(0.6)

    def DrawTrack(graph) :
        graph.SetLineWidth(2)
        cindex = int((abs(Yall[i,0,2])/(meanY[2]+2*sx)) * (palette.GetNColors()-1))
        if cindex>palette.GetNColors()-1 :
            cindex = palette.GetNColors()-1
        graph.SetLineColor(palette.GetColorIndex(cindex))
        graph.SetMarkerSize(0.4)
        graph.SetMarkerStyle(20)
        graph.SetMarkerColor(palette.GetColorIndex(cindex))
        graph.Draw('l')
        
    # --
    frame = TH1F('frame','',10,time[0],time[-1])
    frame.GetYaxis().SetRangeUser(xmin,xmax)
    frame.Draw("axis")
    
    gxvst = []
    for i in plist:
        gxvst.append(TGraph(NT,time,Yall[i,:,2].flatten()))
        DrawTrack(gxvst[-1])        
    DrawFrame()
    C.Print(opath +'/tracks/xvst.pdf')

    # --
    frame2 = TH1F('frame2','',10,xmin,xmax)
    frame2.GetYaxis().SetRangeUser(pxmin,pxmax)
    frame2.Draw("axis")
    
    gpxvsx = []
    for i in plist:
        gpxvsx.append(TGraph(NT,Yall[i,:,2].flatten(),Yall[i,:,3].flatten()))
        DrawTrack(gpxvsx[-1])        
                
    DrawFrame()
    C.Print(opath +'/tracks/pxvsx.pdf')

    # --
    frame3 = TH1F('frame3','',10,zetamin,zetamax)
    frame3.GetYaxis().SetRangeUser(xmin,xmax)
    frame3.Draw("axis")

    gxvszeta = []
    for i in plist:
        gxvszeta.append(TGraph(NT,Yall[i,:,0].flatten()-time[:],Yall[i,:,2].flatten()))
        DrawTrack(gxvszeta[-1])        
                
    DrawFrame()
    C.Print(opath +'/tracks/xvszeta.pdf')

    # --
    frame4 = TH1F('frame4','',10,time[0],time[-1])
    frame4.GetYaxis().SetRangeUser(0.99,1.0)
    frame4.Draw("axis")

    gvzvst = []
    for i in plist:
        vz = Yall[i,:,1].flatten()/Gamma(Yall[i,:,1].flatten(),Yall[i,:,3].flatten())
        gvzvst.append(TGraph(NT,time,vz))
        DrawTrack(gvzvst[-1])        
               
    DrawFrame()
    C.Print(opath +'/tracks/vzvst.pdf')        

    # --
    frame5 = TH1F('frame5','',10,time[0],time[-1])
    frame5.GetYaxis().SetRangeUser(pzmin,pzmax)
    frame5.Draw("axis")

    ggammavst = []
    for i in plist:
        vgamma = Gamma(Yall[i,:,1].flatten(),Yall[i,:,3].flatten())
        ggammavst.append(TGraph(NT,time,vgamma))
        DrawTrack(ggammavst[-1])        
                
    DrawFrame()
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
        DT = int(NT/Nsample)

        tindex = range(0, NT, DT)
        for i in tindex:	
            hname = 'hpzvsz-%i' % i 
            hpzvsz = TH2F(hname,'',100,zetamin,zetamax,100,pzmin,pzmax)
            hname = 'hpxvsx-%i' % i 
            hpxvsx = TH2F(hname,'',100,xmin,xmax,100,pxmin,pxmax)
            hname = 'hxvszeta-%i' % i 
            hxvszeta = TH2F(hname,'',100,zetamin,zetamax,100,xmin,xmax)

            for j in range(NP):
                hpzvsz.Fill(Yall[j,i,0]-time[i],Yall[j,i,1])                  
                hpxvsx.Fill(Yall[j,i,2],Yall[j,i,3])                  
                hxvszeta.Fill(Yall[j,i,0]-time[i],Yall[j,i,2])                  

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
