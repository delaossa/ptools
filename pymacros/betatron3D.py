#!/usr/bin/env python

import scipy.integrate as integrate
import scipy.constants as ct
import numpy as np
import h5py
import os, sys, argparse
import time as clock
import multiprocessing
from functools import partial
import ROOT
from ROOT import PData, PGlobals, PPalette, PUnits, TF1, TH2F, TH1F, TGraph, TCanvas, TLine, TColor, TStyle, gStyle, gPad
import wakefields as wk

def parse_args():

    # Command argument line parser 
    parser = argparse.ArgumentParser(description='Betatron.py: A particle-in-wakefield trajectory solver.')
    parser.add_argument('filename', nargs='?', default='none', help='file name')
    parser.add_argument('opath', nargs='?', default='betaout', help='Output folder')
    parser.add_argument('-NP', type=int, dest='NP', default=1000, help='Number of macroparticles')
    parser.add_argument('-n0', type=float, dest='n0', default=1E18, help='Plasma density in cm^-3')
    parser.add_argument('--ps', action='store_true', default=1, help='Plot phasespace snapshots')
    parser.add_argument('--png', action='store_true', default=0, help='Plot phasespace snapshots in png (slow!)')
    parser.add_argument('--man', action='store_true', default=0, help='Print help')
    
    args = parser.parse_args()

    if ('none' in args.filename) :
        parser.print_help()
        sys.exit(0)

    if args.man :
        parser.print_help()
        sys.exit(0)
        
    return args

    
# Relativistic factor
def Gamma(pz,px,py):
    return np.sqrt( 1 + (pz*pz + px*px + py*py) )


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
    
    # Y is the vector of particle coordinates : Y = (z, pz, x, px, y, py)
    # dY returns the time diferential of that vector.

    kp = np.sqrt(plasma.n(Y[0]))

    zeta = kp*(Y[0]-t)
    r   = kp*np.sqrt(Y[2]**2 + Y[4]**2)
    phi = np.arctan2(Y[4],Y[2])

    return [ Y[1]/Gamma(Y[1],Y[3],Y[5]),
             
             - kp * ( wake.Ez(zeta,r) + ( wake.Ex(zeta,r) - wake.Wx(zeta,r) ) * ( Y[3]*np.cos(phi) + Y[5]*np.sin(phi) ) / Gamma(Y[1],Y[3],Y[5])  ), 
             
             Y[3]/Gamma(Y[1],Y[3],Y[5]),
             
             - kp * ( wake.Wx(zeta,r) * Y[1]/Gamma(Y[1],Y[3],Y[5]) + wake.Ex(zeta,r) * (1 - Y[1]/Gamma(Y[1],Y[3],Y[5]) ) ) * np.cos(phi),

             Y[5]/Gamma(Y[1],Y[3],Y[5]),
             
             - kp * ( wake.Wx(zeta,r) * Y[1]/Gamma(Y[1],Y[3],Y[5]) + wake.Ex(zeta,r) * (1 - Y[1]/Gamma(Y[1],Y[3],Y[5]) ) ) * np.sin(phi) ]


def integra(wake,plasma,time,dY,Y) :
    return integrate.odeint(dY, Y, time, args=(wake,plasma))

def statistics(Y) :
    NP = len(Y[0])

    meanz = [sum(Y[0])/NP, sum(Y[1])/NP]
    meanx = [sum(Y[2])/NP, sum(Y[3])/NP]
    meany = [sum(Y[4])/NP, sum(Y[5])/NP]
    
    covz = np.cov(Y[0],Y[1])
    covx = np.cov(Y[2],Y[3])
    covy = np.cov(Y[4],Y[5])

    return meanz,meanx,meany,covz,covx,covy


def main():

    args = parse_args()

    # Normalization constants
    # n0  = 1.0E23 # m^-3
    um  = ct.micron
    mm  = 1000 * um
    n0  = args.n0 * (1/(10*mm))**3 # m^-3
    kp  = np.sqrt( n0 * ct.e * ct.e / (ct.epsilon_0 * ct.m_e) ) / ct.c
    mc2 = ct.m_e * ct.c * ct.c
    MeV = ct.mega * ct.eV
    GeV = 1000 * MeV
    E0  = mc2 * kp

    # print('Speed of light %f m/s' % ct.c)
    # print('Electron mass  = %f MeV' % (mc2 / MeV))
    print('Density   = %.1e cm^-3' % (1E-6 * n0))
    print('Skindepth = %.2f um' % ((1/kp)/um))
    print('E0        = %.2f GV/m' % (E0/(GeV)))

    # Wakefields definition
    # -------------------------------
    Ez0 = -0.8
    #Ez0 = 0.0
    K = 0.4
    #K = 0.0
    S = 0.28
    #S = 0.0
    Sk = 0.0
    rb = 10.0
    Drho = 1.0
    zeta0 = -np.pi/2
    wake = wk.WakeBlowout(Ez0,K,S,Sk,rb,Drho)

    
    # -------------------------------
    # Bunch definition
    # -------------------------------

    # Bunch definition in SI units    
    
    filename = args.filename

    tclock = clock.time()
    print('Reading file = %s' % (filename))
    if filename.find('.npz')>-1 :
        data = np.load(filename)

        Z0  = data.f.y
        X0  = data.f.x
        Y0  = data.f.z
        PZ0 = data.f.py
        PX0 = data.f.px
        PY0 = data.f.pz
        W0  = data.f.w

        ccenter = 1.403424e-4 # meters


    elif filename.find('.h5')>-1 :
        fileh5 = h5py.File(filename)

        Z0 = fileh5['x1'][:]
        X0 = fileh5['x2'][:]
        Y0 = fileh5['x3'][:]
        PZ0 = fileh5['p1'][:]
        PX0 = fileh5['p2'][:]
        PY0 = fileh5['p3'][:]
        W0 = fileh5['q'][:]

        ccenter = 0.0

        
    NP0   = len(X0)
    print('Total number of particles =  %i' % NP0)
    meanz0 = sum(Z0)/NP0
    meanpz0 = sum(PZ0)/NP0
    #meanpx = sum(PX0)/NP0
    #meanpy = sum(PY0)/NP0
    
    print('Done in %.3f s. ' % ((clock.time() - tclock)))

    # Number of particles
    NP = args.NP
    # NP = 1000
    DP = int(NP0/NP)
    plist = range(0,NP0,DP)
    NP = len(plist)

    # Select a sub-sample and change to polar coodinates (quasi-3D)
    YP0   = np.empty((NP,6))
    YPW0  = np.empty((NP))
    print('Bunch definition : N = %i particles selected from %i' % (NP,NP0))
    
    i = 0
    for index in plist :
        
        if filename.find('.npz')>-1 :
            YP0[i,0] = kp * Z0[index]
            YP0[i,2] = kp * (X0[index] - ccenter)
            YP0[i,4] = kp * (Y0[index] - ccenter)
        elif filename.find('.h5')>-1 :
            YP0[i,0] = Z0[index]
            YP0[i,2] = (X0[index] - ccenter)
            YP0[i,4] = (Y0[index] - ccenter)
            
        YP0[i,1] = PZ0[index]
        YP0[i,3] = PX0[index]
        YP0[i,5] = PY0[index]
        
        # particle weight
        YPW0[i]  = W0[index] * float(NP0/NP)
        i = i+1

    meanz0 = sum(YP0[:,0])/NP
    YP0[:,0] = YP0[:,0] - meanz0 + zeta0
    
    G0 = meanpz0             # Gamma
    K = wake.Kx(zeta0,0)     # Focusing strength (normalized)
    omegaB = np.sqrt(K/G0)   # Betatron frequency
    lambdaB = 2*np.pi/omegaB # Betatron wavelength 

    
    # Plasma profile definition
    # ----------------------------------------- 
    #plasma = Plasma(5 * lambdaB , 0.05 * lambdaB)
    # plasma = Plasma(2400 , 2 * lambdaB)
    plasma = Plasma(0.1*mm*kp, 0.220*mm*kp)
    
    # Time grid
    # --------------------------------
    tmin = 0.0
    dt   = 0.2 / omegaB
    #tmax = 5 * lambdaB
    tmax = 1*mm*kp
    time = np.arange(tmin,tmax,dt)
    NT = len(time)
    # print(' tmax = %f ' % tmax ) 
    print('Time grid : t = (%.2f,%.2f) NT = %i  dt = %.2f' % (time[0],time[-1],NT,dt))   
    
    # ODE solver
    # ----------

    tclock = clock.time()
    print('Solving equations of motion (ODEs) ...')

    # Phase space (Particles and coordinates) 
    # ---------------------------------------
    
    # Matrix storing particle trajectories :
    # (z,pz,x,px,y,py) for each particle at each time. 
    Yall = np.empty((NP,NT,6))

    # Parallel version:
    pintegra = partial(integra,wake,plasma,time,dY)
    pool = multiprocessing.Pool()
    NCPU = multiprocessing.cpu_count()
    print('Number of CPUs = %i' % NCPU)
    
    Ylist = pool.map(pintegra,(YP0[i] for i in range(NP)) )
    pool.close()
    pool.join()

    for i in range(NP):
        Yall[i] = Ylist[i]
        
    # ----------------------------------------------------
        
    print('Done in %.3f s. ' % ((clock.time() - tclock)))
    tclock = clock.time()

    # ---------------------------------
	                
    print('Now analysing...')

    # Units vector
    units = np.array([um,mc2,um,mc2,um,mc2])
    
    # Normalization vector
    norma = np.array([1/kp,mc2,1/kp,mc2,1/kp,mc2])

    # to SI units
    for i in range(NP):
        Yall[i,:] = Yall[i,:] * (norma / units)   

    # Statistics vs time
    # ---------------------------------

    print('Extracting statistical moments...')

    # Statistical moments 
    meanz_all = np.empty((NT,2))
    covz_all  = np.empty((NT,2,2))
    meanx_all = np.empty((NT,2))
    covx_all  = np.empty((NT,2,2))
    meany_all = np.empty((NT,2))
    covy_all  = np.empty((NT,2,2))
    
    # Twiss parameters
    emitx_all = np.empty(NT)
    betax_all = np.empty(NT)
    alphax_all = np.empty(NT)
    emity_all = np.empty(NT)
    betay_all = np.empty(NT)
    alphay_all = np.empty(NT)
    
    # Parallel version:
    pool = multiprocessing.Pool()
    stats = pool.map(statistics,(Yall[:,i].T for i in range(NT)))
    pool.close()
    pool.join()

    # Assign values (quick)
    for i in range(NT):            
        meanz_all[i], meanx_all[i], meany_all[i], covz_all[i], covx_all[i], covy_all[i] = stats[i]

        emitx_all[i] = np.sqrt(covx_all[i,0,0]*covx_all[i,1,1] - covx_all[i,0,1]*covx_all[i,0,1])
        betax_all[i] = meanz_all[i,1] * covx_all[i,0,0] / emitx_all[i]
        alphax_all[i] = - covx_all[i,0,1] / emitx_all[i]

        emity_all[i] = np.sqrt(covy_all[i,0,0]*covy_all[i,1,1] - covy_all[i,0,1]*covy_all[i,0,1])
        betay_all[i] = meanz_all[i,1] * covy_all[i,0,0] / emity_all[i]
        alphay_all[i] = - covy_all[i,0,1] / emity_all[i]

    meanY0 = np.array((meanz_all[0,0],meanz_all[0,1],meanx_all[0,0],meanx_all[0,1],meany_all[0,0],meany_all[0,1]))
    rmsY0  = np.array((covz_all[0,0,0],covz_all[0,1,1],covx_all[0,0,0],covx_all[0,1,1],covy_all[0,0,0],covy_all[0,1,1]))
    rmsY0  = np.sqrt(rmsY0)
    
    time = time * (norma[0]/units[0])

    print('%.3f seconds ' % (clock.time() - tclock))
    tclock = clock.time()

    # --
    
    # Plotting with ROOT
    # ----------------------------------------------------------------- 
    
    ROOT.gROOT.SetBatch(1)

    print('Plotting beam moments...')

    # Ranges
    xmin  = -meanY0[2]-6*rmsY0[2]
    xmax  = meanY0[2]+6*rmsY0[2]
    #pxmin = xmin * np.sqrt(K/meanY0[1]) * meanY0[1]
    #pxmax = xmax * np.sqrt(K/meanY0[1]) * meanY0[1]
    maxspx = np.sqrt(np.max(covx_all[:,1,1]))
    pxmin = 1.0 * np.min(meanx_all[:,1]) - 4 * maxspx 
    pxmax = 1.0 * np.max(meanx_all[:,1]) + 4 * maxspx

    zetamin = meanY0[0] - 10 * np.sqrt(np.max(covz_all[:,0,0]))
    zetamax = meanY0[0] + 4 * rmsY0[0]

    maxspz = np.sqrt(np.max(covz_all[:,1,1]))
    pzmin = 1.0 * np.min(meanz_all[:,1]) - 4 * maxspz 
    pzmax = 1.0 * np.max(meanz_all[:,1]) + 4 * maxspz

    # print(' xmin  = %f  xmax  = %f' % (xmin,xmax) ) 
    # print(' pxmin = %f  pxmax = %f' % (pxmin,pxmax) ) 

    # -----
    
        
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

    opath = args.opath    
    if not os.path.exists(opath):
        os.makedirs(opath)

    # Plot Fields
    if not os.path.exists(opath + '/fields'):
        os.makedirs(opath + '/fields')

    fWx = TF1('Wx',lambda x : wake.Wx(kp*meanY0[0]*units[0],kp*x[0]*units[2]),xmin,xmax)
    fWx.SetNpx(1000)
    maxWx = fWx.GetMaximum()
    minWx = fWx.GetMinimum()

    fEx = TF1('Ex',lambda x : wake.Ex(kp*meanY0[0]*units[0],kp*x[0]*units[2]),xmin,xmax)
    fEx.SetNpx(1000)
    fEx.SetLineStyle(2)
    maxEx = fEx.GetMaximum()
    minEx = fEx.GetMinimum()
    
    hFrame = TH1F('hFrame','',10,xmin,xmax)
    maxV = np.amax([maxWx,maxEx])
    minV = np.amin([minWx,minEx])
    hFrame.GetYaxis().SetRangeUser(minV - 0.2*(maxV-minV),maxV + 0.2*(maxV-minV))
    hFrame.Draw("axis")
    
    fWx.Draw("same")
    fEx.Draw("same")
    C.Print(opath +'/fields/Wx&Ex.pdf')

    fEz = TF1('Ez',lambda x : wake.Ez(kp*x[0]*units[0],0.0),zetamin,zetamax)
    maxEz = fEz.GetMaximum()
    minEz = fEz.GetMinimum()

    fKx = TF1('Kx',lambda x : wake.Kx(kp*x[0]*units[0],0.0),zetamin,zetamax)
    fKx.SetLineStyle(2)
    maxKx = fKx.GetMaximum()
    minKx = fKx.GetMinimum()

    hFrame = TH1F('hFrame','',10,zetamin,zetamax)
    maxV = np.amax([maxEz,maxKx])
    minV = np.amin([minEz,minKx])
    hFrame.GetYaxis().SetRangeUser(minV - 0.2*(maxV-minV),maxV + 0.2*(maxV-minV))
    hFrame.Draw("axis")
    
    fEz.Draw("same")
    fKx.Draw("same")
    
    C.Print(opath +'/fields/Ez&Kx.pdf')

    fProf = TF1('Prof',lambda x : plasma.n(kp*x[0]*units[0]),time[0],time[-1])
    fProf.SetNpx(1000)
    fProf.Draw()
    C.Print(opath +'/fields/DenProf.pdf')

    # Plot bunch parameters
    lcolor = ROOT.kMagenta + 2
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


    # y direction
    femity = TGraph(NT,time,emity_all[:].flatten())
    DrawGraph(femity)
    C.Print(opath +'/emity.pdf')

    fbetay = TGraph(NT,time,betay_all[:].flatten())
    DrawGraph(fbetay)
    C.Print(opath +'/betay.pdf')

    
    print('%.3f seconds ' % (clock.time() - tclock))
    tclock = clock.time()

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
    palette.SetPalette(ROOT.kBird)
    palette.SetAlpha(0.6)

    def DrawTrack(graph) :
        graph.SetLineWidth(2)
        cindex = int((abs(Yall[i,0,2])/(meanY0[2]+2*rmsY0[2])) * (palette.GetNColors()-1))
        if cindex>palette.GetNColors()-1 :
            cindex = palette.GetNColors()-1
        graph.SetLineColor(palette.GetColorIndex(cindex))
        graph.SetMarkerSize(0.4)
        graph.SetMarkerStyle(20)
        graph.SetMarkerColor(palette.GetColorIndex(cindex))
        graph.Draw('l')

    def DrawFrame() :
        gPad.Update()
        gPad.RedrawAxis('g')
        gPad.RedrawAxis()

        l = TLine()
        l.SetLineColor(ROOT.kBlack)
        l.SetLineWidth(2)
        l.DrawLine(gPad.GetUxmin(), gPad.GetUymax(), gPad.GetUxmax(), gPad.GetUymax())
        l.DrawLine(gPad.GetUxmax(), gPad.GetUymin(), gPad.GetUxmax(), gPad.GetUymax())
        l.DrawLine(gPad.GetUxmin(), gPad.GetUymin(), gPad.GetUxmin(), gPad.GetUymax())
        l.DrawLine(gPad.GetUxmin(), gPad.GetUymin(), gPad.GetUxmax(), gPad.GetUymin())
        

        
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
        vz = Yall[i,:,1].flatten()/Gamma(Yall[i,:,1].flatten(),Yall[i,:,3].flatten(),Yall[i,:,5].flatten())
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
        vgamma = Gamma(Yall[i,:,1].flatten(),Yall[i,:,3].flatten(),Yall[i,:,5].flatten())
        ggammavst.append(TGraph(NT,time,vgamma))
        DrawTrack(ggammavst[-1])        
                
    DrawFrame()
    C.Print(opath +'/tracks/gammavst.pdf')
        

    print('%.3f seconds ' % (clock.time() - tclock))
    tclock = clock.time()


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


        print('%.3f seconds ' % (clock.time() - tclock))
        tclock = clock.time()


    if filename.find('.npz')>-1 :
        fileout = filename.replace('.npz','.raw')
    elif filename.find('.h5')>-1 :
        fileout = filename.replace('.h5','.raw')

    fout = open(fileout, 'w')
    for i in range(NP):
        fout.write('%14e  %14e  %14e  %14e  %14e  %14e  %14e\n' %
                   ((Yall[i,NT-1,0]-meanz_all[NT-1,0])*units[0],Yall[i,NT-1,2]*units[2],Yall[i,NT-1,4]*units[4],Yall[i,NT-1,1],Yall[i,NT-1,3],Yall[i,NT-1,5],-(ct.e/ct.pico)*YPW0[i]))

    fout.close()

        
    print('Output at : %s' % opath +'/' )
    print('Ciao.')


if __name__ == '__main__':
    main()
