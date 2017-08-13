#!/usr/bin/env python

import scipy.integrate as integrate
import scipy.constants as ct
import numpy as np
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

    parser.add_argument('name', nargs='?', default='betaout', help='Simulation name (output folder)')
    parser.add_argument('-NP', type=int, dest='NP', default=1000, help='Number of macroparticles')
    parser.add_argument('--ps', action='store_true', default=1, help='Plot phasespace snapshots')
    parser.add_argument('--png', action='store_true', default=0, help='Plot phasespace snapshots in png (slow!)')
    parser.add_argument('--mp', action='store_true', default=0, help='For parallel processing')
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

    # Normalization constants
    n0  = 1.0E23 # m^-3
    kp  = np.sqrt( n0 * ct.e * ct.e / (ct.epsilon_0 * ct.m_e) ) / ct.c
    mc2 = ct.m_e * ct.c * ct.c
    MeV = ct.mega * ct.eV
    GeV = 1000 * MeV
    E0  = mc2 * kp
    um  = ct.micron
    mm  = 1000 * um

    # print('Speed of light %f m/s' % ct.c)
    # print('Electron mass  = %f MeV' % (mc2 / MeV))
    print('Density   = %.1e cm^-3' % (1E-6 * n0))
    print('Skindepth = %.2f um' % ((1/kp)/um))
    print('E0        = %.2f GV/m' % (E0/(GeV)))

    # Wakefields definition
    # -------------------------------
    Ez0 = -0.263
    #Ez0 = 0.0
    K = 0.185
    S = 0.409
    #S = 0.0
    Sk = 0.119
    #Sk = 0.0
    rb = 10.0
    Drho = 1.0
    wake = wk.WakeBlowout(Ez0,K,S,Sk,rb,Drho)
    # -------------------------------

    a0 = 0.8
    L = np.sqrt(2)
    w0 = 2.976
    nb0 = a0**2/2
    sbz = L/np.sqrt(2)
    sbx = w0/2
    #wake = wk.WakeLinearSimple(nb0,sbz,sbx)
    #wake = wk.WakeLinearTimon(a0,w0)

    # -------------------------------------

    # Bunch definition
    # -------------------------------

    # Number of particles
    NP = args.NP
    print('Bunch definition : N = %i particles' % NP)

    # Bunch definition in SI units    

    # Normalization array
    norma = np.array([1/kp,mc2,1/kp,mc2])
    
    # Beam: Initial beam moments (means and covariances)   
    #                [       z,            pz,       x,   px ]
    meanY0 = np.array([  -(np.pi/2)/kp,    1*GeV,  1.0*um,  0.0 ])
    #meanY0 = np.array([  -36*um,    1000*MeV,  0.0*um,  0.0 ])
    rmsY0  = np.array([  3.0*um, 0.0065*meanY0[1], 0.0*um,  0.0 ]) 

    # Normalize
    meanY0 = meanY0 / norma
    rmsY0  = rmsY0  / norma      

    # Gamma
    G0 = meanY0[1]
    
    # Emittance (normalized) 
    emit0 = kp * 0.30 * um

    K = wake.Kx(meanY0[0],0) # Focusing strength (normalized)
    # print(' K = %f' % K)
    omegaB = np.sqrt(K/G0)   # Betatron frequency
    lambdaB = 2*np.pi/omegaB # Betatron wavelength 
    
    # Matched beta 
    betaM = 1/omegaB    

    # Beta
    beta0 = 5.38 * betaM 
    # beta0 = betaM 

    # sigma_x (beam width)
    rmsY0[2] = np.sqrt(beta0 * emit0 / G0)    
    
    # sigma_px
    rmsY0[3]  = emit0 / rmsY0[2]

    # emitM = G0 * np.power(rmsY0[2],2) / betaM  # matched emittance
    # print('Matched emittance = %f um' % (emitM/um) )        
    # emit0 = emitM                
    # beta0 =  G0 * np.power(rmsY0[2],2) / emit0
            
    # Second moments, covariance matrix initialization.
    sz0  = rmsY0[0]
    spz0 = rmsY0[1] 
    sx0  = rmsY0[2]
    spx0 = rmsY0[3] 
        
    covY0  = [[sz0**2, 0.0   ,   0.0,     0.0],
             [  0.0, spz0**2,   0.0,     0.0],
             [  0.0, 0.0    ,sx0**2,     0.0],
             [  0.0, 0.0    ,   0.0, spx0**2]]

    # Initial distribution function f(Y0) = f(z,pz,x,px)
    Y0 = np.random.multivariate_normal(meanY0, covY0, NP)

    print('Beam parameters: ')
    print('Gamma     = %5.1f'  % meanY0[1])
    print('DeltaG    = %.3f %%' % (100*rmsY0[1]/meanY0[1]) )
    print('Emittance = %.3f um (%.3f)' % (emit0/(kp*um),emit0) )
    print('Betax     = %.3f mm (%.2f)' % (beta0/(kp*mm),beta0) )
    print('Sigmax    = %.3f um (%.3f)' % (rmsY0[2]/(kp*um),rmsY0[2]) )
    print('Sigmaz    = %.3f um (%.3f)' % (rmsY0[0]/(kp*um),rmsY0[0]) )

    print('------')
    print('lambdaB   = %.3f mm (%.2f)' % (lambdaB/(kp*mm),lambdaB) )
    print('Beta M    = %.3f mm (%.2f)' % (betaM/(kp*mm),betaM) )

    S = wake.dEz(meanY0[0],0)
    Sk = wake.dKx(meanY0[0],0)
    Ldeco = lambdaB * K / (rmsY0[0]*Sk)
    #Ldeco = (np.pi/(omegaB*rmsY0[0])) * (K/Sk)
    print('L,dc      = %.3f mm (%.2f)' % (Ldeco/(kp*mm),Ldeco))

    print('Ez = %.3f  GV/m      (%.3f)' % ((E0/GeV) * wake.Ez(meanY0[0],0.0),wake.Ez(meanY0[0],0.0)))
    print('S  = %.3f (GV/m)/um  (%.3f)' % ((E0/GeV)*(kp*um)*S,S)) 
    print('K  = %.3f (GV/m)/um  (%.3f)' % ((E0/GeV)*(kp*um)*wake.Kx(meanY0[0],0.0),wake.Kx(meanY0[0],0.0)))
    print('Sk = %.3f (GV/m)/um2 (%.3f)' % ((E0/GeV)*(kp*um)*(kp*um)*Sk,Sk)) 


    # Plasma profile definition
    # ----------------------------------------- 
    plasma = Plasma(5 * lambdaB , 0.05 * lambdaB)
    # plasma = Plasma(2400 , 2 * lambdaB)

    # Time grid
    # --------------------------------
    tmin = 0.0
    dt   = 0.2 / omegaB
    tmax = 5.0 * lambdaB
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


    else :
    # Parallel version:
        pintegra = partial(integra,wake,plasma,time,dY)
        pool = multiprocessing.Pool()
        NCPU = multiprocessing.cpu_count()
        print('Number of CPUs = %i' % NCPU)
    
        Ylist = pool.map(pintegra,(Y0[i] for i in range(NP)) )
        pool.close()
        pool.join()

        for i in range(NP):
            Yall[i] = Ylist[i]

                
    # ----------------------------------------------------
        
    print('Done in %.3f s. ' % ((clock.clock() - tclock)))
    tclock = clock.clock()

    # ---------------------------------
	                
    print('Now analysing...')

    # Units vector
    units = np.array([um,mc2,um,mc2])
    
    # Normalization vector
    norma = np.array([1/kp,mc2,1/kp,mc2])

    # to SI units
    for i in range(NP):
        Yall[i,:] = Yall[i,:] * (norma / units)   

    meanY0 = meanY0 * (norma/units)
    rmsY0  = rmsY0 * (norma/units)    
    time = time * (norma[0]/units[0])

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

    opath = 'betaout'
    if args.name :
        opath = args.name
    
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
