#!/usr/bin/env python

import os, argparse
import numpy as np
import scipy.constants as ct
import ROOT
from ROOT import gStyle, TCanvas, gPad, TFile, TPave, PPalette, PGlobals, TLine, TGaxis, TMarker, TLine, TGraph


def parse_args():

    # Command argument line parser 
    parser = argparse.ArgumentParser(description='')

    parser.add_argument('sim', nargs='?', default=None, help='simulation name')
    parser.add_argument('-index', type=int, dest='index', default=0, help='particle index')
    parser.add_argument('-t', type=int, dest='time', default=0, help='time step dump')
    parser.add_argument('-i', type=int, dest='itime', default=0, help='initial time step dump')
    parser.add_argument('-f', type=int, dest='ftime', default=0, help='final time step dump')
    parser.add_argument('-s', type=int, dest='stime', default=1, help='advance step dump')
    parser.add_argument('-dmax', type=float, dest='dmax', default=1.0, help='maximum spectral 2d density')
    parser.add_argument('-smax', type=float, dest='smax', default=10.0, help='maximum spectral 1d density')
    parser.add_argument('--logz', action='store_true', default=0, help='logaritmic z scale')
    parser.add_argument('opath', nargs='?', default=None, help='Output folder')
    
    args = parser.parse_args()
        
    return args


def DrawFrame() :
    gPad.Update()
    gPad.RedrawAxis('g')
    gPad.RedrawAxis()

    hframe = TPave(gPad.GetUxmin(),gPad.GetUymin(),gPad.GetUxmax(),gPad.GetUymax(),2,'')
    hframe.SetFillStyle(0)
    hframe.SetLineWidth(2)
    hframe.SetLineColor(ROOT.kBlack)
    hframe.SetShadowColor(0)
    hframe.Draw()
    
    return hframe
    
def DrawPaletteFrame(histo) :
    gPad.Update()

    y1 = gPad.GetBottomMargin()
    y2 = 1 - gPad.GetTopMargin()
    x2 = 1 - gPad.GetRightMargin()
    pFrame = TPave((x2 + 0.01),y1,(x2 + 0.03),y2,1,'NDCL')
    pFrame.SetFillStyle(0)
    pFrame.SetLineWidth(2)
    pFrame.SetLineColor(ROOT.kBlack)
    pFrame.SetShadowColor(0)

    palette = histo.GetListOfFunctions().FindObject('palette')
    if palette :
      palette.SetY1NDC(y1)
      palette.SetY2NDC(y2)
      palette.SetX1NDC(x2 + 0.01)
      palette.SetX2NDC(x2 + 0.03)
      palette.SetBorderSize(1)
      palette.SetLineColor(1)
      
      pFrame.Draw()

    return pFrame


def Draw1dhist(hist,ymin,ymax,smin,smax,opts='hist C ax1d same') :
    
    axismax = ymin + (ymax-ymin)/2.0
    slope    = (axismax - ymin)/(smax-smin)
    
    color = 46 
    hclone = hist.Clone('hclone')
    hclone.SetLineColor(color)
    hclone.SetLineWidth(3)

    if opts.find('LF2')>=0 :
        hclone.SetFillColor(ROOT.kRed-10)
        hclone.SetLineWidth(0)

        
    for i in range(0,hclone.GetNbinsX()+1) :
        value = hclone.GetBinContent(i)
        hclone.SetBinContent(i,slope * (value - smin) + ymin)

    xmin = hclone.GetXaxis().GetXmin()
    xmax = hclone.GetXaxis().GetXmax()
    xpos = xmax - (xmax-xmin) * 0.10

    axis = None
    if opts.find('ax1d')>=0 :
        axis = TGaxis(xpos,ymin,xpos,axismax,smin,smax,205,'S+L')
        axis.SetLineWidth(2)
        axis.SetLineColor(color)
        axis.SetLabelColor(color)
        axis.SetLabelFont(43)
        axis.SetLabelSize(22)
        axis.SetLabelOffset(0.005)
        axis.SetTitleColor(color)
        axis.SetTitleFont(43)
        axis.SetTitleSize(22)
        axis.SetTitleOffset(0.4)
        axis.SetTickSize(0.03)
        axis.ChangeLabel(1,-1,0.)
        
        axis.SetTitle('pC/MeV')
        #axis.CenterTitle()
        #axis.ChangeLabel(1,-1,-1,-1,-1,-1,'')
        #axis.SetMaxDigits(2)
        axis.Draw()

    hclone.Draw(opts)

    return hclone, axis


def GetCrossings(hist,baseline=0.0,xmin=None,xmax=None) :

    cross = [0.0]
    extrm = [0.0]
    maxim = 0.0

    if xmin == None :
        xmin = hist.GetXaxis().GetXmin()

    if xmax == None :
        xmax = hist.GetXaxis().GetXmax()

    nbins = hist.GetNbinsX()
    for i in range(nbins,0,-1) :
        x1 = hist.GetBinCenter(i)
        v1 = hist.GetBinContent(i)-baseline
        if x1<xmin or x1>xmax :
            continue

        x2 = hist.GetBinCenter(i-1)
        v2 = hist.GetBinContent(i-1)-baseline

        # print ('(%.2f,%.2f)  (%.2f,%.2f)'%(x1,v1+baseline,x2,v2+baseline))
        
        if v1 * v2 >= 0 :
            if np.absolute(v2+baseline)>=np.absolute(maxim) : 
                maxim = v2+baseline
                extrm[-1] = x2
        else :
            xcross =  -v1 * ( (x2-x1)/(v2-v1) ) + x1
            cross[-1] = xcross
            # print(' CROSS FOUND at %.4f -> value = %.4f' % (cross[-1],hist.GetBinContent(hist.FindBin(cross[-1]))))
            # print(' EXTRM FOUND at %.4f -> value = %.4f' % (extrm[-1],hist.GetBinContent(hist.FindBin(extrm[-1]))))

            maxim = 0.0
            cross.append(0.0)
            extrm.append(0.0)

    ncross = len(cross)
    if (ncross%2>0) :
        del cross[-1]
        del extrm[-1]
    
    return cross, extrm 
            

def GetMoments(hist,xmin=None,xmax=None) :

    if xmin == None :
        xmin = hist.GetXaxis().GetXmin()

    if xmax == None :
        xmax = hist.GetXaxis().GetXmax()

    hcut = hist.Clone('hcut')
    hcut.Reset()
    
    imin = hist.FindBin(xmin)
    imax = hist.FindBin(xmax)
    dx   = hist.GetBinWidth(imin)
    nbins = hist.GetNbinsX()
    norm  = 0.0
    xsum  = 0.0
    xsum2 = 0.0
    for i in range(imin,imax+1,1) :
        norm  = norm  + hist.GetBinContent(i) * dx
        xsum  = xsum  + hist.GetBinCenter(i) * hist.GetBinContent(i) * dx 
        xsum2 = xsum2 + hist.GetBinCenter(i) * hist.GetBinCenter(i) * hist.GetBinContent(i) * dx
        hcut.SetBinContent(i,hist.GetBinContent(i))

    hcut.ResetStats()

    charge = norm
    xmean  = xsum / norm
    xrms   = np.sqrt( (xsum2/norm) - xmean**2 )

    return hcut, charge, xmean, xrms

def distance(timestep) :
    n0=8e24
    kp = np.sqrt( n0 * ct.e * ct.e / (ct.epsilon_0 * ct.m_e) ) / ct.c

    dt   = 0.03270 * 68.0
    time = timestep * dt
    
    zmaxini    = 0.0
    zbeamini   = -4.3
    zplasmaini = 482.17
    
    z =  (time - zmaxini - zbeamini - zplasmaini)/kp
    return z
    
def density(z) :

    n0 = 8e24
    kp = np.sqrt( n0 * ct.e * ct.e / (ct.epsilon_0 * ct.m_e) ) / ct.c
    um = ct.micron

    sigma = 300 * um
    # zstart = 905.90 * um
    zstart = 0.0
    lplateau = 1000.0 * um
    zend = zstart + lplateau 
    
    if z<zstart :
        func = np.exp(-np.power(z-zstart,2.0)/(2*np.power(sigma,2.0)))
    elif z<zend :
        func = 1.0
    else :
        func = np.exp(-np.power(z-zend,2.0)/(2*np.power(sigma,2.0)))
        
    return 0.5 * n0 * func


def SpectralAnalysis(args) : 

    # Canvas setup
    sizex = 1024
    sizey = 320
    C = TCanvas('C','Spectrum',sizex,sizey)
    # C.SetFillStyle(4000)
    #
    gPad.SetLeftMargin(0.10)
    gPad.SetRightMargin(0.14)
    gPad.SetTopMargin(0.05)
    gPad.SetBottomMargin(0.25)
    gPad.SetTickx(1)
    gPad.SetTicky(1)
    if args.logz :
        gPad.SetLogz(1)

    filename = './%s/Plots/Spectrum/Spectrum-%s_%i.root' % (args.sim,args.sim,args.time)
    rfile = TFile(filename,'READ')

    hSpec2D = rfile.Get('hSpectrum2D')
    hSpec2D.SetMaximum(args.dmax)
    
    hSpec2D.Draw('colz')
    pframe = DrawPaletteFrame(hSpec2D)
   
    xaxis = hSpec2D.GetXaxis()
    NxBin = xaxis.GetNbins()
    xMin  = xaxis.GetXmin()
    xMax  = xaxis.GetXmax()
    dx    = (xMax-xMin)/NxBin
    yaxis = hSpec2D.GetYaxis()
    NyBin = yaxis.GetNbins()
    yMin  = yaxis.GetXmin()
    yMax  = yaxis.GetXmax()
    dy    = (yMax-yMin)/NyBin
    
    hSpec1D = hSpec2D.ProjectionX("hEneProj",1,NyBin)
    hSpec1D.Scale(dy)

    maxbin = hSpec1D.GetMaximumBin()
    Maximum = hSpec1D.GetBinContent(maxbin)

    cross, extrm =  GetCrossings(hSpec1D,baseline=Maximum/2.0)

    Np = len(cross)
    # print(' Number of crossings = %i'%Np)
    Npeaks = len(cross)//2
    # print(' Number of peaks     = %i'%Npeaks)
    
    smin = 0.0
    smax = args.smax
    axismax = yMin + (yMax-yMin)/2.0
    slope    = (axismax - yMin)/(smax-smin)

    color = 46
    marker = []
    line = [] 

    moments = {'time' : args.time}
    moments['z'] = distance(args.time)

    if Npeaks>0 :
        i = 0

        # print(' Maximum  = %.4f  value = %.4f' % (extrm[i+1],hSpec1D.GetBinContent(hSpec1D.FindBin(extrm[i+1]))))
        moments['emax'] = extrm[i+1]
        moments['max'] = hSpec1D.GetBinContent(hSpec1D.FindBin(extrm[i+1]))

        # print(' Interval = (%.4f,%.4f)'%(cross[i],cross[i+1]))
        moments['fwhm'] = cross[i]-cross[i+1]

        value = hSpec1D.GetBinContent(hSpec1D.FindBin(cross[i]))
        value = slope * (value - smin) + yMin
        marker.append(TMarker(cross[i],value,20))
        marker[-1].SetMarkerColor(color)
        line.append(TLine(cross[i],yMin,cross[i],value))
        line[-1].SetLineColor(color)
        line[-1].SetLineStyle(2)
        
        value = hSpec1D.GetBinContent(hSpec1D.FindBin(cross[i+1]))
        value = slope * (value - smin) + yMin
        marker.append(TMarker(cross[i+1],value,20))
        marker[-1].SetMarkerColor(color)
        line.append(TLine(cross[i+1],yMin,cross[i+1],value))
        line[-1].SetLineColor(color)
        line[-1].SetLineStyle(2)
        
        value = hSpec1D.GetBinContent(hSpec1D.FindBin(extrm[i+1]))
        value = slope * (value - smin) + yMin
        marker.append(TMarker(extrm[i+1],value,20) )
        marker[-1].SetMarkerColor(color)
        line.append(TLine(extrm[i+1],yMin,extrm[i+1],value))
        line[-1].SetLineColor(color)
        line[-1].SetLineStyle(2)

        hcut, charge, xmean, xrms = GetMoments(hSpec1D,cross[i+1],cross[i])
        # print(' Charge = %.4f  Mean = %.4f  Rms %.4f' % (charge, xmean, xrms))

        moments['charge'] = charge
        moments['emean'] = xmean
        moments['erms'] = xrms

        stats = np.zeros(4,dtype=np.dtype('Float64'))
        hcut.GetStats(stats)
        # print(stats)
        xmeancut = stats[2]/stats[0]
        xrmscut  = np.sqrt( stats[3]/stats[0] - xmeancut * xmeancut )
        chargecut = stats[0] * hcut.GetBinWidth(0)
        # print(' Charge = %.4f  Mean = %.4f  Rms %.4f' % (chargecut,xmeancut,xrmscut))
    

    hclonecut, axiscut = Draw1dhist(hcut,yMin,yMax,smin,smax,'hist LF2 same')
    
    # for i,mark in enumerate(marker) :
    #    mark.Draw()
    
    for i,lin in enumerate(line) :
        lin.Draw()

    hclone, axis = Draw1dhist(hSpec1D,yMin,yMax,smin,smax)

    hframe = DrawFrame()
    
    if args.opath :
        opath = args.opath    
        if not os.path.exists(opath):
            os.makedirs(opath)
        ofile = '%s/Spectrumpy-%s_%i.pdf'%(args.opath,args.sim,args.time)
    else :
        ofile = '%s/Plots/Spectrum/Spectrumpy-%s_%i.pdf'%(args.sim,args.sim,args.time)

    C.Print(ofile)
    C.Clear()
    
    return moments

def main():

    ROOT.gROOT.SetBatch(1)

    home = os.getcwd()    
    args = parse_args()

    # Style
    PGlobals.Initialize()
    
    gStyle.SetOptTitle(0)
    gStyle.SetOptStat(0)
    gStyle.SetOptFit(0)

    gStyle.SetLineWidth(2)
    gStyle.SetLineColor(ROOT.kAzure+2)
    gStyle.SetHistLineWidth(2)
    gStyle.SetMarkerColor(ROOT.kAzure+2)
    gStyle.SetMarkerSize(1)
    
    gStyle.SetNumberContours(255)

    gStyle.SetTickLength(0.03,'x')
    gStyle.SetTitleOffset(0.9,'x')
    gStyle.SetTitleOffset(0.6,'y')
    gStyle.SetTitleSize(32,'y')
    gStyle.SetLabelOffset(0.005,'y')
    
    # color palette
    pal = PPalette('mypal')
    pal.SetPalette('spectrum1')

    listmoms = []
    if args.ftime > args.itime :

        for it in range(args.itime,args.ftime,args.stime) :
            args.time = it
            moments = SpectralAnalysis(args)

            listmoms.append(moments)
            
    else :
        moments = SpectralAnalysis(args)
        listmoms.append(moments)

    # print(listmoms[-1])

    npoints = len(listmoms)
    times = np.zeros(npoints)
    zs = np.zeros(npoints)
    maxs = np.zeros(npoints)
    emaxs = np.zeros(npoints)
    fwhms = np.zeros(npoints)
    charges = np.zeros(npoints)
    emeans = np.zeros(npoints)
    ermss = np.zeros(npoints)
    for i, mom in enumerate(listmoms) :
        times[i] = mom['time']
        zs[i] = mom['z']
        maxs[i] = mom['max']
        emaxs[i] = mom['emax']
        fwhms[i] = mom['fwhm']
        charges[i] = mom['charge']
        emeans[i] = mom['emean']
        ermss[i] = mom['erms']
        
    # sprint(fwhms)

    # Canvas setup
    sizex = 1024
    sizey = 320
    C = TCanvas('C','Evolution',sizex,sizey)
    # C.SetFillStyle(4000)
    #
    gPad.SetLeftMargin(0.14)
    gPad.SetRightMargin(0.10)
    gPad.SetTopMargin(0.05)
    gPad.SetBottomMargin(0.25)
    gPad.SetTickx(1)
    gPad.SetTicky(1)

    if args.opath :
        opath = args.opath    
        if not os.path.exists(opath):
            os.makedirs(opath)
    else :
        args.opath = '%s/Plots/Spectrum' % args.sim

    
    vecden = np.vectorize(density)
    dens = vecden(zs)/1.0e24
    zs = zs/ct.micron


    gStyle.SetNdivisions(510,'x')
    dopt = 'alp'
    
    denvst = TGraph(npoints,zs,dens)
    denvst.GetXaxis().SetTitle('distance [#mum]')
    denvst.GetYaxis().SetTitle('n_{p} [10^{18} cm^{-3}]')
    denvst.Draw('ac')
        
    ofile = '%s/density-%s.pdf'%(args.opath,args.sim)
    C.Print(ofile)
    C.Clear()

    maxvst = TGraph(npoints,zs,maxs)
    maxvst.GetXaxis().SetTitle('distance [#mum]')
    maxvst.GetYaxis().SetTitle('SD_{max} [pC/MeV]')
    maxvst.Draw(dopt)

    ofile = '%s/sdensitypeak-%s.pdf'%(args.opath,args.sim)
    C.Print(ofile)
    C.Clear()

    emaxvst = TGraph(npoints,zs,emaxs)
    emaxvst.GetXaxis().SetTitle('distance [#mum]')
    emaxvst.GetYaxis().SetTitle('E_{max} [MeV]')
    emaxvst.Draw(dopt)

    ofile = '%s/epeak-%s.pdf'%(args.opath,args.sim)
    C.Print(ofile)
    C.Clear()

    # sbrights = charges/(1000 * ermss/emeans)
    sbrights = charges/(1000 * fwhms/emeans)
    sbrightvst = TGraph(npoints,zs,sbrights)
    sbrightvst.GetXaxis().SetTitle('distance [#mum]')
    sbrightvst.GetYaxis().SetTitle('SB_{fwhm} [pC/%0.1]')
    sbrightvst.Draw(dopt)

    ofile = '%s/sbrightfwhm-%s.pdf'%(args.opath,args.sim)
    C.Print(ofile)
    C.Clear()

    fwhmvst = TGraph(npoints,zs,fwhms)
    fwhmvst.GetXaxis().SetTitle('distance [#mum]')
    fwhmvst.GetYaxis().SetTitle('#DeltaE_{fwhm} [MeV]')
    fwhmvst.Draw(dopt)

    ofile = '%s/efwhm-%s.pdf'%(args.opath,args.sim)
    C.Print(ofile)
    C.Clear()
    
    emeanvst = TGraph(npoints,zs,emeans)
    emeanvst.GetXaxis().SetTitle('distance [#mum]')
    emeanvst.GetYaxis().SetTitle('#LTE#GT_{fwhm} [MeV]')
    emeanvst.Draw(dopt)

    ofile = '%s/emeanfwhm-%s.pdf'%(args.opath,args.sim)
    C.Print(ofile)
    C.Clear()

    chargevst = TGraph(npoints,zs,charges)
    chargevst.GetXaxis().SetTitle('distance [#mum]')
    chargevst.GetYaxis().SetTitle('Q_{fwhm} [pC]')
    chargevst.Draw(dopt)

    ofile = '%s/chargefwhm-%s.pdf'%(args.opath,args.sim)
    C.Print(ofile)
    C.Clear()

    
    ofile = '%s/momevol-%s.root'%(args.opath,args.sim)
    rfile = TFile(ofile,'recreate')

    denvst.Write('denvst',ROOT.TObject.kOverwrite)
    chargevst.Write('chargevst',ROOT.TObject.kOverwrite)
    maxvst.Write('maxvst',ROOT.TObject.kOverwrite)
    emaxvst.Write('emaxvst',ROOT.TObject.kOverwrite)
    emeanvst.Write('emeanvst',ROOT.TObject.kOverwrite)
    fwhmvst.Write('efwhmvst',ROOT.TObject.kOverwrite)
    sbrightvst.Write('sbrightvst',ROOT.TObject.kOverwrite)

    rfile.Close()
    

if __name__ == '__main__':
    main()
