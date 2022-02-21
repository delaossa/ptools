#!/usr/bin/env python
import os
from VisualPIC.DataReading.folder_scanners import OpenPMDFolderScanner
from opmd_viewer import OpenPMDTimeSeries
import argparse
import scipy.constants as ct
import numpy as np
import time as clock

import ROOT
from ROOT import PData, PUtils, PPalette
import putils as pu


def parse_args():
    # Command argument line parser
    parser = argparse.ArgumentParser(description='')

    parser.add_argument('sim', nargs='?', default=None, help='simulation name')
    parser.add_argument('spc', nargs='?', default=None, help='species name')

    parser.add_argument('-t', type=int, dest='time', default=0,
                        help='time step dump')
    parser.add_argument('-i', type=int, dest='itime', default=0,
                        help='initial time step dump')
    parser.add_argument('-f', type=int, dest='ftime', default=0,
                        help='final time step dump')
    parser.add_argument('-s', type=int, dest='stime', default=1,
                        help='advance step dump')
    parser.add_argument('-np', type=float, dest='np', default=0.0,
                        help='plasma density')
    parser.add_argument('-zmin', type=float, default=None,
                        help='set min longitudinal value')
    parser.add_argument('-zmax', type=float, default=None,
                        help='set max longitudinal value')
    parser.add_argument('-NS', type=int, default=10,
                        help='number of slices')
    parser.add_argument('-zsmin', type=float, default=None,
                        help='set min (slice) longitudinal value')
    parser.add_argument('-zsmax', type=float, default=None,
                        help='set max (slice) longitudinal value')
    parser.add_argument('-pzmin', type=float, default=None,
                        help='set min longitudinal momentum value')
    parser.add_argument('-pzmax', type=float, default=None,
                        help='set max longitudinal momentum value')
    parser.add_argument('-imax', type=float, default=None,
                        help='set max current value')
    parser.add_argument('-xmax', type=float, default=None,
                        help='set max horizontal value')
    parser.add_argument('-pxmax', type=float, default=None,
                        help='set max horizontal momentum value')
    parser.add_argument('-z0', type=float, default=0.0,
                        help='Set driver''s center position')
    parser.add_argument('-zp', type=float, default=0.0,
                        help='Set the start of plasma plateau')
    parser.add_argument('--comov', action='store_true', default=False,
                        help='comoving frame')
    parser.add_argument('opath', nargs='?', default=None,
                        help='Output folder')
    args = parser.parse_args()

    return args


def main():
    
    args = parse_args()

    # Units
    mc  = ct.m_e * ct.c
    mc2 = ct.m_e * ct.c**2
    MeV = 1E6 * ct.eV
    GeV = 1E9 * ct.eV
    GV  = GeV / ct.e
    um  = ct.micron
    mm  = 1000 * um

    spaunit = um
    sspaunit = '#mum'
    z0 = args.z0 * spaunit
    zp = args.zp * spaunit
    momzunit = GeV / ct.c
    smomzunit = 'GeV/c'
    momxunit = MeV / ct.c
    smomxunit = 'MeV/c'
    
    # plasma constants
    if args.np > 0.0:
        n0 = args.np * 1.E6
        kp  = np.sqrt(n0 * ct.e * ct.e / (ct.epsilon_0 * ct.m_e * ct.c**2))
        E0 = kp * mc2 / ct.e
        print()
        print(' plasma density = %.2e e/cm^3' % (n0 / 1e6))
        print(' E0 = %.2f GV/m' % (E0 / GV))
        
    # access simulation folder
    sim_path = './%s/lab_diags/hdf5' % (args.sim)
    if not os.path.exists(sim_path):
        sim_path = './%s/diags/hdf5' % (args.sim)
        if not os.path.exists(sim_path):
            print(' data not found. Exit.')
            return 0
        
    ts = OpenPMDTimeSeries(sim_path, check_all_files=False)
    it = ts.iterations.tolist().index(args.time)

    time = ts.t[it]
    print(' Time = %.2e seconds (%.2f) ' % (time, (kp * ct.c) * time))

    opmd_fs = OpenPMDFolderScanner()
    splist = opmd_fs.get_list_of_species(sim_path)
    
    def getspecies(spname, sp_list):
        for sp in sp_list:
            if sp.species_name == spname:
                return sp
        return None

    print('\n reading bunch distribution')
    tclock = clock.time()
    if args.spc is not None:
        spname = '%s' % args.spc
    else:
        spname = 'elec_N'

    bunch = getspecies(spname, splist)
    if bunch is None:
        print(' species ''%s'' not found. Exit.' % spname)
        return 0
    
    varlist = ['z', 'pz', 'x', 'px', 'y', 'py']
    varall = varlist[:]
    varall.append('q')
    datalist = bunch.get_data(args.time, varall)
    datarray = []
    for data in datalist:
        array, md = datalist[data]
        # print(md)
        datarray.append(array)

    datarray = np.asarray(datarray)

    Y = datarray[:-1, :]
    NP = Y.shape[1]
    weights = np.abs(datarray[-1, :])
    del datarray
    print(' number of particles = %i' % NP)
    
    # units
    Y[::2] = Y[::2] / spaunit
    Y[1] = Y[1] * mc / momzunit  # GeV/c
    Y[3::2] = Y[3::2] * mc / momxunit  # MeV/c
    weights = weights / 1e-12
    
    print(' .. done in %.3f s. ' % ((clock.time() - tclock)))
    tclock = clock.time()

    print('\n calculating statistical moments (whole bunch)')
    meanz, covz, meanx, covx, meany, covy = statistics(Y, weights)

    grel = meanz[1] * momzunit / mc
    shiftz = meanz[0]
    if args.comov:
        shiftz = ct.c * time / spaunit + z0 / spaunit
        
    Y[0] = Y[0] - shiftz  # Center z position respect to the mean
    meanz[0] = meanz[0] - shiftz
    
    charge = np.sum(weights)
    
    emitx, alphax, betax, gammax = Twisspars(covx, grel)
    emity, alphay, betay, gammay = Twisspars(covy, grel)
    emitx = emitx * momxunit / mc
    emity = emity * momxunit / mc
    
    tau = 2.35 * np.sqrt(covz[0][0]) * spaunit / ct.c / 1e-15
    
    print(' charge = %.3f pC' % (charge))
    print(' energy = %.3f MeV. rel. spread = %.3f %%' % (grel * mc2 / MeV, 100 * np.sqrt(covz[1][1]) * momzunit / mc / grel))
    print(' tau   = %.3f fs. Current (gaussian) = %.3f kA' % (tau, charge / tau))
    print(' emitx = %.3f um. betax = %.3f mm. gammax = %.3f mm-1. alphax = %.3f' % (emitx, betax * spaunit / mm, gammax * mm / spaunit, alphax))
    print(' emity = %.3f um. betay = %.3f mm. gammay = %.3f mm-1. alphay = %.3f' % (emity, betay * spaunit / mm, gammay * mm / spaunit, alphay))

    print(' .. done in %.3f s. ' % ((clock.time() - tclock)))
    tclock = clock.time()

    mins, maxs = getExtremes(Y)
    
    # ROOT style
    pu.SetStyle()
    ROOT.gStyle.SetNumberContours(99)
    ROOT.gStyle.SetHistFillStyle(1001)
    # ROOT.gStyle.SetLineWidth(4)
    ROOT.gStyle.SetHistLineWidth(4)
    ROOT.gStyle.SetFrameLineWidth(2)
    ROOT.gROOT.SetBatch(1)
    
    # scan histograms
    print('\n filling histograms')
    Nbins = 120
    hnames = ['histo_%s' % hname for hname in varlist]
    histos = []
    for i, hname in enumerate(hnames):
        histos.append(ROOT.TH1F(hname, '', Nbins, mins[i], maxs[i]))

    for ih, histo in enumerate(histos):
        # print(NP, Y[ih].shape, weights.shape)
        histo.FillN(NP, Y[ih], weights)
        # PUtils.FillHisto(histo, Y[ih], NP, weights)
        
    for ih, histo in enumerate(histos):
        mini, maxi = findLimits(histo, 0.1)  # This causes problems when mini = maxi
        # mini = ROOT.Double()
        # maxi = ROOT.Double()
        # PUtils.FindLimits(histo, mini, maxi, 0.1)  # This causes problems when mini = maxi
        delta = maxi - mini
        mins[ih] = mini - delta
        maxs[ih] = maxi + delta

    if args.zmin is not None:
        mins[0] = args.zmin

    if args.zmax is not None:
        maxs[0] = args.zmax

    if args.pzmin is not None:
        mins[1] = args.pzmin

    if args.pzmax is not None:
        maxs[1] = args.pzmax

    if args.xmax is not None:
        maxs[2] = np.abs(args.xmax)
        mins[2] = -np.abs(args.xmax)
        maxs[4] = np.abs(args.xmax)
        mins[4] = -np.abs(args.xmax)
        
    if args.pxmax is not None:
        maxs[3] = np.abs(args.pxmax)
        mins[3] = -np.abs(args.pxmax)
        maxs[5] = np.abs(args.pxmax)
        mins[5] = -np.abs(args.pxmax)

    for ih, histo in enumerate(histos):
        histos[ih].Delete()
        histos[ih] = ROOT.TH1F(hnames[ih], '', Nbins, mins[ih], maxs[ih])
        histos[ih].FillN(NP, Y[ih], weights)
        
    histo2D = []
    histo2D.append(ROOT.TH2F('hpzvsz', '', Nbins, mins[0], maxs[0], Nbins, mins[1], maxs[1]))
    histo2D.append(ROOT.TH2F('hpxvsx', '', Nbins, mins[2], maxs[2], Nbins, mins[3], maxs[3]))
    histo2D.append(ROOT.TH2F('hpyvsy', '', Nbins, mins[4], maxs[4], Nbins, mins[5], maxs[5]))
    histo2D.append(ROOT.TH2F('hxvsz', '', Nbins, mins[0], maxs[0], Nbins, mins[2], maxs[2]))
    histo2D.append(ROOT.TH2F('hyvsz', '', Nbins, mins[0], maxs[0], Nbins, mins[4], maxs[4]))
    histo2D.append(ROOT.TH2F('hyvsx', '', Nbins, mins[2], maxs[2], Nbins, mins[4], maxs[4]))

    histo2D[0].FillN(NP, Y[0], Y[1], weights)
    histo2D[1].FillN(NP, Y[2], Y[3], weights)
    histo2D[2].FillN(NP, Y[4], Y[5], weights)
    histo2D[3].FillN(NP, Y[0], Y[2], weights)
    histo2D[4].FillN(NP, Y[0], Y[4], weights)
    histo2D[5].FillN(NP, Y[2], Y[4], weights)

    # Calculate stats and Twiss parameters from histograms
    statsz = statsfromhist(histo2D[0])
    charge = statsz[0]
    pzmean = statsz[4] / charge
    pzrms = np.sqrt(statsz[5] / charge - pzmean**2)
    grel = pzmean * momzunit / mc
    grelrms = pzrms * momzunit / mc
    zmean = statsz[2] / charge
    zrms = np.sqrt(statsz[3] / charge - zmean**2)
    tau = 2.35 * zrms * spaunit / ct.c / 1e-15

    statsx = statsfromhist(histo2D[1])
    chargex = statsx[0]
    pxmean = statsx[4] / chargex
    pxrms = np.sqrt(statsx[5] / chargex - pxmean**2)
    xmean = statsx[2] / chargex
    xrms = np.sqrt(statsx[3] / chargex - xmean**2)
    emitx, alphax, betax, gammax = Twissfromhist(histo2D[1], grel)
    emitx = emitx * momxunit / mc
    alphax = alphax * mc / momxunit
    betax = betax * mc / momxunit
    gammax = gammax * momxunit / mc
    emity, alphay, betay, gammay = Twissfromhist(histo2D[2], grel)
    emity = emity * momxunit / mc
    alphay = alphay * mc / momxunit
    betay = betay * mc / momxunit
    gammay = gammay * momxunit / mc

    print()
    print(' charge = %.3f pC' % (charge))
    print(' energy = %.3f MeV. rel. spread = %.3f %%' % (grel * mc2 / MeV, 100 * grelrms / grel))
    print(' tau   = %.3f fs. Current = %.3f kA' % (tau, charge / tau))

    print(' emitx = %.3f um. betax = %.3f mm. gammax = %.3f mm-1. alphax = %.3f' % (emitx, betax * spaunit / mm, gammax * mm / spaunit, alphax))
    print(' emity = %.3f um. betay = %.3f mm. gammay = %.3f mm-1. alphay = %.3f' % (emity, betay * spaunit / mm, gammay * mm / spaunit, alphay))
    
    print(' .. done in %.3f s. ' % ((clock.time() - tclock)))
    print()

    tclock = clock.time()

    # Slice quantities
    # ----------------
    zsmin, zsmax = findLimits(histos[0], 0.1)
    NS = args.NS
    slidx = np.linspace(zsmin, zsmax, NS + 1)

    YL = []
    WL = []
    h2D_slc = []
    z_slc = []
    charge_slc = []
    grel_slc = []
    grelrms_slc = []
    emitx_slc = []
    emity_slc = []
    for i in range(NS):
        z_slc.append((slidx[i] + slidx[i + 1]) / 2.0)
        
        cond = (Y[0] > slidx[i]) & (Y[0] <= slidx[i + 1])
        YL.append(np.asarray([Y[i, cond] for i in range(Y.shape[0])]))
        WL.append(weights[cond])

        h2D_list = []
        for ih, h2D in enumerate(histo2D):
            h2D_list.append(h2D.Clone('%s_sl%i' % (h2D.GetName(), ih)))
            h2D_list[-1].Reset()
            
        h2D_list[0].FillN(YL[-1].shape[1], YL[-1][0], YL[-1][1], WL[-1])
        h2D_list[1].FillN(YL[-1].shape[1], YL[-1][2], YL[-1][3], WL[-1])
        h2D_list[2].FillN(YL[-1].shape[1], YL[-1][4], YL[-1][5], WL[-1])
        h2D_list[3].FillN(YL[-1].shape[1], YL[-1][0], YL[-1][2], WL[-1])
        h2D_list[4].FillN(YL[-1].shape[1], YL[-1][0], YL[-1][4], WL[-1])
        h2D_list[5].FillN(YL[-1].shape[1], YL[-1][2], YL[-1][4], WL[-1])

        h2D_slc.append(h2D_list)

        statsz = statsfromhist(h2D_slc[-1][0])
        chargez = statsz[0]
        charge_slc.append(chargez)
        pzmean = statsz[4] / chargez
        pzrms = np.sqrt(statsz[5] / chargez - pzmean**2)
        grel_slc.append(pzmean * momzunit / mc)
        grelrms_slc.append(pzrms * momzunit / mc)
        
        semit, _, _, _ = Twissfromhist(h2D_slc[-1][1], grel_slc[-1])
        semit = semit * momxunit / mc
        emitx_slc.append(semit)
        semit, _, _, _ = Twissfromhist(h2D_slc[-1][2], grel_slc[-1])
        semit = semit * momxunit / mc
        emity_slc.append(semit)
        
        '''
        print(' charge = %.3f pC' % (charge_slc[-1]))
        print(' energy = %.3f MeV. rel. spread = %.3f %%' % (grel_slc[-1] * mc2 / MeV, 100 * grelrms_slc[-1] / grel_slc[-1]))
        print(' emitx = %.3f um. betax = %.3f mm. gammax = %.3f mm-1. alphax = %.3f' % (emitx, betax * spaunit / mm, gammax * mm / spaunit, alphax))
        print(' emity = %.3f um. betay = %.3f mm. gammay = %.3f mm-1. alphay = %.3f' % (emity, betay * spaunit / mm, gammay * mm / spaunit, alphay))
        '''

    print()
    print('Total charge within slices = %.3f pC' % (np.sum(charge_slc)))
    print()
    
    # Plots and output data
    # ---------------------
    C = ROOT.TCanvas('C', '')
    C.cd()

    if args.opath:
        opath = args.opath
    else:
        opath = '%s/Plots/Bunch/%s' % (args.sim, args.spc)
    os.makedirs(opath, exist_ok=True)
        
    ofile = ROOT.TFile('%s/histos.root' % opath, 'RECREATE')
    for h in histos:
        C.Clear()
        h.Draw('hist LF')
        frame = pu.DrawFrame(width=4)
        C.Print('%s/%s.pdf' % (opath, h.GetName()))
        
        h.Write(h.GetName(), ROOT.TObject.kOverwrite)

    beampal = PPalette('beampal')
    beampal.SetPalette('elec0')

    for h in histo2D:
        
        C1 = pu.pcanvasPS([h], 'C_%s' % h.GetName())
        C1.cd(0)
        
        pad = ROOT.gROOT.FindObject('pad_0')
        pad.cd()

        textStat = None
        y1 = pad.GetBottomMargin()
        y2 = 1 - pad.GetTopMargin()
        x1 = pad.GetLeftMargin()
        x2 = 1 - pad.GetRightMargin()
        yran = y2 - y1
        xran = x2 - x1

        if 'hpxvsx' in h.GetName():
            h.GetXaxis().SetTitle('x [#mum]')
            h.GetYaxis().SetTitle('p_{x} [%s]' % smomxunit)
            h.GetZaxis().SetTitle('Charge [pC]')

            # textStat = ROOT.TPaveText(x1 + 0.02, y2 - 0.38, x1 + 0.28, y2 - 0.01, 'NDC')
            textStat = ROOT.TPaveText(x1 + 0.02 * xran, y2 - 0.5 * yran, x1 + 0.30 * xran, y2 - 0.02 * yran, 'NDC')
            ROOT.PGlobals.SetPaveTextStyle(textStat, 12)
            textStat.SetTextColor(ROOT.kGray + 3)
            textStat.SetTextFont(42)
            # textStat.SetTextSize(32)
            textStat.AddText('Q = %5.2f pC' % charge)
            textStat.AddText('#sigma_{x} = %5.2f #mum' % xrms)
            textStat.AddText('#sigma_{p_{x}} = %5.2f %s' % (pxrms, smomxunit))
            textStat.AddText('#varepsilon_{n,x} = %5.2f #mum' % emitx)
            textStat.AddText('#beta_{x} = %5.2f mm' % (betax * spaunit / mm))
            textStat.AddText('#gamma_{x} = %5.2f mm^{-1}' % (gammax * mm / spaunit))
            textStat.AddText('#alpha_{x} = %5.2f' % alphax)
                
        h.Draw('colz')

        if textStat is not None:
            textStat.Draw()
        
        frame = pu.DrawFrame(width=4)
        pframe = pu.DrawPaletteFrame(h)
        C1.Print('%s/%s.pdf' % (opath, h.GetName()))
        
        h.Write(h.GetName(), ROOT.TObject.kOverwrite)

    ofile.Close()

    # Side/top spatial view
    C1 = pu.pcanvasPS([histo2D[4], histo2D[3]], 'C_st')
    C1.cd(0)
        
    pad = ROOT.gROOT.FindObject('pad_1')
    pad.cd()
    histo2D[4].Draw('colz')
    pframet = pu.DrawPaletteFrame(histo2D[4])
    framet = pu.DrawFrame(width=4)

    pad = ROOT.gROOT.FindObject('pad_0')
    pad.cd()
    histo2D[3].Draw('colz')
    pframeb = pu.DrawPaletteFrame(histo2D[3])
    frameb = pu.DrawFrame(width=4)

    C1.Print('%s/hyxvsz.pdf' % opath)
        
    
    # ---------
    # longitudinal phasespace
    #
    hcurr = histos[0].Clone('hcurr')
    dt = hcurr.GetXaxis().GetBinWidth(0) * spaunit / ct.c / 1e-15
    hcurr.Scale(1.0 / dt)
    hcurr.GetXaxis().SetTitle('#zeta [#mum]')
    
    hpzvsz = histo2D[0].Clone('hpzvsz_cl')
    hpzvsz.GetYaxis().SetTitle('p_{z} [GeV/c]')
    hpzvsz.GetZaxis().SetTitle('Charge [pC]')

    # Canvas
    C2 = pu.pcanvasPS([hpzvsz, hcurr])

    # Top plot
    C2.cd(0)
    pad = ROOT.gROOT.FindObject('pad_1')
    pad.cd()

    # PhaseSpace
    pal = ROOT.PPalette('pspacepal')
    pal.SetPalette('electron0')
    hpzvsz.Draw('colz0 same')

    # projection
    Nx = histos[1].GetNbinsX()
    xarray = np.zeros(Nx)
    yarray = np.zeros(Nx)
    xMin = mins[0]
    xMax = xMin + (maxs[0] - mins[0]) * 0.2
    yMax = histos[1].GetMaximum()
    for i in range(Nx):
        yarray[i] = histos[1].GetBinCenter(i + 1)
        xarray[i] = xMin + histos[1].GetBinContent(i + 1) * (xMax - xMin) / yMax
    gPz = ROOT.TGraph(Nx, xarray, yarray)
    gPz.SetLineColor(ROOT.PGlobals.elecLine)
    gPz.SetLineWidth(3)
    gPz.SetFillStyle(1001)
    gPz.SetFillColor(ROOT.PGlobals.elecFill)
    gPz.Draw('F')
    gPz.Draw('L')
    
    lpzmean1 = ROOT.TLine(hpzvsz.GetXaxis().GetXmin(), meanz[1],
                          hpzvsz.GetXaxis().GetXmax(), meanz[1])
    lpzmean1.SetLineColor(ROOT.kGray + 2)
    lpzmean1.SetLineStyle(2)
    lpzmean1.Draw()

    lzmean1 = ROOT.TLine(meanz[0], hpzvsz.GetYaxis().GetXmin(),
                         meanz[0], hpzvsz.GetYaxis().GetXmax())
    lzmean1.SetLineColor(ROOT.kGray + 2)
    lzmean1.SetLineStyle(2)
    lzmean1.Draw()

    textTime = ROOT.TPaveText(0.55, 0.76, 0.84, 0.9, 'NDC')
    ROOT.PGlobals.SetPaveTextStyle(textTime, 32)
    textTime.SetTextFont(43)
    textTime.SetTextSize(22)
    textTime.SetTextColor(ROOT.kGray + 2)
    textTime.AddText('z = %5.1f mm' % ((ct.c * time + z0 - zp) / 1e-3))
    textTime.Draw()

    textInfo = ROOT.TPaveText(0.55, 0.25, 0.84, 0.75, 'NDC')
    ROOT.PGlobals.SetPaveTextStyle(textInfo, 32)
    textInfo.SetTextColor(ROOT.kGray + 2)
    textInfo.SetTextFont(43)
    textInfo.SetTextSize(22)

    textInfo.AddText('Q = %5.2f pC' % charge)
    textInfo.AddText('#sigma_{#zeta} = %5.2f #mum' % np.sqrt(covz[0][0]))
    textInfo.AddText('#sigma_{#gamma}/#LT#gamma#GT = %4.1f %%' % (100 * grelrms / grel))
    textInfo.AddText('#varepsilon_{n,x} = %5.2f #mum' % emitx)
    textInfo.AddText('#varepsilon_{n,y} = %5.2f #mum' % emity)

    textInfo.Draw()

    textMom = ROOT.TPaveText(0.55, 0.07, 0.84, 0.17, 'NDC')
    ROOT.PGlobals.SetPaveTextStyle(textMom, 32)
    textMom.SetTextColor(ROOT.kGray + 3)
    textMom.SetTextFont(63)
    textMom.SetTextSize(22)
    textMom.AddText('#LTp_{z}#GT = %5.2f GeV/c' % meanz[1])
    textMom.Draw()
    
    frame1 = pu.DrawFrame(width=4)
    pframe1 = pu.DrawPaletteFrame(hpzvsz)
    
    # Bottom plot
    C2.cd(0)
    pad = ROOT.gROOT.FindObject('pad_0')
    pad.cd()

    imax = 1.2 * hcurr.GetMaximum()
    if args.imax is not None:
        imax = args.imax
        
    hcurr.GetYaxis().SetRangeUser(0.001, imax)
        
    hcurr.Draw('hist LF')

    lzmean0 = ROOT.TLine(meanz[0], 0.0, meanz[0], imax)
    lzmean0.SetLineColor(ROOT.kGray + 2)
    lzmean0.SetLineStyle(2)
    lzmean0.Draw()

    markerSize = 1.0
    lineWidth  = 2
    
    gEmity = ROOT.TGraph(len(z_slc), np.asarray(z_slc), np.asarray(emity_slc))
    gEmity.SetMarkerColor(ROOT.kGray + 2)
    gEmity.SetMarkerSize(markerSize)
    gEmity.SetLineWidth(lineWidth)
    gEmity.SetLineColor(ROOT.kGray + 2)
    gEmity.Draw('PL')

    gEmitx = ROOT.TGraph(len(z_slc), np.asarray(z_slc), np.asarray(emitx_slc))
    gEmitx.SetMarkerColor(ROOT.kGray + 3)
    gEmitx.SetMarkerSize(markerSize)
    gEmitx.SetLineWidth(lineWidth)
    gEmitx.SetLineColor(ROOT.kGray + 3)
    gEmitx.Draw('PL')
        
    gErms = ROOT.TGraph(len(z_slc), np.asarray(z_slc), 100 * np.asarray(grelrms_slc) / np.asarray(grel_slc))
    gErms.SetMarkerColor(ROOT.kOrange + 10)
    gErms.SetMarkerSize(markerSize)
    gErms.SetLineWidth(lineWidth)
    gErms.SetLineColor(ROOT.kOrange + 10)
    gErms.Draw('PL')

    Leg = ROOT.TLegend(0.6, 0.6, 1 - ROOT.gPad.GetRightMargin() - 0.04 - 0.03, 0.95)
    ROOT.PGlobals.SetPaveStyle(Leg)
    Leg.SetTextAlign(12)
    Leg.SetTextColor(ROOT.kGray + 3)
    Leg.SetTextFont(43)
    Leg.SetTextSize(20)
    Leg.SetLineColor(1)
    Leg.SetBorderSize(0)
    Leg.SetFillColor(0)
    Leg.SetFillStyle(1001)
    Leg.SetFillStyle(0)

    Leg.AddEntry(hcurr, 'Current [kA]', 'l')
    Leg.AddEntry(gErms, 'E. spread [%]', 'pl')
    Leg.AddEntry(gEmitx, 'Emitt. x [#mum]', 'pl')
    Leg.AddEntry(gEmity, 'Emitt. y [#mum]', 'pl')

    Leg.Draw()
        
    frame0 = pu.DrawFrame(width=4)
    
    C2.Print('%s/%s.pdf' % (opath, C2.GetName()))
    
    return 0


# Statistics in vectorial form
def statistics(Y, weights=None):
    meanz = [np.average(Y[0], weights=weights), np.average(Y[1], weights=weights)]
    covz = np.cov(Y[0], Y[1], aweights=weights)
    meanx = [np.average(Y[2], weights=weights), np.average(Y[3], weights=weights)]
    covx = np.cov(Y[2], Y[3], aweights=weights)
    meany = [np.average(Y[4], weights=weights), np.average(Y[5], weights=weights)]
    covy = np.cov(Y[4], Y[5], aweights=weights)
    
    return meanz, covz, meanx, covx, meany, covy


def Twisspars(covx, grel):
    emitx = np.sqrt(covx[0, 0] * covx[1, 1] - covx[0, 1] * covx[0, 1])
    alphax = - covx[0, 1] / emitx
    betax = grel * covx[0, 0] / emitx
    gammax = covx[1, 1] / emitx / grel

    return emitx, alphax, betax, gammax


def Twissfromhist(hist, grel):
    stats = statsfromhist(hist)

    meanx = stats[2] / stats[0]
    rmsx2 = stats[3] / stats[0] - meanx**2
    meany = stats[4] / stats[0]
    rmsy2 = stats[5] / stats[0] - meany**2
    rmsxy = stats[6] / stats[0] - stats[2] * stats[4] / stats[0]**2

    emitx = np.sqrt(rmsx2 * rmsy2 - rmsxy**2)

    alphax = -rmsxy / emitx
    betax = grel * rmsx2 / emitx
    gammax = rmsy2 / emitx / grel
    
    return emitx, alphax, betax, gammax


def statsfromhist(hist):
    stats = np.zeros(7, dtype=np.double)
    hist.GetStats(stats)
    return stats


def getExtremes(Y):

    mins = []
    maxs = []
    for X in Y:
        mins.append(np.amin(X))
        maxs.append(np.amax(X))

    return mins, maxs


def findLimits(h, factor=0.1):
    maxval = h.GetBinContent(h.GetMaximumBin())
    for i in range(h.GetNbinsX()):
        binValue = h.GetBinContent(i + 1)
        if binValue > maxval * factor:
            xmin = h.GetBinCenter(i + 1)
            break
  
    for i in range(h.GetNbinsX(), -1, -1):
        binValue = h.GetBinContent(i + 1)
        if binValue > maxval * factor:
            xmax = h.GetBinCenter(i + 1)
            break

    return xmin, xmax


if __name__ == '__main__':
    main()
