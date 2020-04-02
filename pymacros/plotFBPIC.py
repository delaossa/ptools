#!/usr/bin/env python
import os
from visualpic.data_reading.folder_scanners import OpenPMDFolderScanner
from openpmd_viewer import OpenPMDTimeSeries
import argparse
import scipy.constants as ct
import numpy as np
import time as clock

import ROOT
from ROOT import PData, PUtils
import putils as pu


def parse_args():
    # Command argument line parser
    parser = argparse.ArgumentParser(description='')

    parser.add_argument('sim', nargs='?', default=None, help='simulation name')

    parser.add_argument('-t', type=int, dest='time', default=0,
                        help='time step dump')
    parser.add_argument('-i', type=int, dest='itime', default=0,
                        help='initial time step dump')
    parser.add_argument('-f', type=int, dest='ftime', default=0,
                        help='final time step dump')
    parser.add_argument('-s', type=int, dest='stime', default=1,
                        help='advance step dump')
    parser.add_argument('-np', type=float, dest='pden', default=0.0,
                        help='plasma density')
    parser.add_argument('-nj', type=int, dest='nj', default=1,
                        help='joint first <nj> fields')
    parser.add_argument('-zmin', type=float, default=None,
                        help='')
    parser.add_argument('-zmax', type=float, default=None,
                        help='')
    parser.add_argument('-xmin', type=float, default=None,
                        help='')
    parser.add_argument('-xmax', type=float, default=None,
                        help='')
    parser.add_argument('-z0', type=float, default=0.0,
                        help='Set driver''s center position')
    parser.add_argument('-zp', type=float, default=0.0,
                        help='Set the start of plasma plateau')
    parser.add_argument('-zoom', type=float, dest='zoom', default=1.0,
                        help='')
    parser.add_argument('-mode', type=int, dest='mode', default=None,
                        help='')
    parser.add_argument('opath', nargs='?', default=None,
                        help='Output folder')
    parser.add_argument('-al', type=float, default=1.0,
                        help='laser alpha')
    parser.add_argument('--comov', action='store_true', default=False,
                        help='comoving frame')
    parser.add_argument('--logz', action='store_true', default=False,
                        help='logaritmic z scale')
    parser.add_argument('--auto', action='store_true', default=False,
                        help='automatic (ROOT) z range')
    parser.add_argument('--png', action='store_true', default=0,
                        help='png ')
    parser.add_argument('--fix', action='store_true', default=False,
                        help='fix some issue')
    parser.add_argument('--fix2', action='store_true', default=False,
                        help='fix some issue2')
    parser.add_argument('-rfx', type=int, dest='rfx', default=1,
                        help='decrease x resolution by the given factor')
    parser.add_argument('-rfy', type=int, dest='rfy', default=1,
                        help='decrease y resolution by the given factor')
    args = parser.parse_args()

    return args


def main():
    
    args = parse_args()

    # Units
    mc2  = ct.m_e * ct.c**2
    GV  = 1E9 * ct.eV / ct.e
    um  = ct.micron
    k0  = 2 * ct.pi / (0.8 * um)
    n0  = 1.0
    kp  = 1.0
    E0 = 1.0

    spaunit = um
    sspaunit = '#mum'
    z0 = args.z0 * spaunit
    zp = args.zp * spaunit
    
    # plasma constants
    if args.pden > 0.0:
        n0 = args.pden * 1.E6
        kp  = np.sqrt(n0 * ct.e * ct.e / (ct.epsilon_0 * ct.m_e * ct.c**2))
        E0 = kp * mc2 / ct.e
        print(' plasma density = %.2e e/cm^3' % (n0 / 1e6))
        print(' E0 = %.2f GV/m' % (E0 / GV))
        
    # access simulation folder
    sim_path = './%s/lab_diags/hdf5' % (args.sim)
    if not os.path.exists(sim_path):
        sim_path = './%s/diags/hdf5' % (args.sim)
        if not os.path.exists(sim_path):
            print(' data not found. Exit.')
            return 0

    print(' simulation path = %s' % sim_path)
    ts = OpenPMDTimeSeries(sim_path, check_all_files=False)
    it = ts.iterations.tolist().index(args.time)
    
    time = ts.t[it]
    print(' Time = %.2e seconds (%.2f) ' % (time, (kp * ct.c) * time))

    # ROOT
    ROOT.gROOT.SetBatch(1)
    # ROOT style
    ROOT.PGlobals.Initialize()
    pu.SetStyle()
    # ROOT.gStyle.SetPadTopMargin(0.7)

    # from VisualPIC:
    opmd_fs = OpenPMDFolderScanner()
    fld_list = opmd_fs.get_list_of_fields(sim_path)
    
    # ---
    
    # Set fields to plot
    Fields = [pu.Field('rho', -0.01e20, 1e20, label='#rho_{p}/#rho_{0}', norm=-ct.e * n0, texec='beamPalette->cd();', auto=True),
              pu.Field('Ez', -5, 5, label='E_{z} [MV/m]', norm=1e6, texec='fieldPalette->cd();'),
              pu.Field('Ex', -0.5, 0.5, label='E_{x} [GV/m]', norm=1e9, texec='fieldPalette->cd();')]
    #  pu.Field('By', -2.5, 2.5, label='B_{y} [GV/m/c]', norm=1e9 / ct.c, texec='fieldPalette->cd();')]

    if Fields[0].auto is not True:
        pu.SetPlasmaPalette(Fields[0].Min, Fields[0].Max, 0.0, args.logz)
    
    histolist = []
    
    def getfield(fname, fld_list):
        for fld in fld_list:
            if fld.field_name == fname:
                return fld
        return None

    print('\nReading field data and filling histograms')
    tclock = clock.time()
    
    for Field in Fields:
        fname = Field.name
        if fname == 'a0':
            fname = 'Ex'
            
        fld = getfield(fname, fld_list)
        if fld is None:
            print(' field not found: %s' % fname)
            break

        if args.mode is not None:
            data, metadata = fld.get_data(args.time, m=args.mode)
        else:
            data, metadata = fld.get_data(args.time, m='all')
            
        rfx = args.rfx
        rfy = args.rfy
        if rfx > 1 or rfy > 1:
            data = data[::rfy, ::rfx]
                 
        data = data / Field.norm
        
        if 'rho' in Field.name:
            if args.fix:
                data = np.where(data < 1e-10, 0.0, data)
            if args.fix2:
                data = data - 1

        if 'a0' in Field.name:
            data = np.abs(data)
            Field.texec = Field.texec + (' laserPalette->SetAlpha(%.2f);' % args.al)
            
        hname = 'h_%s' % Field.name
        rarray = metadata['axis']['r']['array'] / spaunit
        zarray = metadata['axis']['z']['array'] / spaunit
        nxbins = data.shape[1]
        nybins = data.shape[0]
        zmin = zarray[0]
        zmax = zarray[-1]
        xmin = rarray[0]
        xmax = rarray[-1]
        if args.comov:
            zmin = zmin - ct.c * time / spaunit - z0 / spaunit
            zmax = zmax - ct.c * time / spaunit - z0 / spaunit

        print('Histo name: %s -> (Nx=%i, Ny=%i)' % (hname, nxbins, nybins))
        histolist.append(ROOT.TH2F(hname, '', nxbins, zmin, zmax,
                                   nybins, xmin, xmax))
        # for i in range(nxbins):
        #    for j in range(nybins):
        #        histolist[-1].SetBinContent(i + 1, j + 1, data[j][i])
        PUtils.SetBinsHisto2D(histolist[-1], data.flatten())
        
        if not Field.auto:
            histolist[-1].GetZaxis().SetRangeUser(Field.Min, Field.Max)
        histolist[-1].GetZaxis().SetTitle(Field.label)
        histolist[-1].GetZaxis().CenterTitle()
            
        if args.zmin is not None:
            zmin = args.zmin
        if args.zmax is not None:
            zmax = args.zmax
        histolist[-1].GetXaxis().SetRangeUser(zmin, zmax)
        histolist[-1].GetXaxis().SetTitle('#zeta [%s]' % sspaunit)
        histolist[-1].GetXaxis().CenterTitle()

        if args.xmin is not None:
            xmin = args.xmin
        if args.xmax is not None:
            xmax = args.xmax
        if args.zoom > 1.0:
            xran = (xmax - xmin) / args.zoom
            xmid = (xmin + xmax) / 2.0
            xmin = xmid - xran / 2.0
            xmax = xmid + xran / 2.0
        histolist[-1].GetYaxis().SetRangeUser(xmin, xmax)
        histolist[-1].GetYaxis().SetTitle('x [%s]' % sspaunit)
        histolist[-1].GetYaxis().CenterTitle()
        
        # histolist[-1].ResetStats()
        
    print(' .. done in %.3f s. ' % ((clock.time() - tclock)))
    tclock = clock.time()

    # Save histograms to file
    if args.opath:
        opath = args.opath
    else:
        opath = '%s/Plots/Snapshots' % (args.sim)
    os.makedirs(opath, exist_ok=True)

    ofile = '%s/Snapshot-%s_%i.root' % (opath, args.sim, args.time)
    rfile = ROOT.TFile(ofile, 'RECREATE')
    for histo in histolist:
        if isinstance(histo, list):
            for hist in histo:
                hist.Write(hist.GetName(), ROOT.TObject.kOverwrite)
        else:
            histo.Write(histo.GetName(), ROOT.TObject.kOverwrite)
    rfile.Close()
    
    # Set plotting order and grouping
    def groupfirst(lista, n=1):
        aux = [lista[:n]]
        aux.extend(lista[n:])
        return aux

    fieldname = [Field.name for Field in Fields]
    Min = [Field.Min for Field in Fields]
    Max = [Field.Max for Field in Fields]
    expals = [ROOT.TExec('ex%i' % i, Field.texec) for i, Field in enumerate(Fields)]

    histolist = groupfirst(histolist, args.nj)
    Min = groupfirst(Min, args.nj)
    Max = groupfirst(Max, args.nj)
    expals = groupfirst(expals, args.nj)
    fieldname = groupfirst(fieldname, args.nj)
    
    C = pu.pcanvas(histolist, args)
    hframes = []
    pframes = []
    hists1d = []
    lines = []
    for i, histo in enumerate(histolist):
        padname = 'pad_%i' % (len(histolist) - (i + 1))
        pad = ROOT.gROOT.FindObject(padname)
        pad.cd()
        
        if isinstance(histo, list):
            hframes.append(histo[0].Clone('hframe%i' % i))
        else:
            hframes.append(histo.Clone('hframe%i' % i))
            
        hframes[-1].Reset()
        hframes[-1].Draw('axis')

        if isinstance(histo, list):
            for j, hist in enumerate(histo):
                expals[i][j].Draw()
                if j == 0:
                    hist.Draw('colz same0')
                else:
                    hist.Draw('colz same0')
                pframe = pu.DrawPaletteFrame(hist, len(histo), j)
                pframes.append(pframe)
        else:
            expals[i].Draw()
            if 'a0' in fieldname[i]:
                histo.Draw('colz same0')
            else:
                histo.Draw('colz0 same0')
            pframe = pu.DrawPaletteFrame(histo)
            pframes.append(pframe)

            # Plot 1d outline
            # xmin = histo.GetYaxis().GetXmin()
            # xmax = histo.GetYaxis().GetXmax()
            xmin = ROOT.gPad.GetUymin()
            xmax = ROOT.gPad.GetUymax()
            fmin = histo.GetMinimum()
            fmax = histo.GetMaximum()
            ybin = (xmin + xmax) / 2.0
            bin0 = histo.GetYaxis().FindBin(ybin)
            
            hists1d.append(histo.ProjectionX('%s_px' % histo.GetName(), bin0, bin0))
            hists1d[-1].SetLineColor(ROOT.kOrange + 10)
            hists1d[-1].SetLineWidth(3)
            
            slope = (xmax - xmin) / (fmax - fmin)

            def fx(x, slope=slope, x0=fmin, y0=xmin):
                value = (x - x0) * slope + y0
                if value <= xmin:
                    value = xmin
                elif value >= xmax:
                    value = xmax
                return value

            for j in range(hists1d[-1].GetNbinsX()):
                hists1d[-1].SetBinContent(j + 1, fx(hists1d[-1].GetBinContent(j + 1)))
            lines.append(ROOT.TLine(ROOT.gPad.GetUxmin(), ybin, ROOT.gPad.GetUxmax(), ybin))
            lines[-1].SetLineColor(ROOT.kOrange + 10)
            lines[-1].SetLineStyle(2)
            lines[-1].SetLineWidth(2)
            lines[-1].Draw()
                
            hists1d[-1].Draw('L same')

        if i == 0:  # labels
            xMin = ROOT.gPad.GetUxmin()
            xMax = ROOT.gPad.GetUxmax()
            yMin = ROOT.gPad.GetUymin()
            yMax = ROOT.gPad.GetUymax()
            
            ctext = 'z = %5.2f %s' % ((ct.c * time + z0 - zp) / 1e-3, 'mm')
            texttime = ROOT.TLatex(xMax - (xMax - xMin) / 20., yMax - (yMax - yMin) / 10., ctext)
            texttime.SetTextAlign(32)
            texttime.SetTextFont(43)
            texttime.SetTextSize(28)
            texttime.Draw()

            if args.pden > 0.0:
                ctext = 'n_{0} = %5.2f x %s' % (n0 / 1e6 / 1e18, '10^{18} cm^{-3}')
                textDen = ROOT.TLatex(xMin + (xMax - xMin) / 20., yMax - (yMax - yMin) / 10., ctext)
                textDen.SetTextAlign(12)
                textDen.SetTextFont(43)
                textDen.SetTextColor(ROOT.kGray + 3)
                textDen.SetTextSize(24)
                textDen.Draw()
            
        frame = pu.DrawFrame()
        pframes.append(frame)
    
    # ----------
    ofile = '%s/SnapshotS-%s_%i' % (opath, args.sim, args.time)

    ext = '.pdf'
    if args.png:
        ext = '.png'

    C.cd(0)
    C.Print(ofile + ext)
    
    return 0


if __name__ == '__main__':
    main()
