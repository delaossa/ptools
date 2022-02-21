#!/usr/bin/env python

import numpy as np
import argparse
import ROOT
from ROOT import PData, PGlobals, PPalette, TCanvas, gStyle, gPad, TF1, TF2
import wakefields as wk


def parse_args():

    # Command argument line parser
    parser = argparse.ArgumentParser(description='plotwake.py: plot wakefields')
    parser.add_argument('--TD', dest='TD', action='store_true', default=0, help='Plot 2D wakefields')

    args = parser.parse_args()
        
    return args


def main():

    args = parse_args()

    nb0 = 0.2 / 2  # f0 = a0^2/2 or nb0
    sz = 1.0
    sx = 5.0

    wakeS = wk.WakeLinearSimple(nb0, sz, sx)
    wake = wk.WakeLinear(nb0, sz, sx)
    
    # Style
    PGlobals.SetPlasmaStyle()
    gStyle.SetPadRightMargin(0.05)
    gStyle.SetNumberContours(255)
    gStyle.SetPadGridX(0)
    gStyle.SetPadGridY(0)
    gStyle.SetJoinLinePS(2)
    gStyle.SetFrameLineWidth(2)

    zmin = -4 * np.pi
    zmax = 3 * sz
    xmin = -3 * sx
    xmax = 3 * sx
    
    # Canvas
    C = TCanvas('C', '', 800, 600)

    gPad.SetGridx(1)
    gPad.SetGridy(1)

    fPsi = TF1('fPsi', lambda x: wake.Psi(x[0], 0), zmin, zmax)
    fPsi.SetNpx(100)
    fPsi.SetLineColor(ROOT.kGray + 3)
    fPsi.GetYaxis().SetRangeUser(-2.2 * nb0, 2.2 * nb0)
    fPsi.Draw()
    
    fn1 = TF1('fn1', lambda x: wake.n1(x[0], 0), zmin, zmax)
    fn1.SetNpx(100)
    fn1.SetLineColor(ROOT.kGray + 1)
    fn1.Draw("same")

    fnb = TF1('fnb', lambda x: wake.nb(x[0], 0), zmin, zmax)
    fnb.SetNpx(100)
    fnb.SetLineColor(ROOT.kOrange)
    fnb.Draw("same")

    fEz = TF1('fEz', lambda x: wake.Ez(x[0], 0), zmin, zmax)
    fEz.SetNpx(100)
    fEz.SetLineColor(ROOT.kRed)
    fEz.Draw("same")
    
    fWx = TF1('fWx', lambda x: wake.Wx(x[0], sx), zmin, zmax)
    fWx.SetNpx(100)
    fWx.SetLineColor(ROOT.kAzure + 2)
    fWx.Draw("same")

    fKx = TF1('fKx', lambda x: wake.Kx(x[0], 0.0), zmin, zmax)
    fKx.SetNpx(100)
    fKx.SetLineColor(ROOT.kGreen + 1)
    fKx.Draw("same")


    fsPsi = TF1('fsPsi', lambda x: wakeS.Psi(x[0], 0), zmin, zmax)
    fsPsi.SetNpx(100)
    fsPsi.SetLineColor(ROOT.kGray + 3)
    fsPsi.SetLineStyle(2)
    fsPsi.SetLineWidth(1)
    fsPsi.Draw("same")
    
    fsn1 = TF1('fsn1', lambda x: wakeS.n1(x[0], 0), zmin, zmax)
    fsn1.SetNpx(100)
    fsn1.SetLineColor(ROOT.kGray + 1)
    fsn1.SetLineStyle(2)
    fsn1.SetLineWidth(1)
    fsn1.Draw("same")

    fsEz = TF1('fsEz', lambda x: wakeS.Ez(x[0], 0), zmin, zmax)
    fsEz.SetNpx(100)
    fsEz.SetLineColor(ROOT.kRed)
    fsEz.SetLineStyle(2)
    fsEz.SetLineWidth(1)
    fsEz.Draw("same")
    
    fsWx = TF1('fsWx', lambda x: wakeS.Wx(x[0], sx), zmin, zmax)
    fsWx.SetNpx(100)
    fsWx.SetLineColor(ROOT.kAzure + 2)
    fsWx.SetLineStyle(2)
    fsWx.SetLineWidth(1)
    fsWx.Draw("same")

    fsKx = TF1('fsKx', lambda x: wakeS.Kx(x[0], 0.0), zmin, zmax)
    fsKx.SetNpx(100)
    fsKx.SetLineColor(ROOT.kGreen + 1)
    fsKx.SetLineStyle(2)
    fsKx.SetLineWidth(1)
    fsKx.Draw("same")

    gPad.Print('Fields.pdf')

    if args.TD:

        fKx = TF2('fKx', lambda x: wakeS.Kx(x[0], x[1]), zmin, zmax, xmin, xmax)
        fKx.SetNpx(100)
        fKx.SetNpy(50)
        fKx.Draw("col")
        gPad.Print('Kx2D.pdf')

        fEz = TF2('fEz', lambda x: wakeS.Ez(x[0], x[1]), zmin, zmax, xmin, xmax)
        fEz.SetNpx(100)
        fEz.SetNpy(50)
        fEz.Draw("col")
        gPad.Print('Ez2D.pdf')

        '''
        fEx = TF2('fEx', lambda x: wakeS.Ex(x[0], x[1]), zmin, zmax, xmin, xmax)
        fEx.SetNpx(100)
        fEx.SetNpy(50)
        fEx.Draw("col")
        gPad.Print('Ex2D.pdf')
        '''
        
        fWx = TF2('fWx', lambda x: wakeS.Wx(x[0], x[1]), zmin, zmax, xmin, xmax)
        fWx.SetNpx(100)
        fWx.SetNpy(50)
        fWx.Draw("col")
        gPad.Print('Wx2D.pdf')
    
    
if __name__ == '__main__':
    main()
