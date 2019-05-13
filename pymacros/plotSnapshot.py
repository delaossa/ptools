#!/usr/bin/env python

import os, argparse
import numpy as np
import scipy.constants as ct
import ROOT
from ROOT import gStyle, TCanvas, gPad, TFile, TExec, TPave, TLine, TGaxis, TLatex, TMarker, PData, PDataHiP, PGlobals

def parse_args():

    # Command argument line parser 
    parser = argparse.ArgumentParser(description='')

    parser.add_argument('sim', nargs='?', default=None, help='simulation name')

    parser.add_argument('-t', type=int, dest='time', default=0, help='time step dump')
    parser.add_argument('-i', type=int, dest='itime', default=0, help='initial time step dump')
    parser.add_argument('-f', type=int, dest='ftime', default=0, help='final time step dump')
    parser.add_argument('-s', type=int, dest='stime', default=1, help='advance step dump')
    parser.add_argument('-np', type=float, dest='pden', default=0.0, help='plasma density')
    parser.add_argument('--logz', action='store_true', default=0, help='logaritmic z scale')
    parser.add_argument('--png', action='store_true', default=0, help='png ')
    parser.add_argument('--joint', action='store_true', default=0, help='joints the histograms')
    parser.add_argument('--hp', action='store_true', default=0, help='HiPACE simulation')
    parser.add_argument('opath', nargs='?', default=None, help='Output folder')
    
    args = parser.parse_args()
        
    return args

def DrawFrame() :
    gPad.Update()
    gPad.RedrawAxis('g')
    gPad.RedrawAxis()

    hframe = TPave(gPad.GetUxmin(),gPad.GetUymin(),gPad.GetUxmax(),gPad.GetUymax(),2,'')
    hframe.SetFillStyle(0)
    hframe.SetLineWidth(3)
    hframe.SetLineColor(ROOT.kBlack)
    hframe.SetShadowColor(0)
    hframe.Draw()
    
    return hframe


def DrawPaletteFrame(histo,Ns=1,ip=0) :
    gPad.Update()

    
    y1 = gPad.GetBottomMargin()
    y2 = 1 - gPad.GetTopMargin()
    x1 = gPad.GetLeftMargin()
    x2 = 1 - gPad.GetRightMargin()
    pvsize = (y2-y1)/Ns
    gap = 0.02

    palette = histo.GetListOfFunctions().FindObject('palette')
    if palette :
        if ip==0 :
	    y1b = y1 + ip*pvsize
	else :
	    y1b = y1 + ip*pvsize + gap/2.0
            
	if ip==Ns-1 :  
	    y2b = y1 + (ip+1)*pvsize
	else :
	    y2b = y1 + (ip+1)*pvsize - gap/2.0
        
        palette.SetY1NDC(y1b)
        palette.SetY2NDC(y2b)
        palette.SetX1NDC(x2 + 0.01)
        palette.SetX2NDC(x2 + 0.03)
        palette.SetBorderSize(2)
        palette.SetLineColor(1)
        
        pFrame = TPave((x2 + 0.01),y1b,(x2 + 0.03),y2b,1,'NDCL')
        pFrame.SetFillStyle(0)
        pFrame.SetLineWidth(2)
        pFrame.SetLineColor(ROOT.kBlack)
        pFrame.SetShadowColor(0)
        
        pFrame.Draw()

    return pFrame

def pcanvas(histos,args):
    ypadsize = 500
    ymarginsize = 200
    ysize = ypadsize + ymarginsize
    xsize = 1050
    ar0 = float(ysize)/float(xsize)
    C = TCanvas('C','Snapshot',xsize,ysize)
    C.SetFillStyle(4000)

    bMargin = 0.12 * (950.0/ysize)
    tMargin = 0.04 * (950.0/ysize)
    lMargin = 0.14
    rMargin = 0.15
    mMargin = 0.02 * (950.0/ysize)

    fonttype = 43
    fontsize = 38
    tfontsize = 42
    txsize = tfontsize+6
    lxsize = fontsize+2
    tysize = tfontsize
    lysize = fontsize-2
    tzsize = tfontsize-4
    lzsize = fontsize-6
    txoffset = 2.4 / (950.0/ysize)
    lxoffset = 0.015
    tyoffset = 1.2 / (950.0/ysize)
    lyoffset = 0.01
    tzoffset = 1.1 / (950.0/ysize)
    lzoffset = 0.01
    txlength = 0.02
    tylength = 0.01

    N = len(histos)
    PGlobals.CanvasPartition(C,N,lMargin,rMargin,bMargin,tMargin,mMargin)
    basepad = C.cd(0)

    print(' aspect ratio = ',ar0)
    pad = np.empty(N, dtype=object)
    for i in range(N) :
        ipad = N-1-i
        padname = 'pad_%i'%ipad
        print('pad name = %s'%padname)
        #pad.append(ROOT.gROOT.FindObject(padname))
        pad[i] = ROOT.gROOT.FindObject(padname)
        pad[i].cd()
        pad[i].SetTickx(1)
        pad[i].SetTicky(1)
        if args.logz==True :
            pad[i].SetLogz(1)

        xFactor = pad[0].GetAbsWNDC()/pad[i].GetAbsWNDC()
        yFactor = pad[0].GetAbsHNDC()/pad[i].GetAbsHNDC()
        ar = ar0*xFactor/yFactor
        
        def setaxes(hist,NS=1) :
            hist.GetYaxis().SetLabelFont(fonttype)
            hist.GetYaxis().SetLabelSize(lysize)
            hist.GetYaxis().SetLabelOffset(lyoffset)
            hist.GetYaxis().SetTitleFont(fonttype)
            hist.GetYaxis().SetTitleSize(tysize)
            hist.GetYaxis().SetTitleOffset(tyoffset)
            
            hist.GetYaxis().SetTickLength(ar*tylength)
    
            hist.GetXaxis().SetLabelFont(fonttype)
            hist.GetXaxis().SetLabelSize(lxsize)
            hist.GetXaxis().SetLabelOffset(lxoffset)
            hist.GetXaxis().SetTitleFont(fonttype)
            hist.GetXaxis().SetTitleSize(txsize)
            hist.GetXaxis().SetTitleOffset(txoffset)
            hist.GetXaxis().SetTickLength(txlength/ar)
            
            hist.GetZaxis().SetTitleFont(fonttype)
            hist.GetZaxis().SetTitleOffset(tzoffset)
            hist.GetZaxis().SetTitleSize(tzsize)
            hist.GetZaxis().SetLabelFont(fonttype)
            hist.GetZaxis().SetLabelSize(lzsize)
            hist.GetZaxis().SetTickLength(NS*ar*tylength)     
            if args.logz :
	        hist.GetZaxis().SetLabelOffset(-lyoffset * (250.0/ypadsize))
            else :
                hist.GetZaxis().SetLabelOffset(lyoffset)
            
        if isinstance(histos[i], list) :
            for j,h in enumerate(histos[i]) :
                setaxes(h,len(histos[i]))
                if i!=N-1 :
                    h.GetXaxis().SetLabelSize(0)
                    h.GetXaxis().SetTitleSize(0)
        else :
            setaxes(histos[i])
            if i!=N-1 :
                histos[i].GetXaxis().SetLabelSize(0)
                histos[i].GetXaxis().SetTitleSize(0)

    return C

def diffyhist(hist) :
    NBinsX = hist.GetXaxis().GetNbins()
    NBinsY = hist.GetYaxis().GetNbins()

    histdy = hist.Clone('hdiffy')
    histdy.Reset()
    
    dy = hist.GetYaxis().GetBinWidth(1);
    for i in range(NBinsX):
        for j in range(NBinsY):
	    der = 0.
	    if (j>2 and j<NBinsY-2) : 
	        der =  (4.0 / 3.0) * (hist.GetBinContent(i,j+1) - hist.GetBinContent(i,j-1)) / (2.0 * dy)
	        - (1.0 / 3.0) * (hist.GetBinContent(i,j+2) - hist.GetBinContent(i,j-2)) / (4.0 * dy) 	  
	
	    histdy.SetBinContent(i,j,der)

    return histdy
    

def main():

    ROOT.gROOT.SetBatch(1)
    PGlobals.Initialize()
    gStyle.SetNumberContours(99)
    gStyle.SetLineWidth(2)
    gStyle.SetNdivisions(005,'xyz')
    gStyle.SetNdivisions(010,'x')
    ROOT.gROOT.ForceStyle()
    
    args = parse_args()
    sim = args.sim
    time = args.time
    
    if args.hp :
        pData = PDataHiP.Get(sim)
    else :
        pData = PData.Get(sim)
        
    pData.LoadFileNames(time)
    
    filename = './%s/Plots/Snapshots/Snapshot-%s_%i.root' % (sim,sim,time)
    # print(filename)
    ifile = TFile(filename,'READ')

    NS = pData.NSpecies()
    Min = 0.101 * np.ones(NS)
    Max = 100.0 * np.ones(NS)
        
    hDen2D = []
    for i in range(NS) :
        hName = 'hDen2D_%i'%i
        hDen2D.append(ifile.Get(hName))

        if pData.GetDenMin(i) != -999. : 
	    Min[i] = pData.GetDenMin(i)
        
        if pData.GetDenMax(i) != -999. : 
	    Max[i] = pData.GetDenMax(i)
        
        if args.pden>0.0 :
            hDen2D[-1].Scale(args.pden)
            # hDen2D[-1].GetZaxis().SetTitle('10^{18} cm^{-3}')
            hDen2D[-1].GetZaxis().SetTitle('')

            Min[i] = Min[i] * args.pden
            Max[i] = Max[i] * args.pden    

    Npal = NS
    if args.joint :
        for i in range(1,NS) :
            hDen2D[0].Add(hDen2D[i])
        Npal = 1
        
    hName = 'hE1D_0'
    hEz1D = ifile.Get(hName)

    # Get Focusing field
    hName = 'hW2D_0'
    hW2D = ifile.Get(hName)

    # 
    hK2D = diffyhist(hW2D)
    hK2D.Scale(1E9/ct.c)
    
    histos = []
    histos.append([hDen2D[0],hDen2D[1]])
    histos.append(hW2D)
    histos.append(hK2D)
    
    C = pcanvas(histos,args)

    pad = ROOT.gROOT.FindObject('pad_2')
    pad.cd()

    # Setup plasma palette
    exPlasma = TExec('exPlasma','plasmaPalette->cd();')
    plasmaPalette = ROOT.gROOT.FindObject('plasma')

    baseden = 1.0
    if args.pden>0.0 :
        baseden *= args.pden

    if args.logz :
        a = 1.0/(np.log(Max[0])-np.log(Min[0]))
        b = np.log(Min[0])
        basePos = a*(np.log(baseden) - b)
    else :
        basePos = (1.0/(Max[0]-Min[0]))*(baseden - Min[0])
    
    Stops = np.array( [0.00, basePos, 1.00] )
    Red   = np.array( [0.99, 0.90, 0.00] )
    Green = np.array( [0.99, 0.90, 0.00] )
    Blue  = np.array( [0.99, 0.90, 0.00] )
    NRGBs = 3
    NCont = 99
    plasmaPalette.CreateGradientColorTable(NRGBs, Stops, Red, Green, Blue, NCont, 1.0)
    
    # Get plasma base color
    baseden = 0.5
    if args.pden>0.0 :
        baseden *= args.pden
    ipos = int((baseden-Min[0])*(plasmaPalette.GetNColors()-1)/(Max[0]-Min[0]))
    if args.logz :
        ipos = int((np.log(baseden)-np.log(Min[0]))*(plasmaPalette.GetNColors()-1)/(np.log(Max[0])-np.log(Min[0])))
    
    if ipos<0 :
        ipos = 0
    cindex = plasmaPalette.GetColor(ipos)
    basecolor=ROOT.gROOT.GetColor(cindex)

    # Setup beam palette
    exBeam   = TExec('exBeam','beamPalette->cd();')
    beamPalette = ROOT.gROOT.FindObject('beam')

    Stops = np.array( [0.00, 0.40, 0.50, 0.60, 1.00] )
    Red   = np.array( [basecolor.GetRed(), 0.22, 0.39, 0.70, 1.00] )
    Green = np.array( [basecolor.GetGreen(), 0.34, 0.05, 0.20, 1.00] )
    Blue  = np.array( [basecolor.GetBlue(), 0.58, 0.33, 0.30, 0.20] )
    NRGBs = len(Stops)
    NCont = 99
    beamPalette.CreateGradientColorTable(NRGBs, Stops, Red, Green, Blue, NCont, 1.0)
    # ---

    hFrame = hDen2D[0].Clone('hFrame')
    hFrame.Reset()
    hFrame.Draw('axis')

    exPlasma.Draw()
    hDen2D[0].GetZaxis().SetRangeUser(Min[0],Max[0])
    hDen2D[0].Draw('colz same')
    
    p1Frame = DrawPaletteFrame(hDen2D[0],Npal,0)

    if not args.joint :
        exBeam.Draw()
        hDen2D[1].GetZaxis().SetRangeUser(Min[1],Max[1])
        hDen2D[1].Draw('colz same')
        p2Frame = DrawPaletteFrame(hDen2D[1],Npal,1)
    
    yRange   = (hDen2D[0].GetYaxis().GetXmax() - hDen2D[0].GetYaxis().GetXmin())
    midPoint = (hDen2D[0].GetYaxis().GetXmax() + hDen2D[0].GetYaxis().GetXmin())/2.
    yMin = midPoint-yRange/2
    yMax = midPoint+yRange/2
    xMin = hDen2D[0].GetXaxis().GetXmin()
    xMax = hDen2D[0].GetXaxis().GetXmax()

    Ezmin = -100.0
    Ezmax =  100.0
    y1 = yMin + yRange/20.
    y2 = yMin + yRange/3.0 - yRange/20.
    slope = (y2-y1)/(Ezmax-Ezmin)

    def fez(ez,slope=slope,e1=Ezmin,y1=y1) :
        return (ez-e1)*slope + y1
    
    # color = ROOT.kPink
    color = ROOT.kPink-8
    lineEzero = TLine(xMin,fez(0),xMax,fez(0))
    lineEzero.SetLineColor(color)
    lineEzero.SetLineStyle(2)
    lineEzero.SetLineWidth(1)
    lineEzero.Draw()

    hclone = hEz1D.Clone('hEz1Dclone')
    for i in range(hclone.GetNbinsX()) :
	hclone.SetBinContent(i+1,fez(hclone.GetBinContent(i+1)))
        
    hclone.SetLineStyle(1)
    hclone.SetLineWidth(3)
    hclone.SetLineColor(color)

    hclone.Draw('C same')

    xmin = hclone.GetXaxis().GetXmin()
    xmax = hclone.GetXaxis().GetXmax()
    # xpos = xmax - (xmax-xmin) * 0.10
    # xpos = xmin + (xmax-xmin) * 0.20
    xpos = xmax
    axis = TGaxis(xpos,y1,xpos,y2,Ezmin,Ezmax,203,'S+=R')
    axis.SetLineWidth(3)
    axis.SetLineColor(color)
    axis.SetLabelColor(color)
    axis.SetLabelFont(43)
    axis.SetLabelSize(28)
    axis.SetLabelOffset(0.005)
    #
    axis.SetLabelSize(0.)
    #
    axis.SetTitleColor(color)
    axis.SetTitleFont(43)
    axis.SetTitleSize(22)
    axis.SetTitleOffset(0.8)
    axis.SetTickSize(0.03)
    # axis.ChangeLabel(1,-1,0.)
    axis.ChangeLabel(2,-1,0)
    
    #axis.SetTitle('GV/m')
    axis.SetTitle('')
    #axis.CenterTitle()
    #axis.ChangeLabel(1,-1,-1,-1,-1,-1,'')
    #axis.SetMaxDigits(2)

    padFrame = DrawFrame()

    axis.Draw()

    zpos = -35
    zbin = hEz1D.FindBin(zpos)
    ezvalue = hEz1D.GetBinContent(zbin)
    ezpos = hclone.GetBinContent(zbin)
    point = TMarker(zpos,ezpos,20)
    point.SetMarkerColor(color)
    point.SetMarkerSize(1)
    point.Draw()

    line = TLine(zpos,fez(0),zpos,ezpos)
    line.SetLineColor(color)
    line.SetLineWidth(1)
    line.SetLineStyle(3)
    line.Draw()
    
    ctext = '%.2f GV/m'%ezvalue
    textField = TLatex()
    textField.SetTextAlign(12);
    textField.SetTextFont(43);
    textField.SetTextSize(28)
    textField.SetTextColor(color)
    textField.DrawLatex(zpos+1,fez(-50),ctext)


    # -----------
    pad = ROOT.gROOT.FindObject('pad_1')
    pad.cd()

    exField = TExec('exField','fieldPalette->cd();')
    fieldPalette = ROOT.gROOT.FindObject('field')

    hFrame2 = hW2D.Clone('hFrame')
    hFrame2.Reset()
    hFrame2.Draw('axis')

    hW2D.GetZaxis().SetRangeUser(-50,50)
    hW2D.GetZaxis().SetTitle('W_{x}[GV/m]')
    exField.Draw()
    hW2D.Draw('colz0 same')
    p3Frame = DrawPaletteFrame(hW2D)

    padFrame2 = DrawFrame()
    
    # -----------
    pad = ROOT.gROOT.FindObject('pad_0')
    pad.cd()

    hFrame3 = hK2D.Clone('hFrame3')
    hFrame3.Reset()
    hFrame3.Draw('axis')

    # hK2D.GetZaxis().SetRangeUser(-60,60)
    hK2D.GetZaxis().SetTitle('g[MT/m]')
    exField.Draw()
    hK2D.Draw('colz0 same')
    p4Frame = DrawPaletteFrame(hK2D)

    padFrame3 = DrawFrame()
    
    if args.opath :
        opath = args.opath    
        if not os.path.exists(opath):
            os.makedirs(opath)
        ofile = '%s/SnapshotS-%s_%i'%(opath,sim,time)
    else :
        ofile = '%s/Plots/Snapshots/SnapshotS-%s_%i'%(sim,sim,time)

    ext = 'pdf'
    if args.png :
        ext = 'png'
        
    # ofile = ofile + ext
    PGlobals.imgconv(C,ofile,ext)
    # C.Print(ofile)
    # C.Clear()


if __name__ == '__main__':
    main()
