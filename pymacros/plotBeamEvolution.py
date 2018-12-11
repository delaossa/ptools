#!/usr/bin/env python

import os, argparse
import numpy as np
import ROOT
from ROOT import TCanvas, gStyle, gPad, TFile, TGraph, TH1F, TLine, TGaxis, PData, PGlobals, PPalette


def parse_args():

    # Command argument line parser 
    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-index', type=int, dest='index', default=0, help='particle index')
    parser.add_argument('opath', nargs='?', default='Plots/Evolution', help='Output folder')
    
    args = parser.parse_args()
        
    return args

    
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

def main():

    home = os.getcwd()
    
    args = parse_args()

    if args.index == 0 :
        species = 'plasma'

        xMin = -2.1
        xMax =  2.1

        exMin = 0.01
        exMax = 50.0

        gMin = 0.01
        gMax = 1600.0

        grrMin = 0.01
        grrMax = 50.0

        cMin = 0.0
        cMax = 500000.0

    elif args.index == 1 :
        species = 'He-2'
    
        xMin = -2.0
        xMax = 2.0

        xrmsMin = 0.01
        xrmsMax = 1.00

        exMin = 0.001
        exMax = 5.00
        # exMax = 0.8

        gMin = 0.001
        gMax = 1.599
      
        grrMin = 0.001
        grrMax = 5.0

        cMin = 0.001
        cMax = 299.99

        betaMin = 0.001
        betaMax = 1.0
        
        alphaMin = -0.5
        alphaMax = 0.5
        
        gammatMin = 0.001
        gammatMax = 20.0
        
        
    simlist = []
    colorlist = []
    gXlist = []        
    gXabslist = []        
    gXabsuplist = []        
    gXrmslist = []
    gxPxrmslist = []
    gXuplist = []
    gXdolist = []
    gEmitxlist = []
    gEmitxavglist = []        
    gPxlist = []        
    gPxrmslist = []        
    gBetaxlist = []
    gAlphaxlist = []
    gGammaxlist = []
    gYlist = []        
    gYabslist = []        
    gYabsuplist = []        
    gYrmslist = []
    gyPyrmslist = []
    gYuplist = []
    gYdolist = []
    gEmitylist = []        
    gEmityavglist = []        
    gPylist = []        
    gPyrmslist = []        
    gBetaylist = []
    gAlphaylist = []
    gGammaylist = []
    gGammalist = []        
    gGrmslist = []        
    gGrmsrellist = []        
    gGrmsavglist = []        
    gChargelist = []
    
    opath = args.opath
    opath = opath + ('/%s' % species)
    if not os.path.exists(opath):
        os.makedirs(opath)

    # Style
    PGlobals.SetPlasmaStyle()
    # gStyle.SetPadRightMargin(0.05)
    gStyle.SetPadRightMargin(0.10)
    gStyle.SetPadBottomMargin(0.30)
    gStyle.SetNdivisions(006,'xyz')
    gStyle.SetNdivisions(206,'x')
    gStyle.SetNdivisions(503,'y')
    gStyle.SetTickLength(0.022,'x');
    gStyle.SetTitleOffset(0.5,'y')
    gStyle.SetLabelSize(36,'xyz')
    gStyle.SetTitleSize(38,'xyz')
    gStyle.SetLineWidth(2)
    gStyle.SetNumberContours(255)
    gStyle.SetPadGridX(0)
    gStyle.SetPadGridY(0)
    gStyle.SetPadTickX(1)
    gStyle.SetPadTickY(1)
    gStyle.SetJoinLinePS(1)
    gStyle.SetFrameLineWidth(2)
    gStyle.SetFrameFillStyle(4000)
    gStyle.SetMarkerSize(1)
    gStyle.SetMarkerStyle(20)
    
    ROOT.gROOT.SetBatch(1)

    # Canvas
    C = TCanvas('C',' -- O -- ',1024,380)
    C.SetFillStyle(4000)        

    # Main plot frame
    hFrame = TH1F('hFrame','',10,0,6)

    hFrame.GetXaxis().SetTitle('z [mm]')
    hFrame.GetXaxis().SetNdivisions(510)
    hFrame.GetXaxis().CenterTitle()
    hFrame.GetYaxis().CenterTitle()

    pal = PPalette('pal')
    pal.SetPalette(ROOT.kBird)

    # Select simulations 
    sim = 'rake-v10kA.n8E18.He6.2nd.mrplow.2.3D'
    simlist.append(sim)
            
    Ncurves = len(simlist)
    for sim in simlist :
        # colorlist.append(pal.GetColorIndex((i+1) * (pal.GetNColors()-1) / (Ncurves+1)))
        # colorlist.append(ROOT.kMagenta)
        # colorlist.append(ROOT.kPink-9)
        colorlist.append(ROOT.kOrange+10)

        filename = '%s/Plots/Bunch/%s/Bunch-Evolution-%s.root' % (sim,species,sim)
        rfile = TFile(filename,"READ")

        # Load histos and graphs
        gXlist.append(rfile.Get("gXmeanvsTime")) 
        gXrmslist.append(rfile.Get("gXrmsvsTime")) 
        gPxlist.append(rfile.Get("gPxmeanvsTime")) 
        gPxrmslist.append(rfile.Get("gPxrmsvsTime")) 
        gEmitxlist.append(rfile.Get("gEmitxvsTime")) 
        gEmitxavglist.append(rfile.Get("gEmitxavgvsTime"))
        gxPxrmslist.append(rfile.Get("gxPxrmsvsTime"))
        gYlist.append(rfile.Get("gYmeanvsTime")) 
        gYrmslist.append(rfile.Get("gYrmsvsTime")) 
        gPylist.append(rfile.Get("gPymeanvsTime")) 
        gPyrmslist.append(rfile.Get("gPyrmsvsTime")) 
        gEmitylist.append(rfile.Get("gEmityvsTime")) 
        gEmityavglist.append(rfile.Get("gEmityavgvsTime")) 
        gyPyrmslist.append(rfile.Get("gyPyrmsvsTime"))
        gGammalist.append(rfile.Get("gPzmeanvsTime"))
        gGrmslist.append(rfile.Get("gPzrmsvsTime"))
        gGrmsavglist.append(rfile.Get("gErmsavgvsTime"))
        gChargelist.append(rfile.Get("gChargevsTime"))
        
        rfile.Close()

        # Get arrays from graphs
        NP = gXlist[-1].GetN()
        
        T  = gXlist[-1].GetX()
        X  = gXlist[-1].GetY()
        Xrms  = gXrmslist[-1].GetY()
        xPxrms = gxPxrmslist[-1].GetY()
        Px  = gPxlist[-1].GetY()
        Pxrms  = gPxrmslist[-1].GetY()
        Emitx = gEmitxlist[-1].GetY()
        Y  = gYlist[-1].GetY()
        Yrms  = gYrmslist[-1].GetY()
        Py  = gPylist[-1].GetY()
        Pyrms  = gPyrmslist[-1].GetY()
        yPyrms = gyPyrmslist[-1].GetY()
        Emity = gEmitylist[-1].GetY()
        Gamma = gGammalist[-1].GetY()
        Grms  = gGrmslist[-1].GetY()
        Charge = gChargelist[-1].GetY() 
        
        Grmsrel = np.empty(NP)
        Betax = np.empty(NP)
        Betay = np.empty(NP)
        Alphax = np.empty(NP)
        Alphay = np.empty(NP)
        Gammax = np.empty(NP)
        Gammay = np.empty(NP)
        for j in range(NP) :
            Grmsrel[j] = 100 * Grms[j]/Gamma[j]
            Betax[j]  = 0.001 * (Gamma[j]/0.511E-3) * Xrms[j]**2 / (Emitx[j]*0.1)
            Betay[j]  = 0.001 * (Gamma[j]/0.511E-3) * Yrms[j]**2 / (Emity[j]*0.1)
            Gammax[j] = 1000 * (0.511E-3/Gamma[j]) * (Pxrms[j]/0.511)**2 / (Emitx[j]*0.1)
            Gammay[j] = 1000 * (0.511E-3/Gamma[j]) * (Pyrms[j]/0.511)**2 / (Emity[j]*0.1)
            Alphax[j] = -(xPxrms[j]/0.511) / (Emitx[j]*0.1)
            Alphay[j] = -(yPyrms[j]/0.511) / (Emity[j]*0.1)
            
        gGrmsrellist.append(TGraph(NP,T,Grmsrel))
        gBetaxlist.append(TGraph(NP,T,Betax))
        gBetaylist.append(TGraph(NP,T,Betay))
        gAlphaxlist.append(TGraph(NP,T,Alphax))
        gAlphaylist.append(TGraph(NP,T,Alphay))
        gGammaxlist.append(TGraph(NP,T,Gammax))
        gGammaylist.append(TGraph(NP,T,Gammay))
            
        Xabs = np.empty((NP))
        Xabsup = np.empty((NP))
        Xup = np.empty((NP))
        Xdo = np.empty((NP))
        Yabs = np.empty((NP))
        Yabsup = np.empty((NP))
        Yup = np.empty((NP))
        Ydo = np.empty((NP))
        for j in range(NP) :
            Xup[j] = X[j] + Xrms[j]
            Xdo[j] = X[j] - Xrms[j]
            
            Xabs[j] = np.abs(X[j]) #/X[0]
            Xabsup[j] = Xabs[j] + Xrms[j]
                
            Yup[j] = Y[j] + Yrms[j]
            Ydo[j] = Y[j] - Yrms[j]
            
            Yabs[j] = np.abs(Y[j]) #/Y[0]
            Yabsup[j] = Yabs[j] + Yrms[j]
            '''
            # Auto ranging
            if Xup[j] > xMax : xMax = Xup[j]
            if Xdo[j] < xMin : xMin = Xdo[j]

            if Emitx[j] > exMax : exMax = Emitx[j]
            if Emitx[j] < exMin : exMin = Emitx[j]

            if Gamma[j] > gMax : gMax = Gamma[j]
            if Gamma[j] < gMin : gMin = Gamma[j]

            if Grmsrel[j] > grrMax : grrMax = Grmsrel[j]
            if Grmsrel[j] < grrMin : grrMin = Grmsrel[j]

            if Charge[j] > cMax : cMax = Charge[j]
            if Charge[j] < cMin : cMin = Charge[j]
            '''
                
        gXuplist.append(TGraph(NP,T,Xup))
        gXdolist.append(TGraph(NP,T,Xdo))

        gXuplist[-1].SetLineWidth(1)
        gXuplist[-1].SetLineColor(colorlist[-1])
        gXuplist[-1].SetLineStyle(3)
        
        gXdolist[-1].SetLineWidth(1)
        gXdolist[-1].SetLineColor(colorlist[-1])
        gXdolist[-1].SetLineStyle(3)
        
        gXabslist.append(TGraph(NP,T,Xabs))
        gXabsuplist.append(TGraph(NP,T,Xabsup))
        
        gXabslist[-1].SetLineWidth(3)
        gXabslist[-1].SetLineColor(colorlist[-1])

        gXabsuplist[-1].SetLineWidth(1)
        gXabsuplist[-1].SetLineStyle(1)
        gXabsuplist[-1].SetLineColor(colorlist[-1])

        gXlist[-1].SetLineWidth(3)
        gXlist[-1].SetLineColor(colorlist[-1])

        gXrmslist[-1].SetLineWidth(3)
        gXrmslist[-1].SetLineColor(colorlist[-1])
            
        gEmitxlist[-1].SetLineWidth(3)
        gEmitxlist[-1].SetLineColor(colorlist[-1])

        gEmitxavglist[-1].SetLineWidth(1)
        gEmitxavglist[-1].SetLineColor(colorlist[-1])

        gBetaxlist[-1].SetLineWidth(3)
        gBetaxlist[-1].SetLineColor(colorlist[-1])

        gAlphaxlist[-1].SetLineWidth(3)
        gAlphaxlist[-1].SetLineColor(colorlist[-1])

        gGammaxlist[-1].SetLineWidth(3)
        gGammaxlist[-1].SetLineColor(colorlist[-1])

        gYuplist.append(TGraph(NP,T,Yup))
        gYdolist.append(TGraph(NP,T,Ydo))

        gYuplist[-1].SetLineWidth(1)
        gYuplist[-1].SetLineColor(colorlist[-1])
        gYuplist[-1].SetLineStyle(3)
        
        gYdolist[-1].SetLineWidth(1)
        gYdolist[-1].SetLineColor(colorlist[-1])
        gYdolist[-1].SetLineStyle(3)
        
        gYabslist.append(TGraph(NP,T,Yabs))
        gYabsuplist.append(TGraph(NP,T,Yabsup))
        
        gYabslist[-1].SetLineWidth(3)
        gYabslist[-1].SetLineColor(colorlist[-1])

        gYabsuplist[-1].SetLineWidth(1)
        gYabsuplist[-1].SetLineStyle(1)
        gYabsuplist[-1].SetLineColor(colorlist[-1])

        gYlist[-1].SetLineWidth(3)
        gYlist[-1].SetLineColor(colorlist[-1])

        gYrmslist[-1].SetLineWidth(3)
        gYrmslist[-1].SetLineColor(colorlist[-1])

        gEmitylist[-1].SetLineWidth(3)
        gEmitylist[-1].SetLineColor(colorlist[-1])

        gEmityavglist[-1].SetLineWidth(1)
        gEmityavglist[-1].SetLineColor(colorlist[-1])

        gBetaylist[-1].SetLineWidth(3)
        gBetaylist[-1].SetLineColor(colorlist[-1])

        gAlphaylist[-1].SetLineWidth(3)
        gAlphaylist[-1].SetLineColor(colorlist[-1])

        gGammaylist[-1].SetLineWidth(3)
        gGammaylist[-1].SetLineColor(colorlist[-1])

        gGammalist[-1].SetLineWidth(3)
        gGammalist[-1].SetLineColor(colorlist[-1])

        gGrmslist[-1].SetLineWidth(2)
        gGrmslist[-1].SetLineColor(colorlist[-1])

        gGrmsrellist[-1].SetLineWidth(2)
        gGrmsrellist[-1].SetLineColor(colorlist[-1])

        gGrmsavglist[-1].SetLineWidth(2)
        gGrmsavglist[-1].SetLineColor(colorlist[-1])

        gChargelist[-1].SetLineWidth(2)
        gChargelist[-1].SetLineColor(colorlist[-1])


    hFrame.GetYaxis().SetTitle('#LTX_{b}#GT [#mum]')
    hFrame.GetYaxis().SetRangeUser(xMin, xMax)
    hFrame.Draw("axis")

    for idx,graph in enumerate(gXlist) :
        graph.SetLineColor(ROOT.kGray+3)
        graph.SetLineWidth(3)
        graph.SetMarkerColor(ROOT.kGray+3)
        graph.SetMarkerSize(1)
        graph.Draw('LP')

        gXuplist[idx].SetLineColor(ROOT.kGray+3)
        gXuplist[idx].Draw('L')
        gXdolist[idx].SetLineColor(ROOT.kGray+3)
        gXdolist[idx].Draw('L')

        gYlist[-1].SetLineColor(ROOT.kGray+2)
        gYlist[-1].SetLineWidth(3)
        gYlist[-1].SetLineStyle(1)
        gYlist[-1].SetMarkerColor(ROOT.kGray+2)
        gYlist[-1].SetMarkerSize(1)
        gYlist[-1].Draw('LP')

        gYuplist[idx].SetLineColor(ROOT.kGray+2)
        gYuplist[idx].Draw('L')
        gYdolist[idx].SetLineColor(ROOT.kGray+2)
        gYdolist[idx].Draw('L')

        

    DrawFrame()
        
    C.Print(opath + '/%s-Xmean.pdf' % (species))

    hFrame.GetYaxis().SetTitle('#LT#sigma_{x}#GT [#mum]')
    hFrame.GetYaxis().SetRangeUser(xrmsMin,xrmsMax)

    hFrame.Draw("axis")
    for idx,graph in enumerate(gXrmslist) :
        # if list1[-1] == el1 :
        graph.SetLineColor(ROOT.kGray+3)
        graph.SetLineWidth(3)
        graph.SetMarkerColor(ROOT.kGray+3)
        graph.SetMarkerSize(1)
        graph.Draw('LP')

        gYrmslist[-1].SetLineColor(ROOT.kGray+2)
        gYrmslist[-1].SetLineWidth(3)
        gYrmslist[-1].SetLineStyle(1)
        gYrmslist[-1].SetMarkerColor(ROOT.kGray+2)
        gYrmslist[-1].SetMarkerSize(1)
        gYrmslist[-1].Draw('LP')

        
    DrawFrame()
    C.Print(opath + '/%s-Xrms.pdf' % (species))

    # Emittance        
    # hFrame.GetYaxis().SetTitle('k_{p} #varepsilon_{n,x}')
    hFrame.GetYaxis().SetTitle('#varepsilon_{n} [100 nm]')
    hFrame.GetYaxis().SetRangeUser(exMin,exMax)

    hFrame.Draw("axis")        
    for idx,graph in enumerate(gEmitxlist) :
        # graph.Draw('C')
        
        graph.SetLineColor(ROOT.kGray+3)
        graph.SetLineWidth(3)
        graph.SetMarkerColor(ROOT.kGray+3)
        graph.SetMarkerSize(1)
        graph.Draw('LP')

        gEmitylist[-1].SetLineColor(ROOT.kGray+2)
        gEmitylist[-1].SetLineWidth(3)
        gEmitylist[-1].SetLineStyle(1)
        gEmitylist[-1].SetMarkerColor(ROOT.kGray+2)
        gEmitylist[-1].SetMarkerSize(1)
        gEmitylist[-1].Draw('LP')

    DrawFrame()
    C.Print(opath + '/%s-Emittance.pdf' % (species))

    # Emittance average   
    # hFrame.GetYaxis().SetTitle('k_{p} #varepsilon_{n,x}')
    hFrame.GetYaxis().SetTitle('#LT#varepsilon_{n}#GT [100 nm]')
    hFrame.GetYaxis().SetRangeUser(exMin,exMax)

    hFrame.Draw("axis")        
    for idx,graph in enumerate(gEmitxavglist) :
        # graph.Draw('C')
        
        graph.SetLineColor(ROOT.kGray+3)
        graph.SetLineWidth(3)
        graph.SetMarkerColor(ROOT.kGray+3)
        graph.SetMarkerSize(1)
        graph.Draw('LP')

        gEmityavglist[-1].SetLineColor(ROOT.kGray+2)
        gEmityavglist[-1].SetLineWidth(3)
        gEmityavglist[-1].SetLineStyle(1)
        gEmityavglist[-1].SetMarkerColor(ROOT.kGray+2)
        gEmityavglist[-1].SetMarkerSize(1)
        gEmityavglist[-1].Draw('LP')

    DrawFrame()
    C.Print(opath + '/%s-Emittanceavg.pdf' % (species))

    # Beta        
    # hFrame.GetYaxis().SetTitle('k_{p} #varepsilon_{n,x}')
    hFrame.GetYaxis().SetTitle('#beta_{x,y} [mm]')
    hFrame.GetYaxis().SetRangeUser(betaMin,betaMax)

    hFrame.Draw("axis")        
    for idx,graph in enumerate(gBetaxlist) :
        graph.SetLineColor(ROOT.kGray+3)
        graph.SetLineWidth(3)
        graph.SetMarkerColor(ROOT.kGray+3)
        graph.SetMarkerSize(1)
        graph.Draw('LP')

        gBetaylist[-1].SetLineColor(ROOT.kGray+2)
        gBetaylist[-1].SetLineWidth(3)
        gBetaylist[-1].SetLineStyle(1)
        gBetaylist[-1].SetMarkerColor(ROOT.kGray+2)
        gBetaylist[-1].SetMarkerSize(1)
        gBetaylist[-1].Draw('LP')

    DrawFrame()
    C.Print(opath + '/%s-Beta.pdf' % (species))

    # Alpha        
    # hFrame.GetYaxis().SetTitle('k_{p} #varepsilon_{n,x}')
    hFrame.GetYaxis().SetTitle('#alpha_{x,y}')
    hFrame.GetYaxis().SetRangeUser(alphaMin,alphaMax)

    hFrame.Draw("axis")        
    for idx,graph in enumerate(gAlphaxlist) :
        graph.SetLineColor(ROOT.kGray+3)
        graph.SetLineWidth(3)
        graph.SetMarkerColor(ROOT.kGray+3)
        graph.SetMarkerSize(1)
        graph.Draw('LP')

        gAlphaylist[-1].SetLineColor(ROOT.kGray+2)
        gAlphaylist[-1].SetLineWidth(3)
        gAlphaylist[-1].SetLineStyle(1)
        gAlphaylist[-1].SetMarkerColor(ROOT.kGray+2)
        gAlphaylist[-1].SetMarkerSize(1)
        gAlphaylist[-1].Draw('LP')

    DrawFrame()
    C.Print(opath + '/%s-Alpha.pdf' % (species))

    # Gamma twiss       
    # hFrame.GetYaxis().SetTitle('k_{p} #varepsilon_{n,x}')
    hFrame.GetYaxis().SetTitle('#gamma_{x,y} [mm^{-1}]')
    hFrame.GetYaxis().SetRangeUser(gammatMin,gammatMax)

    hFrame.Draw("axis")        
    for idx,graph in enumerate(gGammaxlist) :
        graph.SetLineColor(ROOT.kGray+3)
        graph.SetLineWidth(3)
        graph.SetMarkerColor(ROOT.kGray+3)
        graph.SetMarkerSize(1)
        graph.Draw('LP')

        gGammaylist[-1].SetLineColor(ROOT.kGray+2)
        gGammaylist[-1].SetLineWidth(3)
        gGammaylist[-1].SetLineStyle(1)
        gGammaylist[-1].SetMarkerColor(ROOT.kGray+2)
        gGammaylist[-1].SetMarkerSize(1)
        gGammaylist[-1].Draw('LP')

    DrawFrame()
    C.Print(opath + '/%s-Gamma-twiss.pdf' % (species))


    # Combined energy plot
    # -----
    gPad.SetTicky(0)
    gPad.SetRightMargin(0.10)
    hFrame.GetYaxis().SetTitle('#LTp_{z}#GT [GeV/c]')
    hFrame.GetYaxis().SetRangeUser(gMin,gMax)
    hFrame.Draw("axis")
    
    hFrame.GetYaxis().SetNdivisions(505)
        
    for idx,graph in enumerate(gGammalist) :
        graph.SetLineColor(ROOT.kGray+3)
        # graph.SetLineColor(PGlobals.elecLine)
        graph.SetLineWidth(3)
        graph.SetMarkerColor(ROOT.kGray+3)
        #graph.SetMarkerColor(PGlobals.elecLine)
        graph.SetMarkerSize(1)
        graph.Draw('LP')

    # color = ROOT.kMagenta+1
    # color = ROOT.kPink-9
    color = ROOT.kOrange+10
    for idx,graph in enumerate(gGrmsrellist) :
        y = graph.GetY()
        yp = np.empty(NP)
        
        graph2 = gGrmsavglist[idx]
        y2 = graph2.GetY()
        yp2 = np.empty(NP)
        
        for j in range(NP) :
            yp[j] = ((gMax-gMin)/(grrMax-grrMin)) * (y[j]-grrMin) + gMin
            y2[j] = y2[j]/10.0
            yp2[j] = ((gMax-gMin)/(grrMax-grrMin)) * (y2[j]-grrMin) + gMin

            
        graphp = TGraph(NP,T,yp)
        graphp.SetLineColor(color)
        graphp.SetLineWidth(3)
        graphp.SetMarkerColor(color)
        graphp.SetMarkerSize(1)
        graphp.Draw('LP')
        
        graphp2 = TGraph(NP,T,yp2)
        graphp2.SetLineColor(color)
        graphp2.SetLineWidth(3)
        graphp2.SetLineStyle(2)
        graphp2.SetMarkerStyle(24)
        graphp2.SetMarkerColor(color)
        graphp2.SetMarkerSize(1)
        graphp2.Draw('L')
        
    DrawFrame()

    axis = TGaxis(gPad.GetUxmax(),gPad.GetUymin(),gPad.GetUxmax(),
			      gPad.GetUymax(),grrMin,grrMax,505,"+LS")

    axis.SetLineWidth(2)
    axis.SetLineColor(color)
    axis.SetLabelColor(color)
    axis.SetLabelFont(43)
    axis.SetLabelSize(36)
    axis.SetTitleFont(43)
    axis.SetTitleSize(38)
    axis.SetTitleOffset(0.4)
    axis.SetTickLength(0.01)
    axis.SetTitle('#sigma_{#gamma}/#bar{#gamma} [%]')
    axis.CenterTitle()
    axis.SetTitleColor(color)
    axis.Draw()
    
    C.Print(opath + '/%s-Gammaandrelspread.pdf' % (species))

    
    hFrame.GetYaxis().SetTitle('Charge')
    hFrame.GetYaxis().SetRangeUser(cMin,cMax)
    hFrame.Draw("axis")
    
    for idx,graph in enumerate(gChargelist) :
        graph.SetLineWidth(3)
        graph.Draw('L')

    DrawFrame()
    C.Print(opath + '/%s-Charge.pdf' % (species))

    
if __name__ == '__main__':
    main()

    
