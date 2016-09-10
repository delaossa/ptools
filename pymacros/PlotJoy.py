#!/usr/bin/env python

from ROOT import PData, PGlobals, gStyle, TObject, TFile, TCanvas, TH1F
import sys, argparse

def main():

    parser = argparse.ArgumentParser(description='Makes a Joy Division style plot.')
    parser.add_argument('sim', help='simulation name')
    parser.add_argument('-b', action='store_true', help='run in batch mode without graphics')
    parser.add_argument('-t', type=int, help='timestep')
    parser.add_argument('-n', type=int, dest="nlines", default=50, help='number of lines')
    parser.add_argument('-max', type=float, dest="max", default=100, help='maximum density')
    parser.add_argument('-f', type=float, dest="factor", default=100, help='vertical factor')

    try:
        args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)
        
    pData = PData(args.sim)

    pData.LoadFileNames(args.t)
    pData.PrintData()

    PGlobals.Initialize()

    gStyle.SetCanvasColor(1)
    gStyle.SetFrameFillColor(1)
    gStyle.SetHistLineColor(0)
    gStyle.SetHistLineWidth(2)

    gStyle.SetPadLeftMargin(0.20)
    gStyle.SetPadRightMargin(0.20)
    gStyle.SetPadTopMargin(0.20)
    gStyle.SetPadBottomMargin(0.20)
    
    ofile = TFile("joy.root","RECREATE")
    NbinsX1 = pData.GetX1N()

    hDen1D = []

    N = args.nlines
    top = args.max
    shift = 0.2
    factor = args.factor
    max = -999
    max2 = factor * (NbinsX1 * shift) 
    
    C = TCanvas("C","",800,1000)

    j =  (pData.GetX3iMin() + pData.GetX3iMax()) / 2

    for i in range(pData.GetX1iMin(),pData.GetX1iMax()) :
        
        hName = "hDen1D_%i" % i
        if pData.Is3D() :
            hDen1D.append(pData.GetH1SliceX3D(pData.GetChargeFileName(0).c_str(),"charge",i,i,j,j,"avg"))

        hDen1D[i].SetName(hName)

        for k in range (1,hDen1D[i].GetNbinsX()+1) :
            
            if(hDen1D[i].GetBinContent(k)<top) :
	        hDen1D[i].SetBinContent(k,hDen1D[i].GetBinContent(k) + i*shift)
            else :
	        hDen1D[i].SetBinContent(k,top + i*shift)
                
        if(max<hDen1D[i].GetMaximum()) :
            max = hDen1D[i].GetMaximum()
                
        #        hDen1D[i].Write(hName,TObject.kOverwrite)
                                
    hFrame = hDen1D[0].Clone("hFrame")
    hFrame.Reset()
    # hFrame.GetYaxis().SetRangeUser(0.0,max2)
    hFrame.GetYaxis().SetRangeUser(0.0,1.05 * max)
                    
    hFrame.Draw("A")
    for i in range(pData.GetX1iMax()-1,pData.GetX1iMin(), -NbinsX1/N):
        hDen1D[i].SetFillColor(1)
        hDen1D[i].Draw("FC same")
                        
    PGlobals.imgconv(C,"joy","pdf")


if __name__ == '__main__':
    main()
