#!/usr/bin/env python

import numpy as np
import scipy.constants as ct
import os, sys, argparse


def parse_args():

    # Command argument line parser 
    parser = argparse.ArgumentParser(description='convertfile.py: ---')
    parser.add_argument('filename', nargs='?', default='none', help='file name')
    parser.add_argument('-NP', type=int, dest='NP', default=-1, help='Number of macroparticles')

    args = parser.parse_args()

    if ('none' in args.filename) :
        parser.print_help()
        sys.exit(0)
    
    return args

def main():

    args = parse_args()

    filename = args.filename
    data = np.load(filename)

    Z0  = data.f.y
    X0  = data.f.x
    Y0  = data.f.z
    PZ0 = data.f.py
    PX0 = data.f.px
    PY0 = data.f.pz
    W0  = data.f.w
    Q0  = -ct.e/ct.pico

    NP    = len(X0)

    ccenter = 1.403424e-4 # meters
 
    # Select a sub-sample
    NPS = args.NP
    if NPS == -1 :
        NPS = NP
    elif NPS >= NP :
        NPS = NP

    DPS = int(NP/NPS)
    plist = range(0,NP,DPS)
    NPS = len(plist)
        
        
    print('N sample = %i' % (NPS))

    fileout = filename.replace('.npz','.txt')
    fout = open(fileout, 'w')
    for index in plist :
        fout.write('%14e  %14e  %14e  %14e  %14e  %14e  %14e\n' % ((Z0[index],X0[index]-ccenter,Y0[index]-ccenter,PZ0[index],PX0[index],PY0[index],W0[index]*Q0*(float(NP/NPS))) ))
    
    fout.close()

if __name__ == '__main__':
    main()
