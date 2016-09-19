#!/usr/bin/env python
from __future__ import print_function
import h5py
from vtk import *
from ROOT import PData, PGlobals
import sys, argparse
import numpy as np


# Command argument line parser 
parser = argparse.ArgumentParser(description='3D renderer using VTK.')
parser.add_argument('sim', nargs='?', default='none', help='simulation name')
parser.add_argument('-b', action='store_true', help='run in batch mode without graphics')
parser.add_argument('-t', type=int, help='simulation time step')
#parser.add_argument('-i', type=int, dest='tstart', help='time start')
#parser.add_argument('-f', type=int, dest='tend', help='time end')
#parser.add_argument('-s', type=int, dest='tdelta', default=1,help='time delta (step between dumps)')
parser.add_argument('-z', type=float, dest='zoom', default=1,help='zoom')
parser.add_argument('-azi', type=float, dest='azimuth', default=0,help='azimuth')
parser.add_argument('-ele', type=float, dest='elevation', default=0,help='elevation')
parser.add_argument('--surf', action='store_true', default=0,help='draw surfaces')
parser.add_argument('--log', action='store_true', default=0,help='log scale')
parser.add_argument('--lden', action='store_true', default=0,help='use local density')
parser.add_argument('--test', action='store_true', default=0,help='run a direct example')
parser.add_argument('--nowin', action='store_true', default=0, help='no windows output (to run in batch)')

try:
    args = parser.parse_args()
except:
    parser.print_help()
    sys.exit(0)
    
if ("none" in args.sim) & (args.test == 0) :
    parser.print_help()
    sys.exit(0)
    
# End of command line setup

# Get data files
hfl = []
if args.test :
    hfl.append(h5py.File('data/charge-plasma-000026.h5','r'))
    hfl.append(h5py.File('data/charge-beam-driver-000026.h5','r'))
    hfl.append(h5py.File('data/charge-He-electrons-000026.h5','r'))
else :
    pData = PData(args.sim)
    pData.LoadFileNames(args.t)
    for i in range(0,pData.NSpecies()) :
        hfl.append(h5py.File(pData.GetChargeFileName(i).c_str(),'r'))

nfiles = len(hfl)
# -----------


# Get data, build multi-component vtkVolume,
# Set colors and opacities for the components... 
data = []
npdata = []
Min = []
Max = []
stype = []
baseden = 1

# Multi component (single) volume 
volumeprop = vtk.vtkVolumeProperty()
#volumeprop.SetIndependentComponents(ncomp)
volumeprop.IndependentComponentsOn()
volumeprop.SetInterpolationTypeToLinear()
#
opacity = []
color = []

# Loop over the files
for i, hf in enumerate(hfl):
    
    data.append(hf.get('charge'))

    axisz = hf.get('AXIS/AXIS1')
    axisy = hf.get('AXIS/AXIS2')
    axisx = hf.get('AXIS/AXIS3')

    dz = (axisz[1]-axisz[0])/data[i].shape[2]
    dy = (axisy[1]-axisy[0])/data[i].shape[1]
    dx = (axisx[1]-axisx[0])/data[i].shape[0]

    print('\nFilename : ',hf.filename)
    print('Axis z range: [%.2f,%.2f]  Nbins = %i  dz = %.4f' % (axisz[0],axisz[1],data[i].shape[2],dz) )
    print('Axis x range: [%.2f,%.2f]  Nbins = %i  dx = %.4f' % (axisx[0],axisx[1],data[i].shape[0],dx) )
    print('Axis y range: [%.2f,%.2f]  Nbins = %i  dy = %.4f' % (axisy[0],axisy[1],data[i].shape[1],dy) )

    data[i] = np.absolute(data[i])    
    maxvalue = np.amax(data[i])
    if maxvalue < 1E-5 : continue
 
    stype.append('default')
    
    npdata.append(np.array(data[i]))
    j = len(npdata) - 1

    denavg = np.average(npdata[j],weights=npdata[j].astype(bool))
    
    Max.append(np.amax(npdata[j]))
    Min.append(0.01*Max[j])
    
    if args.test == 0 :
        if pData.GetDenMax(j) > 0 : Max[j] = pData.GetDenMax(j)
        if pData.GetDenMin(j) > 0 : Min[j] = pData.GetDenMin(j) 
    
    if Min[j] > Max[j] :
        Min[j] = 0.1 * Max[j]
        
    if Max[j] > maxvalue : Max[j] = maxvalue

    if "plasma" in hf.filename :
        Min[j] = 0.99 * baseden
        if args.lden :
            localden = npdata[j][npdata[j].shape[0]-10][npdata[j].shape[1]-10][npdata[j].shape[2]-10]
            print('Local density = %f' % (localden))
            Min[j] = localden
   
    def lscale(x) :
        x = x - Min[j] + 1
        x = np.where(x <= 0.0, 1.0, x)
        x = np.log10(x)
        return x
    
    if args.log :
        npdata[j] = lscale(npdata[j])
        Min[j] = lscale(Min[j])
        Max[j] = lscale(Max[j])
        maxvalue = lscale(maxvalue)
        # stop = lscale(stop)

    # Intermediate points in the density scale
    npts = 101
    stop = np.empty(npts)
    for k in range(stop.size) :
        stop[k] = Min[j] + (1./(npts-1)) * k * (Max[j] - Min[j])

        
    np.set_printoptions(precision=3)
    print('Maximum = %.2f' % maxvalue)
    print('Min = %.2f  Max = %.2f ' % (Min[j],Max[j]),
          '-> scale = ',['{:.2f}'.format(k) for k in stop])
    print('Average = %.2f' % denavg)
    
    # Opacity and color scales
    opacity.append(vtk.vtkPiecewiseFunction())
    color.append(vtk.vtkColorTransferFunction())
    
    if "plasma" in hf.filename :
        stype[j] = 'plasma'
        
        opacity[j].AddPoint(stop[0], 0.0)
        opacity[j].AddPoint(stop[2],0.8)
        opacity[j].AddPoint(stop[10],0.2)
        opacity[j].AddPoint(stop[50],0.8)
        opacity[j].AddPoint(stop[100],0.1)
        opacity[j].AddPoint(maxvalue, 0.0)
        
        color[j].AddRGBPoint(stop[0], 0.7, 0.7, 0.7)
        #color[j].AddRGBPoint(stop[0], 0.188, 0.247, 0.294)
        #color[j].AddRGBPoint(stop[1], 0.8, 0.8, 0.8)
        #color[j].AddRGBPoint(stop[50], 1.0, 1.0, 1.0)
        color[j].AddRGBPoint(stop[100], 1.0, 1.0, 1.0)
        
    elif "beam" in hf.filename :
        stype[j] = 'beam'

        opacity[j].AddPoint(Min[j], 0.0)
        #opacity[j].AddPoint(denavg, 0.5)
        opacity[j].AddPoint(Max[j], 0.9)
        #opacity[j].AddPoint(maxvalue, 0.9)
        
        color[j].AddRGBPoint(stop[0], 0.220, 0.039, 0.235)
        color[j].AddRGBPoint(stop[20], 0.390, 0.050, 0.330)
        color[j].AddRGBPoint(stop[40], 0.700, 0.200, 0.300)
        color[j].AddRGBPoint(stop[100], 1.00, 1.00, 0.20)
        #color[j].AddRGBPoint(maxvalue, 1.00, 1.00, 0.20)
        
    elif "He-electrons" in hf.filename:
        stype[j] = 'witness'
        
        opacity[j].AddPoint(stop[0], 0.0)
        opacity[j].AddPoint(stop[1], 0.4)
        opacity[j].AddPoint(stop[5], 0.6)
        opacity[j].AddPoint(stop[10], 0.8)
        #opacity[j].AddPoint(stop[3], 0.9)
        #opacity[j].AddPoint(stop[4], 0.9)
        opacity[j].AddPoint(stop[100], 0.98)
        opacity[j].AddPoint(stop[100], 0.98)

        color[j].AddRGBPoint(stop[0], 0.220, 0.039, 0.235)
        color[j].AddRGBPoint(stop[1], 0.627, 0.125, 0.235)
        color[j].AddRGBPoint(stop[10], 0.700, 0.200, 0.300)
        #color[j].AddRGBPoint(stop[2], 0.700, 0.200, 0.300)
        #color[j].AddRGBPoint(stop[3], 0.700, 0.200, 0.300)
        #color[j].AddRGBPoint(stop[4], 1.00, 1.00, 0.20)
        color[j].AddRGBPoint(stop[100], 1.00, 1.00, 0.20)


    volumeprop.SetColor(i,color[j])
    volumeprop.SetScalarOpacity(i,opacity[j])
    volumeprop.ShadeOff(i)
    #volumeprop.ShadeOn(i)
    

# Add data components as a 4th dimension 
# npdatamulti = np.stack((npdata[:]),axis=3)
# Alternative way compatible with earlier versions of numpy 
npdatamulti = np.concatenate([aux[...,np.newaxis] for aux in npdata], axis=3)
print('\nShape of the multi-component array: ', npdatamulti.shape,' Type: ',npdatamulti.dtype)

# For VTK to be able to use the data, it must be stored as a VTKimage.
# vtkImageImport imports raw data and stores it in the image.
dataImport = vtk.vtkImageImport()
dataImport.SetImportVoidPointer(npdatamulti)
dataImport.SetDataScalarTypeToFloat()
# Number of scalar components
dataImport.SetNumberOfScalarComponents(len(npdata))
# The following two functions describe how the data is stored
# and the dimensions of the array it is stored in.
dataImport.SetDataExtent(0, npdatamulti.shape[2]-1, 0, npdatamulti.shape[1]-1, 0, npdatamulti.shape[0]-1)
dataImport.SetWholeExtent(0, npdatamulti.shape[2]-1, 0, npdatamulti.shape[1]-1, 0, npdatamulti.shape[0]-1)
dataImport.SetDataSpacing(dz,dy,dx)
dataImport.SetDataOrigin(0.0,axisy[0],axisx[0])
dataImport.Update()

mapper = vtk.vtkGPUVolumeRayCastMapper()
mapper.SetAutoAdjustSampleDistances(1)
#mapper.SetSampleDistance(0.1)
#mapper.SetBlendModeToMaximumIntensity();

# Add data to the mapper
mapper.SetInputConnection(dataImport.GetOutputPort())

# The class vtkVolume is used to pair the previously declared volume
# as well as the properties to be used when rendering that volume.
volume = vtk.vtkVolume()
volume.SetMapper(mapper)
volume.SetProperty(volumeprop)

# Surfaces
threshold = []
dmc = []
mapper2 = []
actor = []

dataImport2 = []


# Loop over the data components
for i in range(0,len(npdata)):
    if (args.surf) & ( ("beam" in stype[i]) ) :

        axisz = hf.get('AXIS/AXIS1')
        axisy = hf.get('AXIS/AXIS2')
        axisx = hf.get('AXIS/AXIS3')

        dz = (axisz[1]-axisz[0])/data[i].shape[2]
        dy = (axisy[1]-axisy[0])/data[i].shape[1]
        dx = (axisx[1]-axisx[0])/data[i].shape[0]
        
        dataImport2.append(vtk.vtkImageImport())
        j = len(dataImport2)-1
        dataImport2[j].SetImportVoidPointer(npdata[i])
        dataImport2[j].SetDataScalarTypeToFloat()
        # Number of scalar components
        dataImport2[j].SetNumberOfScalarComponents(1)
        # The following two functions describe how the data is stored
        # and the dimensions of the array it is stored in.
        dataImport2[j].SetDataExtent(0, npdata[i].shape[2]-1, 0, npdata[i].shape[1]-1, 0, npdata[i].shape[0]-1)
        dataImport2[j].SetWholeExtent(0, npdata[i].shape[2]-1, 0, npdata[i].shape[1]-1, 0, npdata[i].shape[0]-1)
        dataImport2[j].SetDataSpacing(dz,dy,dx)
        dataImport2[j].SetDataOrigin(0.0,axisy[0],axisx[0])
        dataImport2[j].Update()

        maxvalue = np.amax(npdata[i])

        threshold.append(vtkImageThreshold())
        threshold[j].SetInputConnection(dataImport2[j].GetOutputPort())
        threshold[j].ThresholdBetween(0.45*maxvalue,0.55*maxvalue)
        threshold[j].ReplaceInOn()
        threshold[j].SetInValue(0.5*maxvalue)  # set all values in range to 1
        threshold[j].ReplaceOutOn()
        threshold[j].SetOutValue(0)  # set all values out range to 0
        threshold[j].Update()

        dmc.append(vtk.vtkDiscreteMarchingCubes())
        dmc[j].SetInputConnection(threshold[j].GetOutputPort())
        #dmc[j].SetInputConnection(dataImport[i].GetOutputPort())
        dmc[j].GenerateValues(1, 0.5*maxvalue, 0.5*maxvalue)
        dmc[j].Update()

        mapper2.append(vtk.vtkPolyDataMapper())
        mapper2[j].SetInputConnection(dmc[j].GetOutputPort())
        mapper2[j].SetLookupTable(color[i])
        mapper2[j].SetColorModeToMapScalars()
         
        actor.append(vtk.vtkActor())
        actor[j].SetMapper(mapper2[j])
        actor[j].GetProperty().SetOpacity(0.8)

planeClip = vtk.vtkPlane()
planeClip.SetOrigin((axisz[0]+axisz[1])/2.0-axisz[0],0.0,0.0)
planeClip.SetNormal(0.0, 0.0, -1.0)
#mapper.AddClippingPlane(planeClip)

light = vtk.vtkLight()
light.SetColor(1.0, 0.0, 0.0)
light.SwitchOn()
light.SetIntensity(1)
# renderer.AddLight(light)

renderer = vtk.vtkRenderer()
# Set background  
renderer.SetBackground(0,0,0)
#renderer.SetBackground(0.1,0.1,0.1)
# renderer.TexturedBackgroundOn()
# Other colors 
# nc = vtk.vtkNamedColors()
# renderer.SetBackground(nc.GetColor3d('MidnightBlue'))

# Add the volume to the renderer ...
renderer.AddVolume(volume)

for j in range(len(actor)) :
    renderer.AddActor(actor[j])

# Set window
window = vtk.vtkRenderWindow()                    
# ... and set window size.
window.SetSize(1280, 800)
if args.nowin :
    window.SetOffScreenRendering(1)

# add renderer
window.AddRenderer(renderer)

interactor = vtk.vtkRenderWindowInteractor()
#style = vtkInteractorStyleTrackballCamera();
#interactor.SetInteractorStyle(style);
if args.nowin == 0 :
    interactor.SetRenderWindow(window)

# We'll zoom in a little by accessing the camera and invoking a "Zoom"
# method on it.
renderer.ResetCamera()
camera = renderer.GetActiveCamera()
camera.Zoom(args.zoom)
#camera.Roll(45)
camera.Azimuth(args.azimuth)
camera.Elevation(args.elevation)

window.Render()

if args.nowin == 0 :
    interactor.Initialize()
    interactor.Start()

# Output file
if args.test :
    foutname = './snapshot3d.png' 
else :
    foutname = './%s/Plots/Snapshots3D/Snapshot3D-%s_%i.png' % (pData.GetPath(),args.sim,args.t)
PGlobals.mkdir(foutname)

# Write an EPS file.
# exp = vtk.vtkGL2PSExporter()  # Not working with openGL2 yet
# exp.SetRenderWindow(window)
# exp.SetFilePrefix("screenshot")
# exp.DrawBackgroundOn()
# exp.Write()

# Write to PNG file
w2if = vtk.vtkWindowToImageFilter()
w2if.SetInput(window)
w2if.Update();
 
writer = vtk.vtkPNGWriter()
writer.SetFileName(foutname)
writer.SetInputConnection(w2if.GetOutputPort())
writer.Write()

print('\n%s has been created.' % (foutname) )

# END
