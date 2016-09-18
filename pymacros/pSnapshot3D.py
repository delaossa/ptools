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

ncomp = len(hfl)
# -----------


# Get data, build multi-component vtkVolume,
# Set colors and opacities for the components... 
data = []
npdata = []
Min = []
Max = []

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
    if maxvalue < 1E-4 : continue
    
    npdata.append(np.array(data[i]))
    j = len(npdata) - 1
    
    maxvalue = np.amax(npdata[j])
    minvalue = 0.1 * baseden

    if "He-electrons" in hf.filename :
        minvalue = 1E-4 * baseden
    
    if minvalue > maxvalue :
        minvalue = 0.1 * maxvalue
    
    if args.log :
        npdata[j] = np.log10(npdata[j] - minvalue + 1)
        baseden = np.log10(baseden - minvalue + 1)
        minvalue = 0.0
        maxvalue = np.amax(npdata[j])

    Max.append(maxvalue)
    Min.append(minvalue)
        
    print('Minimum value = %.2e  Maximum = %.2e' % (minvalue,maxvalue))
    print('Shape of the array: ', npdata[j].shape,' Type: ',npdata[j].dtype)
    
    # Opacity and color scales
    opacity.append(vtk.vtkPiecewiseFunction())
    color.append(vtk.vtkColorTransferFunction())
    
    if "plasma" in hf.filename:
        localden = npdata[j][npdata[j].shape[0]-10][npdata[j].shape[1]-10][npdata[j].shape[2]-10]
        print('\nLocal density = %f' % (localden))
        if args.lden :
            baseden = localden

        #topden = 10 * baseden
        #if topden > maxvalue : topden = baseden + 0.50 * (maxvalue-baseden)
        
        opacity[j].AddPoint(minvalue, 0.0)
        opacity[j].AddPoint(baseden, 0.0)
        #opacity[j].AddPoint(topden, 0.99)
        opacity[j].AddPoint(baseden + 0.50 * (maxvalue-baseden), 0.8)
        opacity[j].AddPoint(maxvalue, 0.0)
        
        color[j].AddRGBPoint(minvalue, 0.078, 0.078, 0.078)
        color[j].AddRGBPoint(baseden, 0.188, 0.247, 0.294)
        #color[j].AddRGBPoint(topden, 1.0, 1.0, 1.0)
        color[j].AddRGBPoint(baseden + 0.50 * (maxvalue-baseden), 1.0, 1.0, 1.0)
        color[j].AddRGBPoint(maxvalue, 1.0, 1.0, 1.0)
        
    elif "beam" in hf.filename :
        # maxvalue = 10*base       
        opacity[j].AddPoint(minvalue, 0.0)
        # opacity[j].AddPoint(0.2*maxvalue, 0.2)
        # opacity[j].AddPoint(0.8*maxvalue, 1.0)
        opacity[j].AddPoint(minvalue + 1.00 * (maxvalue-minvalue), 1.0)

        color[j].AddRGBPoint(0.0, 0.220, 0.039, 0.235)
        color[j].AddRGBPoint(minvalue + 0.20 * (maxvalue-minvalue), 0.390, 0.050, 0.330)
        color[j].AddRGBPoint(minvalue + 0.40 * (maxvalue-minvalue), 0.700, 0.200, 0.300)
        color[j].AddRGBPoint(maxvalue, 1.00, 1.00, 0.20)
        
    elif "He-electrons" in hf.filename:
        opacity[j].AddPoint(minvalue, 0.0)
        opacity[j].AddPoint(minvalue + 0.01 * (maxvalue-minvalue), 0.6)
        opacity[j].AddPoint(minvalue + 0.40 * (maxvalue-minvalue), 0.8)
        opacity[j].AddPoint(minvalue + 0.60 * (maxvalue-minvalue), 0.9)
        opacity[j].AddPoint(maxvalue, 1.0)

        color[j].AddRGBPoint(minvalue, 0.220, 0.039, 0.235)
        color[j].AddRGBPoint(minvalue + 0.01 * (maxvalue-minvalue), 0.627, 0.125, 0.235)
        color[j].AddRGBPoint(minvalue + 0.40 * (maxvalue-minvalue), 0.700, 0.200, 0.300)
        color[j].AddRGBPoint(minvalue + 1.00 * (maxvalue-minvalue), 1.000, 1.000, 0.200)

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
for i in range(0,ncomp):
    if (args.surf) & ( ("beam" in hfl[i].filename) ) :

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
