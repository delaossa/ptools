#!/usr/bin/env python
from __future__ import print_function
import h5py
from vtk import *
from ROOT import PData
import sys, argparse
import numpy as np


# Command argument line parser 
parser = argparse.ArgumentParser(description='3D renderer using VTK.')
parser.add_argument('sim', help='simulation name')
parser.add_argument('-b', action='store_true', help='run in batch mode without graphics')
parser.add_argument('-t', type=int, help='timestep')
parser.add_argument('--hdf5', action='store_true', help='run a direct hdf5 example')

try:
    args = parser.parse_args()
except:
    parser.print_help()
    sys.exit(0)

hfl = []
if args.hdf5 :
    hfl.append(h5py.File('data/charge-plasma-000026.h5','r'))
    hfl.append(h5py.File('data/charge-beam-driver-000026.h5','r'))
    hfl.append(h5py.File('data/charge-He-electrons-000026.h5','r'))
else :
    pData = PData(args.sim)
    pData.LoadFileNames(args.t)
    for i in range(0,pData.NSpecies()) :
        hfl.append(h5py.File(pData.GetChargeFileName(i).c_str(),'r'))

ncomp = len(hfl)

window = vtk.vtkRenderWindow()
# ... and set window size.
window.SetSize(1280, 800)

renderer = vtk.vtkRenderer()
# Set background  
renderer.SetBackground(0,0,0)
#renderer.TexturedBackgroundOn()
# Other colors 
# nc = vtk.vtkNamedColors()
# renderer.SetBackground(nc.GetColor3d('MidnightBlue'))

data = []
npdata = []
npdatauchar = []
opacity = []
color = []

volumeprop = vtk.vtkVolumeProperty()

#volumeprop.SetIndependentComponents(ncomp)
volumeprop.IndependentComponentsOn()
volumeprop.SetInterpolationTypeToLinear()

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


    # Changing to positive integer types (particle density)
    # it is required by vtkVolumeRayCastMapper
    npdata.append(np.array(np.absolute(data[i])))
    minvalue = np.amin(npdata[i])
    maxvalue = np.amax(npdata[i])
    print('Minimum value = %.2f  Maximum = %.2f' % (minvalue,maxvalue))

    # Rescale data to span from 0 to 255 (unsigned char)
    den1 = 255.0/maxvalue

    # Zooms into low values of the scalar
    # Trick to truncate high scalar values to enhance the resolution of low values
    if "He-electrons" in hf.filename:
        den1 *= 100
    elif "plasma" in hf.filename:
        den1 *= 10
    elif "beam" in hf.filename:
        den1 *= 10

    npdata[i] = np.round(den1 * npdata[i])

    # Change data from float to unsigned char
    npdatauchar.append(np.array(npdata[i], dtype=np.uint8))
    print('Shape of the array: ', npdatauchar[i].shape,' Type: ',npdatauchar[i].dtype)
    minvalue = np.amin(npdatauchar[i])
    maxvalue = np.amax(npdatauchar[i])
    print('Minimum value = %.2f  Maximum = %.2f' % (minvalue,maxvalue))
        
    # Opacity and color scales
    opacity.append(vtk.vtkPiecewiseFunction())
    color.append(vtk.vtkColorTransferFunction())
    if "plasma" in hf.filename:
        opacity[i].AddPoint(0, 0.0)
        opacity[i].AddPoint(den1, 0.01)
        opacity[i].AddPoint(10*den1, 0.8)
        opacity[i].AddPoint(maxvalue, 1.0)

        color[i].AddRGBPoint(0, 0.078, 0.078, 0.078)
        color[i].AddRGBPoint(den1, 0.188, 0.247, 0.294)
        color[i].AddRGBPoint(maxvalue, 1.0, 1.0, 1.0)
        # other palette
        #color[i].AddRGBPoint(0.0, 0.865, 0.865, 0.865)
        #color[i].AddRGBPoint(den1, 0.2313, 0.298, 0.753)
        #color[i].AddRGBPoint(maxvalue, 1.0, 1.0, 1.0)
    elif "beam" in hf.filename :
        opacity[i].AddPoint(0, 0.0)
        opacity[i].AddPoint(maxvalue, 1.0)

        color[i].AddRGBPoint(0.0, 0.220, 0.039, 0.235)
        color[i].AddRGBPoint(0.2*maxvalue, 0.390, 0.050, 0.330)
        color[i].AddRGBPoint(0.4*maxvalue, 0.700, 0.200, 0.300)
        color[i].AddRGBPoint(1.0*maxvalue, 1.000, 1.000, 0.200)
    elif "He-electrons" in hf.filename:
        opacity[i].AddPoint(0.0, 0.0)
        opacity[i].AddPoint(1, 0.3)
        opacity[i].AddPoint(100, 0.8)
        opacity[i].AddPoint(255, 1.0)

        color[i].AddRGBPoint(0.0, 0.220, 0.039, 0.235)
        color[i].AddRGBPoint(0.01*maxvalue, 0.627, 0.125, 0.235)
        color[i].AddRGBPoint(0.10*maxvalue, 0.700, 0.200, 0.300)
        color[i].AddRGBPoint(1.00*maxvalue, 1.000, 1.000, 0.200)

    volumeprop.SetColor(i,color[i])
    volumeprop.SetScalarOpacity(i,opacity[i])
    volumeprop.ShadeOff(i)
    #volumeprop.ShadeOn(i)



# Add data components as a 4th dimension 
# npdatamulti = np.stack((npdatauchar[:]),axis=3)

# Alternative way compatible with earlier versions of numpy 
npdatamulti = np.concatenate([aux[...,np.newaxis] for aux in npdatauchar], axis=3)

print('\nShape of the multi-component array: ', npdatamulti.shape,' Type: ',npdatamulti.dtype)

# For VTK to be able to use the data, it must be stored as a VTK-image.
# This can be done by the vtkImageImport which
# imports raw data and stores it.
dataImport = vtk.vtkImageImport()

dataImport.SetImportVoidPointer(npdatamulti)

dataImport.SetDataScalarTypeToUnsignedChar()

# Number of scalar components
dataImport.SetNumberOfScalarComponents(ncomp)
# The following two functions describe how the data is stored
# and the dimensions of the array it is stored in.
dataImport.SetDataExtent(0, npdatamulti.shape[2]-1, 0, npdatamulti.shape[1]-1, 0, npdatamulti.shape[0]-1)
dataImport.SetWholeExtent(0, npdatamulti.shape[2]-1, 0, npdatamulti.shape[1]-1, 0, npdatamulti.shape[0]-1)
dataImport.SetDataSpacing(dz,dy,dx)
dataImport.SetDataOrigin(0.0,axisy[0],axisx[0])

dataImport.Update()

# Set the mapper
mapper = vtk.vtkGPUVolumeRayCastMapper()
#mapper = vtk.vtkFixedPointVolumeRayCastMapper()
#mapper = vtk.vtkSmartVolumeMapper()
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

planeClip = vtk.vtkPlane()
planeClip.SetOrigin((axisz[0]+axisz[1])/2.0-axisz[0],0.0,0.0)
planeClip.SetNormal(0.0, 0.0, -1.0)
#mapper.AddClippingPlane(planeClip)

light = vtk.vtkLight()
light.SetColor(1.0, 0.0, 0.0)
light.SwitchOn()
light.SetIntensity(1)
#renderer.AddLight(light)

# Add the volume to the renderer ...
renderer.AddVolume(volume)

window.AddRenderer(renderer)

interactor = vtk.vtkRenderWindowInteractor()
interactor.SetRenderWindow(window)
#style = vtkInteractorStyleTrackballCamera();
#interactor.SetInteractorStyle(style);

# We'll zoom in a little by accessing the camera and invoking a "Zoom"
# method on it.
renderer.ResetCamera()
renderer.GetActiveCamera().Zoom(2.0)

window.Render()

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
writer.SetFileName("screenshot.png")
writer.SetInputConnection(w2if.GetOutputPort())
writer.Write()

interactor.Initialize()
interactor.Start()
