#!/usr/bin/env python
import os
import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as pat
import pdb
import eppic_io

##==Declare project path
projPath = 'parametric_wave/run005/'

##==Declare data name and directory
dataName = 'phi'
dataPath = 'data/eppic/'

##==Declare plotting preferences
plotPath = 'python_images'
plotType = 'pdf'
plotName = dataName+'.'+plotType

##==Set up standard path info
homePath = os.path.expanduser('~')
# homePath = '/projectnb/eregion/may/'
# basePath = 'Research'
basePath = 'Documents/BU/research/Projects'
# basePath = 'Stampede_runs'
# fileName = 'parallel004992.h5'
# fileName = 'parallel000000.h5'
# dataFile = os.path.join(homePath,basePath,projPath,dataPath,
#                         'parallel',fileName)
wd = os.path.join(homePath,basePath,projPath,dataPath)

##==Make the image directory, if necessary
try:
    os.mkdir(os.path.join(homePath,basePath,projPath,plotPath))
except OSError:
    pass

##==Read parameter file
#-->Write an eppic_io function to set default values
print("Reading parameter file...")
params = eppic_io.read_parameters(path=wd)

##==Choose time steps to plot
ntMax = eppic_io.calc_timesteps(path=wd)
timeStep = np.arange(ntMax)
nt = len(timeStep)

##==Set up image plane
plane = {'nx':
         params['ny']//params['nout_avg'],
         'ny':
         params['nx']*params['nsubdomains']//params['nout_avg'],
         'dx':
         params['dy'],
         'dy':
         params['dx']}

##==Spatial range of data to extract
pWidth = 64
x0 = int(3*plane['nx']/4 - pWidth)
xf = int(3*plane['nx']/4 + pWidth)
y0 = int(plane['ny']/4 - pWidth)
yf = int(plane['ny']/4 + pWidth)

##==Read HDF5 data
print("Reading "+dataName+"...")
data = np.zeros((xf-x0,yf-y0,nt),
                order='F')
strStep = []
for it,ts in enumerate(timeStep):
    strStep.append('{:06d}'.format(params['nout']*ts))
    fileName = 'parallel'+strStep[it]+'.h5'
    dataFile = os.path.join(homePath,basePath,projPath,dataPath,
                            'parallel',fileName)
    print("Reading",fileName)
    with h5py.File(dataFile,'r') as f:
        tmp = np.array(f['/'+dataName])
        tmp = np.rot90(tmp,k=3)
        tmp = tmp[x0:xf,y0:yf]
        data[:,:,it] = tmp

##==Create axis vectors
xg = plane['dx']*np.linspace(x0,xf,xf-x0)
yg = plane['dy']*np.linspace(y0,yf,yf-y0)

##==Calculate gradient: E = -Grad[phi]
print("Calculating field components...")
Fx = np.zeros_like(data)
Fy = np.zeros_like(data)
for it in np.arange(nt):
    tmp = np.gradient(data[:,:,it])
    Fx[:,:,it] = -tmp[0]
    Fy[:,:,it] = -tmp[1]

##==Calculate magnitude and angle
print("Calculating field magnitude...")
Fr = np.sqrt(np.square(Fx) + np.square(Fy))
print("Calculating field direction...")
Ft = np.arctan2(Fy,Fx)

##==Calculate spatial mean field
patch = np.zeros(nt)
for it in np.arange(nt):
    patch[it] = np.mean(Fr[:,:,it])
# patch = np.convolve(patch,np.ones((5,))/5,mode='same')

##==Create line plot of spatial mean
print("Creating plot...")
yMin = 0
yMax = +6e-3
fig = plt.figure()
ax = fig.gca()
ax.plot(1e3*params['dt']*params['nout']*timeStep,
        patch,'b.')
ax.set_ylim([yMin,yMax])
ax.set_xlabel('Time [ms]')
ax.set_ylabel('$|E|$ [V/m]')
savePath = os.path.join(homePath,basePath,projPath,plotPath,
                        'Er_patch_mean.pdf')
print("Saving",savePath,"...")
plt.savefig(savePath,bbox_inches='tight',dpi=400)

