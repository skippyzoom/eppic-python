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
# basePath = 'Research'
basePath = 'Documents/BU/research/Projects'
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
# timeStep = [1,ntMax-1]
timeStep = np.arange(ntMax)

##==Set up image plane
plane = {'nx':
         params['ny']//params['nout_avg'],
         'ny':
         params['nx']*params['nsubdomains']//params['nout_avg'],
         'dx':
         params['dy'],
         'dy':
         params['dx']}

##==Spatial range for full plot
x0 = 0
xf = plane['nx']
y0 = plane['ny']//2
yf = plane['ny']

##==Read HDF5 data
print("Reading "+dataName+"...")
# with h5py.File(dataFile,'r') as f:
#     data = f['/'+dataName][:]
data = np.zeros((xf-x0,yf-y0,len(timeStep)),
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

##==Adjust data
# data = np.rot90(data,k=3)
# data = data[x0:xf,y0:yf]

##==Create axis vectors
xg = plane['dx']*np.linspace(x0,xf,xf-x0)
yg = plane['dy']*np.linspace(y0,yf,yf-y0)

##==Calculate gradient: E = -Grad[phi]
print("Calculating field components...")
Fx = np.zeros(data.shape)
Fy = np.zeros(data.shape)
for it in timeStep:
    tmp = np.gradient(data[:,:,it])
    Fx[:,:,it] = -tmp[0]
    Fy[:,:,it] = -tmp[1]

##==Calculate magnitude and angle
print("Calculating field magnitude...")
Fr = np.sqrt(Fx*Fx + Fy*Fy)
print("Calculating field direction...")
Ft = np.arctan2(Fy,Fx)

if 0:
    print("Creating plot...")
    plt.plot(xg,Fx[:,plane['ny']//2])
    savePath = os.path.join(homePath,basePath,projPath,plotPath,
                            'Ex_mid-TEST.pdf')
    print("Saving",savePath,"...")
    plt.savefig(savePath,bbox_inches='tight',dpi=400)
    print("Done")

# print("Creating image...")
# maxAbs = np.nanmax(abs(Fy))
# vmin = -maxAbs
# vmax = +maxAbs
# plt.pcolormesh(Fy.T,cmap='seismic',vmin=vmin,vmax=vmax,rasterized=True)
# plt.colorbar()
# print("Done")
# savePath = os.path.join(homePath,basePath,projPath,plotPath,'Ey-TEST.pdf')
# print("Saving",savePath,"...")
# plt.savefig(savePath,bbox_inches='tight',dpi=400)
# print("Done")

# pdb.set_trace()
