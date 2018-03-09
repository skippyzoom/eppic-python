#!/usr/bin/env python
import os
import h5py
import numpy as np
import matplotlib.pyplot as plt
# from matplotlib.pyplot import figure,draw
import pdb
import eppic_io

##==Declare project path
projPath = 'parametric_wave/run005/'

##==Declare data name and directory
dataName = 'den1'
dataPath = 'data/eppic/'

##==Declare plotting preferences
cmap = 'plasma'
plotPath = 'python_images'
plotType = 'pdf'
plotName = dataName+'.'+plotType

##==Set up standard path info
homePath = os.path.expanduser('~')
basePath = 'Documents/BU/research/Projects'
fileName = 'parallel004992.h5'
dataFile = os.path.join(homePath,basePath,projPath,dataPath,
                        'parallel',fileName)
wd = os.path.join(homePath,basePath,projPath,dataPath)
savePath = os.path.join(homePath,basePath,projPath,plotPath,plotName)

##==Read parameter file
print("Reading parameter file...")
params = eppic_io.read_parameters(path=wd)
print("Done")

##==Choose time steps to plot
ntMax = eppic_io.calc_timesteps(path=wd)
timeStep = [1,ntMax-1]

##==Choose spatial range to plot
x0 = 0
xf = params['nx']*params['nsubdomains']//params['nout_avg']
y0 = 0
yf = params['ny']//params['nout_avg']

##==Read data file
print("Reading data...")
with h5py.File(dataFile,'r') as f:
    den = f['/den1'][:]
print("Done")

##==Adjust data
den = den[x0:xf,y0:yf]
den = np.flipud(den)

##==Swap data limits
x0,xf,y0,yf = y0,yf,x0,xf

##==Create axis vectors
xg = params['dx']*np.linspace(x0,xf,xf-x0)
yg = params['dy']*np.linspace(y0,yf,yf-y0)

##==Set data limits
#-->Write a function to create a dict of graphics
#   preferences, given the name of a data quantity?
maxAbs = np.nanmax(abs(den))
vmin = -maxAbs
vmax = +maxAbs

##==Create image
print("Creating image...")
fg = plt.figure()
ax = fg.gca()
hi = ax.pcolormesh(xg,yg,den,
                   vmin=vmin,vmax=vmax,
                   rasterized=True)
ax.set_xlabel('Zonal [m]')
ax.set_ylabel('Vertical [m]')
ax.set_xticks(params['dx']*np.linspace(x0,xf,5))
ax.set_yticks(params['dy']*np.linspace(y0,yf,5))
ax.set_title('Density')
ax.set_aspect('equal')
fg.colorbar(hi).set_label('$\delta n/n_0$')
print("Done")

##==Save image
print("Saving",savePath,"...")
savePath = os.path.join(savePath)
plt.savefig(savePath,bbox_inches='tight')

print("Done")
