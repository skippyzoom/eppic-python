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
#-->Write an eppic_io function to set default values
print("Reading parameter file...")
params = eppic_io.read_parameters(path=wd)
print("Done")

##==Choose time steps to plot
ntMax = eppic_io.calc_timesteps(path=wd)
timeStep = [1,ntMax-1]

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
y0 = 0
yf = plane['ny']

##==Spatial range for left and right call-outs
zWidth = 64
x0Lz = int(xf/4 - zWidth)
xfLz = int(xf/4 + zWidth)
x0Rz = int(3*xf/4 - zWidth)
xfRz = int(3*xf/4 + zWidth)
y0Lz = int(yf/2 - zWidth)
yfLz = int(yf/2 + zWidth)
y0Rz = y0Lz
yfRz = yfLz

##==Read data file
print("Reading data...")
with h5py.File(dataFile,'r') as f:
    den = f['/den1'][:]
print("Done")

##==Adjust data
# den = np.flipud(den)
den = np.rot90(den,k=3)
den = den[x0:xf,y0:yf]

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
fig = plt.figure()
grid = plt.GridSpec(2,2, wspace=0.4, hspace=0.3)
main = fig.add_subplot(grid[0,0:])
botL = fig.add_subplot(grid[1,:1])
botR = fig.add_subplot(grid[1,1:])

im = main.pcolormesh(xg[x0:xf],yg[y0:yf],
                     den[x0:xf,y0:yf].T,
                     vmin=vmin,vmax=vmax,
                     rasterized=True)
main.set_xlabel('Zonal [m]')
main.set_ylabel('Vertical [m]')
main.set_xticks(plane['dx']*np.linspace(x0,xf,5))
main.set_yticks(plane['dy']*np.linspace(y0,yf,5))
main.set_aspect('equal')
fig.suptitle('Density')
fig.subplots_adjust(right=0.80)
cbar_ax = fig.add_axes([0.85,0.15,0.05,0.70])
fig.colorbar(im,cax=cbar_ax).set_label('$\delta n/n0$')
main.add_patch(
    pat.Rectangle(
        (plane['dx']*x0Lz,plane['dy']*y0Lz),
        zWidth,zWidth,
        fill=False,linewidth=1,edgecolor='white'))
main.add_patch(
    pat.Rectangle(
        (plane['dx']*x0Rz,plane['dy']*y0Rz),
        zWidth,zWidth,
        fill=False,linewidth=1,edgecolor='white'))

Lz = botL.pcolormesh(xg[x0Lz:xfLz],yg[y0Lz:yfLz],
                     den[x0Lz:xfLz,y0Lz:yfLz].T,
                     vmin=vmin,vmax=vmax,
                     rasterized=True)
botL.set_xlabel('Zonal [m]')
botL.set_ylabel('Vertical [m]')
botL.set_xticks(plane['dx']*np.linspace(x0Lz,xfLz,5))
botL.set_yticks(plane['dy']*np.linspace(y0Lz,yfLz,5))
botL.set_aspect('equal')

Rz = botR.pcolormesh(xg[x0Rz:xfRz],yg[y0Rz:yfRz],
                     den[x0Rz:xfRz,y0Rz:yfRz].T,
                     vmin=vmin,vmax=vmax,
                     rasterized=True)
botR.set_xlabel('Zonal [m]')
botR.set_ylabel('Vertical [m]')
botR.set_xticks(plane['dx']*np.linspace(x0Rz,xfRz,5))
botR.set_yticks(plane['dy']*np.linspace(y0Rz,yfRz,5))
botR.set_aspect('equal')

print("Done")

##==Save image
print("Saving",savePath,"...")
savePath = os.path.join(savePath)
plt.savefig(savePath,bbox_inches='tight',dpi=400)

print("Done")
