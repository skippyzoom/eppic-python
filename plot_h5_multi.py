#!/usr/bin/env python
import h5py
import numpy as np
import os, sys
from pathlib import Path
import matplotlib.pyplot as plt
import importlib.machinery
import eppic_io

## Declare quantity to plot
# dataName = 'phi'
# colorMap = 'seismic'
dataName = 'den1'
colorMap = 'plasma'
plotType = 'png'

## Set the project directory and import variables
baseDir = '/projectnb/eregion/may/Stampede_runs/'
projDir = 'quasineutral_static_dust/run009/'
path = baseDir+projDir
posixPath = Path(path)
eppic = importlib.machinery.SourceFileLoader('eppic',path+'eppic.py').load_module()

## Set default parameters
## -->This could get bloated if there are a lot
# try:
#     eppic.nz
# except AttributeError:
#     setattr(eppic, 'nz', 1)
if hasattr(eppic, 'nz'):
    pass
else:
    setattr(eppic, 'nz', eppic.nout_avg)

## Declare data plot ranges
## -->Eventually want to import these from a run-specific file
x0 = 0
xf = int(eppic.nx*eppic.nsubdomains/eppic.nout_avg)
y0 = int((eppic.ny/4)/eppic.nout_avg)
yf = int((3*eppic.ny/4)/eppic.nout_avg)

## Choose time steps to plot
ntMax = eppic_io.calc_timesteps(path)
timeStep = [1,ntMax-1]

## Loop over time steps to read
nX = int(eppic.nx*eppic.nsubdomains/eppic.nout_avg)
nY = int(eppic.ny/eppic.nout_avg)
nZ = int(eppic.nz/eppic.nout_avg)
nT = len(timeStep)
# data = np.zeros((nX,nY,nZ,nT))
data = np.zeros((nX,nY,nT),order='F')
strStep = []
# print(data.shape)
for it,ts in enumerate(timeStep):
    # strStep = '{:06d}'.format(eppic.nout*ts)
    strStep.append('{:06d}'.format(eppic.nout*ts))
    fileName = 'parallel'+strStep[it]+'.h5'
    print("Reading",fileName)
    with h5py.File(path+'parallel/'+fileName,'r') as f:
        temp = np.array(f['/'+dataName])
        # temp = np.transpose(temp)
        data[:,:,it] = temp

xg = np.linspace(x0,xf,data[x0:xf,y0:yf,0].shape[0])
yg = np.linspace(y0,yf,data[x0:xf,y0:yf,0].shape[1])
# print(data[x0:xf,y0:yf,0].shape[0])
# print(data[x0:xf,y0:yf,0].shape[1])

dataMaxAbs = np.nanmax(np.absolute(data))
# vmin = data.nanmin()
# vmax = data.nanmax()
vmin = -dataMaxAbs
vmax = dataMaxAbs
# print(vmin,vmax)
cmap = plt.get_cmap(colorMap)
fig, axes = plt.subplots(2, sharex=True)
img = axes[0].pcolormesh(xg,yg,
                         np.transpose(data[x0:xf,y0:yf,0]),
                         cmap=cmap,vmin=vmin,vmax=vmax)
axes[0].set_title(strStep[0])
axes[0].set_aspect('equal')
img = axes[1].pcolormesh(xg,yg,
                         np.transpose(data[x0:xf,y0:yf,1]),
                         cmap=cmap,vmin=vmin,vmax=vmax)
axes[1].set_title(strStep[1])
axes[1].set_aspect('equal')
fig.suptitle('Density')
fig.subplots_adjust(right=0.80)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.70])
cbar = fig.colorbar(img,cax=cbar_ax)
cbar.set_label('$\delta n/n_0$')
plt.show()

print("\n")
