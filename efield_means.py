#!/usr/bin/env python
import os
from pathlib import Path
import h5py
import numpy as np
import matplotlib.pyplot as plt
import pdb
import eppic_io

##==Set options
use_rms = False

##==Declare data name
dataName = 'phi'

##==Declare project 
projName = 'parametric_wave/run005/'

##==Declare rotation, if any
rot_k = 3

##==Declare plotting preferences
plotDir = 'python_images'
plotType = 'pdf'

##==Set up standard path info
paths = eppic_io.set_paths()
basePath = Path(paths[0])
infoPath = basePath / projName / Path(paths[1])
dataPath = infoPath / 'parallel'
savePath = basePath / projName / plotDir

##==Choose time step to plot
ntMax = eppic_io.calc_timesteps(path=str(infoPath))
timeStep = np.linspace(1,ntMax-1,4,dtype='int')
nts = len(timeStep)

##==Make the image directory, if necessary
if not savePath.exists(): savePath.mkdir()

##==Read parameter file
#-->Write an eppic_io function to set default values
print("Reading parameter file...")
params = eppic_io.read_parameters(path=str(infoPath))

##==Create plane dictionary
axes = 'xy'

##==Get image plane and data
fdata = eppic_io.imgplane(dataName,
                          axes=axes,
                          timeStep=timeStep,
                          rotate=rot_k,
                          fft_direction=0,
                          calc_gradient=1,
                          dataPath=dataPath,
                          infoPath=infoPath)
pdb.set_trace() #--> Update below for fdata

##==Calculate gradient for efield
print("Calculating E from phi...")
Ex = np.zeros_like(fdata)
Ey = np.zeros_like(fdata)
tmp = np.gradient(fdata)
Ex = -tmp[0]
Ey = -tmp[1] 
Er = np.sqrt(Ex*Ex + Ey*Ey)
Et = np.arctan2(Ey,Ex)

##==Create axis vectors
xdata = params['dx']*np.linspace(x0,xf,xf-x0)
ydata = params['dy']*np.linspace(y0,yf,yf-y0)

#------------------------#
# Reduce Ex along x axis #
#------------------------#
##==Calculate RMS or mean
if use_rms:
    image = np.sqrt(np.nanmean(np.square(Ex),axis=0))
else:
    image = np.nanmean(Ex,axis=0)

##==Create image
fig = plt.figure()
ax = fig.gca()
for it,ts in enumerate(timeStep):
    label = '{:4.1f}'.format(params['nout']*ts*params['dt']*1e3)+' ms'
    ax.plot(ydata,image[:,it],label=label.lstrip())
ax.legend(loc='best')
ax.set_xlabel('Zonal [m]')
ax.set_ylabel('$\delta E_x$ [V/m]')
plotName = 'Ex_x-mean.'+plotType
savePath = os.path.join(homePath,basePath,projPath,plotPath,plotName)
print("Saving",savePath,"...")
plt.savefig(savePath,bbox_inches='tight',dpi=400)

#------------------------#
# Reduce Ex along y axis #
#------------------------#
##==Calculate RMS or mean
if use_rms:
    image = np.sqrt(np.nanmean(np.square(Ex),axis=1))
else:
    image = np.nanmean(Ex,axis=1)

##==Create image
fig = plt.figure()
ax = fig.gca()
for it,ts in enumerate(timeStep):
    label = '{:4.1f}'.format(params['nout']*ts*params['dt']*1e3)+' ms'
    ax.plot(xdata,image[:,it],label=label.lstrip())
ax.legend(loc='best')
ax.set_xlabel('Zonal [m]')
ax.set_ylabel('$\delta E_x$ [V/m]')
plotName = 'Ex_y-mean.'+plotType
savePath = os.path.join(homePath,basePath,projPath,plotPath,plotName)
print("Saving",savePath,"...")
plt.savefig(savePath,bbox_inches='tight',dpi=400)

#------------------------#
# Reduce Ey along x axis #
#------------------------#
##==Calculate RMS or mean
if use_rms:
    image = np.sqrt(np.nanmean(np.square(Ey),axis=0))
else:
    image = np.nanmean(Ey,axis=0)

##==Create image
fig = plt.figure()
ax = fig.gca()
for it,ts in enumerate(timeStep):
    label = '{:4.1f}'.format(params['nout']*ts*params['dt']*1e3)+' ms'
    ax.plot(ydata,image[:,it],label=label.lstrip())
ax.legend(loc='best')
ax.set_xlabel('Vertical [m]')
ax.set_ylabel('$\delta E_y$ [V/m]')
plotName = 'Ey_x-mean.'+plotType
savePath = os.path.join(homePath,basePath,projPath,plotPath,plotName)
print("Saving",savePath,"...")
plt.savefig(savePath,bbox_inches='tight',dpi=400)

#------------------------#
# Reduce Ey along y axis #
#------------------------#
##==Calculate RMS or mean
if use_rms:
    image = np.sqrt(np.nanmean(np.square(Ey),axis=1))
else:
    image = np.nanmean(Ey,axis=1)

##==Create image
fig = plt.figure()
ax = fig.gca()
for it,ts in enumerate(timeStep):
    label = '{:4.1f}'.format(params['nout']*ts*params['dt']*1e3)+' ms'
    ax.plot(xdata,image[:,it],label=label.lstrip())
ax.legend(loc='best')
ax.set_xlabel('Vertical [m]')
ax.set_ylabel('$\delta E_y$ [V/m]')
plotName = 'Ey_y-mean.'+plotType
savePath = os.path.join(homePath,basePath,projPath,plotPath,plotName)
print("Saving",savePath,"...")
plt.savefig(savePath,bbox_inches='tight',dpi=400)

