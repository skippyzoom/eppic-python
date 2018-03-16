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

##==Declare rotation, if any
k_rot = 3

##==Declare plotting preferences
plotPath = 'python_images'
plotType = 'pdf'
plotName = dataName+'.'+plotType

##==Set up standard path info
homePath = os.path.expanduser('~')
# basePath = 'Research'
basePath = 'Documents/BU/research/Projects'
wd = os.path.join(homePath,basePath,projPath,dataPath)
savePath = os.path.join(homePath,basePath,projPath,plotPath,plotName)

##==Choose time step to plot
ntMax = eppic_io.calc_timesteps(path=wd)
timeStep = ntMax-1

##==Make the image directory, if necessary
try:
    os.mkdir(os.path.join(homePath,basePath,projPath,plotPath))
except OSError:
    pass

##==Read parameter file
#-->Write an eppic_io function to set default values
print("Reading parameter file...")
params = eppic_io.read_parameters(path=wd)

##==Set up image plane
if k_rot % 2 eq 0:
    plane = {'nx':
             params['nx']*params['nsubdomains']//params['nout_avg'],
             'ny':
             params['ny']//params['nout_avg'],
             'dx':
             params['dx']}
             'dy':
             params['dy'],
else:
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

##==Read data file
if dataName[0].casefold() == 'e':
    readName = 'phi'
else:
    readName = dataName
strStep = '{:06d}'.format(params['nout']*timeStep)
fileName = 'parallel'+strStep+'.h5'
dataFile = os.path.join(homePath,basePath,projPath,dataPath,
                        'parallel',fileName)
print("Reading "+readName+"...")
with h5py.File(dataFile,'r') as f:
    fdata = f['/'+readName][:]

##==Adjust data
fdata = np.rot90(fdata,k=rot_k)
fdata = fdata[x0:xf,y0:yf]

##==Calculate gradient for efield
if dataName[0].casefold() == 'e':
    print("Calculating E from phi...")
    Ex = np.zeros_like(fdata)
    Ey = np.zeros_like(fdata)
    tmp = np.gradient(fdata)
    Ex = -tmp[0]
    Ey = -tmp[1] 
    if dataName == 'Er':
        Er = np.sqrt(Ex*Ex + Ey*Ey)
    if dataName == 'Et':
        Et = np.arctan2(Ey,Ex)

##==Create axis vectors
xdata = params['dx']*np.linspace(x0,xf,xf-x0)
ydata = params['dy']*np.linspace(y0,yf,yf-y0)

##==Set data limits
#-->Write a function to create a dict of graphics
#   preferences, given the name of a data quantity?
if dataName[0:3] == 'den':
    vmin = -np.nanmax(abs(fdata))
    vmax = +np.nanmax(abs(fdata))
    cmap = 'PiYG'
    cbar_label = '$\delta n/n0$'
if dataName == 'phi':
    vmin = -np.nanmax(abs(fdata))
    vmax = +np.nanmax(abs(fdata))
    cmap = 'bwr'
    cbar_label = '$\phi$ [V]'
if dataName == 'Ex' or dataName == 'efield_x':
    fdata = Ex
    vmin = -np.nanmax(abs(fdata))
    vmax = +np.nanmax(abs(fdata))
    cmap = 'PiYG'
    cbar_label = '$\delta E_x$ [V/m]'
if dataName == 'Ey' or dataName == 'efield_y':
    fdata = Ey
    vmin = -np.nanmax(abs(fdata))
    vmax = +np.nanmax(abs(fdata))
    cmap = 'PiYG'
    cbar_label = '$\delta E_y$ [V/m]'
if dataName == 'Er':
    fdata = Er
    vmin = 0
    vmax = np.nanmax(fdata)
    cmap = 'gist_heat'
    cbar_label = '$|\delta E|$ [V/m]'
if dataName == 'Et':
    fdata = Et
    vmin = -np.pi
    vmax = +np.pi
    cmap = 'hsv'
    cbar_label = '$tan^{-1}(\delta E_y,\delta E_x)$ [rad.]'

##==Create image
rasterized = True
print("Creating image...")
fig = plt.figure()
ax = fig.gca()
im = ax.pcolormesh(xdata[x0:xf],ydata[y0:yf],
                   fdata[x0:xf,y0:yf].T,
                   vmin=vmin,vmax=vmax,
                   cmap=cmap, 
                   rasterized=rasterized)
ax.set_xlim(plane['dx']*x0,plane['dx']*xf)
ax.set_ylim(plane['dy']*y0,plane['dy']*yf)
ax.set_xlabel('Zonal [m]')
ax.set_ylabel('Vertical [m]')
ax.set_xticks(plane['dx']*np.linspace(x0,xf,5))
ax.set_yticks(plane['dy']*np.linspace(y0,yf,5))
ax.set_aspect('equal')
fig.suptitle('Density')
fig.subplots_adjust(right=0.80)
cbar_ax = fig.add_axes([0.85,0.15,0.05,0.70])
fig.colorbar(im,cax=cbar_ax).set_label(cbar_label)

savePath = os.path.join(homePath,basePath,projPath,plotPath,
                        dataName+'-'+strStep+'.pdf')
print("Saving",savePath,"...")
plt.savefig(savePath,bbox_inches='tight',dpi=400)
print("Done")
