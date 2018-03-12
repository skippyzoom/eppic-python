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
timeStep = [1,ntMax-1]
# timeStep = np.arange(ntMax)

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

##==Calculate FFT
print("Calculating FFT of "+dataName+"...")
maxDim = np.max([data.shape[0],data.shape[1]])
tmp = np.fft.rfft2(data[:,:,0],s=[maxDim,maxDim])
fnx = tmp.shape[0]
fny = tmp.shape[1]
fnt = data.shape[2]
data_fft = np.zeros((fnx,fny,fnt))
dcw = 4
for it in np.arange(len(timeStep)):
    tmp = np.fft.rfft2(data[:,:,it],s=[maxDim,maxDim])
    tmp = tmp.real
    tmp = np.fft.fftshift(tmp,axes=0)
    tmp[fnx//2-dcw:fnx//2+dcw,:] = np.float64(1e-6)
    test = tmp[fnx//2-dcw:fnx//2+dcw,:]
    tmp /= np.nanmax(tmp)
    data_fft[:,:,it] = 10*np.log10(np.square(tmp))
    tmp = None

print("Creating image...")
image = data_fft[:,:,fnt-1]
vmin = -30
vmax = 0
fig = plt.figure()
ax = fig.gca()
im = ax.pcolormesh(image.T,
                   cmap='Spectral_r',
                   vmin=vmin,vmax=vmax,
                   rasterized=True)
ax.set_aspect('equal')
ax.set_position([0.10,0.10,0.70,0.35])
cbar_ax = fig.add_axes([0.85,0.15,0.02,0.35])
fig.colorbar(im,cax=cbar_ax).set_label('Spectral Power [dB]')

savePath = os.path.join(homePath,basePath,projPath,plotPath,
                        dataName+'_fft-TEST.pdf')
print("Saving",savePath,"...")
plt.savefig(savePath,bbox_inches='tight',dpi=400)
print("Done")
