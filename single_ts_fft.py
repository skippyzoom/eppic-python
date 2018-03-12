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
dataPath = ''

##==Declare plotting preferences
plotPath = 'python_images'
plotType = 'pdf'
plotName = dataName+'.'+plotType

##==Set up standard path info
# homePath = os.path.expanduser('~')
homePath = '/projectnb/eregion/may/'
# basePath = 'Research'
# basePath = 'Documents/BU/research/Projects'
basePath = 'Stampede_runs'
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
# pdb.set_trace()
# data_fft = np.zeros((data.shape[0],
#                      data.shape[1]//2 + 1,
#                      data.shape[2]))
maxDim = np.max([data.shape[0],data.shape[1]])
tmp = np.fft.rfft2(data[:,:,0],s=[maxDim,maxDim])
# pdb.set_trace()
data_fft = np.zeros((tmp.shape[0],
                    tmp.shape[1],
                    data.shape[2]))
for it in timeStep:
    tmp = np.fft.rfft2(data[:,:,it],s=[maxDim,maxDim])
    tmp = np.fft.fftshift(tmp)
    tmp /= np.max(tmp)
    data_fft[:,:,it] = 10*np.log10(tmp*tmp)

print("Creating image...")
image = data_fft[:,:,ntMax-1]
maxAbs = np.nanmax(abs(image))
vmin = np.nanmin(image)
vmax = np.nanmax(image)
plt.pcolormesh(image.T,cmap='Spectral',vmin=vmin,vmax=vmax,rasterized=True)
plt.colorbar()
plt.set_aspect('equal')
print("Done")
savePath = os.path.join(homePath,basePath,projPath,plotPath,
                        dataName+'_fft-TEST.pdf')
print("Saving",savePath,"...")
plt.savefig(savePath,bbox_inches='tight',dpi=400)
print("Done")
