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
Fx = np.zeros_like(data)
Fy = np.zeros_like(data)
for it in timeStep:
    tmp = np.gradient(data[:,:,it])
    Fx[:,:,it] = -tmp[0]
    Fy[:,:,it] = -tmp[1]

##==Calculate magnitude and angle
print("Calculating field magnitude...")
Fr = np.sqrt(Fx*Fx + Fy*Fy)
print("Calculating field direction...")
Ft = np.arctan2(Fy,Fx)

##==Calculate FFT
print("Calculating FFT of field magnitude...")
# pdb.set_trace()
Fr_fft = np.zeros((Fr.shape[0],
                   Fr.shape[1]//2 + 1,
                   Fr.shape[2]))
for it in timeStep:
    # tmp = Fr[:,:,it]
    # Fr_fft[:,:,it] = 20*np.log10(np.fft.rfft2(tmp))
    tmp = np.fft.rfft2(Fr[:,:,it])
    tmp /= np.max(tmp)
    Fr_fft[:,:,it] = 10*np.log10(tmp*tmp)
# Fr_fft = 20*np.log10(Fr_fft)
# pdb.set_trace()

print("Creating image...")
image = Fr_fft[:,:,ntMax-1]
maxAbs = np.nanmax(abs(image))
vmin = np.nanmin(image)
vmax = np.nanmax(image)
plt.pcolormesh(image.T,cmap='Spectral',vmin=vmin,vmax=vmax,rasterized=True)
plt.colorbar()
print("Done")
savePath = os.path.join(homePath,basePath,projPath,plotPath,'Er_fft-TEST.pdf')
print("Saving",savePath,"...")
plt.savefig(savePath,bbox_inches='tight',dpi=400)
print("Done")

if 0:
    print("Creating plot...")
    plt.plot(xg,Fx[:,plane['ny']//2])
    savePath = os.path.join(homePath,basePath,projPath,plotPath,
                            'Ex_mid-TEST.pdf')
    print("Saving",savePath,"...")
    plt.savefig(savePath,bbox_inches='tight',dpi=400)
    print("Done")

if 0:
    print("Creating image...")
    maxAbs = np.nanmax(abs(Fy))
    vmin = -maxAbs
    vmax = +maxAbs
    plt.pcolormesh(Fy.T,cmap='seismic',vmin=vmin,vmax=vmax,rasterized=True)
    plt.colorbar()
    print("Done")
    savePath = os.path.join(homePath,basePath,projPath,plotPath,'Ey-TEST.pdf')
    print("Saving",savePath,"...")
    plt.savefig(savePath,bbox_inches='tight',dpi=400)
    print("Done")

# pdb.set_trace()
