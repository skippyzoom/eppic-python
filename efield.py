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

##==Graphics options
fft_images = False

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

timeStep = np.arange(1,ntMax,ntMax//4)
# timeStep = [1,ntMax-1]
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
y0 = 0
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

##==Create axis vectors
xg = plane['dx']*np.linspace(x0,xf,xf-x0)
yg = plane['dy']*np.linspace(y0,yf,yf-y0)

##==Calculate gradient: E = -Grad[phi]
print("Calculating field components...")
Fx = np.zeros_like(data)
Fy = np.zeros_like(data)
for it in np.arange(len(timeStep)):
    tmp = np.gradient(data[:,:,it])
    Fx[:,:,it] = -tmp[0]
    Fy[:,:,it] = -tmp[1]

##==Calculate magnitude and angle
print("Calculating field magnitude...")
Fr = np.sqrt(Fx*Fx + Fy*Fy)
print("Calculating field direction...")
Ft = np.arctan2(Fy,Fx)

##==Create line plot of spatial mean
# image = Fx[:,plane['ny']//2,:]
image = np.mean(Fr,axis=1)
print("Creating plot...")
for it in np.arange(len(timeStep)):
    fig = plt.figure()
    ax = fig.gca()
    ax.plot(xg,image)
    ax.set_ylim([-6e-3,+6e-3])
    ax.set_xlabel('Zonal [m]')
    ax.set_ylabel('$|E|$ [V/m]')
savePath = os.path.join(homePath,basePath,projPath,plotPath,
                        'Er_mean_multi.pdf')
print("Saving",savePath,"...")
plt.savefig(savePath,bbox_inches='tight',dpi=400)

##==Calculate FFT
if fft_images:
    print("Calculating FFT of field magnitude...")
    maxDim = np.max([Fr.shape[0],Fr.shape[1]])
    tmp = np.fft.rfft2(Fr[:,:,0],s=[maxDim,maxDim])
    fnx = tmp.shape[0]
    fny = tmp.shape[1]
    fnt = Fr.shape[2]
    Fr_fft = np.zeros((fnx,fny,fnt))
    dcw = 4
    for it in np.arange(len(timeStep)):
        tmp = np.fft.rfft2(Fr[:,:,it],s=[maxDim,maxDim])
        tmp = tmp.real
        tmp = np.fft.fftshift(tmp,axes=0)
        tmp[fnx//2-dcw:fnx//2+dcw,:] = np.float64(1e-6)
        test = tmp[fnx//2-dcw:fnx//2+dcw,:]
        tmp /= np.nanmax(tmp)
        Fr_fft[:,:,it] = 10*np.log10(np.square(tmp))
        tmp = None
    print("Creating images...")
    for it in np.arange(len(timeStep)):
        image = Fr_fft[:,:,it]
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
                                'Er_fft-'+strStep[it]+'.pdf')
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
