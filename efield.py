#!/usr/bin/env python
import os
import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as pat
import scipy.interpolate as scint
import pdb
import eppic_io

##==Declare project path
projPath = 'parametric_wave/run005/'

##==Declare data name and directory
dataName = 'phi'
# dataPath = 'data/eppic/'
dataPath = ''

##==Graphics options
raw_images = True
calc_fft = True
calc_interp = False
fft_images = True
fft_filter = True
plot_mean_y = False
plot_mean_xy = True
if fft_images and ~calc_fft:
    calc_fft = True
if fft_filter and ~calc_fft:
    calc_fft = True

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
nt = len(timeStep)

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
data = np.zeros((xf-x0,yf-y0,nt),
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
for it in np.arange(nt):
    tmp = np.gradient(data[:,:,it])
    Fx[:,:,it] = -tmp[0]
    Fy[:,:,it] = -tmp[1]

##==Calculate magnitude and angle
print("Calculating field magnitude...")
Fr = np.sqrt(Fx*Fx + Fy*Fy)
print("Calculating field direction...")
Ft = np.arctan2(Fy,Fx)

if raw_images:
    print("Creating images...")
    for it in np.arange(nt):
        image = Fr[:,:,it]
        vmin = 0
        vmax = np.nanmax(Fr)
        fig = plt.figure()
        ax = fig.gca()
        im = ax.pcolormesh(image.T,
                           cmap='gist_heat',
                           vmin=vmin,vmax=vmax,
                           rasterized=True)
        ax.set_aspect('equal')
        ax.set_position([0.10,0.10,0.70,0.35])
        cbar_ax = fig.add_axes([0.85,0.15,0.02,0.35])
        fig.colorbar(im,cax=cbar_ax).set_label('|E| [V/m]')
        savePath = os.path.join(homePath,basePath,projPath,plotPath,
                                'Er-'+strStep[it]+'.pdf')
        print("Saving",savePath,"...")
        plt.savefig(savePath,bbox_inches='tight',dpi=400)


##==Calculate FFT
if calc_fft:
    print("Calculating FFT of field magnitude...")
    maxDim = np.max([Fr.shape[0],Fr.shape[1]])
    # tmp = np.fft.rfft2(Fr[:,:,0],s=[maxDim,maxDim])
    # tmp = np.fft.rfft2(Fr[:,:,0])
    tmp = np.fft.fft2(Fr[:,:,0])
    fnx = tmp.shape[0]
    fny = tmp.shape[1]
    fnt = Fr.shape[2]
    Fr_fft = np.zeros((fnx,fny,fnt),dtype='complex')
    # Fr_fftdB = np.copy(Fr_fft).real
    dcw = 4
    # for it in np.arange(fnt):
    #     tmp = np.fft.rfft2(Fr[:,:,it],s=[maxDim,maxDim])
    #     tmp = tmp.real
    #     tmp = np.fft.fftshift(tmp,axes=0)
    #     tmp[fnx//2-dcw:fnx//2+dcw,:] = np.float64(1e-6)
    #     test = tmp[fnx//2-dcw:fnx//2+dcw,:]
    #     tmp /= np.nanmax(tmp)
    #     Fr_fftdB[:,:,it] = 10*np.log10(np.square(tmp))
    #     tmp = None
    # kx = np.fft.fftfreq(maxDim,plane['dx'])
    # # kx = np.fft.fftshift(kx)
    # ky = np.fft.rfftfreq(maxDim,plane['dy'])
    for it in np.arange(fnt):
        Fr_fft[:,:,it] = np.fft.fft2(Fr[:,:,it])
    kx = np.fft.fftfreq(Fr_fft.shape[0],plane['dx'])
    ky = np.fft.fftfreq(Fr_fft.shape[1],plane['dy'])
    pdb.set_trace()
if fft_images:
    print("Creating images...")
    for it in np.arange(nt):
        # image = Fr_fftdB[:,:,it]
        image = Fr_fft[:,:,it]
        image = image.real
        image = np.fft.fftshift(image,axes=0)
        image[fnx//2-dcw:fnx//2+dcw,:] = np.float64(1e-6)
        image /= np.nanmax(image)
        image = 10*np.log10(np.square(image))
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

if fft_filter:
    print("Filtering FFT")
    lam = 10.0
    # tmpx = (np.where(np.fabs(kx) < 1/lam))[0]
    # tmpy = (np.where(np.fabs(ky) < 1/lam))[0]
    # Fr_fft[tmpx[0]:tmpx[len(tmpx)-1],
    #        tmpy[0]:tmpy[len(tmpy)-1],:] = np.float64(1e-6)
    kxNeg = (np.where(kx[fnx//2:] > -1/lam))[0]
    kxPos = (np.where(kx[:fnx//2] < +1/lam))[0]
    kyNeg = (np.where(ky[fny//2:] > -1/lam))[0]
    kyPos = (np.where(ky[:fny//2] < +1/lam))[0]
    pdb.set_trace()
    print("Performing inverse transform")
    # pdb.set_trace()
    for it in np.arange(fnt):
        Fr[:,:,it] = np.fft.irfft2(Fr_fft[:,:,it],
                                   s=Fr.shape[0:2])

if fft_images:
    print("Creating images...")
    for it in np.arange(nt):
        image = Fr_fft[:,:,it]
        image = image.real
        image = np.fft.fftshift(image,axes=0)
        image[fnx//2-dcw:fnx//2+dcw,:] = np.float64(1e-6)
        image /= np.nanmax(image)
        image = 10*np.log10(np.square(image))
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
                                'Er_fft-filtered-'+strStep[it]+'.pdf')
        print("Saving",savePath,"...")
        plt.savefig(savePath,bbox_inches='tight',dpi=400)

    for it in np.arange(nt):
        image = Fr[:,:,it]
        vmin = 0
        vmax = np.nanmax(Fr)
        fig = plt.figure()
        ax = fig.gca()
        im = ax.pcolormesh(image.T,
                           cmap='gist_heat',
                           vmin=vmin,vmax=vmax,
                           rasterized=True)
        ax.set_aspect('equal')
        ax.set_position([0.10,0.10,0.70,0.35])
        cbar_ax = fig.add_axes([0.85,0.15,0.02,0.35])
        fig.colorbar(im,cax=cbar_ax).set_label('|E| [V/m]')
        savePath = os.path.join(homePath,basePath,projPath,plotPath,
                                'Er-filtered-'+strStep[it]+'.pdf')
        print("Saving",savePath,"...")
        plt.savefig(savePath,bbox_inches='tight',dpi=400)

if plot_mean_y:
    print("Calculating y-axis mean...")
    image = np.mean(Fr,axis=1)
    print("Creating plot...")
    for it in np.arange(nt):
        fig = plt.figure()
        ax = fig.gca()
        ax.plot(xg,image)
        ax.set_ylim([-6e-3,+6e-3])
        ax.set_xlabel('Zonal [m]')
        ax.set_ylabel('$|E|$ [V/m]')
        savePath = os.path.join(homePath,basePath,projPath,plotPath,
                                'Er_y-mean_multi.pdf')
        print("Saving",savePath,"...")
        plt.savefig(savePath,bbox_inches='tight',dpi=400)

if plot_mean_xy:
    print("Calculating full spatial mean...")
    xy_mean = np.zeros(nt)
    for it in np.arange(nt):
        xy_mean[it] = np.mean(Fr[:,:,it])
    print("Creating plot...")
    yMin = 0
    yMax = +9e-3
    fig = plt.figure()
    ax = fig.gca()
    # pdb.set_trace()
    ax.plot(1e3*params['dt']*params['nout']*timeStep,
            xy_mean,'b.')
    # ax.set_ylim([yMin,yMax])
    ax.set_xlabel('Time [ms]')
    ax.set_ylabel('$|E|$ [V/m]')
    savePath = os.path.join(homePath,basePath,projPath,plotPath,
                            'Er-filtered_xy-mean.pdf')
    print("Saving",savePath,"...")
    plt.savefig(savePath,bbox_inches='tight',dpi=400)


if calc_interp:
    Lx = plane['nx']*plane['dx']
    k = np.linspace(0,np.pi/Lx,plane['nx'])
    theta = np.linspace(0,359,360)
    # xi = k*np.cos(theta*(np.pi/180))
    # yi = k*np.sin(theta*(np.pi/180))
    # gd = griddata(xg,yg,data[:,:,nt-1]
    xx, yy = np.meshgrid(k,theta)
    f = scint.interp2d(xx,yy,data[:,:,nt-1],kind='cubic')
    pdb.set_trace()

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
print("Done")
