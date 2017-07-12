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

# ## Loop over time steps and plot
# cmap = plt.get_cmap(colorMap)
# vmin = -0.06
# vmax = 0.06
# for ts in range(len(timeStep)):
#     strStep = '{:06d}'.format(eppic.nout*timeStep[ts])
#     fileName = 'parallel'+strStep+'.h5'
#     plotName = dataName+'_'+strStep+'.'+plotType
#     print("Reading",fileName)
#     with h5py.File(path+'parallel/'+fileName,'r') as f:
#         data = np.rot90(f['/'+dataName][int(x0):int(xf),int(y0):int(yf)],3)
#     print("Plotting",plotName)
#     plt.pcolormesh(data,cmap=cmap,vmin=vmin,vmax=vmax)
#     plt.colorbar()
#     plt.savefig(path+plotName)
#     plt.close()

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
        # print(temp.shape)
        # temp = np.rot90(temp,3)
        # print(temp.shape)
        # data[:,:,:,it] = temp.reshape(nX,nY,nZ)
        data[:,:,it] = temp

xg = np.linspace(x0,xf,data[x0:xf,y0:yf,0].shape[0])
yg = np.linspace(y0,yf,data[x0:xf,y0:yf,0].shape[1])

# plt.pcolormesh(np.transpose(data[x0:xf,y0:yf,1]))
# plt.xlim(x0, xf)
# plt.ylim(y0, yf)
# plt.gca().set_aspect('equal', adjustable='box')
# plt.show()

# f, axarr = plt.subplots(2)
# axarr[0].pcolormesh(data[x0:xf,y0:yf,0])
# axarr[0].set_title(strStep[0])
# axarr[1].pcolormesh(data[x0:xf,y0:yf,1])
# axarr[1].set_title(strStep[1])
# # f.subplots_adjust(hspace=0.3)
# plt.show()

# plt.subplot(2,1,1)
# plt.pcolormesh(xg,yg,data[x0:xf,y0:yf,0])
# plt.axis([x0,xf,y0,yf])
# plt.subplot(2,1,2)
# plt.pcolormesh(xg,yg,data[x0:xf,y0:yf,1])
# plt.axis([x0,xf,y0,yf])
# plt.show()

dataMaxAbs = np.nanmax(np.absolute(data))
# vmin = data.min()
# vmax = data.max()
vmin = -dataMaxAbs
vmax = dataMaxAbs
print(vmin,vmax)
cmap = plt.get_cmap(colorMap)
fig, axes = plt.subplots(2, sharex=True)
img = axes[0].pcolormesh(xg,yg,data[x0:xf,y0:yf,0],
                         cmap=cmap,vmin=vmin,vmax=vmax)
axes[0].set_title(strStep[0])
img = axes[1].pcolormesh(xg,yg,data[x0:xf,y0:yf,1],
                         cmap=cmap,vmin=vmin,vmax=vmax)
axes[1].set_title(strStep[1])
fig.suptitle('Density')
# fig.text(0.50, 0.95, 'Density', transform=fig.transFigure, horizontalalignment='center')
fig.subplots_adjust(right=0.80)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.70])
fig.colorbar(img,cax=cbar_ax)
plt.show()

# ## The following requires from matplotlib.pyplot import figure, show
# xg = np.linspace(x0,xf,data[:,:,0].shape[0])
# yg = np.linspace(y0,yf,data[:,:,0].shape[1])
# fg = figure()
# ax = fg.gca()
# pc = ax.pcolormesh(xg,yg,data[:,:,1])
# ax.set_xlabel('x [m]')
# ax.set_ylabel('z [m]')
# ax.set_title('Density')
# fg.colorbar(pc).set_label('$\delta n/n_0$')


print("\n")
