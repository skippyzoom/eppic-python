#!/usr/bin/env python
import h5py
import numpy as np
import os, sys
from pathlib import Path
import matplotlib.pyplot as plt
import importlib.machinery
import eppic_io

## Declare quantity to plot
dataName = 'phi'
colorMap = 'seismic'
plotType = 'png'

## Set the project directory and import variables
baseDir = '/projectnb/eregion/may/Stampede_runs/'
projDir = 'quasineutral_static_dust/run009/'
path = baseDir+projDir
posixPath = Path(path)
eppic = importlib.machinery.SourceFileLoader('eppic',path+'eppic.py').load_module()

## Declare data plot ranges
## -->Eventually want to import these from a parameter file
x0 = 0
xf = eppic.nx*eppic.nsubdomains/eppic.nout_avg
y0 = (eppic.ny/4)/eppic.nout_avg
yf = (3*eppic.ny/4)/eppic.nout_avg

## Choose time steps to plot
ntMax = eppic_io.calc_timesteps(path)
timeStep = [1,ntMax-1]

## Loop over time steps and plot
for ts in range(len(timeStep)):
    strStep = '{:06d}'.format(eppic.nout*timeStep[ts])
    fileName = 'parallel'+strStep+'.h5'
    plotName = dataName+'_'+strStep+'.'+plotType
    print("Reading",fileName)
    with h5py.File(path+'parallel/'+fileName,'r') as f:
        data = np.rot90(f['/'+dataName][int(x0):int(xf),int(y0):int(yf)],3)
    print("Plotting",plotName)
    plt.pcolormesh(data,cmap=colorMap)
    plt.colorbar()
    plt.savefig(path+plotName)
    plt.close()
