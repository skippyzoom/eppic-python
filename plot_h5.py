#!/usr/bin/env python
import h5py
import numpy as np
import os, sys
from pathlib import Path
import matplotlib.pyplot as plt
import importlib.machinery
import eppic_io

## Declare quantity to plot
dataName = 'den1'
plotType = 'png'

## Set the project directory and import variables
baseDir = '/projectnb/eregion/may/Stampede_runs/'
projDir = 'quasineutral_static_dust/run009/'
path = baseDir+projDir
posixPath = Path(path)
# sys.path.append(baseDir+projDir)
# from eppic import *
eppic = importlib.machinery.SourceFileLoader('eppic',path+'eppic.py').load_module()

## Choose time steps to plot
ntMax = eppic_io.calc_timesteps(path)
timeStep = [1,ntMax-1]

## Load data
# fileName = ['parallel'] * len(timeStep)
# for ts in range(len(timeStep)):
#     fileName[ts] += '{:06d}'.format(timeStep[ts])+'.h5'
# dataPath = path/fileName[0]

# den = []
for ts in range(len(timeStep)):
    # stepName = 'parallel{:06d}'.format(eppic.nout*timeStep[ts])
    # fileName = stepName+'.h5'
    # plotName = stepName+'.pdf'
    strStep = '{:06d}'.format(eppic.nout*timeStep[ts])
    fileName = 'parallel'+strStep+'.h5'
    plotName = dataName+'_'+strStep+'.'+plotType
    print("Reading",fileName)
    with h5py.File(path+'parallel/'+fileName,'r') as f:
        data = np.rot90(f['/'+dataName][:,512-256:512+256],3)
    print(data.shape)
    print("Plotting",plotName)
    plt.pcolormesh(data)
    plt.colorbar()
    plt.savefig(path+plotName)
    plt.close()

# den = np.array(den)
# print(den.shape)
# print(den[0,:,:].shape)
# # # Only plots the last time step
# # plt.pcolormesh(den)
# # plt.colorbar()
# # # plt.show()
# # pdfName = 'den.pdf'
# # print("Saving",pdfName)
# # plt.savefig(baseDir+projDir+pdfName)
# # plt.close()
