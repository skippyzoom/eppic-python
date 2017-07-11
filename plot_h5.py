#!/usr/bin/env python
import h5py
import numpy as np
import os, sys
from pathlib import Path
# from matplotlib.pyplot import figure,show
import matplotlib.pyplot as plt
import importlib.machinery

## Set the project directory and import variables
baseDir = '/projectnb/eregion/may/Stampede_runs/'
projDir = 'quasineutral_static_dust/run009/'
# sys.path.append(baseDir+projDir)
# from eppic import *
eppic = importlib.machinery.SourceFileLoader('eppic',baseDir+projDir+'eppic.py').load_module()
path = Path(baseDir+projDir+'parallel')
# nout = 32
timestep = [256,512]
# fileName = ['parallel'] * len(timestep)
# for ts in range(len(timestep)):
#     fileName[ts] += '{:06d}'.format(timestep[ts])+'.h5'
# dataPath = path/fileName[0]

den = []
for ts in range(len(timestep)):
    stepName = 'parallel{:06d}'.format(nout*timestep[ts])
    fileName = stepName+'.h5'
    plotName = stepName+'.pdf'
    print("Reading",fileName)
    with h5py.File(path/fileName,'r') as f:
        den = np.rot90(f['/den1'][:,512-256:512+255],3)
    print("Plotting",plotName)
    plt.pcolormesh(den)
    plt.colorbar()
    plt.savefig(baseDir+projDir+plotName)
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
