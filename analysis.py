#!/usr/bin/env python
import os
import h5py
import numpy as np
import matplotlib.pyplot as plt
import pdb
import eppic_io

## Declare project path
projPath = 'parametric_wave/run005/'

## Declare data name and directory
dataName = 'den1'
dataPath = 'data/eppic/'

## Declare plotting preferences
cmap = 'plasma'
plotPath = 'python_images'
plotType = 'pdf'
plotName = dataName+'.'+plotType

## Choose time steps to plot
ntMax = eppic_io.calc_timesteps(path)
timeStep = [1,ntMax-1]

## Set up standard path info
homePath = os.path.expanduser('~')
basePath = 'Documents/BU/research/Projects'
# fileName = 'parallel004992.h5'
dataFile = os.path.join(homePath,basePath,projPath,dataPath,
                        'parallel',fileName)

print("Reading data...")
with h5py.File(dataFile,'r') as f:
    den = f['/den1'][:]
print("Done")

print("Creating image...")

plt.pcolormesh(den, rasterized=True)
print("Done")

print("Saving image...")
savePath = os.path.join(homePath,basePath,projPath,plotPath,
                        plotPath,plotName)
plt.savefig(savePath)
print("Done")
