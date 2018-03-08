#!/usr/bin/env python
import os
import h5py
import numpy as np
# from pathlib import Path
# from scipy.interpolate import interp2d
# from matplotlib.pyplot import figure, show
import matplotlib.pyplot as plt
import pdb

## It appears that pathlib.Path doesn't work with h5py.File()
# basePath = Path(os.path.expanduser('~')+
#                '/Documents/BU/research/Projects')
# dataPath = 'parametric_wave/run005/data/eppic/parallel'
# fileName = 'parallel004992.h5'
# dataFile = basePath/dataPath/fileName

homePath = os.path.expanduser('~')
basePath = 'Documents/BU/research/Projects'
dataPath = 'parametric_wave/run005/data/eppic/parallel'
fileName = 'parallel004992.h5'
dataFile = os.path.join(homePath,basePath,dataPath,fileName)

print("Reading data...")
with h5py.File(dataFile,'r') as f:
    den = f['/den1'][:]
print("Done")

print("Creating image...")
# fg = figure()
# ax = fg.gca()
# p = ax.pcolormesh(den)
# ax.set_xlabel('x [m]')
# ax.set_ylabel('y [m]')
# ax.set_title('Relative Perturbed Density')
# fg.colorbar(p).set_label('$\delta n/n_0$')
# fg.savefig('den_test.pdf')

# fg, ax = plt.subplots(nrows=1)
# pc = ax.pcolormesh(den)
# plt.show()

plt.pcolormesh(den, rasterized=True)
# plt.colorbar()
print("Done")

print("Saving image...")
plt.savefig('den_test.pdf')
print("Done")
