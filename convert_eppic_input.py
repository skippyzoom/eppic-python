#!/usr/bin/env python
# Convert eppic.i files to eppic.py files within a directory tree.
# TO DO:
#    Consider only converting an eppic.i file if no local eppic.py exists.
import fnmatch, os, eppic_io

baseDir = '/projectnb/eregion/may/Stampede_runs/'
projDir = 'parametric_wave'
path = os.path.join(baseDir,projDir)
targets = []
for root, dirs, files in os.walk(path):
    for fn in files:
        if fnmatch.fnmatch(fn,"eppic.i"):
            eppic_io.clean_params(root)
