#!/usr/bin/env python
# This file will contain I/O function for processing EPPIC data.

def read_parameters(name='eppic.i',path='./',comment=';'):
    """Read a parameter input file from an EPPIC run.

    This function checks for a file named 'name' in 
    <path> and reads in all lines that do not begin
    with 'comment'.
    """

    import os
    from pathlib import Path
    import pdb

    #-->Check if path is a string or PosixPath
    rfn = os.path.join(str(path),name)
    # rfn = path / Path(name)
    # pdb.set_trace()
    try:
        rf = open(rfn,'r')
    except OSError:
        print("Cannot open",rfn,"for reading.")
        return None
    params = dict()
    for line in rf:
        okay = line.strip(' \t\n\r')
        if okay:
            if okay[0] == comment:
                pass
            else:
                name,value = okay.split('=')
                value = value.lstrip()
                if ' ' in value:
                    pass
                elif '.' in value:
                    value = float(value)
                elif 'e' in value:
                    value = float(value)
                else:
                    value = int(value)
                params[name.strip()] = value
    return params
    rf.close()

def clean_params(path='./'):
    """Clean an EPPIC input file for python use.

    This function checks for a file called 'eppic.i' in <path>.
    If the file exists, this function extracts the variables
    to a new file called 'eppic.py', which python scripts may
    import in order to obtain simulation variable values.
    """

    rfn = path+'/eppic.i'
    wfn = path+'/eppic.py'
    try:
        rf = open(rfn,'r')
    except OSError:
        print("Cannot open",rfn,"for reading.")
        return None
    try:
        wf = open(wfn,'w') 
    except OSError:
       print("Cannot open",wfn,"for writing.")
    for line in rf:
        cleaned = line.strip(' \t\n\r')
        if cleaned:
            if cleaned[0] == ';':
                pass
            else:
                wf.write(cleaned+'\n')
    rf.close()
    wf.close()
    print("Created",wfn)

def calc_timesteps(path='./'):
    """Calculate the actual number of time steps in a simulation run.

    This function checks for the existence of typical informational
    files from an EPPIC run and attempts to calculate the maximum
    available number of time steps for that run. The calculated value
    is different from params.nt in the case that the run timed out.
    """

    import os, glob

    path = str(path)
    fn_list = [os.path.join(path,'moments1.out'),
               os.path.join(path,'moments0.out'),
               os.path.join(path,'domain000','moments1.out'),
               os.path.join(path,'domain000','moments0.out')]
    for i,fn in enumerate(fn_list):
        if os.path.exists(fn):
            with open(fn) as f:
                print("Calculating time steps from",fn)
                return sum(1 for _ in f)-1
    if os.path.exists(path+'/parallel/'):
        print("Calculating time steps from *.h5 files")
        return len(glob.glob(path+'/parallel/*.h5'))

def set_paths():
    """Return machine-specific paths for simulation data analysis.

    This function checks the machine hostname and compares it to a list
    of known names. If it finds a match, it returns the predefined path
    to the project directory as the first value and the relative path to
    data as the second value.
    """

    import os

    hostname = os.uname()[1]
    if hostname.find('scc') >= 0:
        return ['/projectnb/eregion/may/Stampede_runs/', 
                './']
    elif hostname.find('stampede') >= 0:
        return ['/scratch/02994/may/', 
                './']
    elif hostname.find('Matthews') >= 0:
        return ['/Users/matthewyoung/Documents/BU/research/Projects/',
                'data/eppic/']
    elif hostname.find('fluid') >= 0:
        return ['/home/may/Research/',
                'data/eppic/']
    else:
        return ['./',
                './']

def imgplane(dataName,params,
             dataPath='./',
             ranges=[[0,1],[0,1]],timeStep=0,rot_k=0):
    """Return a (2+1)-D array of simulation data for imaging

    This function reads HDF data from an EPPIC run and extracts the
    requested 2-D plane at each time step, and it manipulates it 
    according to given parameters.
    """

    import os
    import h5py
    import numpy as np

    ##==Set up image plane
    if rot_k % 2 == 0:
        plane = {'nx':
                 params['nx']*params['nsubdomains']//params['nout_avg'],
                 'ny':
                 params['ny']//params['nout_avg'],
                 'dx':
                 params['dx'],
                 'dy':
                 params['dy']}
    else:
        plane = {'nx':
                 params['ny']//params['nout_avg'],
                 'ny':
                 params['nx']*params['nsubdomains']//params['nout_avg'],
                 'dx':
                 params['dy'],
                 'dy':
                 params['dx']}

    ##==Set up the data array
    nts = len(timeStep)
    strStep = []
    fdata = np.zeros((params['nx']*params['nsubdomains'],
                      params['ny'],nts))

    ##==Read data at each time step
    for it,ts in enumerate(timeStep):
        strStep.append('{:06d}'.format(params['nout']*ts))
        fileName = 'parallel'+strStep[it]+'.h5'
        # dataFile = os.path.join(homePath,basePath,projPath,dataPath,
        #                         'parallel',fileName)
        dataFile = os.path.join(dataPath,'parallel',fileName)
        print("Reading",dataName,"from",fileName,"...")
        with h5py.File(dataFile,'r') as f:
            fdata[:,:,it] = np.array(f['/'+dataName])

    ##==Manipulate data
    fdata = np.rot90(fdata,k=rot_k)
    x0 = int(plane['nx']*ranges[0][0])
    xf = int(plane['nx']*ranges[0][1])
    y0 = int(plane['ny']*ranges[1][0])
    yf = int(plane['ny']*ranges[1][1])
    fdata = fdata[x0:xf,y0:yf]

    ##==Return
    return fdata,plane
