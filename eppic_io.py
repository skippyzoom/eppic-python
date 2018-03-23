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

def read_ph5(dataName,
             ext='.h5',
             timeStep=0,
             axes='xy',
             dataType='float',
             dataIsFT=0,
             dataPath='./',
             infoPath='./'):

    """Create a (2+1)-D data set from EPPIC HDF data.
    
    This function reads HDF data time steps from an EPPIC run and
    extracts a 2-D plane at each time step, the returns the resultant
    (2+1)-D array.
    """

    import os
    import numpy as np
    import h5py
    import pdb

    ##==Ensure string path for h5py
    dataPath = str(dataPath)

    ##==Read in run parameters
    params = read_parameters(path=infoPath)

    ##==Set up the data array
    nts = len(timeStep)
    strStep = []
    if params['ndim_space'] == 2:
        fdata = np.zeros((params['nx']*params['nsubdomains'],params['ny'],
                          nts),dtype=dataType)
    elif params['ndim_space'] == 3:
        if axes == 'xy' or axes == 'yx':
            fdata = np.zeros((params['nx']*params['nsubdomains'],params['ny'],
                              nts),dtype=dataType)
        elif axes == 'xz' or axes == 'zx':
            fdata = np.zeros((params['nx'],params['nz'],nts),dtype=dataType)
        elif axes == 'yz' or axes == 'zy':
            fdata = np.zeros((params['ny'],params['nz'],nts),dtype=dataType)
        elif axes == 'all':
            fdata = np.zeros((params['nx']*params['nsubdomains'],params['ny'],
                              params['nz'],nts),dtype=dataType)
        else:
            raise IOError('Error in ndim_space')

    ##==Read data at each time step
    for it,ts in enumerate(timeStep):
        strStep.append('{:06d}'.format(params['nout']*ts))
        fileName = 'parallel'+strStep[it]+'.h5'
        dataFile = os.path.join(dataPath,fileName)
        print("Reading",dataName,"from",fileName,"...")
        with h5py.File(dataFile,'r') as f:
            tmp = np.array(f['/'+dataName])
            if tmp.ndim == 2:
                fdata[:,:,it] = tmp
            elif tmp.ndim == 3:
                if axes == 'xy' or axes == 'yx':
                    fdata[:,:,it] = tmp[:,:,0]
                elif axes == 'xz' or axes == 'zx':
                    fdata[:,:,it] = tmp[:,0,:]
                elif axes == 'yz' or axes == 'zy':
                    fdata[:,:,it] = tmp[0,:,:]
                elif axes == 'all':
                    fdata[:,:,:,it] = tmp
            else:
                raise IOError('Error in data dimensions')

    ##==Return data
    return fdata

def imgplane(dataName,
             axes='xy',
             timeStep=0,
             rotate=0,
             fft_direction=0,
             fft_n=[0,0,0],
             calc_gradient=0,
             dataPath='./',
             infoPath='./'):

    """Return a (2+1)-D array of simulation data for imaging

    This function reads HDF data from an EPPIC run and extracts the
    requested 2-D plane at each time step, and it manipulates it 
    according to given parameters.
    """

    import os
    import h5py
    import numpy as np
    import pdb

    ##==Read data at each time step
    fdata = read_ph5(dataName,
                     timeStep=timeStep,
                     # axes=plane['axes'],
                     axes=axes,
                     dataType='float',
                     dataPath=dataPath,
                     infoPath=infoPath)

    ##==Check data dimensions
    try:
        fdata.ndim == 3
    except:
        IOError('Expected fdata to be (2+1)-D')

    ##==Establish full dimensions
    params = read_parameters(path=infoPath)
    nx = fdata.shape[0]
    ny = fdata.shape[1]
    nt = fdata.shape[2]
    dx = params['dx']
    dy = params['dy']

    ##==Create arrays of x- and y-axis data points
    xdata = dx*np.arange(0,nx)
    ydata = dy*np.arange(0,ny)

    ##==Rotate data, if requested
    fdata = np.rot90(fdata,k=rotate)
    if rotate % 2 == 1: xdata,ydata = ydata,xdata

    ##==Set up output container
    data_out = {'x':xdata, 'y':ydata}

    ##==Calculate gradient, if requested
    if calc_gradient:
        Fx = np.zeros_like(fdata)
        Fy = np.zeros_like(fdata)
        for it in np.arange(nt):
            tmp = np.gradient(fdata[:,:,it])
            Fx[:,:,it] = tmp[0]
            Fy[:,:,it] = tmp[1]

    ##==Calculate FFT, if requested
    if fft_direction != 0: 
        if fft_n == [0,0,0]: fft_n = [nx,ny,nt]
        if calc_gradient:
            Fx.astype(complex)
            Fy.astype(complex)
            if fft_direction < 0:
                for it in np.arange(nt):
                    Fx[:,:,it] = np.fft.fft2(Fx[:,:,it],s=fft_n)
                    Fy[:,:,it] = np.fft.fft2(Fy[:,:,it],s=fft_n)
            elif fft_direction > 0:
                for it in np.arange(nt):
                    Fx[:,:,it] = np.fft.ifft2(Fx[:,:,it],s=fft_n)
                    Fy[:,:,it] = np.fft.ifft2(Fy[:,:,it],s=fft_n)
        else:
            fdata.astype(complex)
            if fft_direction < 0:
                for it in np.arange(nt):
                    fdata[:,:,it] = np.fft.fft2(fdata[:,:,it],s=fft_n)
            elif fft_direction > 0:
                for it in np.arange(nt):
                    fdata[:,:,it] = np.fft.ifft2(fdata[:,:,it],s=fft_n)

    ##==Return
    if calc_gradient: data_out['f'] = {'x':Fx, 'y':Fy}
    else: data_out['f'] = fdata
    return data_out
