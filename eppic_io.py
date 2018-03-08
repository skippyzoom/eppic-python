#!/usr/bin/env python
# This file will contain I/O function for processing EPPIC data.
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
    wf = open(wfn,'w')
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
    """Calculate the actual number of time steps in a simulation run."""

    import os, glob

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

# def set_default(default):
#     """Set a variable to a default value if the variable is undefined."""
