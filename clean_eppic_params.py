#!/usr/bin/env python

baseDir = '/projectnb/eregion/may/Stampede_runs/'
projDir = 'quasineutral_static_dust/run009/'
readFile = baseDir+projDir+'eppic.i'
writeFile = baseDir+projDir+'eppic.py'
rf = open(readFile,"r")
wf = open(writeFile,"w")
for line in rf:
    cleaned = line.strip(' \t\n\r')
    if cleaned:
        if cleaned[0] == ';':
            pass
        else:
            wf.write(cleaned+'\n')
rf.close()
wf.close()
