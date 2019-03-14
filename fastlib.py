import os
import glob
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys

import pybra.cmd
try:
    import weio
except:
    import sys
    sys.path.append('../../_libs/pyDatView/')
    import weio
    

FAST_EXE='openfast'
bDEBIAN=False


def readFASTOutAscii(filename):
    # Read with panda
    f = weio.FASTOutFile(filename)
    df = f.toDataFrame()
    #df=pd.read_csv(filename, sep='\t', skiprows=[0,1,2,3,4,5,7])
    #df.rename(columns=lambda x: x.strip(),inplace=True)
    return df


def createStepWind(filename,WSstep=1,WSmin=3,WSmax=25,tstep=100,dt=0.5,tmin=0,tmax=999):
    f = weio.FASTWndFile()
    Steps= np.arange(WSmin,WSmax+WSstep,WSstep)
    print(Steps)
    nCol = len(f.colNames)
    nRow = len(Steps)*2
    M = np.zeros((nRow,nCol));
    M[0,0] = tmin
    M[0,1] = WSmin
    for i,s in enumerate(Steps[:-1]):
        M[2*i+1,0] = tmin + (i+1)*tstep-dt 
        M[2*i+2,0] = tmin + (i+1)*tstep
        M[2*i+1,1] = Steps[i]
        if i<len(Steps)-1:
            M[2*i+2,1] = Steps[i+1]
        else:
            M[2*i+2,1] = Steps[-1]
    M[-1,0]= max(tmax, (len(Steps)+1)*tstep)
    M[-1,1]= WSmax
    f.data=pd.DataFrame(data=M,columns=f.colNames)
    #
    print(f.data)
    f.write(filename)
    #plt.plot(M[:,0],M[:,1])
    #plt.show()

    #print(f.toDataFrame())
    #pass
#createStepWind('test.wnd',tstep=200,WSmax=28)
# createStepWind('test.wnd',tstep=200,WSmin=5,WSmax=7,WSstep=2)

def run_fastfiles(fastfiles, exe=None, parallel=True, ShowOutputs=True, nCores=None):
    if exe is None:
        exe=FAST_EXE
    pybra.cmd.run_cmds(fastfiles, exe, parallel=parallel, ShowOutputs=ShowOutputs, nCores=nCores)

def run_fast(input_file, exe=None, wait=True, ShowOutputs=False, debian=None):
    if exe is None:
        exe=FAST_EXE
    return pybra.cmd.run_cmd(input_file, exe=exe, wait=wait, ShowOutputs=ShowOutputs, debian=debian)

def removeFASTOuputs(workdir):
    # Cleaning folder
    for f in glob.glob(os.path.join(workdir,'*.out')):
        os.remove(f)
    for f in glob.glob(os.path.join(workdir,'*.outb')):
        os.remove(f)
    for f in glob.glob(os.path.join(workdir,'*.ech')):
        os.remove(f)
    for f in glob.glob(os.path.join(workdir,'*.sum')):
        os.remove(f)


# def detectFastFiles(workdir):
#     FstFiles=glob.glob(os.path.join(workdir,'*.fst'))+glob.glob(os.path.join(workdir,'*.FST'))
#     DatFiles=glob.glob(os.path.join(workdir,'*.dat'))+glob.glob(os.path.join(workdir,'*.DAT'))
#     Files=dict()
#     Files['Main']      = FstFiles
#     Files['Inflow']    = None
#     Files['Aero']      = None
#     Files['Tower']     = None
#     Files['Blade']     = None
#     Files['AeroBlade'] = None
#     Files['ServoDyn']  = None
#     for f in DatFiles:
#         b = os.path.basename(f).lower()
#         if b.find('inflow'):
#             Files['Inflow'] = f
#     windfile_ref = 'InflowWind.dat';
#     fastfile_ref = 'Turbine.fst';
#     elasfile_ref = 'ElastoDyn.dat';
#         remove
   
