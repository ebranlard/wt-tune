import os
import glob
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys

import distutils.dir_util
from shutil import copytree, ignore_patterns, rmtree, copyfile

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


def templateReplace(template_dir, PARAMS, workdir=None, main_file=None, name_function=None, RemoveAllowed=False, RemoveRefSubFiles=False):
    """ Replace parameters in a fast folder using a list of dictionaries where the keys are for instance:
        'FAST|DT', 'EDFile|GBRatio', 'ServoFile|GenEff'
    """
    def fileID(s):
        return s.split('|')[0]
    def basename(s):
        return os.path.splitext(os.path.basename(s))[0]
    def rebase(s,sid):
        split = os.path.splitext(os.path.basename(s))
        return os.path.join(workdir,split[0]+sid+split[1])

    # Default value of workdir if not provided
    if template_dir[-1]=='/'  or template_dir[-1]=='\\' :
        template_dir=template_dir[0:-1]
    if workdir is None:
        workdir=template_dir+'_Parametric'

    # Copying template folder to workdir
    print(template_dir, '  ',workdir)
    if os.path.exists(workdir) and RemoveAllowed:
        rmtree(workdir)
    distutils.dir_util.copy_tree(template_dir, workdir)
    #copytree(template_dir, workdir, ignore=ignore_patterns('.git'))
    removeFASTOuputs(workdir)

    print('Generating fast input files...')
    # --- Fast main file use as "master"
    if main_file is None:
        FstFiles=set(glob.glob(os.path.join(template_dir,'*.fst'))+glob.glob(os.path.join(template_dir,'*.FST')))
        if len(FstFiles)>1:
            print(FstFiles)
            raise Exception('More than one fst file found in template folder, provide `main_file` or ensure there is only one `.fst` file') 
        main_file=rebase(FstFiles.pop(),'')
    else:
        main_file=os.path.join(workdir, os.path.basename(main_file))

    # Params need to be a list
    if not isinstance(PARAMS,list):
        PARAMS=[PARAMS]

    fastfiles=[]
    for p in PARAMS:
        if name_function is None:
            raise NotImplementedError('')
        strID =name_function(p)
        FileTypes = set([fileID(k) for k in list(p.keys())])

        # ---Copying main file and reading it
        fst_full = rebase(main_file,strID)
        copyfile(main_file, fst_full )
        Files=dict()
        Files['FAST']=weio.FASTInFile(fst_full)
        # --- Looping through required files and opening them
        for t in FileTypes: 
            # Doing a naive if
            # The reason is that we want to account for more complex file types in the future
            if t=='FAST':
                continue
            filename   = Files['FAST'][t].strip('"')
            fullname   = rebase(filename,'')
            fullrebase = rebase(filename,strID)
            copyfile(fullname, fullrebase)
            Files['FAST'][t] = '"'+os.path.basename(fullrebase)+'"'
            # Reading files
            Files[t]=weio.FASTInFile(fullrebase)
        # --- Replacing in files
        for k,v in p.items():
            t,kk=k.split('|')
            Files[t][kk]=v
            #print(t+'|'+kk+'=',v)
        # --- Rewritting all files
        for t in FileTypes:
            Files[t].write()

        fastfiles.append(fst_full)
    # --- Remove extra files at the end
    if RemoveRefSubFiles:
        FST = weio.FASTInFile(main_file)
        for t in FileTypes:
            if t=='FAST':
                continue
            filename   = FST[t].strip('"')
            fullname   = rebase(filename,'')
            os.remove(fullname)
    os.remove(main_file)

    return fastfiles

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
  


if __name__=='__main__':
    pass
    # --- Test of templateReplace
    def naming(p):
        return '_ws_'+str(p['InflowFile|HWindSpeed'])
    PARAMS                          = {}
    PARAMS['FAST|TMax']             = 10
    PARAMS['FAST|DT']               = 0.01
    PARAMS['FAST|DT_Out']           = 0.1
    PARAMS['EDFile|RotSpeed']       = 100
    PARAMS['EDFile|BlPitch(1)']     = 1
    PARAMS['EDFile|GBoxEff']        = 0.92
    PARAMS['ServoFile|VS_Rgn2K']    = 0.00038245
    PARAMS['ServoFile|GenEff']      = 0.95
    PARAMS['InflowFile|HWindSpeed'] = 8
    templateReplace(ref_dir,PARAMS,name_function=naming,RemoveRefSubFiles=True)

