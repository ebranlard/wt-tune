import copy
import distutils.dir_util
import glob
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import random 
import shutil 
from scipy import interpolate
import subprocess
import time
import pdb
# My libs
from pybra import pandalib
from pybra import cmd
import weio
import fastlib


# def prepare_run_folder(template_dir,sim_dir):
#     # Copying template folder
#     distutils.dir_util.copy_tree(template_dir, sim_dir)
#     # Cleaning folder, just in case
#     fastlib.removeFASTOuputs(sim_dir)
# # 
# def prepare_template_folder(ref_dir,workdir,airfoilFileNames,OPER,FAST):
#     def naming(p):
#         return '_{:02.0f}'.format(p['InflowFile|HWindSpeed'])
# 
#     PARAMS=[]
#     for wsp,rpm,pit in zip(OPER['WS'],OPER['RPM'],OPER['Pitch']):
#         p=dict()
#         if wsp<6:
#             p['FAST|TMax']         = 8
#         elif wsp<9:
#             p['FAST|TMax']         = 6
#         else:
#             p['FAST|TMax']         = 4
#         p['FAST|DT']               = 0.01
#         p['FAST|DT_Out']           = 0.1
#         p['FAST|OutFileFmt']       = 1 # TODO
#         p['EDFile|RotSpeed']       = rpm
#         p['EDFile|BlPitch(1)']     = pit
#         p['EDFile|BlPitch(2)']     = pit
#         p['EDFile|BlPitch(3)']     = pit
#         p['EDFile|GBoxEff']        = 94.
#         p['ServoFile|VS_Rgn2K']    = 0.00038245
#         p['ServoFile|GenEff']      = 94.
#         p['InflowFile|HWindSpeed'] = wsp
#         p['InflowFile|WindType']   = 1 # Setting steady wind
#         p['InflowFile|PLexp']      = 0.353 # 0.209
#         PARAMS.append(p)
# 
#     fastlib.templateReplace(ref_dir,PARAMS,workdir=workdir,name_function=naming,RemoveRefSubFiles=True)
#     # --- Rewriting the airfoil files, purely to reduce diff
#     for f in airfoilFileNames:
#         AF=weio.FASTInFile(os.path.join(workdir,f))
#         AF.write();

def read_airfoils(airfoilFileNames,workdir=''):
    airfoils=[]
    for f in airfoilFileNames:
        AF = weio.FASTInFile(os.path.join(workdir,f))
        af=dict()
        af['name']=os.path.splitext(os.path.basename(f))[0]
        af['polar']=AF['AFCoeff']
        Tab=af['polar']
        # Total range of alpha values, and window where we focus
        aRange=[-180,180]
        aWindow=[-40,40]
        # Finding min/max Cl Cd in windows range of alpha values
        iStart = np.argmin(abs(Tab[:,0] - aWindow[0]))
        iEnd   = np.argmin(abs(Tab[:,0] - aWindow[1]))
        af['Clmin']  = np.min(Tab[iStart:iEnd,1])
        af['Clmax']  = np.max(Tab[iStart:iEnd,1])
        af['Cdmin']  = np.min(Tab[iStart:iEnd,2])
        af['iClmin'] = np.argmin(Tab[iStart:iEnd,1])+iStart
        af['iClmax'] = np.argmax(Tab[iStart:iEnd,1])+iStart
        af['iCdmin'] = np.argmin(Tab[iStart:iEnd,2])+iStart
        # For convenience
        aBefore = [aWindow[0]-2,aWindow[0]-1,aWindow[0]]
        aAfter  = [aWindow[1],aWindow[1]+1,aWindow[1]+2]
        af['aClDelta'] = aBefore + [Tab[af['iClmin'],0],Tab[af['iClmax'],0]] + aAfter
        af['aCdDelta'] = aBefore + [Tab[af['iCdmin'],0]                    ] + aAfter
        af['ClDeltaMax'] = np.array([0,0,0,0 ,1,0,0,0])
        af['ClDeltaMin'] = np.array([0,0,0,-1,0,0,0,0])
        af['CdDelta']    = np.array([0,0,0,1   ,0,0,0])
        airfoils.append(af)
    return airfoils

def interp_spline(x0,y0,x1):
    if len(x0)!=len(y0):
        raise Exception('Arguments of interp_spline must have the same length')
    x1=np.array(x1)
    tck = interpolate.splrep(x0, y0)
    y1  = interpolate.splev(x1, tck)
    y1 [ x1< np.min(x0) ] = y0[0]
    y1 [ x1> np.max(x0) ] = y0[-1]
    return y1

def plot_airfoils(airfoils,fig=None,names=None,**kwargs):
    naf=len(airfoils)
    #naf=1
    if not fig:
        fig=plt.figure()
    if len(fig.axes)<naf:
        ax=[]
        for i in range(naf):
            ax.append(fig.add_subplot(2,naf,i+1))
            ax.append(fig.add_subplot(2,naf,naf+i+1))
    else:
        ax=fig.axes
    for i in range(naf):
        af=airfoils[i]
        Tab=af['polar']
        # Hack
        #ClDelta =  af['ClDeltaMax']*-0.1 + af['ClDeltaMin']*0.05
        #CdDelta =  af['CdDelta']*0.02
        #DeltaCl = interp_spline(af['aClDelta'], ClDelta, Tab[:,0])
        #DeltaCd = interp_spline(af['aCdDelta'], CdDelta, Tab[:,0])
        #ClMut =  Tab[:,1]+DeltaCl
        #CdMut =  Tab[:,2]+DeltaCd
        # Cl
        alpha=Tab[:,0]
        ax[2*i].plot(alpha, Tab[:,1],**kwargs)
        ax[2*i].plot(alpha[af['iClmin']], Tab[af['iClmin'],1],'o')
        ax[2*i].plot(alpha[af['iClmax']], Tab[af['iClmax'],1],'o')
        #ax[2*i].plot(af['aClDelta'], ClDelta,'--')
        #ax[2*i].plot(alpha, DeltaCl,':')
        #ax[2*i].plot(alpha, ClMut,'-')
        if names is not None:
            ax[2*i].set_title(names[i])
        ax[2*i].set_xlim([-50,50])
        # Cd
        ax[2*i+1].plot(alpha, Tab[:,2],**kwargs)
        ax[2*i+1].plot(alpha[af['iCdmin']], Tab[af['iCdmin'],2],'o')
        #ax[2*i+1].plot( af['aCdDelta'], CdDelta,'--')
        #ax[2*i+1].plot(alpha, DeltaCd,':')
        #ax[2*i+1].plot(alpha, CdMut,'-')
        ax[2*i+1].set_xlim([-50,50])
    return fig

def patch_airfoil(gene,af):
    # Genes are assumed to be between 0 and 1, some are converted them to -1 and 1
    pheno=np.zeros(np.size(gene))
    pheno[0]=0.3*(gene[0]*2-1)
    pheno[1]=0.3*(gene[1]*2-1)
    pheno[2]=0.07*(gene[2]) # for drag only increase!
    ClDelta =  af['ClDeltaMax']*pheno[0] + af['ClDeltaMin']*pheno[1]
    CdDelta =  af['CdDelta']*pheno[2]
    DeltaCl = interp_spline(af['aClDelta'], ClDelta, af['polar'][:,0])
    DeltaCd = interp_spline(af['aCdDelta'], CdDelta, af['polar'][:,0])
    af['polar'][:,1] =  af['polar'][:,1]+DeltaCl
    af['polar'][:,2] =  af['polar'][:,2]+DeltaCd
    return af

# def patch_airfoils(chromosome,airfoils_ref):
#     new_airfoils=copy.deepcopy(airfoils_ref)
#     genes = np.array_split(chromosome,len(airfoils_ref))
#     for [af,gn] in zip(new_airfoils,genes):
#         af = patch_airfoil(gn,af)
#     return new_airfoils
# 
# def rewrite_airfoils(airfoils,airfoilFileNames,workdir=''):
#     # Reading the files, chainging the polars and rewriting
#     for af,f in zip(airfoils,airfoilFileNames):
#         #print('Rewriting {}'.format(os.path.join(workdir,f)))
#         AF = weio.FASTInFile(os.path.join(workdir,f))
#         AF['AFCoeff'] = af['polar']
#         AF.write()
# 
# def patch_airfoil_files(chromosome,airfoils_ref,airfoilFileNames,workdir=''):
#     ### --- Create hacked airfoil data
#     airfoils_mut=patch_airfoils(chromosome,airfoils_ref)
#     rewrite_airfoils(airfoils_mut,airfoilFileNames,workdir=workdir)
#     return airfoils_mut

def run_sim(sim_dir,FAST,nSIM,exe=None):
    if FAST==1:
        files=glob.glob(os.path.join(sim_dir,'*.fst'))
        if len(files)==0:
            raise Exception('No fast file found in sim folder: '+sim_dir)
        if len(files)!=nSIM:
            print(files)
            raise Exception('The number of fast files ({}) is not equal to nSIM ({}) in: {}'.format(len(files),nSIM,sim_dir))
        cmd.run_cmds(files, exe, parallel=True, ShowOutputs=False)
    else:
        files=glob.glob(os.path.join(sim_dir,'*.dvr'))
        if len(files)==0:
            raise Exception('No driver file found in sim folder:'+sim_dir)
        #if len(files)!=1:
        #    raise Exception('More than one driver file found in sim folder:'+sim_dir)
        if len(files)!=nSIM:
            raise Exception('Number of driver file ({}) different from number of ws ({}): {}'.format(len(files),nSIM,sim_dir))

        bOK=False
        nTry=0
        while not bOK and nTry<10:
            cmd.run_cmds(files, exe, parallel=True, ShowOutputs=False)
            #for f in files:
            #    run_aerodyn(f,exe)
            #run_aerodyn(files[0])
            outfiles=glob.glob(os.path.join(sim_dir,'*.out'))
            bOK=len(outfiles)==nSIM
            nTry=nTry+1
            if not bOK:
               print('>>>>>>>>>>>>>>>>>>>>>>>>>>>> TRY {}  -  CRASH {}/{}'.format(nTry,len(outfiles),nSIM))
               for f in outfiles:
                   os.remove(f)
               time.sleep(1.00)

def postpro_simdir(sim_dir,TimeAvgWindow,FAST):
    files=glob.glob(os.path.join(sim_dir,'*.out'))
    col_names =  ['WS', 'RPM','Pitch','Pgen', 'Qgen', 'Thrust','FlapM','CPAero','CTAero','Converged']
    #['Time', 'B1Azimuth', 'B1Pitch', 'B1N1AxInd' 'B1N1TnInd', 'B1N1Alpha' 'B1N1Cl','B1N1Cd'
    #'RtSpeed', 'RtTSR', 'RtAeroCp', 'RtAeroCq', 'RtAeroCt', 'RtVAvgxh', 'RtVAvgyh', 'RtVAvgzh', 'RtAeroFxh', 'RtAeroFyh', 'RtAeroFzh', 'RtAeroMxh', 'RtAeroMyh', 'RtAeroMzh', 'RtAeroPwr', 'RtSkew', 'RtArea'],

    #
    #fast_file=glob.glob(os.path.join(sim_dir,'*.fst'))
    #Fst = weio.FASTInFile(fast_file[0]);
    #pdb.set_trace()

    perf=pd.DataFrame(columns = col_names)

    
    for outfilename in files:
        #print(outfilename)
        # Read with panda
        df=pd.read_csv(outfilename, sep='\t', skiprows=[0,1,2,3,4,5,7])
        df.rename(columns=lambda x: x.strip(),inplace=True)
        
        #
        # Start time and end time of window
        iEndTime=df.index[-1]
        endTime=df['Time'][iEndTime]
        if TimeAvgWindow is None:
            Omega=df['RotSpeed'].mean()/60*2*np.pi
            TimeAvgWindow=(2*np.pi/Omega)*2 # averaging about two rotations
            if TimeAvgWindow>endTime:
                print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>Time too short!',TimeAvgWindow,endTime,sim_dir)
            TimeAvgWindow=min(TimeAvgWindow,endTime)
        iStartTime=(df['Time']-(endTime-TimeAvgWindow)).abs().idxmin()
        startTime=df['Time'][df.index[iStartTime]]
        # Absolute and relative differences at wind extremities
        DeltaValuesAbs=(df.iloc[iEndTime]-df.iloc[iStartTime]).abs()
        DeltaValuesRel=(df.iloc[iEndTime]-df.iloc[iStartTime]).abs()/df.iloc[iEndTime]
        EndValues=df.iloc[iEndTime]
        # Stats values during window
        IWindow=df['Time']>startTime
        MeanValues = df[IWindow].mean()
        StdValues  = df[IWindow].std()
        #print(DeltaValuesRel[['GenPwr','RotSpeed']])
        sOK=''
        bOK=True
        if FAST==1 and DeltaValuesRel['RotSpeed']*100>0.5:
            bOK=False
            #print('                                                                    >>> NOT OK')
            sOK='NOT OK'
        #if FAST==1:
        #    print('    WS: {:.1f} - Pitch {:.2f} - RPM {:.1f} - DeltaRPM {:.4f}% - DeltaP {:.4f}%  {}'.format(MeanValues['Wind1VelX'],MeanValues['BldPitch1'],EndValues['RotSpeed'],DeltaValuesRel['RotSpeed']*100,DeltaValuesRel['GenPwr']*100,sOK))    
        #else:
        #    print('    WS: {:.1f} - Pitch {:.2f} - RPM {:.1f}'.format(MeanValues['RtVAvgxh'],MeanValues['B1Pitch'],EndValues['RtSpeed']))

        # --- Plotting
        #wsp=MeanValues['Wind1VelX']
        #dfScaled=df/df.max()
        #fig,ax = plt.subplots();
        #PlotSignals=['RotSpeed','GenPwr','BldPitch1']
        #PlotSignals=['RtSpeed','RtAeroPwr','B1Pitch']
        #dfScaled.plot(x='Time',y=PlotSignals,title='WS={:.1f} - {}'.format(wsp,sOK),ax=ax)
        #ax.legend(['{} - {:.1f}'.format(x,y) for x,y in zip(df[PlotSignals].columns.values,df[PlotSignals].max())])

        # --- Storing
        if FAST==1:
            perf.loc[len(perf)]        = np.nan
            perf.iloc[-1]['WS']        = MeanValues['Wind1VelX']
            perf.iloc[-1]['Pgen']      = MeanValues['GenPwr']
            perf.iloc[-1]['Qgen']      = MeanValues['GenTq']
            perf.iloc[-1]['RPM']       = MeanValues['RotSpeed']
            perf.iloc[-1]['Pitch']     = MeanValues['BldPitch1']
            perf.iloc[-1]['CPAero']    = MeanValues['RtAeroCp']
            perf.iloc[-1]['CTAero']    = MeanValues['RtAeroCt']
            perf.iloc[-1]['FlapM']     = MeanValues['RootMyc1']*1000 # Nm
            perf.iloc[-1]['Converged'] = bOK
        else:
            perf.loc[len(perf)]       = np.nan
            perf.iloc[-1]['WS']       = EndValues['RtVAvgxh']
            perf.iloc[-1]['Pgen' ]    = EndValues['RtAeroPwr'] #<<<<<<<<<<
            perf.iloc[-1]['RPM']      = EndValues['RtSpeed']
            perf.iloc[-1]['Pitch']    = EndValues['B1Pitch']
            perf.iloc[-1]['CPAero']   = EndValues['RtAeroCp']
            perf.iloc[-1]['CTAero']   = EndValues['RtAeroCt']
            perf.iloc[-1]['FlapM'] = np.nan
            perf.iloc[-1]['Converged'] = bOK

    perf.sort_values(['WS'],inplace=True)
    perf.reset_index(drop=True,inplace=True) # very important
    return perf


def get_perf_error(perf_sim,ISel,perf_ref=None,perf_all=None):
    if perf_ref is None:
        # Interpolating to get the reference values where the simulation was done
        ws_sim=perf_sim['WS'].values.astype(float)
        perf_ref=pandalib.pd_interp1('WS',perf_all,ws_sim)
    # Computing relative error
    #print(perf_sim[ISel])
    #print(perf_ref[ISel])
    RelError=abs(perf_sim[ISel]-perf_ref[ISel])/perf_ref[ISel].max()
    AbsError=abs(perf_sim[ISel]-perf_ref[ISel])
    return RelError.mean()*100,perf_ref




