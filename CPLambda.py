import os
import glob
import distutils.dir_util
import pandas as pd
import matplotlib.pyplot as plt
# import compare
import sys
sys.path.append('C:/Work/_libs/pyDatView/weio'); 
sys.path.append('C:/Work/_libs/welib'); 
sys.path.append('C:/Work/_libs/pybra'); 
import welib
import welib.fastlib as fastlib
import weio
import pybra.tictoc
import numpy as np
from shutil import copytree, ignore_patterns, rmtree



def generateParametricInputs(template_dir,workdir=None,main_file=None,OPER=None,RemoveAllowed=False,bStiff=False,bSteadyAero=False,TMax=None):
    if template_dir[-1]=='/'  or template_dir[-1]=='\\' :
        template_dir=template_dir[0:-1]

    if workdir is None:
        workdir=template_dir+'_Parametric'

    # Copying template folder
    print(template_dir, '  ',workdir)
    if os.path.exists(workdir) and RemoveAllowed:
        rmtree(workdir)
    #distutils.dir_util.copy_tree(template_dir, workdir)
    copytree(template_dir, workdir, ignore=ignore_patterns('.git'))
    fastlib.removeFASTOuputs(workdir)

    print('Generating fast input files...')

    # --- Fast main file use as "master"
    if main_file is None:
        FstFiles=set(glob.glob(os.path.join(template_dir,'*.fst'))+glob.glob(os.path.join(template_dir,'*.FST')))
        print(FstFiles)
        if len(FstFiles)>1:
            raise Exception('More than one fst file found in template folder') 
        main_file=FstFiles.pop()
    main_file=os.path.join(workdir, main_file)

    # --- Reading Master File
    print('Reading template file: '+main_file)
    FST = weio.FASTInFile(main_file);
    windfilename_ref = os.path.join(workdir,FST['InflowFile'].strip('"'))
    edfilename_ref   = os.path.join(workdir,FST['EDFile'].strip('"'))
    adfilename_ref   = os.path.join(workdir,FST['AeroFile'].strip('"'))
    sdfilename_ref   = os.path.join(workdir,FST['ServoFile'].strip('"'))
# 
    Wnd = weio.FASTInFile(windfilename_ref);
    ED  = weio.FASTInFile(edfilename_ref);
    AD  = weio.FASTInFile(adfilename_ref);
    windbasename_ref = os.path.basename(windfilename_ref)
    edbasename_ref   = os.path.basename(edfilename_ref)
    adbasename_ref   = os.path.basename(adfilename_ref)
    sdbasename_ref   = os.path.basename(sdfilename_ref)
    fastbasename_ref = os.path.basename(main_file)

    # Rewritting SD file, making sure the controller is off
    SD  = weio.FASTInFile(sdfilename_ref);
    SD['PCMode']=0;
    SD['VSContrl']=0;
    SD['YCMode']=0;
    SD.write (os.path.join(workdir,sdbasename_ref))
# 
    # --- Generating inputs
    fastfiles=[]
    if OPER is None:
        raise Exception('Please provide OPER')
    Keys = list(OPER.keys())
    nValues = len(OPER[Keys[0]])
    print('Number of values ',nValues)
    for i in range(nValues):
        Params = [(k,OPER[k][i]) for k in Keys]
        strID = '_{:03d}'.format(i)
        print(Params)
        for k in Keys:
            val = OPER[k][i]
            if k.lower()=='ws':
                strID += '_ws{:04.1f}'.format(val)
            elif k.lower()=='omega':
                strID += '_om{:04.2f}'.format(val)
            elif k.lower()=='pitch':
                strID += '_pt{:04.2f}'.format(val)
            else:
                raise Exception('Not supported {}'.format(k))
        fastfilename = fastbasename_ref.replace('.fst',strID+'.fst')
        windfilename = windbasename_ref.replace('.dat',strID+'.dat')
        edfilename   = edbasename_ref.replace('.dat',strID+'.dat')
        adfilename   = adbasename_ref.replace('.dat',strID+'.dat')
        sdfilename   = sdbasename_ref.replace('.dat',strID+'.dat')
        for k in Keys:
            val = OPER[k][i]
            if k.lower()=='ws':
                Wnd['WindType']   = 1
                Wnd['HWindSpeed'] = val
                FST['InflowFile'] = '"'+windfilename+'"'
                Wnd.write(os.path.join(workdir,windfilename))
            elif k.lower()=='omega':
                ED['RotSpeed'] = val
            elif k.lower()=='pitch':
                ED['BlPitch(1)'] = val
                ED['BlPitch(2)'] = val
                ED['BlPitch(3)'] = val

        if TMax is not None:
           FST['TMax'] = TMax
        if bSteadyAero:
            AD['AFAeroMod']=1 # remove dynamic effects dynamic
        #
        ED['GenDOF'] = 'False' # important to prescribe rot speed
        if bStiff:
            ED['FlapDOF1']='False'
            ED['FlapDOF2']='False'
            ED['EdgeDOF' ]='False'
            ED['TeetDOF' ]='False'
            ED['DrTrDOF' ]='False'
            ED['YawDOF'  ]='False'
            ED['TwFADOF1']='False'
            ED['TwFADOF2']='False'
            ED['TwSSDOF1']='False'
            ED['TwSSDOF2']='False'
        FST['EDFile'] = '"'+edfilename+'"'
        FST['AeroFile'] = '"'+adfilename+'"'
        #FST['ServoFile'] = '"'+sdfilename+'"'
        ED.write (os.path.join(workdir,edfilename))
        AD.write (os.path.join(workdir,adfilename))
        FST.write(os.path.join(workdir,fastfilename))
        fastfiles.append(os.path.join(workdir,fastfilename))
#     os.remove(os.path.join(workdir,fastbasename_ref))
#     os.remove(edfilename_ref)
    return fastfiles,workdir

def parametricPostPro(workdir,TimeAvgWindow=10,ColKeep=None,ColSort=None):
    OutFiles=glob.glob(os.path.join(workdir,'*.outb'))
    result=None
    if len(OutFiles)<=0:
        raise Exception('No outb files found in {}'.format(workdir))
    for i,f in enumerate(OutFiles):
        print(f)
        df=weio.FASTOutFile(f).toDataFrame()
        if ColKeep is not None:
            df=df[ColKeep]
        ## Start time and end time of window
        iEndTime=df.index[-1]
        endTime=df['Time_[s]'][iEndTime]
        iStartTime=(df['Time_[s]']-(endTime-TimeAvgWindow)).abs().idxmin()
        startTime=df['Time_[s]'][df.index[iStartTime]]
        ## Absolute and relative differences at wind extremities
        DeltaValuesAbs=(df.iloc[iEndTime]-df.iloc[iStartTime]).abs()
        DeltaValuesRel=(df.iloc[iEndTime]-df.iloc[iStartTime]).abs()/df.iloc[iEndTime]
        EndValues=df.iloc[iEndTime]
        ## Stats values during window
        IWindow=df['Time_[s]']>startTime
        MeanValues = pd.DataFrame(df[IWindow].mean()).transpose()
        StdValues  = df[IWindow].std()
        if i==0:
            result = MeanValues.copy()
        else:
            result=result.append(MeanValues, ignore_index=True)

    if ColSort is not None:
        # Sorting 
        result.sort_values([ColSort],inplace=True,ascending=True)
        result.reset_index(drop=True,inplace=True) 

    #result.drop(['Time'],axis=1,inplace=True)
    #result.drop(['Wind1VelY'],axis=1,inplace=True)
    #result.drop(['Wind1VelZ'],axis=1,inplace=True)
    #result.rename(columns=lambda x:'WS'    if x=='Wind1VelX' else x,inplace=True)
    #result.rename(columns=lambda x:'RPM'   if x=='RotSpeed' else x,inplace=True)
    #result.rename(columns=lambda x:'Pitch' if x=='BldPitch1' else x,inplace=True)
    #result.plot(x='Wind1VelX',y=I)
    #print(result)
    #plt.show()
    return result 

# Time      	'Wind1VelX 	'Wind1VelY 	'Wind1VelZ 	'RotSpee 	'BldPitch1 	'GenSpeed  	'NacYaw    	'TwHt3TDxt 	'TwHt3TDyt 	'TwHt3TDzt 	'RootFxc1  	'RootFyc1  	'RootFzc1  	'RootFxb1  	'RootFyb1  	'RootMxc1  	'RootMyc1  	'RootMzc1  	'RootMxb1  	'RootMyb1  	'RotThrust 	'LSShftFya 	'LSShftFza 	'LSShftFys 	'LSShftFzs 	'RotTorq   	'LSSTipMya 	'LSSTipMza 	'LSSTipMys 	'LSSTipMzs 	'HSShftTq  	'HSShftPwr 	'B1Azimuth 	'B1Pitch   	'RtSpeed   	'RtTSR     	'RtAeroCp  	'RtAeroCt  	'RtVAvgxh  	'RtVAvgyh  	'RtVAvgzh  	'RtAeroFxh 	'RtAeroFyh 	'RtAeroFzh 	'RtAeroMxh 	'RtAeroMyh 	'RtAeroMzh 	'RtAeroPwr 	'RtSkew    	'RtArea    	'GenPwr    


# --- Main Parameters
refdir   = '../OF2-SimpleGen/'
workdir  ='Model_Parametric_Focus/'
fastfile = 'SWT-2.3-93OpenFAST2_R2.fst';
TMax=10
ColKeep=['RotSpeed_[rpm]','BldPitch1_[deg]','RtAeroCp_[-]','RtAeroCt_[-]','Wind1VelX_[m/s]','Time_[s]']
#fastlib.FAST_EXE='OpenFAST2_x64d_ebra.exe'
fastlib.FAST_EXE='OpenFAST2_win32d_ebra.exe'




# --- Parametric
fst = weio.FASTInFile(os.path.join(refdir,fastfile))
print(fst['EDFile'])
ed  = weio.FASTInFile(os.path.join(refdir,fst['EDFile'].replace('"','')))
print('R=', ed['TipRad'])

U0=5
R=ed['TipRad']
Lambda = np.linspace(7.5,8.5,11);
Pitch  = np.linspace(-1.0,1.0,22);
Omega  = U0 * Lambda/R*60/(2*np.pi)  # TODO, use more realistic combinations of WS and Omega

OPER=dict()
OPER['WS']    = []
OPER['Pitch'] = []
OPER['Omega'] = []
for p in Pitch:
    for o in Omega:
        OPER['WS'].append(U0)
        OPER['Pitch'].append(p)
        OPER['Omega'].append(o)

print(Lambda)
print(Pitch)
print(Omega)



# --- Generating input files and running them
fastfiles,workdir=generateParametricInputs(template_dir=refdir,workdir=workdir,main_file=fastfile, RemoveAllowed=True, OPER=OPER, bStiff=True, bSteadyAero=True, TMax=TMax)
with open(os.path.join(workdir,'_RUN_ALL.bat'), 'w') as f:
    for l in [fastlib.FAST_EXE + ' '+ os.path.basename(f) for f in fastfiles]:
        f.write("%s\n" % l)
fastlib.run_fastfiles(fastfiles,parallel=True,Outputs=False,nCores=4)



# --- Postpro
result = parametricPostPro(workdir,TimeAvgWindow=5,ColKeep=ColKeep,ColSort='RotSpeed_[rpm]')
print(result)        
result['lambda_[-]'] = result['RotSpeed_[rpm]']*R*2*np.pi/60/result['Wind1VelX_[m/s]']
result.sort_values(['lambda_[-]','BldPitch1_[deg]'],ascending=[True,True],inplace=True)
ColKeep=['lambda_[-]','BldPitch1_[deg]','RtAeroCp_[-]','RtAeroCt_[-]']
result=result[ColKeep]
result.to_csv('CPLambda_Flat_Focus.csv',sep='\t',index=False)


print(result)        
LAMBDA= np.unique(result['lambda_[-]'].values)
PITCH = np.unique(result['BldPitch1_[deg]'].values)
print('>LAMBDA')
print(LAMBDA)
print(Lambda)
print('>PITCH')
print(PITCH)
print(Pitch)
CP = result['RtAeroCp_[-]'].values
MCP=CP.reshape((len(LAMBDA),len(PITCH)))
print(MCP)

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
fig = plt.figure()
ax = fig.gca(projection='3d')

# Make data.
X = LAMBDA
Y = PITCH
X, Y = np.meshgrid(X, Y)
Z = MCP
MCP[MCP<0]=0
print('MAX CP')
print(MCP.argmax())
i,j = np.unravel_index(MCP.argmax(), MCP.shape)
print(i,j)
print('CP:',MCP[i,j],' lambda:',X[j,i],' pitch:',Y[j,i],PITCH[j])

# X = np.arange(-5, 5, 0.25)
# Y = np.arange(-5, 5, 0.25)
# X, Y = np.meshgrid(X, Y)
# R = np.sqrt(X**2 + Y**2)
# Z = np.sin(R)

# Plot the surface.
surf = ax.plot_surface(X, Y, np.transpose(Z), cmap=cm.coolwarm, linewidth=0, antialiased=True,alpha=0.8)
ax.scatter(X[j,i],Y[j,i],MCP[i,j],c='k',marker='o',s=20)

# # Customize the z axis.
# ax.set_zlim(-1.01, 1.01)
# ax.zaxis.set_major_locator(LinearLocator(10))
# ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()


# result1=powerCurvePostPro('../OpenFAST_V27_Power_PowerCurve/')
# result2=powerCurvePostPro('../OpenFAST_V27_Power_PowerCurve/')
# compare.compare(result1,result2)
# print(result1)
# result1.to_csv('OperationalConditions.csv',sep='\t',index=False)

# print('>>>>DONE')







