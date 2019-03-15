import glob
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import numpy.matlib
import os
import random 
import time
import re
import math
import pdb
import distutils.dir_util
import platform

from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import minimize

from deap import base
from deap import creator
from deap import tools

# Mylibraries
from pybra import figlib
from pybra import galib
import weio
import fastlib
from AirfoilTuningTools import *

def genotype_to_FASTphenotype(chromosome):
    """ Given a chromosome, create a FAST simulation folder
        Uses the global variables: RefValues, CH_MAP, ref_dir
    """
    if not hasattr(chromosome,'data'):
        raise NotImplementedError('')

    def naming(p):
        return '_{:02.0f}'.format(p['InflowFile|HWindSpeed'])

    if len(chromosome)!=CH_MAP.nBases:
        raise Exception('Chromosome length ({}) not compatible with CH_MAP length ({})'.format(len(chromosome),CH_MAP.nBases))

    #print('')
    #for gm,gene in zip(CH_MAP, CH_MAP.split(chromosome)):
    #    print(gm.show_full_raw(gene))

    PARAMS=[]
    for wsp,rpm,pit in zip(RefValues['WS'],RefValues['RPM'],RefValues['Pitch']):
        p=dict()
        if wsp<6:
            p['FAST|TMax']         = 28
        elif wsp<9:
            p['FAST|TMax']         = 23
        else:
            p['FAST|TMax']         = 20
        p['FAST|DT']               = 0.01
        p['FAST|DT_Out']           = 0.1
        p['FAST|OutFileFmt']       = 1 # TODO
        p['EDFile|RotSpeed']       = rpm
        p['EDFile|BlPitch(1)']     = pit
        p['EDFile|BlPitch(2)']     = pit
        p['EDFile|BlPitch(3)']     = pit
        p['InflowFile|HWindSpeed'] = wsp
        p['InflowFile|WindType']   = 1 # Setting steady wind
        for gm,gene in zip(CH_MAP, CH_MAP.split(chromosome)):
            if gm.kind=='fast_param':
                p[gm.name]=gm.decode(gene)
            elif gm.kind=='builtin':
                if gm.name=='pitch':
                    p['EDFile|BlPitch(1)']     = gm.decode(gene)
                    p['EDFile|BlPitch(2)']     = gm.decode(gene)
                    p['EDFile|BlPitch(3)']     = gm.decode(gene)
                #print(gm.name, '->',p[gm.name],gene)
    
        PARAMS.append(p)

    sim_dir=chromosome.data['dir']
    fastlib.templateReplace(ref_dir,PARAMS,workdir=sim_dir,name_function=naming,RemoveRefSubFiles=True)

    # --- Patching the airfoil files
    for gm,gene in zip(CH_MAP, CH_MAP.split(chromosome)):
        if gm.kind=='airfoil':
            af_mut = copy.deepcopy(gm.meta)
            af_mut = patch_airfoil(gene,af_mut)
            AF = weio.FASTInFile(os.path.join(sim_dir,gm.name))
            AF['AFCoeff'] = af_mut['polar']
            AF.write()



### --- PARAMETERS
def individualFitness(chromosome,outdir=None,ForceEvaluation=False,stat=''):
    """ 
    Evaluate an individual.
    Global variable used so far: CH_MAP, WS_SIM, GA_DIR, TIMEAVGWINDOW, EXE, ref_dir
    
    """
    # Basic data
    ID      = galib.chromID(chromosome)
    if outdir is not None:
        sim_dir = outdir
    else:
        sim_dir = os.path.normpath(os.path.join(ref_dir,'..',GA_DIR,''+ID))
    chromosome.data = dict()
    chromosome.data['ID']   = ID
    chromosome.data['dir']  = sim_dir 
    #print('--- EVALUATING {} '.format(chromosome.data['ID']),end='')
    print('--- {}{} '.format(stat,CH_MAP.show_full(chromosome,' ')),end='')
    # Checking if all files are present and with non zero size in the sim folder
    outfiles = glob.glob(os.path.join(sim_dir,'*.out'))
    bFilesOK=False
    if (len(outfiles)==len(WS_SIM)):
        bFilesOK=True
        for fo in outfiles:
            bFilesOK=bFilesOK and os.stat(fo).st_size > 0
    if os.path.exists(sim_dir) and bFilesOK and (not ForceEvaluation):
        # No need to rerun
        print(' >>> ', end='') 
    else:
        fast_files = genotype_to_FASTphenotype(chromosome)
        #print('  Running fast in sim folder')
        run_sim(sim_dir,FAST,len(WS_SIM),exe=EXE)
        print(' --- ', end='')

    ## Evaluating performances and fitnesses
    chromosome.data['perf'] = postpro_simdir(sim_dir, TimeAvgWindow = TIMEAVGWINDOW, FAST = FAST)
    perf_err,_ = get_perf_error(chromosome.data['perf'], PerformanceSignals, perf_ref=RefValues)
    fits = [v for v in perf_err.values]
    print(' Fit: ['+','.join(['{:5.2f}'.format(f) for f in fits])+' ]')
    return fits

def evalNeutralChromosome(outdir=None,ForceEvaluation=False):
    print('Neutral chromosome...')
    neutral=galib.Indiv()
    neutral[:]=CH_MAP.neutralChromosome()
    neutral.fitness.values=individualFitness(neutral,outdir=outdir,ForceEvaluation=ForceEvaluation)
    NeutralValues=neutral.data['perf']
    NeutralValues.columns=[c+'_ori' for c in NeutralValues.columns.values]
    return neutral,NeutralValues

def exportBestData(best,best_dir_dest, RefValues=None, NeutralValues=None):
    with open(os.path.join(best_dir_dest,'chromosome.csv'),'w') as f:
        v=[best.data['ID']]+[v for v in best.fitness.values]+best
        sv = ', '.join([str(val) for val in v])
        f.write(sv+'\n')
        f.write(CH_MAP.show_full(best))
    Vals=[]
    if RefValues is not None:
        Vals+=[RefValues]
    Vals+=[best.data['perf']]
    if NeutralValues is not None:
        Vals+=[NeutralValues]
    #df = pd.concat([RefValuesNewCol,best.data['perf']],axis=1)
    df = pd.concat(Vals,axis=1)
    df.to_csv(os.path.join(best_dir_dest,'Results.csv'),index=False)


# --------------------------------------------------------------------------------}
# --- PARAMETER DEFINITIONS 
# --------------------------------------------------------------------------------{
### --- PARAMETERS
RESOLUTION = 1000; # resolution for variable range between min and max
FAST=1
GA_DIR  = '_GA_Parametric'
ref_dir = 'OpenFAST_V27_v2_forGA/'
WS_SIM  = np.array([4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10])
TIMEAVGWINDOW=None
if platform.system()=='Linux':
    EXE     = '_Exe/gcc'
else:
    EXE     = '_Exe/OpenFAST2_x64s_ebra.exe'  ;
RefFile = '_data/swiftData_Half_Binned.csv'

# --- CHOICE OF TUNING PARAMETERS
PerformanceSignals = ['RPM','FlapM','Pgen']

CH_MAP=galib.ChromosomeMap()
# Airfoils
# airfoilFileNames = ['AD_5_63-214_mod.dat','AD_4_63-218_mod.dat','AD_3_63-224_mod.dat','AD_2_63-235_mod.dat']
# airfoilFileNames = [os.path.join('AeroData_AD15',a) for a in airfoilFileNames]
# for af in airfoilFileNames:
#     af_ref=read_airfoils([af],workdir=ref_dir)[0]
#     gene_info=galib.GeneMap(nBases=3, kind='airfoil',name=af, meta=af_ref, protein_ranges=[[0,1],[0,1],[0,1]], protein_neutr=[0.5,0.5,0])
#     CH_MAP.append(gene_info)
# Fast params
CH_MAP.add(galib.GeneMap(nBases=1, kind='fast_param', name='ServoFile|GenEff'  ,protein_ranges=[[90,100]]       , protein_neutr=[94], resolution=RESOLUTION ))
CH_MAP.add(galib.GeneMap(nBases=1, kind='fast_param', name='EDFile|GBoxEff'    ,protein_ranges=[[90,100]]       , protein_neutr=[94], resolution=RESOLUTION ))
CH_MAP.add(galib.GeneMap(nBases=1, kind='fast_param', name='ServoFile|VS_Rgn2K',protein_ranges=[[0.0003,0.0005]], protein_neutr=[0.00038245], resolution=RESOLUTION ))
CH_MAP.add(galib.GeneMap(nBases=1, kind='builtin', name='pitch',protein_ranges=[[-2,3]], protein_neutr=[0.0] , resolution=RESOLUTION))

print('Number of Bases    :',CH_MAP.nBases)
print('Number of Genes    :',CH_MAP.nGenes)
print('Neutral chromosome :',CH_MAP.neutralChromosome())
print('Neutral protein    :',CH_MAP.neutralProtein())
print(CH_MAP)

# --------------------------------------------------------------------------------}
# --- Derived params 
# --------------------------------------------------------------------------------{
# -- Reference operating conditions and performance values (table as function of WS)
RefValues=pd.read_csv(RefFile)
RefValues['Pitch']=RefValues['Pitch']*0+1
RefValues=pandalib.pd_interp1('WS',RefValues,WS_SIM)
RefValuesNewCol=RefValues.copy()
RefValuesNewCol.columns=[c+'_ref' for c in RefValues.columns.values]
print('Simulations:')
print(RefValues)



# --------------------------------------------------------------------------------}
# --- Parametric run and minimization
# --------------------------------------------------------------------------------{
# --- Parametric GA
fits_norm,fits_arr,pop,v,vProt=galib.parameticGA(individualFitness,CH_MAP,[5,5,5,5],len(PerformanceSignals), resolution=RESOLUTION)
bnds     = tuple([(m+1.e-6,M-1e-6) for m,M in CH_MAP.chromosomeBounds()])
print('Neutral chromosome:',CH_MAP.neutralChromosome())
print('v',v)
print('v',bnds)

# --- Find minimum
fInterp = RegularGridInterpolator(tuple(v), fits_norm)
res     = minimize(fInterp    , tuple(CH_MAP.neutralChromosome()), bounds = bnds)
print('Minimum chromosome:',res.x,'Min:',fInterp(res.x))
print('Minimum protein   :',CH_MAP.decode(res.x))
print('Minimum protein   :',CH_MAP.show_full(res.x))

# --- Evaluating Neutral
neutral_dir_dest='_GA_Parametric/_Neutral'
neutral_ori,NeutralValues=evalNeutralChromosome(outdir=neutral_dir_dest,ForceEvaluation=True)

# --- Evaluating best
print('Evaluating best from minimization:')
best_dir_dest='_GA_Parametric/_Best'
best = galib.Indiv(res.x);
best.fitness.values = individualFitness(best,outdir=best_dir_dest,ForceEvaluation=True)
exportBestData(best, best_dir_dest, RefValuesNewCol, NeutralValues)
# print(fits)
# from mpl_toolkits.mplot3d import Axes3D 
# x=v[0]
# y=v[1]
# X, Y = np.meshgrid(x, y)
# Z = fits.transpose()
# 
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.plot_surface(X, Y, Z)
# ax.set_xlabel('x')
# ax.set_ylabel('y')
# plt.show()
# 
