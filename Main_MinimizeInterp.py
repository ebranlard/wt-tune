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
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib
GA_DIR='_Eagle'
fits_norm   = np.load(os.path.join(GA_DIR,'fits_norm.npy'))
fits_arr    = np.load(os.path.join(GA_DIR,'fits_arr.npy') )
vBase       = np.load(os.path.join(GA_DIR,'vBase.npy')    )
vProt       = np.load(os.path.join(GA_DIR,'vProt.npy')    )
NeutProt    = np.load(os.path.join(GA_DIR,'neutProt.npy')    )
NeutChrom   = np.load(os.path.join(GA_DIR,'neutChrom.npy')    )
boundsProt  = np.load(os.path.join(GA_DIR,'boundsProt.npy')    )
boundsChrom = np.load(os.path.join(GA_DIR,'boundsChrom.npy')    )
geneNames   = np.load(os.path.join(GA_DIR,'geneNames.npy' ))
nBases      = np.load(os.path.join(GA_DIR,'nBases.npy' ))
bndsProt     = tuple([(m+1.e-6,M-1e-6) for m,M in boundsProt])

# --- Find minimum
fInterp = RegularGridInterpolator(tuple(vProt), fits_norm)
res     = minimize(fInterp    , tuple(NeutProt), bounds = bndsProt)
# print('Minimum chromosome:',res.x,'Min:',fInterp(res.x))
# print('Minimum protein   :',CH_MAP.decode(res.x))
print('Minimum protein   :',res.x)
# 

geneNames= [g.split('|')[-1] for g in geneNames]

I=[0,1]

colrs=['b','r','g','m']

def myplot3d(ix=0,iy=1,iz=2):
    x=vProt[ix]
    y=vProt[iy]
    z=vProt[iz]
    X, Y = np.meshgrid(x, y)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    if ix==0 and iy==1 and iz==2:
        for k,z in enumerate(z) :
            Z = fits_norm[:,:,k,0].transpose()
            ax.plot_surface(X, Y, Z,color=colrs[k])
            ax.plot(X[0],Y[0],Z[0], linestyle="none", c=colrs[k], marker = 'o',label='{}={}'.format(geneNames[iz],z))
    elif ix==0 and iy==2 and iz==3:
        for k,z in enumerate(z) :
            Z = fits_norm[:,0,:,k].transpose()
            ax.plot_surface(X, Y, Z,color=colrs[k])
            ax.plot(X[0],Y[0],Z[0], linestyle="none", c=colrs[k], marker = 'o',label='{}={}'.format(geneNames[iz],z))
    elif ix==1 and iy==2 and iz==3:
        for k,z in enumerate(z) :
            Z = fits_norm[0,:,:,k].transpose()
            ax.plot_surface(X, Y, Z,color=colrs[k])
            ax.plot(X[0],Y[0],Z[0], linestyle="none", c=colrs[k], marker = 'o',label='{}={}'.format(geneNames[iz],z))
    ax.set_xlabel(geneNames[ix])
    ax.set_ylabel(geneNames[iy])
    ax.legend()

myplot3d(ix=0,iy=1,iz=2)
myplot3d(ix=0,iy=2,iz=3)
myplot3d(ix=1,iy=2,iz=3)

plt.show()

# # 
# 
# # def f(x, y):
# #     return  (x-2)**2 + (y-5)**2 
# # 
# # x = np.linspace(1, 4, 11)
# # y = np.linspace(4, 7, 22)
# # data = f(*np.meshgrid(x, y, indexing='ij', sparse=True))
# # print(data.shape)
# # print(type(data))
# # 
# # # data is now a 3D array with data[i,j,k] = f(x[i], y[j], z[k]). Next, define an interpolating function from this data:
# # 
# # fInterp = RegularGridInterpolator((x, y), data)
# # 
# # # Evaluate the interpolating function at the two points (x,y,z) = (2.1, 6.2, 8.3) and (3.3, 5.2, 7.1):
# # # >>>
# # 
# # 
# # 
# # X, Y = np.meshgrid(x, y)
# # # zs = np.array([fun(x,y) for x,y in zip(np.ravel(X), np.ravel(Y))])
# # zs = np.array([f(x,y) for x,y in zip(np.ravel(X), np.ravel(Y))])
# # print(zs.shape)
# # print(type(zs.shape))
# # Z = zs.reshape(X.shape)
# # print(Z.shape)
# # print(type(Z.shape))
# # 
# # ax.plot_surface(X, Y, Z)
# # plt.show()
# # 
# # 
# # pts = np.array([[2.1, 6.2], [3.3, 5.2]])
# # print(fInterp(pts))
# # bnds = ((1, 4), (4, 7))
# # print('Minimizing')
# # res = minimize(fInterp, (2, 5), bounds=bnds)
# # print(res)
# 

# 
# ### --- PARAMETERS
# def individualFitness(chromosome,outdir=None,ForceEvaluation=False,stat=''):
#     """ 
#     Evaluate an individual.
#     Global variable used so far: CH_MAP, WS_SIM, GA_DIR, TIMEAVGWINDOW, EXE, ref_dir
#     
#     """
#     # Basic data
#     ID      = galib.chromID(chromosome)
#     if outdir is not None:
#         sim_dir = outdir
#     else:
#         sim_dir = os.path.normpath(os.path.join(ref_dir,'..',GA_DIR,''+ID))
#     chromosome.data = dict()
#     chromosome.data['ID']   = ID
#     chromosome.data['dir']  = sim_dir 
#     #print('--- EVALUATING {} '.format(chromosome.data['ID']),end='')
#     print('--- {}{} '.format(stat,CH_MAP.show_full(chromosome,' ')),end='')
#     # Checking if all files are present and with non zero size in the sim folder
#     outfiles = glob.glob(os.path.join(sim_dir,'*.out'))
#     bFilesOK=False
#     if (len(outfiles)==len(WS_SIM)):
#         bFilesOK=True
#         for fo in outfiles:
#             bFilesOK=bFilesOK and os.stat(fo).st_size > 0
#     if os.path.exists(sim_dir) and bFilesOK and (not ForceEvaluation):
#         # No need to rerun
#         print(' >>> ', end='') 
#     else:
#         raise Exception('Simulation not run!')
#     ## Evaluating performances and fitnesses
#     chromosome.data['perf'] = postpro_simdir(sim_dir, TimeAvgWindow = TIMEAVGWINDOW, FAST = FAST)
#     perf_err,_ = get_perf_error(chromosome.data['perf'], PerformanceSignals, perf_ref=RefValues)
#     fits = [v for v in perf_err.values]
#     print(' Fit: ['+','.join(['{:5.2f}'.format(f) for f in fits])+' ]')
#     return fits
# 
# # --------------------------------------------------------------------------------}
# # --- PARAMETER DEFINITIONS 
# # --------------------------------------------------------------------------------{
# ### --- PARAMETERS
# RESOLUTION = None; # resolution for variable range between min and max
# FAST=1
# GA_DIR  = '_Eagle/'
# ref_dir = 'OpenFAST_V27_v2_forGA/'
# WS_SIM  = np.array([4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10])
# TIMEAVGWINDOW=None
# if platform.system()=='Linux':
#     EXE     = '_Exe/gcc'
# else:
#     EXE     = '_Exe/OpenFAST2_x64s_ebra.exe'  ;
# RefFile = '_data/swiftData_Half_Binned.csv'
# 
# # --- CHOICE OF TUNING PARAMETERS
# PerformanceSignals = ['RPM','FlapM','Pgen']
# 
# CH_MAP=galib.ChromosomeMap()
# # Airfoils
# # airfoilFileNames = ['AD_5_63-214_mod.dat','AD_4_63-218_mod.dat','AD_3_63-224_mod.dat','AD_2_63-235_mod.dat']
# # airfoilFileNames = [os.path.join('AeroData_AD15',a) for a in airfoilFileNames]
# # for af in airfoilFileNames:
# #     af_ref=read_airfoils([af],workdir=ref_dir)[0]
# #     gene_info=galib.GeneMap(nBases=3, kind='airfoil',name=af, meta=af_ref, protein_ranges=[[0,1],[0,1],[0,1]], protein_neutr=[0.5,0.5,0])
# #     CH_MAP.append(gene_info)
# # Fast params
# CH_MAP.add(galib.GeneMap(nBases=1, kind='fast_param', name='ServoFile|GenEff'  ,protein_ranges=[[90,100]]       , protein_neutr=[94], resolution=RESOLUTION ))
# CH_MAP.add(galib.GeneMap(nBases=1, kind='fast_param', name='EDFile|GBoxEff'    ,protein_ranges=[[90,100]]       , protein_neutr=[94], resolution=RESOLUTION ))
# CH_MAP.add(galib.GeneMap(nBases=1, kind='fast_param', name='ServoFile|VS_Rgn2K',protein_ranges=[[0.0003,0.0005]], protein_neutr=[0.00038245], resolution=RESOLUTION ))
# CH_MAP.add(galib.GeneMap(nBases=1, kind='builtin', name='pitch',protein_ranges=[[-2,3]], protein_neutr=[0.0] , resolution=RESOLUTION))
# 
# print('Number of Bases    :',CH_MAP.nBases)
# print('Number of Genes    :',CH_MAP.nGenes)
# print('Neutral chromosome :',CH_MAP.neutralChromosome())
# print('Neutral protein    :',CH_MAP.neutralProtein())
# print(CH_MAP)
# 
# # --------------------------------------------------------------------------------}
# # --- Derived params 
# # --------------------------------------------------------------------------------{
# # -- Reference operating conditions and performance values (table as function of WS)
# RefValues=pd.read_csv(RefFile)
# RefValues['Pitch']=RefValues['Pitch']*0+1
# RefValues=pandalib.pd_interp1('WS',RefValues,WS_SIM)
# RefValuesNewCol=RefValues.copy()
# RefValuesNewCol.columns=[c+'_ref' for c in RefValues.columns.values]
# print('Simulations:')
# print(RefValues)
# 
# 
# 
# # --------------------------------------------------------------------------------}
# # --- Parametric run and minimization
# # --------------------------------------------------------------------------------{
# # --- Parametric GA
# fits_norm,fits_arr,pop,v,vProt=galib.parameticGA(individualFitness,CH_MAP,[4,4,4,4],len(PerformanceSignals), resolution=RESOLUTION)
# 
# np.save(os.path.join(GA_DIR,'fits_norm.npy'),fits_norm)
# np.save(os.path.join(GA_DIR,'fits_arr.npy') ,fits_arr)
# np.save(os.path.join(GA_DIR,'vBase.npy')    ,v)
# np.save(os.path.join(GA_DIR,'vProt.npy')    ,vProt)
# np.save(os.path.join(GA_DIR,'neutProt.npy'), CH_MAP.neutralProtein())
# np.save(os.path.join(GA_DIR,'neutChrom.npy'), CH_MAP.neutralChromosome())
# np.save(os.path.join(GA_DIR,'boundsChrom.npy'),CH_MAP.chromosomeBounds  ())
# np.save(os.path.join(GA_DIR,'boundsProt.npy' ),CH_MAP.proteinChainBounds())
# np.save(os.path.join(GA_DIR,'geneNames.npy' ), [gm.name for gm in CH_MAP])
# np.save(os.path.join(GA_DIR,'nBases.npy' ), [gm.nBases for gm in CH_MAP])
# # 
# # 
# # bnds     = tuple([(m+1.e-6,M-1e-6) for m,M in CH_MAP.chromosomeBounds()])
# # print('Neutral chromosome:',CH_MAP.neutralChromosome())
# # print('v',v)
# # print('v',bnds)
# # 
