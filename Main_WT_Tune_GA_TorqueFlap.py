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
        p['FAST|TMax']             = 30
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
    try:
        # If creator exists, we use it
        toolbox = base.Toolbox()
        toolbox.register("attr_float" , random.random)
        toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_float, CH_MAP.nBases) 
        neutral=toolbox.individual()
    except:
        # we a generic implementation (will fail on some of the deap call)
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
        dfRef=RefValues.copy()
        dfRef.columns=[c.replace('_ref','') for c in dfRef.columns]
        dfRef.columns=[c.replace('Loss','') for c in dfRef.columns]
        dfRef=dfRef[['WS', 'RPM', 'Pitch' , 'Pgen', 'Qgen', 'FlapM']]
        dfRef.to_csv(os.path.join(best_dir_dest,'ResultsMeasurements.csv'),index=False)

    Vals+=[best.data['perf']]

    if NeutralValues is not None:
        Vals+=[NeutralValues]
        dfNeut=NeutralValues.copy()
        dfNeut.columns=[c.replace('_ori','') for c in dfNeut.columns]
        dfNeut=dfNeut[['WS', 'RPM', 'Pitch' , 'Pgen', 'Qgen', 'FlapM']]
        dfNeut.to_csv(os.path.join(best_dir_dest,'ResultsPrevious.csv'),index=False)

    dfVals=best.data['perf'].copy()
    dfVals=dfVals[['WS', 'RPM', 'Pitch' , 'Pgen', 'Qgen', 'FlapM']]
    dfVals.to_csv(os.path.join(best_dir_dest,'ResultsNew.csv'),index=False)
    # --- All
    df = pd.concat(Vals,axis=1)
    df.to_csv(os.path.join(best_dir_dest,'ResultsAll.csv'),index=False)

def individualPlot(chromosome,fig=None):
    if not fig:
        fig=plt.figure()
    # Plotting polars
    for a in fig.axes:
        a.clear()
    if len(airfoils_ref)>0:
        # We have to reread the data
        airfoils_mut = read_airfoils(airfoilFileNames,workdir=chromosome.data['dir'])
        fig=plot_airfoils(airfoils_ref,color='k',fig=fig, names=[os.path.splitext(os.path.basename(fn))[0] for fn in airfoilFileNames])
        fig=plot_airfoils(airfoils_mut,fig=fig)

def performancePlot(population,ax=None,kind='',title=''):
    """ Plot the ref signals and the population values for the wanted data
        Typically: Cp,Ct as function of WS
    """
    colrs=['b','r','g','m']
    if kind=='new':
        ax.clear()
    for [ii,c] in enumerate(population):
        data_Scaled=c.data['perf'][['WS']+IPlot]/PerfScale
        if kind=='neutral':
            data_Scaled.plot(x='WS',y=IPlot, ax=ax,color='b',linestyle='--',linewidth=0.5,legend=False)
        elif kind=='new':
            data_Scaled.plot(x='WS',y=IPlot, ax=ax,color=[(0.7,0.7,0.7)],legend=False)
        elif kind=='best':
            lgd=(ii==0) and (title is not None)
            for ip,cl in zip(IPlot,colrs):
                data_Scaled.plot(x='WS',y=ip, ax=ax, color=cl, title=title,label=ip,legend=lgd)
        else:
            raise Exception('Unknown kind {}'.format(kind))

    if kind=='new':
        PlotValues = RefValues[['WS']+IPlot]/PerfScale
        for ip,cl in zip(IPlot,colrs):
            PlotValues.plot(x='WS',y=ip,color=cl,marker='o',linestyle='', ax=ax,legend=False)
        ax.set_xlabel('WS')
        #ax.axis([0,1,0.0,1.2])


def populationPlot(pop,stats=None,fits=None,fig=None,kind='new',title=None):
    """ 
    ax[0]: Performance plot (e.g. signals as function of WS)
    ax[1]: Fitnesses as function of number of generation
    ax[2]: Chromosome plot
    ax[3]: Population fitnesses in X,Y plot (only for 2 fitnesses)
    """
    if not fig:
        fig=plt.figure()
    if len(fig.axes)==0:
        fig.add_subplot(221)
        fig.add_subplot(222)
        fig.add_subplot(223)
        fig.add_subplot(224)
    ax = fig.axes
    nInd=len(pop)
    XYC=np.zeros((nInd*CH_MAP.nBases,3))
    FFC=np.zeros((nInd,3))
    c=0
    if kind=='all':
        SpreadFact=0.5
    else:
        SpreadFact=0.5
    if kind=='offspring':
        # we don't care about fitness
        for j,p in enumerate(pop):
            for i in range(len(p)):
                XYC[c,0]=i
                XYC[c,1]=p[i]
                c += 1
    else:
        # We care about fitness
        if fits is None:
            fits=[]
            for p in pop:
                fits.append(p.fitness.values)

        for p,f,j in zip(pop,fits,range(nInd)):
            for i in range(len(p)):
                XYC[c,0]=i
                XYC[c,1]=p[i]
                XYC[c,2]=1-sum(f)
                c += 1
            FFC[j,0]=f[0]
            FFC[j,1]=f[1]
            FFC[j,2]=1-sum(f)

    if kind=='all':
        for a in ax:
            a.clear()
        ax[2].scatter(XYC[:,0],XYC[:,1],color=[0.6,0.6,0.6],marker='.',s=15)
        ax[3].scatter(FFC[:,0],FFC[:,1],color=[0.6,0.6,0.6],marker='.',s=15)
        ax[3].grid()
        ax[2].axis([-0.1, CH_MAP.nBases+0.1,-0.1, 1.1])
        ax[3].axis([1e-2, 100, 1e-2, 100])
        ax[3].set_yscale('log')
        ax[3].set_xscale('log')
        ax[3].set_xlabel('Fitness '+PerformanceSignals[0])
        ax[3].set_ylabel('Fitness '+PerformanceSignals[1])
        # Separating Genes- TODO
        #for iBase in range(nBasePerGene):
        #    ax[2].plot([iBase*nGenes,iBase*nGenes],[0,1],'k--')
        if title is not None:
            fig.suptitle(title)

    elif kind=='new':
        ax[2].scatter(XYC[:,0],XYC[:,1],c=XYC[:,2],cmap="Reds",vmin=0.2,vmax=1,s=8)
        ax[3].scatter(FFC[:,0],FFC[:,1],c=FFC[:,2],cmap="Reds",vmin=0.2,vmax=1,s=8)
        performancePlot(population=pop,ax=ax[0],kind=kind)
    elif kind=='parents':
        ax[2].scatter(XYC[:,0],XYC[:,1],color=[0.7,0.7,0.7],marker='o',s=80, facecolors='none', edgecolor=[0.7,0.7,0.7])
        ax[3].scatter(FFC[:,0],FFC[:,1],color=[0.7,0.7,0.7],marker='o',s=80, facecolors='none', edgecolor=[0.7,0.7,0.7])
    elif kind=='selected':
        ax[2].scatter(XYC[:,0],XYC[:,1],color='k',marker='o',s=80, facecolors='none', edgecolor='k')
        ax[3].scatter(FFC[:,0],FFC[:,1],color='k',marker='o',s=80, facecolors='none', edgecolor='k')
    elif kind=='offspring':
        ax[2].scatter(XYC[:,0],XYC[:,1],color='k',marker='o',s=4)
    elif kind=='extr':
        ax[2].scatter(XYC[:,0],XYC[:,1],color='g',marker='o',s=7)
        ax[3].scatter(FFC[:,0],FFC[:,1],color='g',marker='o',s=5)
    elif kind=='best':
        performancePlot(pop,ax=ax[0],kind=kind,title=title)

        ax[2].scatter(XYC[:,0],XYC[:,1],color='b',marker='o',s=9)
        ax[3].scatter(FFC[:,0],FFC[:,1],color='b',marker='o',s=9)
    elif kind=='neutral':
        ax[2].scatter(XYC[:,0],XYC[:,1],color='b',marker='+',s=50)
        ax[3].scatter(FFC[:,0],FFC[:,1],color='b',marker='+',s=50)
    else:
        raise Exception('Unknown kind for population Plot `{}`'.format(kind))

    if stats is not None:
        ax[1].set_yscale('log')
        colrs=['b','r','g','m']
        for i,s in enumerate(stats):
            if i==0:
                lgd=True
            else:
                lgd=False
            s.plot(y='Best',color=colrs[i],linestyle='-' ,ax=ax[1],label=PerformanceSignals[i]+' Best',legend=True)
            s.plot(y='Min', color=colrs[i],linestyle='--',ax=ax[1],label=PerformanceSignals[i]+' Min' ,legend=lgd)
            s.plot(y='Max', color=colrs[i],linestyle='--',ax=ax[1],label=PerformanceSignals[i]+' Max' ,legend=lgd)
            s.plot(y='Mean',color=colrs[i],linestyle=':' ,ax=ax[1],label=PerformanceSignals[i]+' Mean',legend=lgd)
        ax[1].set_title('Population Fitnesses')
        ax[1].set_xlabel('Number of generations')
        ax[1].set_ylabel('Fitnesses values')

    
def individualDistance(ind1,ind2):
    return np.linalg.norm(np.array(ind1)-np.array(ind2))

def populationDiversity(pop):
    return [sum([individualDistance(x,y) for x in pop])/(np.sqrt(2)*(len(pop)-1)) for y in pop]



def getRandomPop(n=3):
    toolbox = base.Toolbox()
    toolbox.register("attr_float" , random.random)
    toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_float, CH_MAP.nBases)
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)
    return toolbox.population(n=n)


def mainGA(nBase=2,nInd=10,nIndSelect=10,CXPB=0.3,MUTPB=0.3,nIterMax=100,nPerTournament=2,MutStep=0.1,MutMethod='PolyBound',CXFunction=tools.cxTwoPoint, MutParam1=0.05):
    # Global variables for plotting 
    #random.seed(64)
    random.seed(59)
    toolbox = base.Toolbox()
    toolbox.register("attr_float" , random.random)
    toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_float, nBase) 
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)
    toolbox.register("evaluate", individualFitness)
    toolbox.register("mate"    , CXFunction)
    #toolbox.register("mutate"  , tools.mutFlipBit, indpb=0.05)
    if MutMethod=='PolyBound':
        toolbox.register("mutate", tools.mutPolynomialBounded,eta=MutParam1,low=0,up=1,indpb=0.50)
    elif MutMethod=='GaussianBound':
        toolbox.register("mutate", galib.mutGaussianBounded, mu=0,sigma=MutParam1,low=0,up=1,indpb=0.5)
    elif MutMethod=='UniBound':
        toolbox.register("mutate", galib.mutUniformBounded, low=0,up=1,indpb=0.5)
    else:
        raise Exception('unknown Mut method:'+MutMethod)

    #toolbox.register("select"  , tools.selTournament, tournsize=nPerTournament) #NSGA2, Best
    #toolbox.register("select"  , tools.selBest)
    #toolbox.register("select"  , tools.selNSGA2)
    #toolbox.register("select"  , tools.selSPEA2)
    toolbox.register("select"  , tools.selSPEA2)

    # Create and Evaluate the neutral chromosom
    neutral_dir_dest = os.path.normpath(os.path.join(ref_dir,'..',GA_DIR,'_Neutral'))
    neutral_ori,NeutralValues=evalNeutralChromosome(outdir=neutral_dir_dest)
    # Put neutral in dedicateed directory
    #neutral_dir_src  = os.path.normpath(os.path.join(ref_dir,'..',GA_DIR,''+neutral_ori.data['ID']))
    #neutral_dir_dest = os.path.normpath(os.path.join(ref_dir,'..',GA_DIR,'_Neutral'))
    #distutils.dir_util.copy_tree(neutral_dir_src, neutral_dir_dest)

    # create an initial population
    pop = toolbox.population(n=nInd)
    if bEnforceNeutral:
        pop.append(galib.clone(neutral_ori))
    pop = galib.populationTrimAccuracy(pop,int(np.log10(RESOLUTION)))
    
    # --- EVALUATE the entire population
    for ind in pop:
        ind.fitness.values = toolbox.evaluate(ind)
        #ind.fitness.values = [fit[0]*FitScale +  div*DivScale]
    #galib.populationPrint(pop,1,'INITIAL')
    DB_file      = galib.populationSave(pop,newfile=True,directory=DB_DIR)
    DB_file_Best = galib.populationSave(pop,newfile=True,directory=DB_DIR,basename='GA_DB_BEST_')

    if len(figs)>0:
        populationPlot(pop,fig=figs[0],kind='new')
        plt.pause(0.001)
    best = galib.selBestByNorm(pop,weights=objectiveWeights)
    stats,_ = galib.populationStats(pop,best)
    print('---------------------------------------------------')
    # Variable keeping track of the number of generations
    g = 0
    while g < nIterMax:
        g = g + 1
        # --- SUPER INDIVIDUALS
        # first we select the extremes bests
        best = galib.selBestByNorm     (pop,weights=objectiveWeights)
        extr = galib.selIndependentBest(pop,weights=objectiveWeights)
        enforce=[]
        if bEnforceBests:
            enforce+=[best]+extr
        if bEnforceNeutral:
            enforce+=[neutral_ori]

        # --- PARENTS SELECTION FROM POP 
        parents = list(map(galib.clone, toolbox.select(pop, nIndSelect)  ))
        #populationPrint(offspring,'SELECTED')
        # --- ENFORCING PREVIOUS BEST IN POP
        for e in enforce:
            galib.addIfAbsent(parents,e,'parents')

        if len(figs)>0:
            populationPlot(parents      ,fig = figs[0], kind = 'selected')
            populationPlot([neutral_ori],fig = figs[0], kind = 'neutral' )
            populationPlot([best]       ,fig = figs[0], kind = 'best'    )
            populationPlot(extr         ,fig = figs[0], kind = 'extr'    )
            plt.pause(0.001)

        # --- CROSSOVER of offsprings, with probability CXPB or Mating
        #galib.mateInPlace(parents, toolbox.mate, CPPB)
        offspring=galib.mate(parents, toolbox.mate, nInd)
        nCX=len(offspring)

        # --- ENFORCING BESTS IN NOT-MUTATED OFFSPRING
        for e in enforce:
            galib.addIfAbsent(offspring,e,'clean offspring')

        # --- MUTATION with probability MUTPB 
        nMute = galib.mutateInPlace(offspring,toolbox.mutate,MUTPB)

        # --- ENFORCING PREVIOUS BEST IN OFFSPRING
        for e in enforce:
            galib.addIfAbsent(offspring,e,'mutated offspring')
        offspring = galib.populationTrimAccuracy(offspring,int(np.log10(RESOLUTION)))
        if len(figs)>0:
            populationPlot(offspring,fig=figs[0],kind='offspring')
            plt.pause(0.001)
        # --- EVALUATE FITNESS and DIVERSITY
        # SMART INVALID
        invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
        #diversities = populationDiversity(invalid_ind)
        #diversities = diversities/max(diversities)
        for ind in invalid_ind:
            ind.fitness.values = toolbox.evaluate(ind)
             
        best = galib.selBestByNorm(offspring,weights=objectiveWeights)
        # --- Save Best in dedicated directory
        best_dir_src  = os.path.normpath(os.path.join(ref_dir,'..',GA_DIR,''+best.data['ID']))
        best_dir_dest = os.path.normpath(os.path.join(ref_dir,'..',GA_DIR,'_Best'))
        distutils.dir_util.copy_tree(best_dir_src, best_dir_dest)
        exportBestData(best, best_dir_dest, RefValuesNewCol, NeutralValues)

        # --- STATS
        recent_stats,stats = galib.populationStats(pop,best,stats)
        # --- INFO STRING
        sGen  = "Gen:{0:d} ".format(g)
        sInfo = " Xd:{0:d}  Mut:{1:d}  X+M:{2:d}  Pop:{3:d} ".format(nCX, nMute, len(invalid_ind), len(offspring))
        #sStat = " [{0:.1f} {1:.1f}] Avg: {2:.2f} Std {3:.2f} Div: {5:.2f} {6:.2f} Best: {4:.8f}".format(min(fits), max(fits), mean, std,bestFit,np.mean(diversities),np.max(diversities))
        #print(sGen+sInfo+sStat)
        sBest='Best:'+' '.join(['{:.4f}'.format(x) for x in best.fitness.values])+' '
        print(sGen+sInfo+sBest)
        # --- PLOTTING, STORING
        print('Plotting saving....')
        galib.populationSave(offspring,directory=DB_DIR)
        galib.populationSave([best]   ,directory=DB_DIR,basename='GA_DB_BEST_')
        if len(figs)>0:
            if len(figs)>2:
                individualPlot(best,fig=figs[2])
            #populationPrint(offspring,'NEW')
            dfAll,pop_uniq  = galib.populationLoad(DB_file     ,nFits=len(pop[0].fitness.values))
            dfBest,pop_best = galib.populationLoad(DB_file_Best,nFits=len(pop[0].fitness.values))
            #galib.timelinePlot(df,fig=figs[2])
            galib.timelinePlot(dfBest,fig=figs[1],noFit=True)
            populationPlot(pop_uniq, stats=stats,fig=figs[0],kind='all',title=sGen+sInfo+sBest)
            populationPlot(parents              ,fig=figs[0],kind='parents')
            populationPlot(offspring            ,fig=figs[0],kind='new')
            populationPlot([best]               ,fig=figs[0],kind='best',title=sGen+sBest)
            plt.pause(0.001)

        # --- UPDATE - The population is entirely replaced by the offspring
        pop[:] = offspring
        #time.sleep(1.0)

    if g==nIterMax :
        print(">>>>>>>>>>>>>> Maximum number of iterations reached")

    return pop,best


# --------------------------------------------------------------------------------}
# --- PARAMETER DEFINITIONS 
# --------------------------------------------------------------------------------{
### --- PARAMETERS
RESOLUTION = 1000; # resolution for variable range between min and max
FAST=1
if FAST==1:
    GA_DIR  = '_GA_AF'
    DB_DIR  = '_GA_AF_DB'
    ref_dir = 'OpenFAST_V27_v2_forGA/'
    WS_SIM  = np.array([4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10])

else:
    GA_DIR  = '_GA_Runs_AD'
    DB_DIR  = '_GA_Runs_AD_DB'
    ref_dir = 'AeroDyn_V27_v1_refMulti/'
    #WS_SIM  = np.array([ 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25]);
    WS_SIM  = np.array([ 5, 7, 9, 11, 13, 15, 17, 19]);
TIMEAVGWINDOW=None
if platform.system()=='Linux':
    EXE     = '_Exe/gcc'
else:
	EXE     = '_Exe/OpenFAST2_x64s_ebra.exe'  ;
RefFile = '_data/swiftData_Half_Binned.csv'
# --- GA Specific options
bEnforceBests=True
bEnforceNeutral=True

# --- CHOICE OF TUNING PARAMETERS
PerformanceSignals = ['RPM','FlapM','Pgen']

CH_MAP=galib.ChromosomeMap()
# Airfoils
# airfoilFileNames = ['AD_5_63-214_mod.dat','AD_4_63-218_mod.dat','AD_3_63-224_mod.dat','AD_2_63-235_mod.dat']
# airfoilFileNames = [os.path.join('AeroData_AD15',a) for a in airfoilFileNames]
airfoilFileNames = ['AD_3_63-224_mod.dat','AD_2_63-235_mod.dat']
airfoilFileNames = [os.path.join('AeroData_AD15',a) for a in airfoilFileNames]
airfoils_ref=[]
for af in airfoilFileNames:
    af_ref=read_airfoils([af],workdir=ref_dir)[0]
    gene_info=galib.GeneMap(nBases=3, kind='airfoil',name=af, meta=af_ref, protein_ranges=[[0,1],[0,1],[0,1]], protein_neutr=[0.5,0.5,0])
    airfoils_ref.append(af_ref)
    CH_MAP.append(gene_info)
# Fast params
CH_MAP.add(galib.GeneMap(nBases=1, kind='fast_param', name='ServoFile|GenEff'  ,protein_ranges=[[90,99]]       , protein_neutr=[94], resolution=RESOLUTION ))
CH_MAP.add(galib.GeneMap(nBases=1, kind='fast_param', name='EDFile|GBoxEff'    ,protein_ranges=[[90,99]]       , protein_neutr=[94], resolution=RESOLUTION ))
CH_MAP.add(galib.GeneMap(nBases=1, kind='fast_param', name='ServoFile|VS_Rgn2K',protein_ranges=[[0.0003,0.0005]], protein_neutr=[0.00038245], resolution=RESOLUTION ))
CH_MAP.add(galib.GeneMap(nBases=1, kind='builtin'   , name='pitch'             ,protein_ranges=[[-0.5,2]], protein_neutr=[1.665913124] , resolution=RESOLUTION))

print('Number of Bases    :',CH_MAP.nBases)
print('Number of Genes    :',CH_MAP.nGenes)
print('Neutral chromosome :',CH_MAP.neutralChromosome())
print('Neutral protein    :',CH_MAP.neutralProtein())
print(CH_MAP)

# --------------------------------------------------------------------------------}
# --- Derived params 
# --------------------------------------------------------------------------------{
# --- Derived Params
IPlot              = PerformanceSignals
objectiveWeights=(-1.0,)*len(PerformanceSignals) # negative = minimalization 
# -- Reference operating conditions and performance values (table as function of WS)
RefValues=pd.read_csv(RefFile)
RefValues['Pitch']=RefValues['Pitch']*0+1
RefValues=pandalib.pd_interp1('WS',RefValues,WS_SIM)
PerfScale = RefValues[['WS']+IPlot].max()
PerfScale['WS'] = PerfScale['WS']*0+1
RefValuesNewCol=RefValues.copy()
RefValuesNewCol.columns=[c+'_ref' for c in RefValues.columns.values]
print('Simulations:')
print(RefValues)




# --------------------------------------------------------------------------------}
# --- Full GA
# --------------------------------------------------------------------------------{
### --- INIT
figs=[]
# figs+=figlib.fig_grid(AreaName='Left',ScreenName='RightScreen')
# figs+=figlib.fig_grid(2,1,AreaName='Right',ScreenName='RightScreen')
# figs+=figlib.fig_grid(AreaName='TopRight',ScreenName='RightScreen')
# plt.ion()
# plt.show()
# 

# --- Full GA
creator.create("Fitness", base.Fitness, weights=objectiveWeights)
creator.create("Individual", list, fitness=creator.Fitness)
pop,best_ind=mainGA(nBase=CH_MAP.nBases
          ,nInd=32,nIndSelect=16,CXPB=0.5,MUTPB=0.5
#          ,MutMethod='PolyBound',MutParam1=0.001,nPerTournament=2
#            ,MutMethod='UniBound',MutParam1=np.nan,nPerTournament=2
          ,MutMethod='GaussianBound',MutParam1=0.01,nPerTournament=2
          ,nIterMax=100);

###
# #neutral=getNeutralChromosome();
# #individualPlot(neutral,fig1=figs[0],fig2=figs[1])
##pop,pop_init,best_ind=main(nBase=,nInd=10,CXPB=0.5,MUTPB=0.9,nIterMax=100,nPerTournament=2);
##pop,pop_init,best_ind=main(nBase=,nInd=30,CXPB=0.5,MUTPB=0.5,nIterMax=100,nPerTournament=2);
#individualPlot(neutral,fig1=figs[0],fig2=figs[1])
