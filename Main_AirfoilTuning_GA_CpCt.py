import copy
import distutils.dir_util
import glob
import matplotlib.pyplot as plt
import numpy as np
import numpy.matlib
import nwtcio
import os
import pandas as pd
import random 
import shutil 
import subprocess
import time
import re
import math
import pdb

from deap import base
from deap import creator
from deap import tools

from AirfoilTuningTools import *
from SWIFTWT import *
from pybra import figlib
from pybra import galib

from collections import Sequence
from itertools import repeat

### --- PARAMETERS
def individualFitness(chromosome):
    ID=galib.chromID(chromosome)
    print('--- EVALUATING {} '.format(ID),end='')
    sim_dir = os.path.normpath(os.path.join(ref_dir,'..',GA_DIR,''+ID))
    outfiles=glob.glob(os.path.join(sim_dir,'*.out'))
    bFilesOK=False
    if (len(outfiles)==len(WS_SIM)):
        bFilesOK=True
        for fo in outfiles:
            bFilesOK=bFilesOK and os.stat(fo).st_size > 0
    if os.path.exists(sim_dir) and bFilesOK:
        print(' >>> ', end='')
    else:
        #print('   Preparing run folder')
        prepare_run_folder(template_dir,sim_dir)
        #print('   Changing polars and saving in run folder')
        airfoils_mut  = patch_airfoil_files(chromosome,airfoils_ref,airfoilFileNames,workdir=sim_dir)
        #print('  Running fast in sim folder')
        run_sim(sim_dir,FAST,len(WS_SIM))
        print(' --- ', end='')

    # Storing Data for convenience
    chromosome.data = dict()
    chromosome.data['ID']   = ID
    chromosome.data['perf'] = postpro_simdir(sim_dir, TimeAvgWindow = TIMEAVGWINDOW, FAST = FAST)
    perf_err,_ = get_perf_error(chromosome.data['perf'], PerformanceSignals, perf_ref=perf_ref)
    fits = [perf_err.values[0],perf_err.values[1]]
    print(' Fit: ',np.around(fits,decimals=4))
    #print(type(fits[0]))
    return fits

def individualPlot(chromosome,fig=None):
    if not fig:
        fig=plt.figure()
    # We have to reread the data
    ID=galib.chromID(chromosome)
    sim_dir = os.path.normpath(os.path.join(ref_dir,'..',GA_DIR,''+ID))
    airfoils_mut = read_airfoils(airfoilFileNames,workdir=sim_dir)

    # Plotting polars
    for a in fig.axes:
        a.clear()
    fig=plot_airfoils(airfoils_ref,color='k',fig=fig)
    fig=plot_airfoils(airfoils_mut,fig=fig)

def basesnames():
    return ['Clp{0} Clm{0} Cd{0}'.format(x+1) for x in range(nGenes)]

def fitsnames():
    return ['Clp{0} Clm{0} Cd{0}'.format(x+1) for x in range(nGenes)]


def performancePlot(population,ax=None,kind='',title=''):
    if kind=='new':
        ax.clear()
    for [ii,c] in enumerate(population):
        if kind=='neutral':
            c.data['perf'].plot(x='WS',y=IPlot, ax=ax,color='b',linestyle='--',linewidth=0.5,legend=False)
        elif kind=='new':
            c.data['perf'].plot(x='WS',y=IPlot, ax=ax,color=[(0.7,0.7,0.7)],legend=False)
        elif kind=='best':
            c.data['perf'].plot(x='WS',y=IPlot, ax=ax, color='g', title=title,legend=False)
        else:
            raise Exception('Unknown kind {}'.format(kind))

    if kind=='new':
        perf_ref.plot(x='WS',y=IPlot,color='k',marker='o',linestyle='', ax=ax)
        ax.axis([4,25,0.0,1.2])


def populationPlot(pop,stats=None,fits=None,fig=None,kind='new',title=None):
    if not fig:
        fig=plt.figure()
    if len(fig.axes)==0:
        fig.add_subplot(221)
        fig.add_subplot(222)
        fig.add_subplot(223)
        fig.add_subplot(224)
    ax = fig.axes
    nInd=len(pop)
    XYC=np.zeros((nInd*nGenes*nBasePerGene,3))
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
                iBase=i % 3
                iGene = i/3
                XYC[c,0]=iBase*nGenes+iGene +j*SpreadFact/nInd
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
                iBase=i % 3
                iGene = i/3
                XYC[c,0]=iBase*nGenes+iGene +j*0.2/nInd
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
        ax[2].axis([-0.1, nGenes*nBasePerGene+0.1,-0.1, 1.1])
        ax[3].axis([1e-3, 1.0, 1e-3, 1.0])
        ax[3].set_yscale('log')
        ax[3].set_xscale('log')
        # Separating 
        for iBase in range(nBasePerGene):
            ax[2].plot([iBase*nGenes,iBase*nGenes],[0,1],'k--')
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
        colrs=['b','r']
        for i,s in enumerate(stats):
            s.plot(y='Best',color=colrs[i],linestyle='-' ,ax=ax[1],legend=False)
            s.plot(y='Min', color=colrs[i],linestyle='--',ax=ax[1],legend=False)
            s.plot(y='Max', color=colrs[i],linestyle='--',ax=ax[1],legend=False)
            s.plot(y='Mean',color=colrs[i],linestyle=':' ,ax=ax[1],legend=False)

    
def individualDistance(ind1,ind2):
    return np.linalg.norm(np.array(ind1)-np.array(ind2))

def populationDiversity(pop):
    return [sum([individualDistance(x,y) for x in pop])/(np.sqrt(2)*(len(pop)-1)) for y in pop]



def getRandomPop(n=3):
    toolbox = base.Toolbox()
    toolbox.register("attr_float" , random.random)
    toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_float, nGenes*nBasePerGene) 
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)
    return toolbox.population(n=n)

def getNeutralChromosome():
    print('Neutral chromosome...')
    toolbox = base.Toolbox()
    toolbox.register("attr_float" , random.random)
    toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_float, nGenes*nBasePerGene) 
    neutral=toolbox.individual()
    neutral[:]=np.matlib.repmat([0.5,0.5,0],1,nGenes)[0].tolist()
    neutral.fitness.values=individualFitness(neutral)
    return neutral


def main(nBase=2,nInd=10,nIndSelect=10,CXPB=0.3,MUTPB=0.3,nIterMax=100,nPerTournament=2,MutStep=0.1,MutMethod='PolyBound',CXFunction=tools.cxTwoPoint, MutParam1=0.05):
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

    neutral_ori=getNeutralChromosome()
    # create an initial population
    pop = toolbox.population(n=nInd)
    if bEnforceNeutral:
        pop.append(toolbox.clone(neutral_ori))
    pop = galib.populationTrimAccuracy(pop,DECIMALS)

    #StatsFits.append = np.zeros((nIterMax,nObjectives)
    
    # --- EVALUATE the entire population
    for ind in pop:
        ind.fitness.values = toolbox.evaluate(ind)
        #ind.fitness.values = [fit[0]*FitScale +  div*DivScale]
    #galib.populationPrint(pop,nBasePerGene,'INITIAL')
    DB_file = galib.populationSave(pop,newfile=True,directory=DB_DIR)
    populationPlot(pop,fig=figs[0],kind='new')
    best = galib.selBestByNorm(pop,weights=objectiveWeights)
    stats,_ = galib.populationStats(pop,best)
    plt.pause(0.001)
    print('---------------------------------------------------')
    # Variable keeping track of the number of generations
    g = 0
    while g < nIterMax:
        g = g + 1
        # --- SUPER INDIVIDUALS
        # first we select the extremes bests
        best = galib.selBestByNorm(pop,weights=objectiveWeights)
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

        offspring = galib.populationTrimAccuracy(offspring,DECIMALS)
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

        # --- STATS
        last_stats,stats = galib.populationStats(pop,best,stats)
        # --- INFO STRING
        sGen  = "Gen:{0:d} ".format(g)
        sInfo = " Xd:{0:d}  Mut:{1:d}  X+M:{2:d}  Pop:{3:d} ".format(nCX, nMute, len(invalid_ind), len(offspring))
        #sStat = " [{0:.1f} {1:.1f}] Avg: {2:.2f} Std {3:.2f} Div: {5:.2f} {6:.2f} Best: {4:.8f}".format(min(fits), max(fits), mean, std,bestFit,np.mean(diversities),np.max(diversities))
        #print(sGen+sInfo+sStat)
        sBest='Best:'+' '.join(['{:.3f}'.format(x) for x in best.fitness.values])+' '
        print(sGen+sInfo+sBest)
        # --- PLOTTING, STORING
        print('Plotting saving....')
        galib.populationSave(offspring,directory=DB_DIR)
        individualPlot(best,fig=figs[1])
        #populationPrint(offspring,'NEW')
        df,pop_uniq=galib.populationLoad(DB_file,nFits=len(pop[0].fitness.values))
        galib.timelinePlot(df,fig=figs[2])
        populationPlot(pop_uniq, stats=stats,fig=figs[0],kind='all',title=sGen+sInfo+sBest)
        populationPlot(parents,fig=figs[0],kind='parents')
        populationPlot(offspring,fig=figs[0],kind='new')
        populationPlot([best],fig=figs[0],kind='best',title=sGen+sBest)
        plt.pause(0.001)

        # --- UPDATE - The population is entirely replaced by the offspring
        pop[:] = offspring
        #time.sleep(1.0)

    if g==nIterMax :
        print(">>>>>>>>>>>>>> Maximum number of iterations reached")

    return pop,best


# --------------------------------------------------------------------------------}
# ---  
# --------------------------------------------------------------------------------{
### --- PARAMETERS
PerformanceSignals = ['CPAero','CTAero']
IPlot              = PerformanceSignals
DECIMALS=2
objectiveWeights=(-1.0,)*len(PerformanceSignals) # negative = minimalization 
nObjectives=len(objectiveWeights)
nGenes= len(airfoilFileNames)
nBasePerGene=3
bEnforceBests=True
bEnforceNeutral=False

FAST=0

if FAST==1:
    GA_DIR  = '_GA_Runs'
    DB_DIR  = '../_GA_Runs_DB'
    ref_dir = '../OpenFAST_V27_refForGA/'
    WS_SIM  = np.array([5.0    ,8    ,12   ])
    TIMEAVGWINDOW=1
else:
    GA_DIR  = '_GA_Runs_AD'
    DB_DIR  = '../_GA_Runs_AD_DB'
    ref_dir = '../AeroDyn_V27_v1_refMulti/'
    #WS_SIM  = np.array([ 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25]);
    WS_SIM  = np.array([ 5, 7, 9, 11, 13, 15, 17, 19]);
    TIMEAVGWINDOW=0.5



OPER=pandalib.pd_interp1('WS',perf_all,WS_SIM)
template_dir = os.path.normpath(os.path.join(ref_dir,'..',GA_DIR,'_Template'))
### --- INIT

creator.create("Fitness", base.Fitness, weights=objectiveWeights)
creator.create("Individual", list, fitness=creator.Fitness)

print('- Preparing template folder')
prepare_template_folder(ref_dir,template_dir,airfoilFileNames,OPER,FAST)
print('- Reading reference polars from template')
airfoils_ref=read_airfoils(airfoilFileNames,workdir=template_dir)
# Interpolating to get the reference values where the simulation will be done
perf_ref=pandalib.pd_interp1('WS',perf_all,WS_SIM)
#print(perf_ref)
plt.ion()
figs=[]
figs+=figlib.fig_grid(AreaName='Left',ScreenName='RightScreen')
figs+=figlib.fig_grid(2,1,AreaName='Right',ScreenName='RightScreen')
plt.show()
##
##pop,pop_init,best_ind=main(nBase=nGenes*nBasePerGene,nInd=10,CXPB=0.5,MUTPB=0.9,nIterMax=100,nPerTournament=2);
##pop,pop_init,best_ind=main(nBase=nGenes*nBasePerGene,nInd=30,CXPB=0.5,MUTPB=0.5,nIterMax=100,nPerTournament=2);
pop,best_ind=main(nBase=nGenes*nBasePerGene
          ,nInd=32,nIndSelect=16,CXPB=0.5,MUTPB=0.5
#          ,MutMethod='PolyBound',MutParam1=0.001,nPerTournament=2
           ,MutMethod='UniBound',MutParam1=np.nan,nPerTournament=2
#           ,MutMethod='GaussianBound',MutParam1=0.01,nPerTournament=2
          ,nIterMax=100);
###
#
#neutral=getNeutralChromosome();
#individualPlot(neutral,fig1=figs[0],fig2=figs[1])
# 




#pop=getRandomPop()
#pop[0].fitness.values=[0.1,0.11]
#pop[1].fitness.values=[0.1,0.005]
#pop[2].fitness.values=[0.09999,0.1]
#best_ind = tools.selBest(pop, 1)[0];
#print(best_ind.fitness.values)

# pdb.set_trace()

#plt.pause(0.01)
#input("press enter to continue")

# 
