import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import numpy.matlib
import weio

df=weio.read('swiftData_Half.csv').toDataFrame()

# IKeep   = ['sonic32m_ws_sum','sonic32m_ti','alpha','rpmAvg','GenPwrAvg','GenTqAvg','tsrAvg','BldPitch1Avg','flapM_B1Avg','edgeM_B1Avg']
# IRename = ['WS'             ,'TI'         ,'alpha','RPM'   ,'Pgen'     ,'QgenLoss','TSR'   ,'Pitch'       ,'FlapM'      ,'EdgeM']

IKeep   = ['sonic32m_ws_sum','alpha','rpmAvg','GenPwrAvg','GenTqAvg','tsrAvg','BldPitch1Avg','flapM_B1Avg','edgeM_B1Avg']
IRename = ['WS'             ,'alpha','RPM'   ,'Pgen'     ,'QgenLoss','TSR'   ,'Pitch'       ,'FlapM'      ,'EdgeM']


dfSub=df[IKeep].copy()
dfSub.columns=IRename

# --- Binning
wsMet=df['sonic32m_ws_sum'].values
wsbin=np.floor(wsMet*2)/2.
dfSub['WS_bin'] = wsbin
dfBinStd=dfSub.groupby('WS_bin').std()
dfBinStd.columns=[c+'_std' for  c in dfBinStd.columns.values]

dfBin=dfSub.groupby('WS_bin').mean()
print(dfBin)
# --- Save to file
dfBin.to_csv('swiftData_Half_Binned.csv')

# --- Plot scaled by min max
dfBinScaled=dfBin/dfBin.max()
fig=plt.figure(figsize=(12, 8));
ax=fig.add_subplot(111)
dfBinScaled.plot(y=IRename,ax=ax)
ax.plot(df['sonic32m_ws_sum'],df['rpmAvg'],'.')
# ax.plot(dfBin['sonic32m_ws_sum'],dfBin['rpmAvg'])

# --- Plot with std
fig=plt.figure(figsize=(12, 8));
ax=fig.add_subplot(111)
ax.plot(dfBin['WS'],dfBin['RPM'])
ax.plot(dfBin['WS'],dfBin['RPM']+dfBinStd['RPM_std']/2)
ax.plot(dfBin['WS'],dfBin['RPM']-dfBinStd['RPM_std']/2)

fig=plt.figure(figsize=(12, 8));
ax=fig.add_subplot(111)
ax.plot(dfBin['WS'],dfBin['FlapM'])
ax.plot(dfBin['WS'],dfBin['FlapM']+dfBinStd['FlapM_std']/2)
ax.plot(dfBin['WS'],dfBin['FlapM']-dfBinStd['FlapM_std']/2)

fig=plt.figure(figsize=(12, 8));
ax=fig.add_subplot(111)
ax.plot(dfBin['WS'],dfBin['QgenLoss'])
ax.plot(dfBin['WS'],dfBin['QgenLoss']+dfBinStd['QgenLoss_std']/2)
ax.plot(dfBin['WS'],dfBin['QgenLoss']-dfBinStd['QgenLoss_std']/2)


plt.show()


