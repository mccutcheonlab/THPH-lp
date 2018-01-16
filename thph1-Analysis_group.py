# -*- coding: utf-8 -*-
"""
Created on Mon Jul  3 22:35:33 2017

@author: jaimeHP
"""

import numpy as np
import matplotlib.pyplot as plt
import JM_custom_figs as jmfig

data = np.zeros(Rat.nSessions, dtype = [('Rat',int),
                               ('Session',np.str_,20),
                               ('Drug',np.str_,20),
                               ('Licks',float),
                               ('Intake',float),
                               ('Timestamps',list),
                               ('ILIs',list)])

for i in rats:
    for j in rats[i].sessions:
        x = rats[i].sessions[j]
       
        
cumsumBL = []
cumsumDIS=[]
for i in rats:
    j = 's6'
    x = rats[i].sessions[j]
    cumsumBL.append(x.firstlick)
    j = 's7'
    x = rats[i].sessions[j]
    cumsumDIS.append(x.firstlick)
    
groupFig = plt.figure(figsize=(8.27, 11.69), dpi=100)
gs1 = gridspec.GridSpec(5, 2)
gs1.update(left=0.125, right= 0.9, wspace=0.4, hspace = 0.8)
plt.suptitle('Group data')
    
ax = plt.subplot(gs1[0, 0])
for x in cumsumBL:
    jmfig.cumulativelickFig(ax, x, color='grey')
    avg = [item for rat in cumsumBL for item in rat]
    z = jmfig.cumulativelickFig(ax, avg, color='k')
    ax.set_title('Baseline')
    
    
ax = plt.subplot(gs1[0, 1])
for x in cumsumDIS:
    jmfig.cumulativelickFig(ax, x, color='grey')
    avg = [item for rat in cumsumDIS for item in rat]
    z, m = jmfig.cumulativelickFig(ax, avg, color='k')
    ax.set_title('Distraction Day')
        