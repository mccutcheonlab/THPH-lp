# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 16:36:32 2019

@author: admin
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
import scipy.signal as sig

from scipy import stats

from helper_fx import *

folder = 'C:\\Github\\THPH-lp\\data\\'

#Assemble data
metafilemaker(folder+'simplemetafile.xlsx', folder+'metafile', fileformat='txt')
rows, header=metafilereader(folder+'metafile.txt')

rats=[]
fileArray=[]
lickArray=[]
disArray=[]
for row in rows:
    rats.append(row[0])
    fileArray.append(folder + row[0] + '_distraction')
    lickArray.append(row[1])
    disArray.append(row[2])

change_pdps=[]
nochange_pdps=[]

for idx, file in enumerate(fileArray):
    try:
        output = loadmatfile(file)
        
        fs = output.fs
        data = output.blue
        dataUV = output.uv
        data_filt = correctforbaseline(data, dataUV)
        
        licks = getattr(output, lickArray[idx]).onset
        dis = getattr(output, disArray[idx]).onset
        
        randomevents = makerandomevents(0, max(dis))
        
        pdps = []
        for d in dis:
            try:
                pdp = [l-d for l in licks if l>d][0]
                pdps.append(pdp)
            except IndexError: pass
        
        dis = dis[:len(pdps)] # removes last distractor if not matching pdp
        
        snips = mastersnipper(data_filt, dis, fs, randomevents)
        
        
        
        cutoff=2.3
        avg_change_trials = [trial for trial, peak in zip(snips['snips_z'], snips['peak']) if peak>cutoff]
        avg_nochange_trials = [trial for trial, peak in zip(snips['snips_z'], snips['peak']) if peak<cutoff]
        
        avg_change_pdp = [pdp for peak, pdp in zip(snips['peak'], pdps) if peak>cutoff]
        avg_nochange_pdp = [pdp for peak, pdp in zip(snips['peak'], pdps) if peak<cutoff]
        
        change_pdps.append(avg_change_pdp)
        nochange_pdps.append(avg_nochange_pdp)
        
        #avg_dist_trough = [trough for trough, pdp in zip(snips['trough'], pdps) if pdp>cutoff]
        #avg_nd_trough = [trough for trough, pdp in zip(snips['trough'], pdps) if pdp<cutoff]
        
        f, ax = plt.subplots(figsize=(10,3), ncols=3, gridspec_kw = {'width_ratios':[1,0.5, 0.5]})
        f.subplots_adjust(wspace=0.4)
        shadedError(ax[0], avg_change_trials, linecolor='red')
        shadedError(ax[0], avg_nochange_trials, linecolor='blue')
        
        ax[0].set_ylabel('Z-Score')
        
        barscatter([avg_change_pdp, avg_nochange_pdp], ax=ax[1],
                   barfacecoloroption = 'individual',
                   barfacecolor = ['red', 'blue'])
        
        ax[2].bar([1,2], [np.median(avg_change_pdp), np.median(avg_nochange_pdp)] )
        
        dis_stats = stats.ttest_ind(avg_change_pdp, avg_nochange_pdp)
        print(dis_stats)
    except:
        print('Trouble with rat', rats[idx])


change_pdps = flatten_list(change_pdps)
nochange_pdps = flatten_list(nochange_pdps)

def calculate_pdp_prob(pdps):
    xdata = np.sort(pdps)
    ydata = [1-i/len(xdata) for i,val in enumerate(xdata)]
    
    return xdata, ydata

xdata, ydata = calculate_pdp_prob(change_pdps)
x2, y2 = calculate_pdp_prob(nochange_pdps)

f, ax = plt.subplots()
ax.plot(xdata, ydata, 'r')
ax.set_xscale('log')
ax.plot(x2, y2, 'b')


# Code to look at whether PDP predicts peak 

#cutoff=1
#
#avg_dist_trials = [trial for trial, pdp in zip(snips['snips_z'], pdps) if pdp>cutoff]
#avg_nd_trials = [trial for trial, pdp in zip(snips['snips_z'], pdps) if pdp<cutoff]
#
#avg_dist_peak = [peak for peak, pdp in zip(snips['peak'], pdps) if pdp>cutoff]
#avg_nd_peak = [peak for peak, pdp in zip(snips['peak'], pdps) if pdp<cutoff]
#
#avg_dist_trough = [trough for trough, pdp in zip(snips['trough'], pdps) if pdp>cutoff]
#avg_nd_trough = [trough for trough, pdp in zip(snips['trough'], pdps) if pdp<cutoff]
#
#if len(avg_dist_peak)+len(avg_nd_peak) != len(pdps):
#    print('Check code.')
#    
#f, ax = plt.subplots(figsize=(10,3), ncols=3, gridspec_kw = {'width_ratios':[1,0.5, 0.5]})
#f.subplots_adjust(wspace=0.4)
#shadedError(ax[0], avg_dist_trials, linecolor='red')
#shadedError(ax[0], avg_nd_trials, linecolor='blue')
#
#ax[0].set_ylabel('Z-Score')
#
#barscatter([avg_dist_peak, avg_nd_peak], ax=ax[1],
#           barfacecoloroption = 'individual',
#           barfacecolor = ['red', 'blue'])
#
#barscatter([avg_dist_trough, avg_nd_trough], ax=ax[2],
#           barfacecoloroption = 'individual',
#           barfacecolor = ['red', 'blue'])
#
#dis_stats = stats.ttest_ind(avg_dist_peak, avg_nd_peak)
#print(dis_stats)
#
#dis_stats = stats.ttest_ind(avg_dist_trough, avg_nd_trough)
#print(dis_stats)