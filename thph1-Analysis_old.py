# -*- coding: utf-8 -*-
"""
Created on Wed May 24 11:09:00 2017

MCH6 analysis


To change file set analysed (i.e. blue vs violet vs processed change the following lines:
    in Class Session.__init__
        self.matlabfile = datafolder + self.hrow['rat'] + 's' + self.hrow['session'] + '-p.mat')
    in first input declaration, after classes/functions
        exptsuffix = '-p'
        
@author: James Rig
"""
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import JM_general_functions as jmf
import JM_custom_figs as jmfig

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec

plt.style.use('seaborn-muted')

import os

userhome = os.path.expanduser('~')
datafolder = userhome + '\\Dropbox\\Python\\matlab files\\'

class Rat(object):
    
    nRats = 0
    nSessions = 0
    
    def __init__(self, name):      
        self.rat = name
        self.sessions = {}
        
        Rat.nRats += 1
                
    def loadsession(self, data, header):
        self.session = 's'+str(data[4]) #should reference column of data with session number
        self.sessions[self.session] = Session(data, header, self.rat, self.session)
       
        Rat.nSessions += 1
        
class Session(object):
    
    def __init__(self, data, header, rat, session):
        self.hrow = {}
        for idx, col in enumerate(header):
            self.hrow[col] = data[idx]
        self.matlabfile = datafolder + self.hrow['rat'] + 's' + self.hrow['session'] + '.mat'
        self.medfile = datafolder + self.hrow['medfile']
        self.bottleA = self.hrow['bottleA']
        self.bottleB = self.hrow['bottleB']
        self.rat = str(rat)
        self.session = session
        
        self.bottles = {}
#        self.medfile = currentpath + 'MATLAB\\Experiments\\2015_casein\\behavfiles\\' + data[1]
#        self.lickdata = jmf.medfilereader(self.medfile, ['a', 'b', 'c', 'd'])
    def loadmatfile(self):
        a = sio.loadmat(self.matlabfile, squeeze_me=True, struct_as_record=False) 
        self.output = a['output']
        self.fs = self.output.fs
        self.data = self.output.result
        self.dataUV = self.output.resultUV
        
    def loadUVfile(self):
        a = sio.loadmat(self.uvfile, squeeze_me=True, struct_as_record=False)
        b = a['output']
        self.dataUV = b.result
        
    def time2samples(self):
        tick = self.output.Tick.onset
        maxsamples = len(tick)*int(self.fs)
        if (len(self.data) - maxsamples) > 2*int(self.fs):
            print('Something may be wrong with conversion from time to samples')
            print(str(len(self.data) - maxsamples) + ' samples left over. This is more than double fs.')
            self.t2sMap = np.linspace(min(tick), max(tick), maxsamples)
        else:
            self.t2sMap = np.linspace(min(tick), max(tick), maxsamples)
            
#    def event2sample(self, EOI):
#        idx = (np.abs(self.t2sMap - EOI)).argmin()   
#        return idx
    
    def check4events(self):
        
        if hasattr(self.output, 'LiA'):
            self.licksA = self.output.LiA.onset
            self.offsetA = self.output.LiA.offset
        else:
            self.licksA = []
            self.offsetA = []
            
        if hasattr(self.output, 'LiB'):
            self.licksB = self.output.LiB.onset
            self.offsetB = self.output.LiB.offset
        else:
            self.licksB = []
            self.offsetB = []
            
    def sessionFig(self, ax):
        ax.plot(self.data)
        ax.set_xticks(np.multiply([0, 10, 20, 30, 40, 50, 60],60*self.fs))
        ax.set_xticklabels(['0', '10', '20', '30', '40', '50', '60'])
        ax.set_xlabel('Time (min)')
        ax.set_title('Rat ' + self.rat + ': Session ' + self.session)
        

    def sessionLicksFig(self, ax):
        try:
            ax.vlines(self.lickDataA['licks'], 0.5, 1.5)
        except AttributeError:
            pass
        
        try:
            ax.vlines(self.lickDataB['licks'], 2.5, 3.5)
        except AttributeError:
            pass
        
        ax.set_xlim(0,3600)
        ax.set_ylim(0,5)
        
        ax.set_yticks([])
        ax.set_xlabel('Time (s)')
        ax.set_title('Rat ' + self.rat + ': Session ' + self.session)
        
        ax.text(100, 2, self.bottleA, ha='left',va='center')
        ax.text(100, 4, self.bottleB, ha='left',va='center')
    
    def findlicks1bottle(self):
        if self.bottleB == 'empty':
            self.licks = self.output.LiA.onset
            self.offset = self.output.LiA.offset
        elif self.bottleA == 'empty':
            self.licks = self.output.LiB.onset
            self.offset = self.output.LiB.offset
        else:
            print('No empty bottle. Unable to assign licks. Please update metafile.')

def makeBehavFigs(x):
    # Initialize figure
    behavFig = plt.figure(figsize=(8.27, 11.69), dpi=100)
    gs1 = gridspec.GridSpec(5, 2)
    gs1.update(left=0.10, right= 0.9, wspace=0.5, hspace = 0.7)
    
    # Make figures of lick data from each bottle
    ax = plt.subplot(gs1[0, :])
    x.sessionLicksFig(ax)        
    
    
    if len(x.lickDataA['licks']) > 1:
        ax = plt.subplot(gs1[1, 0])
        jmfig.licklengthFig(ax, x.lickDataA, contents = x.bottleA)
        
        ax = plt.subplot(gs1[2, 0])
        jmfig.iliFig(ax, x.lickDataA, contents = x.bottleA)
        
        ax = plt.subplot(gs1[3, 0])
        jmfig.burstlengthFig(ax, x.lickDataA, contents = x.bottleA)

        ax = plt.subplot(gs1[4, 0])
        jmfig.ibiFig(ax, x.lickDataA, contents = x.bottleA)
    else:
        print('No licks for Bottle A in this session.')
    
    #plt.tight_layout()
    pdf_pages.savefig(behavFig)

def makePhotoFigs(x):
    # Initialize photometry figure
    photoFig = plt.figure(figsize=(8.27, 11.69), dpi=100)
    gs1 = gridspec.GridSpec(4, 2)
    gs1.update(left=0.125, right= 0.9, wspace=0.4, hspace = 0.8)

    ax = plt.subplot(gs1[0, :])
    x.sessionFig(ax)
        
    if len(x.lickDataA['licks']) > 1: 
        ax = plt.subplot(gs1[2, 0])
        jmfig.trialsFig(ax, x.LiATrials, x.pps, eventText = 'LickA')
        
        ax = plt.subplot(gs1[2, 1])
        jmfig.trialsMultShadedFig(ax, [x.LiATrials, x.LiAUVTrials],
                                  x.pps, eventText = 'LickA')

#    plt.savefig(x.rat + '.eps', format='eps', dpi=1000)
    pdf_pages.savefig(photoFig)

def makeHeatmapFigs(x):
    
    heatmapsFig = plt.figure(figsize=(8.27, 11.69), dpi=100)
    gs1 = gridspec.GridSpec(2, 1)
    gs1.update(left=0.125, right= 0.9, wspace=0.4, hspace = 0.8)
    
    ax = plt.subplot(gs1[0, 0])
    jmfig.heatmapFig(ax, x.LiATrials, pps = x.pps)
    data = jmf.nearestevents(x.lickDataA['rStart'], x.lickDataA['licks'])
    nLicks = [len(i) for i in data]
    sortOrder = np.argsort(nLicks)

    ax = plt.subplot(gs1[1, 0])
    jmfig.heatmapFig(ax, x.LiATrials, pps = x.pps, sortEvs = nLicks)
    
    sortedData = []
    for val in sortOrder:
        sortedData.append(data[val])
    
    jmfig.addevent2heatmap(ax, sortedData, pps = x.pps)
   
    pdf_pages.savefig(heatmapsFig)

# Read in metafile

metafileData, metafileHeader = jmf.metafilereader(userhome + '/Dropbox/Python/photometry/thph1-forMatPy.txt')
exptsuffix = ''
includecol = 5

rats = {}

for i in metafileData:
    if int(i[includecol]) == 1:
        rowrat = str(i[2])
        if rowrat not in rats:
            rats[rowrat] = Rat(rowrat)
        rats[rowrat].loadsession(i, metafileHeader)

#for i in ['thph1.3']:
#    pdf_pages = PdfPages(i + '.pdf')
#    for j in ['s1']:
    
for i in rats:
    pdf_pages = PdfPages(userhome + '/Dropbox/Python/photometry/output/' + i + exptsuffix + '.pdf')
    for j in rats[i].sessions:
        
        print('Analysing rat ' + i + ' in session ' + j)
    
        # Load in data from .mat file (convert from Tank first using Matlab script)
        x = rats[i].sessions[j]
        x.loadmatfile()
        # Work out time to samples
        x.time2samples()       
        # Find out which bottles have TTLs/Licks associated with them
        x.check4events()
        x.lickDataA = jmf.lickCalc(x.licksA,
                           offset = x.offsetA)
# Only Bottle A used in this experiment        
#        x.lickDataB = jmf.lickCalc(x.licksB,
#                   offset = x.offsetB)
     
        try:
            x.siTrials, x.pps = jmf.snipper(x.data, x.output.Sir.onset,
                                        t2sMap = x.t2sMap, fs = x.fs, bins=200)   
            x.siUVTrials, x.pps = jmf.snipper(x.dataUV, x.output.Sir.onset,
                                        t2sMap = x.t2sMap, fs = x.fs, bins=200)   
        except TypeError:
            print('Sir event is not present or not an array.')
        
        if len(x.lickDataA['rStart']) > 0:
            x.LiATrials, x.pps = jmf.snipper(x.data, x.lickDataA['rStart'], 
                                        t2sMap = x.t2sMap, fs = x.fs, bins=300)
            x.LiAUVTrials, x.pps = jmf.snipper(x.dataUV, x.lickDataA['rStart'], 
                                        t2sMap = x.t2sMap, fs = x.fs, bins=300)
        
        makeBehavFigs(x)       
        makePhotoFigs(x)
        makeHeatmapFigs(x)

    pdf_pages.close()
"""
test = rats['mch6.6'].sessions['s4'].bottles['BottleB'].iliFig
rats['mch6.6'].sessions['s5'].hrow
http://python-guide-pt-br.readthedocs.io/en/latest/writing/style/

to add:
    histograms on session licks - 4
    individual points on interburstintervals - 5
    colours for different solutions - 3
    loop for behav figs - 2
    photometry data split for different events and coloured - 1
    heatmap for photometry
    lick checker (e.g. remove very very short ilis)
"""

