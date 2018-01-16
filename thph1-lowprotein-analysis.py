# -*- coding: utf-8 -*-
"""
Created on Tue Jul 11 16:03:02 2017

THPH1 Analysis for LOW PROTEIN expt

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
import timeit

tic = timeit.default_timer()

datafolder = 'R:/DA_and_Reward/jem64/Random Datafiles/'

class Rat(object):
    
    nRats = 0
    nSessions = 0
    
    def __init__(self, data):      
        self.rat = data
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
        self.uvfile = datafolder + self.hrow['rat'] + 's' + self.hrow['session'] + '-v.mat'
        self.medfile = datafolder + self.hrow['medfile']
        self.bottleL = self.hrow['bottleA']
        self.bottleR = self.hrow['bottleB']
        self.rat = str(rat)
        self.session = session
        
        self.bottles = {}

    def loadmatfile(self):
        a = sio.loadmat(self.matlabfile, squeeze_me=True, struct_as_record=False) 
        self.output = a['output']
        self.fs = self.output.fs
        self.data = self.output.result
        self.dataUV = self.output.resultUV
        
    def time2samples(self):
        tick = self.output.Tick.onset
        maxsamples = len(tick)*int(self.fs)
        if (len(self.data) - maxsamples) > 2*int(self.fs):
            print('Something may be wrong with conversion from time to samples')
            print(str(len(self.data) - maxsamples) + ' samples left over. This is more than double fs.')
            self.t2sMap = np.linspace(min(tick), max(tick), maxsamples)
        else:
            self.t2sMap = np.linspace(min(tick), max(tick), maxsamples)
            
    def event2sample(self, EOI):
        idx = (np.abs(self.t2sMap - EOI)).argmin()   
        return idx
    
    def check4events(self):        
        if hasattr(self.output, 'LTl'):
            self.leftTrials = True
            self.cuesL = self.output.LTl.onset
            self.cuesL_off = self.output.LTl.offset
            self.licksL = np.array([i for i in self.output.LLk.onset if i<max(self.cuesL_off)])
            self.offsetL = self.output.LLk.offset[:len(self.licksL)]
        else:
            self.leftTrials = False
            self.cuesL = []
            self.cuesL_off = []
            self.licksL = []
            self.offsetL = []
            
        if hasattr(self.output, 'RTl'):
            self.rightTrials = True
            self.cuesR = self.output.RTl.onset
            self.cuesR_off = self.output.RTl.offset
            self.licksR = np.array([i for i in self.output.RLk.onset if i<max(self.cuesR_off)])
            self.offsetR = self.output.RLk.offset[:len(self.licksR)]
        else:
            self.rightTrials = False
            self.cuesR = []
            self.cuesR_off = []
            self.licksR = []
            self.offsetR = []
            
    def removephantomlicks(self):
        if self.leftTrials == True:
            phlicks = jmf.findphantomlicks(self.licksL, self.cuesL, delay=3)
            self.licksL = np.delete(self.licksL, phlicks)
            self.offsetL = np.delete(self.offsetL, phlicks)
    
        if self.rightTrials == True:
            phlicks = jmf.findphantomlicks(self.licksR, self.cuesR, delay=3)
            self.licksR = np.delete(self.licksR, phlicks)
            self.offsetR = np.delete(self.offsetR, phlicks)
                        
    def sessionFig(self, ax):
        ax.plot(self.data, color='blue')
        try:
            ax.plot(self.dataUV, color='m')
        except:
            print('No UV data.')
        ax.set_xticks(np.multiply([0, 10, 20, 30, 40, 50, 60],60*self.fs))
        ax.set_xticklabels(['0', '10', '20', '30', '40', '50', '60'])
        ax.set_xlabel('Time (min)')
        ax.set_title('Rat ' + self.rat + ': Session ' + self.session)
        
    def makephotoTrials(self, bins, events, threshold=10):
        bgMAD = jmf.findnoise(self.data, self.randomevents,
                              t2sMap = self.t2sMap, fs = self.fs, bins=bins,
                              method='sum')          
        blueTrials, self.pps = jmf.snipper(self.data, events,
                                            t2sMap = self.t2sMap, fs = self.fs, bins=bins)        
        UVTrials, self.pps = jmf.snipper(self.dataUV, events,
                                            t2sMap = self.t2sMap, fs = self.fs, bins=bins)
        sigSum = [np.sum(abs(i)) for i in blueTrials]
        sigSD = [np.std(i) for i in blueTrials]
        noiseindex = [i > bgMAD*threshold for i in sigSum]

        return blueTrials, UVTrials, noiseindex
    
        
def makeBehavFigs(x):
    # Initialize figure
    behavFig = plt.figure(figsize=(8.27, 11.69), dpi=100)
    gs1 = gridspec.GridSpec(5, 2)
    gs1.update(left=0.10, right= 0.9, wspace=0.5, hspace = 0.7)
    plt.suptitle('Rat ' + x.rat + ': Session ' + x.session)
    
#    ax = plt.subplot(gs1[0, :])
#    jmfig.sessionlicksFig(ax, x.lickData)

    if x.leftTrials == True:
        behavFigsCol(gs1, 0, x.lickDataL, x.cuesL)
        
    if x.rightTrials == True:
        behavFigsCol(gs1, 1, x.lickDataR, x.cuesR)
        
    ax = plt.subplot(gs1[4, 0])
    jmfig.latencyFig(ax, x)

    pdf_pages.savefig(behavFig)

def behavFigsCol(gs1, col, lickdata, cues):
    ax = plt.subplot(gs1[1, col])
    jmfig.licklengthFig(ax, lickdata)
    
    ax = plt.subplot(gs1[2, col])
    jmfig.iliFig(ax, lickdata)
    
    ax = plt.subplot(gs1[3, col])
    jmfig.cuerasterFig(ax, cues, lickdata['licks'])
    
def makePhotoFigs(x):
    # Initialize photometry figure
    photoFig = plt.figure(figsize=(8.27, 11.69), dpi=100)
    gs1 = gridspec.GridSpec(6, 2)
    gs1.update(left=0.125, right= 0.9, wspace=0.4, hspace = 0.8)
    plt.suptitle('Rat ' + x.rat + ': Session ' + x.session)

    ax = plt.subplot(gs1[0, :])
    x.sessionFig(ax)

    if x.leftTrials == True:
        photoFigsCol(gs1, 0, x.pps,
                     x.cueLTrials, x.cueLUVTrials, x.cueLnoise,
                     x.lickLTrials, x.lickLUVTrials, x.lickLnoise)

    if x.rightTrials == True:
        photoFigsCol(gs1, 1, x.pps,
                     x.cueRTrials, x.cueRUVTrials, x.cueRnoise,
                     x.lickRTrials, x.lickRUVTrials, x.lickRnoise)
        
    if x.leftTrials == True and x.rightTrials == True:
        diffcueL = jmf.findphotodiff(x.cueLTrials, x.cueLUVTrials, x.cueLnoise)
        diffcueR = jmf.findphotodiff(x.cueRTrials, x.cueRUVTrials, x.cueRnoise)

        ax = plt.subplot(gs1[5, 0])
        jmfig.trialsMultShadedFig(ax, [diffcueL, diffcueR], x.pps,
                                  linecolor=['orange', 'g'], eventText = 'Cue')

        difflickL = jmf.findphotodiff(x.lickLTrials, x.lickLUVTrials, x.lickLnoise)
        difflickR = jmf.findphotodiff(x.lickRTrials, x.lickRUVTrials, x.lickRnoise)

        ax = plt.subplot(gs1[5, 1])
        jmfig.trialsMultShadedFig(ax, [difflickL, difflickR], x.pps,
                                  linecolor=['orange', 'g'], eventText = 'Lick')
        
    plt.savefig('R:/DA_and_Reward/jem64/1706_THPH1-lp/output/' + x.rat + '.eps', format='eps', dpi=1000)
    pdf_pages.savefig(photoFig)
    
def photoFigsCol(gs1, col, pps, cues, cuesUV, cuesnoise, licks, licksUV, licksnoise):
    ax = plt.subplot(gs1[1, col])
    jmfig.trialsFig(ax, cues, pps, noiseindex = cuesnoise,
                    eventText = 'Cue',
                    ylabel = 'Delta F / F0')
    
    ax = plt.subplot(gs1[2, col])
    jmfig.trialsMultShadedFig(ax, [cuesUV, cues], pps, noiseindex = cuesnoise,
                              eventText = 'Cue')
    
    ax = plt.subplot(gs1[3, col])
    jmfig.trialsFig(ax, licks, pps, noiseindex = licksnoise,
                    eventText = 'First Lick',
                    ylabel = 'Delta F / F0')
    
    ax = plt.subplot(gs1[4, col])
    jmfig.trialsMultShadedFig(ax, [licksUV, licks], pps, noiseindex = licksnoise,
                              eventText = 'First Lick')

def makeGrantFigs(x):
    if x.leftTrials == True and x.rightTrials == True:       
        
        grantLFig = plt.figure(figsize=(8.27, 11.69), dpi=100)
        gs1 = gridspec.GridSpec(9, 7)
        gs1.update(left=0.125, right= 0.9, wspace=0.4, hspace = 0.8)
        plt.suptitle('Rat ' + x.rat + ': Session ' + x.session)
        jmfig.trialstiledFig(gs1, x.lickLTrials)
        pdf_pages.savefig(grantLFig)
        
        grantRFig = plt.figure(figsize=(8.27, 11.69), dpi=100)
        gs1.update(left=0.125, right= 0.9, wspace=0.4, hspace = 0.8)
        plt.suptitle('Rat ' + x.rat + ': Session ' + x.session)
        jmfig.trialstiledFig(gs1, x.lickRTrials)
        pdf_pages.savefig(grantRFig)
    
    elif x.rightTrials == False:
        grantLFig = plt.figure(figsize=(8.27, 11.69), dpi=100)
        gs1 = gridspec.GridSpec(9, 7)
        gs1.update(left=0.125, right= 0.9, wspace=0.4, hspace = 0.8)
        plt.suptitle('Rat ' + x.rat + ': Session ' + x.session)
        jmfig.trialstiledFig(gs1, x.lickLTrials)
        pdf_pages.savefig(grantLFig)
        
    elif x.leftTrials == False:
        grantRFig = plt.figure(figsize=(8.27, 11.69), dpi=100)
        gs1 = gridspec.GridSpec(9, 7)
        gs1.update(left=0.125, right= 0.9, wspace=0.4, hspace = 0.8)
        plt.suptitle('Rat ' + x.rat + ': Session ' + x.session)
        jmfig.trialstiledFig(gs1, x.lickRTrials)
        pdf_pages.savefig(grantRFig)
        

metafile = 'R:/DA_and_Reward/jem64/1706_THPH1-lp/thph1-lowprotein-metafile.txt'
metafileData, metafileHeader = jmf.metafilereader(metafile)

exptsuffix = '-lp'
includecol = 5

rats = {}

for i in metafileData:
    if int(i[includecol]) == 1:
        rowrat = str(i[2])
        if rowrat not in rats:
            rats[rowrat] = Rat(rowrat)
        rats[rowrat].loadsession(i, metafileHeader)
        
for i in ['thph1.5']:
    pdf_pages = PdfPages('R:/DA_and_Reward/jem64/1706_THPH1-lp/output/' + i + exptsuffix + '.pdf')
    for j in ['s4']:
        print('\nAnalysing rat ' + i + ' in session ' + j)
        
        # Load in data from .mat file (convert from Tank first using Matlab script)
        x = rats[i].sessions[j]
        x.loadmatfile()
        # Work out time to samples
        x.time2samples()       
        # Find out which bottles have TTLs/Licks associated with them     
        x.check4events()
        x.removephantomlicks()
        
        x.lickDataL = jmf.lickCalc(x.licksL,
                          offset = x.offsetL,
                          burstThreshold = 0.50)
        
        x.lickDataR = jmf.lickCalc(x.licksR,
                  offset = x.offsetR,
                  burstThreshold = 0.50)
        
        bins = 300
        
        x.randomevents = jmf.makerandomevents(120, max(x.output.Tick.onset)-120)
        x.bgTrials, x.pps = jmf.snipper(x.data, x.randomevents,
                                        t2sMap = x.t2sMap, fs = x.fs, bins=bins)

        if x.leftTrials == True:
            x.cueLTrials, x.cueLUVTrials, x.cueLnoise = x.makephotoTrials(bins, x.cuesL)
            x.lickLTrials, x.lickLUVTrials, x.lickLnoise = x.makephotoTrials(bins, x.lickDataL['rStart'])
            x.latsL = jmf.latencyCalc(x.lickDataL['licks'], x.cuesL)
        
        if x.rightTrials == True:
            x.cueRTrials, x.cueRUVTrials, x.cueRnoise = x.makephotoTrials(bins, x.cuesR)
            x.lickRTrials, x.lickRUVTrials, x.lickRnoise = x.makephotoTrials(bins, x.lickDataR['rStart'])
            x.latsR = jmf.latencyCalc(x.lickDataR['licks'], x.cuesR)


#        makeBehavFigs(x)
        makePhotoFigs(x)
#        makeGrantFigs(x)
#        
    pdf_pages.close()
    plt.close('all')
    
