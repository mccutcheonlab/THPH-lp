# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 18:16:26 2019

@author: admin
"""
import numpy as np
import matplotlib.pyplot as plt
from itertools import chain, count
from collections import Sequence
import matplotlib.cbook as cbook
import matplotlib.mlab as mlab
import scipy.io as sio
import scipy.signal as sig
import xlrd
import csv



def loadmatfile(file):
    a = sio.loadmat(file, squeeze_me=True, struct_as_record=False) 
    output = a['output']
    
    return output
    # moves data streams and fs into more easily accesible variables

def correctforbaseline(blue, uv):
    pt = len(blue)
    X = np.fft.rfft(uv, pt)
    Y = np.fft.rfft(blue, pt)
    Ynet = Y-X

    datafilt = np.fft.irfft(Ynet)

    datafilt = sig.detrend(datafilt)

    b, a = sig.butter(9, 0.012, 'low', analog=True)
    datafilt = sig.filtfilt(b, a, datafilt)
    
    return datafilt

def distractedOrNot(distractors, licks, delay=1):   
    firstlick = []
    distractedArray = []

    for d in distractors:               
        try:
            firstlick.append([i-d for i in licks if (i > d)][0])
        except IndexError:
            firstlick.append(np.NaN)

    distractedArray = np.array([i>delay for i in firstlick], dtype=bool)
    
    if np.isnan(firstlick)[-1] == 1: 
        distractedArray[-1] = True
    
    return firstlick, distractedArray

def snipper(data, timelock, fs = 1, t2sMap = [], preTrial=10, trialLength=30,
                 adjustBaseline = True,
                 bins = 0):

    if len(timelock) == 0:
        print('No events to analyse! Quitting function.')
        raise Exception('no events')

    pps = int(fs) # points per sample
    pre = int(preTrial*pps) 
#    preABS = preTrial
    length = int(trialLength*pps)
# converts events into sample numbers
    event=[]
    if len(t2sMap) > 1:
        for x in timelock:
            event.append(np.searchsorted(t2sMap, x, side="left"))
    else:
        event = [x*fs for x in timelock]

    new_events = []
    for x in event:
        if int(x-pre) > 0:
            new_events.append(x)
    event = new_events

    nSnips = len(event)
    snips = np.empty([nSnips,length])
    avgBaseline = []

    for i, x in enumerate(event):
        start = int(x) - pre
        avgBaseline.append(np.mean(data[start : start + pre]))
        try:
            snips[i] = data[start : start+length]
        except ValueError: # Deals with recording arrays that do not have a full final trial
            snips = snips[:-1]
            avgBaseline = avgBaseline[:-1]
            nSnips = nSnips-1

    if adjustBaseline == True:
        snips = np.subtract(snips.transpose(), avgBaseline).transpose()
        snips = np.divide(snips.transpose(), avgBaseline).transpose()

    if bins > 0:
        if length % bins != 0:
            snips = snips[:,:-(length % bins)]
        totaltime = snips.shape[1] / int(fs)
        snips = np.mean(snips.reshape(nSnips,bins,-1), axis=2)
        pps = bins/totaltime
              
    return snips, pps

def mastersnipper(signal, events, fs, rand,
                  bins=300,
                  preTrial=10,
                  trialLength=30,
                  threshold=10):
    if len(events) < 1:
        print('Cannot find any events. All outputs will be empty.')
        return
    else:
        snips,_ = snipper(signal, events, fs,
                                   bins=bins,
                                   preTrial=preTrial,
                                   trialLength=trialLength,
                                   adjustBaseline = False)
        
        snips_z = np.asarray(zscore(snips))
        
        bgMAD = findnoise(signal, rand, fs=fs, bins=bins,
                              method='sum')        
        sigSum = [np.sum(abs(i)) for i in snips]
        sigSD = [np.std(i) for i in snips]
        noiseindex = [i > bgMAD*threshold for i in sigSum]
        
        peak_z = [np.max(trial[100:110], axis=0) for trial in snips_z]
        trough_z = [np.min(trial[100:150], axis=0) for trial in snips_z]

        output = {}
        output['snips'] = snips
        output['snips_z'] = snips_z
        output['noise'] = noiseindex
        output['peak'] = peak_z
        output['trough'] = trough_z
        
        return output
    
def zscore(snips, baseline_points=100):
    
    BL_range = range(baseline_points)
    z_snips = []
    for i in snips:
        mean = np.mean(i[BL_range])
        sd = np.std(i[BL_range])
        z_snips.append([(x-mean)/sd for x in i])
        
    return z_snips

def findnoise(data, background, t2sMap = [], fs = 1, bins=0, method='sd'):
    
    bgSnips, _ = snipper(data, background, t2sMap=t2sMap, fs=fs, bins=bins)
    
    if method == 'sum':
        bgSum = [np.sum(abs(i)) for i in bgSnips]
        bgMAD = med_abs_dev(bgSum)
    elif method == 'sd':
        bgSD = [np.std(i) for i in bgSnips]
        bgMAD = med_abs_dev(bgSD)
   
    return(bgMAD)

def removenoise(snipsIn, noiseindex):
    snipsOut = np.array([x for (x,v) in zip(snipsIn, noiseindex) if not v])   
    return snipsOut

def med_abs_dev(data, b=1.4826):
    median = np.median(data)
    devs = [abs(i-median) for i in data]
    mad = np.median(devs)*b
                   
    return mad

def makerandomevents(minTime, maxTime, spacing = 77, n=100):
    events = []
    total = maxTime-minTime
    start = 0
    for i in np.arange(0,n):
        if start > total:
            start = start - total
        events.append(start)
        start = start + spacing
    events = [i+minTime for i in events]
    return events

def barscatter(data, transpose = False, unequal=False,
                groupwidth = .75,
                barwidth = .9,
                paired = False,
                spaced = False,
                yspace = 20,
                xspace = 0.1,
                barfacecoloroption = 'same', # other options 'between' or 'individual'
                barfacecolor = ['white'],
                baredgecoloroption = 'same',
                baredgecolor = ['black'],
                baralpha = 1,
                scatterfacecoloroption = 'same',
                scatterfacecolor = ['white'],
                scatteredgecoloroption = 'same',
                scatteredgecolor = ['grey'],
                scatterlinecolor = 'grey', # Don't put this value in a list
                scattersize = 80,
                scatteralpha = 1,
                spreadscatters = False,
                linewidth=1,
                xlim=[],
                ylim=[],
                ylabel = 'none',
                xlabel = 'none',
                grouplabel = 'auto',
                itemlabel = 'none',
                barlabels = [],
                barlabeloffset=0.025,
                grouplabeloffset=0,
                yaxisparams = 'auto',
                show_legend = 'none',
                legendloc='upper right',
                xfontsize=8,
                ax=[]):

    if unequal == True:
        dims = np.ndim(data)
        data_obj = np.ndarray((np.shape(data)), dtype=np.object)
        if dims == 1:
            for i, dim in enumerate(data):
                data_obj[i] = np.array(dim, dtype=np.object)
            data = data_obj
        elif dims == 2:            
            for i1, dim1 in enumerate(data):
                for i2, dim2 in enumerate(dim1):
                    data_obj[i1][i2] = np.array(dim2, dtype=np.object)
            data = data_obj
        else:
            print('Cannot convert that number of dimensions or data is in wrong format. Attmepting to make graph assuming equal groups.')
    
    if type(data) != np.ndarray or data.dtype != np.object:
        dims = np.shape(data)
        if len(dims) == 2 or len(dims) == 1:
            data = data2obj1D(data)

        elif len(dims) == 3:
            data = data2obj2D(data)
              
        else:
            print('Cannot interpret data shape. Should be 2 or 3 dimensional array. Exiting function.')
            return

    # Check if transpose = True
    if transpose == True:
        data = np.transpose(data)
        
    # Initialize arrays and calculate number of groups, bars, items, and means
    
    barMeans = np.zeros((np.shape(data)))
    items = np.zeros((np.shape(data)))
    
    nGroups = np.shape(data)[0]
    groupx = np.arange(1,nGroups+1)

    if len(np.shape(data)) > 1:
        grouped = True
        barspergroup = np.shape(data)[1]
        barwidth = (barwidth * groupwidth) / barspergroup
        
        for i in range(np.shape(data)[0]):
            for j in range(np.shape(data)[1]):
                barMeans[i][j] = np.mean(data[i][j])
                items[i][j] = len(data[i][j])
        
    else:
        grouped = False
        barspergroup = 1
        
        for i in range(np.shape(data)[0]):
            barMeans[i] = np.mean(data[i])
            items[i] = len(data[i])
    
    # Calculate x values for bars and scatters
    
    xvals = np.zeros((np.shape(data)))
    barallocation = groupwidth / barspergroup
    k = (groupwidth/2) - (barallocation/2)
    
    if grouped == True:
        
        for i in range(np.shape(data)[0]):
            xrange = np.linspace(i+1-k, i+1+k, barspergroup)
            for j in range(barspergroup):
                xvals[i][j] = xrange[j]
    else:
        xvals = groupx
    
    # Set colors for bars and scatters
     
    barfacecolorArray = setcolors(barfacecoloroption, barfacecolor, barspergroup, nGroups, data)
    baredgecolorArray = setcolors(baredgecoloroption, baredgecolor, barspergroup, nGroups, data)
     
    scfacecolorArray = setcolors(scatterfacecoloroption, scatterfacecolor, barspergroup, nGroups, data, paired_scatter = paired)
    scedgecolorArray = setcolors(scatteredgecoloroption, scatteredgecolor, barspergroup, nGroups, data, paired_scatter = paired)
    
    # Initialize figure
    if ax == []:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    
    # Make bars
    barlist = []
    barx = []
    for x, y, bfc, bec in zip(xvals.flatten(), barMeans.flatten(),
                              barfacecolorArray, baredgecolorArray):
        barx.append(x)
        barlist.append(ax.bar(x, y, barwidth,
                         facecolor = bfc, edgecolor = bec,
                         zorder=-1))
    
    # Uncomment these lines to show method for changing bar colors outside of
    # function using barlist properties
    #for i in barlist[2].get_children():
    #    i.set_color('r')
    
    # Make scatters
    sclist = []
    if paired == False:
        for x, Yarray, scf, sce  in zip(xvals.flatten(), data.flatten(),
                                        scfacecolorArray, scedgecolorArray):
            if spaced == True:
                xVals, yVals = xyspacer(ax, x, Yarray, bindist=yspace, space=xspace)
                sclist.append(ax.scatter(xVals, yVals, s = scattersize,
                             c = scf,
                             edgecolors = sce,
                             zorder=20))           
            else:
                for y in Yarray:
                     sclist.append(ax.scatter(x, y, s = scattersize,
                                     c = scf,
                                     edgecolors = sce,
                                     zorder=20))

    elif grouped == True:
        for x, Yarray, scf, sce in zip(xvals, data, scfacecolorArray, scedgecolorArray):
            for y in np.transpose(Yarray.tolist()):
                sclist.append(ax.plot(x, y, '-o', markersize = scattersize/10,
                         color = scatterlinecolor,
                         linewidth=linewidth,
                         markerfacecolor = scf,
                         markeredgecolor = sce,
                         zorder=20))
    elif grouped == False:
        for n,_ in enumerate(data[0]):
            y = [y[n-1] for y in data]
            sclist.append(ax.plot(xvals, y, '-o', markersize = scattersize/10,
                         color = scatterlinecolor,
                         linewidth=linewidth,
                         markerfacecolor = scfacecolorArray[0],
                         markeredgecolor = scedgecolorArray[0],
                         zorder=20))
    
    # Label axes
    if ylabel != 'none':
        ax.set_ylabel(ylabel)
    
    if xlabel != 'none':
        ax.set_xlabel(xlabel)
    
    # Set range and tick values for Y axis
    if yaxisparams != 'auto':
        ax.set_ylim(yaxisparams[0])
        ax.set_yticks(yaxisparams[1])
       
    # X ticks
    ax.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom='off',      # ticks along the bottom edge are off
        top='off') # labels along the bottom edge are off

    if len(xlim) > 0:
        ax.set_xlim(xlim)
    if len(ylim) > 0:
        ax.set_ylim(ylim)
        
    xrange = ax.get_xlim()[1] - ax.get_xlim()[0]
    yrange = ax.get_ylim()[1] - ax.get_ylim()[0]
        
    if grouplabel == 'auto':
        ax.tick_params(labelbottom='off')
    else:
        ax.tick_params(labelbottom='off')

        groupx = np.arange(1, len(grouplabel)+1)
        if len(xlim) > 0:
            groupx = [x for x in groupx]
        xpos = (groupx - ax.get_xlim()[0])/xrange

        for x, label in zip(xpos, grouplabel):
            ax.text(x, -0.05+grouplabeloffset, label, va='top', ha='center', fontsize=xfontsize, transform=ax.transAxes)
        
    if len(barlabels) > 0:
        if len(barlabels) != len(barx):
            print('Wrong number of bar labels for number of bars!')
        else:
            xpos = (barx - ax.get_xlim()[0])/xrange
            ypos = (-ax.get_ylim()[0]/yrange) - barlabeloffset

            for x, label in zip(xpos, barlabels):
                ax.text(x, ypos, label, va='top', ha='center', fontsize=xfontsize, transform=ax.transAxes)

    # Hide the right and top spines and set bottom to zero
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_position('zero')
    
    if show_legend == 'within':
        if len(itemlabel) != barspergroup:
            print('Not enough item labels for legend!')
        else:
            legendbar = []
            legendtext = []
            for i in range(barspergroup):
                legendbar.append(barlist[i])
                legendtext.append(itemlabel[i])
            plt.legend(legendbar, legendtext, loc=legendloc)
    
    return ax, barx, barlist, sclist

#plt.savefig('foo.png')
        
# To do
# check if n's are the same for paired and if not default to unpaired
# add color options for scatters
# add alpha options etc
# add axis options
# remove white background
# work out how to export or save as pdf, tiff, eps etc
# work out how to return handles to scatters so can be altered outside of function
# make help doc
# make html file to show usage using ijupyter

      
def setcolors(coloroption, colors, barspergroup, nGroups, data, paired_scatter = False):
            
    nColors = len(colors)
    
    if (paired_scatter == True) & (coloroption == 'within'):
        print('Not possible to make a Paired scatter plot with Within setting.')
        coloroption = 'same'
        
    if coloroption == 'within':
        if nColors < barspergroup:
            print('Not enough colors for this option! Reverting to one color.')
            coloroption = 'same'
        elif nColors > barspergroup:
            colors = colors[:barspergroup]
        coloroutput = [colors for i in data]
        coloroutput = list(chain(*coloroutput))
        
    if coloroption == 'between':
        if nColors < nGroups:
            print('Not enough colors for this option! Reverting to one color.')
            coloroption = 'same'
        elif nColors > nGroups:
            colors = colors[:nGroups]
        if paired_scatter == False:
            coloroutput = [[c]*barspergroup for c in colors]
            coloroutput = list(chain(*coloroutput))
        else:
            coloroutput = colors
            
    if coloroption == 'individual':
        if nColors < nGroups*barspergroup:
            print('Not enough colors for this color option')
            coloroption = 'same'
        elif nColors > nGroups*barspergroup:
            coloroutput = colors[:nGroups*barspergroup]
        else: 
            coloroutput = colors
    
    if coloroption == 'same':
        coloroutput = [colors[0] for x in range(len(data.flatten()))]

    return coloroutput

def data2obj1D(data):
    obj = np.empty(len(data), dtype=np.object)
    for i,x in enumerate(data):
        obj[i] = np.array(x)  
    return obj

def data2obj2D(data):
    obj = np.empty((np.shape(data)[0], np.shape(data)[1]), dtype=np.object)
    for i,x in enumerate(data):
        for j,y in enumerate(x):
            obj[i][j] = np.array(y)
    return obj

def xyspacer(ax, x, yvals, bindist=20, space=0.1):
    
    histrange=[]
    histrange.append(min(ax.get_ylim()[0], min(yvals)))
    histrange.append(max(ax.get_ylim()[1], max(yvals)))

    yhist = np.histogram(yvals, bins=bindist, range=histrange)
    
    xvals=[]
    for ybin in yhist[0]:
        if ybin == 1:
            xvals.append(x)
        elif ybin > 1:          
            temp_vals = np.linspace(x-space, x+space, num=ybin)
            for val in temp_vals:
                xvals.append(val)
                
    yvals.sort()

    return xvals, yvals

def shadedError(ax, yarray, linecolor='black', errorcolor = 'xkcd:silver', linewidth=1):
    yarray = np.array(yarray)
    y = np.mean(yarray, axis=0)
    yerror = np.std(yarray)/np.sqrt(len(yarray))
    x = np.arange(0, len(y))
    ax.plot(x, y, color=linecolor, linewidth=1)
    ax.fill_between(x, y-yerror, y+yerror, color=errorcolor, alpha=0.4)
    
    return ax

def metafilemaker(xlfile, metafilename, sheetname='metafile', fileformat='csv'):
    with xlrd.open_workbook(xlfile) as wb:
        sh = wb.sheet_by_name(sheetname)  # or wb.sheet_by_name('name_of_the_sheet_here')
        
        if fileformat == 'csv':
            with open(metafilename+'.csv', 'w', newline="") as f:
                c = csv.writer(f)
                for r in range(sh.nrows):
                    c.writerow(sh.row_values(r))
        if fileformat == 'txt':
            with open(metafilename+'.txt', 'w', newline="") as f:
                c = csv.writer(f, delimiter="\t")
                for r in range(sh.nrows):
                    c.writerow(sh.row_values(r))
    
def metafilereader(filename):
    
    f = open(filename, 'r')
    f.seek(0)
    header = f.readlines()[0]
    f.seek(0)
    filerows = f.readlines()[1:]
    
    tablerows = []
    
    for i in filerows:
        tablerows.append(i.split('\t'))
        
    header = header.split('\t')
    # need to find a way to strip end of line \n from last column - work-around is to add extra dummy column at end of metafile
    return tablerows, header

def flatten_list(listoflists):
    try:
        flat_list = [item for sublist in listoflists for item in sublist]
        return flat_list
    except:
        print('Cannot flatten list. Maybe is in the wrong format. Returning empty list.')
        return []
    
    