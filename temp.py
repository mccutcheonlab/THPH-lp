# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 21:52:14 2019

@author: admin
"""
import numpy as np

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
