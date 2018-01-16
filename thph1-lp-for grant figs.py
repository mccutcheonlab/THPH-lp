# -*- coding: utf-8 -*-
"""
Created on Sun Jul 23 07:53:20 2017

@author: jaime
"""
import JM_custom_figs as jmfig



x = rats['thph1.4'].sessions['s5']
EOI = x.lickDataL['rStart'][2]

x = rats['thph1.3'].sessions['s5']
EOI = x.lickDataL['rStart'][3]

shortbins = 300
longlength = 300
longbins = 300


triallicks = [val for val in x.licksL if (val>EOI-longlength/2) and (val<EOI+longlength/2)]
triallicksnorm = triallicks - (EOI-longlength/2)


shortsnip, shortpps = jmf.snipper(x.data, [EOI],
                       t2sMap = x.t2sMap, fs = x.fs, bins=shortbins)

shortsnipUV,_ = jmf.snipper(x.dataUV, [EOI],
                       t2sMap = x.t2sMap, fs = x.fs, bins=shortbins)

longsnip, longpps = jmf.snipper(x.data, [EOI],
                       t2sMap = x.t2sMap, fs = x.fs,
                       bins=longbins,
                       preTrial=longlength/2, 
                       trialLength=longlength)

longsnipUV,_ = jmf.snipper(x.dataUV, [EOI],
                       t2sMap = x.t2sMap, fs = x.fs,
                       bins=longbins,
                       preTrial=longlength/2, 
                       trialLength=longlength)


f = plt.figure(figsize=(8.27, 11.69), dpi=100)
gs1 = gridspec.GridSpec(3, 2)
gs1.update(left=0.125, right= 0.9, wspace=0.4, hspace = 0.8)

ax = plt.subplot(gs1[0, :])
#ax.vlines(triallicks, 1, 2)
#ax.set_xlim([EOI-longlength/2, EOI+longlength/2])
ax.hist(triallicksnorm, range(0, longlength, 1))

ax = plt.subplot(gs1[1, :])

jmfig.trialsMultFig(ax, [longsnipUV, longsnip], longpps, preTrial=longlength/2,
                scale=60,
                eventText='lick')

ax = plt.subplot(gs1[2, 1])
jmfig.trialsMultFig(ax, [shortsnipUV, shortsnip], shortpps)

plt.savefig(userhome + '/Dropbox/Python/photometry/output-thph1-lp/' + x.rat + '.eps', format='eps', dpi=1000)