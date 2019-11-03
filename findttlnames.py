# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 20:42:37 2019

@author: admin
"""

from helper_fx import *



folder = 'C:\\Github\\THPH-lp\\data\\'
filename = 'thph1-1_distraction1'

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
    

for file in files:
    output = loadmatfile(file)
    fields=dir(output)
    for fld in fields[:3]:
        print(file, fld, len(getattr(output, fld).onset))
