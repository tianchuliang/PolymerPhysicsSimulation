# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""



import matplotlib.pyplot as plt
import numpy as np
import sys

filename=sys.argv[1]
#filename2='/users/tianchu/desktop/time-seed-50-out.txt'
f=open(filename,'r')

clns=[]
for line in f:
    line=line.strip()
    column=line.split()
    clns.append(column)
        
info=np.array(clns)
numOfMono=np.zeros(np.size(info[:,0]))
Rg=np.zeros(np.size(info[:,0]))

for x in range(0,np.size(info[:,0])):
    numOfMono[x]=info[x][0]
    Rg[x]=info[x][1]

plt.plot(numOfMono,Rg,'ro')
plt.yscale('log')
plt.xscale('log')
plt.show()

Info2=np.zeros((7,30))

for i in range(0,7):
    for j in range(0,30):
        Info2[i][j]=Rg[i+7*j]
