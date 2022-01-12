#!/usr/bin/env python
# coding: utf-8

import numpy as np
from matplotlib import pyplot as plt
import sys
'''
print('Number of arguments:', len(sys.argv), 'arguments')

print('Argument List:', str(sys.argv))

print((sys.argv[1]))
print(sys.argv[0])
'''
N= int(sys.argv[1])




nr =  int(sys.argv[2])

na =  int(sys.argv[3])


filename= 'run/Protein_N_%d'%N + '_nr%d'%nr +'_na%d'%na +'.txt'
filename2= 'run/gene_N_%d'%N + '_nr%d'%nr +'_na%d'%na +'.txt'


data = np.loadtxt(filename)
data1 = np.loadtxt(filename2)

inp = sys.argv[2]

ncols= range(1,np.shape(data)[1])


b = np.delete(data, 0, axis=1)

t= data[:,0]

genes =np.transpose(np.delete(data1, 0, axis=1))



f, ax = plt.subplots(2)
count= 0;
for i in ncols:
    ax[0].plot(t,b[:,int(i-1)], label = ' Gene_%d'%(count))
    count = count +1;


ax[1].set_xlabel('time', fontsize=15)
ax[1].set_ylabel('Gene expression', fontsize=15)

ax[0].set_ylabel('Protein concentrations', fontsize=15)


ax[1].imshow(genes, interpolation ='None', aspect = 'auto')


plt.show()
