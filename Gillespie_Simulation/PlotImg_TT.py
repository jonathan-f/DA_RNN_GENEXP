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
filename= 'noisy_time_traces' + '_sig%s'%(sys.argv[1]) + '.txt'
filename= '%s'%(sys.argv[1])

norm = '%s'%(sys.argv[2])


#print(filename)

data = np.loadtxt(filename)

b = np.delete(data, 0, axis=1)

if(norm=="norm"):
    bNorm=[]
    for i in b.T:
        i=i/i.max()
        bNorm.append(i)

    b= np.array(bNorm).T



plt.imshow(b.T, aspect= 'auto', interpolation=None)

plt.xlabel('time', fontsize=15)
plt.ylabel('gene expression', fontsize=15)

plt.show()
