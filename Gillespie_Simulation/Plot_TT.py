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
filename= '%s'%(sys.argv[1])

data = np.loadtxt(filename)

input = sys.argv[2]

leg= '%s'%(sys.argv[3])

norm='%s'%(sys.argv[4])

if input == 'All':
	ncols= range(1,np.shape(data)[1])

else:
	ncols=map(float, input.strip('[]').split(','))



#print(filename)



b = np.delete(data, 0, axis=1)


if(norm=="norm"):

    bNorm=[]
    for i in b.T:
        i=i/i.max()
        bNorm.append(i)

    b= np.array(bNorm).T

t= data[:,0]
count= 0;
for i in ncols:
	plt.plot(t,b[:,int(i-1)], label = ' Gene_%d'%(count))
	count = count +1;

if leg=='legend':
	plt.legend(loc='best')

plt.xlabel('time', fontsize=15)
plt.ylabel('gene expression', fontsize=15)

plt.show()
