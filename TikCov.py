#! /usr/bin/python2.7
#----------------------------------------------------------------------------------------
# Name:
#        TickCov.py
#
# Purpose:
#       Create the Covariance for Tikhonov inversion
#
#
#----------------------------------------------------------------------------------------

import numpy as np

nlayer    = 41        #FL0 = 44, MLO = 41, TAB = 47
alpha     = 500
pathOut   = '/data1/ebaumer/mlo/h2co/x.h2co/'

L1=np.zeros((nlayer-1,nlayer),dtype=float)

for i in range(nlayer-1):
    L1[i,i]=1.0
    L1[i,i+1]=-1.0

R=alpha*np.dot(L1.T,L1)

np.savetxt(pathOut+'Tik_'+str(alpha)+'.input',R)