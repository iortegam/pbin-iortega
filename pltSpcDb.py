#!/usr/bin/python
#----------------------------------------------------------------------------------------
# Name:
#        pltSpcDb.py --> 

# Purpose:
#        Read Spec Database and plot
#
# Notes:
#        Initially created to plot trends in temperature
#
# Version History:
#        Created, Jan, 2017  Ivan Ortega (iortega@ucar.edu)
#----------------------------------------------------------------------------------------

import logging
import sys
import os
import getopt
import datetime as dt
import glob
import re
import sfitClasses as sc
import dataOutClass as dc
import shutil
from Layer1Mods import refMkrNCAR, t15ascPrep, errAnalysis    
import matplotlib.pyplot as plt
import itertools
import datetime
import numpy as np
import myfunctions as mf

def main():

    #-----------------------------------------------------------------------------------
    #Input for reading Spec Database
    #-----------------------------------------------------------------------------------
    spcdbFile = '/data/Campaign/TAB/Spectral_DB/CoaddspDB_tab_1999_2014.dat'
    dbData = sc.DbInputFile(spcdbFile)
    dbData.getInputs()

    #-----------------------------------------------------------------------------------
    #Declaring variables to read
    #-----------------------------------------------------------------------------------
    dates = []
    Temp  = []

    for i,(date,time) in enumerate(itertools.izip(dbData.dbInputs['Date'],dbData.dbInputs['Time'])):
        date = str(int(date))
        HH   = time.strip().split(':')[0]
        MM   = time.strip().split(':')[1]
        SS   = time.strip().split(':')[2]
        #-----------------------------------------------------------------------------------
        # Some directories have values of 60 for seconds. This throws an error in date time. 
        # Seconds can only go from 0-59. If this is case reduce seconds by 1. This is a 
        # temp solution until coadd can be fixed.
        #-----------------------------------------------------------------------------------    
        if SS == '60': SS = '59'

        if float(dbData.dbInputs['HouseTemp'][i]) != -9999:

            dates.append(datetime.datetime(int(date[0:4]),int(date[4:6]),int(date[6:]),int(HH),int(MM),int(SS)))
            Temp.append(float(dbData.dbInputs['HouseTemp'][i]))

    #-----------------------------------------------------------------------------------
    #Input for reading 2nd Spec Database
    #-----------------------------------------------------------------------------------
    spcdbFile = '/data/Campaign/TAB/Spectral_DB/HRspDB_tab_2015_2016.dat'
    dbData = sc.DbInputFile(spcdbFile)
    dbData.getInputs()


    for i,(date,time) in enumerate(itertools.izip(dbData.dbInputs['Date'],dbData.dbInputs['Time'])):
        date = str(int(date))
        HH   = time.strip().split(':')[0]
        MM   = time.strip().split(':')[1]
        SS   = time.strip().split(':')[2]
        #-----------------------------------------------------------------------------------
        # Some directories have values of 60 for seconds. This throws an error in date time. 
        # Seconds can only go from 0-59. If this is case reduce seconds by 1. This is a 
        # temp solution until coadd can be fixed.
        #-----------------------------------------------------------------------------------    
        if SS == '60': SS = '59'

        if float(dbData.dbInputs['HouseTemp'][i]) != -9999:

            dates.append(datetime.datetime(int(date[0:4]),int(date[4:6]),int(date[6:]),int(HH),int(MM),int(SS)))
            Temp.append(float(dbData.dbInputs['HouseTemp'][i]))

    dates = np.asarray(dates)
    Temp = np.asarray(Temp)

    #----------------------------
    #Seasonal Temp (by month)
    #----------------------------
    month    = np.array([d.month for d in dates])
    mnthSort = list(set(month))
    mnthMean = np.zeros(len(mnthSort))
    mnthSTD  = np.zeros(len(mnthSort))

    mnthSort = np.asarray(mnthSort)
    
    for i,m in enumerate(mnthSort):
        inds        = np.where(month == m)[0]
        mnthMean[i] = np.mean(Temp[inds])
        mnthSTD[i]  = np.std(Temp[inds])

    mnthVals = mf.mnthlyAvg(Temp,dates,dateAxis=1, meanAxis=0)
    #mnthlyVals['dates'],mnthlyVals['mnthlyAvg']

    #-----------------------------------------------------------------------------------
    #Calculating Temperature anomalies
    #-----------------------------------------------------------------------------------
    TempAnom = []

    for i, m in enumerate(mnthVals['dates']):
        mo = m.month
        indM  = np.where(mnthSort ==mo)[0]
        Anom  = float(mnthVals['mnthlyAvg'][i]) - mnthMean[indM]

        TempAnom.append(Anom)

    TempAnom = np.asarray(TempAnom)

    #-----------------------------------------------------------------------------------
    #Monthly Temp
    #-----------------------------------------------------------------------------------   
    fig,ax1  = plt.subplots(1, figsize=(10,6))
    ax1.plot(mnthSort,mnthMean,'k.',markersize=6)
    ax1.errorbar(mnthSort,mnthMean,yerr=mnthSTD,fmt='k.',markersize=6,ecolor='red')     
    ax1.grid(True,which='both')
    ax1.set_ylabel('Temperature [C]',multialignment='center', fontsize=14)
    ax1.set_xlabel('Month', fontsize=14)
    ax1.set_title('Monthly Mean Temperature at Thule', fontsize=14)
    ax1.set_xlim((0,13))
    ax1.tick_params(labelsize=14)
    ax1.set_xticks(range(1,13))
    
    plt.show(block=False)

    #-----------------------------------------------------------------------------------
    #Temperature anomalies
    #-----------------------------------------------------------------------------------
    ind = np.arange(len(TempAnom)) 
    N = len(ind)
    
    fig , ax = plt.subplots(1, figsize=(10,6))
    
    rects = ax.bar(mnthVals['dates'], TempAnom, color='y')
    ax.set_ylabel('Temperature anomaly [C]', fontsize=14)
    ax.set_title('Temperature anomaly [2000 - 2016]', fontsize=14)
    ax.set_xlabel('Year', fontsize=14)
    #ax.set_xticks(np.arange(0, N, 72))
    ax.tick_params(labelsize=14)
    ax.grid(True)

    plt.show(block=False)

    user_input = raw_input('Press any key to exit >>> ')
    sys.exit()
        


if __name__ == "__main__":
    main()