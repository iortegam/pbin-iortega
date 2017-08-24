#! /usr/bin/python2.7
#----------------------------------------------------------------------------------------
# Name:
#        pltfts_vs_insitu.py
#
# Purpose:
#       The purpose of this program is to plt FTS versus insitu data (initially developed using ICARTT format)
#
# Notes:
#   
#
# Version History:
#       Created, July, 2016  Ivan Ortega (iortega@ucar.edu)

#----------------------------------------------------------------------------------------

                                    #-------------------------#
                                    # Import Standard modules #
                                    #-------------------------#
from scipy.io import netcdf
import os
import datetime as dt
import numpy as np
from numpy import fromstring, vectorize, ndarray, array, genfromtxt
import sys
import glob
import getopt

import matplotlib.dates as md
import matplotlib.dates as md
from matplotlib.dates import DateFormatter, MonthLocator, YearLocator

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.dates import DateFormatter, MonthLocator, YearLocator
from matplotlib.ticker import FormatStrFormatter, MultipleLocator,AutoMinorLocator,ScalarFormatter
from matplotlib.backends.backend_pdf import PdfPages #to save multiple pages in 1 pdf...
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as mplcm
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
import matplotlib.gridspec as gridspec
from itertools import izip
from numpy import *
import myfunctions as mf
from datetime import datetime, timedelta

import dataOutplts as dc

                                    #-------------------------#
                                    # Define helper functions #
                                    #-------------------------#

def ckDir(dirName,logFlg=False,exit=False):
    ''' '''
    if not os.path.exists( dirName ):
        print 'Input Directory %s does not exist' % (dirName)
        if logFlg: logFlg.error('Directory %s does not exist' % dirName)
        if exit: sys.exit()
        return False
    else:
        return True

def ckFile(fName,logFlg=False,exit=False):
    '''Check if a file exists'''
    if not os.path.isfile(fName):
        print 'File %s does not exist' % (fName)
        if logFlg: logFlg.error('Unable to find file: %s' % fName)
        if exit: sys.exit()
        return False
    else:
        return True

                                    #----------------------------#
                                    #                            #
                                    #        --- Main---         #
                                    #                            #
                                    #----------------------------#

def main(argv):

    #----------------
    # Initializations for insitu
    #----------------
    dataDir     = '/data1/ancillary_data/fl0/FRAPPE/BAO/'

    #----------------
    # Initializations for FTIR
    #----------------

    loc        = 'fl0'                  # Name of station location
    gasName    = 'ccl4'                   # Name of gas
    ver        = 'Current_v5'           # Name of retrieval version to process
    ctlF       = 'sfit4_v5.ctl'            # Name of ctl file

    #gasName    = 'nh3'                   # Name of gas
    #ver        = 'Current_v2'           # Name of retrieval version to process
    #ctlF       = 'sfit4_v2.ctl'            # Name of ctl file

    errorFlg   = False                  # Flag to process error data
    fltrFlg    = True                   # Flag to filter the data
    rmsFlg     = True                   # Flag to filter based on max RMS
    dofFlg     = True                   # Flag to filter based on min DOFs

    maxRMS     = 1.5                     # Max Fit RMS to filter data. Data is filtered according to <= maxrms
    minDOF     = 0.5                       # Min DOFs for filtering

    sclfct     = 1.0E9                  # Scale factor to apply to vmr plots (ppmv=1.0E6, ppbv=1.0E9, etc)
    sclfctName = 'ppbv'                 # Name of scale factor for labeling plots

    #----------------------------
    # Partial Columns Bounds [km] 
    #----------------------------
    pCols = [[0.0,8.0]]
    
    #----------------------
    # Date range to process
    #----------------------
    iyear      = 2014
    imnth      = 7
    iday       = 8
    fyear      = 2014
    fmnth      = 8
    fday       = 23

    saveFlg     = False
    pltFile    = '/data/iortega/results/fl0/NH3_FRAPPE.pdf'

    #-------------------------------------------------
    #Reading both standard files
    #-------------------------------------------------

    if saveFlg: pdfsav = PdfPages(pltFile)
    else:       pdfsav = False

    i_date   = dt.date(iyear,imnth,iday)                                                     # Initial Day
    f_date   = dt.date(fyear,fmnth,fday)                                                     # Final Day
    ndays    = (f_date + dt.timedelta(days=1) - i_date).days
    dateList =[i_date + dt.timedelta(days=i) for i in range(0, ndays, 1)]  

    #--------------------------------------------
    # Walk through first level of directories and
    # collect directory names for processing
    #--------------------------------------------

    dirlist = []

    for dd in dateList:
    	# Find year month and day strings
        yrstr   = "{0:02d}".format(dd.year)
        mnthstr = "{0:02d}".format(dd.month)
        daystr  = "{0:02d}".format(dd.day)

        filename = 'frappe-PANGC_GROUND-BAO-TOWER_'+yrstr+mnthstr+daystr+'*.ict'
        #filename = 'FRAPPE-QCTILDAS-NH3_GROUND-BAO-TOWER_'+yrstr+mnthstr+daystr+'*.ict'
        
        icarttfile = glob.glob( dataDir + filename )

        if not icarttfile: continue

        for k in icarttfile:
            dirlist.append(k)

    dirlist.sort()

    insituvmr = []
    insitudate = []

    print 'Reading in-situ Files:'

    for indMain,sngDir in enumerate(dirlist):
    	ckDir(sngDir)
    	keys, data,dateutc, vnames = mf.read_ICARTT(sngDir)
    	insituvmr.append(data[vnames[1]])
    	insitudate.append(dateutc)


    insituvmr = np.asarray(insituvmr, dtype=np.float32)
    insitudate = np.asarray(insitudate)
    insituvmr=np.array(insituvmr).flatten()
    insitudate=np.array(insitudate).flatten()

    index = mf.dbFilterUL(insituvmr, lwrVal=-0.1, uprVal=10000.)

    insituvmr = insituvmr[index]
    insitudate = insitudate[index]


    #---------------------------------------------------------------------------------------
    #
    #                                    READING FTS data
    #
    #---------------------------------------------------------------------------------------

    print '\nReading FTS:'

    retDir = '/data1/ebaumer/'+loc.lower()+'/'+gasName.lower()+'/'+ver+'/'
    ctlFile  = '/data1/ebaumer/'+loc.lower()+'/'+gasName.lower()+'/'+'x.'+gasName.lower()+'/'+ctlF

    #---------------------------------
    # Check for the existance of files 
    # directories from input file
    #---------------------------------
    ckDir(retDir, exit=True)
    ckFile(ctlFile, exit=True)

    if saveFlg:  ckDir(os.path.dirname(os.path.realpath(pltFile)),exit=True)

    #-------------------------
    # Create Instance of Class
    #-------------------------
    gas = dc.PlotData(retDir,ctlFile,iyear=iyear,imnth=imnth,iday=iday,
        fyear=fyear,fmnth=fmnth,fday=fday,outFname=pltFile)

    #----------------------
    # Call to plot profiles
    #----------------------
    ds = gas.pltPrf(fltr=fltrFlg,allGas=True,sclfct=sclfct,sclname=sclfctName,
            errFlg=errorFlg, maxRMS=maxRMS,minTC=0.0, maxTC=1e24, minDOF=minDOF,dofFlg=dofFlg,rmsFlg=rmsFlg,)

    #----------------------
    # Call to plot profiles
    #----------------------
    for gas in ds['localGasList']: print 'Gas included in the Fit: '+str(gas)
    prf     = ds['rPrf']['PAN']
    alt     = ds['alt']
    airmass = ds['airmass']
    vmr     = prf[:,-1]
    vmr2     = prf[:,-1:]
    dates    = ds['dates']

    #---------------------------------------------------------------------------------------
    #
    #                                    PLOTS
    #
    #---------------------------------------------------------------------------------------

    fig, ax = plt.subplots(figsize=(10,6), sharex=True)
    ax.plot(insitudate, insituvmr/1000., 'k', linewidth=2,  label='insitu', alpha=0.5)
    ax.scatter(dates, vmr, facecolors='white', s=70, color='r',label='FTS')
    ax.set_ylabel('VMR [ppbv]', fontsize=16)
    ax.set_xlabel('Time [UT]', fontsize=16)
    ax.set_xlim(i_date, f_date)
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.grid(True)
    ax.legend(prop={'size':12})
    fig.subplots_adjust(left = 0.12, bottom=0.12, top=0.96, right = 0.95)
    fig.autofmt_xdate()

    doy_fts = mf.toYearFraction(dates)
    doy_insitu = mf.toYearFraction(insitudate)

    insituvmr_interp  = np.interp(doy_fts,doy_insitu , insituvmr)
   
    fig2, ax = plt.subplots(figsize=(10,6), sharex=True)
    ax.scatter(dates, insituvmr_interp/1000.,  facecolors='blue', color='blue', s=35, label='insitu')
    ax.scatter(dates, vmr, facecolors='white', s=35, color='r',label='FTS')
    ax.set_ylabel('VMR [ppbv]', fontsize=16)
    ax.set_xlabel('Time [UT]', fontsize=16)
    ax.set_xlim(i_date, f_date)
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.grid(True)
    ax.legend(prop={'size':12})
    fig2.subplots_adjust(left = 0.12, bottom=0.12, top=0.96, right = 0.95)
    fig2.autofmt_xdate()

    if saveFlg: 
       pdfsav.savefig(fig,dpi=200)
       #pdfsav.savefig(fig2,dpi=200)
       
    else:           
        plt.show(block=False)

    if saveFlg: pdfsav.close()

    print('\nFinished Plots.......\n')
       
      #--------------------------------
      # Pause so user can look at plots
      #--------------------------------
    if not saveFlg:
        user_input = raw_input('Press any key to exit >>> ')
        sys.exit()           # Exit program        


if __name__ == "__main__":
    main(sys.argv[1:])