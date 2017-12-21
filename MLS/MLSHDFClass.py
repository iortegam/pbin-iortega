
#----------------------------------------------------------------------------------------
# Name:
#        MLSHDFClass.py
#
# Purpose:
#       This is a collection of classes and functions used for processing and ploting HDF files
#
#
# Notes:
#
#
# Version History:
#       Created, Sep, 2017  Ivan Ortega (iortega@ucar.edu)
#
#
#----------------------------------------------------------------------------------------

import datetime as dt
import time
import math
import sys
import numpy as np
import os
import csv
import itertools
from collections import OrderedDict
import os
from os import listdir
from os.path import isfile, join
import re
#import statsmodels.api as sm
from scipy.integrate import simps
import matplotlib.animation as animation
import matplotlib
# Force matplotlib to not use any Xwindows backend.
#matplotlib.use('Agg')
import hdfBaseRetDat
import matplotlib.dates as md
from matplotlib.dates import DateFormatter, MonthLocator, YearLocator

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import FormatStrFormatter, MultipleLocator,AutoMinorLocator,ScalarFormatter
from matplotlib.backends.backend_pdf import PdfPages #to save multiple pages in 1 pdf...
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as mplcm
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
import matplotlib.gridspec as gridspec
import glob

from pyhdf.SD import SD, SDC
from pyhdf.SD import *
import coda

import h5py
import myfunctions as mf

from netCDF4 import Dataset, Group
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA




                                            #------------------#
 

def ckDir(dirName,logFlg=False,exitFlg=False):
    ''' '''
    if not os.path.exists( dirName ):
        print 'Input Directory %s does not exist' % (dirName)
        if logFlg: logFlg.error('Directory %s does not exist' % dirName)
        if exitFlg: sys.exit()
        return False
    else:
        return True   

def ckFile(fName,logFlg=False,exitFlg=False):
    '''Check if a file exists'''
    if not os.path.isfile(fName):
        print 'File %s does not exist' % (fName)
        if logFlg: logFlg.error('Unable to find file: %s' % fName)
        if exitFlg: sys.exit()
        return False
    else:
        return True   
    
def tryopen(fname):
    ''' Try to open a file and read contents'''
    try:
        with open(fname, 'r' ) as fopen:
            return fopen.readlines()
    except IOError as errmsg:
        print errmsg
        return False


def toYearFraction(dates):
    ''' Convert datetime to year and fraction of year'''

    #def sinceEpoch(date): # returns seconds since epoch
        #return time.mktime(date.timetuple())
    #s = sinceEpoch
    ep_fnc = lambda x: time.mktime(x.timetuple())
    
    retrnDates = np.zeros(len(dates))
    
    for i,sngDate in enumerate(dates):
        year = sngDate.year
        startOfThisYear = dt.datetime(year=year, month=1, day=1)
        startOfNextYear = dt.datetime(year=year+1, month=1, day=1)
    
        yearElapsed = ep_fnc(sngDate) - ep_fnc(startOfThisYear)
        yearDuration = ep_fnc(startOfNextYear) - ep_fnc(startOfThisYear)
        fraction = yearElapsed/yearDuration
        retrnDates[i] = sngDate.year + fraction

    return retrnDates


def dailyAvg(data, dates, dateAxis=1, meanAxis=0, quad=0):
    ''' Creates daily averages of specified quantity'''

    #-----------------------------
    # Initialize output dictionary
    #-----------------------------
    outdata = {}
    dataAvg = []
    cnt     = []
    std     = []

    #-----------------------------------------
    # Convert date time input to strictly date
    #-----------------------------------------
    dates_daily = np.asarray([dt.date(d.year,d.month,d.day) for d in dates])

    #-----------------------------
    # Ensure data is a numpy array
    #-----------------------------
    data = np.asarray(data)

    #--------------------------------------
    # Create a list of days to loop through
    #--------------------------------------
    uniqueDays = list(set(dates_daily))          # Find a list of unique days
    uniqueDays.sort()
    uniqueDays = np.asarray(uniqueDays)    

    #-------------------------
    # Loop through unique days
    #-------------------------
    for indDay in uniqueDays:

        inds = np.where(dates_daily == indDay)[0]                                       # Index for values to use as daily average
       
        if data.ndim == 1:                                                              # For single dimension array   
            if quad == 1: dataAvg.append( np.sqrt(np.sum(data[inds]**2)) / len(inds) )  # Find daily average (sqrt(sum(x^2))/N)
            else:         dataAvg.append( np.mean(data[inds]) )                         # Find daily average (straight mean)
            std.append(np.std(data[inds]))                                              # Find standard deviation
            
        else:                                                                           # For multi-Dimension array
            s           = [slice(None)] * data.ndim
            s[dateAxis] = inds
            tempMat     = data[s]
            if quad == 1: dataAvg.append( np.sqrt(np.sum(tempMat,axis=meanAxis) / len(inds)))
            else:         dataAvg.append( np.mean(tempMat,axis=meanAxis) )
            std.append(np.std(tempMat,axis=meanAxis))                                  # Find standard deviation
        
        cnt.append(len(inds))                                                         # Number of observations used in daily average

    cnt     = np.asarray(cnt)
    dataAvg = np.asarray(dataAvg)
    std     = np.asarray(std)

    outdata['dailyAvg'] = dataAvg
    outdata['dates']    = uniqueDays
    outdata['cnt']      = cnt
    outdata['std']      = std

    return outdata


def mnthlyAvg(data, dates, dateAxis=1, meanAxis=0, quad=0):
    ''' Creates monthly averages of specified quantity'''

    #-----------------------------
    # Initialize output dictionary
    #-----------------------------
    outdata = {}
    dataAvg = []
    cnt     = []
    std     = []

    #-----------------------------------------
    # Convert date time input to strictly date
    #-----------------------------------------
    dates_mnthly = np.asarray([dt.date(d.year,d.month,1) for d in dates])

    #-----------------------------
    # Ensure data is a numpy array
    #-----------------------------
    data = np.asarray(data)

    #----------------------------------------
    # Create a list of months to loop through
    #----------------------------------------
    uniqueMnths = list(set(dates_mnthly))          # Find a list of unique months
    uniqueMnths.sort()
    uniqueMnths = np.asarray(uniqueMnths)    

    #---------------------------
    # Loop through unique months
    #---------------------------
    for indDay in uniqueMnths:

        inds = np.where(dates_mnthly == indDay)[0]                                      # Index for values to use as monhtly average
        
        if data.ndim == 1:                                                              # For single dimension array   
            if quad == 1: dataAvg.append( np.sqrt(np.sum(data[inds]**2)) / len(inds) )  # Find monthly average (sqrt(sum(x^2))/N)
            else:         dataAvg.append( np.mean(data[inds]) )                         # Find monthly average (straight mean)
            std.append(np.std(data[inds]))                                              # Find standard deviation
            
        else:                                                                           # For multi-Dimension array
            s           = [slice(None)] * data.ndim
            s[dateAxis] = inds
            tempMat     = data[s]
            if quad == 1: dataAvg.append( np.sqrt(np.sum(tempMat,axis=meanAxis) / len(inds)))
            else:         dataAvg.append( np.mean(tempMat,axis=meanAxis) )
            std.append(np.std(tempMat,axis=meanAxis))                                  # Find standard deviation    

        cnt.append(len(inds))                                                          # Number of observations used in monhtly average

    cnt     = np.asarray(cnt)
    dataAvg = np.asarray(dataAvg)
    std     = np.asarray(std)

    outdata['mnthlyAvg'] = dataAvg
    outdata['dates']     = uniqueMnths
    outdata['cnt']       = cnt
    outdata['std']       = std

    return outdata


def allMnthAvg(data, dates):
    ''' Creates monthly averages of specified quantity'''

    #-----------------------------
    # Initialize output dictionary
    #-----------------------------
    outdata = {}
    dataAvg = []
    cnt     = []
    std     = []
    mnths   = np.array(range(1,13))

    #-----------------------------
    # Ensure data is a numpy array
    #-----------------------------
    data = np.asarray(data)

    #------------------------------
    # Loop through months (Jan-Dec)
    #------------------------------
    monthlist = np.array([d.month for d in dates])   # Create a month list

    for indMnth in mnths:
        inds = np.where(monthlist == indMnth)[0]     # Index for values to use as monhtly average
        dataAvg.append(np.mean(data[inds]))          # Find monhtly average
        cnt.append(len(inds))                        # Number of observations used in monhtly average
        std.append(np.std(data[inds]))               # Find standard deviation

    cnt     = np.asarray(cnt)
    dataAvg = np.asarray(dataAvg)
    std     = np.asarray(std)

    outdata['mnthlyAvg'] = dataAvg
    outdata['month']     = mnths
    outdata['cnt']       = cnt
    outdata['std']       = std

    return outdata



    

                                                #----------------#
                                                # Define classes #
                                                #----------------#

#------------------------------------------------------------------------------------------------------------------------------    
class _DateRange(object):
    '''
    This is an extension of the datetime module.
    Adds functionality to create a list of days.
    '''
    def __init__(self,iyear,imnth,iday,fyear,fmnth,fday, incr=1):
        self.i_date   = dt.date(iyear,imnth,iday)                                                     # Initial Day
        self.f_date   = dt.date(fyear,fmnth,fday)                                                     # Final Day
        self.dateList =[self.i_date + dt.timedelta(days=i) for i in range(0, self.numDays(), incr)]   # Incremental day list from initial to final day

    def numDays(self):
        '''Counts the number of days between start date and end date'''
        return (self.f_date + dt.timedelta(days=1) - self.i_date).days

    def inRange(self,crntyear,crntmonth,crntday):
        '''Determines if a specified date is within the date ranged initialized'''
        crntdate = dt.date(crntyear,crntmonth,crntday)
        if self.i_date <= crntdate <= self.f_date:
            return True
        else:
            return False

    def nearestDate(self, year, month, day=1, daysList=False):
        ''' Finds the nearest date from a list of days based on a given year, month, and day'''
        testDate = dt.date(year, month, day)
        if not daysList:
            daysList = self.dateList
        return min( daysList, key=lambda x:abs(x-testDate) )

    def yearList(self):
        ''' Gives a list of unique years within DateRange '''
        years = [ singDate.year for singDate in self.dateList]               # Find years for all date entries
        years = list(set(years))                                             # Determine all unique years
        years.sort()
        return years

    def daysInYear(self,year):
        ''' Returns an ordered list of days from DateRange within a specified year '''
        if isinstance(year,int):
            newyears = [inYear for inYear in self.dateList if inYear.year == year]
            return newyears
        else:
            print 'Error!! Year must be type int for daysInYear'
            return False


#------------------------------------------------------------------------------------------------------------------------------ 
#------------------------------------------------------------------------------------------------------------------------------                     
class ReadHDFMLS():

    def __init__(self, dataDir, iyear=False, fyear=False):

        #---------------------------------
        #
        #---------------------------------
        if not( dataDir.endswith('/') ):
                dataDir = dataDir + '/'

        self.dataDir    = dataDir

        self.HDF     = {}
        self.ReadHDFMLSFlg = True

        #-------------------------
        Vars = {}
        Vars.setdefault(self.getGas(),[])
        Vars.setdefault(self.getPressure(),[])
        Vars.setdefault(self.getLatitude(),[])
        Vars.setdefault(self.getLongitude(),[])
        Vars.setdefault(self.getStatus(),[])
        Vars.setdefault(self.getQuality(),[])

        #--------------------------
        if (not iyear) and (not fyear):
            print "Must specify initial and final years"
            exit()

       
        Nyears = fyear - iyear
        years  = [(iyear + i) for i in range(Nyears+1)]

        DirFiles  = []


        for y in years:
            DirFiles.append(glob.glob(dataDir + '*'+ str(y) + '*.he5'))

        DirFiles = np.asarray(DirFiles)
        DirFiles = np.concatenate( DirFiles, axis=0 )
        DirFiles =  np.sort(DirFiles)

  
        for drs in DirFiles:

            print '\nReading HDF File: %s' % (drs)
            ckFile(drs)

            drsSplt  = drs.strip().split('/')[-1]
            yyyy     = int(drsSplt[29:33])
            ddd      = int(drsSplt[34:37])
        
            #--------------------------
            #
            #--------------------------
            try:
                file    = h5py.File(drs, 'r', driver='core', backing_store=False)

                Date =  dt.date(yyyy, 1, 1) + dt.timedelta(ddd - 1)
                self.HDF.setdefault('Date',[]).append(Date)

                for var in Vars.keys():
                
                    data = file[var]

                    isDataset = isinstance(data, h5py.Dataset)

                    #print '(Dataset)', data.name, '    len =', data.shape #, g.dtype
                    self.HDF.setdefault(var,[]).append(np.asarray(data))

                #file.flush()
                file.close()
               
            except Exception as errmsg:
                print errmsg
                continue

        #---------------------------------
        #FLATTENED THE ARRAYS
        #---------------------------------
        # for var in Vars.keys():
        #   try:
        #       self.HDF[var] = list(itertools.chain.from_iterable(self.HDF[var]))
        #       self.HDF[var] = np.asarray(self.HDF[var])
        #   except Exception: continue

        print '\nVariables in HDF Files: {}'.format(Vars.keys())
       
        
    def getGas(self):
        return '/HDFEOS/SWATHS/HCl/Data Fields/HCl'

    def getPressure(self):
        return '/HDFEOS/SWATHS/HCl/Geolocation Fields/Pressure'

    def getLatitude(self):
        return '/HDFEOS/SWATHS/HCl/Geolocation Fields/Latitude'

    def getLongitude(self):
        return '/HDFEOS/SWATHS/HCl/Geolocation Fields/Longitude'

    def getStatus(self):
        return '/HDFEOS/SWATHS/HCl/Data Fields/Status'

    def getQuality(self):
        return '/HDFEOS/SWATHS/HCl/Data Fields/Quality'

class ReadNCMrg():

        def __init__(self, dataDir, iyear=False, fyear=False):

            #---------------------------------
            #
            #---------------------------------
            if not( dataDir.endswith('/') ):
                    dataDir = dataDir + '/'

            self.dataDir    = dataDir

            self.Mrg        = {}

            self.ReadNCMrgFlg  = True

            #--------------------------
            if (not iyear) and (not fyear):
                print "Must specify initial and final years"
                exit()

            
            Nyears = fyear - iyear
            years  = [(iyear + i) for i in range(Nyears+1)]

            DirFiles  = []

            for y in years:
                DirFiles.append(glob.glob(dataDir + '*'+ str(y) + '*.nc4'))

            DirFiles = np.asarray(DirFiles)
            DirFiles = np.concatenate( DirFiles, axis=0 )
            DirFiles =  np.sort(DirFiles)

            for drs in DirFiles:

                print '\nReading Nnc File: %s' % (drs)
                ckFile(drs)

                drsSplt  = drs.strip().split('/')[-1]
                yyyy     = int(drsSplt[26:30])

                #--------------------------
                #
                #--------------------------
                try:
                    nc_fid    = Dataset(drs , 'r')  # Dataset is the class behavior to open the file
                    
                    c_attrs, nc_dims, nc_vars = mf.ncdump(nc_fid, verb=False)
                    vgrp =  str(nc_fid.groups)

                    print vgrp
                    print c_attrs

                    #exit()

                    self.Mrg.setdefault('average',[]).append(np.asarray(nc_fid['/Merged/average']))
                    self.Mrg.setdefault('lat',[]).append(np.asarray(nc_fid['/Merged/lat']))
                    self.Mrg.setdefault('lev',[]).append(np.asarray(nc_fid['/Merged/lev']))
                    self.Mrg.setdefault('time',[]).append(np.asarray(nc_fid['/Merged/time']))

                    self.Mrg.setdefault('data_source',[]).append(np.asarray(nc_fid['/Merged/data_source']))
                    self.Mrg.setdefault('std_dev',[]).append(np.asarray(nc_fid['/Merged/std_dev']))
                    self.Mrg.setdefault('std_error',[]).append(np.asarray(nc_fid['/Merged/std_error']))
                    self.Mrg.setdefault('day_in_month',[]).append(np.asarray(nc_fid['/Merged/day_in_month']))
 
                except Exception as errmsg:
                    print errmsg
                    continue

class ReadAga():

    def __init__(self, dataDir):

        #---------------------------------
        #
        #---------------------------------
        if not( dataDir.endswith('/') ):
                dataDir = dataDir + '/'

        self.dataDir    = dataDir

        self.AGA        = {}

        self.ReadAgaFlg    = True
        #-------------------------

        fileName = glob.glob(dataDir + '*.mon')

        for file in fileName:

            try:
    
                f = open(file, 'r')

                header1 = f.readline()
                header2 = f.readline()

                hsplit = header2.split(',')
                site   = hsplit[1]
   
                lat    = float(hsplit[4])
                lon    = float(hsplit[7])

                self.AGA.setdefault('site',[]).append(site)
                self.AGA.setdefault('lat',[]).append(lat)
                self.AGA.setdefault('lon',[]).append(lon)

                header3 = [f.readline() for i in range(14)]

                time          = []
                MM            = []
                YYYY          = []
                CFC11         = []
                CFC11_e       = []
                
                CFC11         = []
                CFC11_e       = []
                
                CFC12         = []
                CFC12_e       = []
                
                CFC113         = []
                CFC113_e       = []
                
                CCl4         = []
                CCl4_e       = []

                Dates          = []

                for line in f:
                    time        = line.strip()
                    columns     = line.split()

                    time        = float(columns[0])
                    MM          = int(columns[1])
                    YYYY        = int(columns[2])
                    Dates.append(dt.date(YYYY, MM, 15))       
                   
                    CFC11.append(float(columns[4]))
                    CFC11_e.append(float(columns[5]))

                    CFC12.append(float(columns[7]))
                    CFC12_e.append(float(columns[8]))

                    CCl4.append(float(columns[13]))
                    CCl4_e.append(float(columns[14]))

                    CFC113.append(float(columns[19]))
                    CFC113_e.append(float(columns[20]))
                

                self.AGA.setdefault(site+'CFC11',[]).append(CFC11)
                self.AGA.setdefault(site+'CFC11_e',[]).append(CFC11_e)

                self.AGA.setdefault(site+'CFC12',[]).append(CFC12)
                self.AGA.setdefault(site+'CFC12_e',[]).append(CFC12_e)

                self.AGA.setdefault(site+'CFC113',[]).append(CFC113)
                self.AGA.setdefault(site+'CFC113_e',[]).append(CFC113_e)

                self.AGA.setdefault(site+'CCl4',[]).append(CCl4)
                self.AGA.setdefault(site+'CCl4_e',[]).append(CCl4_e)
                
                self.AGA.setdefault(site+'Dates',[]).append(Dates)

            except Exception as errmsg:
                print '\nError: ', errmsg

            


class PlotMLS(ReadHDFMLS):

    def __init__(self,dataDir, iyear=False, fyear=False, saveFlg=False, outFname=''):
       
        #------------------------------------------------------------
        # If outFname is specified, plots will be saved to this file,
        # otherwise plots will be displayed to screen
        #------------------------------------------------------------
        if saveFlg:  self.pdfsav = PdfPages(outFname)
        else:        self.pdfsav = False

        self.dataDir = dataDir

        #---------------
        # ReadOutputData
        #---------------
        ReadHDFMLS.__init__(self,dataDir, iyear=iyear, fyear=fyear)   
       
        
    def closeFig(self):
        self.pdfsav.close()


    def PlotMLSSet(self):

        #----------------------------------------
        #DEFINE VARIABLES
        #----------------------------------------
        try:
            Dates       = np.asarray(self.HDF['Date'])
            Prfs        = np.asarray(self.HDF[self.getGas()])
            Latitude    = np.asarray(self.HDF[self.getLatitude()])
            Longitude   = np.asarray(self.HDF[self.getLongitude()])
            Pressure    = np.asarray(self.HDF[self.getPressure()])

            Status     = np.asarray(self.HDF[self.getStatus()])
            Quality    = np.asarray(self.HDF[self.getQuality()])

        except Exception as errmsg:
            print '\nError: ', errmsg

        Nobs = len(Dates)

        print 'N obs = ', len(Dates)

        #print Prfs.shape
        #print Latitude.shape
        #print Longitude.shape
        #print Pressure.shape


        Prfs2       = []
        Dates2      = []
        Latitude2   = []
        Longitude2  = []
        Pressure2   = []
        

        for i in range(Nobs):
   
            indsQA = np.where( (Status[i] % 2 == 0) & (Quality[i] > 1.2) )[0]

            inds   = np.where( (Latitude[i][indsQA] > 20.0) & (Latitude[i][indsQA] < 40.0) )[0]

            if len(inds) > 1:
                Prfs2.append(np.mean(Prfs[i][indsQA[inds], :], axis=0))
                Dates2.append(Dates[i])
                Pressure2.append(Pressure[i, :])
        
        Pressure2 = np.asarray(Pressure2)
        Pressure2 = np.asarray(Pressure2[0, :])
        Prfs2     = np.asarray(Prfs2)
        Dates2    = np.asarray(Dates2)


        #print Pressure2.shape
        #print Prfs2.shape
        #print Dates2.shape

        indsPres = np.where( (Pressure2 > 0.32) & (Pressure2 < 100.0)  )[0]

        Pressure2 = Pressure2[indsPres]
        Prfs2     = Prfs2[:, indsPres]

        #----------------------------------------------------------------------------
        #
        #----------------------------------------------------------------------------

        years = [ singDate.year for singDate in Dates2]      # Find years for all date entries
        if len(list(set(years))) > 1: yrsFlg = True         # Determine all unique years
        else:                         yrsFlg = False   

        yearsLc      = YearLocator()
        #months       = MonthLocator(bymonth=1,bymonthday=1)
        monthsAll    = MonthLocator()
        months       = MonthLocator()
        if yrsFlg: DateFmt   = DateFormatter('%Y')
        else: DateFmt      = DateFormatter('%m\n%Y')

        POI = [0.46, 1.0, 2.1, 6.8, 10.0, 21.0, 46.0, 68.0]

        #----------------------------------------------------------------------------
        #
        #----------------------------------------------------------------------------

        fig = plt.figure(figsize=(12,10))

        outer_grid = gridspec.GridSpec(4, 2, wspace=0.175, hspace=0.175)

        for i, p in enumerate(POI):

            ind1 = mf.nearestind(p, Pressure2)

            ax = plt.Subplot(fig, outer_grid[i])
            
            ax.plot(Dates2,Prfs2[:,ind1],'k.',markersize=4)

            ax.grid(True,which='both', alpha=0.5)    
            ax.tick_params(which='both',labelsize=12)
            #ax.set_ylim(ymin=0, ymax=3e-9)
            #ax.set_xlim(0,600) 
            ax.set_title(str(p) + ' hPa')

            if yrsFlg:
                ax.xaxis.set_major_locator(yearsLc)
                ax.xaxis.set_minor_locator(months)
                #ax1.xaxis.set_minor_formatter(DateFormatter('%m'))
                ax.xaxis.set_major_formatter(DateFmt)

            else:
                ax.xaxis.set_major_locator(monthsAll)
                ax.xaxis.set_major_formatter(DateFmt)

            
            fig.add_subplot(ax)

        all_axes = fig.get_axes()
        #show only the outside spines
        for ax in all_axes:
            for sp in ax.spines.values():
                sp.set_visible(False)
                plt.setp(ax.get_xticklabels(), visible=False)
            if ax.is_first_row():
                ax.spines['top'].set_visible(True)
            if ax.is_last_row():
                ax.spines['bottom'].set_visible(True)
                plt.setp(ax.get_xticklabels(), visible=True)
            if ax.is_first_col():
                ax.spines['left'].set_visible(True)
            if ax.is_last_col():
                ax.spines['right'].set_visible(True)

        if (2 % 2 == 1): #even

            all_axes[-2].spines['bottom'].set_visible(True)
            plt.setp(all_axes[-2].get_xticklabels(), visible=True)
            all_axes[-2].set_zorder(1)

        all_axes[-1].set_xlabel('Date', fontsize=14)
        all_axes[-2].set_xlabel('Date', fontsize=14)

        fig.text(0.03, 0.5, 'Mixing Ratios ', fontsize=14, va='center', rotation='vertical')
        #plt.suptitle('{} Vertical Profiles'.format(gasName.upper()), fontsize=16  )
        #plt.tight_layout(h_pad=0.25) #w_pad=1.75 pad=1.75,
        fig.autofmt_xdate()
        fig.subplots_adjust(left=0.1, bottom=0.075, right=0.95, top=0.95)
        
        if self.pdfsav: 
            self.pdfsav.savefig(fig,dpi=200)
        else:
            plt.show(block=False)
            #user_input = raw_input('Press any key to exit >>> ')
            #sys.exit()


        #----------------------------------------------------------------------------
        #CONTOUR PLOT: TIME SERIES OF VERTICAL PROFILES
        #---------------------------------------------------------------------------
        # print Pressure2.shape
        # print Dates2.shape
        # print Prfs2.shape
        # #levels1 = np.arange(0., 10,  0.5)        
        # fig,ax1 = plt.subplots(figsize=(10, 6))
        
        # #cax1          = ax1.contourf(Dates2, np.flipud(Pressure2), np.flipud(np.transpose(Prfs2)),  cmap=mplcm.jet)
        
        # cax1          = ax1.contourf(Dates2, Pressure2, np.transpose(Prfs2),  cmap=mplcm.jet)
        
        # divider1      = make_axes_locatable(ax1)
        # cb1           = divider1.append_axes("right",size="5%",pad=0.075)

        # cbar1         = plt.colorbar(cax1,cax=cb1)
        # cbar1.set_label('VMR [ppt]', fontsize=12)
        # cbar1.ax.tick_params(labelsize=12)
       
        # ax1.grid(True)
        # ax1.set_xlabel('Date',  fontsize=14)
        # ax1.set_ylabel('Pressure [hPa]',  fontsize=14)   
        # #ax1.set_ylim((0,35))     
        # ax1.tick_params(labelsize=14)
        # ax1.set_title('Time Series of vertical Profiles', fontsize=14)
        # ax1.set_ylabel('Log Scale VMR')
        # ax1.set_yscale('log')
        

        # if yrsFlg:
        #     ax1.xaxis.set_major_locator(yearsLc)
        #     ax1.xaxis.set_minor_locator(months)
        #     #ax1.xaxis.set_minor_formatter(DateFormatter('%m'))
        #     ax1.xaxis.set_major_formatter(DateFmt)

        # else:
        #     ax1.xaxis.set_major_locator(monthsAll)
        #     ax1.xaxis.set_major_formatter(DateFmt)

        # fig.autofmt_xdate()
        # fig.subplots_adjust(left = 0.1, bottom=0.15, top=0.9, right = 0.92)

        # if self.pdfsav: 
        #     self.pdfsav.savefig(fig,dpi=200)
        # else:
        #     plt.show(block=False)
        #     user_input = raw_input('Press any key to exit >>> ')
        #     sys.exit() 
           
class PlotMrg(ReadNCMrg):

    def __init__(self,dataDir, iyear=False, fyear=False, saveFlg=False, outFname=''):
       
        #------------------------------------------------------------
        # If outFname is specified, plots will be saved to this file,
        # otherwise plots will be displayed to screen
        #------------------------------------------------------------
        if saveFlg:  self.pdfsav = PdfPages(outFname)
        else:        self.pdfsav = False

        self.dataDir = dataDir

        #---------------
        # ReadOutputData
        #---------------
        ReadNCMrg.__init__(self,dataDir, iyear=iyear, fyear=fyear)   
       
        
    def closeFig(self):
        self.pdfsav.close()

    def PlotMrgSet(self):

        #----------------------------------------
        #DEFINE VARIABLES
        #----------------------------------------
        try:
            time        = np.asarray(self.Mrg['time'])
            Prfs        = np.asarray(self.Mrg['average'])
            lev        = np.asarray(self.Mrg['lev'])
            Latitude    = np.asarray(self.Mrg['lat'])
            day         = np.asarray(self.Mrg['day_in_month'])
           
            #Pressure    = np.asarray(self.Mrg['time'])

        except Exception as errmsg:
            print '\nError: ', errmsg

        Nfiles  = time.shape[0]
        NMonths = time.shape[1]
        Nlat    = Latitude.shape[1]
        Nlev    = lev.shape[1]

        print 'N files = ', Nfiles
        print 'N Months = ', NMonths
        print 'N Latitudes = ', Nlat
        print 'N Levels = ', Nlev

        
        Dates = np.asarray([[dt.date(1950, 1,1) + dt.timedelta(days=int(t)) for t in tt] for tt in time])

        Dates    = np.reshape(Dates, (Nfiles*NMonths))
        Prfs     = np.reshape(Prfs, (Nfiles*NMonths, Nlev,Nlat ))
        Pressure = lev[0, :]
        Lat      = Latitude[0, :]

        Prfs[Prfs < 0] = np.nan

        #----------------------------------------------------------------------------
        #
        #----------------------------------------------------------------------------

        years = [ singDate.year for singDate in Dates]      # Find years for all date entries
        if len(list(set(years))) > 1: yrsFlg = True         # Determine all unique years
        else:       yrsFlg = False   


        yearsLc      = YearLocator()
        #months       = MonthLocator(bymonth=1,bymonthday=1)
        monthsAll    = MonthLocator()
        months       = MonthLocator()
        if yrsFlg: DateFmt   = DateFormatter('%Y')
        else: DateFmt      = DateFormatter('%m\n%Y')

        POI = [0.46, 1.0, 2.1, 6.8, 10.0, 21.0, 46.0, 68.0]
        LOI = 45

        indLat = mf.nearestind(LOI, Lat)

        #----------------------------------------------------------------------------
        #
        #----------------------------------------------------------------------------

        fig = plt.figure(figsize=(12,10))

        outer_grid = gridspec.GridSpec(4, 2, wspace=0.175, hspace=0.175)

        for i, p in enumerate(POI):

            ind1 = mf.nearestind(p, Pressure)
            
            ax = plt.Subplot(fig, outer_grid[i])
            
            ax.plot(Dates,Prfs[:,ind1, indLat],'k.',markersize=4)

            ax.grid(True,which='both', alpha=0.5)    
            ax.tick_params(which='both',labelsize=12)
            #ax.set_ylim(ymin=0, ymax=3e-9)
            #ax.set_xlim(0,600) 
            ax.set_title(str(p) + ' hPa')

            if yrsFlg:
                ax.xaxis.set_major_locator(yearsLc)
                ax.xaxis.set_minor_locator(months)
                #ax1.xaxis.set_minor_formatter(DateFormatter('%m'))
                ax.xaxis.set_major_formatter(DateFmt)

            else:
                ax.xaxis.set_major_locator(monthsAll)
                ax.xaxis.set_major_formatter(DateFmt)

            
            fig.add_subplot(ax)

        all_axes = fig.get_axes()
        #show only the outside spines
        for ax in all_axes:
            for sp in ax.spines.values():
                sp.set_visible(False)
                plt.setp(ax.get_xticklabels(), visible=False)
            if ax.is_first_row():
                ax.spines['top'].set_visible(True)
            if ax.is_last_row():
                ax.spines['bottom'].set_visible(True)
                plt.setp(ax.get_xticklabels(), visible=True)
            if ax.is_first_col():
                ax.spines['left'].set_visible(True)
            if ax.is_last_col():
                ax.spines['right'].set_visible(True)

        if (2 % 2 == 1): #even

            all_axes[-2].spines['bottom'].set_visible(True)
            plt.setp(all_axes[-2].get_xticklabels(), visible=True)
            all_axes[-2].set_zorder(1)

        all_axes[-1].set_xlabel('Date', fontsize=14)
        all_axes[-2].set_xlabel('Date', fontsize=14)

        fig.text(0.03, 0.5, 'Mixing Ratios ', fontsize=14, va='center', rotation='vertical')
        #plt.suptitle('{} Vertical Profiles'.format(gasName.upper()), fontsize=16  )
        #plt.tight_layout(h_pad=0.25) #w_pad=1.75 pad=1.75,
        fig.autofmt_xdate()
        fig.subplots_adjust(left=0.1, bottom=0.075, right=0.95, top=0.95)
        
        if self.pdfsav: 
            self.pdfsav.savefig(fig,dpi=200)
        else:
            plt.show(block=False)
            #user_input = raw_input('Press any key to exit >>> ')
            #sys.exit()

class PlotAga(ReadAga):

    def __init__(self,dataDir, iyear=False, fyear=False, saveFlg=False, outFname=''):
       
        #------------------------------------------------------------
        # If outFname is specified, plots will be saved to this file,
        # otherwise plots will be displayed to screen
        #------------------------------------------------------------
        if saveFlg:  self.pdfsav = PdfPages(outFname)
        else:        self.pdfsav = False

        self.dataDir = dataDir

        #---------------
        # ReadOutputData
        #---------------
        ReadAga.__init__(self,dataDir)   
       
        
    def closeFig(self):
        self.pdfsav.close()

    def PloAgaSet(self):

        #----------------------------------------
        #DEFINE VARIABLES
        #----------------------------------------
        Nsites = len(self.AGA['site'])
        print 'N AGAGE sites:' + str(Nsites)

        s      = self.AGA['site']
        Lat   = np.asarray(self.AGA['lat'])


        indLat = mf.nearestind(45., Lat)

        print 'AGAGE sites:' + self.AGA['site'][indLat]
        print 'AGAGE sites Lat:' + str(self.AGA['lat'][indLat])


        Dates    = np.asarray(self.AGA[s[indLat]+'Dates'][0])

        CFC11    = np.asarray(self.AGA[s[indLat]+'CFC11'][0])
        CFC11_e  = np.asarray(self.AGA[s[indLat]+'CFC11_e'][0])

        CFC12    = np.asarray(self.AGA[s[indLat]+'CFC12'][0])
        CFC12_e  = np.asarray(self.AGA[s[indLat]+'CFC12_e'][0])

        CCl4     = np.asarray(self.AGA[s[indLat]+'CCl4'][0])
        CCl4_e   = np.asarray(self.AGA[s[indLat]+'CCl4_e'][0])

        CFC113   = np.asarray(self.AGA[s[indLat]+'CFC113'][0])
        CFC113_e = np.asarray(self.AGA[s[indLat]+'CFC113_e'][0])

        Cl = CFC11 + CFC12 + CCl4 


        #----------------------------------------
        years = [ singDate.year for singDate in Dates]      # Find years for all date entries
        if len(list(set(years))) > 1: yrsFlg = True         # Determine all unique years
        else:                         yrsFlg = False   

        yearsLc      = YearLocator()
        monthsAll    = MonthLocator()
        #months       = MonthLocator()
        months       = MonthLocator(bymonth=1, bymonthday=1)
        if yrsFlg: DateFmt   = DateFormatter('%Y')
        else: DateFmt        = DateFormatter('%m\n%Y')


        #----------------------------------------
        
        fig, ax  = plt.subplots()
        
        ax.plot(Dates,CFC11, color='k',label='CFC11')
        ax.fill_between(Dates,CFC11-CFC11_e,CFC11+CFC11_e,alpha=0.5,color='0.75') 

        ax.plot(Dates,CFC12, color='r',label='CFC12')
        ax.fill_between(Dates,CFC12-CFC12_e,CFC12+CFC12_e,alpha=0.5,color='r') 

        ax.plot(Dates,CCl4, color='b',label='CCl4')
        ax.fill_between(Dates,CCl4-CCl4_e,CCl4+CCl4_e,alpha=0.5,color='b') 

        ax.plot(Dates,CFC113, color='green',label='CFC113')
        ax.fill_between(Dates,CFC113-CFC113_e,CFC113+CFC113_e,alpha=0.5,color='green')  

        ax.plot(Dates,Cl, color='cyan' ,label='Chlorine')
        
        ax.legend(prop={'size':10})

        ax.set_ylabel('VMR [ppt]')
        ax.set_xlabel('Dates')      
        ax.grid(True,which='both', alpha=0.5)  

        ax.set_title(self.AGA['site'][indLat] + ' (Lat: '+str(self.AGA['lat'][indLat]) + ',  Lon:' +str(self.AGA['lon'][indLat])+ ')')   

        if yrsFlg:
            ax.xaxis.set_major_locator(yearsLc)
            ax.xaxis.set_minor_locator(months)
            #ax1.xaxis.set_minor_formatter(DateFormatter('%m'))
            ax.xaxis.set_major_formatter(DateFmt)

        else:
            ax.xaxis.set_major_locator(monthsAll)
            ax.xaxis.set_major_formatter(DateFmt)

        fig.autofmt_xdate()
        fig.subplots_adjust(left=0.1, bottom=0.12, right=0.975, top=0.95)
        
        if self.pdfsav: 
            self.pdfsav.savefig(fig,dpi=200)
        else:
            plt.show(block=False)
            #user_input = raw_input('Press any key to exit >>> ')
            #sys.exit()        
    

class PlotAll(ReadHDFMLS, ReadNCMrg, ReadAga):

    def __init__(self,MLSdataDir, MrgdataDir, AgaDir, iyear=False, fyear=False, saveFlg=False, outFname=''):
       
        #------------------------------------------------------------
        # If outFname is specified, plots will be saved to this file,
        # otherwise plots will be displayed to screen
        #------------------------------------------------------------
        if saveFlg:  self.pdfsav = PdfPages(outFname)
        else:        self.pdfsav = False

        self.ReadHDFMLSFlg = False
        self.ReadNCMrgFlg  = False
        self.ReadAgaFlg    = False

        #---------------
        # ReadOutputData
        #---------------
        #ReadHDFMLS.__init__(self,MLSdataDir, iyear=iyear, fyear=fyear)
        ReadNCMrg.__init__(self,MrgdataDir, iyear=iyear, fyear=fyear)
        ReadAga.__init__(self,AgaDir)   
       
        
    def closeFig(self):
        self.pdfsav.close()

    def PlotAllSet(self):


        #----------------------------------------
        #Define global interes
        #----------------------------------------
        POI   = [0.46, 1.0, 2.1, 6.8, 10.0, 21.0, 46.0, 68.0]   #Level pressure of interest
        LOI   = 45   #Ltitude of interest
        poi_i = 10.
        y4trend = 1997

        

                                        #----------------------------------------
                                        #           Read MLS
                                        #----------------------------------------

        if self.ReadHDFMLSFlg:

            #----------------------------------------
            #DEFINE VARIABLES
            #----------------------------------------
            try:
                Dates       = np.asarray(self.HDF['Date'])
                Prfs        = np.asarray(self.HDF[self.getGas()])*1e9
                Latitude    = np.asarray(self.HDF[self.getLatitude()])
                Longitude   = np.asarray(self.HDF[self.getLongitude()])
                Pressure    = np.asarray(self.HDF[self.getPressure()])

                Status     = np.asarray(self.HDF[self.getStatus()])
                Quality    = np.asarray(self.HDF[self.getQuality()])

            except Exception as errmsg:
                print '\nError: ', errmsg

            Nobs = len(Dates)

            print 'N obs = ', len(Dates)

            PrfsMLS       = []
            StDMLS        = []

            DatesMLS      = []
            LatMLS        = []
            LonMLS        = []
            PresMLS       = []
            

            for i in range(Nobs):
       
                indsQA = np.where( (Status[i] % 2 == 0) & (Quality[i] > 1.2) )[0]

                inds   = np.where( (Latitude[i][indsQA] > LOI - 5.) & (Latitude[i][indsQA] < LOI + 5) )[0]

                if len(inds) > 1:
                    PrfsMLS.append(np.nanmean(Prfs[i][indsQA[inds], :], axis=0))
                    StDMLS.append(np.nanstd(Prfs[i][indsQA[inds], :], axis=0))
                    DatesMLS.append(Dates[i])
                    PresMLS.append(Pressure[i, :])
            
            PresMLS = np.asarray(PresMLS)
            PresMLS = np.asarray(PresMLS[0, :])
            PrfsMLS     = np.asarray(PrfsMLS)
            StDMLS     = np.asarray(StDMLS)
            DatesMLS    = np.asarray(DatesMLS)

            indsPres = np.where( (PresMLS > 0.32) & (PresMLS < 100.0)  )[0]

            PresMLS    = PresMLS[indsPres]
            PrfsMLS    = PrfsMLS[:, indsPres]

                                        #----------------------------------------
                                        #           Read MERGE
                                        #----------------------------------------

        if self.ReadNCMrgFlg:

            #----------------------------------------
            #DEFINE VARIABLES
            #----------------------------------------
            try:
                time        = np.asarray(self.Mrg['time'])
                Prfs        = np.asarray(self.Mrg['average'])*1e9
                lev        = np.asarray(self.Mrg['lev'])
                Latitude    = np.asarray(self.Mrg['lat'])
                day         = np.asarray(self.Mrg['day_in_month'])
                std_dev     = np.asarray(self.Mrg['std_dev'])*1e9
               
                #Pressure    = np.asarray(self.Mrg['time'])

            except Exception as errmsg:
                print '\nError: ', errmsg

            Nfiles  = time.shape[0]
            NMonths = time.shape[1]
            Nlat    = Latitude.shape[1]
            Nlev    = lev.shape[1]

            print 'N files = ', Nfiles
            print 'N Months = ', NMonths
            print 'N Latitudes = ', Nlat
            print 'N Levels = ', Nlev

            
            Dates = np.asarray([[dt.date(1950, 1,1) + dt.timedelta(days=int(t)) for t in tt] for tt in time])

            DatesMrg    = np.reshape(Dates, (Nfiles*NMonths))
            PrfsMrg     = np.reshape(Prfs, (Nfiles*NMonths, Nlev,Nlat ))
            PresMrg     = lev[0, :]
            LatMrg      = Latitude[0, :]
            StDMrg     = np.reshape(std_dev, (Nfiles*NMonths, Nlev,Nlat ))

            #PrfsMrg[PrfsMrg < 0] = np.nan
            #StDMrg[PrfsMrg < 0]  = np.nan

            # rprf_neg = np.asarray(PrfsMrg[:,:,5]) <= 0
            # indsT = np.where( np.sum(rprf_neg,axis=1) > 0 )[0]
            # print ('Total number observations found with negative partial column = {}'.format(len(indsT)))

            indMrg1 = mf.nearestind(poi_i, PresMrg)

            indLatMrg = mf.nearestind(LOI, LatMrg)


            print PrfsMrg.shape
            indsBad =  np.where(PrfsMrg[:, indMrg1, indLatMrg] < 0.0)[0]

            print len(indsBad)
            PrfsMrg       = np.delete(PrfsMrg,indsBad,axis=0)
            StDMrg         = np.delete(StDMrg,indsBad,axis=0)
            DatesMrg       = np.delete(DatesMrg,indsBad,axis=0)

            print PrfsMrg.shape
            # print indsBad

            #exit()



                                        #----------------------------------------
                                        #           Read AGAGE
                                        #----------------------------------------

        if self.ReadAgaFlg:
            #----------------------------------------
            #DEFINE VARIABLES
            #----------------------------------------
            Nsites = len(self.AGA['site'])
            print 'N AGAGE sites:' + str(Nsites)

            s      = self.AGA['site']
            Lat   = np.asarray(self.AGA['lat'])

            indLatAGA = mf.nearestind(53., Lat)  ## This is to read MACE, HEAD
            ##indLatAGA = mf.nearestind(LOI, Lat)  ## This is to read Closest to LOI

            print self.AGA['site'][indLatAGA]
            print Lat[indLatAGA]  


            DatesAGA    = np.asarray(self.AGA[s[indLatAGA]+'Dates'][0])

            CFC11    = np.asarray(self.AGA[s[indLatAGA]+'CFC11'][0])
            CFC11_e  = np.asarray(self.AGA[s[indLatAGA]+'CFC11_e'][0])

            CFC12    = np.asarray(self.AGA[s[indLatAGA]+'CFC12'][0])
            CFC12_e  = np.asarray(self.AGA[s[indLatAGA]+'CFC12_e'][0])

            CCl4     = np.asarray(self.AGA[s[indLatAGA]+'CCl4'][0])
            CCl4_e   = np.asarray(self.AGA[s[indLatAGA]+'CCl4_e'][0])

            CFC113   = np.asarray(self.AGA[s[indLatAGA]+'CFC113'][0])
            CFC113_e = np.asarray(self.AGA[s[indLatAGA]+'CFC113_e'][0])

            Cl = CFC11 + CFC12 + CCl4 
            Cl_e = np.sqrt(CFC11_e**2 + CFC12_e**2 + CCl4_e**2)


                                        #----------------------------------------
                                        #           Read JFJ
                                        #----------------------------------------

        file = '/data1/Campaign/Satellite/MLS/HCl_ClONO2_Cly_Jungfraujoch_FTIR_1983-2016.txt'
        f = open(file, 'r')

        HCl_Date_JFJ       = []
        HCl_JFJ           = []
        CLONO2_Date_JFJ    = []
        CLONO2_JFJ         = []

        Cl_Date_JFJ        = []
        Cl_JFJ             = []

        header1 = f.readline()

        for line in f:
            columns     = line.split()
            atime       = float(columns[0])

            da_i = mf.t2dt(atime)
            if da_i.year >= 1983:
                HCl_Date_JFJ.append(da_i)
                HCl_JFJ.append(float(columns[1])/1e15)

            if len(columns) >=3:

                atime       = float(columns[2])
                da_i = mf.t2dt(atime) 
                if da_i.year >= 1983:
                    CLONO2_Date_JFJ.append(da_i)
                    CLONO2_JFJ.append(float(columns[3])/1e15)

            if len(columns) >= 5:
                atime       = float(columns[4])
                da_i = mf.t2dt(atime)
                if da_i.year >= 1983:
                    Cl_Date_JFJ.append(da_i)
                    Cl_JFJ.append(float(columns[5])/1e15)


        HCl_Date_JFJ = np.asarray(HCl_Date_JFJ)
        HCl_JFJ      = np.asarray(HCl_JFJ)

        Cl_JFJ       = np.asarray(Cl_JFJ)
        Cl_Date_JFJ  = np.asarray(Cl_Date_JFJ)



                                        #----------------------------------------
                                        #           Read MLO
                                        #----------------------------------------

        file = '/data1/ebaumer/mlo/hcl/MLO-FTS-HCL_GROUND-BOULDER.ascii'
        f = open(file, 'r')

        HCl_Date_MLO       = []
        HCl_MLO           = []
        

        headers = [f.readline() for i in range(7)]
        

        for line in f:
    
            columns     = line.split(',')

            yyyy = int(columns[1].strip()[0:4])
            mm   = int(columns[1].strip()[5:7])
            dd   = int(columns[1].strip()[8:10])

            hh   = int(columns[2].strip()[0:2])
            mi   = int(columns[2].strip()[3:5])
            se   = int(columns[2].strip()[6:8])

            da = dt.datetime(yyyy, mm, dd, hh, mi, se)
        
            HCl_Date_MLO.append(da)
            HCl_MLO.append(float(columns[3].strip()))


        HCl_MLO      = np.asarray(HCl_MLO)/1e15
        HCl_Date_MLO = np.asarray(HCl_Date_MLO)

        Avg = mf.mnthlyAvg(HCl_MLO, HCl_Date_MLO, dateAxis=1)
        HCl_Date_MLO        = Avg['dates']
        HCl_MLO             = Avg['mnthlyAvg']
        HCl_std_MLO         = Avg['std']

        
        #----------------------------------------------------------------------------
        #
        #----------------------------------------------------------------------------

        years = [ singDate.year for singDate in DatesMrg]      # Find years for all date entries
        if len(list(set(years))) > 1: yrsFlg = True         # Determine all unique years
        else:                         yrsFlg = False   

        yearsLc      = YearLocator(2)
        monthsAll    = MonthLocator()
        #months       = MonthLocator()
        months       = MonthLocator(bymonth=1, bymonthday=1)
        if yrsFlg: DateFmt   = DateFormatter('%Y')
        else: DateFmt        = DateFormatter('%m\n%Y')

        #----------------------------------------------------------------------------
        #
        #----------------------------------------------------------------------------
        indLatMrg = mf.nearestind(LOI, LatMrg)

        #-----------------------------------------
        #ONE PANEL (THREE Y-AXIS)
        #-----------------------------------------
        fig, ax = plt.subplots(figsize=(12, 7), sharex=True)
       
        axes = [ax, ax.twinx(), ax.twinx()]

        axes[-1].spines['right'].set_position(('axes', 1.125))
        axes[-1].set_frame_on(True)
        axes[-1].patch.set_visible(False)

        #-----------------------------------------
        #PLT GOZCARDS
        #-----------------------------------------
        if self.ReadHDFMLSFlg: ind1MLS = mf.nearestind(poi_i, PresMLS)
        if self.ReadNCMrgFlg:  indMrg1 = mf.nearestind(poi_i, PresMrg)

        indsGood =  np.where(PrfsMrg[:, indMrg1, indLatMrg] > 0.0)[0]

        PrfsMrg2        = PrfsMrg[indsGood, indMrg1, indLatMrg] #   np.delete(PrfsMrg,indsGood,axis=0)
        StDMrg2         = StDMrg[indsGood, indMrg1, indLatMrg]  # np.delete(StDMrg,indsGood,axis=0)
        DatesMrg2       = DatesMrg[indsGood]  #np.delete(DatesM
        
        #ax.plot(DatesMrg,PrfsMrg[:,indMrg1, indLatMrg],'k.',markersize=4)
        #axes[0].plot(DatesMrg2,PrfsMrg2, '.', color='gray', markersize=0)
        
        #axes[0].fill_between(DatesMrg,PrfsMrg2-StDMrg2,PrfsMrg2+StDMrg2, color='gray',alpha=0.5)

        PrfsMrg2_smth     = mf.smooth(PrfsMrg2,  window_len=36, window='flat')
        axes[0].plot(DatesMrg2,PrfsMrg2_smth,  linewidth=5, color='gray')
        
        yearGO = [single.year for single in DatesMrg2]
        yearGO = np.asarray(yearGO)

        inds = np.where( yearGO >= y4trend )[0]

        dateYearFrac = mf.toYearFraction(DatesMrg2[inds])
        weights      = np.ones_like(dateYearFrac)
        res          = mf.fit_driftfourier(dateYearFrac, PrfsMrg2[inds], weights, 2)
        f_drift, f_fourier, f_driftfourier = res[3:6]

        res_b        = mf.cf_driftfourier(dateYearFrac, PrfsMrg2[inds], weights, 2)
        perc, intercept_b, slope_b, pfourier_b = res_b



        #axes[2].plot(DatesAGA[inds],f_drift(dateYearFrac),label='Fitted Anual Trend')
        #axes[2].text(0.02,0.95,"Trends [%/yr]",transform=axes[2].transAxes)
        #axes[2].text(0.02,0.95,"GOZCARDS: {0:.2f} %/yr".format(res[1]/np.nanmean(PrfsMrg2[inds])*100.0),transform=axes[2].transAxes, color='gray', fontsize=12)
        print "Fitted trend GOZCARDS -- slope: {0:.3E} ({1:.3f}%)".format(res[1],res[1]/np.nanmean(PrfsMrg2[inds])*100.0)
        axes[0].scatter(DatesMrg2,PrfsMrg2, s=20, color='gray', edgecolors='gray', alpha=0.75, label='GOZCARDS (10hPa, 40-50$^\circ$ N)\n {0:.2f} $\pm$ {1:.2f} %/yr'.format(res[1]/np.nanmean(PrfsMrg2[inds])*100.0, np.std(slope_b)/np.mean(PrfsMrg2[inds])*100.0))   #facecolors='none'

        #axes[0].grid(True,which='both', alpha=0.5)    
        axes[0].tick_params(which='both',labelsize=14)
        #axes[0].legend(prop={'size':12}, loc=2)
        axes[0].set_ylim(1.7, 3.4)

        if not self.ReadHDFMLSFlg: axes[0].set_ylabel('VMR [ppb], GOZCARDS (HALOE, ACE & MLS) (HCl)', fontsize=14)


        if yrsFlg:
            axes[0].xaxis.set_major_locator(yearsLc)
            axes[0].xaxis.set_minor_locator(months)
            axes[0].xaxis.set_major_formatter(DateFmt)
            #axes[0].xaxis.set_major_formatter(DateFmt)

        else:
            axes[0].xaxis.set_major_locator(monthsAll)
            axes[0].xaxis.set_major_formatter(DateFmt)

        #Data = (DatesMrg2.T,PrfsMrg2.T)
        #np.savetxt('/data1/Campaign/Satellite/MLS/Cl_data.ascii', Data, fmt=['%s', '%1.4e'])
        

        #-----------------------------------------
        #PLT MLS
        #-----------------------------------------
        if self.ReadHDFMLSFlg:
            Avg          = mf.mnthlyAvg(PrfsMLS[:,ind1MLS], DatesMLS, dateAxis=1, meanAxis=0)
            Dates        = Avg['dates']
            AvgData      = Avg['mnthlyAvg']
            #stdData      = Avg['mnthlyAvg']
            axes[0].set_ylabel('VMR [ppb], GOZCARDS & MLS V4.2 (HCl)', fontsize=14)
            #axes[0].set_xlabel('Year', fontsize=14)

            axes[0].plot(Dates,AvgData, 'k.', markersize=0)
            axes[0].scatter(Dates,AvgData, s=20, color='r' , edgecolors='k', label='MLS V4.2')  #facecolors='none'
            Avg          = mf.mnthlyAvg(StDMLS[:,ind1MLS], DatesMLS, dateAxis=1, meanAxis=0)
            stdData      = Avg['mnthlyAvg']
            axes[0].fill_between(Dates,AvgData-stdData,AvgData+stdData, color='r', alpha=0.25)

        #-----------------------------------------
        #PLT AGAGE
        #-----------------------------------------
        
        #axes[1].fill_between(DatesAGA,Cl-Cl_e,Cl+Cl_e,alpha=0.5,color='green') 
        axes[1].set_ylabel('VMR [ppt], AGAGE (CFC11 + CFC12 + CCl$_4$)', color='green', fontsize=14)
        axes[1].tick_params(axis='y', colors='green', labelsize=14)
        #axes[1].legend(prop={'size':12})
        axes[1].set_ylim(380, 1100)

        yearAGA = [single.year for single in DatesAGA]
        yearAGA = np.asarray(yearAGA)

        inds = np.where( yearAGA >= y4trend )[0]

        dateYearFrac = mf.toYearFraction(DatesAGA[inds])
        weights      = np.ones_like(dateYearFrac)
        res          = mf.fit_driftfourier(dateYearFrac, Cl[inds], weights, 2)
        f_drift, f_fourier, f_driftfourier = res[3:6]

        res_b        = mf.cf_driftfourier(dateYearFrac, Cl[inds], weights, 2)
        perc, intercept_b, slope_b, pfourier_b = res_b
        #axes[2].plot(DatesAGA[inds],f_drift(dateYearFrac),label='Fitted Anual Trend')
        #axes[2].text(0.02,0.94,"Fitted trend -- slope: {0:.3E} ({1:.3f}%)".format(res[1],res[1]/np.mean(totClmn)*100.0),transform=ax1.transAxes)
        print "Fitted trend AGAGE -- slope: {0:.3E} ({1:.3f}%)".format(res[1],res[1]/np.mean(Cl[inds])*100.0)
        #axes[1].text(0.02,0.91,"AGAGE: {0:.2f} %/yr".format(res[1]/np.mean(Cl[inds])*100.0), transform=axes[1].transAxes, color='green', fontsize=12)
        #axes[1].plot(DatesAGA,Cl, linewidth=4, color='green',label='AGAGE Trinidad Head, CA (41.04 N): {0:.2f} %/yr'.format(res[1]/np.nanmean(Cl[inds])*100.0))
        axes[1].scatter(DatesAGA,Cl, s=20, color='green', edgecolors='green', alpha=0.75, label='AGAGE Mace Head, Ireland (53.33$^\circ$ N)\n {0:.2f} $\pm$ {1:.2f} %/yr'.format(res[1]/np.nanmean(Cl[inds])*100.0, np.std(slope_b)/np.mean(Cl[inds])*100.0))   #facecolors='none'


        #-----------------------------------------
        #PLT JFJ HRFTIR
        #-----------------------------------------
        
        HCl_JFJ_smth     = mf.smooth(HCl_JFJ,  window_len=24, window='flat')
        axes[2].plot(HCl_Date_JFJ,HCl_JFJ_smth,  linewidth=5, color='blue')

        #axes[2].scatter(Cl_Date_JFJ[len(Cl_Date_JFJ)-len(Cl_JFJ2):],Cl_JFJ2, s=20, color='blue', edgecolors='black', label='HR-FTIR (HCl + ClONO$_2$)')
        #axes[2].scatter(Cl_Date_JFJ,Cl_JFJ_ma, s=20, color='blue', edgecolors='black', label='HR-FTIR (HCl + ClONO$_2$) - Moving average (N=5)')
        #axes[2].scatter(Cl_Date_JFJ,Cl_JFJ_sg, s=20, color='brown', edgecolors='black', label='HR-FTIR (HCl + ClONO$_2$) - Savitzky golay (N=5)')
        #axes[2].scatter(Cl_Date_JFJ[len(Cl_Date_JFJ)-len(Cl_JFJ_Ga):],Cl_JFJ_Ga, s=20, color='blue', edgecolors='black', label='JFJ HR-FTIR (HCl + ClONO$_2$)')
        yearJFJ = [single.year for single in HCl_Date_JFJ]
        yearJFJ = np.asarray(yearJFJ)

        inds = np.where( yearJFJ >= y4trend )[0]

        dateYearFrac = mf.toYearFraction(HCl_Date_JFJ[inds])
        weights      = np.ones_like(dateYearFrac)
        res          = mf.fit_driftfourier(dateYearFrac, HCl_JFJ[inds], weights, 2)
        f_drift, f_fourier, f_driftfourier = res[3:6]

        res_b        = mf.cf_driftfourier(dateYearFrac, HCl_JFJ[inds], weights, 2)
        perc, intercept_b, slope_b, pfourier_b = res_b
        #axes[2].plot(Cl_Date_JFJ[inds],f_drift(dateYearFrac),label='Fitted Anual Trend')
        #axes[2].text(0.02,0.94,"Fitted trend -- slope: {0:.3E} ({1:.3f}%)".format(res[1],res[1]/np.mean(totClmn)*100.0),transform=ax1.transAxes)

        print "Fitted trend JFJ -- slope: {0:.3E} ({1:.3f}%)".format(res[1],res[1]/np.mean(HCl_JFJ[inds])*100.0)
        #axes[2].text(0.02,0.87,"JFJ FTIR: {0:.2f} %/yr".format(res[1]/np.mean(HCl_JFJ[inds])*100.0), transform=axes[2].transAxes, color='blue', fontsize=12)
        axes[2].scatter(HCl_Date_JFJ,HCl_JFJ, s=20, color='blue', edgecolors='blue', alpha=0.75, label='Jungfraujoch NDACC FTIR (46.55$^\circ$ N)\n {0:.2f} $\pm$ {1:.2f} %/yr'.format(res[1]/np.nanmean(HCl_JFJ[inds])*100.0, np.std(slope_b)/np.mean(HCl_JFJ[inds])*100.0))
        
        #-----------------------------------------
        #PLT MLO HRFTIR
        #-----------------------------------------
        
        HCl_MLO_smth     = mf.smooth(HCl_MLO,  window_len=24, window='flat')
        axes[2].plot(HCl_Date_MLO,HCl_MLO_smth,  linewidth=5, color='navy')
        #axes[2].scatter(HCl_Date_MLO, HCl_MLO_smth, s=20, facecolors='none', edgecolors='blue', label='MLO HR-FTIR (HCl)')
        #axes[2].fill_between(HCl_Date_MLO[len(HCl_Date_MLO)-len(HCl_MLO_smth):],HCl_MLO_smth-HCl_std_MLO_smth,HCl_MLO_smth+HCl_std_MLO_smth,alpha=0.5,color='blue')
        #axes[2].errorbar(HCl_Date_MLO[len(HCl_Date_MLO)-len(HCl_MLO_smth):],HCl_MLO_smth,yerr=HCl_std_MLO_smth, ecolor='blue', linestyle='None')


        yearMLO = [single.year for single in HCl_Date_MLO]
        yearMLO = np.asarray(yearMLO)

        inds = np.where( yearMLO >= y4trend )[0]

        dateYearFrac = mf.toYearFraction(HCl_Date_MLO[inds])
        weights      = np.ones_like(dateYearFrac)
        res          = mf.fit_driftfourier(dateYearFrac, HCl_MLO[inds], weights, 2)
        f_drift, f_fourier, f_driftfourier = res[3:6]

        res_b        = mf.cf_driftfourier(dateYearFrac, HCl_JFJ[inds], weights, 2)
        perc, intercept_b, slope_b, pfourier_b = res_b

        #axes[2].plot(Cl_Date_JFJ[inds],f_drift(dateYearFrac),label='Fitted Anual Trend')
        #axes[2].text(0.02,0.94,"Fitted trend -- slope: {0:.3E} ({1:.3f}%)".format(res[1],res[1]/np.mean(totClmn)*100.0),transform=ax1.transAxes)
        print "Fitted trend MLO -- slope: {0:.3E} ({1:.3f}%)".format(res[1],res[1]/np.mean(HCl_MLO[inds])*100.0)
        #axes[2].text(0.02,0.83,"MLO FTIR: {0:.2f} %/yr".format(res[1]/np.mean(HCl_MLO[inds])*100.0), transform=axes[2].transAxes, color='navy', fontsize=12)
        axes[2].scatter(HCl_Date_MLO,HCl_MLO, s=20, facecolors='none', edgecolors='navy', alpha=0.75, label='Mauna Loa Observatory NDAC FTIR (19.54$^\circ$ N)\n {0:.2f} $\pm$ {1:.2f} %/yr'.format(res[1]/np.nanmean(HCl_MLO[inds])*100.0, np.std(slope_b)/np.mean(HCl_MLO[inds])*100.0))

        axes[2].set_ylim(1.2, 5.3)
        axes[2].set_xlim(dt.date(1982, 1, 1), dt.date(2017, 12, 31))
        #axes[2].legend(prop={'size':10})
        #axes[2].set_ylabel('Total Column [x10$^{15}$ molec/cm$^2$], HR-FTIR (HCl + ClONO$_2$)', color='blue', fontsize=14)
        axes[2].set_ylabel('Total Column [x10$^{15}$ molec/cm$^2$], FTIR (HCl)', color='blue', fontsize=14)
        axes[2].tick_params(axis='y', colors='blue', labelsize=14)
        #axes[2].legend(prop={'size':12}, loc=4)
        axes[0].set_xlabel('Year', fontsize=14)

        fig.autofmt_xdate()
        fig.subplots_adjust(left=0.1, bottom=0.125, right=0.8, top=0.94)


        lines, labels   = axes[0].get_legend_handles_labels()
        lines2, labels2 = axes[1].get_legend_handles_labels()
        lines3, labels3 = axes[2].get_legend_handles_labels()
       
        axes[0].legend(lines + lines2 + lines3, labels + labels2 + labels3, prop={'size':10.5}, loc=2, frameon=False, ncol=2)  #'weight':'bold'

        plt.suptitle('NH Inorganic Chlorine 1983 to 2017 (Trends 1997 to 2017)', fontsize=16  )


      
        if self.pdfsav: 
            self.pdfsav.savefig(fig,dpi=200)
            plt.savefig('/data1/Campaign/Satellite/MLS/Cly_all.pdf') # bbox_inches='tight'
        else:
            plt.draw()
            plt.show(block=False)
            user_input = raw_input('Press any key to exit >>> ')
            sys.exit()

        #----------------------------------------------------------------------------
        #
        #----------------------------------------------------------------------------



        # fig = plt.figure(figsize=(14,10))

        # outer_grid = gridspec.GridSpec(4, 2, wspace=0.15, hspace=0.175)

        # for i, p in enumerate(POI):

        #     ax = plt.Subplot(fig, outer_grid[i])

        #     if self.ReadNCMrgFlg:

        #         indMrg1 = mf.nearestind(p, PresMrg)
        #         indsGood =  np.where(PrfsMrg[:, indMrg1, indLatMrg] > 0.0)[0]

        #         PrfsMrg2        = PrfsMrg[indsGood, indMrg1, indLatMrg] #   np.delete(PrfsMrg,indsGood,axis=0)
        #         StDMrg2         = StDMrg[indsGood, indMrg1, indLatMrg]  # np.delete(StDMrg,indsGood,axis=0)
        #         DatesMrg2       = DatesMrg[indsGood]  #np.delete(DatesMrg,indsGood,axis=0)

        #         #ax.plot(DatesMrg,PrfsMrg[:,indMrg1, indLatMrg],'k.',markersize=4)
        #         ax.plot(DatesMrg2,PrfsMrg2, 'k.', markersize=0)
        #         ax.scatter(DatesMrg2,PrfsMrg2, s=20, facecolors='none', edgecolors='k', label='GOZCARDS')

        #         ax.fill_between(DatesMrg2,PrfsMrg2-StDMrg2,PrfsMrg2+StDMrg2,alpha=0.5)

        #     if self.ReadHDFMLSFlg:
        #         ind1MLS = mf.nearestind(p, PresMLS)
        #         Avg          = mf.mnthlyAvg(PrfsMLS[:,ind1MLS], DatesMLS, dateAxis=1, meanAxis=0)
        #         Dates        = Avg['dates']
        #         AvgData      = Avg['mnthlyAvg']
        #         #stdData      = Avg['mnthlyAvg']
                
        #         #ax.plot(Dates,AvgData,'r.',markersize=4)
        #         #ax.plot(Dates,AvgData,'o', mfc='none', color='r')
        #         ax.plot(Dates,AvgData, 'k.', markersize=0)
        #         ax.scatter(Dates,AvgData, s=20, facecolors='none', edgecolors='r', label='MLS V4.2')
        #         Avg          = mf.mnthlyAvg(StDMLS[:,ind1MLS], DatesMLS, dateAxis=1, meanAxis=0)
        #         stdData      = Avg['mnthlyAvg']

        #         ax.fill_between(Dates,AvgData-stdData,AvgData+stdData, color='r', alpha=0.25)

        #     if i == 0: ax.legend(prop={'size':10})

        #     ax.grid(True,which='both', alpha=0.5)    
        #     ax.tick_params(which='both',labelsize=11)
        #     #ax.set_ylim(ymin=0, ymax=3e-9)
        #     #ax.set_xlim(0,600) 
        #     ax.set_title(str(p) + ' hPa')

        #     if yrsFlg:
        #         ax.xaxis.set_major_locator(yearsLc)
        #         ax.xaxis.set_minor_locator(months)
        #         #ax1.xaxis.set_minor_formatter(DateFormatter('%m'))
        #         ax.xaxis.set_major_formatter(DateFmt)

        #     else:
        #         ax.xaxis.set_major_locator(monthsAll)
        #         ax.xaxis.set_major_formatter(DateFmt)

            
        #     fig.add_subplot(ax)

        # all_axes = fig.get_axes()
        # #show only the outside spines
        # for ax in all_axes:
        #     for sp in ax.spines.values():
        #         sp.set_visible(False)
        #         plt.setp(ax.get_xticklabels(), visible=False)
        #     if ax.is_first_row():
        #         ax.spines['top'].set_visible(True)
        #     if ax.is_last_row():
        #         ax.spines['bottom'].set_visible(True)
        #         plt.setp(ax.get_xticklabels(), visible=True)
        #     if ax.is_first_col():
        #         ax.spines['left'].set_visible(True)
        #     if ax.is_last_col():
        #         ax.spines['right'].set_visible(True)

        # if (2 % 2 == 1): #even

        #     all_axes[-2].spines['bottom'].set_visible(True)
        #     plt.setp(all_axes[-2].get_xticklabels(), visible=True)
        #     all_axes[-2].set_zorder(1)

        # all_axes[-1].set_xlabel('Date', fontsize=14)
        # all_axes[-2].set_xlabel('Date', fontsize=14)

        # fig.text(0.03, 0.5, 'Mixing Ratios ', fontsize=14, va='center', rotation='vertical')
        # #plt.suptitle('{} Vertical Profiles'.format(gasName.upper()), fontsize=16  )
        # #plt.tight_layout(h_pad=0.25) #w_pad=1.75 pad=1.75,
        # fig.autofmt_xdate()
        # fig.subplots_adjust(left=0.075, bottom=0.075, right=0.975, top=0.95)
        
        # if self.pdfsav: 
        #     self.pdfsav.savefig(fig,dpi=200)
        # else:
        #     plt.show(block=False)
            #user_input = raw_input('Press any key to exit >>> ')
            #sys.exit()


        

        #-----------
        #ONE PANEL (ONE Y-AXIS IN PERCENTAGE)
        #-----------
        # fig, ax = plt.subplots(figsize=(12, 7), sharex=True, sharey=True)
       
        # ind1MLS = mf.nearestind(poi_i, PresMLS)
        # indMrg1 = mf.nearestind(poi_i, Pressure)

        # y   = PrfsMrg[:,indMrg1, indLatMrg]/np.nanmax(PrfsMrg[:,indMrg1, indLatMrg]) * 100.
        # ye  = StDMrg[:,indMrg1, indLatMrg]/np.nanmax(PrfsMrg[:,indMrg1, indLatMrg]) * 100.
        
        # #ax.plot(DatesMrg,PrfsMrg[:,indMrg1, indLatMrg],'k.',markersize=4)
        # ax.plot(DatesMrg, y, '.', color='gray', markersize=0)
        # ax.scatter(DatesMrg, y, s=20, color='gray', edgecolors='k', label='GOZCARDS')   #facecolors='none'
        # ax.fill_between(DatesMrg,y-ye, y+ye, color='gray',alpha=0.5)

        # Avg          = mf.mnthlyAvg(PrfsMLS[:,ind1MLS], DatesMLS, dateAxis=1, meanAxis=0)
        # Dates        = Avg['dates']
        # AvgData      = Avg['mnthlyAvg']
        # #stdData      = Avg['mnthlyAvg']
        # ax.set_ylabel('Normalized Scale [%]', fontsize=14)
        # ax.set_xlabel('Year', fontsize=14)

        # y  = AvgData/np.nanmax(PrfsMrg[:,indMrg1, indLatMrg]) * 100.0
        
        # ax.plot(Dates, y, 'k.', markersize=0)
        # ax.scatter(Dates, y, s=20, color='r' , edgecolors='k', label='MLS V4.2')  #facecolors='none'

        # Avg          = mf.mnthlyAvg(StDMLS[:,ind1MLS], DatesMLS, dateAxis=1, meanAxis=0)
        # stdData      = Avg['mnthlyAvg']
        # ye  = stdData/np.nanmax(PrfsMrg[:,indMrg1, indLatMrg]) * 100.0
        # ax.fill_between(Dates,y-ye,y+ye, color='r', alpha=0.25)

        # ax.grid(True,which='both', alpha=0.5)    
        # #ax.legend(prop={'size':12})
     
        # if yrsFlg:
        #     ax.xaxis.set_major_locator(yearsLc)
        #     ax.xaxis.set_minor_locator(months)
        #     ax.xaxis.set_major_formatter(DateFmt)
        #     #axes[0].xaxis.set_major_formatter(DateFmt)

        # else:
        #     ax.xaxis.set_major_locator(monthsAll)
        #     ax.xaxis.set_major_formatter(DateFmt)

        # y  =   Cl/np.nanmax(Cl) * 100.
        # ye =   Cl_e/np.nanmax(Cl) * 100. 

        # ax.plot(DatesAGA,y, linewidth=3, color='green',label='AGAGE (CFC11 + CFC12 + CCl$_4$)')
        # ax.fill_between(DatesAGA,y-ye,y+ye,alpha=0.5,color='green') 
        # #ax.set_ylabel('VMR [ppt], AGAGE (CFC11 + CFC12 + CCl4)', color='green', fontsize=14)
        # #ax.tick_params(axis='y', colors='green')
        # #ax.legend(prop={'size':12}, loc=3)
        # #ax.set_ylim(800, 920)

        # y  =   Cl_JFJ/np.max(Cl_JFJ) * 100.

        # ax.plot(Cl_Date_JFJ,y, 'k.', markersize=0)
        # ax.scatter(Cl_Date_JFJ,y, s=20, color='blue', edgecolors='black', label='HR-FTIR (HCl + ClONO$_2$)')
        # ax.set_ylim(50, 105)
        # #axes[2].legend(prop={'size':10})
        # #ax.set_ylabel('Total Column [x10$^{15}$ molec/cm$^2$], HR-FTIR (HCl + ClONO$_2$)', color='blue', fontsize=14)
        # #ax.tick_params(axis='y', colors='blue')
        # ax.legend(prop={'size':12}, loc=4)
        # ax.tick_params(which='both',labelsize=14)

        # fig.autofmt_xdate()
        # fig.subplots_adjust(left=0.1, bottom=0.1, right=0.975, top=0.95)
      
        # if self.pdfsav: 
        #     self.pdfsav.savefig(fig,dpi=200)
        # else:
        #     plt.draw()
        #     plt.show(block=False)
            #user_input = raw_input('Press any key to exit >>> ')
            #sys.exit()                


        #-----------
        #3 panels
        # #-----------
        # fig, axes = plt.subplots(3, figsize=(12, 10), sharex=True)
       
        # ind1MLS = mf.nearestind(poi_i, PresMLS)
        # indMrg1 = mf.nearestind(poi_i, Pressure)
        
        # #ax.plot(DatesMrg,PrfsMrg[:,indMrg1, indLatMrg],'k.',markersize=4)
        # axes[0].plot(DatesMrg,PrfsMrg[:,indMrg1, indLatMrg], '.', color='gray', markersize=0)
        # axes[0].scatter(DatesMrg,PrfsMrg[:,indMrg1, indLatMrg], s=20, color='gray', edgecolors='k', label='GOZCARDS')   #facecolors='none'
        # axes[0].fill_between(DatesMrg,PrfsMrg[:,indMrg1, indLatMrg]-StDMrg[:,indMrg1, indLatMrg],PrfsMrg[:,indMrg1, indLatMrg]+StDMrg[:,indMrg1, indLatMrg], color='gray',alpha=0.5)

        # Avg          = mf.mnthlyAvg(PrfsMLS[:,ind1MLS], DatesMLS, dateAxis=1, meanAxis=0)
        # Dates        = Avg['dates']
        # AvgData      = Avg['mnthlyAvg']
        # #stdData      = Avg['mnthlyAvg']
        # axes[0].set_ylabel('VMR [ppb]', fontsize=14)
        # axes[0].set_title('HCl (GOZCARDS & MLS V4.2)')
        

        # axes[0].plot(Dates,AvgData, 'k.', markersize=0)
        # axes[0].scatter(Dates,AvgData, s=20, color='r' , edgecolors='k', label='MLS V4.2')  #facecolors='none'
        # Avg          = mf.mnthlyAvg(StDMLS[:,ind1MLS], DatesMLS, dateAxis=1, meanAxis=0)
        # stdData      = Avg['mnthlyAvg']
        # axes[0].fill_between(Dates,AvgData-stdData,AvgData+stdData, color='r', alpha=0.25)

        # axes[0].grid(True,which='both', alpha=0.5)    
        # axes[0].tick_params(which='both',labelsize=14)
        # axes[0].legend(prop={'size':12})

     
        
        # axes[1].plot(DatesAGA,Cl, linewidth=3, color='green',label='AGAGE (CFC11 + CFC12 + CCl$_4$)')
        # axes[1].fill_between(DatesAGA,Cl-Cl_e,Cl+Cl_e,alpha=0.5,color='green') 
        # axes[1].set_ylabel('VMR [ppt]', fontsize=14)
        # axes[1].tick_params(which='both', labelsize=14)
        # axes[1].grid(True,which='both', alpha=0.5)  
        # axes[1].set_title('AGAGE (Cl$_y$ = CFC11 + CFC12 + CCl$_4$)')
        # axes[1].set_ylim(800, 920)
        # #axes[1].legend(prop={'size':10})

        # axes[2].plot(Cl_Date_JFJ,Cl_JFJ, 'k.', markersize=0)
        # axes[2].scatter(Cl_Date_JFJ,Cl_JFJ, s=20, color='blue', edgecolors='black', label='HR-FTIR (HCl + ClONO$_2$)')
        # axes[2].set_ylim(1, 6)
        # #axes[2].legend(prop={'size':10})
        # axes[2].set_ylabel('Total Column\n [x10$^{15}$ molec/cm$^2$]', fontsize=14)
        # #axes[2].tick_params(axis='y', colors='blue')
        # axes[2].tick_params(which='both', labelsize=14)
        # axes[2].grid(True,which='both', alpha=0.5)
        # axes[2].set_title('HR-FTIR (Cl$_y$ = HCl + ClONO$_2$)')
        # axes[2].set_xlabel('Year', fontsize=14)

        # if yrsFlg:
        #     axes[0].xaxis.set_major_locator(yearsLc)
        #     axes[0].xaxis.set_minor_locator(months)
        #     axes[0].xaxis.set_major_formatter(DateFmt)
        #     #axes[0].xaxis.set_major_formatter(DateFmt)

        # else:
        #     axes[0].xaxis.set_major_locator(monthsAll)
        #     axes[0].xaxis.set_major_formatter(DateFmt)
        # #axes.legend(prop={'size':10})

        # fig.autofmt_xdate()
        # fig.subplots_adjust(left=0.1, bottom=0.1, right=0.975, top=0.95)
      
        # if self.pdfsav: 
        #     self.pdfsav.savefig(fig,dpi=200)
        # else:
        #     plt.draw()
        #     plt.show(block=False)
        #     #user_input = raw_input('Press any key to exit >>> ')
            #sys.exit()        



    

