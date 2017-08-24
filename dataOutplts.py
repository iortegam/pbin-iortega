#----------------------------------------------------------------------------------------
# Name:
#        dataOutplts.py
#
# Purpose:
#       This is a collection of classes and functions used for processing and ploting 
#       sfit output
#
#
# Notes:
#
#
# Version History:
#       Created, October, 2013  Eric Nussbaumer (ebaumer@ucar.edu)
#
#
# License:
#    Copyright (c) 2013-2014 NDACC/IRWG
#    This file is part of sfit4.
#
#    sfit4 is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    any later version.
#
#    sfit4 is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with sfit4.  If not, see <http://www.gnu.org/licenses/>
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
from scipy import linspace, polyval, polyfit, sqrt, stats, randn
#from scipy.interpolate import interp1d

#MODIFIED (IVAN)
#import matplotlib
# Force matplotlib to not use any Xwindows backend.
#matplotlib.use('Agg')
#

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
import matplotlib.colorbar as colorbar
import glob




                                            #------------------#
                                            # Define functions #
                                            #------------------#

def nearestDate(daysList, year, month, day=1):
    ''' Finds the nearest date from a list of days based on a given year, month, and day'''
    testDate = dt.date(year, month, day)
    return min( daysList, key=lambda x:abs(x-testDate) )

def nearestind(val,array):
    ''' Returns the index in array for closest value to val'''
    return (np.abs(array-val)).argmin()

def daysList(timeList):
    ''' Finds a unique set of days within a listof dates '''
    dateList = [dt.date(x.year,x.month,x.day) for x in timeList]
    dateList = np.unique(dateList)
    
    return np.sort(dateList)

def sortDict(DataDict,keyval,excld=[]):
    ''' Sort all values of dictionary based on values of one key'''
    base = DataDict[keyval]
    for k in DataDict:
        if k not in excld: DataDict[k] = [y for (x,y) in sorted(zip(base,DataDict[k]))]
    return DataDict

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

def convrtD(rhs):
    '''This identifies numbers specified in scientific notation with d 
       for double precision and converts them to floats'''
    if 'd' in rhs.lower():                               # Test and handle strings containing D (ie 1.0D0)
        mant,trsh,expo = rhs.lower().partition('d')
        try:
            rhs = float(mant)*10**int(expo)
        except ValueError:
            pass
    else:                                                # Handle strings => number
        try:
            rhs = float(rhs)
        except ValueError:
            pass        

    return rhs

def dbFilterUL(fltrParam,lwrVal=False,uprVal=False):
    ''' Filter the output data based on upper an lower value and
        return the corresponding indicies'''
    inds = []

    #-------------------------------------------
    # At least one bound must be set. Check this
    #-------------------------------------------
    if (not lwrVal) and (not uprVal):
        print "Must specify at least one bound in dbFilterUL"
        return

    #----------------------------------------
    # Case where lower value is not specified
    #----------------------------------------
    if not lwrVal:
        for ind,val in enumerate(fltrParam):
            if val <= uprVal: inds.append(ind)

    #----------------------------------------
    # Case where upper value is not specified
    #----------------------------------------
    elif not uprVal:
        for ind,val in enumerate(fltrParam):
            if val <= lwrVal: inds.append(ind)

    #---------------------------------------------
    # Case where lower and upper val are specified
    #---------------------------------------------
    else:
        for ind,val in enumerate(fltrParam):
            if ( val >= lwrVal and val <= uprVal ): inds.append(ind)

    return inds


def dbFilterTF(fltrParam,fltrVal):
    ''' Filter the output data based on True False value
        and return the corresponding indicies'''
    inds = []

    #------------------------------
    # Find cases that match fltrVal
    #------------------------------
    for ind,val in enumerate(fltrParam):
        if val == fltrVal: inds.append(ind)

    return inds


def getintrsct(vals1,vals2):
    ''' Return the indicies of vals1 and vals2 for the intersection '''

    intrsctVals = np.intersect1d(vals1,vals2)
    inds1       = np.nonzero( np.in1d( vals1, intrsctVals ) )[0]
    inds2       = np.nonzero( np.in1d( vals2, intrsctVals ) )[0]

    return (inds1,inds2)


def bias(xData,yData):
    ''' Cacluates the mean value of the residuals. 
    Where xData is observation and yData is the model'''

    #--------------------------------
    # Make sure arrays are numpy type
    #--------------------------------
    yData = np.asarray(yData)
    xData = np.asarray(xData)

    #-------------------------------
    # Make sure arrays are same size
    #-------------------------------
    if ( yData.shape[0] != xData.shape[0] ): 
        print 'Data sets must have same size in corrcoef: xData = {}, yData = {}'.format(xData.shape[0],yData.shape[0])
        sys.exit() 

    biasCalc = sum( yData - xData ) / len(yData) 

    return biasCalc    

def rmse(xData, yData):
    '''Calcuates the RMSE and Persons correlation coefficient, 
       xData = observed data, yData = modeled data'''

    #--------------------------------
    # Make sure arrays are numpy type
    #--------------------------------
    yData = np.asarray(yData)
    xData = np.asarray(xData)

    #-------------------------------
    # Make sure arrays are same size
    #-------------------------------
    if ( yData.shape[0] != xData.shape[0] ): 
        print 'Data sets must have same size in corrcoef: xData = {}, yData = {}'.format(xData.shape[0],yData.shape[0])
        sys.exit() 

    #---------------------------------
    # Find the residual sum of squares 
    #---------------------------------
    SS_res = sum( (yData - xData)**2 )

    #------------------------------
    # Root Mean Square Error (RMSE)
    #------------------------------
    rmse = np.sqrt( SS_res / len(yData) )

    return rmse


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


def readParm(fname):
    ''' Function to read parm.output file '''

    #--------------
    # Read all file
    #--------------
    with open(fname, 'r' ) as fopen:
        lines = fopen.readlines()

    params = lines[3].strip().split()
    data   = np.array([[float(x) for x in row.split()[1:]] for row in lines[4:]])

    #-----------------------------------------------------------
    # Create a dictionary where keys are header of parm.out file
    # A list of numpy arrays is created for repeated keys
    #-----------------------------------------------------------
    parm = {}
    for k in set(params):
        inds = [i for i, val in enumerate(params) if val == k]
        parm.setdefault(k,[]).append(data[:,inds])

    #--------------------------------------
    # Un-nest numpy arrays in parm dictionary
    #--------------------------------------
    for k in parm: parm[k] = parm[k][0]   # Check this

    return parm


def wet_to_dry_vmr(f_wet,f_H2O):
    ''' Convert wet VMR to dry VMR '''
    
    #----------------------------------------
    # Convert wet vmr to dry:
    #   --- f_gas(dry)= f_gas(wet)/ (1-f_H2O)
    #----------------------------------------
    f_dry    = f_wet / (1.0 - f_H2O)    
    
    return f_dry


def lyrVMR(lAlt,hAlt,Z,prf,airmass,sclfct):
    ''' Find layer mixing ratio and convert to ppx (based on scale factor sclfct)'''
    
    Zlwr     = (np.abs(Z-lAlt)).argmin()
    Zupr     = (np.abs(Z-hAlt)).argmin()
    part_VMR = np.asarray( np.average(rprf1[Zlwr:Zupr,:], axis=0, weights=Airmass1[Zlwr:Zupr,:] ) ) * sclfct 
    
    return part_VMR

def readCtlF(ctlF):
    ''' Read .ctl file'''
    
    ctl = {}
    lines = tryopen(ctlF)
    
    if lines:
        gas_flg = True
        
        for line in lines:
            line = line.strip()  

            #---------------------------------------------
            # Lines which start with comments and comments
            # embedded within lines
            #---------------------------------------------                
            if line.startswith('#'): continue                # For lines which start with comments
            if '#' in line: line = line.partition('#')[0]    # For comments embedded in lines

            #-------------------
            # Handle empty lines
            #-------------------
            if len(line) == 0: continue

            #--------------------------
            # Populate input dictionary
            #--------------------------
            if '=' in line:
                lhs,_,rhs = line.partition('=')
                lhs       = lhs.strip()
                rhs       = rhs.strip().split()


                ctl[lhs] = [convrtD(val) for val in rhs]
                #for rhs_part in rhs:                          
                    #if 'd' in rhs_part.lower():                               # Test and handle strings containing D (ie 1.0D0)
                        #mant,trsh,expo = rhs_part.lower().partition('d')
                        #try:
                            #rhs_part = float(mant)*10**int(expo)
                        #except ValueError:
                            #pass
                    #else:                                                # Handle strings => number
                        #try:
                            #rhs_part = [float(ind) for ind in rhs_part]
                        #except ValueError:
                            #pass
                    #ctl[lhs] = rhs

            else:                                                        # Handle multiple line assignments
                rhs = line.strip().split()             

                try:
                    rhs = [float(ind) for ind in rhs]
                except ValueError:
                    pass

                ctl[lhs] += rhs          

            #----------------------    
            # Primary Retrieval Gas
            #----------------------
            # Search for primary retrieval gas in gas.column.list or gas.profile.list
            # The first gas listed is considered the primary retrieval gas
            if gas_flg:
                match = re.match(r'\s*gas\.\w+\.list\s*=\s*(\w+)', line)
                if match:
                    PrimaryGas = match.group(1)
                    gas_flg = False

            #----------------   
            # Spectral Region
            #----------------                
            # Find number of bands
            #match = re.match(r'\s*band\s*=\s*(.+)',line)
            #if match:
                #bands = match.group(1).split()
                #bnd_flg = True

            # Find the upper and lower bands for set including all regions    
            #if bnd_flg:
                #for band in bands:
                    #match = re.match(r'\s*band.' + band + '.nu_\w+\s*=\s*(.*)',line)
                    #if match:
                        #nu.append(float(match.group(1)))

    #nu.sort()              # Sort wavenumbers
    #nuUpper = nu[-1]  # Pull off upper wavenumber
    #nuLower = nu[0]   # Pull off lower wavenumber  
    
    if not gas_flg: return (PrimaryGas,ctl)
    else:           return ctl
    

#---------------------------------
# A simple code for finding linear 
# trend with prediciton intervals
#---------------------------------



#-------------------------------------
# Code to do boot strap trend analysis
#-------------------------------------
def fourier_basis(x, degree, half_period):
    """Returns a 2-d array of fourier basis."""
    A = np.ones((x.size, 2 * degree + 1))
    
    for d in range(1, degree + 1):
        A[:, 2*d-1] = np.cos(d * np.pi * x / half_period)
        A[:, 2*d] = np.sin(d * np.pi * x / half_period)
    
    return A


def fit_driftfourier(x, data, weights, degree, half_period=0.5):
    """
    Fit y = f(x - x.min()) to data where f is given by
    fourier series + drift.
    
    Parameters
    ----------
    x : 1-d array
        x-coordinates
    data : 1-d array
        data values
    weights : 1-d array
        weights (>=0)
    degree : int
        degree of fourier series
    half_period : float
        half period
    
    Returns
    -------
    intercept : float
        intercept at x.min()
    slope : float
        slope (drift) for the normalized data
        (x - x.min())
    pfourier : 1-d array
        Fourier series parameters for the
        normalized data
    f_drift : callable
        Can be used to calculate the drift
        given any (non-normalized) x
    f_fourier : callable
        Can be used to calculate fourier series
    f_driftfourier : callable
        Can be used to calculate drift + fourier
    residual_std : float
        estimated standard deviation of residuals
    A : 2-d array
        matrix of "coefficients"
    
    """
    xmin = x.min()
    xnorm = x - xmin
    
    # coefficient matrix
    A = np.ones((x.size, 2 * degree + 2))
    A[:, 1] = xnorm
    A[:, 2:] = fourier_basis(xnorm, degree, half_period)[:, 1:]
    
    # linear weighted least squares
    results = np.linalg.lstsq(A * weights[:, np.newaxis],
                              data * weights)
    
    params = results[0]
    intercept = params[0]
    slope = params[1]
    pfourier = params[2:]
    
    f_drift = lambda t: slope * (t - xmin) + intercept
    f_fourier = lambda t: np.sum(fourier_basis(t - xmin, degree,
                                               half_period)[:, 1:]
                                 * pfourier[np.newaxis, :],
                                 axis=1) + intercept
    f_driftfourier = lambda t: f_drift(t) + f_fourier(t) - intercept
    
    residual_std = np.sqrt(results[1][0] / (x.size - 2 * degree + 2)) 
    
    return (intercept, slope, pfourier,
            f_drift, f_fourier, f_driftfourier,
            residual_std, A)


def cf_driftfourier(x, data, weights, degree,
                    half_period=0.5, nboot=5000,
                    percentiles=(2.5, 50., 97.5)):
    """
    Calculate confidence intervals for the fitted
    parameters from fourier series + drift modelling,
    using bootstrap resampling.
    
    Parameters
    ----------
    nboot : int
        number of bootstrap replicates
    percentiles : sequence of floats
        percentiles of parameter estimate
        distributions to return 
    
    Returns
    -------
    perc : dict
        percentiles for of each parameter
        distribution
    intercept : 1-d array
        intercept estimates from bootstraped
        datasets.
    slope : 1-d array
        slope estimates
    pfourier : 2-d array
        fourier parameters estimates
    
    See Also
    --------
    :func:`fit_driftfourier`
    """
    
    # 1st fit without bootstraping
    results = fit_driftfourier(x, data, weights,
                               degree, half_period)
    f_driftfourier = results[5]
    A = results[7]
    model = f_driftfourier(x)
    residuals = data - model
    
    # generate bootstrap resamples of residuals
    # and new datasets from these resamples
    boot_dataset = np.empty((x.size, nboot))
    for i in range(nboot):
        resample_i = np.floor(np.random.rand(x.size) * x.size).astype(int)
        resample_residuals = residuals[resample_i]
        boot_dataset[:, i] = model + resample_residuals
    
    # fit all bootstrap datasets
    results_boot = np.linalg.lstsq(A * weights[:, np.newaxis],
                                   boot_dataset * weights[:, np.newaxis])
    
    params_boot = results_boot[0]
    
    # compute percentiles
    perc_boot = np.column_stack(np.percentile(params_boot,
                                              percentiles, axis=1))
    
    perc = {'intercept' : perc_boot[0],
            'slope' : perc_boot[1],
            'pfourier' : perc_boot[2:]}
    
    intercept = params_boot[0]
    slope = params_boot[1]
    pfourier = params_boot[2:]
    
    return perc, intercept, slope, pfourier


def weatherout(loc, datadir, dataFileTag, iyear, fyear):
    ''' Read met data (yearly) data from eol at FL0  (initial files are created with the read_FL0_EOL_data.py'''

    tempList       = []
    rhList         = []
    pressList      = []
    tempDPList     = []
    wdirList       = []
    wspdList       = []
    wmaxList       = []
    datet_weather  = []


    nyears = int(fyear) - int(iyear) + 1

    for i in range(nyears):

        #------------------------
        # Open file and read data
        #------------------------
        #fileName = loc + '_met_data_'+str(iyear + i)+'_'+dataFileTag+'.txt'
        fileName = glob.glob(datadir + loc + '_met_data_'+str(iyear + i)+'_'+dataFileTag+'.txt')

        if not fileName:
            continue
        
        fileName = fileName[0]
        print 'reading met data: %s' %fileName

        f = open(fileName, 'r')
        header1 = f.readline()
       
        temp_i          = []
        rh_i            = []
        press_i         = []
        tempdp_i        = []
        wdir_i          = []
        wspd_i          = []
        wmax_i          = []
        datet_weather_i = []
        
        for line in f:
            line = line.strip()
            columns = line.split()
            
            year = int(columns[0])
            month = int(columns[1])
            day = int(columns[2])
            hour = int(columns[3])
            minute = int(columns[4])
            datet_weather_i.append(dt.datetime(year,month,day,hour,minute))
            temp_i.append(float(columns[5]))
            rh_i.append(float(columns[6]))
            press_i.append(float(columns[7]))
            tempdp_i.append(float(columns[8]))
            wdir_i.append(float(columns[9]))
            wspd_i.append(float(columns[10]))
            wmax_i.append(float(columns[11]))

        tempList.extend(temp_i)
        rhList.extend(rh_i)
        pressList.extend(press_i)
        tempDPList.extend(tempdp_i)
        wdirList.extend(wdir_i)
        wspdList.extend(wspd_i)
        wmaxList.extend(wmax_i)
        datet_weather.extend(datet_weather_i)



    # clmap        = 'jet'
    # cm           = plt.get_cmap(clmap) 
    # fig, ax = plt.subplots(sharex=True)
    # ax.scatter(datet_weather,wspdList, facecolors='white', edgecolors='gray',s=35)
    # plt.show(block=False)
    # user_input = raw_input('Press any key to exit >>> ')
    # sys.exit()  

    return wdirList, wspdList, tempList, rhList, datet_weather
             



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
class ReadOutputData(_DateRange):

    def __init__(self,dataDir,primGas='',ctlF='',iyear=False,imnth=False,iday=False,fyear=False,fmnth=False,fday=False,incr=1):
        
        if (not primGas) and (not ctlF):
            print 'Either primGas or ctlF needs to be specify'
            return False
        
        self.PrimaryGas = primGas            
        self.dirLst     = []
        self.fltrFlg    = False
        
        #------------------------------
        # Set flags to indicate whether 
        # data has been read
        #------------------------------
        self.readPbpFlg               = False
        self.readSpectraFlg           = False
        self.readsummaryFlg           = False
        self.readRefPrfFlg            = False
        self.readStateVecFlg          = False
        self.readErrorFlg             = {}
        self.readErrorFlg['totFlg']   = False
        self.readErrorFlg['sysFlg']   = False
        self.readErrorFlg['randFlg']  = False
        self.readErrorFlg['vmrFlg']   = False
        self.readErrorFlg['avkFlg']   = False
        self.readErrorFlg['KbFlg']    = False
        self.readPrfFlgApr            = {}
        self.readPrfFlgRet            = {}

        #---------------------------------
        # Test if date range given. If so 
        # create a list of directories to 
        # read data from
        #---------------------------------
        if all([iyear,imnth,iday,fyear,fmnth,fday]): 
            
            # check if '/' is included at end of path
            if not( dataDir.endswith('/') ):
                dataDir = dataDir + '/'           
                
            # Check if directory exits
            ckDir(dataDir,exitFlg=True)
            
            self.dirDateTime = []
            self.dirFlg      = True
            
            _DateRange.__init__(self, iyear, imnth, iday, fyear, fmnth, fday, incr=1)
            
            #--------------------------------------------
            # Walk through first level of directories and
            # collect directory names for processing
            #--------------------------------------------
            for drs in os.walk(dataDir).next()[1]: 
                
                #-------------------------------------------
                # Test directory to make sure it is a number
                #-------------------------------------------
                try:    int(drs[0:4])
                except: continue

                if _DateRange.inRange(self, int(drs[0:4]), int(drs[4:6]), int(drs[6:8]) ):
                    #-----------------------------------------------------------------------------------
                    # Some directories have values of 60 for seconds. This throws an error in date time. 
                    # Seconds can only go from 0-59. If this is case reduce seconds by 1. This is a 
                    # temp solution until coadd can be fixed.
                    #-----------------------------------------------------------------------------------
                    if drs[13:] == '60': ss = int(drs[13:]) - 1
                    else:                ss = int(drs[13:]) 
                    self.dirDateTime.append(dt.datetime(int(drs[0:4]), int(drs[4:6]), int(drs[6:8]), int(drs[9:11]), int(drs[11:13]), ss ) )
                    self.dirLst.append(dataDir+drs+'/')            
            
            #----------------------------------------------------------------------------
            # This is important in keeping all the gathered data in order (i.e. profiles,
            # summary, etc.). The sorted directory is carried through and dates are
            # collected in the summary and the profile gather
            #----------------------------------------------------------------------------
            self.dirLst.sort()
            
        #-----------------
        # Single Directory
        #-----------------
        else: 
            self.dirLst = [dataDir]
            self.dirFlg = False
            
        #-----------------------
        # Read ctl File if given
        #-----------------------
        if ctlF: 
            (self.PrimaryGas,self.ctl) = readCtlF(ctlF)     
            #-------------------
            # Construct Gas list
            #-------------------
            self.gasList = []
            if 'gas.profile.list' in self.ctl: self.gasList += self.ctl['gas.profile.list'] 
            if 'gas.column.list' in self.ctl:  self.gasList += self.ctl['gas.column.list']
            if not self.gasList: 
                print 'No gases listed in column or profile list....exiting'
                sys.exit()
                
            self.gasList = filter(None,self.gasList)  # Filter out all empty strings
            self.ngas    = len(self.gasList)   
            
            for gas in self.gasList:
                self.readPrfFlgApr[gas.upper()] = False
                self.readPrfFlgRet[gas.upper()] = False
                
        else:
            self.readPrfFlgApr[self.PrimaryGas] = False
            self.readPrfFlgRet[self.PrimaryGas] = False


    def fltrData(self,gasName,mxrms=1.0,minsza=0.0,mxsza=80.0,minDOF=1.0,maxCHI=2.0,minTC=1.0E15,maxTC=1.0E16,mnthFltr=[1,2,3,4,5,6,7,8,9,10,11,12],
                 rmsFlg=True,tcFlg=True,pcFlg=True,cnvrgFlg=True,szaFlg=False,dofFlg=False,chiFlg=False,tcMinMaxFlg=False,mnthFltFlg=False):
        
        #------------------------------------------
        # If filtering has already been done return
        #------------------------------------------
        if self.fltrFlg: return True

        #--------------------------------
        # Filtering has not yet been done
        #--------------------------------
        self.inds = []

        #------------
        # Filter Data
        #------------
        nobs = len(np.asarray(self.summary[gasName+'_FITRMS']))
        print 'Number of total observations before filtering = {}'.format(nobs)

        #---------------------------------
        # Filter based on specified months
        #---------------------------------
        if mnthFltFlg:
            mnthFltr = np.asarray(mnthFltr)
            rminds   = []
            dates    = np.array(self.summary["date"])
            months   = np.array([day.month for day in dates])
            
            for i,month in enumerate(months):
                if month not in mnthFltr: rminds.append(i)
            
            rminds = np.asarray(rminds)
            print ('Total number observations found outside of specified months = {}'.format(len(rminds)))
            self.inds = np.union1d(rminds, self.inds)
        
        #-----------------------------
        # Find total column amount < 0
        #-----------------------------
        if tcFlg:
            if not gasName+'_RetColmn' in self.summary:
                print 'TotColmn values do not exist...exiting..'
                sys.exit()
                
            indsT =  np.where(np.asarray(self.summary[gasName+'_RetColmn']) <= 0.0)[0]
            print ('Total number observations found with negative total column amount = {}'.format(len(indsT)))
            self.inds = np.union1d(indsT, self.inds)
        
        #---------------------------------------------
        # Find total column amount < minTC and > maxTC
        #---------------------------------------------
        if tcMinMaxFlg:
            if not gasName+'_RetColmn' in self.summary:
                print 'TotColmn values do not exist...exiting..'
                sys.exit()
                
            indsT1 = np.where(np.asarray(self.summary[gasName+'_RetColmn']) < minTC)[0]
            indsT2 = np.where(np.asarray(self.summary[gasName+'_RetColmn']) > maxTC)[0]
            indsT  = np.union1d(indsT1,indsT2)
            print "Total number of observations found with total column < minTotalColumn = {}".format(len(indsT1))
            print "Total number of observations found with total column > maxTotalColumn = {}".format(len(indsT2))
            self.inds = np.union1d(indsT, self.inds)        
                
        #-----------------------------
        # Find data with fit RMS > X
        #-----------------------------
        if rmsFlg:
            if not gasName+'_FITRMS' in self.summary:
                print 'RMS values do not exist...exiting..'
                sys.exit()            
                
            indsT = np.where(np.asarray(self.summary[gasName+'_FITRMS']) >= mxrms)[0]
            print ('Total number observations found above max rms value = {}'.format(len(indsT)))
            self.inds = np.union1d(indsT, self.inds)
            
        #------------------------------
        # Find values above max chi_2_y
        #------------------------------
        if chiFlg:
            if not gasName+"_CHI_2_Y" in self.summary:
                print 'CHI_2_Y values do not exist...exiting..'
                sys.exit()            
                
            indsT = np.where(np.asarray(self.summary[gasName+"_CHI_2_Y"]) >= maxCHI)[0]
            print ('Total number observations found above max chi_2_y value = {}'.format(len(indsT)))
            self.inds = np.union1d(indsT, self.inds)            
                    
        #-----------------------------------
        # Find any partial column amount < 0
        #-----------------------------------
        if pcFlg:
            if not gasName in self.rprfs:
                print 'Profile values do not exist...exiting..'
                sys.exit()   
                
            rprf_neg = np.asarray(self.rprfs[gasName]) <= 0
            indsT = np.where( np.sum(rprf_neg,axis=1) > 0 )[0]
            print ('Total number observations found with negative partial column = {}'.format(len(indsT)))
            self.inds = np.union1d(indsT, self.inds)
            
        #-------------------------------------
        # Find observations with SZA > max SZA
        #-------------------------------------
        if szaFlg:
            if 'sza' not in self.pbp:
                print 'SZA not found.....exiting'
                sys.exit()  
            
            sza_inds1 = np.where(self.pbp['sza'] > mxsza)[0]
            sza_inds2 = np.where(self.pbp['sza'] < minsza)[0]
            sza_inds  = np.union1d(sza_inds1, sza_inds2)
            print 'Total number of observations with SZA greater than {0:} = {1:}'.format(mxsza,len(sza_inds1))
            print 'Total number of observations with SZA less than    {0:} = {1:}'.format(minsza,len(sza_inds2))
            self.inds = np.union1d(sza_inds,self.inds)
        
        #--------------------------
        # Filter data based on DOFs
        #--------------------------
        if dofFlg:
            if not gasName+'_DOFS_TRG' in self.summary:
                print 'DOFs values do not exist...exiting..'
                sys.exit() 
                
            indsT = np.where(np.asarray(self.summary[gasName+'_DOFS_TRG']) < minDOF)[0]
            print ('Total number observations found below minimum DOFs = {}'.format(len(indsT)))
            self.inds = np.union1d(indsT, self.inds)      
        
        #------------------------------------
        # Find retrievals where not converged
        #------------------------------------
        if cnvrgFlg:
            if not gasName+'_CONVERGED' in self.summary:
                print 'Converged values do not exist...exiting..'
                sys.exit()
                
            indsT = np.where( np.asarray(self.summary[gasName+'_CONVERGED']) == 'F')[0]
            print ('Total number observations that did not converge = {}'.format(len(indsT)))    
            self.inds = np.union1d(indsT, self.inds)
    
        self.inds = np.array(self.inds)
        print 'Total number of observations filtered = {}'.format(len(self.inds))
        
        self.fltrFlg = True
        
        if nobs == len(self.inds):
            print '!!!! All observations have been filtered....'
            self.empty = True
            return False
        else: self.empty = False
                
            
    def readStatLyrs(self,fname):
        ''' Reads the station layers file '''
        
        #--------------------
        # Check if file exits
        #--------------------
        ckFile(fname,exitFlg=True)

        #--------------------
        # Read data from file
        #--------------------
        with open(fname,'r') as fopen: lines = fopen.readlines()
        
        #---------------
        # Assign heights
        #---------------
        altTemp  = [float(row.strip().split()[1]) for row in lines[3:] ]
        thckTemp = [float(row.strip().split()[2]) for row in lines[3:] ]
        grthTemp = [float(row.strip().split()[3]) for row in lines[3:] ]
        midTemp  = [float(row.strip().split()[4]) for row in lines[3:] ]
        
        self.alt     = np.asarray(altTemp)
        self.thckAlt = np.asarray(thckTemp)
        self.grthAlt = np.asarray(grthTemp)
        self.midAlt  = np.asarray(midTemp)
        

    def readRefPrf(self,fname=''):
            ''' Reads in reference profile, an input file for sfit4 (raytrace) '''
            self.refPrf = {}
                
            parms = ['ALTITUDE','PRESSURE','TEMPERATURE']
            
            if not fname: fname = 'reference.prf'
            
            #------------------------------------
            # Loop through collected directories
            # self.dirLst has already been sorted
            #------------------------------------
            for indMain,sngDir in enumerate(self.dirLst):
                
                #-----------------------------------------
                # Check for the existance of summary file
                # If does not exist then retrieval was not
                # completed...skip 
                #-----------------------------------------
                if not os.path.isfile(sngDir + 'summary'): continue 
    
                try:
                    with open(sngDir + fname,'r') as fopen: lines = fopen.readlines()
                                    
                    #----------------------------------------
                    # Get Altitude, Pressure, and Temperature
                    # from reference.prf file
                    #----------------------------------------
                    nlyrs  = int(lines[0].strip().split()[1])
                    nlines = int(np.ceil(nlyrs/5.0))
                    
                    for ind,line in enumerate(lines):
                        if any(p in line for p in parms):
                            val = [x for x in parms if x in line][0]
                            self.refPrf.setdefault(val,[]).append([float(x[:-1]) for row in lines[ind+1:ind+nlines+1] for x in row.strip().split()])

                except Exception as errmsg:
                    print errmsg
                    continue
    
            self.readRefPrfFlg = True
            
            #------------------------
            # Convert to numpy arrays
            # and sort based on date
            #------------------------
            for k in self.refPrf:
                self.refPrf[k] = np.asarray(self.refPrf[k])

    def readsummary(self,fname=''):
            ''' Reads in summary data from SFIT output '''
            self.summary = {}
                
            if not fname: fname = 'summary'
            
            #-----------------------------------
            # Loop through collected directories
            #-----------------------------------
            for sngDir in self.dirLst:
    
                try:
                    with open(sngDir + fname,'r') as fopen: lines = fopen.readlines()
                                    
                    #--------------------------------
                    # Get retrieved column amount for 
                    # each gas retrieved
                    #--------------------------------
                    ind1       = [ind for ind,line in enumerate(lines) if 'IRET' in line][0]
                    ngas       = int(lines[ind1-1].strip())
                    indGasName = lines[ind1].strip().split().index('GAS_NAME')
                    indRETCOL  = lines[ind1].strip().split().index('RET_COLUMN')
                    indAPRCOL  = lines[ind1].strip().split().index('APR_COLUMN')
    
                    for i in range(ind1+1,ind1+ngas+1):
                        gasname = lines[i].strip().split()[indGasName]
                        self.summary.setdefault(gasname.upper()+'_RetColmn',[]).append(float(lines[i].strip().split()[indRETCOL]))
                        self.summary.setdefault(gasname.upper()+'_AprColmn',[]).append(float(lines[i].strip().split()[indAPRCOL]))
    
                    #---------------------------------------------------------
                    # Get NPTSB, FOVDIA, and INIT_SNR
                    # Currently set up to read SNR from summary file where the
                    # summary file format has INIT_SNR on the line below IBAND 
                    #---------------------------------------------------------
                    ind2     = [ind for ind,line in enumerate(lines) if 'IBAND' in line][0]  
                    self.nbands = int(lines[ind2-1].strip().split()[0])
                    indNPTSB = lines[ind2].strip().split().index('NPTSB')
                    indFOV   = lines[ind2].strip().split().index('FOVDIA')
                    indSNR   = lines[ind2].strip().split().index('INIT_SNR') - 9         # Subtract 9 because INIT_SNR is on seperate line therefore must re-adjust index
                    indFitSNR= lines[ind2].strip().split().index('FIT_SNR') - 9          # Subtract 9 because INIT_SNR is on seperate line therefore must re-adjust index
                    lend     = [ind for ind,line in enumerate(lines) if 'FITRMS' in line][0] - 1
            
                    for lnum in range(ind2+1,lend,2):
                        band = lines[lnum].strip().split()[0] # Get band number
                        self.summary.setdefault('nptsb_'+band,[]).append(  float( lines[lnum].strip().split()[indNPTSB] ) )
                        self.summary.setdefault('FOVDIA_'+band,[]).append( float( lines[lnum].strip().split()[indFOV]   ) )
                        self.summary.setdefault('SNR_'+band,[]).append(    float( lines[lnum+1].strip().split()[indSNR] ) )       # Add 1 to line number because INIT_SNR exists on next line
                        self.summary.setdefault('FIT_SNR_'+band,[]).append(float( lines[lnum+1].strip().split()[indFitSNR] ) )    # Add 1 to line number because FIT_SNR exists on next line
                        
                    #----------------------------------------------------------------
                    # Get fit rms, chi_y^2, degrees of freedom target, converged flag
                    #----------------------------------------------------------------
                    ind2       = [ind for ind,line in enumerate(lines) if 'FITRMS' in line][0]
                    indRMS     = lines[ind2].strip().split().index('FITRMS')
                    indCHIY2   = lines[ind2].strip().split().index('CHI_2_Y')
                    indDOFtrgt = lines[ind2].strip().split().index('DOFS_TRG')
                    indCNVRG   = lines[ind2].strip().split().index('CONVERGED')
    
                    self.summary.setdefault(self.PrimaryGas.upper()+'_FITRMS'   ,[]).append( float( lines[ind2+1].strip().split()[indRMS]     ) )
                    self.summary.setdefault(self.PrimaryGas.upper()+'_CHI_2_Y'  ,[]).append( float( lines[ind2+1].strip().split()[indCHIY2]   ) )
                    self.summary.setdefault(self.PrimaryGas.upper()+'_DOFS_TRG' ,[]).append( float( lines[ind2+1].strip().split()[indDOFtrgt] ) )
                    self.summary.setdefault(self.PrimaryGas.upper()+'_CONVERGED',[]).append(        lines[ind2+1].strip().split()[indCNVRG]     )                        
                    if self.dirFlg: 
                        dirname = os.path.basename(os.path.normpath(sngDir)) 
                        self.summary.setdefault('date',[]).append( dt.datetime(int(dirname[0:4]), int(dirname[4:6]), int(dirname[6:8]), 
                                                                               int(dirname[9:11]), int(dirname[11:13]), int(dirname[13:])))

                except Exception as errmsg:
                    print errmsg
                    continue
    
            self.readsummaryFlg = True
            #------------------------
            # Convert to numpy arrays
            # and sort based on date
            #------------------------
            for k in self.summary:
                self.summary[k] = np.asarray(self.summary[k])
    
            if self.dirFlg: self.summary = sortDict(self.summary, 'date')
            else:           return self.summary    
            
            
            
    def readprfs(self,rtrvGasList,fname='',retapFlg=1):
        ''' Reads in retrieved profile data from SFIT output. Profiles are given as columns. Each row corresponds to
                to an altitude layer [nLayers,nObservations]
                retapFlg determines whether retrieved profiles (=1) or a priori profiles (=0) are read'''
        self.deflt = {}
        #retrvdAll   = ['Z','ZBAR','TEMPERATURE','PRESSURE','AIRMASS','H2O']   # These profiles will always be read
        retrvdAll   = ['Z','ZBAR','TEMPERATURE','PRESSURE','AIRMASS']   # These profiles will always be read

        if not fname: 
            if   retapFlg == 1: fname = 'rprfs.table'
            elif retapFlg == 0: fname = 'aprfs.table'        


        #--------------------------------------
        # Add user specified retrieved gas list 
        # to standard retrievals
        #--------------------------------------
        #orginalRtrvGasList = rtrvGasList
        #rtrvGasList = [g.upper() for g in rtrvGasList if g.upper() != 'H2O']   # Remove water from gas list since this read by default
        retrvdAll.extend(rtrvGasList)
        
        #-----------------------------------
        # Loop through collected directories
        #-----------------------------------
        for sngDir in self.dirLst:       
            
            try:
                with open(sngDir + fname,'r') as fopen:
                    
                    defltLines = fopen.readlines()        
    
                    #--------------------------------
                    # Get Names of profiles retrieved
                    #--------------------------------
                    defltParm = defltLines[3].strip().split()
    
                    #----------------------------------------
                    # Loop through retrieved profiles to read
                    #----------------------------------------
                    for rtrvdSing in retrvdAll:
                        self.deflt.setdefault(rtrvdSing,[]).append([ float(row.strip().split()[defltParm.index(rtrvdSing.upper())]) for row in defltLines[4:] ] )
    
                    #-------------------------------
                    # Get date and time of retrieval
                    #-------------------------------
                    if self.dirFlg: 
                        dirname = os.path.basename(os.path.normpath(sngDir)) 
                        self.deflt.setdefault('date',[]).append( dt.datetime(int(dirname[0:4]), int(dirname[4:6]), int(dirname[6:8]), 
                                                                             int(dirname[9:11]), int(dirname[11:13]), int(dirname[13:])))

            except Exception as errmsg:
                print errmsg
                continue                            

        #----------------------------
        # Test if Dictionary is empty
        #----------------------------
        if not self.deflt: 
            print 'No Profiles exist...exiting'
            self.empty = True
            return  #sys.exit()
        else: self.empty = False
        
        #-----------------------
        # Convert to numpy array 
        #-----------------------
        for rtrvdSing in retrvdAll:
            #self.deflt[rtrvdSing] = np.transpose( np.asarray( self.deflt[rtrvdSing] ) )
            self.deflt[rtrvdSing] = np.asarray( self.deflt[rtrvdSing] )
  
        if self.dirFlg: self.deflt['date'] = np.asarray( self.deflt['date'] )
  
        #--------------------------------------------------------
        # If retrieved profiles is a gas, get total column amount
        #--------------------------------------------------------
        for gas in rtrvGasList:
            self.deflt[gas+'_tot_col'] = np.sum(self.deflt[gas] * self.deflt['AIRMASS'], axis=1)
  
        #-------------------------------------------
        # Assign to aprfs or rprfs according to flag
        #-------------------------------------------
        if retapFlg == 1: 
            try:
                self.rprfs
                gen = (g for g in retrvdAll if g not in self.rprfs)
                for g in gen: self.rprfs[g] = self.deflt[g]
                
            except: self.rprfs = self.deflt
            for gas in retrvdAll: self.readPrfFlgRet[gas.upper()] = True
            
        elif retapFlg == 0: 
            try:
                self.aprfs
                gen = (g for g in retrvdAll if g not in self.aprfs)
                for g in gen: self.aprfs[g] = self.deflt[g]
                
            except: self.aprfs = self.deflt
            for gas in retrvdAll: self.readPrfFlgApr[gas.upper()] = True
            
        if self.dirFlg: del self.deflt
        else:           return self.deflt        
             
    
    def readStateVec(self,fname=""):
        ''' Read retrieved parameters in state vector (not includinng profiles)'''
        
        if not fname: fname = "statevec"
        self.statevec = {}
        
        #-----------------------------------
        # Loop through collected directories
        #-----------------------------------
        for sngDir in self.dirLst:
    
            try:
                with open(sngDir + fname,'r') as fopen: lines = fopen.readlines()
        
                #--------------------------------------------------
                # Find location of non profile retrieved parameters
                # These are burried near the end of the file
                #--------------------------------------------------
                # Determine number of layers
                #---------------------------
                nlyrs  = int(lines[1].strip().split()[0])
                nlines = int(np.ceil(nlyrs/5.0))
                
                #------------------------------------
                # Determine number of retrieved gases
                #------------------------------------
                nskip = nlines*3+6
                ngas  = int(lines[nskip].strip().split()[0])
                
                #--------------------------------------------------------
                # Finally find number of non profile retrieved parameters
                #--------------------------------------------------------
                nskip += 2*ngas*(3+nlines) + 2
                nparms = int(lines[nskip].strip().split()[0])
                
                #----------------------------------------
                # Get retrieved parameters (not a priori)
                #----------------------------------------
                nlines = int(np.ceil(nparms/5.0))
                parms = []
                vals  = []
                
                for lineNum in range(0,nlines):
                    parms.extend(lines[nskip+lineNum+1].strip().split())
                    skipVal = nskip + lineNum + nlines*3 - 1
                    vals.extend(lines[skipVal].strip().split())
                    
                vals = [float(val) for val in vals]
    
                for key,val in zip(*(parms,vals)):
                    self.statevec.setdefault(key,[]).append(val)
                    
            except Exception as errmsg:
                print errmsg
                continue            
    
        #------------------------
        # Convert to numpy arrays
        # and sort based on date
        #------------------------
        for k in self.statevec:
            self.statevec[k] = np.asarray(self.statevec[k])    
            
        self.readStateVecFlg = True
   
    
    def readError(self,totFlg=True,sysFlg=False,randFlg=False,vmrFlg=False,avkFlg=False,KbFlg=False):
        ''' Reads in error analysis data from Layer1 output '''
        self.error   = {}
        self.Kb      = {}
        self.sysErr  = {}
        self.randErr = {}
        self.sysErrDiag  = {}
        self.randErrDiag = {}
            
        #-----------------------------------
        # Loop through collected directories
        #-----------------------------------
        for sngDir in self.dirLst:          
            #-------------------------------
            # Read error summary information
            #-------------------------------
            try:
                with open(sngDir + 'Errorsummary.output','r') as fopen:
                    errSumLines = fopen.readlines()        

                #---------------------------------------------------------
                # Loop through summary information. Does not include units
                #---------------------------------------------------------
                for line in errSumLines[3:]:
                    header = line.strip().split('=')[0].strip()
                    val    = float(line.strip().split('=')[1].split()[0])
                    self.error.setdefault(header,[]).append( val )

                #-------------------------------
                # Get date and time of retrieval
                #-------------------------------
                if self.dirFlg: 
                    dirname = os.path.basename(os.path.normpath(sngDir)) 
                    self.error.setdefault('date',[]).append( dt.datetime(int(dirname[0:4]), int(dirname[4:6]), int(dirname[6:8]), 
                                                                         int(dirname[9:11]), int(dirname[11:13]), int(dirname[13:])))

            except Exception as errmsg:
                print errmsg
                continue                            

            #-------------------------
            # Read in averaging kernel
            #-------------------------
            if avkFlg:
                try:
                    with open(sngDir + 'avk.output','r') as fopen:
                        lines = fopen.readlines() 

                    for ind,line in enumerate(lines):
                        if 'nrows' in line: nrows = int(lines[2].strip().split('=')[1])

                        elif line.strip() == 'AVK_scale_factor':
                            totRand = np.array( [ [ float(x) for x in line.strip().split() ] for line in lines[ind+1:ind+1+nrows] ] )
                            self.error.setdefault('AVK_scale_factor',[]).append(totRand)      
                        elif line.strip() == 'AVK_VMR':
                            totRand = np.array( [ [ float(x) for x in line.strip().split() ] for line in lines[ind+1:ind+1+nrows] ] )
                            self.error.setdefault('AVK_vmr',[]).append(totRand)                            

                except Exception as errmsg:
                    print errmsg
                    continue                  
                
            #------------------
            # Read in Kb matrix
            #------------------
            if KbFlg:                
                #------------------
                # Read in Kb matrix
                #------------------    
                lines    = tryopen(sngDir + 'kb.output')
                Kb_param = lines[2].strip().split()
                
                #---------------------------------------------------------------
                # Some parameter names in the Kb file are appended by either 
                # micro-window or gas name. If they are appended by micro-window
                # strip this and just keep parameter name so that we may group
                # the micro-windows under one key
                #---------------------------------------------------------------
                for ind,val in enumerate(Kb_param):
                    if len(val.split('_')) == 2:
                        pname,appnd = val.split('_')
                        try:               int(appnd); val = pname
                        except ValueError: pass
                    Kb_param[ind] = val
            
                Kb_unsrt = np.array([[float(x) for x in row.split()] for row in lines[3:]])
            
                #----------------------------------
                # Create a dictionary of Kb columns
                # A list of numpy arrays is created
                # for repeated keys
                #----------------------------------
                for k in set(Kb_param):
                    inds = [i for i, val in enumerate(Kb_param) if val == k]
                    self.Kb.setdefault(k,[]).append(Kb_unsrt[:,inds])
            
                #--------------------------------------
                # Un-nest numpy arrays in Kb dictionary
                #--------------------------------------
                #for k in self.Kb: self.Kb[k] = self.Kb[k][0]   # Check this                

            #------------------------
            # Read in error matricies
            #------------------------
            # Total Error: 2 matricies (systematic and random)
            if totFlg:
                try:
                    with open(sngDir + 'Stotal.output','r') as fopen:
                        lines = fopen.readlines()        

                    #-------------------------------------------------
                    # Read total Systematic and Random Error Matricies
                    #-------------------------------------------------
                    for ind,line in enumerate(lines):
                        if 'nrows' in line: nrows = int(lines[2].strip().split('=')[1])

                        elif line.strip() == 'Random':
                            totRand = np.array( [ [ float(x) for x in line.strip().split() ] for line in lines[ind+1:ind+1+nrows] ] )
                            self.error.setdefault('Total_Random_Error',[]).append(totRand)

                        elif line.strip() == 'Systematic':
                            totSys = np.array( [ [ float(x) for x in line.strip().split() ] for line in lines[ind+1:ind+1+nrows] ] )
                            self.error.setdefault('Total_Systematic_Error',[]).append(totSys)                            

                except Exception as errmsg:
                    print errmsg
                    continue                    

                #--------------------------------------------
                # If VMR flag is set to true, read vmr values
                #--------------------------------------------
                if vmrFlg:
                    try:
                        with open(sngDir + 'Stotal.vmr.output','r') as fopen:
                            lines = fopen.readlines()        

                        #-------------------------------------------------
                        # Read total Systematic and Random Error Matricies
                        #-------------------------------------------------
                        for ind,line in enumerate(lines):
                            if 'nrows' in line: nrows = int(lines[2].strip().split('=')[1])

                            elif line.strip() == 'Random':
                                totRand = np.array( [ [ float(x) for x in line.strip().split() ] for line in lines[ind+1:ind+1+nrows] ] )
                                self.error.setdefault('Total_Random_Error_VMR',[]).append(totRand)

                            elif line.strip() == 'Systematic':
                                totSys = np.array( [ [ float(x) for x in line.strip().split() ] for line in lines[ind+1:ind+1+nrows] ] )
                                self.error.setdefault('Total_Systematic_Error_VMR',[]).append(totSys)                            

                    except Exception as errmsg:
                        print errmsg
                        continue                    

            # Systematic Error 
            if sysFlg:
                try:
                    with open(sngDir + 'Ssystematic.output','r') as fopen:
                        lines = fopen.readlines()        

                    #--------------------------------------
                    # Read total Systematic Error Matricies
                    #--------------------------------------
                    for ind,line in enumerate(lines):
                        if 'nrows' in line: nrows = int(lines[2].strip().split('=')[1])

                        elif not line.strip().startswith('#') and len(line.strip().split()) == 1:
                            errLbl  = line.strip().split()[0].strip()
                            errmtrx = np.array( [ [ float(x) for x in line.strip().split() ] for line in lines[ind+1:ind+1+nrows] ] )
                            self.sysErr.setdefault(errLbl+'_Systematic_Error',[]).append(errmtrx)
                            self.sysErrDiag.setdefault(errLbl+'_Systematic_Error',[]).append(np.diag(errmtrx))

                        else: continue

                except Exception as errmsg:
                    print errmsg
                    continue                    

                #--------------------------------------------
                # If VMR flag is set to true, read vmr values
                #--------------------------------------------
                if vmrFlg:
                    try:
                        with open(sngDir + 'Ssystematic.vmr.output','r') as fopen:
                            lines = fopen.readlines()        

                        #--------------------------------------------
                        # Read total Systematic Error Matricies (VMR)
                        #--------------------------------------------
                        for ind,line in enumerate(lines):
                            if 'nrows' in line: nrows = int(lines[2].strip().split('=')[1])

                            elif not line.strip().startswith('#') and len(line.strip().split()) == 1:
                                errLbl  = line.strip().split()[0].strip()
                                errmtrx = np.array( [ [ float(x) for x in line.strip().split() ] for line in lines[ind+1:ind+1+nrows] ] )
                                self.sysErr.setdefault(errLbl+'_Systematic_Error_VMR',[]).append(errmtrx)
                                self.sysErrDiag.setdefault(errLbl+'_Systematic_Error_VMR',[]).append(np.diag(errmtrx))
                                
                            else: continue                        

                    except Exception as errmsg:
                        print errmsg
                        continue                  

            # Systematic Error 
            if randFlg:
                try:
                    with open(sngDir + 'Srandom.output','r') as fopen:
                        lines = fopen.readlines()        

                    #--------------------------------------
                    # Read total Systematic Error Matricies
                    #--------------------------------------
                    for ind,line in enumerate(lines):
                        if 'nrows' in line: nrows = int(lines[2].strip().split('=')[1])

                        elif not line.strip().startswith('#') and len(line.strip().split()) == 1:
                            errLbl  = line.strip().split()[0].strip()
                            errmtrx = np.array( [ [ float(x) for x in line.strip().split() ] for line in lines[ind+1:ind+1+nrows] ] )
                            self.randErr.setdefault(errLbl+'_Random_Error',[]).append(errmtrx)
                            self.randErrDiag.setdefault(errLbl+'_Random_Error',[]).append(np.diag(errmtrx))

                        else: continue

                except Exception as errmsg:
                    print errmsg
                    continue                    

                #--------------------------------------------
                # If VMR flag is set to true, read vmr values
                #--------------------------------------------
                if vmrFlg:
                    try:
                        with open(sngDir + 'Srandom.vmr.output','r') as fopen:
                            lines = fopen.readlines()        

                        #--------------------------------------------
                        # Read total Systematic Error Matricies (VMR)
                        #--------------------------------------------
                        for ind,line in enumerate(lines):
                            if 'nrows' in line: nrows = int(lines[2].strip().split('=')[1])

                            elif not line.strip().startswith('#') and len(line.strip().split()) == 1:
                                errLbl  = line.strip().split()[0].strip()
                                errmtrx = np.array( [ [ float(x) for x in line.strip().split() ] for line in lines[ind+1:ind+1+nrows] ] )
                                self.randErr.setdefault(errLbl+'_Random_Error_VMR',[]).append(errmtrx)
                                self.randErrDiag.setdefault(errLbl+'_Random_Error_VMR',[]).append(np.diag(errmtrx))

                            else: continue                        

                    except Exception as errmsg:
                        print errmsg
                        continue                 

        #----------
        # Set flags
        #----------
        if totFlg:  self.readErrorFlg['totFlg']   = True
        if sysFlg:  self.readErrorFlg['sysFlg']   = True
        if randFlg: self.readErrorFlg['randFlg']  = True
        if vmrFlg:  self.readErrorFlg['vmrFlg']   = True
        if avkFlg:  self.readErrorFlg['avkFlg']   = True
        if KbFlg:   self.readErrorFlg['KbFlg']    = True        

        #-----------------------------------
        # Convert date values to numpy array 
        #-----------------------------------
        for k in self.error:
            if k =='date': self.error[k] = np.asarray(self.error[k])        
    

    def readPbp(self,fname=''):
        ''' Reads data from the pbp file'''
        self.pbp = {}

        if not fname: fname = 'pbpfile'
        
        #-----------------------------------
        # Loop through collected directories
        #-----------------------------------
        for sngDir in self.dirLst:      
            
            #-------------------------------------
            # Catch error for opening/reading file
            #-------------------------------------    
            try:
                with open(sngDir + fname,'r') as fopen: lines = fopen.readlines()   

            except Exception as errmsg:
                print errmsg
                continue       
            
            #--------------------
            # Get Number of bands
            #--------------------
            #nbands = int(lines[1].strip().split()[1])
            nbands = int(lines[1].strip().split()[0])
            
            #-------------------------------------------------------
            # Loop through bands. Header of first band is on line[3]
            #-------------------------------------------------------
            lstart = 3
            for i in range(1,nbands+1):
                #--------------------------
                # Read header for each band
                #--------------------------
                # Only take the SZA from first micro-window
                if i == 1: self.pbp.setdefault('sza',[]).append( float(lines[lstart].strip().split()[0])/ 1000.0 )
                nspac = float(lines[lstart].strip().split()[1])
                npnts = int(lines[lstart].strip().split()[2])
                iWN   = float(lines[lstart].strip().split()[3])
                fWN   = float(lines[lstart].strip().split()[4])

                #--------------------------------------------------------------------
                # Determine the number of lines for each band. From lines 363 and 364 
                # of writeout.f90, 12 spectral points are written for each line
                #--------------------------------------------------------------------
                nlines = int( math.ceil( float(lines[lstart].strip().split()[2])/ 12.0) )
                mw     = lines[lstart].strip().split()[6]
                if not 'wavenumber_'+mw in self.pbp: self.pbp.setdefault('wavenumber_'+mw,[]).append(np.linspace(iWN,fWN,num=npnts))
                
                #---------------------------------------------------------------------
                # Read in observed, fitted, and difference spectra for particular band
                #---------------------------------------------------------------------
                self.pbp.setdefault('Obs_'+mw,[]).append([float(x) for row in lines[lstart+1:(lstart+1+nlines*3):3] for x in row.strip().split()[1:] ])
                self.pbp.setdefault('Fitted_'+mw,[]).append([float(x) for row in lines[lstart+2:(lstart+2+nlines*3):3] for x in row.strip().split() ])
                self.pbp.setdefault('Difference_'+mw,[]).append([float(x) for row in lines[lstart+3:(lstart+3+nlines*3):3] for x in row.strip().split() ])
                
                #----------------------------------------
                # Set new line number start for next band
                #----------------------------------------
                lstart += nlines*3+2
                                      
            #-------------------------------
            # Get date and time of retrieval
            #-------------------------------
            if self.dirFlg: 
                dirname = os.path.basename(os.path.normpath(sngDir)) 
                self.pbp.setdefault('date',[]).append( dt.datetime(int(dirname[0:4]), int(dirname[4:6]), int(dirname[6:8]), 
                                                                   int(dirname[9:11]), int(dirname[11:13]), int(dirname[13:])))                        
    
       # if self.dirFlg: self.pbp = sortDict(self.pbp, 'date')
        
        #------------------------
        # Convert to numpy arrays
        # and sort based on date
        #------------------------
        for k in self.pbp:
            self.pbp[k] = np.asarray(self.pbp[k])    
            
        self.readPbpFlg = True
        
        if not self.dirFlg: return self.pbp
    

#------------------------------------------------------------------------------------------------------------------------------        
class PlotData(ReadOutputData):

    def __init__(self,dataDir,ctlF,iyear=False,imnth=False,iday=False,fyear=False,fmnth=False,fday=False,incr=1,outFname=''):
        primGas = ''
        #------------------------------------------------------------
        # If outFname is specified, plots will be saved to this file,
        # otherwise plots will be displayed to screen
        #------------------------------------------------------------
        if outFname: self.pdfsav = PdfPages(outFname)
        else:        self.pdfsav = False
        
        super(PlotData,self).__init__(dataDir,primGas,ctlF,iyear,imnth,iday,fyear,fmnth,fday,incr)

        
    def closeFig(self):
        self.pdfsav.close()
    
               
    def pltPrf(self,fltr=False,minSZA=0.0,maxSZA=80.0,maxRMS=1.0,minDOF=1.0,maxCHI=2.0,minTC=1.0E15,maxTC=1.0E16,dofFlg=False,rmsFlg=True,tcFlg=True,mnthFltr=[1,2,3,4,5,6,7,8,9,10,11,12],
               pcFlg=True,cnvrgFlg=True,allGas=True,sclfct=1.0,sclname='ppv',pltStats=True,szaFlg=False,errFlg=False,chiFlg=False,tcMMflg=False,mnthFltFlg=False):
        ''' Reading retrieved profiles '''
        
        
        aprPrf          = {}   
        rPrf            = {}                     #Deifine Retrieval profile in VMR
        sys_cmpnts_vmr  = {}
        rand_cmpnts_vmr = {}


        localGasList = [self.PrimaryGas]
        
        #-------------------------------------------------------
        # Get profile, summary for Primary gas.... for filtering
        #-------------------------------------------------------
        if not self.readPrfFlgRet[self.PrimaryGas]: self.readprfs([self.PrimaryGas],retapFlg=1)   # Retrieved Profiles
        if self.empty: return False
        if not self.readPrfFlgApr[self.PrimaryGas]: self.readprfs([self.PrimaryGas],retapFlg=0)   # Apriori Profiles
        aprPrf[self.PrimaryGas] = np.asarray(self.aprfs[self.PrimaryGas][0,:]) * sclfct
        rPrf[self.PrimaryGas]   = np.asarray(self.rprfs[self.PrimaryGas]) * sclfct
        rPrfMol                 = np.asarray(self.rprfs[self.PrimaryGas]) * np.asarray(self.rprfs['AIRMASS'])
        airmass                 = np.asarray(self.rprfs['AIRMASS'])
        if len(self.dirLst) > 1: dates                   = self.rprfs['date'] 

        alt = np.asarray(self.rprfs['Z'][0,:])
        
        if not self.readsummaryFlg: self.readsummary()                                      # Summary File info
        rms     = np.asarray(self.summary[self.PrimaryGas+'_FITRMS'])
        dofs    = np.asarray(self.summary[self.PrimaryGas+'_DOFS_TRG'])
        totClmn = np.asarray(self.summary[self.PrimaryGas.upper()+'_RetColmn'])
        

        if not self.readPbpFlg: self.readPbp()                                              # Pbp file info
        sza   = self.pbp['sza']
                 
        if errFlg:                                    # Error info
            
            if not all((self.readErrorFlg['totFlg'],self.readErrorFlg['sysFlg'],self.readErrorFlg['randFlg'])):
                self.readError(totFlg=True,sysFlg=True,randFlg=True,vmrFlg=False,avkFlg=False,KbFlg=False) 
            
            npnts    = np.shape(self.error['Total_Random_Error'])[0]
            nlvls    = np.shape(alt)[0]
            rand_err = np.zeros((npnts,nlvls))
            sys_err  = np.zeros((npnts,nlvls))

           
            for i in range(npnts):
                rand_err[i,:] = np.diag(self.error['Total_Random_Error'][i][:,:])
                sys_err[i,:]  = np.diag(self.error['Total_Systematic_Error'][i][:,:])
            
            tot_err  = np.sqrt(rand_err + sys_err)          
            rand_err = np.sqrt(rand_err)
            sys_err  = np.sqrt(sys_err)

            
            rand_errvmr = (rand_err/ np.asarray(self.rprfs['AIRMASS'][0:npnts]))*sclfct
            sys_errvmr = (sys_err/ np.asarray(self.rprfs['AIRMASS'][0:npnts]))*sclfct
            tot_errvmr = np.sqrt(rand_errvmr**2 + sys_errvmr**2)

            
            rand_cmpnts = self.randErrDiag
            sys_cmpnts  = self.sysErrDiag

            rand_cmpnts_vmr = dict(self.randErrDiag) #np.array(self.randErrDiag, copy = True) 
            sys_cmpnts_vmr = dict(self.sysErrDiag) #np.array(self.sysErrDiag, copy = True)

            #print np.sqrt(rand_cmpnts['measurement_Random_Error'][0])
            #print rand_err[0]
            #exit()


            for k in sys_cmpnts_vmr:
                 sys_cmpnts_vmr[k] = (np.sqrt(sys_cmpnts_vmr[k])/ np.asarray(self.rprfs['AIRMASS'][0:npnts]))*sclfct
                 
            for k in rand_cmpnts_vmr:
                 rand_cmpnts_vmr[k] = (np.sqrt(rand_cmpnts_vmr[k])/ np.asarray(self.rprfs['AIRMASS'][0:npnts]))*sclfct
    
             
        #--------------------
        # Call to filter data
        #--------------------
        if fltr: self.fltrData(self.PrimaryGas,mxrms=maxRMS,minsza=minSZA,mxsza=maxSZA,minDOF=minDOF,maxCHI=maxCHI,minTC=minTC,maxTC=maxTC,mnthFltr=mnthFltr,
                               dofFlg=dofFlg,rmsFlg=rmsFlg,tcFlg=tcFlg,pcFlg=pcFlg,szaFlg=szaFlg,cnvrgFlg=cnvrgFlg,chiFlg=chiFlg,tcMinMaxFlg=tcMMflg,mnthFltFlg=mnthFltFlg)
        else:    self.inds = np.array([]) 
        
        if self.empty: return False
        
        #--------------------------------------------------------
        # Get profile data for other gases if allGas flag is True
        #--------------------------------------------------------
        if allGas:
            nonPgas = (gas for gas in self.gasList if gas != self.PrimaryGas)

            for gas in nonPgas:
                if not self.readPrfFlgRet[gas.upper()]: self.readprfs([gas.upper()],retapFlg=1)   # Retrieved Profiles
                if not self.readPrfFlgApr[gas.upper()]: self.readprfs([gas.upper()],retapFlg=0)   # Apriori Profiles
                aprPrf[gas.upper()] = np.asarray(self.aprfs[gas.upper()][0,:]) * sclfct
                rPrf[gas.upper()]   = np.asarray(self.rprfs[gas.upper()]) * sclfct
                localGasList.append(gas)

        #----------------------------
        # Remove inds based on filter
        #----------------------------
        nfltr   = len(self.inds)
        rms     = np.delete(rms,self.inds)
        ntot    = len(rms)
        sza     = np.delete(sza,self.inds)
        if len(self.dirLst) > 1: dates   = np.delete(dates,self.inds)
        dofs    = np.delete(dofs,self.inds)
        totClmn = np.delete(totClmn,self.inds)
        rPrfMol = np.delete(rPrfMol,self.inds,axis=0)
        airmass = np.delete(airmass,self.inds,axis=0)
        for gas in rPrf:
            rPrf[gas]  = np.delete(rPrf[gas],self.inds,axis=0)
            
        if errFlg:
            rand_err = np.delete(rand_err,self.inds,axis=0)
            sys_err  = np.delete(sys_err,self.inds,axis=0)  
            tot_err  = np.delete(tot_err,self.inds,axis=0)
            rand_errvmr = np.delete(rand_errvmr,self.inds,axis=0)
            sys_errvmr  = np.delete(sys_errvmr,self.inds,axis=0)  
            tot_errvmr  = np.delete(tot_errvmr,self.inds,axis=0)   

            
            for k in sys_cmpnts:
                sys_cmpnts[k] = np.delete(sys_cmpnts[k],self.inds,axis=0)
                sys_cmpnts_vmr[k] = np.delete(sys_cmpnts_vmr[k],self.inds,axis=0)
                
            for k in rand_cmpnts:
                rand_cmpnts[k] = np.delete(rand_cmpnts[k],self.inds,axis=0)
                rand_cmpnts_vmr[k] = np.delete(rand_cmpnts_vmr[k],self.inds,axis=0)
                
        #---------------------
        # Calculate statistics
        #---------------------
        if len(self.dirLst) > 1:
            prfMean = {gas:np.mean(rPrf[gas],axis=0) for gas in rPrf}
            prfSTD  = {gas:np.std(rPrf[gas],axis=0) for gas in rPrf}
        
        #----------------------------
        # Determine if multiple years
        #----------------------------
        if len(self.dirLst) > 1:
            years = [ singDate.year for singDate in dates]      # Find years for all date entries
            if len(list(set(years))) > 1: yrsFlg = True         # Determine all unique years
            else:                         yrsFlg = False        
        
    
        #------------------------------
        # Plot Errors On (mean) Profile
        #------------------------------
        if errFlg:# and (self.PrimaryGas.lower() == self.PrimaryGas.lower()):   
            if len(self.dirLst) > 1:

                print 'Total Number of significant points = %s ' %(len(rms))
                #-------------------
                # Plot on total mean
                #-------------------
                prf_mean = np.mean(rPrfMol,axis=0)
                prf_std = np.std(rPrfMol,axis=0)

                prfvmr_mean = np.mean(rPrf[self.PrimaryGas],axis=0)
                prfvmr_std = np.std(rPrf[self.PrimaryGas],axis=0)
                
                #----------------------------------
                # Find mean, standar deviation, min and max error components
                #----------------------------------
                rand_mean = np.mean(rand_err,axis=0)
                sys_mean  = np.mean(sys_err,axis=0)
                tot_mean  = np.mean(tot_err,axis=0)

                rand_std = np.std(rand_err,axis=0)
                sys_std  = np.std(sys_err,axis=0)
                tot_std  = np.std(tot_err,axis=0)

                rand_min= np.min(rand_err,axis=0)
                sys_min  = np.min(sys_err,axis=0)
                tot_min  = np.min(tot_err,axis=0)

                rand_max = np.max(rand_err,axis=0)
                sys_max  = np.max(sys_err,axis=0)
                tot_max  = np.max(tot_err,axis=0)

                randvmr_mean = np.mean(rand_errvmr,axis=0)
                sysvmr_mean  = np.mean(sys_errvmr,axis=0)
                totvmr_mean  = np.mean(tot_errvmr,axis=0)

                randvmr_std = np.std(rand_errvmr,axis=0)
                sysvmr_std  = np.std(sys_errvmr,axis=0)
                totvmr_std  = np.std(tot_errvmr,axis=0)

                
                for k in sys_cmpnts:
                    sys_cmpnts_vmr[k] = np.mean(sys_cmpnts_vmr[k], axis=0)
                    sys_cmpnts[k] = np.mean(sys_cmpnts[k], axis=0)

                for k in rand_cmpnts:
                    rand_cmpnts_vmr[k] = np.mean(rand_cmpnts_vmr[k], axis=0)
                    rand_cmpnts[k] = np.mean(rand_cmpnts[k], axis=0)
               
        
            mystr = {'prf_Mean':prf_mean,'prf_std':prf_std,'alt':alt, 
            'rand_mean':rand_mean, 'sys_mean':sys_mean, 'tot_mean':tot_mean,
            'rand_std':rand_std, 'sys_std':sys_std, 'tot_std':tot_std,
            'rand_min':rand_min, 'sys_min':sys_min, 'tot_min':tot_min,
            'rand_max':rand_max, 'sys_max':sys_max, 'tot_max':tot_max,
            'rand_cmpnts':rand_cmpnts, 'rPrfMol':rPrfMol,
            'rand_err':rand_err, 'sys_cmpnts':sys_cmpnts,'sys_err':sys_err,
            'randvmr_mean':randvmr_mean, 'sysvmr_mean':sysvmr_mean, 'totvmr_mean':totvmr_mean,
            'randvmr_std':randvmr_std, 'sysvmr_std':sysvmr_std, 'totvmr_std':totvmr_std,
            'prfvmr_mean':prfvmr_mean, 'prfvmr_std':prfvmr_std,
            'rand_cmpnts_vmr':rand_cmpnts_vmr, 'sys_cmpnts_vmr':sys_cmpnts_vmr, 'airmass':airmass,
            'rPrf':rPrf, 'localGasList':localGasList, 'dates':dates}
            
        else:     
            mystr = {'alt':alt,'airmass':airmass,
            'rPrf':rPrf, 'localGasList':localGasList, 'dates':dates}
                
        return mystr


#---------------------------TOTAL COLUMN INFORMATION----------------------
    def pltTotClmn(self,fltr=False,minSZA=0.0,maxSZA=80.0,maxRMS=1.0,minDOF=1.0,maxCHI=2.0,minTC=1.0E15,maxTC=1.0E16,mnthFltr=[1,2,3,4,5,6,7,8,9,10,11,12],
                   dofFlg=False,errFlg=False,szaFlg=False,sclfct=1.0,sclname='ppv',
                   partialCols=False,cnvrgFlg=True,pcFlg=True,tcFlg=True,rmsFlg=True,chiFlg=False,tcMMflg=False,mnthFltFlg=False,
                   allGas = False):
        ''' Reading Time Series of Total Column '''
        
        print '\nReading Total Column.....\n'

        TCapr          = {}   
        TCret          = {}                     
        localGasList = [self.PrimaryGas]

        
        #------------------------------------------
        # Get Profile and Summary information. Need
        # profile info for filtering
        #------------------------------------------
        if not self.readPrfFlgRet[self.PrimaryGas]: self.readprfs([self.PrimaryGas],retapFlg=1)   # Retrieved Profiles
        if not self.readPrfFlgApr[self.PrimaryGas]: self.readprfs([self.PrimaryGas],retapFlg=0)   # Apriori Profiles
        try:    
            if not self.readPrfFlgApr['H2O']:           self.readprfs(['H2O'],retapFlg=0)             # Apriori H2O Profiles
        except: self.readprfs(['H2O'],retapFlg=0)
        if self.empty: return False
        if not self.readsummaryFlg:                 self.readsummary()                            # Summary File info
        if not self.readPbpFlg:                     self.readPbp()                                # Pbp file info
        sza = self.pbp['sza']
        
        if errFlg:
            if not self.readErrorFlg['totFlg']:
                #-----------------------------------
                # Grab total random error components
                #-----------------------------------
                self.readError(totFlg=True,sysFlg=False,randFlg=False,vmrFlg=False,avkFlg=False,KbFlg=False)
            tempKeys = self.error.keys()
            randErrs = {}
            for k in tempKeys:
                if 'Total random uncertainty' in k: randErrs[k] = np.array(self.error[k])
            
        #--------------------
        # Call to filter data
        #--------------------
        if fltr: self.fltrData(self.PrimaryGas,mxrms=maxRMS,minsza=minSZA,mxsza=maxSZA,minDOF=minDOF,maxCHI=maxCHI,minTC=minTC,maxTC=maxTC,mnthFltr=mnthFltr,
                               dofFlg=dofFlg,rmsFlg=rmsFlg,tcFlg=tcFlg,pcFlg=pcFlg,szaFlg=szaFlg,cnvrgFlg=cnvrgFlg,chiFlg=chiFlg,tcMinMaxFlg=tcMMflg,mnthFltFlg=mnthFltFlg)     
        else:    self.inds = np.array([]) 
        
        if self.empty: return False


        #--------------------------------------------------------
        # Get profile data for other gases if allGas flag is True
        #--------------------------------------------------------
        if allGas:
            nonPgas = (gas for gas in self.gasList if gas != self.PrimaryGas)
            for gas in nonPgas:
                TCapr[gas.upper()]   = np.asarray(self.summary[gas.upper()+'_AprColmn'])
                TCret[gas.upper()]   = np.asarray(self.summary[gas.upper()+'_RetColmn'])
                localGasList.append(gas)
        
        
        #--------------------------
        # Get total column and date
        #--------------------------
        totClmn = np.asarray(self.summary[self.PrimaryGas.upper()+'_RetColmn'])
        totClmnApr = np.asarray(self.summary[self.PrimaryGas.upper()+'_AprColmn'])
        dates   = np.asarray(self.summary['date'])
        rms     = np.asarray(self.summary[self.PrimaryGas.upper()+'_FITRMS'])
        dofs    = np.asarray(self.summary[self.PrimaryGas.upper()+'_DOFS_TRG'])
        chi2y   = np.asarray(self.summary[self.PrimaryGas.upper()+'_CHI_2_Y'])
        nbands  = self.nbands
        Airmass = np.asarray(self.rprfs['AIRMASS'])
        rPrf    = np.asarray(self.rprfs[self.PrimaryGas]) * sclfct
        rPrfDry = np.asarray(self.rprfs[self.PrimaryGas]) / (1.0 - np.asarray(self.aprfs['H2O']))* sclfct
        rPrfMol = np.asarray(self.rprfs[self.PrimaryGas]) * Airmass
        alt     = np.asarray(self.rprfs['Z'][0,:])

        nobstot = len(dates)

        #----------------------------------------- 
        snr    = {}
        fitsnr = {}
        for i in range(1,nbands+1):
            snr[i]     = np.asarray(self.summary['SNR_'+str(i)])
            fitsnr[i]  = np.asarray(self.summary['FIT_SNR_'+str(i)])
        
        #-----------------
        # Get Total Errors 
        #-----------------
        if errFlg:
            tot_rnd = np.array(self.error['Total random uncertainty'])
            tot_sys = np.array(self.error['Total systematic uncertainty'])
            tot_std = np.sqrt(tot_rnd**2 + tot_sys**2)
 
        #---------------------------------
        # Remove data based on filter inds
        #---------------------------------
        totClmn    = np.delete(totClmn,self.inds)
        totClmnApr = np.delete(totClmnApr,self.inds)
        dates      = np.delete(dates,self.inds)
        rms        = np.delete(rms,self.inds)
        sza        = np.delete(sza,self.inds)
        dofs       = np.delete(dofs,self.inds)
        chi2y      = np.delete(chi2y,self.inds)
        rPrf       = np.delete(rPrf,self.inds,axis=0)
        rPrfDry    = np.delete(rPrfDry,self.inds,axis=0)
        rPrfMol    = np.delete(rPrfMol,self.inds,axis=0)
        Airmass    = np.delete(Airmass,self.inds,axis=0)

        #--------------------------------------------------------
        #
        #--------------------------------------------------------
        if allGas:
            nonPgas = (gas for gas in self.gasList if gas != self.PrimaryGas)
            for gas in nonPgas:

                TCapr[gas.upper()]   =  np.delete(TCapr[gas.upper()],self.inds)
                TCret[gas.upper()]   =  np.delete(TCret[gas.upper()],self.inds)

        nobsflt = len(dates)
        
        for i in range(1,nbands+1):
            snr[i]     = np.delete(snr[i],self.inds)
            fitsnr[i]  = np.delete(fitsnr[i],self.inds)
        
        if errFlg: 
            tot_std    = np.delete(tot_std,self.inds)
            tot_rnd    = np.delete(tot_rnd,self.inds)
            tot_sys    = np.delete(tot_sys,self.inds)
        else:
            tot_std = 0.0
            tot_rnd = 0.0
            tot_sys = 0.0

            #for k in randErrs:
            #    randErrs[k] = np.delete(randErrs[k],self.inds)
            #    randErrs[k] = randErrs[k] / totClmn * 100.00  # Convert total random error components to 
            
        #----------------------------
        # Determine if multiple years
        #----------------------------
        years = [ singDate.year for singDate in dates]      # Find years for all date entries
        if len(list(set(years))) > 1: yrsFlg = True         # Determine all unique years
        else:                         yrsFlg = False

        #----------------------------
        # Determine if multiple months
        #----------------------------
        months = [ singDate.month for singDate in dates]      # Find years for all date entries
        if len(list(set(months))) > 1: mntsFlg = True         # Determine all unique years
        else:                          mntsFlg = False

        #----------------------------
        # Determine if multiple days
        #----------------------------
        mdays = [ singDate.day for singDate in dates]      # Find years for all date entries
        if len(list(set(months))) > 1: mdays = True         # Determine all unique years
        else:                          mdays = False

        #RETURN STRUCTURE

        mystr = {'dates':dates,'totClmn':totClmn, 'totClmnApr':totClmnApr,'rms':rms, 
                'sza':sza, 'dofs':dofs, 'chi2y':chi2y,
                'rPrf':rPrf, 'rPrfDry':rPrfDry, 'rPrfMol':rPrfMol,
                'Airmass':Airmass, 'snr':snr, 'fitsnr':fitsnr,
                'tot_std':tot_std, 'tot_rnd':tot_rnd,
                'years':years, 'yrsFlg':yrsFlg, 'alt':alt, 'tot_sys':tot_sys,
                'nobsflt':nobsflt, 'nobstot':nobstot, 'TCapr':TCapr, 'TCret':TCret, 
                'GasList':localGasList, 'mntsFlg':mntsFlg, 'mdays':mdays}
                
        return mystr

#---------------------------TOTAL COLUMN PLOTS OF SEVERAL TRACE GASE----------------------
    def plt_ts_TotClmn(self, ngases, TC, DT, gasname_w, FAT, FATAV, Rate, Rate_e, vmrP, vmrPDry, sumP, sl_fourier, 
                       TC_d, DT_d, FAT_d, FATAV_d, Rate_d, Rate_e_d, TCsd_d, sl_fourier_d,
                       TC_m, DT_m, FAT_m, FATAV_m, Rate_m, Rate_e_m, TCsd_m, sl_fourier_m,
                       vmr_d, FAT_vmr_d, FATAV_vmr_d, Rate_vmr_d, Rate_e_vmr_d, vmrsd_d, sl_fourier_vmr_d,
                       vmr_m, FAT_vmr_m, FATAV_vmr_m, Rate_vmr_m, Rate_e_vmr_m, vmrsd_m, sl_fourier_vmr_m,
                       iyear, imnth, iday, fyear, fmnth,fday, yrsFlg=True, pCols=True, weatherFlag=False,
                       wdir=0.0, wspeed=0.0, temp=0.0, rh=0.0, dt_weather=dt.datetime(2010,1,1)):
        ''' Plot Time Series of Total Column for several gases'''


        #-------------------------------------------------------
        #                          PLOTS
        #-------------------------------------------------------
      
        clmap        = 'jet'
        cm           = plt.get_cmap(clmap)    
        yearsLc      = YearLocator()
        monthsAll    = MonthLocator()
        months       = MonthLocator()
        DateFmt      = DateFormatter('%m\n%Y')
        #DateFmt      = DateFormatter('%Y') 

        print '\nPrinting Plots:\n'

        if ngases >= 3: 
            fsize         = (8, 13)
            bottomadjust  = 0.07
            topadjust     = 0.93
        else:
            fsize         = (8, 9)
            bottomadjust  = 0.09
            topadjust     = 0.91
        
        fig,  ax = plt.subplots(ngases,figsize=fsize, sharex=True)
        fig2, ax2 = plt.subplots(ngases,figsize=fsize, sharex=True)
        fig3, ax3 = plt.subplots(ngases,figsize=fsize, sharex=True)
        fig4, ax4 = plt.subplots(ngases,figsize=fsize, sharex=True)
        fig5, ax5 = plt.subplots(ngases,figsize=fsize, sharex=True)

        doy_weather = toYearFraction(dt_weather)

        for k in range(ngases):

            #-------------------------------------------------------
            #If weather data interpolate to actual data
            #-------------------------------------------------------
            if weatherFlag:
                doy_w = toYearFraction(DT[k])
                wspeed_interp = np.interp(doy_w, doy_weather, wspeed)
                wdir_interp = np.interp(doy_w, doy_weather, wdir)

            #-----------------------------------
            ax[k].plot(DT[k],TC[k], color='k', linestyle='None', marker ='.', markersize=4, label='All data')
            ax[k].scatter(DT_d[k],TC_d[k],facecolors='red', edgecolors='black', s=35, label='Dayly averages')
            #ax[k].scatter(DT[k],TC[k], facecolors='white', edgecolors='gray',s=35)
            #ax[k].errorbar(ds['dates'],ds['totClmn'],yerr=ds['tot_std'],fmt='k.',markersize=4,ecolor='red')
            
            #------Wind right axis----
            #axtw = ax[k].twinx()
            #axtw.plot(DT[k], wdir_interp, 'c')
            #axtw.fill_between(DT[k], 0.0, wdir_interp, color='c', facecolor='c',  alpha=0.35)
            #axtw.set_ylabel('Wind direction [r to N]', color='c')
            #for tl in axtw.get_yticklabels():
            #    tl.set_color('c')
            #--------------------------

            #-------fitted trend
            ax[k].plot(DT_d[k],FAT_d[k],linewidth=2.5)
            ax[k].plot(DT_d[k],FATAV_d[k], linewidth=2.5)
            #ax[k].text(0.02,0.85,"Fitted trend -- slope: {0:.3E} ({1:.2f}%)".format(sl_fourier[k], Rate[k]),transform=ax[k].transAxes, fontsize = 14) 

            
            ax[k].grid(True)
            ax[k].set_ylim([np.min(TC[k])-0.05*np.min(TC[k]), np.max(TC[k])+0.1*np.max(TC[k])])
            if gasname_w[k] == 'CH$_4$': ax3[k].set_ylim(np.min(TC[k])-0.02*np.min(TC[k]), np.max(TC[k])+0.02*np.max(TC[k]))
     
            if yrsFlg:
                ax[k].xaxis.set_major_locator(yearsLc)
                ax[k].xaxis.set_minor_locator(months)
                ax[k].xaxis.set_major_formatter(DateFmt) 
                ax[k].xaxis.set_tick_params(which='major',labelsize=12)
                ax[k].xaxis.set_tick_params(which='minor',labelbottom='off')
                ax[k].set_xlim((dt.date(iyear,imnth, iday), dt.date(fyear,fmnth,fday)))
                if k == ngases-1: ax[k].set_xlabel('Year', fontsize = 16)
            else:
                ax[k].xaxis.set_major_locator(monthsAll)
                ax[k].xaxis.set_major_formatter(DateFmt)
                ax[k].set_xlim((dt.date(iyear,1,1), dt.date(iyear,12,31)))
                ax[k].xaxis.set_minor_locator(AutoMinorLocator())
                if k == ngases-1: ax[k].set_xlabel('Month/Year', fontsize = 16)

            if ngases <=2: 
                ax[k].xaxis.set_tick_params(which='major',labelsize=14)
                ax[k].yaxis.set_tick_params(which='major',labelsize=14)
                ax[k].set_ylabel(gasname_w[k] + ' [ppb$_v$]', fontsize = 16)
                ax[k].xaxis.set_tick_params(which='minor',labelbottom='off')

            else:  
                ax[k].set_title(gasname_w[k])
                fig.text(0.075, 0.5, 'Total Column [molecules$\cdot$cm$^{-2}$]', fontsize = 16,
                          horizontalalignment='right',
                          verticalalignment='center',
                          rotation='vertical',
                          transform=ax[k].transAxes)

            fig.autofmt_xdate()

            fig.suptitle('All data and daily averages', fontsize=18)
            fig.subplots_adjust(bottom=bottomadjust,top=topadjust)
                   
            #-----------------------------------
            ax2[k].plot(DT[k],TC[k], color='k', linestyle='None', marker ='.', markersize=4, label='All data')
            ax2[k].scatter(DT_m[k],TC_m[k],facecolors='red', edgecolors='black', s=35, label= 'Monthly averages')
            
            #-------fitted trend 
            ax2[k].plot(DT_m[k],FAT_m[k],label='Fitted Anual Trend',linewidth=2.5)
            ax2[k].plot(DT_m[k],FATAV_m[k],label='Fitted Anual Trend + intra-annual variability', linewidth=2.5)
            #ax2[k].text(0.02,0.85,"Fitted trend -- slope: {0:.3E} ({1:.2f}%)".format(sl_fourier_m[k], Rate_m[k]),transform=ax2[k].transAxes, fontsize = 14) 
            
            #ax2[k].errorbar(DT_m[k],TC_m[k],yerr=TCsd_m[k],fmt='k.',markersize=0,ecolor='red', capthick=2)
            ax2[k].set_title(gasname_w[k])
            ax2[k].grid(True)
            ax2[k].set_ylim([np.min(TC_m[k])-0.1*np.min(TC_m[k]), np.max(TC_m[k])+0.13*np.max(TC_m[k])])
            ##ax2[k].set_ylabel('Total Column$\cdot$[molecules cm$^{-2}$]')

           #ax2[k].text(0.02,0.9,"Fitted intercept at xmin: {:.3E}".format(res[0]),transform=ax2[k].transAxes)
           # ax2[k].text(0.02,0.86,"STD of residuals: {0:.3E} ({1:.3f}%)".format(res[6],res[6]/np.mean(ds['totClmn'])*100.0),transform=ax2[k].transAxes) 

                
            if yrsFlg:
                #plt.xticks(rotation=45)
                ax2[k].xaxis.set_major_locator(yearsLc)
                ax2[k].xaxis.set_minor_locator(months)
                ax2[k].xaxis.set_major_formatter(DateFmt) 
                ax2[k].xaxis.set_tick_params(which='major',labelsize=12)
                ax2[k].xaxis.set_tick_params(which='minor',labelbottom='off')
                ax2[k].set_xlim((dt.date(iyear,imnth, iday), dt.date(fyear,fmnth,fday)))
                if k == ngases-1: ax2[k].set_xlabel('Year', fontsize = 16)
            else:
                ax2[k].xaxis.set_major_locator(monthsAll)
                ax2[k].xaxis.set_major_formatter(DateFmt)
                ax2[k].set_xlim((dt.date(iyear,1,1), dt.date(iyear,12,31)))
                ax2[k].xaxis.set_minor_locator(AutoMinorLocator())
                if k == ngases-1: ax2[k].set_xlabel('Month/Year', fontsize = 16)

                fig2.autofmt_xdate()

            fig2.text(0.075, 0.5, 'Total Column [molecules$\cdot$cm$^{-2}$]', fontsize = 16,
            horizontalalignment='right',
            verticalalignment='center',
            rotation='vertical',
            transform=ax2[k].transAxes)

            fig2.suptitle('All data and monthly averages $\pm$ standard deviation', fontsize=18)
            fig2.subplots_adjust(bottom=bottomadjust,top=topadjust)

            #------------------------------------
            # Plot time series of partial columns with daily averages
            #------------------------------------
            if pCols:

                for pcol in pCols:
                    ax3[k].plot(DT[k],vmrPDry[k], color='k', linestyle='None', marker ='.', markersize=4, label='All data')
                    ax3[k].scatter(DT_d[k],vmr_d[k],facecolors='red', edgecolors='black', s=35, label='Dayly averages')


                    #-----Fitted trend
                    ax3[k].plot(DT_d[k],FAT_vmr_d[k],linewidth=2.5)
                    ax3[k].plot(DT_d[k],FATAV_vmr_d[k], linewidth=2.5) 
                    
                    ax3[k].grid(True)
                    ax3[k].set_ylim(np.min(vmrPDry[k])-0.1*np.min(vmrPDry[k]), np.max(vmrPDry[k])+0.13*np.max(vmrPDry[k]))
                    if gasname_w[k] == 'CH$_4$': ax3[k].set_ylim(np.min(vmrPDry[k])-0.02*np.min(vmrPDry[k]), np.max(vmrPDry[k])+0.02*np.max(vmrPDry[k]))
                    
                    
                    if yrsFlg:
                    
                        ax3[k].xaxis.set_major_locator(yearsLc)
                        ax3[k].xaxis.set_minor_locator(months)
                        ax3[k].xaxis.set_major_formatter(DateFmt) 
                        ax3[k].xaxis.set_tick_params(which='major',labelsize=12)
                        ax4[k].xaxis.set_tick_params(which='minor',labelbottom='off')
                        ax3[k].set_xlim((dt.date(iyear,imnth, iday), dt.date(fyear,fmnth,fday)))
                        if k == ngases-1: ax3[k].set_xlabel('Year', fontsize = 16)
                    else:
                        ax3[k].xaxis.set_major_locator(monthsAll)
                        ax3[k].xaxis.set_major_formatter(DateFmt)
                        ax3[k].set_xlim((dt.date(iyear,1,1), dt.date(iyear,12,31)))
                        ax3[k].xaxis.set_minor_locator(AutoMinorLocator())
                        if k == ngases-1: ax3[k].set_xlabel('Month/Year', fontsize = 16)

                    if ngases <=2: 
                        ax3[k].xaxis.set_tick_params(which='major',labelsize=14)
                        ax3[k].yaxis.set_tick_params(which='major',labelsize=14)
                        ax3[k].set_ylabel(gasname_w[k] + ' [ppb$_v$]', fontsize = 16)
                        ax3[k].xaxis.set_tick_params(which='minor',labelbottom='off')

                    else:  
                        ax3[k].set_title(gasname_w[k])
                        fig3.text(0.075, 0.5, 'VMR [ppb$_{v}$] Dry Air', fontsize = 16,
                        horizontalalignment='right',
                        verticalalignment='center',
                        rotation='vertical',
                        transform=ax3[k].transAxes)

                    
                    
                    fig3.suptitle('All data and daily averages - weighted VMR ('+str(pCols[0][0]) + '-' + str(pCols[0][1])+ 'km)', fontsize=18)
                    fig3.subplots_adjust(bottom=bottomadjust,top = topadjust)

            #------------------------------------
            # Plot time series of partial columns with monthly averages
            #------------------------------------
            if pCols:

                for pcol in pCols:
                    ax4[k].plot(DT[k],vmrPDry[k], color='k', linestyle='None', marker ='.', markersize=4, label='All data')
                    ax4[k].scatter(DT_m[k],vmr_m[k],facecolors='red', edgecolors='black', s=35, label='Dayly averages')
                    
                    #----fitted trend
                    ax4[k].plot(DT_m[k],FAT_vmr_m[k],label='Fitted Anual Trend',linewidth=2.5)
                    ax4[k].plot(DT_m[k],FATAV_vmr_m[k],label='Fitted Anual Trend + intra-annual variability', linewidth=2.5)

                    ax4[k].set_title(gasname_w[k])
                    ax4[k].grid(True)
                    ax4[k].set_ylim(np.min(vmrPDry[k])-0.1*np.min(vmrPDry[k]), np.max(vmrPDry[k])+0.13*np.max(vmrPDry[k]))
                    if gasname_w[k] == 'CH$_4$': ax4[k].set_ylim(np.min(vmrPDry[k])-0.02*np.min(vmrPDry[k]), np.max(vmrPDry[k])+0.02*np.max(vmrPDry[k]))
                    
                    
                    if yrsFlg:
                    
                        ax4[k].xaxis.set_major_locator(yearsLc)
                        ax4[k].xaxis.set_minor_locator(months)
                        ax4[k].xaxis.set_major_formatter(DateFmt) 
                        ax4[k].xaxis.set_tick_params(which='major',labelsize=12)
                        ax4[k].xaxis.set_tick_params(which='minor',labelbottom='off')
                        ax4[k].set_xlim((dt.date(iyear,imnth, iday), dt.date(fyear,fmnth,fday)))
                        if k == ngases-1: ax4[k].set_xlabel('Year', fontsize = 16)
                    else:
                        ax4[k].xaxis.set_major_locator(monthsAll)
                        ax4[k].xaxis.set_major_formatter(DateFmt)
                        ax4[k].set_xlim((dt.date(iyear,1,1), dt.date(iyear,12,31)))
                        ax4[k].xaxis.set_minor_locator(AutoMinorLocator())
                        if k == ngases-1: ax4[k].set_xlabel('Month/Year', fontsize = 16)

                    if ngases <=2: 
                        ax4[k].xaxis.set_tick_params(which='major',labelsize=14)
                        ax4[k].yaxis.set_tick_params(which='major',labelsize=14)
                        ax4[k].set_ylabel(gasname_w[k] + ' [ppb$_v$]', fontsize = 16)
                        ax4[k].xaxis.set_tick_params(which='minor',labelbottom='off')

                    else:  
                        ax4[k].set_title(gasname_w[k])
                        fig4.text(0.075, 0.5, 'VMR [ppb$_{v}$] Dry Air', fontsize = 16,
                        horizontalalignment='right',
                        verticalalignment='center',
                        rotation='vertical',
                        transform=ax4[k].transAxes)
                        
                    
                    fig4.suptitle('All data and monthly averages - weighted VMR ('+str(pCols[0][0]) + '-' + str(pCols[0][1])+ 'km)', fontsize=18)
                    fig4.subplots_adjust(bottom=bottomadjust,top=topadjust)

            #----------------------------
            # Plot total columns by month
            #----------------------------
            month    = np.array([d.month for d in DT[k]])
            mnthSort = list(set(month))
            mnthMean = np.zeros(len(mnthSort))
            mnthSTD  = np.zeros(len(mnthSort))
        
            for i,m in enumerate(mnthSort):
                inds        = np.where(month == m)[0]
                mnthMean[i] = np.mean(vmrPDry[k][inds])
                mnthSTD[i]  = np.std(vmrPDry[k][inds])   
            
            #ax5[k].plot(mnthSort,mnthMean,'k.',markersize=6)
           
            #ax5[k].errorbar(mnthSort,mnthMean,yerr=mnthSTD, fmt='o', color='white', markersize=7, ecolor='red', capthick=2)
            ax5[k].plot(mnthSort,mnthMean, c='red', marker='o', markersize=10)   # edgecolors='red', facecolors='white',
            ax5[k].fill_between(mnthSort,mnthMean-mnthSTD,mnthMean+mnthSTD,alpha=0.5, facecolor='red', color='0.75')   

            #ax5[k].scatter(mnthSort,mnthMean, facecolors='white', edgecolors='red', alpha=0.75, s=7, label='Dayly averages')     
            ax5[k].grid(True,which='both')
            #ax5[k].set_ylabel('Retrieved Total Column\n[molecules cm$^{-2}$]',multialignment='center')
            if k == ngases-1: ax5[k].set_xlabel('Month', fontsize = 16)
           
            ax5[k].set_xlim((0,13))
            ax5[k].set_xticks(range(1,13))

            if ngases <=2: 
                ax5[k].xaxis.set_tick_params(which='major',labelsize=14)
                ax5[k].yaxis.set_tick_params(which='major',labelsize=14)
                ax5[k].set_ylabel(gasname_w[k] + ' [ppb$_v$]', fontsize = 16)
                ax5[k].xaxis.set_tick_params(which='minor',labelbottom='off')

            else:  
                ax5[k].set_title(gasname_w[k])
                fig5.text(0.065, 0.5, r'VMR [ppb$_{v}$] Dry Air (monthly mean with Standard Deviation)', fontsize = 16,
                horizontalalignment='right',
                verticalalignment='center',
                rotation='vertical',
                transform=ax5[k].transAxes)

            

            fig5.subplots_adjust(bottom=bottomadjust,top=topadjust)



        ind = np.arange(ngases)
        fign, ax = plt.subplots()
        #ax.bar(ind,Rate, width = 0.27, align='center', color = 'r', yerr=Rate_e, ecolor = 'k', label = 'all')
        ax.bar(ind+0.27,Rate_vmr_d, width = 0.27, align='center', color = 'b', yerr=Rate_e_d, ecolor = 'k', label = 'Daily')
        #ax.bar(ind+(0.27*2),Rate_vmr_m, width = 0.27, align='center', color = 'g', yerr=Rate_e_m, ecolor = 'k', label = 'Monthly')
        ax.yaxis.grid(True)
        ax.set_xticks(ind)
        ax.set_ylabel('Annual rate of change (%)', fontsize = 16)
        ax.set_xticklabels(gasname_w)
        ax.set_xticks(ind+0.27)
        ax.set_xlabel('Gas', fontsize = 16)
        ax.axhline(0, color='black', lw=1)
        ax.legend(prop={'size':12})

        if self.pdfsav: 
            self.pdfsav.savefig(fig,dpi=200)
            self.pdfsav.savefig(fig2,dpi=200)
            self.pdfsav.savefig(fig3,dpi=200)
            self.pdfsav.savefig(fig4,dpi=200)
            self.pdfsav.savefig(fig5,dpi=200)
            self.pdfsav.savefig(fign,dpi=200)
           
        else:           
            plt.show(block=False) 

        print('\nFinished Plots.......\n')

    def pltcorr(self, ngases, gasname, TC, DT, gasname_w, iyear, imnth, iday, fyear, fmnth,fday, yrsFlg=True):
        ''' Plot correlation plot of gases versus co'''

            #-----------------------------------
            # Finding the index element of CO
            #-----------------------------------
        for k in range(ngases):
            if gasname[k] == 'co': index_co = k

        if index_co: 
            print 'CO is in the list and its position is:', index_co
            print 'Starting linear correlation with CO'

            #-----------------------------------
            #Starting figure in for loop
            #-----------------------------------

            clmap        = 'jet'
            cm           = plt.get_cmap(clmap)    
            yearsLc      = YearLocator()
            monthsAll    = MonthLocator()
            months       = MonthLocator()
            DateFmt      = DateFormatter('%m\n%Y')
            
            fig, ax = plt.subplots(1,ngases-1,figsize=(18,6), sharex=True)
            fig2, ax2 = plt.subplots(ngases-1,figsize=(8,13), sharex=True)
            p = 0

            #-----------------------------------
            # Finding coincident dates
            #-----------------------------------
            DT_co_set = set(DT[index_co])

            for k in range(ngases):
            
            #-----------------------------------
            #Gases that are not co first then for co
            #-----------------------------------
                if k != index_co:        
                    TC_coi = []
                    DT_coi = []
                
                    for i, item in enumerate(DT[k]):
                        if item in DT_co_set:
                            TC_coi.append(TC[k][i])
                            DT_coi.append(item)

                    DT_set = set(DT[k])
                    TC_co_coi = []
                    for i, item in enumerate(DT[index_co]):
                        if item in DT_set:
                            TC_co_coi.append(TC[index_co][i])

            #-------------------------------------------------------
            #Plot 1
            #-------------------------------------------------------
 
                    ax[p].scatter(TC_co_coi,TC_coi, s =25, c='red', alpha=0.75, edgecolors='k')
                    ax[p].grid(True)
                    if p == 0: ax[p].set_ylabel('Gas Total Column [molecules cm$^{-2}$]',multialignment='center',fontsize = 14)
                    fig.text(0.5, 0.075, 'CO Total Column [molecules$\cdot$cm$^{-2}$]', fontsize = 14,
                    horizontalalignment='center', verticalalignment='center', rotation='horizontal',
                    transform=ax[p].transAxes)
                    ax[p].set_title(gasname_w[k])
                    #if k == 3: ax[k].set_xlabel('CO Total Column\n[molecules cm$^{-2}$]',multialignment='center')
                    #ax[k].set_title('Linear Correlation with CO',multialignment='center')

                    A = np.vstack([TC_co_coi, np.ones(len(TC_co_coi))]).T
                    slope, intercept = np.linalg.lstsq(A, TC_coi)[0]
                    fit = slope*np.array(TC_co_coi) + intercept
                    
                    #ax[p].plot(TC_co_coi, fit, 'k', label='Fitted line')
                    fig.subplots_adjust(bottom=0.15,top=0.9, left=0.05, right=0.97)

                    #-------------------------------------------------------
                    #Plot 2
                    #-------------------------------------------------------
                    
                    ax2[p].plot(DT_coi, np.divide(TC_coi,TC_co_coi), color='k', linestyle='None', marker ='', markersize=4, label='All data')
                    ax2[p].scatter(DT_coi,np.divide(TC_coi,TC_co_coi),facecolors='red', edgecolors='black', s=35, label='Dayly averages')
                    #ax[k].scatter(DT[k],TC[k], facecolors='white', edgecolors='gray',s=35)
                    #ax[k].errorbar(ds['dates'],ds['totClmn'],yerr=ds['tot_std'],fmt='k.',markersize=4,ecolor='red')
            
                    ax2[p].set_title(str(gasname_w[k]+'/CO [ppb$_v$/ppb$_v$]'))
                    ax2[p].axvspan(dt.date(2014, 7,16), dt.date(2014, 8,16), facecolor='g', alpha=0.5)
                    ax2[p].grid(True)
                    #ax[k].set_ylim([np.min(TC[k])-0.05*np.min(TC[k]), np.max(TC[k])+0.1*np.max(TC[k])])
     
                    if yrsFlg:
                        ax2[p].xaxis.set_major_locator(yearsLc)
                        ax2[p].xaxis.set_minor_locator(months)
                        ax2[p].xaxis.set_major_formatter(DateFmt) 
                        ax2[p].xaxis.set_tick_params(which='major',labelsize=12)
                        ax2[p].xaxis.set_tick_params(which='minor',labelbottom='off')
                        ax2[p].set_xlim((dt.date(iyear,imnth, iday), dt.date(fyear,fmnth,fday)))
                        if p == ngases-2: ax2[p].set_xlabel('Year', fontsize = 16)
                    else:
                        ax2[p].xaxis.set_major_locator(monthsAll)
                        ax2[p].xaxis.set_major_formatter(DateFmt)
                        ax2[p].set_xlim((dt.date(iyear,1,1), dt.date(iyear,12,31)))
                        ax2[p].xaxis.set_minor_locator(AutoMinorLocator())
                        if p == ngases-2: ax2[p].set_xlabel('Month/Year', fontsize = 16)

                    fig2.autofmt_xdate()

                    fig2.text(0.035, 0.5, 'Ratios', fontsize = 16,
                    horizontalalignment='right',
                    verticalalignment='center',
                    rotation='vertical',
                    transform=ax2[p].transAxes)
                               
                    fig2.subplots_adjust(bottom=0.07,top=0.93)
                    #fig2.suptitle('Daily ratios ('+str(pCols[0][0]) + '-' + str(pCols[0][1])+ 'km)', fontsize=18)

                    p+=1

                    del(TC_coi)
                    del(DT_coi)
                    del(TC_co_coi)



            
            if self.pdfsav: 
                self.pdfsav.savefig(fig,dpi=200)
                self.pdfsav.savefig(fig2,dpi=200)
            else:           
                plt.show(block=False)  

            print('\nFinished Plots.......\n')

        else:
            'CO is not in the list, exit'
            #DT_common = set(DT_d[0]) & set(DT_d[index_co])

    def pltcorrWind(self, ngases, gasname, TC, DT, gasname_w, weatherFlag=False,
                    wdir=0.0, wspeed=0.0, temp=0.0, rh=0.0, dt_weather=dt.datetime(2010,1,1)):
        ''' Plot correlation plot of gases versus co using wind data'''

            #-----------------------------------
            # Finding the index element of CO
            #-----------------------------------

        for k in range(ngases):
            if gasname[k].lower() == 'co': index_co = k

        if index_co: 
            
            print '\nCO is in the list and its position is:', index_co
            print 'Starting linear correlation with CO'

            #-----------------------------------
            #Starting figure in for loop
            #-----------------------------------
            clmap  = 'jet'
            cm     = plt.get_cmap(clmap)  
            fig, ax = plt.subplots(1,ngases-1, figsize=(18,6), sharex=True)
            if weatherFlag: fig2, ax2 = plt.subplots(4,ngases-1, figsize=(18,12), sharex=True)
            p = 0
            
            wdir_comp = [ [0.0, 90.0],[90.0, 180.0], [180.0, 270.0], [270.0, 360.0]  ]
            wdir_comP_name = ['NE' , 'SE', 'SW', 'NW']
            
            for k in range(ngases):

                Wdir_fix = []
                TC_wdirc = []
                TC_CO_wdirc = []

                if k != index_co:
                    #-------------------------------------------------------
                    #Interpolate CO to gas
                    #-------------------------------------------------------
                    doy_w    = toYearFraction(DT[k])
                    doy_co_w = toYearFraction(DT[index_co])
                    TC_coi  = np.interp(doy_w, doy_co_w, TC[index_co])
                    if weatherFlag:
                        doy_weather   = toYearFraction(dt_weather)
                        wspeed_interp = np.interp(doy_w, doy_weather, wspeed)
                        wdir_interp   = np.interp(doy_w, doy_weather, wdir)
                        
                        #some elements of wind are negative (delete them. they are few points)
                        inds           = np.where(wdir_interp < 0.0)[0]
                        TC_coi2        = np.delete(TC_coi, inds)
                        TC2            = np.delete(TC[k], inds)
                        wdir_interp2   = np.delete(wdir_interp, inds)
                        #print 'elements deleted: %s' %len(inds)

                        for wc in range(len(wdir_comp)): 
                            ind_wdir = np.where((wdir_interp2 >= wdir_comp[wc][0]) & (wdir_interp2 <= wdir_comp[wc][1]))[0]
                            TC_wdirc.append(TC[k][ind_wdir])
                            TC_CO_wdirc.append(TC_coi2[ind_wdir])
                            wdir_interp2[ind_wdir] = (wdir_comp[wc][0] + wdir_comp[wc][1])/2.0


                    #-------------------------------------------------------
                    #Plots all in one plot
                    #-------------------------------------------------------
                    #if weatherFlag: norm = colors.Normalize(vmin=np.nanmin(wdir_interp2), vmax=np.nanmax(wdir_interp2) )
                    bounds = [0.0]
                    bounds.extend([int(wdir_comp[ii][1]) for ii in range(len(wdir_comp))])   #   [wdir_comp[0][0], 90, 180, 270, 360]
                    norm = colors.BoundaryNorm(bounds, cm.N)

                    if weatherFlag: sc = ax[p].scatter(TC_coi2,TC2, s =25, c=wdir_interp2, cmap=cm,norm=norm)
                    else:  ax[p].scatter(TC_coi2,TC2, s =25, c='red', alpha=0.75, edgecolors='k')
                    ax[p].grid(True)
                    if p == 0: ax[p].set_ylabel(r'Gas Total Column [molecules cm$^{-2}$]',multialignment='center',fontsize = 14)
                    #if p == (ngases/2.0): ax[p].set_xlabel(r'CO Total Column [molecules$\cdot$cm$^{-2}$]',multialignment='center',fontsize = 14)
                    
                    fig.text(0.5, 0.075, r'CO Total Column [molecules$\cdot$cm$^{-2}$]', fontsize = 14,
                    horizontalalignment='center', verticalalignment='center', rotation='horizontal', transform=ax[p].transAxes)
                    
                    ax[p].set_title(gasname[k])
                  
                    #-----linear correlation analysis
                    A = np.vstack([TC_coi2, np.ones(len(TC_coi2))]).T
                    slope, intercept = np.linalg.lstsq(A, TC2)[0]
                    fit = slope*np.array(TC_coi2) + intercept
                    #ax[p].plot(TC_coi, fit, 'k', label='Fitted line')

                    fig.subplots_adjust(bottom=0.15,top=0.78, left=0.05, right=0.97)
                    if weatherFlag:
                        cax  = fig.add_axes([0.25, 0.92, 0.5, 0.05])
                        cbar = fig.colorbar(sc, cax=cax, format='%2i', orientation='horizontal',
                        boundaries= bounds )
                        cbar.set_label('Wind direction [relative to North]')
                        #cbar.set_ticklabels(['North', 'East', 'South', 'West'])  # horizontal colorbar 
                    
                    #-------------------------------------------------------
                    #Plots: 4 components of wind
                    #-------------------------------------------------------
                    #if weatherFlag: norm = colors.Normalize(vmin=np.nanmin(wdir_interp2), vmax=np.nanmax(wdir_interp2) )

                    if weatherFlag:
                        for ii in range(len(wdir_comp)): 
                            ax2[ii, p].scatter(TC_CO_wdirc[ii], TC_wdirc[ii], s =25)
                            ax2[ii, p].grid(True)
                            ax2[ii, p].set_xlim( np.min(TC_coi2)-0.15*np.min(TC_coi2), np.max(TC_coi2) + 0.15*np.max(TC_coi2) ) 
                            ax2[ii, p].set_ylim( np.min(TC2) -0.15*np.min(TC2), np.max(TC2) + 0.15*np.max(TC2) )
                            ax2[ii, p].annotate('N =%s'%len(TC_wdirc[ii]), xy=(.95, .89), xycoords='axes fraction', fontsize=20,
                            horizontalalignment='right', verticalalignment='center')
                           
                            #-----linear correlation analysis
                            # A = np.vstack([TC_CO_wdirc[ii], np.ones(len(TC_CO_wdirc[ii]))]).T
                            # slope, intercept = np.linalg.lstsq(A, TC_wdirc[ii])[0]
                            # fit = slope*np.array(TC_CO_wdirc[ii]) + intercept

                            # slope2, intercept2, r_value, p_value, std_err = stats.linregress(TC_CO_wdirc[ii], TC_wdirc[ii])
                            # fit2 = slope2*np.array(TC_CO_wdirc[ii]) + intercept2

                            # ax2[ii, p].plot(TC_CO_wdirc[ii], fit, 'k', label='Fitted line')
                            # ax2[ii, p].plot(TC_CO_wdirc[ii], fit2, 'r', label='Fitted line')

                            #ax3[k].text(0.02,0.85,"Fitted trend -- slope: {0:.3E} ({1:.2f}%)".format(sl_fourier_vmr_d[k], Rate_vmr_d[k]),transform=ax3[k].transAxes, fontsize = 14) 
  

                    
                        ax2[0, p].set_title(gasname_w[k])
                        fig2.text(0.5, 0.04, r'CO Total Column [molecules$\cdot$cm$^{-2}$]', fontsize = 14,
                        horizontalalignment='center', verticalalignment='center', rotation='horizontal')#, transform=ax2[0,p].transAxes)

                        fig2.text(0.04, 0.5, r'Total Column [molecules cm$^{-2}$]', fontsize = 16, horizontalalignment='right',
                        verticalalignment='center', rotation='vertical')#,  transform=ax2.transAxes)
       
                            
                        fig2.subplots_adjust(bottom=0.08,top=0.95, left=0.075, right=0.97)

                    p+=1

                    del(Wdir_fix)
                    del(TC_wdirc)
                    del(TC_CO_wdirc)

              
            if self.pdfsav: 
                self.pdfsav.savefig(fig,dpi=200)
                if weatherFlag: self.pdfsav.savefig(fig2,dpi=200)
            else:           
                plt.show(block=False) 

        else:
            'CO is not in the list, exit'

    def pltqq(self, ngases, gasname, TC, DT,TC_d, DT_d, gasname_w, weatherFlag=False,
                    wdir=0.0, wspeed=0.0, temp=0.0, rh=0.0, dt_weather=dt.datetime(2010,1,1)):
        ''' Plot quntile-quantile plot of gases versus co'''

        #-----------------------------------
        # Finding the index element of CO
        #-----------------------------------

        for k in range(ngases):
            if gasname[k].lower() == 'co': index_co = k

        if index_co: 
            
            print '\nCO is in the list and its position is:', index_co
            print 'Starting linear correlation with CO'

            #-----------------------------------
            #Starting figure in for loop
            #-----------------------------------
            clmap  = 'jet'
            cm     = plt.get_cmap(clmap)  
            fig, ax = plt.subplots(1,ngases-1, figsize=(18,6), sharex=True)
            fig2, ax2 = plt.subplots(1,ngases-1, figsize=(18,6), sharex=True)
            p = 0
            
            wdir_comp = [ [0.0, 90.0], [90.0, 180.0], [180.0, 270.0], [270.0, 360.0]  ]
            wdir_comP_name = ['NE' , 'SE', 'SW', 'NW']
            
            for k in range(ngases):

                Wdir_fix = []
                TC_wdirc = []
                TC_CO_wdirc = []

                if k != index_co:
                
                    #-------------------------------------------------------
                    #Interpolate CO to gas
                    #-------------------------------------------------------
                    doy_w    = toYearFraction(DT[k])
                    doy_co_w = toYearFraction(DT[index_co])
                
                    TC_coi   = np.interp(doy_w, doy_co_w, TC[index_co])

                    doy_d_w    = toYearFraction(DT_d[k])
                    doy_d_co_w = toYearFraction(DT_d[index_co])
                    TC_d_coi   = np.interp(doy_d_w, doy_d_co_w, TC_d[index_co])

                    if weatherFlag:
                        doy_weather   = toYearFraction(dt_weather)
                        wspeed_interp = np.interp(doy_w, doy_weather, wspeed)
                        wdir_interp   = np.interp(doy_w, doy_weather, wdir)
                        
                        #some elements of wind are negative (delete them. they are few points)
                        inds          = np.where(wdir_interp < 0.0)[0]
                        TC_coi        = np.delete(TC_coi, inds)
                        TC_n          = np.delete(TC[k], inds)
                        dt_n          = np.delete(DT[k], inds)
                        wdir_interp   = np.delete(wdir_interp, inds)
                        #print 'elements deleted: %s' %len(inds)

                        for wc in range(len(wdir_comp)): 
                            ind_wdir = np.where((wdir_interp >= wdir_comp[wc][0]) & (wdir_interp <= wdir_comp[wc][1]))[0]
                            TC_wdirc.append(TC[k][ind_wdir])
                            TC_CO_wdirc.append(TC_coi[ind_wdir])
                            wdir_interp[ind_wdir] = (wdir_comp[wc][0] + wdir_comp[wc][1])/2.0
 
                #----------------------------------------------------------------------
                # Substracting mean daily values to all points
                #----------------------------------------------------------------------

                    TC_q    = []
                    TC_co_q = []

                    daysall = []
                    for da in dt_n:
                        daysall.append(dt.date(da.year, da.month, da.day))

                    for i, item in enumerate(DT_d[k]):

                        diff = np.asarray(daysall) - item
                        inds = np.where( diff == dt.timedelta(0) )[0]
                        #print item, i, gasname_w[k], len(inds), len(TC_n)                 
                        TC_q.append( (TC_n[inds] - TC_d[k][i])/TC_d[k][i] *100.0 )
                        TC_co_q.append(  (TC_coi[inds] - TC_d_coi[i])/TC_d_coi[i]*100.0  )

                    TC_q_all    =[]
                    TC_co_q_all = []
                    
                    for r, tt in enumerate(TC_q):
                        TC_q_all = np.concatenate((tt, TC_q_all))
                        TC_co_q_all = np.concatenate((TC_co_q[r], TC_co_q_all))


                    #bounds = [0.0]
                    #bounds.extend([int(wdir_comp[ii][1]) for ii in range(len(wdir_comp))])   #   [wdir_comp[0][0], 90, 180, 270, 360]
                    #norm = colors.BoundaryNorm(bounds, cm.N)

                    tw = zip(TC_q_all, wdir_interp)
                    tw = sorted(tw)
                    tc_sorted = [t for t, w in tw]
                    wd_sorted = [w for t, w in tw]

                    Ra_sorted = np.divide(np.asarray(sorted(TC_q_all)), np.asarray(sorted(TC_co_q_all)))
                    #weights = np.ones_like(lines)/len(lines)


                    #sc = ax[p].scatter(sorted(TC_co_q_all), tc_sorted, s =25, c=wd_sorted, cmap=cm,norm=norm)
                    ax[p].scatter(sorted(TC_co_q_all), sorted(TC_q_all), s =25, c='red', alpha=0.75, edgecolors='k')
                    ax[p].grid(True)
                    if p == 0: ax[p].set_ylabel(r'$\%$ Change',multialignment='center',fontsize = 14)
                    fig.text(0.5, 0.075, r'$\%$ Change in CO', fontsize = 14,
                    horizontalalignment='center', verticalalignment='center', rotation='horizontal', transform=ax[p].transAxes)
                    ax[p].set_title(gasname_w[k])
                    fig.subplots_adjust(bottom=0.15,top=0.78, left=0.05, right=0.97)


                    #ax2[p].hist(Ra_sorted, 20)
                    ax2[p].plot(Ra_sorted, 'k')
                    #ax2[p].hist(Ra_sorted, normed=True, histtype='stepfilled', alpha=0.2)
                    #ax2[p].grid(True)
                    #if p == 0: ax[p].set_ylabel(r'$\%$ Change',multialignment='center',fontsize = 14)
                    #fig2.text(0.5, 0.075, r'$\%$ Change in CO', fontsize = 14,
                    #horizontalalignment='center', verticalalignment='center', rotation='horizontal', transform=ax2[p].transAxes)
                    #ax2[p].set_title(gasname_w[k])
                    fig2.subplots_adjust(bottom=0.15,top=0.78, left=0.05, right=0.97)

                    p+=1


            if self.pdfsav: 
                self.pdfsav.savefig(fig,dpi=200)
            
            else:           
                plt.show(block=False)
  
            #stats.probplot( sorted(TC_q_all), dist="norm", plot=plt)
            #plt.show() 

                

    def pltAvk(self,fltr=False,minSZA=0.0,maxSZA=80.0,maxRMS=1.0,minDOF=1.0,maxCHI=2.0,minTC=1.0E15,maxTC=1.0E16,mnthFltr=[1,2,3,4,5,6,7,8,9,10,11,12],
               dofFlg=False,errFlg=False,szaFlg=False,partialCols=False,cnvrgFlg=True,pcFlg=True,tcFlg=True,rmsFlg=True,chiFlg=False,tcMMflg=False,mnthFltFlg=False,
               sclfct=1.0,sclname=''):
        ''' Plot Averaging Kernel. Only for single retrieval '''
        
        print '\nPlotting Averaging Kernel........\n'


        #-------------------------------------------------------
        # Get profile, summary for Primary gas.... for filtering
        #-------------------------------------------------------
        # if not self.readPrfFlgRet[self.PrimaryGas]: self.readprfs([self.PrimaryGas],retapFlg=1)   # Retrieved Profiles
        # if self.empty: return False
        # if not self.readPrfFlgApr[self.PrimaryGas]: self.readprfs([self.PrimaryGas],retapFlg=0)   # Apriori Profiles
        # aprPrf[self.PrimaryGas] = np.asarray(self.aprfs[self.PrimaryGas][0,:]) * sclfct
        # rPrf[self.PrimaryGas]   = np.asarray(self.rprfs[self.PrimaryGas]) * sclfct
        # rPrfMol                 = np.asarray(self.rprfs[self.PrimaryGas]) * np.asarray(self.rprfs['AIRMASS'])
        # if len(self.dirLst) > 1: dates                   = self.rprfs['date'] 

        # alt = np.asarray(self.rprfs['Z'][0,:])
        
        if not self.readsummaryFlg: self.readsummary()                                      # Summary File info
        rms     = np.asarray(self.summary[self.PrimaryGas+'_FITRMS'])
        dofs    = np.asarray(self.summary[self.PrimaryGas+'_DOFS_TRG'])
        totClmn = np.asarray(self.summary[self.PrimaryGas.upper()+'_RetColmn'])
        

        if not self.readPbpFlg: self.readPbp()                                              # Pbp file info
        sza   = self.pbp['sza']


        
        if fltr: self.fltrData(self.PrimaryGas,mxrms=maxRMS,minsza=minSZA,mxsza=maxSZA,minDOF=minDOF,maxCHI=maxCHI,minTC=minTC,maxTC=maxTC,mnthFltr=mnthFltr,
                               dofFlg=dofFlg,rmsFlg=rmsFlg,tcFlg=tcFlg,pcFlg=pcFlg,szaFlg=szaFlg,cnvrgFlg=cnvrgFlg,chiFlg=chiFlg,tcMinMaxFlg=tcMMflg,mnthFltFlg=mnthFltFlg)
        else:    self.inds = np.array([])
        
        if self.empty: return False    
            
        #-------------------------------------------------
        # Determine if AVK has been created via sfit4 core
        # code or via python error analysis
        #-------------------------------------------------     
        if errFlg:   # Read AVK from error output
            if not any((self.readErrorFlg['vmrFlg'],self.readErrorFlg['avkFlg'])):
                self.readError(totFlg=False,sysFlg=False,randFlg=False,vmrFlg=True,avkFlg=True,KbFlg=False)
            
            if not self.error: 
                print 'No Error output files found for AVK plot...exiting..'
                sys.exit()   
                
            #---------------------
            # Get averaging kernel
            #---------------------   
            if len(self.dirLst) > 1:
                avkSCF  = np.delete(np.asarray(self.error['AVK_scale_factor']),self.inds,axis=0)
                avkVMR  = np.delete(np.asarray(self.error['AVK_vmr']),self.inds,axis=0)                    
                avkSCF  = np.mean(avkSCF,axis=0)    
                avkVMR  = np.mean(avkVMR,axis=0)              
            else:
                avkSCF  = self.error['AVK_scale_factor'][0][:,:]
                avkVMR  = self.error['AVK_vmr'],axis=0[0][:,:]        
                
            dofs    = np.trace(avkSCF)
            dofs_cs = np.cumsum(np.diag(avkSCF)[::-1])[::-1]            
            
        else:        # Read AVK from sfit4 output (only contains scaled AVK)
            avkSCF = []
            for d in self.dirLst:
                lines  = tryopen( d + self.ctl['file.out.ak_matrix'][0])
                if not lines: continue
                avkSCF.append(np.array( [ [ float(x) for x in line.strip().split() ] for line in lines[2:] ] ))
                
            if not self.readPrfFlgApr[self.PrimaryGas]: self.readprfs([self.PrimaryGas],retapFlg=0)   # Apriori Profiles
            
            if len(avkSCF) == 1: 
                avkSCF      = avkSCF[0]
                n_layer     = np.shape(avkSCF)[0]
                Iapriori    = np.zeros((n_layer,n_layer))
                IaprioriInv = np.zeros((n_layer,n_layer))
                np.fill_diagonal(Iapriori,self.aprfs[self.PrimaryGas.upper()])
                np.fill_diagonal(IaprioriInv, 1.0 / (self.aprfs[self.PrimaryGas.upper()]))
                avkVMR      = np.dot(np.dot(Iapriori,avkSCF),IaprioriInv)            
                dofs        = np.trace(avkSCF)
                dofs_cs     = np.cumsum(np.diag(avkSCF)[::-1])[::-1]   
                
            else:                
                avkSCF  = np.asarray(avkSCF)
                nobs    = np.shape(avkSCF)[0]
                n_layer = np.shape(avkSCF)[1]
                avkVMR  = np.zeros((nobs,n_layer,n_layer))
                for obs in range(0,nobs):
                    Iapriori        = np.zeros((n_layer,n_layer))
                    IaprioriInv     = np.zeros((n_layer,n_layer))
                    np.fill_diagonal(Iapriori,self.aprfs[self.PrimaryGas.upper()][obs])
                    np.fill_diagonal(IaprioriInv, 1.0 / (self.aprfs[self.PrimaryGas.upper()][obs]))
                    avkVMR[obs,:,:] = np.dot(np.dot(Iapriori,np.squeeze(avkSCF[obs,:,:])),IaprioriInv)          
                
                avkSCF  = np.delete(avkSCF,self.inds,axis=0)
                avkVMR  = np.delete(avkVMR,self.inds,axis=0)
                avkSCF  = np.mean(avkSCF,axis=0)
                avkVMR  = np.mean(avkVMR,axis=0)
                dofs    = np.trace(avkSCF)
                dofs_cs = np.cumsum(np.diag(avkSCF)[::-1])[::-1]                    

        #-------------
        # Get Altitude
        #-------------
        if not self.readPrfFlgRet[self.PrimaryGas]: self.readprfs([self.PrimaryGas],retapFlg=1)   # Retrieved Profiles
        alt = np.asarray(self.rprfs['Z'][0,:])
        
        #-------------------------------------------------------
        #Return structure with results
        #-------------------------------------------------------

        mystr = {'avkSCF':avkSCF, 'avkVMR':avkVMR,'alt':alt, 'dofs_cs':dofs_cs}
                
        return mystr





        
                  
                
                        