#! /usr/bin/python2.7
#----------------------------------------------------------------------------------------
# Name:
#        pltFRAPPE.py
#
# Purpose:
#       The purpose of this program is to plt FRAPPE
#
# Notes:
#   
#
# Version History:
#       Created, August, 2016  Ivan Ortega (iortega@ucar.edu)

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

from scipy import interpolate

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
from mpl_toolkits.basemap import Basemap
import srtm
#from scipy.spatial import cKDTree
from geographiclib.geodesic import Geodesic
import sys
from dbfread import DBF
from scipy import stats

import PltClass as pc

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

def cmap_discretize(cmap, N):
    """Return a discrete colormap from the continuous colormap cmap.
    
        cmap: colormap instance, eg. cm.jet. 
        N: number of colors.
    
    Example
        x = resize(arange(100), (5,100))
        djet = cmap_discretize(cm.jet, 5)
       imshow(x, cmap=djet)
   """
   
    if type(cmap) == str:
        cmap = plt.get_cmap(cmap)
    colors_i = concatenate((linspace(0, 1., N), (0.,0.,0.,0.)))
    colors_rgba = cmap(colors_i)
    indices = linspace(0, 1., N+1)
    cdict = {}
    for ki,key in enumerate(('red','green','blue')):
        cdict[key] = [ (indices[i], colors_rgba[i-1,ki], colors_rgba[i,ki]) for i in xrange(N+1) ]
    # Return colormap object.
    return colors.LinearSegmentedColormap(cmap.name + "_%d"%N, cdict, 1024)

                                    #----------------------------#
                                    #                            #
                                    #        --- Main---         #
                                    #                            #
                                    #----------------------------#

def main(argv):

    #-------------------------------------------------
    # Initializations for C-130 insitu data
    #-------------------------------------------------
    dataDir     = '/data1/ancillary_data/fl0/FRAPPE/C130/'
    
    #-------------------------------------------------
    # Date range to Read C-130 insitu data (ICART FILES)
    #-------------------------------------------------
    iyear      = 2014
    imnth      = 7
    iday       = 25 
    fyear      = 2014
    fmnth      = 8
    fday       = 20

    pltFile    = '/data/iortega/results/fl0/C130_FRAPPE.pdf'
    #-------------------------------------------------
    #                 FLAGS
    #-------------------------------------------------
    saveFlg         = False
    
    #--------------FLAGS TO READ DATA-----------------
    RFFlag          = True              #FLIGHT INFORMATION (ALWAYS TRUE TO INTERPOLATE LAT AND LON)
    
    C2H6Flag        = False              #C2H6
    COFlag          = False              #CO
    H2COFlag        = False              #H2CO
    NH3Flag         = False              #H2CO

    #--------------FLAGS TO PLOT DATA-----------------
    pltC2H6toCO     = False              #Plt Ratio
    pltNH3toCO      = False              #Plt Ratio
    pltH2COtoCO     = False              #Plt Ratio

    pltMap          = True              #ALWAYS TRUE (IF AT LEAST ONE MAP IS CREATED)
   
    pltMapRF        = True             #Map of all RF
    pltMapC2H6      = False             #Map of all C2H6
    pltMapCO        = False             #Map of all CO
    pltMapH2CO      = False             #Map of all h2co
    pltMapNH3       = False             #Map of all NH3
    
    pltMapC2H6toCO  = False             #Map of the Ratio
    pltMapNH3toCO   = False             #Map of the Ratio
    pltMapH2COtoCO  = False             #Map of the Ratio

    #---------------------------------------------------
    # LOCATIONS for AIR MASSES
    #---------------------------------------------------
    locP = [ [-104.0, 40.5] , [-104.6, 39.2], [-106.75, 40.0]]
    maxD = 75.0          #Maximum Distance in Km
    clr  = ['red', 'blue',  'green' ]
    IDP  = ['O&NG/Feedlot', 'Urban', 'Background']
   
    lonB, latB = -105.245, 40.035      #LOCATION OF NCAR FOOTHILLS LAB in Boulder

    #------------------------------------------------------------------------------------------------------
    #
    #                                               START
    #
    #------------------------------------------------------------------------------------------------------

    #---------------------------------------------------
    # INITIALIZE A CLASS FOR MAPS
    #---------------------------------------------------
    if pltMap: mp = pc.MapClass(origin=(40.0,-105.0) , maxAlt=4600, minAlt=500., DistLeft=250000,  DistRight=350000,  DistTop=250000,  DistBottom=200000, saveFlg=saveFlg)

    #-------------------------------------------------
    #                 START
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
    dirlist      = []
    dirlistCO    = []
    dirlistC2H6  = []
    dirlistH2CO  = []
    dirlistNH3   = []

    for dd in dateList:
    	# Find year month and day strings
        yrstr   = "{0:02d}".format(dd.year)
        mnthstr = "{0:02d}".format(dd.month)
        daystr  = "{0:02d}".format(dd.day)

        filename = 'FRAPPE-NCAR-LRT-NAV_C130_'+yrstr+mnthstr+daystr+'*.ict'
        icarttfile = glob.glob( dataDir + filename )
        #if not icarttfile: continue
        for k in icarttfile:
            dirlist.append(k)

        filename = 'frappe-CO_C130_'+yrstr+mnthstr+daystr+'*.ict'
        icarttfile = glob.glob( dataDir + filename )
        #if not icarttfile: continue
        for k in icarttfile:
            dirlistCO.append(k)

        filename = 'frappe-C2H6_C130_'+yrstr+mnthstr+daystr+'*.ict'
        icarttfile = glob.glob( dataDir + filename )
        #if not icarttfile: continue
        for k in icarttfile:
            dirlistC2H6.append(k)

        filename = 'frappe-CH2O_C130_'+yrstr+mnthstr+daystr+'*.ict'
        icarttfile = glob.glob( dataDir + filename )
        #if not icarttfile: continue
        for k in icarttfile:
            dirlistH2CO.append(k)

        filename = 'FRAPPE-NH3_C130_'+yrstr+mnthstr+daystr+'*.ict'
        icarttfile = glob.glob( dataDir + filename )
        #if not icarttfile: continue
        for k in icarttfile:
            dirlistNH3.append(k)

    dirlist.sort()
    dirlistCO.sort()
    dirlistC2H6.sort()
    dirlistH2CO.sort()
    dirlistNH3.sort()

    if RFFlag:

        print 'Reading C130 in-situ Files:'
        altC130       = []
        LatC130       = []
        LonC130       = []
        DateC130      = []
        RF            = []

        for indMain,sngDir in enumerate(dirlist):
    	    ckDir(sngDir)
    	    keys, data,dateutc, vnames = mf.read_ICARTT(sngDir)

            DateC    = dateutc
    	    altC     = data[vnames[9]]   
            LatC     = data[vnames[11]]
            LonC     = data[vnames[14]]

            altC    = np.array(altC, dtype=np.float)
            LatC     = np.array(LatC, dtype=np.float)
            LonC     = np.array(LonC, dtype=np.float)
            DateC    = np.array(DateC)

            index    = np.where( (altC >= 0.0) & (LatC >= -100.0)  & (LonC >= -3000.0) )[0]

            altC     = altC[index]
            LatC     = LatC[index]
            LonC     = LonC[index]
            DateC    = DateC[index]

            DateC130.append(DateC)
            altC130.append(altC)
            LatC130.append(LatC)
            LonC130.append(LonC)
            RF.append(indMain + 1)

        altC130=np.array(altC130).flatten()
        LatC130=np.array(LatC130).flatten()
        LonC130=np.array(LonC130).flatten()
        RF=np.array(RF).flatten()
        DateC130=np.array(DateC130).flatten()
        doyC130 = [mf.toYearFraction(d) for d in DateC130]

    if COFlag:
        print '\nReading CO from C130 in-situ Files:'
        COC130        = []
        DateC130_CO   = []

        for indMain,sngDir in enumerate(dirlistCO):
            ckDir(sngDir)
            keys, data,dateutc, vnames = mf.read_ICARTT(sngDir)

            DateC    = dateutc
            CO       = data[vnames[1]]   

            CO       = np.array(CO, dtype=np.float)
            DateC    = np.array(DateC)

            index    = np.where( (CO >= 0.0) )[0]

            CO       = CO[index]
            DateC    = DateC[index]

            DateC130_CO.append(DateC)
            COC130.append(CO)

        COC130=np.array(COC130).flatten()
        DateC130_CO=np.array(DateC130_CO).flatten()
        doyC130_CO = [mf.toYearFraction(d) for d in DateC130_CO]

        LatC130_CO_int  = [np.interp(d, doyC130[i], LatC130[i]) for i, d in enumerate(doyC130_CO)]
        LonC130_CO_int  = [np.interp(d, doyC130[i], LonC130[i]) for i, d in enumerate(doyC130_CO)]
    

    if C2H6Flag:
        print '\nReading C2H6 from C130 in-situ Files:'
        C2H6C130        = []
        DateC130_C2H6   = []

        for indMain,sngDir in enumerate(dirlistC2H6):
            ckDir(sngDir)
            keys, data,dateutc, vnames = mf.read_ICARTT(sngDir)

            DateC    = dateutc
            C2H6     = data[vnames[3]]   

            C2H6     = np.array(C2H6, dtype=np.float)/1000.0  #ppt to ppb
            DateC    = np.array(DateC)

            index    = np.where( (C2H6 >= 0.0) )[0]

            C2H6     = C2H6[index]
            DateC    = DateC[index]

            DateC130_C2H6.append(DateC)
            C2H6C130.append(C2H6)

        C2H6C130=np.array(C2H6C130).flatten()
        DateC130_C2H6=np.array(DateC130_C2H6).flatten()
        doyC130_C2H6 = [mf.toYearFraction(d) for d in DateC130_C2H6]

        LatC130_C2H6_int  = [np.interp(d, doyC130[i], LatC130[i]) for i, d in enumerate(doyC130_C2H6)]
        LonC130_C2H6_int  = [np.interp(d, doyC130[i], LonC130[i]) for i, d in enumerate(doyC130_C2H6)]

    if H2COFlag:
        print '\nReading H2CO from C130 in-situ Files:'
        H2COC130        = []
        DateC130_H2CO   = []

        for indMain,sngDir in enumerate(dirlistH2CO):
            ckDir(sngDir)
            keys, data,dateutc, vnames = mf.read_ICARTT(sngDir)

            DateC    = dateutc
            H2CO     = data[vnames[3]]   

            H2CO     = np.array(H2CO, dtype=np.float)/1000.0  #ppt to ppb
            DateC    = np.array(DateC)

            index    = np.where( (H2CO >= 0.0) )[0]

            H2CO     = H2CO[index]
            DateC    = DateC[index]

            DateC130_H2CO.append(DateC)
            H2COC130.append(H2CO)

        H2COC130=np.array(H2COC130).flatten()
        DateC130_H2CO=np.array(DateC130_H2CO).flatten()
        doyC130_H2CO = [mf.toYearFraction(d) for d in DateC130_H2CO]

        LatC130_H2CO_int  = [np.interp(d, doyC130[i], LatC130[i]) for i, d in enumerate(doyC130_H2CO)]
        LonC130_H2CO_int  = [np.interp(d, doyC130[i], LonC130[i]) for i, d in enumerate(doyC130_H2CO)]

    if NH3Flag:
        print '\nReading NH3 from C130 in-situ Files:'
        NH3C130        = []
        DateC130_NH3   = []

        for indMain,sngDir in enumerate(dirlistNH3):
            ckDir(sngDir)
            keys, data,dateutc, vnames = mf.read_ICARTT(sngDir)

            DateC    = dateutc
            NH3     = data[vnames[3]]   

            NH3     = np.array(NH3, dtype=np.float) 
            DateC    = np.array(DateC)

            index    = np.where( (NH3 >= 0.0) )[0]

            NH3     = NH3[index]
            DateC    = DateC[index]

            DateC130_NH3.append(DateC)
            NH3C130.append(NH3)

        NH3C130=np.array(NH3C130).flatten()
        DateC130_NH3=np.array(DateC130_NH3).flatten()
        doyC130_NH3 = [mf.toYearFraction(d) for d in DateC130_NH3]

        LatC130_NH3_int  = [np.interp(d, doyC130[i], LatC130[i]) for i, d in enumerate(doyC130_NH3)]
        LonC130_NH3_int  = [np.interp(d, doyC130[i], LonC130[i]) for i, d in enumerate(doyC130_NH3)]

    #---------------------------------------------------------------------------------------
    #
    #                                    PLOTS
    #
    #---------------------------------------------------------------------------------------
    if pltC2H6toCO:
        x        = np.concatenate(COC130)
        xDate    = np.concatenate(DateC130_CO)
        xdoy     = mf.toYearFraction(xDate)

        y        = np.concatenate(C2H6C130)
        yDate    = np.concatenate(DateC130_C2H6)
        ydoy     = mf.toYearFraction(yDate)

        la       = np.concatenate(LatC130)
        lo       = np.concatenate(LonC130)
        d        = np.concatenate(DateC130)
        ddoy     = mf.toYearFraction(d)

        y_int     = interpolate.interp1d(ydoy, y,  bounds_error=False)(xdoy)
        la_int    = interpolate.interp1d(ddoy, la,  bounds_error=False)(xdoy)
        lo_int    = interpolate.interp1d(ddoy, lo,   bounds_error=False)(xdoy)

        Ra        = np.true_divide(y_int, x)

        fig, ax = plt.subplots(figsize=(10,6), sharex=True)

        #ax.plot(x, y_int, '.k', markersize=1, linewidth=0, alpha=0.15)
        #ax.scatter(x, y_int, facecolors='white', s=50, color='gray', alpha=0.15)

        x     = x*math.exp(-1.6/7.4)
        y_int = y_int*math.exp(-1.6/7.4)


        slope, intercept, r_value, p_value, std_err = stats.linregress(x, y_int)
        print '\nCorrelation of all data' + ' C2H6 to CO'
        print 'Slope = ', str(slope)
        print 'intercept = ', str(intercept)
        print 'r_value = ', str(r_value)
        print 'std_err =', str(std_err)



        for p, loc in enumerate(locP):
            Dist     = []

            for i, k in enumerate(x):
                c = mf.haversine(abs(loc[0]), abs(loc[1]), abs(lo_int[i]), abs(la_int[i]) )
                Dist.append(float(c))

            Dist = np.asarray(Dist)
            inds = np.where(Dist <= maxD)[0]

            x1 = np.asarray(x[inds])
            y1 = np.asarray(y_int[inds])
            
            ax.scatter(x1, y1, facecolors='white', s=70, color=clr[p], label=IDP[p])

            slope, intercept, r_value, p_value, std_err = stats.linregress(x1, y1)
            print '\nCorrelation of '+ IDP[p] + ' C2H6 to CO'
            print 'Slope = ', str(slope)
            print 'intercept = ', str(intercept)
            print 'r_value = ', str(r_value)
            print 'std_err =', str(std_err)
            print 'std_err =', str(std_err*np.sqrt(len(x1)))

            xx = range(0, 500, 2)
            xx = np.asarray(xx)
            Fit = slope*np.asarray(xx) + (intercept)#*math.exp(-1.6/7.4)
            ax.plot(xx, Fit, color=clr[p],  linewidth=2.0, linestyle='--')
            ax.fill_between(xx,Fit - Fit*0.2 ,Fit + Fit*0.2 ,alpha=0.25, color=clr[p])
        
        ax.set_ylabel('C$_2$H$_6$ VMR [ppbv]', fontsize=18)
        ax.set_xlabel('CO VMR [ppbv]', fontsize=18)
        ax.set_ylim(0, 62.0)
        ax.set_xlim(50, 400.0)
        ax.grid(True)
        ax.tick_params(axis='both', which='major', labelsize=18)
        ax.legend(prop={'size':18})
        fig.subplots_adjust(left = 0.12, bottom=0.12, top=0.96, right = 0.95)

        if saveFlg: 
           pdfsav.savefig(fig,dpi=200)
           
        else:           
            plt.show(block=False)

        if pltMapC2H6toCO:
            mp.pltMapZ( LonData=lo_int, LatData=la_int, zData=Ra, zmin=0.01, zmax=0.1, ztitle='Ratio C$_2$H$_6$ to CO [ppb]', 
                        LatID=latB, LonID=lonB, SaveFile='/data/iortega/results/fl0/C130_FRAPPE_MAP_C2H6toCO.pdf')

    if pltNH3toCO:
        x        = np.concatenate(COC130)
        xDate    = np.concatenate(DateC130_CO)
        xdoy     = mf.toYearFraction(xDate)

        y        = np.concatenate(NH3C130)
        yDate    = np.concatenate(DateC130_NH3)
        ydoy     = mf.toYearFraction(yDate)

        la       = np.concatenate(LatC130)
        lo       = np.concatenate(LonC130)
        d        = np.concatenate(DateC130)
        ddoy     = mf.toYearFraction(d)

        y_int     = interpolate.interp1d(ydoy, y,  bounds_error=False)(xdoy)
        la_int    = interpolate.interp1d(ddoy, la,  bounds_error=False)(xdoy)
        lo_int    = interpolate.interp1d(ddoy, lo,   bounds_error=False)(xdoy)

        Ra        = np.true_divide(y_int, x)

        fig, ax = plt.subplots(figsize=(10,6), sharex=True)

        #ax.plot(x, y_int, '.k', markersize=1, linewidth=0, alpha=0.15)
        #ax.scatter(x, y_int, facecolors='white', s=50, color='gray', alpha=0.15)
        slope, intercept, r_value, p_value, std_err = stats.linregress(x, y_int)
        print '\nCorrelation of all data' + ' NH3 to CO'
        print 'Slope = ', str(slope)
        print 'intercept = ', str(intercept)
        print 'r_value = ', str(r_value)
        print 'std_err =', str(std_err)

        for p, loc in enumerate(locP):
            Dist     = []

            for i, k in enumerate(x):
                c = mf.haversine(abs(loc[0]), abs(loc[1]), abs(lo_int[i]), abs(la_int[i]) )
                Dist.append(float(c))

            Dist = np.asarray(Dist)
            inds = np.where(Dist <= maxD)[0]

            x1 = np.asarray(x[inds])
            y1 = np.asarray(y_int[inds])
            
            ax.scatter(x1, y1, facecolors='white', s=70, color=clr[p], label=IDP[p])

            slope, intercept, r_value, p_value, std_err = stats.linregress(x1, y1)
            print '\nCorrelation of '+ IDP[p] + ' NH3 to CO'
            print 'Slope = ', str(slope)
            print 'intercept = ', str(intercept)
            print 'r_value = ', str(r_value)

            xx = range(0, 500, 2)
            xx = np.asarray(xx)
            Fit = slope*np.asarray(xx) + (intercept)#*math.exp(-1.6/7.4)
            ax.plot(xx, Fit, color=clr[p],  linewidth=2.0, linestyle='--')
            ax.fill_between(xx,Fit - Fit*0.2 ,Fit + Fit*0.2 ,alpha=0.25, color=clr[p])
            
        
        ax.set_ylabel('NH$_3$ VMR [ppbv]', fontsize=18)
        ax.set_xlabel('CO VMR [ppbv]', fontsize=18)
        ax.set_ylim(0, 40.0)
        ax.set_xlim(50, 400.0)
        ax.grid(True)
        ax.tick_params(axis='both', which='major', labelsize=18)
        ax.legend(prop={'size':18})
        fig.subplots_adjust(left = 0.12, bottom=0.12, top=0.96, right = 0.95)

        if saveFlg: 
           pdfsav.savefig(fig,dpi=200)
           
        else:           
            plt.show(block=False)

        if pltMapNH3toCO:
            mp.pltMapZ( LonData=lo_int, LatData=la_int, zData=Ra, zmin=0.01, zmax=0.075, ztitle='Ratio NH$_3$ to CO [ppb]', 
                        LatID=latB, LonID=lonB, SaveFile='/data/iortega/results/fl0/C130_FRAPPE_MAP_NH3toCO.pdf')
    
    if pltH2COtoCO:
        x        = np.concatenate(COC130)
        xDate    = np.concatenate(DateC130_CO)
        xdoy     = mf.toYearFraction(xDate)

        y        = np.concatenate(H2COC130)
        yDate    = np.concatenate(DateC130_H2CO)
        ydoy     = mf.toYearFraction(yDate)

        la       = np.concatenate(LatC130)
        lo       = np.concatenate(LonC130)
        d        = np.concatenate(DateC130)
        ddoy     = mf.toYearFraction(d)

        y_int     = interpolate.interp1d(ydoy, y,  bounds_error=False)(xdoy)
        la_int    = interpolate.interp1d(ddoy, la,  bounds_error=False)(xdoy)
        lo_int    = interpolate.interp1d(ddoy, lo,   bounds_error=False)(xdoy)

        Ra        = np.true_divide(y_int, x)

        fig, ax = plt.subplots(figsize=(10,6), sharex=True)

        #ax.plot(x, y_int, '.k', markersize=1, linewidth=0, alpha=0.15)
        #ax.scatter(x, y_int, facecolors='white', s=50, color='gray', alpha=0.15)
        slope, intercept, r_value, p_value, std_err = stats.linregress(x, y_int)
        print '\nCorrelation of all data' + ' H2CO to CO'
        print 'Slope = ', str(slope)
        print 'intercept = ', str(intercept)
        print 'r_value = ', str(r_value)
        print 'std_err =', str(std_err)

        for p, loc in enumerate(locP):
            Dist     = []

            for i, k in enumerate(x):
                c = mf.haversine(abs(loc[0]), abs(loc[1]), abs(lo_int[i]), abs(la_int[i]) )
                Dist.append(float(c))

            Dist = np.asarray(Dist)
            inds = np.where(Dist <= maxD)[0]

            x1 = np.asarray(x[inds])
            y1 = np.asarray(y_int[inds])
            
            ax.scatter(x1, y1, facecolors='white', s=70, color=clr[p], label=IDP[p])

            slope, intercept, r_value, p_value, std_err = stats.linregress(x1, y1)
            print '\nCorrelation of '+ IDP[p] + ' H2CO to CO'
            print 'Slope = ', str(slope)
            print 'intercept = ', str(intercept)
            print 'r_value = ', str(r_value)
            print 'std_err =', str(std_err)

            xx = range(0, 500, 2)
            xx = np.asarray(xx)
            Fit = slope*np.asarray(xx) + (intercept)#*math.exp(-1.6/7.4)
            ax.plot(xx, Fit, color=clr[p],  linewidth=2.0, linestyle='--')
            ax.fill_between(xx,Fit - Fit*0.2 ,Fit + Fit*0.2 ,alpha=0.25, color=clr[p])
            
        
        ax.set_ylabel('H$_2$CO VMR [ppbv]', fontsize=16)
        ax.set_xlabel('CO VMR [ppbv]', fontsize=16)
        ax.set_ylim(0, 7.0)
        ax.set_xlim(50, 400.0)
        ax.grid(True)
        ax.tick_params(axis='both', which='major', labelsize=14)
        ax.legend(prop={'size':12})
        fig.subplots_adjust(left = 0.12, bottom=0.12, top=0.96, right = 0.95)

        if saveFlg: 
           pdfsav.savefig(fig,dpi=200)
           
        else:           
            plt.show(block=False)

        if pltMapH2COtoCO:
            mp.pltMapZ( LonData=lo_int, LatData=la_int, zData=Ra, zmin=0.005, zmax=0.04, ztitle='Ratio H$_2$CO to CO [ppb]', 
                        LatID=latB, LonID=lonB, SaveFile='/data/iortega/results/fl0/C130_FRAPPE_MAP_H2COtoCO.pdf')
   
    #------------------------------------------------------------------
    #                     PLOT MAP OF C2H6  
    #------------------------------------------------------------------
    if pltMapC2H6:

        C2H6C130    = np.concatenate(C2H6C130)
        LonC130_int = np.concatenate(LonC130_C2H6_int)
        LatC130_int = np.concatenate(LatC130_C2H6_int)

        mp.pltMapZ( LonData=LonC130_int, LatData=LatC130_int, zData=C2H6C130, zmin=0.5, zmax=6.5, ztitle='C$_2$H$_6$ [ppb]', 
                    LatID=latB, LonID=lonB, SaveFile='/data/iortega/results/fl0/C130_FRAPPE_MAP_C2H6.pdf')
    
    #------------------------------------------------------------------
    #                     PLOT MAP OF CO  
    #------------------------------------------------------------------
    if pltMapCO:

        COC130      = np.concatenate(COC130)
        LonC130_int = np.concatenate(LonC130_CO_int)
        LatC130_int = np.concatenate(LatC130_CO_int)

        mp.pltMapZ( LonData=LonC130_int, LatData=LatC130_int, zData=COC130, zmin=50, zmax=150, ztitle='CO [ppb]', 
                    LatID=latB, LonID=lonB, SaveFile='/data/iortega/results/fl0/C130_FRAPPE_MAP_CO.pdf')

    #------------------------------------------------------------------
    #                     PLOT MAP OF H2CO  
    #------------------------------------------------------------------
    if pltMapH2CO:

        H2COC130      = np.concatenate(H2COC130)
        LonC130_int = np.concatenate(LonC130_H2CO_int)
        LatC130_int = np.concatenate(LatC130_H2CO_int)

        mp.pltMapZ( LonData=LonC130_int, LatData=LatC130_int, zData=H2COC130, zmin=0.1, zmax=4.0, ztitle='H$_2$CO [ppb]', 
                    LatID=latB, LonID=lonB, SaveFile='/data/iortega/results/fl0/C130_FRAPPE_MAP_H2CO.pdf')

    #------------------------------------------------------------------
    #                     PLOT MAP OF NH3  
    #------------------------------------------------------------------
    if pltMapNH3:

        NH3C130      = np.concatenate(NH3C130)
        LonC130_int = np.concatenate(LonC130_NH3_int)
        LatC130_int = np.concatenate(LatC130_NH3_int)

        mp.pltMapZ( LonData=LonC130_int, LatData=LatC130_int, zData=NH3C130, zmin=0.05, zmax=15.0, ztitle='NH$_3$ [ppb]', 
                    LatID=latB, LonID=lonB, SaveFile='/data/iortega/results/fl0/C130_FRAPPE_MAP_NH3.pdf')
    
    #------------------------------------------------------------------
    #                     PLOT MAP WITH RF TRACKS  
    #------------------------------------------------------------------
    if pltMapRF:
        mp.pltMapRF(LonData=LonC130, LatData=LatC130, zData=RF, ztitle='Research Flight', 
                    LatID=latB, LonID=lonB, locP=locP, SaveFile='/data/iortega/results/fl0/C130_FRAPPE_MAP_RF.pdf')

      #--------------------------------
      # Pause so user can look at plots
      #--------------------------------
    if saveFlg: 
        pdfsav.close()
    else:
        user_input = raw_input('Press any key to exit >>> ')
        sys.exit()           # Exit program        


if __name__ == "__main__":
    main(sys.argv[1:])