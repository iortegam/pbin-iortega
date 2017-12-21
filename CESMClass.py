#----------------------------------------------------------------------------------------
# Name:
#        CESMClass.py
#
# Purpose:
#	     Initial Class to Test to Read NcDF files (Siyuan's)

# Notes:
#
#
# Version History:
#       Created, Oct, 2017  Ivan Ortega (iortega@ucar.edu)
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

from scipy.integrate import simps
from scipy import interpolate

import matplotlib
#
from numpy import fromstring, vectorize, ndarray, array, genfromtxt

import matplotlib.dates as md
from matplotlib.dates import DateFormatter, MonthLocator, YearLocator, DayLocator, WeekdayLocator, MONDAY

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import FormatStrFormatter, MultipleLocator,AutoMinorLocator,ScalarFormatter
from matplotlib.backends.backend_pdf import PdfPages #to save multiple pages in 1 pdf...
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as mplcm
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
import matplotlib.gridspec as gridspec

from netCDF4 import Dataset, Group

from mpl_toolkits.basemap import Basemap
import sfitClasses as sc
import glob
from datetime import datetime, timedelta
from geographiclib.geodesic import Geodesic
from mpl_toolkits.basemap import Basemap
import srtm

from matplotlib.backends.backend_pdf import PdfPages 
from scipy.io               import netcdf

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

def ncdump(nc_fid, verb=True):
    '''
    ncdump outputs dimensions, variables and their attribute information.
    The information is similar to that of NCAR's ncdump utility.
    ncdump requires a valid instance of Dataset.

    Parameters
    ----------
    nc_fid : netCDF4.Dataset
        A netCDF4 dateset object
    verb : Boolean
        whether or not nc_attrs, nc_dims, and nc_vars are printed

    Returns
    -------
    nc_attrs : list
        A Python list of the NetCDF file global attributes
    nc_dims : list
        A Python list of the NetCDF file dimensions
    nc_vars : list
        A Python list of the NetCDF file variables
    '''
    def saveFlg(key):
        """
        Prints the NetCDF file attributes for a given key

        Parameters
        ----------
        key : unicode
            a valid netCDF4.Dataset.variables key
        """
        try:
            print "\t\ttype:", repr(nc_fid.variables[key].dtype)
            for ncattr in nc_fid.variables[key].ncattrs():
                print '\t\t%s:' % ncattr,\
                      repr(nc_fid.variables[key].getncattr(ncattr))
        except KeyError:
            print "\t\tWARNING: %s does not contain variable attributes" % key

    # NetCDF global attributes
    nc_attrs = nc_fid.ncattrs()
    #if verb:
    #    print "NetCDF Global Attributes:"
    #    for nc_attr in nc_attrs:
    #        print '\t%s:' % nc_attr, repr(nc_fid.getncattr(nc_attr))
    nc_dims = [dim for dim in nc_fid.dimensions]  # list of nc dimensions
    # Dimension shape information.
    if verb:
        print "NetCDF dimension information:"
        for dim in nc_dims:
            print "\tName:", dim 
            print "\t\tsize:", len(nc_fid.dimensions[dim])
            saveFlg(dim)
    # Variable information.
    nc_vars = [var for var in nc_fid.variables]  # list of nc variables
    if verb:
        print "NetCDF variable information:"
        for var in nc_vars:
            if var not in nc_dims:
                print '\tName:', var
                print "\t\tdimensions:", nc_fid.variables[var].dimensions
                print "\t\tsize:", nc_fid.variables[var].size
                saveFlg(var)
    return nc_attrs, nc_dims, nc_vars


class ncClass():
    '''
    This class deals with reading nc CESM outputs.
    '''

    def __init__(self,dataDir,  outFname='', saveFlg= False):
        
        # check if '/' is included at end of path
        if not( dataDir.endswith('/') ):
            dataDir = dataDir + '/'           
            
        # Check if directory exits
        ckDir(dataDir, exitFlg=True)
        
        self.dataDir  = dataDir

        files = glob.glob( dataDir + '*.nc')

        files.sort()
                    
        self.files = files

        if saveFlg:
            self.pdfsav = PdfPages(outFname)
            self.saveFlg   = True
        else: self.saveFlg = False

    def ReadNCCESM():

        self.lat     = []
        self.lon     = []
        self.alt     = []

        for indMain, f in enumerate(self.files):

            nc_fid = Dataset(f , 'r')  # Dataset is the class behavior to open the file
            nc_attrs, nc_dims, nc_vars = mf.ncdump(nc_fid, verb=True)

            vgrp =  str(nc_fid.groups)

            lat       = np.asarray(nc_fid.variables['latitude'])
            lon       = np.asarray(nc_fid.variables['longitude'])
            alt       = np.asarray(nc_fid.variables['longitude'])

            self.lat.append(lat)
            self.lon.append(lon)
            self.alt.append(alt)

        self.lat   = np.asarray(self.lat)
        self.lon   = np.asarray(self.lon)
        self.alt   = np.asarray(self.alt)

    def pltCESM():

        #-------------------------------------------
        #PLOTS
        #-------------------------------------------
        latMin   = 70.   
        latMax   = 80.    
        lonMin   = -75.  
        lonMax   = -40.   

        levels = np.arange(0, 10, 0.2)
        extender = "max"

        lons, lats = np.meshgrid(self.lon, self.lat)

        for i, f in enumerate(self.files):

            data = self.data

            fig, ax1  = plt.subplots()

            map = Basemap(llcrnrlat=latMin,urcrnrlat=latMax,
                    llcrnrlon=lonMin,urcrnrlon=lonMax,
                    rsphere=(6378137.00,6356752.3142),
                    resolution='l',area_thresh=1000.,projection='lcc',
                    lat_1=latMin,lon_0=-60)

            x, y = map(lons, lats)

            map.drawcoastlines(color='black')
            map.drawcountries(color='lightgrey')
            map.fillcontinents(color='gray')

            map.drawparallels(np.arange(-80., 81., 5.), labels=[1,0,0,0], alpha=0.5)
            map.drawmeridians(np.arange(-180., 181., 10.), labels=[0,0,0,1],alpha=0.5)

            #CS1 = map.contourf(x,y,chl.T,levels,cmap=cm.get_cmap('RdYlBu_r',len(levels)),extend=extender)#,alpha=0.5)
            CS1 = map.contourf(x,y,chlmatrix[i, :, :], levels)#,alpha=0.5)
            CS1.axis='tight'
            bar = plt.colorbar(CS1,orientation='vertical',extend=extender)#, shrink=0.5)
            bar.set_label('Chlorophyll Concentration [mg m$^{-3}$]')

            #x, y = map(lons[indlon2], lats[indlat2])
            #ax1.plot(x, y, '.k')

            plt.title('Chlorophyll Concentration\n' + '(Year = {},  Month = {})'.format(yyyy, moi))

            if saveFlg:
                pdfsav.savefig(fig,dpi=200)
            else:
                plt.show(block=False)
                #user_input = raw_input('Press any key to exit >>> ')
                #exit()

            