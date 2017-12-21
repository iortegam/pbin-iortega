#----------------------------------------------------------------------------------------
# Name:
#        MainCESM.py
#
# Purpose:
#	     Initial Test to Read NcDF files (Siyuan's)

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
import CESMClass as mf

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

									#----------------------------#
                                    #                            #
                                    #        --- Main---         #
                                    #                            #
                                    #----------------------------#

def main():

    #-------------------------------------
    #INITIALIZATIONS FOR HIPPO
    #-------------------------------------
    dataDir   = '/data1/Campaign/'
    saveFlg   = False
    pltFile   = dataDir+'Test_NC.pdf'


									#----------------------------#
                                    #                            #
                                    #        --- Start---         #
                                    #                            #
                                    #----------------------------#

  
    DataSet   = mf.ncClass(dataDir,  outFname= pltFile, saveFlg= saveFlg)
    Data.ReadNCCESM()


if __name__ == "__main__":
    main()
   