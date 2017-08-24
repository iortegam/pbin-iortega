#! /usr/bin/python2.7
#----------------------------------------------------------------------------------------
# Name:
#         pltmap.py
#
# Purpose:
#         Plot time series of multiple species using HDF GEOMS
#
# Notes:  
#         Initially created for OCS (using trop height) and see trend in straospheric OCS
#   
#
# Version History:
#       Created, Feb, 2017  Ivan Ortega (iortega@ucar.edu)

#----------------------------------------------------------------------------------------

                                    #-------------------------#
                                    # Import Standard modules #
                                    #-------------------------#

import numpy as np
import sys
import glob

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import FormatStrFormatter, MultipleLocator,AutoMinorLocator,ScalarFormatter
from matplotlib.backends.backend_pdf import PdfPages #to save multiple pages in 1 pdf...
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as mplcm
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
import matplotlib.gridspec as gridspec
from itertools import izip
from numpy import *
import srtm
from mpl_toolkits.basemap import Basemap


                                    #----------------------------#
                                    #                            #
                                    #        --- Main---         #
                                    #                            #
                                    #----------------------------#

def main():
    #---------------------------------------------------
    # Order data based on +Lat to -Lat
    #---------------------------------------------------


    pltID        = ['Thule', 'Mauna Loa', 'Boulder']
    Lat          = [76.56, 19.54, 40.04]
    Lon          = [291.26, 204.43, 254.76]
    saveFlg      = False

    
    #---------------------------------------------------
    # Map with Sites
    #---------------------------------------------------

    fig,ax = plt.subplots(figsize=(15,7))
   

    #eq_map = Basemap(projection='robin', resolution = 'l', area_thresh = 1000.0, lat_0=40, lon_0=-130)

    eq_map = Basemap(llcrnrlon=-150,llcrnrlat=1., urcrnrlon=-2.566, urcrnrlat=46.352,\
            rsphere=(6378137.00, 6356752.3142),\
            resolution='l',area_thresh=1000.,projection='lcc',\
            lat_1=50.,lon_0=-107.,  ax=ax)


    eq_map.drawcoastlines()
    eq_map.drawcountries()
    #eq_map.fillcontinents(color = 'gray', alpha=0.5)
    eq_map.drawmapboundary()
    eq_map.drawmapboundary(fill_color='#99ffff')
    eq_map.fillcontinents(color='#cc9966',lake_color='#99ffff')
    #eq_map.drawparallels(np.arange(10,70,20),labels=[1,1,1,1])
    #eq_map.drawmeridians(np.arange(-100,0,20),labels=[0,0,0,1])
    eq_map.drawmeridians(np.arange(0, 360, 20))
    eq_map.drawparallels(np.arange(0, 90, 30))

    for lo, la, idhdf in zip(Lon, Lat, pltID):
        x,y = eq_map(lo, la)
        eq_map.plot(x, y, 'ro', markersize=20.0)
        #plt.text(x+10000,y-800000, idhdf, fontsize=12, color='red')

    #plt.title("HR-FTIR NCAR Sites", fontsize=18)

    if saveFlg: pdfsav.savefig(fig,dpi=200)
    else: 
        plt.show(block=False)
        user_input = raw_input('Press any key to exit >>> ')
        sys.exit()

if __name__ == "__main__":
    main()