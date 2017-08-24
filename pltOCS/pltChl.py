#! /usr/bin/python2.7
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
# Name:
#        pltChl.py
#
# Purpose:
#       Script to plot maps/time series of Chlorophyl A
#
# Sources :
#      - https://oceancolor.gsfc.nasa.gov/
#
#
# Version History:
#       Created, July, 2016  Ivan Ortega (iortega@ucar.edu)
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------

import numpy as np
import os
import subprocess as sp
import sys
import datetime as dt

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from pylab import *
from netCDF4 import Dataset
import datetime, types
import string
from matplotlib.backends.backend_pdf import PdfPages 

from math import acos, asin, atan2, cos, hypot, sin, sqrt
from math import degrees, pi as PI, radians, sin, tan

import pyproj    
import shapely
import shapely.ops as ops
from shapely.geometry.polygon import Polygon
from functools import partial
from matplotlib.patches import Polygon as PatchPolygon

from netCDF4 import Dataset
import myfunctions as mf
import fnmatch
import glob

def findCls2(dataArray, val):
    return min(dataArray, key=lambda p:dist_sq(val, p))

def dist_sq(a, b): # distance squared (don't need the square root)
  return (a[0] - b[0])**2 + (a[1] - b[1])**2

def findCls(dataArray, val):
    ''' Returns the indice and closest value in dataArray to val'''
    return np.argmin(abs(val-dataArray))

def destination(lat, lon, distance, bearing, radius=6371.008771415):
    '''Return the destination LatLon from this point after having
       travelled the given distance on the given initial bearing.
       @param {number} distance - Distance travelled.
       @param {degrees} bearing - Initial bearing in degrees from North.
       @param {number} [radius=R_M] - Mean radius of earth (default.

    Notes: see https://github.com/mrJean1/PyGeodesy/blob/master/geodesy/sphericalTrigonometry.py
    '''
    d        = float(distance)/float(radius) 
    t        = radians(bearing)
    latRad   = radians(lat)
    lonRad   = radians(lon)

    ca1, sa1 = cos(latRad), sin(latRad)
    cd, sd   = cos(d), sin(d) 

    a2 = asin(sa1 * cd + ca1 * sd * cos(t))
    b2 = atan2(sin(t) * sd * ca1, cd - sa1 * sin(a2)) + lonRad

    lat2  = degrees(a2)
    lon2  = degrees(b2)

    return lat2, lon2

def main():

    #----------------------------------------------------------------------------------------
    #----------------------------------INPUTS-----------------------------------------------
    #----------------------------------------------------------------------------------------
    Sensor    = 'MODIS'  # VIIRS
    DataDir   = '/data1/iortega/Chlorophyll/'+Sensor.upper()+'/'
   
    saveFlg   = True
    
    moi       = 6       #MONTHS OF INTEREST FOR MAP

    Loclat    = 76.52   #THULE
    Loclon    = -68.7   
   
    #-------------------------------------------
    #DEFINING THE LATITUDE/LONGITUDE TO SAVE MATRICES
    #--------------------------------------------
    north = 85.
    south = 65.
    west  = -85.
    east  = -45.

    #-------------------------------------------
    #DEFINING THE LATITUDE/LONGITUDE TO AREA TO CALCULATE TIME SERIES
    #--------------------------------------------
    #north2 = 78.
    #south2 = 75.
    #west2 = -78
    #east2 = -70

    north2  = 78.
    south2  = 74.
    west2   = (282 - 360.)
    east2   = (300 - 360.)

    if saveFlg:
        #outFname  = DataDir +'plt_chlorophyll_'+"{0:02d}".format(yoi)+"{0:02d}".format(moi)+'.pdf'
        outFname  ='/data1/iortega/Chlorophyll/plt_chlorophyll_Month'+"{0:02d}".format(moi)+'_'+Sensor.upper()+'.pdf'
        
        pdfsav = PdfPages(outFname)
    else: pdfsav = False

    #-------------------------------------------------
    #Get variables
    #-------------------------------------------------

    Files = glob.glob(DataDir+'/'+str(moi)+'/' + '*.nc')
    Files.sort()

    #-------------------------------------------------
    #Lat and Lon are the same for all Files: Defining LAT and LON. Find Shapes
    #-------------------------------------------------
    nc_fid1 = Dataset(Files[0] , 'r')             
    nc_attrs, nc_dims, nc_vars = mf.ncdump(nc_fid1, verb=False)
    
    lon          = np.squeeze(nc_fid1.variables['lon'][:])
    lat          = np.squeeze(nc_fid1.variables['lat'][:])


    #-------------------------------------------
    #DEFINING THE AREA TO CALCULATE TIME SERIES
    #--------------------------------------------
    indlon = np.where( (lon >= west) & (lon <= east) )[0]   
    indlat = np.where( (lat >= south) & (lat <= north) )[0] 

    lon    = lon[indlon]
    lat    = lat[indlat]

    Nlon  = np.asarray(lon).shape[0]
    Nlat  = np.asarray(lat).shape[0]    

    #-------------------------------------------------
    #Loop through years
    #-------------------------------------------------
    chlmatrix  = np.zeros((len(Files), Nlat, Nlon))  #MATRIX FOR CHLOROPYL
    chlAnomaly = np.zeros((len(Files), Nlat, Nlon))  #MATRIX FOR CHLOROPYL ANOMALY
    years      = []                                  #YEARS
    chlTs      = []                                  #CHL TIME Series
    chlAnoTs   = []                                  #CHL ANOMALY TIME Series
   
    for i, f in enumerate(Files):

    	fsplit = f.split('/')

    	year = int(fsplit[-1][1:5])
    	years.append(year)

    	nc_fid1 = Dataset(f , 'r')             
    	nc_attrs, nc_dims, nc_vars = mf.ncdump(nc_fid1, verb=False)

        #-------------------------------------------------
        #MODIS AND VIIRS DEFINE DIFFERENT CHL
        #-------------------------------------------------
    	if Sensor.upper() == 'VIIRS': chl = np.asarray(nc_fid1.variables['chl_oc3'])
    	
        if Sensor.upper() == 'MODIS': 
    	    if moi == 6:
                if year < 2014: 
                    chl = np.asarray(nc_fid1.variables['chl_ocx'])
                else:
                    chl = np.asarray(nc_fid1.variables['chl_oc3'])
            
            else: 
                if year < 2015: 
                    chl = np.asarray(nc_fid1.variables['chl_ocx'])
                else:
                    chl = np.asarray(nc_fid1.variables['chl_oc3'])

    		

        #-------------------------------------------------
        #TRICK TO SAVE ONLY AREA OF INTEREST (otherwise is to slow)
        #-------------------------------------------------
        chl          = chl[indlat, :]
        chl          = chl[:, indlon]

        chlmatrix[i,:,:]    = chl

    years     = np.asarray(years)
    chlmatrix = np.asarray(chlmatrix)

    #-------------------------------------------------
    #FILTERING NEGATIVE VALUES
    #-------------------------------------------------
    chlmatrix[chlmatrix < 0.0]   = np.nan
    #chlmatrix[chlmatrix > 100.0] = np.nan

    #-------------------------------------------------
    #DEFINING THE AREA OF INTEREST FOR TIME SERIES
    #-------------------------------------------------
    indlon2  = np.logical_and(lon >= west2, lon <= east2)
    indlat2  = np.logical_and(lat >= south2, lat <= north2)

    lat2 = lat[indlat2]
    lon2 = lon[indlon2]

    #-------------------------------------------------
    #POLYGON OF AREA OF INTEREST
    #-------------------------------------------------
    indlon2  = np.logical_and(lon >= west2, lon <= east2)
    indlat2  = np.logical_and(lat >= south2, lat <= north2)

    Inipoly = [ (lon2[0], lat2[0]), (lon2[0], lat2[-1]), (lon2[-1], lat2[-1]), (lon2[-1], lat2[0])]

    poly = Polygon(Inipoly)

    #-------------------------------------------------
    #CALCULATING THE GEOMETRIC AREA 
    #-------------------------------------------------
    indlon2  = np.logical_and(lon >= west2, lon <= east2)
    indlat2  = np.logical_and(lat >= south2, lat <= north2)

    geom_area = ops.transform(
        partial(
            pyproj.transform,
            pyproj.Proj(init='EPSG:4326'),
            pyproj.Proj(
                proj='aea',
                lat1=poly.bounds[1],
                lat2=poly.bounds[3])),
        poly)

    AreaPoly = geom_area.area/1000. /1000.
    print 'Geometric Area in Polygon: {0:.2f}km2'.format(AreaPoly)

    #-------------------------------------------
    #ANOMALIES AND TIME SERIES
    #-------------------------------------------
    chlMean = np.nanmean(chlmatrix, axis=0)

    for i, f in enumerate(Files):

        # chlmatrix2 = chlmatrix[i, indlat2, :]
        # chlmatrix2 = chlmatrix2[:, indlon2]
        # chlTs.append(np.nanmean(chlmatrix2) )

        chlAnomaly[i,:,:] = chlmatrix[i,:,:] - chlMean[:,:]
        chlAnomaly2 = chlAnomaly[i, indlat2, :]
        chlAnomaly2 = chlAnomaly2[:, indlon2]

        chlAnoTs.append(np.nanmean(chlAnomaly2 ) )

    chlAnomaly = np.asarray(chlAnomaly)
    chlAnoTs   = np.asarray(chlAnoTs)
    chlTs      = np.asarray(chlTs)


    #-------------------------------------------
    #PLOTS
    #-------------------------------------------
    latMin   = 70.   # 54.
    latMax   = 80.    # 88.
    lonMin   = -75.   # -65.
    lonMax   = -40.    # 20

    levels = np.arange(0, 10, 0.2)
    extender = "max"

    lons, lats = np.meshgrid(lon, lat)

    if saveFlg:

        for i, yyyy in enumerate(years):

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
    	    
    	    #lons2, lats2 = np.meshgrid(lon2, lat2)
    	    #x2, y2 = map(lons2, lats2)
    	    #ax1.plot(x2, y2, '.k', alpha = 0.05)

    	    if saveFlg:
    	        pdfsav.savefig(fig,dpi=200)
    	    else:
    	        plt.show(block=False)
                #user_input = raw_input('Press any key to exit >>> ')
                #exit()
      
        #-------------------------------------------
        #
        #-------------------------------------------
        levels = np.arange(-5, 5, 0.2)

        for i, yyyy in enumerate(years):

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

    	    CS1 = map.contourf(x,y,chlAnomaly[i, :, :], levels)
    	    CS1.axis='tight'
    	    bar = plt.colorbar(CS1,orientation='vertical',extend=extender)#, shrink=0.5)
    	    bar.set_label('Chlorophyll Concentration [mg m$^{-3}$]')

    	    #ax1.plot(x2, y2, '.k', alpha = 0.05)

    	    plt.title('Chlorophyll Concentration Anomaly' + '({} to {}) \nYear = {},  Month = {}'.format(years[0], years[-1], yyyy, moi))

    	    if saveFlg:
    	        pdfsav.savefig(fig,dpi=200)
    	    else:
    	        plt.show(block=False)

    #-------------------------------------------
    #Area
    #-------------------------------------------
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

    plt.title('Area Used to calculate time series')

    lons2, lats2 = np.meshgrid(lon2, lat2)
    x2, y2 = map(lons2, lats2)
    ax1.plot(x2, y2, '.k', alpha = 0.05)

    if saveFlg:
        pdfsav.savefig(fig,dpi=200)
    else:
        plt.show(block=False)

    #-------------------------------------------
    #PLOT OF TIME SERIES
    #-------------------------------------------
    y_pos = np.arange(len(years))
   
    Colors = []
    for val in chlAnoTs:
        if val <= 0: Colors.append('blue')# chocolate
        elif val > 0.0: Colors.append('red')\

    #fig1 , (ax1, ax2) = plt.subplots(2, figsize=(10,8), sharex=True)
    fig1 , ax2 = plt.subplots(1, figsize=(10,6), sharex=True)

    # ax1.bar(y_pos,chlTs, color='gray', align='center', edgecolor="black")
    # ax1.grid(True)
    # ax1.set_ylabel('Concentration [mg m$^{-3}$]',multialignment='center', fontsize=14)
    # ax1.set_title('Chlorophyll Concentration [Month={}]'.format(moi))
    # ax1.axhline(y=np.nanmean(chlTs),color='r')
    
    ax2.bar(y_pos,chlAnoTs,color=Colors, align='center', edgecolor="black")
    ax2.grid(True)
    ax2.set_ylabel('Anomaly [mg m$^{-3}$]',multialignment='center', fontsize=14)
    ax2.set_title('Chlorophyll Concentration Anomaly [Month={}]'.format(moi))
    ax2.set_xlabel('Year', fontsize=14)
    ax2.set_xlim(y_pos[0]-1,y_pos[-1]+1)
    ax2.set_ylim(-1.0, 1.5)
    ax2.tick_params(which='both',labelsize=14)

    plt.xticks(y_pos, years, rotation=45)
    plt.subplots_adjust(bottom=0.15, left=0.1, right=0.95, top=0.95)


    if saveFlg:
        pdfsav.savefig(fig1,dpi=200)
    else:
	    plt.show(block=False)


    if saveFlg:
        pdfsav.close()
    else:
        user_input = raw_input('Press any key to exit >>> ')
        exit()

if __name__ == "__main__":
    main()