#! /usr/bin/python2.7
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
# Name:
#        pltSSTIce.py
#
# Purpose:
#       Script to plot maps/time series of SST and Ice Conc
#
# External Subprocess Calls:
#      - F2py  -- > Fortran to python 
#      -f2py takes a Fortran subroutine and some additional instructions, compiles the Fortran code and builds a module which can then be imported
#      into Python and used there like a normal function
#      
#       Details: http://websrv.cs.umt.edu/isis/index.php/F2py_example
#
#Steps:
#       f2py -m fib1 -h fib1.pyf fib1.f --> NOTE: For some reason this won't compile right unless the itag arrays is removed 
#       f2py -c fib1.pyf fib1.f
#       run python
#
# Version History:
#       Created, Jan, 2016  Ivan Ortega (iortega@ucar.edu)
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------

import numpy as np
import os
import subprocess as sp
import py_fortran_tools
import sys
sys.path.append('..')
import datetime as dt

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from pylab import *
from netCDF4 import Dataset
import datetime, types
import string
from matplotlib.backends.backend_pdf import PdfPages 

import rsstice
import ldse

from math import acos, asin, atan2, cos, hypot, sin, sqrt
from math import degrees, pi as PI, radians, sin, tan

import pyproj    
import shapely
import shapely.ops as ops
from shapely.geometry.polygon import Polygon
from functools import partial
from matplotlib.patches import Polygon as PatchPolygon

print rsstice.fib.__doc__
print ldse.ldse.__doc__

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

#---------------------------------------------------------------------------------------------
#PLOT ICE ANOMALIES
#---------------------------------------------------------------------------------------------

def pltMap(lon, lat, data, siteloc=False,type=False, saveFlg=False, startYear=False, endYear=False, yoi=False, moi=False, polyFlg=False, globalmapFlg=False, saa=False, xdist=False):
    fig, ax1  = plt.subplots()


    if type == "iceanomaly":
        levels   = np.arange(-45, 50, 5)
        extender = "both"
        barlabel = 'Sea Ice Concentration Anomalies [%]'

    elif type == "sstanomaly":
        levels = np.arange(-5, 5, 0.5)
        extender = "both"
        barlabel = 'Sea Surface Temperature Anomalies [C]'
    
    elif type == "sst":
        levels = np.arange(-10, 10, 1.0)
        extender = "max"
        barlabel = 'Sea Surface Temperature [C]'
    
    elif type == "ice":
        levels = np.arange(0, 100, 5)
        extender = "max"
        barlabel = 'Sea Ice Concentration [%]'
    
    elif type == 'area':
        barlabel = 'Area used to calculate time series'
    
    else:
        print 'Error in type: options --> iceanomaly, sstanomaly, sst, ice'
        exit()
        
    labels=[9,10]

    if globalmapFlg is True:
        #map = Basemap(width=10.e6,height=11.e6,projection='gnom',lat_0=88.,lon_0=-30.)
        map = Basemap(projection='npstere',boundinglat=65,lon_0=270,  resolution='l', ax=ax1)
    else:
        latMin   = 70.   # 54.
        latMax   = 80    # 88.
        lonMin   = -75   # -65.
        lonMax   = -40    # 20
        map = Basemap(llcrnrlat=latMin,urcrnrlat=latMax,
            llcrnrlon=lonMin,urcrnrlon=lonMax,
            rsphere=(6378137.00,6356752.3142),
            resolution='l',area_thresh=1000.,projection='lcc',
            lat_1=latMin,lon_0=-60)

    x, y = map(lon,lat)

    map.drawcoastlines(color='black')
    map.drawcountries(color='lightgrey')
    map.fillcontinents(color='gray')

    map.drawparallels(np.arange(-80., 81., 5.), labels=[1,0,0,0], alpha=0.5)
    map.drawmeridians(np.arange(-180., 181., 10.), labels=[0,0,0,1],alpha=0.5)

    if  polyFlg:
        Lonpoly   = []
        Latpoly   = []

        for a in saa:
            xlat_i, xlon_i = destination(siteloc[1], siteloc[0], xdist, a, radius=6371.008771415)
            Lonpoly.append(xlon_i)
            Latpoly.append(xlat_i)

            map.drawgreatcircle(siteloc[0],siteloc[1],xlon_i,xlat_i,linewidth=2,color='m')
            x2, y2 = map(xlon_i, xlat_i)
            ax1.text(x2, y2, str(a), fontsize=14, color='m')

        Inipoly   = []
        Inipoly.append((siteloc[0], siteloc[1]))

        for (lo,la) in zip(Lonpoly,Latpoly):
            Inipoly.append( (lo,la))

        Inimappoly = [map(xx,yy) for xx, yy in Inipoly]
        mappoly    = PatchPolygon(Inimappoly, color='m', facecolor='none', fill=False, linewidth=2)
        fig.gca().add_patch(mappoly)

    if data.ndim == 2: toplot = np.reshape(data, (data.shape[0],data.shape[1]))
    elif data.ndim == 3: toplot = np.reshape(data, (data.shape[1],data.shape[2]))
    elif data.ndim == 1: toplot = data
    
    try:

        CS1 = map.contourf(x,y,toplot,levels,cmap=cm.get_cmap('RdYlBu_r',len(levels)),extend=extender)#,alpha=0.5)
        CS1.axis='tight'
        bar = plt.colorbar(CS1,orientation='vertical',extend=extender)#, shrink=0.5)
        bar.set_label(barlabel)
    except Exception as errmsg:
        print '\nError: ', errmsg

    if type == 'area':
        ax1.plot(x, y, '.k')
        plt.title(barlabel)
    

    #------
    #for (a,s) in zip(saa2,sza2):
    #   Loclat2, Loclon2 = destination(Loclat, Loclon, xdist, a, radius=6371.008771415)
    #   map.drawgreatcircle(Loclon,Loclat,Loclon2,Loclat2,linewidth=2,color='m')
    #   x2, y2 = map(Loclon2, Loclat2)
    #   ax1.text(x2, y2, str(a), fontsize=14, color='m')
        
    if siteloc:
        x, y = map(siteloc[0], siteloc[1])
        ax1.plot(x, y, marker='D', color='m', markersize=7)
        if (type == "iceanomaly") or (type == "sstanomaly"): plt.title(barlabel + ' ({} to {}) \nYear = {},  Month = {}'.format(startYear, endYear, yoi, moi))
        else: plt.title(barlabel + ' \nYear = {},  Month = {}'.format(yoi, moi))

    if saveFlg:
        print "Saving plot to file: %s"%type
        saveFlg.savefig(fig,dpi=200)
        #pdfsav.close()
    else:
        plt.show(block=False)
        #user_input = raw_input('Press any key to exit >>> ')
        #exit()


def main():

    #----------------------------------------------------------------------------------------
    #----------------------------------INPUTS-----------------------------------------------
    #----------------------------------------------------------------------------------------
    DataDir   = '/data1/iortega/nsidc/datafiles/'
    ftag      = '/data1/iortega/nsidc/datafiles/land.sea.mask.v2.asc'
    saveFlg   = True

    startYear = 1999    #START YEAR
    endYear   = 2017    #END YEAR
    
    yoi       = 2016    #YEAR OF INTEREST FOR MAP
    moi       = 6       #MONTHS OF INTEREST FOR MAP

    Loclat    = 76.52   #THULE
    Loclon    = -68.7   

    polyH     = 16.   #HEIGHT OF POLY IN KM
    avgSZA    = 80.0   #TYPICAL SZA OF MEASUREMENTS LOOKING SOUTH, SOUTH-WEST
    
    saa       = [105., 161., 220]  ##    [11.48, 55.35, 64.45, 108.02, 161.02, 215.08]  ##saa      = [11.48, 12.65, 55.35, 64.45, 108.02, 109.41, 159.58, 161.02, 213.67, 215.08]
    sza       = [61.7, 53.7, 55.3] ##[79.785, 73.88, 71.92, 61.705, 53.734, 55.251]  ###sza      = [79.785, 79.724, 73.88, 71.92, 61.705, 61.403, 53.833, 53.734, 55.083, 55.251]
    
    #-------------------------------------------
    #DEFINING THE LATITUDE/LONGITUDE TO AREA TO CALCULATE TIME SERIES
    #--------------------------------------------
    #north = 76.
    #south = 74.
    #west = 290.
    #east = 300.

    north  = 78.
    south  = 74.
    west   = 282
    east   = 300

    #----------------------------------------------------------------------------------------
    #----------------------------------MAIN-----------------------------------------------
    #----------------------------------------------------------------------------------------
    if saveFlg:
        outFname  = '/data1/iortega/nsidc/plt_IceSsstLOS_Anomaly_'+"{0:02d}".format(yoi)+"{0:02d}".format(moi)+'.pdf'
        pdfsav = PdfPages(outFname)
    else: pdfsav = False

    years     = [year for year in range(startYear,endYear+1)]
    years     = np.asarray(years)
    indyear   = np.where(years == yoi)[0]

    (itag)     = ldse.ldse(ftag)
    itag       = np.asarray(itag)

    #------------------------------------------
    #MATRICES FOR MAP
    #------------------------------------------
    icematrix  = np.zeros((len(years),360,180)) 
    sstmatrix  = np.zeros((len(years),360,180))
    iceAnomaly = np.zeros((len(years),360,180))
    sstAnomaly = np.zeros((len(years),360,180))

    #------------------------------------------
    #VECTORS FOR TIME SERIES USING THE AREA DEFINED
    #------------------------------------------
    iceY       = []
    sstY       = []
    Dates      = []

    #-------------------------------------------
    #DEFINING LATITUDE/LONGITUDE
    #---------------------------------------------------------------------------------------------
    #Defining Lat and Lot arrays: The first gridbox of each array is centered on 0.5E, 89.5S.  The points
    #move eastward to 359.5E, then northward to 89.5N.
    #itagls(1,1)     =   0.5E, 89.5S
    #  itagls(2,1)     =   1.5E, 89.5S
    #  itagls(360,1)   = 359.5E, 89.5S
    #  itagls(1,2)     =   0.5E, 88.5S
    #  itagls(1,180)   =   0.5E, 89.5N
    #  itagls(360,180) = 359.5E, 89.5N
    #---------------------------------------------------------------------------------------------
    Nlon = 360
    Nlat = 180

    lat = np.zeros( (360, 180))
    lon = np.zeros( (360, 180))

    for y in range(Nlat):
        for x in range(Nlon):
            lon[x, y] = x+0.5
            lat[x, y] = -89.5+y


    #-------------------------------------------
    #FINDING THE CLOSEST LATITUDE/LONGITUDE TO THULE
    #--------------------------------------------
    lat_data          = np.squeeze(lat)
    lon_data          = np.squeeze(lon - 360.)

    LatLon = [zip(x,y) for x,y in zip(*(lon_data, lat_data))]
    LatLon = np.squeeze(np.asarray(LatLon))

    coord = (Loclon, Loclat)
    LatLon = np.reshape(LatLon, (360*180,2))

    LatLonclose  = findCls2(LatLon, coord)
    print 'Lon and Lat close to TAB ({}, {}): {}'.format(Loclon, Loclat, LatLonclose)

    #LatLon       = np.reshape(LatLon, (360, 180, 2))
    #indsLoc      = np.where(LatLon == LatLonclose)[0][0]

    #-------------------------------------------
    #DEFINING THE AREA TO CALCULATE TIME SERIES
    #--------------------------------------------
    inside = np.logical_and(np.logical_and(lon >= west,
                                           lon <= east),
                            np.logical_and(lat >= south,
                                           lat <= north))

    lat2 = lat[inside]
    lon2 = lon[inside]
    Inipoly = [ (lon2[0], lat2[0]), (lon2[0], lat2[-1]), (lon2[-1], lat2[-1]), (lon2[-1], lat2[0])]

    poly = Polygon(Inipoly)

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
    #READING YEARLY FILES
    #-------------------------------------------
    for y, year in enumerate(years):

        fname = DataDir + 'oiv2mon.'+str(year)+'.asc'
        (iyrst,imst,idst,iyrnd,imnd,idnd, sst, ice) = rsstice.fib(fname)

        sst = np.asarray(sst)
        ice = np.asarray(ice)
        da = [dt.date(yy, mm, 15) for yy, mm in zip(iyrst,imst) ]
        Dates.append(da)

        for mm in imst:
            iceY.append(np.mean(ice[inside,mm-1]))
            sstY.append(np.mean(sst[inside,mm-1]))

        indMonth = np.where(imst == moi)[0]

        sst      = sst[:,:,indMonth]
        ice      = ice[:,:,indMonth]

        sst = np.reshape(sst, (360,180))
        ice = np.reshape(ice, (360,180))

        icematrix[y,:,:] = ice
        sstmatrix[y,:,:] = sst

    #-------------------------------------------
    #ANOMALIES
    #-------------------------------------------
    sstMean = np.mean(sstmatrix, axis=0)
    iceMean = np.mean(icematrix, axis=0)
    #iceMean[np.where(itag==0)] = float('NAN')

    for y, year in enumerate(years):
        iceAnomaly[y,:,:] = icematrix[y,:,:] - iceMean[:,:]
        sstAnomaly[y,:,:] = sstmatrix[y,:,:] - sstMean[:,:]

    iceY = np.asarray(iceY)
    sstY = np.asarray(sstY)
    Dates = np.asarray(Dates)
    Dates = np.reshape(Dates, (Dates.shape[0]*Dates.shape[1]))

    
    month    = np.array([d.month for d in Dates])
    mnthSort = list(set(month))
    mnthMean = np.zeros(len(mnthSort))


    for i,m in enumerate(mnthSort):
        inds        = np.where(month == m)[0]
        iceY[inds]  = iceY[inds] - np.mean(iceY[inds])
        sstY[inds]  = sstY[inds] - np.mean(sstY[inds])

    #---------------------------------------------------------------------------------------------
    #
    #---------------------------------------------------------------------------------------------
    R_M       = 6371.008771415  # mean (spherical) earth radius in kilo meter
    theta     = radians(90.0 - avgSZA)
    xdist     = polyH/tan(theta)
    print 'Horizontal distance (assuming Trop Height of {0:.1f}km and sza of {1:.2f}deg): {2:.1f}km'.format(polyH, avgSZA, xdist )

    Lonpoly   = []
    Latpoly   = []

    for a in saa:
        xlat_i, xlon_i = destination(Loclat, Loclon, xdist, a, radius=6371.008771415)
        Lonpoly.append(xlon_i)
        Latpoly.append(xlat_i)
    
    Inipoly   = []
    Inipoly.append((Loclon, Loclat))

    for (lo,la) in zip(Lonpoly,Latpoly):
        Inipoly.append( (lo,la))


    poly = Polygon(Inipoly)

    
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


    #-------
    pltMap(lon, lat, icematrix[indyear,:,:], type="ice", siteloc=[Loclon, Loclat], saveFlg=pdfsav, startYear=startYear, endYear=endYear, yoi=yoi, moi=moi, globalmapFlg=False, polyFlg=False, saa=saa, xdist=xdist)
    
    pltMap(lon, lat, sstmatrix[indyear,:,:], type="sst", siteloc=[Loclon, Loclat], saveFlg=pdfsav, startYear=startYear, endYear=endYear, yoi=yoi, moi=moi, globalmapFlg=False, polyFlg=False, saa=saa, xdist=xdist)
    
    pltMap(lon, lat, iceAnomaly[indyear,:,:], type="iceanomaly", siteloc=[Loclon, Loclat], saveFlg=pdfsav, startYear=startYear, endYear=endYear, yoi=yoi, moi=moi, globalmapFlg=False, polyFlg=False, saa=saa, xdist=xdist)
    
    pltMap(lon, lat, sstAnomaly[indyear,:,:], type="sstanomaly",  siteloc=[Loclon, Loclat], saveFlg=pdfsav, startYear=startYear, endYear=endYear, yoi=yoi, moi=moi, globalmapFlg=False, polyFlg=False, saa=saa, xdist=xdist)

    data    = icematrix[indyear,inside]
    pltMap(lon2, lat2, data, type="area", saveFlg=pdfsav, startYear=startYear, endYear=endYear, yoi=yoi, moi=moi, globalmapFlg=False, polyFlg=False, saa=saa, xdist=xdist)

    #pltMap(lon, lat, icematrix, indyear, type="ice", siteloc=[Loclon, Loclat], saveFlg=pdfsav, startYear=startYear, endYear=endYear, yoi=yoi, moi=moi, globalmapFlg=True, polyFlg=False)

    #-------------------------------------------
    #PLOT OF TIME SERIES
    #-------------------------------------------
    #inds2 = np.where( (month >= moi) & (month <=6) )[0]
    inds2 = np.where( month == moi )[0]
    yearsall = [ singDate.year for singDate in Dates[inds2]]

    indI = np.arange(len(inds2)) 
    N = len(indI)


    Colors = []
    for val in sstY[inds2]:
        if val <= 0: Colors.append('blue')# chocolate
        elif val > 0.0: Colors.append('red')

    #fig1 ,(ax1, ax2) = plt.subplots(2, sharex=True)
    fig1 , (ax1, ax2) = plt.subplots(2, figsize=(10,8), sharex=True)

    plt.suptitle('SST and Sea Ice Concentration Time series [Month={}]'.format(moi), fontsize=14)
    sc = ax1.bar(indI,sstY[inds2],color=Colors, align='center', edgecolor="black")
    #ax1.plot(Dates[inds2],sstY[inds2],'k.',markersize=4)
    ax1.grid(True)
    ax1.grid(True)
    ax1.set_ylabel('SST anomaly [C]',multialignment='center', fontsize=14)
    ax1.tick_params(which='both',labelsize=14)
    ax1.set_ylim(-1.5, 1.0)

    Colors = []

    for val in iceY[inds2]:
        if val <= 0: Colors.append('blue')# chocolate
        elif val > 0.0: Colors.append('red')
    
    sc2 = ax2.bar(indI,iceY[inds2], color=Colors, align='center', edgecolor="black")
    ax2.set_xticks(np.arange(0, N, 2))
    ax2.grid(True)
    ax2.set_ylabel('Sea Ice Concentration \n anomaly [%]',multialignment='center', fontsize=14)

    xt = ax2.get_xticks()
    new_xticks = [yearsall[d] for d in xt]
    #ax2.set_xticklabels(new_xticks, rotation=45)
    plt.xticks(indI, years, rotation=45)
    ax2.tick_params(which='both',labelsize=14)

    ax2.set_xlim(indI[0]-1,indI[-1]+1)
    ax2.set_ylim(-30, 30)
    
    ax2.set_xlabel('Year',  fontsize=14)
    #ax2.set_title('Time Series of Retrieved Total Column\n[molecules cm$^{-2}$]',multialignment='center')
    plt.subplots_adjust(bottom=0.12, left=0.1, right=0.95, top=0.95)
    

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

