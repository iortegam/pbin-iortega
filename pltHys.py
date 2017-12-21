#! /usr/bin/python2.7
#----------------------------------------------------------------------------------------
# Name:
#        pltHys.py
#
# Purpose:
#        Plot outputs generated with Hysplit4
#
# Notes:
#        Initial Hysplit was installed in the Mac (Users/iortega/Hysplit4/), apple based hysplit
#        Hysplit is run within the guicode folder usong ./hysplit4.tcl
#
# Version History:
#       Created, Dec, 2016  Ivan Ortega (iortega@ucar.edu)

#----------------------------------------------------------------------------------------
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
import datetime as dt
import sys
import matplotlib.gridspec as gridspec
from matplotlib.dates import DateFormatter, MonthLocator, YearLocator
from matplotlib.ticker import FormatStrFormatter, MultipleLocator,AutoMinorLocator,ScalarFormatter
from matplotlib.backends.backend_pdf import PdfPages #to save multiple pages in 1 pdf...
import os
import glob

from math import acos, asin, atan2, cos, hypot, sin, sqrt
from math import degrees, pi as PI, radians, sin, tan

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

    #-------------------
	#PATHS
    #-------------------
    dataDir    = '/data/iortega/results/tab/Hysplit/tdump/'
    pltDir     = '/data/iortega/results/tab/Hysplit/Fig/'
    saveFlg    = True
    DOI        = ['0525', '0527','0528','0529','0530','0601','0605','0606','0607','0609','0610','0611','0616','0617','0618','0619','0627']   #Days Of Interest
    #DOI        = ['0605']
   
    #-------------------
    #Making sure Directories exist
    #-------------------    
    ckDir(dataDir, exit=True)
    ckDir(pltDir, exit=True)

	#-------------------
    #specify format for numpy recfromtxt routine
    #-------------------
    nam = ['traj', 'gridnb', 'year', 'month', 'day', 'hour', 'minute', 'fhour', 'ageoftraj', 'lat', 'lon', 'height', 'pres']
	#-------------------specify format strings for each column
    form = ['<i4', '<i4', '<i4', '<i4', '<i4', '<i4', '<i4', '<i4', '<f8', '<f8', '<f8', '<f8', '<f8']
	
	#-------------------
    #Creating list of files
    #-------------------
    DirFiles = glob.glob(dataDir + '*'+'tdump*')
    DirFiles.sort()

    drssplt = [drs.split('/')[-1][5:11] for drs in  DirFiles]
    drssplt = np.asarray(drssplt)
    
    if saveFlg:
    	outFname  = pltDir + 'plt_Hysplit4_v2.pdf'
    	pdfsav = PdfPages(outFname)

    lat    = 76.52   #THULE
    lon    = -68.7   
    dist   = 1000.0  #DISTANCE OF THE LOS
    #saa    = [11.48, 12.65, 55.35, 64.45, 108.02, 109.41, 159.58, 161.02, 213.67, 215.08]
    saa    = [11.48, 55.35, 64.45, 108.02, 161.02, 215.08]

    #sza    = [79.785, 79.724, 73.88, 71.92, 61.705, 61.403, 53.833, 53.734, 55.083, 55.251]
    sza    = [79.785, 73.88, 71.92, 61.705, 53.734, 55.251]
    R_M    = 6371.008771415  # mean (spherical) earth radius in kilo meter

    #lat2, lon2 = destination(lat, lon, dist, saa, radius=6371.008771415)
    #h     = tan(radians(sza))*dist

    #-------------------
    #MAP WITH BEARING
    #-------------------

  #   for d in DOI:

  #   	clmap        = 'jet'
  #   	cm           = plt.get_cmap(clmap)   
  #   	fig, ax1  = plt.subplots(1, figsize=(7,6))

  #   	ax1.set_title('Line of Sight in Thule for different SAA')

  #   	map = Basemap(projection='npstere',boundinglat=65,lon_0=270,  resolution='l',  ax=ax1)
  #   	map.fillcontinents(color = 'gray', lake_color='gray', alpha=0.25)
  #   	map.drawcoastlines()
  #   	map.drawparallels(np.arange(-80., 81., 20.), alpha=0.25)
  #   	map.drawmeridians(np.arange(-180., 181., 20.), alpha=0.25)
  #   	map.drawmapboundary(fill_color= 'white')
  #   	#map.drawmapscale(-105., 70, -105.0, 75, 30, barstyle='fancy')


  #   	for (a,s) in zip(saa,sza):

  #   		lat2, lon2 = destination(lat, lon, dist, a, radius=6371.008771415)
  #   		map.drawgreatcircle(lon,lat,lon2,lat2,linewidth=2,color='m')
	 #        x2, y2 = map(lon2, lat2)
	 #        ax1.plot(x2, y2, marker='D', color='m', markersize=7)

	 #        #ax1.annotate(str(a), xytext=(lon2, lat2), textcoords='data')
	 #        ax1.text(x2, y2, str(a), fontsize=14, color='m')



		# x, y = map(lon, lat)
		# ax1.plot(x, y, marker='D', color='b', markersize=7)

		# if saveFlg:
		# 	pdfsav.savefig(fig,dpi=200)
  #           #--------For Single Figures
  #           #outFname  = pltDir + drssplt[-1]+'.pdf'
  #           #plt.savefig(outFname, bbox_inches='tight')
  #           #--------
  #       else:       
  #           plt.show(block=False)
  #           #user_input = raw_input('Press any key to exit >>>')
  #          # sys.exit()
  #           plt.close()

    latMin   = 60.   # 54.
    latMax   = 88.    # 88.
    lonMin   = -90.   # -65.
    lonMax   = -25.    # 20

    #-------------------
    #HYSPLIT: Reading Files
    #-------------------
    #for drs in DirFiles:
    for d in DOI:

    	indS = np.where(drssplt == '16'+d )[0]

    	for dd in indS:

	        hysplit = np.recfromtxt(DirFiles[dd], skip_header=16, dtype={'names':nam, 'formats':form}, delimiter=(6,6,6,6,6,6,6,6,8,9,9,9,9))

	        date = []
	        for y,m,d,h,mi in zip(hysplit.year, hysplit.month, hysplit.day, hysplit.hour, hysplit.minute):
	        	date.append(dt.datetime(y+2000,m,d,h,mi))
	        date = np.array(date)

	        clmap        = 'jet'
	        cm           = plt.get_cmap(clmap)              
	        
	        months       = MonthLocator()
	        DateFmt      = DateFormatter('%H\n%d/%m/%Y') 
	        DateFmt      = DateFormatter('%H\n%d/%m')

	        #-------------------
	        #Plot Map and Height in same Figure
	        #-------------------
	        fig, (ax1, ax2)  = plt.subplots(2) # figsize=(10, 10)

	        gs   = gridspec.GridSpec(2, 1,height_ratios=[3,1], wspace=0.05, hspace=0.05)
	        ax1  = plt.subplot(gs[0], sharex=ax2)

	        ax2  = plt.subplot(gs[1])

	        map = Basemap(llcrnrlat=latMin,urcrnrlat=latMax,
    	            llcrnrlon=lonMin,urcrnrlon=lonMax,
    	            rsphere=(6378137.00,6356752.3142),
    	            resolution='l',area_thresh=1000.,projection='lcc',
    	            lat_1=latMin,lon_0=-60, ax=ax1)

	        #map = Basemap(projection='npstere',boundinglat=40,lon_0=270,  resolution='l',  ax=ax1)


	        #map.drawgreatcircle(lon,lat,lon2,lat2,linewidth=2,color='m')

	        map.fillcontinents(color = 'gray', lake_color='gray', alpha=0.25)

	        map.drawcoastlines(color='black')
	        map.drawcountries(color='lightgrey')

	        map.drawparallels(np.arange(-80., 81., 5.), labels=[1,0,0,0], alpha=0.25)
	        map.drawmeridians(np.arange(-180., 181., 10.),  labels=[0,0,0,1], alpha=0.25)
	        #map.drawmapboundary(fill_color= 'white')
	        #map.drawmapscale(-105., 70, -105.0, 75, 30, barstyle='fancy')

	        # x, y = map(lon, lat)
	        # ax1.plot(x, y, marker='D', color='m', markersize=7)

	        # x2, y2 = map(lon2, lat2)
	        # ax1.plot(x2, y2, marker='D', color='m', markersize=7)
	        
	        #calculate map coordinates and plot all trajectories 1 to 28
	        T = {}
	        H = {}
	        for i in range(1,4):

	         	T[str(i)] = (hysplit.lat[hysplit.traj==i], hysplit.lon[hysplit.traj==i])
	         	x,y = map(T[str(i)][1],T[str(i)][0])
	         	ax1.plot(x, y, '-', linewidth=3)
	         	H[str(i)] = (hysplit.height[hysplit.traj==i], date[hysplit.traj==i])
	        	
	         	ax2.plot(H[str(i)][1], H[str(i)][0]/1000., '-', linewidth=2)
	         	ax2.grid()

	        ax2.xaxis.set_major_formatter(DateFmt)
	        ax2.set_ylabel('Height [km]', fontsize=14)
	        ax2.set_xlabel('Date', fontsize=14)
	        ax1.set_title('Backward trajectories ending at {0:02} UTC on {1:}/{2:}/{3:}'.format(date[0].hour, date[0].day, date[0].month, date[0].year), fontsize=14)
	        #plt.suptitle('Backward trajectories ending at {0:02} UTC on {1:}/{2:}/{3:}'.format(date[0].hour, date[0].day, date[0].month, date[0].year), fontsize=16)
	        ax1.tick_params(which='both',labelsize=14)
	        ax2.tick_params(which='both',labelsize=14)

	        fig.subplots_adjust(bottom=0.12,top=0.92, left=0.1, right=0.95)

	        if saveFlg:
	        	pdfsav.savefig(fig,dpi=200)
	            #--------For Single Figures
	            #outFname  = pltDir + drssplt[-1]+'.pdf'
	            #plt.savefig(outFname, bbox_inches='tight')
	            #--------
	        else:       
	            plt.show(block=False)
	            user_input = raw_input('Press any key to exit >>>')
	            #sys.exit()
	            plt.close()
            #sys.exit()
    if saveFlg: pdfsav.close()

 
if __name__ == "__main__":
    main()