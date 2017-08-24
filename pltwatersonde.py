#! /usr/local/python-2.7/bin/python
#! /usr/bin/python2.7
#----------------------------------------------------------------------------------------
# Name:
#        pltwatersonde.py
#
# Purpose:
#        The purpose of this program is to plot the water profiles obtained with the water sondes profiles
#
#----------------------------------------------------------------------------------------


                                    #-------------------------#
                                    # Import Standard modules #
                                    #-------------------------#
import os
import datetime as dt
import sys
import glob
import DateRange as dr
import SondeReader as sr
import csv
import matplotlib.pyplot as plt
import numpy as np
                                    #-------------------------#
                                    # Define helper functions #
                                    #-------------------------#

def ckDir(dirName):
    '''Check if a directory exists'''
    if not os.path.exists( dirName ):
        print 'Directory %s does not exist' % (dirName)
        return False
    else:
        return True
        
def ckFile(fName):
    '''Check if a file exists'''
    if not os.path.isfile(fName)    :
        print 'File %s does not exist' % (fName)
        sys.exit()


                                    #----------------------------#
                                    #                            #
                                    #        --- Main---         #
                                    #                            #
                                    #----------------------------#

def main():
                                    
    loc     = 'mlo'
    dataDir     = '/data1/ancillary_data/sondes/'+loc.lower()+'/'
    
    #------------------------------
    # Date Range of data to process
    #------------------------------
    # Starting 
    iyear = 2015               # Year
    imnth = 1                  # Month
    iday  = 1                  # Day
    
    # Ending
    fyear = 2015               # Year
    fmnth = 1                 # Month
    fday  = 28                 # Day
    
    #-------------------
    # Call to date class
    #-------------------
    DOI = dr.DateRange(iyear,imnth,iday,fyear,fmnth,fday)
   
    #----------------------------
    # Create list of unique years
    #----------------------------
    years = DOI.yearList()
    
    #------------------
    # Loop through year
    #------------------
    for indvYear in years:

    #--------------------------------------
        # Loop through days/folders within year
        #--------------------------------------
        daysYear = DOI.daysInYear(indvYear)        # Create a list of days within one year
        
         #--------------------------------------
        # Create sonde data dictionary instance
        #--------------------------------------       
        SondeData = sr.WaterReads()
       

        for indvDay in daysYear:
            
            # Find year month and day strings
            yrstr   = "{0:02d}".format(indvDay.year)
            mnthstr = "{0:02d}".format(indvDay.month)
            daystr  = "{0:02d}".format(indvDay.day)  

            #--------------------------------------------
            # Log files changed name after Jan 21st, 2009
            #--------------------------------------------
            if (statstr.lower() == 'mlo'):
            	loc = 'HIH_H2O_'+yrstr+mnthstr+daystr+'.txt'
            elif (statstr.lower() == 'tab'):
            	loc = 'BLD_H2O_'+yrstr+mnthstr+daystr+'.txt'
           
            #------------------------
            # Look for house log file
            #------------------------
            sondefile = glob.glob( dataDir + loc )
            #print sondefile
            
            if not sondefile:
                continue

            sondefile = sondefile[0]

            SondeData.SondeA(sondefile,indvDay.year,indvDay.month,indvDay.day)
            #print SondeData.data['Press_hPa']
            #plt.plot([1,2,3,4], [1,4,9,16], 'ro')
            #plt.show()
            Press = np.array(map(float, SondeData.data['Press_hPa']))
            Height = np.array(map(float, SondeData.data['Alt_km']))
            Temp = np.array(map(float, SondeData.data['Temp_degC']))
            H2Omr = np.array(map(float, SondeData.data['H2Omr_ppmv']))
            H2Osd = np.array(map(float, SondeData.data['H2Osd_ppmv']))
            

            grid = {}
            H2Ovcd = {}
            rho = {}    
            #for i in enumerate(Press):
            rho = (Press*100.00) / ((Temp+273.15)*8.314) * 6.022e23 /100 /100 /100
            
            for i in enumerate(rho):
                if rho == -999: 
                    rho = 'nan'

                #H2Ovcd[i] = 
                #grid[i] = Height[i+1] - Height[i]
   
            plt.plot(rho, Height, 'ro')
            plt.show()


if __name__ == "__main__":
    main()