#! /usr/bin/python2.7
##! /usr/local/python-2.7/bin/python
#----------------------------------------------------------------------------------------
# Name:
#        read_ncdf_WRF.py
#
# Purpose:
#       The purpose of this program is to read netCDF files from Gabi Pfister during FRAPPE
#
# Notes:
#   
#
# Version History:
#       Created, May, 2016  Ivan Ortega (iortega@ucar.edu)

#----------------------------------------------------------------------------------------

                                    #-------------------------#
                                    # Import Standard modules #
                                    #-------------------------#
from scipy.io import netcdf
import os
import datetime as dt
import numpy as np
import sys
import glob
                                    #-------------------------#
                                    # Define helper functions #
                                    #-------------------------#

def ckDir(dirName):
    '''Check if a directory exists'''
    if not os.path.exists( dirName ):
        print 'Input Directory %s does not exist' % (dirName)
        sys.exit()

def ckFile(fName):
    '''Check if a file exists'''
    if not os.path.isfile(fName):
        print 'File %s does not exist' % (fName)
        sys.exit()


                                    #----------------------------#
                                    #                            #
                                    #        --- Main---         #
                                    #                            #
                                    #----------------------------#

def main():

    #----------------
    # Initializations
    #----------------
    dataDir     = '/data1/ancillary_data/fl0/FRAPPE/'
    outDataDir  = '/data1/ancillary_data/fl0/FRAPPE/'
    file1       = 'BoulderFL.nc'
    file2       = 'BoulderFL_loc.nc'


    #-------------------------------------------------
    #Reading both standard files
    #-------------------------------------------------
    indfile1 = dataDir + file1
    indfile2 = dataDir + file2
    
    ckDir(indfile1)
    ckDir(indfile2)

    cdfname1 = netcdf.netcdf_file(indfile1,'r',mmap=False)       # Open netcdf file 1
    cdfname2 = netcdf.netcdf_file(indfile2,'r',mmap=False)       # Open netcdf file 2

    #-------------------------------------------------
    #To see list of variable/dimensions
    #-------------------------------------------------
    print cdfname1.dimensions.keys()
    print cdfname1.variables.keys()

    #print cdfname2.dimensions.keys()
    #print cdfname2.variables.keys()

    #-------------------------------------------------
    #Get variables from file 2
    #-------------------------------------------------
    cosalpha          = cdfname2.variables['COSALPHA'][:]
    sinalpha          = cdfname2.variables['SINALPHA'][:]   

    #-------------------------------------------------
    # Get variables from file 1
    #-------------------------------------------------
    time        = cdfname1.variables['Times'][:]

    a1          = cdfname1.variables['PB'][:]
    a2          = cdfname1.variables['P'][:]   
    pressure    = (a1 + a2)/100.0
    a1          = cdfname1.variables['HGT'][:]
    a2          = cdfname1.variables['PHB'][:]
    a3          = cdfname1.variables['PH'][:]
    a2          = cdfname1.variables['P'][:]

    Ttmp        = cdfname1.variables['T'][:]
    p00         = cdfname1.variables['P00'][:]

    print pressure.shape
    print Ttmp.shape
    print p00.shape

    Temperature = (Ttmp+300.0)* ((pressure)/p00)0.286

    #print len(time[:, 0])
    #npnts       = len(time)
    #print npnts 




    exit()

#for itime=0,ntimes -1 do begin  ; ntimes is number of time entries in file
#for i=0,nz-1 do begin  ; nz is number of vertical levels
#    atemp = (a2[*,*,i, itime]+a3[*,*,i,itime])/9.8/1000.   
#    atemp2 = (a2[*,*,i+1,itime]+a3[*,*,i+1,itime])/9.8/1000.
#   altitude[*,*,i,itime]= (atemp+atemp2)/2.
#endfor

    exit()


    #--------------------------
    # Initialize variable lists
    #--------------------------
    tempList   = []
    timeList   = []
    rhList     = []
    pressList  = []
    tempDPList = []
    rmind      = []

    #--------------------
    # Search through base directory for files
    #--------------------
    files = glob.glob(dataDir + 'BoulderFL_loc*.' + fileExtTag)
    if not files:
        print 'No files found'
        sys.exit()
    else:
        print ' %d files found ' %len(files)

    #-------------------------
    # Loop through found files
    #-------------------------
    for indvfile in files:

        cdfname = netcdf.netcdf_file(indvfile,'r',mmap=False)       # Open netcdf file

        print cdfname.dimensions.keys()
        print cdfname.variables.keys()
        #print cdfname.attributes.keys()
        #print cdfname.variables['hcho']

        # Get variables
        a1          = cdfname.variables['PB'][:]
        a2          = cdfname.variables['P'][:]
        
        
        pressure    = (a1 + a2)/100.0

        print pressure

        exit()


        base_time   = cdfname.variables['base_time']
        time_offset = cdfname.variables['time_offset']
        temp        = cdfname.variables['tdry']
        rh          = cdfname.variables['rh']
        press       = cdfname.variables['pres']
        tempDP      = cdfname.variables['dp']
        wdir      = cdfname.variables['wdir']
        wspd      = cdfname.variables['wspd']
        wmax      = cdfname.variables['wmax']
        wsdev      = cdfname.variables['wsdev']
        cdfname.close()

        #----------------------------------
        # Create an actual time vector from
        # basetime and offset (unix time)
        #----------------------------------
        total_time = base_time[()] + time_offset.data
        #total_time = [dt.datetime.utcfromtimestamp(indtime) for indtime in total_time]
        total_time = total_time.tolist()                # Convert numpy array to list

        #------------------------------------------------------------
        # There seems to be a timing issue in some of the data files.
        # time_offset can have unusually large numbers. Need to check
        # for this and delete observations.
        #------------------------------------------------------------
        for ind, indvtime in enumerate(total_time):
            # Identify erroneous time_offset values
            try:
                total_time[ind] = dt.datetime.utcfromtimestamp(indvtime)
            except ValueError:
                total_time[ind] = -9999
                rmind.append(ind)

        #------------------------------------------------------
        # Remove observations with erroneous time_offset values
        #------------------------------------------------------
        total_time[:] = [item for ind,item in enumerate(total_time) if ind not in rmind]   # List
        temp.data     = np.delete(temp.data,rmind)                                         # Numpy array
        rh.data       = np.delete(rh.data,rmind)                                           # Numpy array
        press.data    = np.delete(press.data,rmind)                                        # Numpy array
        tempDP.data   = np.delete(tempDP.data,rmind)                                       # Numpy array

        #---------------------
        # Append to main lists
        #---------------------
        timeList.extend(total_time)
        tempList.extend(temp.data)
        rhList.extend(rh.data)
        pressList.extend(press.data)
        tempDPList.extend(tempDP.data)

    #------------------------
    # Sort list based on time
    # This returns a tuple
    #------------------------
    timeList, tempList, rhList, pressList, tempDPList = zip(*sorted(zip(timeList, tempList, rhList, pressList, tempDPList)))

    #-----------------------------------
    # Construct a vector of string years
    #-----------------------------------
    years  = ["{0:02d}".format(indtime.year) for indtime in timeList]
    months = ["{0:02d}".format(indtime.month) for indtime in timeList]
    days   = ["{0:02d}".format(indtime.day) for indtime in timeList]
    hours  = ["{0:02d}".format(indtime.hour) for indtime in timeList]
    minutes= ["{0:02d}".format(indtime.minute) for indtime in timeList]


    #--------------------------
    # Write data to output file
    #--------------------------
    with open(outDataDir+'fl0_met_data_'+years[0]+'.txt', 'w') as fopen:
        fopen.write('Year Month Day Hour Minute Temperature[C] RelativeHumidity[%] Pressure[mbars] DewPointTemperature[C]\n')
        for line in zip(years,months,days,hours,minutes,tempList,rhList,pressList,tempDPList):
            fopen.write('%-4s %-5s %-3s %-4s %-6s %-14.1f %-19.1f %-15.1f %-22.6f\n' % line)


if __name__ == "__main__":
    main()