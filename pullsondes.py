#! /usr/local/python-2.7/bin/python

#----------------------------------------------------------------------------------------
# Name:
#        pullsondes.py
#
# Purpose:
#       This program pulls water vapor sondes from NOAA
#           
#
# Notes:
#       1) Command line arguments
#            -d   YYYY       : Specify the start and stop date to get data.
#                                         Default is previous utc day
#
# Usage:
#     pullsondes.py 
#
# Examples:
#    ./pullsondes.py 
#
# Version History:
#  1.0     Created, November, 2013  Ivan Ortega
#
#----------------------------------------------------------------------------------------


                            #-------------------------#
                            # Import Standard modules #
                            #-------------------------#

import datetime as dt
import sys
import os
import itertools
import getopt
import shutil
import subprocess as sp
import logging

                            #-------------------------#
                            # Define helper functions #
                            #-------------------------#
                            
def usage():
    ''' Prints to screen standard program usage'''
    print 'pullsondes.py [-d YYYY]'                            

def subProcRun( sysCall, logF=False, shellFlg=False ):
    '''This runs a system command and directs the stdout and stderr'''
    rtn = sp.Popen( sysCall, stdout=sp.PIPE, stderr=sp.PIPE, shell=shellFlg )
    stdoutInfo, stderrInfo = rtn.communicate()

    if logF:
        if stdoutInfo: logF.info( stdoutInfo )
        if stderrInfo: logF.error( stderrInfo )
               
    return (stdoutInfo,stderrInfo)

def chMod(PrntPath):
    for dirpath, dirnames, filenames in os.walk(PrntPath):
        try:    os.chmod(dirpath,0o777)
        except: pass        
        for filename in filenames:
            path = os.path.join(dirpath, filename)
            try:    os.chmod(path, 0o777)   
            except: pass    

def ckDirMk(dirName,logFlg=False):
    ''' '''
    if not ( os.path.exists(dirName) ):
        os.makedirs( dirName, mode=0777 )
        if logFlg: logFlg.info( 'Created folder {}'.format(dirName))
        return False
    else:
        return True
                            #----------------------------#
                            #                            #
                            #        --- Main---         #
                            #                            #
                            #----------------------------#

def main(argv):

    dateFlg = False
    mloFlg  = True
    tabFlg  = True
                                                #---------------------------------#
                                                # Retrieve command line arguments #
                                                #---------------------------------#
    #------------------------------------------------------------------------------------------------------------#                                             
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'd:')

    except getopt.GetoptError as err:
        print str(err)
        usage()
        sys.exit()
        
    #-----------------------------
    # Parse command line arguments
    #-----------------------------
    for opt, arg in opts:
        
        #-----------
        # Start Date
        #-----------
        if opt == '-d':           
            
            if not arg:
                print 'Date format must be YYYY'
                sys.exit()
                
            dateFlg = True
            dates   = arg.strip()

            iyear   = int(dates[0:4])
            imnth   = int(1)
            iday    = int(1)

       # elif opt == '-s':
       #     print 'site: ' + opt
       #     site = opt   
        #------------------
        # Unhandled options
        #------------------
        else:
            print 'Unhandled option: ' + opt
            usage()
            sys.exit()
    #------------------------------------------------------------------------------------------------------------#    
    
    #-------------------
    # Set some constants
    #-------------------
    # 
    SondeDir = '/data/iortega/test'
    
    # Log File name
    logFname     = '/data/iortega/test.log'
    
    #------------------------------------
    # Get the current date and time (UTC)
    #------------------------------------
    curntDateUTC = dt.datetime.utcnow()
        
    #--------------------------------------------
    # Get date list:
    # -- If command line arguments are used then
    #    datelist corresponds to list of days 
    #    given in command line
    # -- If no command line arguments then create
    #    list starting from previous UTC date
    #    extending back ndaysbck
    #--------------------------------------------
    if dateFlg: curntDateUTC = dt.date(iyear,imnth,iday)    

    #------------------------------------
    # Start log information for this pull
    #------------------------------------
    logFilePull = logging.getLogger('1')
    logFilePull.setLevel(logging.INFO)
    hdlr1   = logging.FileHandler(logFname, mode='a+')
    fmt1    = logging.Formatter('%(asctime)s %(levelname)-8s %(message)s','%a, %d %b %Y %H:%M:%S')
    hdlr1.setFormatter(fmt1)
    logFilePull.addHandler(hdlr1)  
    logFilePull.info('**************** Starting Logging for sonde ***********************')
    logFilePull.info('Current Time (UTC):        ' + str(curntDateUTC) )    

    #----------
    # Pull Data
    #----------
    #-----------------
    # Get Current Year
    #-----------------
    yrstr = '{:4d}'.format(curntDateUTC.year)
    
    # Check if directories exist
    ckDirMk(SondeDir,logFlg=logFilePull)
    
    
    cmnd = ['wget','-r','-l1','-nd','-N','--no-parent','-A*'+yrstr+'*.txt', '-P'+SondeDir, 'ftp://aftp.cmdl.noaa.gov/data/ozwv/WaterVapor/Boulder_New/'  ]
    (stdoutInfo,stderrInfo) = subProcRun( cmnd, logFilePull )
    chMod(SondeDir)    

    cmnd = ['wget','-r','-l1','-nd','-N','--no-parent','-A*'+yrstr+'*.txt', '-P'+SondeDir, 'ftp://aftp.cmdl.noaa.gov/data/ozwv/WaterVapor/Hilo_New/'  ]
    (stdoutInfo,stderrInfo) = subProcRun( cmnd, logFilePull )
    chMod(SondeDir)    


    
if __name__ == "__main__":
    main(sys.argv[1:])
