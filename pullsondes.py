#! /usr/bin/python

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
import glob

                            #-------------------------#
                            # Define helper functions #
                            #-------------------------#
                            
def usage():
    ''' Prints to screen standard program usage'''
    print 'pullsondes.py [-d YYYY -s fl0]'                            

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

def pullFPH():

    dateFlg = False
                                                #---------------------------------#
                                                # Retrieve command line arguments #
                                                #---------------------------------#
    #------------------------------------------------------------------------------------------------------------#                                             
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'd:s:') 

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

        elif opt == '-s':

            if not arg:
                print 'Give a location (fl0 or mlo)'
                sys.exit()

            print 'site: ' + str(arg)
            loc = str(arg)   
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
    #SondeDir = '/data/iortega/test'
    SondeDir     = '/data1/ancillary_data/sondes/'+loc.lower()+'/'

    
    # Log File name
    logFname     = SondeDir+'log_pullsonde.log'

    if ckFile(logFname): os.remove(logFname)
    
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

    if loc.lower() == 'fl0': 
    
        cmnd = ['wget','-r','-l1','-nd','-N','--no-parent','-A*'+yrstr+'*.txt', '-P'+SondeDir, 'ftp://aftp.cmdl.noaa.gov/data/ozwv/WaterVapor/Boulder_New/'  ]
        (stdoutInfo,stderrInfo) = subProcRun( cmnd, logFilePull )
        chMod(SondeDir)

    elif loc.lower()  == 'mlo':  

        cmnd = ['wget','-r','-l1','-nd','-N','--no-parent','-A*'+yrstr+'*.txt', '-P'+SondeDir, 'ftp://aftp.cmdl.noaa.gov/data/ozwv/WaterVapor/Hilo_New/'  ]
        (stdoutInfo,stderrInfo) = subProcRun( cmnd, logFilePull )
        chMod(SondeDir) 

def pullFleOut():

    #------------------------------------------------------------------------------------------------------------#                                             
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'd:s:') 

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

        elif opt == '-s':

            if not arg:
                print 'Give a location (fl0 or mlo)'
                sys.exit()

            print 'site: ' + str(arg)
            loc = str(arg)   
        #------------------
        # Unhandled options
        #------------------
        else:
            print 'Unhandled option: ' + opt
            usage()
            sys.exit()

    SondeDir     = '/data1/ancillary_data/sondes/'+loc.lower()+'/'

    logFname     = SondeDir+'log_pullsonde_fleout.log'

    if ckFile(logFname): os.remove(logFname)

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


    FilesFPH = glob.glob(SondeDir+'/BLD_H2O*.txt')
    FilesFPH.sort()

    for i, f in enumerate(FilesFPH):

        fsplit = f.split('/')

        yearS  = str(fsplit[-1][8:12])
        monthS = str(fsplit[-1][12:14])
        dayS   = str(fsplit[-1][14:16])

        idstr = yearS + '_' + monthS + '_' + dayS


        if loc.lower() == 'fl0': 
    
            cmnd = ['wget', '-nc','-P'+SondeDir, 'ftp://aftp.cmdl.noaa.gov/data/ozwv/Ozonesonde/Boulder, Colorado/Native Resolution (60s, 7s, 1s)/'+'*'+idstr+'_fleout.dat']
            (stdoutInfo,stderrInfo) = subProcRun( cmnd, logFilePull)


        chMod(SondeDir)

        

def main(argv):

    #pullFPH()

    pullFleOut()

    
if __name__ == "__main__":
    main(sys.argv[1:])
