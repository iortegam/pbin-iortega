#! /usr/bin/python2.7
#----------------------------------------------------------------------------------------
# Name:
#        pltsonde.py
#
# Purpose:
#
# Plot sonde vs FTS
#
# Input:
#       input_pltsonde.py
#----------------------------------------------------------------------------------------

#---------------
# Import modules
#---------------
import sys
import os
import getopt
import dataOutplts as dc
import time
import datetime as dt
import matplotlib.dates as md
from matplotlib.dates import DateFormatter, MonthLocator, YearLocator
from matplotlib.ticker import FormatStrFormatter, MultipleLocator,AutoMinorLocator,ScalarFormatter
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm 
import myfunctions as mf
import numpy as np
import sondeclass as rs


def usage():
    ''' Prints to screen standard program usage'''
    print 'pltsonde.py -i <inputfile> -?'
    print '  -i <file> : Run pltsonde.py with specified input file'
    print '  -?        : Show all flags'

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

def main(argv):

    #--------------------------------
    # Retrieve command line arguments
    #--------------------------------
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'i:?')

    except getopt.GetoptError as err:
        print str(err)
        usage()
        sys.exit()

    #-----------------------------
    # Parse command line arguments
    #-----------------------------
    for opt, arg in opts:
        # Check input file flag and path
        if opt == '-i':

            pltInputs = {}

            ckFile(arg,exit=True)

            try:
                execfile(arg, pltInputs)
            except IOError as errmsg:
                print errmsg + ' : ' + arg
                sys.exit()

            if '__builtins__' in pltInputs:
                del pltInputs['__builtins__']               


        # Show all command line flags
        elif opt == '-?':
            usage()
            sys.exit()

        else:
            print 'Unhandled option: ' + opt
            sys.exit()

     
    #---------------------------------
    # Radiosonde Info
    #---------------------------------
    sonde = rs.Pltsonde(pltInputs['retDir'], pltInputs['ctlFile'], pltInputs['sondeDir'], pltInputs['loc'], iyear=pltInputs['iyear'],imnth=pltInputs['imnth'],iday=pltInputs['iday'],
    fyear=pltInputs['fyear'],fmnth=pltInputs['fmnth'],fday=pltInputs['fday'],outFname=pltInputs['pltFile'])

    #----------------------
    # Call to plot profiles
    #----------------------
    sonde.plt1(fltr=pltInputs['fltrFlg'],allGas=False,sclfct=pltInputs['sclfct'],sclname=pltInputs['sclfctName'],mnthFltr=pltInputs["mnths"],mnthFltFlg=pltInputs["mnthFlg"],
               errFlg=pltInputs['errorFlg'],minSZA=pltInputs['minSZA'],maxSZA=pltInputs['maxSZA'],maxRMS=pltInputs['maxRMS'],minTC=pltInputs['minTC'],maxTC=pltInputs['maxTC'],
               minDOF=pltInputs['minDOF'],maxCHI=pltInputs['maxCHI'],dofFlg=pltInputs['dofFlg'],rmsFlg=pltInputs['rmsFlg'],tcFlg=pltInputs['tcNegFlg'],
               pcFlg=pltInputs['pcNegFlg'],szaFlg=pltInputs['szaFlg'],cnvrgFlg=pltInputs['cnvrgFlg'],chiFlg=pltInputs['chiFlg'],tcMMflg=pltInputs['tcMMFlg'], 
               pCols=pltInputs['pCols'], smthFlg=pltInputs['smthFlg'], doiplt=pltInputs['doi'], loc=pltInputs['loc'], latFTS =pltInputs['latFTS'], lonFTS =pltInputs['lonFTS'] )

    #-------------------------------------------------------------------------
    #END PLOTS
    #-------------------------------------------------------------------------
    if pltInputs['saveFlg']: sonde.closeFig()

    print('\nFinished Plots.......\n')
       
 #    #--------------------------------
 #    # Pause so user can look at plots
 #    #--------------------------------
    if not pltInputs['saveFlg']:
        user_input = raw_input('Press any key to exit >>> ')
        sys.exit()           # Exit program       





if __name__ == "__main__":
    main(sys.argv[1:])