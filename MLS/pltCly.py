#! /usr/bin/python2.7
#----------------------------------------------------------------------------------------
# Name:
#        pltCly.py
#
# Purpose:
#        Plot Cly species from GOZCARDS, AGAGE, JFJ, and MLO
#
# Notes:
#   
#
# Version History:
#       Created, Sep, 2017  Ivan Ortega (iortega@ucar.edu)
#----------------------------------------------------------------------------------------

                                    #-------------------------#
                                    # Import Standard modules #
                                    #-------------------------#

import sys
import os
import getopt
import ClyClass as dc
import time
from matplotlib.backends.backend_pdf import PdfPages


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

def main():

	#-------------------------------------
    #INITIALIZATIONS 
    #-------------------------------------
    MrgDir        = '/data1/Campaign/Satellite/MLS/Merge'                                            #MERGE DIR (GOZCARDS)
    #AgaDir        = '/data1/Campaign/Satellite/MLS/Agage_gcmd_gcms.data/gc-md/monthly'               #AGAGE DIR
    AgaDir        = '/data1/Campaign/Satellite/MLS/Agage_gcmd_gcms.data'

    jfjFile       = '/data1/Campaign/Satellite/MLS/HCl_ClONO2_Cly_Jungfraujoch_FTIR_1983-2016.txt'   #JFJ FILE
    mloFile       = '/data1/ebaumer/mlo/hcl/MLO-FTS-HCL.ascii'                                       #MLO FILE

    iyear         = 1991                                                                             #INITIAL YEAR TO READ DATA (Mainly for GOZCARDS, and MLS)
    fyear         = 2017                                                                             #FINAL YEAR TO READ DATA

    saveFlg       = True                                                                             #IF TRUE SAVE PDF (see name below)
    
    pltFile = '/data1/Campaign/Satellite/MLS/Cly_all.pdf'
    
	#---------------------------------
    # Check for the existance of files 
    # directories from input file
    #---------------------------------

    ckDir(MrgDir, exit=True)
    ckDir(AgaDir, exit=True)
    ckFile(jfjFile, exit=True)
    ckFile(mloFile, exit=True)

   
    if saveFlg:  
    	ckDir(os.path.dirname(os.path.realpath(pltFile)), exit=True)

    #-------------------------
    # Create Instance of Class
    #-------------------------
    gas = dc.PlotCly(MrgDir, AgaDir, jfjFile, mloFile, iyear=iyear, fyear=fyear, saveFlg=saveFlg, outFname=pltFile)
    gas.PlotClySet()
    if saveFlg: gas.closeFig()
    
    print('\nFinished Plots.......\n')

    #--------------------------------
    # Pause so user can look at plots
    #--------------------------------
    if not saveFlg:
        user_input = raw_input('Press any key to exit >>> ')
        sys.exit()           # Exit program     


if __name__ == "__main__":
    main()