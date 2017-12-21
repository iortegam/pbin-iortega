#! /usr/bin/python2.7
#----------------------------------------------------------------------------------------
# Name:
#        pltMLS.py
#
# Purpose:
#  
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
import MLSHDFClass as dc
import time
from matplotlib.backends.backend_pdf import PdfPages
import myfunctions as mf


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

    pltMLS        = False
    pltMrg        = False
    pltAga        = False 

    pltAll        = True

    MLSDir        = '/data1/Campaign/Satellite/MLS/hcl'
    MrgDir        = '/data1/Campaign/Satellite/MLS/Merge'
    AgaDir        = '/data1/Campaign/Satellite/MLS/Agage_gcmd_gcms.data/gc-md/monthly'

    iyear         = 1991
    fyear         = 2017

    pltFileMLS    = '/data1/Campaign/Satellite/MLS/HCl_MLS.pdf'
    pltFileMrg    = '/data1/Campaign/Satellite/MLS/HCl_Merge.pdf'
    pltFileAga    = '/data1/Campaign/Satellite/MLS/AGAGE_Cl.pdf'
    
    pltFileAll = '/data1/Campaign/Satellite/MLS/HCl_All.pdf'
    
    saveFlg        = True

	#---------------------------------
    # Check for the existance of files 
    # directories from input file
    #---------------------------------
    ckDir(MLSDir, exit=True)
    ckDir(MrgDir, exit=True)
   
    if saveFlg:  
    	ckDir(os.path.dirname(os.path.realpath(pltFileMLS)),exit=True)
    	ckDir(os.path.dirname(os.path.realpath(pltFileMrg)),exit=True)
        
    #-------------------------
    # Create Instance of Class
    #-------------------------
    if pltMLS: 
        gas = dc.PlotMLS(MLSDir, iyear=iyear, fyear=fyear, saveFlg=saveFlg, outFname=pltFileMLS)
        gas.PlotMLSSet()
        if saveFlg: gas.closeFig()

    if pltMrg: 
        gas = dc.PlotMrg(MrgDir, iyear=iyear, fyear=fyear, saveFlg=saveFlg, outFname=pltFileMrg)
        gas.PlotMrgSet()
        if saveFlg: gas.closeFig()

    if pltAga: 
        gas = dc.PlotAga(AgaDir, saveFlg=saveFlg, outFname=pltFileAga)
        gas.PloAgaSet()
        if saveFlg: gas.closeFig()

    if pltAll:
    	gas = dc.PlotAll(MLSDir, MrgDir, AgaDir, iyear=iyear, fyear=fyear, saveFlg=saveFlg, outFname=pltFileAll)
        gas.PlotAllSet()
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