#----------------------------------------------------------------------------------------
#Name:
#        mklist.py
#
# Purpose:
#       This program is initially created to re-write a a new spec Data base using coincident date of measurements with sondes
#
# Notes:
#       it can be modified to read dates of other instruments such as satellites
#
# Created by :
#       Ivan Ortega, Jan 2018
#
#----------------------------------------------------------------------------------------


                        #-------------------------#
                        # Import Standard modules #
                        #-------------------------#
import sys
import os
import getopt
import datetime as dt
import glob
import shutil  
import matplotlib.pyplot as plt
import numpy as np
import csv
import dataOutClass as dc


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



                            #----------------------------#
                            #                            #
                            #        --- Main---         #
                            #                            #
                            #----------------------------#

def main():

                            #----------------------------#
                            #        --- Inputs---       #
                            #----------------------------#

    #--------------------------------------------
    # FTS inputs
    #--------------------------------------------
    loc          = 'fl0'

    if loc.lower() == 'mlo':
        spcdbFile    = '/data/Campaign/MLO/Spectral_DB/HRspDB_mlo_1995_2017.dat'
    elif loc.lower() == 'fl0':
        spcdbFile    = '/data/Campaign/FL0/Spectral_DB/CoaddspDB_fl0_2010_2016.dat'
    else: 
        print 'Error: Error in spcdbFile location!'
        exit()

    #--------------------------------------------
    # Sonde inputs
    #--------------------------------------------
    sondeDir     = '/data1/ancillary_data/sondes/'+loc.lower()+'/'

    CoiTime      = 60 #Minutes
    
    #--------------------------------------------
    # Range of dates
    #--------------------------------------------
    iyear      = 2010
    imnth      = 1
    iday       = 1
    fyear      = 2016
    fmnth      = 12
    fday       = 31

    #--------------------------------------------
    # Output
    #--------------------------------------------
    if loc.lower() == 'mlo':
        OutspcdbFile = '/data/Campaign/MLO/Spectral_DB/HRspDB_'+loc.lower()+'_'+str(iyear)+'_'+str(fyear)+'_Sonde'+'_'+str(CoiTime)+'Min.dat'
    elif loc.lower() == 'fl0':
        OutspcdbFile = '/data/Campaign/FL0/Spectral_DB/CoaddspDB_'+loc.lower()+'_'+str(iyear)+'_'+str(fyear)+'_Sonde'+'_'+str(CoiTime)+'Min.dat'
    else: 
        print 'Error: Error in Output spcdbFile location!'
        exit()

                            #----------------------------#
                            #        --- Start ---       #
                            #----------------------------#
    
 
    idate      = dt.date(iyear, imnth, iday)
    fdate      = dt.date(fyear, fmnth, fday)


    #--------------------------------------------
    # Read Initial spcd Data base
    #--------------------------------------------
    ckFile(spcdbFile,exit=True)

    dbInputs    = {}
    dbFltInputs = {}

 
    with open(spcdbFile,'rb') as fname:
        # Currently only reads white space delimited file. Need to make general****
        reader = csv.DictReader(fname, delimiter=' ',skipinitialspace=True)                   # Read csv file
        DBfieldNames = reader.fieldnames
        for row in reader:
            #------------------------------------------------------
            # DictReader will add a None key and value if there are
            # extra spaces at the end of database line.
            #------------------------------------------------------
            if None in row: del row[None]                    

            for col,val in row.iteritems():

                dbInputs.setdefault(col,[]).append(val)                                 # Construct input dictionary                  
                if '' in dbInputs: del dbInputs['']

    DatesFTS = []
    dtFTS    = []
    
    for ind,val in enumerate(dbInputs['Date']):

        valstr = val
        datestmp = dt.date(int(valstr[0:4]), int(valstr[4:6]),int(valstr[6:]))

        timestr = dbInputs['Time'][ind]
        dttmp  = dt.datetime(int(valstr[0:4]), int(valstr[4:6]),int(valstr[6:]), int(timestr[0:2]), int(timestr[3:5]),int(timestr[6:8]))

        DatesFTS.append(datestmp)
        dtFTS.append(dttmp)
     
    #DatesFTS = np.asarray(int(valstr[0:4]), int(valstr[4:6]),int(valstr[6:]), )
    DatesFTS = np.asarray(DatesFTS)
    dtFTS    = np.asarray(dtFTS)

    #--------------------------------------------
    # Create list of sonde Dates
    #--------------------------------------------
    dataDirSonde = sondeDir 
    ckDir(dataDirSonde,exit=True)

    if loc.lower() == 'mlo': 
        filename = 'HIH_H2O_*.txt'
    elif loc.lower() == 'fl0':
        filename = 'BLD_H2O_*.txt'
    else:
        print 'Error: Bad location'
        exit()

    sondefiles = glob.glob( dataDirSonde + '/'+ filename )

    sondefiles.sort()

    DatesSonde = []
    dtSonde    = []

    for sf in sondefiles:
        datestmp = sf.strip().split('/')[-1][8:16]

        with open(sf ,'r') as fopen: lines = fopen.readlines()
        info = [ row.strip().split() for row in lines if 'Water Vapor Flight Date' in row]
        info2 = [ row.strip().split() for row in lines if 'Number of header lines' in row] 
        info3 = [ row.strip().split() for row in lines if 'Number of variables' in row]

        nheaders = int(info2[-1][-1])
        nvariable = int(info3[-1][-1])

        lines[:] = [ row.strip().split() for row in lines[nheaders:-1] ] 
        
        #lines[:] = [ row.strip().split() for row in lines if len(row.strip().split()) == 13 and not 'to' in row and not 'Level' in row and not 'Traditionally' in row and not 'Number' in row]              

        npoints = len(lines)

        time = [row[7].split()[0] for row in info]

        time2 = time[0].split(':')
        
        hh = int(time2[0])
        mi = int(time2[1])
        se = int(time2[2])
        
        datestmp2 = dt.date(int(datestmp[0:4]), int(datestmp[4:6]),int(datestmp[6:8]))
        dttmp     = dt.datetime(int(datestmp[0:4]), int(datestmp[4:6]),int(datestmp[6:8]), hh, mi, se)

        if datestmp2 >= idate and datestmp2 < fdate:

            DatesSonde.append(datestmp2)
            dtSonde.append(dttmp)

    DatesSonde = np.asarray(DatesSonde)
    dtSonde    = np.asarray(dtSonde)
    #-------------------------------------------------------
    #FINDING COINCIDENT DATES
    #-------------------------------------------------------
    doy_fts   = dc.toYearFraction(DatesFTS)
    doy_sonde = dc.toYearFraction(DatesSonde)

    intrsctVals = np.intersect1d(doy_fts, doy_sonde, assume_unique=False)
    
    inds1       = np.nonzero( np.in1d( doy_sonde, intrsctVals, assume_unique=False ) )[0]
    inds2       = np.nonzero( np.in1d( doy_fts, intrsctVals, assume_unique=False ) )[0]

    print 'Total Number of coincident dates between FTS and sondes = ' +str(len(intrsctVals))+'\n'

    #-------------------------------------------------------
    #FINDING ONLY DELTA TIMES
    #-------------------------------------------------------
    indsDT = []

    for d, da in enumerate(DatesSonde):
    
        deltaT    = dtSonde[d] -   dtFTS
        deltaTmin =  np.asarray( [k.total_seconds() for k in deltaT])/60.

        indstmp = np.array([x for x, dff in enumerate(deltaTmin) if abs(dff) <= CoiTime])

        if len(indstmp) >= 1:
            indsDT.extend(indstmp)

    indsDT = np.asarray(indsDT)

    #-------------------------------------------------------
    #CREATING NEW DIR LIST FOR COINCIDENT DATES
    #-------------------------------------------------------

    dbOutput = dict((key, [val[i] for i in indsDT]) for (key, val) in dbInputs.iteritems()) 


    order = {DBfieldNames[i]:i for i in range(len(DBfieldNames))}

    with open(OutspcdbFile, 'w') as fopen:
        
        strformat = ['{0:<15}'] + [' {'+str(i)+':<12}' for i in range(1,len(DBfieldNames))]
        strformat = ''.join(strformat).lstrip().rstrip() + '\n'
        
        fopen.write(strformat.format(*[k for k in sorted(dbOutput,key=order.get)]))
        for row in zip(*[dbOutput[k] for k in sorted(dbOutput, key=order.get)]):
            fopen.write(strformat.format(*row))


if __name__ == "__main__":
    main()