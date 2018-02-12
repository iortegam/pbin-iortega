#!/usr/bin/python
##! /usr/local/python-2.7/bin/python
#----------------------------------------------------------------------------------------
# Name:
#      sondeWaterPrf.py
#
# Purpose:
#      This program creates water profiles from sonde Profiles into daily folders that can be used to retrieve gases
#
#
# Input files:
#       1)
#
# Output files:
#       1)
#
#
# Notes:
#       1) Initially created using FPH profiles at FL0
#
# References:
#
#----------------------------------------------------------------------------------------


                        #-------------------------#
                        # Import Standard modules #
                        #-------------------------#
import sys
import os
import datetime                                            as dt
import dataOutClass                                        as dc
import numpy                                               as np
from scipy.interpolate import InterpolatedUnivariateSpline as intrpUniSpl
import matplotlib.pyplot                                   as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.interpolate import InterpolatedUnivariateSpline as intrpUniSpl
from scipy import interpolate
import glob

import classSondeFTS as rs


                        #-------------------------------------#
                        # Define helper functions and classes #
                        #-------------------------------------#

def ckDir(dirName,exitFlg=False):
    ''' '''
    if not os.path.exists( dirName ):
        print 'Input Directory %s does not exist' % (dirName)
        if exitFlg: sys.exit()
        return False
    else:
        return True

def ckFile(fName,exitFlg=False):
    '''Check if a file exists'''
    if not os.path.isfile(fName):
        print 'File %s does not exist' % (fName)
        if exitFlg: sys.exit()
        return False
    else:
        return True

def segmnt(seq,n):
    '''Yeilds successive n-sized segments from seq'''
    for i in xrange(0,len(seq),n): yield seq[i:i+n]

def findCls(dataArray, val):
    ''' Returns the indice and closest value in dataArray to val'''
    return np.argmin(abs(val-dataArray))


                            #----------------------------#
                            #                            #
                            #        --- Main---         #
                            #                            #
                            #----------------------------#

def main():

    #----------------
    # Initializations
    #----------------
    loc          = 'fl0'
    sondeDir     = '/data1/ancillary_data/sondes/'+loc.lower()+'/'
    verW         = 'v77'                 # version 77 is used for individual retrievals at specific times

    pltFlg       = True

    #-------------------
    # Interpolation Oder
    #-------------------
    intrpOrder = 1      
    nSkip      = 1            # Number of points to skip when merging WACCM and NCEP profiles    

    #-----------------------
    # Date Range of interest
    #-----------------------
    iyear          = 2010
    imnth          = 1
    iday           = 1
    fyear          = 2016
    fmnth          = 12
    fday           = 31

    #----------------------
    # Output Data Directory
    #----------------------
    dataDir  = '/data1/'+loc.lower()+'/'

    #-------------------------
    # WACCM monthly means file
    #-------------------------
    WACCMfile = '/data/Campaign/'+loc.upper()+'/waccm/WACCM_pTW-meanV6.'+loc.upper()

    #--------------------------------------------------------------------
    # Read WACCM monthly mean data. WACCM monthly mean file is ascending,
    # adjust so that it is descending. Also units are in km.
    #--------------------------------------------------------------------
    with open(WACCMfile, 'r') as fopen:
        lines = fopen.readlines()
        
    nlyrs = int(lines[0].strip().split()[0])
    s_ind  = 3
    Z      = np.flipud( np.array( [ float(row.strip().split()[0]) for row in lines[s_ind:nlyrs+s_ind] ] ) )
    waccmT = np.flipud( np.array( [ [float(x) for x in line.strip().split()[1:]] for line in lines[s_ind:nlyrs+s_ind] ] ) )
    s_ind  = 3 + nlyrs + 2
    waccmP = np.flipud( np.array( [ [float(x) for x in line.strip().split()[1:]] for line in lines[s_ind:nlyrs+s_ind] ] ) )
    s_ind  = 3 + nlyrs + 2 + nlyrs + 2
    waccmW = np.flipud( np.array( [ [float(x) for x in line.strip().split()[1:]] for line in lines[s_ind:nlyrs+s_ind] ] ) )


    #----------------------------
    # File and directory checking
    #----------------------------
    ckDir(sondeDir,exitFlg=True)
    ckDir(dataDir,exitFlg=True)

    s  = rs.SondeClass(sondeDir, loc, iyear=iyear,imnth=imnth,iday=iday,fyear=fyear,fmnth=fmnth,fday=fday, fleoutFlg=False)

    s.readfilessonde()
    s.sondePrf()

    sonde = s.sonde

    sondedt            = np.asarray(sonde['dt'])
    sondedate          = np.asarray(sonde['date'])
    sondealt           = np.asarray(sonde['Alt_km_a'])
    sondeH2OPrf_a      = np.asarray(sonde['H2Omr_ppmv_a'])*1e-6
    sondeH2OPrfsd_a    = np.asarray(sonde['H2Osd_ppmv_a'])
    sondealt_a         = np.asarray(sonde['Alt_km_a'])
    sondeairmass_a     = np.asarray(sonde['airmass_a'])


    #---------------------------
    # Iterate through sonde Prfs
    #---------------------------
    for i,sngDay in enumerate(sondedt):

        #----------------------------
        # Get corresponding directory
        #----------------------------

        baseDir = dataDir + '{0:04d}{1:02d}{2:02d}/'.format(sngDay.year,sngDay.month,sngDay.day)        

        if not os.path.exists( baseDir ): continue

        #--------------------------------------------------
        # Find month index for monthly WACCM water profiles
        #--------------------------------------------------
        mnthInd = sngDay.month - 1     # -1 because January is in the 0th column

        #---------------------------------
        #
        #---------------------------------
        Z_day = sondealt[i]
        Q_day = sondeH2OPrf_a[i]

        sondetop = Z_day[-1]
        topInd   = np.argmin( abs(Z - sondetop) )   # Where top of NCEP reanalysis fits in WACCM grid height
        
        Zin  = np.concatenate( ( Z[0:(topInd-nSkip)], np.flipud(Z_day)) )
        
        SHin = np.concatenate( ( waccmW[0:(topInd-nSkip),mnthInd], np.flipud(Q_day)) )
        
        #-----------------------------------
        # Interpolate retrieved profile onto
        # sfit input grid
        #-----------------------------------
        H2Oout      = interpolate.interp1d(Zin,SHin , axis=0, fill_value=Q_day[0], bounds_error=False, kind='linear')(Z)
        H2Oout2     = np.flipud( intrpUniSpl( np.flipud(Zin), np.flipud(SHin), k=1 )( np.flipud(Z) ) )

        #------------------------
        # remove v4 for 6-hourly files
        #------------------------
        waterFiles_test = glob.glob(baseDir+'w-120*'+verW)
        waterFiles_test = [i for i in waterFiles_test if len(os.path.basename(i)) > 10]

        if len(waterFiles_test) >= 1:
            for f in waterFiles_test:
                os.remove(f)


        #---------------------
        # Write out water file
        #---------------------
        tstamp = '{0:04d}{1:02d}{2:02d}.{3:02d}{4:02d}{5:02d}'.format(sngDay.year,sngDay.month,sngDay.day,sngDay.hour,sngDay.minute,sngDay.second)

        with open(baseDir+'w-120.'+tstamp+'.'+verW,'w') as fopen:
            fopen.write('    1     H2O from Sonde  \n')

            for row in segmnt(H2Oout,5):
                strformat = ','.join('{:>12.4E}' for i in row) + ', \n'
                fopen.write(strformat.format(*row))

        #--------------------
        # Create plots to pdf
        #--------------------
        if pltFlg:
            pdfsav = PdfPages(baseDir+'Sonde_'+tstamp+'_WaterProfile.pdf')

            fig1,ax1 = plt.subplots()
            ax1.plot(H2Oout,Z,'rx-', label='Interpolated Water Profile')
            ax1.plot(Q_day,Z_day,'bx-',label='Original Sonde Water Profile')
            ax1.grid(True,which='both')
            ax1.legend(prop={'size':9})
            ax1.set_ylabel('Altitude [km]')
            ax1.set_xlabel('VMR')
            ax1.tick_params(axis='x',which='both',labelsize=8)
            ax1.set_ylim((Z[-1],60))
            #ax1.set_xlim((0,np.max((waccmW[-1,mnthInd],dayShum[-1]))))
            ax1.set_title(sngDay)

            pdfsav.savefig(fig1,dpi=250)

            pdfsav.close()

        print 'Finished processing folder: {}'.format(sngDay)


if __name__ == "__main__":
    main()