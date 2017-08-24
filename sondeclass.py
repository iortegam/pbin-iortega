import datetime as dt
import time
import math
import sys
import numpy as np
import os
import csv
import itertools
from collections import OrderedDict
import os
from os import listdir
from os.path import isfile, join
import re
from scipy.integrate import simps
from scipy import linspace, polyval, polyfit, sqrt, stats, randn
import matplotlib.dates as md

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import FormatStrFormatter, MultipleLocator,AutoMinorLocator,ScalarFormatter
from matplotlib.backends.backend_pdf import PdfPages #to save multiple pages in 1 pdf...
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as mplcm
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
import matplotlib.gridspec as gridspec
import matplotlib.colorbar as colorbar
import glob
import myfunctions as mf

import dataOutClass as dc
from scipy import interpolate
from scipy.interpolate import InterpolatedUnivariateSpline as intrpUniSpl
from matplotlib.dates import DateFormatter, MonthLocator, YearLocator, DayLocator, WeekdayLocator, MONDAY
import pylab as P

from math import acos, asin, atan2, cos, hypot
from math import degrees, pi as PI, radians, sin, tan




#------------------------------------------------------------------------------------------------------------------------------    
class _DateRange(object):
    '''
    This is an extension of the datetime module.
    Adds functionality to create a list of days.
    '''
    def __init__(self,iyear,imnth,iday,fyear,fmnth,fday, incr=1):
        self.i_date   = dt.date(iyear,imnth,iday)                                                     # Initial Day
        self.f_date   = dt.date(fyear,fmnth,fday)                                                     # Final Day
        self.dateList =[self.i_date + dt.timedelta(days=i) for i in range(0, self.numDays(), incr)]   # Incremental day list from initial to final day

    def numDays(self):
        '''Counts the number of days between start date and end date'''
        return (self.f_date + dt.timedelta(days=1) - self.i_date).days

    def inRange(self,crntyear,crntmonth,crntday):
        '''Determines if a specified date is within the date ranged initialized'''
        crntdate = dt.date(crntyear,crntmonth,crntday)
        if self.i_date <= crntdate <= self.f_date:
            return True
        else:
            return False

    def nearestDate(self, year, month, day=1, daysList=False):
        ''' Finds the nearest date from a list of days based on a given year, month, and day'''
        testDate = dt.date(year, month, day)
        if not daysList:
            daysList = self.dateList
        return min( daysList, key=lambda x:abs(x-testDate) )

    def yearList(self):
        ''' Gives a list of unique years within DateRange '''
        years = [ singDate.year for singDate in self.dateList]               # Find years for all date entries
        years = list(set(years))                                             # Determine all unique years
        years.sort()
        return years

    def daysInYear(self,year):
        ''' Returns an ordered list of days from DateRange within a specified year '''
        if isinstance(year,int):
            newyears = [inYear for inYear in self.dateList if inYear.year == year]
            return newyears
        else:
            print 'Error!! Year must be type int for daysInYear'
            return False


#------------------------------------------------------------------------------------------------------------------------------                     
class ReadOutputsSonde(_DateRange):

    def __init__(self, dataDirSonde,loc, iyear=False,imnth=False,iday=False,fyear=False,fmnth=False,fday=False,incr=1):       
        
        self.dirLstsonde           = []   #DIRECTORY WITH THE SONDE FILES

        self.loc = loc.lower()
        
        #---------------------------------
        # Test if date range given. If so 
        # create a list of directories to 
        # read data from
        #---------------------------------
        if all([iyear,imnth,iday,fyear,fmnth,fday]):
            
            # check if '/' is included at end of path
            if not( dataDirSonde.endswith('/') ):
                dataDirSonde = dataDirSonde + '/'  

            self.dataDirSonde          = dataDirSonde         
                
            # Check if directory exits
            dc.ckDir(dataDirSonde,exitFlg=True)
            
            _DateRange.__init__(self, iyear, imnth, iday, fyear, fmnth, fday, incr=1)

            #--------------------------------------------
            # Walk through first level of directories and
            # collect directory names for processing
            #--------------------------------------------

            for dd in self.dateList:
                # Find year month and day strings
                yrstr   = "{0:02d}".format(dd.year)
                mnthstr = "{0:02d}".format(dd.month)
                daystr  = "{0:02d}".format(dd.day)

                if loc.lower() == 'mlo':   filename = 'HIH_H2O_'+yrstr+mnthstr+daystr+'*.txt'
                elif loc.lower() == 'fl0': 
                    filename = 'BLD_H2O_'+yrstr+mnthstr+daystr+'*.txt'
                    filenamefleout =  '*'+yrstr+'_'+mnthstr+'_'+daystr+'*'

                
                sondefilefleout = glob.glob( dataDirSonde + filenamefleout )

                sondefile = glob.glob( dataDirSonde + filename )

                if not sondefile: continue

                for k in sondefile:
                    self.dirLstsonde.append(k)

            self.dirLstsonde.sort() 


    def readfilessonde(self):
        ''' Reads in txt files from sondes '''
        self.sonde  = {}
        self.fleout = {}

        #------------------------------------
        # Loop through collected directories
        # self.dirLst has already been sorted
        #------------------------------------
        for indMain,sngDir in enumerate(self.dirLstsonde):

            #-----------------------------------------
            #Reading FPH Sonde
            #-----------------------------------------
            if not os.path.isfile(sngDir): continue

            with open(sngDir ,'r') as fopen: lines = fopen.readlines()
            info = [ row.strip().split() for row in lines if 'Water Vapor Flight Date' in row]
            info2 = [ row.strip().split() for row in lines if 'Number of header lines' in row] 
            info3 = [ row.strip().split() for row in lines if 'Number of variables' in row]

            nheaders = int(info2[-1][-1])
            nvariable = int(info3[-1][-1])

            lines[:] = [ row.strip().split() for row in lines[nheaders:-1] ] 
           
            #lines[:] = [ row.strip().split() for row in lines if len(row.strip().split()) == 13 and not 'to' in row and not 'Level' in row and not 'Traditionally' in row and not 'Number' in row]              

            npoints = len(lines)

            time = [row[7].split()[0] for row in info]
            
            strfile = [x for x in sngDir.split('/')]
            filename = strfile[-1]
            yyyy = filename[8:12]
            mm   = filename[12:14]
            dd   = filename[14:16]
           
            time2 = time[0].split(':')
            
            hh = int(time2[0])
            mi = int(time2[1])
            se = int(time2[2])  

            testQ = [row[7] for row in lines if float(row[7]) > 0.0]
            testQ = np.asarray(testQ)

            if len(testQ) > 2:

                self.sonde.setdefault('Level_Number',[]).append([row[0] for row in lines if float(row[7]) > 0.0])  #and filtering when h2o is negative 
                self.sonde.setdefault('Alt_km',[]).append([row[1] for row in lines if float(row[7]) > 0.0])
                self.sonde.setdefault('Press_hPa',[]).append([row[2] for row in lines if float(row[7]) > 0.0])
                self.sonde.setdefault('Temp_degC',[]).append([row[3] for row in lines if float(row[7]) > 0.0])
                self.sonde.setdefault('Theta_degK',[]).append([row[4] for row in lines if float(row[7]) > 0.0])
                self.sonde.setdefault('RH_FPH',[]).append([row[5] for row in lines if float(row[7]) > 0.0])
                self.sonde.setdefault('RH_RS',[]).append([row[6] for row in lines if float(row[7]) > 0.0])
                self.sonde.setdefault('H2Omr_ppmv',[]).append([row[7] for row in lines if float(row[7]) > 0.0])
                self.sonde.setdefault('H2Osd_ppmv',[]).append([row[8] for row in lines if float(row[7]) > 0.0])
                self.sonde.setdefault('H2Omr_orig_ppmv',[]).append([row[9] for row in lines if float(row[7]) > 0.0])
                self.sonde.setdefault('TFP_FPH_decC',[]).append([row[10] for row in lines if float(row[7]) > 0.0])
                self.sonde.setdefault('TFP_RS_degC',[]).append([row[11] for row in lines if float(row[7]) > 0.0])
                if nvariable > 13: self.sonde.setdefault('O3mr_ppbv',[]).append([row[12] for row in lines if float(row[7]) > 0.0])
                else: self.sonde.setdefault('O3mr_ppbv',[]).append([])
                self.sonde.setdefault('dt',[]).append(dt.datetime( int(yyyy), int(mm), int(dd), hh, mi, se  ))
                self.sonde.setdefault('date',[]).append(dt.date(int(yyyy),  int(mm), int(dd) ))

                if self.loc.lower() == 'fl0':

                    #-----------------------------------------
                    #Reading fleout Sonde
                    #-----------------------------------------
                    filename =  '*'+yyyy+'_'+mm+'_'+dd+'*'

                    sondefilefleout = glob.glob( self.dataDirSonde + filename )

                    if sondefilefleout:

                        with open(sondefilefleout[0] ,'r') as fname: 
                            lines = fname.readlines()
                        
                            info = [ row.strip().split() for row in lines if 'Header lines' in row]
                            info2 = [ row.strip().split() for row in lines if 'Data columns' in row] 

                            nheaders = int(info[-1][-1])
                            nvariable = int(info2[-1][-1])

                            lines[:] = [ row.strip().split(',') for row in lines[nheaders-2:-1] ] 
                            headings = lines[0]

                            headings = [f.strip() for f in headings]
                           
                            indsAlt = [i for i, s in enumerate(headings) if s == 'Alt']#np.where(headings == 'Alt')[0]
                            indsLat = [i for i, s in enumerate(headings) if s == 'GPS lat']
                            indsLon = [i for i, s in enumerate(headings) if s == 'GPS lon']

                            
                            if len(indsLat) == 1:
                                self.fleout.setdefault('Alt',[]).append([row[int(indsAlt[0])] for row in lines[2:-1] if float(row[int(indsLat[0])]) > 0.0])
                                self.fleout.setdefault('GPS lat',[]).append([row[int(indsLat[0])] for row in lines[2:-1] if float(row[int(indsLat[0])]) > 0.0])
                                self.fleout.setdefault('GPS lon',[]).append([row[int(indsLon[0])] for row in lines[2:-1] if float(row[int(indsLat[0])]) > 0.0])

                            else:

                                self.fleout.setdefault('Alt',[]).append([np.nan for row in lines[2:-1]])
                                self.fleout.setdefault('GPS lat',[]).append([np.nan for row in lines[2:-1]])
                                self.fleout.setdefault('GPS lon',[]).append([np.nan for row in lines[2:-1]])

                    else: 

                        self.fleout.setdefault('Alt',[]).append([np.nan for row in lines if float(row[7]) > 0.0])
                        self.fleout.setdefault('GPS lat',[]).append([np.nan for row in lines if float(row[7]) > 0.0])
                        self.fleout.setdefault('GPS lon',[]).append([np.nan for row in lines if float(row[7]) > 0.0])


    def sondePrf(self):
        ''' Post analysis of sonde profiles, i.e., ascent/descent profiles, concentration profiles '''


        for a, alt in enumerate(self.sonde['Alt_km']):


            #------------------------------------------------
            #First, calculate density and concentration profiles in molec/cm2
            #------------------------------------------------
            pres         = np.asarray(self.sonde['Press_hPa'][a], dtype=np.float32) 
            Temp         = np.asarray(self.sonde['Temp_degC'][a], dtype=np.float32)
            H2Omr        = np.asarray(self.sonde['H2Omr_ppmv'][a], dtype=np.float32)
            H2Osd        = np.asarray(self.sonde['H2Osd_ppmv'][a], dtype=np.float32)
            altsonde     = np.asarray(alt, dtype=np.float32)
            Temp         = Temp + 273.15
            rho          = (pres*100.00) / ((Temp)*8.314) * 6.022e23 /100. /100. /100.
            self.sonde.setdefault('rho',[]).append(rho)
            
            
            #------------------------------------------------
            #Midpoints
            #------------------------------------------------
            rho_mid      = np.absolute ((rho[1:] + rho[:-1]) / 2.0)
            altsonde_mid = np.absolute ((altsonde[1:] + altsonde[:-1]) / 2.0)
            H2Omr_mid    = np.absolute ((H2Omr[1:] + H2Omr[:-1]) / 2.0)
            H2Osd_mid    = np.absolute ((H2Osd[1:] + H2Osd[:-1]) / 2.0)

            dz           = np.absolute ((altsonde[1:] - altsonde[:-1])) *1000.*100.
          
            airmass      = rho_mid*dz
            self.sonde.setdefault('airmass',[]).append(airmass)
            CPrf_mid     = (H2Omr_mid*1e-6)*rho_mid
            CPrf_sd_mid  =  (H2Osd_mid*1e-6) *rho_mid
            
            H2OMol_mid   = (H2Omr_mid*1e-6)*airmass
            H2OMolsd_mid = (H2Osd_mid*1e-6)*airmass

            self.sonde.setdefault('H2Omr_ppmv_mid',[]).append(H2Omr_mid)
            self.sonde.setdefault('H2Osd_ppmv_mid',[]).append(H2Osd_mid)
            self.sonde.setdefault('H2OMol_mid',[]).append(H2OMol_mid)
            self.sonde.setdefault('H2OMolsd_mid',[]).append(H2OMolsd_mid)
            self.sonde.setdefault('CPrf_mid',[]).append(CPrf_mid)
            self.sonde.setdefault('CPrf_sd_mid',[]).append(CPrf_sd_mid)
            self.sonde.setdefault('Alt_km_mid',[]).append(altsonde_mid)
            
    
            #------------------------------------------------
            #Calculate ascent and descent profiles
            #----------------------------------------------
            inds = np.where(altsonde_mid == np.amax(altsonde_mid))[0]
            
            altsonde_a    = np.asarray(altsonde_mid[0:inds[0]] , dtype=np.float32)
            altsonde_d    = np.asarray(altsonde_mid[inds[0]+1:] , dtype=np.float32)           

            self.sonde.setdefault('H2Omr_ppmv_a',[]).append(H2Omr_mid[0:inds[0]])
            self.sonde.setdefault('H2Osd_ppmv_a',[]).append(H2Osd_mid[0:inds[0]])
            self.sonde.setdefault('Alt_km_a',[]).append(altsonde_a)
            self.sonde.setdefault('CPrf_a',[]).append(CPrf_mid[0:inds[0]])
            self.sonde.setdefault('rho_a',[]).append(rho_mid[0:inds[0]])
            self.sonde.setdefault('H2OMol_a',[]).append(H2OMol_mid[0:inds[0]])
            self.sonde.setdefault('H2OMolsd_a',[]).append(H2OMolsd_mid[0:inds[0]])
            self.sonde.setdefault('airmass_a',[]).append(airmass[0:inds[0]])


            self.sonde.setdefault('H2Omr_ppmv_d',[]).append(H2Omr_mid[inds[0]+1:])
            self.sonde.setdefault('H2Osd_ppmv_d',[]).append(H2Osd_mid[inds[0]+1:])
            self.sonde.setdefault('Alt_km_d',[]).append(altsonde_d)
            self.sonde.setdefault('CPrf_d',[]).append(CPrf_mid[inds[0]+1:])
            self.sonde.setdefault('rho_d',[]).append(rho_mid[inds[0]+1:])
            self.sonde.setdefault('H2OMol_d',[]).append(H2OMol_mid[inds[0]+1:])
            self.sonde.setdefault('H2OMolsd_d',[]).append(H2OMolsd_mid[inds[0]+1:])
            self.sonde.setdefault('airmass_d',[]).append(airmass[inds[0]+1:])
            
            #------------------------------------------------
            #Interpolate fleout lat/lon to FPH altitude
            #----------------------------------------------
            if self.loc.lower() == 'fl0':

                try: 
                
                    AltFleout    = np.asarray(self.fleout['Alt'][a], dtype=np.float32)
                    latFleout    = np.asarray(self.fleout['GPS lat'][a], dtype=np.float32)  
                    lonFleout    = np.asarray(self.fleout['GPS lon'][a], dtype=np.float32)        

                    latFleout_a       = interpolate.interp1d(AltFleout, latFleout, bounds_error=False)(altsonde_a)
                    latFleout_d       = interpolate.interp1d(AltFleout, latFleout, bounds_error=False)(altsonde_d)

                    self.sonde.setdefault('latFleout_a',[]).append(latFleout_a)
                    self.sonde.setdefault('latFleout_d',[]).append(latFleout_d)

                    lonFleout_a       = interpolate.interp1d(AltFleout, lonFleout, bounds_error=False)(altsonde_a)
                    lonFleout_d       = interpolate.interp1d(AltFleout, lonFleout, bounds_error=False)(altsonde_d)

                    self.sonde.setdefault('lonFleout_a',[]).append(lonFleout_a)
                    self.sonde.setdefault('lonFleout_d',[]).append(lonFleout_d)

                except Exception as errmsg:
                    print errmsg




#------------------------------------------------------------------------------------------------------------------------------        
class Pltsonde(ReadOutputsSonde, dc.ReadOutputData):

    def __init__(self, dataDir, ctlF, dataDirSonde, loc, iyear=False,imnth=False,iday=False,fyear=False,fmnth=False,fday=False,incr=1,outFname=''):\
        
        primGas = ''
        if outFname: self.pdfsav = PdfPages(outFname)
        else:        self.pdfsav = False


        #---------------
        # ReadOutputsSonde
        #---------------
        ReadOutputsSonde.__init__(self,dataDirSonde,loc,iyear,imnth,iday,fyear,fmnth,fday,incr)    

        #------------
        # ReadOutputData
        #------------
        dc.ReadOutputData.__init__(self,dataDir,primGas,ctlF,iyear,imnth,iday,fyear,fmnth,fday,incr) 

   
    def closeFig(self):
        self.pdfsav.close()

    def plt1(self,fltr=False,minSZA=0.0,maxSZA=80.0,maxRMS=10.0,minDOF=1.0,maxCHI=2.0,minTC=1.0E15,maxTC=1.0E16,dofFlg=False,rmsFlg=True,tcFlg=True,mnthFltr=[1,2,3,4,5,6,7,8,9,10,11,12],
               pcFlg=True,cnvrgFlg=True,allGas=False,sclfct=1.0,sclname='ppv',pltStats=True,szaFlg=False,errFlg=False,chiFlg=False,tcMMflg=False,mnthFltFlg=False, pCols=False, smthFlg=False, doiplt=False, loc='', latFTS=False, lonFTS=False):
        
        #-------------------------------------------------------
        #READING ALL THE SONDE FILES
        #-------------------------------------------------------
        self.readfilessonde()

        print 'Total Number of sonde observations = ' +str(len(self.dirLstsonde))

        self.datefts = []
        self.dirDateTime.sort()

        for dd in self.dirDateTime:
            self.datefts.append(dt.date( int(dd.year), int(dd.month), int(dd.day) ) )

        self.datefts = np.asarray(self.datefts)
        self.sonde['date'] = np.asarray(self.sonde['date'])


        doy_fts = dc.toYearFraction(self.datefts)
        doy_sonde = dc.toYearFraction(self.sonde['date'])

        #-------------------------------------------------------
        #FINDING COINCIDENT DATES
        #-------------------------------------------------------
        intrsctVals = np.intersect1d(doy_fts, doy_sonde, assume_unique=False)
    
        inds1       = np.nonzero( np.in1d( doy_sonde, intrsctVals, assume_unique=False ) )[0]
        inds2       = np.nonzero( np.in1d( doy_fts, intrsctVals, assume_unique=False ) )[0]

        print 'Total Number of coincident dates between FTS and sondes = ' +str(len(intrsctVals))+'\n'

        #-------------------------------------------------------
        #CREATING NEW DIR LIST FOR COINCIDENT DATES
        #-------------------------------------------------------
        self.dirLstsonde = np.array(self.dirLstsonde)
        self.dirLst      = np.array(self.dirLst)

        self.dirLstsonde = self.dirLstsonde[inds1]
        self.dirLst      = self.dirLst[inds2]

        #-------------------------------------------------------
        #READING THE SONDE FILES FOR COINCIDENT DATES
        #-------------------------------------------------------   

        self.readfilessonde()
        self.sondePrf()

        #-------------------------------------------------------
        #Test sonde
        #-------------------------------------------------------
        # sondelat    = np.asarray(self.sonde['latFleout_a'])
        # sondelon    = np.asarray(self.sonde['lonFleout_a'])

        # sondealt  = np.asarray(self.sonde['Alt_km_a'])
        # sondeprf  = np.asarray(self.sonde['H2Omr_ppmv_a'])

        # print sondealt.shape, sondeprf.shape

        # opp   = []
        # for ii, r in enumerate(sondeprf):
        #     opp.append(np.average(sondealt[ii], weights= r))

        # opp = np.asarray(opp)

        # print np.mean(opp)
        # exit()
        
        # fig, ax1       = plt.subplots()

        # for i in range(sondelat.shape[0]):
    
        #     ax1.plot(sondelon[i][:], sondelat[i][:], linewidth=0.75)
        
        # ax1.set_ylabel('Latitude')
        # ax1.set_xlabel('Longitude')        
        # ax1.grid(True,which='both')
        # ax1.set_title('Location', fontsize=12)
        # ax1.tick_params(labelsize=12)
        # ax1.set_xlim(-106, -102)
        # ax1.set_ylim(38, 42)

        # if self.pdfsav: 
        #     self.pdfsav.savefig(fig,dpi=200)
        #     plt.close(fig)
        # else:
        #     plt.show(block=False)
        #     user_input = raw_input('Press any key to exit >>> ')
        #     sys.exit()


        #-------------------------------------------------------
        #READING THE FTS DATA ONLY FOR COINCIDENT DATES
        #-------------------------------------------------------
        aprPrf       = {}
        rPrf         = {}
        localGasList = [self.PrimaryGas]
        
        #-------------------------------------------------------
        # Get profile, summary for Primary gas.... for filtering
        #-------------------------------------------------------
        if not self.readPrfFlgRet[self.PrimaryGas]: self.readprfs([self.PrimaryGas],retapFlg=1)   # Retrieved Profiles
        if self.empty: return False
        if not self.readPrfFlgApr[self.PrimaryGas]: self.readprfs([self.PrimaryGas],retapFlg=0)   # Apriori Profiles
        
        aprPrf      = np.asarray(self.aprfs[self.PrimaryGas]) * sclfct
        rPrf        = np.asarray(self.rprfs[self.PrimaryGas]) * sclfct
        rPrfMol     = np.asarray(self.rprfs[self.PrimaryGas]) * np.asarray(self.rprfs['AIRMASS'])
        aprPrffMol  = np.asarray(self.aprfs[self.PrimaryGas]) * np.asarray(self.rprfs['AIRMASS'])
        Airmass     = np.asarray(self.rprfs['AIRMASS'])

        if len(self.dirLst) > 1: 
            dates                   = self.rprfs['date']

        alt = np.asarray(self.rprfs['Z'][0,:])
        
        if not self.readsummaryFlg: self.readsummary()                                      # Summary File info
        rms     = np.asarray(self.summary[self.PrimaryGas+'_FITRMS'])
        dofs    = np.asarray(self.summary[self.PrimaryGas+'_DOFS_TRG'])
        totClmn = np.asarray(self.summary[self.PrimaryGas.upper()+'_RetColmn'])
        
        if not self.readPbpFlg: self.readPbp()                                              # Pbp file info
        sza   = self.pbp['sza']
        saa   = self.pbp['saa']
                 
        if errFlg:                                    # Error info
            
            if not all((self.readErrorFlg['totFlg'],self.readErrorFlg['sysFlg'],self.readErrorFlg['randFlg'])):
                self.readError(totFlg=True,sysFlg=True,randFlg=True,vmrFlg=False,avkFlg=True,KbFlg=False) 
            
            npnts    = np.shape(self.error['Total_Random_Error'])[0]
            nlvls    = np.shape(alt)[0]

            #-------------------------------------------------
            #Error profiles constructed with the diagonal elements of the covariances matrices
            #-------------------------------------------------
            rand_err = np.zeros((npnts,nlvls))
            sys_err  = np.zeros((npnts,nlvls))
            
            for i in range(npnts):
                rand_err[i,:] = np.diag(self.error['Total_Random_Error'][i][:,:])
                sys_err[i,:]  = np.diag(self.error['Total_Systematic_Error'][i][:,:])
                
            tot_err  = np.sqrt(rand_err + sys_err ) #* sclfct         
            rand_err = np.sqrt(rand_err) #* sclfct
            sys_err  = np.sqrt(sys_err)  #* sclfct

            rand_errvmr  = (rand_err/ np.asarray(self.rprfs['AIRMASS'][0:npnts]))*sclfct
            sys_errvmr   = (sys_err/ np.asarray(self.rprfs['AIRMASS'][0:npnts]))*sclfct
            tot_errvmr   = np.sqrt(rand_errvmr**2 + sys_errvmr**2)
            
            #-------------------------------------------------
            #Error profiles of components
            #-------------------------------------------------
            rand_cmpnts      = self.randErrDiag 
            sys_cmpnts       = self.sysErrDiag

            rand_cmpnts_vmr  = self.randErrDiag  
            sys_cmpnts_vmr   = self.sysErrDiag

            for k in sys_cmpnts_vmr:
                 sys_cmpnts_vmr[k] = (np.sqrt(sys_cmpnts_vmr[k])/ np.asarray(self.rprfs['AIRMASS'][0:npnts]))*sclfct
                 
            for k in rand_cmpnts_vmr:
                 rand_cmpnts_vmr[k] = (np.sqrt(rand_cmpnts_vmr[k])/ np.asarray(self.rprfs['AIRMASS'][0:npnts]))*sclfct

            #-------------------------------------------------
            #Error in the summary output. Get Total Errors 
            #-------------------------------------------------
            rand_errTC  = np.asarray(self.error['Total random uncertainty'])
            sys_errTC   = np.asarray(self.error['Total systematic uncertainty'])
            tot_errTC   = np.sqrt(rand_errTC**2 + sys_errTC**2)  
    
        #--------------------
        # Call to filter data
        #--------------------
        if fltr: self.fltrData(self.PrimaryGas,mxrms=maxRMS,minsza=minSZA,mxsza=maxSZA,minDOF=minDOF,maxCHI=maxCHI,minTC=minTC,maxTC=maxTC,mnthFltr=mnthFltr,
                               dofFlg=dofFlg,rmsFlg=rmsFlg,tcFlg=tcFlg,pcFlg=pcFlg,szaFlg=szaFlg,cnvrgFlg=cnvrgFlg,chiFlg=chiFlg,tcMinMaxFlg=tcMMflg,mnthFltFlg=mnthFltFlg)
        else:    self.inds = np.array([]) 
        
        if self.empty: return False
        
        #--------------------------------------------------------
        # Get profile data for other gases if allGas flag is True
        #--------------------------------------------------------
        if allGas:
            nonPgas = (gas for gas in self.gasList if gas != self.PrimaryGas)
            for gas in nonPgas:
                if not self.readPrfFlgRet[gas.upper()]: self.readprfs([gas.upper()],retapFlg=1)   # Retrieved Profiles
                if not self.readPrfFlgApr[gas.upper()]: self.readprfs([gas.upper()],retapFlg=0)   # Apriori Profiles
                aprPrf[gas.upper()] = np.asarray(self.aprfs[gas.upper()]) * sclfct
                rPrf[gas.upper()]   = np.asarray(self.rprfs[gas.upper()]) * sclfct
                localGasList.append(gas)
          
        #----------------------------
        # Remove inds based on filter
        #----------------------------
        nfltr   = len(self.inds)
        rms     = np.delete(rms, self.inds)
        ntot    = len(rms)
        sza     = np.delete(sza, self.inds)
        saa     = np.delete(saa, self.inds)
        if len(self.dirLst) > 1:
            dates   = np.delete(dates, self.inds)
            datefts = []
            for dd in dates:
                datefts.append(dt.date( int(dd.year), int(dd.month), int(dd.day) ) )
            datefts = np.asarray(datefts)
        dofs       = np.delete(dofs,self.inds)
        totClmn    = np.delete(totClmn,self.inds)
        rPrfMol    = np.delete(rPrfMol,self.inds,axis=0)
        aprPrffMol = np.delete(aprPrffMol,self.inds,axis=0)
        Airmass    = np.delete(Airmass,self.inds,axis=0)
        rPrf       = np.delete(rPrf,self.inds,axis=0)
        aprPrf     = np.delete(aprPrf,self.inds,axis=0)

            
        if errFlg:
            rand_err    = np.delete(rand_err,self.inds,axis=0)
            sys_err     = np.delete(sys_err,self.inds,axis=0)  
            tot_err     = np.delete(tot_err,self.inds,axis=0)

            rand_errvmr = np.delete(rand_errvmr,self.inds,axis=0)
            sys_errvmr  = np.delete(sys_errvmr,self.inds,axis=0)  
            tot_errvmr  = np.delete(tot_errvmr,self.inds,axis=0)

            rand_errTC  = np.delete(rand_errTC,self.inds)
            sys_errTC   = np.delete(sys_errTC,self.inds)
            tot_errTC   = np.delete(tot_errTC,self.inds)

            
            for k in sys_cmpnts:
                sys_cmpnts[k]      = np.delete(sys_cmpnts[k],self.inds,axis=0)
                sys_cmpnts_vmr[k]  = np.delete(sys_cmpnts_vmr[k],self.inds,axis=0)
                
            for k in rand_cmpnts:
                rand_cmpnts[k]     = np.delete(rand_cmpnts[k],self.inds,axis=0)
                rand_cmpnts_vmr[k] = np.delete(rand_cmpnts_vmr[k],self.inds,axis=0)

            #---------------------
            # Get averaging kernel
            #---------------------   
            avkSCF      = np.delete(np.asarray(self.error['AVK_scale_factor']),self.inds,axis=0)
            avkVMR      = np.delete(np.asarray(self.error['AVK_vmr']), self.inds,axis=0)   
            avkSCFav    = np.mean(avkSCF,axis=0)    
            avkVMRav    = np.mean(avkVMR,axis=0)

        #--------------------------------------------------------
        # Get location of FTS based on altitude profile weighted from sone, azimuth (bearing) and SZA
        #--------------------------------------------------------
        altCG  = 3.8
        sza    = np.asarray(sza)
        saa    = np.asarray(saa)
        szarad = [radians(s) for s in sza]

        distFTS = [altCG*tan(s) for s in   szarad]

        latFTS2  = []
        lonFTS2  = []

        for i, a in enumerate(saa):
            lat2, lon2 = mf.destination(latFTS, lonFTS, distFTS[i], a, radius=6371.008771415)
            latFTS2.append(lat2)
            lonFTS2.append(lon2)

        latFTS2 = np.asarray(latFTS2)
        lonFTS2 = np.asarray(lonFTS2)
                
        #----------------------------
        # Determine if multiple years
        #----------------------------
        if len(self.dirLst) > 1:
            years = [ singDate.year for singDate in dates]      # Find years for all date entries
            if len(list(set(years))) > 1: yrsFlg = True         # Determine all unique years
            else:                         yrsFlg = False

        #---------------------------------
        # Re - Defining variables (FTIR) using Max/Min altitude
        #---------------------------------
        maxalt      = 20. #km
        indsalt     = np.where(alt<= maxalt)[0]
        alt         = alt[indsalt]
        rPrf        = rPrf[:, indsalt]
        aprPrf      = aprPrf[:, indsalt]
        tot_err     = tot_err[:, indsalt]
        rand_errvmr = rand_errvmr[:, indsalt]
        sys_errvmr  = sys_errvmr[:, indsalt]
        rPrfMol     = rPrfMol[:, indsalt]
        tot_errvmr  = tot_errvmr[:, indsalt]
        aprPrffMol  = aprPrffMol[:, indsalt]
        avkSCFav    = avkSCFav[indsalt[0]:, indsalt[0]:]
        avkVMRav    = avkVMRav[indsalt[0]:, indsalt[0]:]
        avkSCF      = avkSCF[:, indsalt[0]:, indsalt[0]: ]
        avkVMR      = avkVMR[:, indsalt[0]:, indsalt[0]: ]
        Airmass     = Airmass[:, indsalt]


        #---------------------------
        # Define Delta time,  plot Profiles, and calculate columns/VRM(weighted) for both FTS and sonde
        #---------------------------

        prfmean_day            = {}
        prferr_day             = {}
        prfSTD_day             = {}
        aPrf_day               = {}
        Airmassmean            = {}
        NobsFTS                = {}

        latmean                = {}
        lonmean                = {}

        avkSCFday              = {}
        avkSCFVMRday           = {}
        rPrfMolmean_day        = {}
        aPrfMol_day            = {}
        prferrMol_day          = {}

        sondePrf               = {}
        sondePrfsd             = {}
        sondeairmass           = {}
        sondePrfMol            = {}
        sondePrfMolsd          = {}
        sondePrflat            = {}
        sondePrflon            = {}

        Prfsonde_interp        = {}
        Prfsonde_sd_interp     = {}
        PrfsondeMol_interp     = {}
        PrfsondeMol_sd_interp  = {}
        sondeairmass_a_interp  = {}
        Prfsondelat_interp     = {}
        Prfsondelon_interp     = {}
        sondealt               = {}

        PrfDiff                = {}
        PrfDiffApr             = {}
        PrfDiffRel             = {}
        PrfDiffAprRel          = {}
        PrfDist                = {}

        FTSvmrP                = {}
        FTSsumP                = {}
        FTSvmrsdP              = {}
        FTSsumsdP              = {}
        FTSdates               = {}

        FTSAprvmrP             = {}
        FTSAprsumP             = {}

        sondevmrP              = {}
        sondevmrsdP            = {}
        sondesumP              = {}
        sondesumsdP            = {}
        sondedates             = {}

        sondelatP              = {}
        sondelonP              = {}
        distP                  = {}


        #------------------
        diffT      = [3, 10, 15, 30, 60, 90, 120, 180, 240]
        dfValue    = 120

        #---------------------------------
        #
        #---------------------------------

        for df in diffT:

            for d, da in enumerate(self.sonde['date']):

                #---------------------------------
                #Index for same days - FTIR
                #---------------------------------
                dtsonde   = self.sonde['dt'][d]
                deltaT    = dtsonde -   dates 
                deltaTmin =  np.asarray( [k.total_seconds() for k in deltaT])/60.
                
                inds = np.array([x for x, dff in enumerate(deltaTmin) if abs(dff) <= df])
                #inds = np.array([x for x, df in enumerate(datefts) if df == da])

                if len(inds) >= 1:
                    
                    #prfmean_day   = np.mean(rPrf[inds],axis=0)
                    ##WEIGTHED MEAN
                    prfmean_day.setdefault(df, []).append(np.average(rPrf[inds],axis=0, weights=tot_errvmr[inds]) ) # ==> prfmean_day   = np.sum(rPrf[inds]/(tot_errvmr[inds]**2), axis=0 )/ np.sum(1/(tot_err[inds]**2), axis=0); http://www.physics.umd.edu/courses/Phys261/F06/ErrorPropagation.pdf
                    #prfmean_day.setdefault(df, []).append(np.mean(rPrf[inds],axis=0 ))
                    ## uncertainty in the weighted mean can be evaluated using standard error propagation  
                    ##WEIGTHED UNCERTAINTY (DOUBLE CHECK THIS)
                    #prferr_day.setdefault(df, []).append(np.sqrt (1.0/np.sum(1.0/tot_errvmr[inds]**2, axis=0)))
                    ##UNCERTAINTY USING standard error propagation
                    prferr_day.setdefault(df, []).append(np.sqrt (np.sum(tot_errvmr[inds]**2, axis=0))) 

                    prfSTD_day.setdefault(df, []).append(np.std(rPrf[inds], axis=0))
                    
                    aPrf_day.setdefault(df, []).append(np.mean(aprPrf[inds], axis=0))

                    Airmassmean.setdefault(df, []).append(np.mean(Airmass[inds],axis=0))
                    avkSCFday.setdefault(df, []).append(np.mean(avkSCF[inds],axis=0))
                    avkSCFVMRday.setdefault(df, []).append(np.mean(avkVMR[inds],axis=0))
                    
                    #rPrfMolmean_day   = np.mean(rPrfMol[inds],axis=0)
                    rPrfMolmean_day.setdefault(df, []).append(np.average(rPrfMol[inds],axis=0, weights=tot_err[inds]))
                    prferrMol_day.setdefault(df, []).append(np.sqrt (1.0/np.sum(1.0/tot_err[inds]**2, axis=0)))

                    aPrfMol_day.setdefault(df, []).append(np.mean(aprPrffMol[inds],axis=0))

                    NobsFTS.setdefault(df, []).append(len(inds))

                    latmean.setdefault(df, []).append(np.mean(latFTS2[inds]))
                    lonmean.setdefault(df, []).append(np.mean(lonFTS2[inds]))

                    sondeH2OPrf_a      = np.asarray(self.sonde['H2Omr_ppmv_a'][d], dtype=np.float32)
                    sondeH2OPrfsd_a    = np.asarray(self.sonde['H2Osd_ppmv_a'][d], dtype=np.float32)
                    sondeH2OPrfMol_a   = np.asarray(self.sonde['H2OMol_a'][d], dtype=np.float32)
                    sondeH2OPrfMolsd_a = np.asarray(self.sonde['H2OMolsd_a'][d], dtype=np.float32)

                    sondealt_a          = np.asarray(self.sonde['Alt_km_a'][d], dtype=np.float32)
                    sondeairmass_a      = np.asarray(self.sonde['airmass_a'][d], dtype=np.float32)

                    sondelat_a            = np.asarray(self.sonde['latFleout_a'][d])
                    sondelon_a            = np.asarray(self.sonde['lonFleout_a'][d])

                else: continue
                
                #ind1 = mf.nearestind(pcol[0], sondealt_a)
                #ind2 = mf.nearestind(pcol[1], sondealt_a) 

                sondePrf.setdefault(df, []).append(sondeH2OPrf_a)
                sondePrfsd.setdefault(df, []).append(sondeH2OPrfsd_a)
                sondeairmass.setdefault(df, []).append(sondeairmass_a)
                sondePrfMol.setdefault(df, []).append(sondeH2OPrfMol_a)
                sondePrfMolsd.setdefault(df, []).append(sondeH2OPrfMolsd_a)

                sondePrflat.setdefault(df, []).append(sondelat_a)
                sondePrflon.setdefault(df, []).append(sondelon_a)

                sondealt_a  = map(lambda x: float(x),sondealt_a)
                sondealt_a  = np.asarray(sondealt_a)

                sondedates.setdefault(df, []).append(da)
                sondealt.setdefault(df, []).append(sondealt_a)

                day = np.asarray([dt.date(d.year,d.month,d.day) for d in dates[inds]])
                uniqueDay = list(set(day))          # Find a list of unique days   
                FTSdates.setdefault(df, []).append(uniqueDay)


                Prfsonde_interp_i           = interpolate.interp1d(sondealt_a, sondeH2OPrf_a, axis=0, fill_value=sondeH2OPrf_a[0], bounds_error=False)(alt)
                #Prfsonde_interp_i           = interpolate.interp1d(sondealt_a, sondeH2OPrf_a, axis=0, fill_value='extrapolate', bounds_error=False, kind='linear')(alt)
                Prfsonde_sd_interp_i        = interpolate.interp1d(sondealt_a, sondeH2OPrfsd_a, axis=0, fill_value='extrapolate', bounds_error=False, kind='linear')(alt)
                PrfsondeMol_interp_i        = interpolate.interp1d(sondealt_a, sondeH2OPrfMol_a, axis=0, fill_value='extrapolate', bounds_error=False, kind='linear')(alt)
                PrfsondeMol_sd_interp_i     = interpolate.interp1d(sondealt_a, sondeH2OPrfMolsd_a, axis=0, fill_value='extrapolate', bounds_error=False, kind='linear')(alt)
                sondeairmass_a_interp.setdefault(df, []).append(interpolate.interp1d(sondealt_a, sondeairmass_a, axis=0, fill_value='extrapolate', bounds_error=False, kind='linear')(alt))
                
                
                sondePrflon_interp_i            = interpolate.interp1d(sondealt_a, sondelon_a, axis=0,  bounds_error=False, kind='linear')(alt)
                sondePrflat_interp_i            = interpolate.interp1d(sondealt_a, sondelat_a, axis=0,  bounds_error=False, kind='linear')(alt)

                # print sondelon_a
                # print sondePrflon_interp_i
                # print da
                # user_input = raw_input('Press any key to exit >>> ')

                Prfsondelon_interp.setdefault(df, []).append(sondePrflon_interp_i)
                Prfsondelat_interp.setdefault(df, []).append(sondePrflat_interp_i)


                if smthFlg:
                    #---------------------------------
                    # Smoothing using FTIR AK and apriori
                    #---------------------------------
                    Prfsonde_interp.setdefault(df, []).append(np.mean(aprPrf[inds], axis=0)/1e3 + np.dot(np.mean(avkVMR[inds],axis=0), (Prfsonde_interp_i -  np.mean(aprPrf[inds], axis=0)/1e3)))
                    Prfsonde_sd_interp.setdefault(df, []).append(np.mean(aprPrf[inds], axis=0)/1e3 + np.dot(np.mean(avkVMR[inds],axis=0), (Prfsonde_sd_interp_i -  np.mean(aprPrf[inds], axis=0)/1e3)))

                    PrfsondeMol_interp.setdefault(df, []).append(np.mean(aprPrffMol[inds],axis=0)/1e3 + np.dot(np.mean(avkSCF[inds],axis=0), (PrfsondeMol_interp_i -  np.mean(aprPrffMol[inds],axis=0)/1e3)))
                    PrfsondeMol_sd_interp.setdefault(df, []).append(np.mean(aprPrffMol[inds],axis=0)/1e3 + np.dot(np.mean(avkSCF[inds],axis=0), (PrfsondeMol_sd_interp_i -  np.mean(aprPrffMol[inds],axis=0)/1e3)))

                else:
                    Prfsonde_interp.setdefault(df, []).append(Prfsonde_interp_i)
                    Prfsonde_sd_interp.setdefault(df, []).append(Prfsonde_sd_interp_i)
                    PrfsondeMol_interp.setdefault(df, []).append(PrfsondeMol_interp_i)
                    PrfsondeMol_sd_interp.setdefault(df, []).append(PrfsondeMol_sd_interp_i)
                    #sondeairmass_a_interp.setdefault(df, []).append(interpolate.interp1d(sondealt_a, sondeairmass_a, axis=0, fill_value='extrapolate', bounds_error=False, kind='nearest')(alt))
                    

            prfmean_day[df]           = np.asarray(prfmean_day[df])
            prferr_day[df]            = np.asarray(prferr_day[df])
            prfSTD_day[df]            = np.asarray(prfSTD_day[df])
            aPrf_day[df]              = np.asarray(aPrf_day[df])
            Airmassmean[df]           = np.asarray(Airmassmean[df])
            avkSCFday[df]             = np.asarray(avkSCFday[df])
            rPrfMolmean_day[df]       = np.asarray(rPrfMolmean_day[df])
            prferrMol_day[df]         = np.asarray(prferrMol_day[df])
            aPrfMol_day[df]           = np.asarray(aPrfMol_day[df])
            FTSdates[df]              = np.asarray(FTSdates[df])
            NobsFTS[df]               = np.asarray(NobsFTS[df])
            latmean[df]               = np.asarray(latmean[df])
            lonmean[df]               = np.asarray(lonmean[df])

            Prfsonde_interp[df]       = np.asarray(Prfsonde_interp[df])
            Prfsonde_sd_interp[df]    = np.asarray(Prfsonde_sd_interp[df])
            PrfsondeMol_interp[df]    = np.asarray(PrfsondeMol_interp[df])
            PrfsondeMol_sd_interp[df] = np.asarray(PrfsondeMol_sd_interp[df])
            sondeairmass_a_interp[df] = np.asarray(sondeairmass_a_interp[df])
            sondedates[df]            = np.asarray(sondedates[df])
            sondealt[df]              = np.asarray(sondealt[df])

            Prfsondelon_interp[df]    = np.asarray(Prfsondelon_interp[df])
            Prfsondelat_interp[df]    = np.asarray(Prfsondelat_interp[df])

            PrfDiff[df]               = np.asarray(prfmean_day[df] - Prfsonde_interp[df])
            PrfDiffApr[df]            = np.asarray(aPrf_day[df] - Prfsonde_interp[df])

            PrfDiffRel[df]            = np.true_divide(PrfDiff[df], Prfsonde_interp[df])*100.
            PrfDiffAprRel[df]         = np.true_divide(PrfDiffApr[df], Prfsonde_interp[df])*100.

            PrfDist[df]  = [ [mf.haversine(np.nanmean(lonmean[df]), np.nanmean(latmean[df]), lo2, la2 ) for (lo2, la2) in  zip(lo,la)] for (lo, la) in zip(Prfsondelon_interp[df], Prfsondelat_interp[df]) ]
            PrfDist[df]    = np.asarray(PrfDist[df])


            #------------------------------------
            #For ssome reason at ~10k some values does not make sense
            #------------------------------------
            #PrfDist[df][PrfDist[df] > 500.] = np.nan
             
      
            sondeairmass_a_i = Airmassmean[df]

            #------------------------------------
            #
            #------------------------------------
            # for pn, pcol in enumerate(pCols):

            #     inds = np.where( (alt >= pcol[0]) & (alt <=pcol[1])  )[0]

            #     sondevmrP[str(df)+'_'+str(pn)]    = np.ma.average(np.ma.masked_invalid(Prfsonde_interp[df][:,inds]), axis=1,  weights=sondeairmass_a_i[:,inds])
            #     sondevmrsdP[str(df)+'_'+str(pn)]  = np.ma.average(np.ma.masked_invalid(Prfsonde_sd_interp[df][:, inds]), axis=1, weights=sondeairmass_a_i[:,inds])
            #     sondesumP[str(df)+'_'+str(pn)]    = np.nansum(np.ma.masked_invalid(Prfsonde_interp[df][:,inds])*np.ma.masked_invalid(sondeairmass_a_i[:,inds])*1e-6, axis=1)
            #     sondesumsdP[str(df)+'_'+str(pn)]  = np.nansum(np.ma.masked_invalid(Prfsonde_sd_interp[df][:,inds])*np.ma.masked_invalid(sondeairmass_a_i[:,inds])*1e-6, axis=1)
                
            #     FTSvmrP[str(df)+'_'+str(pn)]      = np.ma.average(np.ma.masked_invalid(prfmean_day[df][:,inds]), axis=1,  weights=np.ma.masked_invalid(Airmassmean[df][:,inds]))
            #     FTSvmrsdP[str(df)+'_'+str(pn)]    = np.ma.average(np.ma.masked_invalid(prferr_day[df][:,inds]), axis=1,  weights=np.ma.masked_invalid(Airmassmean[df][:,inds]))
            #     FTSAprvmrP[str(df)+'_'+str(pn)]   = np.ma.average(np.ma.masked_invalid(aPrf_day[df][:,inds]), axis=1,  weights=np.ma.masked_invalid(Airmassmean[df][:,inds]))
            #     FTSsumP[str(df)+'_'+str(pn)]      = np.nansum(np.ma.masked_invalid(rPrfMolmean_day[df][:,inds]), axis=1)
            #     FTSsumsdP[str(df)+'_'+str(pn)]    = np.nansum(np.ma.masked_invalid(prferrMol_day[df][:,inds]), axis=1)
            #     FTSAprsumP[str(df)+'_'+str(pn)]   = np.nansum(np.ma.masked_invalid(aPrfMol_day[df][:,inds]), axis=1)

            #     sondelatP[str(df)+'_'+str(pn)]    = np.ma.mean(np.ma.masked_invalid(Prfsondelat_interp[df][:,inds]), axis=1)
            #     sondelonP[str(df)+'_'+str(pn)]    = np.ma.mean(np.ma.masked_invalid(Prfsondelon_interp[df][:,inds]), axis=1)

            #     #distP[str(df)+'_'+str(pn)]        = [mf.haversine(np.nanmean(lonmean[df]), np.mean(latmean[df]), lo, la ) for (lo, la) in  zip(sondelonP[str(df)+'_'+str(pn)], sondelatP[str(df)+'_'+str(pn)])]
            #     distP[str(df)+'_'+str(pn)]        = np.ma.mean(np.ma.masked_invalid(PrfDist[df][:,inds]), axis=1)
        #--------------------------------
        # Profiles as a function of Month
        #--------------------------------
        month          = np.array([d.month for d in self.sonde['date']])
        fig, ax1       = plt.subplots()
        cm             = plt.get_cmap('jet')
        cNorm          = colors.Normalize( vmin=1, vmax=12 )
        scalarMap      = mplcm.ScalarMappable( norm=cNorm, cmap=cm )
        
        scalarMap.set_array(month)
        
        ax1.set_color_cycle( [scalarMap.to_rgba(x) for x in month] )
        
        for i in range(len(month)):
            ax1.plot(np.asarray(self.sonde['H2Omr_ppmv_a'][i])/1e3, np.asarray(self.sonde['Alt_km_a'][i]), linewidth=0.75)
        
        ax1.set_ylabel('Altitude [km]', fontsize=14)
        ax1.set_xlabel('VMR [x10$^3$ ppm$_v$]', fontsize=15)        
        ax1.grid(True,which='both')
        ax1.set_ylim(1,12)
        ax1.set_title('CFH profiles (2010 - 2016)', fontsize=16)
        ax1.tick_params(labelsize=14)
        
        cbar = fig.colorbar(scalarMap,orientation='vertical')#, format='%2i')
        cbar.set_label('Month')

        labels = np.arange(1, 13 , 1)
        location = labels 
        cbar.set_ticks(location)
        cbar.set_ticklabels(labels)

        if self.pdfsav: 
            self.pdfsav.savefig(fig,dpi=200)
            plt.close(fig)
        else:
            plt.show(block=False)
            #user_input = raw_input('Press any key to exit >>> ')
            #sys.exit()
        
        #---------------------------------
        #Figure: Individual Prfs
        #---------------------------------
        
        nprfs   = prfmean_day[dfValue].shape[0]

        if self.pdfsav:

            for d, da in enumerate(self.sonde['date']):

                #---------------------------------
                #Index for same days - FTIR
                #---------------------------------
                dtsonde   = self.sonde['dt'][d]
                deltaT    = dtsonde -   dates 
                deltaTmin =  np.asarray( [k.total_seconds() for k in deltaT])/60.
                
                inds = np.array([x for x, dff in enumerate(deltaTmin) if abs(dff) <= dfValue])
                #inds = np.array([x for x, df in enumerate(datefts) if df == da])

                if len(inds) >= 1:

                    fig,(ax1,ax2, ax3, ax4) = plt.subplots(1, 4, figsize=(15,7), sharey=True, sharex=False)
                                  
                    for i in inds:
                        #ax1.plot(rPrf[i]/1e3,alt, label='{:%H:%M}'.format(dates[i]),  linewidth=2.0)  #
                        ax1.errorbar(rPrf[i]/1e3,alt, xerr=tot_errvmr[i]/1e3, markersize=35, linestyle='-', label='{:%H:%M}'.format(dates[i]), linewidth=2.0)
                        ax1.scatter(rPrf[i]/1e3,alt,facecolors='white', s=35, color='k')
                      
                    prfmean_i    = np.average(rPrf[inds],axis=0, weights=tot_errvmr[inds])
                    Aprfmean_i   = np.average(aprPrf[inds], axis=0)
                    prferrmean_i = np.sqrt (1.0/np.sum(1.0/tot_errvmr[inds]**2, axis=0))
                    prferrmean_i    = np.sqrt (np.sum(tot_errvmr[inds]**2, axis=0))
                    prfstd_i     = np.std(rPrf[inds], axis=0)


                    ax1.plot(Aprfmean_i/1e3, alt, '-', color='gray', label='a priori', linewidth=2.0)
                    ax1.scatter(Aprfmean_i/1e3, alt,facecolors='white', s=35, color='gray')
                    
                    ax1.plot(np.asarray(self.sonde['H2Omr_ppmv_a'][d])/1e3, np.asarray(self.sonde['Alt_km_a'][d]),color='k', label='sonde\n({:%H:%M})'.format(dtsonde), linewidth=2.0)
                    ax1.scatter(np.asarray(self.sonde['H2Omr_ppmv_a'][d])/1e3, np.asarray(self.sonde['Alt_km_a'][d]),facecolors='white', s=35, color='k')
                    ax1.legend(prop={'size':11})

                    ax2.plot(prfmean_i/1e3,alt, color='b',  linewidth=2.0, label='Mean (Retrieved)')
                    ax2.scatter(prfmean_i/1e3,alt,facecolors='white', s=35, color='b')
                    ax2.fill_betweenx(alt,prfmean_i/1e3-prferrmean_i/1e3,prfmean_i/1e3+prferrmean_i/1e3,alpha=0.25,color='blue')

                    ax2.plot(Aprfmean_i/1e3, alt, '-', color='gray', label='a priori', linewidth=2.0)
                    ax2.scatter(Aprfmean_i/1e3, alt,facecolors='white', s=35, color='gray')

                    Prfsonde_interp_i = interpolate.interp1d(np.asarray(self.sonde['Alt_km_a'][d]), np.asarray(self.sonde['H2Omr_ppmv_a'][d]), axis=0, fill_value=np.asarray(self.sonde['H2Omr_ppmv_a'][d][0]), bounds_error=False)(alt)
                    
                    #if smthFlg:
                    ax2.plot(Prfsonde_interp_i/1e3, alt, color='k',  linewidth=2.0, label='sonde\n({:%H:%M})'.format(dtsonde))
                    ax2.scatter(Prfsonde_interp_i/1e3, alt,facecolors='white', s=35, color='k')
                        
                    rd                 = (Prfsonde_interp_i - prfmean_i)
                    rda                 = (Prfsonde_interp_i - Aprfmean_i)

                    ax3.plot(rd/1e3, alt, '-', color='blue', linewidth=2.0, label='Retrieved')
                    ax3.scatter(rd/1e3, alt,facecolors='white', s=35, color='b')
                    
                    ax3.plot(rda/1e3, alt, '-', color='gray', linewidth=2.0, label='a priori')
                    ax3.scatter(rda/1e3, alt,facecolors='white', s=35, color='gray')
                    ax3.legend(prop={'size':11})
                    ax3.axvline(x=0, linestyle='--', linewidth=2, color='k')

                    rd = rd/Prfsonde_interp_i *100.
                    rda = rda/Prfsonde_interp_i *100.

                    ax4.plot(rd, alt, '-', color='blue', linewidth=2.0, label='Retrieved')
                    ax4.scatter(rd, alt,facecolors='white', s=35, color='b')
                    
                    ax4.plot(rda, alt, '-', color='gray', linewidth=2.0, label='a priori')
                    ax4.scatter(rda, alt,facecolors='white', s=35, color='gray')
                    ax4.legend(prop={'size':11})
                    ax4.axvline(x=0, linestyle='--', linewidth=2, color='k')

                    ax2.legend(prop={'size':12})

                    ax1.set_ylabel('Altitude [km]')
                    if loc.lower() == 'fl0': ax1.set_ylim(1, 15)
                    if loc.lower() == 'mlo': ax1.set_ylim(3, 15)
                    #ax1.set_ylim(1, 15)
                    ax1.set_xlim(xmin=0)

                    ax1.set_xlabel('VMR [x10$^3$ ppm$_v$]')
                    ax2.set_xlabel('VMR [x10$^3$ ppm$_v$]')
                    ax2.set_xlim(xmin=0)
                    ax3.set_xlabel('Residuals [x10$^3$ ppm$_v$]')
                    ax4.set_xlabel('Relative difference [%]')

                    ax1.grid(True,which='both')
                    ax2.grid(True,which='both')
                    ax3.grid(True,which='both')
                    ax4.grid(True,which='both')
                    #print mf.bias(prfmean_day/1e3,Prfsonde_interp)

                    ax1.tick_params(which='both',labelsize=12)
                    ax2.tick_params(which='both',labelsize=12)
                    ax3.tick_params(which='both',labelsize=12)
                    ax4.tick_params(which='both',labelsize=12)

                    plt.suptitle(da, fontsize=16)    

                    fig.subplots_adjust(bottom=0.1,top=0.94, left=0.05, right=0.95)    
                    #ax1.text(-0.1,1.1, 'Number of Obs Filtered        = '+str(nfltr),ha='left',va='center',transform=ax1.transAxes,fontsize=8)
                    #ax1.text(-0.1,1.05,'Number of Obs After Filtering = '+str(ntot), ha='left',va='center',transform=ax1.transAxes,fontsize=8)
               
                    if self.pdfsav: 
                        self.pdfsav.savefig(fig,dpi=200)
                        plt.close(fig)
                    else:
                        plt.show(block=False)
                        user_input = raw_input('Press any key to exit >>> ')
                        #sys.exit()

        #---------------------------------
        #Figure: selected profiles (No smoothed)
        #---------------------------------
        if doiplt:
            doidt = [dt.date(int(d[0:4]), int(d[4:6]), int(d[6:8])) for d in doiplt]
            doidt = np.asarray(doidt)
            
            fig0, ax0 = plt.subplots(2, 5, figsize=(15,10), sharey=True, sharex=False)
            ndoi = 0 

            for d, da in enumerate(sondedates[dfValue]):

                deltadoi = doidt - da

                if ndoi < len(doidt):
                
                    if deltadoi[ndoi] == dt.timedelta(0):

                        if ndoi < len(doidt):

                            if ndoi<=4:

                                if ndoi ==0: ax0[0, ndoi].set_ylabel('Altitude [km]', fontsize=14)

                                ax0[0, ndoi].plot(prfmean_day[dfValue][d]/1e3,alt, color='b',  linewidth=2.0, label='HR-FTIR')
                                ax0[0, ndoi].scatter(prfmean_day[dfValue][d]/1e3,alt,facecolors='white', s=35, color='b')
                                ax0[0, ndoi].fill_betweenx(alt,prfmean_day[dfValue][d]/1e3-prferr_day[dfValue][d]/1e3,prfmean_day[dfValue][d]/1e3+prferr_day[dfValue][d]/1e3,alpha=0.25,color='blue')

                                ax0[0, ndoi].plot(aPrf_day[dfValue][d]/1e3, alt, '-', color='gray', label='a priori', linewidth=2.0)
                                ax0[0, ndoi].scatter(aPrf_day[dfValue][d]/1e3, alt,facecolors='white', s=20, color='gray')
                                
                                ax0[0, ndoi].plot(sondePrf[dfValue][d]/1e3, sondealt[dfValue][d],color='k', label='CFH', linewidth=2.0)
                                ax0[0, ndoi].scatter(sondePrf[dfValue][d]/1e3, sondealt[dfValue][d],facecolors='white', s=20, color='k')

                                if loc.lower() == 'fl0': ax0[0, ndoi].set_ylim(1, 15)
                                if loc.lower() == 'mlo': ax0[0, ndoi].set_ylim(3, 15)
                                
                                ax0[0, ndoi].grid(True,which='both')
                                ax0[0, ndoi].tick_params(which='both',labelsize=14)
                                ax0[0, ndoi].set_title(da, fontsize=14)   

                                ax0[0, ndoi].set_xlim(xmin=0)

                                #ax0[0, ndoi].errorbar(sondePrf[dfValue][d]/1e3, sondealt[dfValue][d], xerr=sondePrf[dfValue][d]*0.1/1e3, color='k', markersize=0, linestyle='-', label='CFH', linewidth=2.0)
                                #ax0[0, ndoi].scatter(sondePrf[dfValue][d]/1e3, sondealt[dfValue][d],facecolors='white', s=35, color='k')

                                ax0[0, ndoi].text(0.05, 0.95, 'N = {}'.format(NobsFTS[dfValue][d]), va='center',transform=ax0[0, ndoi].transAxes,fontsize=14)
                                if ndoi == 0: ax0[0, ndoi].legend(prop={'size':11})

                            if ndoi>=5:

                                ax0[1, (ndoi-5)].plot(prfmean_day[dfValue][d]/1e3,alt, color='b',  linewidth=2.0, label='Mean (Retrieved)')
                                ax0[1, (ndoi-5)].scatter(prfmean_day[dfValue][d]/1e3,alt,facecolors='white', s=35, color='b')
                                ax0[1, (ndoi-5)].fill_betweenx(alt,prfmean_day[dfValue][d]/1e3-prferr_day[dfValue][d]/1e3,prfmean_day[dfValue][d]/1e3+prferr_day[dfValue][d]/1e3,alpha=0.25,color='blue')

                                ax0[1, (ndoi-5)].plot(aPrf_day[dfValue][d]/1e3, alt, '-', color='gray', label='a priori', linewidth=2.0)
                                ax0[1, (ndoi-5)].scatter(aPrf_day[dfValue][d]/1e3, alt,facecolors='white', s=20, color='gray')
                                
                                ax0[1, (ndoi-5)].plot(sondePrf[dfValue][d]/1e3, sondealt[dfValue][d],color='k', label='sonde\n({:%H:%M})'.format(dtsonde), linewidth=2.0)
                                ax0[1, (ndoi-5)].scatter(sondePrf[dfValue][d]/1e3, sondealt[dfValue][d],facecolors='white', s=20, color='k')

                                if ndoi ==5: ax0[1, (ndoi-5)].set_ylabel('Altitude [km]', fontsize=14)
                                #ax0[1, (ndoi-5)].set_ylim(1, 15)
                                if loc.lower() == 'fl0': ax0[1, (ndoi-5)].set_ylim(1, 15)
                                if loc.lower() == 'mlo': ax0[1, (ndoi-5)].set_ylim(3, 15)                                

                                if ndoi ==7 :ax0[1, (ndoi-5)].set_xlabel('VMR [x10$^3$ ppm]', fontsize=14)
                                ax0[1, (ndoi-5)].grid(True,which='both')
                                ax0[1, (ndoi-5)].tick_params(which='both',labelsize=14)
                                ax0[1, (ndoi-5)].set_title(da, fontsize=14)

                                ax0[1, (ndoi-5)].set_xlim(xmin=0)

                                

                                ax0[1, (ndoi-5)].text(0.05, 0.95, 'N = {}'.format(NobsFTS[dfValue][d]), va='center',transform=ax0[1, (ndoi-5)].transAxes,fontsize=14)
                                #ax0[1, (ndoi-5)].errorbar(sondePrf[dfValue][d]/1e3, sondealt[dfValue][d], xerr=sondePrf[dfValue][d]*0.1/1e3, color='k', markersize=0, linestyle='-', label='CFH', linewidth=2.0)
                                #ax0[1, (ndoi-5)].scatter(sondePrf[dfValue][d]/1e3, sondealt[dfValue][d],facecolors='white', s=35, color='k')
                                #ax0[1, (ndoi-5)].legend(prop={'size':11})

                            fig0.subplots_adjust(bottom=0.075,top=0.95, left=0.05, right=0.95) 
                            ndoi+=1

            if self.pdfsav: 
                self.pdfsav.savefig(fig0,dpi=200)
                plt.savefig('/data/iortega/Manuscripts/Fig/PrFSelected_'+self.loc.upper()+'.pdf', bbox_inches='tight')
            else:
                plt.show(block=False)
                #user_input = raw_input('Press any key to exit >>> ')
                #sys.exit()

        #---------------------------------
        #Figure: selected profiles (same grid, smoothed or not)
        #---------------------------------
        if doiplt:
            doidt = [dt.date(int(d[0:4]), int(d[4:6]), int(d[6:8])) for d in doiplt]
            doidt = np.asarray(doidt)
            
            fig0, ax0 = plt.subplots(2, 5, figsize=(15,10), sharey=True, sharex=False)
            ndoi = 0 

            for d, da in enumerate(sondedates[dfValue]):

                deltadoi = doidt - da

                if ndoi < len(doidt):
                
                    if deltadoi[ndoi] == dt.timedelta(0):

                        if ndoi < len(doidt):

                            if ndoi<=4:

                                if ndoi ==0: ax0[0, ndoi].set_ylabel('Altitude [km]', fontsize=14)

                                ax0[0, ndoi].plot(prfmean_day[dfValue][d]/1e3,alt, color='b',  linewidth=2.0, label='HR-FTIR')
                                ax0[0, ndoi].scatter(prfmean_day[dfValue][d]/1e3,alt,facecolors='white', s=35, color='b')
                                ax0[0, ndoi].fill_betweenx(alt,prfmean_day[dfValue][d]/1e3-prferr_day[dfValue][d]/1e3,prfmean_day[dfValue][d]/1e3+prferr_day[dfValue][d]/1e3,alpha=0.25,color='blue')

                                #ax0[0, ndoi].set_ylim(1, 15)
                                if loc.lower() == 'fl0': ax0[0, ndoi].set_ylim(1, 15)
                                if loc.lower() == 'mlo': ax0[0, ndoi].set_ylim(3, 15)
                                
                                ax0[0, ndoi].grid(True,which='both')
                                ax0[0, ndoi].tick_params(which='both',labelsize=14)
                                ax0[0, ndoi].set_title(da, fontsize=14)   

                                ax0[0, ndoi].plot(aPrf_day[dfValue][d]/1e3, alt, '-', color='gray', label='a priori', linewidth=2.0)
                                ax0[0, ndoi].scatter(aPrf_day[dfValue][d]/1e3, alt,facecolors='white', s=20, color='gray')
                                
                                ax0[0, ndoi].plot(Prfsonde_interp[dfValue][d]/1e3, alt,color='k', label='CFH', linewidth=2.0)
                                ax0[0, ndoi].scatter(Prfsonde_interp[dfValue][d]/1e3, alt,facecolors='white', s=20, color='k')

                                #ax0[0, ndoi].errorbar(sondePrf[dfValue][d]/1e3, sondealt[dfValue][d], xerr=sondePrf[dfValue][d]*0.1/1e3, color='k', markersize=0, linestyle='-', label='CFH', linewidth=2.0)
                                #ax0[0, ndoi].scatter(sondePrf[dfValue][d]/1e3, sondealt[dfValue][d],facecolors='white', s=35, color='k')

                                ax0[0, ndoi].text(0.05, 0.95, 'N = {}'.format(NobsFTS[dfValue][d]), va='center',transform=ax0[0, ndoi].transAxes,fontsize=14)

                                ax0[0, ndoi].set_xlim(xmin=0)

                                if ndoi == 0: ax0[0, ndoi].legend(prop={'size':11})

                            if ndoi>=5:

                                ax0[1, (ndoi-5)].plot(prfmean_day[dfValue][d]/1e3,alt, color='b',  linewidth=2.0, label='Mean (Retrieved)')
                                ax0[1, (ndoi-5)].scatter(prfmean_day[dfValue][d]/1e3,alt,facecolors='white', s=35, color='b')
                                ax0[1, (ndoi-5)].fill_betweenx(alt,prfmean_day[dfValue][d]/1e3-prferr_day[dfValue][d]/1e3,prfmean_day[dfValue][d]/1e3+prferr_day[dfValue][d]/1e3,alpha=0.25,color='blue')

                                if ndoi ==5: ax0[1, (ndoi-5)].set_ylabel('Altitude [km]', fontsize=14)
                                #ax0[1, (ndoi-5)].set_ylim(1, 15)
                                if loc.lower() == 'fl0': ax0[1, (ndoi-5)].set_ylim(1, 15)
                                if loc.lower() == 'mlo': ax0[1, (ndoi-5)].set_ylim(3, 15)

                                
                                if ndoi ==7 :ax0[1, (ndoi-5)].set_xlabel('VMR [x10$^3$ ppm]', fontsize=14)
                                ax0[1, (ndoi-5)].grid(True,which='both')
                                ax0[1, (ndoi-5)].tick_params(which='both',labelsize=14)
                                ax0[1, (ndoi-5)].set_title(da, fontsize=14)

                                ax0[1, (ndoi-5)].plot(aPrf_day[dfValue][d]/1e3, alt, '-', color='gray', label='a priori', linewidth=2.0)
                                ax0[1, (ndoi-5)].scatter(aPrf_day[dfValue][d]/1e3, alt,facecolors='white', s=20, color='gray')
                                
                                ax0[1, (ndoi-5)].plot(Prfsonde_interp[dfValue][d]/1e3, alt,color='k', label='sonde\n({:%H:%M})'.format(dtsonde), linewidth=2.0)
                                ax0[1, (ndoi-5)].scatter(Prfsonde_interp[dfValue][d]/1e3, alt,facecolors='white', s=20, color='k')

                                #ax0[1, (ndoi-5)].errorbar(sondePrf[dfValue][d]/1e3, sondealt[dfValue][d], xerr=sondePrf[dfValue][d]*0.1/1e3, color='k', markersize=0, linestyle='-', label='CFH', linewidth=2.0)
                                #ax0[1, (ndoi-5)].scatter(sondePrf[dfValue][d]/1e3, sondealt[dfValue][d],facecolors='white', s=35, color='k')
                                #ax0[1, (ndoi-5)].legend(prop={'size':11})

                                ax0[1, (ndoi-5)].text(0.05, 0.95, 'N = {}'.format(NobsFTS[dfValue][d]), va='center',transform=ax0[1, (ndoi-5)].transAxes,fontsize=14)

                                ax0[1, (ndoi-5)].set_xlim(xmin=0)

                            fig0.subplots_adjust(bottom=0.075,top=0.95, left=0.05, right=0.95) 
                            ndoi+=1

            if self.pdfsav: 
                self.pdfsav.savefig(fig0,dpi=200)
                plt.savefig('/data/iortega/Manuscripts/Fig/PrFSelected_'+self.loc.upper()+'_interp.pdf', bbox_inches='tight')
            else:
                plt.show(block=False)
                #user_input = raw_input('Press any key to exit >>> ')
                #sys.exit()

        #------------------------------------
        #FILTER
        #------------------------------------
        for df in diffT:

            indsBad = []

            for d, da in enumerate(FTSdates[df]):

                iBad = np.where( (np.abs(PrfDiffRel[df][d]) >= 300.0) | (PrfDist[df][d] > 500.0) )[0]
                #iBad = np.where( (np.abs(PrfDiffRel[df][d]) >= 300.0)  )[0]

                if len(iBad) >= 1:

                    indsBad.append(d)

            prfmean_day[df][indsBad,:] = np.nan#                   = np.delete(prfmean_day[df], indsBad, axis=0)
            prferr_day[df][indsBad, :] = np.nan#                      = np.delete(prferr_day[df], indsBad, axis=0)
            prfSTD_day[df][indsBad, :] = np.nan#                      = np.delete(prfSTD_day[df], indsBad, axis=0)
            aPrf_day[df][indsBad,:] = np.nan#                        = np.delete(aPrf_day[df], indsBad, axis=0)
            Airmassmean[df][indsBad,:] = np.nan#                     = np.delete(Airmassmean[df], indsBad, axis=0)
            avkSCFday[df][indsBad,:] = np.nan#                       = np.delete(avkSCFday[df], indsBad, axis=0)
            rPrfMolmean_day[df][indsBad,:] = np.nan#                 = np.delete(rPrfMolmean_day[df], indsBad, axis=0)
            prferrMol_day[df][indsBad,:] = np.nan#                   = np.delete(prferrMol_day[df], indsBad, axis=0)
            aPrfMol_day[df][indsBad, :] = np.nan#                     = np.delete(aPrfMol_day[df], indsBad, axis=0)
            #FTSdates[df][indsBad][:] = np.nan#                        = np.delete(FTSdates[df], indsBad)
            #NobsFTS[df][indsBad][:] = np.nan#                         = np.delete(NobsFTS[df], indsBad)

            Prfsonde_interp[df][indsBad, :] = np.nan#                 = np.delete(Prfsonde_interp[df], indsBad, axis=0)
            Prfsonde_sd_interp[df][indsBad, :] = np.nan#              = np.delete(Prfsonde_sd_interp[df], indsBad, axis=0)
            PrfsondeMol_interp[df][indsBad, :] = np.nan#              = np.delete(PrfsondeMol_interp[df], indsBad, axis=0)
            PrfsondeMol_sd_interp[df][indsBad, :] = np.nan#           = np.delete(PrfsondeMol_sd_interp[df], indsBad, axis=0)
            sondeairmass_a_interp[df][indsBad, :] = np.nan#           = np.delete(sondeairmass_a_interp[df], indsBad, axis=0)
            #sondedates[df][indsBad][:] = np.nan#                      = np.delete(sondedates[df], indsBad)
            Prfsondelon_interp[df][indsBad, :] = np.nan#              = np.delete(Prfsondelon_interp[df], indsBad, axis=0)
            Prfsondelat_interp[df][indsBad, :] = np.nan#              = np.delete(Prfsondelat_interp[df], indsBad, axis=0)
            
            PrfDiff[df][indsBad, :] = np.nan#                         = np.delete(PrfDiff[df], indsBad, axis=0)
            PrfDiffApr[df][indsBad, :] = np.nan#                      = np.delete(PrfDiffApr[df], indsBad, axis=0)
            PrfDiffRel[df][indsBad, :] = np.nan#                      = np.delete(PrfDiffRel[df], indsBad, axis=0)
            #PrfDiffRel[df]                      = np.delete(PrfDiffRel[df], indsBad, axis=0)
            PrfDiffAprRel[df][indsBad, :] = np.nan#                   = np.delete(PrfDiffAprRel[df], indsBad, axis=0)
            PrfDist[df][indsBad, :] = np.nan#                         = np.delete(PrfDist[df], indsBad, axis=0)

            for pn, pcol in enumerate(pCols):

                inds = np.where( (alt >= pcol[0]) & (alt <=pcol[1])  )[0]

                sondevmrP[str(df)+'_'+str(pn)]    = np.ma.average(np.ma.masked_invalid(Prfsonde_interp[df][:,inds]), axis=1,  weights=Airmassmean[df][:,inds])
                sondevmrsdP[str(df)+'_'+str(pn)]  = np.ma.average(np.ma.masked_invalid(Prfsonde_sd_interp[df][:, inds]), axis=1, weights=Airmassmean[df][:,inds])
                sondesumP[str(df)+'_'+str(pn)]    = np.nansum(np.ma.masked_invalid(Prfsonde_interp[df][:,inds])*np.ma.masked_invalid(Airmassmean[df][:,inds])*1e-6, axis=1)
                sondesumsdP[str(df)+'_'+str(pn)]  = np.nansum(np.ma.masked_invalid(Prfsonde_sd_interp[df][:,inds])*np.ma.masked_invalid(Airmassmean[df][:,inds])*1e-6, axis=1)
                
                FTSvmrP[str(df)+'_'+str(pn)]      = np.ma.average(np.ma.masked_invalid(prfmean_day[df][:,inds]), axis=1,  weights=np.ma.masked_invalid(Airmassmean[df][:,inds]))
                FTSvmrsdP[str(df)+'_'+str(pn)]    = np.ma.average(np.ma.masked_invalid(prferr_day[df][:,inds]), axis=1,  weights=np.ma.masked_invalid(Airmassmean[df][:,inds]))
                FTSAprvmrP[str(df)+'_'+str(pn)]   = np.ma.average(np.ma.masked_invalid(aPrf_day[df][:,inds]), axis=1,  weights=np.ma.masked_invalid(Airmassmean[df][:,inds]))
                FTSsumP[str(df)+'_'+str(pn)]      = np.nansum(np.ma.masked_invalid(rPrfMolmean_day[df][:,inds]), axis=1)
                FTSsumsdP[str(df)+'_'+str(pn)]    = np.nansum(np.ma.masked_invalid(prferrMol_day[df][:,inds]), axis=1)
                FTSAprsumP[str(df)+'_'+str(pn)]   = np.nansum(np.ma.masked_invalid(aPrfMol_day[df][:,inds]), axis=1)

                sondelatP[str(df)+'_'+str(pn)]    = np.ma.mean(np.ma.masked_invalid(Prfsondelat_interp[df][:,inds]), axis=1)
                sondelonP[str(df)+'_'+str(pn)]    = np.ma.mean(np.ma.masked_invalid(Prfsondelon_interp[df][:,inds]), axis=1)

                #distP[str(df)+'_'+str(pn)]        = [mf.haversine(np.nanmean(lonmean[df]), np.mean(latmean[df]), lo, la ) for (lo, la) in  zip(sondelonP[str(df)+'_'+str(pn)], sondelatP[str(df)+'_'+str(pn)])]
                distP[str(df)+'_'+str(pn)]        = np.ma.mean(np.ma.masked_invalid(PrfDist[df][:,inds]), axis=1)

                    # for pn, pcol in enumerate(pCols):

                    #     sondevmrP[str(df)+'_'+str(pn)]    = np.delete(sondevmrP[str(df)+'_'+str(pn)], indsBad)
                    #     sondevmrsdP[str(df)+'_'+str(pn)]  = np.delete(sondevmrsdP[str(df)+'_'+str(pn)], indsBad)
                    #     sondesumP[str(df)+'_'+str(pn)]    = np.delete(sondesumP[str(df)+'_'+str(pn)], indsBad)
                    #     sondesumsdP[str(df)+'_'+str(pn)]  = np.delete(sondesumsdP[str(df)+'_'+str(pn)], indsBad)
                        
                    #     FTSvmrP[str(df)+'_'+str(pn)]      = np.delete(FTSvmrP[str(df)+'_'+str(pn)], indsBad)
                    #     FTSvmrsdP[str(df)+'_'+str(pn)]    = np.delete(FTSvmrsdP[str(df)+'_'+str(pn)], indsBad)
                    #     FTSAprvmrP[str(df)+'_'+str(pn)]   = np.delete(FTSAprvmrP[str(df)+'_'+str(pn)], indsBad)
                    #     FTSsumP[str(df)+'_'+str(pn)]      = np.delete(FTSsumP[str(df)+'_'+str(pn)], indsBad)
                    #     FTSsumsdP[str(df)+'_'+str(pn)]    = np.delete(FTSsumsdP[str(df)+'_'+str(pn)], indsBad)
                    #     FTSAprsumP[str(df)+'_'+str(pn)]   = np.delete(FTSAprsumP[str(df)+'_'+str(pn)] , indsBad)
                    #     distP[str(df)+'_'+str(pn)]        = np.delete(distP[str(df)+'_'+str(pn)] , indsBad)

        #-------------------------------------------------
        #Fit: Profile position diference at every coincidence time   
        #-------------------------------------------------
        fig = plt.figure(figsize=(14,9))

        outer_grid = gridspec.GridSpec(3, 3, wspace=0.15, hspace=0.1)
            
        for di, df in enumerate(diffT):
            ax = plt.Subplot(fig, outer_grid[di])

            for pr in PrfDist[df]:
                ax.plot(pr, alt,'k')

            ax.set_xlim(0, 200)
            ax.grid(True, color='gray', alpha=0.5)
            ax.tick_params(which='both',labelsize=10)
            ax.annotate('{} min'.format((df)), xy=(0.9, 0.9), xycoords='axes fraction', fontsize=12, ha='left')

            ax.tick_params(axis = 'both', which = 'minor', labelsize = 0)
      
            ax.xaxis.set_tick_params(which='minor',labelbottom='off')

            fig.add_subplot(ax)

        all_axes = fig.get_axes()
        #show only the outside spines
        for ax in all_axes:
            for sp in ax.spines.values():
                sp.set_visible(False)
                plt.setp(ax.get_xticklabels(), visible=False)
            if ax.is_first_row():
                ax.spines['top'].set_visible(True)
            if ax.is_last_row():
                ax.spines['bottom'].set_visible(True)
                plt.setp(ax.get_xticklabels(), visible=True)
            if ax.is_first_col():
                ax.spines['left'].set_visible(True)
            if ax.is_last_col():
                ax.spines['right'].set_visible(True)

        if (9 % 2 == 1): #even
            all_axes[-2].spines['bottom'].set_visible(True)
            plt.setp(all_axes[-2].get_xticklabels(), visible=True)
            all_axes[-2].set_zorder(1)

        all_axes[-1].set_xlabel('$\Delta$x [km]')
        all_axes[-2].set_xlabel('$\Delta$x [km]')
        all_axes[-3].set_xlabel('$\Delta$x [km]')
        
        #fig.text(0.03, 0.5, ytypeStr +' ['+unitsStr+']', fontsize=16, va='center', rotation='vertical')
        fig.subplots_adjust(left=0.08, bottom=0.05, right=0.975, top=0.95)

        if self.pdfsav:
            pdfsav.savefig(fig,dpi=200)
            
        else:
            plt.show(block=False)
            #user_input = raw_input('Press any key to exit >>> ')
            #sys.exit()
        #-------------------------------------------------------
        #Test sonde
        #-------------------------------------------------------
        fig1, ax1  = plt.subplots()

        for df in diffT:

            PrfDIstMean = np.nanmean(PrfDist[df], axis=0)
            PrfDIstStd  = np.nanstd(PrfDist[df], axis=0)
            #ax1.scatter(diffT, PrfDIstMean, color='k',s=60, label='# of Dates')
            ax1.plot(PrfDIstMean, alt, label=df)
           
        ax1.grid(True) 
        ax1.tick_params(which='both',labelsize=14)
        ax1.set_ylabel('Altitude [km]', fontsize=14)
        ax1.set_ylim(bottom=0)
        ax1.legend(prop={'size':11}, loc=2)
        #ax1.set_xlim(left=0, right=np.max(diffT) + 10)
        ax1.set_xlabel('$\Delta$x [km]', fontsize=14)  

        if self.pdfsav: 
            self.pdfsav.savefig(fig1,dpi=200)
            plt.savefig('/data/iortega/Manuscripts/Fig/StatTime_'+self.loc.upper()+'.pdf', bbox_inches='tight')
        else:
            plt.show(block=False)
            #user_input = raw_input('Press any key to exit >>> ')
            #sys.exit()
        # sondelat    = np.asarray(self.sonde['latFleout_a'])
        # sondelon    = np.asarray(self.sonde['lonFleout_a'])

        # sondealt  = np.asarray(self.sonde['Alt_km_a'])
        # sondeprf  = np.asarray(self.sonde['H2Omr_ppmv_a'])

        # print sondealt.shape, sondeprf.shape

        # opp   = []
        # for ii, r in enumerate(sondeprf):
        #     opp.append(np.average(sondealt[ii], weights= r))

        # opp = np.asarray(opp)

        # print np.mean(opp)
        # exit()
        
        # fig, ax1       = plt.subplots()

        # for i in range(sondelat.shape[0]):
    
        #     ax1.plot(sondelon[i][:], sondelat[i][:], linewidth=0.75)
        
        # ax1.set_ylabel('Latitude')
        # ax1.set_xlabel('Longitude')        
        # ax1.grid(True,which='both')
        # ax1.set_title('Location', fontsize=12)
        # ax1.tick_params(labelsize=12)
        # ax1.set_xlim(-106, -102)
        # ax1.set_ylim(38, 42)

        # if self.pdfsav: 
        #     self.pdfsav.savefig(fig,dpi=200)
        #     plt.close(fig)
        # else:
        #     plt.show(block=False)
        #     user_input = raw_input('Press any key to exit >>> ')
        #     sys.exit()


        #------------------------------------
        #Fig: Stats with Coincident Interval
        #------------------------------------
        fig1, ax1  = plt.subplots(2, figsize=(9,8), sharex=True)
        clr = mf.clrplt()


        for pn, pcol in enumerate(pCols):

            slope_all    =  []
            inter_all    =  []
            r_all        =  []

            slope_odr    = []
            inter_odr    = []
            slopeErr_odr = []
            interErr_odr = []

            stdDF        = []
            Npnts        = []
            Npnts2       = []
            biasDF       = []
            rmseDF       = []

            for df in diffT:

                inds = np.where( (alt > pcol[0]) & (alt <=pcol[1])  )[0]

                x = PrfDiffRel[df][:, inds]
                x = np.asarray(x)
                x = x.flatten()

                biasDF.append(np.nansum( x) / len(x))

                ftsPrf = [i[inds] for i in prfmean_day[df]]
                sonPrf = [i[inds] for i in Prfsonde_interp[df]]

                st = prfSTD_day[df][:, inds]
                st = np.asarray(st)
                st = st.flatten()

                stdDF.append(np.nanmean(st)/np.nanmean(ftsPrf))


                ftsPrf = np.asarray(ftsPrf)
                ftsPrf = ftsPrf.flatten()

                sonPrf = np.asarray(sonPrf)
                sonPrf = sonPrf.flatten()

                # slope, intercept, r_value, p_value, std_err = stats.linregress(sonPrf, ftsPrf)
                # slope_all.append(slope)
                # inter_all.append(intercept)
                # r_all.append(r_value)

                #------------------------------
                # Root Mean Square Error (RMSE)
                #------------------------------
                ss_res = np.nansum( x**2)
                rmseDF.append(np.sqrt( ss_res / len(x) ))

                #stdDF.append(np.std(x))

                #print 'bias (VMR) = {}'.format(np.nansum( x) / len(x))
                slope, intercept, r_value, p_value, std_err = stats.linregress(sondevmrP[str(df)+'_'+str(pn)], FTSvmrP[str(df)+'_'+str(pn)])
                slope_all.append(slope)
                inter_all.append(intercept)
                r_all.append(r_value)

                odr, odrErr = mf.orthoregress(sondevmrP[str(df)+'_'+str(pn)], FTSvmrP[str(df)+'_'+str(pn)], xerr=sondevmrsdP[str(df)+'_'+str(pn)], yerr=FTSvmrsdP[str(df)+'_'+str(pn)], InError=True)
                slope_odr.append(float(odr[0]))
                inter_odr.append(float(odr[1]))
                slopeErr_odr.append(float(odrErr[0]))
                interErr_odr.append(float(odrErr[1]))

                #Npnts.append(len(x))

                Npnts.append(np.sum(NobsFTS[df]))
                Npnts2.append(len(NobsFTS[df]))


            if pn == 0:
                ax1[0].scatter(diffT, Npnts2, color='k',s=60, label='# of Dates')
                ax1[0].plot(diffT, Npnts2,  color='k')
                ax1[0].scatter(diffT, Npnts, color='blue',s=60, label='# of Profiles')
                ax1[0].plot(diffT, Npnts,  color='blue')
                #ax1[0].set_title('(a)', loc='left', fontsize=14)

            #ax1[1].scatter(diffT, r_all, color=clr[pn], s=60)
            #ax1[1].plot(diffT, r_all,   color=clr[pn])
            ax1[1].scatter(diffT, np.asarray(stdDF)*100., color=clr[pn], s=60, label= str(pcol[0])+'-'+str(pcol[1])+' km')
            ax1[1].plot(diffT, np.asarray(stdDF)*100.,   color=clr[pn])
            #ax1[1].set_title('(b)', loc='left', fontsize=14)

            #ax1[2].scatter(diffT, r_all, color=clr[pn], s=60)
            #ax1[2].plot(diffT, r_all,   color=clr[pn])
            #ax1[2].set_title('(c)', loc='left', fontsize=14)


            #ax1[1].errorbar(diffT,slope_odr, yerr=slopeErr_odr, color=clr[pn], markersize=60, linestyle='-')
            #ax1[1].scatter(diffT, slope_odr, s=60,color=clr[pn], label= str(pcol))

            #ax1[2].errorbar(diffT,inter_odr, yerr=interErr_odr,color=clr[pn], markersize=60, linestyle='-', label= str(pcol))
            #ax1[1].scatter(diffT, inter_odr, color=clr[pn],  s=60)

        ax1[0].grid(True) 
        ax1[0].tick_params(which='both',labelsize=14)
        ax1[0].set_ylabel('Number of Dates/Profiles', fontsize=14)
        ax1[0].set_ylim(bottom=0)
        ax1[0].legend(prop={'size':11}, loc=2)

        ax1[1].grid(True)
        ax1[1].set_ylabel('HR-FTIR variability [%]', fontsize=14)
        ax1[1].tick_params(which='both',labelsize=14)
        ax1[1].set_ylim(bottom=0)
        ax1[1].legend(prop={'size':11}, loc=2)
        ax1[1].set_xlim(left=0, right=np.max(diffT) + 10)
        ax1[1].set_xlabel('HR-FTIR/Sonde coincidence interval [min]', fontsize=14)

        #ax1[2].grid(True)
        #ax1[2].set_ylabel('r-value', fontsize=14)            
        #ax1[2].tick_params(which='both',labelsize=14)
        #ax1[2].set_xlim(left=0, right=np.max(diffT) + 10)
        #ax1[2].set_xlabel('HR-FTIR/Sonde coincidence interval [min]', fontsize=14)
        #ax1[2].set_xscale('log')

        if self.pdfsav: 
            self.pdfsav.savefig(fig1,dpi=200)
            plt.savefig('/data/iortega/Manuscripts/Fig/StatTime_'+self.loc.upper()+'.pdf', bbox_inches='tight')
        else:
            plt.show(block=False)
            #user_input = raw_input('Press any key to exit >>> ')
            #sys.exit()

        #---------------------------------
        # Figure: Profile with the difference as function of Time
        #---------------------------------

        clr = mf.clrplt()

        #fig1, (ax1, ax2)  = plt.subplots(1, 2,figsize=(10,8), sharey=True)
        fig1, ax1  = plt.subplots(figsize=(6,7))
        
        for ii, df in enumerate(diffT):

            PrfDiffMean       = np.nanmean(PrfDiff[df], axis=0)
            PrfDiffAprMean    = np.nanmean(PrfDiffApr[df], axis=0)
            PrfDifRelfMean    = np.nanmean(PrfDiffRel[df], axis=0)
            PrfDiffAprRelMean = np.nanmean(PrfDiffAprRel[df], axis=0)

            #ax1.plot(PrfDiffMean/1e3, alt, '-', color=clr[ii], label=str(df) + '-'+str(len(NobsFTS[df])) + '-'+str(np.sum(NobsFTS[df])) , linewidth=3)
            ax1.plot(PrfDiffMean/1e3, alt, '-', color=clr[ii], label=str(df)+' min', linewidth=3)
            ax1.scatter(PrfDiffMean/1e3, alt, facecolors='white', color=clr[ii], s=45)

            #ax0[0, ndoi].scatter(prfmean_day[dfValue][d]/1e3,alt,facecolors='white', s=35, color='b')
            #ax0[0, ndoi].fill_betweenx(alt,prfmean_day[dfValue][d]/1e3-prferr_day[dfValue][d]/1e3,prfmean_day[dfValue][d]/1e3+prferr_day[dfValue][d]/1e3,alpha=0.25,color='blue')

            #ax2.plot(PrfDifRelfMean, alt, '-', color=clr[ii], linewidth=3)
            #ax1.plot(PrfDiffAprMean/1e3, alt, '--', color=clr[ii], linewidth=3)
            #ax2.plot(PrfDiffAprRelMean, alt, '--',  color=clr[ii], linewidth=3)
            
        if loc.lower() == 'fl0': ax1.set_ylim(1, 15)
        if loc.lower() == 'mlo': ax1.set_ylim(3, 15)

        ax1.set_xlim(xmin=-2, xmax=1.5)
        ax1.grid(True,which='both')
        ax1.set_ylabel('Altitude [km]', fontsize=14)
        ax1.set_xlabel('VMR [x10$^3$ ppm$_v$]', fontsize=14)   
        ax1.tick_params(which='both',labelsize=14)
        ax1.axvline(x=0, linestyle='--', linewidth=2, color='k')

        ax1.legend(prop={'size':14}, loc=2)

        xval = range(-2,3)
        for pn, pcol in enumerate(pCols):
            ax1.axhline(y=pcol[0], linewidth=2, color='gray', alpha=0.5)  
            ax1.axhline(y=pcol[1], linewidth=2, color='gray', alpha=0.5)
            ax1.fill_between(xval,pcol[0],pcol[1],alpha=0.15, color='0.5')  


        #ax2.set_ylim(0, 15)
        # if loc.lower() == 'fl0': ax2.set_ylim(1, 15)
        # if loc.lower() == 'mlo': ax2.set_ylim(3, 15)
        # ax2.set_xlim(-50, 100)
        # ax2.grid(True,which='both')
        # ax2.set_xlabel('Relative Difference [%]', fontsize=14)   
        # ax2.tick_params(which='both',labelsize=14)
        # ax2.axvline(x=0, linestyle='--', linewidth=2, color='k')

        if self.pdfsav: 
            self.pdfsav.savefig(fig1,dpi=200)
        else:
            plt.show(block=False)
            #user_input = raw_input('Press any key to exit >>> ')
            #sys.exit()

        #---------------------------------
        # Figure: Profile with the difference of all profiles for fixed time
        #---------------------------------
        # fig1, (ax1, ax2)  = plt.subplots(1, 2,figsize=(10,8), sharey=True)
        
        # PrfDiffMean       = np.nanmean(PrfDiff[dfValue], axis=0)
        # PrfDiffAprMean    = np.nanmean(PrfDiffApr[dfValue], axis=0)
        # PrfDifRelfMean    = np.nanmean(PrfDiffRel[dfValue], axis=0)
        # PrfDiffAprRelMean = np.nanmean(PrfDiffAprRel[dfValue], axis=0)


        # for i, diff in enumerate(PrfDiff[dfValue]):
        #     ax1.plot(diff/1e3, alt, '-', color='gray', linewidth=1, alpha=0.5)
        #     ax2.plot(PrfDiffRel[dfValue][i], alt, '-', color='gray', linewidth=1, alpha=0.5)
            

        # ax1.plot(PrfDiffMean/1e3, alt, '-', color='b', label='Mean (retrieved)', linewidth=3)
        # ax1.plot(PrfDiffAprMean/1e3, alt, '--', color='b',label='Mean (a priori)', linewidth=3)
        # #ax1.plot(tot_err_Mean, alt, '-', color='r', label='a priori', linewidth=3)
        # #ax1.set_ylim(0, 15)
        # if loc.lower() == 'fl0': ax1.set_ylim(1, 15)
        # if loc.lower() == 'mlo': ax1.set_ylim(3, 15)

        # #ax1.set_xlim(xmin=0)
        # ax1.grid(True,which='both')
        # ax1.set_ylabel('Altitude [km]', fontsize=14)
        # ax1.set_xlabel('VMR [x10$^3$ ppm$_v$]', fontsize=14)   
        # ax1.tick_params(which='both',labelsize=14)
        # ax1.axvline(x=0, linestyle='--', linewidth=2, color='k')

        # ax1.legend(prop={'size':14})

        # ax2.plot(PrfDifRelfMean, alt, '-',  color='b',  linewidth=3)
        # ax2.plot(PrfDiffAprRelMean, alt, '--',  color='b', linewidth=3)
        # #ax2.set_ylim(0, 15)
        # if loc.lower() == 'fl0': ax2.set_ylim(1, 15)
        # if loc.lower() == 'mlo': ax2.set_ylim(3, 15)
        # #ax2.set_xlim(-100, 100)
        # ax2.grid(True,which='both')
        # ax2.set_xlabel('Relative Difference [%]', fontsize=14)   
        # ax2.tick_params(which='both',labelsize=14)
        # ax2.axvline(x=0, linestyle='--', linewidth=2, color='k')

        # if self.pdfsav: 
        #     self.pdfsav.savefig(fig1,dpi=200)
        # else:
        #     plt.show(block=False)
        #     user_input = raw_input('Press any key to exit >>> ')
        #    sys.exit()

            
        

        #---------------------------------
        #
        #---------------------------------
        if self.loc == 'mlo': altHist  = pCols[0:]# [ [3, 5], [5, 7], [7, 9], [9, 11], [11, 13] ]
        else: altHist  = pCols[0:-1]# [[0, 3], [3, 5], [5, 7], [7, 9], [9,11] ]
        
        numBins = 20

        if self.loc == 'fl0': dfValue = 15
        elif self.loc == 'mlo': dfValue = 90

        fig = plt.figure(figsize=(10,7.5))

        for df in diffT:

            if df == dfValue:

                #---------------------------------
                # Figure: HISTOGRAMS FOR DIFFERENT ALTITUDE LAYERS
                #---------------------------------
                for p, h in enumerate(altHist):

                    inds = np.where( (alt >= h[0]) & (alt <=h[1])  )[0]

                    #inds = np.where( (alt > pcol[0]) & (alt <=pcol[1])  )[0]

                    x = PrfDiffRel[df][:, inds]
                    x = np.asarray(x)
                    x = x.flatten()

                    x = x[~np.isnan(x)]

                    # print x

                    # ibad = np.where(np.isnan(x)) [0]
                    # print ibad
                    # x = np.delete(x, ibad)

                    # print x


                    biasCalc = np.nansum( x) / len(x)

                    gs1 = gridspec.GridSpec(1, 3)

                    if p == 0: gs1.update(left=0.075, right=0.5, top=0.95, bottom=0.69,  wspace=0.05,  hspace=0.08)
                    if p == 2: gs1.update(left=0.075, right=0.5, top = 0.64, bottom=0.38, wspace=0.05, hspace=0.08)
                    if p == 4: gs1.update(left=0.075, right=0.5, top = 0.33, bottom=0.07, wspace=0.05, hspace=0.08)
                    
                    if p == 1: gs1.update(left=0.555, right=0.98, top=0.95, bottom=0.69,  wspace=0.05,  hspace=0.08)
                    if p == 3: gs1.update(left=0.555, right=0.98, top=0.64, bottom=0.38, wspace=0.05,   hspace=0.08)
                    if p == 5: gs1.update(left=0.555, right=0.98, top=0.33, bottom=0.07, wspace=0.05,   hspace=0.08)

                    ax1 = plt.subplot(gs1[0:1, :])

                    weights = np.ones_like(x)/len(x)
                    #plt.hist(myarray, weights=weights)
                    #normed=True means that the total area under the histogram is equal to 1 but the sum of heights is not equal to 1    
                    n, bins, patches = ax1.hist(x, bins=numBins, normed=True, alpha=0.8 ,histtype='step', fill=True, color='k')  #color='green'
                    P.setp(patches, 'facecolor', 'g', 'alpha', 0.8)
                    mu = np.nanmean(x)
                    sigma = np.nanstd(x)
                    #add a line showing the expected distribution
                    y = P.normpdf( bins, mu, sigma)
                    l = ax1.plot(bins, y, 'k--', linewidth=3)


                    ax1.grid(True,which='both')
                    ax1.tick_params(which='both',labelsize=12, labelbottom='off')
                    #ax1.set_ylabel('Probability', fontsize=12)
                    ax1.set_title(str(h[0])+' - '+str(h[1])+' km',
                              horizontalalignment='center', verticalalignment='baseline', fontsize=14)
                    ax1.set_xlim(-100, 100)

                #------
                    #ax1.set_xlabel('Relative Difference [%]', fontsize=12)
                    #ax1.text(0.55,0.875,"Bias: {0:.3f}$\pm${1:.3}".format(mu, sigma),transform=ax1.transAxes,  fontsize=13)
                   
                    if (p == 4) or (p == 5): 
                        ax1.set_xlabel('Difference [%]', fontsize=14)
                        ax1.tick_params(labelbottom='on')
                    if (p == 0) or (p == 2) or (p == 4): ax1.set_ylabel('Probability', fontsize=14)

                    print '\nAltitude Layer: '+str(h[0])+' - '+str(h[1])+' km'
                    print 'bias (%)        = {0:.2f}'.format(biasCalc)
                    print 'Histogram (%)   = {0:.2f} +/- {1:.2f}'.format(mu, sigma)


        
                if self.pdfsav: 
                    self.pdfsav.savefig(fig,dpi=200)
                    plt.savefig('/data/iortega/Manuscripts/Fig/HistLayers_'+self.loc.upper()+'.pdf', bbox_inches='tight')
                else:
                    plt.show(block=False)
                    #user_input = raw_input('Press any key to exit >>> ')
                    #sys.exit()


        

        #---------------------------------
        #
        #---------------------------------
        fig = plt.figure(figsize=(10,7.5))

        for df in diffT:

            if df == dfValue:

                #---------------------------------
                # Figure: HISTOGRAMS FOR DIFFERENT ALTITUDE LAYERS
                #---------------------------------
                for p, h in enumerate(altHist):

                    inds = np.where( (alt >= h[0]) & (alt <=h[1])  )[0]

                    #x = [i[inds] for i in np.ma.masked_invalid(PrfDiff[df])]
                    #x = np.asarray(x)
                    #x = x.flatten()

                    x = PrfDiff[df][:, inds]
                    x = np.asarray(x)
                    x = x.flatten()

                    x = x[~np.isnan(x)]

                    biasCalc = np.nansum( x) / len(x)

                    gs1 = gridspec.GridSpec(1, 3)

                    if p == 0: gs1.update(left=0.075, right=0.5, top=0.95, bottom=0.69,  wspace=0.05,  hspace=0.08)
                    if p == 2: gs1.update(left=0.075, right=0.5, top = 0.64, bottom=0.38, wspace=0.05, hspace=0.08)
                    if p == 4: gs1.update(left=0.075, right=0.5, top = 0.33, bottom=0.07, wspace=0.05, hspace=0.08)
                    
                    if p == 1: gs1.update(left=0.555, right=0.98, top=0.95, bottom=0.69,  wspace=0.05,  hspace=0.08)
                    if p == 3: gs1.update(left=0.555, right=0.98, top=0.64, bottom=0.38, wspace=0.05,   hspace=0.08)
                    if p == 5: gs1.update(left=0.555, right=0.98, top=0.33, bottom=0.07, wspace=0.05,   hspace=0.08)

                    ax1 = plt.subplot(gs1[0:1, :])

                    weights = np.ones_like(x/1e3)/len(x)
                    #plt.hist(myarray, weights=weights)
                    #normed=True means that the total area under the histogram is equal to 1 but the sum of heights is not equal to 1    
                    n, bins, patches = ax1.hist(x/1e3, bins=numBins, normed=True, alpha=0.8 ,histtype='step', fill=True, color='k')  #color='green'
                    P.setp(patches, 'facecolor', 'g', 'alpha', 0.8)
                    mu = np.nanmean(x/1e3)
                    sigma = np.std(x/1e3)
                    #add a line showing the expected distribution
                    y = P.normpdf( bins, mu, sigma)
                    l = ax1.plot(bins, y, 'k--', linewidth=3)


                    ax1.grid(True,which='both')
                    ax1.tick_params(which='both',labelsize=12, labelbottom='off')
                    #ax1.set_ylabel('Probability', fontsize=12)
                    ax1.set_title(str(h[0])+' - '+str(h[1])+' km',
                              horizontalalignment='center', verticalalignment='baseline', fontsize=14)
                    ax1.set_xlim(-4, 4)

                #------
                    #ax1.set_xlabel('Relative Difference [%]', fontsize=12)
                    #ax1.text(0.55,0.875,"Bias: {0:.3f}$\pm${1:.3}".format(mu, sigma),transform=ax1.transAxes,  fontsize=13)
                   
                    if (p == 4) or (p == 5): 
                        ax1.set_xlabel('VMR [x10$^3$ ppm]', fontsize=14)
                        ax1.tick_params(labelbottom='on')
                    if (p == 0) or (p == 2) or (p == 4): ax1.set_ylabel('Probability', fontsize=14)

                    print '\nAltitude Layer: '+str(h[0])+' - '+str(h[1])+' km'
                    print 'bias (vmr)        = {0:.2f}'.format(biasCalc/1e3)
                    print 'Histogram (vmr)   = {0:.2f} +/- {1:.2f}'.format(mu, sigma)


        
                if self.pdfsav: 
                    self.pdfsav.savefig(fig,dpi=200)
                    plt.savefig('/data/iortega/Manuscripts/Fig/HistLayers_'+self.loc.upper()+'.pdf', bbox_inches='tight')
                else:
                    plt.show(block=False)
                    #user_input = raw_input('Press any key to exit >>> ')
                    #sys.exit()


        #---------------------------------
        #
        #---------------------------------

        fig = plt.figure(figsize=(10,7.5))

        for df in diffT:

            if df == dfValue:

                #---------------------------------
                # Figure: bias as function of Distance
                #---------------------------------
                for p, h in enumerate(altHist):

                    res = np.true_divide( (FTSvmrP[str(df)+'_'+str(p)] - sondevmrP[str(df)+'_'+str(p)]), sondevmrP[str(df)+'_'+str(p)])*100.
                    #res = FTSsumP[str(df)+'_'+str(p)] - sondesumP[str(df)+'_'+str(p)]

                    gs1 = gridspec.GridSpec(1, 3)

                    if p == 0: gs1.update(left=0.075, right=0.5, top=0.95, bottom=0.69,  wspace=0.05,  hspace=0.08)
                    if p == 2: gs1.update(left=0.075, right=0.5, top = 0.64, bottom=0.38, wspace=0.05, hspace=0.08)
                    if p == 4: gs1.update(left=0.075, right=0.5, top = 0.33, bottom=0.07, wspace=0.05, hspace=0.08)
                    
                    if p == 1: gs1.update(left=0.555, right=0.98, top=0.95, bottom=0.69,  wspace=0.05,  hspace=0.08)
                    if p == 3: gs1.update(left=0.555, right=0.98, top=0.64, bottom=0.38, wspace=0.05,   hspace=0.08)
                    if p == 5: gs1.update(left=0.555, right=0.98, top=0.33, bottom=0.07, wspace=0.05,   hspace=0.08)

                    ax1 = plt.subplot(gs1[0:1, :])

                    weights = np.ones_like(x)/len(x)
                    #plt.hist(myarray, weights=weights)
                    #normed=True means that the total area under the histogram is equal to 1 but the sum of heights is not equal to 1  

                    ax1.plot(distP[str(df)+'_'+str(p)], res, '-',  linewidth=3)
                    ax1.scatter(distP[str(df)+'_'+str(p)], res, facecolors='white', s=45)  

                    
                    ax1.grid(True,which='both')
                    ax1.tick_params(which='both',labelsize=12, labelbottom='off')
                    #ax1.set_ylabel('Probability', fontsize=12)
                    ax1.set_title(str(h[0])+' - '+str(h[1])+' km',
                              horizontalalignment='center', verticalalignment='baseline', fontsize=14)
                    ax1.set_xlim(0, 150)
                    ax1.set_ylim(-100, 100)

                #------
                    #ax1.set_xlabel('Relative Difference [%]', fontsize=12)
                    #ax1.text(0.55,0.875,"Bias: {0:.3f}$\pm${1:.3}".format(mu, sigma),transform=ax1.transAxes,  fontsize=13)
                   
                    if (p == 4) or (p == 5): 
                        ax1.set_xlabel('$\Delta$x [km]', fontsize=14)
                        ax1.tick_params(labelbottom='on')
                    if (p == 0) or (p == 2) or (p == 4): ax1.set_ylabel('Probability', fontsize=14)

        
                if self.pdfsav: 
                    self.pdfsav.savefig(fig,dpi=200)
                    plt.savefig('/data/iortega/Manuscripts/Fig/HistLayers_'+self.loc.upper()+'.pdf', bbox_inches='tight')
                else:
                    plt.show(block=False)
                    #user_input = raw_input('Press any key to exit >>> ')
                    #sys.exit()
                    


        #------------------------------------
        # Plot time series of partial columns
        #------------------------------------
        #DateFmt      = DateFormatter('%m/%d/%Y')
        DateFmt      = DateFormatter('%Y')
        dayLc        = DayLocator()
        yearsLc      = YearLocator()
        monthLc      = MonthLocator()

        #------------------------------------
        # Initialize variable for standard linear regression
        #------------------------------------
        slope_all     =  []
        inter_all     =  []
        r_all         =  []

        slopeApr_all  =  []
        interApr_all  =  []
        rApr_all      =  []

        #------------------------------------
        # Initialize variable for Orthogonal Distance Regression (ODR)
        #------------------------------------
        slope_odr     =  []
        inter_odr     =  []
        slopeErr_odr  =  []
        interErr_odr  =  []

        slopeApr_odr  =  []
        interApr_odr  =  []

      
        #fig1, ax1  = plt.subplots(len(pCols),1,figsize=(9,12), sharex=True)
        #fig2, ax2  = plt.subplots(len(pCols),1,figsize=(9,12), sharex=True)

        #fig3, ax3  = plt.subplots(len(pCols),1,figsize=(7,13))
        #fig4, ax4  = plt.subplots(len(pCols),1,figsize=(7,13))

        print '\n'

        for df in diffT:

            if df == dfValue:
       
                for pn, pcol in enumerate(pCols):

                    #-------------------------------------------
                    #
                    #-------------------------------------------

                    slope, intercept, r_value, p_value, std_err = stats.linregress(sondevmrP[str(df)+'_'+str(pn)], FTSvmrP[str(df)+'_'+str(pn)])
                    slope_all.append(slope)
                    inter_all.append(intercept)
                    r_all.append(r_value)

                    mx = np.mean(sondevmrP[str(df)+'_'+str(pn)])
                    sx2 = np.sum(((sondevmrP[str(df)+'_'+str(pn)]-mx)**2))
                    sd_intercept = std_err * np.sqrt(1./len(sondevmrP[str(df)+'_'+str(pn)]) + mx*mx/sx2)
                    sd_slope = std_err * np.sqrt(1./sx2)


                    odr, odrErr = mf.orthoregress(sondevmrP[str(df)+'_'+str(pn)], FTSvmrP[str(df)+'_'+str(pn)], xerr=sondevmrsdP[str(df)+'_'+str(pn)], yerr=FTSvmrsdP[str(df)+'_'+str(pn)], InError=True)
                    slope_odr.append(float(odr[0]))
                    inter_odr.append(float(odr[1]))
                    slopeErr_odr.append(float(odrErr[0]))
                    interErr_odr.append(float(odrErr[1]))

                    print 'ODR (N= {0:4}) - Retrieved (Col:{1:}): slope= {2:}, off= {3:}'.format(len(FTSvmrP[str(df)+'_'+str(pn)]), pcol, odr,odrErr)  
                    
                   
                    slope, intercept, r_value, p_value, std_err = stats.linregress(sondevmrP[str(df)+'_'+str(pn)], FTSAprvmrP[str(df)+'_'+str(pn)])
                    slopeApr_all.append(slope)
                    interApr_all.append(intercept)
                    rApr_all.append(r_value)

                    odr, odrErr = mf.orthoregress(sondevmrP[str(df)+'_'+str(pn)], FTSAprvmrP[str(df)+'_'+str(pn)], xerr=sondevmrsdP[str(df)+'_'+str(pn)], yerr=FTSAprvmrP[str(df)+'_'+str(pn)]*0., InError=True)
                    slopeApr_odr.append(float(odr[0]))
                    interApr_odr.append(float(odr[1]))

                    print 'ODR - Apriori (Col:{0:4}):slope= {1:}, off= {2:}'.format(pcol, odr,odrErr)  

                    #ax3[pn].text(0.025,0.92,"slope: {0:.2f}+/-{1:.2}".format(slope, sd_slope),transform=ax3[pn].transAxes,  fontsize=11)
                    #ax3[pn].text(0.025,0.85,"Intercept: {0:.2f}+/-{1:.2}".format(intercept, sd_intercept),transform=ax3[pn].transAxes,  fontsize=11)
                    #ax3[pn].text(0.025,0.78,"r-value: {0:.2f}".format(r_value),transform=ax3[pn].transAxes,  fontsize=11)        

            #         ax1[pn].scatter(FTSdates[df], FTSAprvmrP[str(df)+'_'+str(pn)],facecolors='white', s=60, color='gray', label='a priori')
            #         ax1[pn].plot(FTSdates[df], FTSAprvmrP[str(df)+'_'+str(pn)], color='gray')
            
            #         ax1[pn].errorbar(FTSdates[df],FTSvmrP[str(df)+'_'+str(pn)],yerr=FTSvmrsdP[str(df)+'_'+str(pn)], markersize=60, linestyle='-', color='k', ecolor='k')
            #         ax1[pn].scatter(FTSdates[df], FTSvmrP[str(df)+'_'+str(pn)],facecolors='white', s=60, color='k', label='Retrieved')
                    
            #         ax1[pn].errorbar(sondedates[df], sondevmrP[str(df)+'_'+str(pn)],yerr=sondevmrsdP[str(df)+'_'+str(pn)],markersize=60, linestyle='-', color='r', ecolor='r')
            #         ax1[pn].scatter(sondedates[df], sondevmrP[str(df)+'_'+str(pn)],facecolors='white', s=60, color='r', label='Sonde')
            #         ax1[0].legend(prop={'size':12}, loc=1)

            #         ax2[pn].scatter(FTSdates[df], FTSAprsumP[str(df)+'_'+str(pn)],facecolors='white', s=60, color='gray', label='a priori')
            #         ax2[pn].plot(FTSdates[df], FTSAprsumP[str(df)+'_'+str(pn)], color='gray')
                    
            #         ax2[pn].errorbar(FTSdates[df],FTSsumP[str(df)+'_'+str(pn)],yerr=FTSsumsdP[str(df)+'_'+str(pn)],markersize=60, linestyle='-', color='k', ecolor='k', label='Retrieved')
            #         ax2[pn].scatter(FTSdates[df], FTSsumP[str(df)+'_'+str(pn)],facecolors='white', s=60, color='k')
            #         ax2[pn].errorbar(sondedates[df], sondesumP[str(df)+'_'+str(pn)],yerr=sondesumsdP[str(df)+'_'+str(pn)],markersize=60, linestyle='-', color='r', ecolor='r', label='Sonde')
            #         ax2[pn].scatter(sondedates[df], sondesumP[str(df)+'_'+str(pn)],facecolors='white', s=60, color='r')
            #         ax2[0].legend(prop={'size':12}, loc=1)


            #         ax1[pn].grid(True)
            #         ax1[pn].set_ylabel('Weighted VMR ['+sclname+']', fontsize=14)
            #         ax1[pn].set_title('Altitude Layer '+str(pcol[0])+' - '+str(pcol[1])+' km',
            #                   multialignment='center',fontsize=14)
            #         ax1[pn].tick_params(which='both',labelsize=14)
            #         ax1[pn].xaxis.set_major_formatter(DateFmt)
            #         ax1[pn].xaxis.set_minor_locator(monthLc)


            #         ax2[pn].grid(True)
            #         ax2[pn].set_ylabel('Partial Column\n[molecules cm$^{-2}$]',multialignment='center', fontsize=14)
            #         ax2[pn].set_title('Altitude Layer '+str(pcol[0])+' - '+str(pcol[1])+' km',
            #                   multialignment='center',fontsize=14)
            #         ax2[pn].tick_params(which='both',labelsize=14)
            #         ax2[pn].xaxis.set_major_formatter(DateFmt)
            #         ax2[pn].xaxis.set_minor_locator(monthLc)
                   
            #         ax1[pn].set_ylim(bottom=0)
            #         ax2[pn].set_ylim(bottom=0)
            #         fig1.autofmt_xdate()
            #         fig2.autofmt_xdate()

            #         fig1.subplots_adjust(bottom=0.075,top=0.95, left=0.12, right=0.95)
            #         fig2.subplots_adjust(bottom=0.075,top=0.95, left=0.12, right=0.95)

                    
            #         ax3[pn].errorbar(sondevmrP[str(df)+'_'+str(pn)], FTSvmrP[str(df)+'_'+str(pn)],xerr=sondevmrsdP[str(df)+'_'+str(pn)], yerr=FTSvmrsdP[str(df)+'_'+str(pn)], markersize=60, linestyle='none', color='k', ecolor='k')
            #         ax3[pn].scatter(sondevmrP[str(df)+'_'+str(pn)], FTSvmrP[str(df)+'_'+str(pn)], facecolors='white', s=60, color='k', label='Retrieved')
            #         ax3[pn].scatter(sondevmrP[str(df)+'_'+str(pn)], FTSAprvmrP[str(df)+'_'+str(pn)], facecolors='white', s=60, color='gray', label='a priori')            
            #         ax3[0].legend(prop={'size':12}, loc=1)

            #         # # o2o = np.arange(0, 1e4, 50)

            #         # # ax1.plot(o2o, o2o, linestyle='--',linewidth=2,  color='gray')
            #         # # ax1.plot(o2o*0.8, o2o, linestyle='--',linewidth=2,  color='gray')
            #         # # ax1.plot(o2o, o2o*0.8, linestyle='--',linewidth=2,  color='gray')

                    
                    
            #         ax3[pn].text(0.65,0.2,"slope: {0:.2f}".format(slope),transform=ax3[pn].transAxes, color='gray', fontsize=14)
            #         ax3[pn].text(0.65,0.13,"Intercept: {:.2f}".format(intercept),transform=ax3[pn].transAxes, color='gray', fontsize=14)
            #         ax3[pn].text(0.65,0.06,"r-value: {0:.2f}".format(r_value),transform=ax3[pn].transAxes, color='gray', fontsize=14)

            #         ax4[pn].errorbar(sondesumP[str(df)+'_'+str(pn)], FTSsumP[str(df)+'_'+str(pn)],yerr=FTSsumsdP[str(df)+'_'+str(pn)], xerr=sondesumsdP[str(df)+'_'+str(pn)], markersize=60, linestyle='none', color='k', ecolor='k')
            #         ax4[pn].scatter( sondesumP[str(df)+'_'+str(pn)], FTSsumP[str(df)+'_'+str(pn)], facecolors='white', s=60, color='k', label='Retrieved')
            #         ax4[pn].scatter( sondesumP[str(df)+'_'+str(pn)], FTSAprsumP[str(df)+'_'+str(pn)],facecolors='white', s=60, color='gray', label='a priori')

            #         ax4[pn].legend(prop={'size':12}, loc=1)

            #         # o2o = np.arange(0, 1e23, 1e22)

            #         # ax2.plot(o2o, o2o, linestyle='--',linewidth=2,  color='gray')
            #         # ax2.plot(o2o*0.8, o2o, linestyle='--',linewidth=2,  color='gray')
            #         # ax2.plot(o2o, o2o*0.8, linestyle='--',linewidth=2,  color='gray')


            #         slope, intercept, r_value, p_value, std_err = stats.linregress(sondesumP[str(df)+'_'+str(pn)], FTSsumP[str(df)+'_'+str(pn)])


            #         ax4[pn].text(0.025,0.92,"slope: {0:.2f}".format(slope),transform=ax4[pn].transAxes,  fontsize=14)
            #         ax4[pn].text(0.025,0.85,"Intercept: {:.2E}".format(intercept),transform=ax3[pn].transAxes,  fontsize=14)
            #         ax4[pn].text(0.025,0.78,"r-value: {0:.2f}".format(r_value),transform=ax4[pn].transAxes,  fontsize=14)

            #         slope, intercept, r_value, p_value, std_err = stats.linregress(sondesumP[str(df)+'_'+str(pn)], FTSAprsumP[str(df)+'_'+str(pn)])
                
            #         ax4[pn].text(0.65,0.2,"slope: {0:.2f}".format(slope),transform=ax4[pn].transAxes, color='gray', fontsize=14)
            #         ax4[pn].text(0.65,0.13,"Intercept: {:.2E}".format(intercept),transform=ax4[pn].transAxes, color='gray', fontsize=14)
            #         ax4[pn].text(0.65,0.06,"r-value: {0:.2f}".format(r_value),transform=ax4[pn].transAxes, color='gray', fontsize=14)

            #         ax3[pn].grid(True)
            #         ax3[pn].grid(True)
                
            #         ax3[pn].set_ylabel('Weighted VMR ['+sclname+']', fontsize=14)
            #         ax3[-1].set_xlabel('Weighted VMR - Sonde ['+sclname+']', fontsize=14)
            #         ax3[pn].tick_params(which='both',labelsize=14)
            #         ax3[pn].set_xlim(xmin=0)
            #         ax3[pn].set_ylim(ymin=0)
            #         ax3[pn].set_title('Altitude Layer '+str(pcol[0])+' - '+str(pcol[1])+' km', multialignment='center',fontsize=14)

            #         ax4[pn].grid(True)
            #         ax4[pn].grid(True)
                
            #         ax4[pn].set_ylabel('partial Column [molecules cm$^{-2}$]',multialignment='center', fontsize=14)
            #         ax4[-1].set_xlabel('Partial Column - sonde [molecules cm$^{-2}$]',multialignment='center', fontsize=14)
            #         ax4[pn].tick_params(which='both',labelsize=14)
            #         ax4[pn].set_xlim(xmin=0)
            #         ax4[pn].set_ylim(ymin=0)
            #         ax4[pn].set_title('Altitude Layer '+str(pcol[0])+' - '+str(pcol[1])+' km', multialignment='center',fontsize=14)


            #         fig3.subplots_adjust(bottom=0.075,top=0.95, left=0.15, right=0.95)
            #         fig4.subplots_adjust(bottom=0.075,top=0.95, left=0.15, right=0.95)

            # if self.pdfsav: 
            #     self.pdfsav.savefig(fig1,dpi=200)
            #     self.pdfsav.savefig(fig2,dpi=200)
            #     self.pdfsav.savefig(fig3,dpi=200)
            #     self.pdfsav.savefig(fig4,dpi=200)

            # else:
            #     plt.show(block=False)



        slope_all     =  np.asarray(slope_all)
        inter_all     =  np.asarray(inter_all)
        r_all         =  np.asarray(r_all)

        slopeApr_all  =  np.asarray(slopeApr_all)
        interApr_all  =  np.asarray(interApr_all)
        rApr_all      =  np.asarray(rApr_all)

        slope_odr     =  np.asarray(slope_odr)
        inter_odr     =  np.asarray(inter_odr)
        slopeErr_odr  =  np.asarray(slopeErr_odr)
        interErr_odr  =  np.asarray(interErr_odr)
       
        slopeApr_odr  =  np.asarray(slopeApr_odr)
        interApr_odr  =  np.asarray(interApr_odr)

        
        fig, (ax, ax2, ax3)  = plt.subplots(3,1,figsize=(7,10), sharex=True)

        ind = np.arange(len(slope_all))                # the x locations for the groups
        width = 0.35                      # the width of the bars
        rects1 = ax.bar(ind, slope_odr, width, color='green', label='Retrieved', yerr=slopeErr_odr, error_kw=dict(elinewidth=1.0, ecolor='k'))
        rects2 = ax.bar(ind+width, slopeApr_odr, width, color='blue', label='A priori')

        ax.set_xlim(-width,len(ind)+width*0.5)
        ax.set_ylabel('Slope', fontsize=14)
        ax.set_ylim(0,1.5)

        #xTickMarks = [str(pcol[0])+'-'+str(pcol[1])+ ' km' for pcol in pCols]
        ax.set_xticks(ind+width)
        #xtickNames = ax.set_xticklabels(xTickMarks)
        #plt.setp(xtickNames, rotation=45, fontsize=12)
        ax.tick_params(which='both',labelsize=14)
        ax.legend(prop={'size':12}, loc=1)

        rects1 = ax2.bar(ind, inter_odr/1e3, width, color='green', yerr=interErr_odr/1e3, error_kw=dict(elinewidth=1.0, ecolor='k'))
        rects2 = ax2.bar(ind+width, interApr_all/1e3, width, color='blue')

        ax2.set_xlim(-width,len(ind)+width*0.5)
        ax2.set_ylabel('Intercept [x10$^3$ ppm$_v$]', fontsize=14)
        #xTickMarks = [str(pcol[0])+'-'+str(pcol[1])+ ' km' for pcol in pCols]
        ax2.set_xticks(ind+width)
        #xtickNames = ax2.set_xticklabels(xTickMarks)
        #plt.setp(xtickNames, rotation=45, fontsize=12)
        ax2.tick_params(which='both',labelsize=14)

        rects1 = ax3.bar(ind, r_all, width, color='green')
        rects2 = ax3.bar(ind+width, rApr_all, width, color='blue')

        ax3.set_xlim(-width,len(ind)+width*0.5)
        ax3.set_ylabel('r-value', fontsize=14)
        xTickMarks = [str(pcol[0])+'-'+str(pcol[1]) for pcol in pCols]
        ax3.set_xticks(ind+width)
        xtickNames = ax3.set_xticklabels(xTickMarks)
        plt.setp(xtickNames, rotation=0, fontsize=11)
        ax3.tick_params(which='both',labelsize=14)
        ax3.set_xlabel('Layer [km]', fontsize=14)

        fig.subplots_adjust(bottom=0.075,top=0.95, left=0.15, right=0.95)
        
        if self.pdfsav: 
            self.pdfsav.savefig(fig,dpi=200)
            plt.savefig('/data/iortega/Manuscripts/Fig/StatsBars_'+self.loc.upper()+'.pdf', bbox_inches='tight')
        else:
            plt.show(block=False)


          #---------------------------------
            # Plot : Averaging Kernel Smoothing Function (row of avk)
            #---------------------------------
        clmap        = 'jet'
        cm           = plt.get_cmap(clmap)
        fig       = plt.figure(figsize=(7,7))
        gs        = gridspec.GridSpec(1,2,width_ratios=[3,1])
        ax        = plt.subplot(gs[0])
        axb       = plt.subplot(gs[1])
        cm        = plt.get_cmap(clmap)
        cNorm     = colors.Normalize(vmin=np.min(alt), vmax=np.max(alt))
        scalarMap = mplcm.ScalarMappable(norm=cNorm,cmap=clmap)
        scalarMap.set_array(alt)
        ax.set_color_cycle([scalarMap.to_rgba(x) for x in alt])
        
        for i in range(len(alt)):
            ax.plot(avkSCFav[i,:],alt)
            
        ax.set_ylabel('Altitude [km]', fontsize=12)
        ax.set_xlabel('Averaging Kernels', fontsize=12)
        ax.grid(True)
        cbar = fig.colorbar(scalarMap, orientation='vertical')
        cbar.set_label('Altitude [km]', fontsize=12)
        #ax.set_title('H2O Averaging Kernels Scale Factor', fontsize=14)
        ax.tick_params(labelsize=12)
        #ax.set_ylim(1, 15)
        if loc.lower() == 'fl0': ax.set_ylim(1, 15)
        if loc.lower() == 'mlo': ax.set_ylim(3, 15)
        
        axb.plot(np.sum(avkSCFav,axis=0), alt,color='k')
        axb.grid(True)
        axb.set_xlabel('Averaging Kernel Area', fontsize=12)
        axb.tick_params(axis='x',which='both',labelsize=12)
        major_ticks = np.arange(0, 3, 1)
        axb.set_xticks(major_ticks) 
        #axb.set_ylim(1, 15)
        if loc.lower() == 'fl0': axb.set_ylim(1, 15)
        if loc.lower() == 'mlo': axb.set_ylim(3, 15)
        #axb.tick_params(labelsize=14) 

        if self.pdfsav: 
            self.pdfsav.savefig(fig,dpi=200)
        else:
            plt.show(block=False)

        clmap        = 'jet'
        cm           = plt.get_cmap(clmap)
        fig       = plt.figure(figsize=(7,7))
        gs        = gridspec.GridSpec(1,2,width_ratios=[3,1])
        ax        = plt.subplot(gs[0])
        axb       = plt.subplot(gs[1])
        cm        = plt.get_cmap(clmap)
        cNorm     = colors.Normalize(vmin=np.min(alt), vmax=np.max(alt))
        scalarMap = mplcm.ScalarMappable(norm=cNorm,cmap=clmap)
        scalarMap.set_array(alt)
        ax.set_color_cycle([scalarMap.to_rgba(x) for x in alt])
        
        for i in range(len(alt)):
            ax.plot(avkVMRav[i,:],alt)
            
        ax.set_ylabel('Altitude [km]', fontsize=14)
        ax.set_xlabel('Averaging Kernels', fontsize=14)
        ax.grid(True)
        cbar = fig.colorbar(scalarMap, orientation='vertical')
        cbar.set_label('Altitude [km]', fontsize=14)
        ax.set_title('H2O Averaging Kernels [VMR]', fontsize=14)
        ax.tick_params(labelsize=14)
        
        axb.plot(np.sum(avkVMRav,axis=0), alt,color='k')
        axb.grid(True)
        axb.set_xlabel('Averaging Kernel Area', fontsize=14)
        axb.tick_params(axis='x',which='both',labelsize=8)
        #axb.tick_params(labelsize=14) 

        if self.pdfsav: 
            self.pdfsav.savefig(fig,dpi=200)
        else:
            plt.show(block=False)

            #---------------------------------

        


       
      
    