
#----------------------------------------------------------------------------------------
# Name:
#        MLSHDFClass.py
#
# Purpose:
#       This is a collection of classes and functions used for processing and ploting HDF files
#
#
# Notes:
#
#
# Version History:
#       Created, Sep, 2017  Ivan Ortega (iortega@ucar.edu)
#
#
#----------------------------------------------------------------------------------------

import datetime as dt
import time
import math
import sys
import numpy as np
import os
import csv

import os
from os import listdir
from os.path import isfile, join
import re

import matplotlib.animation as animation
import matplotlib
import matplotlib.dates as md
from matplotlib.dates import DateFormatter, MonthLocator, YearLocator
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import FormatStrFormatter, MultipleLocator,AutoMinorLocator,ScalarFormatter
from matplotlib.backends.backend_pdf import PdfPages #to save multiple pages in 1 pdf...
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as mplcm
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
import matplotlib.gridspec as gridspec
import glob

import myfunctions as mf

from netCDF4 import Dataset, Group
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA

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
 
                                        
class ReadNCMrg():

        def __init__(self, dataDir, iyear=False, fyear=False):


            if not( dataDir.endswith('/') ):
                    dataDir = dataDir + '/'

            self.dataDir    = dataDir

            self.Mrg        = {}
            self.ReadNCMrgFlg  = True

            if (not iyear) and (not fyear):
                print "Must specify initial and final years"
                exit()
 
            Nyears = fyear - iyear
            years  = [(iyear + i) for i in range(Nyears+1)]

            DirFiles  = []

            for y in years:
                DirFiles.append(glob.glob(dataDir + '*'+ str(y) + '*.nc4'))

            DirFiles = np.asarray(DirFiles)
            DirFiles = np.concatenate( DirFiles, axis=0 )
            DirFiles =  np.sort(DirFiles)

            for drs in DirFiles:

                #print '\nReading Nnc File: %s' % (drs)
                ckFile(drs)

                drsSplt  = drs.strip().split('/')[-1]
                yyyy     = int(drsSplt[26:30])

                try:
                    nc_fid    = Dataset(drs , 'r')  # Dataset is the class behavior to open the file
                    
                    c_attrs, nc_dims, nc_vars = mf.ncdump(nc_fid, verb=False)
                    vgrp =  str(nc_fid.groups)
                    #print vgrp
                    #print c_attrs
                    self.Mrg.setdefault('average',[]).append(np.asarray(nc_fid['/Merged/average']))
                    self.Mrg.setdefault('lat',[]).append(np.asarray(nc_fid['/Merged/lat']))
                    self.Mrg.setdefault('lev',[]).append(np.asarray(nc_fid['/Merged/lev']))
                    self.Mrg.setdefault('time',[]).append(np.asarray(nc_fid['/Merged/time']))

                    self.Mrg.setdefault('data_source',[]).append(np.asarray(nc_fid['/Merged/data_source']))
                    self.Mrg.setdefault('std_dev',[]).append(np.asarray(nc_fid['/Merged/std_dev']))
                    self.Mrg.setdefault('std_error',[]).append(np.asarray(nc_fid['/Merged/std_error']))
 
                except Exception as errmsg:
                    print errmsg
                    continue

            


class ReadAga():

    def __init__(self, dataDir):


        if not( dataDir.endswith('/') ):
                dataDir = dataDir + '/'

        self.dataDir    = dataDir
        self.AGA        = {}
        self.ReadAgaFlg    = True
        #-------------------------

        fileName = glob.glob(dataDir + '*.mon')

        for file in fileName:

            try:
    
                f = open(file, 'r')

                header1 = f.readline()
                header2 = f.readline()

                hsplit = header2.split(',')
                site   = hsplit[1]
   
                lat    = float(hsplit[4])
                lon    = float(hsplit[7])

                self.AGA.setdefault('site',[]).append(site)
                self.AGA.setdefault('lat',[]).append(lat)
                self.AGA.setdefault('lon',[]).append(lon)

                header3 = [f.readline() for i in range(14)]

                time          = []
                MM            = []
                YYYY          = []
                CFC11         = []
                CFC11_e       = []
                
                CFC11         = []
                CFC11_e       = []
                
                CFC12         = []
                CFC12_e       = []
                
                CFC113         = []
                CFC113_e       = []
                
                CCl4         = []
                CCl4_e       = []

                Dates          = []

                for line in f:
                    time        = line.strip()
                    columns     = line.split()

                    time        = float(columns[0])
                    MM          = int(columns[1])
                    YYYY        = int(columns[2])
                    Dates.append(dt.date(YYYY, MM, 15))       
                   
                    CFC11.append(float(columns[4]))
                    CFC11_e.append(float(columns[5]))

                    CFC12.append(float(columns[7]))
                    CFC12_e.append(float(columns[8]))

                    CCl4.append(float(columns[13]))
                    CCl4_e.append(float(columns[14]))

                    CFC113.append(float(columns[19]))
                    CFC113_e.append(float(columns[20]))
                

                self.AGA.setdefault(site+'CFC11',[]).append(CFC11)
                self.AGA.setdefault(site+'CFC11_e',[]).append(CFC11_e)

                self.AGA.setdefault(site+'CFC12',[]).append(CFC12)
                self.AGA.setdefault(site+'CFC12_e',[]).append(CFC12_e)

                self.AGA.setdefault(site+'CFC113',[]).append(CFC113)
                self.AGA.setdefault(site+'CFC113_e',[]).append(CFC113_e)

                self.AGA.setdefault(site+'CCl4',[]).append(CCl4)
                self.AGA.setdefault(site+'CCl4_e',[]).append(CCl4_e)
                
                self.AGA.setdefault(site+'Dates',[]).append(Dates)

            except Exception as errmsg:
                print '\nError in AGAGE file: {}'.format(file)

class ReadAga2():

    def __init__(self, dataDir, fileID=0):


        if not( dataDir.endswith('/') ):
                dataDir = dataDir + '/'

        self.dataDir    = dataDir
        self.AGA        = {}
        self.ReadAgaFlg    = True
        #-------------------------

        if fileID == 0: file = 'MHD-agA-AGAGE.mon'
        elif fileID == 1: file = 'total_cl_agage.out'
        else: 
            'Error: AGAGE File ID incorrect'
            sys.exit()

        fileName = dataDir + file

        ckDir(fileName, exit=True)

        if fileID == 0:

            f = open(fileName, 'r')

            header1 = f.readline()
            header2 = f.readline()

            hsplit = header2.strip().split(',')

            site   ='Mace Head'

            lat    = 53.33
            lon    = 9.9

            self.AGA.setdefault('site',[]).append(site)
            self.AGA.setdefault('lat',[]).append(lat)
            self.AGA.setdefault('lon',[]).append(lon)

            header3 = [f.readline() for i in range(13)]

            time          = []
            MM            = []
            YYYY          = []
            CFC11S         = []
            CFC11S_e       = []
            
            CFC11P         = []
            CFC11P_e       = []

            CFC11         = []
            CFC11_e       = []
            
            CFC12         = []
            CFC12_e       = []
            
            CFC113         = []
            CFC113_e       = []
            
            CCl4         = []
            CCl4_e       = []

            Dates          = []

            for line in f:
                time        = line.strip()
                columns     = line.split()

                time        = float(columns[0])
                MM          = int(columns[1])
                YYYY        = int(columns[2])
                Dates.append(dt.date(YYYY, MM, 15))       
               
                CFC11S.append(float(columns[4]))
                CFC11S_e.append(float(columns[5]))

                CFC11P.append(float(columns[7]))
                CFC11P_e.append(float(columns[8]))

                CFC12.append(float(columns[10]))
                CFC12_e.append(float(columns[11]))

                CCl4.append(float(columns[16]))
                CCl4_e.append(float(columns[17]))

                CFC113.append(float(columns[22]))
                CFC113_e.append(float(columns[23]))

               
                CFC11.append( (float(columns[4]) + float(columns[7])) / 2.)
                CFC11_e.append(np.sqrt( float(columns[5])**2 + float(columns[8])**2))
               
            
            self.AGA.setdefault(site+'CFC11',[]).append(CFC11)
            self.AGA.setdefault(site+'CFC11_e',[]).append(CFC11_e)

            self.AGA.setdefault(site+'CFC12',[]).append(CFC12)
            self.AGA.setdefault(site+'CFC12_e',[]).append(CFC12_e)

            self.AGA.setdefault(site+'CFC113',[]).append(CFC113)
            self.AGA.setdefault(site+'CFC113_e',[]).append(CFC113_e)

            self.AGA.setdefault(site+'CCl4',[]).append(CCl4)
            self.AGA.setdefault(site+'CCl4_e',[]).append(CCl4_e)
            
            self.AGA.setdefault(site+'Dates',[]).append(Dates)

        if fileID == 1:

            f = open(fileName, 'r')

            Dates          = []
            Cly            = []

            site   ='Global'

            for line in f:
                time        = line.strip()
                columns     = line.split()

                time        = float(columns[0])
                col2         = float(columns[1])

                t2dt        = mf.t2dt(time)

                Dates.append(t2dt)
                Cly.append(col2)


            self.AGA.setdefault(site+'Dates',[]).append(Dates)
            self.AGA.setdefault(site+'Cly',[]).append(Cly)

            self.AGA.setdefault('site',[]).append(site)
            self.AGA.setdefault('lat',[]).append(53.33)
            self.AGA.setdefault('lon',[]).append(9.99)





class PlotCly(ReadNCMrg, ReadAga2):

    def __init__(self,MrgdataDir, AgaDir, jfjFile, mloFile, iyear=False, fyear=False, saveFlg=False, outFname=''):
       
        #------------------------------------------------------------
        # If outFname is specified, plots will be saved to this file,
        # otherwise plots will be displayed to screen
        #------------------------------------------------------------
        if saveFlg:  self.pdfsav = PdfPages(outFname)
        else:        self.pdfsav = False

        self.ReadHDFMLSFlg = False
        self.ReadNCMrgFlg  = False
        self.ReadAgaFlg    = False

        self.AgagefileID   = 1

        self.jfjFile = jfjFile
        self.mloFile = mloFile

        #---------------
        # ReadOutputData
        #---------------
        ReadNCMrg.__init__(self,MrgdataDir, iyear=iyear, fyear=fyear)
        #ReadAga.__init__(self,AgaDir) 
        ReadAga2.__init__(self,AgaDir, fileID=self.AgagefileID)   
       
        
    def closeFig(self):
        self.pdfsav.close()

    def PlotClySet(self):
        #----------------------------------------
        #Define global interes
        #----------------------------------------
        loi       = 45            #Ltitude of interest
        poi       = 10.           #Ppressure of interest
        y4trend   = 1997          #initial year if trends


                                        #----------------------------------------
                                        #           GOZCARDS (MERGE OF HALOE, ACE, MLS)
                                        #----------------------------------------

        if self.ReadNCMrgFlg:
            
            time        = np.asarray(self.Mrg['time'])
            Prfs        = np.asarray(self.Mrg['average'])*1e9
            lev        = np.asarray(self.Mrg['lev'])
            Latitude    = np.asarray(self.Mrg['lat'])
            std_dev     = np.asarray(self.Mrg['std_dev'])*1e9

            Nfiles  = time.shape[0]
            NMonths = time.shape[1]
            Nlat    = Latitude.shape[1]
            Nlev    = lev.shape[1]

            print '\nN GOZCARDS files = ', Nfiles
            print 'N GOZCARDS Months = ', NMonths
            print 'N GOZCARDS Latitudes = ', Nlat
            print 'N GOZCARDS Levels = ', Nlev
            
            Dates = np.asarray([[dt.date(1950, 1,1) + dt.timedelta(days=int(t)) for t in tt] for tt in time])

            DatesMrg    = np.reshape(Dates, (Nfiles*NMonths))
            PrfsMrg     = np.reshape(Prfs, (Nfiles*NMonths, Nlev,Nlat ))
            PresMrg     = lev[0, :]
            LatMrg      = Latitude[0, :]
            StDMrg      = np.reshape(std_dev, (Nfiles*NMonths, Nlev,Nlat ))

            #----------------------------------------
            #DELETE POINTS IN THE pressure and latitude of interest (Important because some values are negatives)
            #----------------------------------------
            indMrg1   = mf.nearestind(poi, PresMrg)
            indLatMrg = mf.nearestind(loi, LatMrg)
            indsBad   =  np.where(PrfsMrg[:, indMrg1, indLatMrg] < 0.0)[0]

            print 'N points to remove in GOZCARDS = {}'.format(len(indsBad))
            PrfsMrg       = np.delete(PrfsMrg,indsBad,axis=0)
            StDMrg         = np.delete(StDMrg,indsBad,axis=0)
            DatesMrg       = np.delete(DatesMrg,indsBad,axis=0)


                                        #----------------------------------------
                                        #           Read AGAGE
                                        #----------------------------------------

        if self.ReadAgaFlg:
            #----------------------------------------
            #DEFINE VARIABLES
            #----------------------------------------
            Nsites = len(self.AGA['site'])
            s      = self.AGA['site']
            Lat   = np.asarray(self.AGA['lat'])

            indLatAGA = mf.nearestind(53., Lat)    ## This is to read MACE, HEAD
            ##indLatAGA = mf.nearestind(loi, Lat)  ## This is to read Closest to loi

            print '\nAGAGE site for Plot: {}, --> Latitude: {}'.format(self.AGA['site'][indLatAGA], Lat[indLatAGA]) 


            DatesAGA    = np.asarray(self.AGA[s[indLatAGA]+'Dates'][0])

            if self.AgagefileID == 0:

                CFC11    = np.asarray(self.AGA[s[indLatAGA]+'CFC11'][0])
                CFC11_e  = np.asarray(self.AGA[s[indLatAGA]+'CFC11_e'][0])

                CFC12    = np.asarray(self.AGA[s[indLatAGA]+'CFC12'][0])
                CFC12_e  = np.asarray(self.AGA[s[indLatAGA]+'CFC12_e'][0])

                CCl4     = np.asarray(self.AGA[s[indLatAGA]+'CCl4'][0])
                CCl4_e   = np.asarray(self.AGA[s[indLatAGA]+'CCl4_e'][0])

                CFC113   = np.asarray(self.AGA[s[indLatAGA]+'CFC113'][0])
                CFC113_e = np.asarray(self.AGA[s[indLatAGA]+'CFC113_e'][0])

                ind = np.where( (CFC11 > 0. )  & (CFC12 > 0.)  & (CCl4 > 0.))[0]

                Cl = CFC11[ind] + CFC12[ind] + CCl4[ind] 
                Cl_e = np.sqrt(CFC11_e[ind]**2 + CFC12_e[ind]**2 + CCl4_e[ind]**2)
                DatesAGA = DatesAGA[ind]

            elif self.AgagefileID == 1:

                Cl = np.asarray(self.AGA[s[indLatAGA]+'Cly'][0])



                                        #----------------------------------------
                                        #           Read JFJ
                                        #----------------------------------------

        f = open(self.jfjFile, 'r')

        HCl_Date_JFJ       = []
        HCl_JFJ           = []
        CLONO2_Date_JFJ    = []
        CLONO2_JFJ         = []

        Cl_Date_JFJ        = []
        Cl_JFJ             = []

        header1 = f.readline()

        for line in f:
            columns     = line.split()
            atime       = float(columns[0])

            da_i = mf.t2dt(atime)
            if da_i.year >= 1983:
                HCl_Date_JFJ.append(da_i)
                HCl_JFJ.append(float(columns[1])/1e15)

            if len(columns) >=3:

                atime       = float(columns[2])
                da_i = mf.t2dt(atime) 
                if da_i.year >= 1983:
                    CLONO2_Date_JFJ.append(da_i)
                    CLONO2_JFJ.append(float(columns[3])/1e15)

            if len(columns) >= 5:
                atime       = float(columns[4])
                da_i = mf.t2dt(atime)
                if da_i.year >= 1983:
                    Cl_Date_JFJ.append(da_i)
                    Cl_JFJ.append(float(columns[5])/1e15)


        HCl_Date_JFJ = np.asarray(HCl_Date_JFJ)
        HCl_JFJ      = np.asarray(HCl_JFJ)

        Cl_JFJ       = np.asarray(Cl_JFJ)
        Cl_Date_JFJ  = np.asarray(Cl_Date_JFJ)


                                        #----------------------------------------
                                        #           Read MLO
                                        #----------------------------------------

        f = open(self.mloFile, 'r')

        HCl_Date_MLO       = []
        HCl_MLO            = []

        headers = [f.readline() for i in range(7)]
        
        for line in f:
    
            columns     = line.split(',')

            yyyy = int(columns[1].strip()[0:4])
            mm   = int(columns[1].strip()[5:7])
            dd   = int(columns[1].strip()[8:10])

            hh   = int(columns[2].strip()[0:2])
            mi   = int(columns[2].strip()[3:5])
            se   = int(columns[2].strip()[6:8])

            da = dt.datetime(yyyy, mm, dd, hh, mi, se)
        
            HCl_Date_MLO.append(da)
            HCl_MLO.append(float(columns[3].strip()))


        HCl_MLO      = np.asarray(HCl_MLO)/1e15
        HCl_Date_MLO = np.asarray(HCl_Date_MLO)

        Avg = mf.mnthlyAvg(HCl_MLO, HCl_Date_MLO, dateAxis=1)
        HCl_Date_MLO        = Avg['dates']
        HCl_MLO             = Avg['mnthlyAvg']
        HCl_std_MLO         = Avg['std']

        #----------------------------------------------------------------------------
        #
        #----------------------------------------------------------------------------
        years = [ singDate.year for singDate in DatesMrg]      # Find years for all date entries
        if len(list(set(years))) > 1: yrsFlg = True            # Determine all unique years
        else:                         yrsFlg = False   

        yearsLc      = YearLocator(2)
        monthsAll    = MonthLocator()
        #months       = MonthLocator()
        months       = MonthLocator(bymonth=1, bymonthday=1)
        if yrsFlg: DateFmt   = DateFormatter('%Y')
        else: DateFmt        = DateFormatter('%m\n%Y')

        #----------------------------------------------------------------------------
        #
        #----------------------------------------------------------------------------
        indLatMrg = mf.nearestind(loi, LatMrg)

                                                #-----------------------------------------
                                                #             PLOT (ONE PANEL, THREE Y-AXIS)
                                                #-----------------------------------------
        fig, ax = plt.subplots(figsize=(12, 7), sharex=True)
       
        axes = [ax, ax.twinx(), ax.twinx()]

        axes[-1].spines['right'].set_position(('axes', 1.125))
        axes[-1].set_frame_on(True)
        axes[-1].patch.set_visible(False)

        #-----------------------------------------
        #PLOT GOZCARDS
        #-----------------------------------------
        if self.ReadHDFMLSFlg: ind1MLS = mf.nearestind(poi, PresMLS)
        if self.ReadNCMrgFlg:  indMrg1 = mf.nearestind(poi, PresMrg)

        indsGood =  np.where(PrfsMrg[:, indMrg1, indLatMrg] > 0.0)[0]

        PrfsMrg2        = PrfsMrg[indsGood, indMrg1, indLatMrg] 
        StDMrg2         = StDMrg[indsGood, indMrg1, indLatMrg]  
        DatesMrg2       = DatesMrg[indsGood] 
        
        PrfsMrg2_smth     = mf.smooth(PrfsMrg2,  window_len=24, window='flat')
        axes[0].plot(DatesMrg2,PrfsMrg2_smth,  linewidth=5, color='gray')

        #-----------------------------------------
        #GOZCARDS TRENDS ANALYSIS
        #-----------------------------------------
        yearGO = [single.year for single in DatesMrg2]
        yearGO = np.asarray(yearGO)

        inds = np.where( yearGO >= y4trend )[0]

        dateYearFrac = mf.toYearFraction(DatesMrg2[inds])
        weights      = np.ones_like(dateYearFrac)
        res          = mf.fit_driftfourier(dateYearFrac, PrfsMrg2[inds], weights, 2)
        f_drift, f_fourier, f_driftfourier = res[3:6]

        res_b        = mf.cf_driftfourier(dateYearFrac, PrfsMrg2[inds], weights, 2)
        perc, intercept_b, slope_b, pfourier_b = res_b

        print "\nFitted trend GOZCARDS -- slope: {0:.3E} ({1:.3f}%)".format(res[1],res[1]/np.nanmean(PrfsMrg2[inds])*100.0)
        axes[0].scatter(DatesMrg2,PrfsMrg2, s=20, color='gray', edgecolors='gray', alpha=0.75, label='GOZCARDS (10hPa, 40-50$^\circ$ N)\n {0:.2f} $\pm$ {1:.2f} %/yr'.format(res[1]/np.nanmean(PrfsMrg2[inds])*100.0, np.std(slope_b)/np.mean(PrfsMrg2[inds])*100.0))   #facecolors='none'   
        axes[0].tick_params(which='both',labelsize=14)
        axes[0].set_ylim(1.7, 3.4)
        axes[0].set_ylabel('VMR [ppb], GOZCARDS (HALOE, ACE & MLS) (HCl)', fontsize=14)

        if yrsFlg:
            axes[0].xaxis.set_major_locator(yearsLc)
            axes[0].xaxis.set_minor_locator(months)
            axes[0].xaxis.set_major_formatter(DateFmt)

        else:
            axes[0].xaxis.set_major_locator(monthsAll)
            axes[0].xaxis.set_major_formatter(DateFmt)

        # #-----------------------------------------
        # #SAVE GOZCARDS DATA IN ASCII
        # #----------------------------------------- 
        # YYYYMMDD = np.asarray(['{0:4d}-{1:02d}-{2:02d}'.format(DatesMrg2[i].year, DatesMrg2[i].month, DatesMrg2[i].day)    for i,dum in enumerate(DatesMrg2)])

        # with open('/data1/Campaign/Satellite/MLS/GOZCARD_pltCly.ascii','w') as fopen:
        #     fopen.write('#Hannigan, J.W., Ortega, I\n')
        #     fopen.write('#National Center for Atmospheric Research\n')
        #     fopen.write('#GOZCARDS Monthly data (GOZCARDS (10hPa, 40-50deg N)\n')
        #     fopen.write('#GOZCARDS smoothed data is included using a window of 24 (2 years)\n')
        #     fopen.write('#CONTACT_INFO: Hannigan, Jim, jamesw@ucar.edu, 303-497-1853, NCAR, 3090 Center Green Drive, Boulder, CO 80301\n')
        #     fopen.write('#CONTACT_INFO: Ortega, Ivan, iortega@ucar.edu, 303-497-1861, NCAR, 3090 Center Green Drive, Boulder, CO 80301\n')
  
        #     fopen.write('Index, YYYY-MM-DD, HCl [ppb], HCl_smooth [ppb]\n')
        #     strFormat = '{0:d}, {1:>10s},  {2:.3f}, {3:.3f}\n'

        #     for i,sngTime in enumerate(YYYYMMDD):
        #         fopen.write(strFormat.format((i+1),YYYYMMDD[i], PrfsMrg2[i], PrfsMrg2_smth[i]))

        #-----------------------------------------
        #PLOT AGAGE
        #----------------------------------------- 
        if self.AgagefileID == 0: axes[1].set_ylabel('VMR [ppt], AGAGE (CFC11 + CFC12 + CCl$_4$)', color='green', fontsize=14)
        if self.AgagefileID == 1: axes[1].set_ylabel('Mole Fraction [ppb], AGAGE (Total Cl)', color='green', fontsize=14)
        axes[1].tick_params(axis='y', colors='green', labelsize=14)
        if self.AgagefileID == 0: axes[1].set_ylim(380, 1100)
        if self.AgagefileID == 1: axes[1].set_ylim(1, 4.5)

        yearAGA = [single.year for single in DatesAGA]
        yearAGA = np.asarray(yearAGA)

        #-----------------------------------------
        #AGAGE TRENDS ANALYSIS
        #-----------------------------------------
        inds = np.where( yearAGA >= y4trend )[0]

        dateYearFrac = mf.toYearFraction(DatesAGA[inds])
        weights      = np.ones_like(dateYearFrac)
        res          = mf.fit_driftfourier(dateYearFrac, Cl[inds], weights, 2)
        f_drift, f_fourier, f_driftfourier = res[3:6]

        res_b        = mf.cf_driftfourier(dateYearFrac, Cl[inds], weights, 2)
        perc, intercept_b, slope_b, pfourier_b = res_b

        print "Fitted trend AGAGE -- slope: {0:.3E} ({1:.3f}%)".format(res[1],res[1]/np.mean(Cl[inds])*100.0)
        if self.AgagefileID == 0: axes[1].scatter(DatesAGA,Cl, s=20, color='green', edgecolors='green', alpha=0.75, label='AGAGE Mace Head, Ireland (53.33$^\circ$ N)\n {0:.2f} $\pm$ {1:.2f} %/yr'.format(res[1]/np.nanmean(Cl[inds])*100.0, np.std(slope_b)/np.mean(Cl[inds])*100.0))   #facecolors='none'
        if self.AgagefileID == 1: axes[1].scatter(DatesAGA,Cl, s=20, color='green', edgecolors='green', alpha=0.75, label='AGAGE Global Mean Tropospheric Total Cl\n {0:.2f} $\pm$ {1:.2f} %/yr'.format(res[1]/np.nanmean(Cl[inds])*100.0, np.std(slope_b)/np.mean(Cl[inds])*100.0))   #facecolors='none'

       # #-----------------------------------------
       # #SAVE AGAGE DATA IN ASCII
       # #----------------------------------------- 
       #  YYYYMMDD = np.asarray(['{0:4d}-{1:02d}-{2:02d}'.format(DatesAGA[i].year, DatesAGA[i].month, DatesAGA[i].day)    for i,dum in enumerate(DatesAGA)])

       #  with open('/data1/Campaign/Satellite/MLS/AGAGE_pltCly.ascii','w') as fopen:
       #      fopen.write('#Hannigan, J.W., Ortega, I\n')
       #      fopen.write('#National Center for Atmospheric Research\n')
       #      fopen.write('#AGAGE  Monthly data from Mace Head, Ireland (53.33deg N)\n')
       #      fopen.write('#AGAGE  Cly =  CFC11 + CFC12 + CCl4\n') 
       #      fopen.write('#CONTACT_INFO: Hannigan, Jim, jamesw@ucar.edu, 303-497-1853, NCAR, 3090 Center Green Drive, Boulder, CO 80301\n')
       #      fopen.write('#CONTACT_INFO: Ortega, Ivan, iortega@ucar.edu, 303-497-1861, NCAR, 3090 Center Green Drive, Boulder, CO 80301\n')
       
       #      fopen.write('Index, YYYY-MM-DD, Cly [ppm]\n')
       #      strFormat = '{0:d}, {1:>10s},  {2:.3f}\n'

       #      for i,sngTime in enumerate(YYYYMMDD):
       #          fopen.write(strFormat.format((i+1),YYYYMMDD[i], Cl[i])) 

        #-----------------------------------------
        #PLT JFJ HRFTIR
        #-----------------------------------------
        
        HCl_JFJ_smth     = mf.smooth(HCl_JFJ,  window_len=24, window='flat')
        axes[2].plot(HCl_Date_JFJ,HCl_JFJ_smth,  linewidth=5, color='blue')

        yearJFJ = [single.year for single in HCl_Date_JFJ]
        yearJFJ = np.asarray(yearJFJ)

        inds = np.where( yearJFJ >= y4trend )[0]

        dateYearFrac = mf.toYearFraction(HCl_Date_JFJ[inds])
        weights      = np.ones_like(dateYearFrac)
        res          = mf.fit_driftfourier(dateYearFrac, HCl_JFJ[inds], weights, 2)
        f_drift, f_fourier, f_driftfourier = res[3:6]

        res_b        = mf.cf_driftfourier(dateYearFrac, HCl_JFJ[inds], weights, 2)
        perc, intercept_b, slope_b, pfourier_b = res_b
  

        print "Fitted trend JFJ -- slope: {0:.3E} ({1:.3f}%)".format(res[1],res[1]/np.mean(HCl_JFJ[inds])*100.0)
        axes[2].scatter(HCl_Date_JFJ,HCl_JFJ, s=20, color='blue', edgecolors='blue', alpha=0.75, label='Jungfraujoch NDACC FTIR (46.55$^\circ$ N)\n {0:.2f} $\pm$ {1:.2f} %/yr'.format(res[1]/np.nanmean(HCl_JFJ[inds])*100.0, np.std(slope_b)/np.mean(HCl_JFJ[inds])*100.0))
        
        # #-----------------------------------------
        # #SAVE JFJ DATA IN ASCII
        # #----------------------------------------- 
        # YYYYMMDD = np.asarray(['{0:4d}-{1:02d}-{2:02d}'.format(HCl_Date_JFJ[i].year, HCl_Date_JFJ[i].month, HCl_Date_JFJ[i].day)    for i,dum in enumerate(HCl_Date_JFJ)])

        # with open('/data1/Campaign/Satellite/MLS/JFJ_pltCly.ascii','w') as fopen:
        #     fopen.write('#Hannigan, J.W., Ortega, I\n')
        #     fopen.write('#National Center for Atmospheric Research\n')
        #     fopen.write('#Jungfraujoch NDACC FTIR Monthly Data (46.55deg N)\n')
        #     fopen.write('#Jungfraujoch smoothed data is included using a window of 24 (2 years)\n')
        #     fopen.write('#CONTACT_INFO: Hannigan, Jim, jamesw@ucar.edu, 303-497-1853, NCAR, 3090 Center Green Drive, Boulder, CO 80301\n')
        #     fopen.write('#CONTACT_INFO: Ortega, Ivan, iortega@ucar.edu, 303-497-1861, NCAR, 3090 Center Green Drive, Boulder, CO 80301\n')
  
        #     fopen.write('Index, YYYY-MM-DD, HCl [[x10^15 molec/cm^2]], HCl_smooth [[x10$^15 molec/cm^2]]\n')
        #     strFormat = '{0:d}, {1:>10s},  {2:.3f}, {3:.3f}\n'

        #     for i,sngTime in enumerate(YYYYMMDD):
        #         fopen.write(strFormat.format((i+1),YYYYMMDD[i], HCl_JFJ[i], HCl_JFJ_smth[i]))

        #-----------------------------------------
        #PLT MLO HRFTIR
        #-----------------------------------------
        
        HCl_MLO_smth     = mf.smooth(HCl_MLO,  window_len=24, window='flat')
        axes[2].plot(HCl_Date_MLO,HCl_MLO_smth,  linewidth=5, color='navy')

        yearMLO = [single.year for single in HCl_Date_MLO]
        yearMLO = np.asarray(yearMLO)

        inds = np.where( yearMLO >= y4trend )[0]

        dateYearFrac = mf.toYearFraction(HCl_Date_MLO[inds])
        weights      = np.ones_like(dateYearFrac)
        res          = mf.fit_driftfourier(dateYearFrac, HCl_MLO[inds], weights, 2)
        f_drift, f_fourier, f_driftfourier = res[3:6]

        res_b        = mf.cf_driftfourier(dateYearFrac, HCl_JFJ[inds], weights, 2)
        perc, intercept_b, slope_b, pfourier_b = res_b

        print "Fitted trend MLO -- slope: {0:.3E} ({1:.3f}%)".format(res[1],res[1]/np.mean(HCl_MLO[inds])*100.0)

        axes[2].scatter(HCl_Date_MLO,HCl_MLO, s=20, facecolors='none', edgecolors='navy', alpha=0.75, label='Mauna Loa Observatory NDAC FTIR (19.54$^\circ$ N)\n {0:.2f} $\pm$ {1:.2f} %/yr'.format(res[1]/np.nanmean(HCl_MLO[inds])*100.0, np.std(slope_b)/np.mean(HCl_MLO[inds])*100.0))

        axes[2].set_ylim(1.2, 5.3)
        axes[2].set_xlim(dt.date(1979, 1, 1), dt.date(2017, 12, 31))
       
        axes[2].set_ylabel('Total Column [x10$^{15}$ molec/cm$^2$], FTIR (HCl)', color='blue', fontsize=14)
        axes[2].tick_params(axis='y', colors='blue', labelsize=14)

        axes[0].set_xlabel('Year', fontsize=14)

        fig.autofmt_xdate()
        fig.subplots_adjust(left=0.1, bottom=0.125, right=0.8, top=0.94)

        lines, labels   = axes[0].get_legend_handles_labels()
        lines2, labels2 = axes[1].get_legend_handles_labels()
        lines3, labels3 = axes[2].get_legend_handles_labels()
       
        axes[0].legend(lines + lines2 + lines3, labels + labels2 + labels3, prop={'size':10.5}, loc=2, frameon=False, ncol=2)  #'weight':'bold'

        plt.suptitle('NH Inorganic Chlorine 1983 to 2016 (Trends 1997 to 2016)', fontsize=16  )

        # #-----------------------------------------
        # #SAVE MLO DATA IN ASCII
        # #----------------------------------------- 
        # YYYYMMDD = np.asarray(['{0:4d}-{1:02d}-{2:02d}'.format(HCl_Date_MLO[i].year, HCl_Date_MLO[i].month, HCl_Date_MLO[i].day)    for i,dum in enumerate(HCl_Date_MLO)])

        # with open('/data1/Campaign/Satellite/MLS/MLO_pltCly.ascii','w') as fopen:
        #     fopen.write('#Hannigan, J.W., Ortega, I\n')
        #     fopen.write('#National Center for Atmospheric Research\n')
        #     fopen.write('#Mauna Loa Observatory NDAC FTIR Monthly Data (19.54deg N)\n')
        #     fopen.write('#Mauna Loa smoothed data is included using a window of 24 (2 years)\n')
        #     fopen.write('#CONTACT_INFO: Hannigan, Jim, jamesw@ucar.edu, 303-497-1853, NCAR, 3090 Center Green Drive, Boulder, CO 80301\n')
        #     fopen.write('#CONTACT_INFO: Ortega, Ivan, iortega@ucar.edu, 303-497-1861, NCAR, 3090 Center Green Drive, Boulder, CO 80301\n')
  
        #     fopen.write('Index, YYYY-MM-DD, HCl [[x10^15 molec/cm^2]], HCl_smooth [[x10$^15 molec/cm^2]]\n')
        #     strFormat = '{0:d}, {1:>10s},  {2:.3f}, {3:.3f}\n'

        #     for i,sngTime in enumerate(YYYYMMDD):
        #         fopen.write(strFormat.format((i+1),YYYYMMDD[i], HCl_MLO[i], HCl_MLO_smth[i]))

        if self.pdfsav: 
            self.pdfsav.savefig(fig,dpi=200)
            if self.AgagefileID == 0: plt.savefig('/data1/Campaign/Satellite/MLS/Cly_all_Mace.pdf') # bbox_inches='tight'
            if self.AgagefileID == 1: plt.savefig('/data1/Campaign/Satellite/MLS/Cly_all_Global.pdf') # bbox_inches='tight'
        else:
            plt.draw()
            plt.show(block=False)
            user_input = raw_input('Press any key to exit >>> ')
            sys.exit()



