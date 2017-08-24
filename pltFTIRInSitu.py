#! /usr/bin/python2.7
#------------------------------------------------------
#SCRIPT TO PLOT IN-SITU Vs FTIR
#------------------------------------------------------

import sys
import datetime as dt
import matplotlib
from numpy import fromstring, vectorize, ndarray, array, genfromtxt
import matplotlib.dates as md
from matplotlib.dates import DateFormatter, MonthLocator, YearLocator, DayLocator, WeekdayLocator, MONDAY
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import FormatStrFormatter, MultipleLocator,AutoMinorLocator,ScalarFormatter
from matplotlib.backends.backend_pdf import PdfPages #to save multiple pages in 1 pdf...
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as mplcm
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
import matplotlib.gridspec as gridspec
import numpy as np

from itertools import izip
import myfunctions as mf
import dataModelOutClass as dm
import dataOutClass as dc
from collections                     import OrderedDict
from scipy import linspace, polyval, polyfit, sqrt, stats, randn

from scipy.io import netcdf
import os
import glob

from scipy import interpolate

import matplotlib.dates as md
import matplotlib.dates as md
from matplotlib.dates import DateFormatter, MonthLocator, YearLocator, DayLocator, WeekdayLocator, MONDAY
from scipy.interpolate import InterpolatedUnivariateSpline as intrpUniSpl

import pytz

def month_string_to_number(string):
    m = {
        'jan': 1,
        'feb': 2,
        'mar': 3,
        'apr':4,
         'may':5,
         'jun':6,
         'jul':7,
         'aug':8,
         'sep':9,
         'oct':10,
         'nov':11,
         'dec':12
        }
    s = string.strip()[:3].lower()

    try:
        out = m[s]
        return out
    except:
        raise ValueError('Not a month')

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

def closeFig(self):
    self.pdfsav.close()


def main():

    #-----------------------------------------------------------------------------------------
    #                 Initializations for FTIR
    #-----------------------------------------------------------------------------------------

    loc               = 'mlo'
    gasName           = ['h2o', 'hdo']
    ver               = ['Current_ERA_l', 'Current_l'] 
    ctlF              = ['sfit4.ctl', 'sfit4.ctl'] 
    maxRMS            = [0.5, 0.5]                      
    minDOF            = [1.0, 1.0]

    saveFlg           = True 
    pltFile           =  '/data/iortega/results/'+loc.lower()+'/fig/FTS-inSitu_'+loc.upper()+'_log.pdf'               

    #------
    # Flags
    #------
    errorFlg           = False                   # Flag to process error data
    fltrFlg            = True                   # Flag to filter the data
    byYrFlg            = False                  # Flag to create plots for each individual year in date range
    szaFlg             = True                   # Flag to filter based on min and max SZA
    dofFlg             = True                   # Flag to filter based on min DOFs
    pcNegFlg           = True                  # Flag to filter profiles with negative partial columns
    tcNegFlg           = True                  # Flagsag to filter profiles with negative total columns
    tcMMFlg            = True                  # Flag to filter based on min and max total column amount
    cnvrgFlg           = True                   # Flag to filter profiles that did not converge
    rmsFlg             = True                   # Flag to filter based on max RMS
    chiFlg             = False                  # Flag to filter based on max CHI_2_Y
    mnthFlg            = False                  # Flag to filter based on 
    mnths              = [6,7,8]                # Months to filter on (these are the months to include data)
    minSZA             = 0.0                    # Min SZA for filtering
    maxSZA             = 90.0                   # Max SZA for filtering
    maxCHI             = 2.0                    # Max CHI_y_2 value
    maxTC              = 5.0E24                 # Max Total column amount for filtering
    minTC              = 0.0                    # Min Total column amount for filtering
    sclfct             = 1.0E6                  # Scale factor to apply to vmr plots (ppmv=1.0E6, ppbv=1.0E9, etc)
    sclfctName         = 'ppmv'                 # Name of scale factor for labeling plots

    pCols      = [ [3.0, 5.0], [5.0, 7.0], [7.0, 9.0] ]



    #-------------------------------------------------------
    #READ IN-SITU NDATA - H2O (ppm)
    #-------------------------------------------------------
    File = '/data1/ancillary_data/mlo/mlo_vapor_ppmv.dat'    

    utc = pytz.utc
    Hawaii = pytz.timezone('US/Hawaii')

    with open(File ,'r') as fopen: lines = fopen.readlines()

    Data = [row.strip().split() for row in lines]

    Date_in   = []
    insitu   = []

    for i, row in enumerate(Data):

        ho   =  int(row[1][0:2].strip())
        mi   =  int(row[1][3:5].strip())
        ss   =  int(row[1][6:8].strip())

        #if (ho <= 7) or (ho <= 10) or (ho >= 14):
        
        insitu.append(float(row[2]))
        yyyy = int(row[0][7:11].strip())
        smm  = str(row[0][3:6].strip())
        mm   = month_string_to_number(smm)
        dd   = int(row[0][0:2].strip())

        date_lt  = dt.datetime(yyyy, mm, dd, ho, mi, ss)
        date_lt  = Hawaii.localize(date_lt)

        Date_in.append(date_lt.astimezone(utc))

    insitu    = np.asarray(insitu)
    Date_in   = np.asarray(Date_in)

    #-------------------------------------------------------
    #READ IN-SITU NDATA - dD
    #-------------------------------------------------------
    File = '/data1/ancillary_data/mlo/dD_LGR_MLO.txt'

    with open(File ,'r') as fopen: lines = fopen.readlines()

    Data = [row.strip().split() for row in lines]

    Date_in2   = []
    insitu2    = []

    for i, row in enumerate(Data):

        ho   =  int(row[1][0:2].strip())
        mi   =  int(row[1][3:5].strip())
        ss   =  int(row[1][6:8].strip())

        insitu2.append(float(row[2]))
        yyyy = int(row[0][7:11].strip())
        smm  = str(row[0][3:6].strip())
        mm   = month_string_to_number(smm)
        dd   = int(row[0][0:2].strip())

        date_lt  = dt.datetime(yyyy, mm, dd, ho, mi, ss)
        date_lt  = Hawaii.localize(date_lt)

        Date_in2.append(date_lt.astimezone(utc))

    dDinsitu    = np.asarray(insitu2)
    #Date_in2   = np.asarray(Date_in2)

    #-------------------------------------------------------
    #READ IN-SITU NDATA - HDO
    #-------------------------------------------------------
    File = '/data1/ancillary_data/mlo/mlo_hdo.dat'

    with open(File ,'r') as fopen: lines = fopen.readlines()

    Data = [row.strip().split() for row in lines]

    Date_in3   = []
    insitu3    = []

    for i, row in enumerate(Data):

        ho   =  int(row[1][0:2].strip())
        mi   =  int(row[1][3:5].strip())
        ss   =  int(row[1][6:8].strip())

        insitu3.append(float(row[2]))
        yyyy = int(row[0][7:11].strip())
        smm  = str(row[0][3:6].strip())
        mm   = month_string_to_number(smm)
        dd   = int(row[0][0:2].strip())

        date_lt  = dt.datetime(yyyy, mm, dd, ho, mi, ss)
        date_lt  = Hawaii.localize(date_lt)

        Date_in3.append(date_lt.astimezone(utc))

    HDOinsitu    = np.asarray(insitu3)*2.0   #A CORRECTION OF 2 IS NEEDED

    dDinsitu_r         = (np.divide(HDOinsitu, insitu)/3.1152e-4 - 1.0) *1000.

    #----------------------
    # Date range to process FTIR (depending on insitu dates)
    #----------------------
    iyear              = Date_in[0].year
    imnth              = Date_in[0].month
    iday               = Date_in[0].day
    fyear              = Date_in[-1].year
    fmnth              = Date_in[-1].month
    fday               = Date_in[-1].day
    #fyear              = 2016
    #fmnth              = 6
    #fday               = 30

    #-------------------------------------------------
    #FTS
    #-------------------------------------------------
    retDir = [ '/data1/ebaumer/'+loc.lower()+'/'+g.lower()+'/'+v+'/' for g,v in izip(gasName,ver)] 
    ctlFile  = ['/data1/ebaumer/'+loc.lower()+'/'+g.lower()+'/'+'x.'+g.lower()+'/'+ ctlF[i] for i, g in enumerate(gasName)]

    #---------------------------
    # Check file and directories
    #---------------------------
    for d in retDir:  ckDir(d,exit=True)
    for c in ctlFile: ckFile(c,exit=True)
    ckDir(os.path.dirname(os.path.realpath(pltFile)),exit=True)

    #-------------------------------------
    # Create instance of output data class   
    #-------------------------------------
    statDataCl = OrderedDict()
    for i,gas in enumerate(gasName):
        statDataCl[gas+'_'+ver[i]] = dc.ReadOutputData(retDir[i],'',ctlFile[i],iyear,imnth,iday,fyear,fmnth,fday)

    #--------------
    # Read profiles
    #--------------
    for g, gasVer in enumerate(statDataCl):
        if gasName[g].lower() == 'h2o':
            statDataCl[gasVer].readprfs([statDataCl[gasVer].PrimaryGas, 'HDO'],retapFlg=1)
            statDataCl[gasVer].readprfs([statDataCl[gasVer].PrimaryGas, 'HDO'],retapFlg=0)
        else: 
            statDataCl[gasVer].readprfs([statDataCl[gasVer].PrimaryGas],retapFlg=1)
            statDataCl[gasVer].readprfs([statDataCl[gasVer].PrimaryGas],retapFlg=0)

        if statDataCl[gasVer].empty: 
            print 'No retreivals found for {}. Exiting.......'.format(gasVer)
            sys.exit()

    rPrfVMR      = OrderedDict()
    rPrfMol      = OrderedDict()
    dates        = OrderedDict()
    alt          = OrderedDict()
    Airmass      = OrderedDict()
    waterVMR     = OrderedDict()
    waterMol     = OrderedDict()
    totClmn      = OrderedDict()
    totClmn_e    = OrderedDict()
    TCdryAir     = OrderedDict()
    TCdry        = OrderedDict()
    rPrfDry      = OrderedDict()
    rms          = OrderedDict()
    dofsAvg      = OrderedDict()
    dofsAvg_cs   = OrderedDict()
    LayVMR       = OrderedDict()
    DryLayVMR    = OrderedDict()
    upperHgt     = OrderedDict()
    lowerHgt     = OrderedDict()
    LayThk       = OrderedDict()
    LayDOF       = OrderedDict()
    avkSCF       = OrderedDict()
    avkVMR       = OrderedDict()
    avkSCFav     = OrderedDict()
    avkVMRav     = OrderedDict()
    aPrfVMR      = OrderedDict()
    aPrfMol      = OrderedDict()
    tot_err      = OrderedDict()      
    rand_err     = OrderedDict()
    sys_err      = OrderedDict()
    tot_rnd      = OrderedDict()
    tot_sys      = OrderedDict()
    totClmn_e    = OrderedDict()

    vmrP         = OrderedDict()
    TCp         = OrderedDict()


    for j, gasVer in enumerate(statDataCl):
        rPrfVMR[gasVer]  = np.asarray(statDataCl[gasVer].rprfs[statDataCl[gasVer].PrimaryGas]) * sclfct
        rPrfMol[gasVer]  = np.asarray(statDataCl[gasVer].rprfs[statDataCl[gasVer].PrimaryGas]  * np.asarray(statDataCl[gasVer].rprfs['AIRMASS']))
        dates[gasVer]    = statDataCl[gasVer].rprfs['date']
        alt[gasVer]      = np.asarray(statDataCl[gasVer].rprfs['Z'][0,:])
        Airmass[gasVer]  = np.asarray(statDataCl[gasVer].rprfs['AIRMASS'])
        waterVMR[gasVer] = np.asarray(statDataCl[gasVer].aprfs['H2O'])
        waterMol[gasVer] = np.asarray(statDataCl[gasVer].aprfs['H2O'] * Airmass[gasVer])
        totClmn[gasVer]  = np.sum(rPrfMol[gasVer],axis=1)
        TCdryAir[gasVer] = np.sum(Airmass[gasVer],axis=1) - np.sum(waterMol[gasVer],axis=1)
        TCdry[gasVer]    = (totClmn[gasVer] / TCdryAir[gasVer]) * sclfct
        aPrfVMR[gasVer]  = np.asarray(statDataCl[gasVer].aprfs[statDataCl[gasVer].PrimaryGas]) * sclfct
        aPrfMol[gasVer]  = np.asarray(statDataCl[gasVer].aprfs[statDataCl[gasVer].PrimaryGas]  * np.asarray(statDataCl[gasVer].aprfs['AIRMASS']))

        if gasName[j].lower() == 'h2o': rPrfVMR_2  = np.asarray(statDataCl[gasVer].rprfs['HDO']) * sclfct

        #----------------------------------------
        # This is the mixing ratio for DRY AIR!!!
        #----------------------------------------
        rPrfDry[gasVer] = np.asarray(statDataCl[gasVer].rprfs[statDataCl[gasVer].PrimaryGas]) / (1.0 - waterVMR[gasVer]) * sclfct       
        
        #----------------------------------
        # Read Summary data (For filtering)
        #----------------------------------
        statDataCl[gasVer].readsummary()
        rms[gasVer]     = np.asarray(statDataCl[gasVer].summary[statDataCl[gasVer].PrimaryGas+'_FITRMS'])       
    
        #--------------------
        # Call to filter data
        #--------------------
        if fltrFlg: statDataCl[gasVer].fltrData(statDataCl[gasVer].PrimaryGas,mxrms=maxRMS[j],rmsFlg=rmsFlg, minDOF=minDOF[j],  dofFlg=dofFlg, 
        tcFlg=tcNegFlg, pcFlg=pcNegFlg , cnvrgFlg=True)
        else:       statDataCl[gasVer].inds = np.array([]) 

        #--------------------------------------------
        # Read Error data to get AVK and profile DOFs
        #-------------------------------------------------
        # Determine if AVK has been created via sfit4 core
        # code or via python error analysis
        #-------------------------------------------------     
        if errorFlg:   # Read AVK from error output
            statDataCl[gasVer].readError(totFlg=True,sysFlg=True,randFlg=True,vmrFlg=True,avkFlg=True,KbFlg=False)
            tot_rnd[gasVer] = np.array(statDataCl[gasVer].error['Total random uncertainty'])
            tot_sys[gasVer] = np.array(statDataCl[gasVer].error['Total systematic uncertainty'])
            totClmn_e[gasVer] = np.sqrt(tot_rnd[gasVer]**2 + tot_sys[gasVer]**2)

            npnts    = np.shape(statDataCl[gasVer].error['Total_Random_Error'])[0]
            nlvls    = np.shape(alt[gasVer])[0]
            rand_err[gasVer] = np.zeros((npnts,nlvls))
            sys_err[gasVer]  = np.zeros((npnts,nlvls))

            
            for i in range(npnts):
                rand_err[gasVer][i,:] = np.diag(statDataCl[gasVer].error['Total_Random_Error_VMR'][i][:,:])
                sys_err[gasVer][i,:]  = np.diag(statDataCl[gasVer].error['Total_Systematic_Error_VMR'][i][:,:])

            
            tot_err[gasVer]  = np.sqrt(rand_err[gasVer] + sys_err[gasVer] ) * sclfct         
            rand_err[gasVer] = np.sqrt(rand_err[gasVer]) * sclfct
            sys_err[gasVer]  = np.sqrt(sys_err[gasVer])  * sclfct

            rand_cmpnts = statDataCl[gasVer].randErrDiag
            sys_cmpnts  = statDataCl[gasVer].sysErrDiag
                
            #---------------------
            # Get averaging kernel
            #---------------------   
            avkSCF[gasVer]  = np.delete(np.asarray(statDataCl[gasVer].error['AVK_scale_factor']),statDataCl[gasVer].inds,axis=0)
            avkVMR[gasVer]  = np.delete(np.asarray(statDataCl[gasVer].error['AVK_vmr']),statDataCl[gasVer].inds,axis=0)   
            dofs            = np.diagonal(avkSCF[gasVer],axis1=1,axis2=2)
            avkSCFav[gasVer]  = np.mean(avkSCF[gasVer],axis=0)    
            avkVMRav[gasVer]  = np.mean(avkVMR[gasVer],axis=0)                   
                
            dofsAvg[gasVer]    = np.diag(avkSCFav[gasVer])
            dofsAvg_cs[gasVer] = np.cumsum(np.diag(avkSCFav[gasVer])[::-1])[::-1] 


            
        else:        # Read AVK from sfit4 output (only contains scaled AVK)
            avkSCFi = []
            for d in statDataCl[gasVer].dirLst:
                lines  = dc.tryopen( d + statDataCl[gasVer].ctl['file.out.ak_matrix'][0])
                if not lines: continue
                avkSCFi.append(np.array( [ [ float(x) for x in line.strip().split() ] for line in lines[2:] ] ))
                
            if not statDataCl[gasVer].readPrfFlgApr[statDataCl[gasVer].PrimaryGas]: statDataCl[gasVer].readprfs([statDataCl[gasVer].PrimaryGas],retapFlg=0)   # Apriori Profiles
                      
            avkSCF[gasVer]  = np.asarray(avkSCFi)
            nobs            = np.shape(avkSCF[gasVer])[0]
            n_layer         = np.shape(avkSCF[gasVer])[1]
            avkVMR[gasVer]  = np.zeros((nobs,n_layer,n_layer))
    
            for obs in range(0,nobs):
                Iapriori        = np.zeros((n_layer,n_layer))
                IaprioriInv     = np.zeros((n_layer,n_layer))
                np.fill_diagonal(Iapriori,statDataCl[gasVer].aprfs[statDataCl[gasVer].PrimaryGas.upper()][obs])
                np.fill_diagonal(IaprioriInv, 1.0 / (statDataCl[gasVer].aprfs[statDataCl[gasVer].PrimaryGas.upper()][obs]))
                avkVMR[gasVer][obs,:,:] = np.dot(np.dot(Iapriori,np.squeeze(avkSCF[gasVer][obs,:,:])),IaprioriInv)       
                
            avkSCF[gasVer]     = np.delete(avkSCF[gasVer],statDataCl[gasVer].inds,axis=0)
            avkVMR[gasVer]     = np.delete(avkVMR[gasVer],statDataCl[gasVer].inds,axis=0)
            dofs               = np.diagonal(avkSCF[gasVer],axis1=1,axis2=2)
            avkSCFav[gasVer]     = np.mean(avkSCF[gasVer],axis=0)
            avkVMRav[gasVer]     = np.mean(avkVMR[gasVer],axis=0)
            
            dofsAvg[gasVer]    = np.diag(avkSCFav[gasVer])
            dofsAvg_cs[gasVer] = np.cumsum(np.diag(avkSCFav[gasVer])[::-1])[::-1]
            
        #--------------------------------------
        # Remove retrieval data based on filter
        #--------------------------------------
        nfltr             = len(statDataCl[gasVer].inds)
        rms[gasVer]       = np.delete(rms[gasVer],statDataCl[gasVer].inds)
        ntot              = len(rms[gasVer])
        dates[gasVer]     = np.delete(dates[gasVer],statDataCl[gasVer].inds)
        totClmn[gasVer]   = np.delete(totClmn[gasVer],statDataCl[gasVer].inds)

        rPrfVMR[gasVer]   = np.delete(rPrfVMR[gasVer],statDataCl[gasVer].inds,axis=0)   
        rPrfMol[gasVer]   = np.delete(rPrfMol[gasVer],statDataCl[gasVer].inds,axis=0)   
        rPrfDry[gasVer]   = np.delete(rPrfDry[gasVer],statDataCl[gasVer].inds,axis=0)   
        Airmass[gasVer]   = np.delete(Airmass[gasVer],statDataCl[gasVer].inds,axis=0)
        TCdry[gasVer]     = np.delete(TCdry[gasVer],statDataCl[gasVer].inds)

        aPrfVMR[gasVer]   = np.delete(aPrfVMR[gasVer],statDataCl[gasVer].inds,axis=0)   
        aPrfMol[gasVer]   = np.delete(aPrfMol[gasVer],statDataCl[gasVer].inds,axis=0)

        if gasName[j].lower() == 'h2o':  rPrfVMR_2 = np.delete(rPrfVMR_2,statDataCl[gasVer].inds,axis=0)

        if errorFlg:
            rand_err[gasVer] = np.delete(rand_err[gasVer],statDataCl[gasVer].inds,axis=0)
            sys_err[gasVer]  = np.delete(sys_err[gasVer],statDataCl[gasVer].inds,axis=0)  
            tot_err[gasVer]  = np.delete(tot_err[gasVer],statDataCl[gasVer].inds,axis=0)

            tot_rnd[gasVer] = np.delete(tot_rnd[gasVer],statDataCl[gasVer].inds)
            tot_sys[gasVer] = np.delete(tot_sys[gasVer],statDataCl[gasVer].inds)
            totClmn_e[gasVer] = np.delete(totClmn_e[gasVer],statDataCl[gasVer].inds)

            
            for k in sys_cmpnts:
                sys_cmpnts[k] = np.delete(sys_cmpnts[k],statDataCl[gasVer].inds,axis=0)
                
            for k in rand_cmpnts:
                rand_cmpnts[k] = np.delete(rand_cmpnts[k],statDataCl[gasVer].inds,axis=0)


        #-------------------------------------------------
        # Calculate partial columns and weighted VMR
        #-------------------------------------------------
        for p, pcol in enumerate(pCols):
            
            ind1             = mf.nearestind(pcol[0], alt[gasVer])
            ind2             = mf.nearestind(pcol[1], alt[gasVer])
            #vmrP[gasVer][p]  = np.average(rPrfVMR[gasVer][:,ind2:ind1],axis=1,weights=Airmass[gasVer][:,ind2:ind1])           
            #TCp[gasVer][p]   = np.sum(rPrfMol[gasVer][:,ind2:ind1],axis=1)
            vmrP.setdefault(gasVer+str(p) ,[]).append(np.average(rPrfVMR[gasVer][:,ind2:ind1],axis=1,weights=Airmass[gasVer][:,ind2:ind1]))         
            TCp.setdefault(gasVer+str(p), []).append(np.sum(rPrfMol[gasVer][:,ind2:ind1],axis=1))

        for p, pcol in enumerate(pCols):
            vmrP[gasVer+str(p)] = np.asarray(vmrP[gasVer+str(p)])
            TCp[gasVer+str(p)] = np.asarray(TCp[gasVer+str(p)])
    
    #-------------------------------------------------
   
    #---------------------------------
    # Defining variables (FTIR) - H2O
    #---------------------------------
    maxalt        = 16.0
    h2over        = gasName[0]+'_'+ver[0]
    
    indsalt       = np.where(alt[h2over] <= maxalt)[0]
    hours         = [d.hour for d in dates[h2over]]

    alt_h2o         = alt[h2over][indsalt]
    TC_h2o          = totClmn[h2over]
    Dates_h2o       = dates[h2over]
    Prf_h2o         = rPrfVMR[h2over][:, indsalt]
    avkSCF_h2o      = avkSCF[h2over][:, indsalt[0]:, indsalt[0]:]
    avkSCav_h2o     = np.mean(avkSCF_h2o, axis=0)
    avkVMR_h2o      = avkVMR[h2over][:, indsalt[0]:, indsalt[0]:]
    avkVMRav_h2o     = np.mean(avkVMR_h2o, axis=0)
    aPrf_h2o        = aPrfVMR[h2over][:, indsalt]
    Airmass_h2o     = Airmass[h2over][:, indsalt]
    rPrfMol_h2o     = rPrfMol[h2over][:,indsalt]
    TCpc_h2o        = np.sum(rPrfMol_h2o,axis=1)
    VMR_ns_h2o      = rPrfVMR[h2over][:, -1]
    VMRa_ns_h2o     = aPrfVMR[h2over][:, -1]

    Prf_hdoW        = rPrfVMR_2[:, indsalt]*3.107e-4
    VMR_ns_hdoW     = rPrfVMR_2[:, -1]
    Prf_dDW         = (np.divide(Prf_hdoW, Prf_h2o)/3.1152e-4 - 1) *1000.
    dDftsW          = Prf_dDW[:, -1]

    vmrP_h2o        = []

    for p, pcol in enumerate(pCols):
        vmrP_h2o.append(vmrP[h2over+str(p)])

    if errorFlg: 
        TCe_h2o         = totClmn_e[h2over]
        tot_err_h2o     = tot_err[h2over][:, indsalt]
        VMR_e_ns_h2o    = tot_err[h2over][:, -1]

    doy_h2o  = mf.toYearFraction(Dates_h2o)

    #---------------------------------
    # Defining variables (FTIR) - HDO
    #---------------------------------
    hdover        = gasName[1]+'_'+ver[1]
    indsalt       = np.where(alt[hdover] <= maxalt)[0]
    hours         = [d.hour for d in dates[hdover]]
    
    TC_hdo          = totClmn[hdover]
    Dates_hdo       = dates[hdover]
    Prf_hdo         = rPrfVMR[hdover][:, indsalt]*3.107e-4
    avkSCF_hdo      = avkSCF[hdover][:, indsalt[0]:, indsalt[0]:]
    avkSCav_hdo     = np.mean(avkSCF_hdo, axis=0)
    avkVMR_hdo      = avkVMR[hdover][:, indsalt[0]:, indsalt[0]:]
    avkVMRav_hdo    = np.mean(avkVMR_hdo, axis=0)
    aPrf_hdo        = aPrfVMR[hdover][:, indsalt]*3.107e-4
    Airmass_hdo     = Airmass[hdover][:, indsalt]
    rPrfMol_hdo     = rPrfMol[hdover][:,indsalt]
    TCpc_hdo        = np.sum(rPrfMol_hdo, axis=1)
    VMR_ns_hdo      = Prf_hdo[:, -1]
    VMRa_ns_hdo     = aPrfVMR[hdover][:, -1]*3.107e-4

    vmrP_hdo        = []

    for p, pcol in enumerate(pCols):
        vmrP_hdo.append(vmrP[hdover+str(p)]*3.107e-4)

    
    if errorFlg: 
        TCe_hdo         = totClmn_e[hdover]
        tot_err_hdo     = tot_err[hdover][:, indsalt]
        VMR_e_ns_hdo    = tot_err[hdover][:, -1]

    doy_hhdo  = mf.toYearFraction(Dates_hdo)

    #---------------------------------
    # Defining dD (FTIR)
    #---------------------------------
    Prf_h2o_i      = interpolate.interp1d(doy_h2o, Prf_h2o, axis=0, bounds_error=False)(doy_hhdo)
    aPrf_h2o_i     = interpolate.interp1d(doy_h2o, aPrf_h2o, axis=0, bounds_error=False)(doy_hhdo)
    VMR_ns_h2o_i   = interpolate.interp1d(doy_h2o, VMR_ns_h2o, bounds_error=False)(doy_hhdo)
    VMRa_ns_h2o_i  = interpolate.interp1d(doy_h2o, VMRa_ns_h2o, bounds_error=False)(doy_hhdo)
    avkVMR_h2o_i   = interpolate.interp1d(doy_h2o, avkVMR_h2o, axis=0, bounds_error=False)(doy_hhdo)


    Prf_dD      = (np.divide(Prf_hdo, Prf_h2o_i)/3.1152e-4 - 1.0) *1000.
    TC_dD       = np.sum(Prf_dD, axis=1)
    TC_dD2      = (np.divide(TC_hdo, TC_h2o)/3.1152e-4 - 1.0) *1000.
    dDfts       = (np.divide(VMR_ns_hdo, VMR_ns_h2o_i)/3.1152e-4 - 1.0) *1000.

    vmrP_dD     = []

    for p, pcol in enumerate(pCols):
        vmrP_h2o_i  = interpolate.interp1d(doy_h2o, vmrP_h2o[p], bounds_error=False)(doy_hhdo)

        vmrP_dD.append((np.divide(vmrP_hdo[p], vmrP_h2o_i)/3.1152e-4 - 1.0) *1000.)


    #---------------------------------
    # Interpolation of in-situ to FTS
    #---------------------------------
    doyin_i            = mf.toYearFraction(Date_in)

    hours              =  [float(d.hour) for d in Dates_h2o]
    hours              = np.asarray(hours)
  
    insitu_interp      = interpolate.interp1d(doyin_i, insitu, bounds_error=False)(doy_hhdo)
    HDOinsitu_interp   = interpolate.interp1d(doyin_i, HDOinsitu, bounds_error=False)(doy_hhdo)
    dDinsitu_r_interp  = interpolate.interp1d(doyin_i, dDinsitu_r, bounds_error=False)(doy_hhdo)
    dDinsitu_interp    = interpolate.interp1d(doyin_i, dDinsitu, bounds_error=False)(doy_hhdo)

    insitu_interp_s    = []
    HDOinsitu_interp_s = []

    ##Smoothing using FTS AK
    for itime in range(len(insitu_interp)):
        insitu_interp_s.append(aPrf_h2o_i[itime, -1] + np.dot(avkVMR_h2o_i[itime, -1, -1], (insitu_interp[itime] -  aPrf_h2o_i[itime, -1])) )
        HDOinsitu_interp_s.append(aPrf_hdo[itime, -1] + np.dot(avkVMR_hdo[itime, -1, -1], (HDOinsitu_interp[itime] -  aPrf_hdo[itime, -1])) )

    insitu_interp_s    = np.asarray(insitu_interp_s)
    HDOinsitu_interp_s = np.asarray(HDOinsitu_interp_s)

    dDinsitu_interp_s  = (np.divide(HDOinsitu_interp_s, insitu_interp_s)/3.1152e-4 - 1.0) *1000.

    #---------------------------------
    # Defining Good Data
    #---------------------------------
    inds = np.where( (insitu_interp >= 0.) & (insitu_interp <= 10000.) )[0]
    #inds = np.where( (insitu_interp >= 0.) )[0]

    if saveFlg: pdfsav = PdfPages(pltFile)
    clmap        = 'jet'
    cm           = plt.get_cmap(clmap)
    dayLc        = DayLocator()
    yearsLc      = YearLocator()
    monthLc      = MonthLocator()
    mondays      = WeekdayLocator(MONDAY)
    DateFmt      = DateFormatter('%b %d')

    #--------------------------------
    # FIGURE: H2O Profiles as a function of Month
    #--------------------------------
    month = np.array([d.month for d in Dates_h2o])
    fig,(ax1,ax2)  = plt.subplots(1,2,sharey=True)
    cm             = plt.get_cmap(clmap)
    cNorm          = colors.Normalize( vmin=np.min(month), vmax=np.max(month) )
    scalarMap      = mplcm.ScalarMappable( norm=cNorm, cmap=cm )
    
    scalarMap.set_array(month)
    
    ax1.set_color_cycle( [scalarMap.to_rgba(x) for x in month] )
    ax2.set_color_cycle( [scalarMap.to_rgba(x) for x in month] )
    
    for i in range(len(month)):
        ax1.plot(Prf_h2o[i,:]/1e3,alt_h2o,linewidth=0.75)
        ax2.plot(Prf_h2o[i,:]/1e3,alt_h2o,linewidth=0.75)
    
    #ax1.plot(aprPrf[gas],alt_h2o,'k--',linewidth=4,label='A priori')
    #ax2.plot(aprPrf[gas],alt,'k--',linewidth=4,label='A priori')
    
    ax1.set_ylabel('Altitude [km]', fontsize=14)
    ax1.set_xlabel('VMR [ppm, x10$^3$]', fontsize=14)
    ax2.set_xlabel('Log VMR [ppm, x10$^3$]', fontsize=14)
    
    ax1.grid(True,which='both')
    ax2.grid(True,which='both')
    ax2.set_xscale('log')
    
    cbar = fig.colorbar(scalarMap,orientation='vertical')
    cbar.set_label('Month', fontsize=14)
    
    ax1.legend(prop={'size':12})
    ax2.legend(prop={'size':12})
    
    ax1.tick_params(axis='x',which='both',labelsize=14)
    ax2.tick_params(axis='x',which='both',labelsize=14)  
    plt.suptitle('FTS Profiles: {}'.format('H$_2$O'), fontsize=16)
    
    if saveFlg:     
        pdfsav.savefig(fig,dpi=200)
    else:           
        plt.show(block=False)

    #---------------------------------
    # Plot : Averaging Kernel: H2O
    #---------------------------------
    clmap        = 'jet'
    cm           = plt.get_cmap(clmap)
    fig       = plt.figure(figsize=(7,7))
    gs        = gridspec.GridSpec(1,2,width_ratios=[3,1])
    ax        = plt.subplot(gs[0])
    axb       = plt.subplot(gs[1])
    cm        = plt.get_cmap(clmap)
    cNorm     = colors.Normalize(vmin=np.min(alt_h2o), vmax=np.max(alt_h2o))
    scalarMap = mplcm.ScalarMappable(norm=cNorm,cmap=clmap)
    scalarMap.set_array(alt_h2o)
    ax.set_color_cycle([scalarMap.to_rgba(x) for x in alt_h2o])
    
    for i in range(len(alt_h2o)):
        ax.plot(avkSCav_h2o[i,:], alt_h2o)
        
    ax.set_ylabel('Altitude [km]', fontsize=14)
    ax.set_xlabel('Averaging Kernels', fontsize=14)
    ax.grid(True)

    cbaxes = fig.add_axes([0.52, 0.3, 0.03, 0.5]) 
    cbar = fig.colorbar(scalarMap, orientation='vertical', cax = cbaxes)
    cbar.set_label('Altitude [km]', fontsize=14)
    ax.set_title('H2O Averaging Kernels Scale Factor', fontsize=14)
    ax.tick_params(labelsize=14)
    #ax.set_ylim(1, 15)
    if loc.lower() == 'fl0': ax.set_ylim(1, 15)
    if loc.lower() == 'mlo': ax.set_ylim(3, 15)
    
    axb.plot(np.sum(avkSCav_h2o,axis=0), alt_h2o,color='k')
    axb.grid(True)
    axb.set_xlabel('Averaging Kernel Area', fontsize=14)
    axb.tick_params(axis='x',which='both',labelsize=14)
    major_ticks = np.arange(0, 3, 1)
    axb.set_xticks(major_ticks) 
    #axb.set_ylim(1, 15)
    if loc.lower() == 'fl0': axb.set_ylim(1, 15)
    if loc.lower() == 'mlo': axb.set_ylim(3, 15)
    #axb.tick_params(labelsize=14) 

    if saveFlg: pdfsav.savefig(fig,dpi=200)
    #else:       plt.show(block=False)

    #--------------------------------
    # FIGURE: HDO Profiles as a function of Month
    #--------------------------------
    month = np.array([d.month for d in Dates_hdo])
    fig,(ax1,ax2)  = plt.subplots(1,2,sharey=True)
    cm             = plt.get_cmap(clmap)
    cNorm          = colors.Normalize( vmin=np.min(month), vmax=np.max(month) )
    scalarMap      = mplcm.ScalarMappable( norm=cNorm, cmap=cm )
    
    scalarMap.set_array(month)
    
    ax1.set_color_cycle( [scalarMap.to_rgba(x) for x in month] )
    ax2.set_color_cycle( [scalarMap.to_rgba(x) for x in month] )
    
    for i in range(len(month)):
        ax1.plot(Prf_hdo[i,:],alt_h2o,linewidth=0.75)
        ax2.plot(Prf_hdo[i,:],alt_h2o,linewidth=0.75)
    
    #ax1.plot(aprPrf[gas],alt_h2o,'k--',linewidth=4,label='A priori')
    #ax2.plot(aprPrf[gas],alt,'k--',linewidth=4,label='A priori')
    
    ax1.set_ylabel('Altitude [km]', fontsize=14)
    ax1.set_xlabel('VMR [ppm]', fontsize=14)
    ax2.set_xlabel('Log VMR [ppm]', fontsize=14)
    
    ax1.grid(True,which='both')
    ax2.grid(True,which='both')
    ax2.set_xscale('log')
    
    cbar = fig.colorbar(scalarMap,orientation='vertical')
    cbar.set_label('Month', fontsize=14)
    
    ax1.legend(prop={'size':12})
    ax2.legend(prop={'size':12})
    
    ax1.tick_params(axis='x',which='both',labelsize=14)
    ax2.tick_params(axis='x',which='both',labelsize=14)  
    plt.suptitle('FTS Profiles: {}'.format('HDO'), fontsize=16)
    
    if saveFlg:     
        pdfsav.savefig(fig,dpi=200)
    else:           
        plt.show(block=False)

    #---------------------------------
    # Plot : Averaging Kernel: HDO
    #---------------------------------
    clmap        = 'jet'
    cm           = plt.get_cmap(clmap)
    fig       = plt.figure(figsize=(7,7))
    gs        = gridspec.GridSpec(1,2,width_ratios=[3,1])
    ax        = plt.subplot(gs[0])
    axb       = plt.subplot(gs[1])
    cm        = plt.get_cmap(clmap)
    cNorm     = colors.Normalize(vmin=np.min(alt_h2o), vmax=np.max(alt_h2o))
    scalarMap = mplcm.ScalarMappable(norm=cNorm,cmap=clmap)
    scalarMap.set_array(alt_h2o)
    ax.set_color_cycle([scalarMap.to_rgba(x) for x in alt_h2o])
    
    for i in range(len(alt_h2o)):
        ax.plot(avkSCav_hdo[i,:], alt_h2o)
        
    ax.set_ylabel('Altitude [km]', fontsize=14)
    ax.set_xlabel('Averaging Kernels', fontsize=14)
    ax.grid(True)

    cbaxes = fig.add_axes([0.52, 0.3, 0.03, 0.5]) 
    cbar = fig.colorbar(scalarMap, orientation='vertical', cax = cbaxes)
    cbar.set_label('Altitude [km]', fontsize=14)
    ax.set_title('HDO Averaging Kernels Scale Factor', fontsize=14)
    ax.tick_params(labelsize=14)
    #ax.set_ylim(1, 15)
    if loc.lower() == 'fl0': ax.set_ylim(1, 15)
    if loc.lower() == 'mlo': ax.set_ylim(3, 15)
    
    axb.plot(np.sum(avkSCav_hdo,axis=0), alt_h2o,color='k')
    axb.grid(True)
    axb.set_xlabel('Averaging Kernel Area', fontsize=14)
    axb.tick_params(axis='x',which='both',labelsize=14)
    major_ticks = np.arange(0, 3, 1)
    axb.set_xticks(major_ticks) 
    #axb.set_ylim(1, 15)
    if loc.lower() == 'fl0': axb.set_ylim(1, 15)
    if loc.lower() == 'mlo': axb.set_ylim(3, 15)
    #axb.tick_params(labelsize=14) 

    if saveFlg: pdfsav.savefig(fig,dpi=200)
    #else:       plt.show(block=False)


    #--------------------------------
    # FIGURE: dD Profiles as a function of Month
    #--------------------------------
    month = np.array([d.month for d in Dates_hdo])
    fig,(ax1)  = plt.subplots(1, sharey=True)
    cm             = plt.get_cmap(clmap)
    cNorm          = colors.Normalize( vmin=np.min(month), vmax=np.max(month) )
    scalarMap      = mplcm.ScalarMappable( norm=cNorm, cmap=cm )
    
    scalarMap.set_array(month)
    
    ax1.set_color_cycle( [scalarMap.to_rgba(x) for x in month] )
   
    for i in range(len(month)):
        ax1.plot(Prf_dD[i,:],alt_h2o,linewidth=0.75)
    
    #ax1.plot(aprPrf[gas],alt_h2o,'k--',linewidth=4,label='A priori')
    #ax2.plot(aprPrf[gas],alt,'k--',linewidth=4,label='A priori')
    
    ax1.set_ylabel('Altitude [km]', fontsize=14)
    ax1.set_xlabel('dD [per mil]', fontsize=14)
    ax1.grid(True,which='both')
   
    cbar = fig.colorbar(scalarMap,orientation='vertical')
    cbar.set_label('Month', fontsize=14)
    
    ax1.legend(prop={'size':12})
   
    ax1.tick_params(axis='x',which='both',labelsize=14)
   
    plt.suptitle('FTS Profiles: {}'.format('dD'), fontsize=16)
    
    if saveFlg:     
        pdfsav.savefig(fig,dpi=200)
    else:           
        plt.show(block=False)


    #---------------------------------
    # FIGURE - TIME SERIES OF COLUMNS
    #---------------------------------
    fig, ax = plt.subplots(figsize=(12,6))
    
    ax.scatter(Dates_hdo[inds], TC_hdo[inds]/1e3, facecolors='white', s=40, color='r', label='Retrieved FTIR - near surface')
   
    ax.grid(True)
    ax.set_ylabel('VMR [ppm, x10$^3$]',fontsize=14)
    ax.tick_params(labelsize=14)
    ax.xaxis.set_minor_locator(dayLc)
    ax.xaxis.set_major_formatter(DateFmt)
    ax.legend(prop={'size':14})

    plt.suptitle('Time series of H2O', fontsize=16)

    fig.subplots_adjust(left=0.1, right=0.95)

    if saveFlg:     
        pdfsav.savefig(fig,dpi=200)
    else:           
        plt.show(block=False)
   
    #---------------------------------
    # FIGURE - TIME SERIES OF H2O
    #---------------------------------
    fig, ax = plt.subplots(figsize=(12,6))
    
    ax.scatter(Dates_hdo[inds], insitu_interp[inds]/1e3, facecolors='white', s=40, color='b', label='in-situ (measured)')
    ax.scatter(Dates_hdo[inds], insitu_interp_s[inds]/1e3, facecolors='white', s=40, color='green', label='in-situ (smooth by FTIR AK)')
    
    if errorFlg: ax.errorbar(Dates_hdo[inds],VMR_ns_h2o[inds]/1e3, yerr=VMR_e_ns_h2o[inds]/1e3, fmt='o', markersize=0, color='r', ecolor='r')
    ax.scatter(Dates_hdo[inds], VMR_ns_h2o_i[inds]/1e3, facecolors='white', s=40, color='r', label='Retrieved FTIR - near surface')
    ax.scatter(Dates_hdo[inds], VMRa_ns_h2o_i[inds]/1e3, facecolors='white', s=40, color='gray', label='a priori (ERA) - near surface')
    
    ax.grid(True)
    ax.set_ylabel('VMR [ppm, x10$^3$]',fontsize=14)
    ax.tick_params(labelsize=14)
    ax.xaxis.set_minor_locator(dayLc)
    ax.xaxis.set_major_formatter(DateFmt)
    ax.legend(prop={'size':14})

    plt.suptitle('Time series of H2O', fontsize=16)

    fig.subplots_adjust(left=0.1, right=0.95)

    if saveFlg:     
        pdfsav.savefig(fig,dpi=200)
    else:           
        plt.show(block=False)

    #---------------------------------
    # FIGURE - TIME SERIES OF HDO
    #---------------------------------
    fig, ax = plt.subplots(figsize=(12,6))

    ax.scatter(Dates_hdo[inds], HDOinsitu_interp[inds],   facecolors='white', s=40, color='b', label='in-situ (Calculated)')
    ax.scatter(Dates_hdo[inds], HDOinsitu_interp_s[inds], facecolors='white', s=40, color='green', label='in-situ (smooth by FTIR AK)')
    
    ax.scatter(Dates_hdo[inds], VMR_ns_hdo[inds], facecolors='white', s=40, color='r', label='Retrieved FTIR')
    ax.scatter(Dates_hdo[inds], VMRa_ns_hdo[inds], facecolors='white', s=40, color='gray', label='a priori (ERA)')

    
    ax.grid(True)
    ax.set_ylabel('VMR [ppm]',fontsize=14)
    ax.tick_params(labelsize=14)
    ax.xaxis.set_minor_locator(dayLc)
    ax.xaxis.set_major_formatter(DateFmt)
    ax.legend(prop={'size':14})

    plt.suptitle('Time series of HDO', fontsize=16)

    fig.subplots_adjust(left=0.1, right=0.95)

    if saveFlg:     
        pdfsav.savefig(fig,dpi=200)
    else:           
        plt.show(block=False)

    #---------------------------------
    # FIGURE - TIME SERIES OF dD
    #---------------------------------

    fig, ax = plt.subplots(figsize=(12,6))
   
    ax.scatter(Dates_hdo[inds], dDinsitu_interp[inds], facecolors='white', s=40, color='b', label='in-situ')
    ax.scatter(Dates_hdo[inds], dDinsitu_r_interp[inds], facecolors='white', s=40, color='yellow', label='in-situ (Calculated)')
    ax.scatter(Dates_hdo[inds], dDinsitu_interp_s[inds], facecolors='white', s=40, color='green', label='in-situ (smooth by FTIR AK)')
    ax.scatter(Dates_hdo[inds], dDfts[inds], facecolors='white', s=40, color='r', label='Retrieved FTIR')
   
    ax.grid(True)
    ax.set_ylabel('dD [per mil]',fontsize=14)
    ax.tick_params(labelsize=14)
    ax.set_ylim(-550, 400)
    ax.xaxis.set_minor_locator(dayLc)
    ax.xaxis.set_major_formatter(DateFmt)
    ax.legend(prop={'size':14})

    plt.suptitle('Time series of dD', fontsize=16)

    fig.subplots_adjust(left=0.1, right=0.95)

    if saveFlg:     
        pdfsav.savefig(fig,dpi=200)
    else:           
        plt.show(block=False)

    #---------------------------------
    # FIGURE - CORRELATION (H2O)
    #---------------------------------
    fig, ax = plt.subplots(figsize=(8,6))

    ax.scatter(insitu_interp[inds]/1e3, VMR_ns_h2o_i[inds]/1e3, facecolors='white', s=40, color='b', label='FTIR Vs in-situ')
    ax.scatter(insitu_interp_s[inds]/1e3, VMR_ns_h2o_i[inds]/1e3, facecolors='white', s=40, color='green', label='FTIR Vs in-situ (smooth by FTIR AK)')
    ax.scatter(insitu_interp_s[inds]/1e3, VMRa_ns_h2o_i[inds]/1e3, facecolors='white', s=40, color='gray', label='A priori (ERA) Vs in-situ (smooth by FTIR AK)')
          
    ax.set_ylabel('VMR [ppm, x10$^3$]',fontsize=14)
    ax.set_xlabel('VMR - in-situ [ppm, x10$^3$]',fontsize=14)
    ax.tick_params(labelsize=14)
    ax.grid(True)
    ax.set_xlim(0, 16)
    ax.set_ylim(0, 16)
    ax.legend(prop={'size':12})

    plt.suptitle('Correlation of H2O\n(near-surface FTS Vs in-situ)', fontsize=16)

    slope, intercept, r_value, p_value, std_err = stats.linregress(insitu_interp[inds], VMR_ns_h2o_i[inds])
   
    ax.text(0.025,0.92,"Slope: {0:.2f}".format(slope),transform=ax.transAxes,  fontsize=14, color='blue')
    ax.text(0.025,0.85,"Intercept: {:.2f}".format(intercept/1e3),transform=ax.transAxes,  fontsize=14, color='blue')
    ax.text(0.025,0.78,"r-value: {0:.2f}".format(r_value),transform=ax.transAxes,  fontsize=14, color='blue')

    slope, intercept, r_value, p_value, std_err = stats.linregress(insitu_interp_s[inds], VMR_ns_h2o_i[inds])

    ax.text(0.025,0.7,"Slope: {0:.2f}".format(slope),transform=ax.transAxes,  fontsize=14, color='green')
    ax.text(0.025,0.63,"Intercept: {:.2f}".format(intercept/1e3),transform=ax.transAxes,  fontsize=14, color='green')
    ax.text(0.025,0.56,"r-value: {0:.2f}".format(r_value),transform=ax.transAxes,  fontsize=14, color='green')
    
    bias = mf.bias(insitu_interp_s[inds]/1e3, VMR_ns_h2o_i[inds]/1e3)
    rmse = mf.rmse(insitu_interp_s[inds]/1e3, VMR_ns_h2o_i[inds]/1e3)
 
    print 'bias (H2O) = {}'.format(bias) 
    print 'rmse (H2O) = {}'.format(rmse) 

    slope, intercept, r_value, p_value, std_err = stats.linregress(insitu_interp_s[inds], VMRa_ns_h2o_i[inds])

    ax.text(0.75,0.7,"Slope: {0:.2f}".format(slope),transform=ax.transAxes,  fontsize=14, color='gray')
    ax.text(0.75,0.63,"Intercept: {:.2f}".format(intercept/1e3),transform=ax.transAxes,  fontsize=14, color='gray')
    ax.text(0.75,0.56,"r-value: {0:.2f}".format(r_value),transform=ax.transAxes,  fontsize=14, color='gray')

    if saveFlg:     
        pdfsav.savefig(fig,dpi=200)
    else:           
        plt.show(block=False)

     #---------------------------------
    # FIGURE - CORRELATION (HDO)
    #---------------------------------
    fig, ax = plt.subplots(figsize=(8,6))

    ax.scatter(HDOinsitu_interp[inds], VMR_ns_hdo[inds], facecolors='white', s=40, color='b', label='FTIR Vs in-situ')
    ax.scatter(HDOinsitu_interp_s[inds], VMR_ns_hdo[inds], facecolors='white', s=40, color='green', label='FTIR Vs in-situ (smooth by FTIR AK)')
          
    ax.set_ylabel('VMR [ppm]',fontsize=14)
    ax.set_xlabel('VMR - in-situ [ppm]',fontsize=14)
    ax.tick_params(labelsize=14)
    ax.grid(True)
    #ax.set_xlim(0, 16)
    #ax.set_ylim(0, 16)
    ax.legend(prop={'size':12})

    plt.suptitle('Correlation of HDO\n(near-surface FTS Vs in-situ)', fontsize=16)

    slope, intercept, r_value, p_value, std_err = stats.linregress(HDOinsitu_interp[inds], VMR_ns_hdo[inds])
   
    ax.text(0.025,0.92,"Slope: {0:.2f}".format(slope),transform=ax.transAxes,  fontsize=14, color='blue')
    ax.text(0.025,0.85,"Intercept: {:.3f}".format(intercept),transform=ax.transAxes,  fontsize=14, color='blue')
    ax.text(0.025,0.78,"r-value: {0:.2f}".format(r_value),transform=ax.transAxes,  fontsize=14, color='blue')

    slope, intercept, r_value, p_value, std_err = stats.linregress(HDOinsitu_interp_s[inds], VMR_ns_hdo[inds])

    ax.text(0.025,0.7,"Slope: {0:.2f}".format(slope),transform=ax.transAxes,  fontsize=14, color='green')
    ax.text(0.025,0.63,"Intercept: {:.3f}".format(intercept),transform=ax.transAxes,  fontsize=14, color='green')
    ax.text(0.025,0.56,"r-value: {0:.2f}".format(r_value),transform=ax.transAxes,  fontsize=14, color='green')
    
    bias = mf.bias(HDOinsitu_interp_s[inds], VMR_ns_hdo[inds])
    rmse = mf.rmse(HDOinsitu_interp_s[inds], VMR_ns_hdo[inds])
 
    print 'bias (HDO) = {}'.format(bias) 
    print 'rmse (HDO) = {}'.format(rmse) 

   
    if saveFlg:     
        pdfsav.savefig(fig,dpi=200)
    else:           
        plt.show(block=False)
    
    
    #---------------------------------
    # FIGURE - CORRELATION (dD)
    #---------------------------------
    fig, ax = plt.subplots(figsize=(8,6))

    ax.scatter(dDinsitu_r_interp[inds], dDfts[inds], facecolors='white', s=40, color='b', label='in-situ (Calculated)')
    ax.scatter(dDinsitu_interp[inds],   dDfts[inds], facecolors='white', s=40, color='r', label='in-situ (Measured)')
    ax.scatter(dDinsitu_interp_s[inds], dDfts[inds], facecolors='white', s=40, color='green',label='in-situ smooth (Calculated)')
          
    ax.set_ylabel('dD - FTIR [per mil]',fontsize=14)
    ax.set_xlabel('dD - in-situ [per mil]',fontsize=14)
    ax.tick_params(labelsize=14)
    ax.grid(True)
    ax.set_xlim(-550, 400)
    ax.set_ylim(-550, 400)
    ax.legend(prop={'size':14})

    plt.suptitle('Correlation of dD\n(near-surface FTS Vs in-situ)', fontsize=16)

    slope, intercept, r_value, p_value, std_err = stats.linregress(dDinsitu_interp[inds], dDfts[inds])
    bias = mf.bias(dDinsitu_interp[inds], dDfts[inds])
    rmse = mf.rmse(dDinsitu_interp[inds], dDfts[inds])
   
    ax.text(0.025,0.92,"Slope: {0:.2f}".format(slope),transform=ax.transAxes,  fontsize=14, color='red')
    ax.text(0.025,0.85,"Intercept: {:.2f}".format(intercept/1e3),transform=ax.transAxes,  fontsize=14, color='red')
    ax.text(0.025,0.78,"r-value: {0:.2f}".format(r_value),transform=ax.transAxes,  fontsize=14, color='red')

    ax.text(0.75,0.92,"Bias: {0:.2f}".format(bias),transform=ax.transAxes,  fontsize=14, color='red')
    ax.text(0.75,0.85,"RMSE: {0:.2f}".format(rmse),transform=ax.transAxes,  fontsize=14, color='red')

    print 'bias (dD)= {}'.format(bias) 
    print 'rmse (dD)= {}'.format(rmse)
    weights      = np.ones_like(dDinsitu_interp[inds])

    res    = mf.fit_driftfourier(dDinsitu_interp[inds], dDfts[inds], weights, 2)
    f_drift, f_fourier, f_driftfourier = res[3:6]

    print "Fitted slope: {0:.3E} ({1:.3f}%)".format(res[1],res[1]/np.mean(dDinsitu_interp[inds])*100.0)
    print "Fitted intercept at xmin: {:.3E}".format(res[0])
    print "STD of residuals: {0:.3E} ({1:.3f}%)".format(res[6],res[6]/np.mean(dDinsitu_interp[inds])*100.0)

    slope, intercept, r_value, p_value, std_err = stats.linregress(dDinsitu_interp_s[inds], dDfts[inds])

    ax.text(0.025,0.7,"Slope: {0:.2f}".format(slope),transform=ax.transAxes,  fontsize=14, color='green')
    ax.text(0.025,0.63,"Intercept: {:.3f}".format(intercept),transform=ax.transAxes,  fontsize=14, color='green')
    ax.text(0.025,0.56,"r-value: {0:.2f}".format(r_value),transform=ax.transAxes,  fontsize=14, color='green')


    if saveFlg:     
        pdfsav.savefig(fig,dpi=200)
        pdfsav.close()
    else:           
        plt.show(block=False)
        user_input = raw_input('Press any key to exit >>> ')
        sys.exit()    









if __name__ == "__main__":
    main()



