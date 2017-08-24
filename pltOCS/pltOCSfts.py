#! /usr/bin/python2.7
#----------------------------------------------------------------------------------------
# Name:
#         pltOCSfts.py
#
# Purpose:
#         Plot OCS (OCS Study in Thule)
#
# Notes:
#   
#
# Version History:
#       Created, May, 2016  Ivan Ortega (iortega@ucar.edu)

#----------------------------------------------------------------------------------------

                                    #-------------------------#
                                    # Import Standard modules #
                                    #-------------------------#
from scipy.io import netcdf
import os
import datetime as dt
import numpy as np
import numpy.ma as ma
import sys
import glob

from scipy import interpolate
from scipy.integrate import simps

import matplotlib.dates as md
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
from itertools import izip
from numpy import *
import myfunctions as mf
import dataOutClass as dc
from collections                     import OrderedDict
import PltClass as mp
import sondeclass as rs

                                    #-------------------------#
                                    # Define helper functions #
                                    #-------------------------#

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
                                    #----------------------------#
                                    #                            #
                                    #        --- Main---         #
                                    #                            #
                                    #----------------------------#

def main():
    #-----------------------------------------------------------------------------------------
    #                             Initializations for FTIR
    #-----------------------------------------------------------------------------------------

    loc        = 'tab'

    gasName    = ['ocs']             
    ver        = ['Current_v8']        
    ctlF       = ['sfit4_v8.ctl'] 

    saveFlg    = True 
    pltFile    = '/data1/projects/ocstab/fig/'+loc.lower()+'_Results_ocs.pdf'

    #------
    # Flags
    #------
    errorFlg   = True                  # Flag to process error data
    fltrFlg    = True                   # Flag to filter the data
    byYrFlg    = False                  # Flag to create plots for each individual year in date range
    szaFlg     = True                   # Flag to filter based on min and max SZA
    dofFlg     = True                   # Flag to filter based on min DOFs
    pcNegFlg   = True                   # Flag to filter profiles with negative partial columns
    tcNegFlg   = True                  # Flagsag to filter profiles with negative total columns
    tcMMFlg    = False                  # Flag to filter based on min and max total column amount
    cnvrgFlg   = True                   # Flag to filter profiles that did not converge
    rmsFlg     = True                   # Flag to filter based on max RMS
    chiFlg     = False                   # Flag to filter based on max CHI_2_Y
    mnthFlg    = False                  # Flag to filter based on 

    mnths      = [6,7,8]                # Months to filter on (these are the months to include data)
    maxRMS     = [2.0]                         # Max Fit RMS to filter data. Data is filtered according to <= maxrms
    minDOF     = [0.0]                   # Min DOFs for filtering
    minSZA     = 0.0                   # Min SZA for filtering
    maxSZA     = 90.0                   # Max SZA for filtering
    maxCHI     = 2.0                    # Max CHI_y_2 value
    maxTC      = 5.0E24                # Max Total column amount for filtering
    minTC      = 0.0                 # Min Total column amount for filtering
    sclfct     = 1.0E12                  # Scale factor to apply to vmr plots (ppmv=1.0E6, ppbv=1.0E9, etc)
    sclfctName = 'pptv'                 # Name of scale factor for labeling plots

    #----------------------
    # Date range to process
    #----------------------
    
    iyear      = 1999
    imnth      = 1
    iday       = 1
    fyear      = 2017
    fmnth      = 7
    fday       = 31

    pCols      = [ [0.0, 4.0], [4.0, 10.0], [10.0, 30.0] ]  #meanTpp = 8.78 ; stdTpp  = 1.14

    #----------------------------------------
    if pCols:
        
        try:

            Filein   = '/data1/projects/ocs/tab/Dates_Thule_dNCEP.ascii'
            
            with open(Filein,'r') as fopen:
                lines       = fopen.readlines()
                YYYY        = np.array([int(x[0:4]) for x in lines[5:]])
                MM          = np.array([int(x[4:6]) for x in lines[5:]])
                DD          = np.array([int(x[6:8]) for x in lines[5:]])
                HH          = np.array([int(x[13:15]) for x in lines[5:]])
                MI          = np.array([int(x[16:18]) for x in lines[5:]])
                SS          = np.array([int(x[19:21]) for x in lines[5:]])
                datesdtp    = [dt.datetime(ye, mo, da, ho, mi, se) for ye, mo, da, ho, mi, se in zip(YYYY, MM, DD, HH, MI, SS)  ] 
                dtp         = np.array([float(x[26:33]) for x in lines[5:]]) 

            datesdtp = np.asarray(datesdtp)
            dtp      = np.asarray(dtp)

            dtpDaily        = mf.dailyAvg(dtp, datesdtp, dateAxis=1, meanAxis=0)
            dtpD            = dtpDaily['dailyAvg']
            datesdtpD       = dtpDaily['dates']

        except:
            print 'Missing Dynamical TH file for tab, a dtp will be assumed'

                                    #----------------------------#
                                    #        --- START ---       #
                                    #----------------------------#

    retDir   = [ '/data1/ebaumer/'+loc.lower()+'/'+g.lower()+'/'+v+'/' for g,v in izip(gasName,ver)] 
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
        statDataCl[gas+'_'+ver[i]] = dc.ReadOutputData(retDir[i],gasName[i],ctlFile[i],iyear,imnth,iday,fyear,fmnth,fday)
        
    #--------------
    # Read profiles
    #--------------
    for gasVer in statDataCl:
        statDataCl[gasVer].readprfs([statDataCl[gasVer].PrimaryGas],retapFlg=1)
        statDataCl[gasVer].readprfs([statDataCl[gasVer].PrimaryGas],retapFlg=0)

        if statDataCl[gasVer].empty: 
            print 'No retreivals found for {}. Exiting.......'.format(gasVer)
            sys.exit()


    rPrfVMR          = OrderedDict()
    rPrfMol          = OrderedDict()
    dates            = OrderedDict()
    alt              = OrderedDict()
    Airmass          = OrderedDict()
    waterVMR         = OrderedDict()
    waterMol         = OrderedDict()
    totClmn          = OrderedDict()
    TCdryAir         = OrderedDict()
    TCdry            = OrderedDict()
    rPrfDry          = OrderedDict()
    rms              = OrderedDict()
    dofsAvg          = OrderedDict()
    dofsAvg_cs       = OrderedDict()
    LayVMR           = OrderedDict()
    DryLayVMR        = OrderedDict()
    upperHgt         = OrderedDict()
    lowerHgt         = OrderedDict()
    LayThk           = OrderedDict()
    LayDOF           = OrderedDict()

    vmrP             = {}    #OrderedDict()
    TCp              = {}    #OrderedDict()

    avkSCF           = OrderedDict()
    avkVMR           = OrderedDict()
    avkSCFav         = OrderedDict()
    avkVMRav         = OrderedDict()

    aPrfVMR          = OrderedDict()
    aPrfMol          = OrderedDict()
           
    rand_err         = OrderedDict()
    sys_err          = OrderedDict()
    tot_err          = OrderedDict()

    rand_errvmr      = OrderedDict()
    sys_errvmr       = OrderedDict()
    tot_errvmr       = OrderedDict()

    rand_cmpnts      = OrderedDict() 
    sys_cmpnts       = OrderedDict()

    rand_cmpnts_vmr  = OrderedDict()
    sys_cmpnts_vmr   = OrderedDict()

    tot_rnd          = OrderedDict()
    tot_sys          = OrderedDict()
    tot_std          = OrderedDict()

    err_summary      = OrderedDict()


    for j, gasVer in enumerate(statDataCl):
        rPrfVMR[gasVer]  = np.asarray(statDataCl[gasVer].rprfs[statDataCl[gasVer].PrimaryGas]) * sclfct
        rPrfMol[gasVer]  = np.asarray(statDataCl[gasVer].rprfs[statDataCl[gasVer].PrimaryGas]  * np.asarray(statDataCl[gasVer].rprfs['AIRMASS']))
        dates[gasVer]    = np.asarray(statDataCl[gasVer].rprfs['date'])
        alt[gasVer]      = np.asarray(statDataCl[gasVer].rprfs['Z'][0,:])
        Airmass[gasVer]  = np.asarray(statDataCl[gasVer].rprfs['AIRMASS'])
        waterVMR[gasVer] = np.asarray(statDataCl[gasVer].aprfs['H2O'])
        waterMol[gasVer] = np.asarray(statDataCl[gasVer].aprfs['H2O'] * Airmass[gasVer])
        totClmn[gasVer]  = np.sum(rPrfMol[gasVer],axis=1)
        TCdryAir[gasVer] = np.sum(Airmass[gasVer],axis=1) - np.sum(waterMol[gasVer],axis=1)
        TCdry[gasVer]    = (totClmn[gasVer] / TCdryAir[gasVer]) * sclfct

        aPrfVMR[gasVer]  = np.asarray(statDataCl[gasVer].aprfs[statDataCl[gasVer].PrimaryGas]) * sclfct
        aPrfMol[gasVer]  = np.asarray(statDataCl[gasVer].aprfs[statDataCl[gasVer].PrimaryGas]  * np.asarray(statDataCl[gasVer].aprfs['AIRMASS']))

        #----------------------------------------
        # This is the mixing ratio for DRY AIR!!!
        #----------------------------------------
        rPrfDry[gasVer] = np.asarray(statDataCl[gasVer].rprfs[statDataCl[gasVer].PrimaryGas]) / (1.0 - waterVMR[gasVer]) * sclfct       
        
        #----------------------------------
        # Read Summary data (For filtering)
        #----------------------------------
        statDataCl[gasVer].readsummary()
        rms[gasVer]     = np.asarray(statDataCl[gasVer].summary[statDataCl[gasVer].PrimaryGas+'_FITRMS'])

        #----------------------------------
        # Read readPbp (observed, fitted, and difference spectra)
        #----------------------------------
        mw    = [str(int(x)) for x in statDataCl[gasVer].ctl['band']]     
        numMW = len(mw)
        statDataCl[gasVer].readPbp()

        #----------------------------------
        # Read Spectra for each gas
        #----------------------------------
        statDataCl[gasVer].readSpectra(statDataCl[gasVer].gasList)

        #-------------------- 
        # Call to filter data
        #--------------------
        if fltrFlg: statDataCl[gasVer].fltrData(statDataCl[gasVer].PrimaryGas,mxrms=maxRMS[j],rmsFlg=rmsFlg, minDOF=minDOF[j],  dofFlg=dofFlg, 
        tcFlg=tcNegFlg, pcFlg=pcNegFlg , cnvrgFlg=cnvrgFlg)
        else:       statDataCl[gasVer].inds = np.array([]) 

        #--------------------------------------------
        # Read Error data to get AVK and profile DOFs
        #-------------------------------------------------
        # Determine if AVK has been created via sfit4 core
        # code or via python error analysis
        #-------------------------------------------------     
        if errorFlg:   # Read AVK from error output
            #statDataCl[gasVer].readError(totFlg=False,sysFlg=False,randFlg=False,vmrFlg=True,avkFlg=True,KbFlg=False)
            statDataCl[gasVer].readError(totFlg=True,sysFlg=True,randFlg=True,vmrFlg=False,avkFlg=True,KbFlg=False)

            npnts           = np.shape(statDataCl[gasVer].error['Total_Random_Error'])[0]
            nlvls           = np.shape(alt[gasVer])[0]
            
            #-------------------------------------------------
            #Error profiles constructed with the diagonal elements of the covariances matrices
            #-------------------------------------------------
            rand_err[gasVer]   = np.zeros((npnts,nlvls))
            sys_err[gasVer]    = np.zeros((npnts,nlvls))

            for i in range(npnts):
                rand_err[gasVer][i,:] = np.diag(statDataCl[gasVer].error['Total_Random_Error'][i][:,:])
                sys_err[gasVer][i,:]  = np.diag(statDataCl[gasVer].error['Total_Systematic_Error'][i][:,:])

            tot_err[gasVer]     = np.sqrt(rand_err[gasVer] + sys_err[gasVer])            
            sys_err[gasVer]     = np.sqrt(sys_err[gasVer])
            rand_err[gasVer]    = np.sqrt(rand_err[gasVer])

            sys_errvmr[gasVer]  = sys_err[gasVer]/ np.asarray(statDataCl[gasVer].rprfs['AIRMASS'])*sclfct
            rand_errvmr[gasVer] = rand_err[gasVer]/ np.asarray(statDataCl[gasVer].rprfs['AIRMASS'])*sclfct
            tot_errvmr[gasVer]  = np.sqrt(rand_errvmr[gasVer] + sys_errvmr[gasVer]) 

            #-------------------------------------------------
            #Error profiles of components
            #-------------------------------------------------
            rand_cmpnts[gasVer] = statDataCl[gasVer].randErrDiag
            sys_cmpnts[gasVer]  = statDataCl[gasVer].sysErrDiag

            rand_cmpnts_vmr[gasVer] = statDataCl[gasVer].randErrDiag 
            sys_cmpnts_vmr[gasVer]  = statDataCl[gasVer].sysErrDiag

            for k in sys_cmpnts_vmr[gasVer]:
                
                try:
                    sys_cmpnts_vmr[gasVer][k] = (np.sqrt(sys_cmpnts_vmr[gasVer][k])/ np.asarray(statDataCl[gasVer].rprfs['AIRMASS']))*sclfct
                except Exception as errmsg:
                    print errmsg
              
            for k in rand_cmpnts_vmr[gasVer]:
                rand_cmpnts_vmr[gasVer][k] = (np.sqrt(rand_cmpnts_vmr[gasVer][k])/ np.asarray(statDataCl[gasVer].rprfs['AIRMASS']))*sclfct

            #-------------------------------------------------
            #Error in the summary output. Get Total Errors 
            #-------------------------------------------------
            tot_rnd[gasVer]        = np.array(statDataCl[gasVer].error['Total random uncertainty'])
            tot_sys[gasVer]        = np.array(statDataCl[gasVer].error['Total systematic uncertainty'])
            tot_std[gasVer]        = np.sqrt(tot_rnd[gasVer]**2 + tot_sys[gasVer]**2)
            
            err_summary[gasVer]    = statDataCl[gasVer].error

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
        nfltr           = len(statDataCl[gasVer].inds)
        rms[gasVer]     = np.delete(rms[gasVer],statDataCl[gasVer].inds)
        ntot            = len(rms[gasVer])
        dates[gasVer]   = np.delete(dates[gasVer],statDataCl[gasVer].inds)
        totClmn[gasVer] = np.delete(totClmn[gasVer],statDataCl[gasVer].inds)
        rPrfVMR[gasVer] = np.delete(rPrfVMR[gasVer],statDataCl[gasVer].inds,axis=0)   
        rPrfMol[gasVer] = np.delete(rPrfMol[gasVer],statDataCl[gasVer].inds,axis=0)   
        rPrfDry[gasVer] = np.delete(rPrfDry[gasVer],statDataCl[gasVer].inds,axis=0)   
        Airmass[gasVer] = np.delete(Airmass[gasVer],statDataCl[gasVer].inds,axis=0)
        TCdry[gasVer]   = np.delete(TCdry[gasVer],statDataCl[gasVer].inds)

        aPrfVMR[gasVer] = np.delete(aPrfVMR[gasVer],statDataCl[gasVer].inds,axis=0)   
        aPrfMol[gasVer] = np.delete(aPrfMol[gasVer],statDataCl[gasVer].inds,axis=0)

        if errorFlg:
            rand_err[gasVer] = np.delete(rand_err[gasVer],statDataCl[gasVer].inds,axis=0)
            sys_err[gasVer]  = np.delete(sys_err[gasVer],statDataCl[gasVer].inds,axis=0)  
            tot_err[gasVer]  = np.delete(tot_err[gasVer],statDataCl[gasVer].inds,axis=0)

            rand_errvmr[gasVer] = np.delete(rand_errvmr[gasVer],statDataCl[gasVer].inds,axis=0)
            sys_errvmr[gasVer]  = np.delete(sys_errvmr[gasVer],statDataCl[gasVer].inds,axis=0)  
            tot_errvmr[gasVer]  = np.delete(tot_errvmr[gasVer],statDataCl[gasVer].inds,axis=0)

            tot_rnd[gasVer] = np.delete(tot_rnd[gasVer],statDataCl[gasVer].inds)
            tot_sys[gasVer] = np.delete(tot_sys,statDataCl[gasVer].inds)
            tot_std[gasVer] = np.delete(tot_std,statDataCl[gasVer].inds)

            
            for k in sys_cmpnts[gasVer]:
                sys_cmpnts[gasVer][k]      = np.delete(sys_cmpnts[gasVer][k],statDataCl[gasVer].inds,axis=0)
                sys_cmpnts_vmr[gasVer][k]  = np.delete(sys_cmpnts_vmr[gasVer][k],statDataCl[gasVer].inds,axis=0)
                
            for k in rand_cmpnts[gasVer]:
                rand_cmpnts[gasVer][k]     = np.delete(rand_cmpnts[gasVer][k],statDataCl[gasVer].inds,axis=0)
                rand_cmpnts_vmr[gasVer][k] = np.delete(rand_cmpnts_vmr[gasVer][k],statDataCl[gasVer].inds,axis=0)

            for k in err_summary[gasVer]:
                err_summary[gasVer][k]     = np.delete(err_summary[gasVer][k],statDataCl[gasVer].inds, axis=0)
          

        #--------------------------------------
        # Remove retrieval data based on filter for spectra fit
        #--------------------------------------
        dataSpec  = OrderedDict()
        gasSpec   = OrderedDict()
        gasAbs    = OrderedDict()
        gasAbsSNR = OrderedDict()
        mwList    = OrderedDict()

        for x in mw:  # Loop through micro-windows
            dataSpec['Obs_'+x]        = np.delete(statDataCl[gasVer].pbp['Obs_'+x],statDataCl[gasVer].inds,axis=0)
            dataSpec['Fitted_'+x]     = np.delete(statDataCl[gasVer].pbp['Fitted_'+x],statDataCl[gasVer].inds,axis=0)
            dataSpec['Difference_'+x] = np.delete(statDataCl[gasVer].pbp['Difference_'+x], statDataCl[gasVer].inds, axis=0)
            dataSpec['WaveN_'+x]      = statDataCl[gasVer].spc['MW_'+x]
            dataSpec['All_'+x]        = np.delete(statDataCl[gasVer].spc['All_'+x],statDataCl[gasVer].inds,axis=0)
            if statDataCl[gasVer].solarFlg: dataSpec['Sol_'+x]        = np.delete(statDataCl[gasVer].spc['Solar_'+x],statDataCl[gasVer].inds,axis=0)

        #----------------------------------------
        # Loop through gases and micro-windows
        # Not all gasses are in each micro-window
        # Determine which gasses belong to which
        # micro-window
        #----------------------------------------
                # Use ordered dictionary so that micro-window plots print in order
        for x in mw: mwList[x] = statDataCl[gasVer].ctl['band.'+x+'.gasb']
        
        gasSpec = {gas.upper()+'_'+m:np.delete(statDataCl[gasVer].spc[gas.upper()+'_'+m],statDataCl[gasVer].inds,axis=0) for m in mwList for gas in mwList[m]}

        #---------------------
        # Calculate Statistics
        #---------------------
        for x in mw:  # Loop through micro-windows

            #-------------------------------------------------------------
            # Calculate the total Observed absorption in micro-window
            # This must be done first because below code modifies dataSpec
            #-------------------------------------------------------------
            gasAbs["Total_"+x] = simps(1.0 - dataSpec['Obs_'+x],x=dataSpec['WaveN_'+x],axis=1)               
                        
            if len(dates[gasVer]) > 1:
                dataSpec['Obs_'+x]        = np.mean(dataSpec['Obs_'+x],axis=0)
                dataSpec['Fitted_'+x]     = np.mean(dataSpec['Fitted_'+x],axis=0)
                dataSpec['Difference_'+x] = np.mean(dataSpec['Difference_'+x],axis=0)
                dataSpec['All_'+x]        = np.mean(dataSpec['All_'+x],axis=0)
                if statDataCl[gasVer].solarFlg:dataSpec['Sol_'+x]        = np.mean(dataSpec['Sol_'+x],axis=0)    
                dataSpec['DifSTD_'+x]     = np.std(dataSpec['Difference_'+x],axis=0)
                
            else:
                dataSpec['Obs_'+x]        = dataSpec['Obs_'+x][0]
                dataSpec['Fitted_'+x]     = dataSpec['Fitted_'+x][0]
                dataSpec['Difference_'+x] = dataSpec['Difference_'+x][0]
                dataSpec['All_'+x]        = dataSpec['All_'+x][0]
                if statDataCl[gasVer].solarFlg:dataSpec['Sol_'+x]        = dataSpec['Sol_'+x][0]
                                
            if len(dates[gasVer]) > 1:

                if statDataCl[gasVer].PrimaryGas+"_"+x in gasSpec:  
                #---------------------------------------------------
                # Calculate the integrate absorption for primary gas
                #---------------------------------------------------               
                    gasAbs[statDataCl[gasVer].PrimaryGas+"_"+x] = simps(1.0 - gasSpec[statDataCl[gasVer].PrimaryGas+"_"+x],x=dataSpec['WaveN_'+x],axis=1)       
                
                #-----------------------------------
                # Calculate the peak absorption of 
                # primary gas for each micro-window
                #-----------------------------------
                    gasAbs[statDataCl[gasVer].PrimaryGas+"_trans_"+x] = 1.0 - np.min(gasSpec[statDataCl[gasVer].PrimaryGas+"_"+x],axis=1)
                
                #---------------------------------------------
                # Determine product of SNR and Peak Absorption
                #---------------------------------------------
                    tempSNR  = np.delete(statDataCl[gasVer].summary["SNR_"+x],statDataCl[gasVer].inds)
                    gasAbsSNR[statDataCl[gasVer].PrimaryGas+"_"+x] = gasAbs[statDataCl[gasVer].PrimaryGas+"_"+x] * tempSNR

        if len(dates[gasVer]) > 1:
            gasSpec = {gas.upper()+'_'+x:np.mean(gasSpec[gas.upper()+'_'+x],axis=0) for x in mwList for gas in mwList[x]}   
        else:
            for x in gasSpec: gasSpec[x] = gasSpec[x][0]

        
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
    #                     PLOTS  
    #-------------------------------------------------
    if saveFlg: pdfsav = PdfPages(pltFile)
   
    clmap        = 'jet'
    cm           = plt.get_cmap(clmap)              
    yearsLc      = YearLocator()
    monthsAll    = MonthLocator()
    #months       = MonthLocator(bymonth=1,bymonthday=1)
    monthsFmt     = MonthLocator()
    DateFmt      = DateFormatter('%Y')


    for j, gasVer in enumerate(statDataCl):

        years = np.asarray([ singDate.year for singDate in dates[gasVer]] )     # Find years for all date entries
        if len(list(set(years))) > 1: yrsFlg = True         # Determine all unique years
        else:                         yrsFlg = False

        #--------------------------------
        # VERTICAL PROFILES (CONTOUR PLOTS)
        #--------------------------------
        indsH = np.where(alt[gasVer] <= 40.0)[0]
        #levels = np.arange(np.round(np.min(rPrfVMR[gasVer][:,indsH]),decimals=0),np.round(np.max(rPrfVMR[gasVer][:,indsH]),decimals=0))
        #levels = np.arange(0, np.round(np.max(rPrfVMR[gasVer][:,indsH]), decimals=1), 0.5)
        levels = np.arange(0, 700., 10)
       
        fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12,6), sharey=False)
       
        gs   = gridspec.GridSpec(1, 2, width_ratios=[3,1])
        ax  = plt.subplot(gs[0])
        ax2  = plt.subplot(gs[1])

        cax          = ax.contourf(dates[gasVer],alt[gasVer][indsH],np.transpose(rPrfVMR[gasVer][:,indsH]),levels, cmap=mplcm.jet)

        #--------------------------------
        #Tropopuase
        #--------------------------------
        #ax.plot(datesdtpD, dtpD, color='gray', linewidth=2)
        #--------------------------------

        divider      = make_axes_locatable(ax)
        cb           = divider.append_axes("right",size="3.5%",pad=0.05)
      
        cbar         = plt.colorbar(cax,cax=cb)
        #cbar.ax.tick_params(labelsize=12)
        cbar.set_label('VMR [{}]'.format(sclfctName), fontsize=14)#, fontsize=14)
        
        ax.set_xlabel('Year', fontsize=14)#, fontsize=14)
        ax.set_ylabel('Altitude [km]', fontsize=14)#, fontsize=14)
        ax.set_ylim(0, 30)
        #ax.set_xlim(dt.date(iyear, imnth, iday), dt.date(fyear, fmnth, fday) )
        ax.set_xlim(np.min(dates[gasVer]) , np.max(dates[gasVer]) )


        ax.yaxis.set_tick_params(which='major',labelsize=14)

        ax.xaxis.set_major_locator(yearsLc)
        ax.xaxis.set_major_formatter(DateFmt)
        ax.xaxis.set_minor_locator(monthsFmt)
        
        #ax.xticks(rotation=45)
        ax.set_title('Time Series of Retrieved {} Vertical Profiles'.format(gasName[0].upper()),multialignment='center')
        ax.xaxis.set_tick_params(which='major',labelsize=14)
        ax.xaxis.set_tick_params(which='minor',labelbottom='off')
        ax.tick_params(labelsize=14)


        for tl in ax.get_xticklabels():
            tl.set_rotation(30)
        #----
        months   = np.asarray([d.month for d in dates[gasVer]])  

        for m in list(set(months)):
            if m == 1: lm = 'Jan'
            if m == 2: lm = 'Feb'
            if m == 3: lm = 'Mar'
            if m == 4: lm = 'Apr'
            if m == 5: lm = 'May'
            if m == 6: lm = 'Jun'
            if m == 7: lm = 'Jul'
            if m == 8: lm = 'Aug'
            if m == 9: lm = 'Sep'
            if m == 10: lm = 'Oct'
            if m == 11: lm = 'Nov'
            if m == 12: lm = 'Dec'
            inds     = np.where(months == m)[0]
            mnthMean = np.mean(rPrfVMR[gasVer][inds,:], axis=0)

            ax2.plot(mnthMean[indsH], alt[gasVer][indsH],label=lm)

        ax2.tick_params(labelsize=14)
        ax2.set_xlabel('VMR [{}]'.format(sclfctName), fontsize=12)#, fontsize=14)
        #ax2.axes.get_yaxis().set_visible(False)
        ax2.set_yticklabels([])
        ax2.grid(True)
        ax2.legend(prop={'size':10})
        ax2.set_ylim(0, 30)
        ax2.set_title('Monthly Mean',multialignment='center')
        ax2.set_xlim(0, 600)
        #ax2.yaxis.set_major_formatter(ax2.NullFormatter())
        #ax2.yticks([])
        
        
        fig.subplots_adjust(left=0.075, bottom=0.15,top=0.95, right=0.975)
        
        if saveFlg: 
            pdfsav.savefig(fig,dpi=200)

        else: 
            plt.show(block=False)
            #user_input = raw_input('Press any key to exit >>> ')
            #sys.exit()

        #-------------------------------------------------
        #Time series of columns 
        #-------------------------------------------------

        if yrsFlg:

            fig, ax1     = plt.subplots(len(pCols), figsize=(9,10), sharex=True)

            yoi     = [[2000, 2008], [2008, 2017]]
            slope   = []
            slope_e = []

            yyyy = [sngdate.year for sngdate in dates[gasVer]]
            yyyy = np.asarray(yyyy)

            for p, pcol in enumerate(pCols):

                Data = vmrP[gasVer+str(p)][0]
                
                #Avg            = mf.dailyAvg(vmrP[gasVer+str(p)[0]], dates[gasVer], dateAxis=1, meanAxis=0)
                #AvgDates        = Avg['dates']
                #dateYearFrac    = mf.toYearFraction(Avg['dates'])
                #AvgData         = Avg['dailyAvg']

                dateYearFrac    = mf.toYearFraction(dates[gasVer])
                weights         = np.ones_like(dateYearFrac)
                           
                ax1[p].scatter(dates[gasVer], Data)
            
                ax1[p].grid(True,which='both', alpha=0.5)
                ax1[p].set_ylabel('VMR [{}]'.format(sclfctName), fontsize=14)
                ax1[p].tick_params(which='both',labelsize=14)
                ax1[p].xaxis.set_major_formatter(DateFmt)
                ax1[p].set_title('{} - {} km'.format(pcol[0], pcol[1]),multialignment='center')

                ax1[p].xaxis.set_major_locator(yearsLc)
                ax1[p].xaxis.set_minor_locator(monthsFmt)
                #ax1.xaxis.set_minor_formatter(DateFormatter('%m'))
                ax1[p].xaxis.set_major_formatter(DateFmt) 
                #ax1.xaxis.set_tick_params(which='major', pad=15)

                for y in yoi:
                    indx1  = np.where( (yyyy >= y[0]) &  (yyyy <= y[1]))[0]

                    if len(indx1) > 1:

                        res    = mf.fit_driftfourier(dateYearFrac[indx1], Data[indx1], weights[indx1], 2, half_period=1.0)
                        f_drift, f_fourier, f_driftfourier,  res_std, A = res[3:8]

                        res_b        = mf.cf_driftfourier(dateYearFrac[indx1], Data[indx1], weights[indx1], 2, half_period=1.0)
                        perc, intercept_b, slope_b, pfourier_b = res_b

                        slope.append(res[1]/np.mean(Data[indx1])*100.0)
                        slope_e.append(np.std(slope_b)/np.mean(Data[indx1])*100.0)

                        ax1[p].plot(dates[gasVer][indx1],f_drift(dateYearFrac[indx1]),label='Fitted Anual Trend + intra-annual variability',linewidth=2.0, color='red')
                        ax1[p].plot(dates[gasVer][indx1],f_driftfourier(dateYearFrac[indx1]),label='Fitted Anual Trend + intra-annual variability',linewidth=2.0, color='green')

                        #ax1[p].text(0.02,0.94,"Fitted trend -- : {0:.3f}$\pm${1:.3f}%)".format(res[1]/np.mean(Data[indx1])*100.0,np.std(slope_b)/np.mean(Data[indx1])*100.0),transform=ax1[p].transAxes)
                        
                        print "Fitted trend -- : {0:.3f}$\pm${1:.3f}%) [{2} - {3}]".format(res[1]/np.mean(Data[indx1])*100.0,np.std(slope_b)/np.mean(Data[indx1])*100.0, y[0], y[1])

                    else:
                        print 'Years out of yoi'

            
            for tl in ax1[-1].get_xticklabels():
                tl.set_rotation(30)

            fig.subplots_adjust(right=0.95, top=0.95, left=0.12, bottom=0.075)
            
            if saveFlg: pdfsav.savefig(fig, dpi=200)
            else:           
                plt.show(block=False)
                #user_input = raw_input('Press any key to exit >>> ')
                #sys.exit() 
                

        #-------------------------------------------------
        #Time series of columns 
        #-------------------------------------------------

        if yrsFlg:

            fig, ax1     = plt.subplots(len(pCols), figsize=(9,10), sharex=True)
            
            DateFmt      = DateFormatter('%b') 
           
            DatesFake    = [dt.datetime(2016 ,d.month, d.day, d.hour, d.minute, d.second) for d in dates[gasVer]]
            DatesFake    = np.asarray(DatesFake)

            cm           = plt.get_cmap('jet', len(set(years)) )

            for p, pcol in enumerate(pCols):
                           
                sc1 = ax1[p].scatter(DatesFake, vmrP[gasVer+str(p)][0], c=years, cmap=cm)
            
                ax1[p].grid(True,which='both')
                ax1[p].set_ylabel('VMR [{}]'.format(sclfctName), fontsize=14)
                ax1[p].tick_params(which='both',labelsize=14)
                ax1[p].xaxis.set_major_formatter(DateFmt)
                ax1[p].set_title('{} - {} km'.format(pcol[0], pcol[1]),multialignment='center')

            ax1[-1].set_xlabel('Month', fontsize=14)

            cax  = fig.add_axes([0.2, 0.95, 0.7, 0.03])

            cbar = fig.colorbar(sc1, cax=cax, format='%4i', orientation='horizontal')                           
            #cbar.set_label('Year', y=2, fontsize=14)

            cbar.set_ticks(np.linspace(np.min(years), np.max(years), np.round(len(set(years))/2. )) )
            cbar.set_ticklabels(np.arange(np.min(years), np.max(years), 2))
            #cbar.set_ticklabels(range(np.min(years), np.max(years)), rotation =45)
            
            fig.subplots_adjust(right=0.95, top=0.9, left=0.12, bottom=0.075)
            
            if saveFlg: pdfsav.savefig(fig, dpi=200)
            else:           
                plt.show(block=False)  
                #user_input = raw_input('Press any key to exit >>> ')
                #sys.exit()

            #-------------------------------------------
            #PLOT OF TIME SERIES ANOMALIES
            #-------------------------------------------

            moi = [4, 5, 6, 7]

            for m in moi:
                
                fig, ax1     = plt.subplots(len(pCols), figsize=(9,10), sharex=True)
               

                for p, pcol in enumerate(pCols):

                    mnthlyVals       = mf.mnthlyAvg(vmrP[gasVer+str(p)][0], dates[gasVer], dateAxis=1, meanAxis=0)
                    dateYearFrac     = mf.toYearFraction(mnthlyVals['dates'])
                    weights          = np.ones_like(dateYearFrac)

                    mntMOI           = np.asarray([ singDate.month for singDate in mnthlyVals['dates']])
                    yrsMOI           = np.asarray([ singDate.year for singDate in mnthlyVals['dates']])
                    indsMonth        = np.where(mntMOI == m)[0]

                    yrsMOI            = yrsMOI[indsMonth]
                    y_pos = np.arange(len(yrsMOI))

                    AnolTS           = mnthlyVals['mnthlyAvg'][indsMonth] - np.nanmean(mnthlyVals['mnthlyAvg'][indsMonth])
     
                    Colors = []
                    #AnolTS = 
                    for val in AnolTS:
                        if val <= 0: Colors.append('blue')# chocolate
                        elif val > 0.0: Colors.append('red')
                
                    ax1[p].bar(y_pos, AnolTS,color=Colors, align='center', edgecolor="black")
                    ax1[p].grid(True)
                    ax1[p].set_ylabel('Anomaly VMR [{}]'.format(sclfctName),multialignment='center', fontsize=14)
                    #ax1[p].set_title('OCS Anomaly [Month={}]'.format(5))
                    ax1[p].set_title('{} - {} km'.format(pcol[0], pcol[1]),multialignment='center')
                    ax1[p].set_xlim(y_pos[0]-1,y_pos[-1]+1)
                    #ax1[p].set_ylim(-1.0, 1.5)
                    ax1[p].tick_params(which='both',labelsize=14)

                ax1[-1].set_xlabel('Year', fontsize=14)

                plt.xticks(y_pos, yrsMOI, rotation=45)
                plt.suptitle('OCS Anomalies [Month={}]'.format(m), fontsize=14)
                fig.subplots_adjust(right=0.97, top=0.9, left=0.12, bottom=0.1)


                if saveFlg: pdfsav.savefig(fig, dpi=200)
                else:           
                    plt.show(block=False)  
            
            #user_input = raw_input('Press any key to exit >>> ')
            #sys.exit()


        #-------------------------------------------------
        #Fit: Mean  
        #-------------------------------------------------
        fig = plt.figure(figsize=(16,9))
            
        for p, x in enumerate(mwList):
            gs1 = gridspec.GridSpec(3, 2)

            if p == 0: gs1.update(left=0.06, right=0.48, top=0.95, bottom=0.55,  wspace=0.05,  hspace=0.08)
            if p == 1: gs1.update(left=0.06, right=0.48, top = 0.45, bottom=0.1, wspace=0.05, hspace=0.08)
            if p == 2: gs1.update(left=0.55, right=0.98, top=0.95, bottom=0.55,  wspace=0.05,  hspace=0.08)
            if p == 3: gs1.update(left=0.55, right=0.98, top=0.45, bottom=0.1, wspace=0.05,   hspace=0.08)

            #ax1 = fig.add_subplot(gs1[0])
            #ax2 = fig.add_subplot(gs1[1], sharex=ax1)

            ax1 = plt.subplot(gs1[0:1, :])
            ax2 = plt.subplot(gs1[1:3, :])

           
            ax1.plot(dataSpec['WaveN_'+x],dataSpec['Difference_'+x]*100.0,  color='k')
            ax1.axhline(y=0,color='r')
            if len(dates[gasVer]) > 1: ax1.fill_between(dataSpec['WaveN_'+x],dataSpec['Difference_'+x]*100.0-dataSpec['DifSTD_'+x]*100.0,dataSpec['Difference_'+x]*100.0+dataSpec['DifSTD_'+x]*100.0,alpha=0.5,color='0.75')
            ax1.grid(True)
            if (p == 0) or (p == 1): ax1.set_ylabel('% Difference', fontsize=14,)
            ax1.set_title('Micro-Window '+ x)
            ax1.set_xlim((np.min(dataSpec['WaveN_'+x]),np.max(dataSpec['WaveN_'+x])))
            
            plt.tick_params(which='minor',length=4,color='b')
            
            ax2.plot(dataSpec['WaveN_'+x],dataSpec['Obs_'+x],label='Observed')
            ax2.plot(dataSpec['WaveN_'+x],dataSpec['Fitted_'+x],label='Fitted')
            if statDataCl[gasVer].solarFlg: ax2.plot(dataSpec['WaveN_'+x],dataSpec['Sol_'+x],label='Solar')
            sclfct = 0.0
            for g in mwList[x]: 
                sclfct += 0.02
                
                ax2.plot(dataSpec['WaveN_'+x],gasSpec[g.upper()+'_'+x]+(gasSpec[g.upper()+'_'+x]*sclfct),label=g.upper())
            
            ax2.grid(True)
            if (p == 1) or (p == 3): ax2.set_xlabel('Wavenumber [cm$^{-1}$]', fontsize=14)
            if (p == 0) or (p == 1): ax2.set_ylabel('Transmission', fontsize=14)
            #ax2.set_ylim(bottom=0.0)
            #if (p == 0) or (p == 1):ax2.set_xlim((np.min(dataSpec['WaveN_'+x]),np.max(dataSpec['WaveN_'+x])))
            ax2.set_xlim((np.min(dataSpec['WaveN_'+x]),np.max(dataSpec['WaveN_'+x])))

            ax2.legend(prop={'size':10.5},loc='upper center', bbox_to_anchor=(0.5, 1.065),
                      fancybox=True, ncol=len(mwList[x])+3)  
            
            ax1.tick_params(axis='x',which='both',labelsize=14, bottom='off', labelbottom='off')
            ax1.tick_params(axis='y',which='both',labelsize=14)
            ax2.tick_params(axis='x',which='both',labelsize=14)
            ax2.tick_params(axis='y',which='both',labelsize=14)

            ax1.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
            ax1.xaxis.set_minor_locator(AutoMinorLocator())

            ax2.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
            ax2.xaxis.set_minor_locator(AutoMinorLocator())

            ax1.set_ylim(-2,2)
            ax2.set_ylim(0,1.3)

            #major_ticks = np.arange(-2, 2, 1.0)
            #ax1.set_xticks(major_ticks)

        if saveFlg: 
            pdfsav.savefig(fig,dpi=200)
            #plt.savefig(FigDir+'FitMean1_'+loc.upper()+'.pdf', bbox_inches='tight')
        else:       plt.show(block=False)

        #-------------------------------------------------
        #Fit: Mean   
        #-------------------------------------------------
        fig = plt.figure(figsize=(16,9))
            
        for p, x in enumerate(mwList):
            gs1 = gridspec.GridSpec(1, 2)

            if p == 0: gs1.update(left=0.06, right=0.48, top=0.95, bottom=0.55, wspace=0.05, hspace=0.08)
            if p == 1: gs1.update(left=0.06, right=0.48, top=0.45, bottom=0.1,  wspace=0.05, hspace=0.08)
            if p == 2: gs1.update(left=0.55, right=0.98, top=0.95, bottom=0.55, wspace=0.05, hspace=0.08)
            if p == 3: gs1.update(left=0.55, right=0.98, top=0.45, bottom=0.1,  wspace=0.05, hspace=0.08)

            ax2 = plt.subplot(gs1[0:1, :])

            ax2.plot(dataSpec['WaveN_'+x],dataSpec['Difference_'+x]*10.0,  color='k')
            ax2.axhline(y=0,color='gray', linestyle='-')
            if len(dates[gasVer]) > 1: ax2.fill_between(dataSpec['WaveN_'+x],dataSpec['Difference_'+x]*10.0-dataSpec['DifSTD_'+x]*10.0,dataSpec['Difference_'+x]*10.0+dataSpec['DifSTD_'+x]*10.0,alpha=0.5,color='0.75')
            
            plt.tick_params(which='minor',length=4,color='b')
            
            ax2.plot(dataSpec['WaveN_'+x],dataSpec['Obs_'+x],label='Observed')
            ax2.plot(dataSpec['WaveN_'+x],dataSpec['Fitted_'+x],label='Fitted')
            #if statDataCl[gasVer].solarFlg: ax2.plot(dataSpec['WaveN_'+x],dataSpec['Sol_'+x],label='Solar')
            sclfct = 0.0
            #for g in mwList[x]: 
            #    sclfct += 0.02   
            #    ax2.plot(dataSpec['WaveN_'+x],gasSpec[g.upper()+'_'+x]+(gasSpec[g.upper()+'_'+x]*sclfct),label=g.upper())
            
            ax2.grid(True)
            if (p == 1) or (p == 3): ax2.set_xlabel('Wavenumber [cm$^{-1}$]', fontsize=14)
            if (p == 0) or (p == 1): ax2.set_ylabel('Transmission', fontsize=14)
            if (p == 0) or (p == 1):ax2.set_xlim((np.min(dataSpec['WaveN_'+x]),np.max(dataSpec['WaveN_'+x])))

            #ax2.legend(prop={'size':10.5},loc='upper center', bbox_to_anchor=(0.5, 1.065),
            #          fancybox=True, ncol=len(mwList[x])+3)  
            
            ax2.tick_params(axis='x',which='both',labelsize=14)
            ax2.tick_params(axis='y',which='both',labelsize=14)

            ax2.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
            ax2.xaxis.set_minor_locator(AutoMinorLocator())

            ax2.annotate('x10', xy=(0.65, 0.12),size=16, xycoords='axes fraction', xytext=(0.7, 0.25), textcoords='axes fraction', arrowprops=dict(facecolor='black', shrink=0.05))

            ax2.set_ylim(-0.22,1.075)
            ax2.set_xlim((np.min(dataSpec['WaveN_'+x]),np.max(dataSpec['WaveN_'+x])))

        if saveFlg: 
            pdfsav.savefig(fig,dpi=200)
            #plt.savefig(FigDir+'FitMean2_'+loc.upper()+'.pdf', bbox_inches='tight')
        else:       
            plt.show(block=False)
            #user_input = raw_input('Press any key to exit >>> ')
            #sys.exit()


        #-----------------------------
        #ERRORS:  
        #-----------------------------
        if errorFlg:
            print '\nMean Total Column = {}'.format(np.mean(totClmn[gasVer]))
            #-----------------------------
            #ERRORS: Calculate mean errors
            #-----------------------------
            #fig, ax1 = plt.subplots()
            for k in sorted(err_summary[gasVer].iterkeys()):
                try:
                    cmpnts_error_i = np.mean(err_summary[gasVer][k])
                    print '{:<45s} = {:15.5E} '.format(k, cmpnts_error_i)
                    print '{:<45s} = {:10.4f} %'.format(k, cmpnts_error_i/np.mean(totClmn[gasVer]) * 100.)
                    #ax1.plot(dates[gasVer],cmpnts_error[gasVer][k], label=k)
                    #if k == 'Total systematic uncertainty lineint_h2o':
             #           ax1.scatter(dates[gasVer], err_summary[gasVer][k], color='r', label=k)
                    #elif k == 'Total systematic uncertainty linepair_h2o':
              #          ax1.scatter(dates[gasVer], err_summary[gasVer][k], color='blue', label=k)
                except Exception as errmsg:
                    print errmsg
            print '\n'
            #ax1.legend(prop={'size':8},loc='upper right') 
            #plt.show(block=False)
            #user_input = raw_input('Press any key to exit >>> ')
            #sys.exit()

            #-----------------------------
            #Plot Random error components
            #-----------------------------
            
            fig, ax1 = plt.subplots(1, 2, figsize=(8,7), sharey=True, sharex=True)

            #for k, p in enumerate(list_comp_ran):
            for p in sorted(rand_cmpnts_vmr[gasVer].iterkeys()):

                if npnts > 1:
                    errPlt = np.mean(rand_cmpnts_vmr[gasVer][p],axis=0)
                    retPrf = np.mean(rPrfVMR[gasVer],axis=0)
                else:
                    errPlt = rand_cmpnts_vmr[gasVer][p][0]
                    retPrf = rPrfVMR[gasVer][0]

                errPlt = errPlt / retPrf *100.

                if p == 'measurement_Random_Error': legend = 'Measurement'
                elif p == 'interfering_species_Random_Error': legend = 'Interfering species'
                elif p == 'temperature_Random_Error': legend = 'Temperature'
                elif p == 'sza_Random_Error': legend = 'SZA'
                elif p == 'retrieval_parameters_Random_Error': legend = 'Retrieval parameters'
                elif p == 'curvature_Random_Error': continue
                elif p == 'max_opd_Random_Error': continue     
                else: legend = p
            
                ax1[0].plot(errPlt,alt[gasVer],linewidth=2.0, label=legend)

            ax1[0].plot(np.nanmean(rand_errvmr[gasVer], axis=0)/ retPrf *100.,alt[gasVer],linestyle='--', color='k',linewidth=2.0, label='Total')
            ax1[0].set_ylabel('Altitude [km]', fontsize=14)
            #ax1[0].set_xlabel('VMR [ppm]', fontsize=14)
            ax1[0].set_xlabel('Error [%]', fontsize=14)             
            ax1[0].grid(True,which='both')
            ax1[0].legend(prop={'size':12},loc='upper right')          
            #ax2[k,0].tick_params(axis='both',which='both',labelsize=8) 
            ax1[0].set_title('Random Components')
            ax1[0].set_ylim(0, 60)
            ax1[0].tick_params(which='both',labelsize=14)
            #ax1[0].set_xlim(0, 10)

            #-----------------------------
            #ERRORS:   Plot systematic error components
            #-----------------------------

            for p in sorted(sys_cmpnts_vmr[gasVer].iterkeys()):

                if npnts > 1:
                    errPlt = np.mean(sys_cmpnts_vmr[gasVer][p],axis=0)
                    retPrf = np.mean(rPrfVMR[gasVer],axis=0)
                else:
                    errPlt = sys_cmpnts_vmr[gasVer][p][0]
                    retPrf = rPrfVMR[gasVer][0]

                errPlt = errPlt / retPrf *100.
                
                if p == 'temperature_Systematic_Error': legend = 'Temperature'
                elif p == 'smoothing_Systematic_Error': continue # legend = 'Smoothing'
                #elif p == 'smoothing_Systematic_Error': legend = 'Smoothing'            
                elif p == 'linepair_'+statDataCl[gasVer].PrimaryGas.lower()+'_Systematic_Error': legend = '$\gamma$'
                elif p == 'lineint_'+statDataCl[gasVer].PrimaryGas.lower()+'_Systematic_Error': legend = '$S$'
                elif p == 'linetair_'+statDataCl[gasVer].PrimaryGas.lower()+'_Systematic_Error': legend = '$n$'
                elif p == 'phase_fcn_Systematic_Error': legend = 'Phase fcn'
                elif p == 'apod_fcn_Systematic_Error':  legend = 'Apod fcn'
                else: legend = p

                ax1[1].plot(errPlt, alt[gasVer],linewidth=2.0, label=legend)

            ax1[1].plot(np.nanmean(sys_errvmr[gasVer], axis=0)/ retPrf *100., alt[gasVer],linestyle='--', color='k',linewidth=2.0, label='Total')
            #ax1[1].set_xlabel('VMR [ppm]', fontsize=14)
            ax1[1].set_xlabel('Error [%]', fontsize=14)             
            ax1[1].grid(True,which='both')
            ax1[1].legend(prop={'size':12}, loc='upper right')
            ax1[1].set_title('Systematic Components')
            ax1[1].set_ylim(0, 60)
            ax1[1].tick_params(which='both',labelsize=14)
            ax1[1].set_xlim(0, 15)
                #ax1[1].set_xscale('log')

            fig.subplots_adjust(bottom=0.1, top=0.95, left = 0.1, right = 0.95)
           
            if saveFlg: 
                pdfsav.savefig(fig,dpi=200)
                #plt.savefig(FigDir+'UncertaintyPrf_'+loc.upper()+'.pdf', bbox_inches='tight')
            else:       
                plt.show(block=False)
                #user_input = raw_input('Press any key to exit >>> ')
                #sys.exit()


            


        #---------------------------------
        # Plot : Averaging Kernel Smoothing Function (row of avk)
        #---------------------------------
        clmap        = 'jet'
        cm           = plt.get_cmap(clmap)
        fig       = plt.figure(figsize=(9,7))
        gs        = gridspec.GridSpec(1,3,width_ratios=[3,1,1])
        ax        = plt.subplot(gs[0])
        axb       = plt.subplot(gs[1])
        axc       = plt.subplot(gs[2])
        cm        = plt.get_cmap(clmap)
        cNorm     = colors.Normalize(vmin=np.min(alt[gasVer]), vmax=np.max(alt[gasVer]))
        scalarMap = mplcm.ScalarMappable(norm=cNorm,cmap=clmap)
        scalarMap.set_array(alt[gasVer])

        #---------------------------------
        ax.set_color_cycle([scalarMap.to_rgba(x) for x in alt[gasVer]])
        
        for i in range(len(alt[gasVer])):
            ax.plot(avkSCFav[gasVer][i,:],alt[gasVer])
            
        ax.set_ylabel('Altitude [km]', fontsize=14)
        ax.set_xlabel('AK', fontsize=14)
        ax.grid(True, alpha=0.5)
        ax.set_title('(a)', loc='left', fontsize=14)

        cbaxes = fig.add_axes([0.4, 0.55, 0.02, 0.35]) 
        cbar = fig.colorbar(scalarMap, orientation='vertical', cax = cbaxes)
        cbar.set_label('Altitude [km]', fontsize=14)
        #ax.set_title('H2O Averaging Kernels Scale Factor', fontsize=14)
        ax.tick_params(labelsize=14)
        #ax.set_ylim(1, 15)
        if loc.lower() == 'fl0': ax.set_ylim(1, 60)
        if loc.lower() == 'mlo': ax.set_ylim(3, 60)
        if loc.lower() == 'tab': ax.set_ylim(0, 60)

        #---------------------------------        
        axb.plot(np.sum(avkSCFav[gasVer],axis=0), alt[gasVer],color='k')
        axb.grid(True,  alpha=0.5)
        axb.set_xlabel('AK Area', fontsize=14)
        axb.tick_params(labelsize=14)
        axb.set_title('(b)', loc='left', fontsize=14)

        major_ticks = np.arange(0, 3, 1)
        axb.set_xticks(major_ticks) 
        #axb.set_ylim(1, 15)

        if loc.lower() == 'fl0': axb.set_ylim(1, 60)
        if loc.lower() == 'mlo': axb.set_ylim(3, 60)
        if loc.lower() == 'tab': axb.set_ylim(0, 60) 

        #---------------------------------
        dofs_cs = np.cumsum(np.diag(avkSCFav[gasVer])[::-1])[::-1]
        axc.plot(dofs_cs,alt[gasVer],color='k',label='Cumulative Sum of DOFS (starting at surface)')
        xval = range(0,int(np.ceil(max(dofs_cs)))+2)

        #ind1         = mf.nearestind(Pcol[0], alt[gasVer])
        #ind2         = mf.nearestind(Pcol[1], alt[gasVer])

        #axc.fill_between(xval,alt[gasVer][ind1],alt[gasVer][ind2],alpha=0.5,color='0.75')  
        #axc.axhline(alt[gasVer][ind2],color='k',linestyle='--')
        #dofsPcol = dofs_cs[ind2] - dofs_cs[ind1]
        #axc.text(0.15,(alt[idhdf][ind1]+alt[idhdf][ind2])/2.0, 
        #         'DOFs for layer {0:.2f}-{1:.2f}[km] = {2:.3f}'.format(alt[idhdf][ind1],alt[idhdf][ind2],dofsPcol),
        #         fontsize=9)
        #axc.set_title('DOFs Profile - ' +str(pltID[i]))
        #axc.set_ylabel('Altitude [km]')
        
        axc.set_xlabel('Cumulative\nSum of DOFS', fontsize=14)  
        axc.tick_params(labelsize=14)
        axc.set_title('(c)', loc='left', fontsize=14)
        #axc.set_title('DOFs for layer {0:.1f}-{1:.1f}[km] = {2:.2f}'.format(alt[idhdf][ind1],alt[idhdf][ind2],dofsPcol), fontsize=9)    
        axc.set_ylim((0,120))
        axc.grid(True,which='both',  alpha=0.5)
        if loc.lower() == 'fl0': axc.set_ylim(1, 60)
        if loc.lower() == 'mlo': axc.set_ylim(3, 60)
        if loc.lower() == 'tab': axc.set_ylim(0, 60)

        plt.tight_layout()

        if saveFlg: 
            pdfsav.savefig(fig,dpi=200)
            #plt.savefig(FigDir+'avkSCFavPrf_'+loc.upper()+'.pdf', bbox_inches='tight')
        else:       plt.show(block=False)

        #---------------------------------
        # Plot : VMR Averaging Kernel Smoothing Function (row of avk)
        #---------------------------------
        clmap        = 'jet'
        cm           = plt.get_cmap(clmap)
        fig       = plt.figure(figsize=(7,7))
        gs        = gridspec.GridSpec(1,2,width_ratios=[3,1])
        ax        = plt.subplot(gs[0])
        axb       = plt.subplot(gs[1])
        cm        = plt.get_cmap(clmap)
        cNorm     = colors.Normalize(vmin=np.min(alt[gasVer]), vmax=np.max(alt[gasVer]))
        scalarMap = mplcm.ScalarMappable(norm=cNorm,cmap=clmap)
        scalarMap.set_array(alt[gasVer])
        ax.set_color_cycle([scalarMap.to_rgba(x) for x in alt[gasVer]])
        
        for i in range(len(alt[gasVer])):
            ax.plot(avkVMRav[gasVer][i,:],alt[gasVer])
            
        ax.set_ylabel('Altitude [km]', fontsize=14)
        ax.set_xlabel('Averaging Kernels', fontsize=14)
        ax.grid(True)
        ax.set_title('H2O Averaging Kernels [VMR]', fontsize=14)
        ax.tick_params(labelsize=14)

        if loc.lower() == 'fl0': ax.set_ylim(1, 60)
        if loc.lower() == 'mlo': ax.set_ylim(3, 60)
        if loc.lower() == 'tab': ax.set_ylim(0, 60)
        
        axb.plot(np.sum(avkVMRav[gasVer],axis=0), alt[gasVer],color='k')
        axb.grid(True)
        axb.set_xlabel('Averaging Kernel Area', fontsize=14)
        axb.tick_params(axis='x',which='both',labelsize=8)

        if loc.lower() == 'fl0': axb.set_ylim(1, 60)
        if loc.lower() == 'mlo': axb.set_ylim(3, 60)
        if loc.lower() == 'tab': axb.set_ylim(0, 60)
        #axb.tick_params(labelsize=14) 

        cbar = fig.colorbar(scalarMap, orientation='vertical')
        cbar.set_label('Altitude [km]', fontsize=14)

        if saveFlg: 
            pdfsav.savefig(fig,dpi=200)
            #plt.savefig(FigDir+'avkVMRavPrf_'+loc.upper()+'.pdf', bbox_inches='tight')
        else:       plt.show(block=False)

        #---------------------------------
        # Plot :Partial Columns
        #---------------------------------

        #fig = plt.figure(figsize=(16,12))

        for p, pcol in enumerate(pCols):

            #----------------------------------------------------------------------
            #analysis of time series of DAILY averages
            #----------------------------------------------------------------------
            dailyVals        = mf.dailyAvg(vmrP[gasVer+str(p)][0], dates[gasVer],dateAxis=1, meanAxis=0)
            #mnthlyVals       = mf.mnthlyAvg(vmrP[gasVer+str(p)][0],dates[gasVer],dateAxis=1, meanAxis=0)
            dateYearFrac     = mf.toYearFraction(dailyVals['dates'])
            #dateYearFrac     = mf.toYearFraction(mnthlyVals['dates'])
            weights          = np.ones_like(dateYearFrac)

            res              = mf.fit_driftfourier(dateYearFrac, dailyVals['dailyAvg'], weights, 2)
            #res          = mf.fit_driftfourier(dateYearFrac, mnthlyVals['mnthlyAvg'], weights, 2)
            intercept, slope, pfourier = res[0:3]
            f_drift, f_fourier, f_driftfourier = res[3:6]

            DT_d             = dailyVals['dates']
            #DT_d             = mnthlyVals['dates']
            vmrP_d           = dailyVals['dailyAvg']
            #vmrP_d           = mnthlyVals['mnthlyAvg']
            
            TCsd_d           = dailyVals['std']
            #TCsd_d           = mnthlyVals['std']

            int_fourier_d    = intercept
            sl_fourier_d     = slope
            p_fourier_d      = pfourier
            FAT_d            = f_drift(dateYearFrac)
            FATAV_d          = f_driftfourier(dateYearFrac)

            perc, intercept_boot, slope_boot, pfourier_boot = mf.cf_driftfourier(dateYearFrac, dailyVals['dailyAvg'], weights, 2)
            #perc, intercept_boot, slope_boot, pfourier_boot = mf.cf_driftfourier(dateYearFrac, mnthlyVals['mnthlyAvg'], weights, 2)
            int_boot_d       = intercept_boot
            sl_boot_d        = slope_boot
            p_boot_d         = pfourier_boot

            Rate_d           = np.divide(res[1], np.mean(dailyVals['dailyAvg']))  *100.0 
            Rate_e_d         = np.divide(np.std(slope_boot), np.mean(dailyVals['dailyAvg']))*100.0

            #Rate_d           = np.divide(res[1], np.mean(mnthlyVals['mnthlyAvg']))  *100.0 
            #Rate_e_d         = np.divide(np.std(slope_boot), np.mean(mnthlyVals['mnthlyAvg']))*100.0

            #----------------------------------------------------------------------
            #analysis of time series of Sonde
            #----------------------------------------------------------------------

           #  datesCFH          = np.asarray(sondedates[p])
           #  vmrPCFH           = np.asarray(sondevmrP[p])

           #  if p == 5:
           #      indmax = np.where(vmrPCFH > 150.)[0]
           #      datesCFH       =  np.delete(datesCFH, indmax, axis=0)
           #      vmrPCFH        =  np.delete(vmrPCFH, indmax, axis=0)

           #  mnthlyValsCFH       = mf.mnthlyAvg(vmrPCFH,datesCFH, dateAxis=1, meanAxis=0)
           
           #  dateYearFracCFH   = mf.toYearFraction(mnthlyValsCFH['dates'])
           #  weights           = np.ones_like(dateYearFracCFH)
           #  resCFH            = mf.fit_driftfourier(dateYearFracCFH, mnthlyValsCFH['mnthlyAvg'], weights, 2)

           #  interceptCFH, slopeCFH, pfourierCFH          = resCFH[0:3]
           #  f_driftCFH, f_fourierCFH, f_driftfourierCFH  = resCFH[3:6]

           #  FAT_CFH            = f_drift(dateYearFracCFH)
           #  FATAV_CFH          = f_driftfourier(dateYearFracCFH)

           #  perc, intercept_boot, slope_boot, pfourier_boot = mf.cf_driftfourier(dateYearFracCFH, mnthlyValsCFH['mnthlyAvg'], weights, 2)
           #  int_boot_CFH       = intercept_boot
           #  sl_boot_CFH        = slope_boot
           #  p_boot_CFH         = pfourier_boot

           #  Rate_CFH           = np.divide(resCFH[1], np.mean(mnthlyValsCFH['mnthlyAvg']))  *100.0 
           #  Rate_e_CFH         = np.divide(np.std(slope_boot), np.mean(mnthlyValsCFH['mnthlyAvg']))*100.0

           #  #----------------------------------------------------------------------
           #  gs1 = gridspec.GridSpec(1, 3)

           #  # if p == 0: gs1.update(left=0.075, right=0.5, top=0.95, bottom=0.525,  wspace=0.05,  hspace=0.08)
           #  # if p == 2: gs1.update(left=0.075, right=0.5, top = 0.475, bottom=0.1, wspace=0.05, hspace=0.08)
           
           #  # if p == 1: gs1.update(left=0.555, right=0.98, top=0.95, bottom=0.525,  wspace=0.05,  hspace=0.08)
           #  # if p == 3: gs1.update(left=0.555, right=0.98, top=0.475, bottom=0.1, wspace=0.05,   hspace=0.08)

           #  if p == 0: gs1.update(left=0.075, right=0.5, top=0.95, bottom=0.7,  wspace=0.05,  hspace=0.08)
           #  if p == 2: gs1.update(left=0.075, right=0.5, top = 0.65, bottom=0.4, wspace=0.05, hspace=0.08)
           #  if p == 4: gs1.update(left=0.075, right=0.5, top = 0.35, bottom=0.1, wspace=0.05, hspace=0.08)
           
           #  if p == 1: gs1.update(left=0.555, right=0.98, top=0.95, bottom=0.7,  wspace=0.05,  hspace=0.08)
           #  if p == 3: gs1.update(left=0.555, right=0.98, top=0.65, bottom=0.4, wspace=0.05,   hspace=0.08)
           #  if p == 5: gs1.update(left=0.555, right=0.98, top=0.35, bottom=0.1, wspace=0.05,   hspace=0.08)

           #  ax2 = plt.subplot(gs1[0:1, :])

        
           #  ax2.scatter(DT_d,vmrP_d,facecolors='red', edgecolors='black', s=35, label='Dayly averages')
           #  ax2.plot(DT_d,FAT_d, color='blue', linewidth=2.5)
           #  ax2.plot(DT_d,FATAV_d, color='green', linewidth=2.5)

           #  ax2.scatter(mnthlyValsCFH['dates'],mnthlyValsCFH['mnthlyAvg'],facecolors='blue', edgecolors='black', s=35, label='CFH')


           #  ax2.xaxis.set_major_locator(yearsLc)
           #  ax2.xaxis.set_minor_locator(monthsFmt)
           #  ax2.xaxis.set_major_formatter(DateFmt)
           #  ax2.xaxis.set_tick_params(which='major',labelsize=18)
           #  ax2.xaxis.set_tick_params(which='minor',labelbottom='off')
           
           #  ax2.grid(True)
           #  if (p == 1) or (p == 3): ax2.set_xlabel('Year', fontsize=14)
           #  if (p == 0) or (p == 1): ax2.set_ylabel('VMR [ppm]', fontsize=14)

           #  ax2.set_ylim(ymin=0)
           #  ax2.tick_params(axis='x',which='both',labelsize=14)
           #  ax2.tick_params(axis='y',which='both',labelsize=14)

           #  ax2.text(0.02,0.94,"Layer: {0:.1f} - {1:.1f} km".format(pcol[0], pcol[1]),transform=ax2.transAxes)
           #  ax2.text(0.02,0.88,"Fitted trend -- slope: {0:.3E} ({1:.3f}%)".format(res[1],res[1]/np.mean(dailyVals['dailyAvg'])*100.0),transform=ax2.transAxes)            
           #  #ax2.text(0.02,0.88,"Fitted trend -- slope: {0:.3E} ({1:.3f}%)".format(res[1],res[1]/np.mean(mnthlyVals['mnthlyAvg'])*100.0),transform=ax2.transAxes)
           #  ax2.text(0.02,0.82,"Fitted trend -- slope: {0:.3E} ({1:.3f}%)".format(resCFH[1],resCFH[1]/np.mean(mnthlyValsCFH['mnthlyAvg'])*100.0),transform=ax2.transAxes)
           
           # # ax2.text(0.02,0.89,"Fitted intercept at xmin: {:.3E}".format(res[0]),transform=ax2.transAxes)
           #  #ax2.text(0.02,0.84,"STD of residuals: {0:.3E} ({1:.3f}%)".format(res[6],res[6]/np.mean(dailyVals['dailyAvg'])*100.0),transform=ax2.transAxes) 
            


        #if saveFlg: 
        #    pdfsav.savefig(fig,dpi=200)
        #    plt.savefig(FigDir+'TSPC_'+loc.upper()+'.pdf', bbox_inches='tight')

        #else:       plt.show(block=False)


    if saveFlg: pdfsav.close()
    else:
        plt.show(block=False)
        user_input = raw_input('Press any key to exit >>> ')
        sys.exit()


 
if __name__ == "__main__":
    main()