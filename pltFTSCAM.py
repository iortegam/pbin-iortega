#! /usr/bin/python2.7
#----------------------------------------------------------------------------------------
# Name:
#        pltFTSCAM.py
#
# Purpose:
#       1) Read and plot netCDF files From CAM-CHEM
#       2) Integrate Model outputs and FTS
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

from netCDF4 import Dataset
import dataModelOutClass as dm
import dataOutClass as dc
from collections                     import OrderedDict

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
    #                     Global Inputs independently of FTS and/or Models
    #-----------------------------------------------------------------------------------------

    loc               = 'fl0'
    pCols             = [1.6, 8.0]       #--ALTITUDE TO CALCULATE PARTIAL COLUMNS AND WEIGHTED VMR
    
    pltAnalysis       = True
    saveFlg           = True 
    
    pltFile           =  '/data/iortega/results/'+loc.lower()+'/fig/FTS-Models_'+loc.upper()+'.pdf'


    #-----------------------------------------------------------------------------------------
    #                 Initializations for CAM-CHEM (2000 - 2015)
    #-----------------------------------------------------------------------------------------
    dataDirCAM        = '/data1/ancillary_data/'+loc.lower()+'/'
    fileCAM           = 'CAM_chem_fmerra_FSDSSOA_2deg_2000_2014_extra_Boulder.nc'


    gasNameCAM          = ['c2h6', 'co', 'ch2o', 'o3', 'c2h2', 'ch4', 'hcooh' ]           

    ReadCAM           = True
    pltCAM            = False
    saveCAMFlg        = False
    
    pltCAMFile        = '/data/iortega/results/'+loc.lower()+'/fig/CAM-Chem_FL0.pdf'

    #-----------------------------------------------------------------------------------------
    #                 Initializations for FTIR
    #-----------------------------------------------------------------------------------------
    gasName           = ['co', 'c2h2',  'c2h6', 'h2co', 'hcooh', 'o3']#, 'ch4']                   # Name of gas

    ver               = ['Current_v3', 'Current_v2',   'Current_v2', 'Current_WP_v6',  'Current_v1', 'Current_WP']#, 'Current_WP']          # Name of retrieval version to process
    ctlF              = ['sfit4_v3.ctl', 'sfit4_v2.ctl', 'sfit4_v2.ctl', 'sfit4_v6.ctl', 'sfit4_v1.ctl', 'sfit4.ctl']#, 'sfit4_v3.ctl']            # Name of ctl file

    #------
    # Flags - 1
    #------
    ReadFTS            = True   

    #------
    # Flags - 2
    #------
    errorFlg           = False                   # Flag to process error data
    fltrFlg            = True                   # Flag to filter the data
    byYrFlg            = False                  # Flag to create plots for each individual year in date range
    szaFlg             = True                   # Flag to filter based on min and max SZA
    dofFlg             = True                   # Flag to filter based on min DOFs
    pcNegFlg           = True                  # Flag to filter profiles with negative partial columns
    tcNegFlg           = True                  # Flagsag to filter profiles with negative total columns
    tcMMFlg            = False                  # Flag to filter based on min and max total column amount
    cnvrgFlg           = True                   # Flag to filter profiles that did not converge
    rmsFlg             = True                   # Flag to filter based on max RMS
    chiFlg             = False                  # Flag to filter based on max CHI_2_Y
    mnthFlg            = False                  # Flag to filter based on 

    mnths              = [6,7,8]                # Months to filter on (these are the months to include data)
    maxRMS             = [2.0, 0.65, 1.0, 1.0, 3.0, 4.0, 1.0]                      # Max Fit RMS to filter data. Data is filtered according to <= maxrms
    minDOF             = [1.0, 0.5, 0.5, 0.5, 0.0, 2.0, 0.0]                      # Min DOFs for filtering
    minSZA             = 0.0                    # Min SZA for filtering
    maxSZA             = 90.0                   # Max SZA for filtering
    maxCHI             = 2.0                    # Max CHI_y_2 value
    maxTC              = 5.0E24                 # Max Total column amount for filtering
    minTC              = 0.0                    # Min Total column amount for filtering
    sclfct             = 1.0E9                  # Scale factor to apply to vmr plots (ppmv=1.0E6, ppbv=1.0E9, etc)
    sclfctName         = 'ppbv'                 # Name of scale factor for labeling plots

    #----------------------
    # Date range to process
    #----------------------
    iyear              = 2010
    imnth              = 1
    iday               = 1
    fyear              = 2014
    fmnth              = 12
    fday               = 31

    #-------------------------------
    # Decide how to create layers
    # above and below the tropopause
    #-------------------------------
    redFct  = 5.0
    hgtTrgU = 10.0       # Upper Layer thickness [km]
    hgtTrgL = False       # Lower Layer thickness [km];Set to False if you want surface value

                                    #----------------------------#
                                    #        --- START ---       #
                                    #----------------------------#

    if loc == 'fl0':
        sLat              = 40.4             #--LATITUDE OF BOULDER
        sLon              = -105.24          #--LONGITUDE OF BOULDER
    else:
        print 'Need a location (currently only fl0)'
        sys.exit()

    #-------------------------------------------------
    #CAM-CHEM
    #-------------------------------------------------
    if ReadCAM:
        DataCAM = dm.CAMClass(dataDirCAM, fileCAM,  outFname= pltCAMFile, saveFlg= saveCAMFlg)
        DataCAM.ReadOutputCAM(gasNameCAM, pCols, sLat, sLon)
        if pltCAM:
            DataCAM.PltCAM()

    #-------------------------------------------------
    #FTS
    #-------------------------------------------------
    if ReadFTS:
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
        for gasVer in statDataCl:
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
        avkSCFav       = OrderedDict()
        avkVMRav      = OrderedDict()

        aPrfVMR      = OrderedDict()
        aPrfMol      = OrderedDict()

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
                statDataCl[gasVer].readError(totFlg=False,sysFlg=False,randFlg=False,vmrFlg=True,avkFlg=True,KbFlg=False)
                    
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
        
            #-------------------------------------------------
            # Calculate layer information for each observation
            #-------------------------------------------------
            LayVMR[gasVer]    = np.zeros(len(dates[gasVer]))
            DryLayVMR[gasVer] = np.zeros(len(dates[gasVer]))
            upperHgt[gasVer]  = np.zeros(len(dates[gasVer]))
            lowerHgt[gasVer]  = np.zeros(len(dates[gasVer]))
            LayThk[gasVer]    = np.zeros(len(dates[gasVer]))
            LayDOF[gasVer]    = np.zeros(len(dates[gasVer]))
            
            for i,dateDay in enumerate(dates[gasVer]):        
                indU             = dc.nearestind(hgtTrgU,alt[gasVer])
                if hgtTrgL: indL = dc.nearestind(hgtTrgL,alt[gasVer])
                else:       indL = -1
                
                upperHgt[gasVer][i]  = alt[gasVer][indL]
                lowerHgt[gasVer][i]  = alt[gasVer][indU]
                
                if indL < 0:
                    LayDOF[gasVer][i]    = np.sum(dofs[i,indU+1:])
                    LayVMR[gasVer][i]    = np.average(rPrfVMR[gasVer][i,indU+1:],weights=Airmass[gasVer][i,indU+1:])
                    DryLayVMR[gasVer][i] = np.average(rPrfDry[gasVer][i,indU+1:],weights=Airmass[gasVer][i,indU+1:])
                else:
                    LayDOF[gasVer][i]    = np.sum(dofs[i,indU+1:indL+1])
                    LayVMR[gasVer][i]    = np.average(rPrfVMR[gasVer][i,indU+1:indL+1],weights=Airmass[gasVer][i,indU+1:indL+1])
                    DryLayVMR[gasVer][i] = np.average(rPrfDry[gasVer][i,indU+1:indL+1],weights=Airmass[gasVer][i,indU+1:indL+1])
             
                LayThk[gasVer][i]    = alt[gasVer][indU] - alt[gasVer][indL]
        
        #----------
        # Plot Data
        #----------
    RatioMean = []
    RatioStd  = []
    GasStrall = []
    RateFTS   = []
    RateCAM   = []
    RateFTSe  = []
    RateCAMe  = []


    if pltAnalysis:

        if saveFlg: pdfsav = PdfPages(pltFile)
        clmap        = 'jet'
        cm           = plt.get_cmap(clmap)
        dayLc        = DayLocator()
        yearsLc      = YearLocator()
        monthLc      = MonthLocator()
        mondays      = WeekdayLocator(MONDAY)
        #DateFmt      = DateFormatter('%b %d')
        DateFmt      = DateFormatter('%Y')

        fign,  axn = plt.subplots(len(gasName)/2, 2, figsize=(12, 12), sharex=True)

        for k, gasVer in enumerate(statDataCl):

            gas2plt = gasName[k].lower()
            gasSTR  = dm.GasName(gas2plt)
            GasStrall.append(gasSTR)

            if   gas2plt == 'co':    scFactor = 1.0
            elif gas2plt == 'c2h6':  scFactor = 1.0
            else: scFactor = 1.0

            maxalt = 16.0

            #---------------------------------
            # Defining variables (FTIR)
            #---------------------------------
            indsalt       = np.where(alt[gasVer] <= maxalt)[0]
            alt_i         = alt[gasVer][indsalt]
            TC_i          = totClmn[gasVer]
            Dates_i       = dates[gasVer]
            Prf_i         = rPrfVMR[gasVer][:, indsalt]
            avkSCFavs_i   = avkSCFav[gasVer][indsalt[0]:, indsalt[0]:]
            avkVMRavs_i   = avkVMRav[gasVer][indsalt[0]:, indsalt[0]:]
            
            aPrf_i        = aPrfVMR[gasVer][:, indsalt]
            Airmass_i     = Airmass[gasVer][:, indsalt]
            rPrfMol_i     = rPrfMol[gasVer][:,indsalt]
            TCpc_i        = np.sum(rPrfMol_i,axis=1)

            indh1 = mf.nearestind(pCols[0], alt_i)
            indh2 = mf.nearestind(pCols[1], alt_i)               
            vmrP = np.average(Prf_i[:,indh2:indh1],axis=1,weights=Airmass_i[:,indh2:indh1]) 

            #---------------------------------
            # Defining variables (CAM-CHEM) - Data is monthly mean already
            #---------------------------------
            altCAM      = DataCAM.CAM['midpoints'][0, :, DataCAM.CAM['indsLoc'][0],DataCAM.CAM['indsLoc'][1]]
            indsalt     = np.where(altCAM <= maxalt)[0]
            altCAM      = altCAM[indsalt]

            if gas2plt == 'h2co': gas2plt = 'ch2o'

            DatesCAM    = np.asarray(DataCAM.CAM['dates'])
            PrfCAMi     = np.asarray(DataCAM.CAM['GasPrf_'+gas2plt][:,:,DataCAM.CAM['indsLoc'][0],DataCAM.CAM['indsLoc'][1]])
            TCCAMi      = np.asarray(DataCAM.CAM['GasTC_'+gas2plt][:,DataCAM.CAM['indsLoc'][0],DataCAM.CAM['indsLoc'][1]])
            AirmassCAMi = np.asarray(DataCAM.CAM['AIRMASS_'+gas2plt][:,:,DataCAM.CAM['indsLoc'][0],DataCAM.CAM['indsLoc'][1]])

            PrfCAM      = PrfCAMi[:, indsalt]
            AirmassCAM  = AirmassCAMi[:, indsalt]
            TCCAM       = np.sum( (PrfCAM) *AirmassCAM, axis=1)
 
            #---------------------------------
            # Smoothing CAM-CHEM using FTIR AK and apriori
            #---------------------------------
            PrfCAM_interp     = interpolate.interp1d(altCAM, PrfCAM, axis=1, fill_value='extrapolate', bounds_error=False)(alt_i)*1e9
            AirmassCAM_interp = interpolate.interp1d(altCAM, AirmassCAM, axis=1, fill_value='extrapolate', bounds_error=False)(alt_i)

            aPrfm   = np.mean(aPrf_i, axis=0)

            PrfCAMs = np.zeros( (len(DatesCAM), len(alt_i)) )
           
            for itime in range(len(DatesCAM)):
                #The equations below are equivalent
                PrfCAMs[itime, :]  = aPrfm + np.dot(avkVMRavs_i, (PrfCAM_interp[itime, :] -  aPrfm))         # x_s = Apriori(FTS) + AK(Prf(CAM-CHEM) - Apriori(FTS))
                #PrfCAMs[itime, :]  = np.dot(avkSCFavs, PrfCAM_interp[itime, :]) + np.dot(  np.identity(avkSCFavs.shape[0]) - avkSCFavs,    aPrfm )   # x_s = AK*Prf(CAM-CHEM) + (I-AK)*Apriori

            #BELOW I USE THE MEAN AIRMASS FROM THE FTIR, THE AIRMASS FROM THE MODEL IS NOT ADEQUATE
            AirmassFTSMean = np.mean(Airmass_i, axis=0)
            AirmassFTSstd = np.std(Airmass_i, axis=0)

            TCCAMs  = np.sum( (PrfCAMs/1e9) *AirmassFTSMean, axis=1)
            vmrPCAM = np.average(PrfCAMs[:,indh2:indh1], axis=1, weights=AirmassFTSMean[indh2:indh1]) 
 
            #---------------------------------
            # Plot : Averaging Kernel Smoothing Function (row of avk)
            #---------------------------------
            fig       = plt.figure(figsize=(9,9))
            gs        = gridspec.GridSpec(1,2,width_ratios=[3,1])
            ax        = plt.subplot(gs[0])
            axb       = plt.subplot(gs[1])
            cm        = plt.get_cmap(clmap)
            cNorm     = colors.Normalize(vmin=np.min(alt_i), vmax=np.max(alt_i))
            scalarMap = mplcm.ScalarMappable(norm=cNorm,cmap=clmap)
            scalarMap.set_array(alt_i)
            ax.set_color_cycle([scalarMap.to_rgba(x) for x in alt_i])
            
            for i in range(len(alt_i)):
                ax.plot(avkSCFavs_i[i,:],alt_i)
                
            ax.set_ylabel('Altitude [km]', fontsize=14)
            ax.set_xlabel('Averaging Kernels', fontsize=14)
            ax.grid(True)
            cbar = fig.colorbar(scalarMap, orientation='vertical')
            cbar.set_label('Altitude [km]', fontsize=14)
            ax.set_title(gasSTR + ' Averaging Kernels Scale Factor', fontsize=14)
            ax.tick_params(labelsize=14)
            
            axb.plot(np.sum(avkSCFavs_i,axis=0), alt_i,color='k')
            axb.grid(True)
            axb.set_xlabel('Averaging Kernel Area', fontsize=14)
            axb.tick_params(axis='x',which='both',labelsize=8)
            #axb.tick_params(labelsize=14) 

            if saveFlg:     pdfsav.savefig(fig,dpi=200)
            else:           plt.show(block=False)

            #---------------------------------
            # Plot : Prf Mean
            #---------------------------------

            PrfmeanCAM = np.mean(PrfCAMs, axis=0)
            prfSTDCAM  = np.std(PrfCAMs, axis=0)

            prfMean = np.mean(Prf_i, axis=0)
            prfSTD  = np.std(Prf_i, axis=0)

            PrfmeanCAMi = np.mean(PrfCAM, axis=0)

            fig,ax  = plt.subplots(figsize=(7,9))

            ax.plot(prfMean,alt_i, linewidth = 2.0, color='k', label='FTIR')
            ax.scatter(prfMean,alt_i, facecolors='white', s=60, color='k')
            ax.fill_betweenx(alt_i,prfMean-prfSTD,prfMean+prfSTD, alpha=0.5, color='k') 

            ax.plot(PrfmeanCAMi*1e9*scFactor,altCAM, linewidth = 2.0, color='r', label='CAM-Chem x '+' '+str(scFactor))
            ax.scatter(PrfmeanCAMi*1e9*scFactor,altCAM, facecolors='white', s=30, color='r') 

            ax.plot(PrfmeanCAM*scFactor,alt_i, linewidth = 2.0, color='green', label='CAM-Chem smoothed x '+' '+str(scFactor))
            ax.scatter(PrfmeanCAM*scFactor,alt_i, facecolors='white', s=60, color='green')
            ax.fill_betweenx(alt_i,PrfmeanCAM*scFactor-prfSTDCAM,PrfmeanCAM*scFactor+prfSTDCAM,alpha=0.5,color='green') 

    
            ax.set_title('Mean Profile of '+gasSTR, fontsize=14)
            ax.set_ylabel('Altitude [km]', fontsize=14)
            ax.set_xlabel('VMR [ppb$_v$]', fontsize=14)    
            ax.grid(True,which='both')
            ax.tick_params(labelsize=14)
            ax.legend(prop={'size':12})

            if saveFlg:     pdfsav.savefig(fig,dpi=200)
            else:           plt.show(block=False)   
            
            #---------------------------------
            # Plot : Monthly averages of total columns
            #---------------------------------

            mnthlyVals   = mf.mnthlyAvg(TCpc_i, Dates_i,dateAxis=1, meanAxis=0)
            dateYearFrac = mf.toYearFraction(mnthlyVals['dates'])
            weights      = np.ones_like(dateYearFrac)
            res          = mf.fit_driftfourier(dateYearFrac, mnthlyVals['mnthlyAvg'], weights, 2)
            f_drift, f_fourier, f_driftfourier = res[3:6]

            fig, ax = plt.subplots(1, figsize=(10,6), sharex=True)
            ax.plot(mnthlyVals['dates'],mnthlyVals['mnthlyAvg'], color='k', label='FTIR')
            ax.scatter(mnthlyVals['dates'],mnthlyVals['mnthlyAvg'], facecolors='white', s=60, color='k')
        
            ax.plot(DatesCAM, TCCAM*scFactor, color='r', label='CAM-Chem x '+' '+str(scFactor))
            ax.scatter(DatesCAM, TCCAM*scFactor, facecolors='white', s=60, color='r')

            ax.plot(DatesCAM, TCCAMs*scFactor, color='green', label='CAM-Chem smoothed x '+' '+str(scFactor))
            ax.scatter(DatesCAM, TCCAMs*scFactor, facecolors='white', s=60, color='green')
            
            ax.grid(True)
            ax.set_ylabel('Partial Column\n[molecules$\cdot$cm$^{-2}$]',fontsize=14)
            ax.set_title(gasSTR + ' Partial Column [Monthly average], '+ str(alt_i[-1])+'[km] - '+str(alt_i[0])+'[km]',multialignment='center',fontsize=14)
            ax.tick_params(labelsize=14)
            ax.set_xlabel('Year', fontsize=14)
            ax.xaxis.set_minor_locator(monthLc)
            ax.xaxis.set_major_formatter(DateFmt)
            ax.legend(prop={'size':12})


            if saveFlg:     pdfsav.savefig(fig,dpi=200)
            else:           plt.show(block=False)   

            #---------------------------------
            # Plot : Monthly averages of total columns (Normalized)
            #---------------------------------

            fig, ax = plt.subplots(1, figsize=(10,6), sharex=True)       
             
            ax.plot(mnthlyVals['dates'],mnthlyVals['mnthlyAvg']/np.mean(mnthlyVals['mnthlyAvg']), color='k', label='FTIR')
            ax.scatter(mnthlyVals['dates'],mnthlyVals['mnthlyAvg']/np.mean(mnthlyVals['mnthlyAvg']), facecolors='white', s=60, color='k')

            ax.plot(DatesCAM, TCCAM/np.mean( TCCAM), color='r', label='CAM-Chem x '+' '+str(scFactor))
            ax.scatter(DatesCAM, TCCAM/np.mean( TCCAM), facecolors='white', s=60, color='r')

            ax.plot(DatesCAM, TCCAMs/np.mean( TCCAMs), color='green', label='CAM-Chem smoothed x '+' '+str(scFactor))
            ax.scatter(DatesCAM, TCCAMs/np.mean( TCCAMs), facecolors='white', s=60, color='green')


            ax.grid(True)
            ax.set_ylabel('molecules$\cdot$cm$^{-2}$ / Mean molecules$\cdot$cm$^{-2}$', fontsize=14)
            ax.set_title(gasSTR +' Partial Column Normalized, '+ str(alt_i[-1])+'[km] - '+str(alt_i[0])+'[km]', fontsize=14)
            #box = ax.get_position()
            #ax.set_position([box.x0, box.y0, box.width * 0.90, box.height])
            ax.legend(loc='center left',bbox_to_anchor=(1,0.5),prop={'size':8})
            ax.tick_params(labelsize=14)
            ax.set_xlabel('Year', fontsize=14)
            ax.xaxis.set_minor_locator(monthLc)
            ax.xaxis.set_major_formatter(DateFmt)
            ax.legend(prop={'size':12})


            if saveFlg:     pdfsav.savefig(fig,dpi=200)
            else:           plt.show(block=False)


            #---------------------------------
            # Plot : VMR and partial columns weighted and shown for co-located years/months
            #---------------------------------
            #FInding the same datetime of measurements in diff versions
            doy    = mf.toYearFraction(mnthlyVals['dates'])
            doyCAM = mf.toYearFraction(DatesCAM)

            intrsctVals = np.intersect1d(doy, doyCAM, assume_unique=False)
         
            inds1       = np.nonzero( np.in1d( doyCAM, intrsctVals, assume_unique=False ) )[0]
            inds2       = np.nonzero( np.in1d( doy, intrsctVals, assume_unique=False ) )[0]

            vmrPinters        = vmrP[inds2]
            vmrPCAMinters     = vmrPCAM[inds1]

            TCintersec        = mnthlyVals['mnthlyAvg'][inds2]
            TCCAMintersec     = TCCAMs[inds1]
            TCSDintersec      = mnthlyVals['std'][inds2]

            Datesintersec     = mnthlyVals['dates'][inds2]
            DatesCAMintersec  = DatesCAM[inds1]

            Ratio             = np.divide(TCintersec,TCCAMintersec)

            #---------------------------------
            # Trend Analysis
            #---------------------------------

            #----------------------------------------------------------------------
            #FTIR
            #----------------------------------------------------------------------
            #mnthlyVals       = mf.mnthlyAvg(vmrP[gv], DT[gv],dateAxis=1, meanAxis=0)
            dateYearFrac     = mf.toYearFraction(Datesintersec)
            weights          = np.ones_like(dateYearFrac)
            res              = mf.fit_driftfourier(dateYearFrac,TCintersec, weights, 2)
            intercept, slope, pfourier = res[0:3]
            f_drift, f_fourier, f_driftfourier = res[3:6]

            int_fourier_m    = intercept
            sl_fourier_m     = slope
            p_fourier_m      = pfourier
            FAT_m            = f_drift(dateYearFrac)
            FATAV_m          = f_driftfourier(dateYearFrac)

            #-----------------------------------
            # bootstrap resampling information of Monthly averages
            #-----------------------------------
            perc, intercept_boot, slope_boot, pfourier_boot = mf.cf_driftfourier(dateYearFrac, TCintersec, weights, 2)

            int_boot_m       = intercept_boot
            sl_boot_m        = slope_boot
            p_boot_m         = pfourier_boot
            RateFTS.append( np.divide(res[1], np.mean(TCintersec))*100.0)
            RateFTSe.append( np.divide(np.std(slope_boot), np.mean(TCintersec))*100.0)

            #----------------------------------------------------------------------
            #CAM-CHEM
            #----------------------------------------------------------------------
            res              = mf.fit_driftfourier(dateYearFrac,TCCAMintersec, weights, 2)
            intercept, slope, pfourier = res[0:3]
            f_drift, f_fourier, f_driftfourier = res[3:6]

            int_fourier_m    = intercept
            sl_fourier_m     = slope
            p_fourier_m      = pfourier
            FAT_cam_m            = f_drift(dateYearFrac)
            FATAV_cam_m          = f_driftfourier(dateYearFrac)


            #-----------------------------------
            # bootstrap resampling information of Monthly averages
            #-----------------------------------
            perc, intercept_boot, slope_boot, pfourier_boot = mf.cf_driftfourier(dateYearFrac, TCCAMintersec, weights, 2)

            int_boot_m       = intercept_boot
            sl_boot_m        = slope_boot
            p_boot_m         = pfourier_boot
            RateCAM.append( np.divide(res[1], np.mean(TCCAMintersec))*100.0)
            RateCAMe.append( np.divide(np.std(slope_boot), np.mean(TCCAMintersec))*100.0)

            #---------------------------------
            fig, ax = plt.subplots(1, figsize=(10,6), sharex=True)        
             
            ax.plot(Datesintersec, TCintersec, color='k', label='FTIR')
            #ax.errorbar(Datesintersec, TCintersec,yerr=TCSDintersec,markersize=60, color='k', ecolor='k', label='FTIR')
            ax.fill_between(Datesintersec ,TCintersec-TCSDintersec,TCintersec+TCSDintersec,alpha=0.5,color='0.75')

            ax.scatter(Datesintersec,TCintersec, facecolors='white', s=60, color='k')
            #-------fitted trend 
            ax.plot(Datesintersec,FAT_m,label='Fitted Anual Trend',linewidth=2.5)
            ax.plot(Datesintersec,FATAV_m,label='Fitted Anual Trend + intra-annual variability', linewidth=2.5)

            ax.plot(DatesCAMintersec, TCCAMintersec, color='green', label='CAM-Chem smoothed x '+' '+str(scFactor))
            ax.scatter(DatesCAMintersec, TCCAMintersec, facecolors='white', s=60, color='green')

            ax.plot(Datesintersec,FAT_cam_m,label='Fitted Anual Trend',linewidth=2.5)
            ax.plot(Datesintersec,FATAV_cam_m,label='Fitted Anual Trend + intra-annual variability', linewidth=2.5)

            ax.grid(True)
            ax.set_ylabel('Partial Column\n[molecules$\cdot$cm$^{-2}$]', fontsize=14)
            ax.set_title(gasSTR + ' Partial Column [Monthly average], '+ str(alt_i[-1])+'[km] - '+str(alt_i[0])+'[km]', fontsize=14)
            ax.legend(loc='center left',bbox_to_anchor=(1,0.5),prop={'size':8})
            ax.tick_params(labelsize=14)
            ax.set_xlabel('Year', fontsize=14)
            ax.xaxis.set_minor_locator(monthLc)
            ax.xaxis.set_major_formatter(DateFmt)
            ax.legend(prop={'size':12})

            if saveFlg:     pdfsav.savefig(fig,dpi=200)
            else:           plt.show(block=False)

            #---------------------------------
            # Plot : All Months as function of month
            #---------------------------------

            years     = np.asarray([dt.date(d.year,1,1) for d in mnthlyVals['dates'][inds2]])


            #---------------------------------
            # Plot : Ratio of FTS/CAM-CHEM
            #---------------------------------
            fig, ax = plt.subplots(1, figsize=(10,6), sharex=True)        
             
            ax.plot(Datesintersec, Ratio , color='k')
            ax.scatter(Datesintersec,Ratio, facecolors='white', s=60, color='k')

            ax.grid(True)
            ax.set_ylabel('Partial Column\n[molecules$\cdot$cm$^{-2}$] - FTS/CAM-CHEM', fontsize=14)
            ax.set_title(gasSTR + ' Partial Column [Monthly average], '+ str(alt_i[-1])+'[km] - '+str(alt_i[0])+'[km]\nFTS/CAM-CHEM', fontsize=14)
            ax.legend(loc='center left',bbox_to_anchor=(1,0.5),prop={'size':8})
            ax.tick_params(labelsize=14)
            ax.set_xlabel('Year', fontsize=14)
            ax.xaxis.set_minor_locator(monthLc)
            ax.xaxis.set_major_formatter(DateFmt)
            ax.legend(prop={'size':12})
            ax.text(0.02,0.94, 'Mean Ratio ($\pm$ std) = {0:1.2f} ({1:1.2f})'.format(np.mean(Ratio), np.std(Ratio)), transform=ax.transAxes, fontsize=14)

            RatioMean.append(np.mean(Ratio))
            RatioStd.append(np.std(Ratio))

            if saveFlg:     pdfsav.savefig(fig,dpi=200)
            else:           plt.show(block=False)


            # fig, ax = plt.subplots(1, figsize=(10,6), sharex=True) 

            # ax.plot(Datesintersec, vmrPinters, color='k', label='FTIR')
            # ax.scatter(Datesintersec,vmrPinters, facecolors='white', s=60, color='k')

            # ax.plot(DatesCAMintersec, vmrPCAMinters, color='green', label='CAM-Chem smoothed x '+' '+str(scFactor))
            # ax.scatter(DatesCAMintersec, vmrPCAMinters, facecolors='white', s=60, color='green')

            # ax.grid(True)
            # ax.set_ylabel('VMR [ppb$_v$]', fontsize=14)
            # ax.set_title(gasSTR +' VMR weighted, '+ str(alt_i[indh1])+'[km] - '+str(alt_i[indh2])+'[km]', fontsize=14)
            # ax.set_xlabel('Year')
            # ax.legend(loc='center left',bbox_to_anchor=(1,0.5),prop={'size':8})
            # ax.tick_params(labelsize=14)
            # ax.set_xlabel('Year', fontsize=14)
            # ax.xaxis.set_minor_locator(monthLc)
            # ax.xaxis.set_major_formatter(DateFmt)
            # ax.legend(prop={'size':12})

            # if saveFlg:     pdfsav.savefig(fig,dpi=200)
            # else:           plt.show(block=False)

            #----------------------------
            # Plot total columns by month
            #----------------------------
            #FTS
            month    = np.array([d.month for d in Dates_i])
            mnthSort = list(set(month))
            mnthMean = np.zeros(len(mnthSort))
            mnthSTD  = np.zeros(len(mnthSort))
            
            for i,m in enumerate(mnthSort):
                inds        = np.where(month == m)[0]
                mnthMean[i] = np.mean(TCpc_i[inds])
                mnthSTD[i]  = np.std(TCpc_i[inds])

            #CAM-CHEM
            monthCAM    = np.array([d.month for d in DatesCAM])
            mnthSortCAM = list(set(monthCAM))
            mnthMeanCAM = np.zeros(len(mnthSortCAM))
            mnthSTDCAM  = np.zeros(len(mnthSortCAM))
            
            for i,m in enumerate(mnthSortCAM):
                inds        = np.where(monthCAM == m)[0]
                mnthMeanCAM[i] = np.mean(TCCAMs[inds])
                mnthSTDCAM[i]  = np.std(TCCAMs[inds])

            fig,ax1  = plt.subplots()

            ax1.plot(mnthSort,mnthMean, color='k', label='FTIR')
            ax1.scatter(mnthSort,mnthMean, facecolors='white', s=60, color='k')

            ax1.plot(mnthSortCAM,mnthMeanCAM, color='green', label='CAM-Chem')
            ax1.scatter(mnthSortCAM,mnthMeanCAM, facecolors='white', s=60, color='green')    
            ax1.grid(True,which='both')
            ax1.set_ylabel('molecules$\cdot$cm$^{-2}$',multialignment='center')
            ax1.set_xlabel('Month')
            ax1.set_title('Retrieved Partial Column Monthly Mean')
            ax1.set_xlim((0,13))
            ax1.set_xticks(range(1,13))
            ax1.legend(prop={'size':12})

            if saveFlg:     pdfsav.savefig(fig,dpi=200)
            else:           plt.show(block=False)
                
            fig,ax1  = plt.subplots()

            ax1.plot(mnthSort,mnthMean/np.mean(mnthMean), color='k', label='FTIR')
            ax1.scatter(mnthSort,mnthMean/np.mean(mnthMean), facecolors='white', s=60, color='k')
            #ax1.errorbar(mnthSort,mnthMean/np.mean(mnthMean),yerr=mnthSTD,fmt='k.',markersize=6,ecolor='K')

            ax1.plot(mnthSortCAM,mnthMeanCAM/np.mean(mnthMeanCAM), color='green', label='CAM-Chem')
            ax1.scatter(mnthSortCAM,mnthMeanCAM/np.mean(mnthMeanCAM), facecolors='white', s=60, color='green')
            #ax1.errorbar(mnthSortCAM,mnthMeanCAM/np.mean(mnthMeanCAM),yerr=mnthSTDCAM,fmt='r.',markersize=6,ecolor='red')     
            ax1.grid(True,which='both')
            ax1.set_ylabel('molecules$\cdot$cm$^{-2}$ / Mean molecules$\cdot$cm$^{-2}$',multialignment='center')
            ax1.set_xlabel('Month')
            ax1.set_title('Retrieved Partial Column Monthly Mean')
            ax1.set_xlim((0,13))
            ax1.set_xticks(range(1,13))
            ax1.legend(prop={'size':12})

            if saveFlg:     pdfsav.savefig(fig,dpi=200)
            else:           plt.show(block=False)

            #--------------

            #----------------------------------------
            # Create a list of months to loop through
            #----------------------------------------
            uniqueyears     = list(set(years))          # Find a list of unique months
            uniqueyears.sort()
            uniqueyears     = np.asarray(uniqueyears)

            fig, ax = plt.subplots(1, figsize=(10,6), sharex=True)

            clr = mf.clrplt()

            for ii, yyyy in enumerate(uniqueyears):

                indsY = where(mf.toYearFraction(years) == mf.toYearFraction(uniqueyears)[ii] )[0]
                month = np.array([d.month for d in Datesintersec[indsY]])

                #ax.plot(month, TCintersec[indsY], color='k', label='FTIR')
                #ax.scatter(month,TCintersec[indsY], facecolors='white', s=60, color=clr[ii])
   
                ax.errorbar(month, TCintersec[indsY], yerr=TCSDintersec[indsY],fmt='o', markersize=7.5, color=clr[ii], ecolor=clr[ii], label=yyyy.year)

                
                if k <= ( len(gasName)/2 - 1):
                    axn[k, 0].errorbar(month, TCintersec[indsY], yerr=TCSDintersec[indsY],fmt='o', markersize=7.5, color=clr[ii], ecolor=clr[ii], label=yyyy.year)
                    axn[k, 0].plot(month, TCCAMintersec[indsY], color=clr[ii], alpha=0.5, linewidth = 2.0)
                    axn[k, 0].grid(True)
                    if k == 1: axn[k, 0].set_ylabel('Partial Column [molecules$\cdot$cm$^{-2}$]', fontsize=18)
                    #axn[k, 0].set_title(gasSTR, fontsize=14)
                    axn[k, 0].tick_params(labelsize=18)
                    if k == ( len(gasName)/2 - 1): axn[k, 0].set_xlabel('Month', fontsize=18)
                    axn[k, 0].set_xlim((0,13))
                    if k == 0: axn[k, 0].legend(prop={'size':12})
                else:

                    axn[k-(len(gasName)/2), 1].errorbar(month, TCintersec[indsY], yerr=TCSDintersec[indsY],fmt='o', markersize=7.5, color=clr[ii], ecolor=clr[ii], label=yyyy.year)
                    axn[k-(len(gasName)/2), 1].plot(month, TCCAMintersec[indsY], color=clr[ii], alpha=0.5, linewidth = 2.0)
                    axn[k-(len(gasName)/2), 1].grid(True)
                    #axn[k-(len(gasName)/2), 1].set_title(gasSTR, fontsize=14)
                    axn[k-(len(gasName)/2), 1].tick_params(labelsize=18)
                    if k == ( len(gasName)-1): axn[k-(len(gasName)/2), 1].set_xlabel('Month', fontsize=18)
                    axn[k-(len(gasName)/2), 1].set_xlim((0,13))
                


                #ax.fill_between(Datesintersec ,TCintersec-TCSDintersec,TCintersec+TCSDintersec,alpha=0.5,color='0.75')

                #ax.scatter(month,TCintersec[indsY], facecolors='white', s=60, color='k')

                ax.plot(month, TCCAMintersec[indsY], color=clr[ii], alpha=0.5, linewidth = 2.0)
                #ax.scatter(month, TCCAMintersec[indsY], facecolors='white', s=60, color='green')

            #ax.plot(mnthSortCAM,mnthMeanCAM, color='green', label='CAM-Chem')
            #ax.fill_between(mnthSortCAM ,mnthMeanCAM-mnthSTDCAM,mnthMeanCAM+mnthSTDCAM,alpha=0.5,color='green')

            fign.subplots_adjust(bottom=0.05, top=0.95, left = 0.1, right = 0.97)
            

            ax.grid(True)
            ax.set_ylabel('Partial Column\n[molecules$\cdot$cm$^{-2}$]', fontsize=14)
            ax.set_title(gasSTR + ' Partial Column [Monthly average], '+ str(alt_i[-1])+'[km] - '+str(alt_i[0])+'[km]', fontsize=14)
            ax.legend(loc='center left',bbox_to_anchor=(1,0.5),prop={'size':8})
            ax.tick_params(labelsize=14)
            ax.set_xlabel('Month', fontsize=14)
            ax.set_xlim((0,13))
            #ax.xaxis.set_minor_locator(monthLc)
            #ax.xaxis.set_major_formatter(DateFormatter('%m'))
            ax.legend(prop={'size':12})

            if saveFlg:     pdfsav.savefig(fig,dpi=200)
            else:           plt.show(block=False)

        if saveFlg:     pdfsav.savefig(fign,dpi=200)
        else:           plt.show(block=False)
            
            #---------------------------------
            # Plot :Airmass
            #---------------------------------

            # Prfmairmass_1  = np.mean(AirmassCAM, axis=0)
            # Prfmairmass_2  = np.mean(AirmassCAM_interp, axis=0)        

            # fig,ax  = plt.subplots(figsize=(7,9))
            # ax.plot(Prfmairmass_1 ,altCAM, linewidth = 2.0, color='k')
            # ax.scatter(Prfmairmass_1,altCAM, facecolors='white', s=60, color='k')

            # ax.plot(Prfmairmass_2,alt_i, linewidth = 2.0, color='green')
            # ax.scatter(Prfmairmass_2,alt_i, facecolors='white', s=60, color='green')

            # ax.plot(AirmassFTSMean,alt_i, linewidth = 2.0, color='r')
            # ax.scatter(AirmassFTSMean,alt_i, facecolors='white', s=60, color='r')
            # ax.fill_betweenx(alt_i,AirmassFTSMean-AirmassFTSstd*2.,AirmassFTSMean+AirmassFTSstd*2.,alpha=0.5,color='green') 
        
            # ax.set_title('Mean Profile of '+gasSTR, fontsize=14)
            # ax.set_ylabel('Altitude [km]', fontsize=14)
            # ax.set_xlabel('VMR [ppb$_v$]', fontsize=14)    
            # ax.grid(True,which='both')
            # ax.tick_params(labelsize=14)

            # if saveFlg:     pdfsav.savefig(fig,dpi=200)
            # else:           
            #     plt.show(block=False)
            #     user_input = raw_input('Press any key to exit >>> ')
            #     sys.exit()


            #-----------------------------------------------------------------
            #Plot: Monthly by year (for PAPER)
            #-----------------------------------------------------------------


        #-----------------------------------------------------------------
        #Plot: Ratio of FTS to CAM-Chem in bar plot
        #-----------------------------------------------------------------

        #RatioMean = np.asarray(RatioMean)
        #RateFTS   = np.asarray(RateFTS)
        #RateFTSe  = np.asarray(RateFTSe)
        #GasStrall = np.asarray(GasStrall)

        RatioMean.extend([0.0])
        RatioStd.extend([0.0])
        RateFTS.extend([2.4])
        RateFTSe.extend([3.0])
        GasStrall.extend(['NH$_3$'])

        RateCAM.extend([0.0])
        RateCAMe.extend([0.0])


        N = len(RatioMean)
        ind = np.arange(N)  # the x locations for the groups

        fig, ax = plt.subplots()
        Ra1 = ax.bar(ind, RatioMean, 0.3, color='r', yerr=RatioStd)
        # add some text for labels, title and axes ticks
        ax.set_ylabel('Ratio HR-FTIR/CAM-Chem', fontsize=18)
        #ax.xticks(ind + 0.35/2., gasName)
        ax.set_xlabel('Gas', fontsize=14)
        ax.set_xticks(ind+0.35/2.0)
        ax.set_xticklabels(GasStrall)
        ax.tick_params(labelsize=14)
        ax.set_title('Ratio FTIR/CAM-Chem', fontsize=14)

        if saveFlg:     pdfsav.savefig(fig,dpi=200)
        else:           plt.show(block=False)



        fig, ax = plt.subplots()
        Ra1 = ax.bar(ind, RateFTS, width = 0.3, align='center', color = 'r', yerr=RateFTSe, ecolor = 'k', label = 'HR-FTIR')
        ax.bar(ind+0.3,RateCAM, width = 0.3, align='center', color = 'b', yerr=RateCAMe, ecolor = 'k', label = 'CAM-Chem')
        # add some text for labels, title and axes ticks
        ax.set_ylabel('Annual rate of change (%)', fontsize=18)
        #ax.xticks(ind + 0.35/2., gasName)
        ax.set_xlabel('Gas', fontsize=18)
        ax.set_xticks(ind+0.3/2.0)
        ax.set_xticklabels(GasStrall,  rotation=45)
        ax.tick_params(labelsize=18)
        ax.legend(prop={'size':12}, loc=4)
        ax.axhline(0, color='black', lw=1)
        ax.yaxis.grid(True)
        ax.xaxis.set_tick_params(which='major',labelsize=18)
        ax.yaxis.set_tick_params(which='major',labelsize=18)
        ax.set_xlim(0.5, N-0.5)
        ax.set_ylim(-15, 7)
        #ax.set_title('Annual rate of change (%)', fontsize=14)
        fig.subplots_adjust(left=0.125, bottom=0.15, right=0.95, top=0.95)

        if saveFlg:     pdfsav.savefig(fig,dpi=200)
        else:           plt.show(block=False)


        fig , (ax, ax1) = plt.subplots(2, figsize=(8,8), gridspec_kw = {'height_ratios':[3, 1]}, sharex=True)
        
        ax.bar(ind, RateFTS, width = 0.27, align='center', color = 'r', yerr=RateFTSe, ecolor = 'k', label = 'HR-FTIR')
        ax.bar(ind+0.27,RateCAM, width = 0.27, align='center', color = 'b', yerr=RateCAMe, ecolor = 'k', label = 'CAM-Chem')
        ax.yaxis.grid(True)
        ax.set_xticks(ind)
        ax.set_ylabel('Annual rate of change (%)', fontsize = 18)
        ax.set_xticklabels(GasStrall)
        ax.set_xticks(ind+0.27/2.0)
        ax.axhline(0, color='black', lw=1)
        ax.legend(prop={'size':12}, loc=4)
        ax.xaxis.set_tick_params(which='major',labelsize=18)
        ax.yaxis.set_tick_params(which='major',labelsize=18)

        ax1.bar(ind, RatioMean, 0.27, color='gray', yerr=RatioStd, ecolor = 'k')
        # add some text for labels, title and axes ticks
        ax1.set_ylabel('HR-FTIR/CAM-Chem', fontsize=14)
        ax1.set_xlabel('Gas', fontsize=18)
        ax1.set_xticks(ind+0.27/2.0)
        ax1.set_xticklabels(GasStrall)
        ax1.tick_params(labelsize=18)

        #fig.tight_layout()

        if saveFlg:     pdfsav.savefig(fig,dpi=200)
        else:           plt.show(block=False)

        
    

        #-------------------------------------------------
        #IF SAVING FILES
        #-------------------------------------------------
        #-------------------------
        # Close output figure file
        #-------------------------
 
        if saveFlg: pdfsav.close()
        else:
            user_input = raw_input('Press any key to exit >>> ')
            sys.exit()

    
    
    #if saveFlg: pdfsav.close()
    #else:
    #    user_input = raw_input('Press any key to exit >>> ')
    #    sys.exit()

 
if __name__ == "__main__":
    main()