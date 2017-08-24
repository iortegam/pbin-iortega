#! /usr/bin/python2.7
#----------------------------------------------------------------------------------------
# Name:
#        pltts.py
#
# Purpose:
#
# Plot time series of columns/VMRs for several trace gases defined in pltts_input.py
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
#from sklearn import linear_model, datasets

import numpy as np



def usage():
    ''' Prints to screen standard program usage'''
    print 'pltSet.py -i <inputfile> -?'
    print '  -i <file> : Run pltSet.py with specified input file'
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

    ngases = len(pltInputs['gasName'])

 
    #--------------------------------
    # Define the Retrieval arguments
    #--------------------------------
    TC        =  []         #TOTAL COLUMN
    DT        =  []         #DATE and TIME
    DOF       =  []         #DEGREES OF FREEDOM
    RMS       =  []         #RMS
    TCRE      =  []         #TOTAL COLUMN RANDOM ERROR
    TCSE      =  []         #TOTAL COLUMN SYSTEMATIC ERROR
    TCTE      =  []         #TOTAL COLUMN TOTAL ERROR
    gasname_w =  []         #NAME OF THE GAS

    #--------------------------------
    # Define Additional Parameters 
    #--------------------------------
    int_fourier   =  []    #intercept Fourier analysis
    sl_fourier    =  []    #slope Fourier analysis
    p_fourier     =  []
    FAT           =  []    #Fitted Annual Trend
    FATAV         =  []    #Fitted Annual Trend and annual variability
    int_boot      =  []    #intercept Fourier + boot analysis
    sl_boot       =  []    #slope Fourier + boot analysis
    p_boot        =  []
    Rate          =  []    #Annual Rate of change
    Rate_e        =  []    #Uncertainty of the Annual Rate of change

    #--------------------------------
    # Define the Retrieval arguments (DAILY VALUES)
    #--------------------------------
    TC_d        =  []         #TOTAL COLUMN
    DT_d        =  []         #DATE and TIME
    TCsd_d      =  []         #TOTAL COLUMN DAILY STANDARD DEVIATION
    TCRE_d      =  []         #TOTAL COLUMN RANDOM ERROR
    TCSE_d      =  []         #TOTAL COLUMN SYSTEMATIC ERROR
    TCTE_d      =  []         #TOTAL COLUMN TOTAL ERROR

    #--------------------------------
    # Define Additional Parameters  (DAILY VALUES)
    #--------------------------------
    int_fourier_d   =  []    #intercept Fourier analysis
    sl_fourier_d    =  []    #slope Fourier analysis
    p_fourier_d     =  []
    FAT_d           =  []    #Fitted Annual Trend
    FATAV_d         =  []    #Fitted Annual Trend and annual variability
    int_boot_d      =  []    #intercept Fourier + boot analysis
    sl_boot_d       =  []    #slope Fourier + boot analysis
    p_boot_d        =  []
    Rate_d          =  []    #Annual Rate of change
    Rate_e_d        =  []    #Uncertainty of the Annual Rate of change
   
    #--------------------------------
    # Define the Retrieval arguments (MONTHLY VALUES)
    #--------------------------------
    TC_m        =  []         #TOTAL COLUMN
    DT_m        =  []         #DATE and TIME
    TCsd_m      =  []         #TOTAL COLUMN DAILY STANDARD DEVIATION
    TCRE_m      =  []         #TOTAL COLUMN RANDOM ERROR
    TCSE_m      =  []         #TOTAL COLUMN SYSTEMATIC ERROR
    TCTE_m      =  []         #TOTAL COLUMN TOTAL ERROR

    #--------------------------------
    # Define Additional Parameters  (MONTHLY VALUES)
    #--------------------------------
    int_fourier_m   =  []    #intercept Fourier analysis
    sl_fourier_m    =  []    #slope Fourier analysis
    p_fourier_m     =  []
    FAT_m           =  []    #Fitted Annual Trend
    FATAV_m         =  []    #Fitted Annual Trend and annual variability
    int_boot_m      =  []    #intercept Fourier + boot analysis
    sl_boot_m       =  []    #slope Fourier + boot analysis
    p_boot_m        =  []
    Rate_m          =  []    #Annual Rate of change
    Rate_e_m        =  []    #Uncertainty of the Annual Rate of change

    if pltInputs['pCols']:
        vmrP        =  []           
        vmrPDry     =  []
        sumP        =  []
        #--------------------------------
        # Define the Retrieval arguments (DAILY VALUES)
        #--------------------------------
        vmr_d        =  []         #TOTAL COLUMN
        vmrsd_d      =  []         #TOTAL COLUMN DAILY STANDARD DEVIATION

        #--------------------------------
        # Define Additional Parameters  (DAILY VALUES)
        #--------------------------------
        int_fourier_vmr_d   =  []    #intercept Fourier analysis
        sl_fourier_vmr_d    =  []    #slope Fourier analysis
        p_fourier_vmr_d     =  []
        FAT_vmr_d           =  []    #Fitted Annual Trend
        FATAV_vmr_d         =  []    #Fitted Annual Trend and annual variability
        int_boot_vmr_d      =  []    #intercept Fourier + boot analysis
        sl_boot_vmr_d       =  []    #slope Fourier + boot analysis
        p_boot_vmr_d        =  []
        Rate_vmr_d          =  []    #Annual Rate of change
        Rate_e_vmr_d        =  []    #Uncertainty of the Annual Rate of change

        #--------------------------------
        # Define the Retrieval arguments (MONTHLY VALUES)
        #--------------------------------
        vmr_m        =  []         #TOTAL COLUMN
        vmrsd_m      =  []         #TOTAL COLUMN DAILY STANDARD DEVIATION

        #--------------------------------
        # Define Additional Parameters  (MONTHLY VALUES)
        #--------------------------------
        int_fourier_vmr_m   =  []    #intercept Fourier analysis
        sl_fourier_vmr_m    =  []    #slope Fourier analysis
        p_fourier_vmr_m     =  []
        FAT_vmr_m           =  []    #Fitted Annual Trend
        FATAV_vmr_m         =  []    #Fitted Annual Trend and annual variability
        int_boot_vmr_m      =  []    #intercept Fourier + boot analysis
        sl_boot_vmr_m       =  []    #slope Fourier + boot analysis
        p_boot_vmr_m        =  []
        Rate_vmr_m          =  []    #Annual Rate of change
        Rate_e_vmr_m        =  []    #Uncertainty of the Annual Rate of change



    if pltInputs['weatherFlg']:
        wdir, wspeed, temp, rh, dt_weather = dc.weatherout(pltInputs['loc'], pltInputs['weatherDir'], pltInputs['weatherFileTag'], pltInputs['iyear'], pltInputs['fyear'] )

    #----------------------

    for k in range(ngases):

        retDir = '/data1/ebaumer/'+pltInputs['loc'].lower()+'/'+pltInputs['gasName'][k].lower()+'/'+pltInputs['ver'][k]+'/'
        ctlFile  = '/data1/ebaumer/'+pltInputs['loc'].lower()+'/'+pltInputs['gasName'][k].lower()+'/'+'x.'+pltInputs['gasName'][k].lower()+'/'+pltInputs['ctlF'][k]
        maxRMS = pltInputs['maxRMS'][k]
        minDOF = pltInputs['minDOF'][k]

    #---------------------------------
    # Check for the existance of files 
    # directories from input file
    #---------------------------------
        ckDir(retDir, exit=True)
        ckFile(ctlFile, exit=True)

        if pltInputs['saveFlg']:  ckDir(os.path.dirname(os.path.realpath(pltInputs['pltFile'])),exit=True)

    #-------------------------
    # Create Instance of Class
    #-------------------------
        gas = dc.PlotData(retDir,ctlFile,iyear=pltInputs['iyear'],imnth=pltInputs['imnth'],iday=pltInputs['iday'],
        fyear=pltInputs['fyear'],fmnth=pltInputs['fmnth'],fday=pltInputs['fday'],outFname=pltInputs['pltFile'])

    #----------------------
    # Call to get column data
    #----------------------
        ds = gas.pltTotClmn(fltr=pltInputs['fltrFlg'],sclfct=pltInputs['sclfct'],sclname=pltInputs['sclfctName'],mnthFltr=pltInputs["mnths"],mnthFltFlg=pltInputs["mnthFlg"],
                   partialCols=pltInputs['pCols'],errFlg=pltInputs['errorFlg'],minSZA=pltInputs['minSZA'],minTC=pltInputs['minTC'],maxTC=pltInputs['maxTC'],
                   maxSZA=pltInputs['maxSZA'],maxRMS=maxRMS,minDOF=minDOF,maxCHI=pltInputs['maxCHI'],
                   dofFlg=pltInputs['dofFlg'],rmsFlg=pltInputs['rmsFlg'],tcFlg=pltInputs['tcNegFlg'],
                   pcFlg=pltInputs['pcNegFlg'],szaFlg=pltInputs['szaFlg'],chiFlg=pltInputs['chiFlg'],cnvrgFlg=pltInputs['cnvrgFlg'],tcMMflg=pltInputs['tcMMFlg'])

        
        TC.append(ds['totClmn'])
        DT.append(ds['dates'])
        DOF.append(ds['dofs'])
        RMS.append(ds['rms'])
        TCRE.append(ds['tot_rnd'])
        TCSE.append(ds['tot_sys'])
        TCTE.append(ds['tot_std'])

        if gas.PrimaryGas   == 'C2H6':     gasname = 'C$_{2}$H$_{6}$'
        elif gas.PrimaryGas == 'CH4':      gasname = 'CH$_4$'
        elif gas.PrimaryGas == 'C2H2':     gasname = 'C$_2$H$_2$'
        elif gas.PrimaryGas == 'NH3':      gasname = 'NH$_3$'
        elif gas.PrimaryGas == 'O3':       gasname = 'O$_3$'
        elif gas.PrimaryGas == 'H2CO':     gasname = 'CH$_2$O'
        elif gas.PrimaryGas == 'HNO3':     gasname = 'HNO$_3$'
        elif gas.PrimaryGas == 'N2O':      gasname = 'N$_2$O'
        elif gas.PrimaryGas == 'CLONO2':   gasname = 'ClONO$_2$'
        else: gasname = gas.PrimaryGas



        #-----------------------------------
        # trend analysis of time series of actual data
        #-----------------------------------
        dateYearFrac = dc.toYearFraction(ds['dates'])
        weights      = np.ones_like(dateYearFrac)
        res          = dc.fit_driftfourier(dateYearFrac, ds['totClmn'], weights, 2)
        intercept, slope, pfourier = res[0:3]
        f_drift, f_fourier, f_driftfourier = res[3:6]

        int_fourier.append(intercept)
        sl_fourier.append(slope)
        p_fourier.append(pfourier)
        FAT.append(f_drift(dateYearFrac))
        FATAV.append(f_driftfourier(dateYearFrac))

        #-----------------------------------
        # bootstrap resampling information
        #-----------------------------------
        perc, intercept_boot, slope_boot, pfourier_boot = dc.cf_driftfourier(dateYearFrac, ds['totClmn'], weights, 2)

        int_boot.append(intercept_boot)
        sl_boot.append(slope_boot)
        p_boot.append(pfourier_boot)

        Rate.append(res[1]/np.mean(ds['totClmn'])*100.0)
        Rate_e.append(np.std(slope_boot)/np.mean(ds['totClmn'])*100.0)
        gasname_w.append(gasname)

        #----------------------------------------------------------------------
        #analysis of time series of DAILY averages
        #----------------------------------------------------------------------
        dailyVals    = dc.dailyAvg(ds['totClmn'],ds['dates'],dateAxis=1, meanAxis=0)
        dateYearFrac = dc.toYearFraction(dailyVals['dates'])
        weights      = np.ones_like(dateYearFrac)
        res          = dc.fit_driftfourier(dateYearFrac, dailyVals['dailyAvg'], weights, 2)
        intercept, slope, pfourier = res[0:3]
        f_drift, f_fourier, f_driftfourier = res[3:6]

        DT_d.append(dailyVals['dates'])
        TC_d.append(dailyVals['dailyAvg'])
        TCsd_d.append(dailyVals['std'])
        
        int_fourier_d.append(intercept)
        sl_fourier_d.append(slope)
        p_fourier_d.append(pfourier)
        FAT_d.append(f_drift(dateYearFrac))
        FATAV_d.append(f_driftfourier(dateYearFrac))

        
        #-----------------------------------
        # bootstrap resampling information of DAILY averages
        #-----------------------------------
        perc, intercept_boot, slope_boot, pfourier_boot = dc.cf_driftfourier(dateYearFrac, dailyVals['dailyAvg'], weights, 2)

        int_boot_d.append(intercept_boot)
        sl_boot_d.append(slope_boot)
        p_boot_d.append(pfourier_boot)

        Rate_d.append(res[1]/np.mean(dailyVals['dailyAvg'])*100.0)
        Rate_e_d.append(np.std(slope_boot)/np.mean(dailyVals['dailyAvg'])*100.0)

        #----------------------------------------------------------------------
        #analysis of time series of MONTHLY averages
        #----------------------------------------------------------------------
        mnthlyVals   = dc.mnthlyAvg(ds['totClmn'],ds['dates'],dateAxis=1, meanAxis=0)
        dateYearFrac = dc.toYearFraction(mnthlyVals['dates'])
        weights      = np.ones_like(dateYearFrac)
        res          = dc.fit_driftfourier(dateYearFrac, mnthlyVals['mnthlyAvg'], weights, 2)
        intercept, slope, pfourier = res[0:3]
        f_drift, f_fourier, f_driftfourier = res[3:6]

        DT_m.append(mnthlyVals['dates'])
        TC_m.append(mnthlyVals['mnthlyAvg'])
        TCsd_m.append(mnthlyVals['std'])
        
        int_fourier_m.append(intercept)
        sl_fourier_m.append(slope)
        p_fourier_m.append(pfourier)
        FAT_m.append(f_drift(dateYearFrac))
        FATAV_m.append(f_driftfourier(dateYearFrac))

        #-----------------------------------
        # bootstrap resampling information of MONTHLY averages
        #-----------------------------------
        perc, intercept_boot, slope_boot, pfourier_boot = dc.cf_driftfourier(dateYearFrac, mnthlyVals['mnthlyAvg'], weights, 2)

        int_boot_m.append(intercept_boot)
        sl_boot_m.append(slope_boot)
        p_boot_m.append(pfourier_boot)
        Rate_m.append(res[1]/np.mean(mnthlyVals['mnthlyAvg'])*100.0)
        Rate_e_m.append(np.std(slope_boot)/np.mean(mnthlyVals['mnthlyAvg'])*100.0)

        if pltInputs['pCols']:

            for pcol in pltInputs['pCols']:

                ind1 = dc.nearestind(pcol[0], ds['alt'])
                ind2 = dc.nearestind(pcol[1], ds['alt'])

                vmrP.append(np.average(ds['rPrf'][:,ind2:ind1],axis=1,weights=ds['Airmass'][:,ind2:ind1]) )            
                vmrPDry.append(np.average(ds['rPrfDry'][:,ind2:ind1],axis=1,weights=ds['Airmass'][:,ind2:ind1]) )
                sumP.append(np.sum(ds['rPrfMol'][:,ind2:ind1],axis=1))

               #----------------------------------------------------------------------
               #analysis of time series of DAILY averages of VMRP
               #----------------------------------------------------------------------
                dailyVals    = dc.dailyAvg(vmrPDry[k],ds['dates'],dateAxis=1, meanAxis=0)
                dateYearFrac = dc.toYearFraction(dailyVals['dates'])
                weights      = np.ones_like(dateYearFrac)
                res          = dc.fit_driftfourier(dateYearFrac, dailyVals['dailyAvg'], weights, 2)
                intercept, slope, pfourier = res[0:3]
                f_drift, f_fourier, f_driftfourier = res[3:6]

                vmr_d.append(dailyVals['dailyAvg'])
                vmrsd_d.append(dailyVals['std'])
        
                int_fourier_vmr_d.append(intercept)
                sl_fourier_vmr_d.append(slope)
                p_fourier_vmr_d.append(pfourier)
                FAT_vmr_d.append(f_drift(dateYearFrac))
                FATAV_vmr_d.append(f_driftfourier(dateYearFrac))

                #-----------------------------------
                # bootstrap resampling information of DAILY averages VMR
                #-----------------------------------
                perc, intercept_boot, slope_boot, pfourier_boot = dc.cf_driftfourier(dateYearFrac, dailyVals['dailyAvg'], weights, 2)

                int_boot_vmr_d.append(intercept_boot)
                sl_boot_vmr_d.append(slope_boot)
                p_boot_vmr_d.append(pfourier_boot)

                Rate_vmr_d.append(res[1]/np.mean(dailyVals['dailyAvg'])*100.0)
                Rate_e_vmr_d.append(np.std(slope_boot)/np.mean(dailyVals['dailyAvg'])*100.0)

                #----------------------------------------------------------------------
                #analysis of time series of MONTHLY averages
                #----------------------------------------------------------------------
                mnthlyVals   = dc.mnthlyAvg(vmrPDry[k],ds['dates'],dateAxis=1, meanAxis=0)
                dateYearFrac = dc.toYearFraction(mnthlyVals['dates'])
                weights      = np.ones_like(dateYearFrac)
                res          = dc.fit_driftfourier(dateYearFrac, mnthlyVals['mnthlyAvg'], weights, 2)
                intercept, slope, pfourier = res[0:3]
                f_drift, f_fourier, f_driftfourier = res[3:6]

                vmr_m.append(mnthlyVals['mnthlyAvg'])
                vmrsd_m.append(mnthlyVals['std'])
        
                int_fourier_vmr_m.append(intercept)
                sl_fourier_vmr_m.append(slope)
                p_fourier_vmr_m.append(pfourier)
                FAT_vmr_m.append(f_drift(dateYearFrac))
                FATAV_vmr_m.append(f_driftfourier(dateYearFrac))

                #-----------------------------------
                # bootstrap resampling information of DAILY averages VMR
                #-----------------------------------
                perc, intercept_boot, slope_boot, pfourier_boot = dc.cf_driftfourier(dateYearFrac, mnthlyVals['mnthlyAvg'], weights, 2)

                int_boot_vmr_m.append(intercept_boot)
                sl_boot_vmr_m.append(slope_boot)
                p_boot_vmr_m.append(pfourier_boot)

                Rate_vmr_m.append(res[1]/np.mean(mnthlyVals['mnthlyAvg'])*100.0)
                Rate_e_vmr_m.append(np.std(slope_boot)/np.mean(mnthlyVals['mnthlyAvg'])*100.0)


    #-----------------------------------
    # PLOT TOTAL COLUMN TIME SERIES OF SEVERAL TRACE GASES
    #-----------------------------------
    #gas.plt_ts_TotClmn(ngases, TC, DT, gasname_w, FAT, FATAV, Rate, Rate_e, vmrP, vmrPDry, sumP, sl_fourier,
    #                     TC_d, DT_d, FAT_d, FATAV_d, Rate_d, Rate_e_d, TCsd_d, sl_fourier_d,
    #                     TC_m, DT_m, FAT_m, FATAV_m, Rate_m, Rate_e_m, TCsd_m, sl_fourier_m,
    #                     vmr_d, FAT_vmr_d, FATAV_vmr_d, Rate_vmr_d, Rate_e_vmr_d, vmrsd_d, sl_fourier_vmr_d,
    #                     vmr_m, FAT_vmr_m, FATAV_vmr_m, Rate_vmr_m, Rate_e_vmr_m, vmrsd_m, sl_fourier_vmr_m,
    #                     pltInputs['iyear'],pltInputs['imnth'],pltInputs['iday'], pltInputs['fyear'],pltInputs['fmnth'],pltInputs['fday'],
    #                     yrsFlg=ds['yrsFlg'], pCols=pltInputs['pCols'], weatherFlag=pltInputs['weatherFlg'], wdir=wdir, wspeed=wspeed, temp=temp, rh=rh, dt_weather=dt_weather)

    #-----------------------------------
    # PLOT CORRELATION PLOTS VERSUS CO USING DAILY VALUES
    #-----------------------------------
    #gas.pltcorr(ngases, pltInputs['gasName'], TC_d, DT_d, gasname_w, pltInputs['iyear'],pltInputs['imnth'],pltInputs['iday'], pltInputs['fyear'],pltInputs['fmnth'],pltInputs['fday'],yrsFlg=ds['yrsFlg'])

    #-----------------------------------
    # PLOT CORRELATION PLOTS VERSUS CO USING WIND DATA AND ALL DATA
    #-----------------------------------
    gas.pltcorrWind(ngases, gasname_w, TC, DT, gasname_w, weatherFlag=pltInputs['weatherFlg'], wdir=wdir, wspeed=wspeed, temp=temp, rh=rh, dt_weather=dt_weather)

    #-----------------------------------
    # PLOT Q-Q PLOT VERSUS CO AND USING WIND DATA
    #-----------------------------------
    gas.pltqq(ngases, pltInputs['gasName'], TC, DT, TC_d, DT_d, gasname_w, weatherFlag=pltInputs['weatherFlg'], wdir=wdir, wspeed=wspeed, temp=temp, rh=rh, dt_weather=dt_weather)

    ##----
    if pltInputs['saveFlg']: gas.closeFig()

 #    #--------------------------------
 #    # Pause so user can look at plots
 #    #--------------------------------
    if not pltInputs['saveFlg']:
        user_input = raw_input('Press any key to exit >>> ')
        sys.exit()           # Exit program        


if __name__ == "__main__":
    main(sys.argv[1:])
