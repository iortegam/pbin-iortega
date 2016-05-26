#! /usr/local/python-2.7/bin/python
#----------------------------------------------------------------------------------------
# Name:
#        pltdofth.py
#
# Purpose:
#
# Plot time the relationship of DOF and tropospheric height
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
from matplotlib.pyplot import cm 
import numpy as np
from scipy import stats



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

def toYearFraction(dates):
    ''' Convert datetime to year and fraction of year'''

    #def sinceEpoch(date): # returns seconds since epoch
        #return time.mktime(date.timetuple())
    #s = sinceEpoch
    ep_fnc = lambda x: time.mktime(x.timetuple())
    
    retrnDates = np.zeros(len(dates))
    
    for i,sngDate in enumerate(dates):
        year = sngDate.year
        startOfThisYear = dt.datetime(year=year, month=1, day=1)
        startOfNextYear = dt.datetime(year=year+1, month=1, day=1)
    
        yearElapsed = ep_fnc(sngDate) - ep_fnc(startOfThisYear)
        yearDuration = ep_fnc(startOfNextYear) - ep_fnc(startOfThisYear)
        fraction = yearElapsed/yearDuration
        retrnDates[i] = sngDate.year + fraction

    return retrnDates    

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
    SZA      =  []  #TOTAL COLUMN OF MAIN TRACE GAS
    TH       =  []  #TOTAL COLUMN OF APRIORI OF TRACE GAS
    DT       =  []  #DATE and TIME
    DOF      =  []  #DEGREES OF FREEDOM
    ID       =  []  #ID version or NAME
    DOY      =  []
    GasList  =  []



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

        if gas.PrimaryGas == 'C2H6': gasname = 'C$_2$H$_6$'
        elif gas.PrimaryGas == 'CH4': gasname = 'CH$_4$'
        elif gas.PrimaryGas == 'C2H2': gasname = 'C$_2$H$_2$'
        elif gas.PrimaryGas == 'NH3': gasname = 'NH$_3$'
        elif gas.PrimaryGas == 'O3': gasname = 'O$_3$'
        elif gas.PrimaryGas == 'H2CO': gasname = 'CH$_2$O'
        elif gas.PrimaryGas == 'HCL': gasname = 'HCl'
        elif gas.PrimaryGas == 'HNO3': gasname = 'HNO$_3$'
        elif gas.PrimaryGas == 'N2O': gasname = 'N$_2$O'
        elif gas.PrimaryGas == 'CLONO2': gasname = 'ClONO$_2$'
        elif gas.PrimaryGas == 'H2O': gasname = 'H$_2$O'
        elif gas.PrimaryGas == 'NO2': gasname = 'NO$_2$'
        else: gasname = gas.PrimaryGas


    #----------------------
    # Call to plot profiles
    #----------------------
        ds = gas.pltTotClmn(fltr=pltInputs['fltrFlg'],sclfct=pltInputs['sclfct'],sclname=pltInputs['sclfctName'],mnthFltr=pltInputs["mnths"],mnthFltFlg=pltInputs["mnthFlg"],
                   partialCols=pltInputs['pCols'],errFlg=pltInputs['errorFlg'],minSZA=pltInputs['minSZA'],minTC=pltInputs['minTC'],maxTC=pltInputs['maxTC'],
                   maxSZA=pltInputs['maxSZA'],maxRMS=maxRMS,minDOF=minDOF,maxCHI=pltInputs['maxCHI'],
                   dofFlg=pltInputs['dofFlg'],rmsFlg=pltInputs['rmsFlg'],tcFlg=pltInputs['tcNegFlg'],
                   pcFlg=pltInputs['pcNegFlg'],szaFlg=pltInputs['szaFlg'],chiFlg=pltInputs['chiFlg'],cnvrgFlg=pltInputs['cnvrgFlg'],tcMMflg=pltInputs['tcMMFlg'],
                   allGas = True)

        SZA.append(ds['sza'])
        DOF.append(ds['dofs'])
        TH.append(ds['TH'])
        DT.append(ds['dates'])
        DOY.append(toYearFraction(ds['dates']))

    #----------------------------------------------------------
    #READ THE  DAYLY TROPOSPHERIC HEIGHT USING YEARLY FILES (ERIC CREATED)
    #----------------------------------------------------------
    NCEPhgtDir  = '/data1/ancillary_data/NCEPdata/NCEP_trpp/'

    TH_ncep = []
    TP_ncep = []
    dt_ncep = []
    doy_ncep = []

    nyears = int(pltInputs['fyear']) - int(pltInputs['iyear']) + 1

    for i in range(nyears):

        f2 = open( NCEPhgtDir +  'TropHght_'+pltInputs['loc'].lower()+'_'+str( int(pltInputs['iyear']) +i)+ '.dat', 'r')
        lines = f2.readlines()
        f2.close()

        c1 = []
        c2 = []
        c3 = []
        dt_temp = []
    
        count = 0
        for line in lines:
            if count >=1:
                p = line.split()
                yyyy = int(p[0][0:4])
                mm   = int(p[0][4:6])
                dd   = int(p[0][6:8])        
                c1.append(str(p[0]))
                c2.append(float(p[1]))
                c3.append(float(p[1]))
                dt_temp.append(dt.date(yyyy, mm, dd))
            count+= 1

        TH_ncep.extend(np.asarray(c2))
        TP_ncep.extend(np.asarray(c3))
        dt_ncep.extend(np.asarray(dt_temp))
        doy_ncep.extend(toYearFraction(np.asarray(dt_temp)))

    #----------------------------------------------------------

    clmap        = 'jet'
    cm           = plt.get_cmap(clmap)    
    yearsLc      = YearLocator()
    monthsAll    = MonthLocator()
    months       = MonthLocator()
   
    DateFmt      = DateFormatter('%m\n%Y')

    colors = ['r', 'b', 'g', 'y', 'c', 'k']

    #----------------------------------------------------------
    #SZA vs DOF
    #----------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8,5), sharex=True)
    
    for k in range(ngases):
        ax.scatter(SZA[k], DOF[k],facecolors='white', s=35, color=colors[k])

        #-----linear correlation analysis
        A = np.vstack([SZA[k], np.ones(len(SZA[k]))]).T
        slope, intercept = np.linalg.lstsq(A, DOF[k])[0]
        fit = slope*np.array(SZA[k]) + intercept
        ax.plot(SZA[k], fit, 'k', linewidth=2, label='Fitted line')
        ax.text(0.7,0.17, "slope= {0:4.3f}".format(slope),transform=ax.transAxes, fontsize = 16)
        ax.text(0.7,0.1,  "Intercept= {0:4.2f}".format(intercept),transform=ax.transAxes, fontsize = 16)
        slope2, intercept2, r_value, p_value, std_err = stats.linregress(SZA[k], DOF[k])
        ax.text(0.7,0.03,  "R$^2$= {0:4.2f}".format(r_value**2),transform=ax.transAxes, fontsize = 16)
        ax.annotate(gasname, xy=(.05, .93), xycoords='axes fraction', fontsize=30,
                    horizontalalignment='left', verticalalignment='center')


    ax.set_ylabel('DOF', fontsize=16)
    ax.set_xlabel('SZA [$^\circ$]', fontsize=16)
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.grid(True)
    fig.subplots_adjust(left = 0.12, bottom=0.12, top=0.96, right = 0.95)

    #----------------------------------------------------------
    #PDF of SZAs
    #----------------------------------------------------------
    fig2, ax = plt.subplots(figsize=(8,5), sharex=True)


    for k in range(ngases):

        #ax.hist(SZA[k], color=colors[k], facecolor=colors[k], alpha=0.5)
        ax.hist(SZA[k], color=colors[k], histtype='step')

    ax.set_ylabel('Density per bin', fontsize=16)   
    ax.set_xlabel('SZA [$^\circ$]', fontsize=16)
    ax.tick_params(axis='both',which='both',labelsize=8)
    ax.tick_params(axis='both', which='major', labelsize=14)         
    ax.grid(True)
    ax.annotate(gasname, xy=(.05, .93), xycoords='axes fraction', fontsize=30,
                    horizontalalignment='left', verticalalignment='center')
    fig2.subplots_adjust(left = 0.12, bottom=0.12, top=0.96, right = 0.95)

    #----------------------------------------------------------
    #Tropospheric height vs DOF
    #----------------------------------------------------------
    fig3, ax = plt.subplots(figsize=(8,5), sharex=True)
    
    for k in range(ngases):
        dt_n = DT[k]
        dt_n2 = []

        for p, dd in enumerate(dt_n):
            dt_n2.append( dt.date(dd.year, dd.month, dd.day))

        #dt_n2 = list(set(dt_n2))
        dt_n2 = np.asarray(dt_n2)
        
        TH_coi = []

        for p, dd in enumerate(dt_n2):

            diff = np.asarray(dt_ncep) - dd         
            inds = np.where( diff == dt.timedelta(0)  )[0]
            TH_coi.append(TH_ncep[inds])

        TH_coi = np.asarray(TH_coi)


        #TH_ncep_i  = np.interp(DOY[k], doy_ncep, TH_ncep)
        #ax.scatter( TH_ncep_i[indsSZA],DOF[k][indsSZA], facecolors='white', s=35, color=colors[k])

        ax.scatter( TH_coi,DOF[k], facecolors='white', s=35, color=colors[k])


        #-----linear correlation analysis
        #A = np.vstack([SZA[k], np.ones(len(SZA[k]))]).T
        #slope, intercept = np.linalg.lstsq(A, DOF[k])[0]
        #fit = slope*np.array(SZA[k]) + intercept
        #ax.plot(SZA[k], fit, 'k', linewidth=2, label='Fitted line')
        #ax.text(0.7,0.17, "slope= {0:4.3f}".format(slope),transform=ax.transAxes, fontsize = 16)
        #ax.text(0.7,0.1,  "Intercept= {0:4.2f}".format(intercept),transform=ax.transAxes, fontsize = 16)
        #slope2, intercept2, r_value, p_value, std_err = stats.linregress(SZA[k], DOF[k])
        #ax.text(0.7,0.03,  "R$^2$= {0:4.2f}".format(r_value**2),transform=ax.transAxes, fontsize = 16)
        
    ax.annotate(gasname, xy=(.05, .93), xycoords='axes fraction', fontsize=30,
                    horizontalalignment='left', verticalalignment='center')
    ax.set_ylabel('DOF', fontsize=16)
    ax.set_xlabel('Tropospheric Height [km]', fontsize=16)
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.grid(True)
    ax.set_title('All SZA', fontsize=22)
    fig3.subplots_adjust(left = 0.12, bottom=0.12, top=0.93, right = 0.95)

    #----------------------------------------------------------
    #Tropospheric height vs DOF  - 2
    #----------------------------------------------------------
    fig4, ax = plt.subplots(figsize=(8,5), sharex=True)
    indsSZA = np.where( (SZA[k] >= 65) & (SZA[k] <= 75) )[0]
   
    for k in range(ngases):
        dt_n = DT[k][indsSZA]
        dt_n2 = []

        for p, dd in enumerate(dt_n):
            dt_n2.append( dt.date(dd.year, dd.month, dd.day))

        #dt_n2 = list(set(dt_n2))
        dt_n2 = np.asarray(dt_n2)
        
        TH_coi = []

        for p, dd in enumerate(dt_n2):

            diff = np.asarray(dt_ncep) - dd         
            inds = np.where( diff == dt.timedelta(0)  )[0]
            TH_coi.append(TH_ncep[inds])

        TH_coi = np.asarray(TH_coi)


        #TH_ncep_i  = np.interp(DOY[k], doy_ncep, TH_ncep)
        #ax.scatter( TH_ncep_i[indsSZA],DOF[k][indsSZA], facecolors='white', s=35, color=colors[k])

        ax.scatter( TH_coi,DOF[k][indsSZA], facecolors='white', s=35, color=colors[k])

        #-----linear correlation analysis
        #A = np.vstack([SZA[k], np.ones(len(SZA[k]))]).T
        #slope, intercept = np.linalg.lstsq(A, DOF[k])[0]
        #fit = slope*np.array(SZA[k]) + intercept
        #ax.plot(SZA[k], fit, 'k', linewidth=2, label='Fitted line')
        #ax.text(0.7,0.17, "slope= {0:4.3f}".format(slope),transform=ax.transAxes, fontsize = 16)
        #ax.text(0.7,0.1,  "Intercept= {0:4.2f}".format(intercept),transform=ax.transAxes, fontsize = 16)
        #slope2, intercept2, r_value, p_value, std_err = stats.linregress(SZA[k], DOF[k])
        #ax.text(0.7,0.03,  "R$^2$= {0:4.2f}".format(r_value**2),transform=ax.transAxes, fontsize = 16)
        
    ax.annotate(gasname, xy=(.05, .93), xycoords='axes fraction', fontsize=30,
                    horizontalalignment='left', verticalalignment='center')
    ax.set_ylabel('DOF', fontsize=16)
    ax.set_xlabel('Tropospheric Height [km]', fontsize=16)
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.set_title('{0:2}< SZA <{1:2}'.format(65, 75), fontsize=22)
    ax.grid(True)
    fig4.subplots_adjust(left = 0.12, bottom=0.12, top=0.93, right = 0.95)



 

    #fig.subplots_adjust(left=0.07, bottom=0.1, top=0.95, right= 0.95)        
    

    
    # #Linear Correlation analysis
    # for k in range(ngases):
    #     if pltInputs['ID'][k] == pltInputs['ver_vs']: indT = k
     
    # dtfc_fix     = dc.toYearFraction(DT[indT])
    # tc_fic       = TC[indT]
    # rms_fic      = RMS[indT]
    # DOF_fic      = DOF[indT]
    # TCRE_fic     = TCRE[indT]
    # TCSE_fic     = TCSE[indT]
    # ratio_fic    = TCall[indT]['H2O']/TCallApr[indT]['H2O']

    # one2one = np.arange(np.min(TC[indT]) - 1.0*np.min(TC[indT]), np.max(TC[indT])+1.0*np.max(TC[indT]), np.min(TC[indT])*0.05)

    # fig2, ax = plt.subplots(2,3, figsize=(14,10))

    # for k in range(ngases):

    #     if k != indT:

    #         dtfc = dc.toYearFraction(DT[k])
    #         tci  = np.interp(dtfc_fix, dtfc, TC[k])
    #         tcei = np.interp(dtfc_fix, dtfc, TCTE[k])
    #         rmsi = np.interp(dtfc_fix, dtfc, RMS[k])
    #         dofi = np.interp(dtfc_fix, dtfc, DOF[k])
    #         tcrei = np.interp(dtfc_fix, dtfc, TCRE[k])
    #         tcsei = np.interp(dtfc_fix, dtfc, TCSE[k])
    #         ratioi = np.interp(dtfc_fix, dtfc, TCall[k]['H2O']/TCallApr[k]['H2O'])
           
    #         ax[0,0].scatter(tc_fic,  tci, facecolors='white', s=35, color=colors[k], label = ID[k]) 
    #         #-----linear correlation analysis
    #         A = np.vstack([tc_fic, np.ones(len(tc_fic))]).T
    #         slope, intercept = np.linalg.lstsq(A, tci)[0]
    #         fit = slope*np.array(tc_fic) + intercept
    #         ax[0,0].plot(tc_fic, fit, 'k', linewidth=2, label='Fitted line')
    #         ax[0,0].text(0.4,0.1, "slope: {0:4.2f}".format(slope),transform=ax[0,0].transAxes, fontsize = 14)
    #         ax[0,0].text(0.4,0.04,  "Intercept: {0:.2E}".format(intercept),transform=ax[0,0].transAxes, fontsize = 14) 
    #         #---------------
            
    #         ax[0,1].scatter(rms_fic,  rmsi, facecolors='white', s=35, color=colors[k])
    #         #-----linear correlation analysis
    #         A = np.vstack([rms_fic, np.ones(len(rms_fic))]).T
    #         slope, intercept = np.linalg.lstsq(A, rmsi)[0]
    #         fit = slope*np.array(rms_fic) + intercept
    #         ax[0,1].plot(rms_fic, fit, 'k', linewidth=2, label='Fitted line')
    #         ax[0,1].text(0.4,0.1, "slope: {0:4.2f}".format(slope),transform=ax[0,1].transAxes, fontsize = 14)
    #         ax[0,1].text(0.4,0.04,  "Intercept: {0:4.3f}".format(intercept),transform=ax[0,1].transAxes, fontsize = 14) 
    #         #---------------

    #         ax[0,2].scatter(DOF_fic,  dofi, facecolors='white', s=35, color=colors[k])
    #         #-----linear correlation analysis
    #         A = np.vstack([DOF_fic, np.ones(len(DOF_fic))]).T
    #         slope, intercept = np.linalg.lstsq(A, dofi)[0]
    #         fit = slope*np.array(DOF_fic) + intercept
    #         ax[0,2].plot(DOF_fic, fit, 'k', linewidth=2, label='Fitted line')
    #         ax[0,2].text(0.4,0.1, "slope: {0:4.2f}".format(slope),transform=ax[0,2].transAxes, fontsize = 14)
    #         ax[0,2].text(0.4,0.04,  "Intercept: {0:4.3f}".format(intercept),transform=ax[0,2].transAxes, fontsize = 14) 
    #         #---------------

    #         ax[1,0].scatter(TCRE_fic,  tcrei, facecolors='white', s=35, color=colors[k])
    #         #-----linear correlation analysis
    #         A = np.vstack([TCRE_fic, np.ones(len(TCRE_fic))]).T
    #         slope, intercept = np.linalg.lstsq(A, tcrei)[0]  
    #         fit = slope*np.array(TCRE_fic) + intercept
    #         ax[1,0].plot(TCRE_fic, fit, 'k', linewidth=2, label='Fitted line')
    #         ax[1,0].text(0.4,0.1, "slope: {0:4.2f}".format(slope),transform=ax[1,0].transAxes, fontsize = 14)
    #         ax[1,0].text(0.4,0.04,  "Intercept: {0:.2E}".format(intercept),transform=ax[1,0].transAxes, fontsize = 14) 
            
    #         #---------------
    #         ax[1,1].scatter(TCSE_fic,  tcsei, facecolors='white', s=35, color=colors[k])
    #         #-----linear correlation analysis
    #         A = np.vstack([TCSE_fic, np.ones(len(TCSE_fic))]).T
    #         slope, intercept = np.linalg.lstsq(A, tcsei)[0]  
    #         fit = slope*np.array(TCSE_fic) + intercept
    #         ax[1,1].plot(TCSE_fic, fit, 'k', linewidth=2, label='Fitted line')
    #         ax[1,1].text(0.4,0.1, "slope: {0:4.2f}".format(slope),transform=ax[1,1].transAxes, fontsize = 14)
    #         ax[1,1].text(0.4,0.04,  "Intercept: {0:.2E}".format(intercept),transform=ax[1,1].transAxes, fontsize = 14) 
    #         #---------------
    #         ax[1,2].scatter(ratio_fic,  ratioi, facecolors='white', s=35, color=colors[k])
    #         #-----linear correlation analysis
    #         A = np.vstack([ratio_fic, np.ones(len(ratio_fic))]).T
    #         slope, intercept = np.linalg.lstsq(A, ratioi)[0]  
    #         fit = slope*np.array(ratio_fic) + intercept
    #         ax[1,2].plot(ratio_fic, fit, 'k', linewidth=2, label='Fitted line')
    #         ax[1,2].text(0.4,0.1, "slope: {0:4.2f}".format(slope),transform=ax[1,2].transAxes, fontsize = 14)
    #         ax[1,2].text(0.4,0.04,  "Intercept: {0:4.3f}".format(intercept),transform=ax[1,2].transAxes, fontsize = 14) 
    #         #---------------

    
    # ax[0,0].set_xlabel('Total Column [molec$\cdot$cm$^{-2}$] - %s' %str(pltInputs['ver_vs']))
    # ax[0,0].set_ylabel('Total Column [molec$\cdot$cm$^{-2}$]')    
    # ax[0,0].legend(prop={'size':12})
    # ax[0,0].set_ylim(np.min(TC[indT])-0.5*np.min(TC[indT]), np.max(TC[indT])+0.15*np.max(TC[indT]))
    # ax[0,0].set_xlim(np.min(TC[indT])-0.5*np.min(TC[indT]), np.max(TC[indT])+0.15*np.max(TC[indT]))
    # ax[0,0].plot(one2one, one2one, ls ='--', color='gray', linewidth=2)
    # ax[0,0].plot(one2one, one2one*0.8, ls ='--', color='gray', linewidth=2)
    # ax[0,0].plot(one2one*0.8, one2one, ls ='--', color='gray', linewidth=2) 
    # ax[0,0].grid(True)

    # ax[0,1].set_xlabel('RMS - %s' %str(pltInputs['ver_vs']))
    # ax[0,1].set_ylabel('RMS')    
    # ax[0,1].set_ylim(np.min(RMS[indT])-0.5*np.min(RMS[indT]), np.max(RMS[indT])+0.15*np.max(RMS[indT]))
    # ax[0,1].set_xlim(np.min(RMS[indT])-0.5*np.min(RMS[indT]), np.max(RMS[indT])+0.15*np.max(RMS[indT]))
    # ax[0,1].plot(one2one, one2one, ls ='--', color='gray', linewidth=2)
    # ax[0,1].plot(one2one, one2one*0.8, ls ='--', color='gray', linewidth=2)
    # ax[0,1].plot(one2one*0.8, one2one, ls ='--', color='gray', linewidth=2) 
    # ax[0,1].grid(True)

    # ax[0,2].set_xlabel('DOF - %s' %str(pltInputs['ver_vs']))
    # ax[0,2].set_ylabel('DOF')    
    # ax[0,2].set_ylim(np.min(DOF[indT])-0.5*np.min(DOF[indT]), np.max(DOF[indT])+0.15*np.max(DOF[indT]))
    # ax[0,2].set_xlim(np.min(DOF[indT])-0.5*np.min(DOF[indT]), np.max(DOF[indT])+0.15*np.max(DOF[indT]))
    # ax[0,2].plot(one2one, one2one, ls ='--', color='gray', linewidth=2)
    # ax[0,2].plot(one2one, one2one*0.8, ls ='--', color='gray', linewidth=2)
    # ax[0,2].plot(one2one*0.8, one2one, ls ='--', color='gray', linewidth=2) 
    # ax[0,2].grid(True)

    # ax[1,0].set_xlabel('Random Error - %s' %str(pltInputs['ver_vs']))
    # ax[1,0].set_ylabel('Random Error')    
    # ax[1,0].set_ylim(np.min(TCRE[indT])-0.5*np.min(TCRE[indT]), np.max(TCRE[indT])+0.15*np.max(TCRE[indT]))
    # ax[1,0].set_xlim(np.min(TCRE[indT])-0.5*np.min(TCRE[indT]), np.max(TCRE[indT])+0.15*np.max(TCRE[indT]))
    # ax[1,0].plot(one2one, one2one, ls ='--', color='gray', linewidth=2)
    # ax[1,0].plot(one2one, one2one*0.8, ls ='--', color='gray', linewidth=2)
    # ax[1,0].plot(one2one*0.8, one2one, ls ='--', color='gray', linewidth=2) 
    # ax[1,0].grid(True)

    # ax[1,1].set_xlabel('Systematic Error - %s' %str(pltInputs['ver_vs']))
    # ax[1,1].set_ylabel('Systematic Error')    
    # ax[1,1].set_ylim(np.min(TCSE[indT])-0.5*np.min(TCSE[indT]), np.max(TCSE[indT])+0.15*np.max(TCSE[indT]))
    # ax[1,1].set_xlim(np.min(TCSE[indT])-0.5*np.min(TCSE[indT]), np.max(TCSE[indT])+0.15*np.max(TCSE[indT]))
    # ax[1,1].plot(one2one, one2one, ls ='--', color='gray', linewidth=2)
    # ax[1,1].plot(one2one, one2one*0.8, ls ='--', color='gray', linewidth=2)
    # ax[1,1].plot(one2one*0.8, one2one, ls ='--', color='gray', linewidth=2) 
    # ax[1,1].grid(True)

    # ax[1,2].set_xlabel('H$_2$O (ret)/H$_2$O (apr) - %s' %str(pltInputs['ver_vs']))
    # ax[1,2].set_ylabel('H$_2$O (ret)/H$_2$O (apr)')    
    # ax[1,2].set_ylim(np.min(ratio_fic)-0.5*np.min(ratio_fic), np.max(ratio_fic)+0.15*np.max(ratio_fic))
    # ax[1,2].set_xlim(np.min(ratio_fic)-0.5*np.min(ratio_fic), np.max(ratio_fic)+0.15*np.max(ratio_fic))
    # ax[1,2].plot(one2one, one2one, ls ='--', color='gray', linewidth=2)
    # ax[1,2].plot(one2one, one2one*0.8, ls ='--', color='gray', linewidth=2)
    # ax[1,2].plot(one2one*0.8, one2one, ls ='--', color='gray', linewidth=2) 
    # ax[1,2].grid(True)

    # fig2.subplots_adjust(left=0.075, bottom=0.1, top=0.95, right= 0.95)

    if gas.pdfsav: 
        gas.pdfsav.savefig(fig,dpi=200)
        gas.pdfsav.savefig(fig2,dpi=200)
        gas.pdfsav.savefig(fig3,dpi=200)
        gas.pdfsav.savefig(fig4,dpi=200)
       
    else:           
        plt.show(block=False)  


    #-------------------------------------------------------------------------
    #END PLOTS
    #-------------------------------------------------------------------------
    if pltInputs['saveFlg']: gas.closeFig()

    print('\nFinished Plots.......\n')
       
      #--------------------------------
      # Pause so user can look at plots
      #--------------------------------
    if not pltInputs['saveFlg']:
        user_input = raw_input('Press any key to exit >>> ')
        sys.exit()           # Exit program        


if __name__ == "__main__":
    main(sys.argv[1:])
