#! /usr/local/python-2.7/bin/python
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

    fig, ax = plt.subplots(ngases,figsize=(8,13), sharex=True)
    fig2, ax2 = plt.subplots(ngases,figsize=(8,13), sharex=True)
    fig3, ax3 = plt.subplots(ngases,figsize=(8,13), sharex=True)
    fig4, ax4 = plt.subplots(ngases,figsize=(8,13), sharex=True)

    Rate      =  []
    Rate_e    =  []
    xf        =  []
    gasname_w =  []



    for k in range(ngases):

        retDir = '/data1/ebaumer/'+pltInputs['loc'].lower()+'/'+pltInputs['gasName'][k].lower()+'/'+pltInputs['ver'][k]+'/'
        ctlFile  = '/data1/ebaumer/'+pltInputs['loc'].lower()+'/'+pltInputs['gasName'][k].lower()+'/'+'x.'+pltInputs['gasName'][k].lower()+'/'+pltInputs['ctlF'][k]
        maxRMS = pltInputs['maxRMS'][k]
        minDOF = pltInputs['minDOF'][k]


    #---------------------------------
    # Check for the existance of files 
    # directories from input file
    #---------------------------------
       # ckDir(pltInputs['retDir'], exit=True)
       # ckFile(pltInputs['ctlFile'], exit=True)

        ckDir(retDir, exit=True)
        ckFile(ctlFile, exit=True)

        if pltInputs['saveFlg']:  ckDir(os.path.dirname(os.path.realpath(pltInputs['pltFile'])),exit=True)

     #-------------------------
     # Create Instance of Class
     #-------------------------
        gas = dc.PlotData(retDir,ctlFile,iyear=pltInputs['iyear'],imnth=pltInputs['imnth'],iday=pltInputs['iday'],
        fyear=pltInputs['fyear'],fmnth=pltInputs['fmnth'],fday=pltInputs['fday'],outFname=pltInputs['pltFile'])

 #    #----------------------
 #    # Call to plot profiles
 #    #----------------------
        ds = gas.pltTotClmn(fltr=pltInputs['fltrFlg'],sclfct=pltInputs['sclfct'],sclname=pltInputs['sclfctName'],mnthFltr=pltInputs["mnths"],mnthFltFlg=pltInputs["mnthFlg"],
                   partialCols=pltInputs['pCols'],errFlg=pltInputs['errorFlg'],minSZA=pltInputs['minSZA'],minTC=pltInputs['minTC'],maxTC=pltInputs['maxTC'],
                   maxSZA=pltInputs['maxSZA'],maxRMS=maxRMS,minDOF=minDOF,maxCHI=pltInputs['maxCHI'],
                   dofFlg=pltInputs['dofFlg'],rmsFlg=pltInputs['rmsFlg'],tcFlg=pltInputs['tcNegFlg'],
                   pcFlg=pltInputs['pcNegFlg'],szaFlg=pltInputs['szaFlg'],chiFlg=pltInputs['chiFlg'],cnvrgFlg=pltInputs['cnvrgFlg'],tcMMflg=pltInputs['tcMMFlg'])

        #-------------------------------------------------------------------------
        #PLOTS
        #-------------------------------------------------------------------------
        #raw_input('Continue?')

        print '\nPrinting Plots: %s' %gas.PrimaryGas.lower(), '\n'

        if gas.PrimaryGas == 'C2H6': gasname = 'C$_2$H$_6$'
        elif gas.PrimaryGas == 'CH4': gasname = 'CH$_4$'
        elif gas.PrimaryGas == 'C2H2': gasname = 'C$_2$H$_2$'
        elif gas.PrimaryGas == 'NH3': gasname = 'NH$_3$'
        elif gas.PrimaryGas == 'O3': gasname = 'O$_3$'
        elif gas.PrimaryGas == 'H2CO': gasname = 'CH$_2$O'
        else: gasname = gas.PrimaryGas



        #print ds['sys_cmpnts']
        #exit()
        #-------------------------------------------------------
        # Plot mean systematic and random errors on mean profiles
        #-------------------------------------------------------
  
        clmap        = 'jet'
        cm           = plt.get_cmap(clmap)    
        yearsLc      = YearLocator()
        monthsAll    = MonthLocator()
        months       = MonthLocator()
        DateFmt      = DateFormatter('%m\n%Y')
        #DateFmt      = DateFormatter('%Y')     

        #-----------------------------------
        # Plot trend analysis of time series of actual data
        #-----------------------------------
        dateYearFrac = dc.toYearFraction(ds['dates'])
        weights      = np.ones_like(dateYearFrac)
        res          = dc.fit_driftfourier(dateYearFrac, ds['totClmn'], weights, 2)
        intercept, slope, pfourier = res[0:3]
        f_drift, f_fourier, f_driftfourier = res[3:6]

        print 'intercept: %s' %intercept
        print 'slope: %s' %slope
        print 'pfourier: %s' %pfourier

        #-----------------------------------
        # bootstrap resampling information
        #-----------------------------------
        perc, intercept_boot, slope_boot, pfourier_boot = dc.cf_driftfourier(dateYearFrac, ds['totClmn'], weights, 2)
       
        #-----------------------------------
        ax[k].scatter(ds['dates'],ds['totClmn'],facecolors='white', edgecolors='gray',s=35)
        #ax[k].errorbar(ds['dates'],ds['totClmn'],yerr=ds['tot_std'],fmt='k.',markersize=4,ecolor='red')
        ax[k].plot(ds['dates'],f_drift(dateYearFrac),label='Fitted Anual Trend',linewidth=2.5)
        ax[k].plot(ds['dates'],f_driftfourier(dateYearFrac),label='Fitted Anual Trend + intra-annual variability', linewidth=2.5)
        ax[k].set_title(gasname)
        ax[k].grid(True)
        ax[k].set_ylim([np.min(ds['totClmn'])-0.05*np.min(ds['totClmn']), np.max(ds['totClmn'])+0.1*np.max(ds['totClmn'])])
        ##ax[k].set_ylabel('Total Column$\cdot$[molecules cm$^{-2}$]')
        ax[k].text(0.02,0.85,"Fitted trend -- slope: {0:.3E} ({1:.2f}%)".format(res[1],res[1]/np.mean(ds['totClmn'])*100.0),transform=ax[k].transAxes, fontsize = 14) 
       #ax[k].text(0.02,0.9,"Fitted intercept at xmin: {:.3E}".format(res[0]),transform=ax[k].transAxes)
       # ax[k].text(0.02,0.86,"STD of residuals: {0:.3E} ({1:.3f}%)".format(res[6],res[6]/np.mean(ds['totClmn'])*100.0),transform=ax[k].transAxes) 

        Rate.append(res[1]/np.mean(ds['totClmn'])*100.0)
        Rate_e.append(np.std(slope_boot)/np.mean(ds['totClmn'])*100.0)
        xf.append(k)
        gasname_w.append(gasname)

            
        if ds['yrsFlg']:
            #plt.xticks(rotation=45)
            ax[k].xaxis.set_major_locator(yearsLc)
            ax[k].xaxis.set_minor_locator(months)
            ax[k].xaxis.set_major_formatter(DateFmt) 
            ax[k].xaxis.set_tick_params(which='major',labelsize=12)
            ax[k].xaxis.set_tick_params(which='minor',labelbottom='off')
            ax[k].set_xlim((dt.date(pltInputs['iyear'],pltInputs['imnth'], 1), dt.date(pltInputs['fyear'],pltInputs['fmnth'],pltInputs['fday'])))
            if k == ngases-1: ax[k].set_xlabel('Year', fontsize = 16)
        else:
            ax[k].xaxis.set_major_locator(monthsAll)
            ax[k].xaxis.set_major_formatter(DateFmt)
            ax[k].set_xlim((dt.date(ds['years'][0],1,1), dt.date(ds['years'][0],12,31)))
            ax[k].xaxis.set_minor_locator(AutoMinorLocator())
            if k == ngases-1: ax[k].set_xlabel('Month/Year')

            fig.autofmt_xdate()

        fig.text(0.075, 0.5, 'Total Column [molecules$\cdot$cm$^{-2}$]', fontsize = 16,
        horizontalalignment='right',
        verticalalignment='center',
        rotation='vertical',
        transform=ax[k].transAxes)

        fig.suptitle('All data', fontsize=18)
        fig.subplots_adjust(bottom=0.07,top=0.93)

        #----------------------------------------------------------------------
        # Plot trend analysis of time series of DAILY averages
        #----------------------------------------------------------------------
        dailyVals    = dc.dailyAvg(ds['totClmn'],ds['dates'],dateAxis=1, meanAxis=0)
        dateYearFrac = dc.toYearFraction(dailyVals['dates'])
        weights      = np.ones_like(dateYearFrac)
        res          = dc.fit_driftfourier(dateYearFrac, dailyVals['dailyAvg'], weights, 2)
        intercept, slope, pfourier = res[0:3]
        f_drift, f_fourier, f_driftfourier = res[3:6]
        
        #-----------------------------------
        ax2[k].scatter(dailyVals['dates'],dailyVals['dailyAvg'],facecolors='white', edgecolors='gray',s=35)
        ax2[k].plot(dailyVals['dates'],f_drift(dateYearFrac),label='Fitted Anual Trend',linewidth=2.5)
        ax2[k].plot(dailyVals['dates'],f_driftfourier(dateYearFrac),label='Fitted Anual Trend + intra-annual variability', linewidth=2.5)
        ax2[k].errorbar(dailyVals['dates'],dailyVals['dailyAvg'],yerr=dailyVals['std']*2.,fmt='k.',markersize=0,ecolor='grey')   
        ax2[k].set_title(gasname)
        ax2[k].grid(True)
        ax2[k].set_ylim([np.min(dailyVals['dailyAvg'])-0.1*np.min(dailyVals['dailyAvg']), np.max(dailyVals['dailyAvg'])+0.15*np.max(dailyVals['dailyAvg'])])
        ##ax2[k].set_ylabel('Total Column$\cdot$[molecules cm$^{-2}$]')
        ax2[k].text(0.02,0.85,"Fitted trend -- slope: {0:.3E} ({1:.2f}%)".format(res[1],res[1]/np.mean(dailyVals['dailyAvg'])*100.0),transform=ax2[k].transAxes, fontsize = 14) 
       #ax2[k].text(0.02,0.9,"Fitted intercept at xmin: {:.3E}".format(res[0]),transform=ax2[k].transAxes)
       # ax2[k].text(0.02,0.86,"STD of residuals: {0:.3E} ({1:.3f}%)".format(res[6],res[6]/np.mean(ds['totClmn'])*100.0),transform=ax2[k].transAxes) 

            
        if ds['yrsFlg']:
            #plt.xticks(rotation=45)
            ax2[k].xaxis.set_major_locator(yearsLc)
            ax2[k].xaxis.set_minor_locator(months)
            ax2[k].xaxis.set_major_formatter(DateFmt) 
            ax2[k].xaxis.set_tick_params(which='major',labelsize=12)
            ax2[k].xaxis.set_tick_params(which='minor',labelbottom='off')
            ax2[k].set_xlim((dt.date(pltInputs['iyear'],pltInputs['imnth'], 1), dt.date(pltInputs['fyear'],pltInputs['fmnth'],pltInputs['fday'])))
            if k == ngases-1: ax2[k].set_xlabel('Year', fontsize = 16)
        else:
            ax2[k].xaxis.set_major_locator(monthsAll)
            ax2[k].xaxis.set_major_formatter(DateFmt)
            ax2[k].set_xlim((dt.date(ds['years'][0],1,1), dt.date(ds['years'][0],12,31)))
            ax2[k].xaxis.set_minor_locator(AutoMinorLocator())
            if k == ngases-1: ax2[k].set_xlabel('Month/Year')

            fig2.autofmt_xdate()

        fig2.text(0.075, 0.5, 'Total Column [molecules$\cdot$cm$^{-2}$]', fontsize = 16,
        horizontalalignment='right',
        verticalalignment='center',
        rotation='vertical',
        transform=ax[k].transAxes)
        
        fig2.suptitle('Daily averages $\pm$ standard deviation', fontsize=18)
        fig2.subplots_adjust(bottom=0.07,top=0.93)

        #----------------------------------------------------------------------
        # Plot trend analysis of time series of MONTHLY averages
        #----------------------------------------------------------------------
        mnthlyVals   = dc.mnthlyAvg(ds['totClmn'],ds['dates'],dateAxis=1, meanAxis=0)
        dateYearFrac = dc.toYearFraction(mnthlyVals['dates'])
        weights      = np.ones_like(dateYearFrac)
        res          = dc.fit_driftfourier(dateYearFrac, mnthlyVals['mnthlyAvg'], weights, 2)
        intercept, slope, pfourier = res[0:3]
        f_drift, f_fourier, f_driftfourier = res[3:6]
        
        #-----------------------------------
        ax3[k].scatter(mnthlyVals['dates'],mnthlyVals['mnthlyAvg'],facecolors='white', edgecolors='gray',s=35)
        ax3[k].plot(mnthlyVals['dates'],f_drift(dateYearFrac),label='Fitted Anual Trend',linewidth=2.5)
        ax3[k].plot(mnthlyVals['dates'],f_driftfourier(dateYearFrac),label='Fitted Anual Trend + intra-annual variability', linewidth=2.5)
        ax3[k].errorbar(mnthlyVals['dates'],mnthlyVals['mnthlyAvg'],yerr=mnthlyVals['std'],fmt='k.',markersize=0,ecolor='grey')
        ax3[k].set_title(gasname)
        ax3[k].grid(True)
        ax3[k].set_ylim([np.min(mnthlyVals['mnthlyAvg'])-0.1*np.min(mnthlyVals['mnthlyAvg']), np.max(mnthlyVals['mnthlyAvg'])+0.15*np.max(dailyVals['dailyAvg'])])
        ##ax3[k].set_ylabel('Total Column$\cdot$[molecules cm$^{-2}$]')
        ax3[k].text(0.02,0.85,"Fitted trend -- slope: {0:.3E} ({1:.2f}%)".format(res[1],res[1]/np.mean(mnthlyVals['mnthlyAvg'])*100.0),transform=ax3[k].transAxes, fontsize = 14) 
       #ax3[k].text(0.02,0.9,"Fitted intercept at xmin: {:.3E}".format(res[0]),transform=ax3[k].transAxes)
       # ax3[k].text(0.02,0.86,"STD of residuals: {0:.3E} ({1:.3f}%)".format(res[6],res[6]/np.mean(ds['totClmn'])*100.0),transform=ax2[k].transAxes) 

            
        if ds['yrsFlg']:
            #plt.xticks(rotation=45)
            ax3[k].xaxis.set_major_locator(yearsLc)
            ax3[k].xaxis.set_minor_locator(months)
            ax3[k].xaxis.set_major_formatter(DateFmt) 
            ax3[k].xaxis.set_tick_params(which='major',labelsize=12)
            ax3[k].xaxis.set_tick_params(which='minor',labelbottom='off')
            ax3[k].set_xlim((dt.date(pltInputs['iyear'],pltInputs['imnth'], 1), dt.date(pltInputs['fyear'],pltInputs['fmnth'],pltInputs['fday'])))
            if k == ngases-1: ax3[k].set_xlabel('Year', fontsize = 16)
        else:
            ax3[k].xaxis.set_major_locator(monthsAll)
            ax3[k].xaxis.set_major_formatter(DateFmt)
            ax3[k].set_xlim((dt.date(ds['years'][0],1,1), dt.date(ds['years'][0],12,31)))
            ax3[k].xaxis.set_minor_locator(AutoMinorLocator())
            if k == ngases-1: ax3[k].set_xlabel('Month/Year')

            fig3.autofmt_xdate()

        fig3.text(0.075, 0.5, 'Total Column [molecules$\cdot$cm$^{-2}$]', fontsize = 16,
        horizontalalignment='right',
        verticalalignment='center',
        rotation='vertical',
        transform=ax3[k].transAxes)

        fig3.suptitle('Monthly averages $\pm$ standard deviation', fontsize=18)
        fig3.subplots_adjust(bottom=0.07,top=0.93)

    #------------------------------------
        # Plot time series of partial columns
        #------------------------------------
        if pltInputs['pCols']:

            for pcol in pltInputs['pCols']:

                ind1 = dc.nearestind(pcol[0], ds['alt'])
                ind2 = dc.nearestind(pcol[1], ds['alt'])

                vmrP = np.average(ds['rPrf'][:,ind2:ind1],axis=1,weights=ds['Airmass'][:,ind2:ind1])             
                vmrPDry = np.average(ds['rPrfDry'][:,ind2:ind1],axis=1,weights=ds['Airmass'][:,ind2:ind1]) 
                sumP = np.sum(ds['rPrfMol'][:,ind2:ind1],axis=1)
        
                ax4[k].scatter(ds['dates'],vmrPDry,facecolors='white', edgecolors='gray',s=35)
                ax4[k].set_title(gasname)
                ax4[k].grid(True)
                ax4[k].set_ylim(np.min(vmrPDry)-0.15*np.min(vmrPDry), np.max(vmrPDry)+0.15*np.max(vmrPDry))
                
                if ds['yrsFlg']:
                
                    ax4[k].xaxis.set_major_locator(yearsLc)
                    ax4[k].xaxis.set_minor_locator(months)
                    ax4[k].xaxis.set_major_formatter(DateFmt) 
                    ax4[k].xaxis.set_tick_params(which='major',labelsize=12)
                    ax4[k].xaxis.set_tick_params(which='minor',labelbottom='off')
                    ax4[k].set_xlim((dt.date(pltInputs['iyear'],pltInputs['imnth'], 1), dt.date(pltInputs['fyear'],pltInputs['fmnth'],pltInputs['fday'])))
                    if k == ngases-1: ax4[k].set_xlabel('Year', fontsize = 16)
                else:
                    ax4[k].xaxis.set_major_locator(monthsAll)
                    ax4[k].xaxis.set_major_formatter(DateFmt)
                    ax4[k].set_xlim((dt.date(ds['years'][0],1,1), dt.date(ds['years'][0],12,31)))
                    ax4[k].xaxis.set_minor_locator(AutoMinorLocator())
                    if k == ngases-1: ax4[k].set_xlabel('Month/Year')
                    
                    fig4.autofmt_xdate()

                fig4.text(0.075, 0.5, 'VMR [ppb$_{v}$] Dry Air', fontsize = 16,
                horizontalalignment='right',
                verticalalignment='center',
                rotation='vertical',
                transform=ax4[k].transAxes)
                
                fig4.suptitle('All data - weighted VMR ('+str(pltInputs['pCols'][0][0]) + '-' + str(pltInputs['pCols'][0][1])+ 'km)', fontsize=18)
                fig4.subplots_adjust(bottom=0.07,top=0.93)


   
    fig5, ax = plt.subplots()
    ax.bar(xf,Rate, align='center', color = 'r', yerr=Rate_e, ecolor = 'k')
    ax.set_xticks(xf)
    ax.set_ylabel('Annual rate of change (%)', fontsize = 16)
    ax.set_xticklabels(gasname_w)
    ax.set_xlabel('Gas', fontsize = 16)
    ax.axhline(0, color='black', lw=1)
    
              
    if gas.pdfsav: 
        gas.pdfsav.savefig(fig,dpi=200)
        gas.pdfsav.savefig(fig2,dpi=200)
        gas.pdfsav.savefig(fig3,dpi=200)
        gas.pdfsav.savefig(fig4,dpi=200)
        gas.pdfsav.savefig(fig5,dpi=200)
       
    else:           
        plt.show(block=False)          
                    

    #-------------------------------------------------------------------------
    #END PLOTS
    #-------------------------------------------------------------------------
 # ##----
    if pltInputs['saveFlg']: gas.closeFig()

    print('\nFinished Plots.......\n')
       

 #    #--------------------------------
 #    # Pause so user can look at plots
 #    #--------------------------------
    if not pltInputs['saveFlg']:
        user_input = raw_input('Press any key to exit >>> ')
        sys.exit()           # Exit program        


if __name__ == "__main__":
    main(sys.argv[1:])
