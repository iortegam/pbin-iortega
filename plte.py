#! /usr/local/python-2.7/bin/python
#----------------------------------------------------------------------------------------
#
# Purpose:
#
# Plot retrieval errors for several trace gases defined in plte_input.py
#----------------------------------------------------------------------------------------


#---------------
# Import modules
#---------------
import sys
import os
import getopt
import dataOutplts as dc
import time

import matplotlib.pyplot as plt
from matplotlib.ticker import OldScalarFormatter, ScalarFormatter
import matplotlib.ticker as plticker

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

    if ngases >= 2:
        fig0,  ax0     = plt.subplots(ngases/2,2,figsize=(11,15))
        fig1,  ax1   = plt.subplots(ngases,2,figsize=(10,15))
        fig2,  ax2   = plt.subplots(ngases,2,figsize=(10,15), sharex=True)
        fig3,  ax3   = plt.subplots(ngases,2,figsize=(10,15))

    else:
        fig0,  ax0   = plt.subplots(figsize=(6,7))
        fig1,  ax1   = plt.subplots(ngases, 2, figsize=(10, 7))
        fig2,  ax2   = plt.subplots(ngases, 2,figsize=(10,7))

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
        ds = gas.pltPrf(fltr=pltInputs['fltrFlg'],allGas=False,sclfct=pltInputs['sclfct'],sclname=pltInputs['sclfctName'],mnthFltr=pltInputs["mnths"],mnthFltFlg=pltInputs["mnthFlg"],
                errFlg=pltInputs['errorFlg'],minSZA=pltInputs['minSZA'],maxSZA=pltInputs['maxSZA'],maxRMS=maxRMS,minTC=pltInputs['minTC'],maxTC=pltInputs['maxTC'],
                minDOF=minDOF,maxCHI=pltInputs['maxCHI'],dofFlg=pltInputs['dofFlg'],rmsFlg=pltInputs['rmsFlg'],tcFlg=pltInputs['tcNegFlg'],
                pcFlg=pltInputs['pcNegFlg'],szaFlg=pltInputs['szaFlg'],cnvrgFlg=pltInputs['cnvrgFlg'],chiFlg=pltInputs['chiFlg'],tcMMflg=pltInputs['tcMMFlg'])

        #-------------------------------------------------------------------------
        #PLOTS
        #-------------------------------------------------------------------------

        print '\nPrinting Profile Plots: %s' %gas.PrimaryGas.lower(), '\n'
       
        #-------------------------------------------------------
        # Plot mean systematic and random errors on mean profiles in VMR
        #-------------------------------------------------------

        if ngases >= 2:

            if k <= ((ngases/2) - 1):

                ax0[k, 0].plot(ds['randvmr_mean'],ds['alt'], color='red', label = 'Random')
                #ax0[k, 0].errorbar(ds['randvmr_mean'],ds['alt'],fmt=None,xerr=ds['randvmr_std'],ecolor='r',label='Total Random Error')
                ax0[k, 0].fill_betweenx(ds['alt'],ds['randvmr_mean']-ds['randvmr_std'],ds['randvmr_mean']+ds['randvmr_std'],alpha=0.5,color='0.5', facecolor='red')
                #if k == 0: ax0[k, 0].set_title('Random Error')                    
                ax0[k, 0].set_ylabel('Altitude [km]', fontsize=14)      
                ax0[k, 0].grid(True,which='both')
                ax0[k, 0].annotate(gasname, xy=(.95, .89), xycoords='axes fraction', fontsize=24,
                    horizontalalignment='right', verticalalignment='center')

                ax0[k, 0].set_ylim(0, 50)
                idyrange = np.where(ds['alt'] < 50)[0]
                ax0[k, 0].set_xlim(0, max(ds['totvmr_mean'][idyrange] +ds['totvmr_std'][idyrange]) + 0.5*max(ds['totvmr_std'][idyrange]) )
                if k == (ngases/2)-1: ax0[k, 0].set_xlabel('VMR [ppb$_v$]', fontsize=14)
            
                ax0[k, 0].plot(ds['sysvmr_mean'],ds['alt'], color='blue', label='Systematic')
                ax0[k, 0].fill_betweenx(ds['alt'],ds['sysvmr_mean']-ds['sysvmr_std'],ds['sysvmr_mean']+ds['sysvmr_std'],alpha=0.5,color='0.5', facecolor='blue')
           

                ax0[k, 0].plot(ds['totvmr_mean'],ds['alt'], color='k', linestyle='--', linewidth=2, label = 'Total')

                if k == 0: ax0[k,0].legend(prop={'size':14}, loc='upper left')

            else:
                ax0[k-(ngases/2), 1].plot(ds['randvmr_mean'],ds['alt'], color='red', label = 'Random')
                #ax0[k, 0].errorbar(ds['randvmr_mean'],ds['alt'],fmt=None,xerr=ds['randvmr_std'],ecolor='r',label='Total Random Error')
                ax0[k-(ngases/2), 1].fill_betweenx(ds['alt'],ds['randvmr_mean']-ds['randvmr_std'],ds['randvmr_mean']+ds['randvmr_std'],alpha=0.5,color='0.5', facecolor='red')
                #if k == 0: ax0[k, 1].set_title('Random Error')                    
                #ax0[k-(ngases/2), 1].set_ylabel('Altitude [km]', fontsize=14)      
                ax0[k-(ngases/2), 1].grid(True,which='both')
                ax0[k-(ngases/2), 1].annotate(gasname, xy=(.95, .89), xycoords='axes fraction', fontsize=24,
                    horizontalalignment='right', verticalalignment='center')

                ax0[k-(ngases/2), 1].set_ylim(1, 50)
                idyrange = np.where(ds['alt'] < 50)[0]
                ax0[k-(ngases/2), 1].set_xlim(0, max(ds['totvmr_mean'][idyrange] +ds['totvmr_std'][idyrange]) + 0.5*max(ds['totvmr_std'][idyrange]) )
                if k == ngases-1: ax0[k-(ngases/2), 1].set_xlabel('VMR [ppb$_v$]', fontsize=14)
            
             
                ax0[k-(ngases/2), 1].plot(ds['sysvmr_mean'],ds['alt'], color='blue', label='Systematic')
                ax0[k-(ngases/2), 1].fill_betweenx(ds['alt'],ds['sysvmr_mean']-ds['sysvmr_std'],ds['sysvmr_mean']+ds['sysvmr_std'],alpha=0.5,color='0.5', facecolor='blue')
           
                ax0[k-(ngases/2), 1].plot(ds['totvmr_mean'],ds['alt'],label=gas.PrimaryGas+' Retrieved Monthly Mean', color='k', linestyle='--', linewidth=2)
            
            fig0.subplots_adjust(bottom=0.05, top=0.97, left = 0.07, right = 0.97)
        else:
            if gas.PrimaryGas == 'H2O': 
                mean_r = ds['randvmr_mean']/1000.0
                stdev_r = ds['randvmr_std']/1000.0
                mean_s = ds['sysvmr_mean']/1000.0
                stdev_s = ds['sysvmr_std']/1000.0
                tot     = ds['totvmr_mean']/1000.0
                stdev_t = ds['totvmr_std']/1000.0
            else:
                mean_r = ds['randvmr_mean']
                stdev_r = ds['randvmr_std']
                mean_s = ds['sysvmr_mean']
                stdev_s = ds['sysvmr_std']
                tot     = ds['totvmr_mean']
                stdev_t = ds['totvmr_std']

            #ax0.plot(ds['randvmr_mean'],ds['alt'], color='red', label = 'Random')
            ax0.plot(mean_r,ds['alt'], color='red', label = 'Random')
            ax0.fill_betweenx(ds['alt'],mean_r-stdev_r ,mean_r + stdev_r,alpha=0.5,color='0.5', facecolor='red')                    
            ax0.set_ylabel('Altitude [km]', fontsize=14)      
            ax0.grid(True,which='both')
            ax0.annotate(gasname, xy=(.95, .95), xycoords='axes fraction', fontsize=26, horizontalalignment='right', verticalalignment='center')
            ax0.set_ylim(0, 50)
            idyrange = np.where(ds['alt'] < 50)[0]
            ax0.set_xlim(0, max(tot[idyrange]) + max(tot[idyrange])*0.15  )# + stdev_t[idyrange]))
            if gas.PrimaryGas != 'H2O': ax0.set_xlabel('VMR [ppv$_v$]', fontsize=14)
            if gas.PrimaryGas == 'H2O': ax0.set_xlabel('VMR [ppm$_v$]', fontsize=14)
            
            ax0.plot(mean_s,ds['alt'], color='blue', label='Systematic')
            ax0.fill_betweenx(ds['alt'],mean_s-stdev_s,mean_s+stdev_s,alpha=0.5,color='0.5', facecolor='blue') 
            ax0.plot(tot,ds['alt'], color='k', linestyle='--', linewidth=2, label = 'Total')
            ax0.legend(prop={'size':14}, loc='upper left')

            fig0.subplots_adjust(bottom=0.08, top=0.95, left = 0.1, right = 0.97)

        #------------------------------------------------------------
        #------------------------------------------------------------
        # Plot individual components of error analysis in vmr
        #------------------------------------------------------------
        #------------------------------------------------------------
        list_comp_ran = list(sorted(ds['rand_cmpnts_vmr']))

        if ngases >= 2:
            maxerraxis = []
                    #-----------------------------
                    # Plot Random error components
                    #-----------------------------
            for p in list_comp_ran:

                if p == 'measurement_Random_Error': legend = 'Measurement'
                elif p == 'interfering_species_Random_Error': legend = 'Interfering species'
                elif p == 'temperature_Random_Error': legend = 'Temperature'
                elif p == 'sza_Random_Error': legend = 'SZA'
                elif p == 'retrieval_parameters_Random_Error': legend = 'Retrieval parameters'
                elif p == 'curvature_Random_Error': legend = 'Curvature'
                elif p == 'max_opd_Random_Error': legend = 'Max OPD'     
                else: legend = p 

                errPlt = ds['rand_cmpnts_vmr'][p]
                idyrange = np.where(ds['alt'] < 50)[0]    
                maxerraxis.append(max(errPlt[idyrange]))         
                           
                ax1[k,0].plot(errPlt,ds['alt'],linewidth=1.0, label=legend)            
                ax1[k,0].set_ylabel('Altitude [km]')
                if k == ngases-1: ax1[k,0].set_xlabel('VMR [ppb$_v$]', fontsize=14)             
                ax1[k,0].grid(True,which='both')
                if k == 0: ax1[k,0].legend(prop={'size':9},loc='upper right')          
                #ax2[k,0].tick_params(axis='both',which='both',labelsize=8) 
                if k == 0: ax1[k,0].set_title('Random Error Components')
                ax1[k, 0].set_ylim(0, 50)
                ax1[k, 0].set_xlim(0, max(ds['totvmr_mean'][idyrange] +ds['totvmr_std'][idyrange]) + 0.5*max(ds['totvmr_std'][idyrange]) )

                ax1[k, 0].annotate(gasname, xy=(0.075, 0.83), xycoords='axes fraction', fontsize=22,
                    horizontalalignment='left', verticalalignment='center')

                    #-----------------------------
                    # Plot systematic error components
                    #-----------------------------
            list_comp_sys = list(sorted(ds['sys_cmpnts_vmr']))
            
            for p in list_comp_sys:
                
               
                if p == 'temperature_Systematic_Error': legend = 'Temperature'
                elif p == 'smoothing_Systematic_Error': continue # legend = 'Smoothing'
                #elif p == 'smoothing_Systematic_Error': legend = 'Smoothing'            
                elif p == 'linepair_'+gas.PrimaryGas.lower()+'_Systematic_Error': legend = 'Line Pair'
                elif p == 'lineint_'+gas.PrimaryGas.lower()+'_Systematic_Error': legend = 'Line Int'
                elif p == 'linetair_'+gas.PrimaryGas.lower()+'_Systematic_Error': legend = 'Line Tair'
                elif p == 'phase_fcn_Systematic_Error': legend = 'Phase fcn'
                elif p == 'apod_fcn_Systematic_Error':  legend = 'Apod fcn'
                else: legend = p

                errPlt = ds['sys_cmpnts_vmr'][p]
                maxerraxis.append(max(errPlt[idyrange]))        
                ax1[k,1].plot(errPlt,ds['alt'],linewidth=1.0, label=legend)
        
                #------------------------
                # Plot total random error
                #------------------------
                #sysMean  = np.mean(ds['sys_err'],axis=0) / ds['prf_Mean']#  retPrf
                #ax2[k,1].plot(sysMean,ds['alt'],linewidth=0.75, label='Total Systematic Error')

                #ax2[k,1].set_ylabel('Altitude [km]')
                if k == ngases-1: ax1[k,1].set_xlabel('VMR [ppb$_v$]', fontsize=14)             
                ax1[k,1].grid(True,which='both')
                if k == 0: ax1[k,1].legend(prop={'size':9}, loc='upper right')
                    
                #ax2[k,1].tick_params(axis='both',which='both',labelsize=8) 
                if k == 0: ax1[k,1].set_title('Systematic Error Components')
                #ax2[k,1].locator_params(nbins=6 )
                ax1[k, 1].set_ylim(0, 50)
                ax1[k, 1].set_xlim(0, max(ds['totvmr_mean'][idyrange] +ds['totvmr_std'][idyrange]) + 0.5*max(ds['totvmr_std'][idyrange]) )


            ax1[k, 0].set_xlim(0, max(maxerraxis)+max(maxerraxis)*0.1)
            ax1[k, 1].set_xlim(0, max(maxerraxis)+max(maxerraxis)*0.1)
            del(maxerraxis)

            fig1.subplots_adjust(bottom=0.05, top=0.97, left = 0.07, right = 0.97)


        else:
            maxerraxis = []
                    #-----------------------------
                    # Plot Random error components
                    #-----------------------------
            for p in list_comp_ran:
            
                if p == 'measurement_Random_Error': legend = 'Measurement'
                elif p == 'interfering_species_Random_Error': continue# legend = 'Interfering species'
                elif p == 'temperature_Random_Error': legend = 'Temperature'
                elif p == 'sza_Random_Error': legend = 'SZA'
                elif p == 'curvature_Random_Error': legend = 'Curvature'
                elif p == 'max_opd_Random_Error': legend = 'Max OPD'
                elif p == 'retrieval_parameters_Random_Error': continue # legend = 'Retrieval parameters'
                else: legend = p

                if gas.PrimaryGas == 'H2O': 
                    errPlt = ds['rand_cmpnts_vmr'][p]/1000.0

                else:
                    errPlt = ds['rand_cmpnts_vmr'][p]
             
                idyrange = np.where(ds['alt'] < 50)[0]    
                maxerraxis.append(max(errPlt[idyrange]))        
                           
                ax1[0].plot(errPlt,ds['alt'],linewidth=1.0, label=legend)            
                ax1[0].set_ylabel('Altitude [km]', fontsize =14)
                if gas.PrimaryGas == 'H2O': 
                    ax1[0].set_xlabel('VMR [ppm$_v$]', fontsize=14)
                else:
                     ax1[0].set_xlabel('VMR [ppb$_v$]', fontsize=14)       
                ax1[0].grid(True,which='both')
                ax1[0].legend(prop={'size':12},loc='upper right')                 
                ax1[0].set_title('Random Error Components')
                ax1[0].set_ylim(0, 50)
                ax1[0].annotate(gasname, xy=(0.06, 0.95), xycoords='axes fraction', fontsize=26,
                    horizontalalignment='left', verticalalignment='center')
                if gas.PrimaryGas == 'H2O': 
                    ax1[0].set_xlim(0, max(ds['totvmr_mean'][idyrange]/1000) + max(ds['totvmr_mean'][idyrange]/1000)*0.15)#. +ds['totvmr_std'][idyrange]/1000.) )
                else:
                    ax1[0].set_xlim(0, max(ds['totvmr_mean'][idyrange] +ds['totvmr_std'][idyrange]) )

                    #-----------------------------
                    # Plot systematic error components
                    #-----------------------------
            list_comp_sys = list(sorted(ds['sys_cmpnts_vmr']))
            
            for p in list_comp_sys:        
                
                if p == 'temperature_Systematic_Error': legend = 'Temperature'
                elif p == 'smoothing_Systematic_Error': continue # legend = 'Smoothing'
                #elif p == 'smoothing_Systematic_Error': legend = 'Smoothing'
                elif p == 'linepair_'+gas.PrimaryGas.lower()+'_Systematic_Error': legend = 'Line Pair'
                elif p == 'lineint_'+gas.PrimaryGas.lower()+'_Systematic_Error': legend = 'Line Int'
                elif p == 'linetair_'+gas.PrimaryGas.lower()+'_Systematic_Error': legend = 'Line Tair'
                elif p == 'phase_fcn_Systematic_Error': legend = 'Phase fcn'
                elif p == 'apod_fcn_Systematic_Error':  legend = 'Apod fcn'
                else: legend = p

                if gas.PrimaryGas == 'H2O': 
                    errPlt = ds['sys_cmpnts_vmr'][p]/1000.0

                else:
                    errPlt = ds['sys_cmpnts_vmr'][p]

                
                maxerraxis.append(max(errPlt[idyrange]))         

                ax1[1].plot(errPlt,ds['alt'],linewidth=1.0, label=legend)
                if gas.PrimaryGas == 'H2O': 
                    ax1[1].set_xlabel('VMR [ppm$_v$]', fontsize=16)
                else:
                    ax1[1].set_xlabel('VMR [ppv$_v$]', fontsize=16)

                ax1[1].grid(True,which='both')
                ax1[1].legend(prop={'size':12}, loc='upper right')
                ax1[1].set_title('Systematic Error Components')
                ax1[1].set_ylim(0, 50)
                if gas.PrimaryGas == 'H2O': 
                    ax1[1].set_xlim(0, max(ds['totvmr_mean'][idyrange]/1000.) + max(ds['totvmr_mean'][idyrange]/1000.)*0.15)# +ds['totvmr_std'][idyrange]/1000.) )
                else:
                    ax1[1].set_xlim(0, max(ds['totvmr_mean'][idyrange] +ds['totvmr_std'][idyrange]) )

            #ax1[0].set_xlim(0, max(maxerraxis)+max(maxerraxis)*0.1)
            #ax1[1].set_xlim(0, max(maxerraxis)+max(maxerraxis)*0.1)
            del(maxerraxis)  
                
            fig1.subplots_adjust(bottom=0.08, top=0.95, left = 0.1, right = 0.97)

  
        #------------------------------------------------------------
        #------------------------------------------------------------
        #Plot individual components of error analysis fractional number
        #------------------------------------------------------------
        #------------------------------------------------------------

        if ngases >= 2:
            maxerraxis = [] 
        
            for p in list_comp_ran:

                if p == 'measurement_Random_Error': legend = 'Measurement'
                elif p == 'interfering_species_Random_Error':  legend = 'Interfering species'
                elif p == 'temperature_Random_Error': legend = 'Temperature'
                elif p == 'sza_Random_Error': legend = 'SZA'
                elif p == 'curvature_Random_Error': legend = 'Curvature'
                elif p == 'max_opd_Random_Error': legend = 'Max OPD'
                elif p == 'retrieval_parameters_Random_Error':legend = 'Retrieval parameters'
                else: legend = p

                errPlt = ds['rand_cmpnts_vmr'][p]
                idyrange = np.where(ds['alt'] < 50)[0]
                     
                #-------------------------------------------------
                # Normalize error as fraction of retrieved profile
                #-------------------------------------------------
                errPlt = errPlt /  ds['prfvmr_mean']
                maxerraxis.append(max(errPlt[idyrange]))          
                           
                ax2[k,0].plot(errPlt,ds['alt'],linewidth=1.0, label=legend)
                               

                ax2[k,0].set_ylabel('Altitude [km]')
                if k == ngases-1: ax2[k,0].set_xlabel('Fraction of Retrieved Profile', fontsize=14)             
                ax2[k,0].grid(True,which='both')
                if k == 0: ax2[k,0].legend(prop={'size':9},loc='upper right')          
                #ax2[k,0].tick_params(axis='both',which='both',labelsize=8) 
                if k == 0: ax2[k,0].set_title('Random Error Components')
                ax2[k, 0].set_ylim(0, 50)

                ax2[k, 0].annotate(gasname, xy=(0.075, 0.83), xycoords='axes fraction', fontsize=22,
                    horizontalalignment='left', verticalalignment='center')

                    #-----------------------------
                    # Plot systematic error components
                    #-----------------------------
            
            for p in list_comp_sys:      
                
                if p == 'temperature_Systematic_Error': legend = 'Temperature'
                elif p == 'smoothing_Systematic_Error': continue # legend = 'Smoothing'
                #elif p == 'smoothing_Systematic_Error': legend = 'Smoothing'
                elif p == 'linepair_'+gas.PrimaryGas.lower()+'_Systematic_Error': legend = 'Line Pair'
                elif p == 'lineint_'+gas.PrimaryGas.lower()+'_Systematic_Error': legend = 'Line Int'
                elif p == 'linetair_'+gas.PrimaryGas.lower()+'_Systematic_Error': legend = 'Line Tair'
                elif p == 'phase_fcn_Systematic_Error': legend = 'Phase fcn'
                elif p == 'apod_fcn_Systematic_Error':  legend = 'Apod fcn'
                else: legend = p         

                errPlt = ds['sys_cmpnts_vmr'][p]
                 
               
                #-------------------------------------------------
                # Normalize error as fraction of retrieved profile
                #-------------------------------------------------
                errPlt = errPlt / ds['prfvmr_mean']#  retPrf 
                maxerraxis.append(max(errPlt[idyrange]))   

                ax2[k,1].plot(errPlt,ds['alt'],linewidth=1.0, label=legend)
        
                #------------------------
                # Plot total random error
                #------------------------
                #sysMean  = np.mean(ds['sys_err'],axis=0) / ds['prf_Mean']#  retPrf
                #ax2[k,1].plot(sysMean,ds['alt'],linewidth=0.75, label='Total Systematic Error')

                #ax2[k,1].set_ylabel('Altitude [km]')
                if k == ngases-1: ax2[k,1].set_xlabel('Fraction of Retrieved Profile', fontsize=14)             
                ax2[k,1].grid(True,which='both')
                if k == 0: ax2[k,1].legend(prop={'size':9}, loc='upper right')
                    
                #ax2[k,1].tick_params(axis='both',which='both',labelsize=8) 
                if k == 0: ax2[k,1].set_title('Systematic Error Components')
                #ax2[k,1].locator_params(nbins=6 )
                ax2[k, 1].set_ylim(0, 50)

            if pltInputs['loc'].lower() == 'tab':
                ax2[k, 0].set_xlim(0, 0.5)
                ax2[k, 1].set_xlim(0, 0.5)
            else:
                ax2[k, 0].set_xlim(0, 0.3)
                ax2[k, 1].set_xlim(0, 0.3)


            fig2.subplots_adjust(bottom=0.05, top=0.97, left = 0.07, right = 0.97)

            
        else:
            maxerraxis = []

            for p in list_comp_ran:

                if p == 'measurement_Random_Error': legend = 'Measurement'
                elif p == 'interfering_species_Random_Error': continue # legend = 'Interfering species'
                elif p == 'temperature_Random_Error': legend = 'Temperature'
                elif p == 'sza_Random_Error': legend = 'SZA'
                elif p == 'curvature_Random_Error': legend = 'Curvature'
                elif p == 'max_opd_Random_Error': legend = 'Max OPD'
                elif p == 'retrieval_parameters_Random_Error': continue # legend = 'Retrieval parameters'
                else: legend = p     

                errPlt = ds['rand_cmpnts_vmr'][p]
                idyrange = np.where(ds['alt'] < 50)[0]
                #-------------------------------------------------
                # Normalize error as fraction of retrieved profile
                #-------------------------------------------------
                errPlt = errPlt /  ds['prfvmr_mean']
                maxerraxis.append(max(errPlt[idyrange]))       
                           
                ax2[0].plot(errPlt,ds['alt'],linewidth=1.0, label=legend)  
                ax2[0].set_ylabel('Altitude [km]', fontsize=14)
                ax2[0].set_xlabel('Fraction of Retrieved Profile', fontsize=14)             
                ax2[0].grid(True,which='both')
                ax2[0].legend(prop={'size':12},loc='upper right')          
                #ax2[k,0].tick_params(axis='both',which='both',labelsize=8) 
                ax2[0].set_title('Random Error Components')
                ax2[0].set_ylim(0, 50)

                ax2[0].annotate(gasname, xy=(0.075, 0.95), xycoords='axes fraction', fontsize=26,
                    horizontalalignment='left', verticalalignment='center')

                    #-----------------------------
                    # Plot systematic error components
                    #-----------------------------
            #list_comp_sys = list(sorted(ds['sys_cmpnts_vmr']))
            
            for p in list_comp_sys:
                
                if p == 'temperature_Systematic_Error': legend = 'Temperature'
                elif p == 'smoothing_Systematic_Error': continue #legend = 'Smoothing'
                #elif p == 'smoothing_Systematic_Error': legend = 'Smoothing'
                elif p == 'linepair_'+gas.PrimaryGas.lower()+'_Systematic_Error': legend = 'Line Pair'
                elif p == 'lineint_'+gas.PrimaryGas.lower()+'_Systematic_Error': legend = 'Line Int'
                elif p == 'linetair_'+gas.PrimaryGas.lower()+'_Systematic_Error': legend = 'Line Tair'
                elif p == 'phase_fcn_Systematic_Error': legend = 'Phase fcn'
                elif p == 'apod_fcn_Systematic_Error':  legend = 'Apod fcn'
                else: legend = p

                errPlt = ds['sys_cmpnts_vmr'][p]                   
                #-------------------------------------------------
                # Normalize error as fraction of retrieved profile
                #-------------------------------------------------
                errPlt = errPlt / ds['prfvmr_mean']#  retPrf 
                maxerraxis.append(max(errPlt[idyrange]))        
                   
                ax2[1].plot(errPlt,ds['alt'],linewidth=1.0, label=legend)
        
                #------------------------
                # Plot total random error
                #------------------------
                #sysMean  = np.mean(ds['sys_err'],axis=0) / ds['prf_Mean']#  retPrf
                #ax2[k,1].plot(sysMean,ds['alt'],linewidth=0.75, label='Total Systematic Error')

                #ax2[k,1].set_ylabel('Altitude [km]')
                ax2[1].set_xlabel('Fraction of Retrieved Profile', fontsize=14)             
                ax2[1].grid(True,which='both')
                ax2[1].legend(prop={'size':12}, loc='upper right')
                    
                #ax2[k,1].tick_params(axis='both',which='both',labelsize=8) 
                ax2[1].set_title('Systematic Error Components')
                #ax2[k,1].locator_params(nbins=6 )
                ax2[1].set_ylim(0, 50)
                
                fig2.subplots_adjust(bottom=0.08, top=0.95, left = 0.1, right = 0.97)

            ax2[0].set_xlim(0, max(maxerraxis) + 0.1*max(maxerraxis) )
            ax2[1].set_xlim(0, max(maxerraxis) + 0.1*max(maxerraxis) )


        if ngases >= 2:
            maxerraxis = [] 
        
            for p in list_comp_ran:

                if p == 'measurement_Random_Error': legend = 'Measurement'
                elif p == 'interfering_species_Random_Error': legend = 'Interfering species'
                elif p == 'temperature_Random_Error': legend = 'Temperature'
                elif p == 'sza_Random_Error': legend = 'SZA'
                elif p == 'curvature_Random_Error': legend = 'Curvature'
                elif p == 'max_opd_Random_Error': legend = 'Max OPD'
                elif p == 'retrieval_parameters_Random_Error': legend = 'Retrieval parameters'
                else: legend = p

                errPlt = ds['rand_cmpnts_vmr'][p]
                idyrange = np.where(ds['alt'] < 50)[0]
                     
                #-------------------------------------------------
                # Normalize error as fraction of retrieved profile
                #-------------------------------------------------
                errPlt = errPlt /  ds['prfvmr_mean']
                maxerraxis.append(max(errPlt[idyrange]))          
                           
                ax3[k,0].plot(errPlt,ds['alt'],linewidth=1.0, label=legend)
                               

                ax3[k,0].set_ylabel('Altitude [km]')
                if k == ngases-1: ax3[k,0].set_xlabel('Fraction of Retrieved Profile', fontsize=14)             
                ax3[k,0].grid(True,which='both')
                if k == 0: ax3[k,0].legend(prop={'size':9},loc='upper right')          
                #ax2[k,0].tick_params(axis='both',which='both',labelsize=8) 
                if k == 0: ax3[k,0].set_title('Random Error Components')
                ax3[k, 0].set_ylim(0, 50)

                ax3[k, 0].annotate(gasname, xy=(0.075, 0.83), xycoords='axes fraction', fontsize=22,
                    horizontalalignment='left', verticalalignment='center')

                    #-----------------------------
                    # Plot systematic error components
                    #-----------------------------
            
            for p in list_comp_sys:      
                
                if p == 'temperature_Systematic_Error': legend = 'Temperature'
                elif p == 'smoothing_Systematic_Error': continue # legend = 'Smoothing'
                #elif p == 'smoothing_Systematic_Error': legend = 'Smoothing'
                elif p == 'linepair_'+gas.PrimaryGas.lower()+'_Systematic_Error': legend = 'Line Pair'
                elif p == 'lineint_'+gas.PrimaryGas.lower()+'_Systematic_Error': legend = 'Line Int'
                elif p == 'linetair_'+gas.PrimaryGas.lower()+'_Systematic_Error': legend = 'Line Tair'
                elif p == 'phase_fcn_Systematic_Error': legend = 'Phase fcn'
                elif p == 'apod_fcn_Systematic_Error':  legend = 'Apod fcn'
                else: legend = p         

                errPlt = ds['sys_cmpnts_vmr'][p]
                 
               
                #-------------------------------------------------
                # Normalize error as fraction of retrieved profile
                #-------------------------------------------------
                errPlt = errPlt / ds['prfvmr_mean']#  retPrf 
                maxerraxis.append(max(errPlt[idyrange]))   

                ax3[k,1].plot(errPlt,ds['alt'],linewidth=1.0, label=legend)
        
                #------------------------
                # Plot total random error
               
                if k == ngases-1: ax3[k,1].set_xlabel('Fraction of Retrieved Profile', fontsize=14)             
                ax3[k,1].grid(True,which='both')
                if k == 0: ax3[k,1].legend(prop={'size':9}, loc='upper right')

                if k == 0: ax3[k,1].set_title('Systematic Error Components')
           
                ax3[k, 1].set_ylim(0, 50)


            ax3[k, 0].set_xlim(0, max(maxerraxis)+max(maxerraxis)*0.1)
            ax3[k, 1].set_xlim(0, max(maxerraxis)+max(maxerraxis)*0.1)
            del(maxerraxis)


            fig3.subplots_adjust(bottom=0.05, top=0.97, left = 0.07, right = 0.97)


#         #-------------------------------------------------------
#         # Plot mean systematic and random errors on mean profiles in molec/cm2
#         #-------------------------------------------------------
  
#         ax[k, 0].plot(ds['rand_mean'],ds['alt'],color='k',label=gas.PrimaryGas+' Retrieved Monthly Mean',linewidth=2, color='red')
#         #ax[k, 0].errorbar(ds['rand_mean'],ds['alt'],fmt=None,xerr=ds['rand_std'],ecolor='r',label='Total Random Error')
#         ax[k, 0].fill_betweenx(ds['alt'],ds['rand_mean']-ds['rand_std'],ds['rand_mean']+ds['rand_std'],alpha=0.75,color='0.75', facecolor='red')  
#         if k == 0: ax[k, 0].set_title('Random Error')                    
#         ax[k, 0].set_ylabel('Altitude [km]')      
#         ax[k, 0].grid(True,which='both')
#         #ax[k, 0].text(0.1, 0.8, gasname, fontsize=24, textcoords='axes fraction')
#         #ax[k, 0].annotate(gasname, 0.1, 0.8, fontsize=24, textcoords='axes fraction')
#         ax[k, 0].annotate(gasname, xy=(.05, .89), xycoords='axes fraction', fontsize=24,
#                 horizontalalignment='left', verticalalignment='center')
      
#         ax[k, 0].set_ylim(0, 80)
#         ax[k, 0].set_xlim(0, max(ds['tot_mean'] +ds['tot_std']) + 0.5*max(ds['tot_std'])   )

#         #ax[k, 1]=plt.subplot(temp)           
#         ax[k, 1].plot(ds['sys_mean'],ds['alt'],color='k',label=gas.PrimaryGas+' Retrieved Monthly Mean',linewidth=2)
#         #ax[k, 1].errorbar(ds['sys_mean'],ds['alt'],fmt=None,xerr=ds['sys_std'],ecolor='r',label='Total Systematic Error')
#         ax[k, 1].fill_betweenx(ds['alt'],ds['sys_mean']-ds['sys_std'],ds['sys_mean']+ds['sys_std'],alpha=0.75,color='0.75')     
#         if k == 0: ax[k, 1].set_title('Systematic Error')      
#         if k == ngases-1: ax[k, 1].set_xlabel('molecules cm$^{-2}$')                                         
#         ax[k, 1].grid(True,which='both')
#         ax[k, 1].set_ylim(0, 80)
#         ax[k, 1].set_xlim(0, max(ds['tot_mean'] +ds['tot_std']) + 0.5*max(ds['tot_std'])   ) 
                                   

#         #ax[k, 2]=plt.subplot(temp)                                 
#         ax[k, 2].plot(ds['tot_mean'],ds['alt'],color='k',label=gas.PrimaryGas+' Retrieved Monthly Mean',linewidth=2)
#         #ax[k, 2].errorbar(ds['tot_mean'],ds['alt'],fmt=None,xerr=ds['tot_std'],ecolor='r',label='Total Error')
#         ax[k, 2].fill_betweenx(ds['alt'],ds['tot_mean']-ds['tot_std'],ds['tot_mean']+ds['tot_std'],alpha=0.75,color='0.75')                
#         if k == 0: ax[k, 2].set_title('Total Error')  
#         ax[k, 2].grid(True,which='both')
#         ax[k, 2].set_ylim(0, 80)
#         ax[k, 2].set_xlim(0, max(ds['tot_mean'] +ds['tot_std']) + 0.5*max(ds['tot_std'])   )
       

#         fig.subplots_adjust(bottom=0.07, top=0.95, left = 0.07, right = 0.97)


#         #------------------------------------------------------------
#         # Plot individual components of error analysis in molec/cm2
#         #------------------------------------------------------------

#         maxerraxis = []
#         for p in ds['rand_cmpnts']:
#             if len(gas.dirLst) > 1:
#                 errPlt = np.mean(np.sqrt(ds['rand_cmpnts'][p]),axis=0)
#                 retPrf = np.mean(ds['rPrfMol'],axis=0)

#             else:
#                 errPlt = ds['rand_cmpnts'][p][0]

#             maxerraxis.append(max(errPlt)) 
                    
#             #-------------------------------------------------
#             # Normalize error as fraction of retrieved profile
#             #-------------------------------------------------
#             #errPlt = errPlt /  ds['prf_Mean']#  retPrf

#             if p == 'measurement_Random_Error': legend = 'Measurement'
#             elif p == 'interfering_species_Random_Error': legend = 'Interfering species'
#             elif p == 'temperature_Random_Error': legend = 'Temperature'
#             elif p == 'sza_Random_Error': legend = 'SZA'
#             elif p == 'retrieval_parameters_Random_Error': legend = 'Retrieval parameters'
#             else: legend = p         
                       
#             ax2[k,0].plot(errPlt,ds['alt'],linewidth=1.0, label=legend)
                    
#             #------------------------
#             # Plot total random error
#             #------------------------
#             #randMean = np.mean(ds['rand_err'],axis=0) / ds['prf_Mean']#  retPrf
#             #ax2[k,0].plot(randMean,ds['alt'],linewidth=0.75, label='Total Random Error')                

#             ax2[k,0].set_ylabel('Altitude [km]')
#             if k == ngases-1: ax2[k,0].set_xlabel('molecules cm$^{-2}$')             
#             ax2[k,0].grid(True,which='both')
#             if k == 0: ax2[k,0].legend(prop={'size':9})          
#             #ax2[k,0].tick_params(axis='both',which='both',labelsize=8) 
#             if k == 0: ax2[k,0].set_title('Random Error Components')
#             ax2[k, 0].set_ylim(0, 80)
            

    
#                 #-----------------------------
#                 # Plot random error components
#                 #-----------------------------
#         for p in ds['sys_cmpnts']:
#             if len(gas.dirLst) > 1:
#                 errPlt = np.mean(np.sqrt(ds['sys_cmpnts'][p]),axis=0)
#             else:
#                 errPlt = ds['sys_cmpnts'][p][0]

#             maxerraxis.append(max(errPlt)) 
            
#              #-------------------------------------------------
#             # Normalize error as fraction of retrieved profile
#             #-------------------------------------------------
#             #errPlt = errPlt / ds['prf_Mean']#  retPrf         
            
#             if p == 'temperature_Systematic_Error': legend = 'Temperature'
#             elif p == 'smoothing_Systematic_Error': legend = 'Smoothing'
#             elif p == 'linepair_'+gas.PrimaryGas.lower()+'_Systematic_Error': legend = 'Linepair'
#             elif p == 'lineint_'+gas.PrimaryGas.lower()+'_Systematic_Error': legend = 'Lineint'
#             elif p == 'linetair_'+gas.PrimaryGas.lower()+'_Systematic_Error': legend = 'Measurement'
#             else: legend = p         


#             ax2[k,1].plot(errPlt,ds['alt'],linewidth=1.0, label=legend)
    
#             #------------------------
#             # Plot total random error
#             #------------------------
#             #sysMean  = np.mean(ds['sys_err'],axis=0) / ds['prf_Mean']#  retPrf
#             #ax2[k,1].plot(sysMean,ds['alt'],linewidth=0.75, label='Total Systematic Error')

#             #ax2[k,1].set_ylabel('Altitude [km]')
#             if k == ngases-1: ax2[k,1].set_xlabel('molecules cm$^{-2}$')             
#             ax2[k,1].grid(True,which='both')
#             if k == 0: ax2[k,1].legend(prop={'size':9})
                
#             #ax2[k,1].tick_params(axis='both',which='both',labelsize=8) 
#             if k == 0: ax2[k,1].set_title('Systematic Error Components')
#             #ax2[k,1].locator_params(nbins=6 )
#             ax2[k, 1].set_ylim(0, 80)
#             ax2[k, 0].annotate(gasname, xy=(.05, .89), xycoords='axes fraction', fontsize=24,
#                 horizontalalignment='left', verticalalignment='center')

#             fig2.subplots_adjust(bottom=0.07, top=0.95)

#         ax2[k, 0].set_xlim(0, max(maxerraxis)+max(maxerraxis)*0.1)
#         ax2[k, 1].set_xlim(0, max(maxerraxis)+max(maxerraxis)*0.1)
#         del(maxerraxis) 

#         #------------------------------------------------------------
#         # Plot individual components of error in fractional number
#         #------------------------------------------------------------

#         for p in ds['rand_cmpnts_vmr']:
#             if len(gas.dirLst) > 1:
#                 errPlt = np.mean(ds['rand_cmpnts_vmr'][p],axis=0)

#             else:
#                 errPlt = ds['rand_cmpnts_vmr'][p][0]
                    
#             #-------------------------------------------------
#             # Normalize error as fraction of retrieved profile
#             #-------------------------------------------------
#             errPlt = errPlt /  ds['prfvmr_mean']#  retPrf

#             if p == 'measurement_Random_Error': legend = 'Measurement'
#             elif p == 'interfering_species_Random_Error': legend = 'Interfering species'
#             elif p == 'temperature_Random_Error': legend = 'Temperature'
#             elif p == 'sza_Random_Error': legend = 'SZA'
#             elif p == 'retrieval_parameters_Random_Error': legend = 'Retrieval parameters'
#             else: legend = p         
                       
#             ax23[k,0].plot(errPlt,ds['alt'],linewidth=1.0, label=legend)
                    
#             #------------------------
#             # Plot total random error
#             #------------------------
#             #randMean = np.mean(ds['rand_err'],axis=0) / ds['prf_Mean']#  retPrf
#             #ax2[k,0].plot(randMean,ds['alt'],linewidth=0.75, label='Total Random Error')                

#             ax23[k,0].set_ylabel('Altitude [km]')
#             if k == ngases-1: ax23[k,0].set_xlabel('Fraction of Retrieved Profile')             
#             ax23[k,0].grid(True,which='both')
#             if k == 0: ax23[k,0].legend(prop={'size':9})          
#             #ax2[k,0].tick_params(axis='both',which='both',labelsize=8) 
#             if k == 0: ax23[k,0].set_title('Random Error Components')
#             ax23[k, 0].set_ylim(0, 80)
#             ax23[k, 0].set_xlim(0, 1)  
    
#                 #-----------------------------
#                 # Plot random error components
#                 #-----------------------------
#         for p in ds['sys_cmpnts_vmr']:
#             if len(gas.dirLst) > 1:
#                 errPlt = np.mean(ds['sys_cmpnts_vmr'][p],axis=0)
#             else:
#                 errPlt = ds['sys_cmpnts_vmr'][p][0]
            
#             #-------------------------------------------------
#             # Normalize error as fraction of retrieved profile
#             #-------------------------------------------------
#             errPlt = errPlt / ds['prfvmr_mean']#  retPrf         
            
#             if p == 'temperature_Systematic_Error': legend = 'Temperature'
#             elif p == 'smoothing_Systematic_Error': legend = 'Smoothing'
#             elif p == 'linepair_'+gas.PrimaryGas.lower()+'_Systematic_Error': legend = 'Linepair'
#             elif p == 'lineint_'+gas.PrimaryGas.lower()+'_Systematic_Error': legend = 'Lineint'
#             elif p == 'linetair_'+gas.PrimaryGas.lower()+'_Systematic_Error': legend = 'Measurement'
#             else: legend = p         


#             ax23[k,1].plot(errPlt,ds['alt'],linewidth=1.0, label=legend)
    
#             #------------------------
#             # Plot total random error
#             #------------------------
#             #sysMean  = np.mean(ds['sys_err'],axis=0) / ds['prf_Mean']#  retPrf
#             #ax2[k,1].plot(sysMean,ds['alt'],linewidth=0.75, label='Total Systematic Error')

#             #ax2[k,1].set_ylabel('Altitude [km]')
#             if k == ngases-1: ax23[k,1].set_xlabel('Fraction of Retrieved Profile')             
#             ax23[k,1].grid(True,which='both')
#             if k == 0: ax23[k,1].legend(prop={'size':9})
                
#             #ax2[k,1].tick_params(axis='both',which='both',labelsize=8) 
#             if k == 0: ax23[k,1].set_title('Systematic Error Components')
#             #ax2[k,1].locator_params(nbins=6 )
#             ax23[k, 1].set_ylim(0, 80)
#             ax23[k, 0].annotate(gasname, xy=(.05, .89), xycoords='axes fraction', fontsize=24,
#                 horizontalalignment='left', verticalalignment='center')
#             ax23[k, 1].set_xlim(0, 1) 

#             fig23.subplots_adjust(bottom=0.07, top=0.95) 


#                 #-----------------------------
#                 # Additional Plots in molec/cm2
#                 #-----------------------------

#         ax3[k].plot(ds['prf_Mean'],ds['alt'],color='k')
#         ax3[k].errorbar(ds['prf_Mean'],ds['alt'],fmt=None,xerr=ds['tot_mean'],ecolor='r')
#         #ax[k, 0].fill_betweenx(ds['alt'],ds['rand_mean']-ds['rand_max'],ds['rand_mean']+ds['rand_max'],alpha=0.5,color='0.75')  
#         if k == 0: ax3[k].set_title('Mean Profile $\pm$ total error')                    
#         ax3[k].set_ylabel('Altitude [km]')      
#         ax3[k].grid(True,which='both')
#         ax3[k].annotate(gasname, xy=(.05, .89), xycoords='axes fraction', fontsize=24,
#                 horizontalalignment='left', verticalalignment='center')  
#         ax3[k].set_ylim(0, 80)
#         if k == ngases-1: ax3[k].set_xlabel('molecules cm$^{-2}$')

#         fig3.subplots_adjust(bottom=0.07, top=0.95)  

# # #-------------------------------------------------------------------------------------------------------------------------------------------------------------------
#         #-------------------------------------------------------
#         # Plot mean systematic and random errors on mean profiles in VMR
#         #-------------------------------------------------------
#         #loc = plticker.MultipleLocator(base= np.min(ds['randvmr_mean'])) # this locator puts ticks at regular intervals

#         #if np.mean(ds['randvmr_mean']) < 1.0: 
#              #ax4[k, 0].plot(ds['randvmr_mean']*10.0,ds['alt'],color='k',label=gas.PrimaryGas+' Retrieved Monthly Mean')
#         #else: 
#         ax4[k, 0].plot(ds['randvmr_mean'],ds['alt'],color='k',label=gas.PrimaryGas+' Retrieved Monthly Mean')
#         ax4[k, 0].errorbar(ds['randvmr_mean'],ds['alt'],fmt=None,xerr=ds['randvmr_std'],ecolor='r',label='Total Random Error')
#         #ax4[k, 0].fill_betweenx(ds['alt'],ds['rand_mean']-ds['rand_max'],ds['rand_mean']+ds['rand_max'],alpha=0.5,color='0.75')  
#         if k == 0: ax4[k, 0].set_title('Random Error')                    
#         ax4[k, 0].set_ylabel('Altitude [km]')      
#         ax4[k, 0].grid(True,which='both')
#         ax4[k, 0].annotate(gasname, xy=(.05, .89), xycoords='axes fraction', fontsize=24,
#                 horizontalalignment='left', verticalalignment='center')  
#         ax4[k, 0].set_ylim(0, 80)
#         #ax4[k, 0].xaxis.set_major_locator(loc)  

#         #ax[k, 1]=plt.subplot(temp)           
#         ax4[k, 1].plot(ds['sysvmr_mean'],ds['alt'],color='k',label=gas.PrimaryGas+' Retrieved Monthly Mean')
#         ax4[k, 1].errorbar(ds['sysvmr_mean'],ds['alt'],fmt=None,xerr=ds['sysvmr_std'],ecolor='r',label='Total Systematic Error')
#        # ax[k, 1].fill_betweenx(ds['alt'],ds['sys_mean']-ds['sys_max'],ds['sys_mean']+ds['sys_max'],alpha=0.5,color='0.75')      
#         if k == 0: ax4[k, 1].set_title('Systematic Error')      
#         if k == ngases-1: ax4[k, 1].set_xlabel('VMR [ppb$_v$]')                                         
#         ax4[k, 1].grid(True,which='both')
#         ax4[k, 1].set_ylim(0, 80) 
#         #ax4[k, 1].xaxis.set_major_locator(loc)   
                                   

#         #ax[k, 2]=plt.subplot(temp)                                 
#         ax4[k, 2].plot(ds['totvmr_mean'],ds['alt'],color='k',label=gas.PrimaryGas+' Retrieved Monthly Mean')
#         ax4[k, 2].errorbar(ds['totvmr_mean'],ds['alt'],fmt=None,xerr=ds['totvmr_std'],ecolor='r',label='Total Error')
#         #ax[k, 2].fill_betweenx(ds['alt'],ds['tot_mean']-ds['tot_max'],ds['tot_mean']+ds['tot_max'],alpha=0.5,color='0.75')                
#         if k == 0: ax4[k, 2].set_title('Total Error')  
#         ax4[k, 2].grid(True,which='both')
#         ax4[k, 2].set_ylim(0, 80)
#         #ax4[k, 2].xaxis.set_major_formatter(ScalarFormatter(useMathText=True), interval=2)
#         #ax4[k, 2].xaxis.set_major_locator(loc) 

#         fig4.subplots_adjust(bottom=0.07, top=0.95, left = 0.07, right = 0.97) 

#                  #-----------------------------
#                  # Additional Plots in vmr
#                  #-----------------------------

#         ax5[k].plot(ds['prfvmr_mean'],ds['alt'],color='k')
#         ax5[k].errorbar(ds['prfvmr_mean'],ds['alt'],fmt=None,xerr=ds['totvmr_mean'],ecolor='r')
#         #ax[k, 0].fill_betweenx(ds['alt'],ds['rand_mean']-ds['rand_max'],ds['rand_mean']+ds['rand_max'],alpha=0.5,color='0.75')  
#         if k == 0: ax5[k].set_title('Mean Profile $\pm$ total error')                    
#         ax5[k].set_ylabel('Altitude [km]')      
#         ax5[k].grid(True,which='both') 
#         ax5[k].set_ylim(0, 80)
#         if k == ngases-1: ax5[k].set_xlabel('VMR [ppb$_v$]')
#         ax5[k].annotate(gasname, xy=(.05, .89), xycoords='axes fraction', fontsize=24,
#                 horizontalalignment='left', verticalalignment='center')
#         fig5.subplots_adjust(bottom=0.07, top=0.95)   

                    
    if gas.pdfsav:
        gas.pdfsav.savefig(fig0,dpi=200) 
        gas.pdfsav.savefig(fig1,dpi=200)
        gas.pdfsav.savefig(fig2,dpi=200)
        if ngases >= 2: gas.pdfsav.savefig(fig3,dpi=200)
        # gas.pdfsav.savefig(fig4,dpi=200)
        #gas.pdfsav.savefig(fig22,dpi=200)
        # gas.pdfsav.savefig(fig5,dpi=200)
        # gas.pdfsav.savefig(fig23,dpi=200)
    else:           
        plt.show(block=False)          
                    
    #-------------------------------------------------------------------------
    #Plot individual components of error analysis
    #-------------------------------------------------------------------------                
                    

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
