#!/usr/bin/python
##! /usr/local/python-2.7/bin/python
#----------------------------------------------------------------------------------------
#
# Purpose:
#
# Plot average Averaging kernel of a serie of gases defined in input_pltaks.py
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
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
import matplotlib.cm as mplcm


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
        fig, ax = plt.subplots(ngases/2,2,figsize=(11,15), sharey=True)
        fig1, ax1 = plt.subplots(ngases/2,2,figsize=(11,15), sharey=True)

        if pltInputs['pCols']:
            fig2, ax2 = plt.subplots(ngases/2,2,figsize=(11,15), sharey=True)
            fig22, ax22 = plt.subplots(ngases/2,2,figsize=(11,15), sharey=True)
    
    else:
        fig, ax = plt.subplots(figsize=(6,7))
       

        if pltInputs['pCols']:
            fig2, ax2 = plt.subplots(figsize=(6,7))


    clr = ('red', 'green', 'blue', 'cyan', 'yellow', 'orange')

   
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
        else: gasname = gas.PrimaryGas

        #gasname = gas.PrimaryGas

        # ds = gas.pltPrf(fltr=pltInputs['fltrFlg'],allGas=False,sclfct=pltInputs['sclfct'],sclname=pltInputs['sclfctName'],mnthFltr=pltInputs["mnths"],mnthFltFlg=pltInputs["mnthFlg"],
        #         errFlg=pltInputs['errorFlg'],minSZA=pltInputs['minSZA'],maxSZA=pltInputs['maxSZA'],maxRMS=maxRMS,minTC=pltInputs['minTC'],maxTC=pltInputs['maxTC'],
        #         minDOF=minDOF,maxCHI=pltInputs['maxCHI'],dofFlg=pltInputs['dofFlg'],rmsFlg=pltInputs['rmsFlg'],tcFlg=pltInputs['tcNegFlg'],
        #         pcFlg=pltInputs['pcNegFlg'],szaFlg=pltInputs['szaFlg'],cnvrgFlg=pltInputs['cnvrgFlg'],chiFlg=pltInputs['chiFlg'],tcMMflg=pltInputs['tcMMFlg'])
        # exit()

        ds = gas.pltAvk(fltr=pltInputs['fltrFlg'],errFlg=pltInputs['errorFlg'],partialCols=pltInputs['pCols'],mnthFltr=pltInputs["mnths"],mnthFltFlg=pltInputs["mnthFlg"],
                   minSZA=pltInputs['minSZA'],maxSZA=pltInputs['maxSZA'],maxRMS=maxRMS,minTC=pltInputs['minTC'],maxTC=pltInputs['maxTC'],
                   minDOF=minDOF,maxCHI=pltInputs['maxCHI'],dofFlg=pltInputs['dofFlg'],rmsFlg=pltInputs['rmsFlg'],
                   tcFlg=pltInputs['tcNegFlg'],pcFlg=pltInputs['pcNegFlg'],szaFlg=pltInputs['szaFlg'],
                   chiFlg=pltInputs['chiFlg'],cnvrgFlg=pltInputs['cnvrgFlg'],tcMMflg=pltInputs['tcMMFlg'],sclfct=pltInputs['sclfct'],sclname=pltInputs['sclfctName'])

        #-------------------------------------------------
        # Averaging Kernel Smoothing Function (row of avk)
        #-------------------------------------------------
        clmap     = 'jet'

        gs        = gridspec.GridSpec(1,2,width_ratios=[3,1])
        cm        = plt.get_cmap(clmap)
        cNorm     = colors.Normalize(vmin=np.min(ds['alt']), vmax=np.max(ds['alt']))
        scalarMap = mplcm.ScalarMappable(norm=cNorm,cmap=clmap)
        scalarMap.set_array(ds['alt'])

        #-----------------------------PLOT1---------------------------
        if ngases >= 2:

            if k <= ((ngases/2) - 1):
                for i in range(len(ds['alt'])):
                    im = ax[k, 0].plot(ds['avkSCF'][i,:],ds['alt'])

                ax[k, 0].set_color_cycle([scalarMap.to_rgba(x) for x in ds['alt']])
                ax[k, 0].set_ylabel('Altitude [km]',  fontsize=18)
                if k == (ngases/2)-1: ax[k, 0].set_xlabel('Averaging Kernels', fontsize=18)
                ax[k, 0].grid(True)        
                #ax[k, 0].set_title('Averaging Kernels Scale Factor')
                ax[k, 0].set_ylim(0, 50)
                ax[k, 0].annotate(gasname, xy=(.95, .89), xycoords='axes fraction', fontsize=30,
                    horizontalalignment='right', verticalalignment='center')
                ax[k, 0].tick_params(labelsize=18)
            else:

                for i in range(len(ds['alt'])):
                    im = ax[k-(ngases/2), 1].plot(ds['avkSCF'][i,:],ds['alt'])

                ax[k-(ngases/2), 1].set_color_cycle([scalarMap.to_rgba(x) for x in ds['alt']])
                ax[k-(ngases/2), 1].set_ylabel('Altitude [km]',  fontsize=18)
                if k == ngases-1: ax[k-(ngases/2), 1].set_xlabel('Averaging Kernels',  fontsize=18)
                ax[k-(ngases/2), 1].grid(True)        
                #ax[k-(ngases/2), 1].set_title('Averaging Kernels Scale Factor')
                ax[k-(ngases/2), 1].set_ylim(0, 50)
                ax[k-(ngases/2), 1].annotate(gasname, xy=(.95, .89), xycoords='axes fraction', fontsize=30,
                    horizontalalignment='right', verticalalignment='center')
                ax[k-(ngases/2), 1].tick_params(labelsize=18)

            fig.subplots_adjust(bottom=0.07, top=0.92, left = 0.07, right = 0.97) 

            cbar = fig.add_axes([0.2, 0.96, 0.6, 0.03])
            cbar = fig.colorbar(scalarMap,orientation='horizontal', cax=cbar)
            cbar.set_label('Altitude [km]')

        else:
            for i in range(len(ds['alt'])):
                im = ax.plot(ds['avkSCF'][i,:],ds['alt'])

            ax.set_color_cycle([scalarMap.to_rgba(x) for x in ds['alt']])
            ax.set_ylabel('Altitude [km]', fontsize=14)
            ax.set_xlabel('Averaging Kernels', fontsize=14)
            ax.grid(True)        
            if gas.PrimaryGas == 'H2O': ax.set_ylim(0, 25)
            else: ax.set_ylim(0, 50)
            ax.annotate(gasname, xy=(.95, .95), xycoords='axes fraction', fontsize=26,
                    horizontalalignment='right', verticalalignment='center')

            cbar = fig.add_axes([0.2, 0.96, 0.6, 0.03])
            cbar = fig.colorbar(scalarMap,orientation='horizontal', cax=cbar)
            cbar.set_label('Altitude [km]')

            fig.subplots_adjust(bottom=0.08, top=0.89, left = 0.1, right = 0.96)


        #-----------------------------PLOT2---------------------------
        if ngases >= 2:
            if k <= ((ngases/2) - 1):
            
                for i in range(len(ds['alt'])):
                    im = ax1[k, 0].plot(ds['avkSCF'][i,:],ds['alt'])

                ax1[k, 0].set_color_cycle([scalarMap.to_rgba(x) for x in ds['alt']])
                ax1[k, 0].set_ylabel('Altitude [km]')
                if k == (ngases/2)-1: ax1[k, 0].set_xlabel('Averaging Kernels')
                ax1[k, 0].grid(True)        
                #ax[k, 0].set_title('Averaging Kernels Scale Factor')
                ax1[k, 0].set_ylim(0, 50)
                ax1[k, 0].set_xlim(-0.05, 0.25)
                ax1[k, 0].annotate(gasname, xy=(.95, .89), xycoords='axes fraction', fontsize=24,
                    horizontalalignment='right', verticalalignment='center')
            else:

                for i in range(len(ds['alt'])):
                    im = ax1[k-(ngases/2), 1].plot(ds['avkSCF'][i,:],ds['alt'])

                ax1[k-(ngases/2), 1].set_color_cycle([scalarMap.to_rgba(x) for x in ds['alt']])
                ax1[k-(ngases/2), 1].set_ylabel('Altitude [km]')
                if k == ngases-1: ax1[k-(ngases/2), 1].set_xlabel('Averaging Kernels')
                ax1[k-(ngases/2), 1].grid(True)        
                #ax1[k-(ngases/2), 1].set_title('Averaging Kernels Scale Factor')
                ax1[k-(ngases/2), 1].set_ylim(0, 50)
                ax1[k-(ngases/2), 1].set_xlim(-0.05, 0.25)
                ax1[k-(ngases/2), 1].annotate(gasname, xy=(.95, .89), xycoords='axes fraction', fontsize=24,
                    horizontalalignment='right', verticalalignment='center')

            fig1.subplots_adjust(bottom=0.07, top=0.92, left = 0.07, right = 0.97) 

            cbar = fig1.add_axes([0.2, 0.96, 0.6, 0.03])
            cbar = fig1.colorbar(scalarMap,orientation='horizontal', cax=cbar)
            cbar.set_label('Altitude [km]')

#-----------------------------PLOT3---------------------------
        if pltInputs['pCols']:

            if ngases >= 2:

                Nz = np.shape(ds['avkSCF'])[0]

                for cc, pcol in enumerate(pltInputs['pCols']):
                    ind1 = dc.nearestind(pcol[0], ds['alt'])
                    ind2 = dc.nearestind(pcol[1], ds['alt'])

                    q = np.zeros(Nz)
                    q[ind2:ind1] = 1.0 
                    mat = np.asarray(ds['avkSCF'])
                    smak = np.dot(q, mat)

                    if k <= ((ngases/2) - 1):
                        
                        im = ax2[k, 0].plot(smak,ds['alt'], color= clr[cc], label = str(pcol[0])+' - '+str(pcol[1])+ ' km')
                        ax2[k, 0].set_color_cycle([scalarMap.to_rgba(x) for x in ds['alt']])
                        ax2[k, 0].set_ylabel('Altitude [km]')      
                        if k == (ngases/2)-1: ax2[k, 0].set_xlabel('Summed partial kernels')
                        ax2[k, 0].grid(True)        
                        ax2[k, 0].set_ylim(0, 50)
                        if k == 0: ax2[k,0].legend(prop={'size':14}, loc='upper left')
                        ax2[k, 0].annotate(gasname, xy=(.95, .89), xycoords='axes fraction', fontsize=24,
                        horizontalalignment='right', verticalalignment='center')


                    else:
                        im = ax2[k-(ngases/2), 1].plot(smak,ds['alt'], color= clr[cc], label = str(pcol[0])+' - '+str(pcol[1])+ ' km')
                        ax2[k-(ngases/2), 1].set_color_cycle([scalarMap.to_rgba(x) for x in ds['alt']])
                        ax2[k-(ngases/2), 1].set_ylabel('Altitude [km]')      
                        if k == ngases-1: ax2[k-(ngases/2), 1].set_xlabel('Summed partial kernels')
                        ax2[k-(ngases/2), 1].grid(True)        
                        ax2[k-(ngases/2), 1].set_ylim(0, 50)
                        ax2[k-(ngases/2), 1].annotate(gasname, xy=(.95, .89), xycoords='axes fraction', fontsize=24,
                        horizontalalignment='right', verticalalignment='center')

                fig2.subplots_adjust(bottom=0.07, top=0.97, left = 0.07, right = 0.97)

            else:

                Nz = np.shape(ds['avkSCF'])[0]

                for cc, pcol in enumerate(pltInputs['pCols']):
                    ind1 = dc.nearestind(pcol[0], ds['alt'])
                    ind2 = dc.nearestind(pcol[1], ds['alt'])

                    q = np.zeros(Nz)
                    q[ind2:ind1] = 1.0 
                    mat = np.asarray(ds['avkSCF'])
                    smak = np.dot(q, mat)

                        
                    im = ax2.plot(smak,ds['alt'], color= clr[cc], label = str(pcol[0])+' - '+str(pcol[1])+ ' km')
                    ax2.set_color_cycle([scalarMap.to_rgba(x) for x in ds['alt']])
                    ax2.set_ylabel('Altitude [km]', fontsize=14)      
                    ax2.set_xlabel('Summed partial kernels', fontsize=14)
                    ax2.grid(True)        
                    if gas.PrimaryGas == 'H2O': ax2.set_ylim(0, 25)
                    else: ax2.set_ylim(0, 50)
                    ax2.legend(prop={'size':14}, loc='upper left')
                    ax2.annotate(gasname, xy=(.95, .95), xycoords='axes fraction', fontsize=26,
                    horizontalalignment='right', verticalalignment='center')

                fig2.subplots_adjust(bottom=0.08, top=0.95, left = 0.1, right = 0.97)


#-----------------------------PLOT4---------------------------

        if pltInputs['pCols']:

            if ngases >= 2:

                Nz = np.shape(ds['avkSCF'])[0]

                for cc, pcol in enumerate(pltInputs['pCols']):
                    ind1 = dc.nearestind(pcol[0], ds['alt'])
                    ind2 = dc.nearestind(pcol[1], ds['alt'])

                    q = np.zeros(Nz)
                    q[ind2:ind1] = 1.0 
                    mat = np.asarray(ds['avkSCF'])
                    smak = np.dot(q, mat)

                    if k <= ((ngases/2) - 1):
                        
                        im = ax22[k, 0].plot(smak,ds['alt'], color= clr[cc], label = str(pcol[0])+' - '+str(pcol[1])+ ' km')
                        ax22[k, 0].set_color_cycle([scalarMap.to_rgba(x) for x in ds['alt']])
                        ax22[k, 0].set_ylabel('Altitude [km]')      
                        if k == (ngases/2)-1: ax22[k, 0].set_xlabel('Summed partial kernels')
                        ax22[k, 0].grid(True)        
                        ax22[k, 0].set_ylim(0, 50)
                        ax22[k, 0].set_xlim(-0.2, 1.2)
                        if k == 0: ax22[k,0].legend(prop={'size':14}, loc='upper left')
                        ax22[k, 0].annotate(gasname, xy=(.95, .89), xycoords='axes fraction', fontsize=24,
                        horizontalalignment='right', verticalalignment='center')


                    else:
                        im = ax22[k-(ngases/2), 1].plot(smak,ds['alt'], color= clr[cc], label = str(pcol[0])+' - '+str(pcol[1])+ ' km')
                        ax22[k-(ngases/2), 1].set_color_cycle([scalarMap.to_rgba(x) for x in ds['alt']])
                        ax22[k-(ngases/2), 1].set_ylabel('Altitude [km]')      
                        if k == ngases-1: ax22[k-(ngases/2), 1].set_xlabel('Summed partial kernels')
                        ax22[k-(ngases/2), 1].grid(True)
                        ax22[k-(ngases/2), 1].set_xlim(-0.2, 1.2)        
                        ax22[k-(ngases/2), 1].set_ylim(0, 50)
                        ax22[k-(ngases/2), 1].annotate(gasname, xy=(.95, .89), xycoords='axes fraction', fontsize=24,
                        horizontalalignment='right', verticalalignment='center')

                fig22.subplots_adjust(bottom=0.07, top=0.97, left = 0.07, right = 0.97)

    if gas.pdfsav:
        gas.pdfsav.savefig(fig,dpi=200)
        if ngases >= 2: gas.pdfsav.savefig(fig1,dpi=200)
        gas.pdfsav.savefig(fig2,dpi=200)
        if ngases >= 2: gas.pdfsav.savefig(fig22,dpi=200)  
        
    else:           
        plt.show(block=False)          
                    

    #-------------------------------------------------------------------------
    #END PLOTS
    #-------------------------------------------------------------------------
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
