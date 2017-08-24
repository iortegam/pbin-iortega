#----------------------------------------------------------------------------------------
# Name:
#        input_pltsonde.py --> 

# Purpose:
#        This is the input file for pltsonde.py   --> Plot sonde vs FTS 
#----------------------------------------------------------------------------------------
loc = 'fl0'               
    
if loc.lower() == 'mlo':
   
    gasName    = 'h2o'                                  
    ver        = 'Current_ERA'          
    ctlF       = 'sfit4.ctl'
    doi        = ['20120706', '20120831', '20130116', '20130410', '20130626', '20131031', '20140212', '20141020', '20150326', '20151026']
   
elif loc.lower() == 'fl0':

    gasName    = 'h2o'                                  
    ver        = 'Current_ERA'          
    ctlF       = 'sfit4.ctl'
    doi        = ['20100914', '20101105', '20110728', '20120430', '20130117', '20130816', '20140722', '20150629', '20150805', '20160706']

else:
    print 'You need to input a location!'
    exit()

smthFlg    = False                #Smooth sonde profiles using AK

#------
# Flags
#------
saveFlg    = False	               # Flag to either save data to pdf file (saveFlg=True) or plot to screen (saveFlg=False)
errorFlg   = True                  # Flag to process error data
fltrFlg    = True                  # Flag to filter the data
byYrFlg    = False                  # Flag to create plots for each individual year in date range
szaFlg     = True                   # Flag to filter based on min and max SZA
dofFlg     = True                   # Flag to filter based on min DOFs
pcNegFlg   = True                   # Flag to filter profiles with negative partial columns
tcNegFlg   = True                  # Flagsag to filter profiles with negative total columns
tcMMFlg    = False                  # Flag to filter based on min and max total column amount
cnvrgFlg   = True                   # Flag to filter profiles that did not converge
rmsFlg     = True                   # Flag to filter based on max RMS
chiFlg     = False     	            # Flag to filter based on max CHI_2_Y
mnthFlg    = False                  # Flag to filter based on 

mnths      = [6,7,8]                # Months to filter on (these are the months to include data)

if loc.lower() == 'mlo': 
    maxRMS     = 1.2                      # Max Fit RMS to filter data. Data is filtered according to <= maxrms
    minDOF     = 1.5                   # Min DOFs for filtering

elif loc.lower() == 'fl0':

    maxRMS     = 1.0                      # Max Fit RMS to filter data. Data is filtered according to <= maxrms
    minDOF     = 1.7
    latFTS     = 40.04
    lonFTS     = -105.24

else:
    print 'You need to input the maxRMS and minDOF!'
    exit() 

minSZA     = 0.0                   # Min SZA for filtering
maxSZA     = 90.0                   # Max SZA for filtering
maxCHI     = 2.0                    # Max CHI_y_2 value
maxTC      = 5.0E25                # Max Total column amount for filtering
minTC      = 0.0                 # Min Total column amount for filtering
sclfct     = 1.0E6                  # Scale factor to apply to vmr plots (ppmv=1.0E6, ppbv=1.0E9, etc)
sclfctName = 'ppmv'                 # Name of scale factor for labeling plots

#----------------------
# Date range to process
#----------------------
iyear      = 2010
imnth      = 1
iday       = 1
fyear      = 2016
fmnth      = 12
fday       = 31

#----------------------------
# Partial Columns Bounds [km] 
# [lower bound, upper bound]
# To turn off set to False
#----------------------------
#if loc.lower() == 'mlo': pCols = [[3.8, 7.0], [7.0, 10.0], [10.0, 13.0], [13.0, 16.0]]
#if loc.lower() == 'fl0': pCols = [[1.0, 4.0], [4.0, 8.0], [8.0, 12.0], [12.0, 16.0]]

if loc.lower() == 'mlo': pCols = [ [3.0, 5.5], [5.5, 8], [8, 10.5], [10.5, 13.5], [13.5, 17], [17, 21]]
if loc.lower() == 'fl0': pCols = [ [1.5, 3], [3, 5], [5, 7.5], [7.5, 10], [10,13], [13,16], [16,20] ]


# Directories
#------------ -
retDir = '/data1/ebaumer/'+loc.lower()+'/'+gasName.lower()+'/'+ver+'/'
ctlFile  = '/data1/ebaumer/'+loc.lower()+'/'+gasName.lower()+'/'+'x.'+gasName.lower()+'/'+ctlF

sondeDir     = '/data1/ancillary_data/sondes/'+loc.lower()+'/'

#------
# Files
#-----

if saveFlg: 
    if smthFlg: pltFile  = '/data/iortega/results/'+loc.lower()+'/fig/'+loc.lower()+'_'+gasName.lower()+'_sonde_vs_fts_'+ver+'_smooth.pdf'
    else:       pltFile  = '/data/iortega/results/'+loc.lower()+'/fig/'+loc.lower()+'_'+gasName.lower()+'_sonde_vs_fts_'+ver+'.pdf'
else:
    pltFile  = ''
    
