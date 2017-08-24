#----------------------------------------------------------------------------------------
# Name:
#        input_pltts.py --> input for pltts.py (Plot time series of a series of gases)
#
# Purpose:
#        This is the input file for pltts.py 
#----------------------------------------------------------------------------------------
loc = 'fl0'               # Name of station location

if loc == 'tab':
    gasName    = ['c2h6', 'ch4', 'co', 'hcl', 'hcn']                   # Name of gas
    ver        = ['Current_NewSA_NewScale2', 'Current_WP', 'Current_WP','Current_WP','Current_WP']           # Name of retrieval version to process
    ctlF       = ['sfit4_3.ctl', 'sfit4.ctl','sfit4.ctl', 'sfit4.ctl', 'sfit4.ctl']           # Name of ctl file

elif loc == 'mlo':            
    gasName    = ['c2h6', 'ch4', 'co', 'hcl', 'hcn']                   # Name of gas
    ver        = ['Current_newSA', 'Current_WP', 'Current_WP','Current_WP','Current_WP']           # Name of retrieval version to process
    ctlF       = ['sfit4.ctl', 'sfit4.ctl','sfit4_h2o.ctl', 'sfit4.ctl', 'sfit4_h2o.ctl']           # Name of ctl file

elif loc == 'fl0':
    gasName    = ['c2h6', 'ch4', 'c2h2', 'h2co' ,'co', 'nh3', 'o3']                   # Name of gas
    ver        = ['Current_WP', 'Current_WP', 'Current_WP', 'Current_WP_v6', 'Current_WP','Current_WP', 'Current_WP']           # Name of retrieval version to process
    ctlF       = ['sfit4.ctl', 'sfit4_3.ctl','sfit4.ctl', 'sfit4_v6.ctl', 'sfit4.ctl', 'sfit4_Enrico.ctl', 'sfit4.ctl']           # Name of ctl file

    #gasName    = ['c2h6', 'ch4', 'co', 'nh3']                  # Name of gas
    #ver        = ['Current_WP', 'Current_WP', 'Current_WP', 'Current_WP'  ]#, 'Current_WP', 'Current_WP_v6', 'Current_WP','Current_WP', 'Current_WP']           # Name of retrieval version to process
    #ctlF       = ['sfit4.ctl', 'sfit4_3.ctl', 'sfit4.ctl' , 'sfit4_Enrico.ctl'  ]#,'sfit4.ctl', 'sfit4_v6.ctl', 'sfit4.ctl', 'sfit4_Enrico.ctl', 'sfit4.ctl']           # Name of ctl file


else:
    print 'You need to input a location!'
    exit() 

#----------------------------
# Flags
#----------------------------
saveFlg    = True	                   # Flag to either save data to pdf file (saveFlg=True) or plot to screen (saveFlg=False)
errorFlg   = True                  # Flag to process error data
fltrFlg    = True                   # Flag to filter the data
byYrFlg    = False                  # Flag to create plots for each individual year in date range
szaFlg     = True                   # Flag to filter based on min and max SZA
dofFlg     = True                   # Flag to filter based on min DOFs
pcNegFlg   = False                   # Flag to filter profiles with negative partial columns
tcNegFlg   = False                  # Flagsag to filter profiles with negative total columns
tcMMFlg    = False                  # Flag to filter based on min and max total column amount
cnvrgFlg   = True                   # Flag to filter profiles that did not converge
rmsFlg     = True                   # Flag to filter based on max RMS
chiFlg     = False     	            # Flag to filter based on max CHI_2_Y
mnthFlg    = False                  # Flag to filter based on 

mnths      = [6,7,8]                # Months to filter on (these are the months to include data)

if loc == 'tab':
    maxRMS     = [2.2, 1.5, 2, 1.5, 1.5]                       # Max Fit RMS to filter data. Data is filtered according to <= maxrms
    minDOF     = [0.5, 0.9, 1.5, 0.5, 0.9]                    # Min DOFs for filtering

elif loc == 'mlo': 
    maxRMS     = [1.0, 1.0, 1.5, 1.0, 1.0]                       # Max Fit RMS to filter data. Data is filtered according to <= maxrms
    minDOF     = [0.9, 0.9, 0.9, 0.9, 0.9]                    # Min DOFs for filtering

elif loc == 'fl0':
    maxRMS     = [2.2, 0.3, 0.8, 0.6 , 1.0, 0.6, 4]                       # Max Fit RMS to filter data. Data is filtered according to <= maxrms
    minDOF     = [0.5, 0.9, 0.8, 0.5 , 2,0, 0.9, 3]

    #maxRMS     = [2.2, 0.3, 1.0, 0.6]                        # Max Fit RMS to filter data. Data is filtered according to <= maxrms
    #minDOF     = [0.5, 0.9, 2.0, 0.9]

else:
    print 'You need to input the maxRMS and minDOF!'
    exit() 

minSZA     = 0.0                   # Min SZA for filtering
maxSZA     = 90.0                   # Max SZA for filtering
maxCHI     = 2.0                    # Max CHI_y_2 value
maxTC      = 5.0E25                # Max Total column amount for filtering
minTC      = 0.0                 # Min Total column amount for filtering
sclfct     = 1.0E9                  # Scale factor to apply to vmr plots (ppmv=1.0E6, ppbv=1.0E9, etc)
sclfctName = 'ppbv'                 # Name of scale factor for labeling plots

#----------------------------
# Date range to process
#----------------------------
iyear      = 2010
imnth      = 1
iday       = 1
fyear      = 2016
fmnth      = 1
fday       = 31

#----------------------------
# Partial Columns Bounds [km] 
# [lower bound, upper bound]
# To turn off set to False
#----------------------------
pCols = [[0.0,8.0]]

#--------------
# Directories
#--------------
#retDir = '/Volumes/data1/ebaumer/'+loc.lower()+'/'+gasName.lower()+'/'+ver+'/'  
retDir = '/data1/ebaumer/'+loc[0].lower()+'/'+gasName[0].lower()+'/'+ver[0]+'/'  

#----------------------------
# Files
#----------------------------
#ctlFile  = '/Volumes/data1/ebaumer/'+loc.lower()+'/'+gasName.lower()+'/'+'x.'+gasName.lower()+'/'+ctlF
ctlFile  = '/data1/ebaumer/'+loc.lower()+'/'+gasName[0].lower()+'/'+'x.'+gasName[0].lower()+'/'+ctlF[0]
#pltFile  = '/Volumes/data1/ebaumer/'+loc.lower()+'/' + 'Plots/' + loc + '_' + gasName + '_' + ver + '.pdf'

if saveFlg: 
    pltFile  = '/data/iortega/results/'+loc.lower()+'/fig/'+loc.lower()+'_TG2.pdf'
else:
    pltFile  = ''

#----------------------------
# Initializations for weather (wind, temp) in case is needed
#----------------------------
weatherFlg     = True
weatherDir     = '/data1/ancillary_data/fl0/eol/'
weatherFileTag = 'v2'

