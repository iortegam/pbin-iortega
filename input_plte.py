#----------------------------------------------------------------------------------------
# Name:
#        plte_input.py
#
# Purpose:
#        This is the input file for plte.py   --> Plot error for a serie of target gases and flags defined below
#
#----------------------------------------------------------------------------------------


loc = 'fl0'               # Name of station location

#------------------------------------TAB----------------------------------
if loc == 'tab':

    #gasName    = ['h2o']
    #ver        = ['Current_v1']
    #ctlF       = ['sfit4_v1.ctl']

    
    #gasName    = ['ocs']
    #ver        = ['Current_v2']
    #ctlF       = ['sfit4_v2.ctl']  
      
    gasName    = ['c2h6', 'ch4', 'hcn', 'co', 'hcl', 'hf', 'hno3', 'n2o', 'clono2', 'o3' ]                   
    ver        = ['Current_NewSA_NewScale2', 'Current', 'Current',  'Current','Current', 'Current', 'Current', 'Current','Current', 'Current']           
    ctlF       = ['sfit4_3.ctl', 'sfit4.ctl','sfit4.ctl', 'sfit4.ctl', 'sfit4.ctl','sfit4_h2o.ctl', 'sfit4_h2o.ctl','sfit4_h2o.ctl', 'sfit4_h2o.ctl', 'sfit4.ctl']


elif loc == 'mlo':

    #gasName    = ['h2o']
    #ver        = ['Current_v1']
    #ctlF       = ['sfit4_v1.ctl']    
        
    #gasName    = ['ocs']
    #ver        = ['Current_v2']
    #ctlF       = ['sfit4_v2.ctl']            
    
    gasName    = ['c2h6', 'ch4', 'hcn', 'co', 'hcl', 'hf', 'hno3', 'n2o', 'clono2', 'o3' ]                  
    ver        = ['Current_newSA', 'Current', 'Current','Current','Current', 'Current', 'Current', 'Current','Current', 'Current']           
    ctlF       = ['sfit4.ctl', 'sfit4.ctl','sfit4_h2o.ctl', 'sfit4_h2o.ctl', 'sfit4.ctl', 'sfit4_h2o.ctl', 'sfit4_h2o.ctl','sfit4_h2o.ctl', 'sfit4_h2o.ctl', 'sfit4.ctl' ]           
    
elif loc == 'fl0':

    gasName    = ['c2h6']
    ver        = ['Current']
    ctlF       = ['sfit4.ctl']

    #gasName    = ['no2']
    #ver        = ['Current_v2']
    #ctlF       = ['sfit4_v2.ctl']
    

    #gasName    = ['c2h6', 'ch4', 'hcn', 'c2h2', 'h2co' ,'co', 'nh3' ,'o3']                
    #ver        = ['Current_WP',  'Current_WP',   'Current_WP', 'Current_WP', 'Current_WP_v6', 'Current_WP','Current_WP', 'Current_WP']           #
    #ctlF       = ['sfit4.ctl',  'sfit4_3.ctl', 'sfit4.ctl', 'sfit4.ctl', 'sfit4_v6.ctl', 'sfit4.ctl', 'sfit4_Enrico.ctl', 'sfit4.ctl'] 

else:
    print 'You need to input a name of station location!'
    exit() 

#------
# Flags
#------
saveFlg    = False	                   # Flag to either save data to pdf file (saveFlg=True) or plot to screen (saveFlg=False)
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
   
        maxRMS     = [2.2, 1.5, 2.0, 1.5, 1.5, 1.5, 1.0, 1.5, 1.0, 3.2]                       # Max Fit RMS to filter data. Data is filtered according to <= maxrms
        minDOF     = [0.5, 0.9, 1.5, 0.5, 0.9, 0.5, 0.5, 0.9, 0.0, 2.5] 


elif loc == 'mlo': 
   
        maxRMS     = [1.0, 1.0, 1.0, 1.5, 1.0, 1.5, 1.5, 1.0, 2.0, 2.5]                       # Max Fit RMS to filter data. Data is filtered according to <= maxrms
        minDOF     = [0.9, 0.9, 0.9, 0.9, 0.9, 0.3, 0.5, 0.9, 0.0, 3.0]                    # Min DOFs for filtering


elif loc == 'fl0':
   
        maxRMS     = [2.2, 0.3, 0.6, 0.8, 5.0, 1.0, 1.0, 4.0]                       # Max Fit RMS to filter data. Data is filtered according to <= maxrms
        minDOF     = [0.5, 0.9, 0.9, 0.8, 0.0, 2,0, 0.4, 3.0]
     
        #maxRMS     = [10, 10, 10, 10]                       # Max Fit RMS to filter data. Data is filtered according to <= maxrms
        #minDOF     = [0, 0, 0, 0]
  

else:
    print 'You need to input a location!'
    exit() 

minSZA     = 0.0                   # Min SZA for filtering
maxSZA     = 90.0                   # Max SZA for filtering
maxCHI     = 2.0                    # Max CHI_y_2 value
maxTC      = 5.0E25                # Max Total column amount for filtering
minTC      = 0.0                 # Min Total column amount for filtering
sclfct     = 1.0E9                  # Scale factor to apply to vmr plots (ppmv=1.0E6, ppbv=1.0E9, etc)
sclfctName = 'ppbv'                 # Name of scale factor for labeling plots

#----------------------
# Date range to process
#----------------------
if loc == 'tab': yyyy   = 2012
elif loc == 'mlo': yyyy = 2013
elif loc == 'fl0': yyyy = 2015
else: 
    print 'Input year'
    exit()
 
iyear      = yyyy
imnth      = 1
iday       = 1
fyear      = yyyy
fmnth      = 6
fday       = 30

#----------------------------
# Partial Columns Bounds [km] 
# [lower bound, upper bound]
# To turn off set to False
#----------------------------
pCols = [[0.0,8.0],[8.0,16.0],[16.0,25.0]]
# Directories
#------------ -
#retDir = '/Volumes/data1/ebaumer/'+loc.lower()+'/'+gasName.lower()+'/'+ver+'/'  
retDir = '/data1/ebaumer/'+loc[0].lower()+'/'+gasName[0].lower()+'/'+ver[0]+'/'  

#------
# Files
#-----
#ctlFile  = '/Volumes/data1/ebaumer/'+loc.lower()+'/'+gasName.lower()+'/'+'x.'+gasName.lower()+'/'+ctlF
ctlFile  = '/data1/ebaumer/'+loc.lower()+'/'+gasName[0].lower()+'/'+'x.'+gasName[0].lower()+'/'+ctlF[0]
#pltFile  = '/Volumes/data1/ebaumer/'+loc.lower()+'/' + 'Plots/' + loc + '_' + gasName + '_' + ver + '.pdf'

if saveFlg:
    if len(gasName) == 1: 
        pltFile  = '/data/iortega/results/'+loc.lower()+'/fig/'+loc.lower()+ '_'+ gasName[0] +'_errors_'+str(iyear)+'_old.pdf'
        
    else:
        pltFile  = '/data/iortega/results/'+loc.lower()+'/fig/'+loc.lower()+ '_'+'errors_'+str(iyear)+'_v2.pdf'

else:
    pltFile  = ''

       
