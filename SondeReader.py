#----------------------------------------------------------------------------------------
# Name:
#        SondeReader.py
#
# Purpose:
#       Reads the sonde profiles measured from NOAA 
#
#----------------------------------------------------------------------------------------

import datetime as dt
import re
import csv
from itertools import islice
import sys
import subprocess  


class WaterReads():
    '''  Class for reading water sonde profiles '''    
    def __init__(self):
        ''' Initializations '''
        #-------------------------------------
        # initialize dictionary for house data
        #-------------------------------------
        self.data = {'Level_Number':[],                         # Internallly used for sorting entries
                     'Alt_km':[],                             # 
                     'Press_hPa':[],
                     'Temp_degC':[],
                     'Theta_degK':[],
                     'RH_FPH':[],
                     'RH_RS':[],
                     'H2Omr_ppmv':[],
                     'H2Osd_ppmv':[],
                     'H2Omr_orig_ppmv':[],
                     'TFP_FPH_decC':[],
                     'TFP_RS_degC':[],
                     'O3mr_ppbv':[],
                     'Time':[]}
        
            
    def SondeA(self,fileName,year,month,day):
        ''' '''
        #------------------------
        # Open file and read data
        #------------------------
        print 'Reading file:', fileName
        with open(fileName,'r') as fname:
            try:
                data = fname.readlines()
            except:
                print 'Error in reading file: %s' % fileName
                sys.exit()   
      
        #--------------------------
        # Remove header information 
        #--------------------------
        info = [ row.strip().split() for row in data if 'Water Vapor Flight Date' in row]              
        data[:] = [ row.strip().split() for row in data if len(row.strip().split()) == 13 and not 'to' in row and not 'Level' in row and not 'Traditionally' in row and not 'Number' in row]              

        #---------------------------------
        # Determine number of observations
        #---------------------------------
        npoints = len(data)
       
        time = [row[7].split()[0] for row in info]
        #time[:]  = info[0]

        #ihour   = str(time[0][0:2])
        #imin    = str(time[0][3:2])
        #isec    = str(time[0][6:2]) 
        self.data['Time'].extend([time]*npoints) 
        self.data['Level_Number'].extend([row[0] for row in data])    
        self.data['Alt_km'].extend([row[1] for row in data])  
        self.data['Press_hPa'].extend([row[2] for row in data])   
        self.data['Temp_degC'].extend([row[3] for row in data])   
        self.data['Theta_degK'].extend([row[4] for row in data])   
        self.data['RH_FPH'].extend([row[5] for row in data])   
        self.data['RH_RS'].extend([row[6] for row in data])           
        self.data['H2Omr_ppmv'].extend([row[7] for row in data])
        self.data['H2Osd_ppmv'].extend([row[8] for row in data])
        self.data['H2Omr_orig_ppmv'].extend([row[9] for row in data])
        self.data['TFP_FPH_decC'].extend([row[10] for row in data])
        self.data['TFP_RS_degC'].extend([row[11] for row in data])
        self.data['O3mr_ppbv'].extend([row[12] for row in data])

