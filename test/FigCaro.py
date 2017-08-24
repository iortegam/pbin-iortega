# - *- coding: utf- 8 - *-
###coding=utf-8

#------------------------------------------------------
#SCRIPT TO PLOT THE SPI (CARO'S REQUEST)
#------------------------------------------------------

import sys
import datetime as dt
import matplotlib
from numpy import fromstring, vectorize, ndarray, array, genfromtxt
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
import numpy as np


#-------------------------------------------------------
#READNDATA
#-------------------------------------------------------
#File = '/data/iortega/pbin/test/caro_table2.txt'
File = 'caro_table2.txt'
with open(File ,'r') as fopen: lines = fopen.readlines()

Data = [row.strip().split() for row in lines]

Time   = []
Value  = []

for i, row in enumerate(Data):
    if len(row) == 2:

        if int(len(row[0])) == 8:
        
            Value.append(float(row[1]))
            yyyy = int(row[0][4:8].strip())
            mm   = int(row[0][0:1].strip())
            dd   = int(row[0][2:3].strip())
            Time.append(dt.date(yyyy, mm, dd))

        elif int(len(row[0])) == 9:

            Value.append(float(row[1]))
            yyyy = int(row[0][5:9].strip())
            mm   = int(row[0][0:2].strip())
            dd   = int(row[0][3:4].strip())
            Time.append(dt.date(yyyy, mm, dd))

    else:
        Time.append(dt.date(yyyy+1, 1, 1))
        Value.append(float(row[0]))


Time = np.asarray(Time)
Value = np.asarray(Value)

#-------------------------------------------------------
#DEFINING COLORS FOR PLOT
#-------------------------------------------------------
Colors = []

for val in Value:
    if val <= -1.5: Colors.append('maroon')# chocolate
    elif (val > -1.5) & (val <= 0.0): Colors.append('peru')
    elif (val > 0.0) & (val <= 1.5): Colors.append('lightgreen')    #limegreen
    elif val > 1.5  : Colors.append('darkgreen')


Colors = np.asarray(Colors)
#-------------------------------------------------------
#
#-------------------------------------------------------
yearsall = [ singDate.year for singDate in Time]               # Find years for all date entries
years = list(set(yearsall))                                             # Determine all unique years
years.sort()
Nyears = len(years)

ind = np.arange(len(Time)) 
N = len(ind)

#-------------------------------------------------------
#PLOT
#-------------------------------------------------------

fig, ax = plt.subplots(figsize=(10,6))
sc = ax.bar(ind, Value, color=Colors, edgecolor="none")

ax.set_ylabel('SPI-12', fontsize=14)
ax.set_xticks(np.arange(0, N, 72))
ax.tick_params(labelsize=14)

xt = ax.get_xticks()
new_xticks = [yearsall[d] for d in xt]
ax.set_xticklabels(new_xticks, rotation=45)
#ax.set_xlabel('A\~{n}o', fontsize=13)
ax.set_xlabel(u'Año', fontsize=14)
ax.set_ylim(-3.5, 3)

fig.subplots_adjust(bottom=0.15,top=0.95, left=0.075, right=0.95)

plt.show(block=False)

user_input = raw_input('Press any key to exit >>> ')
sys.exit()

