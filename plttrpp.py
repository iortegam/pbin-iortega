#! /usr/bin/python2.7
#----------------------------------------------------------
#READ THE  DAYLY TROPOSPHERIC HEIGHT USING YEARLY FILES (ERIC CREATED)
#----------------------------------------------------------

#---------------
# Import modules
#---------------
import sys
import time
import datetime as dt
import matplotlib.dates as md
from matplotlib.dates import DateFormatter, MonthLocator, YearLocator
from matplotlib.ticker import FormatStrFormatter, MultipleLocator,AutoMinorLocator,ScalarFormatter
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from matplotlib.backends.backend_pdf import PdfPages 
import numpy as np
from scipy import stats
import myfunctions as mf

from sklearn import linear_model, datasets

def main(argv):
	#----------------------------------------------------------
	#INPUTS
	#----------------------------------------------------------
	
	loc        = 'wlo'

	savePDF    = True
	

	#trppDir  = '/data1/ancillary_data/NCEPdata/NCEP_trpp/'

	if loc == 'eur':
		trppDir  = '/data1/projects/ocs/tor'
		pltFile    = '/data1/projects/ocs/tor/NCEP_trpp/'+loc.upper()+'_TropHeight_DTH-NCEP.pdf'
	else:	
		trppDir  = '/data1/projects/ocs/'+loc.lower()
		pltFile    = '/data1/projects/ocs/'+loc.lower()+'/NCEP_trpp/'+loc.upper()+'_TropHeight_DTH-NCEP.pdf'

	if savePDF: pdfsav = PdfPages(pltFile)

	#----------------------------------------------------------
	#                           MAIN
	#----------------------------------------------------------

	if loc.lower() == 'tab':
		loc2 = 'Thule'
	
	elif loc.lower() == 'mlo':
		loc2 = 'Mauna Loa'

	elif loc.lower() == 'bld':
		loc2 = 'Boulder'

	elif loc.lower() == 'eur':
		loc2 = 'Eureka'

	elif loc.lower() == 'jfj':
		loc2 = 'Jungfraujoch'

	elif loc.lower() == 'tor':
		loc2 = 'Toronto'

	elif loc.lower() == 'wlo':
		loc2 = 'Wollongong'

   #-----------------------------------------------
   #NCEP tropopause
   #-----------------------------------------------
 
	TH_ncep   = []
	date_ncep = []
	dt_ncep   = []

	f2 = open( trppDir+'/Dates_'+loc2+'_dNCEP.ascii', 'r')
	lines = f2.readlines()
	f2.close()

	count = 0
	for line in lines:
		if count >=5:
			p = line.split()
			yyyy = int(p[0][0:4])
			mm   = int(p[0][4:6])
			dd   = int(p[0][6:8])
			hh   = int(p[1][0:2])
			mi   = int(p[1][3:5])
			ss   = int(p[1][6:8])
			date_ncep.append(dt.date(yyyy, mm, dd))
			dt_ncep.append(dt.datetime(yyyy, mm, dd, hh, mi, ss))
			TH_ncep.append(float(p[2]))

		count+= 1

	date_ncep = np.asarray(date_ncep)
	dt_ncep   = np.asarray(dt_ncep)
	TH_ncep   = np.asarray(TH_ncep)

	dailyavg        = mf.dailyAvg(TH_ncep, dt_ncep)
	TH_ncep_da      = dailyavg['dailyAvg']
	date_ncep_da    = dailyavg['dates']

	mnthlyyavg       = mf.mnthlyAvg(TH_ncep, dt_ncep)
	TH_ncep_mth      = mnthlyyavg['mnthlyAvg']
	TH_sdev_ncep_mth = mnthlyyavg['std']
	date_ncep_mth    = mnthlyyavg['dates']


   #-----------------------------------------------
   #Dynamic tropopause
   #-----------------------------------------------
	TH_dtp   = []
	date_dtp = []
	dt_dtp   = []

	f2 = open( trppDir+'/Dates_'+loc2+'_dtp.ascii', 'r')
	lines = f2.readlines()
	f2.close()

	count = 0
	for line in lines:
		if count >=5:
			p = line.split()
			yyyy = int(p[0][0:4])
			mm   = int(p[0][4:6])
			dd   = int(p[0][6:8])
			hh   = int(p[1][0:2])
			mi   = int(p[1][3:5])
			ss   = int(p[1][6:8])
			date_dtp.append(dt.date(yyyy, mm, dd))
			dt_dtp.append(dt.datetime(yyyy, mm, dd, hh, mi, ss))
			TH_dtp.append(float(p[2]))

		count+= 1

	date_dtp = np.asarray(date_dtp)
	dt_dtp   = np.asarray(dt_dtp)
	TH_dtp   = np.asarray(TH_dtp)

	dailyavg       = mf.dailyAvg(TH_dtp, dt_dtp)
	TH_dtp_da      = dailyavg['dailyAvg']
	TH_sdev_dtp_da = dailyavg['std']
	date_dtp_da    = dailyavg['dates']

	mnthlyyavg       = mf.mnthlyAvg(TH_dtp, dt_dtp)
	TH_dtp_mth      = mnthlyyavg['mnthlyAvg']
	TH_sdev_dtp_mth = mnthlyyavg['std']
	date_dtp_mth    = mnthlyyavg['dates']

    #---------------------------------
	date_dtp_da_set = set(date_dtp_da)
	TH_ncep_da_coi=[]

	for i, item in enumerate(date_ncep_da):
		if item in date_dtp_da_set:
			TH_ncep_da_coi.append(TH_ncep_da[i])

	#---------------------------------
	date_dtp_mth_set      = set(date_dtp_mth)
	TH_ncep_mth_coi       = []
	TH_sdev_ncep_mth_coi  = []

	for i, item in enumerate(date_ncep_mth):
		if item in date_dtp_mth_set:
			TH_ncep_mth_coi.append(TH_ncep_mth[i])
			TH_sdev_ncep_mth_coi.append(TH_sdev_ncep_mth[i])


    #-----------------------------------------------
    #plots
    #-----------------------------------------------

	clmap = 'jet'
	cm           = plt.get_cmap(clmap)
	yearsLc      = YearLocator()
	months       = MonthLocator()
	DateFmt      = DateFormatter('%Y')

	#------------------------------
	#plots 2
	#------------------------------

	fig2, ax = plt.subplots(2, figsize=(12,10), sharex=True)
	ax[0].plot(date_dtp_da, TH_dtp_da, color='k', linestyle='None', marker ='', markersize=4)
	ax[0].scatter(date_dtp_da, TH_dtp_da, facecolors='green', s=35, color='green', label = 'Dynamic', edgecolor='k')

	ax[0].plot(date_ncep_da, TH_ncep_da, color='k', linestyle='None', marker ='', markersize=4)
	ax[0].scatter(date_ncep_da, TH_ncep_da, facecolors='r', s=35, color='r', label='NCEP', edgecolor='k')

	ax[0].xaxis.set_minor_locator(months)
	ax[0].xaxis.set_major_formatter(DateFmt) 
	ax[0].xaxis.set_tick_params(which='major',labelsize=12)
	ax[0].xaxis.set_tick_params(which='minor',labelbottom='off')
	#ax[0].set_xlim((dt.date(iyear,imnth, iday), dt.date(fyear,fmnth,fday)))
	#ax[0].set_ylim(4, 14)	
	ax[0].set_ylabel('Tropospheric Height [Km]', fontsize=16)
	ax[0].tick_params(axis='both', which='major', labelsize=14)
	ax[0].grid(True)
	ax[0].set_title('Daily values (Retrieval dates)', fontsize = 16)
	ax[0].legend(prop={'size':14},loc='upper right') 

	ax[1].errorbar(date_dtp_mth,TH_dtp_mth,yerr=TH_sdev_dtp_mth, fmt='o', color='green', markersize=7, ecolor='green', capthick=2, label = 'Dynamic')
	ax[1].errorbar(date_ncep_mth,TH_ncep_mth ,yerr=TH_sdev_ncep_mth, fmt='o', color='r', markersize=7, ecolor='r', capthick=2, label='NCEP')
	ax[1].set_ylabel('Tropospheric Height [Km]', fontsize=16)
	ax[1].tick_params(axis='both', which='major', labelsize=14)
	#ax[1].set_xlim((dt.date(iyear,imnth, iday), dt.date(fyear,fmnth,fday)))
	ax[1].grid(True)
	ax[1].set_ylabel('Tropospheric Height [Km]', fontsize=16)
	ax[1].set_xlabel('Year', fontsize = 16)
	ax[1].set_title('Monthly values (retrieval months)', fontsize = 16)
	#ax[1].set_ylim(4, 14)

	ax[1].legend(prop={'size':14},loc='upper right') 

	fig2.autofmt_xdate()
	fig2.subplots_adjust(left = 0.12, bottom=0.12, top=0.96, right = 0.95)

	#------------------------------
	#plots 3
	#------------------------------
	yyyy = 2006

	fig3, ax = plt.subplots(2, figsize=(12,10), sharex=True)
	ax[0].plot(date_dtp_da, TH_dtp_da, color='k', linestyle='None', marker ='', markersize=4)
	ax[0].scatter(date_dtp_da, TH_dtp_da, facecolors='green', s=35, color='green', label = 'Dynamic', edgecolor='k')

	ax[0].plot(date_ncep_da, TH_ncep_da, color='k', linestyle='None', marker ='', markersize=4)
	ax[0].scatter(date_ncep_da, TH_ncep_da, facecolors='r', s=35, color='r', label='NCEP', edgecolor='k')
	DateFmt      = DateFormatter('%m')
	#DateFmt      = DateFormatter('%Y-%m-%d %H:%M:%S')

	#ax[0].xaxis.set_major_locator(yearsLc)
	ax[0].xaxis.set_minor_locator(months)
	ax[0].xaxis.set_major_formatter(DateFmt)
	ax[0].xaxis.set_tick_params(which='major',labelsize=12)
	ax[0].xaxis.set_tick_params(which='minor',labelbottom='off')
	#ax[0].set_ylim(4, 14)	
	ax[0].set_ylabel('Tropospheric Height [Km]', fontsize=16)
	ax[0].tick_params(axis='both', which='major', labelsize=14)
	ax[0].grid(True)
	ax[0].set_title('Daily values (Retrieval dates)', fontsize = 16)
	ax[0].legend(prop={'size':14},loc='upper right') 
	ax[0].set_xlim((dt.datetime(yyyy,2, 1, 0), dt.datetime(yyyy,11, 1, 0)))

	ax[1].errorbar(date_dtp_mth,TH_dtp_mth,yerr=TH_sdev_dtp_mth, fmt='o', color='green', markersize=7, ecolor='green', capthick=2, label = 'Dynamic')
	ax[1].errorbar(date_ncep_mth,TH_ncep_mth,yerr=TH_sdev_ncep_mth, fmt='o', color='r', markersize=7, ecolor='r', capthick=2, label='NCEP')
	ax[1].set_ylabel('Tropospheric Height [Km]', fontsize=16)
	ax[1].tick_params(axis='both', which='major', labelsize=14)
	ax[1].grid(True)
	ax[1].set_ylabel('Tropospheric Height [Km]', fontsize=16)
	ax[1].set_xlabel('Month in ' + str(yyyy), fontsize = 16)
	ax[1].set_title('Monthly values (retrieval months)', fontsize = 16)
	#ax[1].set_ylim(4, 14)

	ax[1].legend(prop={'size':14},loc='upper right') 

	fig3.autofmt_xdate()
	fig3.subplots_adjust(left = 0.12, bottom=0.12, top=0.96, right = 0.95)

    #------------------------------
	#
	#------------------------------
	fig4, ax = plt.subplots(figsize=(8,6), sharex=True)
	res = TH_dtp_da - TH_ncep_da_coi
	ind = np.where(np.abs(res) < 3.0)[0]
	ind = np.asarray(ind)

	TH_dtp_da      = TH_dtp_da[ind]
	TH_ncep_da_coi = np.asarray(TH_ncep_da_coi)
	TH_ncep_da_coi = TH_ncep_da_coi[ind]
	
	ax.scatter(TH_dtp_da,  TH_ncep_da_coi, facecolors='white', s=20, color='b', alpha=0.5, label='Daily') 
	#ax.errorbar(TH_dtp_mth,  TH_ncep_mth_coi, xerr=TH_sdev_dtp_mth, yerr=TH_sdev_ncep_mth_coi, fmt='o', color='r')
	ax.scatter(TH_dtp_mth,  TH_ncep_mth_coi, facecolors='white', s=45, color='r', label='Monthly')
	ax.set_title(loc2, fontsize=18) 

	odr, odrErr  = mf.orthoregress(np.asarray(TH_dtp_mth),  np.asarray(TH_ncep_mth_coi), xerr=np.asarray(TH_sdev_dtp_mth), yerr=np.asarray(TH_sdev_ncep_mth_coi), InError=True)
	slope_odr    = float(odr[0])
	inter_odr    = float(odr[1])
	slopeErr_odr = float(odrErr[0])
	interErr_odr = float(odrErr[1])


	slope, intercept, r_value, p_value, std_err = stats.linregress(TH_dtp_da,  TH_ncep_da_coi)
	
	#ax.text(0.7,0.32, "slope: {0:4.2f}".format(slope),transform=ax.transAxes, fontsize = 16)
	#ax.text(0.7,0.27,  "Intercept: {0:.2f}".format(intercept),transform=ax.transAxes, fontsize = 16)
	#ax.text(0.7,0.22,  "r: {0:4.2f}".format(r_value),transform=ax.transAxes, fontsize = 16)
	
	slope, intercept, r_value, p_value, std_err = stats.linregress(TH_dtp_mth,  TH_ncep_mth_coi)
	fit = slope*np.array(TH_dtp_mth) + intercept

	ax.plot(TH_dtp_mth, fit, 'k', linewidth=2, color = 'red')
	ax.text(0.7,0.12, "slope: {0:4.2f}".format(slope_odr),transform=ax.transAxes, color='red', fontsize = 16)
	ax.text(0.7,0.07,  "Intercept: {0:.2f}".format(inter_odr),transform=ax.transAxes, color='red',fontsize = 16)
	ax.text(0.7,0.02,  "r: {0:4.2f}".format(r_value),transform=ax.transAxes,color='red', fontsize = 16)

	ax.text(0.82, 0.79, "1:1",transform=ax.transAxes, color='gray', fontsize = 18)
	ax.text(0.61, 0.8,  "x1.2" ,transform=ax.transAxes, color='gray',fontsize = 18)
	ax.text(0.82, 0.58,  "x0.8",transform=ax.transAxes,color='gray', fontsize = 18)
  
	one2one = np.arange(np.min(TH_dtp_mth) - 1.0*np.min(TH_dtp_mth), np.max(TH_dtp_mth)+1.0*np.max(TH_dtp_mth), np.min(TH_dtp_mth)*0.05)
	ax.set_xlabel('Tropopause Height - Dynamic [Km]',  fontsize=14)
	ax.set_ylabel('Tropopause Height - NCEP [Km]',  fontsize=14)
	ax.tick_params(axis='both', which='major', labelsize=14)
	ax.set_ylim(np.min(TH_dtp_da)-0.3*np.min(TH_dtp_da), np.max(TH_dtp_da)+0.3*np.max(TH_dtp_da))
	ax.set_xlim(np.min(TH_dtp_da)-0.3*np.min(TH_dtp_da), np.max(TH_dtp_da)+0.3*np.max(TH_dtp_da))
	ax.plot(one2one, one2one, ls ='--', color='gray', linewidth=2)
	ax.plot(one2one, one2one*0.8, ls ='--', color='gray', linewidth=2)
	ax.plot(one2one*0.8, one2one, ls ='--', color='gray', linewidth=2) 
	ax.grid(True)
	ax.legend(prop={'size':14},loc='upper right') 

	if savePDF: 
	    pdfsav.savefig(fig2,dpi=200)
	    pdfsav.savefig(fig3,dpi=200)
	    pdfsav.savefig(fig4,dpi=200)

	    pdfsav.close() 
	else:           
		plt.show(block=	False)
		user_input = raw_input('Press any key to exit >>> ')
		sys.exit()   

if __name__ == "__main__":
    main(sys.argv[1:])
