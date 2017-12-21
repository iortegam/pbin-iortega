#! /usr/bin/python2.7
#----------------------------------------------------------------------------------------
# Name:
#        syncWPC.py
#
# Purpose:
#  		Sync time in Windows PC using TAB1
#
# Notes:
#   
#
#----------------------------------------------------------------------------------------

                                    #-------------------------#
                                    # Import Standard modules #
                                    #-------------------------#


                            #-------------------------#
                            # Import Standard modules #
                            #-------------------------#

import datetime as dt
from threading import Timer
import sys
import os
import itertools
import getopt
import shutil
import subprocess as sp
import logging

                            #-------------------------#
                            # Define helper functions #
                            #-------------------------#

def usage():
    ''' Prints to screen standard program usage'''
    print 'syncWPC.py'

def printTest():

	print 'Hi there, this is  a test'

def syncwTAB():
	#---------------------------------------
    # Time server failed, sync to TAB1
    # Better if this was once a day...
    #---------------------------------------
    os.system("w32tm /config /update /manualpeerlist:TAB1")
    os.system("w32tm /resync")

   							#----------------------------#
                            #                            #
                            #        --- Main---         #
                            #                            #
                            #----------------------------#

def main():

    while True:

        #---------------------------
        # Find current date and time
        #---------------------------
        crntTime = dt.datetime.utcnow()

        print 'UTC current Time: {}'.format(crntTime)
        printTest()

        #nextTime = crntTime.replace(day=crntTime.day+1, hour=1, minute=0, second=0, microsecond=0)
        nextTime = crntTime.replace(day=0, hour=0, minute=0, second=crntTime.second+5)

        print nextTime

        delta_t = nextTime - crntTime
        secs    = delta_t.seconds+1

        t = Timer(secs, printTest)
        t.start()



if __name__ == "__main__":
    main()