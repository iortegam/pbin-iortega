#! /usr/bin/python2.7
import os
import tarfile
import sys
import getopt
import numpy as np

#----------------------------------------------------------------------------------------
# Name:
#        untarfiles.py
#
# Purpose:
#       1) Untar files
#
# Notes:
#       The Files will be untared in the location where the script is run
#   
#
# Version History:
#       Created, Sep, 2016  Ivan Ortega (iortega@ucar.edu)

#----------------------------------------------------------------------------------------

def usage():
    ''' Prints to screen standard program usage'''
    print 'untarfiles.py -p <path> -?'
    print '  -?             : Show all flags'
    print ' Example: untarfiles.py -p /data1/ancillary_data/fl0/FRAPPE/FRAPPE_FLEXPART_groundsites_backtrajectories/'

def main(argv):
	#--------------------------------
    # Retrieve command line arguments
    #--------------------------------
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'p:?')

    except getopt.GetoptError as err:
        print str(err)
        usage()
        sys.exit()

    #-----------------------------
    # Parse command line arguments
    #-----------------------------
    for opt, arg in opts:
        # Check input file flag and path
        if opt == '-p':
            path = arg
        else:
            print 'Unhandled option: ' + opt
            sys.exit()

    if not( path.endswith('/') ):
            path = path + '/'    
    
    listF = []

    for (dirpath, dirnames, filenames) in os.walk(path):
    	listF.extend(filenames)

    listF = np.asarray(listF)

    for fname in listF:

    	if (fname.endswith("tar.gz")):
    	    tar = tarfile.open(path+fname, "r:gz")
    	    tar.extractall()
    	    tar.close()
    	elif (fname.endswith("tar")):
    	    tar = tarfile.open(path+fname, "r:")
    	    tar.extractall()
    	    tar.close()

 #            tar.close()

if __name__ == "__main__":
    main(sys.argv[1:])