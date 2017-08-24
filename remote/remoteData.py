#! /usr/bin/python
#----------------------------------------------------------------------------------------
# Name:
#        remoteData.py
#
# Purpose:
#       This is the data server for the remote FTIR instrument
#
#
#
# Notes:
#       1) 
#
#
# License:
#    Copyright (c) 2013-2014 NDACC/IRWG
#    This file is part of sfit4.
#
#    sfit4 is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    any later version.
#
#    sfit4 is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with sfit4.  If not, see <http://www.gnu.org/licenses/>
#
#----------------------------------------------------------------------------------------

    #-------------------------#
    # Import Standard modules #
    #-------------------------#


import os
import sys
import socket
import select
import smtplib
import datetime       as     dt
from email.mime.text  import MIMEText
from trackerUtils     import *

 
    
class FTIRdataClient(object):
    
    def __init__(self,TCP_IP="192.168.1.100",TCP_Port=5555,BufferSize=4024):
    
        #-----------------------------------
        # Configuration parameters of server
        #-----------------------------------
        self.TCP_IP          = TCP_IP
        self.TCP_Port        = int(TCP_Port)
        self.RECV_BUFFER     = int(BufferSize) 

    def setParam(self,message):
        try:
            sock = socket.socket(socket.AF_INET,socket.SOCK_STREAM)    
            sock.connect((self.TCP_IP,self.TCP_Port))
            sock.sendall("set "+message)
            
            #-------------------------
            # Loop to recieve all data
            #-------------------------
            incommingTotal = ""
            while True:
                incommingPart = sock.recv(self.RECV_BUFFER)
                if not incommingPart: break
                incommingTotal += incommingPart
            
            sock.close()    
        except:
            print "Unable to connect to data server!!"           
            incommingTotal = False   
            
        return incommingTotal
    
    def getParam(self,message):
        
        try:
            sock = socket.socket(socket.AF_INET,socket.SOCK_STREAM)    
            sock.connect((self.TCP_IP,self.TCP_Port))
            sock.sendall("get "+message)
            
            #-------------------------
            # Loop to recieve all data
            #-------------------------
            incommingTotal = ""
            while True:
                incommingPart = sock.recv(self.RECV_BUFFER)
                if not incommingPart: break
                incommingTotal += incommingPart
            
            sock.close()    
        except:
            print "Unable to connect to data server!!"           
            incommingTotal = False    
            
        return incommingTotal

    def getParamN(self,message):
        
        try:
            sock = socket.socket(socket.AF_INET,socket.SOCK_STREAM)    
            sock.connect((self.TCP_IP,self.TCP_Port))
            sock.sendall("getN "+message)
            
            #-------------------------
            # Loop to recieve all data
            #-------------------------
            incommingTotal = ""
            while True:
                incommingPart = sock.recv(self.RECV_BUFFER)
                if not incommingPart: break
                incommingTotal += incommingPart
            
            sock.close()    
        except:
            print "Unable to connect to data server!!"           
            incommingTotal = False         

        return incommingTotal

    def writeTCP(self,message):
        
        try:
            sock = socket.socket(socket.AF_INET,socket.SOCK_STREAM)    
            sock.connect((self.TCP_IP,self.TCP_Port))
            sock.sendall(message)
            
            #-------------------------
            # Loop to recieve all data
            #-------------------------
            incommingTotal = ""
            while True:
                incommingPart = sock.recv(self.RECV_BUFFER)
                if not incommingPart: break
                incommingTotal += incommingPart
            
            sock.close()    
        except:
            print "Unable to connect to data server!!"           
            incommingTotal = False          
            
        return incommingTotal

    def writeSpectra(self,message):

        try:
            sock = socket.socket(socket.AF_INET,socket.SOCK_STREAM)    
            sock.connect((self.TCP_IP,self.TCP_Port))
            sock.sendall("WRITE_OPUS "+message)
            
            #-------------------------
            # Loop to recieve all data
            #-------------------------
            incommingTotal = ""
            while True:
                incommingPart = sock.recv(self.RECV_BUFFER)
                if not incommingPart: break
                incommingTotal += incommingPart
            
            sock.close()    
        except:
            print "Unable to connect to data server!!"           
            incommingTotal = False   
            
        return incommingTotal
                

class FTIRdataServer(object):
    
    def __init__(self,ctlFvars,inputFileName):

        #-----------------------------------
        # Configuration parameters of server
        #-----------------------------------        
        self.ctlFvars        = ctlFvars
        self.TCP_IP          = ctlFvars["FTS_DataServ_IP"]
        self.TCP_Port        = int(ctlFvars["FTS_DataServ_PORT"])
        self.RECV_BUFFER     = int(ctlFvars["FTS_DATASERV_BSIZE"])
        self.connection_list = []
        self.EWSaddress      = ctlFvars["Bruker_IP"]
        self.dataParams      = {}
        self.dataParamTS     = {}
        self.baseDataDir     = ctlFvars["Dir_baseData"]
        self.inputFname      = inputFileName
        #self.pollLst         = obsList["pollLst"].strip().split(",")
        #self.hdrLst          = obsList["hdrLst"].strip().split(",")
        
    def getTS(self,crntTime):
        
        return "{0:%Y%m%d.%H%M%S}".format(crntTime)
                   
    def runServer(self):
        
        self.server_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.server_socket.bind((self.TCP_IP,self.TCP_Port))
        #self.server_socket.setsockopt(socket.IPPROTO_TCP,socket.TCP_NODELAY,1)
        self.server_socket.listen(10)         
        self.connection_list.append(self.server_socket)

        #-------------------------------------------------------------
        # Set initial parameters upon start up of database
        # TRACKERSTAT -> Initializing -- This says that the tracker is
        #    initializing and not ready for taking a measurement.
        #-------------------------------------------------------------
        crntTime = dt.datetime.utcnow()
        
        self.dataParams["TRACKER_STATUS"]       = ["INITIALIZING"]
        self.dataParamTS["TRACKER_STATUS"]      = self.getTS(crntTime)
        self.dataParams["OPUS_CMND"]            = ["STANDBY"]
        self.dataParamTS["OPUS_CMND"]           = self.getTS(crntTime) 
        self.dataParams["OPUS_LASTSCAN"]        = ["-999"]
        self.dataParamTS["OPUS_LASTSCAN"]       = self.getTS(crntTime)         
        self.dataParams["TRACKER_CMND"]         = ["-999"]
        self.dataParamTS["TRACKER_CMND"]        = self.getTS(crntTime)
        self.dataParams["TRACKER_SUNPIX"]       = ["-999"]
        self.dataParamTS["TRACKER_SUNPIX"]      = self.getTS(crntTime)
        self.dataParams["TRACKER_AZIMUTH"]      = ["-999"]
        self.dataParamTS["TRACKER_AZIMUTH"]     = self.getTS(crntTime)
        self.dataParams["TRACKER_ELEVATION"]    = ["-999"]
        self.dataParamTS["TRACKER_ELEVATION"]   = self.getTS(crntTime)        
        self.dataParams["TRACKER_SUNRAD"]       = ["-999"]
        self.dataParamTS["TRACKER_SUNRAD"]      = self.getTS(crntTime)
        self.dataParams["FILE_DEFAULTS"]        = [self.inputFname]          
        self.dataParamTS["FILE_DEFAULTS"]       = self.getTS(crntTime)
        
        #-------------------------------------
        # Start loop to listen for connections
        #-------------------------------------
        while True:

            #--------------------------
            # Get current date and time
            #--------------------------
            crntTime = dt.datetime.utcnow()
            yrstr    = "{0:04d}".format(crntTime.year)
            mnthstr  = "{0:02d}".format(crntTime.month)
            daystr   = "{0:02d}".format(crntTime.day)  
            datestr  = yrstr + mnthstr + daystr                
            
            self.crntDataDir = self.baseDataDir + datestr + "/"       
                            
            #--------------------
            # Get list of sockets
            #--------------------
            read_sockets,write_sockets,error_sockets = select.select(self.connection_list,[],[],5)
            
            for sock in read_sockets:
        
                #-----------------------
                # Handle new connections
                #-----------------------
                if sock == self.server_socket:
                    
                    #----------------------------------------------
                    # New connection recieved through server_socket
                    #----------------------------------------------
                    sockfd, addr = self.server_socket.accept()
                    
                    self.connection_list.append(sockfd)
                    print "Client (%s, %s) connected" % addr
        
                #-------------------------------------
                # Handle incomming request from client
                #-------------------------------------
                else:
                    
                    #------------------------
                    # Handle data from client
                    #------------------------
                    try:
                        data = sock.recv(self.RECV_BUFFER)
                        
                        #------------------------------------------------
                        # Three types of call to server:
                        #  1) set   -- sets the value of a data parameter
                        #  2) get   -- gets the value of a data parameter
                        #  3) write -- write data to a file
                        #------------------------------------------------
                        splitVals = data.strip().split()
                        
                        if splitVals[0].upper() == 'GET':
                            
                            #-------------------------------------------------
                            # Send value of requested parameter back to client
                            #-------------------------------------------------
                            try:
                                msg = " ".join(self.dataParams[splitVals[1].upper()])
                                sock.sendall(msg)
                            except:
                                sock.sendall("-999")
                                
                        elif splitVals[0].upper() == 'GETTS':
                            
                            #-------------------------------------------------
                            # Send value of requested parameter back to client
                            #-------------------------------------------------
                            try:
                                msg = " ".join(self.dataParamTS[splitVals[1].upper()])
                                sock.sendall(msg)
                            except:
                                sock.sendall("-999")                        
                                
                                
                        elif splitVals[0].upper() == 'GETN':
                            
                            #-------------------------------------------------
                            # Send value of requested parameter back to client
                            #-------------------------------------------------                            
                            valsToGet  = splitVals[1:]
                            valsToSend = []
                            
                            for val in valsToGet:
                                try:                           
                                    strVal = " ".join(self.dataParams[val.upper()])
                                    valsToSend.append(strVal)
                                except:
                                    valsToSend.append("-999")
                                
                            msg = ",".join(valsToSend)
                            sock.sendall(msg)                            
                            
                        elif splitVals[0].upper() == 'SET':
                            
                            #---------------------------------------------------
                            # Set value of parameter sent by client. Send back 1
                            # for success or 0 for failure
                            #---------------------------------------------------
                            self.dataParams[splitVals[1].upper()] = splitVals[2:]
                            
                            crntTime = dt.datetime.utcnow()
                            self.dataParamTS[splitVals[1].upper()] = self.getTS(crntTime)
                            
                            sock.sendall(splitVals[1:])
                            
                        elif splitVals[0].upper() == 'PING':
                            
                            #------------------------
                            # Return status to client
                            #------------------------                                    
                            sock.sendall("PONG!")                     
                            
                        elif splitVals[0].upper() == 'LISTALL':
                            
                            msgLst = []
                            
                            #----------------------------
                            # Create a string of all keys 
                            # and values to send back
                            #----------------------------
                            for k in self.dataParams:
                                msgLst.append("{0:}={1:}".format(k," ".join(self.dataParams[k])))
                            
                            msg = ";".join(msgLst)
                            
                            sock.sendall(msg)
                            
                        elif splitVals[0].upper() == 'LISTALLTS':   # List all time stamps
                            
                            msgLst = []
                            
                            #----------------------------
                            # Create a string of all keys 
                            # and values to send back
                            #----------------------------
                            for k in self.dataParamTS:
                                msgLst.append("{0:}={1:}".format(k,self.dataParamTS[k]))
                            
                            msg = ";".join(msgLst)
                            
                            sock.sendall(msg)                        
                        
                        elif splitVals[0].upper() == 'WRITE_OPUS':
                            #--------------------------------
                            # data[1] = Measurement time
                            # data[2] = Filename
                            # data[2] = SNR [RMS]
                            # data[3] = Peak Amplitude
                            # data[4] = Pre-amp gain
                            # data[5] = Signal gain
                            #--------------------------------
                            
                            #-----------------------------------------
                            # Check if Measurement summary file exists
                            #-----------------------------------------
                            fname = self.crntDataDir+"Measurement.log"
                            if not ckFile(fname):
                                with open(fname,'w') as fopen:
                                    fopen.write("{0:20s} {1:20s} {2:20s} {3:20s} {4:20s} {5:20s}\n".format("Measurement_Time","Filename","SNR_RMS","Peak_Amplitude","Pre_Amp_Gain","Signal_Gain"))
                            
                            else:
                                with open(fname,'a') as fopen:
                                    fopen.write("{0:20s} {1:20s} {2:20s} {3:20s} {4:20s} {5:20s}\n".format(splitVals[1],splitVals[2],splitVals[3],splitVals[4],splitVals[5],splitVals[6]))
                            
                            #---------------------------------------------
                            # Temprorary Store last instance in dictionary
                            #---------------------------------------------
                            crntTime                          = dt.datetime.utcnow()
                            self.dataParams["OPUS_LASTSCAN"]  = "{0:} {1:} {2:} {3:} {4:} {5:}".format(splitVals[1],splitVals[2],splitVals[3],splitVals[4],splitVals[5],splitVals[6])
                            self.dataParamTS["OPUS_LASTSCAN"] = self.getTS(crntTime)
                            #----------------------
                            # All entries are lists
                            #----------------------
                            self.dataParams["OPUS_LASTSCAN"] = self.dataParams["OPUS_LASTSCAN"].split()
                            #------------------------
                            # Return status to client
                            #------------------------                                    
                            sock.sendall("Write Successfull")
                            
                            #--------------------
                            # Get Sun information
                            #--------------------
                            sunPixInfo   = self.dataParams["TRACKER_SUNPIX"]
                            sunAzimuth   = self.dataParams["TRACKER_AZIMUTH"]
                            sunElevation = self.dataParams["TRACKER_ELEVATION"]
                            
                            #----------------------------------
                            # Send email with Measurement Stats
                            #----------------------------------
                            toemails = [onemail for onemail in self.ctlFvars["Email_to"].strip().split(",")]
                            
                            msg = MIMEText("_________________________Measurement__DATA_______________________________\n\n"                                                 + \
                                           "Filename         = {1:}\nSNR              = {2:}\nMeasurement time = {0:}\n".format(splitVals[1],splitVals[2],splitVals[3],)   + \
                                           "Peak Amplitude   = {0:}\nPre-amp Gain     = {1:}\nSignal Gain      = {2:}\n\n".format(splitVals[4],splitVals[5],splitVals[6])  + \
                                           "____________________________SOLAR__DATA__________________________________\n\n"                                                 + \
                                           "Sun Obs Time     = {0:}\nSun Intensity    = {1:}\nArea of Sun      = {2:}\n".format(sunPixInfo[0],sunPixInfo[4],sunPixInfo[5]) + \
                                           "Del_X Pixels     = {0:}\nDel_Y Pixels     = {1:}\nLock On Apt      = {2:}\n".format(sunPixInfo[2],sunPixInfo[3],sunPixInfo[1][:-3]) + \
                                           "Tracker El Angle = {0:}\nEphem El Angle   = {1:}\nTracker Az Angle = {2:}\n".format(sunElevation[1],sunElevation[2],sunAzimuth[1]) + \
                                           "Ephem Az Angle   = {0:}\n\n".format(sunAzimuth[2]) + \
                                           "___________________________Weather__DATA_________________________________\n\n"                                                 + \
                                           "Solar Sensor Nom = {0:}\nSolar Sensor Flip= {1:}\nOutside Temp     = {2:}\n".format(self.dataParams["EXTERN_SOLAR_NOMINAL"][0].strip(),self.dataParams["EXTERN_SOLAR_FLIP"][0].strip(),self.dataParams["ATM_TEMPERATURE"][0].strip()) +\
                                           "Atm Pressure     = {0:}\nRel Humidity     = {1:}\n".format(self.dataParams["ATM_BARO_PRESSURE"][0].strip(),self.dataParams["ATM_REL_HUMIDITY"][0].strip())  +\
                                           "Wind Dir (From S)= {0:}\nWind Speed       = {1:}\n".format(self.dataParams["ATM_WIND_DIR_FROM_S"][0].strip(),self.dataParams["ATM_WIND_SPEED"][0].strip()))
                            msg['Subject']= "Sucessful Measurement taken at Thule at {0:}".format(splitVals[1])
                            msg['From']   = self.ctlFvars['Email_from']
                            msg['To']     = self.ctlFvars['Email_to']            
                        
                            s = smtplib.SMTP('localhost')
                            s.sendmail(self.ctlFvars['Email_from'], toemails, msg.as_string())
                            s.quit()                                
                            
                        else:
                            pass
                        
                        #------------------------
                        # Close socket connection
                        #------------------------
                        sock.close()
                        self.connection_list.remove(sock)   
                    
                    #------------------------------------------------------
                    # Remove client from socket list if client discconnects
                    #------------------------------------------------------
                    except:
                        sock.close()
                        self.connection_list.remove(sock)
                        continue
        
        #-------------
        # Close server
        #-------------
        self.closeServer()
                
    def closeConnect(self,sock):
        sock.close()
        self.connection_list.remove(sock)        

    def closeServer(self):
        ''' Close the TCP data server '''
        self.server_socket.close()
        
        
if __name__ == "__main__":
    
    #--------------------------------
    # Retrieve command line arguments
    #--------------------------------
    try:
        ctlFile  = sys.argv[1]
    except:
        ctlFile  = "/home/tabftir/remote/ops/TAB_Defaults.input"

    #try:
        #obsListFile  = sys.argv[2]
    #except:
        #obsListFile  = "/home/tabftir/remote/ops/TABdataObsList.txt"

    #----------------------------
    # Check existance of ctl file
    #----------------------------
    ckFile(ctlFile,exitFlg=True)
    
    #-------------------------
    # Import ctlFile variables
    #-------------------------
    ctlFvars = mainInputParse(ctlFile)

    #------------------------
    # Import observation list
    #------------------------
    #obsList = mainInputParse(obsListFile)
    
    #-----------
    # Run Server
    #-----------
    server1 = FTIRdataServer(ctlFvars,ctlFile)
    server1.runServer()