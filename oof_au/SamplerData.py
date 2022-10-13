# Copyright (C) 2011 Associated Universities, Inc. Washington DC, USA.
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
# 
# Correspondence concerning GBT software should be addressed as follows:
#       GBT Operations
#       National Radio Astronomy Observatory
#       P. O. Box 2
#       Green Bank, WV 24944-0002 USA

# M&C code to get sampler data
# source /home/sparrow/integration/sparrow.tcsh

import os
import numpy as np
import glob
import sys
import pickle

import astropy.io.fits as pyfits
from astropy.time import Time as DateTime

# set some variables
X  = 0
Y1 = 1
Y2 = 2
# This is not used anywhere
MJD2DATETIME = (2400000.5 - 1721424.5)

class SamplerData:

    def __init__(self, samplerName, debug=0):

        self.debug = debug
        
        self.exprTxt = ['X','Y1', 'Y2']
        self.columnNames = []

        self.roots = ['/home/gbtlogs', '/home/archive/gbtlogs/*/*']
        dir = "/home/gbtlogs/" + samplerName
        self.SetLogDirectory(dir)
        self.ReadColumnInfo()


    def GetPlotData(self, dates, columns, expressions):
        "returns the x, y data, where x is the MJD date, and y is the data from the sampler2log file specified by the column names and their conversions"

        if self.debug:
            print(("getting data for: ", dates, columns, expressions))
       

        # set the FITS file columns whose data we will retrieve
        for col in columns:
            if col >= len(self.columnNames):
                raise IndexError("column indicies must be for: ", self.columnNames)
        if self.debug:        
            print("data for column names: ")
            for i in range(len(columns)):
                print((self.columnNames[columns[i]]))

        # set the conversions used
        for i in range(len(columns)):    
            self.exprTxt[i] = expressions[i]
            
        data     = [[], [], []]
        if len(self.columnNames) == 0:
            return data
        s = dates[0]
        e = dates[1]
        startDateTime = DateTime("{0}-{1}-{2} {3}:{4}:{5}".format(s[0],s[1],s[2],s[3],s[4],s[5]))
        endDateTime   = DateTime("{0}-{1}-{2} {3}:{4}:{5}".format(e[0],e[1],e[2],e[3],e[4],e[5]))
        startMJD      = startDateTime.mjd
        endMJD        = endDateTime.mjd

        # save these off
        self.startDateTime = startDateTime
        self.endDateTime = endDateTime
        
         
        logKeys = self.GetLogFilesInRange(startDateTime, endDateTime)

        i = 0

        numKeys = len(logKeys)
        
        for k in logKeys:
            f = self.logFiles[k]
            i += 1

            #msg = "\rReading Files: %d of %d" % (i, numKeys)
            #print msg,
            #sys.stdout.flush()

            # Note: we aren't worried about compressed files, since the
            # DSS should be reading recent sampler data.
            hdulist = pyfits.open(f)
            hdudata = hdulist[1].data
            numRows = hdudata.size
            dmjds = hdudata.field(0)
            columnData = [ [], [], [] ]
            for i, colNo in enumerate(columns):
                columnData[i].extend( hdudata.field(columns[i]) )
            # Note: we could do this better by using a query in the collection
            # of the data (see Pyfits manual)
            for rowNo in range(numRows):
                dmjd = dmjds[rowNo]
                if dmjd < startMJD or endMJD < dmjd: continue
                for id in range(len(columns)):
                    value = columnData[id][rowNo] 
                    data[id].append(value)
            hdulist.close()

        d2 = [np.array(d) for d in data]
        for d in d2:
            d.shape = (len(d),)

        return self.EvaluateExpr(d2, columns)

    def EvaluateExpr(self, data, columns):
        exprDict = globals()
        exprDict.update(np.__dict__)
        exprDict.update(locals())
        return [eval(self.BuildExpr(id), exprDict) for id in range(len(columns))]

    def BuildExpr(self, id):
        expr = self.GetExpr(id)
        for id in ("X", "Y1", "Y2"):
            expr = expr.replace(id, "data["+id+"]")
        return expr

    def GetExpr(self, id):
        return self.exprTxt[id]

    def GetMax(self, x, y):
        "Only necessary due to problems calling max mulitple times"
        # NOTE: this fails (rets 0) when called the second time in the unit tests
        # might be that np.max is getting called?
        #startIdx = max(0, startIdx - 1)
        return x if y < x else y

    def GetMin(self, x, y):
        "Only necessary due to problems calling min mulitple times"
        # NOTE: same issue as with max
        return x if x < y else y

    def GetLogFilesInRange(self, startDateTime, endDateTime):
        "Returns a list of paths to all the log files between the two times for the currently selected sampler."
        startText = startDateTime.to_datetime().strftime("%Y_%m_%d_%H:%M:%S")
        endText   = endDateTime.to_datetime().strftime("%Y_%m_%d_%H:%M:%S")
        
        length = len(self.logKeys)
        startIdx = 0
        while startIdx < length and \
              self.logKeys[startIdx] < startText:
            startIdx = startIdx + 1
        startIdx = self.GetMax(0, startIdx - 1)
        endIdx   = startIdx
        while endIdx   < length and \
              self.logKeys[endIdx] < endText:
            endIdx   = endIdx   + 1
        endIdxM1 = self.GetMin(length - 1, endIdx)
        if length - 1 < endIdx + 1:
            endIdx = length - 1
        else:
            endIdxM1 = endIdx + 1

        self.noLogsMessage = ''
        if length == 0:
            self.noLogsMessage = 'No logs'
        elif startIdx == endIdxM1 == 0:
            self.noLogsMessage = 'No logs before ' + self.logKeys[startIdx]
        elif startIdx == endIdx == (length - 1):
            self.noLogsMessage = 'No logs after ' + self.logKeys[startIdx]
        else:
            self.noLogsMessage = 'No logs between ' + self.logKeys[startIdx] + " and " + self.logKeys[endIdxM1]

        return self.logKeys[startIdx:endIdx + 1]
        
    def GetStartDateTime(self):
        return self.GetDateTime().GetStart()

    def GetEndDateTime(self):
        return self.GetDateTime().GetEnd()

    def GetDateTime(self):
        return self.dateTime
    
    def ReadColumnInfo(self):
        """
        Selects the last file in the list of all files to read the column
        headers and returns the number of columns.
        """
        
        logFile = self.GetLastLogFile()
        if logFile is None:
            return 0
        hdulist = pyfits.open(logFile)
        columns = hdulist[1].columns
        self.columnNames = columns.names
        hdulist.close()
        return len(self.columnNames) 

    def GetLastLogFile(self):
        "Returns the path to the most recent file for the currently selected sampler."
        self.logKeys, self.logFiles = self.GetAllLogFiles()
        if len(self.logKeys) > 1:
            return self.logFiles[self.logKeys[-1]]
        elif len(self.logKeys) == 1:
            return self.logFiles[self.logKeys[0]]
        else:
            self.noLogsMessage = 'No logs'
            return None

    def GetAllLogFiles(self):
        "Returns a list of paths to all the files for the currently selected sampler."
        paths = { }
        
        for dir in self.GetLogDirectories():
            files = os.listdir(dir)
            for f in files:
                paths[f] = dir + '/' + f
        keys = paths.keys()
        keys = [k for k in keys if k.find(".fits") >= 0]
        keys.sort()
        return keys, paths

    def GetLogDirectories(self):
        "Returns a list of paths to all the directories holding the currently selected sampler."
        return self.directories
    
    def FindDirectories(self, sampler):
        if self.debug:
            print(("Sampler: ", sampler))
        self.logName = sampler
        targetDirectories = [ directory + '/' + sampler
                              for directory in self.roots ]
        targetDirectories.reverse()
        self.directories = []
        for directory in targetDirectories:
            self.directories += glob.glob(directory)

    def SetLogDirectory(self, newDirectory):
        "Sets the list of directories to a single entry for the currently selected sampler."
        #self.directories = [ newDirectory ]
        self.logName = newDirectory[newDirectory.rfind("/")+1:]
        self.FindDirectories(self.logName)

    def GetLogFileName(self, logFilePath):
        "Given a full path to a log file, returns the name to the log file."
        return logFilePath[logFilePath.rfind("/")+1:]


def get_wind(sampler, time_range, cols=(0,1)):
    """
    """
    
    exprs = ('X','Y1')
    mjd, wind = sampler.GetPlotData(time_range, cols, exprs)

    return mjd, wind


def get_temp(sampler, time_range):
    """
    """
    
    cols = (0,)
    exprs = ('X',)
    mjd = sampler.GetPlotData(time_range, cols, exprs)
    temps = np.empty((24, len(mjd[0])), dtype=np.float)
    for i in range(1,25):
        cols = (0,i)
        exprs = ('X','Y1')
        _, temps[i-1] = sampler.GetPlotData(time_range, cols, exprs)
    
    return mjd, temps


def get_acc(sampler, time_range):
    """
    """
    
    conv = np.array([3.78, 0.761, 0.753]) # V / g
    
    cols = (0,)
    exprs = ('X',)
    mjd = sampler.GetPlotData(time_range, cols, exprs)
    acc = np.empty((4, len(mjd[0])), dtype=np.float)
    for i in range(1,4):
        cols = (0,i)
        exprs = ('X','Y1')
        _, acc[i-1] = sampler.GetPlotData(time_range, cols, exprs)
        acc[i-1] *= 1./conv[i-1]
        
    acc[3] = np.sqrt(np.power(acc[:3], 2).sum(axis=0))
        
    return mjd, acc


def get_qd(sampler, time_range):
    """
    """

    cols = (0,)
    exprs = ('X',)
    mjd = sampler.GetPlotData(time_range, cols, exprs)
    qdd = np.empty((4, len(mjd[0])), dtype=np.float)
    for i,c_ in enumerate([12,13,14]):
        cols = (0,c_)
        exprs = ('X','Y1')
        _, qdd[i] = sampler.GetPlotData(time_range, cols, exprs)
    
    qdd[3] = np.sqrt(np.power(qdd[:2], 2).sum(axis=0))
    
    return mjd, qdd


def dt2tuple(dt):
    """
    Converts a datetime object into a (year, month, day, hour, minute, second) tuple.
    """

    return dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second


def scanTime2tuple(time):
    """
    Converts a string representing a date and time to a tuple.

    Parameters
    ----------
    time : str
        Date and time in the format <year>-<month>-<day>T<hour>:<minutes>:<seconds>
    """

    return dt2tuple((DateTime(time, scale='utc')).datetime)
