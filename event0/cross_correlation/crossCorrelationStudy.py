# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 15:04:49 2017

@author: mike
"""

import sys
sys.path.insert(0, '/home/mike/Dropbox/0_firebird_research/microburst_characterization/src/')
from interpolateData import interpolateData

# Plot interesting bounce event. 
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
#from matplotlib import dates
import spacepy.datamodel
import spacepy.time
import numpy as np
import datetime

class CrossCorrelateData:
    """
    
    """
    def __init__(self, hr3, hr4, **kwargs):
        self.times3 = kwargs.get('times3', 
            spacepy.time.Ticktock(hr3['Time']).UTC)
        self.times4 = kwargs.get('times4', 
            spacepy.time.Ticktock(hr4['Time']).UTC)
        self.cadence = kwargs.get('cadence', 18.75E-3)
        self.norm = kwargs.get('norm', True)
        self.mode = kwargs.get('mode', 'valid')
        self.hr3 = hr3
        self.hr4 = hr4
        return
        
    def findCenterIndex(self, cT3, cT4):
        """
        NAME:    findCenterIndex(self, cT3, cT4)
        USE:     Finds the first time index after the input cT3 and cT4.
        INPUTS:  cT3, and cT4 which are datetime objects specifing the center
                 time of interest
        RETURNS: None, but self.centerTime3 and self.centerTime4 are placed in
                 the class namespace.
        MOD:     2017-01-17
        """
        self.cT3 = cT3
        self.cT4 = cT4
        self.centerTime3 = np.where(self.times3 >= cT3)[0][0]
        self.centerTime4 = np.where(self.times4 >= cT4)[0][0]
        return
        
    def crossCorrData(self, window3, window4, data3 = None, data4 = None, \
    channel = 0, hrKey = 'Col_counts'):
        """
        NAME:    crossCorrData(self, window3, window4, channel = 0)
        USE:     Calculates the cross correlation on FIREBIRD data adrressed 
                 with the optional hrKey keyword. 
        INPUTS:  window3 and window4 are window sizes in seconds from which to 
                 create cross-correlation arrays. Arrays centered at 
                 self.centerTime3 and self.centerTime4. Default energy channel
                 is 0 
        RETURNS: cross-correlation array, lag array, peak lag in seconds.
        MOD:     2017-01-17
        """
        self.window3 = window3
        self.window4 = window4
        autoCorrLen3 = int(window3/self.cadence)
        autoCorrLen4 = int(window4/self.cadence)  
        self.dInd3 = list(range(self.centerTime3-autoCorrLen3, \
            self.centerTime3+autoCorrLen3)) 
        self.dInd4 = list(range(self.centerTime4-autoCorrLen4, \
            self.centerTime4+autoCorrLen4))    
        if data3 is None:
            data3 =  self.hr3[hrKey][self.dInd3, channel]
        if data4 is None:
            data4 =  self.hr4[hrKey][self.dInd4, channel]
        
        # Normalize the data before cross-correlation
        if self.norm:
            data3 = np.subtract(data3, np.mean(data3))
            data4 = np.subtract(data4, np.mean(data4))
            self.crossCorr = np.correlate(data3, data4, mode=self.mode)
            self.crossCorr = np.divide(self.crossCorr, \
            np.std(data3)*np.std(data4)*(len(data3)+len(data4))/2)
        else:
            self.crossCorr = np.correlate(data3, data4, mode=self.mode)
            
        # Convert the lag arrags, with a specified cross-correlation mode.
        if self.mode == 'full':
            self.lagArr = self.cadence*np.arange(-autoCorrLen3-autoCorrLen4, \
            autoCorrLen3+autoCorrLen4-1)
        elif self.mode == 'same':
            self.lagArr = self.cadence*np.arange(-np.max([autoCorrLen3, autoCorrLen4]), \
            np.max([autoCorrLen3, autoCorrLen4]))
        elif self.mode == 'valid':
            maxS = np.max([len(data3), len(data4)])//2
            minS = np.min([len(data3), len(data4)])//2
            if len(self.crossCorr) % 2 == 0:
                self.lagArr = self.cadence*np.arange(-maxS + minS, maxS-minS)
            else:
                self.lagArr = self.cadence*np.arange(-maxS + minS, maxS-minS+1)
        # Find the peak in the cross-correlation.
        self.peakLag = self.lagArr[np.argmax(self.crossCorr)]
        return self.crossCorr, self.lagArr, self.peakLag
        
    def plotResults(self, **kwargs):
        """
        
        """
        channel = kwargs.get('channel', 0) # Plot this energy channel.
        
        # Which spacecraft data was interpolataed
        interpSC = kwargs.get('interpSC', None)
        interpCounts = kwargs.get('interpCounts', None)
        interpTimes = kwargs.get('interpTimes', None)
        hrKey3 = kwargs.get('hrKey3', 'Col_counts')
        hrKey4 = kwargs.get('hrKey4', 'Col_counts')
        
        # Set up plotting environment
        fig = plt.figure(figsize=(10, 10), dpi=80)
        gs = gridspec.GridSpec(4, 1)
        plt3 = fig.add_subplot(gs[0, 0])
        plt4 = fig.add_subplot(gs[1, 0], sharex=plt3)
        corrPlt = fig.add_subplot(gs[2, :])
        shiftedPlt = fig.add_subplot(gs[3, :])
        
        # Plot data
        if interpSC is None:
            plt3.plot(self.times3[self.dInd3], self.hr3[hrKey3][self.dInd3, channel],c = 'r')
            plt4.plot(self.times4[self.dInd4], self.hr4[hrKey4][self.dInd4, channel], c = 'g')
        elif interpSC is 3 or interpSC is '3':
            plt3.plot(interpTimes, interpCounts, '*', c = 'r')
            plt4.plot(self.times4[self.dInd4], self.hr4[hrKey3][self.dInd4, channel], c = 'g')
        elif interpSC is 4 or interpSC is '4':
            plt3.plot(self.times3[self.dInd3], self.hr3[hrKey4][self.dInd3, channel], c = 'r')
            plt4.plot(interpTimes, interpCounts, '*', c = 'g')
            
        
        plt3.set_title('a) FU3 | ' + self.times3[0].date().isoformat() + \
        '| window size = ' + str(2*self.window3) + ' s')
        plt4.set_title('b) FU4 | ' + self.times4[0].date().isoformat() + \
        '| window size = ' + str(2*self.window4) + ' s')
        plt3.set_ylabel('Counts/s')
        plt4.set_ylabel('Counts/s')
        plt3.set_xlabel('UTC')
        plt4.set_xlabel('UTC')

        
        # Plot cross-correlation
        corrPlt.plot(self.lagArr+(self.cT3-self.cT4).total_seconds(), self.crossCorr)
        corrPlt.set_title('c) Cross-Correlation')
        # corrPlt.text(0.2, 0.85, 'Total shift (FU3 lags by): {} s'.format(
        #     ((self.cT3-self.cT4).total_seconds()+self.peakLag)), ha='center',
        #     va='center', transform = corrPlt.transAxes, size = 'x-large')
        corrPlt.set_xlabel('Lag (s)')
        corrPlt.set_ylabel('Correlation coefficient')

        # Plot the time shifted results     
        if interpSC is None:
            adjTime4 = [self.times4[self.dInd4][i] + \
            datetime.timedelta(seconds = (self.cT3-self.cT4).total_seconds()+self.peakLag)\
            for i in range(len(self.times4[self.dInd4]))]
            shiftedPlt.plot(self.times3[self.dInd3], \
            self.hr3[hrKey3][self.dInd3, channel], c = 'r')
            
            shiftedPlt.plot(adjTime4, self.hr4[hrKey4][self.dInd4, channel], c = 'g')
            
        elif interpSC is 3 or interpSC is '3':         
            adjTime4 = [self.times4[self.dInd4] + datetime.timedelta(seconds = self.peakLag)\
            for i in range(len(interpTimes))]
            shiftedPlt.plot(interpTimes, interpCounts, c = 'r')
            
            shiftedPlt.plot(adjTime4, self.hr4[hrKey4][self.dInd4, channel], c = 'g')
            
        elif interpSC is 4 or interpSC is '4':
            #adjTime4 = [interpTimes[i] + datetime.timedelta(seconds = self.peakLag)\
            #for i in range(len(interpTimes))]:
            adjTime4 = [interpTimes[i] + \
            datetime.timedelta(seconds = (self.cT3-self.cT4).total_seconds()+self.peakLag)\
            for i in range(len(interpTimes))]
            shiftedPlt.plot(self.times3[self.dInd3], \
            self.hr3[hrKey3][self.dInd3, channel], c = 'r')
            
            shiftedPlt.plot(adjTime4, interpCounts, c = 'g')
            
        shiftedPlt.set_yscale('log')
        shiftedPlt.set_xlabel('UTC')
        shiftedPlt.set_title('d) FU3 shifted by {} s'.format(
            (self.cT3-self.cT4).total_seconds() + self.peakLag))
        shiftedPlt.set_ylabel('Counts/s')
        gs.tight_layout(fig)
        plt.show()
