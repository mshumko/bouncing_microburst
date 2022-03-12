# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 19:18:36 2017
# Lat-Long map of the decay microburst 
@author: mike
"""
import spacepy.datamodel
import spacepy.time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import datetime
    
class GenerateLatLonMap:
    def __init__(self, hr3, hr4, **kwargs):
        self.hr3 = hr3
        self.hr4 = hr4
        self.times3 = kwargs.get('times3', spacepy.time.Ticktock(hr3['Time']).UTC)
        self.times4 = kwargs.get('times4', spacepy.time.Ticktock(hr4['Time']).UTC)
        
    def eventIndicies(self, sc_id, startT, endT):
        """
        # creates a list of indicies for the event taken from startT and endT
        datetimes
        """
        assert (sc_id is 3 or sc_id is '3' or sc_id is 4 or sc_id is '4'), \
        'Error: incorrect spacecraft id!'
        if sc_id is 3 or sc_id is '3':
            self.FU3ind = np.where((times3 > startT) & (times3 < endT))[0]
        elif sc_id is 4 or sc_id is '4':
            self.FU4ind = np.where((times4 > startT) & (times4 < endT))[0]
        return
    
    def createBasicMap(self, ax = None):
        print(hasattr(self, 'self.FU4ind'))
#        assert hasattr(self, 'self.FU3ind') and hasattr(self, 'self.FU4ind'),
#        'Error, a time interval has not been specified. Call self.eventIndicies()'
        
        
        if ax is None:
            fig = plt.figure(figsize=(10, 10), dpi=80, facecolor = 'grey')
            gs = gridspec.GridSpec(1,1)
            self.basicMap = fig.add_subplot(gs[0, 0])
        else:
            self.basicMap = ax
        
        self.basicMap.plot(self.hr3['Lon'][self.FU3ind], \
        self.hr3['Lat'][self.FU3ind], '--', c = 'r', label = 'FU3')
        self.basicMap.plot(self.hr4['Lon'][self.FU4ind], \
        self.hr4['Lat'][self.FU4ind], '--', c = 'b', label = 'FU4')
        self.basicMap.set_xlabel('Longitude')
        self.basicMap.set_ylabel('Latitude')
        self.basicMap.legend(loc = 'best')
        return self.basicMap
        
    def offsetSCSeperation(self, lat, long, d):
        """
        This function uses trig to extract the velocity vector of
        spacecraft in Lat-Long space and then offsets it by a distance d, in
        km. This function returns a new vector of lat-long    
        """
        return
        
    def labelBasicMap(self, sc_id, times, labels, **kwargs):
        """
        Places text labels on data points specified with the times argument.
        """
        ax = kwargs.get('ax', self.basicMap)
        hOffset = kwargs.get('hOffset', 0.05)
        timeOffset = kwargs.get('timeOffset', datetime.timedelta(seconds = 0))

        if type(times) is datetime.datetime:
            print('Single value given')
            times = [times]
            labels = [labels]
            
        for i in range(len(times)):
            if sc_id is 3 or sc_id is '3':
                ind = np.where(self.times3 > times[i] + timeOffset)[0][0]
                ax.plot(self.hr3['Lon'][ind], self.hr3['Lat'][ind], 'o', c = 'r')
                ax.text(hOffset + self.hr3['Lon'][ind], self.hr3['Lat'][ind], labels[i],\
                ha='center', va='center', size = 'x-large', color = 'r')
                
            if sc_id is 4 or sc_id is '4':
                ind = np.where(self.times4 > times[i] + timeOffset)[0][0]
                ax.plot(self.hr4['Lon'][ind], self.hr4['Lat'][ind], 'o', c = 'b')
                ax.text(hOffset + self.hr4['Lon'][ind], self.hr4['Lat'][ind], labels[i],\
                ha='center', va='center', size = 'x-large', color = 'b')
        return
        
if 'hr3' not in globals():
    fDir = '/home/mike/research/firebird/Datafiles/FU_3/hires/level2/'
    fName = 'FU3_Hires_2015-02-02_L2.txt'
    hr3 = spacepy.datamodel.readJSONheadedASCII(fDir + fName)
    times3 = spacepy.time.Ticktock(hr3['Time']).UTC
    fDir = '/home/mike/research/firebird/Datafiles/FU_4/hires/level2/'
    fName = 'FU4_Hires_2015-02-02_L2.txt'
    hr4 = spacepy.datamodel.readJSONheadedASCII(fDir + fName)
    times4 = spacepy.time.Ticktock(hr4['Time']).UTC

# Define peak times and their labels.
peakTimes3 = np.array([datetime.datetime(2015, 2, 2, 6, 12, 55, 800000), \
datetime.datetime(2015, 2, 2, 6, 12, 56, 300000), \
datetime.datetime(2015, 2, 2, 6, 12, 56, 900000), \
datetime.datetime(2015, 2, 2, 6, 12, 57, 478000), \
datetime.datetime(2015, 2, 2, 6, 12, 58, 102860)])
peakLabels3 = np.array(['P1', 'P2', 'P3', 'P4', 'P5'])

peakTimes4 = np.array([datetime.datetime(2015, 2, 2, 6, 12, 53, 500000), \
datetime.datetime(2015, 2, 2, 6, 12, 54, 000000), \
datetime.datetime(2015, 2, 2, 6, 12, 54, 700000), \
datetime.datetime(2015, 2, 2, 6, 12, 55, 200000)])
peakLabels4 = np.array(['P1', 'P2', 'P3', 'P4'])


mapObj = GenerateLatLonMap(hr3, hr4, times3 = times3, times4 = times4)
starT = datetime.datetime(2015, 2, 2, 6, 12, 50, 0)
endT = datetime.datetime(2015, 2, 2, 6, 13, 00, 0)
mapObj.eventIndicies(3, starT, endT)
mapObj.eventIndicies(4, starT, endT)
mapObj.createBasicMap()


mapObj.labelBasicMap(3, peakTimes3, peakLabels3, hOffset = -0.03, \
timeOffset = -datetime.timedelta(seconds = 2.6))
mapObj.labelBasicMap(4, peakTimes4, peakLabels4, hOffset = 0.03)

plt.show()
