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
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
import datetime

plt.rc('font', size=15)


Re = 6371 # km

def distanceCalc(lat1, long1, lat2, long2, alt = 0):
    """
    This function is an implementation of the Haversine Formula to calculate
    a great circle distance between two points
    """
    dlong = np.deg2rad(long2 - long1)
    dlat = np.deg2rad(lat2- lat1)
    a = np.sin(dlat/2)**2 + np.cos(lat1)*np.cos(lat2)*np.sin(dlong/2)**2
    c = 2 * np.arctan(np.sqrt(a)/np.sqrt(1 - a))
    d = (Re + alt)*c
    return d
    
def latLongScaleCalc(lat, alt=0):
    """
    This function calculates the scale sizes of one degree of
    longitude, for a supplied latitude and altitude.    
    
    Note that the latitude scale will always be 122.5 km.
    """
    longScale = (np.pi/180)*(Re + alt)*np.cos(np.deg2rad(lat))
    return longScale
    
class LatLonMap:
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
        assert hasattr(self, 'FU3ind') and hasattr(self, 'FU4ind'), \
        'Error, a time interval has not been specified. Call self.eventIndicies()'
        
        if ax is None:
            fig = plt.figure(figsize=(10, 10), dpi=80, facecolor = 'white')
            gs = gridspec.GridSpec(1,1)
            self.basicMap = fig.add_subplot(gs[0, 0])
        else:
            self.basicMap = ax
        
        self.basicMap.plot(self.hr3['Lon'][self.FU3ind], \
        self.hr3['Lat'][self.FU3ind], '--', c = 'r', label = 'FU3')
        self.basicMap.plot(self.hr4['Lon'][self.FU4ind], \
        self.hr4['Lat'][self.FU4ind], '-', c = 'b', label = 'FU4')
        self.basicMap.set_xlabel('Longitude', fontsize = 20)
        self.basicMap.set_ylabel('Latitude', fontsize = 20)
        self.basicMap.legend(loc = 'upper left')
        gs.tight_layout(fig)
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
        vOffset = kwargs.get('vOffset', 0)
        timeOffset = kwargs.get('timeOffset', datetime.timedelta(seconds = 0))
        DRIFT_FLAG = kwargs.get('DRIFT_FLAG', True)
        DRIFT_CHANNEL = kwargs.get('DRIFT_CHANNEL', 1)
        
        assert sc_id in [3,4], 'Error: Spacecraft id is not 3 or 4.'

        if type(times) is datetime.datetime:
            print('Single value given')
            times = [times]
            labels = [labels]
            
        if DRIFT_FLAG:
            latEvolve = -9999*np.ones((len(times) ,2))
            longEvolve = -9999*np.ones((len(times), 2))
            
        if sc_id is 3:
            self.peakTimes3 = times
        elif sc_id is 4:
            self.peakTimes4 = times
        else:
            print("Error, sc_id is not correct (3 or 4)!")
            
        for i in range(len(times)):
            if sc_id is 3:
                ind = np.where(self.times3 > times[i] + timeOffset)[0][0]
                #print('FU3: ', self.hr3['Lat'][ind], self.hr3['Lon'][ind])
                ax.plot(self.hr3['Lon'][ind], self.hr3['Lat'][ind], 'o', c = 'r')
                ax.text(hOffset + self.hr3['Lon'][ind], \
                vOffset + self.hr3['Lat'][ind], labels[i],\
                ha='center', va='center', size = 'x-large', color = 'r')
                
            elif sc_id is 4:
                ind = np.where(self.times4 > times[i] + timeOffset)[0][0]
                #print('FU4: ', self.hr4['Lat'][ind], self.hr4['Lon'][ind])
                ax.plot(self.hr4['Lon'][ind], self.hr4['Lat'][ind], 'o', c = 'b')
                ax.text(hOffset + self.hr4['Lon'][ind], \
                vOffset + self.hr4['Lat'][ind], labels[i],\
                ha='center', va='center', size = 'x-large', color = 'b')
                
            if DRIFT_FLAG:
                # Calculate the new locations of the bursts.
                latEvolve[i, :], longEvolve[i, :], energy = \
                self.calcLongBounces(sc_id, DRIFT_CHANNEL,
                (times[i] - times[-1]).total_seconds(), ind)
                #print(energy)
        if DRIFT_FLAG:
            for e in range(len(energy)):
                # Print the labels
                if sc_id is 3:
                    ax.plot(longEvolve[0, e], latEvolve[-1, e], '*',
                        c = 'r', ms = 15)
                    ax.text(np.min(longEvolve[0, e]), 
                        np.max(latEvolve[-1, e]) - 0.02,
                        str(energy[e]) + ' keV', ha = 'center', 
                        rotation = 'vertical', fontsize = 15, color = 'r')
                elif sc_id is 4:
                    ax.plot(longEvolve[0, e], latEvolve[-1, e], '*',
                        c = 'b', ms = 15)
                    ax.text(np.min(longEvolve[0, e]), 
                        np.max(latEvolve[-1, e]) - 0.02,
                        str(energy[e]) + ' keV', ha = 'center', 
                        rotation = 'vertical', fontsize = 15, color = 'b')
                
        return
        
    def drawDistanceRef(self, **kwargs):
        ax = kwargs.get('ax', self.basicMap)
        scale = kwargs.get('scale', 2.28)
        orientation = kwargs.get('orientation', 'vertical')
        alt = kwargs.get('alt', 651) # altitude in km
        lat = kwargs.get('lat', 65) # degrees latitude
        scale = kwargs.get('scale', 10) # km
        
        if orientation is 'vertical':
            # Draw the vertical scale.
            ax.annotate('', xy=(15.4, 63.3), xycoords='data', \
            xytext=(15.4, 63.3 - scale/122.5), textcoords = 'data', \
            arrowprops={'arrowstyle': '<->'})
            ax.annotate('d = ' + str(scale) + ' km', xy=(15.4, 63.25), xycoords='data', \
            xytext=(5, 0), textcoords = 'offset points', size = 'x-large')
            
        elif orientation is 'horizontal':
            # Draw the vertical scale.
            ax.annotate('', xy=(15.2, 63.1), xycoords='data', \
            xytext=(15.2 + scale/latLongScaleCalc(lat, alt = alt), \
            63.1), textcoords = 'data', \
            arrowprops={'arrowstyle': '<->'})
            ax.annotate('d = ' + str(scale) + ' km', xy=(15.20, 63.11), xycoords='data', \
            xytext=(0, -5), textcoords = 'offset points', size = 'x-large')
        
        # Now format the x-axis limits to be a ratio of 1:2 to match up to the
        # correct conversion at 64 Lat.
        
        ylims = ax.get_ylim()
        xlims = ax.get_xlim()
        xLowerBound = np.mean(xlims) - \
        0.5*122.5/latLongScaleCalc(lat, alt = alt)*(ylims[1] - ylims[0])
        xUpperBound = np.mean(xlims) + \
        0.5*122.5/latLongScaleCalc(lat, alt = alt)*(ylims[1] - ylims[0])
        ax.set_xlim((xLowerBound, xUpperBound))
        
    def calcLongBounces(self, sc_id, ch, dt, dataInd, **kwargs):
        """
        NAME:    calcLongBounces(self, sc_id, ch, dt, dataInd, **kwargs)
        USE:     Calculates the longitudinal drift of electrons 
        INPUTS:  FIREBIRD spacecraft sc_id (3 or 4), and energy channel ch 
                 (1-3) for a time interval dt. dataInd is the index of the lat 
                 long to pull from HiRes to calculate the new longitude of 
                 where the electrons will be in dt. Positive dt is propegating
                 in the future, negative is in the past.
        RETURNS: lat array and long array of size 3. First index is the lower
                 bound on the energy channel, second is middle energy, and 
                 third is upper energy bound.
        MOD:     2017-01-23
        """
        assert sc_id in [3, 4], 'Error: Spacecraft id is not 3 or 4.'
        assert ch in [1,2,3, 4], 'Error: Channel is not 1, 2, 3 or 4.'
        
        if sc_id is 3:
            lat = self.hr3['Lat'][dataInd]
            long = self.hr3['Lon'][dataInd]
            alt = self.hr3['Alt'][dataInd]
            
            # if ch is 1:
            #     v = np.array([6, 6.9, 7.8])
            #     energy = np.array([231, 265, 299])
            # elif ch is 2:
            #     v = np.array([7.8, 9.2, 10.6])
            #     energy = np.array([299, 353, 408])
            # elif ch is 3:
            #     v = np.array([10.6, 12.6, 14.5])
            #     energy = np.array([408, 481, 555])
            # elif ch is 4:
            #     v = np.array([14.5, 17.3, 20.1])
            #     energy = np.array([555, 663, 771])
            if ch is 1:
                v = np.array([6, 7.8])
                energy = np.array([231, 299])
            elif ch is 2:
                v = np.array([7.8, 10.6])
                energy = np.array([299, 408])
            elif ch is 3:
                v = np.array([10.6, 14.5])
                energy = np.array([408, 555])
            elif ch is 4:
                v = np.array([14.5, 20.1])
                energy = np.array([555, 771])
        elif sc_id is 4:
            lat = self.hr4['Lat'][dataInd]
            long = self.hr4['Lon'][dataInd]
            alt = self.hr4['Alt'][dataInd]
            
            # if ch is 1:
            #     v = np.array([6.1, 6.9, 7.8])
            #     energy = np.array([220, 252, 283])
            # elif ch is 2:
            #     v = np.array([7.8, 9.2, 10.6])
            #     energy = np.array([283, 333, 384])
            # elif ch is 3:
            #     v = np.array([10.6, 12.4, 14.3])
            #     energy = np.array([384, 452, 520])
            # elif ch is 4:
            #     v = np.array([14.3, 17.1, 19.9])
            #     energy = np.array([520, 620, 721])
            if ch is 1:
                v = np.array([6.1,  7.8])
                energy = np.array([220,  283])
            elif ch is 2:
                v = np.array([7.8,  10.6])
                energy = np.array([283, 384])
            elif ch is 3:
                v = np.array([10.6, 14.3])
                energy = np.array([384, 520])
            elif ch is 4:
                v = np.array([14.3,  19.9])
                energy = np.array([520, 721])
        
        driftLong = v*dt/latLongScaleCalc(lat, alt = alt)
        return lat*np.ones(2), long + driftLong, energy
        
    def drawBoxes(self, sc_id, **kwargs):
        ax = kwargs.get('ax', self.basicMap)
        boxesFlag = kwargs.get('boxesFlag', [1, 1])
        timeOffset = kwargs.get('timeOffset', datetime.timedelta(seconds = 0))
        assert sc_id in [3,4], 'Error: sc_id not 3 or 4!'
        # Upper right corner derived from times array
        startInd3 = np.where(self.times3 > self.peakTimes3[0] + timeOffset)[0][0]
        startInd4 = np.where(self.times4 > self.peakTimes4[0])[0][0]
        
        # Last indicies
        endInd3 = np.where(self.times3 > self.peakTimes3[-1] + timeOffset)[0][0]
        endInd4 = np.where(self.times4 > self.peakTimes4[-1])[0][0]
        
        # Now calculate the drift position
        dt3 = (self.peakTimes3[-1] - self.peakTimes3[0]).total_seconds()
        dt4 = (self.peakTimes4[-1] - self.peakTimes4[0]).total_seconds()
        energyCH = 4
        lat3, lon3, E3 = mapObj.calcLongBounces(3, energyCH, -dt3, endInd3)
        lat4, lon4, E4 = mapObj.calcLongBounces(4, energyCH, -dt4, endInd4)
        
        if sc_id is 3:
            # Draw bigger box
            if boxesFlag[1]:
                p = patches.Rectangle((self.hr3['Lon'][startInd3], 
                        self.hr3['Lat'][startInd3]), 
                        lon3[1] - self.hr3['Lon'][endInd3], 
                        self.hr3['Lon'][endInd3] - self.hr3['Lon'][startInd3], 
                        fill = False, linewidth = 3, linestyle = 'dashed', 
                        color = 'r')
                ax.add_patch(p)
            # Draw smaller box
            if boxesFlag[0]:
                p = patches.Rectangle((self.hr3['Lon'][startInd3], 
                        self.hr3['Lat'][startInd3]), 
                        lon3[0] - self.hr3['Lon'][endInd3], 
                        self.hr3['Lon'][endInd3] - self.hr3['Lon'][startInd3], 
                        fill = False, linewidth = 2, linestyle = 'dotted', 
                        color = 'r')
                ax.add_patch(p)
        if sc_id is 4:
            # Draw bigger box
            if boxesFlag[1]:
                p = patches.Rectangle((self.hr4['Lon'][startInd4], 
                    self.hr4['Lat'][startInd4]), lon4[1] - 
                    self.hr4['Lon'][endInd4],
                    self.hr4['Lon'][endInd4] - self.hr4['Lon'][startInd4], 
                    fill = False, linewidth = 3, linestyle = 'solid', 
                    color = 'b')
                ax.add_patch(p)
            # Draw smaller box
            if boxesFlag[0]:
                p = patches.Rectangle((self.hr4['Lon'][startInd4], 
                    self.hr4['Lat'][startInd4]), 
                    lon4[0] - self.hr4['Lon'][endInd4], 
                    self.hr4['Lon'][endInd4] - self.hr4['Lon'][startInd4], 
                    fill = False, linewidth = 1, linestyle = 'solid', color = 'b')
                ax.add_patch(p)
        return

if __name__ == '__main__':        
    #if 'hr3' not in globals():
    fDir = '/home/mike/Dropbox/0_firebird_research/microburst_characterization/data/processed/'
    fName = 'FU3_Hires_2015-02-02_L2_bounce_microburst_filtered.txt'
    hr3 = spacepy.datamodel.readJSONheadedASCII(fDir + fName)
    times3 = spacepy.time.Ticktock(hr3['Time']).UTC
    #fDir = '/home/mike/Dropbox/0_firebird_research/microburst_characterization/data/processed/'
    fName = 'FU4_Hires_2015-02-02_L2_bounce_microburst_filtered.txt'
    hr4 = spacepy.datamodel.readJSONheadedASCII(fDir + fName)
    times4 = spacepy.time.Ticktock(hr4['Time']).UTC
    
    # Define peak times and their labels.
    peakTimes3 = np.array([datetime.datetime(2015, 2, 2, 6, 12, 53, 511493), \
    datetime.datetime(2015, 2, 2, 6, 12, 53, 988309), \
    datetime.datetime(2015, 2, 2, 6, 12, 54, 653055), \
    datetime.datetime(2015, 2, 2, 6, 12, 55, 251562), \
    datetime.datetime(2015, 2, 2, 6, 12, 55, 821201)])
    peakLabels3 = np.array(['P1', 'P2', 'P3', 'P4', 'P5'])
    
    peakTimes4 = np.array([datetime.datetime(2015, 2, 2, 6, 12, 53, 538517), \
    datetime.datetime(2015, 2, 2, 6, 12, 54, 191301), \
    datetime.datetime(2015, 2, 2, 6, 12, 54, 725912), \
    datetime.datetime(2015, 2, 2, 6, 12, 55, 265608)])
    peakLabels4 = np.array(['P1', 'P2', 'P3', 'P4'])
    
    mapObj = LatLonMap(hr3, hr4, times3 = times3, times4 = times4)
    starT = datetime.datetime(2015, 2, 2, 6, 12, 50, 0)
    endT = datetime.datetime(2015, 2, 2, 6, 13, 00, 0)
    mapObj.eventIndicies(3, starT, endT)
    mapObj.eventIndicies(4, starT, endT)
    mapObj.createBasicMap()
    hOffset3 = -0.02
    hOffset4 = 0.02
    vOffset3 = 0.01
    vOffset4 = -0.01
    timeOffset = -datetime.timedelta(seconds = 1.5)
    mapObj.basicMap.set_ylim((63.05, 63.50))
    mapObj.basicMap.set_xlim((14.2, 15.6))
    
    # For each channel, draw the channel labels
    mapObj.labelBasicMap(3, peakTimes3, peakLabels3, hOffset = hOffset3, \
    timeOffset = timeOffset, DRIFT_FLAG = True, vOffset = vOffset3)
    mapObj.labelBasicMap(3, peakTimes3, peakLabels3, hOffset = hOffset3, \
    timeOffset = timeOffset, DRIFT_FLAG = True, \
    DRIFT_CHANNEL = 2, vOffset = vOffset3)
    mapObj.labelBasicMap(3, peakTimes3, peakLabels3, hOffset = hOffset3, \
    timeOffset = timeOffset, DRIFT_FLAG = True, \
    DRIFT_CHANNEL = 3, vOffset = vOffset3)
    mapObj.labelBasicMap(3, peakTimes3, peakLabels3, hOffset = hOffset3, \
    timeOffset = timeOffset, DRIFT_FLAG = True, \
    DRIFT_CHANNEL = 4, vOffset = vOffset3)
    
    mapObj.labelBasicMap(4, peakTimes4, peakLabels4, hOffset = hOffset4, \
    DRIFT_FLAG = True, vOffset = vOffset4)
    mapObj.labelBasicMap(4, peakTimes4, peakLabels4, hOffset = hOffset4, \
    DRIFT_FLAG = True, DRIFT_CHANNEL = 2, vOffset = vOffset4)
    mapObj.labelBasicMap(4, peakTimes4, peakLabels4, hOffset = hOffset4, \
    DRIFT_FLAG = True, DRIFT_CHANNEL = 3, vOffset = vOffset4)
    mapObj.labelBasicMap(4, peakTimes4, peakLabels4, hOffset = hOffset4, \
    DRIFT_FLAG = True, DRIFT_CHANNEL = 4, vOffset = vOffset4)
    
    mapObj.drawDistanceRef(orientation = 'horizontal')
    mapObj.drawBoxes(3, timeOffset=timeOffset, boxesFlag=[0,1])
    mapObj.drawBoxes(4, boxesFlag=[0,1])
    #plt.savefig('/home/mike/Desktop/test.pdf')
    plt.show()
