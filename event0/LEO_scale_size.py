# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 16:49:43 2017

@author: mike
"""

import spacepy.time
import spacepy.datamodel
import numpy as np
import datetime
import time

import lat_long_map

Re = 6371 # km
alt = 651

# Calculate decay microburst scale sizes
if 'hr3' not in globals():
    start_time = time.time()
    fDir = '/home/mike/Dropbox/0_firebird_research/microburst_characterization/data/processed/'
    fName = 'FU3_Hires_2015-02-02_L2_bounce_microburst_filtered.txt'
    hr3 = spacepy.datamodel.readJSONheadedASCII(fDir + fName)
    times3 = spacepy.time.Ticktock(hr3['Time']).UTC
    #fDir = '/home/mike/FIREBIRD/Datafiles/FU_4/hires/level2/'
    fName = 'FU4_Hires_2015-02-02_L2_bounce_microburst_filtered.txt'
    hr4 = spacepy.datamodel.readJSONheadedASCII(fDir + fName)
    times4 = spacepy.time.Ticktock(hr4['Time']).UTC
    end_time = time.time()
    print("Data loaded and times converted in %g seconds" % \
    (end_time - start_time))

# Store peak times for CH0 and CH2.
peakTimes3 = {}
peakTimes4 = {}   

peakTimes3['CH0'] = np.array(
    [datetime.datetime(2015, 2, 2, 6, 12, 53, 511493),
    datetime.datetime(2015, 2, 2, 6, 12, 55, 821201)]
)
peakTimes3['CH3'] = np.array(
    [datetime.datetime(2015, 2, 2, 6, 12, 53, 530072),
    datetime.datetime(2015, 2, 2, 6, 12, 55, 715815)]
)
peakTimes4['CH0'] = np.array(
    [datetime.datetime(2015, 2, 2, 6, 12, 53, 538517),
    datetime.datetime(2015, 2, 2, 6, 12, 55, 265608)]
)
peakTimes4['CH3'] = np.array(
    [datetime.datetime(2015, 2, 2, 6, 12, 53, 539725),
    datetime.datetime(2015, 2, 2, 6, 12, 55, 181064)]
)

# Time correction.
#peakTimes4['CH0'] = [i + datetime.timedelta(seconds = 2.28) for i in 
#peakTimes4['CH0']]
#peakTimes4['CH3'] = [i + datetime.timedelta(seconds = 2.28) for i in 
#peakTimes4['CH3']]
#peakTimes4 = np.array([datetime.datetime(2015, 2, 2, 6, 12, 55, 800000), \
#datetime.datetime(2015, 2, 2, 6, 12, 56, 300000), \
#datetime.datetime(2015, 2, 2, 6, 12, 56, 900000), \
#datetime.datetime(2015, 2, 2, 6, 12, 57, 478000)])

# Calculate the spacecraft seperation, and find the fudge factor that changes 
# the seperation to be ~19.9 km. 
idx3 = np.where(times3 >= peakTimes3['CH0'][0] - datetime.timedelta(seconds = 1.8))[0][0]
idx4 = np.where(times4 >= peakTimes4['CH0'][0])[0][0]
d = lat_long_map.distanceCalc(
        hr3['Lat'][idx3], 
        hr3['Lon'][idx3],
        hr4['Lat'][idx4], 
        hr4['Lon'][idx4], 
        alt=(hr3['Alt'][idx3] + hr4['Alt'][idx4])/2
        )    
print('Spacecraft seperation = ', '%.2f' % d, '$\pm 0.8$ km')

# Calculate dlat, dlong, and distance between P1 on FU3, and P4 on FU4.
idx4f = np.where(times4 >= peakTimes4['CH0'][-1])[0][0]

print('Latitude Scale >= ', '%.2f' % float(122.5*(hr3['Lat'][idx3] - hr4['Lat'][idx4f])), r'$\pm 0.8 \ km$')

# Time FU3 saw the microburst in 3rd energy channel.
dt3Ch2 = (peakTimes3['CH3'][-1] - peakTimes3['CH3'][0]).total_seconds()

# Statistical errors from the fits on the first and last peaks in the 1st
# and 3rd energy channels on FU3.
dt3Ch0Err = 0.02
# dt3Ch3Err = np.sqrt(0.01**2 + 0.03**2) # Assume error is in the peak
dt3Ch3Err = np.sqrt(0.1**2 + 0.08**2) # Assume "error" is coming from the peak width.

# From the seperation analysis, total error is 0.9 km
longSepErr = 0.4 # km
v408 = 10.6; v555 = 14.5; v771 = 20.1 # km/s

# Indicies of the start and end times. 
startDriftL3= np.where(times3 >= peakTimes3['CH3'][0])[0][0]
endDriftL3 = np.where(times3 >= peakTimes3['CH3'][-1])[0][0]

#Longitudinal movement of FU3 which adds to the scale size.
d3L = lat_long_map.latLongScaleCalc(hr3['Lat'][idx3], alt=alt)*\
np.abs(hr3['Lon'][startDriftL3] - hr3['Lon'][endDriftL3])

ch4TimeErr = v555*np.sqrt(0.01**2 + 0.03**2) 


print('Longitude Scale')
# print('(555 keV) >= ', '%.2f' % (d3L + dt3Ch2*v555), r'$\pm $', \
# '%.2f' % (v555*np.sqrt((dt3Ch3Err*v555)**2 + longSepErr**2)), ' km$')
# print('(771 keV) >= ', '%.2f' % (d3L + dt3Ch2*v771), r'$\pm $', \
# '%.2f' % (v771*np.sqrt((dt3Ch3Err*v771)**2 + longSepErr**2)), ' km$' )
print('(555 keV) >= ', '%.2f' % (d3L + dt3Ch2*v555), r'$\pm $', \
'%.2f' % (v555*dt3Ch3Err), ' km$')
print('(771 keV) >= ', '%.2f' % (d3L + dt3Ch2*v771), r'$\pm $', \
'%.2f' % (v771*dt3Ch3Err), ' km$' )

n = 5 # Number of peaks
def lonErr(E, L, dtb):
    """
    Calculate the longitudinal error assuming the error is from the 
    uncertanity in the bounce period.
    L = L shell
    E = Electron energy (MeV)
    dtb = error in bounce time
    """
    return 2*np.pi*n*np.cos(np.deg2rad(hr3['Lat'][idx3]))*(Re+alt)*dtb/(62.7/(L*E))

def lonScale(E, L, tb, dtb):
    """
    Calculate the longitudinal scale size. This function does not assume the spacecraft 
    motion in longitude.

    L = L shell
    E = Electron energy (MeV)
    T_Drift is the Parks approximation for alpha_0 ~ 0
    """
    lonSize = 2*np.pi*n*np.cos(np.deg2rad(hr3['Lat'][idx3]))*(Re+alt)*tb/(62.7*60/(L*E))
    lonErr = 2*np.pi*n*np.cos(np.deg2rad(hr3['Lat'][idx3]))*(Re+alt)*dtb/(62.7*60/(L*E))
    return lonSize, lonErr

L = 4.7
tb = 0.571; dtb = 0.010 # Bounce period and error.
print('Analytic drift scale size =', lonScale(0.555, L, tb, dtb), 
    'km for 555 keV electrons (not assuming any s/c motion)')
print('Analytic drift scale size =', lonScale(0.771, L, tb, dtb), 
    'km for 771 keV electrons (not assuming any s/c motion)')
#print(lonScale(0.5, L, tb), lonErr(0.5, L, dtb))
