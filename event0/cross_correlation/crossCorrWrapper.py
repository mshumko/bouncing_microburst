# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 14:30:13 2017

@author: mike
"""
import sys
import spacepy.datamodel
from datetime import date, datetime, timedelta
import dateutil.parser
import datetime
import os
import numpy as np

from crossCorrelationStudy import CrossCorrelateData

#sys.path.insert(0, '/home/mike/Dropbox/0_firebird_research/microburst_characterization/src/')
#from interpolateData import interpolateData

date = date(2015, 2, 18)# 1:08:06 - 1:14:31 FU4 , 1:19:34 - 1:26 for FU3
fbDir = lambda sc_id: ('/home/mike/research/firebird/Datafiles/FU_{}/'
    'hires/level2'.format(sc_id)) 
fbName = lambda sc_id: 'FU{}_Hires_{}_L2.txt'.format(sc_id, date)

# Read in the data.
hr3 = spacepy.datamodel.readJSONheadedASCII(os.path.join(fbDir(3), fbName(3)))
hr4 = spacepy.datamodel.readJSONheadedASCII(os.path.join(fbDir(4), fbName(4)))

# Convert, fix, and shift times
hr3['Time'] = np.array(list(map(dateutil.parser.parse, hr3['Time'])))
#d3['Time'] = np.array([i + timedelta(seconds=tShift) for i in d3['Time']])
hr4['Time'] = np.array(list(map(dateutil.parser.parse, hr4['Time'])))

# Fix times.
dT3 = timedelta(seconds=hr3['Count_Time_Correction'][0]
    )*np.ones(len(hr3['Time']), dtype=object)
times3 = np.array(hr3['Time']) + dT3
dT4 = timedelta(seconds=hr4['Count_Time_Correction'][0]
    )*np.ones(len(hr4['Time']), dtype=object)
times4 = np.array(hr4['Time']) + dT4

# Potential curtains on February 18th
cT3 = datetime.datetime(2015, 2, 18, 18, 3, 18)
#cT4 = datetime.datetime(2015, 2, 18, 18, 2, 47, int(5E5))
cT4 = datetime.datetime(2015, 2, 18, 18, 2, 48)
w3 = 7
w4 = 1.2

# Find the indicies for the time range
idt3 = np.where((times3 > cT3 - timedelta(seconds=w3)) & 
    (times3 < cT3 + timedelta(seconds=w3)))[0]
idt4 = np.where((times4 > cT4 - timedelta(seconds=w4)) & 
    (times4 < cT4 + timedelta(seconds=w4)))[0]

# The error between the dt of the spatial strcutures is 9.8E-2. This 
# corresponds to 740 m difference... 
channel = 0
cc = CrossCorrelateData(hr3, hr4, times3=times3, times4=times4)
cc.findCenterIndex(cT3, cT4)
cc.crossCorrData(w3, w4, data3=hr3['Sur_counts'][idt3, channel], 
    data4=hr4['Col_counts'][idt4, channel])
cc.plotResults(channel=channel, hrKey3='Sur_counts')

# Now interpolate
#interpData, interpTimes = \
#interpolateData(hr4['Col_counts'][cc.dInd4, 0], times4[cc.dInd4])
#cc.crossCorrData(w3, w4, data4 = interpData)
#cc.plotResults(interpSC = 4, interpCounts = interpData, interpTimes = interpTimes)

#        
#        # SAVE FIGURE 
#        fName = 'spatial_structure_cross_corr_FU3w_' + str(CorrTime3[ind3]*2) + 's_FU4w_'\
#        + str(CorrTime4[ind4]*2) + 's.png'
#        fDir = '/home/mike/Dropbox/0_firebird_research/microburst_characterization/plots/bouncing_packet/cross_correlation/'
#        #plt.savefig(fDir + fName, dpi = 80, facecolor = 'grey')
#        #plt.clf()
##plt.close()
