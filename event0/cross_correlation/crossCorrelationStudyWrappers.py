# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 14:30:13 2017

@author: mike
"""
import sys
sys.path.insert(0, '/home/mike/Dropbox/0_firebird_research/microburst_characterization/src/')
from interpolateData import interpolateData

from crossCorrelationStudy import CrossCorrelateData
import spacepy.datamodel
import spacepy.time
import datetime

if 'hr3' not in globals():
###    fDir = '/home/mike/FIREBIRD/Datafiles/FU_3/hires/level2/'
###    fName = 'FU3_Hires_2015-02-02_L2.txt'
###    hr3 = spacepy.datamodel.readJSONheadedASCII(fDir + fName)
###    times3 = spacepy.time.Ticktock(hr3['Time']).UTC
###    fDir = '/home/mike/FIREBIRD/Datafiles/FU_4/hires/level2/'
###    fName = 'FU4_Hires_2015-02-02_L2.txt'
###    hr4 = spacepy.datamodel.readJSONheadedASCII(fDir + fName)
###    times4 = spacepy.time.Ticktock(hr4['Time']).UTC
    
    fDir = '/home/mike/research/firebird/Datafiles/FU_3/hires/level2/'
    fName = 'FU3_Hires_2015-02-02_L2.txt'
    hr3 = spacepy.datamodel.readJSONheadedASCII(fDir + fName)
    times3 = spacepy.time.Ticktock(hr3['Time']).UTC
    fDir = '/home/mike/research/firebird/Datafiles/FU_4/hires/level2/'
    fName = 'FU4_Hires_2015-02-02_L2.txt'
    hr4 = spacepy.datamodel.readJSONheadedASCII(fDir + fName)
    times4 = spacepy.time.Ticktock(hr4['Time']).UTC
    
# Alex's Smoking Gun spatial event. Lag = 4.97. Don't forget to interpolate!
# cT3 = datetime.datetime(2015, 2, 2, 6, 11, 34, 3000)
# cT4 = datetime.datetime(2015, 2, 2, 6, 11, 29, 000) # Adjustment for better match
# w3 = 2
# w4 = 1
# Another spatial structure that I found. Lag = 4.87
# cT3 = datetime.datetime(2015, 2, 2, 6, 11, 53, 000)
# cT4 = datetime.datetime(2015, 2, 2, 6, 11, 48, 000)
# w3 = 5
# w4 = 2
# Alex's coincident microburst center times. Lag = 2.31
# cT3 = datetime.datetime(2015, 2, 2, 6, 12, 15, 0000)
# cT4 = datetime.datetime(2015, 2, 2, 6, 12, 15, 0000)
# w3 = 15
# w4 = 5
# Alex's three coincident microburst center times. Lag = 2.23
# cT3 = datetime.datetime(2015, 2, 2, 7, 49, 35, 000)
# cT4 = datetime.datetime(2015, 2, 2, 7, 49, 32, 000)
# w3 = 4
# w4 = 1.7
# More coincident microbursts Lag = 2.27
# cT3 = datetime.datetime(2015, 2, 2, 6, 12, 21, 000)
# cT4 = datetime.datetime(2015, 2, 2, 6, 12, 18, 000)
# w3 = 4
# w4 = 1.7
    
# # More coincident microbursts. Lag = 2.25
# cT3 = datetime.datetime(2015, 2, 2, 6, 12, 25, 000)
# cT4 = datetime.datetime(2015, 2, 2, 6, 12, 22, 900)
# w3 = 1
# w4 = 2

# More simultaneous microbursts. Lag = 2.31
# cT3 = datetime.datetime(2015, 2, 2, 7, 49, 32, 000)
# cT4 = datetime.datetime(2015, 2, 2, 7, 49, 29, 00)
# w3 = 8
# w4 = 10

# More simultaneous microbursts. Lag = 2.2
cT3 = datetime.datetime(2015, 2, 2, 9, 26, 18, 000)
cT4 = datetime.datetime(2015, 2, 2, 9, 26, 14, 000)
w3 = 6
w4 = 4
    
# Spatial structure. Lag = 5.08
# cT3 = datetime.datetime(2015, 2, 2, 11, 2, 30, 000)
# cT4 = datetime.datetime(2015, 2, 2, 11, 2, 25, 000)
# w3 = 10
# w4 = 4

# Potential curtains on February 18th
# cT3 = datetime.datetime(2015, 2, 2, 11, 2, 30, 000)
# cT4 = datetime.datetime(2015, 2, 2, 11, 2, 25, 000)
# w3 = 10
# w4 = 4

# The error between the dt of the spatial strcutures is 9.8E-2. This 
# corresponds to 740 m difference... 

cc = CrossCorrelateData(hr3, hr4, times3 = times3, times4 = times4)
cc.findCenterIndex(cT3, cT4)
cc.crossCorrData(w3, w4, channel = 2)
cc.plotResults(channel = 0)

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
