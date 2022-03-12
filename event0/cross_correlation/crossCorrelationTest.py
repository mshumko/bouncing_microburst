# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 14:49:44 2017

@author: mike
"""

# Cross correlation testing scipt
import sys
sys.path.insert(0, '/home/mike/Dropbox/0_firebird_research/microburst_characterization/src/')
sys.path.insert(0, '/home/mike/Dropbox/0_firebird_research/microburst_characterization/src/fitting')
from interpolateData import interpolateData
from adaptive_microburst_fitting import nGaus
from crossCorrelationStudy import CrossCorrelateData

import spacepy.datamodel
#import spacepy.time
import datetime
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates
import dateutil.parser

# Create two data sets
dataOffsetA = datetime.timedelta(seconds = 0)
dataOffsetB = datetime.timedelta(seconds = 5)
eventOffsetA = datetime.timedelta(seconds = 0)
eventOffsetB = datetime.timedelta(seconds = 0)

cadence = 18.75E-3
windowA = 30
windowB = 10

nA = int(windowA//cadence)
nB = int(windowB//cadence)

timeA = np.array([datetime.datetime.min + dataOffsetA + datetime.timedelta(seconds = i*cadence) for i in range(nA)])
timeB = np.array([datetime.datetime.min + dataOffsetB + datetime.timedelta(seconds = i*cadence) for i in range(nB)])

timeA = spacepy.datamodel.dmarray([datetime.datetime(2015, 2, 2, 4, 34, 59, 830000) + dataOffsetA + datetime.timedelta(seconds = i*cadence) for i in range(nA)])
timeB = spacepy.datamodel.dmarray([datetime.datetime(2015, 2, 2, 4, 34, 59, 830000) + dataOffsetB + datetime.timedelta(seconds = i*cadence) for i in range(nB)])
#timeA = [timeA[i].isoformat() for i in range(len(timeA))]
#timeB = [timeB[i].isoformat() for i in range(len(timeB))]
#timeA = [spacepy.time.Ticktock(timeA[i].isoformat()).UTC[0] for i in range(len(timeA))]

# Gaussian parameters
pA = np.array([150, nA/5, 50E-3//cadence])
pB = np.array([130, nB/2, 60E-3//cadence])

countsA = nGaus(np.arange(nA), pA).reshape(nA, 1)
countsB = nGaus(np.arange(nB), pB).reshape(nB, 1)

def addBaselineAndNoise(counts, baseline):
    counts += baseline
    counts += np.random.poisson(counts)
    return counts
# Add baseline and noise 
#bA = 150
#countsA += bA + np.random.poisson(countsA + bA)
#bB = 100
#countsB += bB + np.random.poisson(countsB + bB)

countsA = addBaselineAndNoise(countsA, 150)
countsB = addBaselineAndNoise(countsB, 100)

#dates = matplotlib.dates.date2num(timeA)
#plt.plot_date(dates, countsA, '*', ls = '-')
#dates = matplotlib.dates.date2num(timeB)
#plt.plot_date(dates, countsB, '*', ls = '-')

hrA = {'Col_counts':countsA, 'Time':timeA}
hrB = {'Col_counts':countsB, 'Time':timeB}

cc = CrossCorrelateData(hrA, hrB, times3 = timeA, times4 = timeB)
cc.findCenterIndex(timeA[len(timeA)//5], timeB[len(timeB)//2])
#cc.findCenterIndex(dateutil.parser.parse(timeA[len(timeA)//2]), \
#dateutil.parser.parse(timeB[len(timeB)//2]))
cc.crossCorrData(windowA/2, windowB/2)
cc.plotResults()