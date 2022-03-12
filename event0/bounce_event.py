# -*- coding: utf-8 -*-
"""
Created on Mon Jan  2 13:41:37 2017

@author: mike
"""

# Plot interesting bounce event. 
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import dates
#import spacepy.datamodel
#import spacepy.time
import numpy as np
import datetime

import sys
sys.path.insert(0, '/home/mike/Dropbox/0_firebird_research/microburst_characterization/src/fitting/')
from adaptive_microburst_fitting import adaptiveMicroburstFitting

sys.path.insert(0, '/home/mike/Dropbox/0_firebird_research/microburst_characterization/experimental/')
from microburst_fit_and_bounce_period import bouncePeriodCalculator

fDir = '/home/mike/FIREBIRD/Datafiles/FU_4/hires/level2/'
fName = 'FU4_Hires_2015-02-02_L2.txt'
fitObj = adaptiveMicroburstFitting(fDir, fName)

# Set up plotting environment
fig = plt.figure()
gs = gridspec.GridSpec(5, 1)

data = fig.add_subplot(gs[:, 0])
loc = data.twinx() 

#data = fig.add_subplot(gs[0:3, 0])
#loc = fig.add_subplot(gs[3:5, 0], sharex = data)

# make these tick labels invisible
#plt.setp(data.get_xticklabels(), visible=False)

#Isolate the indicies of the interesting event.
#For FU4
bounceInd = np.where((fitObj.times > datetime.datetime(2015, 2, 2, 6, 12, 30, 0)) 
    & (fitObj.times < datetime.datetime(2015, 2, 2, 6, 12, 58, 0)))[0]
    
#For FU3
#bounceInd = np.where((fitObj.times > datetime.datetime(2015, 2, 2, 6, 12, 55, 0)) 
#    & (fitObj.times < datetime.datetime(2015, 2, 2, 6, 13, 2, 0)))[0]
    
## Coincident microburst times on FU3:
#bounceInd = np.where((fitObj.times > datetime.datetime(2015, 2, 2, 7, 49, 33, 000)) 
#    & (fitObj.times < datetime.datetime(2015, 2, 2, 7, 49, 36, 000)))[0]
#    
## Coincident microburst times on FU4:
#bounceInd = np.where((fitObj.times > datetime.datetime(2015, 2, 2, 7, 49, 30, 000)) 
#    & (fitObj.times < datetime.datetime(2015, 2, 2, 7, 49, 33, 000)))[0]
    
# Now calculate the bounce period for 250 keV and 1 MeV electrons.
calcIndex = bounceInd[len(bounceInd)//2]
LLA = {'dateTime':[fitObj.times[calcIndex].isoformat()], \
    'x1':[fitObj.hr['Alt'][calcIndex]], \
    'x2':[fitObj.hr['Lat'][calcIndex]], \
    'x3':[fitObj.hr['Lon'][calcIndex]]}
kp = [fitObj.hr['kp'][calcIndex]]
            
data.text(0.5, 0.9, r'$T_{b ( 250 keV)} = $' + \
    '%.2f' % bouncePeriodCalculator(LLA, kp, 250, L = fitObj.hr['McllwainL'][calcIndex]) + ' ' + r'$T_{b (500 keV)} = $'\
    + '%.2f' % bouncePeriodCalculator(LLA, kp, 500, L = fitObj.hr['McllwainL'][calcIndex]), ha='center', va='center', \
    transform = data.transAxes, size = 'x-large')

# Plot HiRes data
for i in range(6):
    data.plot(fitObj.times[bounceInd], fitObj.hr['Col_counts'][bounceInd, i], \
        label = 'Ch ' + str(int(i+1)))
hfmt = dates.DateFormatter('%H:%M:%S')

data.xaxis.set_major_locator(dates.SecondLocator())
data.xaxis.set_major_formatter(hfmt)

data.legend(loc = 'best')
data.set_title('FU4 Col HiRes from 2015-02-02')
data.set_ylabel('Counts')
data.minorticks_on()
data.grid(b=True, which='both', color='k')#, linestyle='--')
#data.set_xticklabels(data.xaxis.get_majorticklabels(), rotation=45)

# Plot L, MLT
loc.plot(fitObj.times[bounceInd], fitObj.hr['McllwainL'][bounceInd], \
    label = 'T89 Mcllwain L', lw = 3, ls = '--', c = 'r')
loc.plot(fitObj.times[bounceInd], fitObj.hr['MLT'][bounceInd], \
    label = 'T89 MLT', lw = 3, ls = '--', c = 'b')
loc.set_ylabel('T89 Mcllwain L (Red), T89 MLT (blue)')
    
#loc.legend(loc = 'best')
gs.tight_layout(fig)