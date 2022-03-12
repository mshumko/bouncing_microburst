# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
# Error analysis on the timing and spatial seperation. 

#import uncertainties as unc  
import uncertainties.unumpy as unumpy
import numpy as np
import pandas as pd
import datetime
import matplotlib.pylab as plt
import scipy.optimize

def f(t, m, b):
    return m*t + b

ccErr = 2*18.75E-3 # seconds
spatialEventDt = np.array([4.87, 4.97])
temporalEventDt = np.array([2.31, 2.23, 2.27, 2.25])

spatialEvent = unumpy.uarray(spatialEventDt, ccErr)
temporalEvent = unumpy.uarray(temporalEventDt, ccErr)

flightDt = spatialEvent.mean() - temporalEvent.mean()
D = 7.55*flightDt
print('Mean time error of all points = ', temporalEvent.mean(), 's')
print('Spacecraft Separation (clock drift not accounted for) = ', D, 'km')

dateTimes = [datetime.datetime(2015, 2, 2, 6, 11, 34, 3000), \
datetime.datetime(2015, 2, 2, 6, 11, 53, 000), \
datetime.datetime(2015, 2, 2, 6, 12, 15, 000), \
datetime.datetime(2015, 2, 2, 7, 49, 35, 000), \
datetime.datetime(2015, 2, 2, 6, 12, 21, 000), \
datetime.datetime(2015, 2, 2, 6, 12, 25, 000), \
datetime.datetime(2015, 2, 2, 7, 49, 32, 000), \
datetime.datetime(2015, 2, 2, 9, 26, 18, 000), \
datetime.datetime(2015, 2, 2, 11, 2, 30, 000)]

normDate = datetime.datetime(2015, 2, 2)
secTimes = [(dateTimes[i] - normDate).total_seconds() for i in range(len(dateTimes))]
decayPacketTime = datetime.datetime(2015, 2, 2, 6, 12, 54)
decayPacketSecTime  = (decayPacketTime - normDate).total_seconds()

lags = [4.97, 4.87, 2.31, 2.23, 2.27, 2.25, 2.31,2.2, 5.08]
types = ['s', 's', 't', 't', 't', 't', 't', 't', 's']

events = pd.DataFrame({'dateTime':dateTimes, 'secTimes':secTimes, 'lags':lags, \
'type':types})
tempStructs = events[events['type'] == 't']
spaceStructs = events[events['type'] == 's']
#xlims = [datetime.datetime(2015, 2, 2, 6, 00, 00, 3000), \
#datetime.datetime(2015, 2, 2, 10, 0, 0, 3000)]
xlims = [21000, 40000]

plt.style.use('classic')
tempStructs.plot(x = 'secTimes', y = 'lags', marker="*", \
yerr = ccErr*np.ones(len(tempStructs)), ms = 20, ls = 'None', \
xlim = xlims)

#slope, intercept, r_value, p_value, std_err = \
#stats.linregress(tempStructs['secTimes'], y = tempStructs['lags'])

popt, pcov = scipy.optimize.curve_fit(f, tempStructs['secTimes'], \
tempStructs['lags'], sigma = ccErr*np.ones(len(tempStructs)), absolute_sigma = True)

# Plot the resulting fit
tvals = np.arange(tempStructs['secTimes'].min() - 1000, \
tempStructs['secTimes'].max() + 1000)
fitLine = popt[1] + popt[0]*tvals
plt.plot(tvals, fitLine, label = 'fit')
plt.axvline(decayPacketSecTime, color = 'k', ls = '--')
plt.legend()
plt.title('FB Timing Error for 2015-02-02')
perr = np.sqrt(np.diag(pcov))
print('Fit params: slope = ', '%0.2e' % popt[0], '+/-', '%0.2e' % perr[0], \
's  intercept = ', '%0.2e' % popt[1], '+/-', '%0.2e' % perr[1], 's')


# Calculate the seperation, and it's uncertanty
fitParams = unumpy.uarray(popt, perr)
print('Time lag of decay event using fit = ', (fitParams[1] + fitParams[0]*decayPacketSecTime), 's')

Dfit = 7.55*(spatialEvent.mean() - (fitParams[1] + fitParams[0]*decayPacketSecTime))
print('Separation using fit', Dfit, 'km')

plt.show()
