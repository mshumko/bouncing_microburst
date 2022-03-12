# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 10:04:13 2017

@author: mike
"""
import numpy as np

# Decay event timing error and seperation, derived from Arlo's work.

# Seperation analysis. This takes the cross-spacecraft seperation calculated by
# propegating TLEs back from TLEDate arrays to the decay microburst event with
# STK, SGP-4 Algorithm.

TLEDate3 = ['datetime.datetime(2015, 2, 4, 11, 59)', \
'datetime.datetime(2015, 2, 5, 2, 21)','datetime.datetime(2015, 2, 5, 11, 55)',\
'datetime.datetime(2015, 2, 5, 23, 6)', 'datetime.datetime(2015, 2, 6, 11, 52)',\
'datetime.datetime(2015, 2, 7, 11, 48)', 'datetime.datetime(2015, 2, 7, 19, 47)']
TLEDate4 = ['datetime.datetime(2015, 2, 4, 10, 23)', \
'datetime.datetime(2015, 2, 5, 13, 34)','datetime.datetime(2015, 2, 5, 11, 55)',\
'datetime.datetime(2015, 2, 5, 23, 6)', 'datetime.datetime(2015, 2, 6, 11, 52)'\
'datetime.datetime(2015, 2, 7, 11, 48)', 'datetime.datetime(2015, 2, 7, 19, 47)']
seperation = np.array([21.20, 17.15, 17.23, 17.54, 17.80, 17.68, 20.06])
propegateType = ['b', 'b', 'b', 'b', 'b', 'b', 'b']

print('Mean seperation = ', np.mean(seperation), '+/-', np.std(seperation))

# Timing analysis
center = np.array([-19.55, -17.10])
lowerTotalErr = np.array([-19.47, -16.91])
upperTotalErr = np.array([-19.58, -17.50])
centerTime = np.abs(center[0] - center[1])

relativeLowerErr = -center + lowerTotalErr
relativeUpperErr = -center + upperTotalErr

lowerError =  centerTime*np.sqrt(relativeLowerErr[0]**2 + relativeLowerErr[1]**2)
upperError =   centerTime*np.sqrt(relativeUpperErr[0]**2 + relativeUpperErr[1]**2)

print('Timing Error = ', centerTime, '^{+', "%0.2f" % upperError, \
'}_-{', "%0.2f" % lowerError, '}')