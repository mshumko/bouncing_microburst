# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 17:36:23 2017

@author: mike
"""
import uncertainties.unumpy as unumpy
import numpy as np
import dateutil.parser
import matplotlib.pylab as plt
from datetime import datetime
# Fitting time differences analysis.

def loadParameters(sc_id, energyCH, CONFUSING_PEAK = False):
    if sc_id is 3:
        if energyCH is 0:
            peakT = np.array(['2015-02-02T06:12:55.799905', \
            '2015-02-02T06:12:56.318852', '2015-02-02T06:12:56.951249', \
            '2015-02-02T06:12:57.523854', '2015-02-02T06:12:58.102860'])
            peakErr = np.array([0, 0.01, 0.01, 0.01, 0.02])
            if CONFUSING_PEAK:
                peakT = np.append(peakT, '2015-02-02T06:12:58.363441')
                peakErr = np.append(peakErr, 0.02)
        elif energyCH is 1:
            peakT = np.array(['2015-02-02T06:12:55.785175', \
            '2015-02-02T06:12:56.308691', '2015-02-02T06:12:56.934323',\
            '2015-02-02T06:12:57.505072', '2015-02-02T06:12:58.103163'])
            peakErr = np.array([0, 0.01, 0.02, 0.01, 0.03])
            if CONFUSING_PEAK:
                peakT = np.append(peakT, '2015-02-02T06:12:58.312358')
                peakErr = np.append(peakErr, 0.02)
        elif energyCH is 2:
            peakT = np.array(['2015-02-02T06:12:55.774251', \
            '2015-02-02T06:12:56.321166', '2015-02-02T06:12:56.950942',\
            '2015-02-02T06:12:57.529786', '2015-02-02T06:12:58.066881'])
            peakErr = np.array([0.01, 0.01, 0.02, 0.02, 0.03])
            # Potentially include this peak
            #'2015-02-02T06:12:58.697039' +/- 0.07
    if sc_id is 4:
        if energyCH is 0:
            peakT = np.array(['2015-02-02T06:12:53.526902', \
            '2015-02-02T06:12:54.037012', \
            '2015-02-02T06:12:54.709651', \
            '2015-02-02T06:12:55.325733'])
            peakErr = np.array([0, 0.01, 0.02, 0.03])
            
        elif energyCH is 1:
            peakT = np.array(['2015-02-02T06:12:53.523730', \
            '2015-02-02T06:12:54.025800', \
            '2015-02-02T06:12:54.644637', \
            '2015-02-02T06:12:55.232891'])
            peakErr = np.array([0, 0.01, 0.02, 0.03])
        elif energyCH is 2:
            peakT = np.array(['2015-02-02T06:12:53.524320', \
            '2015-02-02T06:12:54.029855', \
            '2015-02-02T06:12:54.550979', \
            '2015-02-02T06:12:55.115897'])
            peakErr = np.array([0, 0.01, 0.02, 0.02])
    return peakT, peakErr
    
sc_id = 4

CALCULATE_BOUNCE_PERIOD = True
CALCULATE_DISPERSION = False

if CALCULATE_BOUNCE_PERIOD:
    for energyCH in range(3):
        peakT, peakErr = loadParameters(sc_id, energyCH)
        # Normalize peak times.
        peakT = [(dateutil.parser.parse(i) - \
        datetime(2015, 2, 2, 6, 12, 50)).total_seconds() for i in peakT]
        peakTimes = unumpy.uarray(peakT, peakErr)
        
        # Ignore the last dt, since its probably not from the microburst
        dt = peakTimes[:-1] - peakTimes[1:]
        #print('dt =  \n', dt, '\n')
        print('Mean bounce period on CH = ', str(energyCH), ', FU', str(sc_id),\
        ' = ', dt.mean(), 's')

if CALCULATE_DISPERSION:
    print('Now calculating the rate dispersion.')
    peakT = {}
    peakErr = {}
    peakTimes = {}
    
    for energyCh in range(3):
        # Load in the data
        peakT[energyCh], peakErr[energyCh] = loadParameters(sc_id, energyCh)
        # Convert datetime strings to objects for analysis.
        peakT[energyCh] = [(dateutil.parser.parse(i) - \
        datetime(2015, 2, 2, 6, 12, 50)).total_seconds() for i in peakT[energyCh]]
        peakTimes[energyCh] = unumpy.uarray(peakT[energyCh], peakErr[energyCh])
        
    # Dispersion between the energy channels.
    ch1 = 0
    ch2 = 2
    values = unumpy.nominal_values(peakTimes[ch1] - peakTimes[ch2])
    errors = unumpy.std_devs(peakTimes[ch1] - peakTimes[ch2])
    plt.errorbar(range(1, len(values)+1), values, yerr = errors, \
    label = 'CH' + str(ch1) + '-' + 'CH' + str(ch2))
    plt.xlim([0, len(values)+1])
    plt.xlabel('Peak')
    plt.ylabel(r'$dt$')
    plt.legend()
    