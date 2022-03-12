# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 14:57:28 2017

@author: mike
"""
# Calculate Microburst flux from the data 

import pickle
import dateutil.parser
import matplotlib.pylab as plt
import spacepy.datamodel
import numpy as np

import sys
sys.path.insert(0, '/home/mike/Dropbox/0_firebird_research/microburst_characterization/src')
from calcMicroburstFlux import CalcMicroburstFlux
from energy_spectra import fit_exponent

sc_id = 4
fDir = '/home/mike/FIREBIRD/Datafiles/FU_' + str(sc_id) + '/hires/level2/'
fName = 'FU' + str(sc_id) + '_Hires_2015-02-02_L2.txt'

if 'hr' not in globals():
    # Load HiRes data
    hr = spacepy.datamodel.readJSONheadedASCII(fDir + fName)
    times = np.array([dateutil.parser.parse(i) for i in hr['Time']])

saveDir = '/home/mike/FIREBIRD/bouncing_packet_microburst/'
with open(saveDir + 'FU' + str(sc_id) + '_fit_dictionary.pickle', 
'rb') as handle:
    fitParams = pickle.load(handle)

flux = -9999*np.ones((len(fitParams['CH0']['vals'][:, 2]), 3))
fluxErr = -9999*np.ones((len(fitParams['CH0']['vals'][:, 2]), 3))
ch = 0

for ch in range(3):
    peakTimes = fitParams['CH' + str(ch)]['vals'][:, 2]
    peakTimes = np.array([dateutil.parser.parse(i) for i in peakTimes])
    fluxObj = CalcMicroburstFlux(times, hr['Col_flux'][:, ch], peakTimes, 
                                 counts = hr['Col_counts'][:, ch])
    sigma = np.abs(np.array(list(map(float, fitParams['CH' + str(ch)]['vals'][:, 1]))))
    fluxObj.sigmaFitBounds(sigma)
    fluxObj.calculateAvgFlux(baseline_subtract = True, baselineTimeWidth = 0.5)
    flux[:, ch] = fluxObj.flux
    fluxErr[:, ch] = fluxObj.fluxErr

if sc_id is 3:
    energy_middles = np.array([265, 353, 481])
elif sc_id is 4:
    energy_middles = np.array([251, 333, 452])
    
# Now fit the exponent
E0Arr = -9999*np.ones(flux.shape[0])
J0Arr = -9999*np.ones_like(E0Arr)
E0errArr = -9999*np.ones_like(E0Arr)
J0errArr = -9999*np.ones_like(E0Arr)

for i in range(len(E0Arr)):
    fit = fit_exponent(energy_middles, flux[i, :], sigma = fluxErr[i, :])
    E0Arr[i] = fit['E0']
    J0Arr[i] = fit['J0']
    E0errArr[i] = fit['E0err']
    J0errArr[i] = fit['J0err']

c = ['r', 'b', 'g', 'c', 'k']
fitEnergy = np.linspace(energy_middles[0] - 50, energy_middles[-1] + 50, 50)

for i in range(len(E0Arr)):
    plt.errorbar(energy_middles, flux[i, :], c = c[i], yerr = fluxErr[i, :],
                 label = 'Peak ' + str(int(i)), marker='o', ls = 'None', ms = 2)
    plt.plot(fitEnergy, J0Arr[i]*np.exp(-fitEnergy/E0Arr[i]), c = c[i],
             label = 'Peak ' + str(int(i)) + ' fit')
    
plt.title('FU' + str(sc_id) + ' Bouncing Packet Microburst Flux')
plt.ylabel('Average Flux from the data')
plt.xlabel('Energy (keV)')
plt.legend(fontsize = 8)
plt.yscale('log')

# Plot evolution of E0
plt.figure(2)
plt.errorbar(range(1, len(E0Arr)+1), E0Arr, yerr = E0errArr)
plt.title('FU' + str(sc_id) + r' Bouncing Packet Evolution of $E_0$')
plt.xlabel('Peak')
plt.ylabel(r'$E_0$ (keV)')
plt.xticks(range(1, len(E0Arr)+1))