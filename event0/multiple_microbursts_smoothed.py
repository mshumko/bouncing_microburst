# -*- coding: utf-8 -*-
"""
Created on Mon Jan  2 15:56:54 2017

This script compares the FIREBIRD HiRes data of the bouncing packet event. Also adjusts for the timing difference in the two clocks with the "secOffset" variable.

Plots FU3 col and/or sur flux and FU4 col flux. FU3 surface flux can to toggled by "plotSur" variable

@author: mike
"""

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

# Define script parameters
secOffset = -2.28
plotSur = False
COl_G = 9 # cm^2 sr
SUR_G = 23 # cm^2 sr
cadence = 18.75E-3 # seconds
plotErr = True
downsampleErr = 10
roll_avg = 8


# Alex's detector calibrations
FU3ColWidths = np.array([68.7, 107.9, 147.2, 215.9, 284.5])
FU3ColBounds = np.array([231.0, 299.7, 407.6, 554.8, 770.7, 1055.2, 2330])
FU3SurWidths = np.array([ 64.5, 101.3, 138.1, 202.7, 267.1])
FU3SurBounds = np.array([176.4, 240.9, 342.2, 480.3, 683.0, 950.1, 2147.3])
FU4ColWidths = np.array([63.7 , 100.2 , 136.7 , 200.4 , 264.25])
FU4ColBounds = np.array([219.7, 283.4, 383.6, 520.3, 720.7, 984.95, 2169.3])

if 'fitObj3' not in globals():
    print('Reading data')
    fDir = '/home/mike/FIREBIRD/Datafiles/FU_3/hires/level2/'
    fName = 'FU3_Hires_2015-02-02_L2.txt'
    fitObj3 = adaptiveMicroburstFitting(fDir, fName)
    
    fDir = '/home/mike/FIREBIRD/Datafiles/FU_4/hires/level2/'
    fName = 'FU4_Hires_2015-02-02_L2.txt'
    fitObj4 = adaptiveMicroburstFitting(fDir, fName)

# Set time ranges.
startTime = datetime.datetime(2015, 2, 2, 7, 49, 10, 0)
endTime = datetime.datetime(2015, 2, 2, 7, 49, 50, 0)

bounceInd3 = np.where((fitObj3.times > startTime) 
    & (fitObj3.times < endTime))[0]
bounceInd4 = np.where((fitObj4.times > startTime) 
    & (fitObj4.times < endTime))[0]

assert (len(bounceInd3) > 0 and len(bounceInd4) > 0), 'Error, no data found'+\
' in the time range specified!'

################ PLOTTING STUFF HERE ###############
plt.rcParams.update({'font.size': 15})
fig = plt.figure(figsize=(20, 10), dpi=80, facecolor = 'white')

if plotSur:
    gs = gridspec.GridSpec(3, 1)
else:
    gs = gridspec.GridSpec(2, 1)

if plotSur:
    surPlt3 = fig.add_subplot(gs[0, 0])
    colPlt3 = fig.add_subplot(gs[1, 0], sharex = surPlt3)
    colPlt4 = fig.add_subplot(gs[2, 0], sharex = colPlt3)
else:
    colPlt3 = fig.add_subplot(gs[0, 0])
    colPlt4 = fig.add_subplot(gs[1, 0], sharex = colPlt3)

colFlux3 = lambda i : fitObj3.hr['Col_counts'][bounceInd3, i]/(
    COl_G*FU3ColWidths[i]*cadence)
colFlux4 = lambda i : fitObj4.hr['Col_counts'][bounceInd4, i]/(
    COl_G*FU4ColWidths[i]*cadence)

colFlux3err = lambda i : np.sqrt(fitObj3.hr['Col_counts'][bounceInd3, i])/(
    COl_G*FU3ColWidths[i]*cadence)
colFlux4err = lambda i : np.sqrt(fitObj4.hr['Col_counts'][bounceInd4, i])/(
    COl_G*FU4ColWidths[i]*cadence)

def smoothedFlux3(i, roll_avg):
    flux = colFlux3(i)
    return np.convolve(np.ones(roll_avg)/roll_avg, flux, mode = 'same')

def smoothedFlux4(i, roll_avg):
    flux = colFlux4(i)
    return np.convolve(np.ones(roll_avg)/roll_avg, flux, mode = 'same')

c = ['r', 'b', 'g', 'c', 'k']

# PLOT flux
for i in range(5):
    # Plot FU3 sur and col Flux (convert here as well)
    colPlt3.plot(fitObj3.times[bounceInd3] +
    datetime.timedelta(seconds = secOffset),
    smoothedFlux3(i, roll_avg), c = c[i],
    label = ('%.0f' % FU3ColBounds[i] + ' < E < ' + 
    '%.0f' % FU3ColBounds[i+1]))

    if plotSur:   
        surPlt3.plot(fitObj3.times[bounceInd3] +
        datetime.timedelta(seconds = secOffset),
        fitObj3.hr['Sur_counts'][bounceInd3, i]/(
        SUR_G*FU3SurWidths[i]*cadence), c = c[i],
        label = ('%.0f' % FU3SurBounds[i] + ' < E < ' + 
        '%.0f' % FU3SurBounds[i+1]))

    # Plot FU4 col Flux
    colPlt4.plot(fitObj4.times[bounceInd4],
    smoothedFlux4(i, roll_avg), label = ('%.0f' %
    FU4ColBounds[i] + ' < E < ' + '%.0f' % 
    FU4ColBounds[i+1]), c = c[i])

    if plotErr:
        # Plot FU3 sur and col Flux (convert here as well)
        colPlt3.errorbar(fitObj3.times[bounceInd3][::downsampleErr] +
        datetime.timedelta(seconds = secOffset),
        smoothedFlux3(i, roll_avg)[::downsampleErr], fmt = 'none',
        yerr = colFlux3err(i)[::downsampleErr],
        ecolor = c[i])

        # Plot FU4 col Flux
        colPlt4.errorbar(fitObj4.times[bounceInd4][::downsampleErr],
        smoothedFlux4(i, roll_avg)[::downsampleErr], fmt = 'none', yerr = colFlux4err(i)[::downsampleErr],
        ecolor = c[i])

colPlt3.set_yscale('log')
colPlt4.set_yscale('log')
if plotSur:
    surPlt3.set_yscale('log')

# Format the time stamps 
hfmt = dates.DateFormatter('%H:%M:%S')
colPlt4.xaxis.set_major_formatter(hfmt)

# Legend
colPlt3.legend(loc = 'upper left', fontsize = 12)
colPlt4.legend(loc = 'upper left', fontsize = 12)
if plotSur:
    surPlt3.legend(loc = 'upper left', fontsize = 12)

# Set the title, and axies labels.
if plotSur:
    surPlt3.set(title = 'FU3', ylabel = r'Flux $(s \ cm^2 \ sr \ keV)^{-1}$')
    colPlt3.set(title = 'b) FU3', ylabel = r'Flux $(s \ cm^2 \ sr \ keV)^{-1}$')
    colPlt4.set(title = 'c) FU4', ylabel = r'Flux $(s \ cm^2 \ sr \ keV)^{-1}$', xlabel = 'UTC')
else:
    colPlt3.set(title = 'FU3', ylabel = r'Flux $(s \ cm^2 \ sr \ keV)^{-1}$')
    colPlt4.set(title = 'FU4', ylabel = r'Flux $(s \ cm^2 \ sr \ keV)^{-1}$', xlabel = 'UTC')

# Set plot x and y axis limits.
colPlt3.set_ylim([0.01, 100])
colPlt4.set_ylim([0.01, 100])
#if plotSur:
#    surPlt3.set_ylim([0, 30])
#colPlt3.set_xlim([datetime.datetime(2015, 2, 2, 6, 
#12, 50, 0), datetime.datetime(2015, 2, 2, 6, 12, 
#56, 500000)])

# Grid settings
colPlt3.minorticks_on()
colPlt4.minorticks_on()
if plotSur:
    surPlt3.minorticks_on()

colPlt3.grid(True, linewidth=2)
colPlt4.grid(True, linewidth=2)

# Turn off tick labels for all but buttom plot.
plt.setp(colPlt3.get_xticklabels(), visible=False)
if plotSur:
    plt.setp(surPlt3.get_xticklabels(), visible=False)
gs.tight_layout(fig)
plt.show()
