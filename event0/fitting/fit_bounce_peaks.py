# -*- coding: utf-8 -*-
"""
Created on Sun Jan 15 20:44:56 2017

@author: mike
"""
import sys
import numpy as np
import scipy.optimize
import scipy.signal
import datetime
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pprint
import pickle

sys.path.insert(0, '/home/mike/Dropbox/0_firebird_research/' + \
'microburst_characterization/src/fitting/')
from adaptive_microburst_fitting import adaptiveMicroburstFitting

sc_id = 3
fitDict = {'sc_id':sc_id}
PLOT_RELATIVE_TIME = True

# Open data and determine which indicies to fit
fDir = '/home/mike/FIREBIRD/Datafiles/FU_' + str(sc_id) + '/hires/level2/'
fName = 'FU' + str(sc_id) + '_Hires_2015-02-02_L2.txt'

if 'fitObj' not in globals():
    fitObj = adaptiveMicroburstFitting(fDir, fName)

# Find the time indicies, and create initial conditions for the peaks in each
# energy channel.
c = np.array(['r', 'b', 'g'])
fig = plt.figure(figsize=(20, 5), dpi=80, facecolor = 'grey')
gs = gridspec.GridSpec(1,1)
dataPlt = fig.add_subplot(gs[0, 0])

for fitEnergy in range(3):
    if sc_id == 3:
        startTime = datetime.datetime(2015, 2, 2, 6, 12, 55, 0000)
        # FU3 time range is 06:12:55 to 06:12:59
        eventInd = np.where(\
        (fitObj.times >= datetime.datetime(2015, 2, 2, 6, 12, 55, 0000)) & \
        (fitObj.times <= datetime.datetime(2015, 2, 2, 6, 12, 59, 0000)))[0]
        if fitEnergy == 0:
            p0 = np.array([229, 0.8, 100E-3, 154, 1.288, 150E-3, 92, 1.93, 150E-3,\
            81, 2.5, 100E-3, 80, 3.07, 100E-3, 81, 3.34, 100E-3, 81, 3.34, 1])
        elif fitEnergy == 1:
            p0 = np.array([139, 0.8, 100E-3, 105, 1.29, 100E-3, 57, 1.89, 150E-3,\
            57, 2.49, 100E-3, 80, 3.11, 100E-3, 54, 3.31, 150E-3, 20, 3.2, 0.5])
        elif fitEnergy == 2:
            p0 = np.array([54, 0.8, 100E-3, 26, 1.34, 150E-3, 21, 1.92, 200E-3, \
            20, 2.52, 150E-3, 22, 3.1, 150E-3, 20, 3.2, 0.5])
        # Constrain paramaters. If inf, then no constain
        ubound = np.inf*np.ones(len(p0))
        lbound = -np.inf*np.ones(len(p0))
    elif sc_id == 4:
        # FU4 time range is 06:12:53 to 06:12:58
        startTime = datetime.datetime(2015, 2, 2, 6, 12, 53, 0000)
        eventInd = np.where(\
        (fitObj.times >= startTime) & \
        (fitObj.times <= datetime.datetime(2015, 2, 2, 6, 12, 57, 0000)))[0]
        
        
        # Give a good guess on the paramaters of the Gaussians.
        if fitEnergy == 0:
            p0 = np.array([225, 0.5, 100E-3, 80, 0.98, 200E-3, 30, 1.72, 200E-3, \
            25, 2.33, 200E-3, 25, 2.5, 2])
        elif fitEnergy == 1:
            p0 = np.array([100, 0.52, 60E-3, 45, 1.03, 130E-3, 20, 1.64, 150E-3, \
            10, 2.27, 150E-3, 13, 2.24, 100E-3])
        elif fitEnergy == 2:
            p0 = np.array([46, 0.52, 75E-3, 16, 1.03, 100E-3, 8.5, 1.61, 100E-3, \
            6, 2.20, 100E-3])
        # Constrain paramaters. If inf, then no constain
        ubound = np.inf*np.ones(len(p0))
        lbound = -np.inf*np.ones(len(p0))
    
    # Detrend Data
    #fitObj.hr['Col_counts'][eventInd, fitEnergy] = scipy.signal.detrend(\
    #fitObj.hr['Col_counts'][eventInd, fitEnergy], type = 'constant')
            
    # NORMALIZE TIME, where t = 0 at the left end of the event array.
    fitObj.tSec -= fitObj.tSec[eventInd[0]]
    print(fitEnergy)
    # Fit the event
    try:
        # I added 1 to sigma for the algorigm to converge. 
        popt, pcov = scipy.optimize.curve_fit(fitObj.nGaus, fitObj.tSec[eventInd], \
        fitObj.hr['Col_counts'][eventInd, fitEnergy], p0 = p0, \
        sigma =  1 + np.sqrt(fitObj.hr['Col_counts'][eventInd, fitEnergy]), \
        bounds = (lbound, ubound), method = 'lm', absolute_sigma = False)
        #print('Parameter estimates \n ', popt)
        
        # Calculate error from the covariance matrix
        perr = np.sqrt(np.diag(pcov))
        
        # Print info to screen and write to file.
#        f = open('FU' + str(sc_id) + '_CH_' + str(fitEnergy) + '_' + 'fit_params.txt', 'w')
#        print('FU', str(sc_id), ', Channel = ', str(fitEnergy), '\n')
#        f.write('FU' + str(sc_id) +  ', Channel = ' + str(fitEnergy) +'\n')
#        for i in range(len(popt)//3):
#            if popt[2 + 3*i] > 1:
#                # If the peak width is large (did not converge to the correct peak)
#                # skip it.
#                continue
#            print('Peak ', str(i), ':')
#            print('t = ', (startTime + datetime.timedelta(seconds = popt[1+i*3])).isoformat(),\
#            ' $ \pm $', "%.2f" % perr[1+3*i], 's')
#            print('A = ', popt[0 + 3*i], '$\pm $', perr[0 + 3*i], ' counts')
#            print('sigma = ', popt[2 + 3*i], '$\pm $', perr[2 + 3*i], ' s \n')
#            f.write('Peak ' + str(i) + ': \n')
#            f.write('t = ' + (startTime + datetime.timedelta(seconds = popt[1+i*3])).isoformat()+\
#            ' $ \pm $' + "%.2f" % perr[1+3*i] + 's \n')
#            f.write('A = ' + "%.2f" % popt[0 + 3*i] + ' $\pm $' + "%.2f" % perr[0 + 3*i] + ' counts \n')
#            f.write('sigma = ' + "%.2f" % popt[2 + 3*i] + '$\pm $' + "%.2f" % perr[2 + 3*i] + ' s \n')
#            
#        f.close()
        n_peaks = len(popt)//3
        
        # Write to dictionary
        fitDict['CH' + str(fitEnergy)] = {'vals':-9999*np.ones((n_peaks, 3), dtype = object), 
        'std':-9999*np.ones((n_peaks, 3))}
        
        for i in range(n_peaks):
            # Saves in format: amplitude, sigma, t0
            fitDict['CH' + str(fitEnergy)]['vals'][i, :] = np.array(
            [popt[0 + 3*i], popt[2 + 3*i], 
            (startTime + datetime.timedelta(seconds = popt[1+i*3])).isoformat()])
            fitDict['CH' + str(fitEnergy)]['std'][i, :] = np.array([
            perr[0 + 3*i], perr[2 + 3*i], perr[1+3*i]])
        
        pprint.pprint(fitDict)
        # Write the dictionary to a pickle file.
        saveDir = '/home/mike/FIREBIRD/bouncing_packet_microburst/'
        saveName = 'FU' + str(sc_id) + '_fit_dictionary.pickle'
        
        # Save picked object
        with open(saveDir + saveName, 'wb') as toFile:
            pickle.dump(fitDict, toFile, protocol=pickle.HIGHEST_PROTOCOL)
            
        # Calculate chi squared goodness of fit. 
        data = np.array(fitObj.hr['Col_counts'][eventInd, fitEnergy])
        model = np.array(fitObj.nGaus(fitObj.tSec[eventInd], *popt))
        goodnessOfFit = scipy.stats.chisquare(data, f_exp = model)[0]/(len(data) - 3*n_peaks)
        print('Goodness of fit: ', goodnessOfFit, ' NEED TO FIGURE THIS OUT!')
    except RuntimeError:
        print("Error - curve_fit failed")
    
    # PLOT RESULTS    
    if PLOT_RELATIVE_TIME:    
        dataPlt.errorbar(fitObj.tSec[eventInd], \
        fitObj.hr['Col_counts'][eventInd, fitEnergy], label = \
        'Channel ' + str(fitEnergy) + ' counts', lw = 1, fmt = 'o', markersize = 3,\
        yerr = np.sqrt(fitObj.hr['Col_counts'][eventInd, fitEnergy]), c = c[fitEnergy])
        dataPlt.plot(fitObj.tSec[eventInd], \
        fitObj.nGaus(fitObj.tSec[eventInd], *popt),\
        label = 'Channel ' + str(fitEnergy) + ' fit', lw = 2, ls = '-', \
        c = c[fitEnergy])
    else:
        dataPlt.errorbar(fitObj.times[eventInd], \
        fitObj.hr['Col_counts'][eventInd, fitEnergy], label = \
        'Channel ' + str(fitEnergy) + ' counts', lw = 1, fmt = 'o', markersize = 3,\
        yerr = np.sqrt(fitObj.hr['Col_counts'][eventInd, fitEnergy]), c = c[fitEnergy])
        dataPlt.plot(fitObj.times[eventInd], \
        fitObj.nGaus(fitObj.tSec[eventInd], *popt),\
        label = 'Channel ' + str(fitEnergy) + ' fit', lw = 2, ls = '-', \
        c = c[fitEnergy])
dataPlt.legend()
dataPlt.set_title('FU' + str(sc_id) + ' microburst decay data and fits')
dataPlt.set_xlabel('Time (UTC)')
dataPlt.set_ylabel('Counts')
dataPlt.grid(b=True, which='both', color='k')
gs.tight_layout(fig)
