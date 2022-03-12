# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 12:59:55 2017

@author: mike
"""
import scipy.integrate as integrate
import scipy.stats as stats
import scipy.optimize
import scipy.special as special
import numpy as np
import matplotlib.pylab as plt
from matplotlib.ticker import MaxNLocator

CADENCE = 18.75E-3
G_FACTOR = 9

# The J bar, or average flux integral (it is surprisingly not analytic!)
fGaus = lambda t, t0, sigma: np.exp(-np.power(t - t0, 2)/(2*np.power(sigma, 2)))

def jBarSigmaNorm(j0, t0, sigma):
    val = j0*np.sqrt(np.pi/2)*special.erf(np.sqrt(1/2))
    return val

def jSigmaNorm(j0, t0, sigma):
    val = np.sqrt(2)*sigma*j0*np.sqrt(np.pi/2)*special.erf(np.sqrt(1/2))
    return val
    
def fit_exponent(xData, yData, sigma = None, absolute_sigma = True):
    """
    NAME: fit_exponent(self, xData, yData)
    USE:  Fits an exponent to xData and yData using linear regression after
          the exponent has been transformed into a linear function.
    RETURNS: A dictionary of E0 and J0 for that fit. NEED TO INCLUDE ERROR
          OUTPUT!!!!
    MOD: 2016-10-21
    AUTHOR: Mykhaylo Shumko
    """
    
    wData = np.log(yData)
    popt, pcov = scipy.optimize.curve_fit(lambda x,m,b: m*x + b, xData, wData,
    sigma = sigma, absolute_sigma = absolute_sigma)
    slope = popt[0]
    intercept = popt[1]
    fitParams = {'J0': np.exp(intercept), 'E0':(-1/slope), 
    'J0err':np.sqrt(np.diag(pcov))[1], 'E0err':np.sqrt(np.diag(pcov))[0]/slope**2}
    return fitParams
    
if __name__ == '__main__':
    # Energy spectra
    # Need to organize by spacecraft, energy, peak number, paramaters
    # Dictionary keys are spacecraft and within it is a dictionary of  energy
    # channels. The value in the dict will be a 2D array of peak, and values
    # j0, t0, sigma from fit data. 
    j = {}
    j['FU3'] = {}
    j['FU4'] = {}
    j['FU3']['CH1'] = np.array([[186.26, 55.799, 0.1], [112.08, 56.319, 0.13],
                                [43.05, 56.951, 0.13], [31.37, 57.524, 0.07],
                                [15.25, 58.103, 0.04]])
    j['FU3']['CH2'] = np.array([[107.54, 55.785, 0.11], [73.28, 56.309, 0.12],
                                [27.30, 56.934, 0.16], [29.89, 57.505, 0.08],
                                [13.45, 58.103, 0.07]])
    j['FU3']['CH3'] = np.array([[45.01, 55.774, 0.1], [24.12, 56.321, 0.17],
                                [19.75, 56.951, 0.18], [17.48, 57.530, 0.17],
                                [18.61, 58.067, 0.15]])
    j['FU3']['energyMiddle'] = np.array([265, 353, 481])
    j['FU3']['energyWidth'] = np.array([69, 108, 216])
    
    j['FU4']['CH1'] = np.array([[160.26, 53.527, 0.07], [56.77, 54.037, 0.15],
                                [17.58, 54.710, 0.16], [10.68, 55.327, 0.14]])
    j['FU4']['CH2'] = np.array([[87.94, 53.524, 0.07], [36.02, 54.024, 0.15],
                                [15.08, 54.636, 0.17], [8.56, 55.233, 0.13]])
    j['FU4']['CH3'] = np.array([[30.57, 53.524, 0.08], [14.93, 54.029, 0.1],
                                [6.83, 54.551, 0.13], [4.32, 55.116, 0.13]])                            
    j['FU4']['energyMiddle'] = np.array([251, 333, 452])
    j['FU4']['energyWidth'] = np.array([64, 100, 137])
                                
                                
    # Calculate average flux.
    jBar3 = -9999*np.ones((5, 3))
    jBar4 = -9999*np.ones((4, 3))
    E0Arr3 = -9999*np.ones(5)
    E0errArr3 = -9999*np.ones(5)
    E0Arr4 = -9999*np.ones(4)
    E0errArr4 = -9999*np.ones(4)
    
    for ch in range(3):
        for peak in range(5):
            jBar3[peak, ch] = jBarSigmaNorm(*j['FU3']['CH' + str(ch+1)][peak])
        for peak in range(4):
            jBar4[peak, ch] = jBarSigmaNorm(*j['FU4']['CH' + str(ch+1)][peak])
    
    # Do the fits.
    for i in range(5):
        fitParams3 = fit_exponent(j['FU3']['energyMiddle'], 
                                  np.divide(jBar3[i, :], 
                                j['FU3']['energyWidth']*CADENCE*G_FACTOR),
                                sigma = np.divide(np.sqrt(jBar3[i, :]), 
                                j['FU3']['energyWidth']*CADENCE*G_FACTOR))
        E0Arr3[i] = fitParams3['E0']
        E0errArr3[i] = fitParams3['E0err']
        
        if i == 4: # Last peak on FU3 that FU4 did not see
            continue
        fitParams4 = fit_exponent(j['FU4']['energyMiddle'], 
                                  np.divide(jBar4[i, :], 
                                j['FU4']['energyWidth']*CADENCE*G_FACTOR), 
                                sigma = np.divide(np.sqrt(jBar4[i, :]), 
                                j['FU4']['energyWidth']*CADENCE*G_FACTOR))
        
        E0Arr4[i] = fitParams4['E0']
        E0errArr4[i] = fitParams4['E0err']
    
    # Plot flux
    c = ['r', 'b', 'g', 'k', 'm']
    plt.scatter(j['FU3']['energyMiddle'], np.divide(jBar3[-1, :], j['FU3']['energyWidth']*CADENCE*G_FACTOR)
    , c = c[-1], s = 50)
    for i in range(4):
        plt.scatter(j['FU3']['energyMiddle'], np.divide(jBar3[i, :], 
                    j['FU3']['energyWidth']*CADENCE*G_FACTOR)
    , c = c[i], s = 50)
        plt.scatter(j['FU4']['energyMiddle'], 
        np.divide(jBar4[i, :], j['FU4']['energyWidth']*CADENCE*G_FACTOR), c = c[i], s = 50)
        
    plt.ylabel(r'Flux $(s \ cm^2 \ sr \ keV)^{-1}$')
    plt.title(r'Bouncing packet integrated $\sigma$ flux')
    plt.xlabel('Energy (keV)')
    plt.yscale('log')
    
    # Plot evolution of E0
    plt.figure(2)
    ax = plt.subplot(111)
    ya = ax.axes.get_yaxis()
    #ax = plt.figure().gca()
    ya.set_major_locator(MaxNLocator(integer=True))
    #plt.figure(2)
    ax.errorbar(range(1, 6), E0Arr3[:], yerr = E0errArr3[:], lw = 2, capsize=10, label = 'FU3 spectra')
    ax.errorbar(range(1, 5), E0Arr4, yerr = E0errArr4, lw = 2, capsize=10, label = 'FU4 spectra')
    ax.set_xlim((0.5, 5.5))
    ax.set_ylabel(r'$E_0$')
    ax.set_xlabel('Peak')
    ax.set_title(r'Evolution of $E_0$ of the bouncing packet')
    ax.legend(loc = 'best')