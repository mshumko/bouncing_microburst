import sys, pickle, pprint
import numpy as np
import scipy.optimize, scipy.signal, scipy.stats
from datetime import datetime
from datetime import timedelta
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from uncertainties import ufloat
import spacepy.datamodel
import dateutil.parser

# Some fitting parameters
sc_id = 4
smooth_plt = True
#energy = 2
fitDict = {'sc_id':sc_id}
PLOT_RELATIVE_TIME = True
startTime = datetime(2015, 2, 2, 6, 12, 53, 0000)
endTime = datetime(2015, 2, 2, 6, 12, 58, 0000)
fDir = '/home/mike/Dropbox/0_firebird_research/microburst_characterization/data/processed/'
fName = 'FU' + str(sc_id) + '_Hires_2015-02-02_L2_bounce_microburst_filtered.txt'

def nGaus(t, *args):
    """
    NAME:    nGaus(self, t, *args)
    USE:     N-peaked Gaussian function
    INPUTS:  independent variable t, and *args are the Gaussian parameters
             with the format p11, p12, p13, p21, p22, p23... where the 
             first indicie is the peak number, the second indicie is the 
             amplitude, peak center, and sigma for each Gaussian.
    RETURNS: Numpy array of amplitude values.
    MOD:     2016-12-22
    """
    # If args are passed as an array.
    if len(args) == 1:
        args = args[0]
    assert len(args) % 3 == 0, 'Invalid number of arguments. ' + \
                                'Need 3 arguments for each Gaussian.'
    val = 0.0
    args = np.array(args)
    
    #  This for loop sums over all of the gaussians.
    for i in range(0,len(args),3):
        x = np.divide(np.power((t-args[i+1]), 2),(2*np.power(args[i+2], 2)))
        val += args[i]*np.exp(-x.astype(float))
    return val
    
    
def redchisqg(ydata,ymod,deg=2,sd=None):  
    """  
    Returns the reduced chi-square error statistic for an arbitrary model,   
    chisq/nu, where nu is the number of degrees of freedom. If individual   
    standard deviations (array sd) are supplied, then the chi-square error   
    statistic is computed as the sum of squared errors divided by the standard   
    deviations. See http://en.wikipedia.org/wiki/Goodness_of_fit for reference.  
     
    ydata,ymod,sd assumed to be Numpy arrays. deg integer.  
   
    Usage:  
    >>> chisq=redchisqg(ydata,ymod,n,sd)  
    where  
    ydata : data  
    ymod : model evaluated at the same x points as ydata  
    n : number of free parameters in the model  
    sd : uncertainties in ydata  
   
    Rodrigo Nemmen  
    http://goo.gl/8S1Oo  
    """  
    # Chi-square statistic  
    if sd is None:  
        chisq=np.sum((ydata-ymod)**2)  
    else:  
        chisq=np.sum( ((ydata-ymod)/sd)**2 )  
             
    # Number of degrees of freedom assuming 2 free parameters  
    nu=ydata.size-1-deg  
    return chisq/nu    


# Read in HiRes data.
hr = spacepy.datamodel.readJSONheadedASCII(fDir + fName)
t = np.array([dateutil.parser.parse(i) for i in hr['Time']])

# Filter data by start/end times
eventInd = np.where((t >= startTime) & (t <= endTime))[0]
t = t[eventInd]
hr['Col_counts'] = hr['Col_counts'][eventInd, :]
hr['Col_counts_baseline'] = hr['Col_counts_baseline'][eventInd, :]

# Normalize time
tSec = np.array([(i - t[0]).total_seconds() for i in t])

# Set initial parameters
chArr = np.arange(4)
c = ['r', 'b', 'g', 'c']

# Prepare plot environment
c = np.array(['r', 'b', 'g', 'c'])
fig = plt.figure(figsize=(20, 5), dpi=80, facecolor = 'grey')
gs = gridspec.GridSpec(1,1)
dataPlt = fig.add_subplot(gs[0, 0])

for energy in chArr:
    if energy != 1:
        continue
    if sc_id == 3:
        if energy == 2:
            p0 = np.array([100, 0.5, 100E-3, 
            50, 1.3, 25E-3, 50, 1.6, 50E-3, 32, 2.25, 50E-3,
            20, 2.8, 25E-3, 23, 3.1, 25E-3])
#        elif energy == 3:
#            p0 = np.array([100, 0.5, 100E-3, 100, 1, 10E-3, 
#            50, 1.3, 25E-3, 50, 1.6, 50E-3, 32, 2.25, 50E-3,
#            20, 2.8, 25E-3])
        else:
            p0 = np.array([100, 0.5, 100E-3, 100, 1, 10E-3, 
            50, 1.3, 25E-3, 50, 1.6, 50E-3, 32, 2.25, 50E-3,
            20, 2.8, 25E-3, 23, 3.1, 25E-3])

    elif sc_id == 4:
        if energy == 3:
            p0 = np.array([100, 0.5, 100E-3, 4, 1.1, 200E-3,
            2, 1.6, 50E-3, 10, 2.1, 10E-3])
            p0 = np.array([100, 0.5, 100E-3, 4, 1.1, 200E-3,
            2, 1.6, 50E-3, 2, 2.2, 50E-3, 2, 2.5, 100E-3])
        else:
            p0 = np.array([100, 0.5, 100E-3, 50, .9, 50E-3,
            50, 1.25, 50E-3, 25, 1.6, 50E-3, 25, 2.2, 25E-3,
            25, 3, 50E-3])

    # Do the fitting and calculate the errors from the diagonal terms 
    # of the covariance matrix. 
    popt, pcov = scipy.optimize.curve_fit(nGaus, tSec, hr['Col_counts'][:, energy] - hr['Col_counts_baseline'][:, energy], p0 = p0, sigma =  1 + np.sqrt(hr['Col_counts'][:, energy] + hr['Col_counts_baseline'][:, energy]), method = 'lm', absolute_sigma = False)
    perr = np.sqrt(np.diag(pcov))
    
    # Calculate and print the reduced ChiSquared
    data = np.array(hr['Col_counts'][:, energy])
    model = np.array(nGaus(tSec, *popt))
    print('Reduced Chi Squared', redchisqg(data, model, deg = len(popt), sd = 1 + np.sqrt(hr['Col_counts'][:, energy])))
    #goodnessOfFit = scipy.stats.chisquare(data, f_exp = model)[0]/(len(data) - len(popt))
    #print(scipy.stats.chisquare(data, f_exp = model))
    #print(goodnessOfFit)
    
    # Plot the data
    dataPlt.errorbar(tSec, hr['Col_counts'][:, energy] - 
    hr['Col_counts_baseline'][:, energy], c = c[energy],
    label = 'Channel ' + str(energy) + ' counts',
    lw = 1, fmt = 'o', markersize = 3,
    yerr = np.sqrt(hr['Col_counts'][:, energy]))
    if smooth_plt:
        kernelWidth = 8
        kernel = np.ones(kernelWidth)/kernelWidth 
        dataPlt.plot(tSec, np.convolve(kernel, (hr['Col_counts'][:, energy] - 
        hr['Col_counts_baseline'][:, energy]), mode = 'same'), '--', 
        c = c[energy], lw = 3)
    # Plot the fit superposition of Gaussians        
    dataPlt.plot(tSec, nGaus(tSec, *popt),
    label = 'Channel ' + str(energy) + ' fit',
    lw = 2, ls = '-', c = c[energy])
    
    # Print the peak times
    for i in range(len(popt)//3):
        print(energy, 'Peak ', i, t[0] + timedelta(seconds = popt[3*i + 1]), perr[3*i + 1])
    
    # Calculate and print the average bounce period
    t0 = ufloat(popt[1], perr[1])
    tend = ufloat(popt[len(popt)- 5], perr[len(popt)- 5])
    if sc_id == 3:
        print('Bounce period', (tend - t0)/4)
    else:
        print('Bounce period', (tend - t0)/3)
    print(popt[1], popt[len(popt)- 5])
    print(t[0] + timedelta(seconds = popt[1]), t[0] + 
    timedelta(seconds = popt[len(popt)- 5]))

dataPlt.set_yscale('log')
dataPlt.set_ylim(bottom = 10**0)
plt.show()
