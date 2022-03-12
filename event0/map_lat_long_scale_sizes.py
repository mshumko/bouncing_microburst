# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 11:46:28 2017

@author: mike
"""

from map_to_magnetic_equator import lineOfSight
from lat_long_map import latLongScaleCalc
import copy

import sys
sys.path.insert(0, '/home/mike/FIREBIRD/irbem-code/python')
from IRBEM import IRBEM
import numpy as np


# Map scale sizes map.
latkm = 1/122.5 #degrees
longkm = lambda lat, alt: 1/latLongScaleCalc(lat, alt) # degrees
    
def LEOCoords(X, kmscale = 1):
    """
    Given initial lat, long, alt coordinates in X in IRBEM format, will output
    the same format of X, but the coordinates will be adjusted by half of the
    km scale. latTop = X + lat(kmscale)/2, latBottom = X - lat(kmscale)/2, 
    longLeft = X - long(kmscale)/2, longRight = X + long(kmscale)/2
    """
    latTop, latBottom, longLeft, longRight = \
    [copy.deepcopy(X) for i in range(4)]
    latTop['x2'] += kmscale*latkm/2
    latBottom['x2'] -= kmscale*latkm/2
    longLeft['x3'] -= kmscale*longkm(X['x2'], X['x1'])/2
    longRight['x3'] += kmscale*longkm(X['x2'], X['x1'])/2
    return latTop, latBottom, longLeft, longRight
    
def mapScaleSizes(latArr, long, alt = 650, kp = 40):
    """
    This function calculates the line of sight latitudinal and longitudinal
    scale sizes for a 1 km^2 patch at LEO mapped to the magnetic equator. 
    
    INPUT: latArr, a single value or an array of latitudes, which will
    be transformed into L shell values. long is the longitude value. Alt and kp
    are optional kwargs, and are set to 650 km and 40 by default.
    """
    if (type(latArr) is not list) and (type(latArr) is not np.ndarray):
        latArr = np.array([latArr])
        
    mapLatScales, mapLongScales, L, MLT= \
                [-9999*np.ones(len(latArr)) for i in range(4)]
        
    model = IRBEM()
    maginput = {'Kp':kp}
    for i in range(len(latArr)):
        X = {'x1':alt, 'x2':latArr[i], 'x3':long, \
        'dateTime':'2015-02-02T06:12:53'}
        latTop, latBottom, longLeft, longRight = LEOCoords(X)
        magGeoCoords = model.make_lstar(X, maginput)
        L[i] = np.abs(magGeoCoords['Lm'][0])
        MLT[i] = magGeoCoords['MLT'][0]
        mappedLatTop = model.find_magequator(latTop, maginput)
        mappedLatBottom = model.find_magequator(latBottom, maginput)
        mappedLongLeft = model.find_magequator(longLeft, maginput)
        mappedLongRight = model.find_magequator(longRight, maginput)  
        
        mapLatScales[i] = lineOfSight(mappedLatTop['XGEO'], mappedLatBottom['XGEO'])
        mapLongScales[i] = lineOfSight(mappedLongLeft['XGEO'], mappedLongRight['XGEO'])
    
    return mapLatScales, mapLongScales, L, MLT
    
if __name__ == '__main__':

    fig = plt.figure(figsize=(8, 8), dpi=80, facecolor = 'grey')
    gs = gridspec.GridSpec(1,1)
    absScales = fig.add_subplot(gs[0, 0])
    ratioScale = absScales.twinx()
    
    kpArr = [0, 40]
    longArr = [-10, -95, -170, -270]    
    
    for kp in kpArr:
        for long in longArr:
            mapLatScales, mapLongScales, L, MLT = \
            mapScaleSizes(np.arange(1, 70).astype(float), long, kp = kp)
            validInd = np.where((L != 1.00000000e+31))[0]
    
            # Plot the scale sizes
            latPlt = absScales.plot(L[validInd], mapLatScales[validInd], \
            c = 'r',  label = 'Raidal scale')
            longPlt = absScales.plot(L[validInd], mapLongScales[validInd], \
            c = 'b',  label = 'Azimuthal scale')
            
            # Plot the ratios
            rPlt = ratioScale.plot(L[validInd], \
            mapLatScales[validInd]/mapLongScales[validInd], ls = '--', c = 'k', \
            label = 'lat/long')
            ratioScale.set_ylabel('lat/long')
            
            # Beautify the legend and set axis limits.
            plots = latPlt + longPlt + rPlt
            lns = [l.get_label() for l in plots]
            plt.legend(plots, lns, loc=2)
            absScales.set_xlim([1, 8])
            absScales.set_ylim([0, 50])
            ratioScale.set_ylim([0, 4])
            absScales.set_xlabel('Lm')
            absScales.set_ylabel('Size (km)')
            absScales.set_title('Size of a ' + r'$ 1 \ km^2$' + \
            ' patch at LEO, mapped to mag equator. MLT = ' + \
            str('%.1f' % np.mean(MLT)) + ' kp = ' + str(kp//10))
            gs.tight_layout(fig)
            
            saveName = 'one_sq_km_LEO_mapped_to_mag_eq_MLT_' + \
            str('%.0f' % np.mean(MLT)) + '_kp_' + str(kp//10)
            saveDir = '/home/mike/Dropbox/0_firebird_research/'+ \
            'microburst_characterization/plots/IRBEM_mapping/'
            fig.savefig(saveDir + saveName, dpi=fig.dpi)
            
            plt.sca(absScales)
            plt.cla()
            plt.sca(ratioScale)
            plt.cla()
    plt.close('all')