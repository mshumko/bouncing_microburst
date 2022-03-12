# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 21:32:58 2017

@author: mike
"""
import sys
sys.path.insert(0, '/home/mike/FIREBIRD/irbem-code/python')
from IRBEM import IRBEM
import numpy as np

from lat_long_map import latLongScaleCalc
# Map scale sizes map.
latkm = 1/122.5 #degrees
longkm = lambda lat, alt: 1/latLongScaleCalc(lat, alt) # degrees

def lineOfSight(X1, X2):
    """
    Calculates the line of sight between two vectors in GEO coordinates. The
    function returns the line of sight distance in km.
    """
    Re = 6371 # km
    D = Re*np.sqrt((X1[0] - X2[0])**2 + (X1[1] - X2[1])**2 + (X1[2] - X2[2])**2)
    return D

if __name__ == '__main__':
    # Map the microburst decay event to the equatorial region using IRBEM.
    # alti, lati, East longi, time

    # North end of the event
    latS = 28.81; latErr = 0.8
    
    startLatX = {'x1':652, 'x2':63.35, 'x3':15.18, 'dateTime':'2015-02-02T06:12:53'}
    # South end of the event
    endLatX = {'x1':656, 'x2':63.12, 'x3':15.18, 'dateTime':'2015-02-02T06:12:59'}
    
    # Min west scale # 31.49 pm 3.80 km
    # Max west scale # 40.43 pm 8.82 km
    minLonS = 38.54; minLonErr = 8.82 
    maxLonS = 50.78; maxLonErr = 11.40
    
    startLongXsmall = {'x1':651, 'x2':63.22, 'x3':15.18 - minLonS*longkm(63.35, 650), 'dateTime':'2015-02-02T06:12:53'}
    startLongXlarge = {'x1':651, 'x2':63.22, 'x3':15.18 - maxLonS*longkm(63.35, 650), 'dateTime':'2015-02-02T06:12:53'}
    # futhest seen eastward. 
    endLongX = {'x1':651, 'x2':63.36, 'x3':15.29, 'dateTime':'2015-02-02T06:12:59'}
    
    # Run the model.
    model = IRBEM(options = [0,0,0,0,0])
    maginput = {'Kp':40.0}
    
    # LATITUDINAL SCALE SIZE AT MAG EQUATOR CALCULATION:
    startLatGEO = model.find_magequator(startLatX, maginput)
    endLatGEO = model.find_magequator(endLatX, maginput)
    magRadScale = float(lineOfSight(startLatGEO['XGEO'], endLatGEO['XGEO']))
    print('Radial scale size at the magnetic equator >=', \
    '%.2f' % magRadScale, '$\pm$ %.2f' % (latErr*magRadScale/latS), ' km')
    
    # LONGIUDINAL SCALE SIZE AT MAG EQUATOR CALCULATION:
    startLongGEOsmall = model.find_magequator(startLongXsmall, maginput)
    startLongGEOlarge = model.find_magequator(startLongXlarge, maginput)
    endLongGEO = model.find_magequator(endLongX, maginput)
    smallLongScale = float(lineOfSight(startLongGEOsmall['XGEO'], \
    endLongGEO['XGEO']))
    largeLongScale = float(lineOfSight(startLongGEOlarge['XGEO'],\
    endLongGEO['XGEO']))
    print('408 keV azimuthal scale size at the magnetic equator >=', \
    '%.2f' % smallLongScale, '$\pm$ %.2f' % (minLonErr*smallLongScale/minLonS)\
    , ' km')
    print('555 keV azimuthal scale size at the magnetic equator >=', \
    '%.2f' % largeLongScale, '$\pm$ %.2f' % (maxLonErr*largeLongScale/maxLonS)\
    , ' km')
