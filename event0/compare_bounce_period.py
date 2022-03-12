# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 10:13:20 2017

@author: mike
"""

import sys
#sys.path.insert(0, '/home/mike/Dropbox/0_firebird_research/IRBEM/scripts')
#sys.path.insert(0, '/home/mike/FIREBIRD/irbem-code/python/')
import numpy as np
import spacepy.datamodel
import datetime
import dateutil.parser
from IRBEM import MagFields
import matplotlib.pylab as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
plt.rcParams.update({'font.size': 18})

#import matplotlib
#font = {'family' : 'normal',
#        'size'   : 18}
#matplotlib.rc('font', **font)

PLOT_MODELS = True

plt.rc('font', family='serif', size=20)

import copy

# Schulz and Lanzerotti Bounce period equation.
Tsl = lambda L, alpha0, v: 4*6.371E6*np.divide(L, v) * \
       (1.3802 - 0.3198*(np.sin(np.deg2rad(alpha0)) + \
       np.sqrt(np.sin(np.deg2rad(alpha0)))))
beta = lambda Ek: np.sqrt(1-((Ek/511)+1)**(-2))
c = 3.0E8 # m/s
EkArr = np.arange(200, 800)

##### OBSERVED BOUNCE PERIOD DATA ########
Tobs = {}
Tobs['FU3'] = np.array([0.577, 0.576, 0.571, 0.543])
Tobs['FU3err'] = np.array([0.005, 0.010, 0.005, 0.017])

Tobs['FU4'] = np.array([0.576, 0.580, 0.554, 0.547])
Tobs['FU4err'] = np.array([0.006, 0.011, 0.009, 0.003])

Tobs['E_FU3'] = [np.arange(231, 300), np.arange(300, 408), np.arange(408, 555), 
            np.arange(555, 771)]    
Tobs['E_FU4'] = [np.arange(220, 283), np.arange(283, 383), np.arange(383, 520), 
            np.arange(520, 720)]
        
           
Tobs['Oleksiy_FU3_max'] = np.array([0.66,0.67])
Tobs['Oleksiy_FU3_min'] = np.array([0.64,0.62])            
#Tobs['Oleksiy_FU3_max'] = np.array([0.66,0.67, 0.66])
#Tobs['Oleksiy_FU3_min'] = np.array([0.64,0.62, 0.5])

if PLOT_MODELS:
    # Run the IRBEM bounce time calculator for T89
    X = {}
    X['x1'] = 651
    X['x2'] = 63
    X['x3'] = 15.9
    X['dateTime'] = '2015-02-02T06:12:43'
    maginput = {'Kp':40.0}
    model = MagFields(options = [0,0,0,0,0])
    TbT89 = model.bounce_period(X, maginput, EkArr)

    # Run the bounce time calculator for T05.
    model = MagFields(options = [0,0,0,0,0], kext = 11)
    omniLoc = '/home/mike/.spacepy/data/omnidata.h5'
    omniData = spacepy.datamodel.fromHDF5(omniLoc)
    omniT = np.array([dateutil.parser.parse(i.decode()) for i in omniData['UTC']])
    t = dateutil.parser.parse(X['dateTime'])
    idx = np.where(t >= omniT)[0][-1]

    # Prepare the magnetic field inputs
    T05Keys = ['Dst', 'Pdyn', 'ByIMF', 'BzIMF', 'W1', 'W2', 'W3', 'W4', 'W5', 'W6']
    maginput = {}
    for i in T05Keys:
        maginput[i] = float(omniData[i][idx])
    TbT05 = np.array(model.bounce_period(X, maginput, EkArr)) 

    # Run bounce period calculator on Olson & Pfitzer quiet.
    model = MagFields(options = [0,0,0,0,0], kext = 5)
    maginput = None
    TbOP = np.array(model.bounce_period(X, maginput, EkArr))

    # Calculate L shell for the Shulz and Lanzerotti model.
    model.make_lstar(X, maginput)
    L = np.abs(model.make_lstar_output['Lm'][0])
    TbSL = np.array(Tsl(L, 3.7, c*beta(EkArr)))

# Plot theoretical bounce times
fig = plt.figure(figsize=(13, 13), dpi=80, facecolor = 'white')
gs = gridspec.GridSpec(1,1)
tbPlt = fig.add_subplot(gs[0, 0])
if PLOT_MODELS:
    T89Tb = tbPlt.plot(EkArr, TbT89, 'r:', lw = 3, label = r'T89')
    T05Tb = tbPlt.plot(EkArr, TbT05, 'c--', lw = 3, label = r'T04')
    OPTb = tbPlt.plot(EkArr, TbOP, 'b-.', lw = 3, label = r'Olson & Pfitzer Quiet')
    SLPlt = plt.plot(EkArr, TbSL, 'k-', lw = 2, label = 'Dipole')
plt.xlabel('Energy (keV)')
plt.ylabel(r'$t_b$ (s)')

for i in range(len(Tobs['Oleksiy_FU3_max'])):
    FU3_oleksiy_patch = tbPlt.add_patch(patches.Rectangle(
    (Tobs['E_FU3'][i][0], Tobs['Oleksiy_FU3_min'][i]), len(Tobs['E_FU3'][i]), 
        Tobs['Oleksiy_FU3_max'][i] - Tobs['Oleksiy_FU3_min'][i], label = r'FU3 minima', alpha = 0.5, hatch = "o"))  

# Plot the observed bounce times as boxes.
#for i in range(len(Tobs['FU4err'])):
#    FU4_patch = tbPlt.add_patch(patches.Rectangle((Tobs['E'][i][0], 
#       Tobs['FU4'][i] - Tobs['FU4err'][i]), len(Tobs['E'][i]),
#         2*Tobs['FU4err'][i], color = 'b', label = r'FU4', alpha = 0.5, hatch = ' | '))
#    FU3_patch = tbPlt.add_patch(patches.Rectangle((Tobs['E'][i][0], 
#       Tobs['FU3'][i] - Tobs['FU3err'][i]), len(Tobs['E'][i]), 
#        2*Tobs['FU3err'][i], color = 'g', label = r'FU3', alpha = 0.5, hatch = '-'))
for i in range(len(Tobs['FU4err'])):
    FU4_patch = tbPlt.add_patch(patches.Rectangle((Tobs['E_FU4'][i][0], 
       Tobs['FU4'][i] - Tobs['FU4err'][i]), len(Tobs['E_FU4'][i]),
         2*Tobs['FU4err'][i], color = 'b', label = r'FU4 fit', alpha = 0.5, hatch = ' / '))
    FU3_patch = tbPlt.add_patch(patches.Rectangle((Tobs['E_FU3'][i][0], 
       Tobs['FU3'][i] - Tobs['FU3err'][i]), len(Tobs['E_FU3'][i]), 
        2*Tobs['FU3err'][i], color = 'g', label = r'FU3 fit', alpha = 0.5, hatch = "\ "))      

if PLOT_MODELS:
    tbPlt.legend(handles = [FU3_patch, FU4_patch, FU3_oleksiy_patch, T89Tb[0], T05Tb[0], OPTb[0], SLPlt[0]], fontsize = 20)
else:
    tbPlt.legend(handles = [FU3_patch, FU3_oleksiy_patch, FU4_patch], fontsize = 20)

tbPlt.set_xlim((Tobs['E_FU3'][0][0] - 50, Tobs['E_FU3'][-1][-1] + 50))
plt.title(r'Theoretical and observed $t_b$ during bouncing packet microburst')
tbPlt.set_ylim((0.51, 0.7))
gs.tight_layout(fig)
if PLOT_MODELS:
    plt.savefig('/home/mike/Desktop/obs_model_bounce_period.pdf')
else:
    plt.savefig('/home/mike/Desktop/obs_bounce_period.pdf')
#plt.show()
