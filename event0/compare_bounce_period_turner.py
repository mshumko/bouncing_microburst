# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 10:13:20 2017

@author: mike
"""

import sys
#sys.path.insert(0, '/home/mike/Dropbox/0_firebird_research/IRBEM/scripts')
sys.path.insert(0, '/home/mike/FIREBIRD/irbem-code/python/')
import numpy as np
import spacepy.datamodel
import datetime
import dateutil.parser
#from IRBEM_bounce_period import find_bounceperiod
from IRBEM import IRBEM
import matplotlib.pylab as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
plt.rcParams.update({'font.size': 15})
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
#Tobs['FU3'] = np.array([0.576, 0.579, 0.573])
#Tobs['FU4'] = np.array([0.60, 0.57, 0.531])
## Include error bars
#Tobs['FU3err'] = np.array([0.005, 0.007, 0.008])
#Tobs['FU4err'] = np.array([0.010, 0.010, 0.007])

# Fit data
Tobs['FU3'] = np.array([0.577, 0.577, 0.570, 0.546])
Tobs['FU3err'] = np.array([0.004, 0.009, 0.005, 0.010])

Tobs['FU4'] = np.array([0.576, 0.581, 0.551, 0.547])
Tobs['FU4err'] = np.array([0.005, 0.010, 0.008, 0.003])

Tobs['E'] = [np.arange(220, 300), np.arange(300, 408), np.arange(408, 555), 
            np.arange(555, 771)]

# Run the IRBEM bounce time calculator for T89
X = {}
X['x1'] = 651
X['x2'] = 62
X['x3'] = 15.9
X['dateTime'] = '2015-02-02T06:12:43'
#maginput = {'Kp':40.0}
#model = IRBEM(options = [0,0,0,0,0])
#TbT89 = model.bounce_period(X, maginput, EkArr)

# Run the bounce time calculator for T05.
model = IRBEM(options = [0,0,0,0,0], kext = 11)
omniLoc = '/home/mike/.spacepy/data/omnidata.h5'
omniData = spacepy.datamodel.fromHDF5(omniLoc)
omniT = np.array([dateutil.parser.parse(i.decode()) for i in omniData['UTC']])
t = dateutil.parser.parse(X['dateTime'])
idx = np.where(t >= omniT)[0][-1]
print(omniT[idx])

# Prepare the magnetic field inputs
T05Keys = ['Dst', 'Pdyn', 'ByIMF', 'BzIMF', 'W1', 'W2', 'W3', 'W4', 'W5', 'W6']
maginput = {}
for i in T05Keys:
    maginput[i] = float(omniData[i][idx])
    
model.make_lstar(X, maginput, STATUS_FLAG = False)
print(np.abs(model.lstar1_output['Lm'][0]))
L = np.abs(model.lstar1_output['Lm'][0])
TbT05 = model.bounce_period(X, maginput, EkArr)

# Run bounce period calculator on Olson & Pfitzer quiet.
#model = IRBEM(options = [0,0,0,0,0], kext = 5)
#maginput = None
#TbOP = model.bounce_period(X, maginput, EkArr)

# Calculate L shell for the Shulz and Lanzerotti model.
#model.make_lstar(X, maginput, STATUS_FLAG = False)
#L = np.abs(model.lstar1_output['Lm'][0])
#TbSL = Tsl(L, 3.7, c*beta(EkArr))

# Plot theoretical bounce times
fig = plt.figure(figsize=(10, 10), dpi=80, facecolor = 'white')
gs = gridspec.GridSpec(1,1)
tbPlt = fig.add_subplot(gs[0, 0])
#T89Tb = tbPlt.plot(EkArr, TbT89, 'r:', label = r'$t_b$ from T89 (Kp = 4)')
T05Tb = tbPlt.plot(EkArr, TbT05, 'c--', label = (r'Fudged $t_b$ from T05, L = %.2f' % L))
#OPTb = tbPlt.plot(EkArr, TbOP, 'b-.', label = r'$t_b$ from Olson & Pfitzer quiet')
#SLPlt = plt.plot(EkArr, TbSL, 'k-', label = 'Schulz and Lanzerotti 1974')
plt.xlabel('Energy (keV)')
plt.ylabel(r'$t_b$ (s)')


# Plot the observed bounce times as boxes.
for i in range(len(Tobs['FU4err'])):
    FU4_patch = tbPlt.add_patch(patches.Rectangle((Tobs['E'][i][0], 
       Tobs['FU4'][i] - Tobs['FU4err'][i]), len(Tobs['E'][i]),
         2*Tobs['FU4err'][i], color = 'b', label = r'Observed FU4 $t_b$', alpha = 0.5, hatch = ' | '))
    FU3_patch = tbPlt.add_patch(patches.Rectangle((Tobs['E'][i][0], 
       Tobs['FU3'][i] - Tobs['FU3err'][i]), len(Tobs['E'][i]), 
        2*Tobs['FU3err'][i], color = 'g', label = r'Observed FU3 $t_b$', alpha = 0.5, hatch = '-'))
    
tbPlt.legend(handles = [T05Tb[0], FU3_patch, FU4_patch], fontsize = 20)

plt.title(r'Theoretical and observed $t_b$ during bouncing packet microburst')
tbPlt.set_ylim((0.5, 0.7))
gs.tight_layout(fig)
plt.show()
