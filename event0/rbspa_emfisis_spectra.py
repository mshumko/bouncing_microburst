import numpy as np
import spacepy.datamodel
import dateutil.parser
from datetime import datetime, timedelta
from matplotlib import gridspec, pylab as plt

import sys
sys.path.insert(0, '/home/mike/Dropbox/0_aerospace_Internship')
from plot_emfisis_spectrogram import plotWFRSpectra

# Load in the magEpehem data
magEphem = spacepy.datamodel.readJSONheadedASCII('/home/mike/Dropbox/0_firebird_research/microburst_characterization/scipts/bounce_event/rbspa_def_MagEphem_OP77Q_20150202_v1.0.0.txt')
#print('EDMAG_MLT' in magEphem.keys())

magT = np.array([dateutil.parser.parse(i) for i in magEphem['DateTime']])

date = datetime(2015, 2, 2)
tBounds = [datetime(2015,2,2,2), datetime(2015,2,2,10)]
fig = plt.figure(figsize=(10, 8), dpi=80, facecolor = 'white')
gs = gridspec.GridSpec(2,1)
ax = fig.add_subplot(gs[0, 0], facecolor='k')
bx = fig.add_subplot(gs[1, 0], sharex = ax)
validInd = np.where((np.array(magEphem['L']) > 0) & (np.array(magEphem['L']) < 10))[0]
bx.plot(magT[validInd], magEphem['L'][validInd], 'k', label = 'L')
bx.plot(magT[validInd], magEphem['EDMAG_MLT'][validInd], '--k', label = 'MLT')
bx.set_xlim((tBounds[0], tBounds[1]))
bx.set_ylim((0, 25))
plotWFRSpectra('A', date, tBounds = tBounds, cLevels = .1, ax = ax)
gs.tight_layout(fig)
plt.show()
