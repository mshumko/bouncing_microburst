import matplotlib.pyplot as plt
import spacepy.datamodel
import matplotlib.gridspec as gridspec
import dateutil.parser
import numpy as np
from datetime import datetime

fDir = '/home/mike/FIREBIRD/Datafiles/FU_3/hires/level2/'
fName = 'FU3_Hires_2015-02-02_L2.txt'
hr = spacepy.datamodel.readJSONheadedASCII(fDir + fName)
t = np.array([dateutil.parser.parse(i) for i in hr['Time']])
startTime = datetime(2015, 2, 2, 6, 12, 49, 0)
endTime = datetime(2015, 2, 2, 6, 13, 00, 0)
ind = np.where((t > startTime) 
    & (t < endTime))[0]

fig = plt.figure(figsize=(15, 10), dpi=80, facecolor = 'white')
gs = gridspec.GridSpec(2, 1)
surPlt = fig.add_subplot(gs[0,0])
colPlt = fig.add_subplot(gs[1, 0])

for i in range(6):
     surPlt.plot(t[ind], hr['Sur_counts'][ind, i])
     colPlt.plot(t[ind], hr['Col_counts'][ind, i])

surPlt.set(title = 'FU3 HiRes data from 2015-02-02', ylabel = 'Surface Counts', yscale = 'log')
colPlt.set(ylabel = 'Collimated Counts', xlabel = 'UTC', yscale = 'log')
plt.savefig('/home/mike/Dropbox/0_firebird_research/microburst_characterization/plots/bouncing_packet/FU3_compare_sur_col.png')
plt.show()

