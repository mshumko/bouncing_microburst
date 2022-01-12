# Calculate the bounce period using IRBEM.
from datetime import datetime

import numpy as np
import IRBEM

from event1.load_hires import load_hires

time = datetime(2015, 5, 25, 18, 18, 52)
hr = load_hires(4, time)
idt = np.where(hr['Time'] > time)[0][0]

X = {'Time':hr['Time'][idt], 'x1':hr['Alt'][idt], 'x2':hr['Lat'][idt], 'x3':hr['Lon'][idt]}
maginput = {'Kp':hr['kp'][idt]}

m = IRBEM.MagFields(kext='T89')
tb = m.bounce_period(X, maginput, 200)
print(f'The bounce period is {tb} seconds.')