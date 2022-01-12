# Calculate the footprint in the Sothern hemisphere IRBEM.
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
n_footprint = m.find_foot_point(X, maginput, 100, 1)
s_footprint = m.find_foot_point(X, maginput, 100, -1)
print(f"FIREBIRD location (alt, lat, lon) = ({hr['Alt'][idt]}, {hr['Lat'][idt]}, {hr['Lon'][idt]})")
print(f'northern footprint (alt, lat, lon) =  {n_footprint["XFOOT"]}')
print(f'southern footprint (alt, lat, lon) = {s_footprint["XFOOT"]}')