# Calculate the bounce period using IRBEM.
from datetime import datetime

import numpy as np
import IRBEM

from event1.load_hires import load_hires

def bounce_period(sc_id, time, energy=1000):
    hr = load_hires(sc_id, time)
    idt = np.where(hr['Time'] > time)[0][0]

    X = {'Time':hr['Time'][idt], 'x1':500, 'x2':hr['Lat'][idt], 'x3':hr['Lon'][idt]}
    maginput = {'Kp':hr['kp'][idt]}

    m = IRBEM.MagFields(kext='T89')
    print(m.make_lstar(X, maginput))
    return m.bounce_period(X, maginput, energy, alpha=64)

if __name__ == '__main__':
    tb = bounce_period(3, datetime(2016, 1, 21, 22, 47, 11))
    print(f'The bounce period is {tb} seconds.')