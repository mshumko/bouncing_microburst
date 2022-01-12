from datetime import datetime, timedelta

import numpy as np
import scipy.signal
import matplotlib.pyplot as plt
import matplotlib.ticker
import matplotlib.dates

from event1.load_hires import load_hires
from event1.bounce_period import bounce_period

time_range = [
    datetime(2015, 5, 25, 18, 18, 47),
    datetime(2015, 5, 25, 18, 18, 58)
]
sc_id = 4
tb = bounce_period(sc_id, time_range[0])

hr = load_hires(sc_id, time_range[0])
idt = np.where(
    (hr['Time'] >= time_range[0]) &
    (hr['Time'] <= time_range[1])
    )[0]

for key in hr:
    if len(hr[key].shape) == 1:
        hr[key] = hr[key][idt]
hr['Col_counts'] = hr['Col_counts'][idt, 0]
fig, ax = plt.subplots(2, sharex=True)

ax[0].plot(hr['Time'], hr['Col_counts'], 'k')
s = (
    f'L={round(hr["McIlwainL"][0], 1)}\n'
    f'MLT={round(hr["MLT"][0], 1)}\n'
    f'(lat,lon)=({round(hr["Lat"][0], 1)}, {round(hr["Lon"][0], 1)})'
    )
ax[0].text(0.7, 0.99, s, va='top', transform=ax[0].transAxes)
ax[0].set(
    title=f'FU{sc_id} | {time_range[0].date()} | bouncing packet',
    ylabel=f'CH0 [counts/{hr.attrs["CADENCE"]*1000} ms]'
)

freqs, times, spectrogram = scipy.signal.spectrogram(
    hr['Col_counts'], fs=1/hr.attrs['CADENCE'], 
    nperseg=50, noverlap=25)
spectrogram_times = [hr['Time'][0] + timedelta(seconds=dt) for dt in times]
p = ax[1].pcolormesh(
    spectrogram_times, freqs, spectrogram[:-1, :-1], 
    norm=matplotlib.colors.LogNorm(vmin=10**0, vmax=10**3)
    )
# plt.colorbar(p)
ax[1].text(0, 1, f'Modeled $f_b$={round(1/tb, 1)} Hz ($t_b$={round(tb, 3)} s)', 
    c='white', fontsize=13, va='top', transform=ax[1].transAxes) 
ax[1].text(spectrogram_times[0]+timedelta(seconds=0.5), 1/tb+0.5, f'$f_b$', 
    c='white', fontsize=13, ha='right')    
ax[1].axhline(1/tb, c='w', ls=':')
ax[1].text(spectrogram_times[0]+timedelta(seconds=0.5), 2/tb+0.5, f'2$f_b$', 
    c='white', fontsize=13, ha='right')
ax[1].axhline(2/tb, c='w', ls=':')
ax[1].text(spectrogram_times[0]+timedelta(seconds=0.5), 3/tb+0.5, f'3$f_b$', 
    c='white', fontsize=13, ha='right')
ax[1].axhline(3/tb, c='w', ls=':')

# locator=matplotlib.ticker.MaxNLocator(nbins=5)
# ax[-1].xaxis.set_major_locator(locator)
ax[-1].xaxis.set_major_locator(matplotlib.dates.SecondLocator(interval=2))
fmt = matplotlib.dates.DateFormatter('%H:%M:%S')
ax[-1].xaxis.set_major_formatter(fmt)
ax[-1].set(ylabel=f'Frequency [Hz]', xlabel='Time')
ax[-1].xaxis.set_minor_locator(matplotlib.dates.SecondLocator())

plt.tight_layout()
plt.show()