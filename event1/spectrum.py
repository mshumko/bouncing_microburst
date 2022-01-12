from datetime import datetime, timedelta

# import numpy as np
import scipy.signal
import matplotlib.pyplot as plt
import matplotlib.ticker
import matplotlib.dates

from event1.load_hires import load_hires

time_range = [
    datetime(2015, 5, 25, 18, 18, 48),
    datetime(2015, 5, 25, 18, 18, 58)
]
hr = load_hires(4, time_range[0])

fig, ax = plt.subplots(2, sharex=True)

ax[0].plot(hr['Time'], hr['Col_counts'][:, 0], 'k')
# ax[0].set_xlim(time_range)
locator=matplotlib.ticker.MaxNLocator(nbins=5)
ax[-1].xaxis.set_major_locator(locator)
fmt = matplotlib.dates.DateFormatter('%H:%M:%S')
ax[-1].xaxis.set_major_formatter(fmt)

freqs, times, spectrogram = scipy.signal.spectrogram(
    hr['Col_counts'][:, 0], fs=1/hr.attrs['CADENCE'])
spectrogram_times = [hr['Time'][0] + timedelta(seconds=dt) for dt in times]
p = ax[1].pcolormesh(
    spectrogram_times, freqs, spectrogram[:-1, :-1], 
    norm=matplotlib.colors.LogNorm(vmin=None, vmax=None)
    )
plt.tight_layout()
plt.show()