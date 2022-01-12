from datetime import datetime

import matplotlib.pyplot as plt
import matplotlib.ticker
import matplotlib.dates

from event1.load_hires import load_hires

time_range = [
    datetime(2015, 5, 25, 18, 18, 48),
    datetime(2015, 5, 25, 18, 18, 58)
]
hr = load_hires(4, time_range[0])

fig, ax = plt.subplots()

ax.plot(hr['Time'], hr['Col_counts'][:, 0], 'k')
ax.set_xlim(time_range)
locator=matplotlib.ticker.MaxNLocator(nbins=5)
ax.xaxis.set_major_locator(locator)
fmt = matplotlib.dates.DateFormatter('%H:%M:%S')
ax.xaxis.set_major_formatter(fmt)
plt.tight_layout()
plt.show()