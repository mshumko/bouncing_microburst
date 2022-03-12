import IRBEM
import numpy as np
import matplotlib.pylab as plt
import spacepy.datamodel
import dateutil.parser
import copy
from mpl_toolkits.mplot3d import Axes3D

X = {}
X['x1'] = 651
X['x2'] = 63
X['x3'] = 15.9
X['dateTime'] = '2015-02-02T06:12:43'

X2 = copy.deepcopy(X)

# Get TS04 field inputs
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


def azimuthalFieldLineVisualization(ax, lat = 50, dLon = 30, pltDensity = 10):
    """
    This function draws the megnetic field lines defined by the lat argument,
    at different longitudes from 0 to 360, in dLon angle steps. pltDensity 
    defines at what inerval to plot the field lines, since it is very 
    computationaly expensive to plot. 
    """
    model2 = IRBEM.IRBEM(options = [0,0,0,0,0], kext = 'T89', verbose = False)
    startLon = 0
    endLon = 360
    
    N = (endLon - startLon)//dLon
    # We will have to append since we can't tell how big the output will be
    xGEO = np.array([])
    yGEO = np.array([])
    zGEO = np.array([])
    
    for i in range(N):#np.arange(startLon, endLon, dLon):
        #maginput = {'Kp':0.0}
        X['x3'] = i*dLon
        print(X)
        out = model2.trace_field_line(X, maginput)
        # pltDensity is to plot every pltDensity location of the field line,
        # to ease the graphical visualization. 
        if len(out['POSIT']) == 0:
            continue
        xGEO = np.append(xGEO, out['POSIT'][::pltDensity, 0])
        yGEO = np.append(yGEO, out['POSIT'][::pltDensity, 1])
        zGEO = np.append(zGEO, out['POSIT'][::pltDensity, 2])
    
    # Now plot the field line
    ax.plot(xGEO, yGEO, zGEO)

"""
Test function to plot a fieldline and a sphere.
"""
model = IRBEM.IRBEM(options = [0,0,0,0,0], kext = 'T89', verbose = False)
#maginput = {'Kp':40.0}
pltDensity = 10
out = model.trace_field_line(X, maginput)

# Now plot the field lines
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
xGEO = out['POSIT'][::pltDensity, 0] 
yGEO = out['POSIT'][::pltDensity, 1] 
zGEO = out['POSIT'][::pltDensity, 2] 

ax.plot(xGEO, yGEO, zGEO, linewidth = 5)

# Draw sphere    
u, v = np.mgrid[0:2*np.pi:40j, 0:np.pi:20j]
x=np.cos(u)*np.sin(v)
y=np.sin(u)*np.sin(v)
z=np.cos(v)
ax.plot_wireframe(x, y, z, color="k")
ax.set_ylim([-10, 10])
ax.set_xlim([-10, 10])
ax.set_zlim([-5, 5])
ax.set_xlabel('x GEO')
ax.set_ylabel('y GEO')
ax.set_zlabel('z GEO')

azimuthalFieldLineVisualization(ax)

plt.show()
