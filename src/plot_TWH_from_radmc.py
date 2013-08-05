import numpy as np
from matplotlib.mlab import griddata
import matplotlib.pyplot as plt

data = np.loadtxt('TWH_from_radmc.dat', comments='!')

xs0 = data[:,0]/1.5e13
ys0 = data[:,1]/1.5e13
zs0 = data[:,3]

N = 70j
extent = (0.0, 8.0, 0.0, 9.0)

xs,ys = np.mgrid[extent[0]:extent[1]:N, extent[2]:extent[3]:N]

resampled = griddata(xs0, ys0, zs0, xs, ys)

plt.imshow(np.log(resampled.T), extent=extent, origin='lower')
#plt.plot(xs0, ys0, "r.")
#plt.plot(xs, ys, "b.")
#plt.title("imshow for irregularly spaced data using griddata")
plt.show()
