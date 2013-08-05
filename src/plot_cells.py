import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from matplotlib.patches import Polygon
from matplotlib.mlab import griddata
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.collections import PolyCollection
matplotlib.rcParams['axes.linewidth'] = 0.001

pp = PdfPages('cells.pdf')

small_number = 1E-6

xRange = (0E0, 57.0)
yRange = (0E0, 57.0)

data_file = 'cells.dat'
dtmp = np.loadtxt(data_file, comments='!')

#minval = 1.0E-2
maxval = np.max(dtmp[:,1])
minval = np.min(dtmp[:,1])
log_max = np.log(maxval)
log_min = np.log(minval)

print maxval, minval

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111,
  aspect='equal',
  xlabel='r', ylabel='z',
  autoscalex_on=False, autoscaley_on=False,
  #xscale='log', yscale='log',
  xlim=xRange, ylim=yRange)
#plt.axis('off')
#ax.axhline(linewidth=0.001)

colormap = cm.rainbow

poly = np.empty((dtmp.shape[0], 4, 2))
poly[:,0,:] = dtmp[:,[2,3]]
poly[:,1,:] = dtmp[:,[4,3]]
poly[:,2,:] = dtmp[:,[4,5]]
poly[:,3,:] = dtmp[:,[2,5]]
sca_col = (np.log(dtmp[:,1]+small_number)-log_min) / (log_max-log_min)
poly_coll = PolyCollection(poly, array=sca_col, cmap=colormap, edgecolors='black', closed=True, linewidth=0.002)
ax.add_collection(poly_coll)

ax.add_artist(plt.Circle((0,0), radius=0.005, linewidth=0.002, color='red', fill=True, zorder=10, clip_on=False))

plt.savefig(pp, format='pdf')

pp.close()

#plt.show()
