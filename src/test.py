import numpy as np
import matplotlib.pyplot as pl
nlevel = 50
xRange = (0, 10)
yRange = (0, 10)
vRange = (-5, 5)
levels = np.linspace(vRange[0], vRange[1], nlevel)

d = np.loadtxt('tmp1')
nx = 11
ny = 11
xx = np.reshape(d[:,0], (nx, ny))
yy = np.reshape(d[:,1], (nx, ny))
zz = np.reshape(d[:,2], (nx, ny))

pl.figure(1)
pl.contourf(xx, yy, zz, levels=levels, vmin=vRange[0], vmax=vRange[1])
pl.xlim(xRange)
pl.ylim(yRange)

d = np.loadtxt('tmp2')
nx = 201
ny = 201
xx = np.reshape(d[:,0], (nx, ny))
yy = np.reshape(d[:,1], (nx, ny))
zz = np.reshape(d[:,2], (nx, ny))

pl.figure(2)
pl.contourf(xx, yy, zz, levels=levels, vmin=vRange[0], vmax=vRange[1])
pl.xlim(xRange)
pl.ylim(yRange)
pl.show()
