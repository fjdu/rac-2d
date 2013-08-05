import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.patches import Polygon
from matplotlib.mlab import griddata
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.collections import PolyCollection
import os
import os.path
import my_script as my

r_pos = 10.0

data_dir = './results/results_20130618_nirgal_c/'
filename_save_results =  os.path.join(data_dir, 'iter_0008.dat')

data = np.loadtxt(filename_save_results, comments='!')

filelen = len(data[:,0])

f = open(filename_save_results, 'r')
str_comment = f.readline()[1:].split()
f.close()

dic = {}
for i in xrange(len(str_comment)):
  dic.update({str_comment[i]: data[:, i]})

#dic.update({'Pressure': dic['Tgas']*dic['n_gas']})

pp = PdfPages(os.path.join(data_dir, 'plot_column.pdf'))

name_list = \
  [ \
    {'name': 'Tgas'    , 'scale': 'log', 'cmap': cm.rainbow,  'vr': (50, 1e4)},
    {'name': 'Tgrain'  , 'scale': 'linear', 'cmap': cm.rainbow},#  'vr': (50, 5e2)},
    {'name': 'n_gas'   , 'scale': 'log', 'cmap': cm.rainbow},#  'vr': (40, 2.5e10)},
    {'name': 'H2O'     , 'scale': 'log', 'cmap': cm.rainbow,  'vr': (1e-10, 1.7e-4)},
    {'name': 'OH'      , 'scale': 'log', 'cmap': cm.rainbow,  'vr': (1e-10, 3e-5)},
    {'name': 'H2'      , 'scale': 'log', 'cmap': cm.rainbow,  'vr': (1e-6, 0.5)},
    {'name': 'H'       , 'scale': 'log', 'cmap': cm.rainbow,  'vr': (1e-6, 1)},
    {'name': 'E-'      , 'scale': 'log', 'cmap': cm.rainbow,  'vr': (1e-9, 7e-4)},
    {'name': 'C'       , 'scale': 'log', 'cmap': cm.rainbow,  'vr': (1e-10, 7e-5)},
    {'name': 'C+'      , 'scale': 'log', 'cmap': cm.rainbow,  'vr': (1e-10, 7e-5)},
    {'name': 'CO'      , 'scale': 'log', 'cmap': cm.rainbow,  'vr': (1e-10, 7.3e-5)},
    {'name': 'f_H2'    , 'scale': 'log', 'cmap': cm.rainbow,  'vr': (1e-6, 1e0)},
    {'name': 'f_H2O'   , 'scale': 'log', 'cmap': cm.rainbow},#  'vr': (1e-3, 1e0)},
  ]

idx = []
for i in xrange(filelen):
  rmin = dic['rmin'][i]
  rmax = dic['rmax'][i]
  if (rmin <= r_pos) and (rmax > r_pos):
    idx.append(i)

nidx = len(idx)
z = 0.5 * (dic['zmin'][idx] + dic['zmax'][idx])
minz = np.min(z)
maxz = np.max(z)

for item in name_list:
  name = item['name']
  v = dic[name][idx]

  maxval = np.nanmax(v)
  minval = np.nanmin(v)
  minval_nonzero = np.nanmin(v[np.where(v > 0.0)])
  print name, ' Max, min ', maxval, minval
  if 'vr' in item:
    if len(item['vr']) >= 2:
      minval = item['vr'][0]
      minval_nonzero = item['vr'][0]
      maxval = item['vr'][1]

  #xRange = (minz/1.5, maxz*1.1)
  xRange = (0, 6)
  yRange = (minval/1.5, maxval*1.5)
  
  xlen = xRange[1] - xRange[0]
  ylen = xlen #yRange[1] - yRange[0]
  figsize_x = 10
  x_start = 0.13
  y_start = 0.15
  x_end = 0.98
  y_end = 0.90
  x_sep = 0.01
  y_sep = 0.08
  del_x = x_end - x_start - x_sep
  del_x_1 = 0.85 * del_x
  del_x_2 = 0.06 * del_x
  del_y = y_end - y_start

  ntick_cbar = 5
  
  figsize_y = figsize_x * (ylen / xlen) * del_x_1/del_y
  figsize = (figsize_x, figsize_y)
  
  fig = plt.figure(figsize=figsize)
  ax = fig.add_axes([x_start, y_start, del_x_1, del_y],
    xlabel='z (AU)', ylabel=name,
    autoscalex_on=False, autoscaley_on=False,
    xscale='linear', yscale=item['scale'],
    xlim=xRange, ylim=yRange)
  ax.xaxis.label.set_fontsize(25)
  ax.yaxis.label.set_fontsize(25)
  for label in ax.get_xticklabels():
    label.set_fontsize(20)
  for label in ax.get_yticklabels():
    label.set_fontsize(20)

  ax.plot(z, v,
      linestyle='-',
      color='blue', linewidth=5)
  
  plt.savefig(pp, format='pdf')

pp.close()
