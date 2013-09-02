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
from scipy import interpolate
from matplotlib.ticker import AutoMinorLocator

r_pos = 0.55

data_dir = '/Users/fdu/work/protoplanetary_disk/res/results_20130830_gondolin_a_4/'
filename_save_results =  os.path.join(data_dir, 'iter_0004.dat')

data = np.loadtxt(filename_save_results, comments='!')

filelen = len(data[:,0])

f = open(filename_save_results, 'r')
str_comment = f.readline()[1:].split()
f.close()

dic = {}
for i in xrange(len(str_comment)):
  dic.update({str_comment[i]: data[:, i]})

#dic.update({'Pressure': dic['Tgas']*dic['n_gas']})

color_list = ['blue', 'red', 'green', 'magenta',
    (0.7,0.5,0.3), (0.8,0.85,0.3), (0.2,0.2,0.2), (0.5,0.5,0.5)]

pp = PdfPages(os.path.join(data_dir, 'plot_column.pdf'))

name_list = \
  [ \
    {'name': 'Tgas'    , 'scale': 'log', 'cmap': cm.rainbow,
        'vr': (50, 1e4), 'hold_on':True},
    {'name': 'Tdust'   , 'scale': 'linear', 'cmap': cm.rainbow,
        'hold_on':False},
    {'name': 'O+'       , 'scale': 'log', 'cmap': cm.rainbow,
        'vr': (1e-10, 2e-4), 'hold_on': True},
    {'name': 'O'       , 'scale': 'log', 'cmap': cm.rainbow,
        'vr': (1e-10, 2e-4), 'hold_on': True},
    {'name': 'O2'      , 'scale': 'log', 'cmap': cm.rainbow,
        'vr': (1e-10, 2e-4), 'hold_on': True},
    {'name': 'H2O'     , 'scale': 'log', 'cmap': cm.rainbow,
        'vr': (1e-10, 1.7e-4), 'hold_on': True},
    {'name': 'CO'      , 'scale': 'log', 'cmap': cm.rainbow,
        'vr': (1e-10, 7.3e-5), 'hold_on':False},
    {'name': 'C'       , 'scale': 'log', 'cmap': cm.rainbow,
        'vr': (1e-10, 7e-5), 'hold_on':True},
    {'name': 'C+'      , 'scale': 'log', 'cmap': cm.rainbow,
        'vr': (1e-10, 7e-5), 'hold_on':True},
    {'name': 'CO'      , 'scale': 'log', 'cmap': cm.rainbow,
        'vr': (1e-10, 7.3e-5), 'hold_on':True},
    {'name': 'CO2'     , 'scale': 'log', 'cmap': cm.rainbow,
        'vr': (1e-10, 7.3e-5), 'hold_on':True},
    {'name': 'HC3N'    , 'scale': 'log', 'cmap': cm.rainbow,
        'vr': (1e-10, 7.3e-5), 'hold_on':True},
    {'name': 'C9N'     , 'scale': 'log', 'cmap': cm.rainbow,
        'vr': (1e-10, 7.3e-5), 'hold_on':True},
    {'name': 'gCO2'    , 'scale': 'log', 'cmap': cm.rainbow,
        'vr': (1e-10, 7.3e-5), 'hold_on':True},
    {'name': 'gCO'     , 'scale': 'log', 'cmap': cm.rainbow,
        'vr': (1e-10, 7.3e-5), 'hold_on':False},
    {'name': 'H2O'     , 'scale': 'log', 'cmap': cm.rainbow,
        'vr': (1e-18, 1.7e-4), 'hold_on': True},
    {'name': 'gH2O'    , 'scale': 'log', 'cmap': cm.rainbow,
        'vr': (1e-18, 1.7e-4), 'hold_on': True},
    {'name': 'OH'      , 'scale': 'log', 'cmap': cm.rainbow,
        'vr': (1e-10, 3e-5)},
    {'name': 'H2'      , 'scale': 'log', 'cmap': cm.rainbow,
        'vr': (1e-6, 1.0), 'hold_on': True},
    {'name': 'H'       , 'scale': 'log', 'cmap': cm.rainbow,
        'vr': (1e-6, 1), 'hold_on': False},
    {'name': 'H+'      , 'scale': 'log', 'cmap': cm.rainbow,
        'vr': (1e-16, 1e-1), 'hold_on': True},
    {'name': 'H2+'     , 'scale': 'log', 'cmap': cm.rainbow,
        'vr': (1e-6, 1.0), 'hold_on': True},
    {'name': 'H3+'     , 'scale': 'log', 'cmap': cm.rainbow,
        'vr': (1e-6, 1.0), 'hold_on': True},
    {'name': 'E-'      , 'scale': 'log', 'cmap': cm.rainbow,
        'vr': (1e-9, 7e-4)},
    {'name': 'Ncol'    , 'scale': 'log', 'cmap': cm.rainbow},
    {'name': 'Av'      , 'scale': 'log', 'cmap': cm.rainbow},
    {'name': 'n_gas'   , 'scale': 'log', 'cmap': cm.rainbow},
    {'name': 'f_H2'    , 'scale': 'log', 'cmap': cm.rainbow,
        'vr': (1e-6, 1e0)},
    {'name': 'f_H2O'   , 'scale': 'log', 'cmap': cm.rainbow},
    {'name': 'f_CO'    , 'scale': 'log', 'cmap': cm.rainbow},
  ]

idx = []
for i in xrange(filelen):
  rmin = dic['rmin'][i]
  rmax = dic['rmax'][i]
  if (rmin <= r_pos) and (rmax > r_pos):
    idx.append(i)

nidx = len(idx)

xaxisName = 'z (AU)'
#xaxisName = 'Ncol'
z = 0.5 * (dic['zmin'][idx] + dic['zmax'][idx])
#z = np.log10(dic['Ncol'][idx])
minz = np.min(z)
maxz = np.max(z)
xRange = (0, min(maxz*1.1, r_pos*0.6))
#xRange = (20, 22)#maxz)

hold_on_prev = False
icolor = 0

for item in name_list:
  name = item['name']
  if name in dic:
    v = dic[name][idx]
  else:
    print 'Not exist.'
    v = z * 1e-100

  maxval = np.nanmax(v)
  minval = np.nanmin(v)
  minval_nonzero = np.nanmin(v[np.where(v > 0.0)])
  print name, ' Max, min ', maxval, minval
  if 'vr' in item:
    if len(item['vr']) >= 2:
      minval = item['vr'][0]
      minval_nonzero = item['vr'][0]
      maxval = item['vr'][1]

  yRange = (minval_nonzero/1.5, maxval*2.5)
  
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

  if not hold_on_prev:
    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes([x_start, y_start, del_x_1, del_y],
      xlabel=xaxisName, ylabel='', title='r = {0:4.1f} AU'.format(r_pos),
      autoscalex_on=False, autoscaley_on=False,
      xscale='linear', yscale=item['scale'],
      xlim=xRange, ylim=yRange)
    ax.xaxis.label.set_fontsize(25)
    ax.yaxis.label.set_fontsize(25)
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.xaxis.grid(which='major', color=(0.7,0.8,0.7), linewidth=1,   linestyle='-')
    ax.xaxis.grid(which='minor', color=(0.7,0.9,0.7), linewidth=0.2, linestyle='-')
    ax.yaxis.grid(which='major', color=(0.7,0.8,0.7), linewidth=1,   linestyle='-')
    ax.yaxis.grid(which='minor', color=(0.7,0.9,0.7), linewidth=0.2, linestyle='-')
    ax.set_axisbelow(True)
    for label in ax.get_xticklabels():
      label.set_fontsize(20)
    for label in ax.get_yticklabels():
      label.set_fontsize(20)

  zsorted, vsorted = zip(*sorted(zip(z, v)))
  f = interpolate.interp1d(zsorted, vsorted)
  znew = np.linspace(minz, maxz, 50)
  ax.plot(znew, f(znew),
      linestyle='-', label=name,
      color=color_list[icolor%len(color_list)], linewidth=5)
  #ax.plot(z, v,
  #    linestyle='None', marker='o', markersize=10,
  #    color=color_list[icolor%len(color_list)], linewidth=1)

  if not 'hold_on' in item:
    hold_on_this = False
  else:
    hold_on_this = item['hold_on']
  
  if not hold_on_this:
    lgd = ax.legend(loc='lower left', bbox_to_anchor=(0.7, 0.65), prop={'size':15},
      fancybox=False, shadow=False, ncol=1)
    plt.savefig(pp, format='pdf')
    icolor = 0
  else:
    icolor += 1

  hold_on_prev = hold_on_this

pp.close()
