import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.patches import Polygon
#from matplotlib.mlab import griddata
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.collections import PolyCollection
import os
import os.path
import my_script as my
from scipy.interpolate import griddata

data_dir = './results/results_20130726_gondolin_a_1/'
filename_save_results =  os.path.join(data_dir, 'iter_0004.dat')

data = np.loadtxt(filename_save_results, comments='!')

filelen = len(data[:,0])

f = open(filename_save_results, 'r')
str_comment = f.readline()[1:].split()
f.close()

dic = {}
for i in xrange(len(str_comment)):
  dic.update({str_comment[i]: data[:, i]})

dic.update({'Pressure': dic['Tgas']*dic['n_gas']})

#pp = PdfPages(os.path.join(data_dir, 'plots.pdf'))

name_list = \
  [ \
    {'name': 'Tgas'    , 'linearscale': False, 'cmap': cm.rainbow,  'vr': (50, 1e4)},
    {'name': 'Tdust'   , 'linearscale': False, 'cmap': cm.rainbow,  'vr': (50, 5e2)},
    {'name': 'n_gas'   , 'linearscale': False, 'cmap': cm.rainbow,  'vr': (40, 1.5e10)},
    {'name': 'H2O'     , 'linearscale': False, 'cmap': cm.rainbow,  'vr': (1e-14, 1.7e-4)},
    {'name': 'OH'      , 'linearscale': False, 'cmap': cm.rainbow,  'vr': (1e-10, 3e-5)},
    #{'name': 'H2'      , 'linearscale': False, 'cmap': cm.rainbow,  'vr': (1e-6, 0.5)},
    #{'name': 'H'       , 'linearscale': False, 'cmap': cm.rainbow,  'vr': (1e-6, 1)},
    #{'name': 'E-'      , 'linearscale': False, 'cmap': cm.rainbow,  'vr': (1e-9, 7e-4)},
    #{'name': 'C'       , 'linearscale': False, 'cmap': cm.rainbow,  'vr': (1e-10, 7e-5)},
    #{'name': 'C+'      , 'linearscale': False, 'cmap': cm.rainbow,  'vr': (1e-10, 7e-5)},
    #{'name': 'CO'      , 'linearscale': False, 'cmap': cm.rainbow,  'vr': (1e-10, 7.3e-5)},
    #{'name': 'f_H2'    , 'linearscale': False, 'cmap': cm.rainbow,  'vr': (1e-6, 1e0)},
    #{'name': 'f_H2O'   , 'linearscale': False, 'cmap': cm.rainbow,  'vr': (1e-1, 1e0)},
  ]

#for i in xrange(2, len(name_list)):
#  dic[name_list[i]['name']] *= dic['n_gas']

font = {'family' : 'normal',
        #'weight' : 'bold',
        'size'   : 25}
matplotlib.rc('font', **font)

x = 0.5 * (dic['rmin'] + dic['rmax'])
y = 0.5 * (dic['zmin'] + dic['zmax'])

for item in name_list:
  name = item['name']
  #use_linear_scale = True
  use_linear_scale = item['linearscale']

  minr = np.min(dic['rmin'])
  maxr = np.max(dic['rmax'])
  minz = np.min(dic['zmin'])
  maxz = np.max(dic['zmax'])
  xRange = (minr/1.5, min(maxr*1.1, 30))
  yRange = (minz/1.5, min(maxz*1.5, xRange[1]*0.8))
  
  colormap = item['cmap']
  
  maxval = np.nanmax(dic[name])
  minval = np.nanmin(dic[name])
  minval_nonzero = np.nanmin(dic[name][np.where(dic[name] > 0.0)])
  print name, ' Max, min ', maxval, minval
  if 'vr' in item:
    if len(item['vr']) >= 2:
      minval = item['vr'][0]
      minval_nonzero = item['vr'][0]
      maxval = item['vr'][1]
  log_max = np.log10(maxval)
  log_min = np.log10(max(minval, minval_nonzero, maxval/1E10))
  #minval = np.exp(log_min)
  
  xlen = xRange[1] - xRange[0]
  ylen = yRange[1] - yRange[0]
  figsize_x = 10
  x_start = 0.1
  y_start = 0.12
  x_end = 0.9
  y_end = 0.90
  x_sep = 0.04
  y_sep = 0.08
  del_x = x_end - x_start - x_sep
  del_x_1 = 0.85 * del_x
  del_x_2 = 0.04 * del_x
  del_y = y_end - y_start

  ntick_cbar = 5
  
  figsize_y = figsize_x * (ylen / xlen) * del_x_1/del_y
  figsize = (figsize_x, figsize_y)
  
  fig = plt.figure(figsize=figsize)
  ax = fig.add_axes([x_start, y_start, del_x_1, del_y],
    xlabel='r (AU)', ylabel='z (AU)', title=name,
    autoscalex_on=False, autoscaley_on=False,
    xscale='linear', yscale='linear',
    xlim=xRange, ylim=yRange)

  x_s = 0.5 * (dic['rmin'] + dic['rmax'])
  y_s = 0.5 * (dic['zmin'] + dic['zmax'])
  #z_s = np.log(dic[name])
  z_s = dic[name]
  
  xx, yy = np.mgrid[xRange[0]:xRange[1]:400j, yRange[0]:yRange[1]:400j]
  zz = griddata((x_s, y_s), z_s, (xx, yy))

  plt.imshow(np.log10(zz.T), extent=xRange + yRange, vmin=log_min, vmax=log_max, cmap=colormap, origin='lower', interpolation='bilinear')
  #plt.contourf(xx, yy, np.log10(zz), 160, vmin=log_min, vmax=log_max, cmap = colormap, extend='both', origin='lower')
  #plt.pcolormesh(xx, yy, zz)
  
  #xi = np.linspace(xRange[0], xRange[1], 100)
  #yi = np.linspace(yRange[0], yRange[1], 100)
  #xi, yi = np.meshgrid(xi, yi)
  #z = np.log(dic[name])
  #zi = griddata((x, y), z, (xi, yi), method='linear')
  #zi[yi/xi > 0.75] = np.nan
  #CS = ax.contour(xi, yi, zi, 6, colors='white', linestyle='-')

  ax1 = fig.add_axes([x_start+del_x_1+x_sep, y_start, del_x_2, del_y])
  ticks = np.linspace(0.0, 1.0, ntick_cbar)
  if use_linear_scale:
    tickvals = np.linspace(minval, maxval, num=ntick_cbar)
  else:
    tickvals = np.logspace(log_min, log_max, num=ntick_cbar)
  ticklabels = [my.num2str4fig(tickval) for tickval in tickvals]
  cbar = matplotlib.colorbar.ColorbarBase(ax1, ticks=ticks, cmap=colormap, orientation='vertical')
  cbar.set_ticklabels(ticklabels)

  plt.savefig(os.path.join(data_dir, name + '_imshow.eps'))
  #plt.savefig(pp, format='pdf')

#pp.close()
