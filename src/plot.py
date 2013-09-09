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
from scipy.interpolate import griddata

data_dir = '/Users/fdu/work/protoplanetary_disk/res/results_20130905_gondolin_a_2a/'
filename_save_results =  os.path.join(data_dir, 'iter_0002.dat')

data = np.loadtxt(filename_save_results, comments='!')

filelen = len(data[:,0])

f = open(filename_save_results, 'r')
str_comment = f.readline()[1:].split()
f.close()

dic = {}
for i in xrange(len(str_comment)):
  dic.update({str_comment[i]: data[:, i]})

dic.update({'Pressure': dic['Tgas']*dic['n_gas']})
dic.update({'Tgas_Tdust': dic['Tgas']/dic['Tdust']})

#pp = PdfPages(os.path.join(data_dir, 'plots.pdf'))

draw_rectangles = True

name_list = \
  [ \
    {'name': 'Tgas'    , 'linearscale': False, 'cmap': cm.rainbow},
    {'name': 'Tdust'   , 'linearscale': False, 'cmap': cm.rainbow},
    {'name': 'Tgas_Tdust', 'linearscale': False, 'cmap': cm.rainbow, 'levels': (0.1,1.0,10.0)},
    {'name': 'n_gas'   , 'linearscale': False, 'cmap': cm.rainbow},
    {'name': 'H2O'     , 'linearscale': False, 'cmap': cm.rainbow, 'levels': (1e-8, 1e-7, 1e-6, 1e-5)},
    {'name': 'OH'      , 'linearscale': False, 'cmap': cm.rainbow},
    {'name': 'H2'      , 'linearscale': False, 'cmap': cm.rainbow,  'vr': (1e-6, 0.5)},
    #{'name': 'H'       , 'linearscale': False, 'cmap': cm.rainbow,  'vr': (1e-6, 1)},
    #{'name': 'E-'      , 'linearscale': False, 'cmap': cm.rainbow,  'vr': (1e-9, 7e-4)},
    #{'name': 'C'       , 'linearscale': False, 'cmap': cm.rainbow,  'vr': (1e-10, 7e-5)},
    #{'name': 'C+'      , 'linearscale': False, 'cmap': cm.rainbow,  'vr': (1e-10, 7e-5)},
    {'name': 'CO'      , 'linearscale': False, 'cmap': cm.rainbow,  'vr': (1e-10, 7.3e-5)},
    {'name': 'gCO'     , 'linearscale': False, 'cmap': cm.rainbow,  'vr': (1e-10, 7.3e-5)},
    {'name': 'gH2O', 'linearscale': False, 'cmap': cm.rainbow, 'levels': (1e-8, 1e-7, 1e-6, 1e-5)},
    #{'name': 'f_H2'    , 'linearscale': False, 'cmap': cm.rainbow,  'vr': (1e-6, 1e0)},
    #{'name': 'f_H2O'   , 'linearscale': False, 'cmap': cm.rainbow,  'vr': (1e-1, 1e0)},
  ]

#name_list = \
#  [ \
#    {'name': 'gH2O', 'linearscale': False, 'cmap': cm.rainbow, 'levels': (1e-8, 1e-7, 1e-6, 1e-5)},
#  ]

#for i in xrange(2, len(name_list)):
#  dic[name_list[i]['name']] *= dic['n_gas']

font = {'family' : 'normal',
        #'weight' : 'bold',
        'size'   : 25}
matplotlib.rc('font', **font)

x = dic['rmin'] #0.5 * (dic['rmin'] + dic['rmax'])
y = dic['zmin'] #0.5 * (dic['zmin'] + dic['zmax'])

for item in name_list:
  name = item['name']
  #use_linear_scale = True
  use_linear_scale = item['linearscale']

  minr = np.min(dic['rmin'])
  maxr = np.max(dic['rmax'])
  minz = np.min(dic['zmin'])
  maxz = np.max(dic['zmax'])
  #xRange = (minr/1.5, min(maxr*1.1, 50))
  #yRange = (minz/1.5, min(maxz*1.5, xRange[1]*0.8))
  xRange = (0.5, 50)
  yRange = (0, 30)
  
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
  minval = np.exp(log_min)
  
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
  #ax.xaxis.label.set_fontsize(30)
  #ax.yaxis.label.set_fontsize(30)
  matplotlib.rcParams['axes.linewidth'] = 0.001
  
  if draw_rectangles:
    for i in xrange(filelen):
      x1 = dic['rmin'][i]
      y1 = dic['zmin'][i]
      x2 = dic['rmax'][i]
      y2 = dic['zmax'][i]
      ##if y2/x2 > 0.75:
      ##  continue
      if x1 > xRange[1] or y1 > yRange[1]:
        continue
      if x2 < xRange[0] or y2 < yRange[0]:
        continue
      pxy = np.zeros((5,2))
      pxy[0, :] = [x1, y1]
      pxy[1, :] = [x2, y1]
      pxy[2, :] = [x2, y2]
      pxy[3, :] = [x1, y2]
      pxy[4, :] = [x1, y1]
      val = dic[name][i]
      if use_linear_scale:
        sca_col = (val - minval) / (maxval - minval)
      else:
        if val <= 0.0:
          sca_col = np.nan
        else:
          sca_col = (np.log10(val) - log_min) / (log_max-log_min)
      thiscolor = colormap(sca_col)
      poly = Polygon(pxy, closed=True,
        facecolor=thiscolor, edgecolor=thiscolor, linewidth=0.000001)
      ax.add_patch(poly)
  else:
    xi = np.linspace(xRange[0], xRange[1], 200)
    yi = np.linspace(yRange[0], yRange[1], 200)
    xi, yi = np.meshgrid(xi, yi)
    z = np.log10(dic[name])
    zi = griddata((x, y), z, (xi, yi), method='linear')
    #zi[yi/xi > 0.75] = np.nan
    ax.contourf(xi, yi, zi, 100, cmap=colormap, vmin=log_min, vmax=log_max)
    #if 'levels' in item:
    #  CS = ax.contour(xi, yi, zi, levels=np.log10(item['levels']), linewidth=1,
    #      colors='white', vmin=log_min, vmax=log_max)
    #else:
    #  CS = ax.contour(xi, yi, zi, 6, linewidth=1, colors='white', vmin=log_min, vmax=log_max)
    #plt.clabel(CS, inline=1, fontsize=5)

  ax1 = fig.add_axes([x_start+del_x_1+x_sep, y_start, del_x_2, del_y])
  ticks = np.linspace(0.0, 1.0, ntick_cbar)
  if use_linear_scale:
    tickvals = np.linspace(minval, maxval, num=ntick_cbar)
  else:
    tickvals = np.logspace(log_min, log_max, num=ntick_cbar)
  ticklabels = [my.num2str4fig(tickval) for tickval in tickvals]
  cbar = matplotlib.colorbar.ColorbarBase(ax1, ticks=ticks, cmap=colormap, orientation='vertical')
  cbar.set_ticklabels(ticklabels)

  plt.savefig(os.path.join(data_dir, name + '.eps'))
  #plt.savefig(pp, format='pdf')

#pp.close()
