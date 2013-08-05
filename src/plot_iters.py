import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.patches import Polygon
from matplotlib.mlab import griddata
from matplotlib.backends.backend_pdf import PdfPages
import os, os.path

data_dir = 'results/results_20130617_nirgal_f/'

name_list = \
  [ \
    {'name': 'Tgas'  , 'linearscale': False, 'cmap': cm.coolwarm},#, 'vr': (1e2,   2e4)},
    #{'name': 'n_gas' , 'linearscale': False, 'cmap': cm.BuPu  },
    {'name': 'H2O'   , 'linearscale': False, 'cmap': cm.Blues_r },#, 'vr': (1e-10, 1e-4)},
    {'name': 'OH'    , 'linearscale': False, 'cmap': cm.Blues},
    #{'name': 'H2'    , 'linearscale': False, 'cmap': cm.Greys},
    #{'name': 'H'     , 'linearscale': False, 'cmap': cm.Greys},
    #{'name': 'E-'    , 'linearscale': False, 'cmap': cm.Greys},
    #{'name': 'C'     , 'linearscale': False, 'cmap': cm.Greys},
    #{'name': 'C+'    , 'linearscale': False, 'cmap': cm.Greys},
    #{'name': 'CO'    , 'linearscale': False, 'cmap': cm.Greys},
    {'name': 'UV_G0'  , 'linearscale': False, 'cmap': cm.Greys},
    {'name': 'f_H2'  , 'linearscale': False, 'cmap': cm.BuGn},
    {'name': 'f_H2O' , 'linearscale': False, 'cmap': cm.BuGn},
  ]

n_files = len([name for name in os.listdir(data_dir) if name[0:4] == 'iter' and name[-3:] == 'dat'])
#os.path.isfile(os.path.join(data_dir, name))])
print n_files

for item in name_list:

  name = item['name']

  pp = PdfPages(os.path.join(data_dir, 'plots_iter_' + name + '.pdf'))

  for i in xrange(n_files):

    filename_save_results = data_dir + 'iter_{0:04d}.dat'.format(i+1)

    data = np.loadtxt(filename_save_results, comments='!')

    filelen = len(data[:,0])

    f = open(filename_save_results, 'r')
    str_comment = f.readline()[1:].split()
    f.close()

    dic = {}
    for j in xrange(len(str_comment)):
      dic.update({str_comment[j]: data[:, j]})

    use_linear_scale = item['linearscale']

    minr = np.min(dic['rmin'])
    maxr = np.max(dic['rmax'])
    minz = np.min(dic['zmin'])
    maxz = np.max(dic['zmax'])
    xRange = (minr/1.2, maxr*1.1)
    yRange = (minz/1.2, maxz*1.5)

    colormap = item['cmap']

    maxval = np.nanmax(dic[name])
    minval = np.nanmin(dic[name])
    minval_nonzero = np.nanmin(dic[name][np.where(dic[name] > 0.0)])
    print name, i, ' Max, min ', maxval, minval
    if 'vr' in item:
      if len(item['vr']) >= 2:
        minval = item['vr'][0]
        minval_nonzero = item['vr'][0]
        maxval = item['vr'][1]
    log_max = np.log(maxval)
    log_min = np.log(max(minval, minval_nonzero, maxval/1E10))
    minval = np.exp(log_min)

    xlen = xRange[1] - xRange[0]
    ylen = yRange[1] - yRange[0]
    figsize_x = 8
    x_start = 0.13
    y_start = 0.1
    x_end = 0.9
    y_end = 0.93
    x_sep = 0.06
    y_sep = 0.08
    del_x = x_end - x_start - x_sep
    del_x_1 = 0.85 * del_x
    del_x_2 = 0.06 * del_x
    del_y = y_end - y_start

    figsize_y = figsize_x * (ylen / xlen) * del_x_1/del_y
    figsize = (figsize_x, figsize_y)

    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes([x_start, y_start, del_x_1, del_y],
      xlabel='r (AU)', ylabel='z (AU)', title=name + ' iter_{0:04d}'.format(i+1),
      autoscalex_on=False, autoscaley_on=False,
      xlim=xRange, ylim=yRange)
    ax.xaxis.label.set_fontsize(15)
    ax.yaxis.label.set_fontsize(15)

    for j in xrange(filelen):
      x1 = dic['rmin'][j]
      y1 = dic['zmin'][j]
      x2 = dic['rmax'][j]
      y2 = dic['zmax'][j]
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
      val = dic[name][j]
      if use_linear_scale:
        sca_col = (val - minval) / (maxval - minval)
      else:
        if val <= 0.0:
          sca_col = np.nan
        else:
          sca_col = (np.log(val) - log_min) / (log_max-log_min)
      thiscolor = colormap(sca_col)
      poly = Polygon(pxy, closed=True,
        facecolor=thiscolor, edgecolor=thiscolor, linewidth=0.000001)
      ax.add_patch(poly)

    ax1 = fig.add_axes([x_start+del_x_1+x_sep, y_start, del_x_2, del_y])
    cbar = matplotlib.colorbar.ColorbarBase(ax1, ticks=[0, 0.5, 1], cmap=colormap, orientation='vertical')
    if use_linear_scale:
      mid_val = (maxval+minval)*0.5
    else:
      mid_val = np.sqrt(maxval*minval)
    cbar.set_ticklabels(['{0:7.1e}'.format(minval),
      '{0:7.1e}'.format(mid_val), '{0:7.1e}'.format(maxval)])

    plt.savefig(pp, format='pdf')

  print 'Saving pdf...'
  pp.close()


