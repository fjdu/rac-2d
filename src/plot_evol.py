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

data_dir = './'
filename_save_results =  os.path.join(data_dir, 'chem_evol_tmp.dat')

data = np.loadtxt(filename_save_results, comments='!')

filelen = len(data[:,0])

f = open(filename_save_results, 'r')
str_comment = f.readline()[1:].split()
f.close()

dic = {}
for i in xrange(len(str_comment)):
  dic.update({str_comment[i]: data[:, i]})

pp = PdfPages(os.path.join(data_dir, 'plot_evol.pdf'))

name_list = [{'name': s, 'scale': 'log', 'cmap': cm.rainbow} for s in dic.keys()]

z = dic['Time']
minz = np.min(z)
maxz = np.max(z)

i = 1
for item in name_list:
  name = item['name']
  v = dic[name]

  maxval = np.nanmax(v)
  minval = np.nanmin(v)
  minval_nonzero = np.nanmin(v[np.where(v > 0.0)])
  print name, ' Max, min ', maxval, minval
  if 'vr' in item:
    if len(item['vr']) >= 2:
      minval = item['vr'][0]
      minval_nonzero = item['vr'][0]
      maxval = item['vr'][1]

  xRange = (maxz/100.0, maxz)
  yRange = (minval/1.5, maxval*1.5)
  
  xlen = xRange[1] - xRange[0]
  ylen = xlen*0.62#yRange[1] - yRange[0]
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

  figsize_y = figsize_x * (ylen / xlen) * del_x_1/del_y
  figsize = (figsize_x, figsize_y)
  
  fig = plt.figure(figsize=figsize)
  ax = fig.add_axes([x_start, y_start, del_x_1, del_y],
    xlabel='Time (yr)', ylabel=name,
    autoscalex_on=False, autoscaley_on=False,
    xscale='log', yscale=item['scale'],
    xlim=xRange, ylim=yRange)
  ax.xaxis.label.set_fontsize(25)
  ax.yaxis.label.set_fontsize(25)
  for label in ax.get_xticklabels():
    label.set_fontsize(20)
  for label in ax.get_yticklabels():
    label.set_fontsize(20)

  ax.plot(z, v,
      linestyle='-',
      marker = 'd',
      color='blue', linewidth=5)
  
  plt.savefig(pp, format='pdf')

pp.close()
