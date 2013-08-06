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

font = {'family' : 'normal',
        #'weight' : 'bold',
        'size'   : 25}
matplotlib.rc('font', **font)

r_pos = 8.0

data_dir = '/Users/fdu/work/protoplanetary_disk/res/results_20130805_gondolin_a_3/'
filename_save_results =  os.path.join(data_dir, 'iter_0001.dat')

data = np.loadtxt(filename_save_results, comments='!')

filelen = len(data[:,0])

f = open(filename_save_results, 'r')
str_comment = f.readline()[1:].split()
f.close()

dic = {}
for i in xrange(len(str_comment)):
  dic.update({str_comment[i]: data[:, i]})

h_gg_co_name = 'h_gg_co'
str_comment.append(h_gg_co_name)
dic.update({h_gg_co_name: -1.0*dic['c_gg_co']})

nlines = 5

#pp = PdfPages(os.path.join(data_dir, 'plot_column_h_c.pdf'))

name_list = []

names_heating = []
names_cooling = []

for item in str_comment:
  if item[0:2] == 'h_':
    names_heating.append({'name': item})
  if item[0:2] == 'c_':
    names_cooling.append({'name': item})

name_lists = [names_heating, names_cooling]

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

color_list = ['blue', 'red', 'green', 'magenta', (0.7,0.5,0.3), (0.8,0.85,0.3)]

h_c_names_dic = \
{
'h_ph_gr': 'photo-elec',
'h_fo_H2': 'H2-form',
'h_cosmi': 'cosm-ray',
'h_vi_H2': 'H2-vib',
'h_io_CI': 'C-ioniz',
'h_ph_H2': 'phdis-H2',
'h_ph_wa': 'phdis-H2O',
'h_ph_OH': 'phdis-OH',
'h_Xray': 'Xray',
'h_visco': 'viscosity',
'h_gg_co': 'gas-dust',
'c_el_gr': 'elec-recom',
'c_vi_H2': 'H2-vib',
'c_gg_co': 'gas-dust',
'c_OI': 'OI',
'c_CII': 'CII',
'c_wa_ro': 'H2O-rot',
'c_wa_vi': 'H2O-vib',
'c_CO_ro': 'CO-rot',
'c_CO_vi': 'CO-vib',
'c_H2_ro': 'H2-rot',
'c_LyAlp': 'LyA',
'c_fb': 'free-b',
'c_ff': 'free-f',
}
  
for name_list in name_lists:
  maxvals = []
  for item in name_list:
    name = item['name']
    item.update({'max': np.max(dic[name][idx])})
    maxvals.append(item['max'])
  maxvals.sort(reverse=True)
  
  nth = min(len(maxvals)-1, nlines)
  xRange = (0, r_pos*0.6)
  yRange = (maxvals[nth]/100, maxvals[0]*1.1)
  
  xlen = xRange[1] - xRange[0]
  ylen = xlen * 0.9
  figsize_x = 10
  x_start = 0.17
  y_start = 0.1
  x_end = 0.98
  y_end = 0.90
  x_sep = 0.01
  y_sep = 0.08
  del_x = x_end - x_start - x_sep
  del_x_1 = 0.95 * del_x
  del_x_2 = 0.06 * del_x
  del_y = y_end - y_start
  
  figsize_y = figsize_x * (ylen / xlen) * del_x_1/del_y
  figsize = (figsize_x, figsize_y)
  
  fig = plt.figure(figsize=figsize)
  ax = fig.add_axes([x_start, y_start, del_x_1, del_y],
    xlabel='z (AU)', ylabel='Rate (erg cm$^{-3}$ s$^{-1}$)',
    title = 'r = {0:3.1f} AU'.format(r_pos),
    autoscalex_on=False, autoscaley_on=False,
    xscale='linear', yscale='log',
    xlim=xRange, ylim=yRange)
  #fig.suptitle('test title', fontsize=20)
  #ax.xaxis.label.set_fontsize(25)
  #ax.yaxis.label.set_fontsize(25)
  #ax.title.set_fontsize(20)
  #for label in ax.get_xticklabels():
  #  label.set_fontsize(20)
  #for label in ax.get_yticklabels():
  #  label.set_fontsize(20)
  
  i = 0
  for item in name_list:
    name = item['name']
    if item['max'] < maxvals[nth]:
      continue
  
    v = dic[name][idx]
    f = interpolate.interp1d(z, v)
    znew = np.linspace(z[0], z[-1], 100)
  
    ax.plot(z, v, linestyle='-', label=name,
    #ax.plot(znew, f(znew), linestyle='-', label=h_c_names_dic[name],
        color=color_list[i%(nlines+1)], linewidth=5)
    i += 1 
  
  lgd = ax.legend(loc='lower left', bbox_to_anchor=(0.0, 0.69), prop={'size':15},
    fancybox=False, shadow=False, ncol=1)

  if name[0] == 'h':
    filename = 'heatings.eps'
  if name[0] == 'c':
    filename = 'coolings.eps'
  plt.savefig(os.path.join(data_dir, filename))
  #plt.savefig(pp, format='pdf')

#pp.close()
