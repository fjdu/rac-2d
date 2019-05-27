from matplotlib import *
use('Agg')
from matplotlib.pyplot import *
from numpy import *

from os.path import join as opj
from glob import glob

from parse_ana import *
from long_function_definitions import *
from my_script import *

import drawing

mycm = make_my_colormap(c_list=[(0.0, (0.5, 0.0, 0.5)),
                                (0.2, (0.0, 0.0, 1.0)),
                                (0.4, (0.0, 0.8, 1.0)),
                                (0.6, (0.0, 0.8, 0.0)),
                                (0.8, (1.0, 0.8, 0.0)),
                                (1.0, (1.0, 0.0, 0.0))])

rcParams['axes.color_cycle'] = mycolors

from scipy.interpolate import griddata
seterr(divide='ignore')

def plot_contribution_function(c, fig_fname,
                               nn=200, xmin=0, xmax=20, ymin=-10, ymax=10, nlev=100,
                               thr=1e-5, normalize=True,
                               return_data = False,
                               xscale='linear', yscale='linear',
                               figsize=(8,8)):
    #
    r = sqrt(c[:,0]**2 + c[:,1]**2)
    #
    xi = linspace(xmin, xmax, nn)
    yi = linspace(ymin, ymax, nn)
    #
    xi, yi = np.meshgrid(xi, yi)
    maxval = c[:,3].max()
    c[c[:,3] < maxval*thr, 3] = maxval*thr #0.0
    #
    zi = griddata((r, c[:,2]), log10(c[:,3]), (xi, yi), method='linear')
    if normalize:
        zi = zi - np.nanmax(zi)
    #
    ##xi = 1.15 * xi
    #
    f = figure(figsize=figsize)
    pos = (0.15, 0.15, 0.8, 0.4)
    ax = f.add_axes(pos,
          xlabel='r (AU)',
          ylabel='z (AU)',
          autoscalex_on=False, autoscaley_on=False,
          xscale=xscale, yscale=yscale,
          xlim=(xmin, xmax), ylim=(ymin, ymax))
    #
    C = ax.contourf(xi, yi, zi, nlev, cmap=mycm, extend='neither')
    for _ in C.collections:
        _.set_rasterized(True)
    #colorbar(C)
    #set_axis_format(ax, graygrid=True, majorgridon=True, minorgridon=False)
    pos = (0.15, 0.6, 0.8, 0.03)
    cax = f.add_axes(pos)
    drawing.add_colorbar(cax, np.nanmin(zi), np.nanmax(zi), mycm,
                 label='', title='log$_{10}$ (Relative contributions)',
                 orientation='horizontal', scale='linear')
    #
    flx = xi * 10**zi
    idx = logical_not(isfinite(zi))
    flx[idx] = 0.0
    flx = sum(flx, axis=0)
    flx_t = sum(flx)
    flx = flx / flx_t
    pos = (0.15, 0.68, 0.8, 0.33)
    #maxval = flx.max()
    #maxtik = ceil(maxval/0.02)*0.02
    #yticks = linspace(0,maxtik,num=int(maxtik/0.02)+1)
    ax = f.add_axes(pos,
          xlabel='', ylabel='',
          xticklabels=[],
          #yticks=yticks,
          autoscalex_on=False, autoscaley_on=False,
          xscale=xscale, yscale=yscale,
          xlim=(xmin, xmax), ylim=(-0.1,1.1))
    x_, y_ = loglin_interpol(xi[0,:], flx, num=100, method='linear',
                            dosmooth=True, winwidth=3)
    y_ = y_ / np.nanmax(y_)
    ax.plot(x_, y_, lw=2, color='red', label='Contrib.')
    fsum = cumsum(y_)
    fsum = fsum/np.nanmax(fsum)
    #set_axis_format(ax, majorgridon=False, minorgridon=False)
    #
    ax.plot(x_, fsum, lw=1, color='blue', label='Accum. Contrib.')
    ax.legend(loc="center right", fontsize=15)
    #
    savefig(fig_fname, bbox_inches='tight')
    #
    if return_data:
        return xi, yi, zi
    else:
        return

def plot_contribution_together(c_s, labels, fig_fname,
                               nn=200, xmin=0, xmax=20, ymin=-10, ymax=10, nlev=100,
                               thr=1e-5, normalize=True,
                               return_data = False,
                               xscale='linear', yscale='linear',
                               figsize=(8,5)):
    #
    f = figure(figsize=figsize)
    pos = (0.15, 0.15, 0.8, 0.8)
    #
    for icount, c in enumerate(c_s):
        r = sqrt(c[:,0]**2 + c[:,1]**2)
        #
        xi = linspace(xmin, xmax, nn)
        yi = linspace(ymin, ymax, nn)
        #
        xi, yi = np.meshgrid(xi, yi)
        maxval = c[:,3].max()
        c[c[:,3] < maxval*thr, 3] = 0.0
        if normalize:
            c[:,3] = c[:,3] / maxval
        #
        zi = griddata((r, c[:,2]), log10(c[:,3]), (xi, yi), method='linear')
        #
        idx = logical_not(isfinite(zi))
        zi[idx] = -99.0
        flx = sum(xi * 10**zi, axis=0)
        flx_t = sum(flx)
        flx = flx / flx_t
        fsum = cumsum(flx)
        #maxval = flx.max()
        #maxtik = ceil(maxval/0.1)*0.1
        #yticks = linspace(0,maxtik,num=int(maxtik/0.1)+1)
        if icount == 0:
            ax = f.add_axes(pos,
                  xlabel='r (AU)',
                  ylabel='Fraction (within r)',
                  autoscalex_on=False, autoscaley_on=False,
                  xscale=xscale, yscale=yscale,
                  xlim=(xmin, xmax), ylim=(0, 1.05))
        ax.plot(xi[0], fsum, lw=2, label=labels[icount])
    #
    legend(loc="upper left", fontsize=20)
    set_axis_format(ax, majorgridon=False, minorgridon=False)
    #
    savefig(fig_fname, bbox_inches='tight')
    #
    return

if __name__ == '__main__':
    fig_dir = '/n/Users/fdu/now/'
    #res_dir = '/n/Fdu1/work/grid_20150703_redo/results/'
    res_dir = '/n/Users/fdu/res/'
    #
   #models = ['run_011']
    models = [
                #'20150609_C2H_5a_4a_2d_8cl_B6l_41_h',
                '20160714_TWH_d',
                '20160714_TWH_e',
                '20160714_TWH_f',
                #'20151117_a_dep1a',
                #'20151117_a_dep2',
             ]
    #
    lines_contri = [ \
        {
         'dir': 'C2H_2_NLTE',
         'postfix': 'C2H_2_NLTE_contri',
         'file': 'line_00003_00032_00001_3.49339E+11_7.00.fits_contri.dat'},
         #'dir': 'C18O_contri',
         #'postfix': 'C18O_contri',
         #'file': 'line_00003_00003_00001_3.45798E+11_7.00.fits_contri.dat'},
        ]
    #
    for model_name in models:
        for line in lines_contri:
            line_dir = line['dir']
            postfix = line['postfix']
            data_file = line['file']
            #
            model_dir = opj(res_dir, model_name)
            data_dir = opj(model_dir, line_dir, 'images')
            fpath = opj(data_dir, data_file)
            if not os.path.exists(fpath):
                continue

            print 'Running ', fpath
            #
            c = loadtxt(fpath)
            fig_fname = opj(fig_dir, model_name + '_' + line_dir + postfix + '.pdf')
            plot_contribution_function(c, fig_fname, thr=1e-15,
                                       xmin=-2, xmax=200, ymin=-80, ymax=80)
            print fig_fname, 'saved.'
#       c_s = []
#       for line in lines_contri:
#           line_dir = line['dir']
#           postfix = line['postfix']
#           data_file = line['file']
#           #
#           model_dir = opj(res_dir, model_name)
#           data_dir = opj(model_dir, line_dir, 'images')
#           fpath = opj(data_dir, data_file)
#           print 'Running ', fpath
#           #
#           c = loadtxt(fpath)
#           c_s.append(c)
#       fig_fname = opj(fig_dir, model_name + '_accu_all.pdf')
#       plot_contribution_together(c_s,
#                                  [r'$1_{10}\rightarrow 1_{01}$', r'$3_{12}\rightarrow 2_{21}$'],
#                                  fig_fname, thr=1e-18,
#                                  xmin=-2, xmax=400, ymin=-100, ymax=200)
