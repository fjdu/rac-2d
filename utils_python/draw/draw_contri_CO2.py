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

from cycler import cycler
rcParams['axes.prop_cycle'] = cycler(color=mycolors)

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
                               show_scaling=False,
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
    if show_scaling:
        coeff = fsum[-1] / (xi[0][-1])**2
        ax.plot(xi[0], coeff * xi[0]**2, lw=2, linestyle='--', label='$r^2$ scaling')
    #
    legend(loc="upper left", fontsize=20)
    set_axis_format(ax, majorgridon=False, minorgridon=False)
    #
    savefig(fig_fname, bbox_inches='tight')
    #
    return

if __name__ == '__main__':
    fig_dir = '/Users/fjdu/work_local/upgrade20180626/figures/'
    #
    #c = loadtxt('/Users/fjdu/work_local/upgrade20180626/storage/cjmerch/20150307_dep_p1/CO2_0deg/images/line_00001_10941_00001_1.93049E+13_0.00.fits_contri.dat')
    #c = loadtxt('/Users/fjdu/work_local/upgrade20180626/storage/cjmerch/20150307_dep_p1/CO2/images/line_00001_10941_00001_1.93049E+13_7.00.fits_contri.dat')
    #c = loadtxt('/Users/fjdu/work_local/upgrade20180626/storage/cjmerch/20150307_dep_p1/C2H/images/line_00001_00018_00001_2.62004E+11_7.00.fits_contri.dat')
    c = loadtxt('/Users/fjdu/work_local/upgrade20180626/storage/cjmerch/20150307_dep_p1/CO2_14.97/images/line_00010_13905_00001_2.00310E+13_7.00.fits_contri.dat')
    fig_fname = opj(fig_dir, 'contri_CO2_14.97.pdf')
    plot_contribution_function(c, fig_fname, xmin=-2, xmax=100, ymin=-100, ymax=100, thr=1e-10)
