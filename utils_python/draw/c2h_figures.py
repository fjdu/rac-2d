import os
import sys
sys.path.append('/Users/fjdu/Dropbox/myCodes/my_python_lib/')
from draw.long_function_definitions import set_axis_format
from draw.fits_image import cube, combine_two_cubes
from scipy.ndimage.interpolation import rotate

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import astropy
from astropy.io import fits
from astropy.wcs import WCS

fits_dir = '/n/Users/fdu/now/C2H_obs_data/'
fits_dir_CO = '/n/Users/fdu/now/CO_obs_data/'
res_dir = '/n/Users/fdu/now/res/'


def load_fits_image(fname):
    hdulist = fits.open(fname)
    hdu_using = hdulist[0]
    header = hdu_using.header

    w = WCS(fname)
    ra_s, dec_s, _, _ = w.all_pix2world(np.arange(header['naxis1']),
                                        np.arange(header['naxis2']),
                                        0, 0, 0)
    d = hdu_using.data.squeeze()

    try:
        bmaj, bmin, bpa = header['bmaj'], header['bmin'], header['bpa']
    except:
        bmaj, bmin, bpa = None, None, None

    hdulist.close()
    return ra_s, dec_s, d, bmaj, bmin, bpa


def draw_beam(ax, lx, ly, pa, loc=None, edgecolor='none', facecolor='white'):
    if loc == None:
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        loc = (xlim[0] + (xlim[1] - xlim[0]) * 0.1,
               ylim[0] + (ylim[1] - ylim[0]) * 0.1)
    ax.add_patch(mpl.patches.Ellipse(loc, lx, ly,
                                     angle=pa, edgecolor=edgecolor, facecolor=facecolor))
    return


def astro_pos_angle_to_image_pos_angle(pos_angle):
    return -(90 + pos_angle)


def astro_pos_angle_to_image_pos_angle_new(pos_angle):
    return pos_angle - 90


def convolve_with_gaussian2d(x, y, d, lx, ly, pa):
    from scipy.ndimage.filters import convolve
    import numpy as np
    
    def gausian2d(lx, ly, t):
        sint = np.sin(t)
        cost = np.cos(t)
        lx, ly = lx/np.sqrt(4*np.log(2)), ly/np.sqrt(4*np.log(2))
        def f(x, y):
            xx = x * cost - y * sint
            yy = x * sint + y * cost
            return np.exp(-((xx/lx)**2+(yy/ly)**2))
        return f
    
    pa_radian = pa * np.pi / 180.0
    lx_coor = lx * np.cos(pa_radian)
    ly_coor = ly * np.sin(pa_radian)
    
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    
    nx = 9 + 2*int(abs(lx_coor / dx))
    ny = 9 + 2*int(abs(ly_coor / dy))
    
    kx, ky = np.meshgrid(np.arange(-nx, nx+1) * dx,
                         np.arange(-ny, ny+1) * dy)
    
    kk = gausian2d(lx, ly, pa_radian)(kx, ky)
    #print kx.shape
    return convolve(d, kk, mode='constant', cval=0.0, origin=0)


def add_colorbar(ax, vmin, vmax, cmap, label='', title='',
                 orientation='horizontal', scale='linear'):
    
    import draw.my_script as my
    from matplotlib.ticker import AutoLocator, LogLocator
    
    ax.set_xscale(scale)
    ax.set_xlim((vmin, vmax))
    
    tickvals = ax.xaxis.get_ticklocs()
    
    ax_inv_transform = ax.transAxes.inverted().transform
    ax_transData = ax.transData.transform
    
    ticks = [ax_inv_transform(ax_transData(np.array([(tic, 0.1)])))[0][0] for tic in tickvals]
    
    #ticklabels = [my.num2str4fig(tv, scientific=False) for tv in tickvals]
    ticklabels = [str(tv) for tv in tickvals]
    
    ticks, ticklabels = zip(*filter(lambda x: 0.0<=x[0]<=1.0, zip(ticks, ticklabels)))
    
    if len(ticks) < 10:
        if scale == 'linear':
            L = AutoLocator()
        else:
            L = LogLocator(base=10, subs=np.linspace(1.0,10.0,10))
        tkvals = L.tick_values(vmin, vmax)
        minticks = filter(lambda x: 0<=x<=1,
                          [ax_inv_transform(ax_transData(np.array([(tic,0.1)])))[0][0] for tic in tkvals])

    ax.set_xscale('linear')

    cbar = mpl.colorbar.ColorbarBase(ax, ticklocation='auto', cmap=cmap, ticks=ticks,
                                     label=label, 
                                     norm=mpl.colors.Normalize(vmin=0, vmax=1, clip=False),
                                     orientation=orientation)

    cbar.solids.set_rasterized(True)
    cbar.set_ticklabels(ticklabels)
    
    if orientation == 'horizontal':
        ax_working = cbar.ax.xaxis
    else:
        ax_working = cbar.ax.yaxis
    
    if len(ticks) < 10:
        ax_working.set_ticks(minticks, minor=True)

    ax_working.set_tick_params(which='minor', length=4)
    ax_working.set_tick_params(which='major', length=7)
    
    ax.set_title(title)

    return


def xy2radec(x, y, ra, dec):
    """
    x, y: relative coordinates in arcsec
    ra, dec: in deg
    """
    ra_cen = (ra[0] + ra[-1]) * 0.5
    dec_cen = (dec[0] + dec[-1]) * 0.5
    return (x/(3600.0*np.cos(dec_cen * np.pi/180.0)) + ra_cen,
            y/3600.0 + dec_cen)


def draw_fits_image(ra, dec, d, fig=None, ax=None, to_arcsec=True,
                    xlim=None, ylim=None, cmap=None, checkDim=False,
                    nan_2_min=True, vmin=None, vmax=None,
                    beam_x=0.37, beam_y=0.56, beam_pa=28.29, drawbeam=True,
                    relative_coordinate=True, ra_cen=None, dec_cen=None):
    if checkDim and (d.shape != (len(dec), len(ra))):
        raise Exception("Dimension does not match!", "d.shape, len(dec), len(ra):", d.shape, len(dec), len(ra))
    if relative_coordinate:
        if ra_cen == None:
            ra_cen = (ra[0] + ra[-1]) * 0.5
        if dec_cen == None:
            dec_cen = (dec[0] + dec[-1]) * 0.5
        x = (ra - ra_cen) * np.cos(dec_cen * np.pi/180.0)
        y = dec - dec_cen
    else:
        x = ra
        y = dec
    if to_arcsec:
        x = x * 3600.0
        y = y * 3600.0
    xlabel = r'$\Delta\theta_{\rm RA}$ (")'
    ylabel = r'$\Delta\theta_{\rm Dec}$ (")'
    if fig == None:
        fig = figure(figsize=(7,7))
    if xlim == None:
        xlim = (x[0], x[-1])
    if ylim == None:
        ylim = (y[0], y[-1])
    if ax == None:
        ax = fig.add_axes((0.1,0.1, 0.89, 0.89),
                          xlim=xlim, ylim=ylim,
                          xlabel=xlabel, ylabel=ylabel,
                          aspect='equal',
                         )
    if cmap == None:
        cmap = mpl.cm.Spectral_r
    if nan_2_min:
        d[np.isnan(d)] = np.nanmin(d)
    ###d = convolve_with_gaussian2d(x, y, d, 2.1, 0.5, 110.0)
    if vmin == None or vmax == None:
        vmin, vmax = np.nanmin(d), np.nanmax(d)
    ax.imshow(d, extent=(x[0], x[-1], y[0], y[-1]), origin='lower',
              vmin=vmin, vmax=vmax,
              cmap=cmap)
    #ax.contourf(x, y, d, 100, cmap=cmap)
    set_axis_format(ax, onlyfont=True)
    if drawbeam:
        draw_beam(ax, beam_x, beam_y, beam_pa)
    return


def my_axes_grid(nrow, ncol, fig=None, figsize=(8,8), xpad=0, ypad=0, aspect=None,
                 xlim=None, ylim=None, xlabel='', ylabel=''):
    if fig == None:
        fig = plt.figure(figsize=figsize)
    xMin, yMin = 0.1, 0.1
    W, H = 0.88, 0.88
    xP, yP = W*xpad, H*ypad
    x_pan = (W - xP*(ncol-1)) / ncol
    y_pan = (H - yP*(nrow-1)) / nrow
    if aspect != None:
        fx, fy = fig.get_size_inches()
        fy_new = x_pan * fx / y_pan * aspect
        fig.set_size_inches((fx, fy_new))
    
    if xlim == None:
        xlim = (0,1)
    if ylim == None:
        ylim = (0,1)
    ax_s = {}
    for i in xrange(nrow):
        ax_s[i] = {}
        for j in xrange(ncol):
            pos_this = (xMin + j * (x_pan + xP),
                        yMin + i * (y_pan + yP),
                        x_pan, y_pan
                       )
            ax = fig.add_axes(pos_this,
                              xlim=xlim, ylim=ylim,
                              xlabel=xlabel, ylabel=ylabel)
            if j > 0:
                ax.set_ylabel('')
                ax.set_yticklabels('')
            if i > 0:
                ax.set_xlabel('')
                ax.set_xticklabels('')
            ax_s[i][j] = ax
    return ax_s


def overlay_contour(ax, x, y, z, levels=None, colors=None, linewidth=None, fmt='%1.2f', fontsize=10, draw_label=True):
    CS = ax.contour(x, y, z, levels, colors=colors, linewidth=linewidth)
    if draw_label:
        plt.clabel(CS, inline=1, colors=colors, fontsize=fontsize, fmt=fmt)



def task1(model_pair, RT_pair):
    model_TWH = os.path.join(res_dir, model_pair['TWH'], RT_pair['TWH'])
    model_DMT = os.path.join(res_dir, model_pair['DMT'], RT_pair['DMT'])
    
    pdfname = '/n/Users/fdu/now/fits_' + \
                model_pair['TWH'] + '_' + RT_pair['TWH'] + '_' + \
                model_pair['DMT'] + '_' + RT_pair['DMT'] +  '.pdf'
    
    fits_data = [
        {'beam': (0.37, 0.56, 28.3), 'label': 'DMT: N=3-2,J=7/2-5/2,F=4-3&3-2',
         'position_angle': 157.0,
         'fits_name': 'dmtau_selfcal_comb_N3_2_J3p5_2p5_ms1.integrated.fits',
         'model_files': [
                os.path.join(model_DMT, 'images/line_00001_00018_00001_2.62004E+11_32.00.fits'),
                os.path.join(model_DMT, 'images/line_00002_00021_00001_2.62004E+11_32.00.fits')]
        },
        {'beam': (0.37, 0.56, 28.3), 'label': 'DMT: N=3-2,J=5/2-3/2,F=3-2&2-1',
         'position_angle': 157.0,
         'fits_name': 'dmtau_selfcal_comb_N3_2_J2p5_1p5_ms1.integrated.fits',
         'model_files': [
                os.path.join(model_DMT, 'images/line_00003_00024_00001_2.62064E+11_32.00.fits'),
                os.path.join(model_DMT, 'images/line_00004_00028_00001_2.62067E+11_32.00.fits')]
        },
        {'beam': (0.41, 0.35, 47.0), 'label': 'TWH: N=4-3,J=7/2-5/2,F=4-3&3-2',
         'position_angle': 150.0,
         'fits_name': 'twhya_selfcal_ext_N4_3_J4p5_3p5_F43_F32.mom8.integrated.fits',
         'model_files': [
                os.path.join(model_TWH, 'images/line_00004_00035_00001_3.49399E+11_7.00.fits'),
                os.path.join(model_TWH, 'images/line_00005_00039_00001_3.49402E+11_7.00.fits')]
        },
        {'beam': (0.41, 0.35, 47.0), 'label': 'TWH: N=4-3,J=9/2-7/2,F=5-4&4-3',
         'position_angle': 150.0,
         'fits_name': 'twhya_selfcal_ext_N4_3_J4p5_3p5_F54_F43.mom8.integrated.fits',
         'model_files': [
                os.path.join(model_TWH, 'images/line_00003_00032_00001_3.49339E+11_7.00.fits'),
                os.path.join(model_TWH, 'images/line_00002_00029_00001_3.49339E+11_7.00.fits')]
        },
    ]
    
    cmap = mpl.cm.Spectral_r
    
    xlim = (5, -5)
    ylim = (-5, 5)
    xlabel = r'$\Delta\theta_{\rm RA}$ (")'
    ylabel = r'$\Delta\theta_{\rm Dec}$ (")'
    
    figsize = (9,18)
    fig = plt.figure(figsize = figsize)
    agrid = my_axes_grid(4, 2, fig=fig,
                         xlim=xlim, ylim=ylim,
                         xlabel=xlabel, ylabel=ylabel,
                         xpad=0.2, ypad=0.01)
    
    for i in xrange(len(fits_data)):
        f_info = fits_data[i]
        ax = agrid[i][0]
    
        fits_name, beam, label = f_info['fits_name'], f_info['beam'], f_info['label']
        #
        ra, dec, d, bmaj, bmin, bpa = load_fits_image(os.path.join(fits_dir, fits_name))
        #
        if bmaj != None:
            beam = (bmaj*3600, bmin*3600, astro_pos_angle_to_image_pos_angle(bpa))
        #
        draw_fits_image(ra, dec, d, fig=fig, ax=ax,
                        beam_x=beam[0], beam_y=beam[1], beam_pa=beam[2])
        vmin, vmax = np.nanmin(d), np.nanmax(d)
    
        xl, yl = ax.get_xlim(), ax.get_ylim()
        ax.text(xl[0] + (xl[1]-xl[0]) * 0.02,
                yl[0] + (yl[1]-yl[0]) * 0.93,
                label, color='white', verticalalignment='bottom')
    
        ax_R = agrid[i][1]
    
        model_files = f_info['model_files']
        
        cb = combine_two_cubes(cube(model_files[0]),
                               cube(model_files[1]))
    
        cb.remove_cube_baseline(n_edgechannel=20)
    
        d = cb.data.sum(axis=0) * abs(cb.v[1]-cb.v[0])
        d = rotate(d, astro_pos_angle_to_image_pos_angle(f_info['position_angle']), reshape=False)
        x_s, y_s = cb.x/cb.dist, cb.y/cb.dist
    
        d = convolve_with_gaussian2d(x_s, y_s, d, beam[0], beam[1], beam[2])
        draw_fits_image(ra, dec, d, fig=fig, ax=ax_R, drawbeam=False, vmin=vmin, vmax=vmax)
    
        xL, yL, wL, hL = agrid[i][0].get_position().bounds
        xR, yR, wR, hR = agrid[i][1].get_position().bounds
    
        x_cen = (xL + wL + xR) / 2
        w_bar = wL * 0.06
        x_sep = xR - xL - wL
        h_cbar = hL * 0.7
        s_cbar = (hL - h_cbar) * 0.5
    
        cbar_pos = (xL + wL + x_sep * 0.2, yL + s_cbar, w_bar, h_cbar)
        ax_cbar = fig.add_axes(cbar_pos)
        add_colorbar(ax_cbar, vmin, vmax, cmap, label='Jy/beam km/s',
                     orientation='vertical', scale='linear')
    
        set_axis_format(agrid[i][1], onlyfont=True)
    
    plt.savefig(pdfname, bbox_inches='tight')


def task2(model_pair, RT_pair, vrange_pair=None):
    model_TWH = os.path.join(res_dir, model_pair['TWH'], RT_pair['TWH'])
    model_DMT = os.path.join(res_dir, model_pair['DMT'], RT_pair['DMT'])
    
    pdfname = '/n/Users/fdu/now/fits_' + \
                model_pair['TWH'] + '_' + RT_pair['TWH'] + '_' + \
                model_pair['DMT'] + '_' + RT_pair['DMT'] +  '.pdf'
    
    fits_data = [
        {'beam': (0.41, 0.35, 47.0), 'label': 'TWH: N=4-3,J=9/2-7/2,F=5-4&4-3',
         'position_angle': 150.0,
         'fits_name': 'twhya_selfcal_ext_N4_3_J4p5_3p5_F54_F43.mom8.integrated.fits',
         'model_files': [
                os.path.join(model_TWH, 'images/line_00003_00032_00001_3.49339E+11_7.00.fits'),
                os.path.join(model_TWH, 'images/line_00002_00029_00001_3.49339E+11_7.00.fits')]
        },
        {'beam': (0.41, 0.35, 47.0), 'label': 'TWH: N=4-3,J=7/2-5/2,F=4-3&3-2',
         'position_angle': 150.0,
         'fits_name': 'twhya_selfcal_ext_N4_3_J4p5_3p5_F43_F32.mom8.integrated.fits',
         'model_files': [
                os.path.join(model_TWH, 'images/line_00004_00035_00001_3.49399E+11_7.00.fits'),
                os.path.join(model_TWH, 'images/line_00005_00039_00001_3.49402E+11_7.00.fits')]
        },
        {'beam': (0.37, 0.56, 28.3), 'label': 'DMT: N=3-2,J=7/2-5/2,F=4-3&3-2',
         'position_angle': 157.0,
         'fits_name': 'dmtau_selfcal_comb_N3_2_J3p5_2p5_ms1.integrated.fits',
         'model_files': [
                os.path.join(model_DMT, 'images/line_00001_00018_00001_2.62004E+11_32.00.fits'),
                os.path.join(model_DMT, 'images/line_00002_00021_00001_2.62004E+11_32.00.fits')]
        },
        {'beam': (0.37, 0.56, 28.3), 'label': 'DMT: N=3-2,J=5/2-3/2,F=3-2&2-1',
         'position_angle': 157.0,
         'fits_name': 'dmtau_selfcal_comb_N3_2_J2p5_1p5_ms1.integrated.fits',
         'model_files': [
                os.path.join(model_DMT, 'images/line_00003_00024_00001_2.62064E+11_32.00.fits'),
                os.path.join(model_DMT, 'images/line_00004_00028_00001_2.62067E+11_32.00.fits')]
        },
    ]
    
    cmap = mpl.cm.Spectral_r
    
    xlim = (5, -5)
    ylim = (-5, 5)
    xlabel = r'$\Delta\theta_{\rm RA}$ (")'
    ylabel = r'$\Delta\theta_{\rm Dec}$ (")'
    
    figsize = (18,9)
    fig = plt.figure(figsize = figsize)
    agrid = my_axes_grid(2, 4, fig=fig,
                         xlim=xlim, ylim=ylim,
                         xlabel=xlabel, ylabel=ylabel,
                         xpad=0.005, ypad=0.09, aspect=1.0)
    
    for i in xrange(len(fits_data)):
        f_info = fits_data[i]
        ax = agrid[1][i]
    
        fits_name, beam, label = f_info['fits_name'], f_info['beam'], f_info['label']
        #
        ra, dec, d, bmaj, bmin, bpa = load_fits_image(os.path.join(fits_dir, fits_name))
        #
        if bmaj != None:
            beam = (bmaj*3600, bmin*3600, astro_pos_angle_to_image_pos_angle(bpa))
        #
        if vrange_pair == None:
            vmin, vmax = np.nanmin(d), np.nanmax(d)
        else:
            if 'TWH' in f_info['label']:
                vmin, vmax = vrange_pair['TWH']
            elif 'DMT' in f_info['label']:
                vmin, vmax = vrange_pair['DMT']
            else:
                raise 'No range provided!'
        #
        draw_fits_image(ra, dec, d, fig=fig, ax=ax,
                        beam_x=beam[0], beam_y=beam[1], beam_pa=beam[2], vmin=vmin, vmax=vmax)
    
        xl, yl = ax.get_xlim(), ax.get_ylim()
        ax.text(xl[0] + (xl[1]-xl[0]) * 0.02,
                yl[0] + (yl[1]-yl[0]) * 0.93,
                label, color='white', verticalalignment='bottom')
    
        ax_R = agrid[0][i]
    
        model_files = f_info['model_files']
        
        cb = combine_two_cubes(cube(model_files[0]),
                               cube(model_files[1]))
    
        cb.remove_cube_baseline(n_edgechannel=20)
    
        d = cb.data.sum(axis=0) * abs(cb.v[1]-cb.v[0])
        d = rotate(d, astro_pos_angle_to_image_pos_angle(f_info['position_angle']), reshape=False)
        x_s, y_s = cb.x/cb.dist, cb.y/cb.dist
    
        d = convolve_with_gaussian2d(x_s, y_s, d, beam[0], beam[1], beam[2])
        draw_fits_image(ra, dec, d, fig=fig, ax=ax_R, drawbeam=False, vmin=vmin, vmax=vmax)
    
        xL, yL, wL, hL = agrid[1][i].get_position().bounds

        cbar_pos = (xL, yL - hL * 0.06, wL, hL * 0.04)
        ax_cbar = fig.add_axes(cbar_pos)
        add_colorbar(ax_cbar, vmin, vmax, cmap, label='Jy/beam km/s',
                     orientation='horizontal', scale='linear')
    
        set_axis_format(agrid[0][i], onlyfont=True)
    
    plt.savefig(pdfname, bbox_inches='tight')



def task3(model_pair, RT_pair, vrange_pair=None):
    model_TWH = os.path.join(res_dir, model_pair['TWH'], RT_pair['TWH'])
    model_DMT = os.path.join(res_dir, model_pair['DMT'], RT_pair['DMT'])
    
    pdfname = '/n/Users/fdu/now/fits_' + \
                model_pair['TWH'] + '_' + RT_pair['TWH'] + '_' + \
                model_pair['DMT'] + '_' + RT_pair['DMT'] +  '.pdf'
    
    fits_data = [
        {'beam': (0.41, 0.35, 47.0), 'label': 'TWH: N=4-3,J=9/2-7/2,F=5-4&4-3',
         'position_angle': 150.0,
         'fits_name': 'twhya_selfcal_ext_N4_3_J4p5_3p5_F54_F43.mom8.integrated.fits',
         'model_files': [
                os.path.join(model_TWH, 'images/line_00003_00032_00001_3.49339E+11_7.00.fits'),
                os.path.join(model_TWH, 'images/line_00002_00029_00001_3.49339E+11_7.00.fits')]
        },
        {'beam': (0.41, 0.35, 47.0), 'label': 'TWH: N=4-3,J=7/2-5/2,F=4-3&3-2',
         'position_angle': 150.0,
         'fits_name': 'twhya_selfcal_ext_N4_3_J4p5_3p5_F43_F32.mom8.integrated.fits',
         'model_files': [
                os.path.join(model_TWH, 'images/line_00004_00035_00001_3.49399E+11_7.00.fits'),
                os.path.join(model_TWH, 'images/line_00005_00039_00001_3.49402E+11_7.00.fits')]
        },
        {'beam': (0.37, 0.56, 28.3), 'label': 'DMT: N=3-2,J=7/2-5/2,F=4-3&3-2',
         'position_angle': 157.0,
         'fits_name': 'dmtau_selfcal_comb_N3_2_J3p5_2p5_ms1.integrated.fits',
         'model_files': [
                os.path.join(model_DMT, 'images/line_00001_00018_00001_2.62004E+11_32.00.fits'),
                os.path.join(model_DMT, 'images/line_00002_00021_00001_2.62004E+11_32.00.fits')]
        },
        {'beam': (0.37, 0.56, 28.3), 'label': 'DMT: N=3-2,J=5/2-3/2,F=3-2&2-1',
         'position_angle': 157.0,
         'fits_name': 'dmtau_selfcal_comb_N3_2_J2p5_1p5_ms1.integrated.fits',
         'model_files': [
                os.path.join(model_DMT, 'images/line_00003_00024_00001_2.62064E+11_32.00.fits'),
                os.path.join(model_DMT, 'images/line_00004_00028_00001_2.62067E+11_32.00.fits')]
        },
    ]
    
    cmap = mpl.cm.Spectral_r
    
    xlim = (4, -4)
    ylim = (-4, 4)
    xlabel = r'$\Delta\theta_{\rm RA}$ (")'
    ylabel = r'$\Delta\theta_{\rm Dec}$ (")'
    
    figsize = (18,9)
    fig = plt.figure(figsize = figsize)
    agrid = my_axes_grid(2, 4, fig=fig,
                         xlim=xlim, ylim=ylim,
                         xlabel=xlabel, ylabel=ylabel,
                         xpad=0.005, ypad=0.09, aspect=1.0)
    
    #CO_data = load_fits_image(os.path.join(fits_dir_CO, 'TWHya_C18O_65_natural_2sig_newvis.mom0.fits'))
    CO_data = load_fits_image(os.path.join(fits_dir_CO, 'TW_Hya_C18O_32_briggs0.5_LSRK_2sig_new.mom0.fits'))

    for i in xrange(len(fits_data)):
        f_info = fits_data[i]
        ax = agrid[1][i]
    
        fits_name, beam, label = f_info['fits_name'], f_info['beam'], f_info['label']
        #
        ra, dec, d, bmaj, bmin, bpa = load_fits_image(os.path.join(fits_dir, fits_name))
        #
        if bmaj != None:
            beam = (bmaj*3600, bmin*3600, astro_pos_angle_to_image_pos_angle_new(bpa))

        print "Data min,max: ", np.nanmin(d), np.nanmax(d)

        if vrange_pair == None:
            vmin, vmax = np.nanmin(d), np.nanmax(d)
        else:
            if 'TWH' in f_info['label']:
                vmin, vmax = vrange_pair['TWH']
            elif 'DMT' in f_info['label']:
                vmin, vmax = vrange_pair['DMT']
            else:
                raise 'No range provided!'
        #
        draw_fits_image(ra, dec, d, fig=fig, ax=ax,
                        beam_x=beam[0], beam_y=beam[1], beam_pa=astro_pos_angle_to_image_pos_angle(beam[2]+90), vmin=vmin, vmax=vmax)
    
        xl, yl = ax.get_xlim(), ax.get_ylim()
        ax.text(xl[0] + (xl[1]-xl[0]) * 0.02,
                yl[0] + (yl[1]-yl[0]) * 0.93,
                label, color='white', verticalalignment='bottom')
    
        ax_R = agrid[0][i]
    
        model_files = f_info['model_files']
        
        cb = combine_two_cubes(cube(model_files[0]),
                               cube(model_files[1]))
    
        cb.remove_cube_baseline(n_edgechannel=40)
    
        d = cb.data.sum(axis=0) * abs(cb.v[1]-cb.v[0])
        d = rotate(d, astro_pos_angle_to_image_pos_angle_new(f_info['position_angle']), reshape=False)
        x_s, y_s = cb.x/cb.dist, cb.y/cb.dist
    
        d = convolve_with_gaussian2d(x_s, y_s, d, beam[0], beam[1], beam[2])
        #draw_fits_image(ra, dec, d, fig=fig, ax=ax_R, drawbeam=False, vmin=vmin, vmax=vmax)
        x_s, y_s = xy2radec(x_s, y_s, ra, dec)
        draw_fits_image(x_s, y_s, d, fig=fig, ax=ax_R, drawbeam=False, vmin=vmin, vmax=vmax)

        print "Model min,max: ", np.nanmin(d), np.nanmax(d)
    
        xL, yL, wL, hL = agrid[1][i].get_position().bounds

        cbar_pos = (xL, yL - hL * 0.06, wL, hL * 0.04)
        ax_cbar = fig.add_axes(cbar_pos)
        add_colorbar(ax_cbar, vmin, vmax, cmap, label='Jy/beam km/s',
                     orientation='horizontal', scale='linear')
    
        set_axis_format(agrid[0][i], onlyfont=True)
    #
    ax = agrid[1][0]
    linewidth = 1
    #
    ra, dec, d, bmaj, bmin, bpa = CO_data

    print "CO Data min,max: ", np.nanmin(d), np.nanmax(d)

    ra_cen  = 0.5 * (ra[0]  + ra[-1])
    dec_cen = 0.5 * (dec[0] + dec[-1])
    x = (ra - ra_cen) * np.cos(dec_cen * np.pi/180.0)
    y = dec - dec_cen
    x = x * 3600.0
    y = y * 3600.0
    xx, yy = np.meshgrid(x, y)
    r = np.sqrt(xx*xx + yy*yy)
    d[r>3] = 0
    vmax, vmin = np.nanmax(d), np.nanmin(d)
    d = d/vmax
    levels = [0.1,0.3,0.5,0.7,0.9]
    #levels = [0.2,0.4,0.6,0.8,0.95]
    overlay_contour(ax, x, y, d, levels=levels, colors='black', linewidth=0.3, draw_label=False)
    #
    plt.savefig(pdfname, bbox_inches='tight')



def task4(model_pair, RT_pair, vrange_pair=None):
    model_TWH = os.path.join(res_dir, model_pair['TWH'], RT_pair['TWH'])
    model_DMT = os.path.join(res_dir, model_pair['DMT'], RT_pair['DMT'])
    
    pdfname = '/n/Users/fdu/now/fits_' + \
                model_pair['TWH'] + '_' + RT_pair['TWH'] + '_' + \
                model_pair['DMT'] + '_' + RT_pair['DMT'] +  '.pdf'
    
    fits_data = [
        {'beam': (0.41, 0.35, 47.0), 'label': 'TWH: N=4-3,J=9/2-7/2,F=5-4&4-3',
         'position_angle': 150.0,
         'fits_name': 'twhya_selfcal_ext_N4_3_J4p5_3p5_F54_F43.mom8.integrated.fits',
         'model_files': [
                os.path.join(model_TWH, 'images/line_00003_00032_00001_3.49339E+11_7.00.fits'),
                os.path.join(model_TWH, 'images/line_00002_00029_00001_3.49339E+11_7.00.fits')]
        },
        {'beam': (0.41, 0.35, 47.0), 'label': 'TWH: N=4-3,J=7/2-5/2,F=4-3&3-2',
         'position_angle': 150.0,
         'fits_name': 'twhya_selfcal_ext_N4_3_J4p5_3p5_F43_F32.mom8.integrated.fits',
         'model_files': [
                os.path.join(model_TWH, 'images/line_00004_00035_00001_3.49399E+11_7.00.fits'),
                os.path.join(model_TWH, 'images/line_00005_00039_00001_3.49402E+11_7.00.fits')]
        },
        {'beam': (0.37, 0.56, 28.3), 'label': 'DMT: N=3-2,J=7/2-5/2,F=4-3&3-2',
         'position_angle': 157.0,
         'fits_name': 'dmtau_selfcal_comb_N3_2_J3p5_2p5_ms1.integrated.fits',
         'model_files': [
                os.path.join(model_DMT, 'images/line_00001_00018_00001_2.62004E+11_32.00.fits'),
                os.path.join(model_DMT, 'images/line_00002_00021_00001_2.62004E+11_32.00.fits')]
        },
        {'beam': (0.37, 0.56, 28.3), 'label': 'DMT: N=3-2,J=5/2-3/2,F=3-2&2-1',
         'position_angle': 157.0,
         'fits_name': 'dmtau_selfcal_comb_N3_2_J2p5_1p5_ms1.integrated.fits',
         'model_files': [
                os.path.join(model_DMT, 'images/line_00003_00024_00001_2.62064E+11_32.00.fits'),
                os.path.join(model_DMT, 'images/line_00004_00028_00001_2.62067E+11_32.00.fits')]
        },
    ]
    
    cmap = mpl.cm.Spectral_r
    
    xlim = (5, -5)
    ylim = (-5, 5)
    xlabel = r'$\Delta\theta_{\rm RA}$ (")'
    ylabel = r'$\Delta\theta_{\rm Dec}$ (")'
    
    figsize = (18,9)
    fig = plt.figure(figsize = figsize)
    agrid = my_axes_grid(2, 4, fig=fig,
                         xlim=xlim, ylim=ylim,
                         xlabel=xlabel, ylabel=ylabel,
                         xpad=0.005, ypad=0.09, aspect=1.0)
    
    #CO_data = load_fits_image(os.path.join(fits_dir_CO, 'TWHya_C18O_65_natural_2sig_newvis.mom0.fits'))
    CO_data = load_fits_image(os.path.join(fits_dir_CO, 'TW_Hya_C18O_32_briggs0.5_LSRK_2sig_new.mom0.fits'))
    #
    CO_model_file = os.path.join(res_dir, '20160714_TWH_flatten/C18O/images/line_00003_00003_00001_3.45798E+11_7.00.fits')

    for i in xrange(len(fits_data)):
        f_info = fits_data[i]
        ax = agrid[1][i]
    
        fits_name, beam, label = f_info['fits_name'], f_info['beam'], f_info['label']
        #
        ra, dec, d, bmaj, bmin, bpa = load_fits_image(os.path.join(fits_dir, fits_name))
        #
        if bmaj != None:
            beam = (bmaj*3600, bmin*3600, astro_pos_angle_to_image_pos_angle(bpa))

        print "Data min,max: ", np.nanmin(d), np.nanmax(d)

        if vrange_pair == None:
            vmin, vmax = np.nanmin(d), np.nanmax(d)
        else:
            if 'TWH' in f_info['label']:
                vmin, vmax = vrange_pair['TWH']
            elif 'DMT' in f_info['label']:
                vmin, vmax = vrange_pair['DMT']
            else:
                raise 'No range provided!'
        #
        draw_fits_image(ra, dec, d, fig=fig, ax=ax,
                        beam_x=beam[0], beam_y=beam[1], beam_pa=beam[2], vmin=vmin, vmax=vmax)
    
        xl, yl = ax.get_xlim(), ax.get_ylim()
        ax.text(xl[0] + (xl[1]-xl[0]) * 0.02,
                yl[0] + (yl[1]-yl[0]) * 0.93,
                label, color='white', verticalalignment='bottom')
    
        ax_R = agrid[0][i]
    
        model_files = f_info['model_files']
        
        cb = combine_two_cubes(cube(model_files[0]),
                               cube(model_files[1]))
    
        cb.remove_cube_baseline(n_edgechannel=40)
    
        d = cb.data.sum(axis=0) * abs(cb.v[1]-cb.v[0])
        d = rotate(d, astro_pos_angle_to_image_pos_angle(f_info['position_angle']), reshape=False)
        x_s, y_s = cb.x/cb.dist, cb.y/cb.dist
    
        d = convolve_with_gaussian2d(x_s, y_s, d, beam[0], beam[1], beam[2])
        draw_fits_image(ra, dec, d, fig=fig, ax=ax_R, drawbeam=False, vmin=vmin, vmax=vmax)
        #x_s, y_s = xy2radec(x_s, y_s, ra, dec)
        #draw_fits_image(x_s, y_s, d, fig=fig, ax=ax_R, drawbeam=False, vmin=vmin, vmax=vmax)

        print "Model min,max: ", np.nanmin(d), np.nanmax(d)
    
        xL, yL, wL, hL = agrid[1][i].get_position().bounds

        cbar_pos = (xL, yL - hL * 0.06, wL, hL * 0.04)
        ax_cbar = fig.add_axes(cbar_pos)
        add_colorbar(ax_cbar, vmin, vmax, cmap, label='Jy/beam km/s',
                     orientation='horizontal', scale='linear')
    
        set_axis_format(agrid[0][i], onlyfont=True)
    #
    ax = agrid[1][0]
    linewidth = 1
    #
    ra, dec, d, bmaj, bmin, bpa = CO_data

    print "CO Data min,max: ", np.nanmin(d), np.nanmax(d)

    ra_cen  = 0.5 * (ra[0]  + ra[-1])
    dec_cen = 0.5 * (dec[0] + dec[-1])
    x = (ra - ra_cen) * np.cos(dec_cen * np.pi/180.0)
    y = dec - dec_cen
    x = x * 3600.0
    y = y * 3600.0
    xx, yy = np.meshgrid(x, y)
    r = np.sqrt(xx*xx + yy*yy)
    d[r>3] = 0
    vmax, vmin = np.nanmax(d), np.nanmin(d)
    d = d/vmax
    levels = [0.1,0.3,0.5,0.7,0.9]
    overlay_contour(ax, x, y, d, levels=levels, colors='black', linewidth=0.3, draw_label=False)
    #
    def overlay_contour_cb(ax, cb, ra=None, dec=None, levels=None):
        cb.remove_cube_baseline(n_edgechannel=40)
        d = cb.data.sum(axis=0) * abs(cb.v[0] - cb.v[1])
        d = rotate(d, astro_pos_angle_to_image_pos_angle(f_info['position_angle']), reshape=False)
        x_s, y_s = cb.x/cb.dist, cb.y/cb.dist
        d = convolve_with_gaussian2d(x_s, y_s, d, 0.4, 0.4, 0.0)
        #x_s, y_s = xy2radec(x_s, y_s, ra, dec)
        vmax, vmin = np.nanmax(d), np.nanmin(d)
        d = d/vmax
        overlay_contour(ax, x_s, y_s, d, levels=levels, colors='black', linewidth=0.3, draw_label=False)
    #
    ax = agrid[0][0]
    cb_CO = cube(CO_model_file)
    overlay_contour_cb(ax, cb_CO, ra=ra, dec=dec, levels=levels)
    #
    plt.savefig(pdfname, bbox_inches='tight')



if __name__ == '__main__':
    
    model_pair_s = [\
                    #{'TWH': '20150609_C2H_5a_4a_2d_8cl_B0',
                    # 'DMT': '20150609_DMT_C2H_5x_4_5a_3a_rerun24'},
                    {'TWH': '20160714_TWH_d_s_noDep_moreOuterMass1',
                     'DMT': '20160830_DMT_d'},
                    #{'TWH': '20160714_TWH_i',
                    # 'DMT': '20150609_DMT_C2H_5x_4_5a_3a_rerun24'},
                    #{'TWH': '20160714_TWH_j',
                    # 'DMT': '20150609_DMT_C2H_5x_4_5a_3a_rerun24'},
                   ]
    
    RT_pair_s = [\
                    #{'TWH': 'C2H_2_NLTE_redu1',
                    # 'DMT': 'C2H_1_NLTE'},
                    {'TWH': 'C2H_2_NLTE',
                     'DMT': 'C2H_1_NLTE'},
                ]*5
    
    vranges = [\
                    #None,
                    {'TWH': (0,0.29),
                     'DMT': (0,0.11)},
              ]*5
    
    
    for model_pair, RT_pair, v_pair in zip(model_pair_s, RT_pair_s, vranges):
        task3(model_pair, RT_pair, vrange_pair = v_pair)
