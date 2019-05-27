import numpy as np
import matplotlib as mpl

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


def add_colorbar(ax, vmin, vmax, cmap, label='', title='',
                 orientation='horizontal', scale='linear'):
    
    #import draw.my_script as my
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
        ax_working.set_ticks([_ for _ in minticks], minor=True)

    ax_working.set_tick_params(which='minor', length=4)
    ax_working.set_tick_params(which='major', length=7)
    
    ax.set_title(title)

    return


def draw_fits_image(ra, dec, d, fig=None, ax=None, to_arcsec=True,
                    xlim=None, ylim=None, cmap=None,
                    nan_2_min=True, vmin=None, vmax=None,
                    beam_x=0.37, beam_y=0.56, beam_pa=28.29, drawbeam=True,
                    relative_coordinate=True, ra_cen=None, dec_cen=None):
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
    for i in range(nrow):
        ax_s[i] = {}
        for j in range(ncol):
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
