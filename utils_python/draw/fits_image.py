import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation

from os.path import join as opj

from astropy.io import fits
import numpy as np
import consts.physical_constants as consts

def least_square_fit_linear(x, y):
    """
    y = a * x + b
    """
    n = x.size
    sx = x.sum()
    sxx = (x*x).sum()
    sy = y.sum()
    syy = (y*y).sum()
    sxy = (x*y).sum()
    a = (n * sxy - sx * sy) / (n*sxx - sx*sx)
    b = (sxx * sy - sx * sxy) / (n*sxx - sx*sx)
    return (a, b)


def least_square_fit_cube(cb, x = None, deg=1):
    nx, ny, nz = cb.shape
    if x is None:
        x = np.arange(0, nx)
    a = [np.empty((ny, nz)) for i in xrange(deg+1)]
    for i in xrange(ny):
        for j in xrange(nz):
            #a[i,j], b[i,j] = least_square_fit_linear(x, cb[:, i, j])
            res = np.polyfit(x, cb[:, i, j], deg)
            for k in xrange(deg+1):
                a[k][i, j] = res[k]
    return a



def number_rounding(val):
    sign = 1
    if val < 0:
        sign = -1
        val = -val
    expo = np.floor(np.log10(val))
    ex = 10.0**expo
    coef = np.ceil(val / ex*10)/10
    #print locals()
    return sign * coef * ex



class cube():
    def __init__(self, fname, saveTauMap=False):
        self.h = fits.open(fname)
        ext = self.h[0]
        #
        hd = ext.header
        if hd['EXTNAME'] != 'LineCube':
            raise 'Not a line cube!'
        #
        self.fname = fname
        self.f0 = hd['F0']
        self.E_up = hd['EUP']
        self.E_low = hd['ELOW']
        self.Aul = hd['AUL']
        self.Bul = hd['BUL']
        self.Blu = hd['BLU']
        self.intfluxl = hd['INTFLUXL']
        self.intflux = hd['INTFLUX']
        self.maxflux = hd['MAXFLUX']
        self.maxtau = hd['MAXTAU']
        try:
            self.molname = hd['MOL-DB']
            self.molname_disp = hd['MOL-DSP']
            if self.molname_disp.strip() == '':
                self.molname_disp = self.molname
        except:
            self.molname = hd['MOL']
            self.molname_disp = hd['MOL']
        self.dist = hd['DIST']
        self.theta = hd['THETA']
        self.lam0 = hd['LAM0']
        if 'QNUM' in hd.keys():
            self.qnum = hd['QNUM']
        else:
            self.qnum = ''
        #
        self.nx  = hd['NAXIS1']
        self.ny  = hd['NAXIS2']
        self.nf  = hd['NAXIS3']
        self.dx  = hd['CDELT1']
        self.dy  = hd['CDELT2']
        self.df  = hd['CDELT3']
        self.k0x = hd['CRPIX1']
        self.k0y = hd['CRPIX2']
        self.k0f = hd['CRPIX3']
        self.xr  = hd['CRVAL1']
        self.yr  = hd['CRVAL2']
        self.fr  = hd['CRVAL3']
        #
        self.x = (np.arange(1, self.nx+1) - self.k0x) * self.dx + self.xr
        self.y = (np.arange(1, self.ny+1) - self.k0y) * self.dy + self.yr
        self.f = (np.arange(1, self.nf+1) - self.k0f) * self.df + self.fr
        #
        self.v = (self.f0 - self.f) / self.f0 * consts.phy_SpeedOfLight_SI / 1e3

        tmp_coeff = self.dx * self.dy * (consts.phy_AU2cm / (self.dist*consts.phy_pc2cm))**2 / \
                    consts.phy_jansky2CGS
        #
        # Jansky per pixel
        self.data = ext.data * tmp_coeff
        #
        if saveTauMap:
            self.tauMapData = self.h[1].data

    def remove_cube_baseline(self, n_edgechannel=2):
        """
        self.data: nx * ny * nz
        """
        nx, ny, nz = self.data.shape
        idx = np.hstack((np.arange(0,n_edgechannel), np.arange(-n_edgechannel, 0)))
        x = np.hstack((np.arange(0, n_edgechannel), np.arange(nx-n_edgechannel, nx)))
        a = least_square_fit_cube(self.data[idx, :, :], x = x, deg = 2)
        a = a[::-1]
        xnew = np.arange(0, nx)
        for i in xrange(ny):
            for j in xrange(nz):
                self.data[:, i, j] = self.data[:, i, j] - \
                    np.polynomial.polynomial.polyval(xnew, np.array([co[i, j] for co in a]))

    def trim_tiny_value(self, threshold_relative = 1e-20):
        smallest = np.median(self.data) * threshold_relative
        self.data[self.data < smallest] = smallest

    def convol(self, FWHM=1.0, sumbeam=True, tau_also=False):
        from scipy.ndimage.filters import gaussian_filter as gsfl
        from scipy.ndimage.filters import generic_filter as grfl
        pixel_arcsec = self.dx/self.dist
        nx,_,_ = self.data.shape
        sigma_in_pixel = np.ceil(FWHM/np.sqrt(8*np.log(2))/pixel_arcsec)
        self.FWHM = FWHM
        if np.isscalar(FWHM) or len(FWHM) == 1:
            for i in xrange(nx):
                self.data[i, :, :] = gsfl(self.data[i, :, :],
                                          sigma_in_pixel,
                                          mode = 'constant', cval = 0.0)
        #else:
        #    def func_gauss_rot():
        #    for i in xrange(nx):
        #        self.data[i, :, :] = grfl(self.data[i, :, :],
        #                                  func_gauss_rot,
        #                                  mode = 'constant', cval = 0.0)
        if tau_also and hasattr(self, 'tauMapData'):
            self.tauMapData = gsfl(self.tauMapData,
                                   sigma_in_pixel,
                                   mode = 'constant', cval = 0.0)
        if sumbeam:
            # Jansky per beam
            self.data = self.data *  (2.0 * np.pi * sigma_in_pixel**2)


    def draw_a_channel(self, ich, pdfname, axisUnit='"',
                       xRange=None, yRange=None):
        #draw_an_image(np.log10(self.data[ich, :, :]), self.x, self.y, pdfname,
        draw_an_image(self.data[ich, :, :],
                      self.x / (1.0 if axisUnit == 'AU' else self.dist),
                      self.y / (1.0 if axisUnit == 'AU' else self.dist),
                      pdfname,
                      axisUnit=axisUnit,
                      xRange=xRange, yRange=yRange,
                      tag = '$\lambda$ = {0:.2e} $\mu$m, peak = {1:.2e} Jy/beam'.\
                        format(self.lam0/1e4, np.nanmax(self.data[ich, :, :])))

    def draw_peak_channel(self, pdfname):
        ich,ix,iy = np.unravel_index(self.data.argmax(), self.data.shape)
        fig = draw_an_image(self.data[ich, :, :], self.x, self.y, pdfname,
                      tag = '{0:s}: peak = {1:.2e} Jy/beam at V = {2:.2f} km/s'.\
                        format(self.molname_disp, np.nanmax(self.data[ich, :, :]),
                               self.v[ich]),
                      returnfig = True,
                      savepdf=False)
        draw_a_spectrum(self.v, self.data[:,ix,iy], fig, xRange=(-10,10))
        plt.savefig(pdfname, bbox_inches='tight')


    def fit_gauss_cube(self):
        from misc.gauss import fit_gauss, get_gauss_area
        nx, ny, nv = self.nx, self.ny, self.nf
        d = np.zeros((nx,ny))
        for i in xrange(nx):
            for j in xrange(ny):
                print self.v[0], self.v[-1]
                d[i,j], _ = get_gauss_area(*fit_gauss(self.v, self.data[:, i,j]))
                print i, j, d[i,j]
        return d


    def draw_sum_channels(self, pdfname, returnfig=False, returnax=False, savepdf=True,
                logscale=False, notext=False,
                drawPeak=False,
                pos = (0.1,0.1,0.85,0.85), fig=None, ax=None, vmin=None, vmax=None,
                noColorbar=False, ytextoff=False,
                cax=None, orientation='horizontal', cbar_ticks=None, pos_angle=None,
                axisUnit='AU', xRange=None, yRange=None, rescaleTo=None, cbarcolor='white',
                gaussianfit=False,
            ):
        if logscale:
            d = np.log10(self.data.sum(axis=0) * abs(self.v[1]-self.v[0]))
            tagstr = '{0:s}\tUnit: log10(Jy km s$^{{-1}}$/beam)\nTotal: {1:.2e} (W m$^{{-2}}$)'.\
                            format(self.molname_disp, self.intfluxl)
            if rescaleTo != None:
                _vmax = np.nanmax(d)
                d = d * (rescaleTo / _vmax)
            return draw_an_image(d,
                      self.x / (1.0 if axisUnit == 'AU' else self.dist),
                      self.y / (1.0 if axisUnit == 'AU' else self.dist),
                      pdfname,
                      axisUnit = axisUnit,
                      tag = tagstr if (not notext) else '',
                      vmax = vmax, vmin = vmin,
                      returnfig = returnfig, returnax=returnax,
                      savepdf = savepdf, xRange=xRange, yRange=yRange,
                      pos=pos, fig=fig, ax=ax, noColorbar=noColorbar, ytextoff=ytextoff, cbarcolor=cbarcolor,
                      cax=cax, orientation=orientation, cbar_ticks=cbar_ticks, pos_angle=pos_angle)
        elif drawPeak:
            #d = self.data.max(axis=0)  # The unit is Jy
            d = np.nanmax(self.data, axis=0)  # The unit is Jy
            #d = self.data[self.nf/2,:,:]
            jy_per_beam_to_K = 1e-26/(self.FWHM*np.pi/3600/180)**2 * (self.lam0*1e-10)**2 / (2*1.38e-23)
            tagstr = '{0:s}\nF: {1:.6e} Hz  $A_{{\\rm ul}}$: {2:.1e} s$^{{-1}}$\n$E_{{\\rm up}}$: {3:.1f} K  $E_{{\\rm low}}$: {4:.1f} K\nTotal: {5:.2e} (W m$^{{-2}}$)\nBeam: {6:.2f}$^{{\prime\prime}}$\nUnit: Jy/beam\nJy/beam to K: {7:.1e}'.\
                            format(self.molname_disp + '  ' + format_qnum(self.qnum),
                            self.f0, self.Aul, self.E_up, self.E_low,
                            self.intfluxl, self.FWHM, jy_per_beam_to_K)
            return draw_an_image(d,
                      self.x / (1.0 if axisUnit == 'AU' else self.dist),
                      self.y / (1.0 if axisUnit == 'AU' else self.dist),
                      pdfname,
                      axisUnit = axisUnit,
                      tag = tagstr if (not notext) else '',
                      returnfig = returnfig, returnax=returnax,
                      savepdf = savepdf, xRange=xRange, yRange=yRange,
                      pos=pos, fig=fig, ax=ax, vmin=vmin, vmax=vmax, noColorbar=noColorbar, ytextoff=ytextoff, cbarcolor=cbarcolor,
                      cax=cax, orientation=orientation, cbar_ticks=cbar_ticks, pos_angle=pos_angle)
        else:
            if gaussianfit:
                d = self.fit_gauss_cube()
            else:
                d = self.data.sum(axis=0) * abs(self.v[1]-self.v[0])
            print self.intfluxl
            if rescaleTo != None:
                _vmax = np.nanmax(d)
                d = d * (rescaleTo / _vmax)
            jy_km_s_per_beam_to_K_km_s = 1e-26/(self.FWHM*np.pi/3600/180)**2 * (self.lam0*1e-10)**2 / (2*1.38e-23)
            tagstr = '{0:s}\nF: {1:.6e} Hz  $A_{{\\rm ul}}$: {2:.1e} s$^{{-1}}$\n$E_{{\\rm up}}$: {3:.1f} K  $E_{{\\rm low}}$: {4:.1f} K\nTotal: {5:.2e} (W m$^{{-2}}$)\nBeam: {6:.2f}$^{{\prime\prime}}$\nUnit: Jy km s$^{{-1}}$/beam\nJy/beam to K: {7:.1e}'.\
                            format(self.molname_disp + '  ' + format_qnum(self.qnum),
                            self.f0, self.Aul, self.E_up, self.E_low,
                            self.intfluxl, self.FWHM, jy_km_s_per_beam_to_K_km_s)
            return draw_an_image(d,
                      self.x / (1.0 if axisUnit == 'AU' else self.dist),
                      self.y / (1.0 if axisUnit == 'AU' else self.dist),
                      pdfname,
                      axisUnit = axisUnit,
                      tag = tagstr if (not notext) else '',
                      returnfig = returnfig, returnax=returnax,
                      savepdf = savepdf, xRange=xRange, yRange=yRange,
                      pos=pos, fig=fig, ax=ax, vmin=vmin, vmax=vmax, noColorbar=noColorbar, ytextoff=ytextoff, cbarcolor=cbarcolor,
                      cax=cax, orientation=orientation, cbar_ticks=cbar_ticks, pos_angle=pos_angle)

    def toArcsec(self):
        self.x = self.x / self.dist
        self.y = self.y / self.dist
        self.dx = self.dx / self.dist
        self.dy = self.dy / self.dist


def draw_an_image(data, x, y, pdfname, tag='', tagcolor='white',
                  axisUnit='AU', xRange=None, yRange=None,
                  figsize=(10,8), pos=(0.1,0.1,0.85,0.85),
                  noColorbar=False, fontsize=10,
                  vmin=None, vmax=None,
                  savepdf=True,
                  returnfig = False, returnax=False,
                  cmap = None, badcolor=None, cax=None, orientation='horizontal', cbar_ticks=None, cbarcolor='white',
                  fig=None, ax=None, ytextoff=False,
                  pos_angle=None
                  ):
    from draw.long_function_definitions import set_axis_format
    from scipy.ndimage.interpolation import rotate
    if fig == None:
        fig = plt.figure(figsize = figsize)
    if ax == None:
        ax = fig.add_axes(pos)
    if vmin == None:
        vmin = np.nanmin(data)
    if vmax == None:
        vmax = np.nanmax(data)
    if cmap == None:
        cmap = mpl.cm.Spectral_r
        #cmap = mpl.cm.Greys_r
    if badcolor==None:
        badcolor = cmap(0)
    cmap.set_bad(badcolor, 1.)
    if pos_angle != None:
        data = rotate(data, pos_angle, reshape=False)
    C = ax.imshow(data,
                  interpolation='none', origin='lower',
                  vmin=vmin, vmax=vmax,
                  extent=(x.min(), x.max(),
                          y.min(), y.max()),
                  cmap=cmap)
    if tag != '':
        ax.text(0.05, 0.98, tag, color=tagcolor,
                transform=ax.transAxes,
                fontsize=fontsize,
                fontdict={'horizontalalignment': 'left',
                          'verticalalignment': 'top'})
    ax.set_xlabel('X (' + axisUnit + ')')
    if not ytextoff:
        ax.set_ylabel('Y (' + axisUnit + ')')
    if xRange != None:
        ax.set_xlim(xRange)
    if yRange != None:
        ax.set_ylim(yRange)
    #
    if (not noColorbar) and cax:
        cb = plt.colorbar(C, cax=cax, orientation=orientation, ticks=cbar_ticks)
        set_axis_color(cax, color=cbarcolor)
        cb.outline.set_edgecolor(cbarcolor)
    ax.set_aspect('equal')
    set_axis_format(ax, xscale='linear', yscale='linear', labelfontsize=fontsize, tickfontsize=fontsize*0.8,
                    majorgridon=False, minorgridon=False, graygrid=True)
    if ytextoff:
        ax.yaxis.set_ticklabels([])
    if savepdf:
        plt.savefig(pdfname, bbox_inches='tight')
        print 'Pdf saved:', pdfname
    if returnax:
        return fig, ax
    if returnfig:
        return fig


def set_axis_color(ax, color='black'):
    import matplotlib
    ax.title.set_color(color)
    ax.tick_params(labelcolor=color, color=color)
    plt.setp(ax.spines.values(), color=color)
    plt.setp([ax.get_xticklines(), ax.get_yticklines()], color=color)
    ax.spines['bottom'].set_color(color)
    ax.spines['top'].set_color(color) 
    ax.spines['right'].set_color(color)
    ax.spines['left'].set_color(color)
    for child in ax.get_children():
        if isinstance(child, matplotlib.spines.Spine):
                child.set_color(color)


def format_qnum(qnum_str):
    s = qnum_str.split()
    idx = s.index('->')
    return s[idx-1] + ' -> ' + s[-1]


def draw_a_spectrum(v, y, fig, ax=None, xRange=None, yRange=None,
                    xlab='$V$ (km s$^{{-1}}$)', ylab='Jy/beam',
                    linestyle='-', color='white', pos=(0.25,0.2,0.2,0.2),
                    returnax=False):
    from draw.long_function_definitions import set_axis_format
    if ax == None:
        ax = fig.add_axes(pos)
        ax.patch.set_alpha(0)
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
        if xRange != None:
            ax.set_xlim(xRange)
        if yRange != None:
            ax.set_ylim(yRange)
        set_axis_format(ax, xscale='linear', yscale='linear',
                        labelfontsize=15, tickfontsize=10,
                        majorgridon=False, minorgridon=False)
    ax.plot(v, y, lw=1, color=color, linestyle=linestyle)
    if returnax:
        return ax


def make_RGB(d1, d2, d3, x, y, fname, tag='', tagcolor='none',
             axisUnit='AU', xRange=None, yRange=None,
             thrsh=None, scale='linear',
             contourLines=None,
             contourLevels=None,
             ):
    if contourLines != None:
        _d1 = d1.copy()
        _d2 = d2.copy()
        _d3 = d3.copy()
    if thrsh != None:
        ds = [d1,d2,d3]
        for i in xrange(3):
            d = ds[i]
            mx = np.nanmax(d) * thrsh[i]
            d[d >= mx] = mx
    if scale == 'linear':
        draw_an_image( \
            np.dstack((d1/np.nanmax(d1),
                       d2/np.nanmax(d2),
                       d3/np.nanmax(d3))),
            x, y, fname, tag=tag, tagcolor=tagcolor,
            axisUnit=axisUnit, xRange=xRange, yRange=xRange,
            noColorbar=True, savepdf=False)
    elif scale == 'log':
        for d in [d1,d2,d3]:
            d = np.log10(d)
        d1mx, d1mn = np.nanmax(d1), np.nanmin(d1)
        d2mx, d2mn = np.nanmax(d2), np.nanmin(d2)
        d3mx, d3mn = np.nanmax(d3), np.nanmin(d3)
        draw_an_image( \
            np.dstack(((d1-d1mn)/(d1mx-d1mn),
                       (d2-d2mn)/(d2mx-d2mn),
                       (d3-d3mn)/(d3mx-d3mn))),
            x, y, fname, tag=tag, tagcolor=tagcolor,
            axisUnit=axisUnit, xRange=xRange, yRange=xRange,
            noColorbar=True, savepdf=False)
    else:
        return
    if contourLines != None:
        ax = plt.gca()
        ds = [_d1, _d2, _d3]
        for i in xrange(3):
            if contourLines[i] != '':
                ax.contour(x, y, ds[i],
                           levels=[np.nanmax(ds[i])*contourLevels[i]],
                           linestyles=contourLines[i],
                           linewidths=2,
                           colors='white',
                           )
    plt.savefig(fname, bbox_inches='tight')


def make_a_video(cb, fname, ifrRange=None,
                 figsize=(8,8), pos=(0.15,0.1,0.75,0.75),
                 shapeType=1,
                 labelfmtstr='{0:.2f} km s$^{{-1}}$'):
    from draw.long_function_definitions import set_axis_format
    #
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=24, metadata=dict(artist='Me'), bitrate=1000)
    #
    cbar_y = pos[1]*1.5 + pos[3]
    pos_colorbar = (pos[0], cbar_y, pos[2], min(0.1, max(0.95 - cbar_y, 0.02)))
    fig = plt.figure(figsize = figsize)
    ax = fig.add_axes(pos, xlabel='X (AU)', ylabel='Y (AU)')
    ax_colorbar = fig.add_axes(pos_colorbar)
    ax_colorbar.set_xlabel('')
    #
    mn, mx = np.nanmin(cb.data), np.nanmax(cb.data)
    mn = max(mx*1e-4, mn)
    #
    x = cb.x
    y = cb.y
    nx = cb.data.shape[shapeType-1]
    ims = []
    if ifrRange == None:
        imin, imax = 0, nx
    else:
        imin, imax = ifrRange[0], ifrRange[1]
    for i in xrange(imin, imax):
        print 'Frame', i-imin
        vel = cb.v[i]
        D = selectSliceIofAxisJ(cb.data, i, shapeType-1)
        C = ax.imshow(D,
                      interpolation='none', origin='lower',
                      vmin=mn, vmax=mx,
                      extent=(x.min(), x.max(),
                              y.min(), y.max()),
                      cmap=mpl.cm.Spectral_r)
        T = ax.text(0.1, 0.87, labelfmtstr.format(vel),
                    transform=ax.transAxes,
                    fontdict={'horizontalalignment': 'left',
                              'verticalalignment': 'bottom'})
        if i == imin:
            fig.colorbar(C, cax=ax_colorbar, orientation='horizontal')
        #
        ims.append((C, T))
        ax.set_aspect('equal')
        set_axis_format(ax, xscale='linear', yscale='linear',
                        majorgridon=False, minorgridon=False)
    #
    print 'Converting video...'
    im_ani = animation.ArtistAnimation(fig, ims, interval=40, 
                repeat=True, repeat_delay=2000, blit=True)
    print 'Saving video...'
    im_ani.save(fname, writer=writer)


def selectSliceIofAxisJ(a, i, j):
    n = a.shape[j]
    return np.squeeze(np.compress(np.arange(n) == i, a, axis=j), axis=j)


def reverse_idx(v, val):
    '''
    v must be monotonic
    '''
    from scipy import interpolate
    x = np.arange(0, len(v))
    f = interpolate.interp1d(v, x)
    return f(val)


def eval_at_fractional_idx(v, idx):
    from scipy import interpolate
    x = np.arange(0, len(v))
    f = interpolate.interp1d(x, v)
    return f(idx)


def combine_two_cubes(cb1, cb2):
   cb1.data = cb1.data + cb2.data
   return cb1

def task1():
    FWHM = 0.1
    #
    model_dir = '20150418_btfl'
    res_dir = '/u/Moria2/fdu/now/res/'
    #
    SO_dir = 'SO'
    CS_dir = 'CS'
    #
    linefname_SO = 'line_00002_00042_00001_2.15221E+11_89.50.fits'
    linefname_CS = 'line_00001_00005_00001_2.44936E+11_89.50.fits'
    #
    print 'Loading data'
    fname_SO = opj(res_dir, model_dir, SO_dir, 'images', linefname_SO)
    cb = cube(fname_SO)
    #
    print 'Convolving'
    cb.convol(FWHM=FWHM)
    #
    print 'Saving a continuum channel'
    cb.draw_a_channel(10, fname_SO + '_' + model_dir + '_cont.pdf')
    d_cont = cb.data[10,:,:].copy()
    #
    print 'Removing baseline'
    cb.remove_cube_baseline(n_edgechannel=20)
    #
    print 'Saving the peak line channel'
    cb.draw_peak_channel(fname_SO + '_' + model_dir + '_peak.pdf')
    #
    #cb.trim_tiny_value()
    print 'Drawing total flux'
    cb.draw_sum_channels(fname_SO + '_' + model_dir + '_sum.pdf')
    d_SO_sum = cb.data.sum(axis=0).copy()
    #
    print 'Making video'
    #make_a_video(cb, fname_SO + '_' + model_dir + '.mp4', ifrRange=(50,150))

    print 'Loading data'
    fname_CS = opj(res_dir, model_dir, CS_dir, 'images', linefname_CS)
    cb1 = cube(fname_CS)
    #
    print 'Convolving'
    cb1.convol(FWHM=FWHM)
    #
    print 'Saving a continuum channel'
    cb1.draw_a_channel(10, fname_CS + '_' + model_dir + '_cont.pdf')
    #
    print 'Removing baseline'
    cb1.remove_cube_baseline(n_edgechannel=20)
    #
    print 'Saving the peak line channel'
    cb1.draw_peak_channel(fname_CS + '_' + model_dir + '_peak.pdf')
    #
    #cb1.trim_tiny_value()
    print 'Drawing total flux'
    cb1.draw_sum_channels(fname_CS + '_' + model_dir + '_sum.pdf')
    d_CS_sum = cb1.data.sum(axis=0).copy()
    #
    print 'Draing RGB map'
    make_RGB(d_SO_sum, d_CS_sum, d_cont, cb.x/cb.dist, cb.y/cb.dist,
             opj(res_dir, model_dir, 'SO_CS_cont' + '_' + model_dir + '.pdf'),
             tag = 'Red: SO 215 GHz\nGreen: CS 245 GHz\nBlue: 1.4 mm continuum',
             tagcolor = 'white',
             axisUnit='"', xRange=(-1,1), yRange=(-1,1),
             thrsh = [0.5, 1.0, 1.0], scale='linear',
             contourLevels=[0.7,0.7,0.7], contourLines=['solid', 'dashed', 'dotted'],
            )
    #
    print 'Making video'
    #make_a_video(cb1, fname_CS + '_' + model_dir + '.mp4', ifrRange=(50,150))

    print 'Drawing ratio of sum'
    cbsum = cb.data.sum(axis=0)
    ratio = np.log10(cbsum / cb1.data.sum(axis=0))
    ratio[cbsum < cbsum.max()*1e-1] = np.nan
    draw_an_image(ratio, cb.x, cb.y,
                  opj(res_dir, model_dir, 'SO_dvd_CS_sum' + '_' + model_dir + '.pdf'),
                  tag = 'log10(SO/CS)')

    #print 'Making video for channel ratios'
    #nx,_,_ = cb.data.shape
    #for i in xrange(nx):
    #    maxval = np.nanmax(cb.data[i,:,:])
    #    cb.data[i,:,:][cb.data[i,:,:] < maxval*0.1] = np.nan
    #    cb.data[i,:,:] = cb.data[i,:,:] / cb1.data[i,:,:]
    #    cb.data[i,:,:] = cb.data[i,:,:] / np.nanmax(cb.data[i,:,:])
    #make_a_video(cb, opj(res_dir, model_dir, 'SO_dvd_CS_channel.mp4'),
    #             ifrRange=(50,150))


def task2():
    FWHM = 0.6
    #
    print 'Loading data'
    cb = cube('/u/Moria2/fdu/now/res/20150415_btfl/SO/images/line_00001_00038_00001_2.06176E+11_89.50.fits')
    #
    print 'Convolving'
    cb.convol(FWHM=FWHM)
    #
    print 'Saving a continuum channel'
    cb.draw_a_channel(10, '/n/Users/fdu/now/figures/SO_cont_206GHz.pdf')


def task3():
    import glob
    import os
    FWHM = 0.1
    #
    image_dir = '/u/Moria2/fdu/now/res/20150418_btfl/CS_angles*/images/'
    linefnames = glob.glob(os.path.join(image_dir, '*.fits'))
    linefnames.sort()
    #
    slices = []
    thetas = []
    for fname in linefnames:
        print 'Loading ', fname
        cb = cube(fname)
        cb.convol(FWHM=FWHM)
        cb.remove_cube_baseline(n_edgechannel=10)
        #
        d_SO_sum = cb.data.sum(axis=0).copy()
        #mx = np.nanmax(d_SO_sum)
        #d_SO_sum[d_SO_sum > mx * 0.7] = mx*0.7
        slices.append(d_SO_sum)
        thetas.append(cb.theta)
    print 'Making videos'
    slices, thetas = zip(*sorted(zip(slices, thetas), key=lambda x: x[1]))
    cb.data = np.dstack(slices)
    cb.v = np.array(thetas)
    make_a_video(cb, '/n/Users/fdu/now/figures/views.mp4',
                 shapeType=3,
                 labelfmtstr='View angle = {0:.2f} deg')


def task4():
    FWHM = 0.1
    #
    model_dir = '20150418_btfl'
    res_dir = '/u/Moria2/fdu/now/res/'
    #
    C18O_dir = 'C18O'
    #
    linefname_C18O = 'line_00001_00002_00001_2.30537E+11_89.50.fits'
    #
    print 'Loading data'
    fname_C18O = opj(res_dir, model_dir, C18O_dir, 'images', linefname_C18O)
    cb = cube(fname_C18O)
    #
    print 'Convolving'
    cb.convol(FWHM=FWHM)
    #
    print 'Saving a continuum channel'
    cb.draw_a_channel(10, fname_C18O + '_' + model_dir + '_cont.pdf')
    d_cont = cb.data[10,:,:].copy()
    #
    print 'Removing baseline'
    cb.remove_cube_baseline(n_edgechannel=20)
    #
    #cb.trim_tiny_value()
    print 'Drawing total flux'
    cb.draw_sum_channels(fname_C18O + '_' + model_dir + '_sum.pdf')
    d_C18O_sum = cb.data.sum(axis=0).copy()
    #
    #for _ in [d_C18O_sum]:
    #    mx = np.nanmax(_)
    #    _[_ > mx * 0.5] = mx*0.5
    ##
    print 'Making video'
    make_a_video(cb, fname_C18O + '_' + model_dir + '.mp4', ifrRange=(50,150))


def task4a(model_dir = None,
           molecule_dir = None,
           linefname_molecule = None,
           FWHM = 0.5,
        ):
    #
    res_dir = '/n/Users/fdu/now/res/'
    #
    #molecule_dir = 'CO'
    #
    #linefname_molecule = 'line_00001_00003_00001_3.45798E+11_7.00.fits'
    #linefname_molecule = 'line_00002_00005_00001_5.76266E+11_7.00.fits'
    #
    print 'Loading data'
    fname_molecule = opj(res_dir, model_dir, molecule_dir, 'images', linefname_molecule)
    cb = cube(fname_molecule)
    #
    print 'Convolving'
    cb.convol(FWHM=FWHM)
    #
    print 'Saving a continuum channel'
    cb.draw_a_channel(10, fname_molecule + '_' + model_dir + '_cont.pdf')
    d_cont = cb.data[10,:,:].copy()
    #
    print 'Removing baseline'
    cb.remove_cube_baseline(n_edgechannel=4)
    #
    #cb.trim_tiny_value()
    print 'Drawing total flux'
    cb.draw_sum_channels(fname_molecule + '_' + model_dir + '_' + molecule_dir + '_sum.pdf', xRange=(4,-4), yRange=(-4,4), axisUnit='"', gaussianfit=False)
    #cb.draw_peak_channel(fname_molecule + '_' + model_dir + '_peak.pdf')


def task5():
    from draw.long_function_definitions import customaxis
    def do_a_spectrum(fname, model_dir, img_dir, nbin=100):
        print 'Loading data', fname
        cb = cube(fname)
        #
        print 'Removing baseline'
        cb.remove_cube_baseline(n_edgechannel=20)
        #
        print 'Drawing total flux'
        pdfname = fname + '_' + model_dir + '_' + img_dir + '_sum.pdf'
        fig = cb.draw_sum_channels(pdfname, savepdf=False, returnfig=True, logscale=True)
        #
        cb_sum = cb.data.sum(axis=0).reshape(-1)
        #
        xx, yy = np.meshgrid(cb.x, cb.y)
        r = np.sqrt(xx*xx + yy*yy).reshape(-1)
        rmin, rmax = 0.0, cb.x.max()
        rs = np.linspace(rmin, rmax, num=nbin)
        #
        idx = np.digitize(r, rs)
        vals = np.zeros_like(rs)
        for i in np.ndindex(*r.shape):
            ind = idx[i]
            if ind >= len(vals):
                continue
            vals[ind] += cb_sum[i]
        vals /= vals.max()
        accs = np.cumsum(vals)
        accs /= accs.max()
        #
        ax = plt.gca()
        xcen, ycen = fig.transFigure.inverted().transform(ax.transData.transform((0, 0)))
        xmax, _    = fig.transFigure.inverted().transform(ax.transData.transform((cb.x.max(), 0)))
        #
        ax = draw_a_spectrum(rs, vals, fig, xlab='', ylab='', pos=(xcen, ycen, xmax-xcen, 0.2), returnax=True)
        draw_a_spectrum(rs, accs, fig, xlab='', ylab='', ax=ax, color='yellow')
        for perct in [0.2, 0.5, 0.8]:
            i_half = reverse_idx(accs, perct)
            r_half = eval_at_fractional_idx(rs, i_half)
            ax.plot([r_half, r_half], [ax.get_ylim()[0], perct], color='black', alpha=0.5, linestyle=':')
            ax.plot([ax.get_xlim()[0], r_half], [perct, perct],  color='black', alpha=0.5, linestyle=':')
        ax.get_xaxis().set_ticklabels([])
        customaxis(ax, lw=0.5)
        #for child in ax.get_children():
        #    if isinstance(child, mpl.spines.Spine):
        #        child.set_color('#dddddd')
        plt.savefig(pdfname, bbox_inches='tight')
        print 'pdf file saved: ', pdfname
    #
    import glob
    model_dir = '20150609_C2H_5a_4a_2d_8cf'
    #model_dir = '20150307_dep_s'
    res_dir = '/n/Users/fdu/now/res/'
    img_dir_s = ['C18O']
    for img_dir in img_dir_s:
        for fname in glob.glob(opj(res_dir, model_dir, img_dir, 'images', '*.fits')):
            do_a_spectrum(fname, model_dir, img_dir)


def task6():
    FWHM = 1.0
    #
    model_dir = '20150307_dep_s'
    res_dir = '/n/Users/fdu/now/res/'
    #
    C18O_dir = 'CO_img'
    #
    linefname_C18O = 'line_00002_00003_00001_3.45798E+11_0.00.fits'
    #
    print 'Loading data'
    fname_C18O = opj(res_dir, model_dir, C18O_dir, 'images', linefname_C18O)
    cb = cube(fname_C18O)
    #
    print 'Convolving'
    cb.convol(FWHM=FWHM)
    #
    print 'Saving a continuum channel'
    cb.draw_a_channel(10, fname_C18O + '_' + model_dir + '_cont.pdf')
    #
    print 'Removing baseline'
    cb.remove_cube_baseline(n_edgechannel=20)
    #
    #cb.trim_tiny_value()
    #print 'Drawing total flux'
    #cb.draw_sum_channels(fname_C18O + '_' + model_dir + '_tot.pdf')
    #d_C18O_sum = cb.data.sum(axis=0).copy()
    #
    #for _ in [d_C18O_sum]:
    #    mx = np.nanmax(_)
    #    _[_ > mx * 0.5] = mx*0.5
    ##
    print 'Making video'
    make_a_video(cb, fname_C18O + '_' + model_dir + '.mp4', ifrRange=(50,150))


def task7():
    from draw.long_function_definitions import customaxis
    from matplotlib.lines import Line2D
    #
    res_dir = '/n/Users/fdu/now/res/'
    #linefile_s = ['CO_img/images/line_00002_00003_00001_3.45798E+11_0.00.fits',
    #              'oH2O_A_img/images/line_00001_00001_00001_5.56936E+11_0.00.fits',
    #              'HD_img/images/line_00001_00002_00001_2.67499E+12_0.00.fits']
    linefile_s = ['oH2O_A_img/images/line_00001_00001_00001_5.56936E+11_0.00.fits',
                  'HD_img/images/line_00001_00002_00001_2.67499E+12_0.00.fits']
    #linestyles = ['-', '--', ':']
    linestyles = ['-', '--']
    #labels = ['CO ($3-2$)', 'H$_2$O ($1_{{10}}-1_{{01}}$)', 'HD ($1-0$)']
    labels = ['H$_2$O ($1_{{10}}-1_{{01}}$)', 'HD ($1-0$)']
    model_dir_s = ['20150307_dep_s',
                   '20150307_noD_d']
    colors = ['blue', 'red']
    #
    nbin = 100
    winwidth = 3
    figsize = (5,2)
    pos = (0.1,0.1,0.85,0.85)
    fig = plt.figure(figsize = figsize)
    #ax = fig.add_axes(pos, xlabel='r (AU)', ylabel='FLux within r (normalized)')
    ax = fig.add_axes(pos, xlabel='r (AU)', ylabel='dF/dr (normalized)')
    pdfname = opj('/n/Users/fdu/now/', 'cmp_dFdr_vs_radius.pdf')
    #
    xtmp = np.arange(winwidth*4)
    window = np.exp(-0.5*((xtmp-0.5*(xtmp.min()+xtmp.max()))/winwidth)**2)
    window = window / window.sum()
    #
    for iL in xrange(len(linefile_s)):
        for iM in xrange(len(model_dir_s)):
            linefile = opj(res_dir, model_dir_s[iM], linefile_s[iL])
            print 'Working on ', iL, iM, linefile
            cb = cube(linefile)
            print 'Removing baseline'
            cb.remove_cube_baseline(n_edgechannel=20)
            cb_sum = cb.data.sum(axis=0).reshape(-1)
            #
            xx, yy = np.meshgrid(cb.x, cb.y)
            r = np.sqrt(xx*xx + yy*yy).reshape(-1)
            rmin, rmax = 0.0, cb.x.max()
            rs = np.linspace(rmin, rmax, num=nbin)
            #
            idx = np.digitize(r, rs)
            vals = np.zeros_like(rs)
            for i in np.ndindex(*r.shape):
                ind = idx[i]
                if ind >= len(vals):
                    continue
                vals[ind] += cb_sum[i]
            vals = np.convolve(vals, window, mode='same')
            vals /= vals.max()
            #accs = np.cumsum(vals)
            #accs /= accs.max()
            draw_a_spectrum(rs, vals, fig, xlab='', ylab='', ax=ax, color=colors[iM], linestyle=linestyles[iL])
    p = [Line2D([],[], linestyle=ls, color='black') for ls in linestyles] + [Line2D([],[],linestyle='none')] * 2
    L = labels + ['Red: full'] + ['Blue: depleted']
    lgd = ax.legend(p, L, bbox_to_anchor=(0.1, 1.14, 0.7, 0.15), frameon=False,
                    ncol=len(labels), numpoints=1, fontsize=10)
    ax.set_ylim(0, 1.02)
    plt.savefig(pdfname, bbox_inches='tight')
    print 'File saved: ', pdfname


def task8(model_dir=None, FWHM=1.0, C2H_dir='C2H', linefname_C2H_combine = None):
    #
    res_dir = '/n/Users/fdu/now/res/'
    #
    fname_C2H = ''
    cbs = []
    for linename_C2H in linefname_C2H_combine:
        fname_C2H = opj(res_dir, model_dir, C2H_dir, 'images', linename_C2H)
        print 'Loading', fname_C2H
        cbs.append(cube(fname_C2H))
    cb = combine_two_cubes(cbs[0], cbs[1])
    #
    print 'Convolving'
    cb.convol(FWHM=FWHM)
    #
    #print 'Saving a continuum channel'
    #cb.draw_a_channel(10, fname_C2H + '_' + model_dir + '_cont.pdf')
    #
    print 'Removing baseline'
    cb.remove_cube_baseline(n_edgechannel=20)
    #
    #cb.trim_tiny_value()
    print 'Drawing total flux'
    cb.draw_sum_channels(fname_C2H + '_' + model_dir + '_combined_tot.pdf')
    #d_C18O_sum = cb.data.sum(axis=0).copy()
    #
    #for _ in [d_C18O_sum]:
    #    mx = np.nanmax(_)
    #    _[_ > mx * 0.5] = mx*0.5
    ##
    #print 'Making video'
    #make_a_video(cb, fname_C2H + '_' + model_dir + '.mp4', ifrRange=(50,150))


def task9(res_dir = '/n/Users/fdu/now/res/',
          model_dir=None,
          FWHM=1.0,
          fits_dirs=None,
          fname_marker='',
          drawPeak=False,
          drawCont=False, cbarcolor='black',
          fits_filename_patterns=None):
    import os
    import glob
    #
    linefnames = []
    for fits_dir in fits_dirs:
        for pattern in fits_filename_patterns:
            linefnames += glob.glob(os.path.join(res_dir,
                                             model_dir,
                                             fits_dir,
                                             'images',
                                             pattern))
        linefnames = list(set(linefnames))
    #
    iL = 1
    nL = len(linefnames)
    for linefname in linefnames:
        print 'Loading', iL, 'of', nL, linefname
        iL += 1
        cb = cube(linefname)
        #
        print 'Convolving'
        cb.convol(FWHM=FWHM)
        #
        if drawCont:
            print 'Saving a continuum channel'
            cb.draw_a_channel(10, linefname + '_' + model_dir + '_' + fname_marker + '_cont.pdf')
        #
        print 'Removing baseline'
        cb.remove_cube_baseline(n_edgechannel=40)
        #
        fig = plt.figure(figsize = (10, 10))
        #
        fontsize = 10
        cax1 = fig.add_axes([0.14,0.74,0.32,0.03], zorder=10)
        #
        cbar_ori = 'horizontal'
        n_cbar_ticks = 5
        #
        print 'Drawing total flux'
        if drawPeak:
            d = cb.data.max(axis=0)
            cax1.set_title('Jy/beam', fontdict={'fontsize': fontsize})
        else:
            d = cb.data.sum(axis=0) * abs(cb.v[1]-cb.v[0])
            cax1.set_title('Jy/beam km s$^{-1}$', fontdict={'fontsize': fontsize})
        minval1, maxval1 = 0, number_rounding(np.nanmax(d)*1.01)
        #
        cb.draw_sum_channels(linefname + '_' + model_dir + '_' + fname_marker + '_tot.pdf',
            notext=False, noColorbar=False, cbarcolor=cbarcolor,
            vmin=minval1, vmax=maxval1, axisUnit='"',
            cax=cax1, orientation=cbar_ori, cbar_ticks=np.linspace(minval1, maxval1, num=n_cbar_ticks),
            drawPeak=drawPeak,
            returnax=True,
            savepdf=True, fig=fig, pos=(0.1, 0.1, 0.4, 0.8), pos_angle=0.0)


def task10(res_dir = '/n/Users/fdu/now/res/',
          model_dir=None,
          FWHM=0.35,
          fits_dirs=None,
          fits_filename_pattern=None):
    import os
    import glob
    from draw.long_function_definitions import set_axis_format
    #
    #line1_fname = 'line_00001_00018_00001_2.62004E+11_7.00.fits'
    #line2_fname = 'line_00002_00021_00001_2.62004E+11_7.00.fits'
    line1_fname = 'line_00004_00035_00001_3.49399E+11_7.00.fits'
    line2_fname = 'line_00005_00039_00001_3.49402E+11_7.00.fits'
    #
   #linedir_1 = '/n/Users/fdu/now/res/20150307_noD_d_1comp_b/C2H_2_NonLTE/images/'
   #linedir_2 = '/n/Users/fdu/now/res/20150307_noD_d/C2H_2_NonLTE/images/'
 #  linedir_1 = '/n/Users/fdu/now/res/20150307_noD_d_1comp_b/C2H_2_NonLTE_r/images/'
 #  linedir_2 = '/n/Users/fdu/now/res/20150307_noD_d/C2H_2_NonLTE_r/images/'
    linedir_1 = '/n/Users/fdu/now/res/20150307_noD_d_1comp_rerun/C2H_2_NLTE/images/'
  # linedir_2 = '/n/Users/fdu/now/res/20150307_noD_d_rerun2a/C2H_2_NLTE/images/'
    linedir_2 = '/n/Users/fdu/now/res/20151117_a/C2H_2_r/images/'
    #
    fig = plt.figure(figsize = (10,5))
    #
    cbs = []
    for linedir in [linedir_1, linedir_2]:
        cb_tmp = []
        for linename in [line1_fname, line2_fname]:
            linefname = os.path.join(linedir, linename)
            print 'Loading', linefname
            cb_tmp.append(cube(linefname))
        #
        cb = combine_two_cubes(cb_tmp[0], cb_tmp[1])
        #
        print 'Convolving'
        cb.convol(FWHM=FWHM)
        #
        print 'Removing baseline'
        cb.remove_cube_baseline(n_edgechannel=40)
        #
        cbs.append(cb)
    #
    #vmin = min(np.nanmin([cbs[0].data.sum(axis=0) * abs(cbs[0].v[1]-cbs[0].v[0]))
    #           np.nanmin([cbs[1].data.sum(axis=0) * abs(cbs[1].v[1]-cbs[1].v[0])))
    #vmax = max(np.nanmax([cbs[0].data.sum(axis=0) * abs(cbs[0].v[1]-cbs[0].v[0]))
    #           np.nanmax([cbs[1].data.sum(axis=0) * abs(cbs[1].v[1]-cbs[1].v[0])))
    transition_name = 'N=4-3,J=7/2-5/2,F=4-3&3-2'
    fontsize = 13.5
    cax1 = fig.add_axes([0.14,0.74,0.32,0.03], zorder=10)
    cax1.set_title('Jy/beam km s$^{-1}$', fontdict={'fontsize': fontsize})
    cax2 = None #fig.add_axes([0.59,0.74,0.32,0.03], zorder=10)
    #cax2.set_title('Jy/beam km s$^{-1}$', fontdict={'fontsize': fontsize})
    cbar_ori = 'horizontal'
    maxval1, minval1 = 0.1, 0.0
    maxval2, minval2 = 0.1, 0.0
    fig, ax = cbs[0].draw_sum_channels('', notext=True, noColorbar=False, vmin=minval1, vmax=maxval1,
        cax=cax1, orientation=cbar_ori, cbar_ticks=[0.0,0.05,0.1], returnax=True, cbarcolor='white',
        rescaleTo=0.1,
        savepdf=False, fig=fig, pos=(0.1, 0.1, 0.4, 0.8))
    ax.text(0.05, 0.93, 'Smooth disk', transform=ax.transAxes, color='white')
    fig, ax = cbs[1].draw_sum_channels('', notext=True, noColorbar=True, ytextoff=True, vmin=minval2, vmax=maxval2,
        cax=cax2, orientation=cbar_ori, cbar_ticks=[], returnax=True,
        savepdf=False, fig=fig, pos=(0.55, 0.1, 0.4, 0.8))
    ax.text(0.05, 0.93, 'With sharp transition', transform=ax.transAxes, fontsize=fontsize, color='white')
    set_axis_format(cax1, onlyfont=True, labelfontsize=fontsize, tickfontsize=fontsize)
    #set_axis_format(cax2, onlyfont=True, labelfontsize=fontsize, tickfontsize=fontsize)
    set_axis_color(cax1, color='white')
    plt.savefig('/n/Users/fdu/now/smooth_vs_transition_new_rerun1.pdf', bbox_inches='tight')
    

def task11(res_dir = '/n/Users/fdu/now/res/',
          model_dir=None,
          FWHM=0.35,
          fits_dirs=None,
          fits_filename_pattern=None):
    import os
    import glob
    from draw.long_function_definitions import set_axis_format
    #
    #line1_fname = 'line_00001_00018_00001_2.62004E+11_7.00.fits'
    #line2_fname = 'line_00002_00021_00001_2.62004E+11_7.00.fits'
    line1_fname = 'line_00004_00035_00001_3.49399E+11_7.00.fits'
    line2_fname = 'line_00005_00039_00001_3.49402E+11_7.00.fits'
    #
    linedir_s = [\
        '/n/Users/fdu/now/res/20150307_noD_d_highC2H_1E3yr/C2H_2_NonLTE_r/images/',
        '/n/Users/fdu/now/res/20150307_noD_d_highC2H_1E4yr/C2H_2_NonLTE_r/images/',
        '/n/Users/fdu/now/res/20150307_noD_d_highC2H_1E5yr/C2H_2_NonLTE_r/images/',
        ]
    #
    fig = plt.figure(figsize = (12,5))
    #
    cbs = []
    for linedir in linedir_s:
        cb_tmp = []
        for linename in [line1_fname, line2_fname]:
            linefname = os.path.join(linedir, linename)
            print 'Loading', linefname
            cb_tmp.append(cube(linefname))
        #
        cb = combine_two_cubes(cb_tmp[0], cb_tmp[1])
        #
        print 'Convolving'
        cb.convol(FWHM=FWHM)
        #
        print 'Removing baseline'
        cb.remove_cube_baseline(n_edgechannel=20)
        #
        cbs.append(cb)
    #
    #vmin = min(np.nanmin([cbs[0].data.sum(axis=0) * abs(cbs[0].v[1]-cbs[0].v[0]))
    #           np.nanmin([cbs[1].data.sum(axis=0) * abs(cbs[1].v[1]-cbs[1].v[0])))
    #vmax = max(np.nanmax([cbs[0].data.sum(axis=0) * abs(cbs[0].v[1]-cbs[0].v[0]))
    #           np.nanmax([cbs[1].data.sum(axis=0) * abs(cbs[1].v[1]-cbs[1].v[0])))
    transition_name = 'N=4-3,J=7/2-5/2,F=4-3&3-2'
    fontsize = 12
    cax = fig.add_axes([0.42,0.67,0.22,0.03], zorder=10)
    cax.set_title('Jy/beam km s$^{-1}$', fontdict={'fontsize': fontsize})
    cbar_ori = 'horizontal'
    maxval, minval = 0.4, 0.0
    fig, ax = cbs[0].draw_sum_channels('', notext=True, noColorbar=False, vmin=minval, vmax=maxval,
        cax=cax, orientation=cbar_ori, cbar_ticks=[0,0.2,0.4], returnax=True, cbarcolor='white',
        savepdf=False, fig=fig, pos=(0.1, 0.1, 0.26, 0.8))
    ax.text(0.05, 0.93, 't = $10^3$ yr', transform=ax.transAxes, color='white')
    fig, ax = cbs[1].draw_sum_channels('', notext=True, noColorbar=True, ytextoff=True, vmin=minval, vmax=maxval, returnax=True,
        savepdf=False, fig=fig, pos=(0.4, 0.1, 0.26, 0.8))
    ax.text(0.05, 0.93, 't = $10^4$ yr', transform=ax.transAxes, color='white')
    fig, ax = cbs[2].draw_sum_channels('', notext=True, noColorbar=True, ytextoff=True,
        vmin=minval, vmax=maxval, returnax=True,
        savepdf=False,  fig=fig, pos=(0.7, 0.1, 0.26, 0.8))
    ax.text(0.05, 0.93, 't = $10^5$ yr', transform=ax.transAxes, fontsize=fontsize, color='white')
    set_axis_format(cax, onlyfont=True, labelfontsize=fontsize, tickfontsize=fontsize)
    set_axis_color(cax, color='white')
    plt.savefig('/n/Users/fdu/now/fig1_time1E31E41E5.pdf', bbox_inches='tight')



def task11d(res_dir = '/n/Users/fdu/now/res/',
          model_dir=None,
          FWHM=0.35,
          fits_dirs=None,
          fits_filename_pattern=None):
    import os
    import glob
    from draw.long_function_definitions import set_axis_format
    #
    #line1_fname = 'line_00001_00018_00001_2.62004E+11_7.00.fits'
    #line2_fname = 'line_00002_00021_00001_2.62004E+11_7.00.fits'
    line1_fname = 'line_00004_00035_00001_3.49399E+11_7.00.fits'
    line2_fname = 'line_00005_00039_00001_3.49402E+11_7.00.fits'
    #
    linedir_s = [\
        '/n/Users/fdu/now/res/20150307_noD_d_highC2H_1E3yr/C2H_2_NonLTE_r/images/',
        '/n/Users/fdu/now/res/20150307_noD_d_highC2H_1E4yr/C2H_2_NonLTE_r/images/',
        '/n/Users/fdu/now/res/20150307_noD_d_highC2H_1E5yr/C2H_2_NonLTE_r/images/',
        ]
    #
    fig = plt.figure(figsize = (12,5))
    #
    cbs = []
    for linedir in linedir_s:
        cb_tmp = []
        for linename in [line1_fname, line2_fname]:
            linefname = os.path.join(linedir, linename)
            print 'Loading', linefname
            cb_tmp.append(cube(linefname))
        #
        cb = combine_two_cubes(cb_tmp[0], cb_tmp[1])
        #
        print 'Convolving'
        cb.convol(FWHM=FWHM)
        #
        print 'Removing baseline'
        cb.remove_cube_baseline(n_edgechannel=2)
        #
        cbs.append(cb)
    #
    transition_name = 'N=4-3,J=7/2-5/2,F=4-3&3-2'
    fontsize = 12

    cax = fig.add_axes([0.42,0.67,0.22,0.03], zorder=10)
    cax.set_title('Jy/beam', fontdict={'fontsize': fontsize})
    set_axis_format(cax, onlyfont=True, labelfontsize=fontsize, tickfontsize=fontsize)
    set_axis_color(cax, color='white')

    cbar_ori = 'horizontal'
    maxval, minval = 0.5, 0.0

    fig, ax = cbs[0].draw_sum_channels('', notext=True, noColorbar=False, vmin=minval, vmax=maxval, drawPeak=True,
        cax=cax, orientation=cbar_ori, cbar_ticks=[0,0.25,0.5], returnax=True, cbarcolor='white',
        savepdf=False, fig=fig, pos=(0.1, 0.1, 0.26, 0.8))
    ax.text(0.05, 0.93, 't = $10^3$ yr', transform=ax.transAxes, color='white')

    fig, ax = cbs[1].draw_sum_channels('', notext=True, noColorbar=True, ytextoff=True, vmin=minval, vmax=maxval, returnax=True, drawPeak=True,
        savepdf=False, fig=fig, pos=(0.4, 0.1, 0.26, 0.8))
    ax.text(0.05, 0.93, 't = $10^4$ yr', transform=ax.transAxes, color='white')

    fig, ax = cbs[2].draw_sum_channels('', notext=True, noColorbar=True, ytextoff=True, drawPeak=True,
        vmin=minval, vmax=maxval, returnax=True,
        savepdf=False,  fig=fig, pos=(0.7, 0.1, 0.26, 0.8))
    ax.text(0.05, 0.93, 't = $10^5$ yr', transform=ax.transAxes, fontsize=fontsize, color='white')

    plt.savefig('/n/Users/fdu/now/fig1_time1E31E41E5_peak.pdf', bbox_inches='tight')




def get_radial_profile(cb, normalize=True):
    d = cb.data.sum(axis=0)


def task11a(res_dir = '/n/Users/fdu/now/res/',
          model_dir=None,
          FWHM=0.35,
          fits_dirs=None,
          fits_filename_pattern=None):
    import os
    import glob
    from draw.long_function_definitions import set_axis_format
    #
    #line1_fname = 'line_00001_00030_00001_3.49312E+11_7.00.fits'
    #line2_fname = 'line_00002_00029_00001_3.49339E+11_7.00.fits'
    line1_fname = 'line_00004_00035_00001_3.49399E+11_7.00.fits'
    line2_fname = 'line_00005_00039_00001_3.49402E+11_7.00.fits'
    #
    linedir_s = [\
        '/n/Users/fdu/now/res/20151117_a_dep2/C2H_2_r/images/',
        '/n/Users/fdu/now/res/20151117_a_dep1a/C2H_2_r/images/',
        '/n/Users/fdu/now/res/20151117_a/C2H_2_r/images/',
        #'/n/Users/fdu/now/res/20151117_a_iniC/C2H_2/images/',
        #'/n/Users/fdu/now/res/20150609_C2H_5a_2a/C2H_2/images/',
        #'/n/Users/fdu/now/res/20150609_C2H_5a_4a/C2H_2/images/',
        #'/n/Users/fdu/now/res/20150609_C2H_5a_nodep/C2H_2/images/',
        ]
    #
    fig = plt.figure(figsize = (12,5))
    #
    cbs = []
    for linedir in linedir_s:
        cb_tmp = []
        for linename in [line1_fname, line2_fname]:
            linefname = os.path.join(linedir, linename)
            print 'Loading', linefname
            cb_tmp.append(cube(linefname))
        #
        cb = combine_two_cubes(cb_tmp[0], cb_tmp[1])
        #
        print 'Convolving'
        cb.convol(FWHM=FWHM)
        #
        print 'Removing baseline'
        cb.remove_cube_baseline(n_edgechannel=20)
        #
        cbs.append(cb)
    #
    transition_name = 'N=4-3,J=7/2-5/2,F=4-3&3-2'
    fontsize = 12
    cax = fig.add_axes([0.42,0.67,0.22,0.03], zorder=10)
    cax.set_title('Jy/beam km s$^{-1}$', fontdict={'fontsize': fontsize})
    set_axis_format(cax, onlyfont=True, labelfontsize=fontsize, tickfontsize=fontsize)
    set_axis_color(cax, color='white')
    cbar_ori = 'horizontal'
    maxval, minval = 0.3, 0.0

    fig, ax = cbs[0].draw_sum_channels('', notext=True, noColorbar=False, vmin=minval, vmax=maxval,
        cax=cax, orientation=cbar_ori, cbar_ticks=[0,0.1,0.2,0.3], returnax=True, cbarcolor='white',
        savepdf=False, fig=fig, pos=(0.7, 0.1, 0.26, 0.8), ytextoff=True)
    ax.text(0.05, 0.93, 'C/O = 2', transform=ax.transAxes, color='white')

    fig, ax = cbs[1].draw_sum_channels('', notext=True, noColorbar=True, ytextoff=True, vmin=minval, vmax=maxval, returnax=True,
        savepdf=False, fig=fig, pos=(0.4, 0.1, 0.26, 0.8))
    ax.text(0.05, 0.93, 'C/O = 1', transform=ax.transAxes, color='white')

    fig, ax = cbs[2].draw_sum_channels('', notext=True, noColorbar=True, ytextoff=False,
        vmin=minval, vmax=maxval, returnax=True,
        savepdf=False,  fig=fig, pos=(0.1, 0.1, 0.26, 0.8))
    ax.text(0.05, 0.93, 'C/O = 0.4', transform=ax.transAxes, fontsize=fontsize, color='white')

    plt.savefig('/n/Users/fdu/now/figDep_uniform_4332_r_pos.pdf', bbox_inches='tight')


def task11b(res_dir = '/n/Users/fdu/now/res/',
          FWHM=0.7,
          fits_dirs=None,
          fits_filename_pattern=None):
    import os
    import glob
    from draw.long_function_definitions import set_axis_format
    from draw.long_function_definitions import calc_colden
    from draw.long_function_definitions import load_data_as_dic
    from draw.long_function_definitions import plot_model_results_radial
    #
    line1_fname = 'line_00004_00035_00001_3.49399E+11_7.00.fits'
    model_s = [\
                ('/n/Users/fdu/now/res/20151117_a_dep2/',
                 '/n/Users/fdu/now/res/20151117_a_dep2/C2H_2_r/images/'),
                ('/n/Users/fdu/now/res/20151117_a_dep1a/',
                 '/n/Users/fdu/now/res/20151117_a_dep1a/C2H_2_r/images/'),
                ('/n/Users/fdu/now/res/20151117_a/',
                 '/n/Users/fdu/now/res/20151117_a/C2H_2_r/images/')
                #'/n/Users/fdu/now/res/20151117_a_iniC/',
                #'/n/Users/fdu/now/res/20150609_C2H_5a_2a/',
                #'/n/Users/fdu/now/res/20150609_C2H_5a_4a/',
                #'/n/Users/fdu/now/res/20150609_C2H_5a_nodep/',
        ]
    #
    fig = plt.figure(figsize = (8,5))
    #
    d_s = []
    rp_s = []
    for model, linedir in model_s:
        fname = os.path.join(model, 'iter_0001.dat')
        print 'Loading ', fname
        d = load_data_as_dic(fname)
        calc_colden(d, 'C2H')
        #
        linefile = os.path.join(linedir, line1_fname)
        cb = cube(linefile)
        rprof = get_radial_profile(cb, normalize)
        #
        d_s.append(d)
        rp_s.append(rprof)
    #
    fontsize = 12
    zpos = 0.0
    lis = [[{'name': 'N(C2H)', 'vr': (1e10, 1e17)}]]
    
    ax = fig.add_axes([0.1,0.1,0.86,0.95], zorder=10,
                      xlabel='r (AU)', ylabel='Column density (cm$^{-2}$)',
                      xscale='linear', yscale='log', ylim=(1e10,1e17))
    color_s = ['blue', 'red', 'green']
    label_s = ['C/O = 2', 'C/O = 1', 'C/O = 0.4']
    for i in xrange(3):
        plot_model_results_radial(d_s[i], lis, zpos, ax = ax,
            returndata=False, noLegend=False, noTitle=True,
            color=color_s[i], label=label_s[i], display_dic=None,
            returnaxis=False, linestyle='-',lw=3, markersize=5,
            do_interpol=True, interpolNum=200,
            dosmooth=True, winwidth=10,
            whichx='rmin')
    set_axis_format(ax, labelfontsize=fontsize, tickfontsize=fontsize, majorgridon=False, minorgridon=False)

    plt.savefig('/n/Users/fdu/now/figDep_Ncol_uniform_a.pdf', bbox_inches='tight')


def task11c(res_dir = '/n/Users/fdu/now/res/',
          FWHM=0.7,
          fits_dirs=None,
          fits_filename_pattern=None):
    import os
    import glob
    from draw.long_function_definitions import set_axis_format
    from draw.long_function_definitions import calc_colden
    from draw.long_function_definitions import load_data_as_dic
    from draw.long_function_definitions import plot_model_results_radial
    #
    model_s = [\
                '/n/Users/fdu/now/res/20151117_a/',
        ]
    #
    fig = plt.figure(figsize = (8,5))
    #
    d_s = []
    for model in model_s:
        fname = os.path.join(model, 'iter_0001.dat')
        print 'Loading ', fname
        d = load_data_as_dic(fname)
        calc_colden(d, 'rhodus_1', is_abundance=False)
        calc_colden(d, 'rhodus_2', is_abundance=False)
        calc_colden(d, 'rhodus_3', is_abundance=False)
        #
        d_s.append(d)
    #
    fontsize = 12
    zpos = 0.0
    #
    ax = fig.add_axes([0.1,0.1,0.86,0.95], zorder=10,
                      xlabel='r (AU)', ylabel='Column density (g cm$^{-2}$)',
                      xscale='linear', yscale='log', ylim=(1e-6, 1e2))
    color_s = ['blue', 'red', 'green']
    label_s = ['mm dust', '$\mu$m dust (outer)', '$\mu$m dust (inner)']
    name_s = ['N(rhodus_1)', 'N(rhodus_2)', 'N(rhodus_3)']
    for i in xrange(3):
        lis = [[{'name': name_s[i], 'vr': (1e-6, 1e6)}]]
        plot_model_results_radial(d_s[0], lis, zpos, ax = ax,
            returndata=False, noLegend=False, noTitle=True,
            color=color_s[i], label=label_s[i], display_dic=None,
            returnaxis=False, linestyle='-',lw=3, markersize=5,
            do_interpol=True, interpolNum=200,
            dosmooth=False, winwidth=3,
            whichx='rmin')
    set_axis_format(ax, labelfontsize=fontsize, tickfontsize=fontsize, majorgridon=False, minorgridon=False)
    #
    plt.savefig('/n/Users/fdu/now/fig_dust_sigma.pdf', bbox_inches='tight')


def task12(linedir_s,
          res_dir = '/n/Users/fdu/now/res/',
          pos_angle = (None, -45.0),
          #FWHM=(0.76, 0.46),
          FWHM=(0.76, 0.38),
          #FWHM=(0.76, 0.2),
          fits_dirs=None,
          fits_filename_pattern=None,
          pdfname = 'fig'):
    import os
    import glob
    from draw.long_function_definitions import set_axis_format
    #
    line_fname_s = [[\
                     'line_00004_00035_00001_3.49399E+11_7.00.fits', # 4_3.5_4 -> 3_2.5_3
                     'line_00005_00039_00001_3.49402E+11_7.00.fits', # 4_3.5_3 -> 3_2.5_2
                    ],
                    [\
                     'line_00001_00018_00001_2.62004E+11_32.00.fits', # 3_3.5_3 -> 2_2.5_2
                     'line_00002_00021_00001_2.62004E+11_32.00.fits', # 3_3.5_4 -> 2_2.5_3
                    ]]
    #
    fig = plt.figure(figsize = (10,5))
    #
    cbs = []
    for i in xrange(2):
        cb_tmp = []
        linedir = linedir_s[i]
        for linename in line_fname_s[i]:
            linefname = os.path.join(linedir, linename)
            print 'Loading', linefname
            cb_tmp.append(cube(linefname))
        #
        cb = combine_two_cubes(cb_tmp[0], cb_tmp[1])
        #
        print 'Convolving'
        cb.convol(FWHM=FWHM[i])
        #
        cb.draw_a_channel(10, '/n/Users/fdu/now/' + pdfname + '_continuum.pdf')
        #
        print 'Removing baseline'
        cb.remove_cube_baseline(n_edgechannel=20)
        #
        cbs.append(cb)
    #
    #vmin = min(np.nanmin([cbs[0].data.sum(axis=0) * abs(cbs[0].v[1]-cbs[0].v[0]))
    #           np.nanmin([cbs[1].data.sum(axis=0) * abs(cbs[1].v[1]-cbs[1].v[0])))
    #vmax = max(np.nanmax([cbs[0].data.sum(axis=0) * abs(cbs[0].v[1]-cbs[0].v[0]))
    #           np.nanmax([cbs[1].data.sum(axis=0) * abs(cbs[1].v[1]-cbs[1].v[0])))
    #
    fig = plt.figure(figsize = (10,5))
    fontsize = 10
    cax1 = fig.add_axes([0.14,0.74,0.32,0.03], zorder=10)
    cax1.set_title('Jy/beam km s$^{-1}$', fontdict={'fontsize': fontsize})
    cax2 = fig.add_axes([0.59,0.74,0.32,0.03], zorder=10)
    cax2.set_title('Jy/beam km s$^{-1}$', fontdict={'fontsize': fontsize})
    cbar_ori = 'horizontal'
    maxval1, minval1 = 1.0, 0.0
    maxval2, minval2 = 0.1, 0.0
    n_cbar_ticks = 5
    #
    fig, ax = cbs[0].draw_sum_channels('', notext=True, noColorbar=False, vmin=minval1, vmax=maxval1, axisUnit='"',
        cax=cax1, orientation=cbar_ori, cbar_ticks=np.linspace(minval1, maxval1, num=n_cbar_ticks),
        returnax=True, #[0,0.1,0.2,0.3,0.4,0.5]
        savepdf=False, fig=fig, pos=(0.1, 0.1, 0.4, 0.8), pos_angle=pos_angle[0])
    ax.text(0.05, 0.93, 'TW Hya', transform=ax.transAxes, fontsize=fontsize)
    #
    fig, ax = cbs[1].draw_sum_channels('', notext=True, noColorbar=False, ytextoff=False, vmin=minval2, vmax=maxval2, axisUnit='"',
        cax=cax2, orientation=cbar_ori, cbar_ticks=np.linspace(minval2, maxval2, num=n_cbar_ticks),
        returnax=True, savepdf=False, fig=fig, pos=(0.55, 0.1, 0.4, 0.8), pos_angle=pos_angle[1])
    ax.set_ylabel('')
    ax.text(0.05, 0.93, 'DM Tau', transform=ax.transAxes, fontsize=fontsize)
    #
    set_axis_format(cax1, onlyfont=True, labelfontsize=fontsize, tickfontsize=fontsize)
    set_axis_format(cax2, onlyfont=True, labelfontsize=fontsize, tickfontsize=fontsize)
    plt.savefig('/n/Users/fdu/now/' + pdfname + '.pdf', bbox_inches='tight')



def task12a(linefname,
          pdfname = '/n/Users/fdu/now/fig_continuum.pdf',
          xRange=None, yRange=None,
          pos_angle = None,
          FWHM=None):
    from scipy.ndimage.interpolation import rotate
    #
    fig = plt.figure(figsize = (6,6))
    #
    cb = cube(linefname)
    #
    #if pos_angle != None:
    #    cb.data = rotate(cb.data, pos_angle, reshape=False)
    #
    print 'Convolving'
    cb.convol(FWHM=FWHM)
    #
    cb.draw_a_channel(10, pdfname, xRange=xRange, yRange=yRange)
    print 'File saved: ', pdfname



def plot_SED_match(f_s, ax = None):
    from draw.long_function_definitions import jy2nuFnu
    for f in f_s:
        d = np.loadtxt(f['fname'])
        ax.plot(d[:,0]*1e-4, jy2nuFnu(d[:,1], d[:, 0]*1e-4),
                lw=2, label=f['label'])
    return

def plot_SED_adhoc(f_s = None, fsize = (16,12), save_pdf = False,
                   SED_cont = None, SED_points = None, SED_star = None, dist_pc = None,
                   xscale='log', yscale='log',
                   xlim=(1e-1,1e4), ylim=(1e-14,1e-8)):

    from draw.long_function_definitions import set_axis_format

    fig = plt.figure(figsize=fsize)
    ax = fig.add_axes((0.1,0.1,0.85,0.85))

    if SED_points != None:
        ax.plot(SED_points[:, 0], SED_points[:, 1], lw=3, color=(0.0,0.0,0.0), marker='o')

    if SED_cont != None:
        ax.plot(SED_cont[:, 0], SED_cont[:, 2], lw=3, color=(0.0,0.0,0.0))

    if SED_star != None:
        pc = 3.1e18
        dist = dist_pc * pc
        factor = 1.0/(4*pi*dist**2)
        ax.plot(SED_star[:,1]*1e-4, SED_star[:,1]*SED_star[:,2]*factor, color=(0.0,0.0,0.0))

    plot_SED_match(f_s, ax = ax)

    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    ax.set_xlabel(r'$\lambda$ ($\mu$m)')
    ax.set_ylabel(r'$\nu F_\nu$ (erg s$^{-1}$ cm$^{-2}$)')
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    
    set_axis_format(ax, graygrid=True)
    ax.legend(fontsize=20, framealpha=0.6, columnspacing=1, handletextpad=0, loc='best')

    if save_pdf:
        s = f_s[0]['fname']
        pdfname = s[:s.rfind('.')] + '.pdf'
        print pdfname
        plt.savefig(pdfname, bbox_inches='tight')
    return


def task13(f_s = None, save_pdf = False):
    #f_s = [{'fname':'/Users/fdu/now/res/20140924_a0/mc/spec_ang_4.00_10.00_iter_0001.dat',
    #        'label':'0924a0(1)'},
    #       {'fname':'/Users/fdu/now/res/20140929_a1_K_n2/mc/spec_ang_4.00_10.00_iter_0001.dat',
    #        'label':'20140929a1Kn2(1)'},
    #       ]
    plot_SED_adhoc(f_s, fsize=(10, 6), xlim=(3e-1, 2e3), ylim=(1e-13, 3e-9), save_pdf = save_pdf)


def rescale(y, vrange):
    ymin, ymax = np.nanmin(y), np.nanmax(y)
    return (y-ymin) * ((vrange[1] - vrange[0]) / (ymax-ymin)) + vrange[0]


def draw_single_spec(ax, x, y, removeBase=False, rescaleTo=None):
    #xmin, xmax = np.nanmin(x), np.nanmax(x)
    #ymin, ymax = np.nanmin(y), np.nanmax(y)
    #ax.set_xlim((xmin, xmax))
    #ax.set_ylim((ymin, ymax))
    ax.set_xlim((-6,6))
    ax.set_ylim((-5e-4,2e-4))
    if removeBase:
        k = (y[-1] - y[0]) / (x[-1] - x[0])
        ybase = k * (x - x[0]) + y[0]
        ynew = y - ybase
    else:
        ynew = y
    if rescaleTo != None:
        ynew = rescale(ynew, rescaleTo)
    ax.plot(x, ynew)
    return (np.nanmin(ynew), np.nanmax(ynew))


def iter_len(it):
    i = 0
    for _ in it:
        i += 1
    return i


def get_idx_range(vx, x_range):
    return [i for i in xrange(len(vx)) if x_range[0] <= vx[i] <= x_range[1]]


def task14(modelname = '20150609_C2H_5a_4a_2d_8cg',
           line_dir = 'C2H_2_LTE',
           linefname = 'line_00001_00030_00001_3.49312E+11_7.00.fits',
           res_dir = '/n/Users/fdu/now/res/',
           xRange=None, yRange=None,
           pos_angle = None,
           x_range = None,
           y_range = None,
           v_range = None,
           FWHM=None,
           removeBase=True,
           rescaleTo=None,
           pdfname=None,
          ):
    from scipy.ndimage.interpolation import rotate
    from draw.c2h_figures import my_axes_grid
    #
    line_path = opj(res_dir, modelname, line_dir, 'images', linefname)
    cb = cube(line_path)
    #
    print 'Convolving'
    if FWHM != None:
        cb.convol(FWHM=FWHM)
    #
    print 'Getting ranges'
    x_indices = get_idx_range(cb.x, x_range)
    y_indices = get_idx_range(cb.y, y_range)
    #
    nx = len(x_indices)
    ny = len(y_indices)
    #
    print nx, ny
    print 'Making axis grid'
    ag = my_axes_grid(ny, nx, fig=None, figsize=(30,30), xpad=0, ypad=0, aspect=None,
                      xlim=None, ylim=None, xlabel='', ylabel='')
    print 'Plotting'
    ymin, ymax = 1e99, -1e99
    for i in xrange(nx):
        ii = x_indices[i]
        for j in xrange(ny):
            jj = y_indices[j]
            print i, "/", nx, j, "/", ny
            ax = ag[i][j]
            #ax.set_yscale('log')
            thismin, thismax = draw_single_spec(ax, cb.v, cb.data[:, ii, jj], removeBase=removeBase, rescaleTo=rescaleTo)
            ymin, ymax = min(ymin, thismin), max(ymax, thismax)
    for i in xrange(nx):
        ii = x_indices[i]
        for j in xrange(ny):
            jj = y_indices[j]
            ax = ag[i][j]
            if v_range == None:
                ax.set_ylim((ymin, ymax))
            else:
                ax.set_ylim(v_range)
    #
    plt.savefig(pdfname, bbox_inches='tight')
    print 'Pdf saved:', pdfname

def task15():
    FWHM = 0.5
    #
    print 'Loading data'
    cb = cube('/n/Users/fdu/now/res/20150609_C2H_5a_4a_2d_8cl_B6k/C2H_1_NLTE/images/line_00001_00018_00001_2.62004E+11_7.00.fits')
    #
    print 'Convolving'
    cb.convol(FWHM=FWHM)
    cb.remove_cube_baseline(n_edgechannel=10)
    #
    fig = plt.figure(figsize = (10, 10))
    #
    fontsize = 10
    cax1 = fig.add_axes([0.14,0.74,0.32,0.03], zorder=10)
    cax1.set_title('Jy/beam km s$^{-1}$', fontdict={'fontsize': fontsize})
    #
    minval1, maxval1 = 0.0, 0.1
    cbar_ori = 'horizontal'
    n_cbar_ticks = 5
    #
    cb.draw_sum_channels('/n/Users/fdu/now/C2H_1_1.pdf', notext=True, noColorbar=False,
        vmin=minval1, vmax=maxval1, axisUnit='"',
        cax=cax1, orientation=cbar_ori, cbar_ticks=np.linspace(minval1, maxval1, num=n_cbar_ticks),
        returnax=True, #[0,0.1,0.2,0.3,0.4,0.5]
        savepdf=True, fig=fig, pos=(0.1, 0.1, 0.4, 0.8), pos_angle=0.0)


def task16(FWHM='',
           fits_fname='',
           pdfname='',
           axisUnit='"',
           tagstr='',
          ):
    cb = cube(fits_fname, saveTauMap=True)
    cb.convol(FWHM=FWHM, tau_also=False)
    #
    fig = plt.figure(figsize = (10, 10))
    #
    fontsize = 10
    cax1 = fig.add_axes([0.14,0.74,0.32,0.03], zorder=10)
    cax1.set_title(r'$\tau$', fontdict={'fontsize': fontsize})
    #
    minval1, maxval1 = 0.0, 500.0
    cbar_ori = 'horizontal'
    n_cbar_ticks = 5
    #
    draw_an_image(cb.tauMapData,
                  cb.x / (1.0 if axisUnit == 'AU' else cb.dist),
                  cb.y / (1.0 if axisUnit == 'AU' else cb.dist),
                  pdfname,
                  axisUnit = axisUnit,
                  tag = tagstr,
                  returnfig=False, returnax=False,
                  savepdf=True,
                  pos=(0.1,0.1,0.85,0.85), fig=fig, ax=None,
                  vmin=minval1, vmax=maxval1,
                  noColorbar=False, ytextoff=False, cbarcolor='black',
                  cax=cax1, orientation='horizontal',
                  cbar_ticks=np.linspace(minval1, maxval1, num=n_cbar_ticks),
                  pos_angle=None)


if __name__ == '__main__':
   #task16(
   #        FWHM=0.1,
   #        fits_fname='/n/Users/fdu/now/res/20150609_C2H_5a_4a_2d_8cl_B0/13CO_redu_1/images/line_00001_00001_00001_1.15271E+11_7.00.fits',
   #        #fits_fname='/n/Users/fdu/now/res/20150609_C2H_5a_4a_2d_8cl_B0/C18O_redu_1/images/line_00001_00001_00001_1.15271E+11_7.00.fits',
   #        pdfname='/n/Users/fdu/now/13CO_tau.pdf',
   #        #pdfname='/n/Users/fdu/now/C18O_tau.pdf',
   #        axisUnit='"',
   #        tagstr=''
   #      )
    ##for L in [
    ##            {
    ##             'modelname': '20150609_C2H_5a_4a_2d_8cl_B6l_41_g',
    ##             'line_dir':  'C18O_contri',
    ##             'linefname': 'line_00003_00003_00001_3.45798E+11_7.00.fits',
    ##             'pdfname':   '/n/Users/fdu/now/specGrid_C18O_20150609_C2H_5a_4a_2d_8cl_B6l_41_g_3.45E11.pdf',
    ##            },
    ##         ]:
    ##    task14(x_range=(-30,30), y_range=(-30,30), **L)
    ##raise Exception("Stop")
    #######################################
    #
 #  task11a()
 #  task11d()
    #
    for model in [
                    '20160714_TWH_d_s_noDep_moreOuterMass1',
                 ]:
        task9(model_dir=model,
              FWHM=0.4,
              #fits_dirs=['C2H_1_NLTE'],
              #fits_dirs=['C18O', 'oC3H2_NLTE', 'C2H_2_NLTE', ],
              #fits_dirs=['C18O', 'C18O_rerun', 'C18O_rerun1', 'C2H_2_NLTE_contri', 'C2H_2_NLTE'],
              #fits_dirs=['C18O'],
              fits_dirs=['C18O'],
              fname_marker='',
              drawPeak=False,
              drawCont=False,
              #fits_filename_patterns=['*1.15271E+11*.fits','*2.30537E+11*.fits'])
              fits_filename_patterns=['*3.45798E+11*.fits','*6.91475E+11*.fits'])
              #fits_filename_patterns=['*.fits'])
    #task5()
 #  for mdir in [\
 #                  '20150609_C2H_5a_4a_2d_8cl_B6k',
 #              ]:
 #      task4a(model_dir = mdir,
 #             molecule_dir = 'C18O_1',
 #             linefname_molecule = 'line_00001_00001_00001_1.15271E+11_7.00.fits',
 #             #linefname_molecule = 'line_00001_00003_00001_3.45798E+11_32.00.fits',
 #             FWHM = 0.1,
 #           )
    ####
 #  for mdir in [ \
 #                  '20150609_C2H_5a_4a_2d_8cl_B6k',
 #              ]:
 #      linedir_s = ['/n/Users/fdu/now/res/20150609_C2H_5a/C2H_2/images/',
 #                   opj('/n/Users/fdu/now/res/', mdir, 'C2H_1/images/')]
 #      task12(linedir_s, pdfname = mdir)

    #    mc_file = opj('/n/Users/fdu/now/res/', mdir, 'mc', 'spec_ang_27.00_37.00_iter_0001.dat')

    #    task13(f_s = [
    #            {'fname': mc_file,
    #             'label': mdir},
    #           ], save_pdf = True)
    ####
    #task12a('/n/Users/fdu/now/res/20150609_C2H_5a_4a_2d_5/C2H_2/images/line_00001_00030_00001_3.49312E+11_7.00.fits',
    #        xRange=(-2, 2), yRange=(-2,2),
    #        pos_angle=0, FWHM=0.41)
    #
#   task10()
#   task11()
    #
    #task11b()
    #task11c()
#   for mdir in [ \
#               '20150609_DMT_C2H',
#               '20150609_DMT_C2H_1',
#               '20150609_DMT_C2H_2',
#               '20150609_DMT_C2H_3',
#               '20150609_DMT_C2H_4',
#           ]:
#       task8(FWHM = 0.5,
#           model_dir = mdir,
#           linefname_C2H_combine = ('line_00001_00018_00001_2.62004E+11_32.00.fits',
#                                    'line_00002_00021_00001_2.62004E+11_32.00.fits'))
    # 23.694, 23.72
