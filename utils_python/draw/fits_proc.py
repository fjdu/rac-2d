import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt



class specline:
    __fig__ = None
    def __init__(self, fname):
        from astropy.io import fits
        from numpy import arange
        h = fits.open(fname)
        n_ext = len(h)
        c = 2.99792458e8
        self.found = False
        iext = 0
        self.headerfound = False
        self.obs_data = None
        for ext in h:
            hd = ext.header
            if hd['EXTNAME'] == 'FluxSpec' or hd['EXTNAME'] == 'LineCube':
                if not self.headerfound:
                    self.headerfound = True
                    self.iext = iext
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
                    #
                    try:
                        self.molname = hd['MOL']
                        self.molname_disp = hd['MOL']
                    except:
                        self.molname = hd['MOL-DB']
                        self.molname_disp = hd['MOL-DSP']
                    if self.molname_disp.strip() == '':
                        self.molname_disp = self.molname
                    #
                    self.dist = hd['DIST']
                    self.theta = hd['THETA']
                    self.lam0 = hd['LAM0']
                    if 'QNUM' in hd.keys():
                        self.qnum = hd['QNUM']
                    else:
                        self.qnum = ''
                #
                if hd['EXTNAME'] == 'FluxSpec':
                    self.spec = ext.data[0]
                    self.n = hd['NAXIS1']
                    self.k0 = hd['CRPIX1']
                    if hd['CTYPE1'] == 'V':
                        self.dv = hd['CDELT1']
                        self.v0 = hd['CRVAL1']
                        self.v = (arange(1, n+1) - self.k0) * dv + v0
                        self.f = f0 * (1 - self.v/c)
                        self.df = abs(self.f[1] - self.f[0])
                    elif hd['CTYPE1'] == 'F':
                        self.fmin = hd['CRVAL1']
                        self.df = hd['CDELT1']
                        self.f = (arange(1, self.n+1) - self.k0) * self.df + self.fmin
                    else:
                        raise "Unsupported CTYPE1: " + hd['CTYPE1']
                    self.lam = c / self.f * 1e6
                    if not hasattr(self, 'v'):
                        self.v = (self.f0 - self.f) * (c / self.f0)
                    self.found = True
                #break
            else:
                iext += 1
        h.close()
        if not self.found:
            raise 'No spectrum data found in the fits file!\nFilename: ' + fname + '\n'
        return

    def qnum_reformat(self):
        import re
        if '%' in self.qnum:
            idx = self.qnum.index('%')
            s = self.qnum[idx+1:].strip()
        else:
            s = self.qnum.strip()
        s = re.sub("\s*\d+\.\d{4,}\s*", "", s).strip()
        s = re.sub("\s+", " ", s)
        s1, s2 = s.split('->')
        return \
            '$' + ''.join(s1.strip().split()[1:]) + '$' + \
            r' $\rightarrow$ ' + \
            '$' + ''.join(s2.strip().split()[1:]) + '$'
        #return s



    def format_info(self, s = '\n', brief=False):
        if brief:
            return s.join( \
                [r'Name: {0:s}'.format(self.molname if self.obs_data == None else \
                                       self.obs_data['molecule']) + \
                 ' '*4 +
                 r'Qnum: ' + (self.qnum_reformat() if self.obs_data == None else \
                              self.obs_data['transition']),
                 #
                 r'Lam: {0:.3f} '.format(self.lam0/1e4) + '$\mu$m' + ' '*4 +
                 r'F: {0:.6e} Hz'.format(self.f0),
                 #
                 r'$E_{\rm up}$: ' + '{0:.1f} K'.format(self.E_up) + ' '*4 +
                 r'$E_{\rm low}$: ' + '{0:.1f} K'.format(self.E_low) +
                 ' '*4 +
                 r'$A_{\rm ul}$: ' + '{0:.2e}'.format(self.Aul) + ' s$^{{-1}}$' + ' '*4,
                ])
        else:
            return s.join( \
                [r'Name: {0:s}'.format(self.molname if self.obs_data == None else \
                                       self.obs_data['molecule']) + \
                 ' '*4 +
                 r'Qnum: ' + (self.qnum_reformat() if self.obs_data == None else \
                              self.obs_data['transition']),
                 #
                 r'Lam: {0:.3f} '.format(self.lam0/1e4) + '$\mu$m' + ' '*4 +
                 r'F: {0:.6e} Hz'.format(self.f0),
                 #
                 r'$E_{\rm up}$: ' + '{0:.1f} K'.format(self.E_up) + ' '*4 +
                 r'$E_{\rm low}$: ' + '{0:.1f} K'.format(self.E_low) +
                 ' '*4 +
                 r'$A_{\rm ul}$: ' + '{0:.2e}'.format(self.Aul) + ' s$^{{-1}}$' + ' '*4 +
                    r'$B_{\rm ul}$: ' + '{0:.2e}'.format(self.Bul) + ' '*4 +
                    r'$B_{\rm lu}$: ' + '{0:.2e}'.format(self.Blu),
                 #
                 r'Dist: ' + '{0:.1f} pc'.format(self.dist) + ' '*4 +
                 r'Incl: ' + '{0:.1f} deg'.format(self.theta),
                 #
                 r'IntLineFlux: {0:.2e}'.format(self.intfluxl/1e-18) + \
                 ' ($10^{-18}$ W m$^{{-2}}$)' +
                 ('' if (self.obs_data == None or self.obs_data == []) else \
                 ' '*4 +
                 r'ObsLineFlux: {0:s}'.format(self.obs_data['obs'].ss) + \
                 ' ($10^{-18}$ W m$^{{-2}}$)')
                ])


    def saveas_pdf(self, pdfname, figsize=(5,2), pos=(0.1,0.1,0.85,0.75),
                   unit='lam', vr=None, brief=False):
        #
        import numpy as np
        if self.__fig__ == None:
            self.__fig__ = plt.figure(figsize = figsize)
        else:
            self.__fig__.set_size_inches(figsize[0], figsize[1])
        #
        ax = self.__fig__.add_axes(pos)
        #
        if unit == 'lam':
            ax.plot(self.lam, self.spec)
            ax.set_xlabel(r'$\lambda$ ($\mu$m)')
            ax.set_ylabel(r'Flux (jy)')
        else:
            ax.plot(self.v/1e3, self.spec, drawstyle='steps-mid')
            ax.set_xlabel(r'$V$ (km s$^{-1}$)')
            ax.set_ylabel(r'Flux (jy)')
        if vr != None:
            ax.set_xlim(vr)
        mx,mn = np.nanmax(self.spec), np.nanmin(self.spec)
        if (mx-mn) <= mn:
            ylim = (mn*0.9, mx*1.1)
        else:
            ylim = (0, mx*1.1)
        ax.set_ylim(ylim)
        #
        plt.figtext(0.1, 0.87, self.format_info(brief=brief),
                    fontdict={'horizontalalignment': 'left', 'verticalalignment': 'bottom'})
        plt.figtext(0.1, 0.60, r'$\int F_\nu d\nu$ = {0:.1e} W m$^{{-2}}$'.format(self.intfluxl),
                    fontdict={'horizontalalignment': 'left', 'verticalalignment': 'bottom'})
        plt.savefig(pdfname, bbox_inches='tight')
        self.__fig__.clf()
        return


    def saveas_pdf_swift(self, unit='lam', vr=None, brief=False):
        fname_pdf = self.fname[:-4] + 'pdf'
        print('Saving pdf: ', fname_pdf)
        self.saveas_pdf(fname_pdf, unit=unit, vr=vr, brief=brief)
        return


class batch_proc:
    def __init__(self, fnames, _obs_data=None, noObs=False):
        self.specs = [specline(f) for f in fnames]
        if _obs_data != None:
            for spec in self.specs:
                self.match_obs(spec, _obs_data)
                self.wash_matched_obs(spec)
        for spec in self.specs:
            if spec.obs_data == None:
                print('No observational data:')
                print('\t', spec.molname, spec.fname)
                if noObs:
                    return
                else:
                    print('This line will not be included.')
            #else:
            #    print spec.lam0, spec.obs_data['lam']
        self.specs = filter(lambda s: s.obs_data != None, self.specs)
        def spcmp(x, y):
            cx = (x.obs_data['molecule'], x.obs_data['transition'])
            cy = (y.obs_data['molecule'], y.obs_data['transition'])
            if cx < cy:
                return -1
            elif cx > cy:
                return 1
            else:
                return 0
        self.specs.sort(cmp = spcmp)
        return
    
    def match_obs(self, spec, _obs_data):
        spec.obs_data = \
            [od for od in _obs_data \
                if (abs(od['lam'] - spec.lam0*1e-4) < 0.02 * (od['lam'] + spec.lam0*1e-4) and \
                    abs(od['Eup'] - spec.E_up)      < 0.02 * (od['Eup'] + spec.E_up))
            ]
        return

    def wash_matched_obs(self, spec):
        from math import log
        if spec.obs_data == None:
            return
        if   len(spec.obs_data) == 0:
            spec.obs_data = None
        elif len(spec.obs_data) == 1:
            spec.obs_data = spec.obs_data[0]
        else:
            print(len(spec.obs_data), spec.molname, spec.lam0/1e4, spec.fname)
            tmp = 1e99
            idx = 0
            for i in xrange(len(spec.obs_data)):
                print('\t', spec.obs_data[i]['molecule'],
                    spec.obs_data[i]['transition'], spec.obs_data[i]['lam'],
                    spec.lam0, spec.intfluxl, spec.obs_data[i]['obs'].val)
                if spec.obs_data[i]['molecule'] == spec.molname_disp and \
                   spec.molname_disp != spec.molname:
                    idx = i
                    break
                tmp1 = abs(log(spec.intfluxl / spec.obs_data[i]['obs'].val))
                if tmp > tmp1 and spec.obs_data[i]['molecule'] == spec.molname_disp:
                    idx = i
                    tmp = tmp1
            print('\tChoose:')
            print('\t', spec.obs_data[idx]['molecule'],
                spec.obs_data[idx]['transition'], spec.obs_data[idx]['lam'],
                spec.lam0, spec.intfluxl, spec.obs_data[idx]['obs'].val)
            spec.obs_data = spec.obs_data[idx]
        return

    def all_to_pdf(self, unit='lam', vr=None, brief=False):
        for s in self.specs:
            s.saveas_pdf_swift(unit=unit, vr=vr, brief=brief)
        return

    def all_to_html(self, fname):
        with open(fname, 'w') as f:
            write_html_header(f)
            col_names = ['Molecule', 'Line', '&lambda; (&mu;m)',
                         'E<sub>up</sub> (K)', 'Model flux (10<sup>-18</sup> W m<sup>-2</sup>)',
                         'Obs flux (10<sup>-18</sup> W m<sup>-2</sup>)', 'M/O']
            write_table_header(f, col_names=col_names)
            for s in self.specs:
                tab_items = [s.obs_data['molecule'],
                             s.obs_data['transition'],
                             s.lam0/1e4,
                             s.E_up, s.intfluxl/1e-18, s.obs_data['obs'].ss,
                             s.intfluxl/s.obs_data['obs'].val]
                write_table_row(f, tab_items)
            write_table_tail(f)
            write_html_tail(f)

    def all_to_html_cmp(self, mcmp, fname):
        with open(fname, 'w') as f:
            write_html_header(f)
            col_names = ['Molecule', 'Line', '&lambda; (&mu;m)',
                         'E<sub>up</sub> (K)',
                         'Model flux (10<sup>-18</sup> W m<sup>-2</sup>)',
                         'CMP flux (10<sup>-18</sup> W m<sup>-2</sup>)',
                         'Obs flux (10<sup>-18</sup> W m<sup>-2</sup>)',
                         'M/O',
                         'C/O']
            write_table_header(f, col_names=col_names)
            for i in xrange(len(self.specs)):
                s = self.specs[i]
                scmp = mcmp.specs[i]
                tab_items = [s.obs_data['molecule'],
                             s.obs_data['transition'],
                             s.lam0/1e4,
                             s.E_up,
                             s.intfluxl/1e-18,
                             scmp.intfluxl/1e-18,
                             s.obs_data['obs'].ss,
                             s.intfluxl/s.obs_data['obs'].val,
                             scmp.intfluxl/s.obs_data['obs'].val]
                write_table_row(f, tab_items)
            write_table_tail(f)
            write_html_tail(f)

    def draw_overview(self, fname, spec_cmp):
        cpr = [{'name': spec.obs_data['molecule'] + ' ' + spec.obs_data['transition'],
                'observed': spec.obs_data['obs'].val,
                'full': scmp.intfluxl,
                'depleted': spec.intfluxl,
                'is_upperlim': spec.obs_data['obs'].ul != None} \
               for spec, scmp in zip(self.specs, spec_cmp)]
        draw_cpr(cpr, fname)



def draw_cpr(cpr, fname, labelfontsize=12, k0='observed', k1='full', k2='depleted'):
    from draw.long_function_definitions import set_axis_format
    from matplotlib.lines import Line2D
    nitem = len(cpr)
    #
    mx, mn = -1e99, 1e99
    for i in xrange(len(cpr)):
        item = cpr[i]
        f2o = item[k1]/item[k0]
        d2o = item[k2]/item[k0]
        mx = max(mx, max(f2o, d2o))
        mn = min(mn, min(f2o, d2o))
    #
    mn = mn / 1.4
    mx = mx * 1.4
    fig = plt.figure(figsize=(8, 4.5))
    ax = fig.add_axes([0.1,0.4,0.85,0.55])
    ax.set_xlim((0,1))
    ax.set_xscale('linear')
    ax.set_xticks([])
    ax.set_xticklabels([])
    ax.set_ylim([mn, mx])
    ax.set_yscale('log')
    ax.set_ylabel('Ratios')
    ax.axhline(y=1, color='black', linestyle='--')
    set_axis_format(ax, labelfontsize=labelfontsize, tickfontsize=labelfontsize, onlyfont=True)

    dx = 1.0 / (nitem + 1)
    ls1 = {'marker': 'o', 'linestyle': 'none', 'markersize': 10,
           'markerfacecolor': 'red', 'markeredgecolor': 'none'}
    ls2 = {'marker': 's', 'linestyle': 'none', 'markersize': 8,
           'markerfacecolor': 'blue', 'markeredgecolor': 'none', 'alpha': 0.7}
    
    for i in xrange(len(cpr)):
        item = cpr[i]
        f2o = item[k1]/item[k0]
        d2o = item[k2]/item[k0]
        x = dx * (i+1)
        if item['is_upperlim']:
            ax.arrow(x, 1, 0, -0.7, 
                     width=0.008, head_width=0.02, head_length=0.15, overhang=0.15,
                     facecolor='red', edgecolor='none', alpha=0.2)
        ax.text(x, mn/1.1, item['name'], rotation='vertical', fontsize=labelfontsize,
                horizontalalignment='center', verticalalignment='top')
        ax.plot([x], [f2o], **ls1)
        ax.plot([x], [d2o], **ls2)
    p = [Line2D([],[], **ls1), Line2D([],[], **ls2)]
    L = [k1, k2]
    lgd = ax.legend(p, L, bbox_to_anchor=(0.1, 1.0, 0.5, 0.15), frameon=False,
                    ncol=2, numpoints=1, fontsize=labelfontsize)
    plt.savefig(fname, bbox_inches='tight')
    return




def write_html_header(f):
    f.write('''
    <!DOCTYPE html>
    <html>
    <title>Integrated line intensities</title>
    <body>
    ''')
    return


def write_html_tail(f):
    f.write('''
    </body>
    </html>
    ''')
    return


def write_table_header(f, col_names = None):
    f.write('''
    <table>
    ''')
    if col_names != None:
        for c in col_names:
            f.write('<td>' + c + '</td>')
    return


def write_table_tail(f):
    f.write('''
    </table>
    ''')
    return


def sub_html_symbol(s):
    import re
    s = re.sub("&", "&amp", s).strip()
    s = re.sub("<-+", "&larr;", s)
    s = re.sub("-+>", "&rarr;", s)
    s = re.sub("<", "&lt", s)
    s = re.sub(">", "&gt", s)
    return s


def write_table_row(f, items):
    f.write('<tr>')
    for it in items:
        if type(it) == str:
            f.write('<td>' + sub_html_symbol(it) + '</td>')
        else:
            f.write('<td>' + it.__str__() + '</td>')
    f.write('</tr>')
    return


class obs_value():
    def __init__(self, val=None, pe=None, ne=None, ul=None, ll=None):
        self.val = val
        self.pe  = pe
        self.ne  = ne
        self.ul  = ul
        self.ll  = ll


def get_obs_value_from_latex_str(t, scale=1e-18):
    v = obs_value()
    s_lt = '$<$'
    s_gt = '$>$'
    s_pm = r'$\pm$'
    v.ss = t
    if   t.startswith(s_lt):
        try:
            v.ul = float(t[len(s_lt):]) * scale
        except:
            raise ValueError('Fail to convert ' + t[len(s_lt):] + ' to float.')
        v.val = v.ul
    elif t.startswith(s_gt):
        try:
            v.ll = float(t[len(s_gt):]) * scale
        except:
            raise ValueError('Fail to convert ' + t[len(s_gt):] + ' to float.')
        v.val = v.ll
    elif s_pm in t:
        _val, _pe = t.split(s_pm)
        try:
            v.val = float(_val) * scale
            v.pe  = float(_pe)  * scale
            v.ne  = float(_pe)  * scale
        except:
            raise ValueError('Fail to convert ' + t + ' to float.')
    else:
        try:
            v.val = float(t) * scale
        except:
            raise ValueError('Fail to convert ' + t + ' to float.')
    return v


def parse_latex_txt_table(fname):
    from load_latex_tab import load_txt_tab
    lines = [line.split('&') for line in load_txt_tab(fname)]
    d = [{'molecule': line[0].strip(),
          'transition': line[1].strip(),
          'lam': float(line[2].strip()),
          'Eup': float(line[3].strip()),
          'obs': get_obs_value_from_latex_str(line[4].strip())} for line in lines]
    return d


def glob_fnames(parent_dir, sub_dirs, filetype = 'images/*.fits'):
    from glob import glob
    from os.path import join as opj
    fs = []
    for subdir in sub_dirs:
        fs = fs + glob(opj(parent_dir, subdir, filetype))
    fs.sort()
    return fs


def fits_job1(model_names=None, model_cmps=None):
    from os.path import join as opj
    from itertools import product
    #
    obs_datfile = '/n/Users/fdu/now/tab_obs.dat'
    res_dir = '/n/Users/fdu/now/res/'
    out_dir = '/n/Users/fdu/'
    dir_lines = ['CO', '13CO', 'C18O', 'O', 'C', 'C+', 'oH2O_A', 'pH2O_A', 'HD']
    #
    for model_name, model_cmp in product(model_names, model_cmps):
        model_dir = opj(res_dir, model_name)
        fname_html = opj(out_dir, 'tab_{0:s}.html'.format(model_name.replace('/', '')))
        fname_overview = opj(out_dir, 'cmp_{0:s}.pdf'.format(model_name.replace('/', '')))
        #
        obsd = parse_latex_txt_table(obs_datfile)
        #
        fnames = glob_fnames(model_dir, dir_lines)
        b = batch_proc(fnames, _obs_data = obsd)
        #
        cmp_dir = opj(res_dir, model_cmp)
        fnames_cmp = glob_fnames(cmp_dir, dir_lines)
        bcmp = batch_proc(fnames_cmp, _obs_data = obsd)
        #
        b.all_to_html_cmp(bcmp, fname_html)
        print('Html file generated: ', fname_html)
        #
        b.draw_overview(fname_overview, bcmp.specs)
        print('Pdf file generated: ', fname_overview)


def fits_job2(res_dir = '/n/Users/fdu/now/res/',
              dir_lines = None,
              model_names=None):
    from os.path import join as opj
    from itertools import product
    #
    #
    for model_name in model_names:
        model_dir = opj(res_dir, model_name)
        #
        fnames = glob_fnames(model_dir, dir_lines)
        b = batch_proc(fnames, noObs=True)
        b.all_to_pdf(unit='vel', vr=(-10,10), brief=True)


def fits_job3(model_path = None,
              pdfname = None,
              H2O_pattern = '[op]H2O_B*',
              OH_pattern = 'OH_A*',
              obsspec_dir = '/n/Users/fdu/now/inp/twhya_resi_lh_from_KZhang.dat'):
    #
    import glob
    import os
    import numpy as np
    from draw.long_function_definitions import set_axis_format, plot_spectra_dir, stepify
    from matplotlib.lines import Line2D
    #
    tw_kzh = np.loadtxt(obsspec_dir, skiprows=4)
    
    # Line window 1
    
    sdata = tw_kzh
    figsize = (8,6)
    xlim = (27.95,35)
    ylim = (-0.1, 0.5)
    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes((0.1,0.1,0.85,0.4),
                      xlabel=r'$\lambda$ ($\mu$m)',
                      ylabel='Flux (Jy)',
                      autoscalex_on=False, autoscaley_on=False,
                      xscale='linear', yscale='linear',
                      xlim=xlim, ylim=ylim)
    
    ax.fill_between(*stepify(sdata[:,0], sdata[:,2]), edgecolor='gray', facecolor='gray', lw=1)
    
    pt = os.path.join(model_path, H2O_pattern, 'images/*.fits')
    files = glob.glob(pt)
    #print files
    
    if len(files) > 0:
        plot_spectra_dir(files, set_ax=False,
                     yRange=ylim, printout=False, hide_comp=True, lamRange=xlim,
                     ax=ax, fig=fig, color='red', lw=0.5, color_coadd='red', lw_coadd=1.0,
                     removebase=True, resolvingpower=600, coadd=True)
    else:
        print('No file found with pattern', pt)

    pt = os.path.join(model_path, OH_pattern, 'images/*.fits')
    files = glob.glob(pt)
    
    if len(files) > 0:
        plot_spectra_dir(files, set_ax=False,
                     yRange=ylim, printout=False, hide_comp=True, lamRange=xlim,
                     ax=ax, fig=fig, color='blue', lw=0.5, color_coadd='blue', lw_coadd=1.0,
                     removebase=True, resolvingpower=600, coadd=True)
    else:
        print('No file found with pattern', pt)

    set_axis_format(ax, graygrid=True, majorgridon=False, minorgridon=False)

    # Line window 2
    
    sdata = tw_kzh
    xlim = (21,28.05)
    ylim = (-0.1, 0.5)
    ax = fig.add_axes((0.1,0.57,0.85,0.4),
                      xlabel='',
                      ylabel='Flux (Jy)',
                      autoscalex_on=False, autoscaley_on=False,
                      xscale='linear', yscale='linear',
                      xlim=xlim, ylim=ylim)
    
    ax.fill_between(*stepify(sdata[:,0], sdata[:,2]), edgecolor='gray', facecolor='gray', lw=1)
    
    files = glob.glob(os.path.join(model_path, H2O_pattern, 'images/*.fits'))
    
    if len(files) > 0:
        plot_spectra_dir(files, set_ax=False,
                     lamRange=xlim, yRange=ylim, printout=False, hide_comp=True,
                     ax=ax, fig=fig, color='red', lw=0.5, color_coadd='red', lw_coadd=1.0,
                     removebase=True, resolvingpower=600, coadd=True)
    
    files = glob.glob(os.path.join(model_path, OH_pattern, 'images/*.fits'))
    
    if len(files) > 0:
        plot_spectra_dir(files, set_ax=False,
                     lamRange=xlim, yRange=ylim, printout=False, hide_comp=True,
                     ax=ax, fig=fig, color='blue', lw=0.5, color_coadd='blue', lw_coadd=1.0,
                     removebase=True, resolvingpower=600, coadd=True)
    
    set_axis_format(ax, graygrid=True, majorgridon=False, minorgridon=False)

    _ = plt.legend([Line2D([], [], color='red'),
                    Line2D([], [], color='blue')],
                   ['H$_2$O', 'OH'],
                   ncol=2, frameon=False,
                   fontsize=20, framealpha=0.3, loc=(0.1, 0.75))
    
    plt.savefig(pdfname, bbox_inches='tight')


if __name__ == '__main__':
    fits_job2(\
              res_dir = '/n/Fdu1/work/grid_20150703_redo/results',
              model_names = ['run_133', 'run_134'],
              dir_lines = ['H2O_ground', 'H2O_1'])
    #fits_job1( \
    #    model_names = ['20150307_dep_s17b'],
    #    model_cmps  = ['20150307_noD_d']
    #)
  # fits_job2(\
  #           #res_dir = '/u/Moria2/fdu/now/res',
  #           res_dir = '/n/Users/fdu/now/res',
  #           model_names = ['20160714_TWH_a'],
  #           dir_lines = ['C2H_2_NLTE', 'C2H_2_NLTE_vwidth'])
    #fits_job2(model_names = ['20150609_DMT_C2H_1'], dir_lines = ['HD'])
    #fits_job2(model_names = ['20150307_dep_s_NH3'], dir_lines = ['oNH3', 'pNH3'])
    #fits_job2(model_names = ['20150307_noD_d', '20150307_dep_s', '20150307_dep_s_NH3'], dir_lines = ['oNH3'])
    #fits_job3(model_path = '/n/Users/fdu/now/res/20150307_dep_s17b/', 
    #          pdfname = '/n/Users/fdu/now/H2O_OH_Spitzer_s17b.pdf')
    #fits_job3(model_path = '/n/Users/fdu/now/res/20140929_a1_K/', 
    #          H2O_pattern = 'H2O_[op]_iter1',
    #          OH_pattern = 'OH_IRS_hitran_iter1_NLTE',
    #          pdfname = '/n/Users/fdu/now/H2O_OH_Spitzer_K.pdf')
