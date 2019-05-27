import numpy as np
from astropy.io import fits
from astropy.io.fits import getval
import glob
import os
from load_latex_tab import load_latex_tab
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as pl
import hashlib

def name2num(name):
    return int(hashlib.md5(name).hexdigest()[-6:], 16)%1000/1000.0

#def get_val_from_key(h, key, ExtName = 'FluxSpec'):
#    for ext in h:
#        hd = ext.header
#        if hd['EXTNAME'] == ExtName and key in hd.keys():
#            return hd[key]
#    return None


def lutstr_2_dic(s):
    s = s.replace(',', ' ')
    s = s.replace(':', ' ')
    s = s.replace('=', ' ')
    a = s.split()
    d = {'name': a[0]}
    for i in xrange(1,len(a),2):
        try:
            d[a[i]] = float(a[i+1])
        except:
            d[a[i]] = a[i+1]
    return d


def flux2Kelvin(flx, nu, dV, D=3.5, alpha=1.0):
    """
    flx: flux in W m-2
    nu: frequency in Hz
    dV: velocity width in km s-1
    D: telescope diameter in m
    alpha: dtheta = lambda/D * alpha
    """
    from consts.physical_constants import phy_kBoltzmann_SI, phy_Pi, phy_SpeedOfLight_SI
    return D*D * flx / \
           (2.0 * phy_kBoltzmann_SI * phy_Pi * alpha**2 * nu * dV * 1e3 / phy_SpeedOfLight_SI)

def de_latexmath(s):
    a = s.replace('$', '')
    a = a.replace('{', '')
    return a.replace('}', '')


def load_lut(res_dir, RT_dir, images_dir, fRT_pattern, LUT_file, dV=1.0, key_flux = 'INTFLUXL', toKelvin=True, nrows=None):
    with open(LUT_file, 'r') as f:
        LUT = f.readlines()
    LUT.pop(0)

    if nrows != None:
        LUT = LUT[0:nrows]

    LUT[:] = [lutstr_2_dic(lut) for lut in LUT]

    for lut in LUT:
        fpattern = os.path.join(res_dir, 'run_' + lut['name'], RT_dir, images_dir, fRT_pattern)
        fname = glob.glob(fpattern)
        if len(fname) > 1:
            raise Exception('Multiple fname match the pattern:', fpattern)
        if len(fname) == 0:
            raise Exception('No file matches the pattern:', fpattern)
        fname = fname[0]
        nu = getval(fname, 'F0', 'FluxSpec')
        if toKelvin:
            lut[key_flux] = flux2Kelvin(getval(fname, key_flux, 'FluxSpec'), nu, dV)
        else:
            lut[key_flux] = getval(fname, key_flux, 'FluxSpec')

        lut['nu'] = nu
        lut['lam'] = getval(fname, 'LAM0', 'FluxSpec')
        lut['qnum'] = getval(fname, 'QNUM', 'FluxSpec')
        lut['Eup'] = getval(fname, 'EUP', 'FluxSpec')
        lut['molname'] = getval(fname, 'MOL-DB', 'FluxSpec')
    return LUT


def spectralType2num(spt):
    #spectral_to_temperature = \
    #    { 'O': 41.00e3,
    #      'B': 31.00e3,
    #      'A':  9.50e3,
    #      'F':  7.24e3,
    #      'G':  5.92e3,
    #      'K':  5.30e3,
    #      'M':  3.85e3,
    #    }
    from batchproc.functions import spectral_to_temperature
    return spectral_to_temperature[spt]

def f_cmp1(x, y):
    if x['rin'] < y['rin']:
        return -1
    elif x['rin'] == y['rin'] and x['mdust'] < y['mdust']:
        return -1
    elif x['rin'] == y['rin'] and x['mdust'] == y['mdust'] and \
         spectralType2num(x['spectralType']) < spectralType2num(y['spectralType']):
        return -1
    elif x['rin'] == y['rin'] and x['mdust'] == y['mdust'] and \
         x['spectralType'] == y['spectralType']:
        return 0
    else:
        return 1

def f_cmp2(x, y):
    if spectralType2num(x[0]) < spectralType2num(y[0]):
        return -1
    elif spectralType2num(x[0]) == spectralType2num(y[0]) and x[1] == y[1]:
        return 0
    else:
        return 1

def get_len_1(L, key):
    for i in xrange(1,len(L)):
        if L[i][key] == L[0][key]:
            return i
    return None

def get_len_2(L, len1, key):
    for i in xrange(1,len(L)/len1):
        if L[i*len1][key] == L[0][key]:
            return i
    return None



def make_sub_axis(fig, npanx, npany, kx, ky,
        pansepxfrac=0.2, pansepyfrac=0.06,
        xmarginleft = 0.1, ymarginlower = 0.1,
        labelfontsize=25, tickfontsize=20,
        xtitle = '', ytitle = '',
        xscale = 'linear', yscale = 'linear',
        xticks = None, yticks = None,
        xticklabels = None, yticklabels = None,
        xRange = None, yRange = None,
        less_xlabel=False, less_ylabel=False,
        ):
    #
    import long_function_definitions as lfd
    #
    panwx = (0.99 - xmarginleft)  / npanx
    panwy = (0.99 - ymarginlower) / npany
    pansepx = panwx * pansepxfrac
    pansepy = panwy * pansepyfrac
    panwx -= pansepx
    panwy -= pansepy
    #
    xleft  = xmarginleft  + panwx * kx + pansepx * kx
    ylower = ymarginlower + panwy * ky + pansepy * ky
    pos = (xleft, ylower, panwx, panwy)

    ax = fig.add_axes(pos,
        xlabel=xtitle, ylabel=ytitle,
        autoscalex_on=False, autoscaley_on=False,
        xscale=xscale, yscale=yscale,
        xlim=xRange, ylim=yRange)

    if xticks != None:
        ax.set_xticks(xticks)
        if xticklabels != None:
            ax.set_xticklabels(xticklabels)
    if yticks != None:
        ax.set_yticks(yticks)
        if yticklabels != None:
            ax.set_yticklabels(yticklabels)
    if kx > 0:
        ax.set_ylabel('')
        ax.set_yticklabels([])
    if ky > 0:
        ax.set_xlabel('')
        ax.set_xticklabels([])
    if less_xlabel and kx > 0:
        ax.set_xlabel('')
    if less_ylabel and ky > 0:
        ax.set_ylabel('')
    lfd.set_axis_format(ax, labelfontsize=labelfontsize, tickfontsize=tickfontsize, onlyfont=True)
    #ax.set_axisbelow(False)
    return ax


def split_LUT(LUT, keys):
    '''
    keys: ['k1', 'k2', ...]
    Return: 
    '''
    vks = []
    for lut in LUT:
        vks.append((lut[k] for k in keys))
    vks = list(set(vks))


def get_distinct_vals(LUT, key):
    return np.array(list(set([lut[key] for lut in LUT])))


def select_LUT(LUT, kv):
    LUT_new = []
    for lut in LUT:
        match = True
        for k in kv.keys():
            if lut[k] != kv[k]:
                match = False
                break
        if match:
            d = {}
            for k in lut.keys():
                if k not in kv.keys():
                    d.update({k: lut[k]})
            LUT_new.append(d)
    return LUT_new



def extract_xy(LUT, kx, ky):
    pairs = [(lut[kx], lut[ky]) for lut in LUT]
    return pairs



def make_data(LUT, key_flux = 'INTFLUXL'):
    coords_1 = list(set([(lut['rin'], lut['mdust']) for lut in LUT]))
    coords_1.sort()
    dict_1 = {c: [] for c in coords_1}
    minflux, maxflux = 1e99, -1e99
    for lut in LUT:
        dict_1[(lut['rin'], lut['mdust'])].append((lut['spectralType'], lut[key_flux]))
        minflux = min(minflux, lut[key_flux])
        maxflux = max(maxflux, lut[key_flux])
    for c in coords_1:
        dict_1[c].sort(cmp=f_cmp2)
        #print c
        #print zip(*dict_1[c])
    return coords_1, dict_1, minflux, maxflux


def mass2spectral_idx(m):
    from batchproc.functions import spectral_to_mass
    a = sorted(spectral_to_mass[k] for k in ['M', 'K', 'G', 'F', 'A', 'B'])
    if m <= a[0]:
        return 0
    for i in xrange(len(a)-1):
        if a[i] <= m and m < a[i+1]:
            return i + (m - a[i]) / (a[i+1] - a[i])
    return len(a) - 1


def getAbbrevSourceName(srcName):
    parts = srcName.split('~')
    try:
        first = parts[0][0]
        second = parts[1][:2]
    except:
        print parts
        raise
    return first + second


def collect_data(dobs, sources, trans, co_this, mdust_s):
    g2d = 1e2
    dist_model = 1e2
    radius_model = 4e2
    dV_assume = 1.0 # km s-1
    d_collect = []
    for line in dobs:
        if line[1] == trans:
            sourcename = line[0]
            #rms        = float(line[2])
            rms_integrated = float(line[2])
            if len(line[3].strip()) == 0:
                val_integrated = 0.0
            else:
                val_integrated = float(line[3])
            #dV_channel = float(line[4])
            Mstar = float(sources[sourcename][r'$M_\text{star}$'])
            Mdisk = float(sources[sourcename][r'$M_\text{disk}$'])
            dist  = float(sources[sourcename][r'$d$'])
            r_out  = float(sources[sourcename][r'$R_\text{out}$'])
            Mdust = Mdisk / g2d
            #rms_integrated = rms * np.sqrt(dV_channel * dV_assume)
            coeff_convert = (radius_model/r_out)**2
            rms_cmpwithmodel = 1e-3 * rms_integrated * (dist/dist_model)**2 * coeff_convert
            val_cmpwithmodel = 1e-3 * val_integrated * (dist/dist_model)**2 * coeff_convert
            spec_idx = mass2spectral_idx(Mstar)
            in_this_plot = np.min(abs(Mdust - mdust_s)) == abs(Mdust - co_this[1])
            if in_this_plot:
                d_collect.append((spec_idx, rms_cmpwithmodel, val_cmpwithmodel, getAbbrevSourceName(sourcename)))
    return d_collect



def make_plot(coords_1, LUT, dict_1, minflux, maxflux, fig_fname, figsize=(16,16),
              dobs = None, sources = None, trans = None, mdust_s = None):
    l_x = get_len_1(coords_1, 1)
    l_y = len(coords_1) / l_x
    print l_x, l_y
    
    lut0 = LUT[0]
    
    textsize = 15
    labelfontsize = textsize*1.2
    
    fig = pl.figure(figsize = figsize)
    pl.figtext(0.1, 1,
            'Molecule: ' + lut0['molname'] + \
                '\t\tTransition: ' + lut0['qnum'].replace('->', r'$\rightarrow$') + '\n' + 
            r'$\nu$: ' + '{0:.6e}'.format(lut0['nu']) + ' Hz\t\t' + \
                r'$\lambda$: ' + '{0:.6e}'.format(lut0['lam']/1e4) + ' $\mu$m\t\t' + \
            r'$E_{{\rm up}}$: ' + '{0:.2f}'.format(lut0['Eup']) + ' K',
        fontsize = textsize,
        fontdict={'horizontalalignment': 'left', 'verticalalignment': 'bottom'})
    
    maxy, miny = maxflux, minflux
    miny = max(maxy/1e4, miny)
    
    if not (dobs is None or sources is None or trans is None or mdust_s is None):
        for co in coords_1:
            d_collect = collect_data(dobs, sources, trans, co, mdust_s)
            if len(d_collect) > 0:
                d_rms_s = [_[1] for _ in d_collect]
                miny = min(miny, min(d_rms_s)/2.0)
                maxy = max(maxy, max(d_rms_s))

    iax = 1
    ax_s = []
    for i in xrange(l_y):
        for j in xrange(l_x):
            ax = make_sub_axis(fig, l_x, l_y, j, i,
                    pansepxfrac=0.04, pansepyfrac=0.0,
                    labelfontsize=labelfontsize, tickfontsize=textsize,
                    xtitle = 'Stellar type', ytitle = r'$\int T_{} dV$ (K km s$^{{-1}}$)',
                    xticks = [0,1,2,3,4,5],
                    xticklabels = ['M', 'K', 'G', 'F', 'A', 'B'],
                    yscale = 'log' if maxy/miny>30 else 'linear',
                    xRange = (-0.5, 5.5),
                    yRange = (miny, maxy),
                    less_xlabel=True, less_ylabel=True,
                 )
            co = coords_1[iax-1]
            x, y = zip(*dict_1[co])
            ax.plot(y, marker='o')
            ax.set_title(r'($r_{{\rm in}}$={0:.1e}, $m_{{\rm dust}}$={1:.1e})'.format(co[0], co[1]),
                         fontsize=textsize,
                         fontdict={'verticalalignment': 'top', 'y':0.95})
            ax_s.append(ax)
            iax += 1
            #
            if not (dobs is None or sources is None or trans is None or mdust_s is None):
                d_collect = collect_data(dobs, sources, trans, co, mdust_s)
                if len(d_collect) > 0:
                    for d_c in d_collect:
                        ax.plot([d_c[0]-0.1, d_c[0]+0.1], [d_c[1], d_c[1]], color='red')
                        ax.arrow(d_c[0], d_c[1], 0.0, -d_c[1]*0.25,
                                 width = 0.01, overhang=0.5,
                                 head_length = d_c[1]*0.2, length_includes_head=False,
                                 edgecolor = 'red', facecolor = 'red')
                        #
                        ax.plot([d_c[0]-0.1, d_c[0]+0.1], [d_c[1]*16, d_c[1]*16], color='magenta')
                        ax.arrow(d_c[0], d_c[1]*16, 0.0, -d_c[1]*16*0.5,
                                 width = 0.01, overhang=0.5,
                                 head_length = d_c[1]*16*0.2, length_includes_head=False,
                                 edgecolor = 'magenta', facecolor = 'magenta')
    
    pl.savefig(fig_fname, bbox_inches='tight')


def task2(remove_d2gs=None, remove_mdusts=None, figsize = (16,16)):
    grid_dir = '/n/Fdu1/work/grid_20150703_redo/'
    res_dir = os.path.join(grid_dir, 'results')
    fig_dir = os.path.join(grid_dir, 'figures')
    images_dir = 'images'
    LUT_file = os.path.join(grid_dir, 'config_files/LUT.dat')
    dV = 1.0
    colors = ['black', 'red', 'green', 'blue', 'magenta', 'Brown']
    arrow_color = (0.1,0.1,0.1)
    errorbar_color = (0.6,0.6,0.6)
    textsize = 15
    labelfontsize = textsize*1.2

    groundstate_dir = 'H2O_ground_new'
    upperstate_dir = 'H2O_1_new'
    fname_postfix = '_20170428.pdf'

    fPatterns = {'$1_{10}-1_{01}$': [groundstate_dir, 'line_00001_00212_00001_5.56936E+11_0.00.fits'],
                 '$1_{11}-0_{00}$': [groundstate_dir, 'line_00002_00437_00001_1.11334E+12_0.00.fits'],
                 '$3_{12}-2_{21}$': [upperstate_dir,  'line_00002_00456_00001_1.15313E+12_0.00.fits'],
                 '$3_{12}-3_{03}$': [upperstate_dir,  'line_00001_00428_00001_1.09736E+12_0.00.fits']}
    
    detections_aux = {'$1_{10}-1_{01}$': [{'val': 25.2, 'sigma': 1.1, 'Source': 'TWH'},
                                          {'val': 15.5, 'sigma': 2.5, 'Source': 'DMT'}],
                      '$1_{11}-0_{00}$': [{'val': 44.6, 'sigma': 2.8, 'Source': 'TWH'}],
                     }

    if not os.path.exists(fig_dir):
        os.mkdir(fig_dir)

    #tab_file = 'tab_for_python_Aug22.tab'
    tab_file = 'tab_for_python_Apr28.tab'
    #dobs, sources, transitions = load_latex_tab(os.path.join(grid_dir, 'tab_for_python_new.tab'))
    dobs, sources, transitions = load_latex_tab(os.path.join(grid_dir, tab_file))
    
    maxy, miny = -1e99, 1e99

    LUT_s = {}
    for trans in fPatterns.keys():
        LUT = load_lut(res_dir, fPatterns[trans][0], images_dir, fPatterns[trans][1], LUT_file,
                       dV = dV, key_flux = 'INTFLUXL')
        LUT_s.update({trans: LUT})

        rins     = get_distinct_vals(LUT, 'rin')
        routs    = get_distinct_vals(LUT, 'rout')
        mdusts   = get_distinct_vals(LUT, 'mdust')
        d2gs     = get_distinct_vals(LUT, 'd2g')
        f_o_deps = get_distinct_vals(LUT, 'f_o_dep')
        f_c_deps = get_distinct_vals(LUT, 'f_c_dep')
        #
        if remove_d2gs:
            d2gs = np.array(list(set(d2gs) - set(remove_d2gs)))
        if remove_mdusts:
            mdusts = np.array(list(set(mdusts) - set(remove_mdusts)))
        #
        rins.sort()
        routs.sort()
        mdusts.sort()
        d2gs.sort()
        f_o_deps.sort()
        f_c_deps.sort()

        nd2g = len(d2gs)
        nmdust = len(mdusts)
        nf_o_dep = len(f_o_deps)

        for id2g, d2g in enumerate(d2gs):
            for imdust, mdust in enumerate(mdusts):
                for if_o_dep, f_o_dep in enumerate(f_o_deps):
                    dic_select = {'rin': rins[0],
                                  'rout': routs[0],
                                  'd2g': d2g,
                                  'f_c_dep': 1.0,
                                  'f_o_dep': f_o_dep,
                                  'mdust': mdust, 
                                 }
                    LUT_this = select_LUT(LUT, dic_select)

                    if len(LUT_this) < 2:
                        print LUT_this, dic_select
                        raise Exception('Error in select LUT')
                    spflx = extract_xy(LUT_this, 'spectralType', 'INTFLUXL')
                    spflx.sort(cmp=f_cmp2)
                    sp, flx = zip(*spflx)
                    maxy, miny = max(maxy, max(flx)), min(miny, min(flx))

    for trans in fPatterns.keys():
        for mdust in mdusts:
            d_collect = collect_data(dobs, sources, trans, (None, mdust), mdusts)
            if len(d_collect) > 0:
                d_rms_s = [_[1] for _ in d_collect]
                maxy = max(maxy, max(d_rms_s))

    maxy *= 2
    miny = max(maxy/1e7, miny)

    for trans in fPatterns.keys():
        LUT = LUT_s[trans]
        lut0 = LUT[0]

        fig = pl.figure(figsize = figsize)
        pl.figtext(0.1, 1,
                'Molecule: ' + lut0['molname'] + \
                    '\t\tTransition: ' + \
                    trans.replace('-', r'\rightarrow') + \
                    #lut0['qnum'].replace('->', r'$\rightarrow$') + \
                    '\n' + 
                r'$\nu$: ' + '{0:.6e}'.format(lut0['nu']) + ' Hz\t\t' + \
                    r'$\lambda$: ' + '{0:.6e}'.format(lut0['lam']/1e4) + ' $\mu$m\t\t' + \
                r'$E_{{\rm up}}$: ' + '{0:.2f}'.format(lut0['Eup']) + ' K',
            fontsize = textsize,
            fontdict={'horizontalalignment': 'left', 'verticalalignment': 'bottom'})

        for id2g, d2g in enumerate(d2gs):
            for imdust, mdust in enumerate(mdusts):
                ax = make_sub_axis(fig, nmdust, nd2g, imdust, id2g,
                        pansepxfrac=0.04, pansepyfrac=0.0,
                        labelfontsize=labelfontsize, tickfontsize=textsize,
                        xtitle = 'Stellar type', ytitle = r'$\int T_{} dV$ (K km s$^{{-1}}$)',
                        xticks = [0,1,2,3,4,5],
                        xticklabels = ['M', 'K', 'G', 'F', 'A', 'B'],
                        yscale = 'log' if maxy/miny>30 else 'linear',
                        xRange = (-0.5, 5.5),
                        yRange = (miny, maxy),
                        less_xlabel=True, less_ylabel=True,
                     )

                ax.set_title(r'(d2g={0:.1e}, $m_{{\rm dust}}$={1:.1e})'.format(d2g, mdust),
                             fontsize=textsize,
                             fontdict={'verticalalignment': 'top', 'y':0.95})

                for if_o_dep, f_o_dep in enumerate(f_o_deps):
                    LUT_this = select_LUT(LUT,
                                     {'rin': rins[0],
                                      'rout': routs[0],
                                      'd2g': d2g,
                                      'f_c_dep': 1.0,
                                      'f_o_dep': f_o_dep,
                                      'mdust': mdust, 
                                      })

                    if len(LUT_this) < 2:
                        print 'Error in select LUT'
                        print LUT_this
                    spflx = extract_xy(LUT_this, 'spectralType', 'INTFLUXL')
                    spflx.sort(cmp=f_cmp2)
                    sp, flx = zip(*spflx)
                    flx = np.array(flx)
                    print sp, flx
                    flx[np.where(flx < 0.0)] = np.nan
                    ax.plot(flx, marker='o', markersize=8, markeredgecolor='none', color=colors[if_o_dep])

                if not (dobs is None or sources is None or trans is None or mdusts is None):
                    d_collect = collect_data(dobs, sources, trans, (None, mdust), mdusts)
                    if len(d_collect) > 0:
                        for d_c in d_collect:
                            #print d_c[2] == 0.0, d_c[2]
                            text_hor_shift = 0.4 * name2num(d_c[3])
                            if d_c[2] == 0.0: # Only upper limit
                                ax.plot([d_c[0]-0.15, d_c[0]+0.15], [d_c[1], d_c[1]], color=arrow_color)
                                ax.arrow(d_c[0], d_c[1], 0.0, -d_c[1]*(2.0/3.0),
                                         width = 0.02, overhang=0.5, zorder=100,
                                         head_length = d_c[1]*0.1, length_includes_head=True,
                                         edgecolor =arrow_color, facecolor =arrow_color)
                                ax.text(d_c[0]+0.1+text_hor_shift, d_c[1]/3, d_c[3], rotation="vertical", size="small", ha="left", va="center")
                            else:
                                #ax.plot([d_c[0]], [d_c[2]],
                                #        marker='s', markersize=2, markeredgecolor='none', color=errorbar_color)
                                ax.errorbar([d_c[0]], [d_c[2]], yerr=[d_c[1]],
                                            fmt='none', markersize=8,
                                            color=errorbar_color, ecolor=errorbar_color,
                                            capsize=6, capthick=1, elinewidth=4)
                                ax.text(d_c[0]+0.1+text_hor_shift, d_c[2], d_c[3], rotation="vertical", size="small", ha="left", va="center")

        fig_fname = de_latexmath(trans) + fname_postfix
        pl.savefig(os.path.join(fig_dir, fig_fname), bbox_inches='tight')
    
    return




def groupby(f, x):
    """
    Given the function f, group the list x into a dictionary based on the value of f evaluated for each item of x.

    Let's say x = [1,2,3,4,5], and f = lambda x: x % 2, then the result will be
        {0: [2,4], 1: [1,3,5]}

    Return the dictionary.
    """
    res = {}
    for item in x:
        key = f(item)
        if key in res.keys():
            res[key].append(item)
        else:
            res[key] = [item]
    return res


def orderby(f, x, cmpfunc=None):
    """
    Given the function f, sort the list x based on the value of f evaluated for each item of x.

    Let's say x = [1,2,3,4,5], and f = lambda x: x % 2, then the result will be
        [2,4,1,3,5]
    """
    if cmpfunc == None:
        cmpfunc = cmp
    return sorted(x, lambda x,y: cmpfunc(f(x), f(y)))




if __name__ == '__main__':
    task2(remove_d2gs=[1.0], remove_mdusts=[5e-5, 1e-4], figsize=(8,8))
