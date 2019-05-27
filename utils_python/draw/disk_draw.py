if __name__ == '__main__':
  from matplotlib import *
  use('Agg')
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
from numpy import *

from os.path import join as opj
from glob import glob

from parse_ana import *
from long_function_definitions import *
from my_script import *



def format_a_ana_reaction(s):
    def combine_spe_to_side(ss):
        sp = ss[0].strip()
        for p in ss[1:]:
            if p[0] != ' ':
                sp =  sp + ' + ' + p.strip()
        return sp
    r1 = s[0:12]
    r2 = s[12:24]
    p1 = s[24:36]
    p2 = s[36:48]
    p3 = s[48:60]
    p4 = s[60:72]
    sr = combine_spe_to_side([r1, r2])
    sp = combine_spe_to_side([p1, p2, p3, p4])
    return sr + r' $\rightarrow$ ' + sp


def get_ana_info(fpath, species=None, time = 1e6):
    c = contri(fpath)
    d = c.get_info(species, time)
    d1 = [format_a_ana_reaction(c.lines[i][54:125]) for i in range(*d['productions'])]
    d1 = d1[0:min(5,len(d1))]
    d['s_prod'] = '\n'.join(d1)
    d2 = [format_a_ana_reaction(c.lines[i][54:125]) for i in range(*d['destructions'])]
    d2 = d2[0:min(5,len(d2))]
    d['s_dest'] = '\n'.join(d2)
    return d


def draw_name_list(d, name_list, figsize = (12, 9), xs = 'linear', ys = 'linear',
                   model_name='', now_dir='',
                   draw_rectangles = True,
                   panelwidth = 5,
                   panelheight = 5,
                   overlay_list = None,
                   xr = (0,800), yr = (0,150), savefig_dir = None):
    #
    npanx = ceil(sqrt(len(name_list)))
    npany = ceil(float(len(name_list))/npanx)
    figsize = (npanx * panelwidth, npany * panelheight)
    fig, ax_s = plot_model_results_map(d, name_list, display_dic=display_dic, all_in_one_fig=True,
                                   majorgridon=False, xmingrid=False, minorgridon=False, graygrid=True,
                                   xRange=xr, yRange=yr, xscale=xs, yscale=ys, linthreshy=1,
                                   draw_rectangles=draw_rectangles, nx=100, ny=100, figsize=figsize)
    #
    if overlay_list != None:
        for top, bot, levels in overlay_list:
            ax = ax_s[bot]
            overlap_contour(fig, ax, d, top, levels, nxy=100, linewidth=1, color='black', fontsize=12,
                            method='nearest', fmt='%1.3f', manual=False, nwidth=1, whichx='rmin', whichy='zmin')
    #
    fig_fname = model_name.replace('/', '') + '_' + make_fname_from_namelist(name_list) + '.pdf'
    if savefig_dir == None:
        savefig_dir = now_dir
    savefig(opj(savefig_dir, fig_fname), bbox_inches='tight')
    print 'File saved at ', opj(savefig_dir, fig_fname)
    return


def draw_radial(dic, name_list, zpos, savefig_dir=None,
                xRange=None, xscale='linear', model_name='', now_dir=''):
    #
    plot_model_results_radial(dic, name_list, zpos, xRange=xRange, figsize = None,
        xscale=xscale, yscale='log', returndata=True, noLegend=False, noTitle=False,
        color=None, label=None, display_dic=None,
        ax=None, returnaxis=False, linestyle='-',lw=3,markersize=5,
        do_interpol=False, interpolNum=200,
        dosmooth=False, winwidth=None,
        whichx='rmin')
    fig_fname = model_name.replace('/', '') + '_' + make_fname_from_namelist(name_list) + '_radial.pdf'
    if savefig_dir == None:
        savefig_dir = opj(now_dir, '')
    savefig(opj(savefig_dir, fig_fname), bbox_inches='tight')
    print 'File saved at ', opj(savefig_dir, fig_fname)
    return


def de_weird(d, name, thrsh=10.0, minv=0.0):
    import numpy as np
    dc = d[name]
    n = len(dc)
    for i in range(1,n-1):
        if ((d['rmin'][i+1] == d['rmin'][i-1]) and
            (min(dc[i-1], dc[i]) > minv) and
            ((max(dc[i]/dc[i-1], dc[i-1]/dc[i]) > thrsh) or
             (max(dc[i]/dc[i+1], dc[i+1]/dc[i]) > thrsh))):
            #
            dc[i] = np.sqrt(dc[i-1] * dc[i+1])



mycm = make_my_colormap(c_list=[(0.0, (0.5, 0.0, 0.5)),
                                (0.2, (0.0, 0.0, 1.0)),
                                (0.4, (0.0, 0.8, 1.0)),
                                (0.6, (0.0, 0.8, 0.0)),
                                (0.8, (1.0, 0.8, 0.0)),
                                (1.0, (1.0, 0.0, 0.0))])

rcParams['axes.color_cycle'] = mycolors


def task1(now_dir = '/n/Users/fdu/now/',
          model_name = '20150307_dep_s',
          iter_fname = 'iter_0001.dat',
          xr = None, yr = None,
          justRadial = False,
          ):

    savefig_dir = opj(now_dir, 'figures/')
    res_dir = opj(now_dir, 'res/')
    model_dir = opj(res_dir, model_name)

    data_fname = opj(model_dir, iter_fname)

    d = load_data_as_dic(data_fname)
    update_stuff(d, data_fname)

    calc_colden(d, 'C')
    calc_colden(d, 'CO')
    calc_colden(d, 'gCO')
    calc_colden(d, 'H2O')
    calc_colden(d, 'gH2O')
    calc_colden(d, 'CO2')
    calc_colden(d, 'gCO2')
    calc_colden(d, 'CH4')
    calc_colden(d, 'C2H')
    calc_colden(d, 'C3H2')
    calc_colden(d, 'C2H3')
    #
    draw_radial(d,
        [
         [
            {'name': 'Tgas',  'vr': (10, 1e4)},
            #{'name': 'Tdust',  'vr': (10, 2e3)},
          # {'name': 'Ncol_I',  'vr': (1e20, 1e27)},
         ],
        #[\
        #   {'name': 'Ncol_I',  'vr': (1e9, 1e27)},
        #  #{'name': 'N(C)',    'vr': (1e4, 1e22)},
        #  #{'name': 'N(CO)',   'vr': (1e4, 1e22)},
        #  #{'name': 'N(H2O)',  'vr': (1e4, 1e22)},
        #   {'name': 'N(C2H)',  'vr': (1e9, 1e19)},
        #  #{'name': 'N(C3H2)', 'vr': (1e4, 1e19)},
        #  #{'name': 'N(C2H3)', 'vr': (1e4, 1e19)},
        #],
        ],
        2.0,
        model_name=model_name, now_dir=now_dir, xscale='log', xRange=(1,200))

    if justRadial:
        return

    d['C/O'] = d['X[C]'] / d['X[O]']
    d['G0_UV'] = d['G0_UV'] + d['UV_G0_I'] * exp(-d['Av_ISM'])
    draw_name_list(d,
        model_name=model_name,
        now_dir=now_dir,
        name_list = \
            [\
             [{'name': 'n_gas', 'linearscale': False, 'cmap': mycm, 'vr': (1e3, 1e13),
              'name_disp': r'$n_{\rm gas}$ (cm$^{-3}$)'}],
             [{'name': 'Tgas', 'linearscale': False, 'cmap': mycm, 'vr': (1e1, 3e3),
              'name_disp': r'$T_{\rm gas}$ (K)'}],
             [{'name': 'Tdust', 'linearscale': False, 'cmap': mycm, 'vr': (1e1, 3e3),
              'name_disp': r'$T_{\rm dust}$ (K)'}],
             [{'name': 'G0_UV', 'linearscale': False, 'cmap': mycm, 'vr': (1e-3, 1e7),
              'name_disp': r'$G_{{\rm 0,UV}}$'}],
             [{'name': 'Av_ISM', 'linearscale': False, 'cmap': mycm, 'vr': (1e-2, 1e5),
              'name_disp': r'Av'}],
             [{'name': 'O', 'linearscale': False, 'cmap': mycm, 'vr': (1e-10, 3e-4),
              'name_disp': r'O'}],
             [{'name': 'C', 'linearscale': False, 'cmap': mycm, 'vr': (1e-10, 3e-4),
              'name_disp': r'C'}],
             [{'name': 'X[O]', 'linearscale': False, 'cmap': mycm, 'vr': (1e-8, 3e-4),
              'name_disp': r'$X[{\rm O}]$'}],
             [{'name': 'X[C]', 'linearscale': False, 'cmap': mycm, 'vr': (1e-8, 3e-4),
              'name_disp': r'$X[{\rm C}]$'}],
             [{'name': 'C/O', 'linearscale': False, 'cmap': mycm, 'vr': (1e-1, 1e5),
              'name_disp': r'C/O'}],
             [{'name': 'CO', 'linearscale': False, 'cmap': mycm, 'vr': (1e-14, 1e-4),
              'name_disp': r'CO'}],
             [{'name': 'N_CO_I', 'linearscale': False, 'cmap': mycm, 'vr': (1e8, 1e19),
              'name_disp': r'N(CO)'}],
             [{'name': 'gCO', 'linearscale': False, 'cmap': mycm, 'vr': (1e-14, 1e-4),
              'name_disp': r'gCO'}],
             [{'name': 'N(gCO)', 'linearscale': False, 'cmap': mycm, 'vr': (1e8, 1e19),
              'name_disp': r'N(gCO)'}],
             [{'name': 'CO2', 'linearscale': False, 'cmap': mycm, 'vr': (1e-15, 1e-5),
              'name_disp': r'CO2'}],
             [{'name': 'N(CO2)', 'linearscale': False, 'cmap': mycm, 'vr': (1e7, 1e17),
              'name_disp': r'N(CO2)'}],
             [{'name': 'gCO2', 'linearscale': False, 'cmap': mycm, 'vr': (1e-15, 1e-7),
              'name_disp': r'gCO2'}],
             [{'name': 'N(gCO2)', 'linearscale': False, 'cmap': mycm, 'vr': (1e7, 1e17),
              'name_disp': r'N(gCO2)'}],
             [{'name': 'H2O', 'linearscale': False, 'cmap': mycm, 'vr': (1e-15, 1e-6),
              'name_disp': r'H2O'}],
             [{'name': 'N_H2O_I', 'linearscale': False, 'cmap': mycm, 'vr': (1e8, 1e18),
              'name_disp': r'N(H2O)'}],
             [{'name': 'gH2O', 'linearscale': False, 'cmap': mycm, 'vr': (1e-15, 1e-6),
              'name_disp': r'gH2O'}],
             [{'name': 'N(gH2O)', 'linearscale': False, 'cmap': mycm, 'vr': (1e8, 1e18),
              'name_disp': r'N(gH2O)'}],
             [{'name': 'C2H', 'linearscale': False, 'cmap': mycm, 'vr': (1e-15, 1e-7),
              'name_disp': r'C2H'}],
             [{'name': 'N(C2H)', 'linearscale': False, 'cmap': mycm, 'vr': (1e6, 1e16),
              'name_disp': r'N(C2H)'}],
             [{'name': 'C3H2', 'linearscale': False, 'cmap': mycm, 'vr': (1e-15, 1e-7),
              'name_disp': r'C3H2'}],
             [{'name': 'N(C3H2)', 'linearscale': False, 'cmap': mycm, 'vr': (1e6, 1e16),
              'name_disp': r'N(C3H2)'}],
             [{'name': 'N2', 'linearscale': False, 'cmap': mycm, 'vr': (1e-12, 1e-4),
              'name_disp': r'N2'}],
             [{'name': 'gN2', 'linearscale': False, 'cmap': mycm, 'vr': (1e-12, 1e-4),
              'name_disp': r'gN2'}],
             [{'name': 'NH3', 'linearscale': False, 'cmap': mycm, 'vr': (1e-12, 1e-4),
              'name_disp': r'NH3'}],
             [{'name': 'gNH3', 'linearscale': False, 'cmap': mycm, 'vr': (1e-12, 1e-4),
              'name_disp': r'gNH3'}],
            ],
        #name_list = \
        #    [\
        #     [{'name': 'n_gas', 'linearscale': False, 'cmap': mycm, 'vr': (1e3, 1e13),
        #      'name_disp': r'$n_{\rm gas}$ (cm$^{-3}$)'}],
        #     [{'name': 'Tgas', 'linearscale': False, 'cmap': mycm, 'vr': (1e1, 3e3),
        #      'name_disp': r'$T_{\rm gas}$ (K)'}],
        #     [{'name': 'Tdust', 'linearscale': False, 'cmap': mycm, 'vr': (1e1, 3e3),
        #      'name_disp': r'$T_{\rm dust}$ (K)'}],
        #     [{'name': 'G0_UV', 'linearscale': False, 'cmap': mycm, 'vr': (1e-3, 1e7),
        #      'name_disp': r'$G_{{\rm 0,UV}}$'}],
        #     [{'name': 'Av_ISM', 'linearscale': False, 'cmap': mycm, 'vr': (1e-2, 1e5),
        #      'name_disp': r'Av'}],
        #     [{'name': 'C2H', 'linearscale': False, 'cmap': mycm, 'vr': (1e-15, 1e-7),
        #      'name_disp': r'C2H'}],
        #    ],
        xr = xr, yr = yr,
        )



def task1a(now_dir = '/n/Users/fdu/now/',
          model_name = '20150307_dep_s',
          iter_fname = 'iter_0001.dat',
          xr = None, yr = None,
          C_O_ratio_r_thr=None,
          C_O_ratio_zr_thr=None,
          freq_rough = None,
          ):

    savefig_dir = opj(now_dir, 'figures/')
    res_dir = opj(now_dir, 'res/')
    model_dir = opj(res_dir, model_name)

    data_fname = opj(model_dir, iter_fname)

    d = load_data_as_dic(data_fname)
    update_stuff(d, data_fname)
    #
    all_ions = get_all_ions(d.keys())
    maj_ions = get_major_ions(d, all_ions)
    d['major ions'] = [ion[0] for ion in maj_ions]
    ion_counts = counts_ions_weighted(maj_ions)
    cdict_ions = make_colors_for_ion(ion_counts)
    def draw_ion_legend(ax):
        return draw_rect_color_legend(ax, cdict_ions, ion_counts, nlim=15)
    #
    maj_C = get_major_species_for_element('C', d, data_fname, all_elements)
    d['major C'] = [c[0] for c in maj_C]
    C_counts = counts_ions_weighted(maj_C)
    cdict_C = make_colors_for_ion(C_counts)
    def draw_C_legend(ax):
        return draw_rect_color_legend(ax, cdict_C, C_counts, nlim=15)
    #
    maj_N = get_major_species_for_element('N', d, data_fname, all_elements)
    d['major N'] = [c[0] for c in maj_N]
    N_counts = counts_ions_weighted(maj_N)
    cdict_N = make_colors_for_ion(N_counts)
    def draw_N_legend(ax):
        return draw_rect_color_legend(ax, cdict_N, N_counts, nlim=15)
    #
    maj_O = get_major_species_for_element('O', d, data_fname, all_elements)
    d['major O'] = [c[0] for c in maj_O]
    O_counts = counts_ions_weighted(maj_O)
    cdict_O = make_colors_for_ion(O_counts)
    def draw_O_legend(ax):
        return draw_rect_color_legend(ax, cdict_O, O_counts, nlim=15)
    #
    calc_colden(d, 'C2H')
    calc_colden(d, 'CO')
    calc_colden(d, 'rhodus_1', is_abundance=False)
    calc_colden(d, 'rhodus_2', is_abundance=False)
    calc_colden(d, 'rhodus_3', is_abundance=False)
    d['n(C2H)'] = d['C2H'] * d['n_gas']
    #
    #de_weird(d, 'CO', thrsh=1e2, minv=1e-12)
    #
    def opti_str_2_kappa(st):
        s = map(float, st.split())
        ab_s = [s[0], s[3], s[6]]
        sc_s = [s[1], s[4], s[7]]
        g_s  = [s[2], s[5], s[8]]
        return ab_s, sc_s, g_s

    if freq_rough == 262:
        ab_s, sc_s, _ = opti_str_2_kappa('''2.63555E+000   7.87484E+000  -2.86927E-007
                    2.89751E-001   6.75577E-007  -2.78497E-007   2.98063E-001   1.26782E-005   3.00201E-005''')
        tau_freq_str = ' (262 GHz)'
    elif freq_rough == 349:
        ab_s, sc_s, _ = opti_str_2_kappa('''4.10727E+000   1.38004E+001  -1.52808E-006
                    4.84421E-001   2.13443E-006  -1.51790E-006   4.96661E-001   4.00530E-005   4.23444E-005''')
        tau_freq_str = ' (349 GHz)'
    elif freq_rough == 115:
        ab_s, sc_s, _ = opti_str_2_kappa('''4.69258E-001   2.49644E-001  -2.86927E-007
                    5.15899E-002   2.14168E-008  -2.78497E-007   5.30699E-002   4.01919E-007   3.00201E-005''')
        tau_freq_str = ' (115 GHz)'
    elif freq_rough == 230:
        ab_s, sc_s, _ = opti_str_2_kappa('''1.86629E+000   3.94873E+000  -2.86927E-007
                    2.05179E-001   3.38759E-007  -2.78497E-007   2.11065E-001   6.35732E-006   3.00201E-005''')
        tau_freq_str = ' (230 GHz)'
    else:
        raise "Error freq_rough"
    #print ab_s, sc_s
    ex_s = np.array(ab_s) + np.array(sc_s)
    d['tau_1'] = d['N(rhodus_1)'] * ex_s[0]
    d['tau_2'] = d['N(rhodus_2)'] * ex_s[1]
    d['tau_3'] = d['N(rhodus_3)'] * ex_s[2]
    d['tau_tot'] = d['tau_1'] + d['tau_2'] + d['tau_3']
    #
    update_stuff(d, data_fname)
    d['C/O'] = d['X[C]'] / d['X[O]']
    if C_O_ratio_zr_thr != None and C_O_ratio_r_thr != None:
        d['C/O'][logical_and(d['zmax']/d['rmax'] > C_O_ratio_zr_thr,
                             d['rmin']>C_O_ratio_r_thr)] = d['C/O'].max()
    d['gCO/gH2O'] = d['gCO'] / d['gH2O']
    #d['G0_UV'] = d['G0_UV'] + d['UV_G0_I'] * exp(-d['Av_ISM'])
    phy_UVext2Av = 2.6
    d['G0_UV'] = d['UV_G0_S'] * exp(-phy_UVext2Av*d['Av_Star']) + \
                 d['UV_G0_I'] * exp(-phy_UVext2Av*d['Av_ISM'])
    #
    #   def min_nonezero(x):
    #       m = x[0]
    #       for t in x[1:]:
    #           if t > 1 and t < m:
    #               m = t
    #       return m
    #   #
    #   d['delT'] = d['Tgas'] - np.array([min_nonezero(x) for x in zip(d['Tdust1'],  d['Tdust2'],  d['Tdust3'],  d['Tdust4'])])
    #
    draw_name_list(d,
        model_name=model_name,
        now_dir=now_dir,
        draw_rectangles = True,
        name_list = \
            [\
                [{'name': 'n_gas', 'linearscale': False, 'cmap': mycm, 'vr': (1e3, 1e13),
                 'name_disp': r'$n_{\rm gas}$ (cm$^{-3}$)'}],
                [{'name': 'rhodus_1', 'linearscale': False, 'cmap': mycm, 'vr': (1e-23, 1e-13),
                 'name_disp': r'$\rho_{\rm mm\ dust}$ (g cm$^{-3}$)'}],
                [{'name': 'rhodus_2', 'linearscale': False, 'cmap': mycm, 'vr': (1e-23, 1e-13),
                 'name_disp': r'$\rho_{\mu{\rm m}{\rm\ dust}}$ (g cm$^{-3}$)'}],
          #     [{'name': 'tau_tot', 'linearscale': False, 'cmap': mycm, 'vr': (1e-2, 1e2),
          #      'name_disp': r'$\tau_{\rm dust}$' + tau_freq_str}],
                [{'name': 'C2H', 'linearscale': False, 'cmap': mycm, 'vr': (1e-15, 1e-7),
                 'name_disp': r'X[C2H]'}],
              # [{'name': 'n(C2H)', 'linearscale': False, 'cmap': mycm, 'vr': (1e-7, 1e4),
              #  'name_disp': r'n(C2H)'}],
              # [{'name': 'N(C2H)', 'linearscale': False, 'cmap': mycm, 'vr': (1e6, 1e16),
              #  'name_disp': r'N(C2H)'}],
               [{'name': 'Ncol_I', 'linearscale': False, 'cmap': mycm, 'vr': (1e16, 1e26),
                'name_disp': r'$N_{{\rm H}}$ cm$^{-2}$'}],
                [{'name': 'gCO', 'linearscale': False, 'cmap': mycm, 'vr': (1e-12, 1e-3),
                 'name_disp': r'gCO'}],
                [{'name': 'gCO2', 'linearscale': False, 'cmap': mycm, 'vr': (1e-12, 1e-4),
                 'name_disp': r'gCO2'}],
                [{'name': 'gH2O', 'linearscale': False, 'cmap': mycm, 'vr': (1e-12, 1e-4),
                 'name_disp': r'gH2O'}],
                [{'name': 'CO', 'linearscale': False, 'cmap': mycm, 'vr': (1e-10, 1e-4),
                 'name_disp': r'CO'}],
                [{'name': 'CO2', 'linearscale': False, 'cmap': mycm, 'vr': (1e-12, 1e-4),
                 'name_disp': r'CO2'}],
                [{'name': 'H2O', 'linearscale': False, 'cmap': mycm, 'vr': (1e-12, 1e-4),
                 'name_disp': r'H2O'}],
                [{'name': 'C', 'linearscale': False, 'cmap': mycm, 'vr': (1e-13, 1e-4),
                 'name_disp': r'C'}],
                [{'name': 'N(CO)', 'linearscale': False, 'cmap': mycm, 'vr': (1e13, 1e23),
                 'name_disp': r'N(CO)'}],
                [{'name': 'C/O', 'linearscale': False, 'cmap': mycm, 'vr': (1e-4, 1e4),
                 'name_disp': r'C/O'}],
                [{'name': 'X[C]', 'linearscale': False, 'cmap': mycm, 'vr': (1e-12, 2e-3),
                 'name_disp': r'X[C]'}],
                [{'name': 'X[O]', 'linearscale': False, 'cmap': mycm, 'vr': (1e-12, 3e-3),
                 'name_disp': r'X[O]'}],
              # [{'name': 'gCO/gH2O', 'linearscale': False, 'cmap': mycm, 'vr': (1e-10, 1e-4),
              #  'name_disp': r'gCO/gH2O'}],
                [{'name': 'Tgas', 'linearscale': False, 'cmap': mycm, 'vr': (1e1, 1000),
                 'name_disp': r'$T_{\rm gas}$ (K)'}],
                [{'name': 'Tdust', 'linearscale': False, 'cmap': mycm, 'vr': (1e1, 1000),
                 'name_disp': r'$T_{\rm dust}$ (K)'}],
                #[{'name': 'delT', 'linearscale': True, 'cmap': mycm, 'vr': (-50,50),
                # 'name_disp': r'$T_{\rm gas}-T_{\rm dust}$ (K)'}],
                [{'name': 'G0_UV', 'linearscale': False, 'cmap': mycm, 'vr': (1e-3, 1e7),
                 'name_disp': r'$G_{{\rm 0,UV}}$'}],
                [{'name': 'E-', 'linearscale': False, 'cmap': mycm, 'vr': (1e-12, 1e-3),
                 'name_disp': r'E-'}],
                [{'name': 'major ions', 'cdict': cdict_ions, 'showColorbar': False,
                  'legendFunc': draw_ion_legend}],
                [{'name': 'major C', 'cdict': cdict_C, 'showColorbar': False,
                  'legendFunc': draw_C_legend}],
                #[{'name': 'major N', 'cdict': cdict_N, 'showColorbar': False,
                #  'legendFunc': draw_N_legend}],
                [{'name': 'major O', 'cdict': cdict_O, 'showColorbar': False,
                  'legendFunc': draw_O_legend}],
                [{'name': 'zeta_X', 'linearscale': False, 'cmap': mycm, 'vr': (1e-17, 1e-7),
                 'name_disp': r'$\zeta_{{\rm X-ray}}$ (s$^{-1}$)'}],
            ],
        xr = xr, yr = yr,
        panelwidth = 5,
        panelheight = 5 * (yr[1]-yr[0]) / (xr[1]-xr[0]),
        #overlay_list = [\
        #        ['tau_tot', 'N(C2H)', (0.1, 0.5, 2)],
        #        ['tau_tot', 'N(CO)',  (0.1, 0.5, 2)],
        #    ],
        )

def task2(now_dir = '/n/Users/fdu/now/',
          model_name = '20150307_noD_d_1comp',
          iter_fname = 'iter_0001.dat',
          xr = None, yr = None,
          ):

    import os

    res_dir = opj(now_dir, 'res/')
    model_dir = opj(res_dir, model_name)
    savefig_dir = opj(model_dir, 'figures/')
    if not os.path.exists(savefig_dir):
        os.mkdir(savefig_dir)

    data_fname = opj(model_dir, iter_fname)

    d = load_data_as_dic(data_fname)
    update_stuff(d, data_fname)

    d['C/O'] = d['X[C]'] / d['X[O]']
    d['n(C2H)'] = d['C2H'] * d['n_gas']
    #draw_name_list(d,
    #    model_name=model_name,
    #    now_dir=now_dir,
    #    savefig_dir = savefig_dir,
    #    xr = xr, yr = yr,
    #    draw_rectangles = True,
    #    name_list = \
    #        [\
    #            [{'name': 'n_gas', 'linearscale': False, 'cmap': mycm, 'vr': (1e3, 1e13),
    #             'name_disp': r'$n_{\rm gas}$ (cm$^{-3}$)'}],
    #            #[{'name': 'Tgas', 'linearscale': False, 'cmap': mycm, 'vr': (1e1, 3e3),
    #            # 'name_disp': r'$T_{\rm gas}$ (K)'}],
    #            [{'name': 'C/O', 'linearscale': False, 'cmap': mycm, 'vr': (1e-1, 1e3),
    #             'name_disp': r'C/O'}],
    #            #[{'name': 'Tdust', 'linearscale': False, 'cmap': mycm, 'vr': (1e1, 3e3),
    #            # 'name_disp': r'$T_{\rm dust}$ (K)'}],
    #            #[{'name': 'C2H', 'linearscale': False, 'cmap': mycm, 'vr': (1e-12, 1e-7),
    #            # 'name_disp': r'C2H'}],
    #            [{'name': 'X[C]', 'linearscale': False, 'cmap': mycm, 'vr': (1e-10, 1e-3),
    #             'name_disp': r'X[C]'}],
    #            [{'name': 'X[O]', 'linearscale': False, 'cmap': mycm, 'vr': (1e-10, 1e-3),
    #             'name_disp': r'X[O]'}],
    #            [{'name': 'n(C2H)', 'linearscale': False, 'cmap': mycm, 'vr': (1e-5, 1e3),
    #             'name_disp': r'n(C2H)'}],
    #            [{'name': 'G0_UV', 'linearscale': False, 'cmap': mycm, 'vr': (1e-3, 1e7),
    #             'name_disp': r'$G_{{\rm 0,UV}}$'}],
    #        ],
    #    )
    draw_name_list(d,
        model_name=model_name,
        now_dir=now_dir,
        savefig_dir = savefig_dir,
        xr = xr, yr = yr,
        draw_rectangles = True,
        name_list = \
            [\
                [{'name': 'n_gas', 'linearscale': False, 'cmap': mycm, 'vr': (1e3, 1e13),
                 'name_disp': r'$n_{\rm gas}$ (cm$^{-3}$)'}],
                [{'name': 'Tgas', 'linearscale': False, 'cmap': mycm, 'vr': (1e1, 3e3),
                 'name_disp': r'$T_{\rm gas}$ (K)'}],
                [{'name': 'Tdust', 'linearscale': False, 'cmap': mycm, 'vr': (1e1, 3e3),
                 'name_disp': r'$T_{\rm dust}$ (K)'}],
                [{'name': 'gH2O', 'linearscale': False, 'cmap': mycm, 'vr': (1e-12, 5e-4),
                 'name_disp': r'gH2O'}],
                [{'name': 'H2O', 'linearscale': False, 'cmap': mycm, 'vr': (1e-12, 5e-4),
                 'name_disp': r'H2O'}],
                [{'name': 'CO', 'linearscale': False, 'cmap': mycm, 'vr': (1e-12, 5e-4),
                 'name_disp': r'CO'}],
                [{'name': 'gCO', 'linearscale': False, 'cmap': mycm, 'vr': (1e-16, 5e-4),
                 'name_disp': r'gCO'}],
                [{'name': 'CO2', 'linearscale': False, 'cmap': mycm, 'vr': (1e-12, 5e-4),
                 'name_disp': r'CO2'}],
                [{'name': 'gCO2', 'linearscale': False, 'cmap': mycm, 'vr': (1e-16, 5e-4),
                 'name_disp': r'gCO2'}],
                [{'name': 'NO2', 'linearscale': False, 'cmap': mycm, 'vr': (1e-16, 5e-4),
                 'name_disp': r'NO2'}],
                [{'name': 'gNO2', 'linearscale': False, 'cmap': mycm, 'vr': (1e-16, 5e-4),
                 'name_disp': r'gNO2'}],
                [{'name': 'NO', 'linearscale': False, 'cmap': mycm, 'vr': (1e-16, 5e-4),
                 'name_disp': r'NO'}],
                [{'name': 'gNO', 'linearscale': False, 'cmap': mycm, 'vr': (1e-20, 1e-10),
                 'name_disp': r'gNO'}],
                [{'name': 'OH', 'linearscale': False, 'cmap': mycm, 'vr': (1e-16, 1e-6),
                 'name_disp': r'OH'}],
                [{'name': 'gOH', 'linearscale': False, 'cmap': mycm, 'vr': (1e-16, 1e-10),
                 'name_disp': r'gOH'}],
                [{'name': 'gH', 'linearscale': False, 'cmap': mycm, 'vr': (1e-22, 1e-8),
                 'name_disp': r'gH'}],
                [{'name': 'H', 'linearscale': False, 'cmap': mycm, 'vr': (1e-12, 5e-1),
                 'name_disp': r'H'}],
                [{'name': 'gH2', 'linearscale': False, 'cmap': mycm, 'vr': (1e-16, 1e-6),
                 'name_disp': r'gH2'}],
                [{'name': 'H2', 'linearscale': False, 'cmap': mycm, 'vr': (1e-12, 5e-1),
                 'name_disp': r'H2'}],
                [{'name': 'C', 'linearscale': False, 'cmap': mycm, 'vr': (1e-10, 1e-3),
                 'name_disp': r'C'}],
                [{'name': 'C2H', 'linearscale': False, 'cmap': mycm, 'vr': (1e-15, 1e-7),
                 'name_disp': r'C2H'}],
                [{'name': 'O', 'linearscale': False, 'cmap': mycm, 'vr': (1e-10, 1e-3),
                 'name_disp': r'O'}],
                [{'name': 'gO', 'linearscale': False, 'cmap': mycm, 'vr': (1e-20, 1e-10),
                 'name_disp': r'gO'}],
                [{'name': 'G0_UV', 'linearscale': False, 'cmap': mycm, 'vr': (1e-3, 1e7),
                 'name_disp': r'$G_{{\rm 0,UV}}$'}],
                [{'name': 'flx_Xray', 'linearscale': False, 'cmap': mycm, 'vr': (1e-8, 1e4),
                 'name_disp': r'$F_{{\rm X}}$'}],
            ],
        )


def task3():
    ds = []
    for mname in ['20150307_noD_d_rerun2a', '20150307_noD_d_1comp_rerun']:
        fname = opj('/n/Users/fdu/now/res', mname, 'iter_0001.dat')
        d = load_data_as_dic(fname)
        update_stuff(d, fname)
        d['n(C2H)'] = d['C2H'] * d['n_gas']
        ds.append(d)
    name_list = \
            [\
             [{'name': 'n_gas', 'linearscale': False, 'cmap': mycm, 'vr': (1e3, 1e13),
              'name_disp': r'$n_{\rm gas}$ (cm$^{-3}$)'}],
             [{'name': 'G0_UV', 'linearscale': False, 'cmap': mycm, 'vr': (1e-3, 1e7),
               'name_disp': r'$G_{{\rm 0,UV}}$'}],
             [{'name': 'n(C2H)', 'linearscale': False, 'cmap': mycm, 'vr': (1e-7, 1e3),
              'name_disp': r'$n_{\rm C_{2}H}$ (cm$^{-3}$)'}],
            ]
    xr = (0,100)
    yr = (0,60)
    xs, ys = 'linear', 'linear'
    draw_rectangles = True
    figsize = (15,6)
    fig = plt.figure(figsize = figsize)
    p_w = 0.7 / 3
    p_h = 0.8 / 2
    sep_hor = 0.05
    sep_ver = 0.05
    for i in range(2):
        d = ds[i]
        for j in range(3):
            pos = (0.1+(p_w+sep_hor)*j, 0.1+(p_h+sep_ver)*i, p_w, p_h)
            if i > 0:
                xlabel = ''
                xticks = []
            else:
                xlabel = 'r (AU)'
            if j > 0:
                ylabel = ''
                yticks = []
            else:
                ylabel = 'z (AU)'
            name_list_this = [name_list[j]]
            ax = fig.add_axes(pos,
                  xlabel=xlabel,
                  ylabel=ylabel,
                  autoscalex_on=False, autoscaley_on=False,
                  xscale=xs, yscale=ys, xlim=xr, ylim=yr)
            plot_model_results_map(d, name_list_this,
                       fig_prev = fig, ax_prev = ax, panel_pos = pos,
                       display_dic=display_dic, all_in_one_fig=False,
                       majorgridon=False, xmingrid=False, minorgridon=False, graygrid=True,
                       xRange=xr, yRange=yr, xscale=xs, yscale=ys, linthreshy=1,
                       draw_rectangles=draw_rectangles, nx=200, ny=200, figsize=figsize)
            if i > 0:
                ax.set_xticklabels(xticks)
            if j > 0:
                ax.set_yticklabels(yticks)
    plt.savefig('/n/Users/fdu/now/fig4_cont_rerun2a_n.pdf', bbox_inches='tight')



def task4(now_dir = '/n/Users/fdu/now/',
          model_name = '20150307_dep_s',
          iter_fname = 'iter_0001.dat',
          xr = None, yr = None,
          ):

    savefig_dir = opj(now_dir, 'figures/')
    res_dir = opj(now_dir, 'res/')
    model_dir = opj(res_dir, model_name)

    data_fname = opj(model_dir, iter_fname)

    d = load_data_as_dic(data_fname)
    update_stuff(d, data_fname)

    Q_3micron = 1.4e-2

    d['n(H2O)'] = d['H2O'] * d['n_gas']
    d['n(gH2O)'] = d['gH2O'] * d['n_gas']
    d['tau'] = d['Ncol_I'] * d['d2gnum'] * d['sigd_av'] * Q_3micron

    draw_name_list(d, panelwidth=6, panelheight=5,
        model_name=model_name,
        now_dir=now_dir,
        name_list = \
            [\
             [{'name': 'n_gas', 'linearscale': False, 'cmap': mycm, 'vr': (1e3, 1e13),
              'name_disp': r'$n_{\rm gas}$ (cm$^{-3}$)'}],
             [{'name': 'n(H2O)', 'linearscale': False, 'cmap': mycm, 'vr': (1e-9, 1e1),
              'name_disp': r'n(water vapor; cm$^{-3}$)'}],
             [{'name': 'n(gH2O)', 'linearscale': False, 'cmap': mycm, 'vr': (1e-9, 1e1),
              'name_disp': r'n(water ice; cm$^{-3}$)'}],
             [{'name': 'tau', 'linearscale': False, 'cmap': mycm, 'vr': (1e-3, 1e1),
              'name_disp': r'$\tau$ (3 $\mu$m; $r=0.1\mu$m)'}],
            ],
            xr = xr, yr = yr,
            )



def task4a(now_dir = '/n/Fdu1/work/grid_20150703_redo/',
          model_name = 'run_004',
          savefig_dir=None,
          iter_fname = 'iter_0001.dat',
          xr = None, yr = None,
          ):

    if not savefig_dir:
        savefig_dir = opj(now_dir, 'figures/')
    res_dir = opj(now_dir, 'results/')
    model_dir = opj(res_dir, model_name)

    data_fname = opj(model_dir, iter_fname)

    d = load_data_as_dic(data_fname)
    update_stuff(d, data_fname)

    draw_name_list(d, panelwidth=6, panelheight=5,
        model_name=model_name,
        now_dir=now_dir,
        draw_rectangles=True,
        name_list = \
            [\
             [{'name': 'n_gas', 'linearscale': False, 'cmap': mycm, 'vr': (1e5, 1e12),
              'name_disp': r'$n_{\rm gas}$ (cm$^{-3}$)'}],
            #[{'name': 'Tgas', 'linearscale': False, 'cmap': mycm, 'vr': (1e1, 1000),
            # 'name_disp': r'$T_{\rm gas}$ (K)'}],
             [{'name': 'Tdust', 'linearscale': False, 'cmap': mycm, 'vr': (1e1, 1000),
              'name_disp': r'$T_{\rm dust}$ (K)'}],
             [{'name': 'G0_UV', 'linearscale': False, 'cmap': mycm, 'vr': (1e-3, 1e7),
              'name_disp': r'$G_{{\rm 0,UV}}$'}],
             [{'name': 'H2O', 'linearscale': False, 'cmap': mycm, 'vr': (1e-12, 1e-6),
              'name_disp': r'H2O'}],
            #[{'name': 'gH2O', 'linearscale': False, 'cmap': mycm, 'vr': (1e-12, 5e-4),
            #  'name_disp': r'gH2O'}],
            ],
            xr = xr, yr = yr,
            )


def task5(data_fname):
    items, d = load_data_as_dic(data_fname, returnOriginalKeys=True)
    nitem = len(items)
    panel_wh = (5, 3)
    xr = (0, 200)
    yr = (0, 120)

    def choose_items(header, data, n):
        nitem = len(header)
        for i in range(0, nitem, n):
            names = header[i : min(i + n, nitem)]
            lis = []
            for name in names:
                vmax, vmin, vmed = np.nanmax(d[name]), np.nanmin(d[name]), np.median(d[name])
                if vmax <= 0:
                    linscale = True
                elif vmed / vmax > 0.1:
                    linscale = True
                else:
                    linscale = False
                    if vmin / vmax < 1e-15:
                        vmin = vmax * 1e-15
                lis.append([{'name': name}])
            d
            yield lis


    from matplotlib.backends.backend_pdf import PdfPages
    pdf_pages = PdfPages('distri_all.pdf')

    for lis in choose_items(items, d):
        #
        page_wh = ()
        fig = plt.figure(figsize = page_wh)
        #
        plot_model_results_map(\
                                d, lis,
                                panel_pos=(0.10,0.15,0.75,0.8), bgcolor='none',
                                hcmode=False, hc_list=None, hc_colors=None, useSecond=False,
                                lgd_alpha=0.5, lgd_bbox=None, showlgd=True,
                                draw_rectangles = True, fig_prev=fig, ax_prev=None,
                                ax_cbar=None, cbar_orien=u'vertical', cbar_tick_loc='auto',
                                xRange = xr, yRange = yr,
                                xscale='linear', yscale='linear',
                                xtitle='r (AU)', ytitle='z (AU)',
                                showColorbar=True, rasterized=True,
                                figsize = None, display_dic=None, show_name=True, name_position=(0.5,0.98),
                                labelfontsize=25, tickfontsize=20, textfontsize=25,
                                majorgridon=True, minorgridon=True,
                                xmingrid=False, xmajgrid=False, ymingrid=False, ymajgrid=False,
                                cbarlabel='', useoldstylecbar=False,
                                linthreshy=1e-3, graygrid=False,
                                all_in_one_fig=False, pansepxfrac=0.2, pansepyfrac=0.09,
                                edge_color=None,
                                edge_width=None,
                                draw_ray_to=None,
                                draw_ray_slope_range=None,
                                draw_ray_r_range=None,
                                draw_box_only=False,
                              )
        pdf_pages.savefig(fig)
    pdf_pages.close(fig)



def task6(now_dir = '/n/Users/fdu/now/',
          model_name = '20150307_dep_s',
          iter_fname = 'iter_0001.dat',
          xr = None, yr = None,
          C_O_ratio_r_thr=None,
          C_O_ratio_zr_thr=None,
          freq_rough = None,
          ):

    savefig_dir = opj(now_dir, 'figures/')
    res_dir = opj(now_dir, 'res/')
    model_dir = opj(res_dir, model_name)

    data_fname = opj(model_dir, iter_fname)

    d = load_data_as_dic(data_fname)
    molecule_names = get_molecule_names_from_file(data_fname)
    all_hydrocarbons = get_hydrocarbons(molecule_names)
    all_nitrogens = get_nitrogens(molecule_names)
    hc_abundances = [(hc, np.nanmax(d[hc])) for hc in all_hydrocarbons]
    ni_abundances = [(ni, np.nanmax(d[ni])) for ni in all_nitrogens]
    hc_abundances = sorted(hc_abundances, key = lambda k: k[1])
    ni_abundances = sorted(ni_abundances, key = lambda k: k[1])
    for item in hc_abundances:
        if item[1] > 1e-10:
            print '{0:12s}:  {1:.4e}'.format(item[0], item[1])
    for item in ni_abundances:
        if item[1] > 1e-10:
            print '{0:12s}:  {1:.4e}'.format(item[0], item[1])


def get_major_species_for_element(ele, d, fname, all_elements):
    from numpy import zeros
    molecule_names = get_molecule_names_from_file(fname)
    n = len(d['H'])
    m = len(molecule_names)
    X_ele_s = zeros((n,m))
    for imol, mol in enumerate(molecule_names):
        c = get_ele_counts(ele, mol, all_elements)
        if c != 0:
            try:
                X_ele_s[:, imol]= c * d[mol]
            except:
                raise KeyError(mol)
    imaxmol = np.argmax(X_ele_s, axis=1)
    ele_max = []
    for i in range(n):
        imol = imaxmol[i]
        this_area = (d['rmax'][i] - d['rmin'][i]) * (d['zmax'][i] - d['zmin'][i])
        ele_max.append((molecule_names[imol], X_ele_s[i, imol], this_area))
    return ele_max



def get_all_ions(all_keys):
    return [k for k in all_keys if k.endswith('+')]


def get_major_ions(d_dic, list_all_ions):
    nlen = len(d_dic[d_dic.keys()[0]])
    ion_max = []
    for i in range(nlen):
        list_ion_abun = [d_dic[ion][i] for ion in list_all_ions]
        i_max = np.argmax(list_ion_abun)
        this_area = (d_dic['rmax'][i] - d_dic['rmin'][i]) * (d_dic['zmax'][i] - d_dic['zmin'][i])
        ion_max.append((list_all_ions[i_max], list_ion_abun[i_max], this_area))
    return ion_max

def counts_ions(ion_abun):
    from collections import Counter
    ions = [i_[0] for i_ in ion_abun]
    return Counter(ions)

def counts_ions_weighted(ion_abun):
    from collections import Counter
    ions = {ion: 0.0 for ion in set(i_[0] for i_ in ion_abun)}
    for i_ in ion_abun:
        ion = i_[0]
        ions[ion] += i_[2]
    return sorted([(ion, ions[ion]) for ion in ions.keys()], key=lambda x: -x[1])


def make_colors_for_ion(ion_abun):
    from matplotlib import colors as CLS
    colors_ = ['red', 'green', 'blue', 'CadetBlue', 'magenta', 'BurlyWood',
           'Brown', 'DarkOliveGreen', 'DarkSlateBlue', 'DodgerBlue', 'BlueViolet', 'Chocolate',
           'Tomato', 'Olive', 'SteelBlue', 'SlateBlue', 'Orchid', 'Orange', 'black']
    colors_ = [CLS.colorConverter.to_rgb(c) for c in colors_]
    def color_transf(c):
        def do_1(v):
            if v < 0.5:
                return (v + 0.05) * 1.8
            else:
                return (v - 0.1) / 1.8
        def do_2(v):
            if v < 0.4:
                return (v + 0.15) * 1.8
            else:
                return (v - 0.2) / 1.2
        return (do_1(c[0]), do_2(c[1]), do_2(c[2]))
    colors_ = colors_ + [color_transf(c) for c in colors_]
    n = len(ion_abun)
    # print n, len(colors_)
    ion_cmap = {ion_abun[i][0]: colors_[i] for i in range(n)}
    return ion_cmap


def draw_rect_color_legend(ax, c_dict, ions, startx=0.05, starty=None,
                           height=None, width=0.15, frac_sep=0.3, nlim=None):
    n = len(ions)
    if nlim == None:
        nlim = n
    else:
        nlim = min(n, nlim)
    if height == None:
        height = min(0.04, 0.7 / (frac_sep + 1.0) / nlim)
    if starty == None:
        starty = 0.95- height
    y = starty
    for i in range(nlim):
        ion = ions[i][0]
        rect = Rectangle((startx, starty), width, height, facecolor=c_dict[ion],
                  edgecolor=c_dict[ion], transform=ax.transAxes)
        #print ion, (startx, starty), width, height
        ax.add_patch(rect)
        ax.text(startx, starty, ion, transform=ax.transAxes, fontsize=8)
        starty -= height * (1.0 + frac_sep)
    return


if __name__ == '__main__':
    #
    task4a(savefig_dir='/n/Users/fdu/', xr=(0,400), yr=(0,150))
 #  task3()
    #
#   for mname in [
#                   '20150609_C2H_5a_4a_2d_8cl_B6l_41',
#                   '20150609_C2H_5a_4a_2d_8cl_B6l_41_flatten',
#                ]:
#       task1(now_dir = '/n/Users/fdu/now/',
#           model_name = mname,
#           iter_fname = 'iter_0001.dat',
#           xr=(0,100), yr=(0,80), justRadial=True,
#           )
#   import sys
#   sys.exit()
    #
#   for mname, freq_rough in \
#       [\
#           #('20150609_DMT_C2H_5x_4_5a_3a_rerun24', 262),
#           #('20150609_DMT_C2H_5x_4_5a_3a_rerun24_lowX', 262),
#           #('20150609_DMT_C2H_5x_4_5a_3a_rerun24_highX', 262),
#           #('20150609_C2H_5a_4a_2d_8cl_B6l_41', 115),
#           #('20150609_C2H_5a_4a_2d_8cl_B6l_41_b', 115),
#           #('20160714_TWH_d_rerun', 115),
#           #('20150609_C2H_5a_4a_2d_8cl_B6l_41_rerun', 115),
#           ('20160714_TWH_d_s_noDep_moreOuterMass', 115),
#           ('20160714_TWH_d_s_noDep_uniformD2G', 115),
#           ('20160714_TWH_d_s', 115),
#           #('20160714_TWH_g', 115),
#           #('20160714_TWH_h', 115),
#           #('20160714_TWH_j', 115),
#           #('20150609_C2H_5a_4a_2d_8cl_B6l_41_longT_lessCosmic', 115),
#           #('20150609_C2H_5a_4a_2d_8cl_B0', 349),
#           #('20160714_TWH_a', 115),
#       ]:
#       task1a(now_dir = '/n/Users/fdu/now/',
#       #task6(now_dir = '/n/Users/fdu/now/',
#             model_name = mname,
#             iter_fname = 'iter_0001.dat',
#             xr=(0,150), yr=(0,100),
#             #xr=(0,30), yr=(0,20),
#             #xr=(0,500), yr=(0,300),
#             #xr=(0,80), yr=(0,50),
#             #xr=(0,1), yr=(0,1),
#             #xr=(0,6), yr=(0,4),
#             C_O_ratio_r_thr=70,
#             C_O_ratio_zr_thr=0.8,
#             freq_rough=freq_rough,
#             )
    #task4(now_dir = '/n/Users/fdu/now/',
    #      #model_name = '20150307_noD_d',
    #      model_name = '20150609_C2H_5a',
    #      iter_fname = 'iter_0001.dat',
    #      xr=(0,200), yr=(0,100),
    #      )
    #task1(now_dir = '/n/Users/fdu/now/',
    #      model_name = '20150307_dep_s',
    #      iter_fname = 'iter_0001.dat',
    #      xr=(0,200), yr=(0, 50),
    #      )
    #task1(now_dir = '/n/Users/fdu/now/',
    #      model_name = '20150307_dep_s_NH3',
    #      iter_fname = 'iter_0001.dat',
    #      xr=(0,200), yr=(0, 50),
    #      )
#   task2(now_dir = '/n/Users/fdu/now/',
#         model_name = '20150609_DMT_C2H_5x_4_5a_3a',
#         iter_fname = 'iter_0001.dat',
#         xr = (0,500), yr = (0,300),
#         )
    #task2(now_dir = '/n/Users/fdu/now/',
    #      model_name = '20150609_DMT_C2H_1',
    #      iter_fname = 'iter_0001.dat',
    #      xr = (0,800), yr = (0,200),
    #      )
    #task2(now_dir = '/n/Users/fdu/now/',
    #      model_name = '20150609_DMT_C2H_reac',
    #      iter_fname = 'iter_0001.dat',
    #      xr = (0,800), yr = (0,200),
    #      )
    #task2(now_dir = '/n/Users/fdu/now/',
    #      model_name = '20150609_C2H_2',
    #      iter_fname = 'iter_0001.dat',
    #      xr = (0,200), yr = (0,100),
    #      )
    #task2(now_dir = '/n/Users/fdu/now/',
    #      model_name = '20150307_noD_d_1comp',
    #      iter_fname = 'iter_0001.dat',
    #      xr = (0,200), yr = (0,100),
    #      )
    #task2(now_dir = '/n/Users/fdu/now/',
    #      model_name = '20150307_noD_d',
    #      iter_fname = 'iter_0001.dat',
    #      xr = (0,200), yr = (0,100),
    #      )
    #for modname in \
    #        ['20150307_noD_d_highC2H_1E4yr',
    #         '20150307_noD_d_highC2H_5E4yr',
    #         '20150307_noD_d_highC2H_1E5yr',
    #         '20150307_noD_d_highC2H_2E5yr',
    #         '20150307_noD_d_highC2H_4E5yr',
    #         '20150307_noD_d_highC2H_1E6yr',
    #        ]:
    #    task2(now_dir = '/n/Users/fdu/now/',
    #          model_name = modname,
    #          iter_fname = 'iter_0001.dat',
    #          xr = (0,200), yr = (0,100),
    #          )
