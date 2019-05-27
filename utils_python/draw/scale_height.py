from matplotlib import *
use('Agg')
from matplotlib.pyplot import *
from numpy import *

from os.path import join as opj
from glob import glob

from long_function_definitions import *
from my_script import *

import consts.physical_constants as consts

model_name = '20150307_dep_p1'
iter_fname = 'iter_0001.dat'

now_dir = '/n/Users/fdu/now/'
#now_dir = '/u/Moria2/fdu/now/'
savefig_dir = opj(now_dir, 'figures/')
res_dir = opj(now_dir, 'res/')
model_dir = opj(res_dir, model_name)

data_fname = opj(model_dir, iter_fname)


def update_stokes_num(d, a0_grain_CGS, beta=-0.3, rho_grain_CGS=2.0):
    a_grain_CGS = a0_grain_CGS * ((d['rmin']+d['rmax'])*0.5)**beta
    d['StokesNum'] = d['w_Kep'] * a_grain_CGS * rho_grain_CGS / \
                     (d['c_sound'] * d['n_gas'] * consts.phy_mProton_CGS)

def update_scale_height(d, alpha_scaling=1.0):
    d['ScaleHeight'] = sqrt(alpha_scaling * d['alpha'] / (minimum(d['StokesNum'], 0.5) * (1.0 + d['StokesNum'])))


def draw_radial(dic, name_list, zpos, savefig_dir=None):
    plot_model_results_radial(dic, name_list, zpos, xRange = None, figsize = None,
        xscale='linear', yscale='log', returndata=True, noLegend=False, noTitle=False,
        color=None, label=None, display_dic=None,
        ax=None, returnaxis=False, linestyle='-',lw=3,markersize=5,
        do_interpol=False, interpolNum=200,
        dosmooth=False, winwidth=None,
        whichx='rmin')
    fig_fname = model_name.replace('/', '') + '_' + make_fname_from_namelist(name_list) + '_radial.pdf'
    if savefig_dir == None:
        savefig_dir = opj(now_dir, 'figures/')
    savefig(opj(savefig_dir, fig_fname), bbox_inches='tight')
    print 'File saved at ', opj(savefig_dir, fig_fname)
    return

mycm = make_my_colormap(c_list=[(0.0, (0.5, 0.0, 0.5)),
                                (0.2, (0.0, 0.0, 1.0)),
                                (0.4, (0.0, 0.8, 1.0)),
                                (0.6, (0.0, 0.8, 0.0)),
                                (0.8, (1.0, 0.8, 0.0)),
                                (1.0, (1.0, 0.0, 0.0))])

rcParams['axes.color_cycle'] = mycolors

d = load_data_as_dic(data_fname)
update_stuff(d, data_fname)
update_stokes_num(d, a0_grain_CGS=0.01, rho_grain_CGS=2.0)
update_scale_height(d, alpha_scaling=100.0)

draw_radial(d, [[{'name': 'ScaleHeight', 'vr': (1e-4, 1e2)}, {'name': 'StokesNum'}]], 0.0)
raise SystemExit

draw_radial(d, [[{'name': 'c_sound', 'vr': (1e4, 1e14)}, {'name': 'n_gas'}]], 0.0)


d['n(SO)'] = d['SO'] * d['n_gas']
d['n(CS)'] = d['CS'] * d['n_gas']
d['SO/H2O'] = d['SO'] / d['H2O']
calc_colden(d, 'SO')
calc_colden(d, 'CS')
d['N(SO)/N(CS)'] = d['N(SO)'] / d['N(CS)']

draw_name_list(d,
    name_list = \
        [
         [{'name': 'n_gas', 'linearscale': False, 'cmap': mycm, 'vr': (1e4, 1e14),
          'name_disp': r'$n_{\rm gas}$ (cm$^{-3}$)'}],
         [{'name': 'Tgas', 'linearscale': False, 'cmap': mycm, 'vr': (1e1, 2e3),
          'name_disp': r'$T_{\rm gas}$ (K)'}],
         [{'name': 'Tdust', 'linearscale': False, 'cmap': mycm, 'vr': (1e1, 2e3),
          'name_disp': r'$T_{\rm dust}$ (K)'}],
         [{'name': 'SO', 'linearscale': False, 'cmap': mycm, 'vr': (1e-16, 1e-8),
          'name_disp': r'SO'}],
         [{'name': 'CS', 'linearscale': False, 'cmap': mycm, 'vr': (1e-16, 1e-8),
          'name_disp': r'CS'}],
         [{'name': 'n(SO)', 'linearscale': False, 'cmap': mycm, 'vr': (1e-4, 1e6),
          'name_disp': r'n(SO)'}],
         [{'name': 'n(CS)', 'linearscale': False, 'cmap': mycm, 'vr': (1e-4, 1e6),
          'name_disp': r'n(CS)'}],
         [{'name': 'SO/H2O', 'linearscale': False, 'cmap': mycm, 'vr': (1e-6, 1e2),
          'name_disp': r'SO/H2O'}],
         [{'name': 'OH', 'linearscale': False, 'cmap': mycm, 'vr': (1e-12, 1e-7),
          'name_disp': r'OH'}],
         [{'name': 'H2O', 'linearscale': False, 'cmap': mycm, 'vr': (1e-12, 1e-7),
          'name_disp': r'H2O'}],
         [{'name': 'O', 'linearscale': False, 'cmap': mycm, 'vr': (1e-14, 1e-4),
          'name_disp': r'O'}],
         [{'name': 'C', 'linearscale': False, 'cmap': mycm, 'vr': (1e-14, 1e-4),
          'name_disp': r'C'}],
         [{'name': 'CO', 'linearscale': False, 'cmap': mycm, 'vr': (1e-14, 1e-4),
          'name_disp': r'CO'}],
         [{'name': 'NO', 'linearscale': False, 'cmap': mycm, 'vr': (1e-12, 1e-7),
          'name_disp': r'NO'}],
         [{'name': 'CN', 'linearscale': False, 'cmap': mycm, 'vr': (1e-12, 1e-7),
          'name_disp': r'CN'}],
         [{'name': 'G0_UV', 'linearscale': False, 'cmap': mycm, 'vr': (1e-3, 1e6),
          'name_disp': r'$G_{{\rm 0,UV}}$'}],
        ])

raise SystemExit

d['Tgas-Tdust'] = d['Tgas'] - d['Tdust']

draw_name_list(d,
    name_list = \
        [[{'name': 'n_gas', 'linearscale': False, 'cmap': mycm, 'vr': (1e3, 1e13),
          'name_disp': r'$n_{\rm gas}$ (cm$^{-3}$)'}],
         [{'name': 'Tgas', 'linearscale': False, 'cmap': mycm, 'vr': (1e1, 5e2),
          'name_disp': r'$T_{\rm gas}$ (K)'}],
         [{'name': 'Tdust', 'linearscale': False, 'cmap': mycm, 'vr': (1e1, 5e2),
          'name_disp': r'$T_{\rm dust}$ (K)'}],
         [{'name': 'Tgas-Tdust', 'linearscale': True, 'cmap': mycm, 'vr': (0, 20.01),
          'name_disp': r'$T_{\rm gas} - T_{\rm dust}$'}],
        ])

draw_name_list(d,
    name_list = \
        [[{'name': 'n_gas', 'linearscale': False, 'cmap': mycm, 'vr': (1e3, 1e13),
          'name_disp': r'$n_{\rm gas}$ (cm$^{-3}$)'}],
         [{'name': 'X[O]', 'linearscale': False, 'cmap': mycm, 'vr': (1e-8, 3e-4),
          'name_disp': r'$X[{\rm O}]$'}],
         [{'name': 'Tgas', 'linearscale': False, 'cmap': mycm, 'vr': (1e1, 3e3),
          'name_disp': r'$T_{\rm gas}$ (K)'}],
         [{'name': 'Tdust', 'linearscale': False, 'cmap': mycm, 'vr': (1e1, 3e3),
          'name_disp': r'$T_{\rm dust}$ (K)'}],
         [{'name': 'X[C]', 'linearscale': False, 'cmap': mycm, 'vr': (1e-8, 3e-4),
          'name_disp': r'$X[{\rm C}]$'}],
         [{'name': 'N_CO_I', 'linearscale': False, 'cmap': mycm, 'vr': (1e8, 1e19),
          'name_disp': r'$N_{\rm CO}$'}],
         [{'name': 'N_H2O_I', 'linearscale': False, 'cmap': mycm, 'vr': (1e8, 1e19),
          'name_disp': r'$N_{\rm H2O}$'}],
         [{'name': 'G0_UV', 'linearscale': False, 'cmap': mycm, 'vr': (1e-3, 1e7),
          'name_disp': r'$G_{{\rm 0,UV}}$'}],
         [{'name': 'Av_ISM', 'linearscale': False, 'cmap': mycm, 'vr': (1e-2, 1e5),
          'name_disp': r'Av'}],
        ])


draw_name_list(d,
    name_list = \
        [[{'name': 'G0_UV', 'linearscale': False, 'cmap': mycm, 'vr': (1e-3, 1e7),
          'name_disp': r'$G_{{\rm 0,UV}}$'}],
        ])

calc_colden(d, 'C2H')
calc_colden(d, 'C2H3')
draw_name_list(d,
    name_list = \
        [[{'name': 'C2H', 'linearscale': False, 'cmap': mycm, 'vr': (1e-15, 1e-7),
          'name_disp': r'C2H'}],
         [{'name': 'C2H3', 'linearscale': False, 'cmap': mycm, 'vr': (1e-15, 1e-7),
          'name_disp': r'C2H3'}],
         [{'name': 'N(C2H)', 'linearscale': False, 'cmap': mycm, 'vr': (1e7, 1e17),
          'name_disp': r'N(C2H)'}],
         [{'name': 'N(C2H3)', 'linearscale': False, 'cmap': mycm, 'vr': (1e7, 1e17),
          'name_disp': r'N(C2H3)'}],
        ])

draw_name_list(d,
    name_list = \
        [[{'name': 'NH2+', 'linearscale': False, 'cmap': mycm, 'vr': (1e-15, 1e-7),
          'name_disp': r'NH2+'}],
         [{'name': 'HF', 'linearscale': False, 'cmap': mycm, 'vr': (1e-15, 1e-7),
          'name_disp': r'HF'}],
         [{'name': 'OH+', 'linearscale': False, 'cmap': mycm, 'vr': (1e-15, 1e-7),
          'name_disp': r'OH+'}],
         [{'name': 'H2O+', 'linearscale': False, 'cmap': mycm, 'vr': (1e-15, 1e-7),
          'name_disp': r'H2O+'}],
         [{'name': 'OH', 'linearscale': False, 'cmap': mycm, 'vr': (1e-12, 1e-6),
          'name_disp': r'OH'}],
         [{'name': 'H2O', 'linearscale': False, 'cmap': mycm, 'vr': (1e-12, 1e-6),
          'name_disp': r'H2O'}],
         [{'name': 'O', 'linearscale': False, 'cmap': mycm, 'vr': (1e-14, 1e-4),
          'name_disp': r'O'}],
         [{'name': 'G0_UV', 'linearscale': False, 'cmap': mycm, 'vr': (1e-3, 1e7),
          'name_disp': r'$G_{{\rm 0,UV}}$'}],
         [{'name': 'Av_ISM', 'linearscale': False, 'cmap': mycm, 'vr': (1e-2, 1e5),
          'name_disp': r'Av'}],
        ])
