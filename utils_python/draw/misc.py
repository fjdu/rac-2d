def load_data_as_dic(filepath, comments='!', returnOriginalKeys=False):
    import numpy as np
    
    data = np.loadtxt(filepath, comments=comments)

    ftmp = open(filepath, 'r')
    str_comment = ftmp.readline()[1:].split()
    ftmp.close()

    dic = {}
    for i in range(len(str_comment)):
        dic.update({str_comment[i]: data[:, i]})

    del data

    if returnOriginalKeys:
        return str_comment, dic
    else:
        return dic



def to_spherical(d, r_grid, theta_grid, phi_grid, items):
    """Export the model data for selected items using spherical coordinates.

    d: model data in the form of a dictionary; obtained by calling load_data_as_dic.

    r_grid: radial grid boundary points; unit: AU

    theta_grid: theta grid boundary points; unit: radian

    phi_grid: phi grid boundary points; unit: radian

    items: a dictionary of items to be exported in the following format; note
           that the filename is the output filename (to be used as input for RADMC 3D)
           to save data for each item:
        items = {
                    'rhodus_1': {'filename': 'rhodus_1.inp'},
                    'rhodus_2': {'filename': 'rhodus_2.inp'},
                    'rhodus_3': {'filename': 'rhodus_3.inp'},
                    'Tdust1': {'filename': 'Tdust1.inp'},
                    'Tdust2': {'filename': 'Tdust2.inp'},
                    'Tdust3': {'filename': 'Tdust3.inp'},
                    'C2H':  {'filename': 'C2H.inp'},
                    'C3H2': {'filename': 'C3H2.inp'}
                }
    """
    import numpy as np
    nr, nt, nphi = len(r_grid), len(theta_grid), len(phi_grid)
    files = {}

    for key in items:
        files.update({key: open(items[key]['filename'], 'w')})

    state = query_state()

    for i in range(nphi-1):
        phi = 0.5 * (phi_grid[i] + phi_grid[i+1])
        for j in range(nt-1):
            theta = 0.5 * (theta_grid[j] + theta_grid[j+1])
            for k in range(nr-1):
                r = 0.5 * (r_grid[k] + r_grid[k+1])
                rho = r * np.sin(theta)
                z = r * np.cos(theta)
                for key in items:
                    val = state.query(d, rho, z, key)
                    files[key].write('{0:.6e}\n'.format(val))

    for key in items:
        files[key].close()



def num2str4fig(num, scientific=True):
  import re
  s1 = '{0:16.1f}'.format(num).strip()
  s2 = [re.sub(r'^([\-+]*)(0+)([1-9]+)', r'\1\3', item) for item in \
        [item.replace('+','') for item in '{0:16.1e}'.format(num).strip().split('e')]]
  len1 = len(s1)
  len2 = 1
  for item in s2:
    len2 += len(item)
  choose = 1
  if (num != 0.0):
    if float(s1) == 0.0:
      choose = 2
  int_0 = int(float(s1))
  if abs(int_0 - float(s1)) < 0.01:
    s1 = '{:d}'.format(int_0)
  if (not scientific) or (len1 <= len2 and choose==1):
    return r'${' + s1 + r'}$'
  else:
    if float(s2[1]) == 0.0 or float(s2[1]) == 1.0:
      return r'${' + s1 + r'}$'
    else:
      int_ = int(float(s2[0]))
      if abs(float(s2[0])-int_) < 0.01:
        if int_ == 1:
          return r'$10^{' + s2[1] + r'}$'
        else:
          return r'${' + '{:d}'.format(int_) + r'}{\times}10^{' + s2[1] + r'}$'
      else:
        return r'${' + s2[0] + r'}{\times}10^{' + s2[1] + r'}$'



def clip_fraction(num_str, th=0.01):
  num = float(num_str)
  num_round = round(num)
  if abs(num_round-num)/num < th:
    return num_str.split('.')[0]
  else:
    return num_str


all_elements = ['H', 'D', 'He', 'C', 'N', 'O', 'Si', 'S', 'Fe', 'Na', 'Mg', 'Cl', 'P', 'F', 'Ne', 'Ar', 'K']


def get_all_matches(el, word):
    m = []
    len_el = len(el)
    len_w = len(word)
    i = 0
    while(1):
        j = i + len_el
        if j > len_w:
            break
        if word[i:j] == el:
            m.append(i)
        i += 1
    return m

def get_leading_digits(word):
    s = ''
    for L in word:
        if '0' <= L and L <= '9':
            s += L
        else:
            break
    return s
    
def get_ele_counts(el, word, all_el):
    len_el = len(el)
    len_w = len(word)
    el_c = [c for c in all_el if (c.startswith(el) and c!=el)]
    idx = get_all_matches(el, word)
    ct = 0
    for i in idx:
        if len([_ for _ in filter(word[i:].startswith, el_c)]) != 0:
            continue
        ii = i+len_el
        if i+len_el == len_w:
            ct += 1
            continue
        Ldig = get_leading_digits(word[ii:])
        if Ldig == '':
            ct += 1
            continue
        else:
            ct += int(Ldig)
    return ct

def get_molecule_names_from_file(fname, col_start=None):
    #col_start = 133
    with open(fname, 'r') as f:
        str_comment = f.readline()[1:].split()
    if col_start == None:
        col_start = str_comment.index('H')
    return str_comment[col_start:]

def update_elemental_abundance(ele, d, fname, all_elements):
    from numpy import zeros
    molecule_names = get_molecule_names_from_file(fname)
    n = len(d['H'])
    X_ele = zeros(n)
    for k in molecule_names:
        c = get_ele_counts(ele, k, all_elements)
        if c != 0:
            try:
                X_ele += c * d[k]
            except:
                raise KeyError(k)
    d['X['+ele+']'] = X_ele
    return

def update_density_of_X(X, d):
    d['n('+X+')'] = d[X] * d['n_gas']
    return


def update_stuff(ca, filename):
    update_density_of_X('H2O', ca)
    update_density_of_X('OH', ca)
    update_density_of_X('CO', ca)
    update_density_of_X('O', ca)
    update_elemental_abundance('O', ca, filename, all_elements)
    update_elemental_abundance('C', ca, filename, all_elements)
    ca['Tgas-Tdust'] = ca['Tgas'] - ca['Tdust']


def is_hydrocarbon(name):
    c_H = get_ele_counts('H', name, all_elements)
    c_C = get_ele_counts('C', name, all_elements)
    if c_H == 0 or c_C == 0:
        return False
    for el in all_elements:
        if el == 'H' or el == 'C':
            continue
        c_el = get_ele_counts(el, name, all_elements)
        if c_el > 0:
            return False
    return True


def has_nitrogen(name):
    c_N = get_ele_counts('N', name, all_elements)
    if c_N > 0:
        return True
    else:   
        return False


def get_hydrocarbons(all_names):
    hc = []
    for name in all_names:
        if is_hydrocarbon(name):
            hc.append(name)
    return hc


def get_nitrogens(all_names):
    ni = []
    for name in all_names:
        if has_nitrogen(name):
            ni.append(name)
    return ni
