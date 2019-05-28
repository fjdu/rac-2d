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
