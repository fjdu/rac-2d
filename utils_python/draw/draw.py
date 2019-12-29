import matplotlib as mpl
import numpy as np
import ipywidgets as wdg


def draw_rect(dic, name, ax, xRange=None, yRange=None, cmap=None,
              norm=None, edge_color=None, edge_width=None,
              draw_box_only=False, rasterized=False):
  from  matplotlib.collections import PolyCollection
  import numpy as np
  poly_collec = []
  facecolors = []
  edgecolors = []
  linewidths = []
  nlen = dic['rmin'].size

  maxval = np.nanmax(dic[name])
  minval = np.nanmin(dic[name])

  all_colors = cmap(norm(dic[name]))

  def make_rect(x1, y1, x2, y2):
    return np.array([[x1, y1], [x2, y1], [x2, y2], [x1, y2], [x1, y1]])

  poly_collec = [make_rect(dic['rmin'][i], dic['zmin'][i], dic['rmax'][i], dic['zmax'][i]) for i in range(nlen)]
  if draw_box_only:
    facecolors = ['none'] * nlen
    edgecolors = ['black'] * nlen
    linewidths = [0.1] * nlen
  else:
    facecolors = all_colors
    edgecolors = all_colors if edge_color == None else [edge_color] * nlen
    linewidths = [0.00001] * nlen if edge_width == None else [edge_width] * nlen

  #for i in range(nlen):
  #    # Make rectangle
  #    x1 = dic['rmin'][i]
  #    y1 = dic['zmin'][i]
  #    x2 = dic['rmax'][i]
  #    y2 = dic['zmax'][i]
  #    if x1 > xRange[1] or y1 > yRange[1]:
  #        continue
  #    if x2 < xRange[0] or y2 < yRange[0]:
  #        continue
  #    pxy = np.zeros((5,2))
  #    pxy[0, :] = [x1, y1]
  #    pxy[1, :] = [x2, y1]
  #    pxy[2, :] = [x2, y2]
  #    pxy[3, :] = [x1, y2]
  #    pxy[4, :] = [x1, y1]
  #    #
  #    thiscolor = all_colors[i]
  #    #
  #    # Draw the rectangle filled with color
  #    if draw_box_only:
  #      facecolor = 'none'
  #      edgecolor = 'black'
  #      linewidth = 0.1
  #    else:
  #      facecolor = thiscolor
  #      if edge_color != None or edge_width != None:
  #          edgecolor = edge_color if edge_color != None else 'black'
  #          edgewidth = edge_width if edge_width != None else 0.1
  #      else:
  #          edgecolor = thiscolor
  #          edgewidth = 0.00001
  #    poly_collec.append(pxy)
  #    facecolors.append(facecolor)
  #    edgecolors.append(edgecolor)
  #    linewidths.append(edgewidth)

  ax.add_collection(PolyCollection(poly_collec, edgecolors=edgecolors,
      facecolors=facecolors, linewidths=linewidths, rasterized=rasterized))


def get_tick_values(vmin, vmax, scale='linear'):
    """Return nice-looking tick values."""
    if scale == 'log':
        L = mpl.ticker.LogLocator(base=10, subs=np.linspace(1.0,9.0,9))
    else:
        L = mpl.ticker.AutoLocator()
    L.create_dummy_axis()
    return [_ for _ in L.tick_values(vmin, vmax) if vmin <= _ <= vmax]


def get_color_norm(vmin, vmax, scale='linear', clip=True):
    """Normalize values, with norm(vmin) = 0 and norm(vmax) = 1."""
    if scale == 'linear':
        norm = mpl.colors.Normalize
    else:
        norm = mpl.colors.LogNorm
    return norm(vmin=vmin, vmax=vmax, clip=clip)


def add_colorbar(ax=None, ax_cbar=None, cmap=None, wfrac=0.03, norm=None):
    if not ax_cbar:
        fig = ax.get_figure()
        pos = ax.get_position()
        cb_width = (pos.x1 - pos.x0) * wfrac
        ax_cbar = fig.add_axes([pos.x1 + cb_width, pos.y0, cb_width, pos.y1-pos.y0])

    cbar = mpl.colorbar.ColorbarBase(ax_cbar, cmap=cmap,
                                 orientation='vertical', norm=norm)



def get_idx_from_xy(d, x, y):
    idx = np.where(np.logical_and.reduce(
        (d['rmin'] <= x, d['rmax'] >= x,
         d['zmin'] <= y, d['zmax'] >= y)))
    if idx[0].size != 0:
        return idx[0][0]
    return None


def onclick_pixval(d, name=None, d_axes=None):
    def onclick(event):
        inaxes = event.inaxes
        if d_axes:
            name = d_axes[inaxes]
        x, y = event.xdata, event.ydata
        idx = get_idx_from_xy(d, x, y)
        if idx != None:
            val = d[name][idx]
            txt.value = '(x,y,{0:})=({1:.4g}, {2:.4g}, {3:.4g})'.format(name, x, y, val)
        else:
            txt.value = 'NULL'
    
    txt = wdg.Textarea(
        value='',
        layout={'width': '100%', 'height': '40px'},
        placeholder='',
        description='Val:',
        disabled=True)
    display(txt)
    
    return onclick



def get_range(vec, logrange=1e10):
    vmin, vmax = np.nanmin(vec), np.nanmax(vec)
    if logrange > 0:
      if vmin*logrange < vmax:
          vmin = vmax / logrange
    return vmin, vmax



def draw_one_species(ax, d, d_axes, name, xRange, yRange, cmap,
                     logrange=1e6, scale='log', vmin=None, vmax=None, hidextick=False, hideytick=False):
    if (vmin is not None) and (vmax is not None):
      vRange = (vmin, vmax)
    elif scale == 'log':
      vRange = get_range(d[name], logrange=logrange)
    else:
      vRange = get_range(d[name], logrange=0)
    norm = get_color_norm(*vRange, scale=scale, clip=True)

    draw_rect(d, name, ax, xRange=xRange, yRange=yRange, cmap=cmap, norm=norm)
    add_colorbar(ax=ax, cmap=cmap, norm=norm)
    ax.text(0.08, 0.87, name, transform=ax.transAxes)
    if hidextick:
        ax.set_xticklabels([])
    if hideytick:
        ax.set_yticklabels([])
    d_axes[ax] = name
    return


def get_nx_ny(nitems):
    npanx = np.ceil(np.sqrt(nitems))
    npany = np.ceil(float(nitems)/npanx)
    return npanx, npany


def draw_multi_species(fig, d, d_axes, items=None, xRange=None, yRange=None, cmap=None,
                       xtitle='x (au)', ytitle='y (au)',
                       xscale='linear', yscale='linear',
                       scale='log', logrange=1e8,
                       vmin=None, vmax=None,
                       xmarginleft = 0.1, ymarginlower = 0.15,
                       pansepxfrac=0.3, pansepyfrac=0.09):    
    nitems = len(items)
    npanx, npany = get_nx_ny(nitems)
    
    panwx = (0.99 - xmarginleft)  / npanx
    panwy = (0.99 - ymarginlower) / npany
    pansepx = panwx * pansepxfrac
    pansepy = panwy * pansepyfrac
    panwx -= pansepx
    panwy -= pansepy
    
    for ii in range(nitems):
        print(items[ii]['name'])
        kx = np.mod(ii, npanx)
        ky = npany - np.ceil((ii+1)/npanx)
        xleft  = xmarginleft  + panwx * kx + pansepx * kx
        ylower = ymarginlower + panwy * ky + pansepy * ky
        pos = (xleft, ylower, panwx, panwy)
        ax = fig.add_axes(pos,
              xlabel=xtitle if ky == 0 else '',
              ylabel=ytitle if kx == 0 else '',
              autoscalex_on=False, autoscaley_on=False,
              xscale=items[ii].get('xscale') or xscale,
              yscale=items[ii].get('yscale') or yscale,
              xlim=items[ii].get('xRange') or xRange,
              ylim=items[ii].get('yRange') or yRange)
        draw_one_species(ax, d, d_axes, items[ii]['name'],
                         items[ii].get('xRange') or xRange,
                         items[ii].get('yRange') or yRange,
                         cmap, logrange=items[ii].get('logrange') or logrange,
                         scale=items[ii].get('scale') or scale,
                         scale=items[ii].get('vmin') or vmin,
                         scale=items[ii].get('vmax') or vmax,
                         hidextick=False if ky == 0 else True,
                         hideytick=False if kx == 0 else True)
    return
