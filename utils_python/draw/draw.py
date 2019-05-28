def draw_rect(dic, name, ax, xRange=None, yRange=None, colormap=None,
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

  for i in range(nlen):
      # Make rectangle
      x1 = dic['rmin'][i]
      y1 = dic['zmin'][i]
      x2 = dic['rmax'][i]
      y2 = dic['zmax'][i]
      if x1 > xRange[1] or y1 > yRange[1]:
          continue
      if x2 < xRange[0] or y2 < yRange[0]:
          continue
      pxy = np.zeros((5,2))
      pxy[0, :] = [x1, y1]
      pxy[1, :] = [x2, y1]
      pxy[2, :] = [x2, y2]
      pxy[3, :] = [x1, y2]
      pxy[4, :] = [x1, y1]
      #
      # Calculate the color
      thiscolor = colormap(norm(dic[name][i]))
      #
      # Draw the rectangle filled with color
      if draw_box_only:
        facecolor = 'none'
        edgecolor = 'black'
        linewidth = 0.1
      else:
        facecolor = thiscolor
        if edge_color != None or edge_width != None:
            edgecolor = edge_color if edge_color != None else 'black'
            edgewidth = edge_width if edge_width != None else 0.1
        else:
            edgecolor = thiscolor
            edgewidth = 0.00001
      poly_collec.append(pxy)
      facecolors.append(facecolor)
      edgecolors.append(edgecolor)
      linewidths.append(edgewidth)

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


def add_colorbar(ax, minval, maxval, ax_cbar=None, wfrac=0.03, norm=None):
    if not ax_cbar:
        fig = ax.get_figure()
        pos = ax.get_position()
        cb_width = (pos.x1 - pos.x0) * wfrac
        ax_cbar = fig.add_axes([pos.x1 + cb_width, pos.y0, cb_width, pos.y1-pos.y0])

    cbar = colorbar.ColorbarBase(ax_cbar, cmap=cm.rainbow,
                                 orientation='vertical', norm=norm)
